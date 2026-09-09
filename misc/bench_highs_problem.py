"""Task 0 benchmark: drive HighsProblem directly with a synthetic flow-like model.

Usage: python misc/bench_highs_problem.py N [--no-solve] [--time-limit S] [--chunk K] [--threads T]

Columns (N total, approximately):
  T targets, C = T/2 cobras.  Each target has 1..5 arcs to cobras (avg 3)
  plus one sink arc  ->  ~4T columns.  A handful of (0, None) overflow columns
  exercise the None-bound path.
Rows (M ~ 0.6 N):
  collision pairs      x_a + x_b <= 1            (T rows, 2 nnz)   added "lazy"
  cobra capacity       sum(arcs into c) <= 1     (C rows, ~6 nnz)
  target conservation  sum(arcs of t) + sink == 1 (T rows, ~4 nnz)
Objective: prob.cost += var * coef for every column.
"""
import argparse
import sys
import time

import numpy as np

import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from ets_fiber_assigner.netflow import HighsProblem  # noqa: E402


def build_and_time(N, cls, seed=1, solve=True, time_limit=600.0, chunk=20000,
                   threads=None, keep_prob=False):
    rng = np.random.default_rng(seed)
    T = max(1, N // 4)
    C = max(1, T // 2)
    opts = {"time_limit": float(time_limit)}
    if threads is not None:
        opts["threads"] = int(threads)
    prob = cls(extraOptions=opts)
    res = {"N_target": N}

    # ---- variables + objective accumulation ------------------------------
    t0 = time.perf_counter()
    arcs_by_cobra = [[] for _ in range(C)]
    arcs_by_target = []
    all_arcs = []
    for t in range(T):
        k = int(rng.integers(1, 6))
        cobras = rng.choice(C, size=min(k, C), replace=False)
        tarcs = []
        for c in cobras:
            f = prob.addVar(f"Tv_Cv_{t}_{c}", 0, 1)
            prob.cost += f * float(rng.random())
            arcs_by_cobra[c].append(f)
            tarcs.append(f)
            all_arcs.append(f)
        s = prob.addVar(f"ST_sink_{t}", 0, 1)
        prob.cost += s * 10.0
        tarcs.append(s)
        arcs_by_target.append(tarcs)
    for j in range(8):
        f = prob.addVar(f"STC_sink_{j}", 0, None)
        prob.cost += f * 1.0
    t1 = time.perf_counter()
    res["t_addVar_loop"] = t1 - t0          # Python-side, no HiGHS call yet
    prob.update()                            # forces _flush (addCols [+names])
    t2 = time.perf_counter()
    res["t_flush"] = t2 - t1
    res["t_addVar_total"] = t2 - t0
    res["ncols"] = prob._prob.getNumCol()

    # ---- constraints -----------------------------------------------------
    # Build all expressions first so expression-building cost is separated
    # from the add_constraint cost.
    t3 = time.perf_counter()
    exprs = []
    narcs = len(all_arcs)
    for p in range(T):
        i, j = rng.integers(0, narcs, size=2)
        if i == j:
            continue
        exprs.append((f"Coll_{p}", prob.sum([all_arcs[i], all_arcs[j]]) <= 1))
    ncoll = len(exprs)
    for c in range(C):
        if arcs_by_cobra[c]:
            exprs.append((f"Cvlim_{c}", prob.sum(arcs_by_cobra[c]) <= 1))
    for t in range(T):
        exprs.append((f"TvIO_{t}", prob.sum(arcs_by_target[t]) == 1))
    t4 = time.perf_counter()
    res["t_expr_build"] = t4 - t3
    res["nrows_planned"] = len(exprs)

    chunk_times = []
    tc = time.perf_counter()
    for k, (name, e) in enumerate(exprs):
        if k < ncoll:
            prob.add_lazy_constraint(name, e)
        else:
            prob.add_constraint(name, e)
        if (k + 1) % chunk == 0:
            now = time.perf_counter()
            chunk_times.append(now - tc)
            tc = now
    t5 = time.perf_counter()
    res["t_add_constraint"] = t5 - t4
    res["chunk_times"] = chunk_times
    # HighsProblem may buffer rows until something needs the full model;
    # update() forces that, so the row flush is timed on its own here.
    prob.update()
    t5b = time.perf_counter()
    res["t_flush_rows"] = t5b - t5
    res["t_constraints_total"] = t5b - t4
    res["nrows"] = prob._prob.getNumRow()
    res["nnz"] = prob._prob.getNumNz()
    res["nrow_flushes"] = getattr(prob, "_nrow_flushes", None)

    if solve:
        t6 = time.perf_counter()
        prob.solve()
        t7 = time.perf_counter()
        res["t_solve"] = t7 - t6
        res["status"] = prob._prob.modelStatusToString(prob._prob.getModelStatus())
        res["objective"] = prob._prob.getObjectiveValue()
        t8 = time.perf_counter()
        vals = [prob.value(v) for v in prob._vardict.values()]
        t9 = time.perf_counter()
        res["t_value_all"] = t9 - t8
        res["sum_values"] = float(sum(vals))
    if keep_prob:
        res["prob"] = prob
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=float)
    ap.add_argument("--no-solve", action="store_true")
    ap.add_argument("--time-limit", type=float, default=600.0)
    ap.add_argument("--chunk", type=int, default=20000)
    ap.add_argument("--threads", type=int, default=None)
    a = ap.parse_args()
    cls = HighsProblem
    res = build_and_time(int(a.N), cls, solve=not a.no_solve,
                         time_limit=a.time_limit, chunk=a.chunk,
                         threads=a.threads)
    print(f"=== {cls.__name__}  N={int(a.N):,} ===")
    for k, v in res.items():
        if k == "chunk_times":
            print(f"  {k:18s} " + " ".join(f"{x:.2f}" for x in v))
        elif isinstance(v, float):
            print(f"  {k:18s} {v:.4f}")
        else:
            print(f"  {k:18s} {v}")


if __name__ == "__main__":
    main()
