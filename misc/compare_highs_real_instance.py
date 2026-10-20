"""Compare the HiGHS and Gurobi backends of buildProblem() on a real instance.

Builds one buildProblem() instance from a target list and a few pointings
the way pfs_target_uploader's PPP does (same Bench, classdict, cobraMoveCost
and solver options), once with solver="gurobi" and once with solver="highs",
and compares the two models by content: column bounds and integrality, rows
as (lower, upper, coefficients keyed by column name) after normalising the
sign convention, and the objective as accumulated on prob.cost. With --solve
it also compares objective value, status, solve time and the set of Tv_Cv_*
arcs that ended up at 1.

buildProblem() draws from numpy's global random state (RandomTargetSelector
in _get_vis_and_elbow), so the order in which it creates variables changes
from one call to the next unless the state is reset; np.random.seed(--seed)
is called before each build. The comparison itself is by name and does not
depend on that order.

Usage:
  python misc/compare_highs_real_instance.py INPUT_DIR [--nvisit 4] [--solve]
        [--gap 0.0] [--time-limit 600] [--out DIR] [--max-targets N] [--seed 20]
INPUT_DIR holds target_<id>.ecsv (ob_code, ra, dec, exptime, priority) and
ppc_<id>.ecsv (ppc_ra, ppc_dec, ppc_pa). Needs gurobipy, highspy, cobraOps,
cobraCharmer and pfs.instdata importable; the spt_target_uploader venv has
them all. gurobipy's pip wheel comes with a size-limited licence (2000
variables / constraints), enough for --max-targets 300 --nvisit 1; larger
instances need a full Gurobi licence.

With --gap 0 both solvers return a proven optimum and the assignments should
coincide whenever the optimum is unique; with a positive gap, only the
objective values are expected to agree, to within twice the gap.
"""
import argparse
import glob
import logging
import os
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(HERE)
sys.path.insert(0, REPO_ROOT)

import ets_fiber_assigner.netflow as nf  # noqa: E402

BIG = 1e30   # anything beyond this is treated as infinite (Gurobi uses 1e100)


def make_bench():
    from pfs.instdata import setup_envvar
    from ics.cobraOps.Bench import Bench
    from ics.cobraCharmer.cobraCoach.cobraCoach import CobraCoach
    setup_envvar()
    with tempfile.TemporaryDirectory() as d:
        cc = CobraCoach(loadModel=True, trajectoryMode=True, rootDir=d)
        return Bench(cobraCoach=cc, blackDotsMargin=1.65)


def classdict_like_ppp():
    # NetflowPreparation() in pfs_target_uploader/utils/ppp.py
    return {f"sci_P{p}": {"nonObservationCost": 100 - 10 * p,
                          "partialObservationCost": 200, "calib": False}
            for p in range(10)}


def observation_time(ra, dec):
    try:
        from pfs_target_uploader.utils.ppp import set_observation_time
        return set_observation_time(ra, dec=dec)
    except Exception as e:  # noqa: BLE001
        print(f"set_observation_time unavailable ({e}); using a fixed time")
        return "2026-10-01T10:00:00Z"


def load_inputs(input_dir, nvisit, max_targets):
    from astropy.table import Table
    tfile = glob.glob(os.path.join(input_dir, "target_*.ecsv"))[0]
    pfile = glob.glob(os.path.join(input_dir, "ppc_*.ecsv"))[0]
    tab = Table.read(tfile)
    if max_targets and len(tab) > max_targets:
        tab = tab[:max_targets]
    ppc = Table.read(pfile)
    seen, tel = set(), []
    for row in ppc:
        key = (float(row["ppc_ra"]), float(row["ppc_dec"]), float(row["ppc_pa"]))
        if key in seen:
            continue
        seen.add(key)
        tel.append(key)
        if len(tel) == nvisit:
            break
    tgt = [nf.ScienceTarget(r["ob_code"], r["ra"], r["dec"], r["exptime"],
                            r["priority"], "sci") for r in tab]
    return tgt, tel


# ---------------------------------------------------------------------------
# canonical model: backend-independent description by column name
# ---------------------------------------------------------------------------
def _inf(x):
    x = float(x)
    if x >= BIG:
        return float("inf")
    if x <= -BIG:
        return float("-inf")
    return x


def _canonical_row(lo, hi, terms):
    """terms: dict name -> coef. Fix the sign so the alphabetically first
    column has a positive coefficient; flip and swap the bounds if not."""
    terms = dict(terms)
    if terms:
        first = min(terms)
        if terms[first] < 0:
            terms = {k: -v for k, v in terms.items()}
            lo, hi = -hi, -lo
    return (round(_inf(lo), 9), round(_inf(hi), 9),
            tuple(sorted((k, round(v, 9)) for k, v in terms.items())))


def canonical_highs(prob):
    import highspy
    prob.update()
    prob._passNames()
    lp = prob._prob.getLp()
    cnames = list(lp.col_names_)
    integ = list(lp.integrality_) if len(lp.integrality_) else [None] * lp.num_col_
    cols = {n: (_inf(lo), _inf(hi), it == highspy.HighsVarType.kInteger)
            for n, lo, hi, it in zip(cnames, lp.col_lower_, lp.col_upper_, integ)}
    start = np.asarray(lp.a_matrix_.start_)
    index = np.asarray(lp.a_matrix_.index_)
    value = np.asarray(lp.a_matrix_.value_)
    entries = {r: {} for r in range(lp.num_row_)}
    if lp.a_matrix_.format_ == highspy.MatrixFormat.kRowwise:
        for r in range(lp.num_row_):
            for k in range(start[r], start[r + 1]):
                entries[r][cnames[index[k]]] = float(value[k])
    else:
        for c in range(lp.num_col_):
            for k in range(start[c], start[c + 1]):
                entries[index[k]][cnames[c]] = float(value[k])
    rows = sorted(_canonical_row(lo, hi, entries[r])
                  for r, (lo, hi) in enumerate(zip(lp.row_lower_, lp.row_upper_)))
    obj = {}
    cost = prob.cost
    if hasattr(cost, "idxs"):
        for i, v in zip(cost.idxs, cost.vals):
            obj[cnames[i]] = obj.get(cnames[i], 0.0) + float(v)
        const = cost.constant or 0.0
    else:
        obj[cnames[cost.index]] = 1.0
        const = 0.0
    obj = {k: round(v, 9) for k, v in obj.items()}
    return cols, rows, obj, float(const)


def canonical_gurobi(prob):
    import gurobipy as gbp
    m = prob._prob
    m.update()
    cols = {v.VarName: (_inf(v.LB), _inf(v.UB), v.VType in (gbp.GRB.BINARY, gbp.GRB.INTEGER))
            for v in m.getVars()}
    rows = []
    for c in m.getConstrs():
        row = m.getRow(c)
        terms = {}
        for i in range(row.size()):
            n = row.getVar(i).VarName
            terms[n] = terms.get(n, 0.0) + float(row.getCoeff(i))
        rhs = float(c.RHS)
        if c.Sense == gbp.GRB.LESS_EQUAL:
            lo, hi = float("-inf"), rhs
        elif c.Sense == gbp.GRB.GREATER_EQUAL:
            lo, hi = rhs, float("inf")
        else:
            lo, hi = rhs, rhs
        rows.append(_canonical_row(lo, hi, terms))
    rows.sort()
    obj = {}
    cost = prob.cost
    if isinstance(cost, gbp.Var):
        obj[cost.VarName] = 1.0
        const = 0.0
    else:
        for i in range(cost.size()):
            n = cost.getVar(i).VarName
            obj[n] = obj.get(n, 0.0) + float(cost.getCoeff(i))
        const = float(cost.getConstant())
    obj = {k: round(v, 9) for k, v in obj.items()}
    return cols, rows, obj, const


def compare_canonical(a, b):
    """Return a list of human-readable differences (empty when identical)."""
    diffs = []
    for what, da, db in (("columns", a[0], b[0]), ("objective", a[2], b[2])):
        if set(da) != set(db):
            diffs.append(f"{what}: name sets differ ({len(set(da) ^ set(db))} names)")
            continue
        bad = [k for k in da if da[k] != db[k]]
        if bad:
            diffs.append(f"{what}: {len(bad)} entries differ, e.g. {bad[0]}: {da[bad[0]]} vs {db[bad[0]]}")
    if a[1] != b[1]:
        sa, sb = set(a[1]), set(b[1])
        ex = next(iter(sa ^ sb), None)
        diffs.append(f"rows: {len(a[1])} vs {len(b[1])} rows, {len(sa ^ sb)} differ, e.g. {ex}")
    if a[3] != b[3]:
        diffs.append(f"objective constant {a[3]} vs {b[3]}")
    return diffs


# ---------------------------------------------------------------------------
def build(solver, bench, tgt, tpos, classdict, nvisit, gap, time_limit):
    if solver == "gurobi":
        opts = {"MIPGap": gap, "Seed": 0, "OutputFlag": 0, "TimeLimit": float(time_limit)}
    else:
        opts = {"mip_rel_gap": gap, "random_seed": 0, "output_flag": False,
                "time_limit": float(time_limit)}
    t0 = time.perf_counter()
    prob = nf.buildProblem(
        bench, tgt, tpos, classdict, 900.0, [0] * nvisit,
        cobraMoveCost=lambda d: 0.1 * d,
        collision_distance=2.0, elbow_collisions=True,
        solver=solver, solverOptions=opts,
        alreadyObserved={}, forbiddenPairs=[[] for _ in range(nvisit)],
        avoidFiducials=False, brokenCobrasMargin=0.0)
    return prob, time.perf_counter() - t0


def sizes(prob):
    if isinstance(prob, nf.GurobiProblem):
        m = prob._prob
        m.update()
        return m.NumVars, m.NumConstrs, m.NumNZs
    prob.update()
    h = prob._prob
    return h.getNumCol(), h.getNumRow(), h.getNumNz()


def status_and_objective(prob):
    if isinstance(prob, nf.GurobiProblem):
        import gurobipy as gbp
        names = {getattr(gbp.GRB, n): n for n in ("OPTIMAL", "TIME_LIMIT", "INFEASIBLE",
                                                   "UNBOUNDED", "INTERRUPTED", "SUBOPTIMAL")}
        return names.get(prob._prob.Status, str(prob._prob.Status)), prob._prob.ObjVal
    h = prob._prob
    return h.modelStatusToString(h.getModelStatus()), h.getObjectiveValue()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input_dir")
    ap.add_argument("--nvisit", type=int, default=4)
    ap.add_argument("--solve", action="store_true")
    ap.add_argument("--gap", type=float, default=0.0,
                    help="MIP relative gap for both solvers (default 0: proven optimum)")
    ap.add_argument("--time-limit", type=float, default=600.0)
    ap.add_argument("--out", default=None)
    ap.add_argument("--max-targets", type=int, default=0)
    ap.add_argument("--seed", type=int, default=20)
    a = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    for noisy in ("cobraCoach", "butler", "root"):
        logging.getLogger(noisy).setLevel(logging.WARNING)

    out = a.out or tempfile.mkdtemp(prefix="highs_vs_gurobi_")
    os.makedirs(out, exist_ok=True)

    tgt, tel = load_inputs(a.input_dir, a.nvisit, a.max_targets)
    otime = observation_time(tel[0][0], tel[0][1])
    telescopes = [nf.Telescope(ra, dec, pa, otime) for ra, dec, pa in tel]
    print(f"{len(tgt)} targets, {len(telescopes)} pointings, otime {otime}")
    bench = make_bench()
    tpos = [t.get_fp_positions(tgt) for t in telescopes]
    classdict = classdict_like_ppp()

    results = {}
    for solver in ("gurobi", "highs"):
        np.random.seed(a.seed)   # buildProblem's RandomTargetSelector uses np.random
        prob, t_build = build(solver, bench, tgt, tpos, classdict, len(telescopes),
                              a.gap, a.time_limit)
        r = {"t_build": t_build}
        r["ncols"], r["nrows"], r["nnz"] = sizes(prob)
        r["canon"] = canonical_gurobi(prob) if solver == "gurobi" else canonical_highs(prob)
        prob.dump(os.path.join(out, f"{solver}.lp"))
        if a.solve:
            t0 = time.perf_counter()
            prob.solve()
            r["t_solve"] = time.perf_counter() - t0
            r["status"], r["objective"] = status_and_objective(prob)
            r["assigned"] = {k for k, v in prob._vardict.items()
                             if k.startswith("Tv_Cv_") and prob.value(v) > 0.5}
        results[solver] = r
        print(f"[{solver:6s}] build {t_build:.2f}s  cols {r['ncols']} rows {r['nrows']} nnz {r['nnz']}"
              + (f"  solve {r['t_solve']:.1f}s {r['status']} obj {r['objective']:.6f} "
                 f"assigned {len(r['assigned'])}" if a.solve else ""))
        del prob

    G, H = results["gurobi"], results["highs"]
    ok = (G["ncols"], G["nrows"], G["nnz"]) == (H["ncols"], H["nrows"], H["nnz"])
    print(f"sizes identical: {ok}")
    diffs = compare_canonical(G["canon"], H["canon"])
    print(f"model identical by content (columns, rows, objective): {not diffs}")
    for d in diffs:
        print("   ", d)
    ok &= not diffs
    if a.solve:
        og, oh = G["objective"], H["objective"]
        tol = (2.0 * a.gap + 1e-9) * max(1.0, abs(og))
        same_obj = abs(og - oh) <= tol
        same_asg = G["assigned"] == H["assigned"]
        print(f"objective agree within tolerance {tol:.3g}: {same_obj} (|diff| {abs(og - oh):.3g})  "
              f"status gurobi {G['status']} / highs {H['status']}  "
              f"assignment identical: {same_asg} (sym. diff {len(G['assigned'] ^ H['assigned'])})")
        ok &= same_obj
        if a.gap == 0.0 and not same_asg:
            print("    note: gap 0 but assignments differ -> the optimum is degenerate "
                  "(check objective agreement above)")
    print(f"build time gurobi {G['t_build']:.2f}s / highs {H['t_build']:.2f}s"
          + (f"; solve gurobi {G['t_solve']:.1f}s / highs {H['t_solve']:.1f}s" if a.solve else "")
          + f"; files in {out}")
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
