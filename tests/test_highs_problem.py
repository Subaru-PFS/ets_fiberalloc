"""Tests for HighsProblem, the HiGHS backend in ets_fiber_assigner/netflow.py.

Three kinds of checks:

1. HighsProblem's own behaviour: when value() may be read, what dump()
   writes before and after solve(), names, row flushes, error cases.
2. Model identity: HighsProblem buffers rows and hands them to HiGHS in one
   addRows call. The model that produces must be the one highspy's own
   one-row-at-a-time API (Highs.addVariable / Highs.addConstr) builds from
   the same expressions -- matrix arrays and the written .lp/.mps files are
   compared byte for byte.
3. Solutions: small models with a known optimum, and agreement with
   GurobiProblem on the same models (skipped when gurobipy is not
   importable; the pip wheel's size-limited licence covers these models).

Run with:  python -m pytest tests -v
"""
import filecmp
import os
import sys

import numpy as np
import pytest

highspy = pytest.importorskip("highspy")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from ets_fiber_assigner.netflow import HighsProblem, GurobiProblem  # noqa: E402

try:
    import gurobipy  # noqa: F401
    HAVE_GUROBI = True
except ImportError:
    HAVE_GUROBI = False

needs_gurobi = pytest.mark.skipif(not HAVE_GUROBI, reason="gurobipy not installed")


# ---------------------------------------------------------------------------
# Reference: the same LPProblem-style interface on top of highspy's own
# convenience API, one column and one row at a time. This is the path
# HighsProblem replaced with bulk addCols/addRows, so the two must build
# identical models.
# ---------------------------------------------------------------------------
class ReferenceHighs(object):
    def __init__(self):
        self._prob = highspy.Highs()
        self._prob.setOptionValue("output_flag", False)
        self.cost = self._prob.addVariable(0.0, highspy.kHighsInf, name="cost")
        self.sum = self._prob.qsum
        self._vardict = {}
        self._constraintdict = {}
        self._bounds = {}

    def addVar(self, name, lo, hi):
        lo = -highspy.kHighsInf if lo is None else lo
        hi = highspy.kHighsInf if hi is None else hi
        var = self._prob.addIntegral(lo, hi, name=name)
        self._vardict[name] = var
        self._bounds[var.index] = (float(lo), float(hi))
        return var

    def add_constraint(self, name, constraint):
        self._constraintdict[name] = constraint
        self._prob.addConstr(constraint, name=name)

    add_lazy_constraint = add_constraint

    def update(self):
        pass

    def dump(self, filename):
        self._prob.writeModel(filename)

    def solve(self):
        self._prob.minimize(self.cost)

    def value(self, var):
        return self._prob.val(var)

    def varBounds(self, var):
        return self._bounds[var.index]

    def changeVarBounds(self, var, lower=None, upper=None):
        lb, ub = self._bounds[var.index]
        lb = lb if lower is None else lower
        ub = ub if upper is None else upper
        self._bounds[var.index] = (float(lb), float(ub))
        self._prob.changeColBounds(var.index, lb, ub)


def make_gurobi():
    return GurobiProblem(extraOptions={"OutputFlag": 0, "MIPGap": 0.0})


def make_highs():
    return HighsProblem(extraOptions={"mip_rel_gap": 0.0})


# ---------------------------------------------------------------------------
# Build sequences, written against the LPProblem interface so the same code
# drives HighsProblem, ReferenceHighs and GurobiProblem. Each returns the
# handles it created.
# ---------------------------------------------------------------------------
def scen_mixed(prob):
    x = prob.addVar("x", 0, 1)
    y = prob.addVar("y", 0, 1)
    z = prob.addVar("z", 0, None)
    w = prob.addVar("w", None, None)
    u = prob.addVar("u", 2, 5)
    a = [prob.addVar(f"a_{i}", 0, 1) for i in range(4)]
    b = [prob.addVar(f"b_{i}", 0, 1) for i in range(3)]
    for i, v in enumerate([x, y, z, w, u] + a + b):
        prob.cost += v * float(i + 1)
    prob.cost += w * 3.0                       # same column twice in the objective
    prob.add_constraint("le", prob.sum([x, y]) <= 1)
    prob.add_constraint("ge", prob.sum([x, y, z]) >= 1)
    prob.add_constraint("eq", prob.sum(a) == 2)
    prob.add_constraint("lhs_const", x + 2 <= 5)
    prob.add_constraint("both_sides", a[0] >= 0.25 * prob.sum(b))
    prob.add_constraint("both_sides2", prob.sum(a) <= prob.sum(b) + 1)
    prob.add_constraint("dup", prob.sum([a[2], a[3], a[2]]) <= 1)
    prob.add_constraint("cancel", x - x + y == 0)
    prob.add_constraint("empty", prob.sum([]) <= 3)
    prob.add_constraint("free_lo", w >= -7)
    prob.add_constraint("free_hi", w - z <= 4)
    prob.add_constraint("neg", prob.sum([a[1]] + [-v for v in b]) == 0)
    prob.add_lazy_constraint("lazy", prob.sum([a[2], b[2]]) <= 1)
    prob.add_constraint("scaled",
                        prob.sum([v * t for v, t in zip(b, [900.0, 450.0, 1800.0])])
                        >= 900.0 * y)
    prob.add_constraint("eq_var", x == 1)
    return dict(x=x, y=y, z=z, w=w, u=u, a=a, b=b)


def scen_no_constraints(prob):
    vs = [prob.addVar(f"v_{i}", 0, 1) for i in range(3)]
    for i, v in enumerate(vs):
        prob.cost += v * float(i + 1)
    return dict(vs=vs)


def scen_single_var(prob):
    x = prob.addVar("x", 0, 1)
    prob.cost += x * 2.0
    prob.add_constraint("c", x >= 1)
    return dict(x=x)


def scen_objective_untouched(prob):
    # cost stays the bare variable handed out in __init__
    x = prob.addVar("x", 0, 1)
    prob.add_constraint("c", x <= 1)
    return dict(x=x)


def scen_bounds_before_solve(prob):
    h = scen_mixed(prob)
    # changeVarBounds in the middle of the build: the reference applies it to
    # HiGHS at once, HighsProblem creates the still-pending column with the
    # new bounds. Both must end up with the same model.
    prob.changeVarBounds(h["u"], lower=3)
    prob.add_constraint("late", prob.sum([h["a"][3], h["b"][0]]) <= 1)
    prob.changeVarBounds(h["b"][1], upper=0)
    return h


def scen_assignment_small(prob):
    """Three targets, two cobras, hand-checkable optimum.

    Arc costs: t0-c0 1, t0-c1 5, t1-c0 2, t2-c1 3; leaving a target
    unobserved (its sink arc) costs 10. Best: t0 on c0, t2 on c1, t1 unobserved
    -> 1 + 3 + 10 = 14, and no other assignment reaches 14.
    """
    arcs = {}
    for (t, c, cost) in [(0, 0, 1.0), (0, 1, 5.0), (1, 0, 2.0), (2, 1, 3.0)]:
        f = prob.addVar(f"Tv_Cv_{t}_{c}", 0, 1)
        prob.cost += f * cost
        arcs[(t, c)] = f
    sinks = []
    for t in range(3):
        s = prob.addVar(f"ST_sink_{t}", 0, 1)
        prob.cost += s * 10.0
        sinks.append(s)
    for t in range(3):
        flows = [f for (tt, c), f in arcs.items() if tt == t] + [sinks[t]]
        prob.add_constraint(f"TvIO_{t}", prob.sum(flows) == 1)
    for c in range(2):
        flows = [f for (t, cc), f in arcs.items() if cc == c]
        prob.add_constraint(f"Cvlim_{c}", prob.sum(flows) <= 1)
    return dict(arcs=arcs, sinks=sinks, optimum=14.0,
                assigned={"Tv_Cv_0_0", "Tv_Cv_2_1"})


def scen_assignment_random(prob, ntargets=60, ncobras=25, seed=5):
    """Flow-shaped model like buildProblem's: targets, cobras, collision pairs.

    Random continuous costs make the optimum unique in practice, so two exact
    solvers must agree on the assignment, not only on the objective.
    """
    rng = np.random.default_rng(seed)
    arcs, by_cobra, by_target, all_arcs = {}, [[] for _ in range(ncobras)], [], []
    for t in range(ntargets):
        cobras = rng.choice(ncobras, size=int(rng.integers(1, 4)), replace=False)
        tarcs = []
        for c in cobras:
            f = prob.addVar(f"Tv_Cv_{t}_{c}", 0, 1)
            prob.cost += f * float(rng.random())
            arcs[(t, int(c))] = f
            by_cobra[c].append(f)
            tarcs.append(f)
            all_arcs.append(f)
        s = prob.addVar(f"ST_sink_{t}", 0, 1)
        prob.cost += s * 10.0
        tarcs.append(s)
        by_target.append(tarcs)
    for p in range(ntargets // 2):
        i, j = rng.integers(0, len(all_arcs), size=2)
        if i != j:
            prob.add_lazy_constraint(f"Coll_{p}", prob.sum([all_arcs[i], all_arcs[j]]) <= 1)
    for c in range(ncobras):
        if by_cobra[c]:
            prob.add_constraint(f"Cvlim_{c}", prob.sum(by_cobra[c]) <= 1)
    for t in range(ntargets):
        prob.add_constraint(f"TvIO_{t}", prob.sum(by_target[t]) == 1)
    return dict(arcs=arcs)


def scen_random_rows(prob, n=300, seed=7):
    """Random rows of every shape, duplicates included; exercises the sort."""
    rng = np.random.default_rng(seed)
    vs = [prob.addVar(f"v_{i}", 0, 1) for i in range(n)]
    for v in vs:
        prob.cost += v * float(rng.random())
    for k in range(n):
        m = int(rng.integers(1, 6))
        idx = rng.integers(0, n, size=m)            # duplicates allowed
        expr = prob.sum([vs[i] * float(rng.integers(-2, 3)) for i in idx])
        kind = k % 3
        if kind == 0:
            prob.add_constraint(f"r_{k}", expr <= float(rng.integers(0, 3)))
        elif kind == 1:
            prob.add_constraint(f"r_{k}", expr >= float(-rng.integers(0, 3)))
        else:
            prob.add_lazy_constraint(f"r_{k}", expr <= prob.sum(vs[:2]) + 1)
    return dict(vs=vs)


SCENARIOS = {
    "mixed": scen_mixed,
    "no_constraints": scen_no_constraints,
    "single_var": scen_single_var,
    "objective_untouched": scen_objective_untouched,
    "bounds_before_solve": scen_bounds_before_solve,
    "assignment_small": scen_assignment_small,
    "assignment_random": scen_assignment_random,
    "random_rows": scen_random_rows,
}


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def assert_same_highs_model(ref, new):
    """Element-by-element comparison of the two HiGHS-side models."""
    a, b = ref._prob, new._prob
    assert (a.getNumCol(), a.getNumRow(), a.getNumNz()) == \
        (b.getNumCol(), b.getNumRow(), b.getNumNz())
    la, lb = a.getLp(), b.getLp()
    for attr in ("col_cost_", "col_lower_", "col_upper_", "row_lower_", "row_upper_"):
        assert np.array_equal(np.asarray(getattr(la, attr)),
                              np.asarray(getattr(lb, attr))), attr
    assert la.a_matrix_.format_ == lb.a_matrix_.format_
    for attr in ("start_", "index_", "value_"):
        assert np.array_equal(np.asarray(getattr(la.a_matrix_, attr)),
                              np.asarray(getattr(lb.a_matrix_, attr))), attr
    assert list(la.integrality_) == list(lb.integrality_)
    assert list(la.col_names_) == list(lb.col_names_)
    assert list(la.row_names_) == list(lb.row_names_)


def assert_same_dump(ref, new, tmp_path, tag):
    for ext in ("lp", "mps"):
        fa, fb = tmp_path / f"{tag}_ref.{ext}", tmp_path / f"{tag}_new.{ext}"
        ref.dump(str(fa))
        new.dump(str(fb))
        assert filecmp.cmp(fa, fb, shallow=False), f".{ext} differs: {fa} {fb}"


def objective(prob):
    if isinstance(prob, GurobiProblem):
        return prob._prob.ObjVal
    return prob._prob.getObjectiveValue()


def assigned(prob):
    return {k for k, v in prob._vardict.items()
            if k.startswith("Tv_Cv_") and prob.value(v) > 0.5}


# ---------------------------------------------------------------------------
# 2. model identity against highspy's own addConstr path
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name", sorted(SCENARIOS))
def test_model_matches_highspy_addconstr_path(name, tmp_path):
    ref, new = ReferenceHighs(), HighsProblem()
    SCENARIOS[name](ref)
    SCENARIOS[name](new)
    assert_same_dump(ref, new, tmp_path, name)      # before solve: no objective
    assert_same_highs_model(ref, new)               # names reached HiGHS in dump()
    ref.solve()
    new.solve()
    assert_same_dump(ref, new, tmp_path, name + "_solved")
    assert (ref._prob.modelStatusToString(ref._prob.getModelStatus())
            == new._prob.modelStatusToString(new._prob.getModelStatus()))
    assert abs(objective(ref) - objective(new)) <= 1e-9 * max(1.0, abs(objective(new)))


def test_dump_twice_and_incremental(tmp_path):
    ref, new = ReferenceHighs(), HighsProblem()
    hr, hn = scen_mixed(ref), scen_mixed(new)
    assert_same_dump(ref, new, tmp_path, "d1")
    assert_same_dump(ref, new, tmp_path, "d2")
    assert filecmp.cmp(tmp_path / "d1_new.lp", tmp_path / "d2_new.lp", shallow=False)
    # add more columns and rows after a dump, then dump again
    for prob, h in ((ref, hr), (new, hn)):
        q = prob.addVar("q", 0, 3)
        prob.cost += q * 0.5
        prob.add_constraint("after_dump", prob.sum([q, h["x"]]) >= 2)
    assert_same_dump(ref, new, tmp_path, "d3")
    new.update()
    assert_same_highs_model(ref, new)
    # and the incremental build equals a one-shot build of the same model
    fresh = HighsProblem()
    h = scen_mixed(fresh)
    q = fresh.addVar("q", 0, 3)
    fresh.cost += q * 0.5
    fresh.add_constraint("after_dump", fresh.sum([q, h["x"]]) >= 2)
    fresh.dump(str(tmp_path / "fresh.lp"))
    assert filecmp.cmp(tmp_path / "d3_new.lp", tmp_path / "fresh.lp", shallow=False)


def test_change_bounds_after_solve(tmp_path):
    ref, new = ReferenceHighs(), HighsProblem()
    hr, hn = scen_mixed(ref), scen_mixed(new)
    ref.solve()
    new.solve()
    obj1 = objective(new)
    for prob, h in ((ref, hr), (new, hn)):
        prob.changeVarBounds(h["u"], lower=4)
        prob.changeVarBounds(h["a"][3], upper=0)
        prob.changeVarBounds(h["w"], lower=-2, upper=10)
        prob.solve()
    assert_same_dump(ref, new, tmp_path, "cb")
    assert abs(objective(ref) - objective(new)) <= 1e-9 * max(1.0, abs(objective(new)))
    assert objective(new) != obj1
    for x, y in zip(ref._vardict.values(), new._vardict.values()):
        assert ref.varBounds(x) == new.varBounds(y)


# ---------------------------------------------------------------------------
# 3. solutions: known optimum, and agreement with Gurobi
# ---------------------------------------------------------------------------
def test_known_optimum():
    new = make_highs()
    h = scen_assignment_small(new)
    new.solve()
    assert abs(objective(new) - h["optimum"]) < 1e-9
    assert assigned(new) == h["assigned"]
    assert new.value(h["sinks"][1]) == 1.0
    assert new.value(h["sinks"][0]) == 0.0 and new.value(h["sinks"][2]) == 0.0


SOLVABLE = ["mixed", "no_constraints", "single_var", "objective_untouched",
            "bounds_before_solve", "assignment_small", "assignment_random"]


@needs_gurobi
@pytest.mark.parametrize("name", SOLVABLE)
def test_agrees_with_gurobi(name):
    g, h = make_gurobi(), make_highs()
    SCENARIOS[name](g)
    SCENARIOS[name](h)
    g.solve()
    h.solve()
    og, oh = objective(g), objective(h)
    assert abs(og - oh) <= 1e-9 * max(1.0, abs(og)), (og, oh)
    assert list(g._vardict) == list(h._vardict)
    for a, b in zip(g._vardict.values(), h._vardict.values()):
        assert g.varBounds(a) == h.varBounds(b)
    if name.startswith("assignment"):
        # random continuous costs: the optimum is unique, so the solvers must
        # pick the same arcs, not just reach the same value
        assert assigned(g) == assigned(h)


@needs_gurobi
def test_agrees_with_gurobi_after_bound_changes():
    g, h = make_gurobi(), make_highs()
    scen_assignment_random(g)
    scen_assignment_random(h)
    g.solve()
    h.solve()
    # forbid two arcs of the first solution and re-solve
    forbid = sorted(assigned(h))[:2]
    for prob in (g, h):
        for name in forbid:
            prob.changeVarBounds(prob.varByName(name), upper=0)
        prob.solve()
    assert abs(objective(g) - objective(h)) <= 1e-9 * max(1.0, abs(objective(g)))
    assert assigned(g) == assigned(h)
    assert not (assigned(h) & set(forbid))


# ---------------------------------------------------------------------------
# 1. HighsProblem behaviour
# ---------------------------------------------------------------------------
def test_infeasible_raises():
    new = HighsProblem()
    x = new.addVar("x", 0, 1)
    y = new.addVar("y", 0, 1)
    new.cost += x * 1.0 + y * 1.0
    new.add_constraint("c", new.sum([x, y]) >= 3)
    with pytest.raises(RuntimeError):
        new.solve()
    with pytest.raises(RuntimeError):       # no all-zero "solution" either
        new.value(x)


def test_unbounded_expression_rejected():
    new = HighsProblem()
    x = new.addVar("x", 0, 1)
    y = new.addVar("y", 0, 1)
    with pytest.raises(Exception):
        new.add_constraint("bad", new.sum([x, y]))
    # nothing half-added: the model still builds and solves
    new.add_constraint("ok", new.sum([x, y]) <= 1)
    new.solve()
    assert new._prob.getNumRow() == 1
    assert "bad" not in new._constraintdict


def test_names_reach_highs_in_order(tmp_path):
    new = HighsProblem()
    scen_mixed(new)
    new.solve()          # so the objective is in the file as well
    new.dump(str(tmp_path / "n.lp"))
    lp = new._prob.getLp()
    assert list(lp.col_names_) == ["cost"] + list(new._vardict)
    assert list(lp.row_names_) == list(new._constraintdict)
    text = (tmp_path / "n.lp").read_text()
    for name in ("both_sides2", "eq_var", "a_3"):
        assert name in text
    obj_line = text.split("obj:", 1)[1].split("\n", 1)[0]
    assert "x" in obj_line and "w" in obj_line, obj_line


def test_dump_objective_only_after_solve(tmp_path):
    """The objective reaches HiGHS in solve(), so a dump written before
    solve() has an empty objective -- exactly what GurobiProblem.dump()
    writes before solve() ("Minimize 0 cost").
    """
    new = HighsProblem()
    x = new.addVar("x", 0, 1)
    y = new.addVar("y", 0, 1)
    new.cost += x * 2.0 + y * 3.0
    new.cost += x * 0.5
    new.add_constraint("c", new.sum([x, y]) >= 1)
    new.dump(str(tmp_path / "pre.lp"))
    pre = (tmp_path / "pre.lp").read_text()
    obj_line = pre.split("obj:", 1)[1].split("\n", 1)[0]
    assert obj_line.strip() == "", obj_line
    new.solve()
    new.dump(str(tmp_path / "post.lp"))
    post = (tmp_path / "post.lp").read_text()
    # "+1 cost": the objective is accumulated onto the free `cost` column, so
    # that column itself carries coefficient 1 -- same as with Gurobi.
    assert "obj: +1 cost +2.5 x +3 y" in post, post
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.readModel(str(tmp_path / "post.lp"))
    h.run()
    assert abs(h.getObjectiveValue() - objective(new)) < 1e-9


def test_row_flush_count():
    new = HighsProblem()
    scen_mixed(new)
    assert new._nrow_flushes == 0
    assert new._prob.getNumRow() == 0          # rows still buffered
    new.solve()
    assert new._nrow_flushes == 1
    assert new._prob.getNumRow() == len(new._constraintdict)
    new.solve()                                # nothing pending: no new flush
    new.update()
    assert new._nrow_flushes == 1


def test_lookup_by_name():
    ref, new = ReferenceHighs(), HighsProblem()
    scen_mixed(ref)
    scen_mixed(new)
    for name in ("x", "w", "a_2"):
        assert new.varByName(name) is new._vardict[name]
        assert new.varByName(name).index == ref._vardict[name].index
    for name in ("dup", "cancel", "empty"):
        ca, cb = ref._constraintdict[name], new.constraintByName(name)
        assert (ca.idxs, ca.vals, ca.bounds) == (cb.idxs, cb.vals, cb.bounds)


def test_value_without_solution_raises():
    """value() refuses whenever the cached solution is invalid.

    Invalidation follows Gurobi's timing: a change that has reached the
    solver (changeVarBounds, or a flush through update()/dump()/solve())
    discards the solution; a change that is still buffered leaves the old
    solution readable, just as a pending Gurobi modification leaves var.X.
    """
    new = HighsProblem()
    x = new.addVar("x", 0, 1)
    y = new.addVar("y", 0, 1)
    new.cost += x * 1.0 + y * 2.0
    new.add_constraint("c", new.sum([x, y]) >= 1)
    with pytest.raises(RuntimeError):        # nothing flushed yet
        new.value(x)
    new.update()
    with pytest.raises(RuntimeError):        # flushed, not solved
        new.value(x)
    new.solve()
    assert new.value(x) == 1.0 and new.value(y) == 0.0
    new.changeVarBounds(x, upper=0)          # pending: old solution readable
    assert new.value(x) == 1.0
    new.update()                             # applied: stale
    with pytest.raises(RuntimeError):
        new.value(x)
    new.solve()
    assert new.value(x) == 0.0 and new.value(y) == 1.0
    new.add_constraint("late", new.sum([x, y]) <= 5)
    assert new.value(y) == 1.0               # row still buffered: readable
    new.update()                             # row reached HiGHS: stale
    with pytest.raises(RuntimeError):
        new.value(y)
    new.solve()
    assert new.value(y) == 1.0
    z = new.addVar("z", 0, 1)                # column still buffered
    new.cost += z * 1.0
    assert new.value(y) == 1.0               # old column: readable
    with pytest.raises(RuntimeError):        # new column: no value yet
        new.value(z)
    new.update()                             # column reached HiGHS: stale
    with pytest.raises(RuntimeError):
        new.value(y)
    new.solve()
    assert new.value(z) == 0.0


@needs_gurobi
def test_value_rule_matches_gurobi():
    """Same sequence on both backends: value() is readable or raises together."""
    def probe(prob, var):
        try:
            return ("value", float(prob.value(var)))
        except Exception:            # GurobiError / RuntimeError
            return ("raises",)

    steps = []
    for prob in (make_gurobi(), make_highs()):
        x = prob.addVar("x", 0, 1)
        y = prob.addVar("y", 0, 1)
        prob.cost += x * 1.0 + y * 2.0
        prob.add_constraint("c", prob.sum([x, y]) >= 1)
        prob.update()
        seq = [probe(prob, x)[0]]                       # before solve
        prob.solve()
        seq.append(probe(prob, x))                       # after solve
        prob.changeVarBounds(x, upper=0)
        seq.append(probe(prob, x)[0])                    # pending change
        prob.update()
        seq.append(probe(prob, x)[0])                    # applied change
        prob.solve()
        seq.append(probe(prob, x))                       # re-solved
        prob.add_constraint("late", prob.sum([x, y]) <= 5)
        seq.append(probe(prob, y)[0])                    # pending row
        prob.update()
        seq.append(probe(prob, y)[0])                    # applied row
        steps.append(seq)
    assert steps[0] == steps[1], steps
