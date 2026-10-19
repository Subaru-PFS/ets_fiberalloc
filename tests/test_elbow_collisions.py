import tracemalloc
from collections import defaultdict

import numpy as np

from ics.cobraOps.CollisionSimulator import CollisionSimulator
from ets_fiber_assigner.netflow import _get_elbow_collisions


def _reference(bench, tpos, vis, dist):
    """_get_elbow_collisions as of v3.9, one call per Cobra-target pair."""
    tpos = np.array(tpos)
    ivis = defaultdict(list)
    epos = defaultdict(list)
    for tidx, cidx_elbow in vis.items():
        for (cidx, elbowpos) in cidx_elbow:
            ivis[cidx].append(tidx)
            epos[cidx].append((tidx, elbowpos))
    res = defaultdict(list)
    for cidx, thing in epos.items():
        nb = bench.getCobraNeighbors(cidx)
        tmp = [ivis[j] for j in nb if j in ivis]
        if tmp != []:
            i2 = np.concatenate([ivis[j] for j in nb if j in ivis])
            i2 = np.unique(i2).astype(int)
            for tidx, elbowpos in thing:
                ebp = np.full(len(i2), elbowpos)
                tp = np.full(len(i2), tpos[tidx])
                d = CollisionSimulator.distancesToLineSegments(tpos[i2], tp, ebp)
                res[(cidx, tidx)] += list(i2[d < dist])
    return res


class _LineBench:
    """Cobras on a line, 8 mm apart; neighbours are the adjacent Cobras."""
    def __init__(self, n):
        self.n = n
        self.centers = 8.0 * np.arange(n) + 0j

    def getCobraNeighbors(self, cidx):
        return np.array([j for j in (cidx - 1, cidx + 1) if 0 <= j < self.n])


def _random_case(seed, ncobras=60, ntargets=400):
    rng = np.random.default_rng(seed)
    bench = _LineBench(ncobras)
    tpos = list(8.0 * ncobras * rng.random(ntargets)
                + 6j * (rng.random(ntargets) - 0.5))
    vis = defaultdict(list)
    for tidx, p in enumerate(tpos):
        for cidx in np.flatnonzero(np.abs(bench.centers - p) < 5.0):
            elbow = bench.centers[cidx] + 2.4 * np.exp(1j * rng.uniform(0, 2 * np.pi))
            vis[tidx].append((int(cidx), elbow))
    return bench, tpos, vis


def _as_lists(res):
    return [(k, [int(x) for x in v]) for k, v in res.items()]


def test_matches_reference():
    for seed in range(5):
        bench, tpos, vis = _random_case(seed)
        assert _as_lists(_get_elbow_collisions(bench, tpos, vis, 2.0)) == \
            _as_lists(_reference(bench, tpos, vis, 2.0))


def test_keys_are_python_ints():
    bench, tpos, vis = _random_case(0)
    res = _get_elbow_collisions(bench, tpos, vis, 2.0)
    assert res and all(type(c) is int and type(t) is int for c, t in res)


def test_target_on_elbow_segment_end_and_duplicates():
    bench = _LineBench(2)
    tpos = [1 + 0j, 1 + 0j, 9 + 0j, 3 + 0j]
    vis = {0: [(0, 3 + 0j)], 1: [(0, 3 + 0j)], 2: [(1, 7 + 0j)],
           3: [(0, 1 + 0j), (1, 5 + 0j)]}
    assert _as_lists(_get_elbow_collisions(bench, tpos, vis, 2.0)) == \
        _as_lists(_reference(bench, tpos, vis, 2.0))


def test_crowded_cobras_keep_memory_bounded():
    # 3000 targets around two neighbouring Cobras: one call over all
    # (target, neighbour target) combinations would need ~9e6-element arrays
    rng = np.random.default_rng(1)
    bench = _LineBench(2)
    tpos = list(4.0 + 3.0 * (rng.random(3000) - 0.5)
                + 3j * (rng.random(3000) - 0.5))
    vis = {tidx: [(cidx, bench.centers[cidx] + 2.4 * np.exp(1j * rng.uniform(0, 2 * np.pi)))
                  for cidx in (0, 1)] for tidx in range(len(tpos))}
    # a small collision distance keeps the result itself small, so the peak
    # measures the temporary arrays
    tracemalloc.start()
    try:
        res = _get_elbow_collisions(bench, tpos, vis, 0.01)
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    assert peak < 64 * 2**20
    assert _as_lists(res) == _as_lists(_reference(bench, tpos, vis, 0.01))


def test_no_neighbours_with_targets():
    bench = _LineBench(3)
    tpos = [0j]
    vis = {0: [(0, 2 + 0j)]}
    assert dict(_get_elbow_collisions(bench, tpos, vis, 2.0)) == {}
