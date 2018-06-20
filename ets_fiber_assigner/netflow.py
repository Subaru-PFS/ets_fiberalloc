from __future__ import print_function
import numpy as np
import pycconv
from collections import defaultdict


def _get_colliding_pairs(bench, tpos, vis, dist):
    tpos = np.array(tpos)
    ivis = defaultdict(list)
    for tidx, thing in vis.items():
        for (cidx, _) in thing:
            ivis[cidx].append(tidx)

    pairs = set()
    for cidx, i1 in ivis.items():
        # determine target indices visible by this cobra and its neighbors
        nb = bench.getCobraNeighbors(cidx)
        i2 = np.concatenate([ivis[j] for j in nb if j in ivis])
        i2 = np.concatenate((i1, i2))
        i2 = np.unique(i2).astype(np.int)
        d = np.abs(np.subtract.outer(tpos[i1], tpos[i2]))
        for m in range(d.shape[0]):
            for n in range(d.shape[1]):
                if d[m][n] < dist:
                    if i1[m] < i2[n]:
                        pairs.add((i1[m], i2[n]))
    return pairs


def _get_vis_and_elbow(bench, tpos):
    from ics.cobraOps.TargetGroup import TargetGroup
    from ics.cobraOps.TargetSelector import TargetSelector
    tgroup = TargetGroup(np.array(tpos))
    tselect = TargetSelector(bench, tgroup)
    tselect.calculateAccessibleTargets()
    tmp = tselect.accessibleTargetIndices
    elb = tselect.accessibleTargetElbows
    res = defaultdict(list)

    for cbr in range(tmp.shape[0]):
        for i, tidx in enumerate(tmp[cbr, :]):
            if tidx >= 0:
                res[tidx].append((cbr, elb[cbr, i]))
    return res


def _get_elbow_collisions(bench, tpos, vis, dist):
    tpos = np.array(tpos)
    ivis = defaultdict(list)
    epos = defaultdict(list)
    for tidx, cidx_elbow in vis.items():
        for (cidx, elbowpos) in cidx_elbow:
            ivis[cidx].append(tidx)
            epos[cidx].append((tidx, elbowpos))

    res = defaultdict(list)
    for cidx, thing in epos.items():
        # determine target indices visible by neighbors of this cobra
        nb = bench.getCobraNeighbors(cidx)
        i2 = np.concatenate([ivis[j] for j in nb if j in ivis])
        i2 = np.unique(i2).astype(np.int)
        # for each target visible by this cobra and the corresponding elbow
        # position, find all targets which are too close to the "upper arm"
        # of the cobra
        for tidx, elbowpos in thing:
            ebp = np.full(len(i2), elbowpos)
            tp = np.full(len(i2), tpos[tidx])
            ti2 = tpos[i2]
            d = bench.distancesToLineSegments(ti2, tp, ebp)
            res[(cidx, tidx)] += list(i2[d<dist])
    return res


def observeWithNetflow(bench, targets, tpos, classdict, tvisit, vis_cost=None,
                       cobraMoveCost=None, collision_distance=0.,
                       elbow_collisions=True, gurobi=False):
    Cv_i = defaultdict(list)  # Cobra visit inflows
    Tv_o = defaultdict(list)  # Target visit outflows
    Tv_i = defaultdict(list)  # Target visit inflows
    T_o = defaultdict(list)  # Target outflows (only science targets)
    T_i = defaultdict(list)  # Target inflows (only science targets)
    CTCv_o = defaultdict(list)  # Calibration Target class visit outflows
    STC_o = defaultdict(list)  # Science Target outflows
    if gurobi:
        import gurobipy as gbp
        prob = gbp.Model("problem")
        cost = prob.addVar(vtype=gbp.GRB.CONTINUOUS, name="cost")
        lpSum = gbp.quicksum

        def add_constraint(problem, constraint):
            problem.addConstr(constraint)

        def varValue(var):
            return var.X
    else:
        import pulp
        prob = pulp.LpProblem("problem", pulp.LpMinimize)
        cost = pulp.LpVariable("cost", 0)
        lpSum = pulp.lpSum

        def add_constraint(problem, constraint):
            problem += constraint

        def varValue(var):
            return pulp.value(var)

    nreqvisit = []
    for t in targets:
        if isinstance(t, ScienceTarget):
            nreqvisit.append(int(t.obs_time/tvisit))
        else:
            nreqvisit.append(0)

    nvisits = len(tpos)

    if vis_cost is None:
        vis_cost = [0.]*nvisits

    def newvar(lo, hi):
        newvar._varcount += 1
        if gurobi:
            if lo is None:
                lo = -gbp.GRB.INFINITY
            if hi is None:
                hi = gbp.GRB.INFINITY
            if lo == 0 and hi == 1:
                return prob.addVar(vtype=gbp.GRB.BINARY,
                                   name="v{}".format(newvar._varcount))
            else:
                return prob.addVar(vtype=gbp.GRB.INTEGER,
                                   name="v{}".format(newvar._varcount),
                                   lb=lo, ub=hi)
        else:
            return pulp.LpVariable("v{}".format(newvar._varcount), lo, hi,
                                   cat=pulp.LpInteger)
    newvar._varcount = 0

    # define LP variables

    print("Creating network topology")
    # create nodes for every visit and calibration target class
    for key, value in classdict.items():
        if value["calib"]:
            for ivis in range(nvisits):
                f = newvar(0, None)
                CTCv_o[(key, ivis)].append(f)
                cost += f*value["nonObservationCost"]

    constr = []
    for ivis in range(nvisits):
        print("  exposure {}".format(ivis+1))
        print("Calculating visibilities")
        vis = _get_vis_and_elbow(bench, tpos[ivis])
        for tidx, thing in vis.items():
            tgt = targets[tidx]
            TC = tgt.targetclass
            Class = classdict[TC]
            if isinstance(tgt, ScienceTarget):
                # Target node to target visit node
                f = newvar(0, 1)
                T_o[tidx].append(f)
                Tv_i[(tidx, ivis)].append(f)
                if len(T_o[tidx]) == 1:  # freshly created
                    # Science Target class node to target node
                    f = newvar(0, 1)
                    T_i[tidx].append(f)
                    STC_o[TC].append(f)
                    if len(STC_o[TC]) == 1:  # freshly created
                        # Science Target class node to sink
                        f = newvar(0, None)
                        STC_o[TC].append(f)
                        cost += f*Class["nonObservationCost"]
                    # Science Target node to sink
                    f = newvar(0, None)
                    T_o[tidx].append(f)
                    cost += f*Class["partialObservationCost"]
            elif isinstance(tgt, CalibTarget):
                # Calibration Target class node to target visit node
                f = newvar(0, 1)
                Tv_i[(tidx, ivis)].append(f)
                CTCv_o[(TC, ivis)].append(f)
            for (cidx, _) in thing:
                # target visit node to cobra visit node
                f = newvar(0, 1)
                Cv_i[(cidx, ivis)].append(f)
                Tv_o[(tidx, ivis)].append((f, cidx))
                tcost = vis_cost[ivis]
                if cobraMoveCost is not None:
                    dist = np.abs(bench.cobras.centers[cidx]-tpos[ivis][tidx])
                    tcost += cobraMoveCost(dist)
                cost += f*tcost

        # Constraints
        print("adding constraints")
        # avoid endpoint collisions
        if collision_distance > 0.:
            print("adding collision constraints")
            if not elbow_collisions:
                cpair = _get_colliding_pairs(bench, tpos[ivis], vis,
                                             collision_distance)
                keys = Tv_o.keys()
                keys = set(key[0] for key in keys if key[1] == ivis)
                for p in cpair:
                    if p[0] in keys and p[1] in keys:
                        flows = [v[0] for v in
                                 Tv_o[(p[0], ivis)] + Tv_o[(p[1], ivis)]]
                        constr.append(lpSum(flows) <= 1)
            else:
                elbowcoll = _get_elbow_collisions(bench, tpos[ivis], vis,
                                                  collision_distance)
                for (cidx, tidx1), tidx2 in elbowcoll.items():
                    for f2, cidx2 in Tv_o[(tidx1, ivis)]:
                        if cidx2 == cidx:
                            flow0 = f2
                    for idx2 in tidx2:
                        if True:#idx2 != tidx1:
                            flows = [flow0]
                            flows += [f2 for f2, cidx2 in Tv_o[(idx2, ivis)]
                                      if cidx2 != cidx]
                            constr.append(lpSum(flows) <= 1)

    # NOTE: at least in Pulp, it is important to set the objective funcation
    # before any of the constraints!
    if gurobi:
        prob.setObjective(cost)
    else:
        prob += cost

    for c in constr:
        add_constraint(prob, c)

    # every Cobra can observe at most one target per visit
    for inflow in Cv_i.values():
        add_constraint(prob, lpSum([f for f in inflow]) <= 1)

    # every calibration target class must be observed a minimum number of times
    # every visit
    for key, value in CTCv_o.items():
        add_constraint(prob, lpSum([v for v in value]) >=
                       classdict[key[0]]["numRequired"])

    # inflow and outflow at every Tv node must be balanced
    for key, ival in Tv_i.items():
        oval = Tv_o[key]
        add_constraint(prob,
                       lpSum([v for v in ival]+[-v[0] for v in oval]) == 0)

    # inflow and outflow at every T node must be balanced
    for key, ival in T_i.items():
        oval = T_o[key]
        nvis = nreqvisit[key]
        add_constraint(prob,
                       lpSum([nvis*v for v in ival]+[-v for v in oval]) == 0)

    # Science targets must be either observed or go to the sink
    for key, val in STC_o.items():
        add_constraint(prob, lpSum([v for v in val]) == len(val)-1)

    print("solving the problem")
    if gurobi:
        prob.ModelSense = 1  # minimize
        prob.optimize()
    else:
        status = prob.solve(pulp.COIN_CMD(msg=1, keepFiles=0, maxSeconds=100,
                                          threads=1, dual=10.))

    res = [{} for _ in range(nvisits)]
    for k1, v1 in Tv_o.items():
        for i2 in v1:
            visited = varValue(i2[0]) > 0
            if visited:
                tidx, ivis = k1
                cidx = i2[1]
                res[ivis][tidx] = cidx
    return res


class Telescope(object):
    """An object describing a telescope configuration to be used for observing
    a target field. Includes a list of Cobras, telescope RA/Dec and position
    angle and an observation time according to ISO8601 UTC.
    The minimum allowed distance between two Cobra tips is also stored here.
    Based on this information, and given a list of targets, this object can
    compute assignment strategies using different algorithms (currently network
    flow and ETs approaches).
    """
    def __init__(self, ra, dec, posang, time):
        self._ra = float(ra)
        self._dec = float(dec)
        self._posang = float(posang)
        self._time = str(time)

    def subtract_obs_time(self, tgt, obs, tvisit):
        res = []
        for k in obs.keys():
            tgt[k].reduce_time(tvisit)
        for t in tgt:
            if t.obs_time > 0.:
                res.append(t)
        return res

    def get_fp_positions(self, tgt):
        return [t.fp_position(self._ra, self._dec, self._posang, self._time)
                for t in tgt]


class Target(object):
    """Base class for all target types observable with PFS. All targets are
    initialized with RA/Dec and an ID string. From RA/Dec, the target can
    determine its position on the focal plane, once also the telescope attitude
    is known. """
    def __init__(self, ID, ra, dec, targetclass):
        self._ID = str(ID)
        self._ra = float(ra)
        self._dec = float(dec)
        self._targetclass = targetclass

    @property
    def ID(self):
        """ID of the object : str"""
        return self._ID

    @property
    def ra(self):
        """the rectascension : float"""
        return self._ra

    @property
    def dec(self):
        """the declination : float"""
        return self._dec

    def fp_position(self, raTel, decTel, posang, time):
        return pycconv.cconv([self._ra], [self._dec],
                             raTel, decTel, posang, time)[0]

    @property
    def targetclass(self):
        """string representation of the target's class"""
        return self._targetclass


class ScienceTarget(Target):
    """Derived from the Target class, with the additional attributes priority
    and observation time.

    All different types of ScienceTarget need to be derived from this class."
    """
    def __init__(self, ID, ra, dec, obs_time, pri, prefix):
        super(ScienceTarget, self).__init__(ID, ra, dec,
                                            "{}_P{}".format(prefix, pri))
        self._obs_time = float(obs_time)
        self._pri = int(pri)

    @property
    def obs_time(self):
        """required observation time : float"""
        return self._obs_time

    def reduce_time(self, dt):
        self._obs_time -= float(dt)


class CalibTarget(Target):
    """Derived from the Target class."""


def telescopeRaDecFromFile(file):
    with open(file) as f:
        ras = []
        decs = []
        ll = f.readlines()
        for l in ll[1:]:
            if not l.startswith("#"):
                tt = l.split()
                ra, dec = (float(tt[1]), float(tt[2]))
                ras.append(ra)
                decs.append(dec)
    return float(np.average(ras)), float(np.average(decs))


def readScientificFromFile(file, prefix):
    with open(file) as f:
        res = []
        ll = f.readlines()
        for l in ll[1:]:
            if not l.startswith("#"):
                tt = l.split()
                id_, ra, dec, tm, pri = (
                    str(tt[0]), float(tt[1]), float(tt[2]),
                    float(tt[3]), int(tt[4]))
                res.append(ScienceTarget(id_, ra, dec, tm, pri, prefix))
    return res


def readCalibrationFromFile(file, targetclass):
    with open(file) as f:
        res = []
        ll = f.readlines()
        for l in ll[1:]:
            if not l.startswith("#"):
                tt = l.split()
                id_, ra, dec = (str(tt[0]), float(tt[1]), float(tt[2]))
                res.append(CalibTarget(id_, ra, dec, targetclass))
    return res
