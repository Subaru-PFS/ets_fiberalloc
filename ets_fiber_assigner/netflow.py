import numpy as np
from collections import defaultdict
from astropy.table import Table


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
        i2 = np.unique(i2).astype(int)
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

    class DummyTargetSelector(TargetSelector):
        def run(self):
            return

        def selectTargets(self):
            return

    tgroup = TargetGroup(np.array(tpos))
    tselect = DummyTargetSelector(bench, tgroup)
    tselect.calculateAccessibleTargets()
    tmp = tselect.accessibleTargetIndices
    elb = tselect.accessibleTargetElbows
    res = defaultdict(list)

#    observable_targets = set()
    for cbr in range(tmp.shape[0]):
        #        observable_targets = observable_targets | set(tmp[cbr, :])
        for i, tidx in enumerate(tmp[cbr, :]):
            if tidx >= 0:
                res[tidx].append((cbr, elb[cbr, i]))

#    nonobservable_targets = set(range(len(tpos))).difference(observable_targets)
#    return res, nonobservable_targets
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
        tmp = [ivis[j] for j in nb if j in ivis]
        if tmp != []:
            i2 = np.concatenate([ivis[j] for j in nb if j in ivis])
            i2 = np.unique(i2).astype(int)
            # for each target visible by this cobra and the corresponding elbow
            # position, find all targets which are too close to the "upper arm"
            # of the cobra
            for tidx, elbowpos in thing:
                ebp = np.full(len(i2), elbowpos)
                tp = np.full(len(i2), tpos[tidx])
                ti2 = tpos[i2]
                d = bench.distancesToLineSegments(ti2, tp, ebp)
                res[(cidx, tidx)] += list(i2[d < dist])
    return res


def getClosestDots(bench):
    res = []
    for cidx in range(len(bench.cobras.centers)):
        nb = bench.getCobraNeighbors(cidx)
        res.append(bench.blackDots.centers[[cidx] + list(nb)])
    return res


class LPProblem(object):
    def __init__(self):
        self._vardict = {}
        self._constraintdict = {}

    def varByName(self, name):
        return self._vardict[name]

    def constraintByName(self, name):
        return self._constraintdict[name]


class GurobiProblem(LPProblem):
    def __init__(self, name="problem", extraOptions=None):
        LPProblem.__init__(self)
        import gurobipy as gbp
        self._prob = gbp.Model(name)
        self.cost = self._prob.addVar(name="cost", vtype=gbp.GRB.CONTINUOUS)
        self._prob.ModelSense = 1  # minimize
        if extraOptions is not None:
            for key, value in extraOptions.items():
                self._prob.setParam(key, value)
        self._prob.update()
        self.sum = gbp.quicksum

    def addVar(self, name, lo, hi):
        import gurobipy as gbp
        if lo is None:
            lo = -gbp.GRB.INFINITY
        if hi is None:
            hi = gbp.GRB.INFINITY
        if lo == 0 and hi == 1:
            var = self._prob.addVar(name=name, vtype=gbp.GRB.BINARY)
        else:
            var = self._prob.addVar(lb=lo, ub=hi, name=name,
                                    vtype=gbp.GRB.INTEGER)
        self._vardict[name] = var
        return var

    def add_constraint(self, name, constraint):
        self._constraintdict[name] = constraint
        self._prob.addConstr(constraint)

    def add_lazy_constraint(self, name, constraint):
        constraint.Lazy = 1
        self.add_constraint(name, constraint)

    @staticmethod
    def value(var):
        return var.X

    def solve(self):
        self._prob.setObjective(self.cost)
        self._prob.optimize()

    def update(self):
        self._prob.update()

    def dump(self, filename):
        self._prob.write(filename)

    @staticmethod
    def varBounds(var):
        return var.lb, var.ub

    @staticmethod
    def changeVarBounds(var, lower=None, upper=None):
        if lower is not None:
            var.lb = lower
        if upper is not None:
            var.ub = upper


class PulpProblem(LPProblem):
    def __init__(self, name="problem"):
        LPProblem.__init__(self)
        import pulp
        self._prob = pulp.LpProblem("problem", pulp.LpMinimize)
        self.cost = pulp.LpVariable("cost", 0)
        self.sum = pulp.lpSum
        self._constr = []

    def addVar(self, name, lo, hi):
        import pulp
        if lo == 0 and hi == 1:
            var = pulp.LpVariable(name, cat=pulp.LpBinary)
        else:
            var = pulp.LpVariable(name, lo, hi, cat=pulp.LpInteger)
        self._vardict[name] = var
        return var

    def add_constraint(self, name, constraint):
        self._constraintdict[name] = constraint
        self._constr.append(constraint)

    def add_lazy_constraint(self, name, constraint):
        self.add_constraint(name, constraint)

    @staticmethod
    def value(var):
        import pulp
        return pulp.value(var)

    def solve(self):
        import pulp
        self._prob += self.cost
        for i in self._constr:
            self._prob += i
        self._prob.solve(pulp.COIN_CMD(msg=1, keepFiles=0, maxSeconds=100,
                                       threads=1))

    def update(self):
        pass

    def dump(self, filename):
        self._prob.writeLP(filename)

    @staticmethod
    def varBounds(var):
        return var.lowBound, var.upBound

    @staticmethod
    def changeVarBounds(var, lower=None, upper=None):
        if lower is not None:
            var.lowBound = lower
        if upper is not None:
            var.upBound = upper


def makeName(*stuff):
    return "_".join([str(x) for x in stuff])


def buildProblem(bench, targets, tpos, classdict, tvisit, vis_cost=None,
                 cobraMoveCost=None, collision_distance=0.,
                 elbow_collisions=True, gurobi=True, gurobiOptions=None,
                 alreadyObserved=None, forbiddenPairs=None,
                 cobraLocationGroup=None, minSkyTargetsPerLocation=None,
                 locationGroupPenalty=None,
                 cobraInstrumentRegion=None, minSkyTargetsPerInstrumentRegion=None,
                 instrumentRegionPenalty=None, blackDotPenalty=None):
    """Build the ILP problem for a given observation task

    Parameters
    ==========
    bench : ics.cobraOps.Bench.Bench
        description of the focal plane
    targets : list of Target objects
        all targets available for observation,
        i.e. science, sky and calibration targets
    tpos : list of list of complex
        focal plane positions for all targets (inner list) and all visits
        (outer list)
    classdict : dict(string : dict(string : <param>))
        properties of target classes
        Each target class *must* have the following entries:
         - "calib" : bool
           whether the target class is a calibration target class
           (otherwise it is a science target class)
         - "nonObservationCost" : float
           the penalty for each object in that class which should be observed,
           but isn't
        Each calibration target class has the following additional entries:
         - "numRequired" : int
           the number of targets from that class which must be observed
           *during every exposure*
        Each science target class has the following additional entries:
         - "partialObservationCost" : float
           the penalty for every target of that class which is observed, but
           not to its full required time.
           This *must* be higher than "nonObservationCost".
         - "nobs_max" : integer, optional
           The number of members of this class that should be observed (default:
           all of them).
    tvisit : float
        duration of a single visit in seconds
    vis_cost : list of float (nvisit entries)
        cost for observing a target in a given visit. This can be used to make
        the solver prefer earlier over later observations
    cobraMoveCost : function taking a float and returning a float
        extra cost when observing a target at a given distance from the center
        of the observing cobra. Can be used to make Cobras prefer targets closer
        to the center of their patrol region.
    collision_distance : float
        collision distance between cobra tips
    elbow_collisions : bool
        if True, avoid elbow collisions in the endpoint configuration
        (increases the number of constraints, especially for long target lists)
    gurobi : bool
        if True, use the Gurobi optimizer, otherwise use PuLP
    gurobiOptions : dict(string : <param>)
        optional additional parameters for the Gurobi solver
    alreadyObserved : None or dict{string: int}
        if not None, this is a dictionary containing IDs of science targets
        and the number of visits they have already been observed
    forbiddenPairs : None or list(list(tuple of 2 ints))
        Pairs of targets that cannot be both observed during the same visit,
        because this would lead to trajectory collisions
    cobraLocationGroup : integer, array-like
        if provided, this must be indexable with the Cobra indices from "bench"
        and return the "location group index" of the respective Cobra.
        As an example, the focal plane could be divided into 7 hexagonal
        sub-areas with indices 0 to 6.
    minSkyTargetsPerLocation : integer
        how many sky targets have to be observed in every location group
    locationGroupPenalty : float
        how much to increase the cost function for every "missing" sky target
        in a focal plane region
    cobraInstrumentRegion : integer, array-like
        if provided, this must be indexable with the Cobra indices from "bench"
        and return the "instrument region index" of the respective Cobra.
        Instrument regions will typically be used to describe continuous regions
        along a spectrograph slit.
    minSkyTargetsPerInstrumentRegion : integer
        how many sky targets have to be observed in every instrument region
    instrumentRegionPenalty : float
        how much to increase the cost function for every "missing" sky target
        in an instrument region
    blackDotPenalty : function taking a float and returning a float
        extra cost when observing a target at a given distance from the center
        of a black dot. Can be used to make Cobras prefer targets farther away
        from a black dot to reduce vignetting.
        The distance parameter is expected in millimeters.

    Returns
    =======
    LPProblem : the ILP problem object
    """
    Cv_i = defaultdict(list)  # Cobra visit inflows
    Tv_o = defaultdict(list)  # Target visit outflows
    Tv_i = defaultdict(list)  # Target visit inflows
    T_o = defaultdict(list)  # Target outflows (only science targets)
    T_i = defaultdict(list)  # Target inflows (only science targets)
    CTCv_o = defaultdict(list)  # Calibration Target class visit outflows
    STC_o = defaultdict(list)  # Science Target outflows

    if gurobi:
        prob = GurobiProblem(extraOptions=gurobiOptions)
    else:
        prob = PulpProblem()

    if blackDotPenalty is not None:
        closestDotsList = getClosestDots(bench)

    nreqvisit = []
    ndone = []
    for t in targets:
        if isinstance(t, ScienceTarget):
            nreqvisit.append(int(t.obs_time/tvisit))  # FIXME: maybe round up?
            tmp = 0
            if alreadyObserved is not None:
                if t.ID in alreadyObserved:
                    tmp = alreadyObserved[t.ID]
            ndone.append(tmp)
        else:
            nreqvisit.append(0)
            ndone.append(0)

    nvisits = len(tpos)

    # sanity check for science targets: make sure that partialObservationCost
    # is larger than nonObservationCost
    for key, val in classdict.items():
        if not val["calib"]:
            if val["partialObservationCost"] < val["nonObservationCost"]:
                raise ValueError(
                    "found a target class where partialObservationCost "
                    "is smaller than nonObservationCost")

    if cobraLocationGroup is not None:
        maxLocGroup = max(cobraLocationGroup)
        locationVars = [[[] for _ in range(maxLocGroup+1)] for _ in range(nvisits)]
        # add overflow arcs.
        for i in range(nvisits):
            for j in range(maxLocGroup+1):
                f = prob.addVar(makeName("locgroup_sink", i, j), 0, None)
                prob.cost += f*locationGroupPenalty
                locationVars[i][j].append(f)

    if cobraInstrumentRegion is not None:
        maxInstRegion = max(cobraInstrumentRegion)
        regionVars = [[[] for _ in range(maxInstRegion+1)] for _ in range(nvisits)]
        # add overflow arcs.
        for i in range(nvisits):
            for j in range(maxInstRegion+1):
                f = prob.addVar(makeName("instregion_sink", i, j), 0, None)
                prob.cost += f*instrumentRegionPenalty
                regionVars[i][j].append(f)

    if vis_cost is None:
        vis_cost = [0.] * nvisits

    # define LP variables

    print("Creating network topology")
    # create nodes for every visit and calibration target class
    for key, value in classdict.items():
        if value["calib"]:
            for ivis in range(nvisits):
                f = prob.addVar(makeName("CTCv_sink", key, ivis), 0, None)
                CTCv_o[(key, ivis)].append(f)
                prob.cost += f*value["nonObservationCost"]

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
                f = prob.addVar(makeName("T_Tv", tgt.ID, ivis), 0, 1)
                T_o[tidx].append(f)
                Tv_i[(tidx, ivis)].append(f)
                if len(T_o[tidx]) == 1:  # freshly created
                    # Science Target class node to target node
                    tmp = 1 if ndone[tidx] > 0 else 0
                    f = prob.addVar(makeName("STC_T", TC, tgt.ID), tmp, 1)
                    T_i[tidx].append(f)
                    STC_o[TC].append(f)
                    if len(STC_o[TC]) == 1:  # freshly created
                        # Science Target class node to sink
                        f = prob.addVar(makeName("STC_sink", TC), 0, None)
                        STC_o[TC].append(f)
                        prob.cost += f*Class["nonObservationCost"]
                    # Science Target node to sink
                    f = prob.addVar(makeName("ST_sink", tgt.ID), 0, None)
                    T_o[tidx].append(f)
                    prob.cost += f*Class["partialObservationCost"]
            elif isinstance(tgt, CalibTarget):
                # Calibration Target class node to target visit node
                f = prob.addVar(makeName("CTCv_Tv", TC, tgt.ID, ivis), 0, 1)
                if tgt._penalty != 0:
                    prob.cost += f*tgt._penalty
                Tv_i[(tidx, ivis)].append(f)
                CTCv_o[(TC, ivis)].append(f)
            for (cidx, _) in thing:
                # target visit node to cobra visit node
                f = prob.addVar(makeName("Tv_Cv", tidx, cidx, ivis), 0, 1)
                Cv_i[(cidx, ivis)].append(f)
                Tv_o[(tidx, ivis)].append((f, cidx))
                if cobraLocationGroup is not None \
                        and isinstance(tgt, CalibTarget) \
                        and tgt.targetclass == "sky":
                    locationVars[ivis][cobraLocationGroup[cidx]].append(f)
                if cobraInstrumentRegion is not None \
                        and isinstance(tgt, CalibTarget) \
                        and tgt.targetclass == "sky":
                    regionVars[ivis][cobraInstrumentRegion[cidx]].append(f)
                tcost = vis_cost[ivis]
                if cobraMoveCost is not None:
                    dist = np.abs(bench.cobras.centers[cidx]-tpos[ivis][tidx])
                    tcost += cobraMoveCost(dist)
                prob.cost += f*tcost
                if blackDotPenalty is not None:
                    dist = np.min(np.abs(closestDotsList[cidx]-tpos[ivis][tidx]))
                    tcost += blackDotPenalty(dist)
                prob.cost += f*tcost

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
                        tname0 = targets[p[0]].ID
                        tname1 = targets[p[1]].ID
                        constr.append([makeName("Coll_", tname0, tname1, ivis),
                                       prob.sum(flows) <= 1])
            else:
                elbowcoll = _get_elbow_collisions(bench, tpos[ivis], vis,
                                                  collision_distance)
                for (cidx, tidx1), tidx2 in elbowcoll.items():
                    for f2, cidx2 in Tv_o[(tidx1, ivis)]:
                        if cidx2 == cidx:
                            flow0 = f2
                    for idx2 in tidx2:
                        if True:  # idx2 != tidx1:
                            flows = [flow0]
                            flows += [f2 for f2, cidx2 in Tv_o[(idx2, ivis)]
                                      if cidx2 != cidx]
                            constr.append([makeName("Coll_", cidx, cidx2, ivis),
                                           prob.sum(flows) <= 1])

        # add constraints for forbidden pairs of targets
        if forbiddenPairs is not None:
            print("adding forbidden pair constraints")
            keys = Tv_o.keys()
            keys = set(key[0] for key in keys if key[1] == ivis)
            for p in forbiddenPairs[ivis]:
                if len(p) == 2:
                    if p[0] in keys and p[1] in keys:
                        flows = [v[0] for v in
                                 Tv_o[(p[0], ivis)] + Tv_o[(p[1], ivis)]]
                        tname0 = targets[p[0]].ID
                        tname1 = targets[p[1]].ID
                        constr.append([makeName("forbiddenPair_", tname0, tname1, ivis),
                                       prob.sum(flows) <= 1])
                elif len(p) == 1:
                    if p[0] in keys:
                        flows = [v[0] for v in Tv_o[(p[0], ivis)]]
                        tname0 = targets[p[0]].ID
                        constr.append([makeName("forbiddenPair_", tname0, ivis),
                                       prob.sum(flows) == 0])
                else:
                    raise RuntimeError("oops")

    for c in constr:
        # We add the collision constraints as lazy in the hope that this will
        # speed up the solution
        prob.add_lazy_constraint(c[0], c[1])

    # every Cobra can observe at most one target per visit
    for key, inflow in Cv_i.items():
        prob.add_constraint(makeName("Cvlim", key[0], key[1]),
                            prob.sum([f for f in inflow]) <= 1)

    # every calibration target class must be observed a minimum number of times
    # every visit
    for key, value in CTCv_o.items():
        prob.add_constraint(makeName("CTCv_min", key[0], key[1]),
                            prob.sum([v for v in value]) >= classdict[key[0]]["numRequired"])

    # inflow and outflow at every Tv node must be balanced
    for key, ival in Tv_i.items():
        oval = Tv_o[key]
        prob.add_constraint(makeName("TvIO", targets[key[0]].ID, key[1]),
                            prob.sum([v for v in ival]+[-v[0] for v in oval]) == 0)

    # inflow and outflow at every T node must be balanced
    for key, ival in T_i.items():
        oval = T_o[key]
        nvis = nreqvisit[key]-ndone[key]
        if nvis < 0:
            nvis = 0
        prob.add_constraint(makeName("TIO", targets[key].ID),
                            prob.sum([nvis*v for v in ival]+[-v for v in oval]) == 0)

    # Science targets must be either observed or go to the sink
    for key, val in STC_o.items():
        n_obs = len(val)-1
        if "nobs_max" in classdict[key]:
            n_obs = classdict[key]["nobs_max"]
        prob.add_constraint(makeName("ST", key[0], key[1]),
                            prob.sum([v for v in val]) == n_obs)

    # Make sure that there are enough sky targets in every Cobra location group
    if cobraLocationGroup is not None:
        for ivis in range(nvisits):
            for i in range(maxLocGroup+1):
                prob.add_constraint(makeName("LocGrp", ivis, i),
                                    prob.sum([v for v in locationVars[ivis][i]]) >= minSkyTargetsPerLocation)

    # Make sure that there are enough sky targets in every Cobra instrument region
    if cobraInstrumentRegion is not None:
        for ivis in range(nvisits):
            for i in range(maxInstRegion+1):
                prob.add_constraint(makeName("InstReg", ivis, i),
                                    prob.sum([v for v in regionVars[ivis][i]]) >= minSkyTargetsPerInstrumentRegion)

    return prob


class Telescope(object):
    """An object describing a telescope configuration to be used for observing
    a target field.

    Parameters
    ==========
    ra : float
        right ascension of the focal plane center
    dec : float
        declination of the focal plane center
    posang : float
        position angle of the focal plane
    time : string
        time of observation in ISO8601 UTC format
    """

    def __init__(self, ra, dec, posang, time):
        self._ra = float(ra)
        self._dec = float(dec)
        self._posang = float(posang)
        self._time = str(time)

    def get_fp_positions(self, tgt):
        """Returns positions on the focal plane for a list of targets

        Parameters
        ==========
        tgt : list of Target objects
            the targets whose focal plane position is requested

        Returns
        =======
        list of complex numbers
            the focal plane positions encoded as complex numbers
        """
        from pfs.utils.coordinates.CoordTransp import CoordinateTransform as ctrans
        tmp = np.zeros((2, len(tgt)), dtype=np.float64)
        for i, t in enumerate(tgt):
            tmp[0, i], tmp[1, i] = t.ra, t.dec
        tmp = ctrans(
            xyin=tmp,
            mode="sky_pfi",
            pa=self._posang,
            cent=np.array([self._ra, self._dec]).reshape((2, 1)),
            pm=np.stack([[t.pmra for t in tgt], [t.pmdec for t in tgt]], axis=0),
            par=np.array([t.parallax for t in tgt]),
            time=self._time,
            # epoch=np.array([t.epoch for t in tgt]),
            epoch=2000.0,
        )

        return tmp[0, :] + 1j*tmp[1, :]


class Target(object):
    """Base class for all target types observable with PFS. All targets are
    initialized with RA/Dec and an ID string. From RA/Dec, the target can
    determine its position on the focal plane, once also the telescope attitude
    is known. """

    def __init__(self, ID, ra, dec, targetclass, pmra=0, pmdec=0, parallax=0, epoch=2000.0):
        self._ID = str(ID)
        self._ra = float(ra)
        self._dec = float(dec)
        self._pmra = float(pmra)
        self._pmdec = float(pmdec)
        self._parallax = float(parallax)
        self._epoch = float(epoch)
        self._targetclass = targetclass

    @property
    def ID(self):
        """ID of the object : str"""
        return self._ID

    @property
    def ra(self):
        """the rectascension in degrees : float"""
        return self._ra

    @property
    def dec(self):
        """the declination in degrees : float"""
        return self._dec

    @property
    def pmra(self):
        """the rectascension component of proper motion in mas/yr : float"""
        return self._pmra

    @property
    def pmdec(self):
        """the declination component of proper motion in mas/yr : float"""
        return self._pmdec

    @property
    def parallax(self):
        """the parallax in mas : float"""
        return self._parallax

    @property
    def epoch(self):
        """epoch : float"""
        return self._epoch

    @property
    def targetclass(self):
        """string representation of the target's class"""
        return self._targetclass


class ScienceTarget(Target):
    """Derived from the Target class, with the additional attributes priority
    and observation time.

    All different types of ScienceTarget need to be derived from this class."
    """

    def __init__(self, ID, ra, dec, obs_time, pri, prefix, pmra=0, pmdec=0, parallax=0, epoch=2000.):
        super(ScienceTarget, self).__init__(ID, ra, dec,
                                            "{}_P{}".format(prefix, pri),
                                            pmra=pmra, pmdec=pmdec, parallax=parallax, epoch=epoch
                                            )
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

    def __init__(self, ID, ra, dec, targetclass, penalty=0., pmra=0, pmdec=0, parallax=0, epoch=2000.):
        super(CalibTarget, self).__init__(ID, ra, dec, targetclass,
                                          pmra=pmra, pmdec=pmdec, parallax=parallax, epoch=epoch
                                          )
        self._penalty = penalty


def telescopeRaDecFromFile(file):
    """Obtain a heuristic telescope pointing from a file with science targets

    Parameters
    ==========
    file : string
        the name os the text file containing the target information

    Returns
    =======
    float, float : the approximate center of the targets in the file

    Notes
    =====
    For testing/demonstration purposes only. DO NOT USE FOR SERIOUS WORK.
    """
    try:
        # first try reading as ecsv format
        t = Table.read(file, format="ascii.ecsv")
        return float(np.average(t["R.A."])), float(np.average(t['Dec.']))
    except:
        pass

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
    """Read a set of scientific targets from a file

    Parameters
    ==========
    file : string
        the name os the text file containing the target information
    prefix : string
        the name of the target class to which these targets will belong

    Returns
    =======
    list of ScienceTarget : the created ScienceTarget objects
    """

    # FIXME: add a routine that can read proper motion and parallax information
    #        from the file, as soon as the file format has been defined

    try:
        # first try reading as ecsv format
        t = Table.read(file, format="ascii.ecsv")
        res = []
        for r in t:
            res.append(ScienceTarget(r["ID"], r["R.A."], r["Dec."],
                       r["Exposure Time"], r["Priority"], prefix))
        return res
    except:
        pass

    # try legacy format
    with open(file) as f:
        res = []
        ll = f.readlines()
        for l in ll[1:]:
            if not l.startswith("#"):
                tt = l.split()
                id_, ra, dec, tm, pri = (
                    str(tt[0])[1:], float(tt[1]), float(tt[2]),
                    float(tt[3]), int(tt[4]))
                res.append(ScienceTarget(id_, ra, dec, tm, pri, prefix))
    return res


def readCalibrationFromFile(file, targetclass):
    """Read a set of calibration targets from a file

    Parameters
    ==========
    file : string
        the name os the text file containing the target information
    targetclass : string
        the name of the target class to which these targets will belong

    Returns
    =======
    list of CalibrationTarget : the created CalibrationTarget objects
    """
    try:
        # first try reading as ecsv format
        t = Table.read(file, format="ascii.ecsv")
        res = []
        for r in t:
            res.append(CalibTarget(r["ID"], r["R.A."], r["Dec."], targetclass))
        return res
    except:
        pass

    # try legacy format
    with open(file) as f:
        res = []
        ll = f.readlines()
        for l in ll[1:]:
            if not l.startswith("#"):
                tt = l.split()
                id_, ra, dec = (str(tt[0])[1:], float(tt[1]), float(tt[2]))
                res.append(CalibTarget(id_, ra, dec, targetclass))
    return res


def readCalibrationWithPenaltyFromFile(file, targetclass):
    """Read a set of calibration targets from a file

    Parameters
    ==========
    file : string
        the name os the text file containing the target information
    targetclass : string
        the name of the target class to which these targets will belong

    Returns
    =======
    list of CalibrationTarget : the created CalibrationTarget objects
    """
    t = Table.read(file, format="ascii.ecsv")
    res = []
    for r in t:
        res.append(CalibTarget(r["ID"], r["R.A."], r["Dec."], targetclass, r["penalty"]))
    return res
