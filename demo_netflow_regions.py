import ets_fiber_assigner.netflow as nf
import numpy as np
from collections import defaultdict
from getBench import getBench

# make runs reproducible
np.random.seed(20)

# define locations of the input files
catalog_path = "data/"
fsky_pos = catalog_path+"pfs_preliminary_target_cosmology_sky.dat"

tgt = nf.readCalibrationFromFile(fsky_pos, "sky")

bench = getBench()

# point the telescope at the center of all science targets
raTel, decTel = nf.telescopeRaDecFromFile(fsky_pos)
posang = 0.
otime = "2016-04-03T08:00:00Z"
telescopes = []

# number of distinct observations
nvisit = 6

# generate randomly jittered telescope pointings for every observation
for _ in range(nvisit):
    telescopes.append(nf.Telescope(raTel+np.random.normal()*1e-2,
                      decTel+np.random.normal()*1e-2, posang, otime))

# get focal plane positions for all targets and all visits
tpos = [tele.get_fp_positions(tgt) for tele in telescopes]

# create the dictionary containing the costs and constraints for all classes
# of targets
classdict = {}
classdict["sky"] = {"numRequired": 240,
                    "nonObservationCost": 1e6, "calib": True}

# optional: slightly increase the cost for later observations,
# to observe as early as possible
vis_cost = [i*10. for i in range(nvisit)]


# optional: penalize assignments where the cobra has to move far out
def cobraMoveCost(dist):
    return 0.1*dist


# duration of one observation in seconds
t_obs = 300.

gurobiOptions = dict(seed=0, presolve=1, method=4, degenmoves=0,
                     heuristics=0.8, mipfocus=0, mipgap=1.0e-04)

# let's pretend that most targets have already been completely observed,
# and that the rest has been partially observed
alreadyObserved={}
for t in tgt:
    alreadyObserved[t.ID] = 3
for t in tgt[::10]:
    alreadyObserved[t.ID] = 1

# define two location regions arbitrarily
ncobras = bench.cobras.nCobras
cobraRegions=np.zeros(ncobras,dtype=np.int32)
cobraRegions[ncobras//2:] = 1

# compute observation strategy
prob = nf.buildProblem(bench, tgt, tpos, classdict, t_obs,
                       vis_cost, cobraMoveCost=cobraMoveCost,
                       collision_distance=2., elbow_collisions=True,
                       gurobi=False, gurobiOptions=gurobiOptions,
                       alreadyObserved=alreadyObserved,
                       cobraLocationGroup=cobraRegions,
                       minSkyTargetsPerLocation=40,
                       locationGroupPenalty=1e9)

# print("writing problem to file ", mpsName)
# prob.dump(mpsName)

print("solving the problem")
prob.solve()

# extract solution
res = [{} for _ in range(nvisit)]
for k1, v1 in prob._vardict.items():
    if k1.startswith("Tv_Cv_"):
        visited = prob.value(v1) > 0
        if visited:
            _, _, tidx, cidx, ivis = k1.split("_")
            res[int(ivis)][int(tidx)] = int(cidx)

# write output file
with open("output.txt", "w") as f:
    for i, (vis, tp, tel) in enumerate(zip(res, tpos, telescopes)):
        print("exposure {}:".format(i))
        print("  assigned Cobras: {}".format(len(vis)))
        tdict = defaultdict(int)
        f.write("# Exposure {}: duration {}s, RA: {}, Dec: {}, PA: {}\n".
                format(i+1, t_obs, tel._ra, tel._dec, tel._posang))
        f.write("# Target    Fiber          X          Y         "
                "RA        DEC\n")
        for tidx, cidx in vis.items():
            tdict[tgt[tidx].targetclass] += 1
            f.write("{:} {:6d} {:10.5f} {:10.5f} {:10.5f} {:10.5f}\n"
                    .format(tgt[tidx].ID, cidx+1, tp[tidx].real, tp[tidx].imag,
                            tgt[tidx].ra, tgt[tidx].dec))
        for cls, num in tdict.items():
            print("   {}: {}".format(cls, num))
