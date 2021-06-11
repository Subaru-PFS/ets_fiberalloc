import numpy as np
from collections import defaultdict

import ets_fiber_assigner.netflow as nf
import ets_fiber_assigner.io_helpers

from ics.cobraOps.Bench import Bench
from ics.cobraOps.TargetGroup import TargetGroup
from ics.cobraOps.CobrasCalibrationProduct import CobrasCalibrationProduct
from ics.cobraOps.CollisionSimulator2 import CollisionSimulator2
from procedures.moduleTest.cobraCoach import CobraCoach
from ics.cobraOps.cobraConstants import NULL_TARGET_POSITION, NULL_TARGET_ID
from ics.cobraOps import plotUtils


def getBench():
    import os
    from procedures.moduleTest.cobraCoach import CobraCoach
    os.environ["PFS_INSTDATA_DIR"] = "/home/martin/codes/pfs_instdata"
    cobraCoach = CobraCoach(
        "fpga", loadModel=False, trajectoryMode=True,
        rootDir="/home/martin/codes/efa/")
    cobraCoach.loadModel(version="ALL", moduleVersion="final_20210512")
    
    # Get the calibration product
    calibrationProduct = cobraCoach.calibModel
    
    # Set some dummy center positions and phi angles for those cobras that have
    # zero centers
    zeroCenters = calibrationProduct.centers == 0
    calibrationProduct.centers[zeroCenters] = np.arange(np.sum(zeroCenters)) * 300j
    calibrationProduct.phiIn[zeroCenters] = -np.pi
    calibrationProduct.phiOut[zeroCenters] = 0
    print("Cobras with zero centers: %i" % np.sum(zeroCenters))
    
    # Transform the calibration product cobra centers and link lengths units from
    # pixels to millimeters
    calibrationProduct.centers -= 5048.0 + 3597.0j
    calibrationProduct.centers *= np.exp(1j * np.deg2rad(1.0)) / 13.02
    calibrationProduct.L1 /= 13.02
    calibrationProduct.L2 /= 13.02
    
    # Use the median value link lengths in those cobras with zero link lengths
    zeroLinkLengths = np.logical_or(
        calibrationProduct.L1 == 0, calibrationProduct.L2 == 0)
    calibrationProduct.L1[zeroLinkLengths] = np.median(
        calibrationProduct.L1[~zeroLinkLengths])
    calibrationProduct.L2[zeroLinkLengths] = np.median(
        calibrationProduct.L2[~zeroLinkLengths])
    print("Cobras with zero link lenghts: %i" % np.sum(zeroLinkLengths))
    
    # Use the median value link lengths in those cobras with too long link lengths
    tooLongLinkLengths = np.logical_or(
        calibrationProduct.L1 > 100, calibrationProduct.L2 > 100)
    calibrationProduct.L1[tooLongLinkLengths] = np.median(
        calibrationProduct.L1[~tooLongLinkLengths])
    calibrationProduct.L2[tooLongLinkLengths] = np.median(
        calibrationProduct.L2[~tooLongLinkLengths])
    print("Cobras with too long link lenghts: %i" % np.sum(tooLongLinkLengths))
   
    # Create the bench instance
    bench = Bench(layout="calibration", calibrationProduct=calibrationProduct)
    print("Number of cobras:", bench.cobras.nCobras)

    return cobraCoach, bench


# make runs reproducible
np.random.seed(20)

# define locations of the input files
catalog_path = "data/"
fscience_targets = catalog_path+"test_sci.dat"
# So far, we only have test data for targets.
# Once we have files for calibration stars and sky locations, we can add them
# here.
fcal_stars = catalog_path+"test_cal.dat"
fsky_pos = catalog_path+"test_sky.dat"

# read all targets into a single list, giving them their proper types
tgt = nf.readScientificFromFile(fscience_targets, "sci")
# add calibration targets
tgt += nf.readCalibrationFromFile(fcal_stars, "cal")
tgt += nf.readCalibrationFromFile(fsky_pos, "sky")

# get a complete, idealized focal plane configuration
cobraCoach, bench = getBench()

#telescope = inputParamsFromPfsDesign("", ".")

# point the telescope at the center of all science targets
raTel, decTel = 0.0, 0.0
posang = 0.
otime = "2016-04-03T08:00:00Z"
telescopes = []

# number of distinct observations
nvisit = 1

# generate randomly jittered telescope pointings for every observation
for _ in range(nvisit):
    telescopes.append(nf.Telescope(raTel,
                      decTel, posang, otime))

# get focal plane positions for all targets and all visits
tpos = [tele.get_fp_positions(tgt) for tele in telescopes]

# create the dictionary containing the costs and constraints for all classes
# of targets
classdict = {}
classdict["sci_P1"] = {"nonObservationCost": 100,
                       "partialObservationCost": 1000, "calib": False}
classdict["sky"] = {"numRequired": 240,
                    "nonObservationCost": 1e9, "calib": True}
classdict["cal"] = {"numRequired": 40,
                    "nonObservationCost": 1e9, "calib": True}

tclassdict = {'sci_P1' : 1, 'sky' : 2, 'cal' : 3}

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

forbiddenPairs = []
for i in range(nvisit):
    forbiddenPairs.append([])

done = False
while not done:
    # compute observation strategy
    prob = nf.buildProblem(bench, tgt, tpos, classdict, t_obs,
                           vis_cost, cobraMoveCost=cobraMoveCost,
                           collision_distance=2., elbow_collisions=True,
                           gurobi=False, gurobiOptions=gurobiOptions,
                           alreadyObserved=alreadyObserved,
                           forbiddenPairs=forbiddenPairs)

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

    print("Checking for trajectory collisions")
    ncoll = 0
    for ivis, (vis, tp) in enumerate(zip(res, tpos)):
        selectedTargets = np.full(len(bench.cobras.centers), NULL_TARGET_POSITION)
        ids = np.full(len(bench.cobras.centers), NULL_TARGET_ID)
        for tidx, cidx in vis.items():
            selectedTargets[cidx] = tp[tidx]
            ids[cidx] = ""
        for i in range(selectedTargets.size):
            if selectedTargets[i] != NULL_TARGET_POSITION:
                dist = np.abs(selectedTargets[i]-bench.cobras.centers[i])

        simulator = CollisionSimulator2(bench, cobraCoach, TargetGroup(selectedTargets, ids))
        simulator.run()
        if np.any(simulator.endPointCollisions):
            print("ERROR: detected end point collision, which should be impossible")
        coll_tidx = []
        for tidx, cidx in vis.items():
            if simulator.collisions[cidx]:
                coll_tidx.append(tidx)
        ncoll += len(coll_tidx)
        for i1 in range(0,len(coll_tidx)):
            for i2 in range(i1+1,len(coll_tidx)):
                if np.abs(tp[coll_tidx[i1]]-tp[coll_tidx[i2]])<10:
                    forbiddenPairs[ivis].append((coll_tidx[i1],coll_tidx[i2]))

    print("trajectory collisions found:", ncoll)
    done = ncoll == 0

for vis, tp in zip(res, tpos):
    selectedTargets = np.full(len(bench.cobras.centers), NULL_TARGET_POSITION)
    ids = np.full(len(bench.cobras.centers), NULL_TARGET_ID)
    for tidx, cidx in vis.items():
        selectedTargets[cidx] = tp[tidx]
        ids[cidx] = ""
    for i in range(selectedTargets.size):
        if selectedTargets[i] != NULL_TARGET_POSITION:
            dist = np.abs(selectedTargets[i]-bench.cobras.centers[i])

    simulator = CollisionSimulator2(bench, cobraCoach, TargetGroup(selectedTargets, ids))
    simulator.run()
    simulator.plotResults(paintFootprints=False)
    plotUtils.pauseExecution()

# write output file
with open("output.txt", "w") as f:
    for i, (vis, tp, tel) in enumerate(zip(res, tpos, telescopes)):
        # Write legacy output.txt
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

        ets_fiber_assigner.io_helpers.writePfsDesign(
            pfsDesignId=i,
            pfsDesignDirectory='.',
            vis=vis,
            tp=tp,
            tel=tel,
            tgt=tgt,
            classdict=tclassdict)

