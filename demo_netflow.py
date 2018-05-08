from __future__ import print_function
import ets_fiber_assigner.netflow as nf
import numpy as np
from cobraOps.Bench import Bench
from cobraOps.TargetGroup import TargetGroup
from cobraOps.CobrasCalibrationProduct import CobrasCalibrationProduct
from cobraOps.CollisionSimulator import CollisionSimulator
from cobraOps.cobraConstants import NULL_TARGET_POSITION, NULL_TARGET_ID
from cobraOps import plotUtils
from collections import defaultdict

# make runs reproducible
np.random.seed(20)

# define locations of the input files
catalog_path = "data/"
fscience_targets = catalog_path+"pfs_preliminary_target_cosmology.dat"
# So far, we only have test data for targets.
# Once we have files for calibration stars and sky locations, we can add them
# here.
fcal_stars       = catalog_path+"pfs_preliminary_target_cosmology_fcstars.dat"
fsky_pos         = catalog_path+"pfs_preliminary_target_cosmology_sky.dat"

# read all targets into a single list, giving them their proper types
tgt = nf.readScientificFromFile(fscience_targets, "sci")
# add calibration targets
tgt += nf.readCalibrationFromFile(fcal_stars, "cal")
tgt += nf.readCalibrationFromFile(fsky_pos, "sky")

# get a complete, idealized focal plane configuration
bench = Bench(layout="full")
# if you have the XML file, you can also generate a more realistic focal plane
#bench = Bench(calibrationProduct=CobrasCalibrationProduct(
#    "../ics_cobraOps/python/ics/demos/updatedMaps6.xml"))

# point the telescope at the center of all science targets
raTel, decTel = nf.telescopeRaDecFromFile(fscience_targets)
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
classdict["sci_P1"] = {"nonObservationCost": 100,
                       "partialObservationCost": 1e9, "calib": False}
classdict["sci_P2"] = {"nonObservationCost": 90,
                       "partialObservationCost": 1e9, "calib": False}
classdict["sci_P3"] = {"nonObservationCost": 80,
                       "partialObservationCost": 1e9, "calib": False}
classdict["sci_P4"] = {"nonObservationCost": 70,
                       "partialObservationCost": 1e9, "calib": False}
classdict["sci_P5"] = {"nonObservationCost": 60,
                       "partialObservationCost": 1e9, "calib": False}
classdict["sci_P6"] = {"nonObservationCost": 50,
                       "partialObservationCost": 1e9, "calib": False}
classdict["sci_P7"] = {"nonObservationCost": 40,
                       "partialObservationCost": 1e9, "calib": False}
classdict["sky"] = {"numRequired": 20,
                    "nonObservationCost": 1000, "calib": True}
classdict["cal"] = {"numRequired": 5,
                    "nonObservationCost": 1000, "calib": True}

# optional: slightly increase the cost for later observations,
# to observe as early as possible
vis_cost = [0.1 + 0.1*i for i in range(nvisit)]


# optional: penalize assignments where the cobra has to move far out
def cobraMoveCost(dist):
    return 5.*dist


# duration of one observation in seconds
t_obs = 900.

# compute observation strategy
res = nf.observeWithNetflow(bench, tgt, tpos, classdict, t_obs,
                            vis_cost, cobraMoveCost=None,#cobraMoveCost,
                            collision_distance=2., elbow_collisions=True,
                            gurobi=True)


# write output file
with open("output.txt", "w") as f:
    for i, (vis, tp, tel) in enumerate(zip(res, tpos, telescopes)):
        print("exposure {}:".format(i))
        print("  assigned Cobras: {}".format(len(vis)))
        tdict = defaultdict(int)
        f.write("# Exposure {}: duration {}s, RA: {}, Dec: {}, PA: {}\n".
                format(i+1, t_obs, tel._ra, tel._dec, tel._posang))
        f.write("# Target    Fiber          X          Y         RA        DEC\n")
        for tidx, cidx in vis.items():
            tdict[tgt[tidx].targetclass] += 1
            f.write("{:} {:6d} {:10.5f} {:10.5f} {:10.5f} {:10.5f}\n"
                .format(tgt[tidx].ID, cidx+1, tp[tidx].real, tp[tidx].imag,
                        tgt[tidx].ra, tgt[tidx].dec))
        for cls, num in tdict.items():
            print("   {}: {}".format(cls, num))

for vis, tp in zip(res, tpos):
    selectedTargets = np.full(len(bench.cobras.centers), NULL_TARGET_POSITION)
    ids = np.full(len(bench.cobras.centers), NULL_TARGET_ID)
    for tidx, cidx in vis.items():
        selectedTargets[cidx] = tp[tidx]
        ids[cidx] = ""
    for i in range(selectedTargets.size):
        if selectedTargets[i] != NULL_TARGET_POSITION:
            dist = np.abs(selectedTargets[i]-bench.cobras.centers[i])

    simulator = CollisionSimulator(bench, TargetGroup(selectedTargets, ids))
    simulator.run()
    simulator.plotResults(paintFootprints=False)
    plotUtils.pauseExecution()

    # Animate the trajectory collisions
    (problematicCobras,) = np.where(np.logical_and(
        simulator.collisions, ~simulator.endPointCollisions))
    for cbr in problematicCobras:
        simulator.animateCobraTrajectory(cbr)
        plotUtils.pauseExecution()
