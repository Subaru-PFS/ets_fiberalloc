import ets_fiber_assigner.netflow as nf
import numpy as np

# make runs reproducible
np.random.seed(20)

# define locations of the input files
catalog_path = "data/"
fscience_targets = catalog_path+"ets_test_data.dat"
# So far, we only have test data for targets.
# Once we have files for calibration stars and sky locations, we can add them
# here.
#fcal_stars       = catalog_path+"/ets_test_data_fcstars.dat"
#fsky_pos         = catalog_path+"/ets_test_data_sky.dat"

# read all targets into a sigle list, giving them their proper types
tgt = nf.readScientificFromFile(fscience_targets, "sci")
# add calibration targets
#tgt += dm.readCalibrationFromFile(fcal_stars, "cal")
#tgt += dm.readCalibrationFromFile(fsky_pos, "sky")

# get a complete, idealized focal plane configuration
cobras = nf.getFullFocalPlane()

Rmax = 100.  # radius of the reduced focal plane

# build reduced Cobra list to speed up calculation
cobras = [c for c in cobras if abs(c.center) < Rmax]


# point the telescope at the center of all science targets
raTel, decTel = nf.telescopeRaDecFromFile(fscience_targets)
posang = 0.
otime = "2016-04-03T08:00:00Z"
telescopes = []

# number of distinct observations
nvisit = 3

# generate randomly jittered telescope pointings for every observation
for _ in range(nvisit):
    telescopes.append(nf.Telescope(cobras, 1., raTel+np.random.normal()*1e-2,
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
#classdict["sky"] = {"numRequired": 2,
#                    "nonObservationCost": 1000, "calib": True}
#classdict["cal"] = {"numRequired": 1,
#                    "nonObservationCost": 1000, "calib": True}

# optional: slightly increase the cost for later observations,
# to observe as early as possible
vis_cost = [0.1 + 0.1*i for i in range(nvisit)]


# optional: penalize assignments where the cobra has to move far out
def cobraMoveCost(dist):
    return 5.*dist

# compute observation strategy
res = nf.observeWithNetflow(telescopes[0].Cobras, tgt, tpos, classdict, 900.,
                            vis_cost, cobraMoveCost=cobraMoveCost,
                            collision_distance=2., gurobi=False)

# plot assignments for every observation
for vis, tp in zip(res, tpos):
    nf.plot_assignment(cobras, tgt, tp, vis)
