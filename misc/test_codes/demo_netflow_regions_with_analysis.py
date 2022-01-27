import ets_fiber_assigner.netflow as nf
import numpy as np
from ics.cobraOps.Bench import Bench
from ics.cobraOps.TargetGroup import TargetGroup
from ics.cobraOps.CobrasCalibrationProduct import CobrasCalibrationProduct
from ics.cobraOps.CollisionSimulator import CollisionSimulator
from ics.cobraOps.cobraConstants import NULL_TARGET_POSITION, NULL_TARGET_ID
from ics.cobraOps import plotUtils
from collections import defaultdict

# make runs reproducible
np.random.seed(20)

# define locations of the input files
catalog_path = "data/"
fsky_pos = catalog_path+"pfs_preliminary_target_cosmology_sky.dat"


tgt = nf.readCalibrationFromFile(fsky_pos, "sky")

# get a complete, idealized focal plane configuration
bench = Bench(layout="full")
# if you have the XML file, you can also generate a more realistic focal plane
# bench = Bench(calibrationProduct=CobrasCalibrationProduct(
#     "../ics_cobraOps/python/ics/demos/updatedMaps6.xml"))

# point the telescope at the center of all science targets
raTel, decTel = nf.telescopeRaDecFromFile(fsky_pos)
posang = 0.
otime = "2016-04-03T08:00:00Z"
telescopes = []



#########################################################################
# code to generate target and sky target distribution for
# testing of focal plane regions and instrumental regions.
def save_targets(filename, tgtra, tgtdec):
    s = """
#  1  ID
#  2  RA            [deg.]
#  3  DEC           [deg.]
#  4  EXP_TIME      [sec.]
#  5  Priority      [1(highest)-15(lowest)]
#  6  Magnitude     [AB mag] (g-band,HSC-CModel)
#  7  Redshift
#  8  ObjectType
    """
    for i, (a,d) in enumerate(zip(tgtra, tgtdec)):
        s += f"K{i:06d} {a:.5f}  {d:.5f}   900.0   1  99.99  N/A  sky(cosmology)\n"
    with open(filename, 'w') as fout:
        fout.write(s)

fsci_pos = catalog_path+"mock_uniform_sci.dat"
fsky_pos = catalog_path+"mock_gradient_sky.dat"

# generate sky targets with strong east west gradient in density
rr = np.linspace( raTel-.5, raTel+.5, 1000)
pp = rr*100 - np.min(rr*100) + 1 # setup propability distribution, 100 times more
                               # more targets on the right
pp=pp/np.sum(pp) # normalize

M,N = 10000,10000
tgtra  = np.random.choice(rr, p=pp, size=M)
tgtdec = np.random.uniform(low=decTel-.5,high=decTel+.5, size=M)
save_targets(fsky_pos, tgtra, tgtdec)

# generate sci objects with uniform distribution
tgtra = np.random.uniform(low=raTel-.5,high=raTel+.5, size=N)
tgtdec = np.random.uniform(low=decTel-.5,high=decTel+.5, size=N)
save_targets(fsci_pos, tgtra, tgtdec)

tgt_sci = nf.readScientificFromFile(fsci_pos, "sci")
tgt_sky = nf.readCalibrationFromFile(fsky_pos, "sky")

from matplotlib import pyplot as plt


xx_sci = [ t.ra for t in tgt_sci ]
yy_sci = [ t.dec for t in tgt_sci ]
xx_sky = [ t.ra for t in tgt_sky ]
yy_sky = [ t.dec for t in tgt_sky ]

f = plt.figure()
plt.plot(xx_sci, yy_sci, '.', label='sci', alpha=.3)
plt.plot(xx_sky, yy_sky, '.', label='sky', alpha=.3)
plt.xlabel("ra [deg]")
plt.ylabel("dec [deg]")
plt.legend()
f.tight_layout()


tgt = tgt_sci + tgt_sky
#tgt = tgt_sci 
#tgt = tgt_sky

#1/0
##########################################################################


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
classdict["sci_P1"] = {"numRequired" : np.inf, "nonObservationCost": 100,\
                    "partialObservationCost": 1000, "calib": False}


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


# define two location regions arbitrarily
ncobras = bench.cobras.nCobras
cobraRegions=np.zeros(ncobras,dtype=np.int32)
cobraRegions[ncobras//2:] = 1

##########################################################################
# code to generate focal plane region definitions
if False:
# mF: bin focal plan into grid and assign corresponding region ids
    _ = np.imag(bench.cobras.centers); miny=np.min(_); maxy=np.max(_)
    _ = np.real(bench.cobras.centers); minx=np.min(_); maxx=np.max(_)
    N = 7 # bin size in mm
    minSkyTargetsPerLocation = 10
    xx=np.linspace(minx,maxx,N)
    yy=np.linspace(minx,maxy,N)
    rid = 0
    for i in range(len(xx)-1):
        ii  = np.imag(bench.cobras.centers) >= xx[i]
        ii *= np.imag(bench.cobras.centers) <  xx[i+1]
        for j in range(len(yy)-1):
            jj  = np.real(bench.cobras.centers) >= yy[j]
            jj *= np.real(bench.cobras.centers)  < yy[j+1]
            jj *= ii
            cobraRegions[ii * jj] = rid
            rid += 1
if True:
# voronoi based binning
# requires: pip install vorbin
    from vorbin import voronoi_2d_binning
    M = 75 # bin size in no of fibers
    minSkyTargetsPerLocation = 7

    x = np.real(bench.cobras.centers)
    y = np.imag(bench.cobras.centers)
    signal = np.ones_like(x)
    noise = signal
    targetSN = np.sqrt(M)
    binNum, xBin, yBin, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning.voronoi_2d_binning(x, y, signal, noise, targetSN)
    cobraRegions = binNum

    
    #from matplotlib import pyplot as plt
    f = plt.figure()
    plt.ion()
    for rid in np.unique(cobraRegions):
        ii = cobraRegions == rid
        plt.plot(np.real(bench.cobras.centers[ii]), np.imag(bench.cobras.centers[ii]), '.')
    plt.show()

# switch off region based distribution of calibrators
if False:
    cobraRegions[:] = 0

#1/0
##########################################################################

##########################################################################
# code to generate instrumental region definitions


minSkyTargetsPerInstrumentRegion = 10

from astropy.io import ascii
from cycler import cycler
import numpy as np

# get colormap
cmap=plt.cm.turbo
cmap=plt.cm.nipy_spectral
# cobraInstrumentRegion cycler with 5 equally spaced colors from that colormap
c = cycler('color', cmap(np.linspace(0,1,20)) )
# supply cycler to the rcParam
plt.rcParams["axes.prop_cycle"] = c

t = ascii.read("/Users/mxhf/ownCloudRZG/work/MPE/pfs/src/pfs_utils/data/fiberids/grandfibermap.20210314.txt")

# filter for actual science fibers on cobras
ii = t['\cob'] != '-'
t = t[ii]


spectrographs = np.unique( t['sp'] )

instrumRegions = np.zeros(ncobras,dtype=np.int32)

reg_count = 0

f = plt.figure()

cc = np.array([int(c)-1 for c in t['\cob']]) # cobra number for each fiber in gand table


for i,s in enumerate(spectrographs):
    print("spectrograph: ", s)

    jj = t['sp'] == s # all fibers that go to this spectrograph
    ff = t[jj]['fh']  # fiber position in slit

    for j in range(5):
        ii = ff >= j*120 # define five regions a 120 fibers along slit
        ii*= ff < (j+1)*120

        instrumRegions[cc[jj][ii]] = reg_count

        reg_count += 1

        plt.plot(ff[ii], np.zeros_like(ff[ii]) + i ,'o')
        #print(list( ff[ii] ))

plt.yticks([])
plt.xlabel("fiber number")


f = plt.figure()


cobra_x = np.real(bench.cobras.centers)
cobra_y = np.imag(bench.cobras.centers)
rr = np.unique(instrumRegions)
for r in rr:
    ii = instrumRegions == r
    plt.plot(cobra_x[ii], cobra_y[ii], '.')

print(s)
plt.xlabel('x [mm]')
plt.ylabel('y [mm]')
f.tight_layout()


cobraInstrumentRegion = instrumRegions

#1/0

##########################################################################



done = False
while not done:
    # compute observation strategy
    if False:
        print("Constraining distribution accross instrumental regions.")
        prob = nf.buildProblem(bench, tgt, tpos, classdict, t_obs,
                               vis_cost, cobraMoveCost=cobraMoveCost,
                               collision_distance=2., elbow_collisions=True,
                               gurobi=True, gurobiOptions=gurobiOptions,
                               alreadyObserved=alreadyObserved,
                               forbiddenPairs=forbiddenPairs,
                               cobraLocationGroup=cobraRegions,
                               minSkyTargetsPerLocation=minSkyTargetsPerLocation,
                               locationGroupPenalty=1e9,
                               cobraInstrumentRegion=cobraInstrumentRegion,
                               minSkyTargetsPerInstrumentRegion=minSkyTargetsPerInstrumentRegion,
                               instrumentRegionPenalty=1e9)
    else:
        print("NOT constraining distribution accross instrumental regions.")
        prob = nf.buildProblem(bench, tgt, tpos, classdict, t_obs,
                               vis_cost, cobraMoveCost=cobraMoveCost,
                               collision_distance=2., elbow_collisions=True,
                               gurobi=True, gurobiOptions=gurobiOptions,
                               alreadyObserved=alreadyObserved,
                               forbiddenPairs=forbiddenPairs,
                               cobraLocationGroup=cobraRegions,
                               minSkyTargetsPerLocation=minSkyTargetsPerLocation,
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

        simulator = CollisionSimulator(bench, TargetGroup(selectedTargets, ids))
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


###########################################################################
# analyse coverage of focal plane regions
# with calibrators

regionCoverage = {}
for i, (vis, tp, tel) in enumerate(zip(res, tpos, telescopes)):
    regionCoverage[i] = np.zeros(np.max(cobraRegions)+1, dtype=int )

    for tidx, cidx in vis.items():
        if tgt[tidx].targetclass == 'sky':
            regionCoverage[i] [cobraRegions[cidx]] +=1


N = len(tpos)
f=plt.figure(figsize=[10,10])
for i, (vis, tp, tel) in enumerate(zip(res, tpos, telescopes)):

    xx_sky = []
    yy_sky = []
    for tidx, cidx in vis.items():

        if tgt[tidx].targetclass == 'sky':
            xx_sky += [tp[tidx].real]
            yy_sky += [tp[tidx].imag]


    plt.subplot( int(np.ceil(np.sqrt(N))), int(np.ceil(np.sqrt(N))), i+1 )
    plt.plot(xx_sky, yy_sky, '.')
    plt.xlabel("x [mm]")
    plt.ylabel("y [mm]")

f.tight_layout()

def entropy(a):
    return - np.sum(a * np.log(a))

a = np.array( [regionCoverage[k] for k in regionCoverage] )
a = a.flatten()


f = plt.figure()
plt.hist(a)
plt.ylabel("number of regions")
plt.xlabel("number of calibrators per region")
plt.xlim([0,25])
plt.ylim([0,130])
f.tight_layout()

mean_cals = np.mean(a)
min_cals = np.min(a)
max_cals = np.max(a)
std_cals = np.std(a)
entropy_cals = entropy(a)



# compute lowest possible entropy
# (=all one, expect one, that contains rest)
b = np.zeros_like(a)
b[0] = np.sum(a) - (len(a)-1)
b[1:] = 1
min_entropy = entropy(b)

# compute largest possible entropy (= all equal) 
c = np.ones_like(a) * mean_cals
max_entropy = entropy(c)

entropy_poisson = np.nanmean( [entropy( np.random.poisson(mean_cals,size=len(a)) ) for i in range(1000) ] )

nentropy_cals = (entropy_cals-min_entropy)/(max_entropy-min_entropy)


print(f"Mean number of calibrators in focal plane regions:   {mean_cals:.2f}")
print(f"Min number of calibrators in focal plane regions:   {min_cals:.2f}")
print(f"Max number of calibrators in focal plane regions:   {max_cals:.2f}")
print(f"STD of number of calibrators in focal plane regions: {std_cals:.2f}")
print(f"Entropy:                    {entropy_cals:.1f}")
print(f"Norm. entropy:              {nentropy_cals:.3f}")
print(f"Entropy of all in one:      {min_entropy:.1f}")
print(f"Entropy of all equal:       {max_entropy:.1f}")
print(f"Entropy of poisson process: {entropy_poisson:.1f}")

#############################################################################

#############################################################################
# analyse coverage of instrumental regions
instRegionCoverage = {}
for i, (vis, tp, tel) in enumerate(zip(res, tpos, telescopes)):
    instRegionCoverage[i] = np.zeros(np.max(instrumRegions)+1, dtype=int )

    for tidx, cidx in vis.items():
        if tgt[tidx].targetclass == 'sky':
            instRegionCoverage[i] [instrumRegions[cidx]] +=1

a = np.array( [instRegionCoverage[k] for k in instRegionCoverage] )
a = a.flatten()

inst_mean_cals = np.mean(a)
inst_min_cals = np.min(a)
inst_max_cals = np.max(a)
inst_std_cals = np.std(a)
inst_entropy_cals = entropy(a)

# compute lowest possible entropy
# (=all one, expect one, that contains rest)
b = np.zeros_like(a)
b[0] = np.sum(a) - (len(a)-1)
b[1:] = 1
inst_min_entropy = entropy(b)

# compute largest possible entropy (= all equal) 
c = np.ones_like(a) * inst_mean_cals
inst_max_entropy = entropy(c)

# compute entropy from a poisson process
inst_entropy_poisson = np.nanmean( [entropy( np.random.poisson(inst_mean_cals,size=len(a)) ) for i in range(1000) ] )
inst_nentropy_cals = (inst_entropy_cals - inst_min_entropy)/(inst_max_entropy - inst_min_entropy)



fiber_slit_position = np.zeros(ncobras,dtype=np.int32)
fiber_spectrograph   = np.zeros(ncobras,dtype=np.int32)
for c,s,h in zip(cc, t['sp'], t['fh']):
    #print(f"cobra {c} spectrograph {s} fiberhole {h}")
    fiber_slit_position[c] = h
    fiber_spectrograph[ c ] = s

f = plt.figure()
plt.hist(a)
plt.ylabel("number of regions")
plt.xlabel("number of calibrators per instrum region")
plt.xlim([0,50])
plt.ylim([0,100])
f.tight_layout()

f = plt.figure()

N = len(tpos)
ax = None
for i, (vis, tp, tel) in enumerate(zip(res, tpos, telescopes)):
    ax = plt.subplot(N,1,i+1, sharex=ax)
    plt.plot(fiber_slit_position, fiber_spectrograph, '.' , alpha=.5)
    for tidx, cidx in vis.items():
        if tgt[tidx].targetclass == 'sky':
            plt.plot([fiber_slit_position[cidx]],[fiber_spectrograph[cidx]],'bx')
    if i == N-1:
        plt.xlabel('fibernummer')
    else:
        plt.xticks([])
    plt.ylabel('SU')
f.tight_layout()

print(f"Mean number of calibrators in instrumental regions:   {inst_mean_cals:.2f}")
print(f"Min number of calibrators in instrumental regions:   {inst_min_cals:.2f}")
print(f"Max number of calibrators in instrumental regions:   {inst_max_cals:.2f}")
print(f"STD of number of calibrators in instrumental regions: {inst_std_cals:.2f}")
print(f"Entropy:                    {inst_entropy_cals:.1f}")
print(f"Norm. entropy:              {inst_nentropy_cals:.3f}")
print(f"Entropy of all in one:      {inst_min_entropy:.1f}")
print(f"Entropy of all equal:       {inst_max_entropy:.1f}")
print(f"Entropy of poisson process: {inst_entropy_poisson:.1f}")
#############################################################################

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

1/0

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
