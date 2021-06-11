import numpy as np
from collections import OrderedDict
import argparse

import matplotlib
from astropy.table import Table, Column

from ics.cobraOps.Bench import Bench
from ics.cobraOps.RandomTargetSelector import RandomTargetSelector
from ics.cobraOps.TargetGroup import TargetGroup
from ics.cobraOps.CobrasCalibrationProduct import CobrasCalibrationProduct
from ics.cobraOps.CollisionSimulator import CollisionSimulator
from ics.cobraOps.cobraConstants import NULL_TARGET_POSITION, NULL_TARGET_ID, NULL_TARGET_INDEX
from pfs import datamodel

DIST_TYPES = OrderedDict()
DIST_TYPES["hom"] = "Homogeneous distribution."
DIST_TYPES["fib"] = "One per cobra."
DIST_TYPES["exp"] = "Exponential dropoff."


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


def dist_hom(args):
    """ implements homogeneous target distribution """
    s = args.radius * 2.
    size = int( args.rho0*s**2.)
    xx,yy = np.random.uniform(size=size) * s - s/2., np.random.uniform(size=size) * s - s/2.
    dsq = xx**2. + yy**2.
    ii = dsq <= (args.radius)**2.
    xx = xx[ii]
    yy = yy[ii]
    xx = xx / np.cos(np.deg2rad(args.dec)) + args.ra
    yy = yy + args.dec
    N = len(xx)
    ids = np.arange( N , dtype = int) + args.start_id
    tt = [args.obj_type] * N
    mm = np.random.uniform(size=N) * (args.mag_max - args.mag_min) + args.mag_min
    zz = np.random.uniform(size=N) * (args.redshift_max - args.redshift_min) + args.redshift_min
    pri = [args.prioritiy] * N
    ee = [args.exp_time] * N
    cc = [Column(ids, name="ID", dtype=int, format="010d"), \
            Column(xx, name="R.A.", unit="deg.", dtype=float, format="11.6f"),\
            Column(yy, name="Dec.", unit="deg.", dtype=float, format="11.6f"), \
            Column(pri, name="Priority", description="Priority class",dtype=str, format="10s"), \
            Column(ee,  name="Exposure Time", unit="sec.",dtype=float, format="6.1f"), \
            Column(mm, name="Magnitude", unit="AB mag", description="(g-band,HSC-CModel)", dtype=float, format="6.3f"), \
            Column(zz, name="Redshift", dtype=float, format="6.3f"), \
            Column(tt, name="Object Type", dtype=str, format="20s")]
    return Table( cc )

def filterTargets(args, bench):
    import ets_fiber_assigner.netflow as nf

    # compute focal plane positions
    tgt = nf.readScientificFromFile(args.out_file_name, 'sci')
    telescope = nf.Telescope(args.tel_ra, args.tel_dec, args.tel_pos_ang, args.tel_obs_time)
    tpos = telescope.get_fp_positions(tgt)
    print("read", tpos.shape[0], "targets")
    sel = RandomTargetSelector(bench,TargetGroup(tpos))
    sel.calculateAccessibleTargets()
    tgtidx = sel.accessibleTargetIndices.copy().flatten()
    tgtidx = tgtidx[tgtidx!=NULL_TARGET_INDEX]
    bc = np.bincount(tgtidx, minlength=tpos.shape[0])
    return (bc==1)
    
def write_to_pfsdesign(args, t):
    import ets_fiber_assigner.netflow as nf

    # compute focal plane positions
    tgt = nf.readScientificFromFile(args.out_file_name, 'sci')
    tel = nf.Telescope(args.tel_ra, args.tel_dec, args.tel_pos_ang, args.tel_obs_time)
    tpos = tel.get_fp_positions(tgt)
    pfiNominal = [[tp.real, tp.imag] for tp in tpos]

    N = t["R.A."].shape[0]
    d = dict(pfsDesignId = 0,
            raBoresight=tel._ra,
            decBoresight=tel._dec,
            posAng=tel._posang,
            fiberId=[0] * N,
            tract=[np.nan] * N,
            patch=["nan,np.nan"] * N,
            ra=t["R.A."],
            dec=t["Dec."],
            catId=[np.nan] * N,
            objId=np.arange(N),
            targetType=[1]*N,
            fiberStatus=[1] * N,
            fiberFlux=[[np.nan]] * N,
            psfFlux=[[np.nan]] * N,
            totalFlux=[[np.nan]] * N,
            fiberFluxErr=[[np.nan]] * N,
            psfFluxErr=[[np.nan]] * N,
            totalFluxErr=[[np.nan]] * N,
            filterNames=[['g']] * N,
            pfiNominal=pfiNominal,
            arms=None,
            guideStars=None
            )
    pfsDesign = datamodel.PfsDesign(**d)
    pfsDesign.write(dirName='.', fileName="pfsdesign_exp{:03d}.fits".format(0))


parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")

parser.add_argument("-a", "--ra", type=float, default=0.,
                    help="Central RA [Deg.].")
parser.add_argument("-d", "--dec", type=float, default=0.,
                    help="Central Dec [Deg.].")
parser.add_argument("-R", "--radius", type=float, default=0.7,
                    help="Maximum radius [Deg].")
parser.add_argument("-t", "--type", type=str, default=list(DIST_TYPES.keys())[0],
                    help="Type of target distribution [{}].".format(", ".join(DIST_TYPES.keys()) ))
parser.add_argument("-r", "--rho0", type=float, default=2394,
                    help="Density of targets per sq. deg. Refers to central density of exponential model.")
parser.add_argument("-s", "--start_id", type=int, default=1e6,
                    help="Start ID. All target will be numbered consecutively starting with start_id.")
parser.add_argument("-T", "--obj_type", type=str, default="cosmology",
                    help="Object type (will be put into \"Object Type\" column).")
parser.add_argument("-m", "--mag_min", type=float, default=22.86,
                    help="Minimun magnitude (default = 22.86).")
parser.add_argument("-M", "--mag_max", type=float, default=24.29,
                    help="Maxmimum magnitude (default = 24.29).")
parser.add_argument("-e", "--exp_time", type=float, default=900.,
                    help="Exposure time (will be put into \"Exposure\" column).")
parser.add_argument("-P", "--prioritiy", type=int, default=1,
                    help="Target priority class (will be put into \"Priority\" column).")
parser.add_argument("-z", "--redshift_min", type=float, default=1.,
                    help="Minimum redshift.")
parser.add_argument("-Z", "--redshift_max", type=float, default=3.,
                    help="Maximum redshift.")
parser.add_argument("-o", "--out_file_name", type=str, default="test.dat",
                    help="Output file name.")

parser.add_argument("--plot", 
                    help="Visualize targets.", action="store_true")
parser.add_argument("--tel_ra", type=float, default=None,
                    help="Only used for plotting. Telescop RA. Will default to \"--ra\" if not set [Deg.].")
parser.add_argument("--tel_dec", type=float, default=None,
                    help="Only used for plotting. Telescop Dec. Will default to \"--ra\" if not set [Deg.].")
parser.add_argument("--tel_pos_ang", type=float, default=0.,
                    help="Only used for plotting. Telescop position angle [Deg.].")
parser.add_argument("--tel_obs_time", type=str, default= "2016-04-03T08:00:00Z",
                    help="Only used for plotting. Telescope observation time [ \"2016-04-03T08:00:00Z\"].")
parser.add_argument("--select-single-cobra", 
                    help="only select targets that can be seen from exactly one Cobra", action="store_true")

args = parser.parse_args()

if not args.type in DIST_TYPES.keys():
  raise argparse.ArgumentTypeError(\
          'Unknow distribution type "-t", expected is either of: {}.'.format(", ".join(DIST_TYPES.keys()) ))
  sys.exit(1)

if args.tel_ra == None:
    print("Option tel_ra not set, setting to {:.6f}".format(args.ra))
    args.tel_ra = args.ra
if args.tel_dec == None:
    print("Option tel_dec not set, setting to {:.6f}".format(args.dec))
    args.tel_dec = args.dec

cobraCoach, bench = getBench()

if args.type == "hom":
    t = dist_hom(args)
    t.write(args.out_file_name, format="ascii.ecsv", overwrite=True)

if args.select_single_cobra:
    # select targets that are accessible from exactly one Cobra
    selection = filterTargets(args, bench)
    t = t[selection]
    t.write(args.out_file_name, format="ascii.ecsv", overwrite=True)

#write_to_pfsdesign(args, t)

if args.plot:
    from matplotlib import pyplot as plt
    import ets_fiber_assigner.netflow as nf

    # compute focal plane positions
    tgt = nf.readScientificFromFile(args.out_file_name, 'sci')
    telescope = nf.Telescope(args.tel_ra, args.tel_dec, args.tel_pos_ang, args.tel_obs_time)
    tpos = telescope.get_fp_positions(tgt)

    f = plt.figure(figsize=[10,5])
    ax1 = plt.subplot(1,2,1)
    plt.plot(t["R.A."], t["Dec."], 'kx',markersize=.5, alpha=0.8)
    plt.ylabel('RA [Deg.]')
    plt.xlabel('Dec [Deg.]')
    plt.axis('equal')

    ax2 = plt.subplot(1,2,2)
    plt.plot(np.real(bench.cobras.centers), np.imag(bench.cobras.centers), 'bo', markersize=1, alpha=0.8)
    plt.plot(np.real( tpos ), np.imag( tpos ), 'kx',markersize=.5, alpha=0.8)
    plt.ylabel('y [mm]')
    plt.xlabel('x [mm]')
    plt.axis('equal')

    f.tight_layout()
    plt.show()

