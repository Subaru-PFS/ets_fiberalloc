import numpy as np
from collections import OrderedDict
import argparse

import matplotlib
from astropy.table import Table, Column

from ics.cobraOps.Bench import Bench
from ics.cobraOps.CobrasCalibrationProduct import CobrasCalibrationProduct
from ics.cobraOps.CollisionSimulator import CollisionSimulator

DIST_TYPES = OrderedDict()
DIST_TYPES["hom"] = "Homogeneous distribution."
DIST_TYPES["fib"] = "One per cobra."
DIST_TYPES["exp"] = "Exponential dropoff."

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
parser.add_argument("-x", "--bench_xml", type=str, default="",
                    help="Cobra bench XML.")

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

args = parser.parse_args()

if not args.type in DIST_TYPES.keys():
  raise argparse.ArgumentTypeError(\
          'Unknow distribution type "-t", expected is either of: {}.'.format(", ".join(DIST_TYPES.keys()) ))
  sys.exit(1)

if args.type == "hom":
    t = dist_hom(args)
    t.write(args.out_file_name, format="ascii.ecsv", overwrite=True)


if args.plot:
    from matplotlib import pyplot as plt
    import ets_fiber_assigner.netflow as nf

    if args.tel_ra == None:
        print("Option tel_ra not set, setting to {:.6f}".format(args.ra))
        args.tel_ra = args.ra
    if args.tel_dec == None:
        print("Option tel_dec not set, setting to {:.6f}".format(args.dec))
        args.tel_dec = args.dec

    if args.bench_xml == "":
        print("WARNING: Using default bench for visualisation.")
        bench = Bench(layout="full")
    else:
        bench = Bench(calibrationProduct=CobrasCalibrationProduct( args.bench_xml ))

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

