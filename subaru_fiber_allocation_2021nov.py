# Script for commissioning runs 09/2021

# Necessary preparations for running:
#
# This script depends on several other modules from https://github.com/Subaru-PFS
# All of them were at the HEAD of the respective master branches, with the
# exception of "ets_fiber_assigner" (must be on branch "commissioning_demo").
#
# Also the "pulp" Python package (version 1.6!) is required to solve the fiber assignment
# problem.
#
# Also, the environment variable PFS_INSTDATA_DIR must be set correctly.

import argparse
import configparser
import os
import tempfile
import time

import ets_fiber_assigner.netflow as nf
import matplotlib.path as mppath
import numpy as np
import pandas as pd
import pfs.datamodel
import psycopg2
import psycopg2.extras
from astropy import units as u
from astropy.table import Table
from astropy.time import Time
from ets_shuffle import query_utils
from ets_shuffle.convenience import flag_close_pairs
from ets_shuffle.convenience import guidecam_geometry
from ets_shuffle.convenience import update_coords_for_proper_motion
from ics.cobraOps.Bench import Bench
from ics.cobraOps.BlackDotsCalibrationProduct import BlackDotsCalibrationProduct
from ics.cobraOps.cobraConstants import NULL_TARGET_ID
from ics.cobraOps.cobraConstants import NULL_TARGET_POSITION
from ics.cobraOps.CollisionSimulator2 import CollisionSimulator2
from ics.cobraOps.TargetGroup import TargetGroup
from pfs.utils.coordinates.CoordTransp import CoordinateTransform as ctrans
from pfs.utils.coordinates.CoordTransp import ag_pfimm_to_pixel
from pfs.utils.pfsDesignUtils import makePfsDesign
from procedures.moduleTest.cobraCoach import CobraCoach
from targetdb import targetdb


# This was needed for fixing some issues with the XML files.
# Can probably be simplified. Javier?
def getBench(args):

    os.environ["PFS_INSTDATA_DIR"] = args.pfs_instdata_dir
    cobraCoach = CobraCoach(
        "fpga", loadModel=False, trajectoryMode=True, rootDir=args.cobra_coach_dir
    )

    cobraCoach.loadModel(version="ALL", moduleVersion=args.cobra_coach_module_version)

    # Get the calibration product
    calibrationProduct = cobraCoach.calibModel

    # Set some dummy center positions and phi angles for those cobras that have
    # zero centers
    zeroCenters = calibrationProduct.centers == 0
    calibrationProduct.centers[zeroCenters] = np.arange(np.sum(zeroCenters)) * 300j
    calibrationProduct.phiIn[zeroCenters] = -np.pi
    calibrationProduct.phiOut[zeroCenters] = 0
    print("Cobras with zero centers: %i" % np.sum(zeroCenters))

    # Use the median value link lengths in those cobras with zero link lengths
    zeroLinkLengths = np.logical_or(
        calibrationProduct.L1 == 0, calibrationProduct.L2 == 0
    )
    calibrationProduct.L1[zeroLinkLengths] = np.median(
        calibrationProduct.L1[~zeroLinkLengths]
    )
    calibrationProduct.L2[zeroLinkLengths] = np.median(
        calibrationProduct.L2[~zeroLinkLengths]
    )
    print("Cobras with zero link lenghts: %i" % np.sum(zeroLinkLengths))

    # Use the median value link lengths in those cobras with too long link lengths
    tooLongLinkLengths = np.logical_or(
        calibrationProduct.L1 > 100, calibrationProduct.L2 > 100
    )
    calibrationProduct.L1[tooLongLinkLengths] = np.median(
        calibrationProduct.L1[~tooLongLinkLengths]
    )
    calibrationProduct.L2[tooLongLinkLengths] = np.median(
        calibrationProduct.L2[~tooLongLinkLengths]
    )
    print("Cobras with too long link lenghts: %i" % np.sum(tooLongLinkLengths))

    calibrationFileName = os.path.join(
        os.environ["PFS_INSTDATA_DIR"], "data/pfi/dot", "black_dots_mm.csv"
    )
    blackDotsCalibrationProduct = BlackDotsCalibrationProduct(calibrationFileName)

    # Create the bench instance
    bench = Bench(
        layout="calibration",
        calibrationProduct=calibrationProduct,
        blackDotsCalibrationProduct=blackDotsCalibrationProduct,
    )
    print("Number of cobras:", bench.cobras.nCobras)

    return cobraCoach, bench


def get_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--ra",
        type=float,
        default=0.0,
        help="Telescope center RA [degrees] (default: 0.0)",
    )
    parser.add_argument(
        "--dec",
        type=float,
        default=0.0,
        help="Telescope center Dec [degrees] (default: 0.0)",
    )
    parser.add_argument(
        "--pa",
        type=float,
        default=0.0,
        help="Telescope position angle [degrees] (default: 0.0)",
    )
    parser.add_argument(
        "--observation_time",
        type=str,
        default="2021-11-20T15:00:00Z",
        help="planned time of observation in UTC (default: 2021-11-20 15:00:00)",
    )
    parser.add_argument(
        "--lim_target_mag",
        type=float,
        default="19.",
        help="magnitude of the faintest targets (obsolete) (default:19)",
    )

    parser.add_argument(
        "--design_dir",
        type=str,
        default=".",
        help="directory for storing PFS designs (default: .)",
    )

    parser.add_argument(
        "--guidestar_mag_max",
        type=float,
        default=19.0,
        help="maximum magnitude for guide star candidates (default: 19.)",
    )
    parser.add_argument(
        "--guidestar_neighbor_mag_min",
        type=float,
        default=21.0,
        help="minimum magnitude for objects in the vicinity of guide star candidates (default: 21.)",
    )
    parser.add_argument(
        "--guidestar_minsep_deg",
        type=float,
        default=1.0 / 3600,
        help="radius of guide star candidate vicinity (default: 1/3600)",
    )

    parser.add_argument(
        "--use_gurobi",
        action="store_true",
        help="use Gurobi",
    )
    parser.add_argument(
        "--cobra_coach_dir",
        type=str,
        default=".",
        help="path for temporary cobraCoach files (default: .)",
    )

    parser.add_argument(
        "--cobra_coach_module_version",
        type=str,
        default="final_20210920_mm",
        help="version of the bench decription file (default: final_20210920_mm)",
    )

    parser.add_argument(
        "--targetdb_conf",
        type=str,
        default="targetdb_config.ini",
        help="Config file for targetDB (default: targetdb_config.ini)",
    )
    parser.add_argument(
        "--gaiadb_conf",
        type=str,
        default="gaiadb_config_hilo.ini",
        help="Config file for Subaru's Gaia DB (default: gaiadb_config_hilo.ini",
    )
    parser.add_argument(
        "--target_mag_max",
        type=float,
        default=19.0,
        help="Maximum (faintest) magnitude for stars in fibers (default: 19.)",
    )
    parser.add_argument(
        "--target_mag_min",
        type=float,
        default=0.0,
        help="Minimum (brightest) magnitude for stars in fibers (default: 0)",
    )
    parser.add_argument(
        "--target_mag_filter",
        type=str,
        default="g",
        help="Photometric band (grizyj of PS1) to apply magnitude cuts (default: g)",
    )
    parser.add_argument(
        "--fluxstd_min_prob_f_star",
        type=float,
        default=0.0,
        help="Minimum acceptable prob_f_star (default: 0)",
    )
    parser.add_argument(
        "--telescope_elevation",
        type=float,
        default=60.0,
        help="Telescope elevation in degree (default: 60)",
    )
    parser.add_argument(
        "--n_fluxstd",
        type=int,
        default=50,
        help="Number of FLUXSTD stars to be allocated. (default: 50)",
    )
    parser.add_argument(
        "--pfs_instdata_dir",
        type=str,
        default="/Users/monodera/Dropbox/NAOJ/PFS/Subaru-PFS/pfs_instdata/",
        help="Location of pfs_instdata (default: /Users/monodera/Dropbox/NAOJ/PFS/Subaru-PFS/pfs_instdata/)",
    )

    args = parser.parse_args()

    if args.observation_time.lower() == "now":
        print("converting to the current time")
        args.observation_time = (
            Time.now().iso
        )  # astropy.time.Time.now() uses datetime.utcnow()

    return args


def connect_subaru_gaiadb(conf=None):
    config = configparser.ConfigParser()
    config.read(conf)
    conn = psycopg2.connect(**dict(config["dbinfo"]))
    return conn


def gen_target_list_from_targetdb(args):
    def connect_db(conf=None):
        config = configparser.ConfigParser()
        config.read(conf)
        db = targetdb.TargetDB(**dict(config["dbinfo"]))
        db.connect()
        return db

    def generate_query_simple_boxsearch(
        ra1, ra2, dec1, dec2, mag_min, mag_max, mag_filter, min_prob_f_star
    ):
        # FIXME: I know this is too simple and stupid,
        #        but should be enough for the November 2021 commissioning run.
        query_target = """SELECT
    obj_id,
    ra,
    dec,
    epoch,
    priority,
    effective_exptime,
    psf_flux_g,
    psf_flux_r,
    psf_flux_i,
    psf_flux_z,
    psf_flux_y,
    target_type_id,
    input_catalog_id
    FROM target
    WHERE ra >= {:f} AND ra < {:f}
    AND dec >= {:f} AND dec < {:f}
    AND psf_mag_{:s} BETWEEN {:f} AND {:f}
    AND prob_f_star > {:f};
    """.format(
            ra1,
            ra2,
            dec1,
            dec2,
            mag_filter,
            mag_min,
            mag_max,
            min_prob_f_star,
        )
        return query_target

    db = connect_db(args.targetdb_conf)

    fp_rad_deg = 260.0 * 10.2 / 3600
    fp_fudge_factor = 1.5

    dw = fp_rad_deg * fp_fudge_factor

    cos_term = 1.0 / np.cos(args.dec * u.deg)

    dw_ra = dw * cos_term

    dec1, dec2 = args.dec - dw, args.dec + dw

    if args.ra - dw_ra < 0.0:
        ra1, ra2 = 0.0, args.ra + dw_ra
        q1 = generate_query_simple_boxsearch(
            ra1,
            ra2,
            dec1,
            dec2,
            args.target_mag_min,
            args.target_mag_max,
            args.target_mag_filter,
            args.args.fluxstd_min_prob_f_star,
        )
        ra1, ra2 = args.ra - dw_ra + 360.0, 360.0
        q2 = generate_query_simple_boxsearch(
            ra1,
            ra2,
            dec1,
            dec2,
            args.target_mag_min,
            args.target_mag_max,
            args.target_mag_filter,
            args.fluxstd_min_prob_f_star,
        )
        qlist = [q1, q2]
    elif args.ra + dw_ra >= 360.0:
        ra1, ra2 = 0.0, args.ra + dw_ra - 360.0
        q1 = generate_query_simple_boxsearch(
            ra1,
            ra2,
            dec1,
            dec2,
            args.target_mag_min,
            args.target_mag_max,
            args.target_mag_filter,
            args.fluxstd_min_prob_f_star,
        )
        ra1, ra2 = args.ra - dw_ra, 360.0
        q2 = generate_query_simple_boxsearch(
            ra1,
            ra2,
            dec1,
            dec2,
            args.target_mag_min,
            args.target_mag_max,
            args.target_mag_filter,
            args.fluxstd_min_prob_f_star,
        )
        qlist = [q1, q2]
    else:
        ra1, ra2 = args.ra - dw_ra, args.ra + dw_ra
        q1 = generate_query_simple_boxsearch(
            ra1,
            ra2,
            dec1,
            dec2,
            args.target_mag_min,
            args.target_mag_max,
            args.target_mag_filter,
            args.fluxstd_min_prob_f_star,
        )
        qlist = [q1]

    df = pd.DataFrame(
        columns=[
            "obj_id",
            "ra",
            "dec",
            "epoch",
            "priority",
            "effective_exptime",
            "psf_flux_g",
            "psf_flux_r",
            "psf_flux_i",
            "psf_flux_z",
            "psf_flux_y",
            "target_type_id",
            "input_catalog_id",
        ]
    )

    for q in qlist:
        print(q)
        t_begin = time.time()
        df_tmp = db.fetch_query(q)
        t_end = time.time()
        print("Time spent for querying: {:f}".format(t_end - t_begin))
        df = df.append(df_tmp, ignore_index=True)

    print(df)

    tbl_tmp = Table.from_pandas(df)

    tbl = Table()
    tbl["ID"] = np.array(tbl_tmp["obj_id"], dtype=np.int64)
    tbl["R.A."] = tbl_tmp["ra"]
    tbl["Dec."] = tbl_tmp["dec"]
    tbl["Epoch"] = tbl_tmp["epoch"]
    tbl["Exposure Time"] = tbl_tmp["effective_exptime"]
    tbl["Priority"] = np.array(tbl_tmp["priority"], dtype=int)

    # FIXME: I think it is worth putting the table file in a non-tmp directory
    with tempfile.NamedTemporaryFile(dir="/tmp", delete=False) as tmpfile:
        outfile = tmpfile.name
    tbl.write(outfile, format="ascii.ecsv", overwrite=True)

    tbl["psfFlux"] = [
        np.array(
            [
                tbl_tmp["psf_flux_g"][i],
                tbl_tmp["psf_flux_r"][i],
                tbl_tmp["psf_flux_i"][i],
            ]
        )
        for i in range(len(tbl["ID"]))
    ]
    tbl["filterNames"] = [["g_ps1", "r_ps1", "i_ps1"]] * len(tbl["ID"])
    tbl["target_type_id"] = tbl_tmp["target_type_id"]
    tbl["input_catalog_id"] = tbl_tmp["input_catalog_id"]

    db.close()

    return outfile, tbl


def gen_target_list_from_gaiadb(args):

    fp_rad_deg = 260.0 * 10.2 / 3600
    fp_fudge_factor = 1.2

    conn = connect_subaru_gaiadb(args.gaiadb_conf)
    cur = conn.cursor()

    query_string = """SELECT source_id,ref_epoch,ra,dec,pmra,pmdec,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag
    FROM gaia
    WHERE q3c_radial_query(ra, dec, {:}, {:}, {:})
    AND {:s} BETWEEN {:} AND {:}
    ;
    """.format(
        args.ra,
        args.dec,
        fp_rad_deg * fp_fudge_factor,
        "phot_g_mean_mag",
        args.target_mag_min,
        args.target_mag_max,
    )

    cur.execute(query_string)

    df_res = pd.DataFrame(
        cur.fetchall(),
        columns=[
            "source_id",
            "ref_epoch",
            "ra",
            "dec",
            "pmra",
            "pmdec",
            "phot_g_mean_mag",
            "phot_bp_mean_mag",
            "phot_rp_mean_mag",
        ],
    )

    cur.close()
    conn.close()

    tbl_tmp = Table.from_pandas(df_res)

    # ZPs are taken from Weiler (2018, A&A, 617, A138)
    tbl_tmp["g_mag_ab"] = (tbl_tmp["phot_g_mean_mag"] + (25.7455 - 25.6409)) * u.ABmag
    tbl_tmp["bp_mag_ab"] = (tbl_tmp["phot_bp_mean_mag"] + (25.3603 - 25.3423)) * u.ABmag
    tbl_tmp["rp_mag_ab"] = (tbl_tmp["phot_rp_mean_mag"] + (25.1185 - 24.7600)) * u.ABmag

    tbl_tmp["g_flux_njy"] = tbl_tmp["g_mag_ab"].to("nJy")
    tbl_tmp["bp_flux_njy"] = tbl_tmp["bp_mag_ab"].to("nJy")
    tbl_tmp["rp_flux_njy"] = tbl_tmp["rp_mag_ab"].to("nJy")

    n_target = tbl_tmp["source_id"].size

    tbl = Table()
    tbl["ID"] = tbl_tmp["source_id"]
    tbl["R.A."] = tbl_tmp["ra"]
    tbl["Dec."] = tbl_tmp["dec"]
    tbl["Epoch"] = tbl_tmp["ref_epoch"]
    tbl["Exposure Time"] = np.full(n_target, 900.0)
    tbl["Priority"] = np.full(n_target, 1, dtype=int)

    filternames = [["g_gaia", "bp_gaia", "rp_gaia"]] * n_target
    totalfluxes = np.empty(n_target, dtype=object)

    for i in range(n_target):
        totalfluxes[i] = np.array(
            [
                tbl_tmp["g_flux_njy"][i],
                tbl_tmp["bp_flux_njy"][i],
                tbl_tmp["rp_flux_njy"][i],
            ]
        )

    # FIXME: I think it is worth putting the table file in a non-tmp directory
    with tempfile.NamedTemporaryFile(dir="/tmp", delete=False) as tmpfile:
        outfile = tmpfile.name
    tbl.write(outfile, format="ascii.ecsv", overwrite=True)

    tbl["totalFlux"] = totalfluxes
    tbl["filterNames"] = filternames
    tbl["target_type_id"] = np.full(n_target, 1)  # 1: SCIENCE
    tbl["input_catalog_id"] = np.full(n_target, 2)  # 2: gaia_dr2

    return outfile, tbl


def gen_assignment(args, listname_targets, listname_fluxstds):
    tgt = nf.readScientificFromFile(listname_targets, "sci")
    tgt += nf.readCalibrationFromFile(listname_fluxstds, "cal")
    # tgt += nf.readCalibrationFromFile(listname_sky, "sky")  # need a list of sky positions, which looks very hard.
    cobraCoach, bench = getBench(args)
    telescopes = [nf.Telescope(args.ra, args.dec, args.pa, args.observation_time)]

    # get focal plane positions for all targets and all visits
    tpos = [tele.get_fp_positions(tgt) for tele in telescopes]

    # create the dictionary containing the costs and constraints for all classes
    # of targets
    # For the purpose of this demonstration we assume that all targets are
    # scientific targets with priority 1.
    classdict = {
        "sci_P1": {
            "nonObservationCost": 100,
            "partialObservationCost": 1000,
            "calib": False,
        },
        "cal": {
            "numRequired": args.n_fluxstd,
            "nonObservationCost": 1e5,
            "calib": True,
        },
        "sky": {
            "numRequired": 100,
            "nonObservationCost": 1e6,
            "calib": True,
        },
    }
    tclassdict = {"sci_P1": 1, "sky": 2, "cal": 3}

    t_obs = 900.0

    alreadyObserved = {}
    forbiddenPairs = []
    for i in range(1):
        forbiddenPairs.append([])

    # We penalize targets near the edge of a patrol region slightly to reduce
    # the chance of endpoint collisions with unllocated Cobras
    # (see note below).
    def cobraMoveCost(dist):
        return 0.1 * dist

    gurobiOptions = dict(
        seed=0,
        presolve=1,
        method=4,
        degenmoves=0,
        heuristics=0.8,
        mipfocus=0,
        mipgap=1.0e-04,
    )

    done = False
    while not done:
        # compute observation strategy
        prob = nf.buildProblem(
            bench,
            tgt,
            tpos,
            classdict,
            t_obs,
            None,
            cobraMoveCost=cobraMoveCost,
            collision_distance=2.0,
            elbow_collisions=True,
            gurobi=args.use_gurobi,
            gurobiOptions=gurobiOptions if args.use_gurobi else None,
            alreadyObserved=alreadyObserved,
            forbiddenPairs=forbiddenPairs,
        )

        print("solving the problem")
        prob.solve()

        # extract solution
        res = [{} for _ in range(1)]
        for k1, v1 in prob._vardict.items():
            if k1.startswith("Tv_Cv_"):
                visited = prob.value(v1) > 0
                if visited:
                    _, _, tidx, cidx, ivis = k1.split("_")
                    res[int(ivis)][int(tidx)] = int(cidx)

        # NOTE: the following block would normally be used to "fix" the trajectory
        # collisions detected by the collision simulator.
        # However, this does not work currently, since the current version of
        # cobraCharmer does not actively move unassigned Cobras out of the way of
        # assigned ones, which can result in endpoint collisions which the fiber
        # assigner itself cannot avoid (since it does not know anything about the
        # positioning of unassigned Cobras).
        # So we skip this for now, hoping that it will become possible again with future
        # releases of cobraCharmer.

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
                    dist = np.abs(selectedTargets[i] - bench.cobras.centers[i])

            simulator = CollisionSimulator2(
                bench, cobraCoach, TargetGroup(selectedTargets, ids)
            )
            simulator.run()
            # If you want to see the result of the collision simulator, uncomment the next three lines
            #            from ics.cobraOps import plotUtils
            #            simulator.plotResults(paintFootprints=False)
            #            plotUtils.pauseExecution()
            #
            #            if np.any(simulator.endPointCollisions):
            #                print("ERROR: detected end point collision, which should be impossible")
            #                raise RuntimeError()
            coll_tidx = []
            for tidx, cidx in vis.items():
                if simulator.collisions[cidx]:
                    coll_tidx.append(tidx)
            ncoll += len(coll_tidx)
            for i1 in range(0, len(coll_tidx)):
                found = False
                for i2 in range(i1 + 1, len(coll_tidx)):
                    if np.abs(tp[coll_tidx[i1]] - tp[coll_tidx[i2]]) < 10:
                        forbiddenPairs[ivis].append((coll_tidx[i1], coll_tidx[i2]))
                        found = True
                if not found:  # not a collision between two active Cobras
                    forbiddenPairs[ivis].append((coll_tidx[i1],))

        print("trajectory collisions found:", ncoll)
        done = ncoll == 0

    return res[0], tpos[0], telescopes[0], tgt, tclassdict


def generate_pfs_design(vis, tp, tel, tgt, classdict, tbl_targets, tbl_fluxstds):

    n_fiber = 2394
    fiber_id = np.arange(n_fiber, dtype=int) + 1  # fiberID starts with 0 or 1?

    idx_array = np.arange(n_fiber)

    ra = np.full(n_fiber, np.nan)
    dec = np.full(n_fiber, np.nan)
    pfiNominal = np.full((n_fiber, 2), [np.nan, np.nan])
    catId = np.full(n_fiber, -1, dtype=int)
    objId = np.full(n_fiber, -1, dtype=np.int64)
    targetType = np.full(n_fiber, 4, dtype=int)  # filled as unassigned number

    totalFlux = [np.array([np.nan, np.nan, np.nan])] * n_fiber
    psfFlux = [np.array([np.nan, np.nan, np.nan])] * n_fiber
    filterNames = [["none", "none", "none"]] * n_fiber

    for tidx, cidx in vis.items():

        idx_fiber = fiber_id == (cidx + 1)
        i_fiber = idx_array[idx_fiber][0]

        ra[idx_fiber] = tgt[tidx].ra
        dec[idx_fiber] = tgt[tidx].dec
        # netflow's Target class convert object IDs to string.
        objId[idx_fiber] = np.int64(tgt[tidx].ID)
        pfiNominal[idx_fiber] = [tp[tidx].real, tp[tidx].imag]
        targetType[idx_fiber] = classdict[tgt[tidx].targetclass]

        idx_target = np.logical_and(
            tbl_targets["ID"] == np.int64(tgt[tidx].ID),
            tbl_targets["target_type_id"] == classdict[tgt[tidx].targetclass],
        )
        idx_fluxstd = np.logical_and(
            tbl_fluxstds["ID"] == np.int64(tgt[tidx].ID),
            tbl_fluxstds["target_type_id"] == classdict[tgt[tidx].targetclass],
        )

        if np.any(idx_target):
            catId[i_fiber] = tbl_targets["input_catalog_id"][idx_target][0]
            totalFlux[i_fiber] = tbl_targets["totalFlux"][idx_target][0]
            filterNames[i_fiber] = tbl_targets["filterNames"][idx_target][0].tolist()
        if np.any(idx_fluxstd):
            catId[i_fiber] = tbl_fluxstds["input_catalog_id"][idx_fluxstd][0]
            psfFlux[i_fiber] = tbl_fluxstds["psfFlux"][idx_fluxstd][0]
            filterNames[i_fiber] = tbl_fluxstds["filterNames"][idx_fluxstd][0].tolist()

    design = makePfsDesign(
        pfiNominal,
        ra,
        dec,
        raBoresight=tel._ra,
        decBoresight=tel._dec,
        posAng=tel._posang,
        # arms="br",
        # tract=1,
        # patch="1,1",
        catId=catId,
        objId=objId,
        targetType=targetType,
        # fiberStatus=FiberStatus.GOOD,
        # fiberFlux=np.NaN,
        psfFlux=psfFlux,
        totalFlux=totalFlux,
        # fiberFluxErr=np.NaN,
        # psfFluxErr=np.NaN,
        # totalFluxErr=np.NaN,
        filterNames=filterNames,
        # guideStars=None,
        # designName=None,
    )

    return design


def create_guidestars_from_gaiadb(args):
    # Get ra, dec and position angle from input arguments
    raTel_deg, decTel_deg, pa_deg = args.ra, args.dec, args.pa

    # this should come from the pfsDesign as well, but is not yet in there
    # (DAMD-101)
    obs_time = args.observation_time

    guidestar_mag_max = args.guidestar_mag_max
    guidestar_neighbor_mag_min = args.guidestar_neighbor_mag_min
    guidestar_minsep_deg = args.guidestar_minsep_deg

    # guide star cam geometries
    agcoord = guidecam_geometry()

    # internal, technical parameters
    # set focal plane radius
    fp_rad_deg = 260.0 * 10.2 / 3600
    fp_fudge_factor = 1.2

    # Find guide star candidates
    conn = connect_subaru_gaiadb(args.gaiadb_conf)
    # cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cur = conn.cursor()
    coldict = {
        "id": "source_id",
        "ra": "ra",
        "dec": "dec",
        "parallax": "parallax",
        "pmra": "pmra",
        "pmdec": "pmdec",
        "mag": "phot_g_mean_mag",
        "color": "bp_rp",
    }
    racol, deccol = coldict["ra"], coldict["dec"]
    # req_columns = [
    #     coldict["id"],
    #     racol,
    #     deccol,
    #     coldict["pmra"],
    #     coldict["pmdec"],
    #     "phot_g_mean_mag",
    # ]

    query_string = """SELECT source_id,ra,dec,parallax,pmra,pmdec,phot_g_mean_mag,bp_rp
    FROM gaia
    WHERE q3c_radial_query(ra, dec, {:}, {:}, {:})
    AND {:s} IS NOT NULL AND {:s} IS NOT NULL
    AND {:s} BETWEEN {:} AND {:}
    ;
    """.format(
        raTel_deg,
        decTel_deg,
        fp_rad_deg * fp_fudge_factor,
        coldict["pmra"],
        coldict["pmdec"],
        "phot_g_mean_mag",
        0.0,
        guidestar_neighbor_mag_min,
    )

    cur.execute(query_string)

    df_res = pd.DataFrame(
        cur.fetchall(),
        columns=[
            "source_id",
            "ra",
            "dec",
            "parallax",
            "pmra",
            "pmdec",
            "phot_g_mean_mag",
            "bp_rp",
        ],
    )

    res = {}
    for col in df_res.columns:
        res[col] = df_res[col].to_numpy()

    cur.close()
    conn.close()

    # # FIXME: run similar query, but without the PM requirement, to get a list of
    # # potentially too-bright neighbours

    # adjust for proper motion
    epoch = Time(args.observation_time).jyear
    res[racol], res[deccol] = update_coords_for_proper_motion(
        res[racol],
        res[deccol],
        res[coldict["pmra"]],
        res[coldict["pmdec"]],
        2015.5,  # Gaia DR2 uses 2015.5
        epoch,
    )

    # compute PFI coordinates
    tmp = np.array([res[racol], res[deccol]])
    tmp = ctrans(
        xyin=tmp,
        za=0.0,
        mode="sky_pfi",
        inr=0.0,
        pa=pa_deg,
        cent=np.array([raTel_deg, decTel_deg]),
        time=obs_time,
    )

    res["xypos"] = np.array([tmp[0, :], tmp[1, :]]).T

    # determine the subset of sources falling within the guide cam FOVs
    # For the moment I'm using matplotlib's path functionality for this task
    # Once the "pfi_sky" transformation direction is available in
    # pfs_utils.coordinates, we can do a direct polygon query for every camera,
    # which should be more efficient.
    targets = {}
    tgtcam = []

    for i in range(agcoord.shape[0]):

        p = mppath.Path(agcoord[i])

        # find all targets in the slighty enlarged FOV
        tmp = p.contains_points(res["xypos"], radius=1.0)  # 1mm more
        tdict = {}

        for key, val in res.items():
            tdict[key] = val[tmp]

        # eliminate close neighbors
        flags = flag_close_pairs(tdict[racol], tdict[deccol], guidestar_minsep_deg)

        for key, val in tdict.items():
            tdict[key] = val[np.invert(flags)]

        # eliminate all targets which are not bright enough to be guide stars
        flags = tdict["phot_g_mean_mag"] < guidestar_mag_max

        for key, val in tdict.items():
            tdict[key] = val[flags]

        # eliminate all targets which are not really in the camera's FOV
        flags = p.contains_points(tdict["xypos"])  # 1mm more

        for key, val in tdict.items():
            tdict[key] = val[flags]

        # add AG camera ID
        tdict["agid"] = [i] * len(tdict[coldict["id"]])

        # compute and add pixel coordinates
        tmp = []
        for pos in tdict["xypos"]:
            tmp.append(ag_pfimm_to_pixel(i, pos[0], pos[1]))
        tdict["agpix_x"] = np.array([x[0] for x in tmp])
        tdict["agpix_y"] = np.array([x[1] for x in tmp])

        # append the results for this camera to the full list
        tgtcam.append(tdict)
        for key, val in tdict.items():
            if key not in targets:
                targets[key] = val
            else:
                targets[key] = np.concatenate((targets[key], val))

    # Write the results to a new pfsDesign file. Data fields are according to
    # DAMD-101.
    # required data:
    # ra/dec of guide star candidates: in racol, deccol
    # PM information: in pmra, pmdec
    # parallax: currently N/A
    # flux: currently N/A
    # AgId: trivial to obtain from data structure
    # AgX, AgY (pixel coordinates): only computable with access to the full
    #   AG camera geometry
    # output_design = input_design

    ntgt = len(targets[coldict["id"]])
    guidestars = pfs.datamodel.guideStars.GuideStars(
        targets[coldict["id"]],
        np.full(ntgt, "J{:.1f}".format(epoch)),  # convert float epoch to string
        targets[coldict["ra"]],
        targets[coldict["dec"]],
        targets[coldict["pmra"]],
        targets[coldict["pmdec"]],
        targets[coldict["parallax"]],
        targets[coldict["mag"]],
        np.full(ntgt, "g_gaia"),  # passband
        targets[coldict["color"]],  # color
        targets["agid"],  # AG camera ID
        targets["agpix_x"],  # AG x pixel coordinate
        targets["agpix_y"],  # AG y pixel coordinate
        args.telescope_elevation,  # telescope elevation, don't know how to obtain,
        2,  # numerical ID assigned to the GAIA catalogue
    )

    return guidestars


def main():

    args = get_arguments()

    print(args)
    # exit()

    for d in [args.design_dir, args.cobra_coach_dir]:
        try:
            os.makedirs(d, exist_ok=False)
        except:
            pass

    listname_fluxstds, tbl_fluxstds = gen_target_list_from_targetdb(args)
    listname_targets, tbl_targets = gen_target_list_from_gaiadb(args)
    vis, tp, tel, tgt, classdict = gen_assignment(
        args, listname_targets, listname_fluxstds
    )
    design = generate_pfs_design(
        vis, tp, tel, tgt, classdict, tbl_targets, tbl_fluxstds
    )

    guidestars = create_guidestars_from_gaiadb(args)
    design.guideStars = guidestars

    design.write(dirName=args.design_dir, fileName=design.filename)

    print(
        "pfsDesign file {:s} is created in the {:s} directory.".format(
            design.filename, args.design_dir
        )
    )


if __name__ == "__main__":
    main()

# Example:
# python ./commissioning_demo_2_monodera.py --use_gurobi --design_dir="design2" --cobra_coach_dir="cobracoach" --ra=150 --dec=2 --targetdb_conf ../../../database_configs/targetdb_config_pfsa-db01-gb.ini --gaiadb_conf ../../../database_configs/gaiadb_config_hilo.ini
