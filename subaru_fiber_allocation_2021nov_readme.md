## Introduction

The fiber allocation script `subaru_fiber_allocation_2021nov.py` is based on the script made by the MPA/MPE team. You should follow their [instruction](https://github.com/Subaru-PFS/ets_fiberalloc/blob/commissioning_demo_2/Commissioning.md) in the [`commissioning_demo_2` branch in the `ets_fiberalloc` repository](https://github.com/Subaru-PFS/ets_fiberalloc/blob/commissioning_demo_2/) about the basic idea.

Just to use it, you can look at the [Running](#running-the-script-to-make-a-pfsdesign-file) section.

## Updates from the MPA/MPE version

- Targets are selected from gaiaDB at Hilo.  This will be updated to use targetDB in the future.
- Flux standard stars are selected from targetDB at Hilo.
- Guide stars are selected from gaiaDB at Hilo.
- Fiber IDs are updated from science fiber ID to real fiber ID including empty and engineering ones.
- Add proper target types and input catalogs (need to be discussed).
- You can pass `--observation_time now` to make a `pfsDesign` file at the time of running the script.
- Use `pfs.utils.pfsDesignUtils.makePfsDesign()` to generate a `pfsDesign` instance.

## Installation

### Create an environment

Move to a shared work directory and create a work directory.
```
cd /work/monodera/
mkdir fiber_allocation_commissioning_2021nov
cd fiber_allocation_commissioning_2021nov
```

Create a virtual environment and activate it (not sure if it conflicts with the above setup).

```
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
```

### Install prerequisites.

```
pip install numpy
pip install scipy
pip install matplotlib
pip install astropy
pip install opencv-python
pip install psycopg2-binary
pip install astroquery
pip install pybind11
pip install pyyaml
pip install sep
pip install pandas
pip install sqlalchemy
pip install astroplan
pip install bokeh
```

Further prerequisites for `targetDB` (thye are actually optional).

```
pip install sqlalchemy-utils
pip install tabulate
pip install logzero
```

You also need a Gurobi package and setup required environment variables.


### Install PFS packages

```
cd ../
mkdir Subaru-PSF

git clone https://github.com/Subaru-PFS/ics_cobraOps.git
git clone https://github.com/Subaru-PFS/ics_fpsActor.git
git clone https://github.com/Subaru-PFS/spt_operational_database.git
git clone https://github.com/Subaru-PFS/datamodel.git
git clone https://github.com/Subaru-PFS/ics_cobraCharmer
git clone https://github.com/Subaru-PFS/pfs_instdata
git clone https://github.com/Subaru-PFS/ets_shuffle.git
git clone https://github.com/Subaru-PFS/pfs_utils.git
git clone https://github.com/Subaru-PFS/ets_fiberalloc.git
```

```
cd ets_fiberalloc
git checkout commissioning_demo_2
```

Link packages.

```
cd datamodel/
pip install -e .
cd ../ics_cobraCharmer
pip install -e .
cd ../ics_cobraOps
pip install -e .
cd ../pfs_utils
pip install -e .
```

I created a `subaru_pfs.pth` containing the following lines and put it in `/work/monodera/fiber_allocation_commissioning_2021nov/venv/lib/python3.6/site-packages`.

```
/work/monodera/Subaru-PFS/spt_operational_database/python
/work/monodera/Subaru-PFS/ics_fpsActor/python/
/work/monodera/Subaru-PFS/ets_shuffle/
/work/monodera/Subaru-PFS/ets_fiberalloc/
/work/monodera/gurobi/gurobi912/linux64/lib/python3.6_utf32/
```

### Install targetDB package

```
cd /work/monodera/Subaru-PFS/
git clone https://github.com/Subaru-PFS/ets_target_database.git
cd ets_target_database
git checkout commissioning_2021nov
pip install -e .
```


## Running the script to make a pfsDesign file

Before you run the script you need a couple of configuration files for gaiaDB and targetDB at Subaru.  Please contact M. Onodera for these files.  A template file is included in `examples/commissioning_2021nov`.


### Some additional setup

You need to setup some additional environment variables.

```
EUPS_DIR="/work/stack/eups/current/"
source "${EUPS_DIR}/bin/setups.zsh"
export EUPS_PKGROOT='https://eups.lsst.codes/stack/redhat/el7/devtoolset-6/miniconda3-4.5.12-1172c30|https://eups.lsst.codes/stack/src'
export PFS_UTILS_DIR=/work/monodera/Subaru-PFS/pfs_utils/

GUROBI_DIR=/work/monodera/gurobi

export GUROBI_HOME="${GUROBI_DIR}/gurobi912/linux64"
export PATH="${GUROBI_HOME}/bin":$PATH
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib
export GRB_LICENSE_FILE="${GUROBI_DIR}/gurobi_`hostname -s`.lic"
```

Move to a working directory and copy the script.
```shell
cd /work/monodera/fiber_allocation_commissioning_2021nov
. ./venv/bin/activate
mkdir work
cd work
# cp ../../Subaru-PFS/ets_target_database/examples/commissioning_2021nov/subaru_fiber_allocation_2021nov.py .
```


### Run the script

```shell
python /work/monodera/Subaru-PFS/ets_target_database/examples/commissioning_2021nov/subaru_fiber_allocation_2021nov.py \
  --use_gurobi \
  --design_dir "design" \
  --cobra_coach_dir "cobracoach" \
  --ra 150 \
  --dec 2 \
  --targetdb_conf ../../Subaru-PFS/database_configs/targetdb_config_pfsa-db01-gb.ini \
  --gaiadb_conf ../../Subaru-PFS/database_configs/gaiadb_config_hilo.ini \
  --pfs_instdata_dir /work/monodera/Subaru-PFS/pfs_instdata

```

The output is saved as `design/pfsDesign-0x6a119411f0b96b28.fits`.

The full options of the script can be found as follows.

```shell
$ python ./subaru_fiber_allocation_2021nov.py -h
usage: subaru_fiber_allocation_2021nov.py [-h] [--ra RA] [--dec DEC] [--pa PA]
                                          [--observation_time OBSERVATION_TIME]
                                          [--lim_target_mag LIM_TARGET_MAG]
                                          [--design_dir DESIGN_DIR]
                                          [--guidestar_mag_max GUIDESTAR_MAG_MAX]
                                          [--guidestar_neighbor_mag_min GUIDESTAR_NEIGHBOR_MAG_MIN]
                                          [--guidestar_minsep_deg GUIDESTAR_MINSEP_DEG]
                                          [--use_gurobi]
                                          [--cobra_coach_dir COBRA_COACH_DIR]
                                          [--cobra_coach_module_version COBRA_COACH_MODULE_VERSION]
                                          [--targetdb_conf TARGETDB_CONF]
                                          [--gaiadb_conf GAIADB_CONF]
                                          [--target_mag_max TARGET_MAG_MAX]
                                          [--target_mag_min TARGET_MAG_MIN]
                                          [--target_mag_filter TARGET_MAG_FILTER]
                                          [--fluxstd_min_prob_f_star FLUXSTD_MIN_PROB_F_STAR]
                                          [--telescope_elevation TELESCOPE_ELEVATION]
                                          [--n_fluxstd N_FLUXSTD]
                                          [--pfs_instdata_dir PFS_INSTDATA_DIR]

optional arguments:
  -h, --help            show this help message and exit
  --ra RA               Telescope center RA [degrees] (default: 0.0)
  --dec DEC             Telescope center Dec [degrees] (default: 0.0)
  --pa PA               Telescope position angle [degrees] (default: 0.0)
  --observation_time OBSERVATION_TIME
                        planned time of observation in UTC (default:
                        2021-11-20 15:00:00)
  --lim_target_mag LIM_TARGET_MAG
                        magnitude of the faintest targets (obsolete)
                        (default:19)
  --design_dir DESIGN_DIR
                        directory for storing PFS designs (default: .)
  --guidestar_mag_max GUIDESTAR_MAG_MAX
                        maximum magnitude for guide star candidates (default:
                        19.)
  --guidestar_neighbor_mag_min GUIDESTAR_NEIGHBOR_MAG_MIN
                        minimum magnitude for objects in the vicinity of guide
                        star candidates (default: 21.)
  --guidestar_minsep_deg GUIDESTAR_MINSEP_DEG
                        radius of guide star candidate vicinity (default:
                        1/3600)
  --use_gurobi          use Gurobi
  --cobra_coach_dir COBRA_COACH_DIR
                        path for temporary cobraCoach files (default: .)
  --cobra_coach_module_version COBRA_COACH_MODULE_VERSION
                        version of the bench decription file (default:
                        final_20210920_mm)
  --targetdb_conf TARGETDB_CONF
                        Config file for targetDB (default:
                        targetdb_config.ini)
  --gaiadb_conf GAIADB_CONF
                        Config file for Subaru's Gaia DB (default:
                        gaiadb_config_hilo.ini
  --target_mag_max TARGET_MAG_MAX
                        Maximum (faintest) magnitude for stars in fibers
                        (default: 19.)
  --target_mag_min TARGET_MAG_MIN
                        Minimum (brightest) magnitude for stars in fibers
                        (default: 0)
  --target_mag_filter TARGET_MAG_FILTER
                        Photometric band (grizyj of PS1) to apply magnitude
                        cuts (default: g)
  --fluxstd_min_prob_f_star FLUXSTD_MIN_PROB_F_STAR
                        Minimum acceptable prob_f_star (default: 0)
  --telescope_elevation TELESCOPE_ELEVATION
                        Telescope elevation in degree (default: 60)
  --n_fluxstd N_FLUXSTD
                        Number of FLUXSTD stars to be allocated. (default: 50)
  --pfs_instdata_dir PFS_INSTDATA_DIR
                        Location of pfs_instdata (default:
                        /Users/monodera/Dropbox/NAOJ/PFS/Subaru-
                        PFS/pfs_instdata/)
```