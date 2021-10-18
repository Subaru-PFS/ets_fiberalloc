# Introduction
For the 2021 commissioning runs MPA/MPE have created a commissioning script
which will for a give RA/Dec/PA query for stars from the Gaia public archive
and then run netflow and shuffle. The result is a pfs_design file.

# Caveats
A number of workarounds were necessary to deal with current shortcomings:
* By default, the code reads the cobra information from the file
ALL_final_20210920_mm.xml stored in pfs_instdata.
This can be overridden by setting the command line parameter `--cobra_coach_module_version` to something different from its default, which is `final_20210920_mm`.
* The data stored in ALL_final_20210920_mm.xml is not correct for some of
the cobras and cobraops does adapts it in order to be able to run the collision simulations:
   * The cobra centers for bad cobras is set to zero --> Cobraops
sets the positions to be outside the PFI field of view, so colliding with those cobras is
impossible. That's not true and should be fixed in the XML.
   * Link lengths are zero for some cobras --> Cobraops sets those
lengths to the median value measured for the rest of the cobras.
   * Link lengths are too large for some cobras, probably associated
to bad measurements --> Cobraops sets those lengths to the median value measured for the rest of the cobras.
* Cobraops sets the position of those cobras not assigned to targets to
their home position. Maybe a different strategy should be use to avoid fiber-to-elbow collisions.
Thus, unless also a sufficiently large number of sky positions are give (which this
example will not take care of) there will be unassigned fibers in the
pfs_design file. Endpoint and trajectory collisions are only handled (and
avoided) for the allocated fibers.  Appropriate measures in the down stream
processing and cobra commanding must be takes to avoid collisions.
* Cobraops uses cobracharmer to calculate the cobras trajectories using
the following method:
https://github.com/Subaru-PFS/ics_cobraOps/blob/master/python/ics/cobraOps/CollisionSimulator2.py#L107
We should make sure that this is the same method/strategy used to
move the cobras.
* The position of the black spots is hard-coded  in cobraops and is
probably not correct. We should find a method to read those positions from the database or pfs_instdata.
* It is possible that a cobra can reach below the black dot of a neighboring
cobra. This is currently not caught and no light will be observed from the
associated object.


# Installation
## Choice of optimiser
You will have to install either the Gurobi Optimiser OR the COIN optimiser.
COIN is open source and available without a license file but by far not as powerful as gurobi. 
For gurobi you will have to create an academic account (one-step email verification) and
create a license file. We found this to be rather straightforward. Be aware though that the
license is only good for one year. You can renew it thereafter.

For COIN:
	https://projects.coin-or.org/Cbc
	brew tap coin-or-tools/coinor
	brew install cbc
Do also install the python bindings:
	pip install pulp==1.6.0

For Gurobi:
Navigate to 
https://www.gurobi.com/academia/academic-program-and-licenses/
get an academic account, install gurobi and get an academic license.
Do not forget to the grbgetkey command as indicated when you request the license and
then install the python bindings.
pip install gurobipy

## Prerequisites
(Optional) You may want to consider to work inside of a virtual environment:

	python3 -m venv netflow
	source netflow/bin/activate
	pip install --upgrade pip

Create a work directory:

	mkdir fiber_allocation
	cd fiber_allocation
	
Download the fiber allocation code:
	git clone -b commissioning_demo git@github.com:Subaru-PFS/ets_fiberalloc.git

Create a subdirectory to host the actual pfs_designs:

	mkdir commissioning
	cd commissioning
	cp ../ets_fiberalloc/commissioning.py .
	cd ..
  
Install pip-installable prerequisites:
Do

	pip install pulp==1.6.0
	pip install opencv-python
	pip install psycopg2
	pip install astroquery
	pip install scipy 
	pip install pybind11
	pip install pyyaml
	pip install sep
	pip install pandas
	pip install sqlalchemy
	
OR just

	python -m pip install -r ets_fiberalloc/requirements.txt


Install PFS specific prerequisites:

	git clone https://github.com/Subaru-PFS/ics_cobraOps.git
	git clone https://github.com/Subaru-PFS/ics_fpsActor.git
	git clone https://github.com/Subaru-PFS/spt_operational_database.git
	git clone https://github.com/Subaru-PFS/datamodel.git
	git clone https://github.com/Subaru-PFS/ics_cobraCharmer
	git clone https://github.com/Subaru-PFS/pfs_instdata
	git clone git@github.com:Subaru-PFS/ets_shuffle.git
	git clone git@github.com:Subaru-PFS/pfs_utils.git
  
	cd datamodel/
	pip install -e .
	cd ../ics_cobraCharmer
	pip install -e .
	cd ../ics_cobraOps
	pip install -e .
	cd ..

Set the enviroment:

	export PYTHONPATH=$PYTHONPATH:../spt_operational_database/python
	export PYTHONPATH=$PYTHONPATH:../ics_fpsActor/python
	export PYTHONPATH=$PYTHONPATH:../ets_shuffle/
	export PYTHONPATH=$PYTHONPATH:../ets_fiberalloc
	export PFS_INSTDATA_DIR=../pfs_instdata

# Execution
Minimal at RA = 0 Deg and Dec = 0 Deg.:

	cd commissioning
	python commissioning.py

More realistic example with gurobi:

	python commissioning.py --use_gurobi True --ra 56.745 --dec 24.113 --pa 0

(CAUTION: due to some weird logic in the "argparse" package, you have to
completely omit the "--use_gurobi" flag if you don't want to use Gurobi.
Specifying "--use_gurobi False" will NOT do the trick ...)

Help arguments:

	python commissioning.py --help
