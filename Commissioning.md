# Introduction
For the 2021 commissioning runs MPA/MPE have created a commissioning script
which will for a give RA/Dec/PA query for stars from the Gaia public archive
and then run netflow and shuffle. The result is a pfs_design file.

# Caveats
A number of workarounds were necessary to deal with current shortcomings:
* Some link lengths have unexpected values, we currently set them to a median value
* Broken cobras have no link lengths assigned with them, we remove them by placing them outside of the 
focal plane
* The black dot positions are currently not taken from an actual measurement but a nominal position is used.
* It is possible that a cobra can reach below the black dot of a neighboring
cobra. This is currently not caught and no light will be observed from the
associated object.
* Unless also a sufficiently large number of sky positions are give (which this
example will not take care of) there will be unassigned fibers in the
pfs_design file.  Endpoint and trajectory collisions are only handled (and
avoided) for the allocated fibers.  Appropriate measures in the down stream
processing and cobra commanding must be takes to avoid collisions.

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

Create a subdirectory to host the actual pfs_designs:

	mkdir commissioning
	cd commissioning
	cp WHEREVER_YOU_DOWNLOADED_IT/commissioning.py .
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

	cp WHEREVER_YOU_DOWNLOADED_IT/requirements.txt .
	python -m pip install -r requirements.txt


Install PFS specific prerequisites:

	git clone https://github.com/Subaru-PFS/ics_cobraOps.git
	git clone https://github.com/Subaru-PFS/ics_fpsActor.git
	git clone https://github.com/Subaru-PFS/spt_operational_database.git
	git clone https://github.com/Subaru-PFS/datamodel.git
	git clone https://github.com/Subaru-PFS/ics_cobraCharmer
	git clone https://github.com/Subaru-PFS/pfs_instdata
	git clone -b tickets/FIBERALLOC-28 git@github.com:Subaru-PFS/ets_shuffle.git
	git clone -b commissioning_demo git@github.com:Subaru-PFS/ets_fiberalloc.git
	git clone -b tickets/INSTRM-1037 git@github.com:Subaru-PFS/pfs_utils.git
  
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

Help arguments:

	python commissioning.py --help
