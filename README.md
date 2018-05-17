## Overview over the files:

Current development:

- `ets_fiber_assigner/netflow.py`:
  fiber assignment tool based on a network flow algorithm


Mostly of historical interest:

- `src/pyETS.cc`:
  wrapper file describing the functionality exported to Python

- `src/ets.*, src/ets_helpers.h`:
  C++ implementation of the assigners and associated functionality

- Directory `src/external/`:
  C and C++ sources that were originally developed for the Planck simulation
  pipeline and can be re-used for ETS

## Compiling/Installing:

### Prerequisites

As far as possible, the package will install its dependencies automatically.
However, it also depends on the "cobraOps" Python package, which currently
has to be installed manually; see https://github.com/Subaru-PFS/ics_cobraOps/
for details.

The package allows to choose between the PULP package and the commercial
(but free for academic use) Gurobi package for solving the network flow
problem. One of those two needs to be installed and the appropriate flag needs
to be set when calling the network solving routine `observeWithNetflow()`.

### Package installation

Simply do a `python setup.py install` or similar. `pip install .` should also
work.


## Demo code:

The file demo_netflow.py contains a very brief overview over the typical
workflow, i.e. reading in target data from a file, converting target positions
to x/y coordinates on the focal plane, determining which targets can be seen by
which cobras, defining a cost function and planning several assignments

The input files containing the targets are ASCII and contain one target per
line. The columns within a line contain the following information:

 1. ID
 2. R.A. [deg.]
 3. Dec. [deg.]
 4. Exposure Time [sec.]
 5. Priority [1(highest) - 15(lowest)]
 6. Magnitude [ABmag]. (optional)
 7. Redshift (optional)
 8. Object Type (optional)

An example file can be found in the `data` subdirectory.

In case of any questions, please don't hesitate to contact me
(martin@mpa-garching.mpg.de)!

Martin Reinecke
