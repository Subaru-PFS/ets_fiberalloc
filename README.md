## Overview over the files:

Current development:

- `ets_fiber_assigner/netflow.py`:
  fiber assignment tool based on a network flow algorithm

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


### Partial recalculation

It is possible to do a partial recalculation of an assignment after some visits
have already been observed. This can become necessary for a large variety of
reasons, a few examples being
- a fiber was discovered to be broken during the first visits
- some targets were not observed due to fiber collisions
- new, very high priority targets are added to the target lists
- it is discovered that some targets need more observation time than
  anticipated

In such a situation, the data structure must be set up in the same way as it
would be for the full assignment calculation, but in addition a dictionary is
passed to the assigner, which simply contains target IDs of already observed
targets and the number of visits they have already been observed. This
information is sufficient to compute the optimal assignment strategy for the
remaining visits given the new parameters.

Example: An assignment is computed for a given target list and 5 visits.
After three visits have been carried out, it is noticed that Cobra #0007 does
not work.
To compute the optimal strategy for the remaining two visits:
- make a dictionary of all science targets and the number of times they have
  been observed during the first three visits (excluding the one observed by
  Cobra #0007).
- remove Cobra #0007 from the list of active cobras
- run a fiber assignment task for 2 visits, passing the updated Cobra list and
  the above dictionary.


In case of any questions, please don't hesitate to contact me
(martin@mpa-garching.mpg.de)!

Martin Reinecke
