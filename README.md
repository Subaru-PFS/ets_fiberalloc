## Overview over the files:

src/pyETS.cc:
wrapper file describing the functionality exported to Python

src/ets.*, src/ets_helpers.h:
C++ implementation of the assigners and associated functionality

Directory src/external/
C and C++ sources that were originally developed for the Planck simulation
pipeline and can be re-used for ETS

## Compiling/Installing:

Simply do a "python setup.py install" or similar. "pip install ." should also
work.

## Demo code:

The file demo.py contains a very brief overview over the typical workflow, i.e.
reading in target data from a file, converting target positions to x/y
coordinates on the focal plane, determining which targets can be seen by which
cobras, and performing a single assignment.

The input file containing the targets is ASCII and contains one target per line.
The columns within a line contain the following information:

 1. ID
 2. R.A. [deg.]
 3. Dec. [deg.]
 4. Exposure Time [sec.]
 5. Priority [1(highest) - 15(lowest)]
 6. Magnitude [ABmag]. (optional)
 7. Redshift (optional)
 8. Object Type (optional)

### Short description of the assigner algorithms:

1. Naive:
  - For each fiber, assign the target with the highest priority visible with this
    fiber. If there is more than one valid target, choose one randomly.

2. Draining:
  - While targets are still visible with any fiber:
    - find the fiber with the lowest number of visible targets, and assign the
      visible target with the highest prioity to it. If there is more than one
      valid target, choose one randomly.

2. Draining_closest:
  - While targets are still visible with any fiber:
    - find the fiber with the lowest number of visible targets, and assign the
      visible target with the highest prioity to it. If there is more than one
      valid target, choose the one closest to the Cobra center.

3. New:
  - While targets are still visible with any fiber:
    - compute an importance function for every target, which depends on number and
      proximity of nearby targets and on the target's remaining observation time.
    - from the list of visible targets with the highest priority, assign the target
      with the highest importance to a fiber. If multiple fibers can see the target,
      choose the fiber with the lowest number of potential targets.
    - the importance function is designed to measure the 'crowdedness' at the
      location of a given source. For target i, it is given as
      I_i = \sum_j t_i t_j K(r_ij)
      where t_i and t_j are the remaining observation times for targets i and j, and
      r_ij is the distance between the targets. The kernel function K should be 1 for
      r_ij=0 and drop off to zero for radii around the fiber patrol radius.
      This function ensures that each fiber is assigned to more 'crowded' areas of the
      target field first, with the goal of homogenizing the distribution of targets.

### Comments on the choice of algorithm:

It is a justified question why ETS does not simply make use of an assigning
algorithm based on the maximum-flow graph approach (Goldberg 1997, J. Algorithms
22, 1), similar to SDSS (Blanton et al. 2003, AJ 125, 2276). This algorithm is
proven to be optimal and runs in polynomial time.

Unfortunately the characteristics of both the instruments and the surveys differ
very much between SDSS and PFS. In SDSS
- all fibers could be placed almost anywhere on the whole focal plane, as long
  as they were not too close to other fibers
- the target distribution contained only very few target groups that could lead
  to fiber collisions.

For PFS
- every fiber is only movable within a circular patrol area much smaller than
  the focal plane
- the test target catalogs are (for a large part) so spatially dense that almost
  all targets are in a group with many potentially colliding neighbors.

The reduced patrol area of individual fibers implies that the optimum
observation strategy depends much more on telescope pointing and orientation
than it does for SDSS. However, selecting optimal pointing and orientation is
not covered by the network flow algorithm and is a NP-hard problem.

SDSS had to adopt a special (probably not optimal) approach for potentially
colliding targets; for its target distribution this was not a big problem. In
the case of the planned PFS surveys, the fact that practically all targets are
potentially colliding defeats the purpose of the network-flow approach.

It appears that other groups developing assignment algorithms for instruments
similar to PFS have come to the same conclusion (see, e.g., Morales et al. 2012,
MNRAS 419, 1187).


In case of any questions, please don't hesitate to contact me
(martin@mpa-garching.mpg.de)!

Martin Reinecke
