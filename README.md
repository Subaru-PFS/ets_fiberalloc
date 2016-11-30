## Overview over the files:

src/ets_demo.cc:
Main file containing the experimental fiber assignment code.

src/astronomy.*:
Code for computing the RA/Dec -> alt/AZ conversion.

Directory src/external/
C and C++ sources that were originally developed for the Planck simulation
pipeline and can be re-used for ETS

## Compiling the code:

Simply type "make". This requires a fairly recent version of GNU g++ (tested
with version 5.3 and above, but any 5.x will probably work).

## Running the demo:

Try, for example:
`./ets_demo assigner=naive input=data/ets_test_data.dat n_exposures=10 output=output.txt time=2016-04-03T08:00:00Z`

Supported values for "assigner" are "naive", "draining" and "new".

The algorithms are documented in the source code.

- "n_exposures" is the number of individual observations performed (with
  potential repointing of the telescope in between).

- "time" is an ISO 8601 string containing the date and time of the observation.
  This is needed to calculate the telescope elevation, the exact azimuth and altitude
  of the targets as seen from Subaru, and the distortion of target positions in the
  instrument's focal plane.
  NOTE: Currenty, only strings in the exact format "yyyy-mm-ddThh:mm:ssZ" are
  accepted. If necessary, the parser can be made more flexible.

- "output" is an optional parameter. If present, it is the name of a file into
which detailed fiber assignment data is written, including target IDs and the
fiber IDs to which they are assigned, as well as focal plane coordinates and
RA/DEC of the targets.

- "ra" and "dec" are optional parameters specifying the approximate telescope
pointing (in degrees). If not specified, the code will use the geometrical
center of the input data.

- "posang" is an optional parameter specifying the desired position angle (in
degrees) of the PFS. If unspecified, it is assumed to be 0.

- "dposang" is an optional parameter describing the maximum deviation (in degrees)
from the specified position angle when searching for the optimal telescope
orientation. Default is 4 degrees.

- "nposang" is the number of different position angles tried (in the range
posang+-dposang) when searching for the optimal telescope orientation.
Default is 5.

- "dptg" is an optional parameter describing the maximum deviation (in millimeters
on the PFI plane) from the specified RA/DEC when searching for the optimal
telescope pointing. Default is 4mm.

- "nptg" is the number of different pointings tried (in both directions)
when searching for the optimal telescope pointing. Default is 5.

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

The code first loads the required data set, and then runs the fiber assignment
algorithm assuming that the telescope points at the specified position
position.

For each exposure, the code
- determines which targets are observed with which fibers
  (this depends on the concrete assigner selected by the user)
  This is repeated nptg x nptg x nposang times with the telescope shifted by small amounts in
  x and y directions, and rotating the PFS by small amounts, and the attempt
  with the most assigned fibers is chosen.
  The exposure time is the minimum remaining observation time of all assigned
  targets.
  These lists are written to the output file, if specified.
- subtracts the exposure time from the planned exposure times of the observed
  targets. If their time becomes <=0, they are removed from the target list.
- repeats the above steps until the requested number of exposures is reached.

At each step, the fraction of allocated fibers and the accumulated fraction of
total observation time is printed.

Target priorities are taken into account by the assignment algorithms.

### Short description of the assigner algorithms:

1. Naive:
  - For each fiber, assign the target with the highest priority visible with this
    fiber.

2. Draining:
  - While targets are still visible with any fiber:
    - find the fiber with the lowest number of visible targets, and assign the
      visible target with the highest prioity to it

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
