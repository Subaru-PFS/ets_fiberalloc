## Overview over the files:

src/ets_demo.cc
Main file containing the experimental fiber assignment code.

Directory src/external/
C and C++ sources that were originally developed for the Planck simulation
pipeline and can be re-used for ETS

## Compiling the code:

Simply type "make". This requires a fairly recent version of GNU g++ (tested
with version 5.3 and above, but any 5.x will probably work).

## Running the demo:

Try, for example:
`./ets_demo assigner=naive input=<path/to>/Halo.dat fract=0.95 output=output.txt`

Supported values for "assigner" are "naive", "draining" and "new".

The algorithms are documented in the source code.

- "fract" is the fraction of coverage that must be reached before the algorithm
stops. It is computed like this:

  total_time := sum over all sources to be observed times their planned
                observation time

  obs_time   := sum over all sources observed so far times their observation time
                so far

  fract := obs_time/total_time

- "output" is an optional parameter. If present, it is the name of a file into
which detailed fiber assignment data is written.

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
telescope pointing. Default is 4 degrees.

- "nptg" is the number of different pointings tried (in both directions)
when searching for the optimal telescope pointing. Default is 5.

The code first loads the required data set, and then runs the fiber assignment
algorithm assuming that the telescope points at the specified position
position.

Then all sources with a distance >190mm from the center are discarded. This is
done for testing purposes only - to make sure all targets can be observed
easily.

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
- repeats the above steps until a fraction of "fract" of the necessary total
  observation time has been reached.

At each step, the fraction of allocated fibers and the fraction of obervations
done so far and the necessary observations is printed.

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

In case of any questions, please don't hesitate to contact me
(martin@mpa-garching.mpg.de)!

Martin Reinecke