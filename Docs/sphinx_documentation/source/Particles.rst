Initializing the Particles
==========================

Particles are initialized from an ASCII file, identified in the IAMR inputs file:

**particles.particle\_init\_file =**\ *particle\_file*

| Here *particle\_file* is the user-specified name of the file. The first line in this file is
  assumed to contain the number of particles. Each line after that contains the position of the particle as
| x y z

Output Format
=============

Checkpoint Files
----------------

| The particle positions and velocities are stored in a binary file in each checkpoint directory.
  This format is designed for being read by the code at restart rather than for diagnostics.

Plot Files
----------

If **particles.write\_in\_plotfile =** 1 in the inputs file
then the particle positions and velocities will be written in a binary file in each plotfile directory.

| In addition, we can also
  visualize the particle locations as represented on the grid. The “derived quantity”
  **particle\_count** represents the number of particles in a grid cell.
  To add it to plotfiles, set
| **amr.derive\_plot\_vars = particle\_count**
| in the inputs file

ASCII Particle Files
--------------------

| To generate an ASCII file containing the particle positions and velocities,
  one needs to restart from a checkpoint file from a run with particles, but one doesn’t need to run any steps.
  For example, if chk00350 exists, then one can set:
| **amr.restart = chk00350**
| **max\_step = 350**
| **particles.particle\_output\_file =** *particle\_output*
| which would tell the code to restart from chk00350, not to take any further time steps, and to write an ASCII-format
  file called *particle\_output*.
| This file has the same format as the ASCII input file:
| number of particles
| x y z

Run-time Screen Output
----------------------

The verbosity of the particle-related sections of code can be controled with

| **particles.pverbose** 
