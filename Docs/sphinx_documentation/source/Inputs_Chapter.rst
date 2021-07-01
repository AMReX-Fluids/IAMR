.. _Chap:Inputs:

Run-time Inputs
===============

The IAMR executable uses an inputs file at runtime to set and alter the
behavior of the algorithm and initial conditions.
The inputs file, typically named inputsXXX is used to
set AMReX parameters and the control flow in
the IAMR code. Each parameter here has a namespace (like amr.\ *optionname* or ns.\ *optionname*). Parameters
set here can be read using the AMReX ParmParse class
infrastructure.

The inputs file is specified on the commandline. Options can also be specified directly on the
command line in addition to being used in an inputs file.


.. toctree::
   :maxdepth: 1

   InputsProblemDefinition
   InputsTimeStepping 
   InputsInitialization 
   InputsLoadBalancing 
   InputsMultigrid
   InputsPlotFiles 
   InputsCheckpoint
   InputsVerbosity
