Introduction
===================

IAMR is a parallel, adaptive mesh refinement (AMR) code that 
solves the variable-density incompressible Navier-Stokes equations.

Key software and algorithmic features of IAMR include:

- 2- and 3-D support
  
- Optional subcycling in time

- The representation of the complex geometry uses the embedded boundary, or cut-cell, approach

- Second-order projection methodology for enforcing the incompressibility constraint

- Higher-order Godunov integration schemes for advection using an intermediate MAC projection for face-centered advection velocities
  
- Implicit viscosity

- Support for particles
  
- Hybrid parallelization via MPI+X where X = OpenMP for multicore machines, and CUDA/HIP/DCP++ for CPU/GPU systems
  
- Parallel I/O using AMReX native I/O or HDF5.

- Plotfile format supported by AmrVis, VisIt, ParaView, and yt.
  

The IAMR source code can be found at https://github.com/AMReX-Codes/IAMR/.
The core libraries for managing the subcycling AMR
grids and communication are found in the AMReX library
(see https://github.com/AMReX-Codes/amrex/).
