Introduction
===================

IAMR is a massively parallel, adaptive mesh refinement (AMR) code that 
solves the variable-density incompressible Navier-Stokes equations in 2-D or 3-D with the
option for an embedded boundary (cut cell) representation of complex geometries. 

It is built on top of AMReX, a publicly available software framework designed for building
massively parallel block-structured adaptive mesh refinement (AMR)
applications.

Another AMReX-based code, `incflo <https://amrex-codes.github.io/incflo/>`_, also solves the variable-density incompressible 
Navier-Stokes equations in 2-D or 3-D but does not support subcycling in time.

Key software and algorithmic features of IAMR include:

 + Fluid velocity, density and tracers are defined at cell centroids; pressure is defined at nodes.

 + Possible advection algorithms: a Method-Of-Lines (MOL) approach and a Godunov-method algorithm. Both use an intermediate MAC projection for face-centered advection velocities.

 + Incompressibility of the fluid is imposed through the use of an approximate projection at the end of the time step.

 + Implicit or explicit discretization of viscous terms with variable viscosity.

 + The representation of the complex geometry uses the embedded boundary, or cut-cell, approach.

 + Hybrid parallelization via MPI+X where X = OpenMP for multicore machines, and CUDA/HIP/DCP++ for CPU/GPU systems.

 + Parallel I/O using AMReX native I/O or HDF5.

 + Plotfile format supported by AmrVis, VisIt, ParaView, and yt.
  
 + Optional subcycling in time.

 + Support for particles.

The IAMR source code can be found at https://github.com/AMReX-Codes/IAMR/.
IAMR heavily leverages AMReX (see https://amrex-codes.github.io/) which is supported by
ECP as part of the AMReX Co-Design Center.

Active development in IAMR is ongoing in the development branch. 
Changes are merged into the main branch at the beginning of each month.
