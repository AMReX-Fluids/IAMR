# IAMR
AMReX Support for Low-Mach Number Algorithms

IAMR is a parallel, adaptive mesh refinement (AMR) code that solves the variable-density incompressible Navier-Stokes equations.

The core libraries for managing the subcycling AMR grids and communication are found in the AMReX library (see https://github.com/AMReX-Codes/amrex).

The algorithm is described in the following paper (and references therein):

-A Conservative Adaptive Projection Method for the Variable Density Incompressible Navier-Stokes Equations, A. S. Almgren, J. B. Bell, P. Colella, L. H. Howell, and M. L. Welcome,
J. Comp. Phys., 142, pp. 1-46, 1998.

Key software and algorithmic features of IAMR include:

-Written in cpp and Fortran.

-2- and 3-D support

-Optional subcycling in time

-Support for particles

-Parallelization via flat MPI, hybrid MPI/OpenMP, or MPI/MPI

-Parallel I/O

-Plotfile format supported by VisIt, yt, and AmrVis

-Second-order projection methodology for enforcing the incompressibility constraint

-Several higher-order Godunov integration schemes for advection.

-Implicit viscosity

Refer to the guide in IAMR/UsersGuide for the most current information.
