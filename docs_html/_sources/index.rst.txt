.. IAMR documentation master file, created by
   sphinx-quickstart on Sun Mar  7 15:21:57 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to IAMR's documentation!
================================

IAMR is a parallel, adaptive mesh refinement (AMR) code that uses
subcycling in time and solves the variable-density incompressible Navier-Stokes equations 
in the presence of complex geometries.
It is built on top of AMReX, a publicly available software framework designed for 
building massively parallel block-structured adaptive mesh refinement (AMR)
applications.

The IAMR source code is available at 
https://amrex-codes.github.io/IAMR/

For an AMReX-based incompressible flow code without subcycling in time,
see incflo (https://amrex-codes.github.io/incflo/)

.. note::
   **This documentation is a work in progress. Please check back later for more content.**

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Introduction_Chapter
   Getting_Started
   Units
   Software_Chapter
   Visualization_Chapter
   Inputs_Chapter
   ManagingGridHierarchy_Chapter
   AMR
   Parallel
   Fluids_Chapter
   EB
   Particles_Chapter
   Debugging
   
.. toctree::
   :caption: References

   references

The copyright notice and license agreement is included in the
IAMR home directory in OpenSource.txt.
