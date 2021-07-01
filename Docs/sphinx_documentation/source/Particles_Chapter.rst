.. role:: cpp(code)
   :language: c++

Particles
=========

IAMR also has the ability to include data-parallel particle simulations. 
Our particles can interact with data defined on a (possibly adaptive)
block-structured hierarchy of meshes. Example applications include
Particle-in-Cell (PIC) simulations, Lagrangian tracers, or particles that exert
drag forces onto a fluid, such as in multi-phase flow calculations.

Within IAMR, we provide an example of passively advected tracer particles in
``IAMR/Exec/run_2d_particles``.

Here we provide a brief introduction to using particles in IAMR. For more detailed information
on particles, see AMReX's documentation: :ref:`amrex:Chap:Particles`.


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   Particles
