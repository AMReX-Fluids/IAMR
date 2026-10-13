
Brief Overview of Low Speed Approximations
==========================================

The IAMR code can be used as a basis for more
general low Mach number flow algorithms (e.g., low Mach number combustion,
low Mach number astrophysics, porous media flow)
There are many low speed formulations of the equations of hydrodynamics
in use, each with their own applications. All of these methods share in
common a constraint equation on the velocity field that augments the
equations of motion.

The simplest low Mach number approximation is incompressible
hydrodynamics. This approximation is formally the :math:`M \rightarrow 0`
limit of the Navier-Stokes equations. In incompressible hydrodynamics,
the velocity satisfies a constraint equation:

.. math:: \nabla \cdot {{\bf U}}= 0

which acts to instantaneously equilibrate the flow, thereby filtering
out soundwaves. The constraint equation implies that

.. math:: D\rho/Dt = 0

(through the continuity equation) which says that the density is
constant along particle paths. This means that there are no
compressibility effects modeled in this approximation.

IAMR uses a constraint of the form

.. math:: \nabla \cdot {{\bf U}} = S

which filters sound waves while capturing compressibilty effects due to thermal
diffusion.
Projection methodology is used to enforce the constraint.
To achieve second order accuracy, IAMR includes two projections per timestep.
The first (the ‘MAC’ projection :cite:`bellColellaHowell:1991`)
operates on the half-time, edge-centered advective velocities, making
sure that they satisfy the divergence constraint. These advective
velocities are used to construct the fluxes through the interfaces to
advance the solution to the new time. The second/final projection
operates on the cell-centered velocities at the new time, again
enforcing the divergence constraint. Additional information on the projections
is in AMReX-Hydro's documentation: :ref:`hydro:projections`.

