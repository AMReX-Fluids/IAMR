
The Variable-Density Incompressible Navier-Stokes Equations
===========================================================

IAMR solves the variable-density incompressible Navier-Stokes equations:

.. math::
  
   \begin{aligned}
   {{\bf U}}_t + ({{\bf U}}\cdot\nabla){{\bf U}}&=& \frac{1}{\rho}(-\nabla p + \mu\nabla^2 {{\bf U}}+ {\bf H}_U), \\
   \rho_t + \nabla\cdot(\rho {{\bf U}}) &=& 0, \\
   c_t + ({{\bf U}}\cdot\nabla)c &=& k\nabla^2 c + H_c, \\
   \nabla\cdot {{\bf U}}&=& 0,\end{aligned}

where :math:`{{\bf U}}= (u, v, w), \rho, c`, and :math:`p` represent the velocity, density, concentration of an
advected scalar, and pressure, respectively, and :math:`{\bf H}_U = (H_x , H_y , H_z )` represents any external
forces. Here :math:`\mu` is the dynamic viscosity coefficient, :math:`k` is the diffusive coefficient for :math:`c`, and
:math:`H_c` is the source term for :math:`c`. In general one could advect an arbitrary number of scalars,
either passively or conservatively.

Brief Overview of Low Speed Approximations
==========================================

The IAMR code can also used as a basis for more
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
  
Projection Methods 101
======================

We include a brief discussion of projection methodology for
incompressible and low Mach number flow.
The compressible Navier-Stokes equations can be written in the form:

.. math:: {{\bf U}}_t + \nabla \cdot F({{\bf U}}) = S

where :math:`{{\bf U}}` is a vector of conserved quantities, :math:`{{\bf U}}= (\rho, \rho u,
\rho E)`, with :math:`\rho` the density, :math:`u` the velocity, :math:`E` the total
energy per unit mass, and :math:`S` are source terms. This system
can be expressed as a coupled set of advection/diffusion equations:

.. math:: {\bf q}_t + A({\bf q}) \nabla {\bf q} + D = {\cal S}

where :math:`{\bf q}` are called the primitive variables, :math:`A` is the advective
flux Jacobian, :math:`A \equiv \partial F / \partial U`, :math:`D` are diffusion terms,
and :math:`{\cal S}` are the transformed sources. The eigenvalues of the
matrix :math:`A` are the characteristic speeds—the real-valued speeds at which
information propagates in the system, :math:`u` and :math:`u
\pm c`, where :math:`c` is the sound speed. Solution methods for the
compressible equations that are strictly conservative make use of this wave-nature to compute advective fluxes
at the interfaces of grid cells. Diffusive fluxes can be computed
either implicit or explicit in time, and are added to the advective fluxes,
and used, along with the source terms to update the state in time. An
excellent introduction to these methods is provided by LeVeque’s book
:cite:`leveque`. The timestep for these methods is limited by all three processes
and their numerical implementation. Typically, advection terms are treated
time-explicitly, and the time step will be constrained by the time
it takes for the maximum characteristic speed to traverse one grid cell.
However, in low speed flow applications, it can be shown the acoustics
transport very little energy in the system. As a result, the time-step
restrictions arising from numerical treatement of the advection terms
can be unnecessarily limited, even if A-stable methods are used to incorporate
the diffusion and source terms.

In contrast, solving incompressible or low Mach number systems
typically involves a stage where one or more
advection-like equations are solved (representing, e.g. conservation of mass and
momentum), and coupling that advance with a divergence constraint on the velocity field.
For example, the equations of invicid constant-density incompressible flow
are:

.. math::

   \begin{aligned}
   {{\bf U}}_t &=& -{{\bf U}}\cdot \nabla {{\bf U}}- \frac{1}{\rho}\nabla p \label{eq:incompressible_u} \\
   \nabla \cdot {{\bf U}}&=& 0\end{aligned}

Here, :math:`{{\bf U}}` represents the velocity vector
and :math:`p` is the dynamical pressure. The time-evolution equation for
the velocity (Eq. [eq:incompressible\_u]) can be solved using
techniques similar to those developed for compressible hydrodynamics,
updating the old velocity, :math:`{{\bf U}}^n`, to the new time-level, :math:`{{\bf U}}^\star`.
Here the “:math:`^\star`” indicates that the updated velocity does not, in
general, satisfy the divergence constraint. A projection method will
take this updated velocity and force it to obey the constraint
equation. The basic idea follows from the fact that any vector
field can be expressed as the sum of a divergence-free quantity and
the gradient of a scalar. For the velocity, we can write:

.. math:: {{\bf U}}^\star = {{\bf U}}^d + \nabla \phi \label{eq:decomposition}

where :math:`{{\bf U}}^d` is the divergence free portion of the velocity vector,
:math:`{{\bf U}}^\star`, and :math:`\phi` is a scalar. Taking the divergence of
Eq. ([eq:decomposition]), we have

.. math:: \nabla^2 \phi = \nabla \cdot {{\bf U}}^\star

(where we used :math:`\nabla \cdot {{\bf U}}^d = 0`).
With appropriate boundary conditions, this Poisson equation can be
solved for :math:`\phi`, and the final, divergence-free velocity can
be computed as

.. math:: {{\bf U}}^{n+1} = {{\bf U}}^\star - \nabla \phi

Because soundwaves are filtered, the timestep constraint now depends only
on :math:`|{{\bf U}}|`.

Extensions to variable-density incompressible
flows :cite:`bellMarcus:1992b` involve a slightly different
decomposition of the velocity field and, as a result, a slightly
different, variable-coefficient Poisson equation.
There are also a variety of different ways
to express what is being projected :cite:`almgren:bell:crutchfield`,
and different discretizations of the divergence and gradient operators
lead to slightly different mathematical properties of the methods
(leading to “approximate
projections” :cite:`almgrenBellSzymczak:1996`). Finally, for
second-order methods, two projections are typically done per timestep.
The first (the ‘MAC’ projection :cite:`bellColellaHowell:1991`)
operates on the half-time, edge-centered advective velocities, making
sure that they satisfy the divergence constraint. These advective
velocities are used to construct the fluxes through the interfaces to
advance the solution to the new time. The second/final projection
operates on the cell-centered velocities at the new time, again
enforcing the divergence constraint. The IAMR algorithm performs
both of these projections.

