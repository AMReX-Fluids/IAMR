.. _Chap:Tutorials:

Tutorials
=========

These tutorials provide problem setups to help familiarize you with running
IAMR. First we provide some useful pointers to relevant sections of this guide
where you can get more detailed information on how IAMR works,
then we give problem descriptions.

There are a large number of options which can be specified from the inputs.
Most of the options which can be specified from the
inputs files are left to their default values in the sample calculations.
For more thorough information on available options see :ref:`Chap:RuntimeOptions`.

As a starting point for code changes: the initial data are specified
in ``/Source/prob/prob_init.cpp``.  We have included several different
subroutines for defining different initial data.  The variable "probtype"
set in the inputs file selects between these subroutines.  You may also,
of course, write your own by modifying ``prob_init.cpp``.
For information on how to set up your own problem, see :ref:`Chap:ProblemSetup`.

For AMR, the criteria used for error estimation used in tagging can be
specified in the inputs file for some of the most common choices or more specialized
choices can be defined in ``NS_error.cpp`` (see :ref:`sec:tagging` Section for more information).
The estimation can depend on any or
all of the state variables or derived quantities, which are specified in
``NS_setup.cpp`` and defined in ``NS_derive.cpp``.

This code is a research code, and is continually being modified and
improved  as our needs evolve.  Because the list of options is so
extensive, and the updates relatively frequent, we heartily encourage
you to contact us directly by opening an issue on github
if you would like help modifying the code supplied here for your
own calculations.  There is extensive but undocumented capability.

That said, we welcome your comments, suggestions, and other feedback.
Again, please open an issue on github.


Problem Descriptions
--------------------

Each problem is its own directory within ``/Tutorials``. Many contain both 2D
and 3D inputs files. Note that running a 2D inputs file requires building
and using a 2D executable. Similarly, 3D inputs files require a 3D executable.
Dimensionality is set in the GNUmakefile.


Embedded Boundaries:
**************************

* **DoubleShearLayer**:
  This test case is a fully periodic double shear layer in a constant
  density fluid. A blob of tracer is passively advected with the flow.
  Contains an embedded boundary with AMR around the tracer and EB.


* **FlowPastCylinder**:
  Constant density flow around a cylinder. A blob of tracer is passively
  advected with the flow.
  ``inputs.2d.flow_past_cylinder-x`` and ``inputs.3d.flow_past_cylinder-x`` use AMR
  around the tracer and EB.


Non-EB:
**************************

* **Bubble**: This test case is a falling drop in a closed box with a density of
  twice the surrounding medium.  The calculation allows two levels of factor
  of 2 refinement. Here the refinement criteria are vorticity and the
  presence of the tracer, which coincides initially with the heavy drop.
  Also contains inputs file for using RZ coordinates.


* **Hotspot**: A hot bubble rising in closed box. Evolves a temperature field and uses
  a low Mach number constraint in place of incompressible. AMR refinement criteria is
  based on the temperature, but only 2D uses AMR by default. The 3D setup features an
  open top (outflow BC) and demonstates how to allow refinement at the outflow (which
  is turned off in IAMR by default).


* **RayleighTaylor**: Rayleigh-Taylor instability; heavy fluid on top of light fluid with gravity.
  AMR refinement is based on vorticity.


* **ConvectedVortex**: Euler vortex in isentropic flow. Analytic solution is translation of the
  initial conditions based on propagation speed and simulation time.
  There are several references, for example
  Spiegel, Seth & Huynh, H.T. & DeBonis, James. (2015). A Survey of the Isentropic Euler Vortex Problem using High-Order Methods. 10.2514/6.2015-2444.


* **LidDrivenCavity**: Lid-driven cavity ia a popular test case for incompressible, viscous flow.
  No-slip conditions are enforced on all walls, but the top ("lid") has a
  prescribed, constant velocity. Velocity and density are normalised so that
  changing the viscosity coefficient :math:`\mu` alters the Reynolds number according to
  :math:`Re = 1 / \mu,`
  where :math:`\mu` is the viscosity coefficient.
  [Reference: Ghia et al., "High-Re solutions for incompressible flow using the
  Navier-Stokes equations and a multigrid method", J. Comp. Phys. (1982)]


* **Poiseuille**: Simple Poiseuille flow in a square duct. The constant pressure gradient
  :math:`p_0` is enforced by setting the gravity parameter. The analytical solution is

  .. math::
     u = p_0 y (L - y) / (2 \mu).

  We use :math:`p_0 = \mu = L = 1`.


* **Euler**: The test case is a "vortex tube" in a constant density fluid
  in a triply periodic geometry.  The refinement criteria are the
  presence of a tracer and the magnitude of vorticity.


* **TaylorGreen**: This case is an unsteady viscous benchmark for which the
  exact solution in 2D is

  .. math::
    u(x,y,t) &=  && V_0 Sin(2\pi x) Cos(2\pi y) Cos(2\pi z) \exp(-2 (2\pi)^2 \nu t) \\
    v(x,y,t) &= -&& V_0 Cos(2\pi x) Sin(2\pi y) Cos(2\pi z) \exp(-2 (2\pi)^2 \nu t) \\
    p(x,y,t) &= -&& \rho_0 V_0^2 \{Cos(4 \pi x) + Cos(4 \pi y)\} \exp(-4 (2\pi)^2 \nu t) / 4

  In ``TaylorGreen/benchmarks``, there is a tool, ViscBench2d.cpp, that reads a plot file and compares the solution against this exact solution. This benchmark was originally derived by G.I. Taylor (Phil. Mag., Vol. 46, No. 274, pp. 671-674, 1923) and Ethier & Steinman (Intl. J. Num. Meth. Fluids, Vol. 19, pp. 369-375, 1994) give the pressure field.

  In 3D, the problem is initialized with

  .. math::
    u(x,y,z) &= &&  V_0 Sin(2\pi x) Cos(2\pi y) Cos(2\pi z) \\
    v(x,y,z) &= -&& V_0 Cos(2\pi x) Sin(2\pi y) Cos(2\pi z) \\
    w(x,y,z) &= &&  0.0 \\
    p(x,y,t) &= -&& \rho_0 V_0^2 \{2 + Cos(4 \pi z)\}\{Cos(4 \pi x) + Cos(4 \pi y)\} \exp(-4 (2\pi)^2 \nu t) / 16


* **HIT**: Homogeneous isentropic forced turbulence with constant density.
  This demonstrates defining a new forcing function by using a local edited
  version of ``NS_getForce.cpp``. IAMR's make system is automatically configured
  to select any local versions of files and ignore the corresponding verions in
  ``IAMR/Source``. This problem is 3D only.


* **Particles**: Particles in a double shear layer. Uses 2 levels of refinement
  and fixed grids. With fixed grids, a grid file (called ``fixed_grids_ml`` here) is
  used to define the grids for levels >= 1.


