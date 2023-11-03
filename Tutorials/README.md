These tutorials provide problem setups to help familiarize you with running
IAMR. Below we provide a quick overview. For more information,
please see IAMR's documentation at https://amrex-fluids.github.io/IAMR/
or open an issue on github.

There are a large number of options which can be specified from the inputs,
as well as significant changes which can be made with relatively small
changes to the program files, especially /Source/prob/prob_init.cpp
(see IAMR's documenation for more thorough information).
Most of the options which can be specified from the
inputs files are left to their default values in the sample calculations.

As a starting point for code changes: the initial data are specified
in /Source/prob/prob_init.cpp.  We have included several different
subroutines for defining different initial data.  The variable "probtype"
set in the inputs file selects between these subroutines.  You may also,
of course, write your own by modifying prob_init.cpp.

For AMR, the criteria used for error estimation used in tagging can be
specified in the inputs file for some of the most common choices
(see IAMR's documentation), or more specialized
choices can be defined in NS_error.cpp. The estimation can depend on any or
all of the state variables or derived quantities, which are specified in
NS_setup.cpp and defined in NS_derive.cpp.

This code is a research code, and is continually being modified and
improved  as our needs evolve.  Because the list of options is so
extensive, and the updates relatively frequent, we heartily encourage
you to contact us directly (open an issue on github or email to
ASAlmgren@lbl.gov) if you want to modify the code supplied here for your
own calculations.  There is extensive but undocumented capability.

That said, we welcome your comments, suggestions, and other feedback.
Again, you can open an issue on github or email to ASALmgren@lbl.gov.


PROBLEM DESCRIPTIONS

Each problem is its own directory within /Tutorials. Many contain both 2D
and 3D inputs files. Note that running a 2D inputs file requires building
and using a 2D executable. Similarly, 3D inputs files require a 3D executable.
Dimensionality is set in the GNUmakefile.


EMBEDDED BOUNDARIES:

**************************
Double Shear Layer

This test case is a fully periodic double shear layer in a constant
density fluid. A blob of tracer is passively advected with the flow.
Contains an embedded boundary with AMR around the tracer and EB.


**************************
Flow past cylinder

Constant density flow around a cylinder. A blob of tracer is passively
advected with the flow.
inputs.2d.flow_past_cylinder-x and inputs.3d.flow_past_cylinder-x use AMR.


NON-EB:

**************************
Bubble

This test case is a falling drop in a closed box with a density of
twice the surrounding medium.  The calculation allows two levels of factor
of 2 refinement. Here the refinement criteria are vorticity and the
presence of the tracer, which coincides initially with the heavy drop.
Also contains inputs file for using RZ coordinates.


**************************
Hotspot

A hot bubble rising in closed box. Uses a low Mach number constraint
in place of incompressible.

Average? LES?


**************************
Rayleigh-Taylor

Rayleigh-Taylor instability; heavy fluid on top of light fluid with gravity.


**************************
Convected Vortex

Euler vortex in isentropic flow. Analytic solution is translation of the
initial conditions based on propagation speed and simulation time.
There are several references, for example
Spiegel, Seth & Huynh, H.T. & DeBonis, James. (2015). A Survey of the Isentropic Euler Vortex Problem using High-Order Methods. 10.2514/6.2015-2444.


**************************
Lid Driven Cavity

Lid-driven cavity ia a popular test case for incompressible, viscous flow.
No-slip conditions are enforced on all walls, but the top ("lid") has a
prescribed, constant velocity. Velocity and density are normalised so that
changing the viscosity coefficient mu alters the Reynolds number according to
        Re = 1 / mu,
where mu is the viscosity coefficient.

Reference: Ghia et al., "High-Re solutions for incompressible flow using the
Navier-Stokes equations and a multigrid method", J. Comp. Phys. (1982)


**************************
Poiseuille

Simple Poiseuille flow in a square duct. The constant pressure gradient p0 is
enforced by setting the gravity parameter.

Analytical solution:
        u = p0 y (L - y) / (2 mu)

We use p0 = mu = L = 1.


**************************
Euler

The test case is a "vortex tube" in a constant density fluid
in a triply periodic geometry.  The refinement criteria are the
presence of a tracer and the magnitude of vorticity.


**************************
Taylor-Green vortex

This case is an unsteady viscous benchmark for which the
exact solution in 2D is
    u(x,y,t) =   V_0 Sin(2Pi x) Cos(2Pi y) Cos(2Pi z) Exp(-2 (2Pi)^2 Nu t)
    v(x,y,t) = - V_0 Cos(2Pi x) Sin(2Pi y) Cos(2Pi z) Exp(-2 (2Pi)^2 Nu t)
    p(x,y,t) = - rho_0 V_0^2 {Cos(4 Pi x) + Cos(4 Pi y)} Exp(-4 (2Pi)^2 Nu t) / 4
In TaylorGreen/benchmarks, there is a tool, ViscBench2d.cpp, that reads a
plot file and compares the solution against this exact solution.
This benchmark was originally derived by G.I. Taylor (Phil. Mag.,
Vol. 46, No. 274, pp. 671-674, 1923) and Ethier and Steinman
(Intl. J. Num. Meth. Fluids, Vol. 19, pp. 369-375, 1994) give
the pressure field.

In 3D, the problem is initialized with
    u(x,y,z) =   V_0 Sin(2Pi x) Cos(2Pi y) Cos(2Pi z)
    v(x,y,z) = - V_0 Cos(2Pi x) Sin(2Pi y) Cos(2Pi z)
    w(x,y,z) =   0.0
    p(x,y,t) = - rho_0 V_0^2 {2 + Cos(4 Pi z)}{Cos(4 Pi x) + Cos(4 Pi y)} Exp(-4 (2Pi)^2 Nu t) / 16


**************************
HIT

Homogeneous isentropic forced turbulence with constant density.


**************************
Particles

Put particles in a double shear layer. Uses 2 levels of refinement
and fixed grids.


