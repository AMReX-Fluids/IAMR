
Included in the run directories /run2d, /run3d, /eb_run2d, /eb_run3d,
/run_2d_particles
are sample inputs files that supply problem-dependent paramters for the
calculations.
In the directories /eb_run2d and /eb_run3d are problems with embedded
boundaries, and /run_2d_particles provides an example of using pasively
advected particles. 
Files starting with "regtest" are inputs files for IAMR's
regression tests. IAMR/Test/README.md has more information about the
regression tests, including how to set them up for your own IAMR fork.
Note that the regression tests are designed to test the code, and
therefore might not always use the more physically meaningful parameters.

There are a large number of options which can be specified from the inputs,
as well as significant changes which can be made with relatively small
changes to the program files, especially /Source/prob/prob_init.cpp
(see IAMR/UsersGuide for more thurough information).
Most of the options which can be specified from the
inputs files are left to their default values in the sample calculations.

As a starting point for code changes: the initial data are specified 
in /Source/prob/prob_init.cpp.  We have included several different
subroutines for defining different initial data.  The variable "probtype,"
set in the inputs file selects between these subroutines.  You may also,
of course, write your own by modifying prob_init.cpp.

The criteria used for error estimation can be specified in the inputs file
for some of the most common choices (see IAMR/UsersGuide), or more specialized
choices can be defined in NS_error.cpp. The estimation can depend on any or
all of the state variables or derived quantities, which are specified in
NS_setup.cpp and defined in NS_derive.cpp.

This code is a research code, and is continually being modified and 
improved  as our needs evolve.   Because the list of options is so 
extensive, and the updates relatively frequent, we heartily encourage 
you to contact us directly (email to ASAlmgren@lbl.gov) if you want to
modify the code supplied here for your own calculations.   There is
extensive but undocumented capability.

That said, we welcome your comments, suggestions, and other feedback.
Again, email to ASALmgren@lbl.gov.



PROBLEM DESCRIPTIONS

**************************
Double Shear Layer

This test case is a fully periodic double shear layer in a constant
density fluid. A blob of tracer is passively advected with the flow.


**************************
Bubble

This test case is a falling drop in a closed box with a density of 
twice the surrounding medium.  The calculation allows two levels of factor
of 2 refinement. Here the refinement criteria are vorticity and the
presence of the tracer, which coincides initially with the heavy drop.  


**************************
Rayleigh-Taylor

Rayleigh-Taylor instability; heavy fluid on top of light fluid with gravity.


**************************
Tracer Advect

Advecting a blob of tracer material with constant velocity field.


**************************
Lid Driven Cavity

Lid-driven cavity, popular test case for incompressible, viscous flow. 
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
Flow past cylinder

Constant density flow around a cylinder. A blob of tracer is passively
advected with the flow.


**************************
Hotspot

Simple hot bubble rising in closed box.  


**************************
Euler

The test case is a "vortex tube" in a constant density fluid
in a triply periodic geometry.  The refinement criteria are the
presence of a tracer and the magnitude of vorticity.


**************************
Taylor-Green vortex

This case is an unsteady viscous benchmark for which the
exact solution in 2D is
    u(x,y,t) = - V_0 Cos(Pi x) Sin(Pi y) Exp(-2 Pi^2 Nu t)
    v(x,y,t) =   V_0 Sin(Pi x) Cos(Pi y) Exp(-2 Pi^2 Nu t)
    p(x,y,t) = - rho_0 V_0^2 {Cos(2 Pi x) + Cos(2 Pi y)} Exp(-4 Pi^2 Nu t) / 4
In Exec/benchmarks, there is a tool ViscBench2d.cpp that reads a
plot file and compares the solution against this exact solution.
This benchmark was originally derived by G.I. Taylor (Phil. Mag.,
Vol. 46, No. 274, pp. 671-674, 1923) and Ethier and Steinman
(Intl. J. Num. Meth. Fluids, Vol. 19, pp. 369-375, 1994) give
the pressure field.

In 3D the problem is initialized with
    u(x,y,z) =   V_0 Sin(2Pi x) Cos(2Pi y) Cos(2Pi z)
    v(x,y,z) = - V_0 Cos(2Pi x) Sin(2Pi y) Cos(2Pi z)
    w(x,y,z) =   0.0
