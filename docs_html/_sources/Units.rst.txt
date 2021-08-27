Units
=====

For incompressible flow, IAMR supports any self-consistent units, as long as the length, time,
and mass are consistent with the viscosity and diffusivity.

For low Mach number flow with temperature variations, the temperature solve requires a specific heat capacity,
:math:`c_p`, which is currently hard-coded as :math:`1004.6 J kg^-1 K^-1`
(the specific heat capacity of dry air). Users can either keep this value, or change the specific heat capacity
in their own version of the source code.

