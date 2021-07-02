Problem Definition
==================

Resolution
----------

The grid resolution is specified by defining the resolution at the
coarsest level (level 0), the number of refinement levels, and
factor of refinement between levels. 
The following inputs must be preceded by "amr."

+-------------------+-------------------------------------------------------------------------+-------------+-----------+
|                   | Description                                                             |   Type      | Default   |
+===================+=========================================================================+=============+===========+
| n_cell            | Number of cells at level 0 in each coordinate direction                 | Int Int Int | None      |
+-------------------+-------------------------------------------------------------------------+-------------+-----------+
| max_level         | Maximum level of refinement allowed (0 when single-level)               |    Int      | None      |
+-------------------+-------------------------------------------------------------------------+-------------+-----------+
| ref_ratio         | Ratio of coarse to fine grid spacing between subsequent levels (2 or 4) |    Int      | None      |
+-------------------+-------------------------------------------------------------------------+-------------+-----------+

Some examples:

::

    amr.n_cell = 32 64 64

would define the domain to have 32 cells in the :math:`x`-direction, 64 cells
in the :math:`y`-direction, and 64 cells in the :math:`z`-direction *at the
coarsest level*. (If this line appears in a 2D inputs file then the
final number will be ignored.)

::

    amr.max_level = 2 

would allow a maximum of 2 refined levels in addition to the coarse
level. Note that these additional levels will only be created only if
the tagging criteria are such that cells are flagged as needing
refinement. The number of refined levels in a calculation must be
:math:`\leq` amr.max_level, but can change in time and need not
always be equal to amr.max_level.

::

    amr.ref_ratio = 2 4 

would set factor of 2 refinement between levels 0 and 1, and factor of 4
refinement between levels 1 and 2. Note that you must have at least
amr.max_level values of amr.ref_ratio (Additional values
may appear in that line and they will be ignored).

   
Problem Geometry
----------------

The following inputs must be preceded by "geometry."

+-----------------+-----------------------------------------------------------------------+-------------+-----------+
|                 | Description                                                           |   Type      | Default   |
+=================+=======================================================================+=============+===========+
| coord_sys       | 0 for Cartesian                                                       |   Int       |   0       |
+-----------------+-----------------------------------------------------------------------+-------------+-----------+
| is_periodic     | 1 for true, 0 for false (one value for each coordinate direction)     |   Ints      | 0 0 0     |
+-----------------+-----------------------------------------------------------------------+-------------+-----------+
| prob_lo         | Low corner of physical domain (physical not index space)              |   Reals     | None      |
+-----------------+-----------------------------------------------------------------------+-------------+-----------+
| prob_hi         | High corner of physical domain (physical not index space)             |   Reals     | None      |
+-----------------+-----------------------------------------------------------------------+-------------+-----------+

As an example, the following:

::

    geometry.prob_lo = 0 0 0
    geometry.prob_hi = 10.0 10.0 10.0
    geometry.coord_sys = 0 
    geometry.is_periodic = 0 1 0 

defines the domain to run from :math:`(0,0,0)` at the lower left to
:math:`(10,10,10)` at the upper right in physical space, specifies a
Cartesian geometry, and makes the domain periodic in the :math:`y`-direction
only.

The following inputs must be preceded by "ns."

+----------------------+-------------------------------------------------------------------------+----------+-----------+
|                      | Description                                                             |   Type   | Default   |
+======================+=========================================================================+==========+===========+
| gravity              | Gravity, taken to be in the -y direction for 2d and -z direction in 3d  |  Reals   |  0        |
+----------------------+-------------------------------------------------------------------------+----------+-----------+



Domain Boundary Conditions
--------------------------

IAMR provides two ways to specify boundary condition types: integer keys or keywords.
An inputs file must choose one style or the other, they cannot be mixed.
Here we provide examples of both styles. We then discuss how to provide Dirichlet
boundary conditions.

Boundary conditions with integer keys
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To specify boundary conditions with integer keys we use

-  ns.lo_bc: boundary type of each low face

-  ns.hi_bc: boundary type of each high face

The valid boundary types are:

+---------------------------+--------------------+
| 0 – Interior / Periodic   | 3 – Symmetry       |
+---------------------------+--------------------+
| 1 – Inflow                | 4 – Slip Wall      |
+---------------------------+--------------------+
| 2 – Outflow               | 5 – No Slip Wall   |
+---------------------------+--------------------+

Note: ns.lo_bc and ns.hi_bc must be consistent with
geometry.is_periodic—if the domain is periodic in a particular
direction then the low and high bc’s must be set to 0 for that direction.

As an example, the following:

::

    ns.lo_bc = 1 4 0 
    ns.hi_bc = 2 4 0 

    geometry.is_periodic = 0 0 1

would define a problem with inflow (1) in the low- :math:`x` direction,
outflow (2) in the high- :math:`x` direction, slip wall (4) on
the low and high :math:`y`-faces, and periodic in the :math:`z`-direction.

Boundary conditions with keywords
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To specify boundary conditions with keywords, we use the following options
preceded by “xlo”, “xhi”, “ylo”, “yhi”, “zlo”, and “zhi”:

+--------------------+---------------------------------------------------------------------------+-------------+-----------+
|                    | Description                                                               |   Type      | Default   |
+====================+===========================================================================+=============+===========+
| type               | Used to define boundary type. Available options include:                  |  String     |  None     |
|                    |                                                                           |             |           |
|                    | * 'po'  or 'pressure_outflow'                                             |             |           |
|                    | * 'mi'  or 'mass_inflow'                                                  |             |           |
|                    | * 'sw'  or 'slip_wall'                                                    |             |           |
|                    | * 'nsw' or 'no_slip_wall'                                                 |             |           |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+

To use the same example problem as above, the following:

::

    xlo.type = mi
    xhi.type = po
    ylo.type = sw
    yhi.type = sw

    geometry.is_periodic = 0 0 1

would define a problem with inflow in the low-\ :math:`x` direction,
outflow in the high-\ :math:`x` direction, slip wall on
the low and high :math:`y`-faces, and periodic in the :math:`z`-direction.
Note that no keyword is needed for a periodic boundary, here only the
specification in geometry.is\_periodic is needed.

Dirichlet Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IAMR provides the ability to specify constant Dirichlet BCs
in the inputs file. We use the following options
preceded by “xlo”, “xhi”, “ylo”, “yhi”, “zlo”, and “zhi”:

+--------------------+---------------------------------------------------------------------------+-------------+-----------+
|                    | Description                                                               |   Type      | Default   |
+====================+===========================================================================+=============+===========+
| velocity           | Sets boundary velocity for mass inflows                                   |    Real     |  None     |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| density            | Sets boundary density for mass inflows                                    |    Real     |  None     |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| tracer             | Sets boundary tracer concentration for mass inflows                       |    Real     |  None     |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| temp               | Sets temperature for mass inflows                                         |    Real     |  None     |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| pressure           | Sets boundary pressure for pressure inflows, outflows and mass inflows    |    Real     |  None     |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+


As an example,

::

    xlo.type                =   "mass_inflow"
    xlo.velocity            =   1.  0.  0.
    xlo.density             =   1.
    xlo.tracer              =   0.
    xlo.temp                =   1.

sets the boundary condtion type at the low x face to be an inflow with
xlo.type = “mass_inflow”.
Then xlo.velocity = 1. 0. 0. sets the inflow velocity,
xlo.density = 1. sets the inflow density,
xlo.tracer = 0. sets the inflow tracer value, and
xlo.temp = 1. sets the inflow temperature.
Another example, from the lid driven cavity problem setup, is

::

    ns.lo_bc                =  4 4 5
    ns.hi_bc                =  5 5 5

    # 0 = Interior/Periodic  3 = Symmetry
    # 1 = Inflow             4 = SlipWall
    # 2 = Outflow            5 = NoSlipWall

    # Boundary condition
    zhi.velocity            =   1.  0.  0.

Here, ns.hi_bc = 5 5 5 sets the boundary conditions on all high faces to
no slip walls.
zhi.velocity = 1. 0. 0. sets the wall at the high z face to be moving in the
x-direction.
Note that IAMR allows walls to move tangentially, but not in the normal direction.

Users can create more complex Dirichlet boundary condtions by writing
their own fill function in NS_bcfill.H, then using that function to create
an amrex::StateDescriptor::BndryFunc object and specifying which variables
will use it in NS_setup.cpp. More information on boundary conditions is in
section :ref:`sec:physicalBCs`.
