.. _Chap:InputsLoadBalancing:

Gridding and Load Balancing Inputs
==================================

The details of the regridding strategy are described in other sections
(see :ref:`sec:gridCreation` and :ref:`chap:Amr`);
here we cover how the input parameters can control the gridding.

As described later, the user defines error estimation functions to tag individual
cells at a given level if they need refinement. This list of tagged cells is
sent to a grid generation routine, which uses the Berger-Rigoutsos algorithm :cite:`br-refine`
to create rectangular grids that contain the tagged cells.

The following inputs must be preceded by "amr" and determine how we create the grids and how often we regrid.
(Additional information can also be found in AMReX documentation at :ref:`amrex:ss:amrcore`.)

+----------------------+-----------------------------------------------------------------------+-------------+-----------+
|                      | Description                                                           |   Type      | Default   |
+======================+=======================================================================+=============+===========+
| regrid_int           | How often to regrid (in number of steps at level 0)                   |   Int       |    -1     |
|                      | if regrid_int = -1 then no regridding will occur                      |             |           |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| regrid_on_restart    | Should we regrid immediately after restarting?                        |    Int      |  0        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| max_grid_size_x      | Maximum number of cells in each grid in x-direction, for all levels   |    Int      | 32        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| max_grid_size_y      | Maximum number of cells in each grid in y-direction, for all levels   |    Int      | 32        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| max_grid_size_z      | Maximum number of cells in each grid in z-direction, for all levels   |    Int      | 32        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| max_grid_size        | Maximum number of cells in each grid in all directions.               |    Int      | 32        |
|                      | Specify multiple values to give levels a different max_grid_size      |             |           |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| blocking_factor_x    | Each grid must be divisible by blocking_factor_x in x-direction       |    Int      |  8        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| blocking_factor_y    | Each grid must be divisible by blocking_factor_y in y-direction       |    Int      |  8        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| blocking_factor_z    | Each grid must be divisible by blocking_factor_z in z-direction       |    Int      |  8        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| blocking_factor      | Each grid must be divisible by blocking_factor in all directions.     |    Int      |  8        |
|                      | Specify multiple values to give levels a different blocking_factor    |             |           |
|                      | Must be a power of 2 at every level and the domain size must be a     |             |           |
|                      | multiple of blocking_factor at level 0.                               |             |           |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| grid_eff             | grid efficiency (must be between 0 and 1)                             |    Real     |  0.7      |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| n_error_buf          | radius of additional tagging around already tagged cells              |    Int      |  1        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| refine_grid_layout   | refine grids more if # of processors :math:`>` # of grids             |    Int      |  1        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| regrid_file          | Name of file from which to read the grids, if specifying fixed grids  |    Text     |  None     |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+

Note that if regrid_file is set (e.g. amr.regrid_file = fixed_grids), then the
list of grids at each fine level are read in from this file during the gridding
procedure. These grids must not violate the amr.max_grid_size criterion. The rest of the gridding procedure
will not occur if amr.regrid_file is set.

Note also that amr.ref_ratio, amr.n_error_buf, amr.max_grid_size and
amr.blocking_factor can be read in as a single value which is
assigned to every level, or as multiple values, one for each level.

As an example, consider:

::

    amr.regrid_int = 2 2
    amr.grid_eff = 0.9
    amr.max_grid_size = 64 
    amr.blocking_factor = 32

tells the code to regrid every 2 steps. Thus in this example, new
level 1 grids will be created every 2 level-0 time steps, and new
level 2 grids will be created every 2 level-1 time steps.
If amr.regrid_int :math:`<` 0 for any level, then regridding starting at that
level will be disabled. If amr.regrid_int = -1 only, then we
never regrid for any level. Note that this is not compatible with amr.regrid_on_restart = 1.

The grid efficiency, amr.grid_eff, means that during the grid
creation process, at least 90% of the cells in each grid at the level
at which the grid creation occurs must be tagged cells. A higher
grid efficiency means fewer cells at higher levels, but may result
in the production of lots of small grids, which have inefficient cache
and OpenMP performance and higher communication costs.

The amr.max_grid_size parameter means that the final grids
will be no longer than 64 cells on a side at every level.
Alternately, we could specify a value for each level of refinement as:
amr.max_grid_size = 64 32 16, in which case our final grids
will be no longer than 64 cells on a side at level 0, 32 cells on a
side at level 1, and 16 cells on a side at level 2. The amr.blocking_factor
means that all of the final grids will be multiples of 32 at all levels.
Again, this can be specified on a level-by-level basis, like
amr.blocking_factor = 32 16 8, in which case the
dimensions of all the final grids will be multiples of 32
at level 0, multiples of 16 at level 1, and multiples of 8 at level 2.



The following inputs must be preceded by "fabarray_mfiter" and determine how we create the logical tiles
(for more information on tiling see :ref:`sec:tiling`):

+----------------------+-----------------------------------------------------------------------+----------+-------------+
|                      | Description                                                           | Type     | Default     |
+======================+=======================================================================+==========+=============+
| tile_size            | Maximum number of cells in each direction for (logical) tiles         | IntVect  | 1024000     |
|                      |        (3D CPU-only)                                                  |          | 1024000,8,8 |
+----------------------+-----------------------------------------------------------------------+----------+-------------+



Getting good performance
~~~~~~~~~~~~~~~~~~~~~~~~

These parameters can have a large impact on the performance
of IAMR, so taking the time to experiment with is worth the effort.
Having grids that are large enough to coarsen multiple levels in a
V-cycle is essential for good multigrid performance.

How grids are created
~~~~~~~~~~~~~~~~~~~~~

The gridding algorithm proceeds in this order:

#. Grids are created using the Berger-Rigoutsos clustering algorithm
   modified to ensure that all new fine grids are divisible by amr.blocking_factor.

#. Next, the grid list is chopped up if any grids are larger than max_grid_size.
   Note that because amr.max_grid_size is a multiple of amr.blocking_factor the amr.blocking_factor criterion is
   still satisfied.

#. Next, if amr.refine_grid_layout = 1 and there are more processors than grids, and
   if amr.max_grid_size / 2 is a multiple of amr.blocking_factor,
   then the grids will be redefined, at each level independently, so that
   the maximum length of a grid at level :math:`\ell`, in any dimension, is
   amr.max_grid_size[:math:`\ell`] / 2.

#. Finally, if amr.refine_grid_layout = 1, and there are still more processors
   than grids, and if amr.max_grid_size / 4 is a multiple of amr.blocking_factor, then the grids will be redefined, at each level
   independently, so that the maximum length of a grid at level :math:`\ell`,
   in any dimension, is amr.max_grid_size[:math:`\ell`] / 4.
