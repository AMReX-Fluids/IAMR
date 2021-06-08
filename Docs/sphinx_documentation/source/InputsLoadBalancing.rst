.. _Chap:InputsLoadBalancing:

Gridding and Load Balancing
===========================

The following inputs must be preceded by "amr" and determine how we create the grids and how often we regrid.

+----------------------+-----------------------------------------------------------------------+-------------+-----------+
|                      | Description                                                           |   Type      | Default   |
+======================+=======================================================================+=============+===========+
| regrid_int           | How often to regrid (in number of steps at level 0)                   |   Int       |    -1     |
|                      | if regrid_int = -1 then no regridding will occur                      |             |           |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| max_grid_size_x      | Maximum number of cells at level 0 in each grid in x-direction        |    Int      | 32        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| max_grid_size_y      | Maximum number of cells at level 0 in each grid in y-direction        |    Int      | 32        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| max_grid_size_z      | Maximum number of cells at level 0 in each grid in z-direction        |    Int      | 32        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| blocking_factor_x    | Each grid must be divisible by blocking_factor_x in x-direction       |    Int      |  8        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| blocking_factor_y    | Each grid must be divisible by blocking_factor_y in y-direction       |    Int      |  8        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+
| blocking_factor_z    | Each grid must be divisible by blocking_factor_z in z-direction       |    Int      |  8        |
+----------------------+-----------------------------------------------------------------------+-------------+-----------+

The following inputs must be preceded by "fabarray_mfiter" and determine how we create the logical tiles:

+----------------------+-----------------------------------------------------------------------+----------+-------------+
|                      | Description                                                           | Type     | Default     |
+======================+=======================================================================+==========+=============+
| tile_size            | Maximum number of cells in each direction for (logical) tiles         | IntVect  | 1024000     |
|                      |        (3D CPU-only)                                                  |          | 1024000,8,8 |
+----------------------+-----------------------------------------------------------------------+----------+-------------+

The following inputs must be preceded by "incflo" and determine how we load balance:

+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|                      | Description                                                           |   Type      | Default      |
+======================+=======================================================================+=============+==============+
| load_balance_type    | What strategy to use for load balancing                               |  String     | KnapSack     |
|                      | Options are "KnapSack"or "SFC"                                        |             |              |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| knapsack_weight_type | What weighting function to use if using Knapsack load balancing       |  String     | RunTimeCosts |
|                      | Options are "RunTimeCosts"                                            |             |              |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| knapsack_nmax        | Maximum number of grids per MPI process if using knapsack algorithm   |  Int        | 128          | 
+----------------------+-----------------------------------------------------------------------+-------------+--------------+

