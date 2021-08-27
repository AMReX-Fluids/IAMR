Running Simulations
===================

Physics Parameters
-------------------

The following inputs must be preceded by "ns." For information on units, see :ref:`sec:units`.

+----------------------+-------------------------------------------------------------------------+----------+-----------+
|                      | Description                                                             |   Type   | Default   |
+======================+=========================================================================+==========+===========+
| gravity              | Gravity, taken to be in the -y direction for 2d and -z direction in 3d  |  Real    |  0        |
+----------------------+-------------------------------------------------------------------------+----------+-----------+


.. _sec:InputsTimeStepping:

Time Stepping
-------------

The first three inputs below do not take a prefix.  Note that the first two are both specified, both criteria
are used and the simulation still stop when the first criterion is hit.  
The simulation will stop when either the number of steps reaches max_step or time reaches stop_time.

+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|                      | Description                                                           |   Type      | Default      |
+======================+=======================================================================+=============+==============+
| max_step             | Maximum number of time steps to take                                  |    Int      |  -1          |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| stop_time            | Maximum time to reach                                                 |    Real     | -1.0         |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+

The inputs below must be preceded by "ns."  

+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|                      | Description                                                           |   Type      | Default      |
+======================+=======================================================================+=============+==============+
| fixed_dt             | Value of fixed dt if > 0                                              |    Real     |   -1.        |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| cfl                  | CFL constraint (dt < cfl * dx / u) if fixed_dt not > 0                |    Real     |   0.5        |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+

  * If you want to fix the dt, simply set :cpp:`ns.fixed_dt = XXX` and the fluid time
    step will always be that number. 

  * If you want to let the code determine the appropriate time step using the advective CFL
    condition, then set :cpp:`ns.cfl = 0.7` for example, and the fluid time step will
    be computed to be dt = 0.7 * dx / max(vel).

  * Note that the cfl defaults to 0.5 so it does not have to be set in the inputs file. If neither
    :cpp:`ns.cfl` nor :cpp:`fixed_dt` is set, then default value of cfl will be used.
    If :cpp:`ns.fixed_dt` is set, then it will override the cfl option whether 
    :cpp:`ns.cfl` is set or not.


Output Options
--------------
	 
.. _Chap:InputsPlotfiles:

Plotfiles
~~~~~~~~~

The following inputs must be preceded by "amr." and control frequency and naming of plotfile generation as well
as whether the EB geometry should be written out.

+---------------------+-----------------------------------------------------------------------+-------------+-----------+
|                     | Description                                                           |   Type      | Default   |
+=====================+=======================================================================+=============+===========+
| plot_int            | Frequency of plotfile output;                                         |    Int      | -1        |
|                     | if -1 then no plotfiles will be written at this frequency             |             |           |
| plot_per            | Time period of plotfile output (approximate); does not modify dt      |    Real     | -1        |
|                     | if -1 then no plotfiles will be written at this frequency             |             |           |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plotfile_on_restart | Should we write a plotfile when we restart (only used if plot_int>0)  |   Bool      | False     |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plot_file           | Prefix to use for plotfile output                                     |  String     | plt       |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plot_vars           | State variables to include in plotfile                                |  String     | ALL       |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| derive_plot_vars    | Derived variables to include in plotfile                              |  String     | NONE      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+

	 
.. _Chap:InputsCheckpoint:

Checkpointing and Restarting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following inputs must be preceded by "amr." and control checkpoint/restart.

+------------------+-----------------------------------------------------------------------+-------------+-----------+
|                  | Description                                                           |   Type      | Default   |
+==================+=======================================================================+=============+===========+
| restart          | If present, then the name of file to restart from                     |    String   | None      |
+------------------+-----------------------------------------------------------------------+-------------+-----------+
| check_int        | Frequency of checkpoint output;                                       |    Int      | -1        |
|                  | if -1 then no checkpoints will be written                             |             |           |
+------------------+-----------------------------------------------------------------------+-------------+-----------+
| check_file       | Prefix to use for checkpoint output                                   |  String     | chk       |
+------------------+-----------------------------------------------------------------------+-------------+-----------+


Particles Output
~~~~~~~~~~~~~~~~

Checkpoint Files
^^^^^^^^^^^^^^^^

The particle positions and velocities are stored in a binary file in each checkpoint directory.
This format is designed for being read by the code at restart rather than for diagnostics.

Plot Files
^^^^^^^^^^

If ``particles.write\_in\_plotfile = 1`` in the inputs file
then the particle positions and velocities will be written in a binary file in each plotfile directory.

In addition, we can also
visualize the particle locations as represented on the grid. The “derived quantity”
``particle\_count`` represents the number of particles in a grid cell.
To add it to plotfiles, set
``amr.derive\_plot\_vars = particle\_count``
in the inputs file

ASCII Particle Files
^^^^^^^^^^^^^^^^^^^^

To generate an ASCII file containing the particle positions and velocities,
one needs to restart from a checkpoint file from a run with particles, but one doesn’t need to run any steps.
For example, if chk00350 exists, then one can set:

::
   
   amr.restart = chk00350
   max\_step = 350
   particles.particle\_output\_file = *particle\_output*

This tells the code to restart from chk00350, not to take any further time steps, and to write an ASCII-format
file called *particle\_output*.
This file has the same format as the ASCII input file:

::
   
   number of particles
   x y z


.. _Chap:InputsLoadBalancing:

Gridding and Load Balancing Inputs
----------------------------------

The details of the regridding strategy are described in :ref:`sec:gridCreation`;
here we cover how the input parameters can control the gridding.

These parameters can have a large impact on the performance
of IAMR, so taking the time to experiment with is worth the effort.
Having grids that are large enough to coarsen multiple levels in a
V-cycle is essential for good multigrid performance.


Gridding
~~~~~~~~

The following inputs must be preceded by "amr." and determine how we create the grids and how often we regrid.
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

Note that if regrid_file is set (e.g. ``amr.regrid_file = fixed_grids``), then the
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


.. _sec:tilingInputs:

Tiling
~~~~~~

For details on IAMR's approach to tiling see :ref:`Chap:Parallel`.

The following inputs determine how we create the logical tiles and must be preceded by "fabarray_mfiter." :

+----------------------+-----------------------------------------------------------------------+----------+-------------+
|                      | Description                                                           | Type     | Default     |
+======================+=======================================================================+==========+=============+
| tile_size            | Maximum number of cells in each direction for (logical) tiles         | IntVect  | 1024000     |
|                      |        (3D CPU-only)                                                  |          | 1024000,8,8 |
+----------------------+-----------------------------------------------------------------------+----------+-------------+


.. _Chap:InputsVerbosity:

Verbosity
---------

Different classes control their own verbosity. In some cases, values > 1 will generate additional verbosity.
Here is some of the more frequently used options:

+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|                      | Description                                                           |   Type      | Default      |
+======================+=======================================================================+=============+==============+
| ns.v                 |  Verbosity in IAMR routines                                           |    Int      |   0          |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| particles.pverbose   |  Verbosity in particle routines                                       |    Int      |   0          |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| nodal_proj.verbose   |  Verbosity in nodal projection                                        |    Int      |   0          |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| mac_proj.verbose     |  Verbosity in MAC projection                                          |    Int      |   0          |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+

.. _Chap:InputsMultigrid:

Multigrid Inputs
----------------

IAMR uses AMReX's multigrid functionality to perform the nodal projection (which enures the cell-centered
velocity field obeys the constraint), the MAC projection (which ensures that the edge-based velocity field
used in advection obeys the constraint), and the diffusive solves.

Here we go over some inputs parameters that can be used to control these solves. For more information
on AMReX's linear solvers, see :ref:`amrex:Chap:LinearSolvers`


Nodal Projection
~~~~~~~~~~~~~~~~

These control the nodal projection and must be preceded by "nodal_proj.": 

+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
|                         |  Description                                                          |   Type      | Default      |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| verbose                 |  Verbosity of multigrid solver in nodal projection                    |    Int      |   0          |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| bottom_verbose          |  Verbosity of BiCGStab solver in nodal projection                     |    Int      |   0          |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| rtol                    |  Relative tolerance in nodal projection                               |    Real     |   1.e-11     |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| atol                    |  Absolute tolerance in nodal projection                               |    Real     |   1.e-14     |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| maxiter                 |  Maximum number of iterations in the nodal projection                 |    Int      |   100        |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| bottom_maxiter          |  Maximum number of iterations in the nodal projection                 |    Int      |   100        |
|                         |  bottom solver if using bicg, cg, bicgcg or cgbicg                    |             |              |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| mg_max_coarsening_level |  Maximum number of coarser levels to allowin the nodal projection     |    Int      |   100        |
|                         |  If set to 0, the bottom solver will be called at the current level   |             |              |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| bottom_solver           |  Which bottom solver to use in the nodal projection                   |  String     |   bicgcg     |
|                         |  Options are bicgcg, bicgstab, cg, cgbicg, smoother or hypre          |             |              |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+

MAC Projection
~~~~~~~~~~~~~~

These control the MAC projection and must be preceded by "mac_proj.":

+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
|                         | Description                                                           |   Type      | Default      |
+=========================+=======================================================================+=============+==============+
| verbose                 |  Verbosity of multigrid solver in MAC projection                      |    Int      |   0          |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| bottom_verbose          |  Verbosity of BiCGStab solver in MAC projection                       |    Int      |   0          |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| rtol                    |  Relative tolerance in MAC projection                                 |    Real     |   1.e-11     |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| atol                    |  Absolute tolerance in MAC projection                                 |    Real     |   1.e-14     |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| maxiter                 |  Maximum number of iterations in the MAC projection                   |    Int      |   200        |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| bottom_maxiter          |  Maximum number of iterations in the MAC projection                   |    Int      |   200        |
|                         |  bottom solver if using bicg, cg, bicgcg or cgbicg                    |             |              |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| mg_max_coarsening_level |  Maximum number of coarser levels to allow in the MAC projection      |    Int      |   100        |
|                         |  If set to 0, the bottom solver will be called at the current level   |             |              |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| bottom_solver           |  Which bottom solver to use in the MAC projection                     |  String     |   bicgcg     |
|                         |  Options are bicgcg, bicgstab, cg, cgbicg, smoother or hypre          |             |              |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+

Viscous and Diffusive Solve
~~~~~~~~~~~~~~~~~~~~~~~~~~~

These control the diffusion solver and must be preceded by "diffusion.":

+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
|                         | Description                                                           |   Type      | Default      |
+=========================+=======================================================================+=============+==============+
| v                       |  Verbosity of linear solver for diffusion solve                       |    Int      |   0          |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+


.. _Chap:InputsInitialization:

Initializing the Calculation
----------------------------

These options determine how we initialize the data for the calculation. The data initialization process
ensures that the initial state is consistent with the constraint, and if applicable, the various AMR levels
are consistent with eachother. 

The following inputs must be preceded by "ns." 

+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|                      | Description                                                           |   Type      | Default      |
+======================+=======================================================================+=============+==============+
| do_init_proj         | Do the initial projections? False is primarily for debugging.         |    Bool     |  True        |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| init_iter            | How many pressure iterations before starting the first timestep       |  Int        |    3         |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| init_vel_iters       | How many projection iterations to ensure the velocity satisfies the   |  Int        |    3         |
|                      |  constraint. Set = 0 to skip this part of the initialization.         |             |              |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+


