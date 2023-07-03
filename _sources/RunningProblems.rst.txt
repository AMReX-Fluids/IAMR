.. role:: cpp(code)

.. _Chap:RuntimeOptions:

Runtime Options
===================

.. _sec:InputsTimeStepping:

Time Stepping
-------------

The inputs below do not take a prefix.  Note that if the first two are both specified, both criteria
are used and the simulation still stop when the first criterion is hit.
The simulation will stop when either the number of steps reaches :cpp:`max_step` or time reaches :cpp:`stop_time`.
Note that if the code reaches :cpp:`stop_time` then the final time
step will be shortened so as to end exactly at :cpp:`stop_time`, not
past it.


+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|                      | Description                                                           |   Type      | Default      |
+======================+=======================================================================+=============+==============+
| max_step             | Maximum number of time steps to take                                  |    Int      |  -1          |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| stop_time            | Maximum time to reach                                                 |    Real     | -1.0         |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| stop_when_steady     | Stop when steady state is reached                                     |    Bool     | false        |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| steady_tol           | Specify tolerance to define steady state                              |    Real     | 1e-10        |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+

The inputs below must be preceded by "ns."

+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|                      | Description                                                           |   Type      | Default      |
+======================+=======================================================================+=============+==============+
| fixed_dt             | Value of fixed dt if > 0                                              |    Real     |   -1.        |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| cfl                  | CFL constraint (dt < cfl * dx / u) if fixed_dt not > 0                |    Real     |   0.5        |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| init_shrink          | Factor by which to shrink the initial time step                       |    Real     |   1.0        |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| max_change           | Factor by which time step can grow in subsequent steps                |    Real     |   1.1        |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| dt_cutoff            | Time step below which the simulation will abort                       |    Real     |   0.0        |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+

  * If you want to fix the dt, simply set :cpp:`ns.fixed_dt = XXX` to set the fluid time
    step at level 0. Note that :cpp:`init_shrink` can be used
    in conjunction with :cpp:`fixed_dt` and the initial time steps will be reduced.

  * If you want to let the code determine the appropriate time step using the advective CFL
    condition, then set :cpp:`ns.cfl = 0.7` for example, and the fluid time step will
    be computed to be dt = 0.7 * dx / max(vel).

  * Note that the cfl defaults to 0.5 so it does not have to be set in the inputs file. If neither
    :cpp:`ns.cfl` nor :cpp:`fixed_dt` is set, then default value of cfl will be used.
    If :cpp:`ns.fixed_dt` is set, then it will override the cfl option whether
    :cpp:`ns.cfl` is set or not.

As an example, consider:

.. code-block:: c++

    ns.cfl = 0.9
    ns.init_shrink = 0.01
    ns.change_max = 1.1
    ns.dt_cutoff = 1.e-20

This defines the :cpp:`cfl` parameter to be 0.9,
but sets (via :cpp:`init_shrink`) the first timestep we take
to be 1% of what it would be otherwise. This allows us to
ramp up to the hydrodynamic timestep at the start of a simulation.
The :cpp:`change_max` parameter restricts the timestep from increasing
by more than 10% over a coarse timestep. Note that the time step
can shrink by any factor; this only controls the extent to which it can grow.
The :cpp:`dt_cutoff` parameter will force the code to abort if the
timestep ever drops below :math:`10^{-20}`. This is a safety feature—if the
code hits such a small value, then something likely went wrong in the
simulation, and by aborting, you won’t burn through your entire allocation
before noticing that there is an issue.



Output Options
--------------

Here we give an overview of some of the output options for IAMR. For more detailed information on I/O, please see
:ref:`sec:ParallelIO` and AMReX's documentation :ref:`amrex:Chap:IO`

.. _sec:InputsPlotfiles:

Plotfiles
~~~~~~~~~

The following inputs must be preceded by "amr." and control frequency, naming, and content of the plotfiles.

+---------------------+-----------------------------------------------------------------------+-------------+-----------+
|                     | Description                                                           |   Type      | Default   |
+=====================+=======================================================================+=============+===========+
| plot_int            | Frequency of plotfile output;                                         |    Int      | -1        |
|                     | if -1 then no plotfiles will be written at this frequency             |             |           |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plot_per            | Time period of plotfile output (approximate); does not modify dt.     |    Real     | -1        |
|                     | If -1 then no plotfiles will be written at this frequency             |             |           |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plotfile_on_restart | Should we write a plotfile when we restart (only used if plot_int>0)  |   Bool      | False     |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plot_file           | Prefix to use for plotfile output                                     |  String     | plt       |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plot_vars           | State variables to include in plotfile                                |  String     | ALL       |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| derive_plot_vars    | Derived variables to include in plotfile                              |  String     | NONE      |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plot_nfiles         | Max number of files to use in parallel output if the number of        |    Int      |  64       |
|                     | processors exceeds this number, else the two are equal.               |             |           |
|                     | If -1, number of files is always set equal to number of processors.   |             |           |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+

The following inputs must be preceded by "ns."

+--------------------------+-----------------------------------------------------------------------+-------------+-----------+
|                          | Description                                                           |   Type      | Default   |
+==========================+=======================================================================+=============+===========+
| set_plot_coveredCell_val | For EB, set all covered cells? If 1, sets covered cells to zero.      |    Int      |  1        |
|                          | If 0, show covered cells without setting a value for the plotfile;    |             |           |
|                          | for debugging, as covered cells are not used in the algorithm.        |             |           |
+--------------------------+-----------------------------------------------------------------------+-------------+-----------+

Note:
 * ``amr.plot_per`` will write a plotfile at the first
   timestep whose ending time is past an integer multiple of this interval.
   The timestep is **not** modified to match this interval, so
   you won’t get a plotfile at exactly the time you requested, unless you've also set up ``ns.fixed_dt``.

 * You can specify both ``amr.plot_int`` or ``amr.plot_per``
   if you so desire; the code will print a warning in case you did this
   unintentionally. If both are non-negative values, you will get plotfiles
   at integer multiples of ``amr.plot_int`` timesteps and at
   multiples of ``amr.plot_per`` simulation time intervals.

 * All the options for ``amr.derive_plot_vars`` are kept in ``derive_lst`` in ``NS_setup.cpp``.
   Feel free to look at it and see what’s there.

As an example:

::

    amr.plot_file = plt_run
    amr.plot_int = 10

means that plotfiles (really directories) starting with the prefix
``plt_run`` will be generated every 10 level-0 time steps. The
directory names will be plt_run00000, plt_run00010, plt_run00020, etc.

If instead you specify

::

    amr.plot_file = plt_run
    amr.plot_per = 0.1

then plotfiles (really directories) starting with the prefix
``plt_run`` will be generated every 0.1 units of simulation time. The
directory names may be something like plt_run00000, plt_run00043, plt_run00061, etc, where
43 level-0 steps is the first point when simulation time :math:`>= 0.1`,
and 61 level-0 steps is the first point when simulation time :math:`>=0.2`, etc.



.. _sec:InputsCheckpoint:

Checkpointing and Restarting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following inputs must be preceded by "amr." and control checkpoint/restart.

+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
|                         | Description                                                           |   Type      | Default   |
+=========================+=======================================================================+=============+===========+
| restart                 | If present, then the name of file to restart from                     |    String   | None      |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
| check_int               | Frequency of checkpoint output by level 0 time steps;                 |    Int      | -1        |
|                         | if -1 then no checkpoints will be written                             |             |           |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
| check_per               | Frequency of checkpoint output by (approximate) simulation time;      |    Real     | -1.0      |
|                         | if -1 then no checkpoints will be written. (See note below)           |             |           |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
| checkpoint_files_output | Should we write checkpoint files? (Takes precident over check_int     |    Int      |  1        |
|                         | and check_per)                                                        |             |           |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
| check_file              | Prefix to use for checkpoint output                                   |  String     | chk       |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
| checkpoint_on_restart   | Write a checkpoint immediately after restarting? (If 1 yes, if 0 no)  |    Int      |  0        |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
| checkpoint_nfiles       | Max number of files to use in parallel output if the number of        |    Int      |  64       |
|                         | processors exceeds this number, else the two are equal.               |             |           |
|                         | If -1, number of files is always set equal to number of processors.   |             |           |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+

Note:
 * ``amr.check_per`` will write a checkpoint at the first
   timestep whose ending time is past an integer multiple of this interval.
   The timestep is **not** modified to match this interval, so
   you won’t get a checkpoint at exactly the time you requested, unless you've also set up ``ns.fixed_dt``.

 * You can specify both ``amr.check_int`` or ``amr.check_per``
   if you so desire; the code will print a warning in case you did this
   unintentionally. If both are non-negative values, you will get checkpoints
   at integer multiples of ``amr.check_int`` timesteps and at
   multiples of ``amr.check_per`` simulation time intervals.

As an example,

::

    amr.check_file = chk_run
    amr.check_int = 10

means that restart files (really directories) starting with the prefix
``chk_run`` will be generated every 10 level-0 time steps. The
directory names will be chk_run00000, chk_run00010, chk_run00020, etc.


To restart from chk_run00061, for example, then set

::

    amr.restart = chk_run00061



Particles Output
~~~~~~~~~~~~~~~~

Checkpoint Files
^^^^^^^^^^^^^^^^

The particle positions and velocities are stored in a binary file in each checkpoint directory.
This format is designed for being read by the code at restart rather than for diagnostics.

Plot Files
^^^^^^^^^^

The “derived quantity” ``particle_count`` represents the number of particles in a grid cell.
To visualize the particle locations as represented on the grid, add ``particle_count`` to the list
of derived quanties in ``amr.derive_plot_vars =`` in the inputs file.

If ``particles.write_in_plotfile = 1`` in the inputs file,
then the particle positions and velocities will be written in a binary file in each plotfile directory.
This allows the use of the AMReX tools such as the particle comparison tool found in ``amrex/Tools/Postprocessing/C_Src/``,
and/or ``amrex/Tools/Py_util/amrex_particles_to_vtp`` to generate a vtp file you can open with ParaView.


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


Log Files
~~~~~~~~~~~~~~~~

The following inputs must be preceeded by "amr."

+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
|                         | Description                                                           |   Type      | Default   |
+=========================+=======================================================================+=============+===========+
| grid_log                | If present, then the name of file to write grids; otherwise nothing   |    String   |           |
|                         | written.                                                              |             |           |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
| run_log                 | If present, then the name of file to write run output; otherwise      |    String   |           |
|                         | nothing written.                                                      |             |           |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
| run_log_terse           | If present, then the name of file to write number of steps, time,     |    String   |           |
|                         | and dt; otherwise nothing written.                                    |             |           |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+

As an example:

::

    amr.grid_log = grdlog
    amr.run_log = runlog

Every time the code regrids it prints a list of grids at all relevant
levels. Here the code will write these grids lists into the file grdlog. Additionally, every time step the code prints certain
statements to the screen, such as (if amr.v = 1):

::

    STEP = 1 TIME = 1.91717746 DT = 1.91717746
    PLOTFILE: file = plt00001

The run_log option will write these statements into the file *runlog* as well.

Terser output can be obtained via:

::

    amr.run_log_terse = runlogterse

This file, runlogterse differs from runlog, in that it
only contains lines of the form

::

    10  0.2  0.005

in which “10” is the number of steps taken, “0.2” is the
simulation time, and “0.005” is the level-0 time step. This file
can be plotted very easily to monitor the time step.



.. _sec:InputsLoadBalancing:

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

+----------------------+-----------------------------------------------------------------------+----------+---------------+
|                      | Description                                                           | Type     | Default       |
+======================+=======================================================================+==========+===============+
| tile_size            | Maximum number of cells in each direction for (logical) tiles.        | IntVect  | 1024000       |
|                      | (Default for 3D CPU-only)                                             |          | (1024000,8,8) |
+----------------------+-----------------------------------------------------------------------+----------+---------------+



.. _sec:InputsVerbosity:

Verbosity
---------

Different classes control their own verbosity. In some cases, values > 1 will generate additional verbosity.
Here is some of the more frequently used options:

+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|                      | Description                                                           |   Type      | Default      |
+======================+=======================================================================+=============+==============+
| ns.v                 |  Verbosity in IAMR routines                                           |    Int      |   0          |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| ns.amr               |  Verbosity in AMR routines                                            |    Int      |   0          |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| particles.verbose    |  Verbosity in particle routines                                       |    Int      |   0          |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| nodal_proj.verbose   |  Verbosity in nodal projection                                        |    Int      |   0          |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| mac_proj.verbose     |  Verbosity in MAC projection                                          |    Int      |   0          |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+

.. _sec:InputsMultigrid:

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
| verbose                 |  Verbosity in nodal projection                                        |    Int      |   0          |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| bottom_verbose          |  Verbosity of the bottom solver in nodal projection                   |    Int      |   0          |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| proj_tol                |  Relative tolerance in nodal projection                               |    Real     |   1.e-12     |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| sync_tol                |  Relative tolerance in sync projection                                |    Real     |   1.e-10     |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| proj_abs_tol            |  Absolute tolerance in nodal & sync projections                       |    Real     |   1.e-16     |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| maxiter                 |  Maximum number of iterations in the nodal projection                 |    Int      |   100        |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| bottom_maxiter          |  Maximum number of iterations in the nodal projection                 |    Int      |   100        |
|                         |  bottom solver if using bicg, cg, bicgcg or cgbicg                    |             |              |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| mg_max_coarsening_level |  Maximum number of coarser levels to allow in the nodal projection    |    Int      |    30        |
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
| bottom_verbose          |  Verbosity of bottom solver in MAC projection                         |    Int      |   0          |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| mac_tol                 |  Relative tolerance in MAC projection                                 |    Real     |   1.e-12     |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| mac_sync_tol            |  Relative tolerance in MAC sync projection                            |    Real     |   1.e-10     |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| mac_abs_tol             |  Absolute tolerance in MAC & sync projection                          |    Real     |   1.e-16     |
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


.. _sec:InputsInitialization:

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
| init_iter            | How many pressure iterations before starting the first timestep.      |  Int        |    3         |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| init_vel_iters       | How many projection iterations to ensure the velocity satisfies the   |  Int        |    3         |
|                      | constraint. Set = 0 to skip this part of the initialization.          |             |              |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+


Checking Conservation
---------------------

IAMR has the option to sum up the total mass (and a few other quantites) at a user-defined interval.

+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
|                         | Description                                                           |   Type      | Default      |
+=========================+=======================================================================+=============+==============+
| ns.sum_interval         | How often (in level-0 time steps) to compute and print integral       |    Int      |   -1         |
|                         | quantities. If <= 0, do nothing.                                      |             |              |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+

The integral quantities include total mass, tracer, kinetic energy and magnitude of vorticity.
They are summed over the entire the domain every ns.sum_interval level-0 steps.
The print statements have the form

::

   TIME= 1.91717746 MASS= 1.792410279e+02

for example.

Note that by default the tracer not conservative. To conservatively advect the tracer,
that option must be set in the inputs (see :ref:`sec:conserv`).

