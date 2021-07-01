
Common inputs Options
=====================

**Important**: because the inputs file is handled by the C++,
any quantities you specify in scientific notation, must take the
form 1.e5 and not 1.d5—the ‘d’ specifier is not recognized.

Additionally, note that in IAMR, all quantities are in MKS units.






Simulation Time
---------------

There are three paramters that can define when a simulation ends:

-  max\_step: maximum number of level 0 time steps (Integer :math:`\geq 0`; default: -1)

-  stop\_time: final simulation time (Real :math:`\geq 0`; default: -1.0)

-  stop\_when\_steady: final simulation time (bool; default: false)

stop\_when\_steady has an optional parameter to specify steady state:

-  steady\_tol: final simulation time (Real :math:`\geq 0`; default: 1.0e-10)

To control the number of time steps, you can limit by the maximum
number of level 0 time steps (max\_step) or by the final
simulation time (stop\_time), or both. The code will stop at
whichever criterion comes first.

Note that if the code reaches stop\_time then the final time
step will be shortened so as to end exactly at stop\_time, not
past it.

As an example:

::

    max_step  = 1000
    stop_time  = 1.0

will end the calculation when either the simulation time reaches 1.0 or
the number of level 0 steps taken equals 1000, whichever comes first.

Time Step
---------

The following parameters affect the timestep choice:

-  ns.cfl: CFL number (Real :math:`> 0` and :math:`\leq 1`; default: 0.8)

-  ns.init\_shrink: factor by which to shrink the initial
   time step (Real :math:`> 0` and :math:`\leq 1`; default: 1.0)

-  ns.change\_max: factor by which the time step can grow in
   subsequent steps (Real :math:`\geq 1`; default: 1.1)

-  ns.fixed\_dt: level 0 time step regardless of cfl or other settings
   (Real :math:`> 0`; unused if not set)

-  ns.dt\_cutoff: time step below which calculation will abort
   (Real :math:`> 0`; default: 0.0)

As an example, consider:

::

    ns.cfl = 0.9 
    ns.init_shrink = 0.01 
    ns.change_max = 1.1
    ns.dt_cutoff = 1.e-20

This defines the :math:`\mathtt{cfl}` parameter in Eq. [eq:cfl] to be 0.9,
but sets (via init\_shrink) the first timestep we take
to be 1% of what it would be otherwise. This allows us to
ramp up to the hydrodynamic timestep at the start of a simulation.
The change\_max parameter restricts the timestep from increasing
by more than 10% over a coarse timestep. Note that the time step
can shrink by any factor; this only controls the extent to which it can grow.
The dt\_cutoff parameter will force the code to abort if the
timestep ever drops below :math:`10^{-20}`. This is a safety feature—if the
code hits such a small value, then something likely went wrong in the
simulation, and by aborting, you won’t burn through your entire allocation
before noticing that there is an issue.

If we know what we are doing, then we can force a particular timestep:

::

    ns.fixed_dt = 1.e-4

sets the level 0 time step to be 1.e-4 for the entire simulation,
ignoring the other timestep controls. Note that if ns.init\_shrink :math:`\neq 1` then the first time step will in fact be
ns.init\_shrink :math:`\cdot` ns.fixed\_dt.

Restart Capability
------------------

IAMR has a standard sort of checkpointing and restarting capability.
In the inputs file, the following options control the generation of
checkpoint files (which are really directories):

-  amr.check\_file: prefix for restart files (text; default: chk)

-  amr.check\_int: how often (by level 0 time steps) to write
   restart files (Integer :math:`> 0`; default: -1)

-  amr.check\_per: how often (by simulation time) to
   write restart files (Real :math:`> 0`; default: -1.0)

   Note that amr.check\_per will write a checkpoint at the first
   timestep whose ending time is past an integer multiple of this interval.
   In particular, the timestep is not modified to match this interval, so
   you won’t get a checkpoint at exactly the time you requested.

-  amr.restart: name of the file (directory) from which to restart
   (Text; not used if not set)

-  amr.checkpoint\_files\_output: should we write checkpoint files? (0 or 1; default: 1)

   If you are doing a scaling study then set amr.checkpoint\_files\_output = 0 so you can test scaling of the
   algorithm without I/O.

-  amr.check\_nfiles: how parallel is the writing of the checkpoint files?
   (Integer :math:`\geq 1`; default: 64)

   See the Software Section for more details on parallel I/O and the
   amr.check\_nfiles parameter.

-  amr.checkpoint\_on\_restart: should we write a checkpoint immediately after restarting?
   (0 or 1; default: 0)

Note:

-  You can specify both amr.check\_int or amr.check\_per,
   if you so desire; the code will print a warning in case you did this
   unintentionally. It will work as you would expect – you will get checkpoints
   at integer multiples of amr.check\_int timesteps and at integer
   multiples of amr.check\_per simulation time intervals.

-  amr.plotfile\_on\_restart and amr.checkpoint\_on\_restart
   only take effect if amr.regrid\_on\_restart is in effect.

As an example,

::

    amr.check_file = chk_run
    amr.check_int = 10

means that restart files (really directories) starting with the prefix
“chk\_run” will be generated every 10 level-0 time steps. The
directory names will be chk\_run00000, chk\_run00010, chk\_run00020, etc.

If instead you specify

::

    amr.check_file = chk_run
    amr.check_per = 0.5

then restart files (really directories) starting with the prefix
“chk\_run” will be generated every 0.1 units of
simulation time. The directory names will be chk\_run00000,
chk\_run00043, chk\_run00061, etc, where :math:`t = 0.1` after
43 level-0 steps, :math:`t = 0.2` after 61 level-0 steps, etc.

To restart from chk\_run00061, for example, then set

::

    amr.restart = chk_run00061

Controlling Plotfile Generation
-------------------------------

The main output from IAMR is in the form of plotfiles (which are
really directories). The following options in the inputs file control
the generation of plotfiles:

-  amr.plot\_file: prefix for plotfiles (text; default:
   “plt”)

-  amr.plot\_int: how often (by level-0 time steps) to write
   plot files (Integer :math:`> 0`; default: -1)

-  amr.plot\_per: how often (by simulation time) to write
   plot files (Real :math:`> 0`; default: -1.0)

   Note that amr.plot\_per will write a plotfile at the first
   timestep whose ending time is past an integer multiple of this interval.
   In particular, the timestep is not modified to match this interval, so
   you won’t get a checkpoint at exactly the time you requested.

-  amr.plot\_vars: name of state variables to include in
   plotfiles (valid options: ALL, NONE or a list; default:
   ALL)

-  amr.derive\_plot\_vars: name of derived variables to
   include in plotfiles (valid options: ALL, NONE or a
   list; default: NONE

-  amr.plot\_files\_output: should we write plot files? (0 or
   1; default: 1)

   If you are doing a scaling study then set amr.plot\_files\_output
   = 0 so you can test scaling of the algorithm without I/O.

-  amr.plotfile\_on\_restart: should we write a plotfile
   immediately after restarting? (0 or 1; default: 0)

-  amr.plot\_nfiles: how parallel is the writing of the
   plotfiles? (Integer :math:`\geq 1`; default: 64)

   See the Software Section for more details on parallel I/O and the amr.plot\_nfiles parameter.

All the options for amr.derive\_plot\_vars are kept in
``derive_lst`` in Iamr\_setup.cpp. Feel free to look at
it and see what’s there.

Some notes:

-  You can specify both amr.plot\_int or amr.plot\_per,
   if you so desire; the code will print a warning in case you did this
   unintentionally. It will work as you would expect – you will get plotfiles
   at integer multiples of amr.plot\_int timesteps and at integer
   multiples of amr.plot\_per simulation time intervals.

As an example:

::

    amr.plot_file = plt_run
    amr.plot_int = 10

means that plot files (really directories) starting with the prefix
“plt\_run” will be generated every 10 level-0 time steps. The
directory names will be plt\_run00000, plt\_run00010, plt\_run00020, etc.

If instead you specify

::

    amr.plot_file = plt_run
    amr.plot_per = 0.5

then restart files (really directories) starting with the prefix
“plt\_run” will be generated every 0.1 units of simulation time. The
directory names will be plt\_run00000, plt\_run00043, plt\_run00061, etc, where :math:`t = 0.1` after 43 level-0 steps, :math:`t =
0.2` after 61 level-0 steps, etc.

Screen Output
-------------

There are several options that set how much output is written to the
screen as IAMR runs:

-  amr.v: verbosity of Amr.cpp (0 or 1; default: 0)

-  ns.v: verbosity of NavierStokesBase.cpp (0 or 1; default: 0)

-  diffusion.v: verbosity of Diffusion.cpp (0 or 1; default: 0)

-  mg.v: verbosity of multigrid solver (for gravity) (allow
   values: 0,1,2,3,4; default: 0)

-  amr.grid\_log: name of the file to which the grids are
   written (text; not used if not set)

-  amr.run\_log: name of the file to which certain output is
   written (text; not used if not set)

-  amr.run\_log\_terse: name of the file to which certain
   (terser) output is written (text; not used if not set)

-  amr.sum\_interval: if :math:`> 0`, how often (in level-0 time
   steps) to compute and print integral quantities (Integer; default: -1)

   The integral quantities include total mass, momentum and energy in
   the domain every ns.sum\_interval level-0 steps.
   The print statements have the form

   ::

           TIME= 1.91717746 MASS= 1.792410279e+34
         

   for example. If this line is commented out then
   it will not compute and print these quanitities.

As an example:

::

    amr.grid_log = grdlog
    amr.run_log = runlog 

Every time the code regrids it prints a list of grids at all relevant
levels. Here the code will write these grids lists into the file grdlog. Additionally, every time step the code prints certain
statements to the screen (if amr.v = 1), such as:

::

    STEP = 1 TIME = 1.91717746 DT = 1.91717746 
    PLOTFILE: file = plt00001 

The run\_log option will output these statements into *runlog* as well.

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

Other parameters
----------------

There are a large number of solver-specific runtime parameters. We describe these
together with the discussion of the physics solvers in later chapters.
