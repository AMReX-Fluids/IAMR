.. sec:InputsTimeStepping:

Time Stepping
=============

The following inputs must be preceded by "incflo."   Note that if both are specified, both criteria
are used and the simulation still stop when the first criterion is hit.  In the case of unsteady flow,
the simulation will stop when either the number of steps reaches max_step or time reaches stop_time.
In the case of unsteady flow, the simulation will stop when either the tolerance (difference between
subsequent steps) is reached or the number of iterations reaches the maximum number specified.

+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|                      | Description                                                           |   Type      | Default      |
+======================+=======================================================================+=============+==============+
| max_step             | Maximum number of time steps to take                                  |    Int      |  -1          |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| stop_time            | Maximum time to reach                                                 |    Real     | -1.0         |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| fixed_dt             | Should we use a fixed timestep?                                       |    Int      |   0          |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| cfl                  | CFL constraint (dt < cfl * dx / u) if fixed_dt not 1                  |    Real     |   0.5        |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+

Setting the Time Step 
---------------------

  * If you want to fix the dt, simply set :cpp:`incflo.fixed_dt = XXX` and the fluid time
    step will always be that number. 

  * If you want to let the code determine the appropriate time step using the advective CFL
    condition, then set :cpp:`incflo.cfl = 0.7` for example, and the fluid time step will
    be computed to be dt = 0.5 * dx / max(vel).

  * Note that the cfl defaults to 0.5 so it does not have to be set in the inputs file. If neither
    :cpp:`incflo.cfl` nor :cpp:`fixed_dt` is set, then default value of cfl will be used.
    If :cpp:`incflo.fixed_dt` is set, then it will override the cfl option whether 
    :cpp:`incflo.cfl` is set or not.
