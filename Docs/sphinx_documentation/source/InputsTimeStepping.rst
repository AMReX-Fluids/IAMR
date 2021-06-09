.. sec:InputsTimeStepping:

Time Stepping
=============

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

The inputs below must be preceded by "ns".  

+----------------------+-----------------------------------------------------------------------+-------------+--------------+
|                      | Description                                                           |   Type      | Default      |
+======================+=======================================================================+=============+==============+
| fixed_dt             | Value of fixed dt if > 0                                              |    Real     |   -1.        |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+
| cfl                  | CFL constraint (dt < cfl * dx / u) if fixed_dt not > 0                |    Real     |   0.5        |
+----------------------+-----------------------------------------------------------------------+-------------+--------------+

Setting the Time Step 
---------------------

  * If you want to fix the dt, simply set :cpp:`ns.fixed_dt = XXX` and the fluid time
    step will always be that number. 

  * If you want to let the code determine the appropriate time step using the advective CFL
    condition, then set :cpp:`ns.cfl = 0.7` for example, and the fluid time step will
    be computed to be dt = 0.7 * dx / max(vel).

  * Note that the cfl defaults to 0.5 so it does not have to be set in the inputs file. If neither
    :cpp:`ns.cfl` nor :cpp:`fixed_dt` is set, then default value of cfl will be used.
    If :cpp:`ns.fixed_dt` is set, then it will override the cfl option whether 
    :cpp:`ns.cfl` is set or not.
