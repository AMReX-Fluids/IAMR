.. _Chap:InputsPlotfiles:

Plotfiles and Other Output
==========================

The following inputs must be preceded by "amr" and control frequency and naming of plotfile generation as well
as whether the EB geometry should be written out.

+---------------------+-----------------------------------------------------------------------+-------------+-----------+
|                     | Description                                                           |   Type      | Default   |
+=====================+=======================================================================+=============+===========+
| plot_int            | Frequency of plotfile output;                                         |    Int      | -1        |
|                     | if -1 then no plotfiles will be written at this frequency             |             |           |
| plot_per_exact      | Time period of plotfile output (exact); will modify dt                |    Real     | -1        |
|                     | if -1 then no plotfiles will be written at this frequency             |             |           |
| plot_per_approx     | Time period of plotfile output (approximate); does not modify dt      |    Real     | -1        |
|                     | if -1 then no plotfiles will be written at this frequency             |             |           |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plotfile_on_restart | Should we write a plotfile when we restart (only used if plot_int>0)  |   Bool      | False     |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plot_file           | Prefix to use for plotfile output                                     |  String     | plt       |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| write_eb_surface    | Should we write out the EB geometry in vtp format                     |   Bool      | False     |
|                     | If true, it will only be written once,after initialization or restart |             |           |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+

The following inputs must be preceded by "amr" and control what variables will be written in plotfiles.

+---------------------+-----------------------------------------------------------------------+-------------+-----------+
|                     | Description                                                           |   Type      | Default   |
+=====================+=======================================================================+=============+===========+
| plt_ccse_regtest    | Save all variables to plot file (overrides all other IO flags)        |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_velx            | Save x-velocity to plot file                                          |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_vely            | Save y-velocity to plot file                                          |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_velz            | Save z-velocity to plot file                                          |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_p               | Save pressure to plot file                                            |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_rho             | Save density to plot file                                             |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_tracer          | Save tracer to plot file                                              |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_eta             | Save viscosity to plot file                                           |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_gpx             | Save dp/dx to plot file                                               |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_gpy             | Save dp/dy to plot file                                               |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_gpz             | Save dp/dz to plot file                                               |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_vort            | Save vorticity to plot file                                           |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_forcing         | Save forcing term for velocity to plot file                           |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_strainrate      | Save strain rate to plot file                                         |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_divu            | Save velocity divergence to plot file                                 |    Int      | 0         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
| plt_vfrac           | Save EB volume fraction to plot file                                  |    Int      | 1         |
+---------------------+-----------------------------------------------------------------------+-------------+-----------+
