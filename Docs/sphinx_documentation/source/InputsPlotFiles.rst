.. _Chap:InputsPlotfiles:

Plotfiles and Other Output
==========================

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
