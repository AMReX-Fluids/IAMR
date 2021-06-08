.. _Chap:InputsMultigrid:

Multigrid Inputs
================

The following inputs can be set directly in the AMReX solver classes but we 
set them via the incflo routines because we may want different inputs for the 
different solvers called by incflo
NOTE: the nodal solver settings are read in directly by AMReX, 
the MAC and diffusion settings by incflo.

These control the nodal projection and must be preceded by "nodal_proj": 

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

These control the MAC projection and must be preceded by "mac_proj":

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

These control the diffusion solver and must be preceded by "diffusion":

+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
|                         | Description                                                           |   Type      | Default      |
+=========================+=======================================================================+=============+==============+
| verbose                 |  Verbosity of linear solver for diffusion solve                       |    Int      |   0          |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| bottom_verbose          |  Verbosity of BiCGStab solver in diffusion solve                      |    Int      |   0          |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| rtol                    |  Relative tolerance in diffusion solve                                |    Real     |   1.e-11     | 
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| atol                    |  Absolute tolerance in diffusion solve                                |    Real     |   1.e-14     | 
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| maxiter                 |  Maximum number of iterations in diffusion solve                      |    Int      |   100        |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| bottom_maxiter          |  Maximum number of iterations in diffusion solve                      |    Int      |   100        |
|                         |  bottom solver if using bicg, cg, bicgcg or cgbicg                    |             |              |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| mg_max_coarsening_level |  Maximum number of coarser levels to allow in diffusion solve         |    Int      |   100        |
|                         |  If set to 0, the bottom solver will be called at the current level   |             |              |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| bottom_solver           |  Which bottom solver to use in the diffusion solve                    |  String     |   bicgcg     |
|                         |  Options are bicgcg, bicgstab, cg, cgbicg, smoother or hypre          |             |              | 
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
