
.. _Chap:AlgorithmOptions:

Algorithm Options
=================

.. _sec:conserv:

Conservative vs. Non-conservative
---------------------------------

The following must be preceded by "ns."

+-------------------------+------------------------------------------------------------------------------+-------------+--------------+
|                         | Description                                                                  |   Type      | Default      |
+=========================+==============================================================================+=============+==============+
| do_mom_diff             | If 0, solve velocity equation in convective form, else use conservation form |    Int      |   0          |
+-------------------------+------------------------------------------------------------------------------+-------------+--------------+
| do_cons_trac            | If 0, solve for a passively advected tracer, else advect conservatively      |    Int      |   0          |
+-------------------------+------------------------------------------------------------------------------+-------------+--------------+
| do_cons_trac2           | If 0, solve for a passively advected 2nd tracer, else advect conservatively  |    Int      |   0          |
+-------------------------+------------------------------------------------------------------------------+-------------+--------------+

Note that Temperature is only non-conservative. For more details, see :ref:`sec:FluidEquations`.


Advection
---------

IAMR computes the advective terms with an unsplit Godunov scheme (piecewise linear or piecewise
parabolic reconstruction) or with the Bell-Dawson-Shubin (BDS) scheme. The following must be
preceded by "ns."

+-------------------------+-------------------------------------------------------------------------+-------------+--------------+
|                         | Description                                                             |   Type      | Default      |
+=========================+=========================================================================+=============+==============+
| advection_scheme        | Godunov_PLM, Godunov_PPM or BDS.  Godunov_PPM and BDS are not           |   String    | Godunov_PLM  |
|                         | available with embedded boundaries.                                     |             |              |
+-------------------------+-------------------------------------------------------------------------+-------------+--------------+

Note that the old ``ns.use_godunov`` key and the MOL scheme have been removed; setting either
aborts the run.


For problems without embedded boundaries, there is an additional option for the Godunov method. The following must
be preceded by "godunov."

+-------------------------+-------------------------------------------------------------------------+-------------+--------------+
|                         | Description                                                             |   Type      | Default      |
+=========================+=========================================================================+=============+==============+
| use_forces_in_trans     | Use external forcing terms in constructing transverse derivatives       |    bool     |   false      |
+-------------------------+-------------------------------------------------------------------------+-------------+--------------+


Diffusion
---------

The following must be preceded by "ns."

+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
|                         | Description                                                           |   Type      | Default      |
+=========================+=======================================================================+=============+==============+
| be_cn_theta             | Diffusion solve fully implicit (1.0) or semi-implicit (<1 && >0.5)    |   Real      |   0.5        |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+

Note the default value of ``ns.be_cn_theta = 0.5`` corresponds to the Crank-Nicolson method.


.. _sec:LES:

Large Eddy Simulation
---------------------

IAMR can add a subgrid-scale eddy viscosity to the viscous terms.  The following must be
preceded by "ns."

+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
|                         | Description                                                           |   Type      | Default      |
+=========================+=======================================================================+=============+==============+
| do_LES                  | Add a subgrid-scale eddy viscosity to the molecular viscosity         |    Int      |   0          |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| LES_model               | Which model to use: Smagorinsky or Sigma.  Any other value aborts.    |  String     | Smagorinsky  |
|                         | Sigma is 3D only.                                                     |             |              |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| smago_Cs_cst            | Model constant, used only when LES_model = Smagorinsky                |   Real      |   0.18       |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| sigma_Cs_cst            | Model constant, used only when LES_model = Sigma                      |   Real      |   1.5        |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| getLESVerbose           | Print the model and constant in use from the LES routine              |    Int      |   0          |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+

Note that each model reads its own constant, so changing ``ns.smago_Cs_cst`` has no effect
when ``ns.LES_model = Sigma``, and vice versa.  Both models compute a kinematic eddy
viscosity :math:`\nu_t`; since IAMR's viscosity is dynamic, :math:`\mu_t = \rho \nu_t`
(with :math:`\rho` averaged to faces) is what is added to it.  The Sigma model is described in
Nicoud et al., *Using singular values to build a subgrid-scale model for large eddy
simulations*, Phys. Fluids 23, 085106 (2011).


.. _sec:EBOptions:

Embedded Boundaries
-------------------

These apply to builds with ``USE_EB=TRUE``; see :ref:`sec:EB-basics` for how the geometry
itself is constructed.  The following must be preceded by "ns."

+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
|                         | Description                                                           |   Type      | Default      |
+=========================+=======================================================================+=============+==============+
| redistribution_type     | How the advective update of a small cut cell is redistributed to its  |  String     | StateRedist  |
|                         | neighbours: StateRedist, FluxRedist or NoRedist.  Any other value     |             |              |
|                         | aborts.                                                               |             |              |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
| refine_cutcells         | Tag every cut cell for refinement, so that the embedded boundary      |    Int      |   1          |
|                         | never crosses a coarse/fine boundary.  Setting 0 allows a partially   |             |              |
|                         | refined EB, which is still under development and issues a warning.    |             |              |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+

Note that ``ns.advection_scheme = Godunov_PPM`` and ``BDS`` are not available with embedded
boundaries.
