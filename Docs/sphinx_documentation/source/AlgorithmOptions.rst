
.. _Chap:AlgorithmOptions:

Algortihm Options
=================

.. _sec:conserv:

Conservative vs. Non-conservative
---------------------------------

The following must be preceeded by "ns."

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

IAMR has the option to use a Method of Lines (MOL) or Godunov scheme to compute the advective terms.

+-------------------------+-------------------------------------------------------------------------+-------------+--------------+
|                         | Description                                                             |   Type      | Default      |
+=========================+=========================================================================+=============+==============+
| ns.use_godunov          | If true, use Godunov, else use MOL.                                     |    bool     |   true       |
+-------------------------+-------------------------------------------------------------------------+-------------+--------------+


For problems without embedded boundaries, there are additional options when using the Godunov method. The following must
be preceeded by "godunov."

+-------------------------+-------------------------------------------------------------------------+-------------+--------------+
|                         | Description                                                             |   Type      | Default      |
+=========================+=========================================================================+=============+==============+
| use_ppm                 | Use the Piecewise Parabolic Method to construct edge states             |    bool     |   false      |
+-------------------------+-------------------------------------------------------------------------+-------------+--------------+
| use_forces_in_trans     | Use external forcing terms in constructing transverse derivatives       |    bool     |   false      |
+-------------------------+-------------------------------------------------------------------------+-------------+--------------+


Diffusion
---------

The following must be preceeded by "ns."

+-------------------------+-----------------------------------------------------------------------+-------------+--------------+
|                         | Description                                                           |   Type      | Default      |
+=========================+=======================================================================+=============+==============+
| be_cn_theta             | Diffusion solve fully implicit (1.0) or semi-implicit (<1 && >0.5)    |   Real      |   0.5        |
+-------------------------+-----------------------------------------------------------------------+-------------+--------------+

Note the default value of ``ns.be_cn_theta = 0.5`` corresponds to the Crank-Nicolson method.
