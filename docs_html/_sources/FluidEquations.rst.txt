

Fluid Variables
===============

   +-----------------------+--------------------------------------------------+
   | Variable              | Definition                                       |
   +=======================+==================================================+
   | :math:`\rho`          | Fluid density                                    |
   +-----------------------+--------------------------------------------------+
   | :math:`U`             | Fluid velocity                                   |
   +-----------------------+--------------------------------------------------+
   | :math:`\tau`          | Viscous stress tensor                            |
   +-----------------------+--------------------------------------------------+
   | :math:`{\bf H}_U`     | :math:`= (H_x , H_y , H_z )`, External Forces    |
   +-----------------------+--------------------------------------------------+
   | :math:`H_s`           | External sources                                 |
   +-----------------------+--------------------------------------------------+
   

.. _sec:FluidEquations:
   
Fluid Equations
===============

Conservation of fluid mass:

.. math:: \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho U)  = 0

Conservation of fluid momentum:

.. math:: \frac{ \partial (\rho U)}{\partial t} 
   + \nabla \cdot (\rho U U) + \nabla p = \nabla \cdot \tau + {\bf H}_U

Incompressibility constraint:

.. math:: \nabla \cdot U = 0

Tracer(s):

.. math:: \frac{\partial \rho s}{\partial t} + \nabla \cdot (\rho s U)  = \nabla \cdot \beta \nabla s + \rho H_s

for conservatively advected scalars and 

.. math:: \frac{\partial s}{\partial t} + U \cdot \nabla s  = \nabla \cdot \beta \nabla s + H_s

for passively advected scalars. In general, one could advect an arbitrary number of scalars.

IAMR has the ability to incorporate general, user-defined external forcing and source terms. The default behaviour is that 
:math:`H_c=0`, and :math:`{\bf H}_U` represents gravitational forces, with :math:`{\bf H}_U= (0 , 0 , -\rho g )` in 3d and
:math:`{\bf H}_U= (0 , -\rho g )` in 2d, where :math:`g` is the magnitude of the gravitational acceleration. However, since
by default, :math:`g=0`, :math:`{\bf H}_U = 0` unless ``ns.gravity`` is set (for more info see :ref:`sec:PhysicsParams`).

By default, IAMR solves the momentum equation in convective form. The inputs parameter ``ns.do_mom_diff`` is used to
switch to conservation form. Tracers are passively advected by default. The inputs parameter ``ns.do_cons_trac = 1``
switches the first tracer to conservative. A second tracer can be included with ``ns.do_trac2 = 1``, and it can be
conservatively advected with ``ns.do_cons_trac2``.

IAMR also has the option to solve for temperature, along with a modified divergence constraint on the velocity field:

.. math:: \rho c_p \left( \frac{\partial T}{\partial t} + U \cdot \nabla T \right)  = \nabla \cdot \lambda \nabla T + H_T

	  \nabla \cdot U = \frac{1}{\rho c_p T} \nabla \cdot \lambda \nabla T 

To enable the temperature solve, use ``ns.do_temp = 1`` and set ``ns.temp_cond_coef`` to represent :math:`\lambda / c_p`,
which is taken to be constant. More sophiticated treatments are possible; if interested, please open an issue on github:
https://github.com/AMReX-Codes/IAMR/issues
