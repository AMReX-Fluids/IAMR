

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
   
   
   
Fluid Equations
===============

Conservation of fluid mass:

.. math:: \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho U)  = 0

Conservation of fluid momentum:

.. math:: \frac{ \partial (\rho U)}{\partial t} 
   + \nabla \cdot (\rho U) + \nabla p = \nabla \cdot \tau + {\bf H}_U)

Incompressibility constraint:

.. math:: \nabla \cdot U = 0

Tracer(s):

.. math:: \frac{\partial \rho s}{\partial t} + \nabla \cdot (\rho U s)  = \nabla \cdot \beta \nabla \frac{s}{\rho} + H_s

for conservatively advected scalars and 

.. math:: \frac{\partial s}{\partial t} + U \cdot \nabla s  = \nabla \cdot \beta \nabla s + H_s

for passively advected scalars. In general, one could advect an arbitrary number of scalars,
either passively or conservatively.

IAMR has the ability to incorporate general, user-defined external forcing and source terms. The default behaviour is that 
:math:`H_c=0`, and :math:`{\bf H}_U` represents gravitational forces, with :math:`{\bf H}_U= (0 , 0 , -\rho g )` in 3d and
:math:`{\bf H}_U= (0 , -\rho g )` in 2d, where :math:`g` is the magnitude of the gravitational acceleration. However, since
by default, :math:`g=0`, :math:`{\bf H}_U = 0` unless ``ns.gravity`` is set (for more info see :ref:`sec:Inputs` and ).

