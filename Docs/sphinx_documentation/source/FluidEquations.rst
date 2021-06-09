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
   | :math:`g`             | Gravitational acceleration                       |
   +-----------------------+--------------------------------------------------+

Fluid Equations
===============

Conservation of fluid mass:

.. math:: \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho U)  = 0

Conservation of fluid momentum:

.. math:: \frac{ \partial (\rho U)}{\partial t} 
   + \nabla \cdot (\rho U g) + \nabla p = \nabla \cdot \tau + \rho g

Incompressibility constraint:

.. math:: \nabla \cdot U = 0

Tracer(s) advection:

.. math:: \frac{\partial \rho s}{\partial t} + \nabla \cdot (\rho U s)  = 0
