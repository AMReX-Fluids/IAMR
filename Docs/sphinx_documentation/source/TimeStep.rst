
To learn how the convective terms are constructed, see `AMReX-Hydro <https://amrex-codes.github.io/amrex/hydro_html>`_

Time Step -- MOL
~~~~~~~~~~~~~~~~

In the predictor

-  Define :math:`U^{MAC,n}`, the face-centered (staggered) MAC velocity which is used for advection, using :math:`U^n`

-  Define an approximation to the new-time state, :math:`(\rho U)^{\ast}` by setting 

.. math:: (\rho U)^{\ast} &= (\rho U)^n -  
           \Delta t \left( \nabla \cdot (\rho U^{MAC} U) + \nabla {p}^{n-1/2} \right) \\ &+ 
           \Delta t \left( \nabla \cdot \tau^n + \sum_p \beta_p (V_p - {U}^{\ast}) + \rho g \right)

-  Project :math:`U^{\ast}` by solving

.. math:: \nabla \cdot \frac{1}{\rho} \nabla \phi = \nabla \cdot \left( \frac{1}{\Delta t} 
          U^{\ast}+ \frac{1}{\rho} \nabla {p}^{n-1/2} \right)

then defining

.. math:: U^{\ast \ast} = U^{\ast} - \frac{\Delta t}{\rho} \nabla \phi

and 

.. math:: {p}^{n+1/2, \ast} = \phi


In the corrector

-  Define :math:`U^{MAC,\ast \ast}` at the "new" time using :math:`U^{\ast \ast}`

-  Define a new approximation to the new-time state, :math:`(\rho U)^{\ast \ast \ast}` by setting  

.. math:: (\rho U)^{\ast \ast \ast} &= (\rho U)^n - \frac{\Delta t}{2} \left( \nabla \cdot (\rho U^{MAC} U)^n + \nabla \cdot (\rho U^{MAC} U)^{\ast \ast}\right) + \\ &+ \frac{\Delta t}{2} \left( \nabla \cdot \tau^n + \nabla \cdot \tau^{\ast \ast \ast} \right) + \Delta t \left( - \nabla {p}^{n+1/2,\ast} + \sum_p \beta_p (V_p - {U}^{\ast \ast \ast}) + \rho g \right)

-  Project :math:`U^{\ast \ast \ast}` by solving

.. math:: \nabla \cdot \frac{1}{\rho} \nabla \phi = \nabla \cdot \left( \frac{1}{\Delta t} U^{\ast \ast \ast} + \frac{1}{\rho} \nabla {p}^{n+1/2,\ast} \right)

then defining

.. math:: U^{n+1} = U^{\ast \ast \ast} - \frac{\Delta t}{\rho} \nabla \phi

and 

.. math:: {p}^{n+1/2} = \phi

Time Step -- Godunov
~~~~~~~~~~~~~~~~~~~~

When we use the time-centered Godunov advection, we no longer need the predictor and corrector steps.

-  Define the time-centered face-centered (staggered) MAC velocity which is used for advection: :math:`U^{MAC,n+1/2}`

-  Define the new-time density, :math:`\rho^{n+1} = \rho^n - \Delta t (\rho^{n+1/2,pred} U^{MAC,n+1/2})` by setting 

-  Define an approximation to the new-time state, :math:`(\rho U)^{\ast}` by setting 

   .. math:: (\rho^{n+1} U^{\ast}) &= (\rho^n U^n) -  
             \Delta t \nabla \cdot (\rho U^{MAC} U) + \Delta t \nabla {p}^{n-1/2}  \\ &+ 
             \frac{\Delta t}{2}  (\nabla \cdot \tau^n + \nabla \cdot \tau^\ast) +
             \Delta t \rho g 

   (for implicit diffusion, which is the current default)

-  Project :math:`U^{\ast}` by solving

.. math:: \nabla \cdot \frac{1}{\rho} \nabla \phi = \nabla \cdot \left( \frac{1}{\Delta t} 
          U^{\ast}+ \frac{1}{\rho} \nabla {p}^{n-1/2} \right)

then defining

.. math:: U^{n+1} = U^{\ast} - \frac{\Delta t}{\rho} \nabla \phi

and 

.. math:: {p}^{n+1/2} = \phi


The algorithm is further described in the following paper (and references therein):

-  *A Conservative Adaptive Projection Method for the Variable Density Incompressible Navier-Stokes Equations*,
   A. S. Almgren, J. B. Bell, P. Colella, L. H. Howell, and M. L. Welcome,
   J. Comp. Phys., 142, pp. 1-46, 1998.
   http://www.sciencedirect.com/science/article/pii/S0021999198958909 :cite:`IAMR`


