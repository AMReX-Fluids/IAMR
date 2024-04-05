.. include:: CustomCommands.rst

To learn how the convective terms are constructed, see `AMReX-Hydro <https://amrex-fluids.github.io/amrex-hydro/docs_html>`_

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


Time Step -- BDS
~~~~~~~~~~~~~~~~

When we use the Bell-Dawson-Shubin (BDS) algorithm, we advance the solution in time using a
three step procedure described below. In the notation,
:math:`s` is a scalar field of the form :math:`s=s(x,y,z,t)`
and :math:`{\bf u}=(u,v,w)` represents a known velocity field. :math:`s^n_{ijk}` represents
the average value of :math:`s` over the cell with index :math:`(ijk)` at time :math:`t^n`.
At each face the normal velocity (e.g., :math:`u_{i+1/2,j,k}`) is assumed constant
over the time step.

- **Step 1**: Construct a limited piecewise trilinear (bilinear in 2D) representation of the solution in
  each grid cell of the form,

.. math::
    \begin{eqnarray}
    s_{ijk}(x,y,z) &=& s_{ijk} + s_{x,ijk}\cdot(x-x_i) + s_{y,ijk}\cdot(y-y_j) + s_{z,ijk}\cdot(z-z_k) \nonumber \\
    && + s_{xy,ijk}\cdot(x-x_i)(y-y_j) + s_{xz,ijk}\cdot(x-x_i)(z-z_k) \nonumber \\
    && + s_{yz,ijk}\cdot(y-y_j)(z-z_k) + s_{xyz,ijk}\cdot(x-x_i)(y-y_j)(z-z_k).
    \end{eqnarray}

- **Step 2**: Construct edge states :math:`s_{i+1/2,j,k}`, etc. by integrating limited
  piecewise trilinear (bilinear in 2D) profiles over the space-time region determined by the characteristic
  domain of dependence on the face.

- **Step 3**: Advance the solution in time using the conservative update equation,

.. math::
    \begin{eqnarray}
    s_{ijk}^{n+1} = s_{ijk}^n &&
    - \frac{\dt}{\Delta x}(u_{i+\half,j,k}s_{i+\half,j,k} - u_{i-\half,j,k}s_{i-\half,j,k}) \nonumber \\
    && - \frac{\dt}{\Delta y}(v_{i,j+\half,k}s_{i,j+\half,k} - v_{i,j-\half,k}s_{i,j-\half,k}) \nonumber \\
    && - \frac{\dt}{\Delta z}(w_{i,j,k+\half}s_{i,j,k+\half} - w_{i,j,k-\half}s_{i,j,k-\half}).
    \end{eqnarray}


Addition details are located in the :ref:`BDS section of the AMReX-Hydro docs <hydro:BDS:BDS Algorithm>` and in the following paper:

- *A Three-Dimensional, Unsplut Godunov Method For Scalar Conservation Laws*,
  A. Nonaka, S. May, A. S. Almgren, and J. B. Bell,
  SIAM Journal of Scientific Computation, Vol. 33, No.4, pp. 2039-2062.
  https://ccse.lbl.gov/Publications/nonaka/BDS_3d.pdf :cite:`BDS_3d`
