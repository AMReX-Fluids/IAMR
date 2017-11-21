c
c ::: This module contains helper functions for computing derived
c ::: quantities related to viscoplastic fluids. 
c
      module viscoplasticity_module

#include <NSCOMM_F.H>

      contains 
 
      REAL_T function strnrt_fun(ux,vx,uy,vy)
         ! Magnitude of rate-of-strain tensor (2D)
         REAL_T, intent(in) :: ux,vx,uy,vy
         strnrt_fun = sqrt(ux**2+half*(vx+uy)**2+vy**2)
      end function strnrt_fun

      REAL_T function regularisation(strnrt)
         ! Function for computing the regularised term in the viscosity
         ! When the strain rate is small, Taylor expand the exponential
         ! to avoid floating point errors.
         REAL_T, intent(in) :: strnrt
         if (strnrt > SMALL) then
            regularisation = (one-exp(-strnrt/eps))/strnrt
         else
            regularisation = one/eps*
     &                (one
     &                -half*strnrt/eps
     &                +sixth*(strnrt/eps)**2
     &                -sixth*fourth*(strnrt/eps)**3
     &                +sixth*fourth*fifth*(strnrt/eps)**4)
          end if
      end function regularisation

      REAL_T function visc_fun(strnrt)
         ! Function for computing the regularised Herschel-Bulkley
         ! viscosity depending on the strain rate and power-law index
         REAL_T, intent(in)  :: strnrt
         REAL_T              :: reg
         ! Calculate visc_fun based on flow index value
         reg = regularisation(strnrt)
         if (n .lt. one) then
            ! Pseudoplastic (shear-thinning)
            ! Regularise power-law term and yield stress term
            visc_fun = (mu*two**((n-one)/two)*strnrt**n+half*tau)*reg
         else if (n .eq. one) then
            ! Bingham fluid
            ! Regularise yield stress term only
            visc_fun = mu+half*tau*reg
         else 
            ! Dilatant (shear-thickening)
            ! Regularise yield stress term only
            visc_fun = mu*two**((n-one)/two)*strnrt**(n-one)+half*tau*reg
         end if
      end function visc_fun

      ! Function for computing the regularised Herschel-Bulkley
      ! stress depending on the strain rate and power-law index
      REAL_T function stress_fun(strnrt)
         REAL_T, intent(in)  :: strnrt
         REAL_T              :: reg
         if (varvisc .eq. 0) then
            ! If non-variable, stress equals 2 mu strnrt everywhere
            stress_fun = 2*mu*strnrt
         else
            ! Calculate viscosity based on flow index value
            reg = regularisation(strnrt)
            if (n .lt. one) then
               ! Pseudoplastic (shear-thinning)
               ! Regularise power-law term and yield stress_fun term
               stress_fun = (mu*two**((n+one)/two)*strnrt**n+tau)*(one-exp(-strnrt/eps))
            else if (n .eq. one) then
               ! Bingham fluid
               ! Regularise yield stress_fun term only
               stress_fun = mu*two*strnrt+tau*(one-exp(-strnrt/eps))
            else 
               ! Dilatant (shear-thickening)
               ! Regularise yield stress_fun term only
               stress_fun = mu*two**((n+one)/two)*strnrt**n+tau*(one-exp(-strnrt/eps))
            end if
         end if
      end function stress_fun

      end module viscoplasticity_module

