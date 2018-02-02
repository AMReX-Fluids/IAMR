#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H> 
#include <AMReX_CONSTANTS.H> 
#include <AMReX_BC_TYPES.H> 
#include <AMReX_ArrayLim.H>

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define SMALL 1.0e-10
#else
#define SMALL 1.0d-10
#endif

c
c ::: This module contains helper functions for computing derived
c ::: quantities related to viscoplastic fluids. 
c
      module viscoplasticity_module
#include <NSCOMM_F.H>

      contains 

      REAL_T function strnrt_fun_2d(ux,vx,uy,vy)
         ! Magnitude of rate-of-strain tensor (2D)
         REAL_T, intent(in) :: ux,vx,uy,vy
         strnrt_fun_2d = sqrt(ux**2+half*(vx+uy)**2+vy**2)
      end function strnrt_fun_2d

      REAL_T function strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
         ! Magnitude of rate-of-strain tensor (3D)
         REAL_T, intent(in) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
         strnrt_fun_3d = sqrt(ux**2+vy**2+wz**2
     &              +half*(uy+vx)**2+half*(uz+wx)**2+half*(vz+wy)**2)
      end function strnrt_fun_3d

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
         ! Regularised Herschel-Bulkley viscosity 
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

      REAL_T function stress_fun(strnrt)
      ! Regularised Herschel-Bulkley stress
         REAL_T, intent(in)  :: strnrt
         REAL_T              :: visc
         if (varvisc .eq. 0) then
            ! If non-variable, stress equals 2 mu strnrt everywhere
            stress_fun = 2*mu*strnrt
         else
            ! Calculate viscosity based on flow index value, use this to
            ! find stress
            visc = visc_fun(strnrt)
            stress_fun = 2*visc*strnrt
         end if
      end function stress_fun

      end module viscoplasticity_module

