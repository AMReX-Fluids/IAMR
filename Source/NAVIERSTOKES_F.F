#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <NAVIERSTOKES_F.H>
#include <AMReX_ArrayLim.H>

c=========================================================

      subroutine FORT_SET_NS_PARAMS(vel_visc_coef,yield_stress,reg_param)
      implicit none
      REAL_T vel_visc_coef,yield_stress,reg_param

#include <NSCOMM_F.H>

      mu = vel_visc_coef
      tau = yield_stress
      eps = reg_param

      if (mu .lt. 0) then
        write(6,*)'FORT_SET_PARAMS : illegal value of mu = ',mu
        call bl_abort(" ")
      end if

      if (tau .lt. 0) then
        write(6,*)'FORT_SET_PARAMS : illegal value of tau = ',tau
        call bl_abort(" ")
      end if

      if (eps .lt. 0) then
        write(6,*)'FORT_SET_PARAMS : illegal value of eps = ',eps
        call bl_abort(" ")
      end if

      end
