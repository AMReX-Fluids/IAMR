
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <GODUNOV_F.H>
#include <AMReX_ArrayLim.H>

module godunov_module

  implicit none

  private

  public :: set_params

contains

!c=========================================================

      subroutine set_params(slope_order_in,use_unlim_in) bind(C,name="set_params")
      implicit none
      integer slope_order_in,use_unlim_in

#include <GODCOMM_F.H>

      slope_order = slope_order_in

      use_unlimited_slopes = (use_unlim_in .eq. 1)

      if(slope_order.ne.1.and.&
#if (BL_SPACEDIM==2)
        slope_order.ne.2.and.&
#endif
        slope_order.ne.4)then
        write(6,*)'set_params : illegal value of slope_order = ',&
                 slope_order

        call bl_abort(" ")
      end if

    end subroutine set_params
  end module godunov_module
