
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>
#include <CONSTANTS.H>
#include <BC_TYPES.H>
#include <GODUNOV_F.H>
#include <ArrayLim.H>

c=========================================================

      subroutine FORT_SET_PARAMS(slope_order_in,use_unlim_in)
      implicit none
      integer slope_order_in,use_unlim_in

#include <GODCOMM_F.H>

      slope_order = slope_order_in

      use_unlimited_slopes = (use_unlim_in .eq. 1)

      if(slope_order.ne.1.and.
#if (BL_SPACEDIM==2)
     &   slope_order.ne.2.and.
#endif
     &   slope_order.ne.4)then
        write(6,*)'FORT_SET_PARAMS : illegal value of slope_order = ',
     &            slope_order

        call bl_abort(" ")
      end if

      end
