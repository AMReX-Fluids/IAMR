c
c $Id: GODUNOV_F.F,v 1.3 1997-12-12 22:44:01 lijewski Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "GODUNOV_F.H"
#include "ArrayLim.H"

c=========================================================

      subroutine FORT_SET_PARAMS(slope_order_in)
      implicit none
      integer slope_order_in

#include "GODUNOV.FH"

      slope_order = slope_order_in

      if(slope_order.ne.1.and.
#if (BL_SPACEDIM==2)
     &   slope_order.ne.2.and.
#endif
     &   slope_order.ne.4)then
        write(6,*)'FORT_SET_PARAMS : illegal value of slope_order = ',
     &            slope_order
#ifndef	BL_NO_FORT_FLUSH
        call flush(6)
#endif
        stop
      endif
      return
      end
