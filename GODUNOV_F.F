c $Id: GODUNOV_F.F,v 1.2 1997-07-22 21:57:02 car Exp $
#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "GODUNOV_F.H"
#include "DIMS.H"

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
