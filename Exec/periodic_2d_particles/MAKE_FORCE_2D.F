
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "PROB_NS_F.H"
#include "ArrayLim.H"

c
c     This routine computes the forcing terms that will be added to the momentum equation
c
c     subroutine FORT_MAKEFORCE(force,
      subroutine makeforce(force,
     &                          scal,
     &                          DIMS(force),
     &                          DIMS(scal),
     &                          dx,xlo,xhi,scomp,ncomp)

      implicit none

      integer    DIMDEC(force)
      integer    DIMDEC(scal)
      integer    scomp, ncomp
      REAL_T     dx(2), xlo(2), xhi(2)
      REAL_T     force  (DIMV(force),scomp:scomp+ncomp-1)
      REAL_T     scal   (DIMV(scal))

      integer    i,j

c     Here the  force array has Xvel in the first component, Yvel in the second

      do j = force_l2, force_h2
      do i = force_l1, force_h1
          force(i,j,scomp  ) = 0.d0
          force(i,j,scomp+1) = 0.d0
      end do
      end do

      end
