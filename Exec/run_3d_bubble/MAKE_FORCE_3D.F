
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "PROB_NS_F.H"
#include "ArrayLim.H"

#define SDIM 3
#ifdef BOUSSINESQ
c
c     This routine computes the forcing terms that will be added to the momentum equation
c
      subroutine FORT_MAKEFORCE(force,
     &                          scal,
     &                          DIMS(force),
     &                          DIMS(scal),
     &                          dx,xlo,xhi,scomp,ncomp,
     &                          nscal)

      implicit none

      integer    DIMDEC(force)
      integer    DIMDEC(scal)
      integer    scomp, ncomp, nscal
      REAL_T     dx(SDIM), xlo(SDIM), xhi(SDIM)
      REAL_T     force  (DIMV(force),scomp:scomp+ncomp-1)
      REAL_T     scal   (DIMV(scal),0:nscal-1)

      integer    i,j,k

      print *,'SCOMP NCOMP ',scomp,ncomp

c     Here the scalar array has rho in the first component and 
c                               phi in the second component 

#include "probdata.H"

      do k = force_l3, force_h3
      do j = force_l2, force_h2
      do i = force_l1, force_h1
          force(i,j,k,SDIM) = -scal(i,j,k,1) * 1.d0
      end do
      end do
      end do

      end
#endif
