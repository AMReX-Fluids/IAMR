
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PROB_NS_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 2

module make_force_2d_moldule

  implicit none

  private

  public :: FORT_MAKEFORCE
  
contains

!c
!    This routine computes the forcing terms that will be added to the momentum equation
!c
      subroutine FORT_MAKEFORCE(time,force, &
                               vel, &
                               scal, &
                               DIMS(force), &
                               DIMS(vel), &
                               DIMS(scal), &
                               dx,xlo,xhi,gravity,scomp,ncomp, &
                               nscal,getForceVerbose &
     )bind(C, name="FORT_MAKEFORCE")

      implicit none

      integer    DIMDEC(force)
      integer    DIMDEC(scal)
      integer    scomp, ncomp
      REAL_T     time, dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     force  (DIMV(force),scomp:scomp+ncomp-1)
      REAL_T     gravity
      integer    DIMDEC(vel)
      integer    getForceVerbose, nscal
      REAL_T     vel    (DIMV(vel),0:SDIM-1)
      REAL_T     scal   (DIMV(scal),0:nscal-1)

      integer    i,j,n
      integer ilo, jlo
      integer ihi, jhi
      integer nXvel, nYvel, nRho, nTrac, nTrac2
      integer nRhoScal, nTracScal, nTrac2Scal

      ilo = force_l1
      jlo = force_l2
      ihi = force_h1
      jhi = force_h2

!c     Assumes components are in the following order
      nXvel = 0
      nYvel = 1
      nRho  = 2
      nTrac = 3
      nTrac2= 4

      nRhoScal   = nRho-SDIM
      nTracScal  = nTrac-SDIM
      nTrac2Scal = nTrac2-SDIM

      if (scomp.eq.0) then
!c
!c     Do velocity forcing
!c
         do j = jlo, jhi
            do i = ilo, ihi
               force(i,j,nXvel) = zero
               force(i,j,nYvel) = zero
            enddo
         enddo
!c     End of velocity forcing
      endif

      if ((scomp+ncomp).gt.BL_SPACEDIM) then
!c
!c     Scalar forcing
!c
         do n = max(scomp,nRho), scomp+ncomp-1
            if (n.eq.nRho) then
!c
!c     Density
!c
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,n) = zero
                  enddo
               enddo
            else if (n.eq.nTrac) then
!c
!c     Tracer
!c
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,n) = zero
                  enddo
               enddo
            else
!c
!c     Other scalar
!c
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,n) = zero
                  enddo
               enddo
            endif
         enddo
      endif


    end subroutine FORT_MAKEFORCE

  end module make_force_2d_moldule
