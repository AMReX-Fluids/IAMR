
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PROB_NS_F.H>
#include <AMReX_ArrayLim.H>

module make_force_2d_moldule

  implicit none

  private

  public :: makeforce
  
contains

!c
!    This routine computes the forcing terms that will be added to the momentum equation
!c
!    subroutine FORT_MAKEFORCE(force,
      subroutine makeforce(force,&
                               scal,&
                               DIMS(force),&
                               DIMS(scal),&
                               dx,xlo,xhi,scomp,ncomp) bind(c, name="makeforce")

      implicit none

      integer    DIMDEC(force)
      integer    DIMDEC(scal)
      integer    scomp, ncomp
      REAL_T     dx(2), xlo(2), xhi(2)
      REAL_T     force  (DIMV(force),scomp:scomp+ncomp-1)
      REAL_T     scal   (DIMV(scal))

      integer    i,j

!    Here the  force array has Xvel in the first component, Yvel in the second

      do j = force_l2, force_h2
      do i = force_l1, force_h1
          force(i,j,scomp  ) = 0.d0
          force(i,j,scomp+1) = 0.d0
      end do
      end do

    end subroutine makeforce

  end module make_force_2d_moldule
