
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PROB_NS_F.H>
#include <AMReX_ArrayLim.H>

module make_force_3d_module

  implicit none

  private

  public :: FORT_MAKEFORCE

contains
!c
!    This routine computes the forcing terms that will be added to the momentum equation
!c
    subroutine FORT_MAKEFORCE(force, &
                               scal,&
                               DIMS(force),&
                               DIMS(scal),&
                               dx,xlo,xhi,scomp,ncomp) &
                               bind(c,name="FORT_MAKEFORCE")

      implicit none

      integer    DIMDEC(force)
      integer    DIMDEC(scal)
      integer    scomp, ncomp
      REAL_T     dx(3), xlo(3), xhi(3)
      REAL_T     force  (DIMV(force),scomp:scomp+ncomp-1)
      REAL_T     scal   (DIMV(scal))

      integer    i,j,k
      REAL_T     beta

      beta = 5.d0

!    Here the  force array has Xvel in the first component, Yvel in the second, Zvel in the third

      do k = force_l3, force_h3
      do j = force_l2, force_h2
      do i = force_l1, force_h1
          force(i,j,k,scomp  ) = 0.d0
          force(i,j,k,scomp+1) = 0.d0
          force(i,j,k,scomp+2) = 0.d0
          force(i,j,k,scomp+2) = scal(i,j,k) * beta
      end do
      end do
      end do
      call flush(6)

    end subroutine FORT_MAKEFORCE
  end module make_force_3d_module
