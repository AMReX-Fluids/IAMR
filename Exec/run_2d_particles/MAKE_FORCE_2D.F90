
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

  use amrex_fort_module, only : dim=>amrex_spacedim

  implicit none

  private

  public :: FORT_MAKEFORCE
  
contains

!
! ::: -----------------------------------------------------------
!     This routine add the forcing terms to the momentum equation

   subroutine FORT_MAKEFORCE( time, &
                              force, f_lo, f_hi,&
                              vel, v_lo, v_hi,&
                              scal, s_lo, s_hi,&
                              dx, xlo, xhi, gravity, scomp, ncomp, &
                              nscal, getForceVerbose ) &
                              bind(C, name="FORT_MAKEFORCE")

      implicit none

! In/Out
      integer :: f_lo(3), f_hi(3)
      integer :: v_lo(3), v_hi(3)
      integer :: s_lo(3), s_hi(3)
      integer :: scomp, ncomp
      integer :: nscal, getForceVerbose
      REAL_T  :: time, dx(3)
      REAL_T  :: xlo(3), xhi(3)
      REAL_T  :: gravity
      REAL_T, dimension(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),scomp:scomp+ncomp-1) :: force
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),0:dim-1) :: vel
      REAL_T, dimension(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),0:nscal-1) :: scal

! Local
      integer :: nXvel, nYvel, nZvel, nRho, nTrac, nRhoScal
      integer :: nTrac2, nTracScal, nTrac2Scal

      integer :: i, j, k, n

!     Assumes components are in the following order
      nXvel = 0
      nYvel = 1
      nRho  = 2
      nTrac = 3
      nTrac2= 4

      nRhoScal   = nRho-dim
      nTracScal  = nTrac-dim
      nTrac2Scal = nTrac2-dim

      if (scomp.eq.0) then
!
!     Do velocity forcing
!
         do k = f_lo(3), f_hi(3)
            do j = f_lo(2), f_hi(2)
               do i = f_lo(1), f_hi(1)
                  force(i,j,k,nXvel) = zero
                  force(i,j,k,nYvel) = zero
#if ( AMREX_SPACEDIM == 3 )
                  force(i,j,k,nZvel) = zero
#endif
               enddo
            enddo
         enddo
!     End of velocity forcing
      endif

      if ((scomp+ncomp).gt.AMREX_SPACEDIM) then
!        Scalar forcing
         do n = max(scomp,nRho), scomp+ncomp-1
            if (n.eq.nRho) then
               !     Density
               do k = f_lo(3), f_hi(3)
                  do j = f_lo(2), f_hi(2)
                     do i = f_lo(1), f_hi(1)
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            else if (n.eq.nTrac) then
               !     Tracer
               do k = f_lo(3), f_hi(3)
                  do j = f_lo(2), f_hi(2)
                     do i = f_lo(1), f_hi(1)
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            else
               !     Other scalar
               do k = f_lo(3), f_hi(3)
                  do j = f_lo(2), f_hi(2)
                     do i = f_lo(1), f_hi(1)
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            endif
         enddo
      endif

    end subroutine FORT_MAKEFORCE

  end module make_force_2d_moldule
