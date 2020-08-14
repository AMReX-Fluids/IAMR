#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PROB_NS_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 3

module MakeForce_3D_module

  implicit none

  private

  public :: FORT_MAKEFORCE

contains

! ::: -----------------------------------------------------------
!
!     This is a dummy routine to add the forcing terms to the momentum equation
!     It should be overwritten in the local run folder.  
!
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
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),0:SDIM-1) :: vel
      REAL_T, dimension(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),0:nscal-1) :: scal

   end subroutine FORT_MAKEFORCE

end module MakeForce_3D_module
