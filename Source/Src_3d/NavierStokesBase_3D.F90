
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <NAVIERSTOKESBASE_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 3


module navierstokesbase_3d_module
  
  implicit none

  private
  
#ifdef AMREX_USE_EB
  public fort_set_body_state
#endif
  
contains


#ifdef AMREX_USE_EB

  subroutine fort_set_body_state(lo, hi, S, Slo, Shi, mask, mlo, mhi, b, nc, bval) &
                                         bind(C,name="fort_set_body_state")

    integer,          intent(in   ) :: nc, bval
    integer,          intent(in   ) :: lo(1:3),hi(1:3)
    integer,          intent(in   ) :: Slo(1:3),Shi(1:3)
    integer,          intent(in   ) :: mlo(1:3),mhi(1:3)
    integer,          intent(in   ) :: mask(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
    REAL_T, intent(inout) :: S(Slo(1):Shi(1),Slo(2):Shi(2),Slo(3):Shi(3),1:nc)
    REAL_T, intent(in   ) :: b(1:nc)
    integer :: i,j,k,n

    do n=1,nc
      do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             if (mask(i,j,k).eq.bval) S(i,j,k,n)=b(n)
          enddo
       enddo
      enddo
    enddo

  end subroutine fort_set_body_state

#endif

end module navierstokesbase_3d_module
