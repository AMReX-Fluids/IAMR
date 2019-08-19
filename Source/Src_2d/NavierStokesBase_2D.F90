
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <NAVIERSTOKESBASE_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 2


module navierstokesbase_2d_module
  
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
    integer,          intent(in   ) :: lo(1:2),hi(1:2)
    integer,          intent(in   ) :: Slo(1:2),Shi(1:2)
    integer,          intent(in   ) :: mlo(1:2),mhi(1:2)
    integer,          intent(in   ) :: mask(mlo(1):mhi(1),mlo(2):mhi(2))
    REAL_T, intent(inout) :: S(Slo(1):Shi(1),Slo(2):Shi(2),1:nc)
    REAL_T, intent(in   ) :: b(1:nc)
    integer :: i,j,n

    do n=1,nc
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             if (mask(i,j).eq.bval) S(i,j,n)=b(n)
          enddo
       enddo
    enddo

  end subroutine fort_set_body_state

#endif

end module navierstokesbase_2d_module
