
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PROJECTION_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 3

module projection_3d_module

  implicit none

  private

  public ::  anelcoeffmpy

contains

      subroutine anelcoeffmpy (lo,hi,a,DIMS(a),domlo,domhi, &
                               anel_coeff,anel_lo,anel_hi,bogus_value,mult) &
                               bind(C, name="anelcoeffmpy")
!c 
!c     multiply A by anel_coeff
!c
!c 
!c     NOTE: THIS ROUTINE HAS BEEN MODIFIED SO THAT ALL VALUES
!c           OUTSIDE THE DOMAIN ARE SET TO BOGUS VALUE
      !c-- only sets boudnary cells at top and bottom to Bogus Val,
      !c   boundary cells on the sides are left alone

      implicit none
      integer    lo(SDIM),hi(SDIM)
      integer    DIMDEC(a)
      integer    domlo(3), domhi(3)
      REAL_T     a(DIMV(a))
      integer    mult, anel_lo, anel_hi
      REAL_T     anel_coeff(anel_lo:anel_hi)
      REAL_T     bogus_value

      integer i, j, k
      integer klo,khi

      klo = lo(3)
      khi = hi(3)
      
      if (lo(3) .lt. domlo(3)) then
         klo = domlo(3)
         do k = lo(3), domlo(3)-1
         do j = lo(2), hi(2)
         do i = lo(1),hi(1)
           a(i,j,k) = bogus_value
         end do
         end do
         end do
      end if
      if (hi(3) .gt. domhi(3)) then
         khi = domhi(3)
         do k = domhi(3)+1, hi(3)
         do j = lo(2), hi(2)
         do i = lo(1),hi(1)
           a(i,j,k) = bogus_value
         end do
         end do
         end do
      end if

      if (mult .eq. 1) then
         do k = klo,khi
         do j = lo(2),hi(2)
         do i = lo(1),hi(1)
           a(i,j,k) = a(i,j,k) * anel_coeff(k)
         end do
         end do
         end do
      else if (mult .eq. 0) then
         do k = klo,khi
         do j = lo(2),hi(2)
         do i = lo(1),hi(1)
           a(i,j,k) = a(i,j,k) / anel_coeff(k)
         end do
         end do
         end do
      else 
         print *,'BOGUS MULT IN ANELCOEFFMULT ',mult
         stop
      end if

      end subroutine anelcoeffmpy

end module projection_3d_module
