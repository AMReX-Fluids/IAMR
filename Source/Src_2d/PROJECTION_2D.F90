
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PROJECTION_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 2

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define SMALL 1.0e-10
#else
#define SMALL 1.0d-10
#endif


module projection_2d_module
  
  implicit none

  private 

  public :: radmpyscal, radmpyvel, fort_raddiv
  
contains

    subroutine radmpyscal(lo,hi,a,DIMS(a),domlo,domhi,r,rlo,rhi) &
         bind(C,name="radmpyscal")
!c 
!c     multiply A by Radius r
!c
!c 
!c     NOTE: THIS ROUTINE HAS BEEN MODIFIED SO THAT ALL VALUES
!c           OUTSIDE THE DOMAIN ARE SET TO ZERO
!c
      implicit none
      integer    rlo,rhi
      integer    DIMDEC(a)
      integer    domlo(2), domhi(2)
      integer    lo(SDIM),hi(SDIM)
      REAL_T     a(DIMV(a))
      REAL_T     r(rlo:rhi)
      REAL_T     bogus_value

      integer i, j

      do j = lo(2),hi(2)
         do i = lo(1),min(domhi(1),hi(1))
           a(i,j) = r(i)*a(i,j)
         end do
      end do

      !c     NOTE: We used to set these to bogus_value to be sure that we
      !c           didn't use them. But now in the divu routine in the F90
      !c           solvers we need to include these values in the stencil
      !c           because they might contain inflow values, for example, 
      !c           and the only test is on the BC for the pressure solve,
      !c           which doesn't differentiate between inflow, reflecting
      !c           and symmetry.

      if (lo(1) .lt. domlo(1)) then
         do j = lo(2),hi(2)
         do i = lo(1), domlo(1)-1
           a(i,j) = 0.d0
         end do
         end do
      end if

      if (hi(1) .gt. domhi(1)) then
         do j = lo(2),hi(2)
         do i = domhi(1)+1, hi(1)
           a(i,j) = 0.d0
         end do
         end do
      end if

    end subroutine radmpyscal

    subroutine radmpyvel(lo,hi,a,DIMS(a),domlo,domhi,r,rlo,rhi,ndim) &
         bind(C,name="radmpyvel")
!c 
!c     multiply A by Radius r
!c
!c 
!c     NOTE: THIS ROUTINE HAS BEEN MODIFIED SO THAT ALL VALUES
!c           OUTSIDE THE DOMAIN ARE SET TO ZERO
!c
      implicit none
      integer    ndim
      integer    rlo,rhi
      integer    DIMDEC(a)
      integer    domlo(2), domhi(2)
      integer    lo(SDIM),hi(SDIM)
      REAL_T     a(DIMV(a))
      REAL_T     r(rlo:rhi)

      integer    i,j 
      REAL_T     dr

      do j = lo(2),hi(2)
         do i = lo(1), min(hi(1),domhi(1))
           a(i,j) = r(i)*a(i,j)
         end do
      end do

!c     NOTE: We used to set these to bogus_value to be sure that we didn't use them.
!c           But now in the divu routine in the F90 solvers we need to include these
!c           values in the stencil because they might contain inflow values, for
!c           example, and the only test is on the BC for the pressure solve, which 
!c           doesn't differentiate between inflow, reflecting and symmetry.

      if (lo(1) .lt. domlo(1)) then
         do j = lo(2),hi(2)
         do i = lo(1), domlo(1)-1
           a(i,j) = 0.d0
         end do
         end do
      end if

!c     Here we only multiply a possibly inflow x-velocity from the hi-r side 
      if (ndim .eq. 0) then
         if (hi(1) .gt. domhi(1)) then
            dr = r(hi(1)) - r(hi(1)-1)
            do j = lo(2),hi(2)
            do i = domhi(1)+1, hi(1)
              a(i,j) = (r(domhi(1)) + (i-domhi(1))*dr - 0.5*dr)*a(i,j)
            end do
            end do
         end if
      else
         if (hi(1) .gt. domhi(1)) then
            do j = lo(2),hi(2)
            do i = domhi(1)+1, hi(1)
              a(i,j) = 0.d0
            end do
            end do
         end if
      end if

    end subroutine radmpyvel

    subroutine fort_raddiv(lo,hi,a,DIMS(a),domlo,domhi,r,rlo,rhi,bogus_value)&
         bind(C,name="fort_raddiv")
!c 
!c     divide A by Radius r
!c
!c 
!c     NOTE: THIS ROUTINE HAS BEEN MODIFIED SO THAT ALL VALUES
!c           OUTSIDE THE DOMAIN ARE SET TO BOGUS VALUE
!c
      implicit none
      integer    rlo,rhi
      integer    DIMDEC(a)
      integer    domlo(2), domhi(2)
      integer    lo(SDIM),hi(SDIM)
      REAL_T     a(DIMV(a))
      REAL_T     r(rlo:rhi)
      REAL_T     bogus_value

      integer i, j

      do j = lo(2), hi(2)
         do i = lo(1),min(domhi(1),hi(1))
           a(i,j) = a(i,j)/r(i)
         end do
      end do

      if (lo(1) .lt. domlo(1)) then
         do j = lo(2),hi(2)
         do i = lo(1), domlo(1)-1
           a(i,j) = bogus_value
         end do
         end do
      end if

      if (hi(1) .gt. domhi(1)) then
         do j = lo(2),hi(2)
         do i = domhi(1)+1, hi(1)
           a(i,j) = bogus_value
         end do
         end do
      end if

    end subroutine fort_raddiv

  end module projection_2d_module
