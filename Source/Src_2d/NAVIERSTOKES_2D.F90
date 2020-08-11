
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif
  
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <NAVIERSTOKES_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 2


module navierstokes_2d_module
  
  implicit none

  private 

  public fort_maxval, cen2edg
  
contains

!c :: ----------------------------------------------------------
!c :: MAXVAL
!c ::             maxval = max{ rho(i,j) }
!c ::
!c :: ----------------------------------------------------------
!c ::
     subroutine fort_maxval(rho,DIMS(rho),DIMS(grid),mxval) &
          bind(C,name="fort_maxval")

       implicit none
       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  rho(DIMV(rho))
       REAL_T  mxval

       integer i,j

       mxval = -Huge(zero)

       do j = ARG_L2(grid), ARG_H2(grid)
          do i = ARG_L1(grid), ARG_H1(grid)
             mxval = max(mxval, rho(i,j))
          end do
       end do

     end subroutine fort_maxval

!c ::
!c :: ----------------------------------------------------------
!c :: This routine fills an edge-centered fab from a cell-centered
!c :: fab using simple linear interpolation.
!c ::
!c :: INPUTS / OUTPUTS:
!c ::  lo,hi      => index limits of the of the cell-centered fab
!c ::  DIMS(cfab) => index limits of the cell-centered fab
!c ::  cfab       => cell-centered data
!c ::  DIMS(efab) => index limits of the edge-centered fab
!c ::  efab       => edge-centered fab to fill
!c ::  n!c         => Number of components in the fab to fill
!c ::  dir        => direction data needs to be shifted to get to edges
!c :: ----------------------------------------------------------
!c ::
      subroutine cen2edg(lo, hi, &
          DIMS(cfab), cfab,&
          DIMS(efab), efab, nc, dir,&
          isharm) bind(C,name="cen2edg")
      implicit none
      integer lo(SDIM), hi(SDIM), nc, dir, isharm
      integer DIMDEC(cfab)
      integer DIMDEC(efab)
      REAL_T  cfab(DIMV(cfab), nc)
      REAL_T  efab(DIMV(efab), nc)

      integer i,j,n

      if ( isharm .eq. 0 ) then
         if (dir .EQ. 0) then
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     efab(i,j,n) = half*(cfab(i,j,n) + cfab(i-1,j,n))
                  end do
               end do
            end do
         else
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     efab(i,j,n) = half*(cfab(i,j,n) + cfab(i,j-1,n))
                  end do
               end do
            end do
         end if
      else
         if (dir .EQ. 0) then
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if((cfab(i,j,n) * cfab(i-1,j,n)).gt.zero)then
                        efab(i,j,n)&
                            = 2*(cfab(i,j,n) * cfab(i-1,j,n))/&
                            (cfab(i,j,n) + cfab(i-1,j,n))
                     else
                        efab(i,j,n)=zero
                     endif
                  end do
               end do
            end do
         else
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if((cfab(i,j,n) * cfab(i,j-1,n)).gt.zero)then
                        efab(i,j,n)&
                            = 2*(cfab(i,j,n) * cfab(i,j-1,n))/&
                            (cfab(i,j,n) + cfab(i,j-1,n))
                     else
                        efab(i,j,n)=zero
                     endif
                  end do
               end do
            end do
         end if
      end if
    end subroutine cen2edg

  end module navierstokes_2d_module
