#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <MACPROJ_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 2


module macproj_2d_module
  
  implicit none

  private

  public macdiv, fort_scalearea

contains

! c :: ----------------------------------------------------------
! c :: MACDIV:  compute the MAC divergence in generalized coordinates
! c ::
! c :: INPUTS / OUTPUTS:
! c ::  dmac        <=  MAC divergence (cell centered)
! c ::  DIMS(dmac)   => index limits for dmac
! c ::  lo,hi        => index limits of grid interior
! c ::  ux           => X edge velocity
! c ::  DIMS(ux)     => index limits for ux
! c ::  uy           => Y edge velocity
! c ::  DIMS(uy)     => index limits for uy
! c ::  xarea        => area of cell faces in X dircetion
! c ::  DIMS(ax)     => index limits for ax
! c ::  yarea        => area of cell faces in Y dircetion
! c ::  DIMS(ay)     => index limits for ay
! c ::  vol          => cell volume
! c ::  DIMS(vol)    => index limits for vol
! c :: ----------------------------------------------------------
! c ::
       subroutine macdiv (dmac,DIMS(dmac),lo,hi, &
                              ux,DIMS(ux),uy,DIMS(uy), &
                              xarea,DIMS(ax),yarea,DIMS(ay), &
                              vol,DIMS(vol)) bind(C,name="macdiv")

       implicit none
       integer DIMDEC(dmac)
       integer lo(SDIM), hi(SDIM)
       integer DIMDEC(ux)
       integer DIMDEC(uy)
       integer DIMDEC(ax)
       integer DIMDEC(ay)
       integer DIMDEC(vol)
       REAL_T  dmac(DIMV(dmac))
       REAL_T    ux(DIMV(ux))
       REAL_T    uy(DIMV(uy))
       REAL_T xarea(DIMV(ax))
       REAL_T yarea(DIMV(ay))
       REAL_T   vol(DIMV(vol))

       integer i, j

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dmac(i,j) = ( xarea(i+1,j)*ux(i+1,j) - xarea(i,j)*ux(i,j)&
                      +   yarea(i,j+1)*uy(i,j+1) - yarea(i,j)*uy(i,j)&
                 )/vol(i,j)

          end do
       end do

     end subroutine macdiv

! c :: ----------------------------------------------------------
! c :: SCALEAREA
! c ::          area = area * anel_coeff
! c ::                 OR 
! c ::          area = area / anel_coeff
! c :: ----------------------------------------------------------

     subroutine fort_scalearea (lo,hi,vbxhi, &
            xarea,DIMS(ax),yarea,DIMS(ay), &
            anel_coeff,anel_coeff_lo,anel_coeff_hi, &
            mult) bind(C,name="fort_scalearea")

       implicit none
       integer lo(SDIM), hi(SDIM), vbxhi(SDIM)
       integer anel_coeff_lo,anel_coeff_hi
       integer DIMDEC(ax)
       integer DIMDEC(ay)
       REAL_T xarea(DIMV(ax))
       REAL_T yarea(DIMV(ay))
       REAL_T anel_coeff(anel_coeff_lo:anel_coeff_hi)
       integer mult

       integer i,j
       integer lo2, hi1, hi2

       ! check to see if we're at a cc box boundary
       ! if so, we need to include 1 more point at high end because
       ! area is nodal in one dim
       if (hi(1) .eq. vbxhi(1)) then
          hi1 = hi(1)+1
       else
          hi1=hi(1)
       endif
       if (hi(2) .eq. vbxhi(2)) then
          hi2 = hi(2)+1
       else
          hi2=hi(2)
       endif

       ! do something different at the very bottom
       if (lo(2) .eq. 0) then
          lo2 = 1
       else
          lo2 = lo(2)
       end if


       if (mult .eq. 1) then

          do j = lo(2), hi(2)
          do i = lo(1), hi1
             xarea(i,j) =  xarea(i,j) * anel_coeff(j)
          end do
          end do

          do j = lo2, hi2
          do i = lo(1), hi(1)
             yarea(i,j) =  yarea(i,j) * 0.5d0 * (anel_coeff(j)+anel_coeff(j-1))
          end do
          end do

          if (lo(2) .eq. 0) then
             j = lo(2)
             do i = lo(1), hi(1)
                yarea(i,j) =  yarea(i,j) * anel_coeff(j-1)
             end do
          end if

       else if (mult .eq. -1) then

          do j = lo(2), hi(2)
          do i = lo(1), hi1
             xarea(i,j) =  xarea(i,j) / anel_coeff(j)
          end do
          end do

          do j = lo2, hi2
          do i = lo(1), hi(1)
             yarea(i,j) =  yarea(i,j) / (0.5d0 * (anel_coeff(j)+anel_coeff(j-1)))
          end do
          end do

          if (lo(2) .eq. 0) then
             j = lo(2)
             do i = lo(1), hi(1)
                yarea(i,j) =  yarea(i,j) / anel_coeff(j-1) 
             end do
          end if

       else 
 
          print *,'BOGUS MULT IN SCALEAREA '
          stop

       end if

     end subroutine fort_scalearea

   end module macproj_2d_module
