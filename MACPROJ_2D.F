c
c $Id: MACPROJ_2D.F,v 1.16 2008-02-08 18:39:48 lijewski Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "MACPROJ_F.H"
#include "ArrayLim.H"

#define SDIM 2

c :: ----------------------------------------------------------
c :: MACDIV:  compute the MAC divergence in generalized coordinates
c ::
c :: INPUTS / OUTPUTS:
c ::  dmac        <=  MAC divergence (cell centered)
c ::  DIMS(dmac)   => index limits for dmac
c ::  lo,hi        => index limits of grid interior
c ::  ux           => X edge velocity
c ::  DIMS(ux)     => index limits for ux
c ::  uy           => Y edge velocity
c ::  DIMS(uy)     => index limits for uy
c ::  xarea        => area of cell faces in X dircetion
c ::  DIMS(ax)     => index limits for ax
c ::  yarea        => area of cell faces in Y dircetion
c ::  DIMS(ay)     => index limits for ay
c ::  vol          => cell volume
c ::  DIMS(vol)    => index limits for vol
c :: ----------------------------------------------------------
c ::
       subroutine FORT_MACDIV (dmac,DIMS(dmac),lo,hi,
     &                         ux,DIMS(ux),uy,DIMS(uy),
     &                         xarea,DIMS(ax),yarea,DIMS(ay),
     &                         vol,DIMS(vol))

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
             dmac(i,j) = ( xarea(i+1,j)*ux(i+1,j) - xarea(i,j)*ux(i,j)
     &                 +   yarea(i,j+1)*uy(i,j+1) - yarea(i,j)*uy(i,j)
     &            )/vol(i,j)

          end do
       end do

       end


c :: ----------------------------------------------------------
c :: UCORUPDT
c ::          Snew = Snew + mult*Sold*dmac
c ::
c :: INPUTS / OUTPUTS:
c ::  snew        <=  state vector at new time
c ::  sold         => state vector at old time
c ::  DIMS(s)      => index limits for sold and snew
c ::  lo,hi        => index limits of grid interior
c ::  nvar         => number of components in state vector
c ::  dmac         => MAC divergence of Ucorr
c ::  DIMS(dmac)   => index limits for dmac
c ::  mult         => scalar multiplier
c :: ----------------------------------------------------------
c ::
       subroutine FORT_UCORUPDT (snew,sold,DIMS(s),nvar,lo,hi,
     &                           dmac,DIMS(dmac),mult)

       implicit none
       integer DIMDEC(s)
       integer DIMDEC(dmac)
       integer nvar
       integer lo(SDIM), hi(SDIM)
       REAL_T  snew(DIMV(s),nvar)
       REAL_T  sold(DIMV(s),nvar)
       REAL_T  dmac(DIMV(dmac))
       REAL_T  mult

       integer i, j, n

       do n = 1, nvar
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                snew(i,j,n) = snew(i,j,n) + mult*dmac(i,j)*sold(i,j,n)
             end do
          end do
       end do

       end

c :: ----------------------------------------------------------
c :: SCALEAREA
c ::          area = area * anel_coeff
c ::                 OR 
c ::          area = area / anel_coeff
c :: ----------------------------------------------------------

       subroutine FORT_SCALEAREA (xarea,DIMS(ax),yarea,DIMS(ay),
     &                            anel_coeff,anel_coeff_lo,anel_coeff_hi,lo,hi,mult)

       implicit none
       integer lo(SDIM), hi(SDIM)
       integer anel_coeff_lo,anel_coeff_hi
       integer DIMDEC(ax)
       integer DIMDEC(ay)
       REAL_T xarea(DIMV(ax))
       REAL_T yarea(DIMV(ay))
       REAL_T anel_coeff(anel_coeff_lo:anel_coeff_hi)
       integer mult

       integer i,j

       if (mult .eq. 1) then

          do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             xarea(i,j) =  xarea(i,j) * anel_coeff(j)
          end do
          end do

          do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             yarea(i,j) =  yarea(i,j) * 0.5d0 * (anel_coeff(j)+anel_coeff(j-1))
          end do
          end do

          if (lo(2) .eq. 0) then
             j = lo(2)
             do i = lo(1), hi(1)
                yarea(i,j) =  yarea(i,j) * anel_coeff(j-1) / 
     $                        (0.5d0 * (anel_coeff(j)+anel_coeff(j-1)))
             end do
          end if

       else if (mult .eq. -1) then

          do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             xarea(i,j) =  xarea(i,j) / anel_coeff(j)
          end do
          end do

          do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             yarea(i,j) =  yarea(i,j) / (0.5d0 * (anel_coeff(j)+anel_coeff(j-1)))
          end do
          end do

          if (lo(2) .eq. 0) then
             j = lo(2)
             do i = lo(1), hi(1)
                yarea(i,j) =  yarea(i,j) / anel_coeff(j-1) *
     $                        (0.5d0 * (anel_coeff(j)+anel_coeff(j-1)))
             end do
          end if

       else 
 
          print *,'BOGUS MULT IN SCALEAREA '
          stop

       end if

       end
