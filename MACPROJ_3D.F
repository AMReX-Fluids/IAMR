#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "MACPROJ_F.H"
#include "DIMS.H"

#define SDIM 3

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
c ::  uz           => Y edge velocity
c ::  DIMS(uz)     => index limits for uz
c ::  xarea        => area of cell faces in X dircetion
c ::  DIMS(ax)     => index limits for xarea
c ::  yarea        => area of cell faces in Y dircetion
c ::  DIMS(ay)     => index limits for yarea
c ::  zarea        => area of cell faces in Z dircetion
c ::  DIMS(az)     => index limits for zarea
c ::  vol          => cell volume
c ::  DIMS(vol)    => index limits for vol
c :: ----------------------------------------------------------
c ::
       subroutine FORT_MACDIV (dmac,DIMS(dmac),lo,hi,
     &                         ux,DIMS(ux),uy,DIMS(uy),uz,DIMS(uz),
     &                         xarea,DIMS(ax),yarea,DIMS(ay),
     &                         zarea,DIMS(az),vol,DIMS(vol))

       integer DIMDEC(dmac)
       integer lo(SDIM), hi(SDIM)
       integer DIMDEC(ux)
       integer DIMDEC(uy)
       integer DIMDEC(uz)
       integer DIMDEC(ax)
       integer DIMDEC(ay)
       integer DIMDEC(az)
       integer DIMDEC(vol)
       REAL_T  dmac(DIMV(dmac))
       REAL_T    ux(DIMV(ux))
       REAL_T    uy(DIMV(uy))
       REAL_T    uz(DIMV(uz))
       REAL_T xarea(DIMV(ax))
       REAL_T yarea(DIMV(ay))
       REAL_T zarea(DIMV(az))
       REAL_T   vol(DIMV(vol))

       integer i, j, k

       do k = lo(3), hi(3)
       do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          dmac(i,j,k) = ( xarea(i+1,j,k)*ux(i+1,j,k) - xarea(i,j,k)*ux(i,j,k)
     &              +     yarea(i,j+1,k)*uy(i,j+1,k) - yarea(i,j,k)*uy(i,j,k)
     &              +     zarea(i,j,k+1)*uz(i,j,k+1) - zarea(i,j,k)*uz(i,j,k)
     &                )/vol(i,j,k)
       enddo
       enddo
       enddo

       return
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

       integer DIMDEC(s)
       integer lo(SDIM), hi(SDIM)
       integer DIMDEC(dmac)
       integer nvar
       REAL_T  mult
       REAL_T  snew(DIMV(s),nvar)
       REAL_T  sold(DIMV(s),nvar)
       REAL_T  dmac(DIMV(dmac))

       integer i, j, k, n, lenx, leny, lenz

       lenx = hi(1)-lo(1)+1
       leny = hi(2)-lo(2)+1
       lenz = hi(3)-lo(3)+1

       if (lenx .ge. leny .and. lenx .ge. lenz) then
          do n = 1, nvar
          do k = lo(3), hi(3)
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             snew(i,j,k,n) = snew(i,j,k,n) + mult*dmac(i,j,k)*sold(i,j,k,n)
          enddo
          enddo
          enddo
          enddo
       else if (leny .ge. lenx .and. leny .ge. lenz) then
          do n = 1, nvar
          do k = lo(3), hi(3)
          do i = lo(1), hi(1)
          do j = lo(2), hi(2)
             snew(i,j,k,n) = snew(i,j,k,n) + mult*dmac(i,j,k)*sold(i,j,k,n)
          enddo
          enddo
          enddo
          enddo
       else 
          do n = 1, nvar
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
          do k = lo(3), hi(3)
             snew(i,j,k,n) = snew(i,j,k,n) + mult*dmac(i,j,k)*sold(i,j,k,n)
          enddo
          enddo
          enddo
          enddo
       endif
c
c     end of routine
c
       return
       end

