c $Id: MACPROJ_2D.F,v 1.2 1997-07-22 21:57:05 car Exp $
#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "MACPROJ_F.H"
#include "DIMS.H"

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

       integer i, j, lenx, leny

       lenx = hi(1)-lo(1)+1
       leny = hi(2)-lo(2)+1
       if (lenx .ge. leny) then
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dmac(i,j) = ( xarea(i+1,j)*ux(i+1,j) - xarea(i,j)*ux(i,j)
     &                 +   yarea(i,j+1)*uy(i,j+1) - yarea(i,j)*uy(i,j)
     &                   )/vol(i,j)

          enddo
          enddo
       else
          do i = lo(1), hi(1)
          do j = lo(2), hi(2)
             dmac(i,j) = ( xarea(i+1,j)*ux(i+1,j) - xarea(i,j)*ux(i,j)
     &                 +   yarea(i,j+1)*uy(i,j+1) - yarea(i,j)*uy(i,j)
     &                   )/vol(i,j)
          enddo
          enddo
       endif

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
       integer DIMDEC(dmac)
       integer nvar
       integer lo(SDIM), hi(SDIM)
       REAL_T  snew(DIMV(s),nvar)
       REAL_T  sold(DIMV(s),nvar)
       REAL_T  dmac(DIMV(dmac))
       REAL_T  mult

       integer i, j, n, lenx, leny

       lenx = hi(1)-lo(1)+1
       leny = hi(2)-lo(2)+1
       if (lenx .ge. leny) then
          do n = 1, nvar
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             snew(i,j,n) = snew(i,j,n) + mult*dmac(i,j)*sold(i,j,n)
          enddo
          enddo
          enddo
       else
          do n = 1, nvar
          do i = lo(1), hi(1)
          do j = lo(2), hi(2)
             snew(i,j,n) = snew(i,j,n) + mult*dmac(i,j)*sold(i,j,n)
          enddo
          enddo
          enddo
       endif

       return
       end



c=========================================================


      subroutine FORT_MACPHIBC(
     &     DIMS(u),uhalfx,DIMS(div),divu,DIMS(rho),rho,
     &     DIMS(rcen),rcen,DIMS(rhalf),rhalf,hx,
     &     DIMS(phi),phi)
c
c ... compute the value of phi for mac to be used at an 
c     upper, y-outflow face
c
      implicit none

      integer DIMDEC(u)
      integer DIMDEC(div)
      integer DIMDEC(rho)
      integer DIMDEC(rcen)
      integer DIMDEC(rhalf)
      integer DIMDEC(phi)
      REAL_T    phi(DIMV(phi))
      REAL_T uhalfx(DIMV(u))
      REAL_T   divu(DIMV(div))
      REAL_T   rho(DIMV(rho))
      REAL_T  rcen(DIM1(rcen))
      REAL_T rhalf(DIM1(rhalf))
      REAL_T hx
c
      REAL_T fst,fstm
      REAL_T rhoedge, divu_ave, vtot
      integer i,j,is,ie,js,je
      REAL_T max_divu, min_divu
c
#if 0
      return
#endif
      is = ARG_L1(phi)+1
      ie = ARG_H1(phi)-1
      js = ARG_L2(phi)
      je = ARG_H2(phi)


        divu_ave = zero
        vtot = zero
        max_divu = zero
        min_divu = zero
        do j=je-1,je-1
        do i=is,ie 
           divu_ave = divu_ave+rcen(i)*divu(i,j)
           vtot = vtot+rcen(i)
           max_divu = max(max_divu,divu(i,j))
           min_divu = min(min_divu,divu(i,j))
        enddo
        enddo
        divu_ave = divu_ave/vtot

      if(max_divu.eq.zero.and.min_divu.eq.zero) then
        return
      endif

      fstm = zero
      phi(is,je) = zero

      do i = is+1,ie 
         fst = (rhalf(i-1)*fstm + 
     $          hx*rcen(i-1)*(divu(i-1,je-1)-divu_ave))/ rhalf(i)
         rhoedge = (rho(i,je-1)+rho(i-1,je-1))*half
         phi(i,je) = phi(i-1,je) + hx * rhoedge * 
     &                             (uhalfx(i,je-1) - fst)
         fstm = fst
      enddo

      return
      end


