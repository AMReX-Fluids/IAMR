c
c $Id: PROJECTION_2D.F,v 1.4 1997-12-12 22:44:09 lijewski Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "PROJECTION_F.H"
#include "ArrayLim.H"

#define SDIM 2


c=========================================================

       subroutine FORT_FILTERP (
     &     pnew,scratch,DIMS(p),
     &     lo,hi,dx,fac)
c
c     ----------------------------------------------------------
c     filter the pressure coming out of the level projection
c     
c     pnew := pnew - fac * lap_5pt(pnew)
c     
c     INPUTS / OUTPUTS:
c     pnew       => 
c     pnew      <=> 
c     p_lo,p_hi  => index limits of pressure array
c     lo,hi      => index limits of cell centered interior
c     fac        => fudge factor for filtering
c     
c     ----------------------------------------------------------
c 
       integer    DIMDEC(p)
       REAL_T     pnew(DIMV(p))
       REAL_T     scratch(DIMV(p))
       integer    lo(2), hi(2)
       REAL_T     dx(2)
       REAL_T     fac
       integer i, j
c
       if ( abs(dx(1)-dx(2)) .gt. 1.e-8 .and. fac .gt. 1.e-4 )then
          print *, 'is filter properly formatted for dx != dy ?'
       endif
c
       do j = lo(2)+1, hi(2)
          do i = lo(1)+1, hi(1)
             scratch(i,j) = four * pnew(i,j) - 
     $         (pnew(i-1,j) + pnew(i+1,j) + pnew(i,j-1) + pnew(i,j+1) )
             scratch(i,j) = - eighth * scratch(i,j) * fac
          enddo
       enddo
c
       do j = lo(2)+1, hi(2)
          do i = lo(1)+1, hi(1)
             pnew(i,j) = pnew(i,j) + scratch(i,j)
          enddo
       enddo
c
c     end of routine
c
       return 
       end


c=========================================================

       subroutine FORT_ACCEL_TO_VEL( lo, hi, 
     &     uold,DIMS(uold),
     &     dt,
     &     unew,DIMS(unew) )
c
c     This function converts unew into a velocity via
c     Unew = Uold + alpha*Unew
c
       integer    lo(SDIM), hi(SDIM)
       REAL_T     dt
       integer    DIMDEC(uold),DIMDEC(unew)
       REAL_T     uold(DIMV(uold),SDIM)
       REAL_T     unew(DIMV(unew),SDIM)
c
       integer i, j, n
c
       do n = 1, SDIM
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1
                unew(i,j,n) = uold(i,j,n) + dt*unew(i,j,n)
             enddo
          enddo
       enddo
c
c     end of routine
c
       return
       end



c=========================================================

      subroutine FORT_VEL_TO_ACCEL( lo, hi, 
     &     unew,DIMS(unew),
     &     uold,DIMS(uold),
     &     dt )
c     
c     This function converts unew into an acceleration
c     
      integer    lo(SDIM), hi(SDIM)
      REAL_T     dt
      integer    DIMDEC(unew),DIMDEC(uold)
      REAL_T     uold(DIMV(uold),SDIM)
      REAL_T     unew(DIMV(unew),SDIM)
c     
      integer i, j, n
c
      do n = 1, SDIM
         do j = lo(2)-1, hi(2)+1
            do i = lo(1)-1, hi(1)+1
               unew(i,j,n) = (unew(i,j,n)-uold(i,j,n))/dt
            enddo
         enddo
      enddo
c
c     end of routine
c
      return
      end


c=========================================================

      subroutine FORT_PROJ_UPDATE(
     &     boxlo, boxhi, nvar, ngrow,
     &     un, DIMS(un),
     &     alpha,
     &     uo, DIMS(uo) )
c     
c     This function updates un via un = un + alpha*uo
c     The loop bounds are determined in the C++
c     
      integer    boxlo(SDIM), boxhi(SDIM), nvar, ngrow
      REAL_T     alpha
      integer    DIMDEC(un),DIMDEC(uo)
      REAL_T     un(DIMV(un),nvar)
      REAL_T     uo(DIMV(uo),nvar)
c     
      integer i, j, n
c     
      do n = 1, nvar
         do j = boxlo(2), boxhi(2)
            do i = boxlo(1), boxhi(1)
               un(i,j,n) = un(i,j,n) + alpha*uo(i,j,n)
            enddo
         enddo
      enddo
c
c     end of routine
c
      return
      end
      


c=========================================================

       subroutine FORT_RADMPY(a,DIMS(grid),ng,r,nr,n)
c 
c     multiply A by a Radius r
c
       integer    ng,nr
       integer    n(SDIM)
       integer    DIMDEC(grid)
       REAL_T     a(ARG_L1(grid)-ng:ARG_H1(grid)+ng, 
     $              ARG_L2(grid)-ng:ARG_H2(grid)+ng)
       REAL_T     r(ARG_L1(grid)-nr:ARG_H1(grid)+nr)
       integer i, j
c
       do j = ARG_L2(grid)-n(2), ARG_H2(grid)+n(2)
          do i = ARG_L1(grid)-n(1), ARG_H1(grid)+n(1)
	     a(i,j) = r(i)*a(i,j)
          enddo
       enddo
c
c     end of routine
c
       return
       end


c=========================================================


       subroutine FORT_RADDIV(a,DIMS(grid),ng,r,nr,n)
c 
c     divide A by a Radius r
c
       integer    ng,nr
       integer    n(SDIM)
       integer    DIMDEC(grid)
       REAL_T     a(ARG_L1(grid)-ng:ARG_H1(grid)+ng, 
     $              ARG_L2(grid)-ng:ARG_H2(grid)+ng)
       REAL_T     r(ARG_L1(grid)-nr:ARG_H1(grid)+nr)
c
       integer i, j
c
       do j = ARG_L2(grid)-n(2), ARG_H2(grid)+n(2)
          do i = ARG_L1(grid)-n(1), ARG_H1(grid)+n(1)
	    a(i,j) = a(i,j)/r(i)
          enddo
       enddo
c
c     end of routine
c
       return
       end

c=========================================================


      subroutine FORT_HGN2C(
     &     isrz,lrweighted, DIMS(nodedat), nodedat,
     &     DIMS(ccdat), lo, hi, ccdat)

c     ----------------------------------------------------------
c     HGN2C
c     averages node centered data to cell centers for use in 
c     holy grail projection

      implicit none
      integer isrz,lrweighted
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(ccdat)
      integer DIMDEC(nodedat)
      REAL_T  nodedat(DIMV(nodedat))
      REAL_T  ccdat(DIMV(ccdat))

      integer i,j

c ---- begin ----------

      if (ARG_H1(ccdat)   .lt. lo(1) .or. 
     &    ARG_L1(ccdat)   .gt. hi(1) .or. 
     &    ARG_H2(ccdat)   .lt. lo(2) .or. 
     &    ARG_L2(ccdat)   .gt. hi(2) .or. 
     &    ARG_H1(nodedat) .lt. lo(1)+1 .or. 
     &    ARG_L1(nodedat) .gt. hi(1) .or. 
     &    ARG_H2(nodedat) .lt. lo(2)+1 .or. 
     &    ARG_L2(nodedat) .gt. hi(2) ) then 
        write(6,*) "FORT_HG_CELL_TO_NODE: bad index limits"
        stop
      endif

      if(isrz.eq.1.and.lrweighted.ne.1)then
        write(6,*)
     &    'FORT_HGN2C : isrz=1 and lrweighted!=1 not yet implemented'
#ifndef	BL_NO_FORT_FLUSH
        call flush(6)
#endif
        stop
      endif

      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          ccdat(i,j) = fourth*(nodedat(i,j)+nodedat(i+1,j)+
     &                         nodedat(i,j+1)+nodedat(i+1,j+1))
        enddo
      enddo

      return
      end

c=========================================================

      subroutine FORT_HGC2N(
     &     nghost, DIMS(dat), dat, rcen,
     &     DIMS(rhs), rhs,
     &     domlo, domhi, lowfix, hifix, dr,
     &     extrap_edges, extrap_corners, is_rz) 
c
c     ----------------------------------------------------------
c     HGC2N
c     averages cell centered data to nodes for use in 
c     holy grail projection
c     
c     INPUTS / OUTPUTS:
c     nghost      => indicates buffer of rhs that does not need values
c     dat         => cell centered array to be averaged
c     DIMS(dat)   => index limits of dat
c     rcen        => r-coordinate cell centers if geoem is r-z; 
c     otherwise, should be 1
c     rhslo,rhshi => index extents of rhs
c     rhs         <= node centered array with results
c     ----------------------------------------------------------
c 
      implicit none
      integer nghost, is_rz, lowfix, hifix
      integer domlo(SDIM), domhi(SDIM)
      integer DIMDEC(dat)
      integer DIMDEC(rhs)
      integer extrap_edges, extrap_corners
      REAL_T  dr
      REAL_T  rcen(DIM1(dat))
      REAL_T  dat(DIMV(dat))
      REAL_T  rhs(DIMV(rhs))

      integer i, j
      REAL_T  factor, rhi, rlo

      if (extrap_edges .eq. 1) then

c ... fix up bcs at domain edges

        if(ARG_L1(rhs)+nghost .eq. domlo(1)) then
          i = domlo(1)
          do j = ARG_L2(rhs)+nghost, ARG_H2(rhs)-nghost
             dat(i-1,j) = dat(i,j)
          enddo
        endif

        if(ARG_L1(rhs)-nghost .eq. domhi(1)+1) then
          i = domhi(1)
          do j = ARG_L2(rhs)+nghost, ARG_H2(rhs)-nghost
            dat(i+1,j) = dat(i,j)
          enddo
        endif

        if(ARG_L2(rhs)+nghost .eq. domlo(2)) then
          j = domlo(2)
          do i = ARG_L1(rhs)+nghost, ARG_H1(rhs)-nghost
            dat(i,j-1) = dat(i,j)
          enddo
        endif

        if(ARG_H2(rhs)-nghost .eq. domhi(2)+1) then
          j = domhi(2)
          do i = ARG_L1(rhs)+nghost, ARG_H1(rhs)-nghost
             dat(i,j+1) = dat(i,j)
          enddo
        endif
      endif

      if (extrap_corners .eq. 1) then

        if(ARG_L1(rhs)+nghost .eq. domlo(1) .and.
     &     ARG_L2(rhs)+nghost .eq. domlo(2)) then
          i = domlo(1)
          j = domlo(2)
          dat(i-1,j)   = dat(i,j)
          dat(i,j-1)   = dat(i,j)
          dat(i-1,j-1) = dat(i,j)
        endif

        if(ARG_H1(rhs)-nghost .eq. domhi(1)+1 .and.
     &     ARG_H2(rhs)+nghost .eq. domlo(2)) then
          i = domhi(1)
          j = domlo(2)
          dat(i+1,j)   = dat(i,j)
          dat(i,j-1)   = dat(i,j)
          dat(i+1,j-1) = dat(i,j)
        endif

        if(ARG_L1(rhs)+nghost .eq. domlo(1) .and.
     &     ARG_L2(rhs)-nghost .eq. domhi(2)+1) then
          i = domlo(1)
          j = domhi(2)
          dat(i-1,j)   = dat(i,j)
          dat(i,j+1)   = dat(i,j)
          dat(i-1,j+1) = dat(i,j)
        endif

        if(ARG_H1(rhs)-nghost .eq. domhi(1)+1 .and.
     &     ARG_H2(rhs)-nghost .eq. domhi(2)+1) then
          i = domhi(1)
          j = domhi(2)
          dat(i+1,j)   = dat(i,j)
          dat(i,j+1)   = dat(i,j)
          dat(i+1,j+1) = dat(i,j)
        endif
      endif

      if (ARG_L1(rhs)+1 .lt. ARG_L1(dat) .or. 
     $    ARG_H1(rhs)-1 .gt. ARG_H1(dat) .or.
     &    ARG_L2(rhs)+1 .lt. ARG_L2(dat) .or. 
     $    ARG_H2(rhs)-1 .gt. ARG_H2(dat)) then
        write(6,*) "FORT_HG_CELL_TO_NODE: bad index limits"
        stop
      endif

      if (is_rz.ne.1) then
        do j=ARG_L2(rhs)+nghost,ARG_H2(rhs)-nghost
          do i=ARG_L1(rhs)+nghost,ARG_H1(rhs)-nghost
            rhs(i,j) = fourth*(dat(i-1,j-1)+dat(i-1,j)+
     &                         dat(i,j-1)+dat(i,j))
          enddo
        enddo
      else
        do j=ARG_L2(rhs)+nghost,ARG_H2(rhs)-nghost
          do i=ARG_L1(rhs)+nghost,ARG_H1(rhs)-nghost
            if (hifix.eq.1 .and. ARG_H1(rhs)-nghost .eq. domhi(1)+1 .and.
     &          i .eq. ARG_H1(rhs)-nghost ) then
              rhi = rcen(i-1)
            else 
              rhi = rcen(i)
            endif
            if(lowfix.eq.1 .and. ARG_L1(rhs)+nghost .eq. domlo(1) .and.
     &          i .eq. ARG_L1(rhs)+nghost ) then
              rlo = rcen(i)
            else
              rlo = rcen(i-1)
            endif
            rhs(i,j) = fourth*(rlo * (dat(i-1,j-1) + dat(i-1,j)) +
     &                         rhi * (dat(i,j-1)   + dat(i,j)))
          enddo
        enddo
        factor = dr/24.0
        do j=ARG_L2(rhs)+nghost,ARG_H2(rhs)-nghost
          do i=ARG_L1(rhs)+nghost,ARG_H1(rhs)-nghost
            if (hifix.eq.1 .and. ARG_H1(rhs)-nghost .eq. domhi(1)+1 .and.
     &          i .eq. ARG_H1(rhs)-nghost ) then
              rhi = -one
            else 
              rhi = one
            endif
            if(lowfix.eq.1 .and. ARG_L1(rhs)+nghost .eq. domlo(1) .and.
     &          i .eq. ARG_L1(rhs)+nghost ) then
              rlo = -one
            else
              rlo = one
            endif
            rhs(i,j) = rhs(i,j) + factor *
     &                    (rlo * (dat(i-1,j-1) + dat(i-1,j)) -
     &                     rhi * (dat(i  ,j-1) + dat(i  ,j)))
          enddo
        enddo
      endif

c
c     end of routine
c
      return
      end



c=========================================================


      subroutine FORT_HGPHIBC(DIMS(u),u,DIMS(div),divcor,DIMS(rho),rho,
     &                        DIMS(r),r,hx,DIMS(phi),phi)
c
c ... compute the value of phi for hgproj to be used at an 
c     upper, y-outflow face
c
      implicit none

      integer DIMDEC(phi)
      integer DIMDEC(u)
      integer DIMDEC(div)
      integer DIMDEC(rho)
      integer DIMDEC(r)
      REAL_T    phi(DIMV(phi))
      REAL_T      u(DIMV(u),SDIM)
      REAL_T divcor(DIMV(div))
      REAL_T    rho(DIMV(rho))
      REAL_T      r(DIM1(r))
      REAL_T hx,rhs,fst,fstm,sigmax
      integer i,is,ie,js,je,iend,j
      integer ivp
      REAL_T  divcor_ave, vtot
      REAL_T  div_edge, rho_edge, div_edge_m1, rho_edge_m1, div_edge_m2, rho_edge_m2

      REAL_T max_divu, min_divu

#if 0
      return
#endif

      is = ARG_L1(phi)+1
      ie = ARG_H1(phi)-1
      js = ARG_L2(phi)+1
      je = ARG_H2(phi)-1
      iend =ie

#if 0
        divcor_ave = zero
        vtot = zero
        do i=is,ie-1 
#if 0
           div_edge = half*(divcor(i,je-1)+divcor(i,je))
#endif
#if 1
           div_edge = divcor(i,je)
#endif
           divcor_ave = divcor_ave+r(i)*div_edge
           vtot = vtot+r(i)
        enddo
        divcor_ave = divcor_ave/vtot

        fst = zero
#if 0
        div_edge = half*(divcor(is,je-1)+divcor(is,je))
        rho_edge = half*(rho(is,je-1)+rho(is,je))
#endif
#if 1
        div_edge = divcor(is,je)
        rho_edge = rho(is,je)
#endif
        fstm = (div_edge-divcor_ave)*hx
        phi(is,je)   = zero
        phi(is+1,je) = rho_edge*hx*(u(is,je-1,1) -
     $      half*(div_edge-divcor_ave)*hx)

        do i=is+2,ie 
c          write(6,*)'rho,u,divcor=',rho(i-1,je-1),u(i-1,je-1,1),divcor(i-1,je-1)
#if 0
          div_edge_m1 = half*(divcor(i-1,je-1)+divcor(i-1,je))
          div_edge_m2 = half*(divcor(i-2,je-1)+divcor(i-2,je))
          rho_edge_m1 = half*(rho(i-1,je-1)+rho(i-1,je))
          rho_edge_m2 = half*(rho(i-2,je-1)+rho(i-2,je))
#endif
#if 1
          div_edge_m1 = divcor(i-1,je)
          div_edge_m2 = divcor(i-2,je)
          rho_edge_m1 = rho(i-1,je)
          rho_edge_m2 = rho(i-2,je)
#endif
          fst = (r(i-2)*fstm + half*(r(i-1)*(div_edge_m1-divcor_ave) +
     $         r(i-2)*(div_edge_m2-divcor_ave))*hx)/r(i-1)
          rhs = r(i-2)*(phi(i-1,je)-phi(i-2,je))/(rho_edge_m2*hx) +
     $          r(i-1)*(u(i-1,je-1,1) - fst) - 
     &          r(i-2)*(u(i-2,je-1,1) - fstm)
          phi(i,je) = phi(i-1,je) + hx*rho_edge_m1*rhs / r(i-1)
          fstm = fst
        enddo
#endif
#if 1
        divcor_ave = zero
        vtot = zero
        max_divu = zero
        min_divu = zero
        do j=je-1,je-1
        do i=is,ie-1 
           divcor_ave = divcor_ave+r(i)*divcor(i,j)
           vtot = vtot+r(i)
           max_divu = max(max_divu,divcor(i,j))
           min_divu = min(min_divu,divcor(i,j))
        enddo
        enddo
        divcor_ave = divcor_ave/vtot

        if(max_divu.eq.zero.and.min_divu.eq.zero) then
          return
        endif

c        write(6,*)"divcor_ave",divcor_ave
c        call flush(6)


        fst = zero
        fstm = (divcor(is,je-1)-divcor_ave)*hx
        phi(is,je)   = zero
c        write(6,*)'rho,u,divcor=',rho(is,je-1),u(is,je-1,1),divcor(is,je-1)
        phi(is+1,je) = rho(is,je-1)*hx*(u(is,je-1,1) -
     $      half*(divcor(is,je-1)-divcor_ave)*hx)

        do i=is+2,ie 
          fst = (r(i-2)*fstm + half*(r(i-1)*(divcor(i-1,je-1)-divcor_ave) +
     $         r(i-2)*(divcor(i-2,je-1)-divcor_ave))*hx)/r(i-1)
          rhs = r(i-2)*(phi(i-1,je)-phi(i-2,je))/(rho(i-2,je-1)*hx) +
     $          r(i-1)*(u(i-1,je-1,1) - fst) - 
     &          r(i-2)*(u(i-2,je-1,1) - fstm)
          phi(i,je) = phi(i-1,je) + hx*rho(i-1,je-1)*rhs / r(i-1)
          fstm = fst
        enddo
#endif
c
c     end of routine
c
      return
      end


