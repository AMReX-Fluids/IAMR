c
c $Id: MACPROJ_2D.F,v 1.6 1998-11-05 18:44:25 lijewski Exp $
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
     &               +   yarea(i,j+1)*uy(i,j+1) - yarea(i,j)*uy(i,j)
     &               )/vol(i,j)

             end do
          end do
       else
          do i = lo(1), hi(1)
             do j = lo(2), hi(2)
                dmac(i,j) = ( xarea(i+1,j)*ux(i+1,j) - xarea(i,j)*ux(i,j)
     &               +   yarea(i,j+1)*uy(i,j+1) - yarea(i,j)*uy(i,j)
     &               )/vol(i,j)
             end do
          end do
       end if

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
                end do
             end do
          end do
       else
          do n = 1, nvar
             do i = lo(1), hi(1)
                do j = lo(2), hi(2)
                   snew(i,j,n) = snew(i,j,n) + mult*dmac(i,j)*sold(i,j,n)
                end do
             end do
          end do
       end if

       end

      subroutine FORT_MACPHIBC(DIMS(u),u,DIMS(div),divcor,DIMS(rho),rho,
     &                         DIMS(r),r,hx,DIMS(phi),phi,per)
c
c     Compute the value of phi for macproj to be used at an  upper, y-outflow face,
c     assuming that the tangential velocity on the edges of the outflow boundary
c     are either zero or periodic.
c
      implicit none

      integer DIMDEC(phi)
      integer DIMDEC(u)
      integer DIMDEC(div)
      integer DIMDEC(rho)
      integer DIMDEC(r)
      integer per
      REAL_T    phi(DIMV(phi))
      REAL_T      u(DIMV(u))
      REAL_T divcor(DIMV(div))
      REAL_T    rho(DIMV(rho))
      REAL_T      r(DIM1(r))
      REAL_T hx
      
c     Local variables
      integer NstripMAX
      parameter (NstripMAX = 500)
      REAL_T small_udiff, small_pert
      parameter (small_udiff = 1.d-10)
      parameter ( small_pert = 1.d-10)
      integer i, ii, j, neq
      REAL_T a(NstripMAX), b(NstripMAX), c(NstripMAX), s(NstripMAX)
      REAL_T x(NstripMAX), alpha, beta, sVal
      logical rNormed
      REAL_T divcor_ave, vtot
      REAL_T rL, rR, diff
      REAL_T max_divu, min_divu, max_pert

c     NOTE: Assumes that rho at edge between i, i-1 = half*(rho(i)+rho(i-1))
c             (1) Linear fit of rho between nodes
c             (2) rho, divcor on same boxes (box), and r(ARG_L1(rho),ARG_L2(rho))
c             (3) phi is on box, shifted up one
c             (4) u is edge-based, on surroundingNodes(box)
#ifndef NDEBUG
      if (ARG_H1(phi) - ARG_L1(phi) + 1 .gt. NstripMAX) then
         
         write(6,*) 'FORT_MACPHIBC: recompile with NstripMAX >',
     &        ARG_H1(phi) - ARG_L1(phi) + 1
         stop
         
      end if
      
c     Check assumptions on input data
      if ( (ARG_L1(rho).ne.ARG_L1(div))   .or.
     &     (ARG_L1(rho).ne.ARG_L1(u))     .or.
     &     (ARG_L1(rho).ne.ARG_L1(r))     .or.
     &     (ARG_L1(rho).ne.ARG_L1(phi))   .or.
     &     (ARG_H1(rho).ne.ARG_H1(div))   .or.
     &     (ARG_H1(rho).ne.ARG_H1(u)-1)   .or.
     &     (ARG_H1(rho).ne.ARG_H1(r))     .or.
     &     (ARG_H1(rho).ne.ARG_H1(phi))   .or.
     &     (ARG_L2(rho).ne.ARG_L2(div))   .or.
     &     (ARG_L2(rho).ne.ARG_L2(phi)-1) .or.
     &     (ARG_L2(rho).ne.ARG_L2(u))     .or.
     &     (ARG_H2(rho).ne.ARG_H2(div))   .or.
     &     (ARG_H2(rho).ne.ARG_H2(u))     .or.
     &     (ARG_H2(rho).ne.ARG_H2(phi)-1)  ) then
         write(6,*) 'FORT_MACPHIBC: bad size for input data'
         stop
      endif
      
#endif

c     Solve d/dx( 1/rho d/dx( phi ) ) = dU/dx - (S - S_ave) [S = divu] w/periodic or
c     Neumann BC's, using a tridiagonal solve which detects, and deals with,
c     the singular equations.  In the Neumann case, arbitrarily set the upper right
c     corner to zero to pin the solution.  Note that the RHS of this equation
c     satisfies the solvability constraint that Int[RHS.dV] = 0 by construction.
c     This implies that the normal component takes up the slack:
c     
c                        d/dy( 1/rho d/dy( phi ) ) = dV/dy - S_ave
c     
c     This information should be used to construct the normal gradient of the
c     normal velocity, for the advective/diffusive step, for example.
c     In this implementation, use that d/dy == 0 at top, so y-edge centered values
c     come directly from cell-centers just inside domain
c      
c     Compute average of divu over outflow bc.  Set trivial solution if average
c     is zero, or if divu is constant
      divcor_ave = zero
      vtot = zero
      j = ARG_H2(rho)
      max_divu = divcor(ARG_L1(rho),j)
      min_divu = max_divu
      do i = ARG_L1(rho), ARG_H1(rho)
         divcor_ave = divcor_ave+r(i)*divcor(i,j)
         vtot = vtot+r(i)
         max_divu = max(max_divu,divcor(i,j))
         min_divu = min(min_divu,divcor(i,j))
      end do
      divcor_ave = divcor_ave/vtot
      max_pert = ABS(divcor(ARG_L1(rho),j)-divcor_ave)
      do i = ARG_L1(rho)+1, ARG_H1(rho)
         max_pert = MAX(max_pert,ABS(divcor(i,j)-divcor_ave))
      end do
      max_pert = max_pert/(ABS(divcor_ave)+small_pert)
      if ((max_divu.eq.zero.and.min_divu.eq.zero)
     &     .or.(max_pert.le.small_pert)) then
         do j = ARG_L2(phi),ARG_H2(phi)
            do i = ARG_L1(phi),ARG_H1(phi)
               phi(i,j) = zero
            end do
         end do
         return
      end if

c     Carry out non-trivial solve.  First set interior equations, then do BC's
      neq = ARG_H1(rho) - ARG_L1(rho) + 1
      do i = 2,neq-1
         ii = i + ARG_L1(rho) - 1
         rL = half*(r(ii)+r(ii-1))
         rR = half*(r(ii)+r(ii+1))
         a(i)=two*rL/(rho(ii,j)+rho(ii-1,j))
         c(i)=two*rR/(rho(ii,j)+rho(ii+1,j))
         b(i)=- a(i) - c(i)
         s(i)=hx*(rR*u(ii+1,j)-rL*u(ii,j)-hx*(divcor(ii,j)-divcor_ave))
      end do

      if (per .eq. 1) then

c     Make sure u_mac is periodic
         diff = u(ARG_L1(u),j) - u(ARG_H1(u),j)
         if (ABS(diff) .gt. small_udiff) then
            write(6,*) 'MACPHIBC: periodic umac screwy.  diff =',diff
            stop
         endif
         
c     Do left-side periodic BC (keep r in there to guarantee correct scaling)
         ii = ARG_L1(rho)
         rL = half*(r(ii)+r(ARG_H1(rho)))
         rR = half*(r(ii)+r(ii+1))
         beta=two*rL/(rho(ii,j)+rho(ARG_H1(rho),j))
         c(1)=two*rR/(rho(ii,j)+rho(ii+1,j))
         b(1)=- beta - c(1)
         s(1)=hx*(rR*u(ii+1,j)-rL*u(ii,j)-hx*(divcor(ii,j)-divcor_ave))
         
c     Do right-side periodic
         ii = ARG_H1(rho)
         rL = half*(r(ii)+r(ii-1))
         rR = half*(r(ARG_L1(rho))+r(ii))
         a(neq)=two*rL/(rho(ii,j)+rho(ii-1,j))
         alpha=two*rR/(rho(ii,j)+rho(ARG_L1(rho),j))
         b(neq)=- a(neq) - alpha
         s(neq)=hx*(rR*u(ii+1,j)-rL*u(ii,j)-hx*(divcor(ii,j)-divcor_ave))
         
c     Solve the equations
         call cyclic(a,b,c,alpha,beta,s,x,neq)
         
      else

c     Solid walls, Nuemann conditions (dphi/dx=u=0)
         ii = ARG_L1(rho)
         rR = half*(r(ii)+r(ii+1))
         c(1) = two*rR/(rho(ii,j)+rho(ii+1,j))
         b(1) = - c(1)
         s(1) = hx*(rR*u(ii+1,j)-hx*(divcor(ii,j)-divcor_ave))

         ii = ARG_H1(rho)
         rL = half*(r(ii)+r(ii-1))
         a(neq) = two*rL/(rho(ii,j)+rho(ii-1,j))
         b(neq) = - a(neq)
         s(neq) = hx*(-rL*u(ii,j)-hx*(divcor(ii,j)-divcor_ave))
c     Solve the equations (we know they're singular, pass the arbitrary value, and
c     a flag that we've already normalized the rhs, in the sense that
c                          Int[dU/dx - (S-S_ave)] == 0
         sVal = zero
         rNormed = .true.
         call tridag_sing(a,b,c,s,x,neq,sVal,rNormed)
         
      end if
      
c     Write the solution into the arg
      j = ARG_H2(phi)
      do i = 1,neq
         phi(ARG_L1(phi)-1+i,j) = x(i)
      end do
      
      end
