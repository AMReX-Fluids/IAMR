c
c $Id: MACPROJ_2D.F,v 1.12 1999-07-12 23:16:03 propp Exp $
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

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define SMALL 1.0e-10
#else
#define SMALL 1.0d-10
#endif

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
     &                         r_lo,r_hi,r,hx,DIMS(phi),phi,face,per)
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
      integer r_lo,r_hi
      integer per
      integer face
      REAL_T    phi(DIMV(phi))
      REAL_T      u(DIMV(u))
      REAL_T divcor(DIMV(div))
      REAL_T    rho(DIMV(rho))
      REAL_T      r(r_lo:r_hi)
      REAL_T hx
      
c     Local variables
      integer NstripMAX
      parameter (NstripMAX = 1000)
      REAL_T small_udiff, small_pert
      parameter (small_udiff = SMALL)
      parameter ( small_pert = SMALL)
      integer i, ii, j, jj, neq,n
      REAL_T a(NstripMAX), b(NstripMAX), c(NstripMAX), s(NstripMAX)
      REAL_T x(NstripMAX), alpha, beta, sVal
      logical rNormed
      REAL_T divcor_ave, vtot
      REAL_T rL, rR, diff
      REAL_T max_divu, min_divu, max_pert, phitot
      integer rstart
      integer rho_start,rho_end
#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3
c     NOTE: Assumes that rho at edge between i, i-1 = half*(rho(i)+rho(i-1))
c             (1) Linear fit of rho between nodes
c             (2) rho, divcor on same boxes (box), and r(ARG_L1(rho),ARG_L2(rho))
c             (3) phi is on box, shifted up one
c             (4) u is edge-based, on surroundingNodes(box)
#ifndef NDEBUG
#if 0
      if (ARG_H1(phi) - ARG_L1(phi) + 1 .gt. NstripMAX) then
         
         write(6,*) 'FORT_MACPHIBC: recompile with NstripMAX >',
     &        ARG_H1(phi) - ARG_L1(phi) + 1
         stop
         
      end if
      
c     Check assumptions on input data
      if ( (ARG_L1(rho).ne.ARG_L1(div))   .or.
     &     (ARG_L1(rho).ne.ARG_L1(u))     .or.
     &     (ARG_L1(rho).ne.r_lo)     .or.
     &     (ARG_L1(rho).ne.ARG_L1(phi))   .or.
     &     (ARG_H1(rho).ne.ARG_H1(div))   .or.
     &     (ARG_H1(rho).ne.ARG_H1(u)-1)   .or.
     &     (ARG_H1(rho).ne.r_hi)     .or.
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

      if (face .eq. XLO .or. face .eq. XHI) then
         i = ARG_H1(rho)
         max_divu = divcor(i,ARG_L2(rho))
         min_divu = max_divu
         do j = ARG_L2(rho), ARG_H2(rho)
            divcor_ave = divcor_ave+r(j)*divcor(i,j)
            vtot = vtot+r(j)
            max_divu = max(max_divu,divcor(i,j))
            min_divu = min(min_divu,divcor(i,j))
         end do
         divcor_ave = divcor_ave/vtot
         max_pert = ABS(divcor(i,ARG_L2(rho))-divcor_ave)
         do j = ARG_L2(rho)+1, ARG_H2(rho)
            max_pert = MAX(max_pert,ABS(divcor(i,j)-divcor_ave))
         end do
         max_pert = max_pert/(ABS(divcor_ave)+small_pert)
      else if (face .eq. YLO .or. face .eq. YHI) then
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
      endif
      
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
      if (face .eq. XLO .or. face .eq. XHI) then
         rho_start = ARG_L2(rho)
         rho_end   = ARG_H2(rho)
         neq = rho_end - rho_start + 1
         i = ARG_H1(rho)
         do n = 2,neq-1
            j = n + rho_start - 1
            rL = half*(r(j)+r(j-1))
            rR = half*(r(j)+r(j+1))
            a(n)=two*rL/(rho(i,j)+rho(i,j-1))
            c(n)=two*rR/(rho(i,j)+rho(i,j+1))
            b(n)=- a(n) - c(n)
            s(n)=(rR*u(i,j+1)-rL*u(i,j))*hx - r(j)*(divcor(i,j)-divcor_ave)*hx*hx
         end do
         
         if (per .eq. 1) then
            
c     Make sure u_mac is periodic
            diff = u(i,ARG_L2(u)) - u(i,ARG_H2(u))
            if (ABS(diff) .gt. small_udiff) then
               write(6,*) 'MACPHIBC: periodic umac screwy.  diff =',diff
               stop
            endif
            
c     Do left-side periodic BC (keep r in there to guarantee correct scaling)
            j = rho_start
            rL = half*(r(j)+r(rho_end))
            rR = half*(r(j)+r(j+1))
            beta=two*rL/(rho(i,j)+rho(i,rho_end))
            c(1)=two*rR/(rho(i,j)+rho(i,j+1))
            b(1)=- beta - c(1)
            s(1)=(rR*u(i,j+1)-rL*u(i,j))*hx - r(j)*(divcor(i,j)-divcor_ave)*hx*hx
            
c     Do right-side periodic
            j = rho_end
            rL = half*(r(j)+r(j-1))
            rR = half*(r(rho_start)+r(j))
            a(neq)=two*rL/(rho(i,j)+rho(i,j-1))
            alpha=two*rR/(rho(i,j)+rho(i,rho_start))
            b(neq)=- a(neq) - alpha
            s(neq)=(rR*u(i,j+1)-rL*u(i,j))*hx - r(j)*(divcor(i,j)-divcor_ave)*hx*hx
            
c     Solve the equations
            call cyclic(a,b,c,alpha,beta,s,x,neq)
            
         else
            
c     Solid walls, Nuemann conditions (dphi/dx=u=0)
            j = rho_start
            rR = half*(r(j)+r(j+1))
            c(1) = two*rR/(rho(i,j)+rho(i,j+1))
            b(1) = - c(1)
            s(1)= rR*u(i,j+1)*hx - r(j)*(divcor(i,j)-divcor_ave)*hx*hx
            
            
            j = rho_end
            rL = half*(r(j)+r(j-1))
            a(neq) = two*rL/(rho(i,j)+rho(i,j-1))
            b(neq) = - a(neq)
            s(neq)=-rL*u(i,j)*hx - r(j)*(divcor(i,j)-divcor_ave)*hx*hx
            
c     Solve the equations (we know they're singular, pass the arbitrary value, and
c     a flag that we've already normalized the rhs, in the sense that
c     Int[dU/dx - (S-S_ave)] == 0
            sVal = zero
            rNormed = .true.
            call tridag_sing(a,b,c,s,x,neq,sVal,rNormed)
            
         end if
      else if (face .eq. YLO .or. face .eq. YHI) then
         rho_start = ARG_L1(rho)
         rho_end   = ARG_H1(rho)
         neq = rho_end - rho_start + 1
         j = ARG_H2(rho)
         do n = 2,neq-1
            i = n + rho_start - 1
            rL = half*(r(i)+r(i-1))
            rR = half*(r(i)+r(i+1))
            a(n)=two*rL/(rho(i,j)+rho(i-1,j))
            c(n)=two*rR/(rho(i,j)+rho(i+1,j))
            b(n)=- a(n) - c(n)
            s(n)=(rR*u(i+1,j)-rL*u(i,j))*hx - r(i)*(divcor(i,j)-divcor_ave)*hx*hx
         end do
         
         if (per .eq. 1) then
            
c     Make sure u_mac is periodic
            diff = u(ARG_L1(u),j) - u(ARG_H1(u),j)
            if (ABS(diff) .gt. small_udiff) then
               write(6,*) 'MACPHIBC: periodic umac screwy.  diff =',diff
               stop
            endif
            
c     Do left-side periodic BC (keep r in there to guarantee correct scaling)
            i = rho_start
            rL = half*(r(i)+r(rho_end))
            rR = half*(r(i)+r(i+1))
            beta=two*rL/(rho(i,j)+rho(rho_end,j))
            c(1)=two*rR/(rho(i,j)+rho(i+1,j))
            b(1)=- beta - c(1)
            s(1)=(rR*u(i+1,j)-rL*u(i,j))*hx - r(i)*(divcor(i,j)-divcor_ave)*hx*hx
            
c     Do right-side periodic
            i = rho_end
            rL = half*(r(i)+r(i-1))
            rR = half*(r(rho_start)+r(i))
            a(neq)=two*rL/(rho(i,j)+rho(i-1,j))
            alpha=two*rR/(rho(i,j)+rho(rho_start,j))
            b(neq)=- a(neq) - alpha
            s(neq)=(rR*u(i+1,j)-rL*u(i,j))*hx - r(i)*(divcor(i,j)-divcor_ave)*hx*hx
            
c     Solve the equations
            call cyclic(a,b,c,alpha,beta,s,x,neq)
            
         else
            
c     Solid walls, Nuemann conditions (dphi/dx=u=0)
            i = rho_start
            rR = half*(r(i)+r(i+1))
            c(1) = two*rR/(rho(i,j)+rho(i+1,j))
            b(1) = - c(1)
            s(1)= rR*u(i+1,j)*hx - r(i)*(divcor(i,j)-divcor_ave)*hx*hx
            
            
            i = rho_end
            rL = half*(r(i)+r(i-1))
            a(neq) = two*rL/(rho(i,j)+rho(i-1,j))
            b(neq) = - a(neq)
            s(neq)=-rL*u(i,j)*hx - r(i)*(divcor(i,j)-divcor_ave)*hx*hx
            
c     Solve the equations (we know they're singular, pass the arbitrary value, and
c     a flag that we've already normalized the rhs, in the sense that
c     Int[dU/dx - (S-S_ave)] == 0
            sVal = zero
            rNormed = .true.
            call tridag_sing(a,b,c,s,x,neq,sVal,rNormed)
            
         end if
      end if

c     Try normalizing phi to average to zero
      phitot = zero
      vtot = zero
      if (face .eq. XLO .or. face .eq. XHI) then
         rstart = ARG_L2(rho)
      else if (face .eq. YLO .or. face .eq. YHI) then
         rstart = ARG_L1(rho)
      endif
      do n = 1,neq
         phitot = phitot + x(n)*r(rstart+n-1)
         vtot = vtot + r(rstart+n-1)
      end do
      phitot = phitot / vtot
      do n = 1,neq
         x(n) = x(n) - phitot
      end do

c     Write the solution into the arg
      if (face .eq. XLO .or. face .eq. XHI) then
         do n = 1,neq
            phi(ARG_H1(phi),ARG_L2(phi)-1+n) = x(n)
         end do
      else if (face .eq. YLO .or. face .eq. YHI) then
         do n = 1,neq
            phi(ARG_L1(phi)-1+n,ARG_H2(phi)) = x(n)
         end do
      endif
#undef XLO
#undef YLO
#undef XHI
#undef YHI
      end
