c
c $Id: PROJECTION_2D.F,v 1.17 1999-02-26 22:22:06 almgren Exp $
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

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define SMALL 1.0e-10
#else
#define SMALL 1.0d-10
#endif

       subroutine FORT_FILTERP (
     &     pnew,scratch,DIMS(p),
     &     lo,hi,dx,fac,
     &     wrap_around_x,wrap_around_y)
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
c     wrap_around_x => is this a single periodic box in the x-direction
c     wrap_around_y => is this a single periodic box in the y-direction
c     
c     ----------------------------------------------------------
c 
       integer    DIMDEC(p)
       REAL_T     pnew(DIMV(p))
       REAL_T     scratch(DIMV(p))
       integer    lo(2), hi(2)
       REAL_T     dx(2)
       REAL_T     fac
       integer    wrap_around_x
       integer    wrap_around_y
       integer i, j

       if ( abs(dx(1)-dx(2)) .gt. 1.e-8 .and. fac .gt. 1.e-4 )then
          print *, 'is filter properly formatted for dx != dy ?'
       end if

       do j = lo(2)+1, hi(2)
          do i = lo(1)+1, hi(1)
             scratch(i,j) = four * pnew(i,j) - 
     $         (pnew(i-1,j) + pnew(i+1,j) + pnew(i,j-1) + pnew(i,j+1) )
             scratch(i,j) = - eighth * scratch(i,j) * fac
          end do
       end do
 
       if (wrap_around_x .eq. 1) then
         i = lo(1)
         do j = lo(2)+1, hi(2)
             scratch(i,j) = four * pnew(i,j) - 
     $         (pnew(hi(1),j) + pnew(i+1,j) + pnew(i,j-1) + pnew(i,j+1) )
             scratch(i,j) = - eighth * scratch(i,j) * fac
         end do

         i = hi(1)+1
         do j = lo(2)+1, hi(2)
             scratch(i,j) = four * pnew(i,j) - 
     $         (pnew(i-1,j) + pnew(lo(1)+1,j) + pnew(i,j-1) + pnew(i,j+1) )
             scratch(i,j) = - eighth * scratch(i,j) * fac
         end do
       end if
 
       if (wrap_around_y .eq. 1) then
         j = lo(2)
         do i = lo(1)+1, hi(1)
             scratch(i,j) = four * pnew(i,j) - 
     $         (pnew(i-1,j) + pnew(i+1,j) + pnew(i,hi(2)) + pnew(i,j+1) )
             scratch(i,j) = - eighth * scratch(i,j) * fac
         end do

         j = hi(2)+1
         do i = lo(1)+1, hi(1)
             scratch(i,j) = four * pnew(i,j) - 
     $         (pnew(i-1,j) + pnew(i+1,j) + pnew(i,j-1) + pnew(i,lo(2)+1) )
             scratch(i,j) = - eighth * scratch(i,j) * fac
         end do
       end if
 
       if (wrap_around_x .eq. 1 .and. wrap_around_y .eq. 1) then
         i = lo(1)
         j = lo(2)
           scratch(i,j) = four * pnew(i,j) - 
     $       (pnew(hi(1),j) + pnew(i+1,j) + pnew(i,hi(2)) + pnew(i,j+1) )
           scratch(i,j) = - eighth * scratch(i,j) * fac

         i = hi(1)+1
         j = lo(2)
           scratch(i,j) = four * pnew(i,j) - 
     $       (pnew(i-1,j) + pnew(lo(1)+1,j) + pnew(i,hi(2)) + pnew(i,j+1) )
           scratch(i,j) = - eighth * scratch(i,j) * fac

         i = lo(1)
         j = hi(2)+1
           scratch(i,j) = four * pnew(i,j) - 
     $       (pnew(hi(1),j) + pnew(i+1,j) + pnew(i,j-1) + pnew(i,lo(2)+1) )
           scratch(i,j) = - eighth * scratch(i,j) * fac

         i = hi(1)+1
         j = hi(2)+1
           scratch(i,j) = four * pnew(i,j) - 
     $       (pnew(i-1,j) + pnew(lo(1)+1,j) + pnew(i,j-1) + pnew(i,lo(2)+1) )
           scratch(i,j) = - eighth * scratch(i,j) * fac
       end if

       do j = lo(2)+1, hi(2)
          do i = lo(1)+1, hi(1)
             pnew(i,j) = pnew(i,j) + scratch(i,j)
          end do
       end do

       if (wrap_around_x .eq. 1) then
         do j = lo(2)+1, hi(2)
           pnew(lo(1)  ,j) = pnew(lo(1)  ,j) + scratch(lo(1)  ,j)
           pnew(hi(1)+1,j) = pnew(hi(1)+1,j) + scratch(hi(1)+1,j)
         end do
       end if

       if (wrap_around_y .eq. 1) then
         do i = lo(1)+1, hi(1)
           pnew(i,lo(2)  ) = pnew(i,lo(2)  ) + scratch(i,lo(2)  )
           pnew(i,hi(2)+1) = pnew(i,hi(2)+1) + scratch(i,hi(2)+1)
         end do
       end if

       if (wrap_around_x .eq. 1 .and. wrap_around_y .eq. 1) then
         pnew(lo(1)  ,lo(2)  ) = pnew(lo(1)  ,lo(2)  ) + scratch(lo(1)  ,lo(2)  )
         pnew(hi(1)+1,lo(2)  ) = pnew(hi(1)+1,lo(2)  ) + scratch(hi(1)+1,lo(2)  )
         pnew(lo(1)  ,hi(2)+1) = pnew(lo(1)  ,hi(2)+1) + scratch(lo(1)  ,hi(2)+1)
         pnew(hi(1)+1,hi(2)+1) = pnew(hi(1)+1,hi(2)+1) + scratch(hi(1)+1,hi(2)+1)
       end if

       end

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

       integer i, j, n

       do n = 1, SDIM
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1
                unew(i,j,n) = uold(i,j,n) + dt*unew(i,j,n)
             end do
          end do
       end do

       end

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

      integer i, j, n

      do n = 1, SDIM
         do j = lo(2)-1, hi(2)+1
            do i = lo(1)-1, hi(1)+1
               unew(i,j,n) = (unew(i,j,n)-uold(i,j,n))/dt
            end do
         end do
      end do

      end

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

      integer i, j, n

      do n = 1, nvar
         do j = boxlo(2), boxhi(2)
            do i = boxlo(1), boxhi(1)
               un(i,j,n) = un(i,j,n) + alpha*uo(i,j,n)
            end do
         end do
      end do

      end

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

       do j = ARG_L2(grid)-n(2), ARG_H2(grid)+n(2)
          do i = ARG_L1(grid)-n(1), ARG_H1(grid)+n(1)
	     a(i,j) = r(i)*a(i,j)
          end do
       end do

       end

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

       integer i, j

       do j = ARG_L2(grid)-n(2), ARG_H2(grid)+n(2)
          do i = ARG_L1(grid)-n(1), ARG_H1(grid)+n(1)
	    a(i,j) = a(i,j)/r(i)
          end do
       end do

       end

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
      end if

      if(isrz.eq.1.and.lrweighted.ne.1)then
        write(6,*)
     &    'FORT_HGN2C : isrz=1 and lrweighted!=1 not yet implemented'
#ifndef	BL_NO_FORT_FLUSH
        call flush(6)
#endif
        stop
      end if

      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          ccdat(i,j) = fourth*(nodedat(i,j)+nodedat(i+1,j)+
     &                         nodedat(i,j+1)+nodedat(i+1,j+1))
        end do
      end do

      end

      subroutine FORT_HGC2N(
     &     nghost, DIMS(dat), dat, rcen,
     &     DIMS(rhs), rhs,
     &     domlo, domhi, lowfix, hifix, dr, is_rz) 
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
      REAL_T  dr
      REAL_T  rcen(DIM1(dat))
      REAL_T  dat(DIMV(dat))
      REAL_T  rhs(DIMV(rhs))

      integer i, j
      REAL_T  factor, rhi, rlo

      if (ARG_L1(rhs)+1 .lt. ARG_L1(dat) .or. 
     $     ARG_H1(rhs)-1 .gt. ARG_H1(dat) .or.
     &     ARG_L2(rhs)+1 .lt. ARG_L2(dat) .or. 
     $     ARG_H2(rhs)-1 .gt. ARG_H2(dat)) then
         write(6,*) "FORT_HG_CELL_TO_NODE: bad index limits"
         stop
      end if

      if (is_rz.ne.1) then
         do j=ARG_L2(rhs)+nghost,ARG_H2(rhs)-nghost
            do i=ARG_L1(rhs)+nghost,ARG_H1(rhs)-nghost
               rhs(i,j) = fourth*(dat(i-1,j-1)+dat(i-1,j)+
     &              dat(i,j-1)+dat(i,j))
            end do
         end do
      else
         do j=ARG_L2(rhs)+nghost,ARG_H2(rhs)-nghost
            do i=ARG_L1(rhs)+nghost,ARG_H1(rhs)-nghost
               if (hifix.eq.1 .and. ARG_H1(rhs)-nghost .eq. domhi(1)+1 .and.
     &              i .eq. ARG_H1(rhs)-nghost ) then
                  rhi = rcen(i-1)
               else 
                  rhi = rcen(i)
               end if
               if(lowfix.eq.1 .and. ARG_L1(rhs)+nghost .eq. domlo(1) .and.
     &              i .eq. ARG_L1(rhs)+nghost ) then
                  rlo = rcen(i)
               else
                  rlo = rcen(i-1)
               end if
               rhs(i,j) = fourth*(rlo * (dat(i-1,j-1) + dat(i-1,j)) +
     &              rhi * (dat(i,j-1)   + dat(i,j)))
            end do
         end do
         factor = dr/24.0
         do j=ARG_L2(rhs)+nghost,ARG_H2(rhs)-nghost
            do i=ARG_L1(rhs)+nghost,ARG_H1(rhs)-nghost
               if (hifix.eq.1 .and. ARG_H1(rhs)-nghost .eq. domhi(1)+1 .and.
     &              i .eq. ARG_H1(rhs)-nghost ) then
                  rhi = -one
               else 
                  rhi = one
               end if
               if(lowfix.eq.1 .and. ARG_L1(rhs)+nghost .eq. domlo(1) .and.
     &              i .eq. ARG_L1(rhs)+nghost ) then
                  rlo = -one
               else
                  rlo = one
               end if
               rhs(i,j) = rhs(i,j) + factor *
     &              (rlo * (dat(i-1,j-1) + dat(i-1,j)) -
     &              rhi * (dat(i  ,j-1) + dat(i  ,j)))
            end do
         end do
      end if

      end

      subroutine FORT_HGPHIBC(DIMS(u),u,DIMS(div),divcor,DIMS(rho),rho,
     &                        DIMS(r),r,hx,DIMS(phi),phi,per)
c
c     Compute the value of phi for hgproj to be used at an  upper, y-outflow face,
c     assuming that the tangential velocity on the edges of the outflow boundary
c     are either zero or periodic.  Algorithm incorrect for stress-free side walls,
c     or if there is outflow on any other boundary.
c
      implicit none

      integer DIMDEC(phi)
      integer DIMDEC(u)
      integer DIMDEC(div)
      integer DIMDEC(rho)
      integer DIMDEC(r)
      integer per
      REAL_T    phi(DIMV(phi))
      REAL_T      u(DIMV(u),SDIM)
      REAL_T divcor(DIMV(div))
      REAL_T    rho(DIMV(rho))
      REAL_T      r(DIM1(r))
      REAL_T hx
      
c     Local variables
      integer NstripMAX
      parameter (NstripMAX = 500)
      integer ics, ice, jcs, jce, ins, ine, jns, jne
      integer ni, nj, i, j, icL, icR, neq
      REAL_T a(NstripMAX), b(NstripMAX), c(NstripMAX), s(NstripMAX)
      REAL_T x(NstripMAX), alpha, beta, sVal
      REAL_T max_pert, small_pert
      parameter ( small_pert = SMALL)
      logical rNormed
      REAL_T divcor_ave, vtot
      REAL_T div_edge, rho_edge
      REAL_T div_edge_m1, rho_edge_m1, div_edge_m2, rho_edge_m2
      REAL_T max_divu, min_divu
      REAL_T uExt(0:NstripMAX-1)

      ins = ARG_L1(phi)
      ine = ARG_H1(phi)
      jns = ARG_L2(phi)
      jne = ARG_H2(phi)
      
      ics = ARG_L1(u)
      ice = ARG_H1(u)
      jcs = ARG_L2(u)
      jce = ARG_H2(u)
      
c     Get average and extremes of divu.  If == 0, then boundary nodes should be set to
c     phi = 0.  Assume that cell-centered values are given on a 3-wide strip across the
c     entire domain--the top row of the box are cell-centered "ghost" cells.  The phi
c     values are given along a 1-wide strip along the top of the domain.
#ifndef NDEBUG
      if ( (ins .ne. ics)   .or.
     &     (ine .ne. ice+1) .or.
     &     (jns .ne. jcs+2) .or.
     &     (jne .ne. jce)   .or.
     &     (jns .ne. jne) ) then
         
         write(6,*) 'FORT_HGPHIBC: bad bounds'
         stop
         
      else if (ine - ins + 1 .gt. NstripMAX) then
         
         write(6,*) 'FORT_HGPHIBC: recompile with NstripMAX >',ine-ins+1
         stop
         
      end if
#endif
      
c     In this implementation, use that d/dy == 0 at top, so y-edge centered values come
c     directly from cell-centers just inside domain
      
      divcor_ave = zero
      vtot = zero
      j = jcs + 1
      max_divu = divcor(ics,j)
      min_divu = max_divu
      do i = ics, ice
         divcor_ave = divcor_ave+r(i)*divcor(i,j)
         vtot = vtot+r(i)
         max_divu = max(max_divu,divcor(i,j))
         min_divu = min(min_divu,divcor(i,j))
         uExt(i-ics) = half*(three*u(i,jcs+1,1) - u(i,jcs,1))
      end do
      divcor_ave = divcor_ave/vtot
      max_pert = ABS(divcor(ics,j)-divcor_ave)
      do i = ics+1, ice
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
      
c     Solve d/dx( 1/rho d/dx( phi ) ) = dU/dx - (S - S_ave) [S = divu if U is 
c     zero, S = d/dt(divu) if U = (ustar - uold)/dt] with periodic or Neumann
c     boundary conditions, using a tridiagonal solve which detects, and deals with,
c     the singular equations.  In the Neumann case, arbitrarily set the upper right
c     corner to zero to pin the solution.  Note that the RHS of this equation
c     satisfies the solvability constraint that Int[RHS.dV] = 0 by construction.
c     This implies that the normal component takes up the slack:
c     
c                        d/dy( 1/rho d/dy( phi ) ) = dV/dy - S_ave
c     
c     This information should be used to construct the normal gradient of the
c     normal velocity, for the advective/diffusive step, for example.
      
      j = jns
      do i = 2,ine-ins
         icR = ins + i - 1
         icL = icR - 1
         a(i) = r(icL) / rho(icL,j)
         c(i) = r(icR) / rho(icR,j)
         b(i) = - a(i) - c(i)
         s(i) = hx*
     &        (r(icR)*( uExt(icR-ics)-half*hx*(divcor(icR,j)-divcor_ave))
     &        +r(icL)*(-uExt(icL-ics)-half*hx*(divcor(icL,j)-divcor_ave)))
      end do

      if (per .eq. 1) then
         
c     Do left-side periodic BC (since first/last node coincide, use first node only
c     (retain r-stuff here, just to be sure scaling is not destroyed)
         neq = ine - ins
         icL = ice
         icR = ics
         beta = r(icL) / rho(icL,j)
         c(1) = r(icR) / rho(icR,j)
         b(1) = - beta - c(1)
         s(1) = hx*
     &        (r(icR)*( uExt(icR-ics)-half*hx*(divcor(icR,j)-divcor_ave))
     &        +r(icL)*(-uExt(icL-ics)-half*hx*(divcor(icL,j)-divcor_ave)))
         
c     Do right-side periodic on penultimate node
         icL = ice - 1
         icR = ice
         a(neq) = r(icL) / rho(icL,j)
         alpha  = r(icR) / rho(icR,j)
         b(neq) = - a(neq) - alpha
         s(neq) = hx*
     &        (r(icR)*( uExt(icR-ics)-half*hx*(divcor(icR,j)-divcor_ave))
     &        +r(icL)*(-uExt(icL-ics)-half*hx*(divcor(icL,j)-divcor_ave)))
         
c     Solve the equations
         call cyclic(a,b,c,alpha,beta,s,x,neq)
         
      else

c     Solid walls, Nuemann conditions
         neq = ine - ins + 1
         icR = ics
         c(1) = r(icR) / rho(icR,j)
         b(1) = - c(1)
         s(1) = hx*
     &        r(icR)*( uExt(icR-ics)-half*hx*(divcor(icR,j)-divcor_ave))
         
         icL = ice
         a(neq) = r(icL) / rho(icL,j)
         b(neq) = - a(ine-1)
         s(neq) = hx*
     &        r(icL)*(-uExt(icL-ics)-half*hx*(divcor(icL,j)-divcor_ave))
         
c     Solve the equations (we know they're singular, pass the arbitrary value, and
c     a flag that we've already normalized the rhs, in the sense that
c                          Int[dU/dx - (S-S_ave)] == 0
         sVal = zero
         rNormed = .true.
         call tridag_sing(a,b,c,s,x,neq,sVal,rNormed)
         
      end if
      
c     Write the solution into the arg
      do i = 1,neq
         phi(i+ins-1,jns) = x(i)
      end do
      if (per .eq. 1) then
         phi(ine,jns) = phi(ins,jns)
      end if
      
      end
      
      SUBROUTINE tridag_sing(a,b,c,r,u,n,sVal,rNormed)
      INTEGER n,NMAX
      REAL_T a(n),b(n),c(n),r(n),u(n), sVal
      PARAMETER (NMAX=500)
      INTEGER j
      REAL_T bet,gam(NMAX)
      logical rNormed, singular
      singular = .false.
      if(b(1).eq.zero)then
         write(6,*) 'tridag: rewrite equations'
         stop
      end if
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.) then
           if (rNormed .and.  (j .eq. n)) then
              singular = .true.
           else
              pause 'tridag failed'
           end if
        end if
        if (singular) then
           u(j) = sVal
        else
           u(j)=(r(j)-a(j)*u(j-1))/bet
        end if
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END
      
      SUBROUTINE tridag(a,b,c,r,u,n)
      INTEGER n,NMAX
      REAL_T a(n),b(n),c(n),r(n),u(n)
      PARAMETER (NMAX=500)
      INTEGER j
      REAL_T bet,gam(NMAX)
      if(n.gt.NMAX)then
         write(6,*) 'NMAX too small in tridag'
      end if
      if(b(1).eq.0.)pause 'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.)pause 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END

      SUBROUTINE cyclic(a,b,c,alpha,beta,r,x,n)
      INTEGER n,NMAX
      REAL_T alpha,beta,a(n),b(n),c(n),r(n),x(n)
      PARAMETER (NMAX=500)
      INTEGER i
      REAL_T fact,gamma,bb(NMAX),u(NMAX),z(NMAX),den
      if(n.le.2)then
         write(6,*) 'n too small in cyclic'
         stop
      end if
      if(n.gt.NMAX)then
         write(6,*) 'NMAX too small in cyclic'
      end if
      gamma=-b(1)
      bb(1)=b(1)-gamma
      bb(n)=b(n)-alpha*beta/gamma
      do 11 i=2,n-1
        bb(i)=b(i)
11    continue
      call tridag(a,bb,c,r,x,n)
      u(1)=gamma
      u(n)=alpha
      do 12 i=2,n-1
        u(i)=0.
12    continue
      call tridag(a,bb,c,u,z,n)
      den=one+z(1)+beta*z(n)/gamma
      fact=cvmgt(zero,(x(1)+beta*x(n)/gamma)/den,den.eq.zero)
      do 13 i=1,n
        x(i)=x(i)-fact*z(i)
13    continue
      return
      END

c ::: -----------------------------------------------------------
c ::: This routine will compute D (U / dt + sigma G^perp phi)
c :::
c ::: NOTE: phi is node-centered
c :::       sigma is cell-centered      
c :::       U is cell-centered      
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: phi        => node centered data
c ::: DIMS(phi)  => index limits for phi
c ::: sigma      => cell-centered coefficient array
c ::: DIMS(sigma)=> index limits for sigma      
c ::: u          => cell-centered velocity array
c ::: DIMS(u)    => index limits for u      
c ::: divgp     <=  node-centered divergence of gradient
c ::: DIMS(divgp)=> index limits for divgp
c ::: lo,hi      => node-centered subregion of divgp to define
c ::: dx         => cell sizes      
c ::: mult       => scalar multiplicative factor      
c ::: -----------------------------------------------------------

      subroutine FORT_FILTRHS(phi,DIMS(phi),sigma,DIMS(sigma),
     &                        u,DIMS(u),divgp,DIMS(divgp),
     &                        lo,hi,dx,mult,is_rz)

      integer    DIMDEC(phi)
      integer    DIMDEC(sigma)
      integer    DIMDEC(u)
      integer    DIMDEC(divgp)
      integer    lo(SDIM), hi(SDIM)
      REAL_T     dx(SDIM), mult
      REAL_T     phi(DIMV(phi))
      REAL_T     sigma(DIMV(sigma))
      REAL_T     u(DIMV(u),SDIM)
      REAL_T     divgp(DIMV(divgp))
      integer    is_rz
      
      integer    i, j
      REAL_T     p1, p2, p3, p4, denom
      REAL_T     rminv,rpinv,fac
      REAL_T     gpy_ul,gpy_ur,gpy_ll,gpy_lr
      REAL_T     pext,vext,ux,vy

#if BL_PRVERSION==9
c     FOR NINE-POINT STENCIL
      denom = one/(twelve*dx(1)**2) + one/(twelve*dx(2)**2) 
#elif BL_PRVERSION==5
c     FOR FIVE-POINT STENCIL
      denom = one/(four*dx(1)**2) + one/(four*dx(2)**2)
#else
#error "Unknown PRVERSION setting"
#endif

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            p1 = sigma(i-1,j-1)*(-phi(i-1,j-1) + phi(i  ,j-1)
     &                           +phi(i-1,j  ) - phi(i  ,j  ))
            p2 = sigma(i-1,j  )*(-phi(i-1,j+1) + phi(i  ,j+1)
     &                           +phi(i-1,j)   - phi(i  ,j  ))
            p3 = sigma(i  ,j-1)*(-phi(i+1,j-1) + phi(i  ,j-1)
     &                           +phi(i+1,j)   - phi(i  ,j  ))
            p4 = sigma(i  ,j  )*(-phi(i+1,j+1) + phi(i  ,j+1)
     &                           +phi(i+1,j)   - phi(i  ,j  ))
            divgp(i,j) = denom*(p1+p2+p3+p4)
         enddo
      enddo

      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        ux = half*(u(i,j,1)-u(i-1,j,1)+u(i,j-1,1)-u(i-1,j-1,1))/dx(1)
        vy = half*(u(i,j,2)-u(i,j-1,2)+u(i-1,j,2)-u(i-1,j-1,2))/dx(2)
        divgp(i,j) = divgp(i,j) + ux + vy
      enddo
      enddo

      if (is_rz .eq. 1) then

c      extra term added to DGphi for cylindrical coordinates
       do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            rpinv = half / (dfloat(i)+half)
            rminv = half / (dfloat(i)-half)

            pext = ((sigma(i-1,j-1) * (phi(i,j-1) - phi(i,j)) +
     $               sigma(i-1,j  ) * (phi(i,j+1) - phi(i,j))) * rminv -
     $              (sigma(i  ,j-1) * (phi(i,j-1) - phi(i,j)) +
     $               sigma(i  ,j  ) * (phi(i,j+1) - phi(i,j))) * rpinv)

            divgp(i,j) = divgp(i,j) + denom*pext
         enddo
       enddo

c      extra term added to DGbarphi for cylindrical coordinates --
c        here we subtract it since we want D G^perp phi
       fac = one / (twelve * dx(2)**2)
       do j = lo(2), hi(2)
       do i = lo(1), hi(1)
           rminv = fac / (dfloat(i)-half)
           rpinv = fac / (dfloat(i)+half)
   
	   gpy_ur = half * sigma(i  ,j  ) * ( phi(i+1,j+1) + phi(i,j+1) 
     $                                       -phi(i+1,j  ) - phi(i,j  ))
	   gpy_ul = half * sigma(i-1,j  ) * ( phi(i-1,j+1) + phi(i,j+1) 
     $                                       -phi(i-1,j  ) - phi(i,j  ))
	   gpy_lr = half * sigma(i  ,j-1) * ( phi(i+1,j  ) + phi(i,j  ) 
     $                                       -phi(i+1,j-1) - phi(i,j-1))
	   gpy_ll = half * sigma(i-1,j-1) * ( phi(i-1,j  ) + phi(i,j  ) 
     $                                       -phi(i-1,j-1) - phi(i,j-1))
           vext =  (gpy_lr - gpy_ur)*rpinv 
     $            +(gpy_ul - gpy_ll)*rminv
           divgp(i,j) = divgp(i,j) - vext
        enddo
        enddo

        fac = one / (twelve * dx(2))
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          rminv = fac / (dfloat(i)-half)
          rpinv = fac / (dfloat(i)+half)

          vext = (u(i  ,j-1,2) - u(i  ,j  ,2))*rpinv 
     $          +(u(i-1,j  ,2) - u(i-1,j-1,2))*rminv
   
          divgp(i,j) = divgp(i,j) + vext
        enddo
        enddo

      endif

      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
           divgp(i,j) = divgp(i,j) * mult
      enddo
      enddo

      end
