c
c $Id: PROJECTION_2D.F,v 1.32 1999-07-13 22:59:35 marc Exp $
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
      REAL_T dt_inv
   
      dt_inv = one/dt

      do n = 1, SDIM
         do j = lo(2)-1, hi(2)+1
            do i = lo(1)-1, hi(1)+1
               unew(i,j,n) = (unew(i,j,n)-uold(i,j,n))*dt_inv
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

      subroutine FORT_RADMPY(a,DIMS(grid),domlo,domhi,ng,r,nr)
c 
c     multiply A by Radius r
c
c 
c     NOTE: THIS ROUTINE HAS BEEN MODIFIED SO THAT R IS REFLECTED AROUND
c           R=0 AND R=RMAX.
c
      integer    ng,nr
      integer    DIMDEC(grid)
      integer    domlo(2), domhi(2)
      REAL_T     a(ARG_L1(grid)-ng:ARG_H1(grid)+ng, 
     $             ARG_L2(grid)-ng:ARG_H2(grid)+ng)
      REAL_T     r(ARG_L1(grid)-nr:ARG_H1(grid)+nr)

      integer i, j

      do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do i = ARG_L1(grid)-ng, ARG_H1(grid)+ng
           a(i,j) = r(i)*a(i,j)
         end do
      end do

      if (ARG_L1(grid)-ng .lt. domlo(1)) then
         do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do i = ARG_L1(grid)-ng, domlo(1)-1
           a(i,j) = -r(i)*a(i,j)
         end do
         end do
      end if

      if (ARG_H1(grid)+ng .gt. domhi(1)) then
         do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do i = domhi(1)+1, ARG_H1(grid)+ng
           a(i,j) = r(domhi(1))*a(i,j)
         end do
         end do
      end if

      end

      subroutine FORT_RADDIV(a,DIMS(grid),domlo,domhi,ng,r,nr,bogus_value)
c 
c     divide A by Radius r
c
c     NOTE: THIS ROUTINE HAS BEEN MODIFIED SO THAT R IS REFLECTED AROUND
c           R=0 AND R=RMAX.
c
      integer    ng,nr
      integer    DIMDEC(grid)
      integer    domlo(2), domhi(2)
      REAL_T     a(ARG_L1(grid)-ng:ARG_H1(grid)+ng, 
     $             ARG_L2(grid)-ng:ARG_H2(grid)+ng)
      REAL_T     r(ARG_L1(grid)-nr:ARG_H1(grid)+nr)
      REAL_T     bogus_value

      integer i, j

      do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do i = ARG_L1(grid)-ng, ARG_H1(grid)+ng
           a(i,j) = a(i,j)/r(i)
         end do
      end do

      if (ARG_L1(grid)-ng .lt. domlo(1)) then
         do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do i = ARG_L1(grid)-ng, domlo(1)-1
           a(i,j) = -a(i,j)/r(i)
         end do
         end do
      end if

      if (ARG_H1(grid)+ng .gt. domhi(1)) then
         do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do i = domhi(1)+1, ARG_H1(grid)+ng
           a(i,j) = a(i,j)/r(domhi(1))
         end do
         end do
      end if

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
c        call flush(6)
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
     &     domlo, domhi, dr, is_rz, imax) 
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
      integer nghost 
      integer domlo(SDIM), domhi(SDIM)
      integer DIMDEC(dat)
      integer DIMDEC(rhs)
      REAL_T  dr
      REAL_T  rcen(DIM1(dat))
      REAL_T  dat(DIMV(dat))
      REAL_T  rhs(DIMV(rhs))
      integer is_rz, imax

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
     &                            dat(i  ,j-1)+dat(i  ,j) )
            end do
         end do

      else

#if BL_PRVERSION == 9
         do j=ARG_L2(rhs)+nghost,ARG_H2(rhs)-nghost
            do i=ARG_L1(rhs)+nghost,ARG_H1(rhs)-nghost
               if (i .eq. imax) then
                  rhi = rcen(i-1)
               else 
                  rhi = rcen(i)
               end if
               if (i .eq. 0) then
                  rlo = rcen(i)
               else
                  rlo = rcen(i-1)
               end if

               rhs(i,j) = fourth*(rlo * (dat(i-1,j-1) + dat(i-1,j)) +
     &                            rhi * (dat(i  ,j-1) + dat(i  ,j)))
            end do
         end do
         factor = dr/24.0
         do j=ARG_L2(rhs)+nghost,ARG_H2(rhs)-nghost
            do i=ARG_L1(rhs)+nghost,ARG_H1(rhs)-nghost
               if (i .eq. imax) then
                  rhi = -one
               else 
                  rhi = one
               end if
               if (i .eq. 0) then
                  rlo = -one
               else
                  rlo = one
               end if
               rhs(i,j) = rhs(i,j) + factor *
     &              (rlo * (dat(i-1,j-1) + dat(i-1,j)) -
     &               rhi * (dat(i  ,j-1) + dat(i  ,j)))
            end do
         end do

#elif BL_PRVERSION == 5
         do j=ARG_L2(rhs)+nghost,ARG_H2(rhs)-nghost
            do i=ARG_L1(rhs)+nghost,ARG_H1(rhs)-nghost
               if (i .eq. imax) then
                  rhi = rcen(i-1)
               else 
                  rhi = rcen(i)
               end if
               if (i .eq. 0) then
                  rlo = rcen(i)
               else 
                  rlo = rcen(i-1)
               end if
               rhs(i,j) = fourth*(rlo * (dat(i-1,j-1) + dat(i-1,j)) +
     &                            rhi * (dat(i  ,j-1) + dat(i  ,j)))
               if (i .eq. 0) rhs(i,j) = half * rhs(i,j)
            end do
         end do
#endif
      end if

      end

      subroutine FORT_EXTRAP_PROJ(DIMS(u),u,DIMS(divu),divu,DIMS(rho),rho,
     &     r_lo,r_hi,r,DIMS(uExt),uExt,DIMS(divuExt),divuExt,
     &     DIMS(rhoExt),rhoExt,face, zeroIt)
      implicit none

      integer DIMDEC(u)
      integer DIMDEC(divu)
      integer DIMDEC(rho)
      integer DIMDEC(uExt)
      integer DIMDEC(divuExt)
      integer DIMDEC(rhoExt)
      integer r_lo,r_hi
      integer face
      REAL_T      u(DIMV(u),SDIM)
      REAL_T   divu(DIMV(divu))
      REAL_T    rho(DIMV(rho))
      REAL_T   uExt(DIMV(uExt))
      REAL_T   divuExt(DIMV(divuExt))
      REAL_T   rhoExt(DIMV(rhoExt))
      REAL_T   r(r_lo:r_hi)
      integer  zeroIt

c local variables
      integer ics,ice,jcs,jce
      integer ifs,ife,jfs,jfe
      integer if,jf
      REAL_T divu_ave, vtot
      REAL_T max_divu, min_divu
      REAL_T max_pert, small_pert
      parameter ( small_pert = SMALL)
      integer i,j

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3

      ics = ARG_L1(u)
      ice = ARG_H1(u)
      jcs = ARG_L2(u)
      jce = ARG_H2(u)

      ifs = ARG_L1(uExt)
      ife = ARG_H1(uExt)
      jfs = ARG_L2(uExt)
      jfe = ARG_H2(uExt)

      zeroIt = 0
      divu_ave = zero
      vtot = zero

      if (face .eq. XLO) then
         if=ife
         max_divu = divu(ice-1,jcs)
         min_divu = max_divu
         do j = jcs, jce
            uExt(if,j)    = half*(three*u(ice-1,j,2)    - u(ice,j,2))
            divuExt(if,j) = half*(three*divu(ice-1,j) - divu(ice,j))
            rhoExt(if,j)  = half*(three*rho(ice-1,j)    - rho(ice,j))
            
            divu_ave = divu_ave+divuExt(if,j)
            vtot = vtot+one
            max_divu = max(max_divu,divuExt(if,j))
            min_divu = min(min_divu,divuExt(if,j))
         end do
         divu_ave = divu_ave/vtot

         max_pert = ABS(divuExt(if,jcs)-divu_ave)
         do j = jcs+1, jce
            max_pert = MAX(max_pert,ABS(divuExt(if,j)-divu_ave))
         end do
         max_pert = max_pert/(ABS(divu_ave)+small_pert)
         if ((max_divu.eq.zero.and.min_divu.eq.zero)
     &        .or.(max_pert.le.small_pert)) then
            zeroIt = 1
            return
         end if
         do j = jcs, jce
            divuExt(if,j) = divuExt(if,j) - divu_ave
         end do
      else if (face .eq. YLO) then
         jf = jfe
         max_divu = divu(ics,jce-1)
         min_divu = max_divu
         do i = ics, ice
            uExt(i,jf)    = half*(three*u(i,jce-1,1)    - u(i,jce,1))
            divuExt(i,jf) = half*(three*divu(i,jce-1) - divu(i,jce))
            rhoExt(i,jf)  = half*(three*rho(i,jce-1)    - rho(i,jce))
            
            divu_ave = divu_ave+r(i)*divuExt(i,jf)
            vtot = vtot+r(i)
            max_divu = max(max_divu,divuExt(i,jf))
            min_divu = min(min_divu,divuExt(i,jf))
         end do
         divu_ave = divu_ave/vtot

         max_pert = ABS(divuExt(ics,jf)-divu_ave)
         do i = ics+1, ice
            max_pert = MAX(max_pert,ABS(divuExt(i,jf)-divu_ave))
         end do
         max_pert = max_pert/(ABS(divu_ave)+small_pert)
         if ((max_divu.eq.zero.and.min_divu.eq.zero)
     &        .or.(max_pert.le.small_pert)) then
            zeroIt = 1
            return
         end if
         do i = ics, ice
            divuExt(i,jf) = divuExt(i,jf) - divu_ave
         end do
      else if (face .eq. XHI) then
         if = ife
         max_divu = divu(ics+1,jcs)
         min_divu = max_divu
         do j = jcs, jce
            uExt(if,j)    = half*(three*u(ics+1,j,2)    - u(ics,j,2))
            divuExt(if,j) = half*(three*divu(ics+1,j) - divu(ics,j))
            rhoExt(if,j)  = half*(three*rho(ics+1,j)    - rho(ics,j))
            
            divu_ave = divu_ave+divuExt(if,j)
            vtot = vtot+one
            max_divu = max(max_divu,divuExt(if,j))
            min_divu = min(min_divu,divuExt(if,j))
         end do
         divu_ave = divu_ave/vtot

         max_pert = ABS(divuExt(if,jcs)-divu_ave)
         do j = jcs+1, jce
            max_pert = MAX(max_pert,ABS(divuExt(if,j)-divu_ave))
         end do
         max_pert = max_pert/(ABS(divu_ave)+small_pert)
         if ((max_divu.eq.zero.and.min_divu.eq.zero)
     &        .or.(max_pert.le.small_pert)) then
            zeroIt = 1
            return
         end if
         do j = jcs, jce
            divuExt(if,j) = divuExt(if,j) - divu_ave
         end do
      else if (face .eq. YHI) then
         jf = jfe
         max_divu = divu(ics,jcs+1)
         min_divu = max_divu
         do i = ics, ice
            uExt(i,jf)    = half*(three*u(i,jcs+1,1)    - u(i,jcs,1))
            divuExt(i,jf) = half*(three*divu(i,jcs+1) - divu(i,jcs))
            rhoExt(i,jf)  = half*(three*rho(i,jcs+1)    - rho(i,jcs))
            
            divu_ave = divu_ave+r(i)*divuExt(i,jf)
            vtot = vtot+r(i)
            max_divu = max(max_divu,divuExt(i,jf))
            min_divu = min(min_divu,divuExt(i,jf))
         end do
         divu_ave = divu_ave/vtot

         max_pert = ABS(divuExt(ics,jf)-divu_ave)
         do i = ics+1, ice
            max_pert = MAX(max_pert,ABS(divuExt(i,jf)-divu_ave))
         end do
         max_pert = max_pert/(ABS(divu_ave)+small_pert)
         if ((max_divu.eq.zero.and.min_divu.eq.zero)
     &        .or.(max_pert.le.small_pert)) then
            zeroIt = 1
            return
         end if
         do i = ics, ice
            divuExt(i,jf) = divuExt(i,jf) - divu_ave
         end do
      endif
#undef XLO
#undef YLO
#undef XHI
#undef YHI

      end


      subroutine FORT_HGPHIBC(hx,r,uExt, divuExt,rhoExt,length,
     &                        DIMS(phi),phi,face,per)
c
c     Compute the value of phi for hgproj to be used at an  upper, y-outflow face,
c     assuming that the tangential velocity on the edges of the outflow boundary
c     are either zero or periodic.  Algorithm incorrect for stress-free side walls,
c     or if there is outflow on any other boundary.
c
      implicit none

      integer DIMDEC(phi)
      integer length
      integer per
      integer face
      REAL_T    phi(DIMV(phi))
      REAL_T      r(0:length-1)
      REAL_T      uExt(0:length-1)
      REAL_T      divuExt(0:length-1)
      REAL_T      rhoExt(0:length-1)
      REAL_T hx
      
c     Local variables
      integer NstripMAX
      parameter (NstripMAX = 1000)
      integer ins, ine, ics, ice
      integer i, j, icL, icR, neq
      REAL_T a(NstripMAX), b(NstripMAX), c(NstripMAX), s(NstripMAX)
      REAL_T x(NstripMAX), alpha, beta, sVal
      logical rNormed
      REAL_T vtot, phitot, rnode, hdr

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3
      
      ics = 0
      ice = length-1

      if (face .eq. XLO .or. face .eq. XHI) then
         ins = ARG_L2(phi)
         ine = ARG_H2(phi)
      else if (face .eq. YLO .or. face .eq. YHI) then
         ins = ARG_L1(phi)
         ine = ARG_H1(phi)
      endif

      
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
      do i = 2,ine-ins
         icR = ics + i - 1
         icL = ics + i - 2
         a(i) = r(icL) / rhoExt(icL)
         c(i) = r(icR) / rhoExt(icR)
         b(i) = - a(i) - c(i)
         rnode = half*(r(icL)+r(icR))
         s(i) = (r(icR)*uExt(icR)-r(icL)*uExt(icL))*hx
     &        -  rnode*half*(divuExt(icL)+divuExt(icR))*hx*hx
      end do

      if (per .eq. 1) then
         
c     Do left-side periodic BC (since first/last node coincide, use first node only
c     (retain r-stuff here, just to be sure scaling is not destroyed)
         neq = ine - ins
         icL = ice
         icR = ics
         beta = r(icL) / rhoExt(icL)
         c(1) = r(icR) / rhoExt(icR)
         b(1) = - beta - c(1)
         rnode = half*(r(icL)+r(icR))
         s(1) = (r(icR)*uExt(icR)-r(icL)*uExt(icL))*hx
     &        -  rnode*half*(divuExt(icL)+divuExt(icR))*hx*hx
         
c     Do right-side periodic on penultimate node
         icL = ice - 1
         icR = ice
         a(neq) = r(icL) / rhoExt(icL)
         alpha  = r(icR) / rhoExt(icR)
         b(neq) = - a(neq) - alpha
         s(neq) = (r(icR)*uExt(icR)-r(icL)*uExt(icL))*hx
     &        -  rnode*half*(divuExt(icL)+divuExt(icR))*hx*hx
         
c     Solve the equations
         call cyclic(a,b,c,alpha,beta,s,x,neq)
         
      else

c     Solid walls, Nuemann conditions
c
c     Note (because I always forget...):: 
c        Can derive these expressions by integrating the FE equations directly
c     (ignoring r-z):
c             Int[1/rho grad(phi)]dx  = Int[du/dx.dpsi/dx dx]- Int[S.psi.dx]
c     But, Int[du/dx.dpsi/dx.dx] = Int[du/dx.psi.dx] - Int[u.dpsi/dx]
c                                = - u Int[dpsi/dx]
c     And, since psi=0 on boundary of FE, and (cell-centered) u constant over FE
c     and since S is constant, Int[psi.dx] = dx/2 (for psi=x/dx), and finally
c     Int[dpsi/dx.dx] = 1/dx, then RHS = -u - S.dx/2...then mult all by dx
c         
         hdr = half*(r(ics+1) - r(ics))
         neq = ine - ins + 1
         icR = ics
         c(1) = r(icR) / rhoExt(icR-ics)
         b(1) = - c(1)
         rnode = r(ics) - hdr
         s(1) = r(icR)*uExt(icR)*hx
     &        -  rnode*half*divuExt(icR)*hx*hx
        
         icL = ice
         a(neq) = r(icL) / rhoExt(icL)
         b(neq) = - a(neq)
         rnode = r(ice) + hdr
         s(neq) = -r(icL)*uExt(icL)*hx
     &        -  rnode*half*divuExt(icL)*hx*hx
         
         
c     Solve the equations (we know they're singular, pass the arbitrary value, and
c     a flag that we've already normalized the rhs, in the sense that
c                          Int[dU/dx - (S-S_ave)] == 0
         sVal = zero
         rNormed = .true.
         call tridag_sing(a,b,c,s,x,neq,sVal,rNormed)
         
      end if
      
c     Try normalizing phi to average to zero
      phitot = zero
      vtot = zero
      do i = 1,ine-ins-1
         rnode = r(ics+i-1) - hdr
         phitot = phitot + x(i)*rnode
         vtot = vtot + rnode
      end do
      do i = ine-ins,neq
         rnode = r(ics+i-2) + hdr
         phitot = phitot + x(i)*rnode
         vtot = vtot + rnode
      end do
      phitot = phitot / vtot
      do i = 1,ine-ins-1
         x(i) = x(i) - phitot
      end do
      do i = ine-ins,neq
         x(i) = x(i) - phitot
      end do
      
c     Write the solution into the arg
      if (face .eq. XLO .or. face .eq. XHI) then
         do i = 1,neq
            phi(ARG_L1(phi),i+ARG_L2(phi)-1) = x(i)
         end do
         if (per .eq. 1) then
            phi(ARG_L1(phi),ARG_H2(phi)) = phi(ARG_L1(phi),ARG_L2(phi))
         end if
      else if (face .eq. YLO .or. face .eq. YHI) then
         do i = 1,neq
            phi(i+ARG_L1(phi)-1,ARG_L2(phi)) = x(i)
         end do
         if (per .eq. 1) then
            phi(ARG_H1(phi),ARG_L2(phi)) = phi(ARG_L1(phi),ARG_L2(phi))
         end if
      endif
#undef XLO
#undef YLO
#undef XHI
#undef YHI
      
      end
      SUBROUTINE tridag_sing(a,b,c,r,u,n,sVal,rNormed)
      INTEGER n,NMAX
      REAL_T a(n),b(n),c(n),r(n),u(n), sVal
      REAL_T eps
      PARAMETER (NMAX=1000)
      PARAMETER (eps=1.e-12)
      INTEGER j
      REAL_T bet,gam(NMAX)
      logical rNormed, singular
      singular = .false.
      if((b(1).eq.zero).or.(ABS(b(n)).le.eps))then
         write(6,*) 'tridag: rewrite equations'
         stop
      end if
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(ABS(bet).le.eps) then
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
      PARAMETER (NMAX=1000)
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
      PARAMETER (NMAX=1000)
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
     &                        divgp,DIMS(divgp),
     &                        lo,hi,dx,mult,is_rz)
      implicit none

      integer    DIMDEC(phi)
      integer    DIMDEC(sigma)
      integer    DIMDEC(divgp)
      integer    lo(SDIM), hi(SDIM)
      REAL_T     dx(SDIM), mult
      REAL_T     phi(DIMV(phi))
      REAL_T     sigma(DIMV(sigma))
      REAL_T     divgp(DIMV(divgp))
      integer    is_rz
      
      integer    i, j
      REAL_T     p1, p2, p3, p4, denom
      REAL_T     rminv,rpinv,fac
      REAL_T     gpy_ul,gpy_ur,gpy_ll,gpy_lr
      REAL_T     pext, vext

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

      endif

      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
           divgp(i,j) = divgp(i,j) * mult
      enddo
      enddo

      end
