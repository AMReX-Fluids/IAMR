#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "MACOUTFLOWBC_F.H"
#include "ArrayLim.H"

#define SDIM 2

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define SMALL 1.0e-10
#else
#define SMALL 1.0d-10
#endif


c *************************************************************************
c ** EXTRAP_MAC
c *************************************************************************

      subroutine FORT_EXTRAP_MAC(DIMS(u0),u0,DIMS(u1),u1,DIMS(div),divu,DIMS(rho),rho,
     &                         r_len,redge,DIMS(divuExt),divuExt,
     &                         DIMS(rhoExt),rhoExt,dx,lo,hi,face,per,zeroIt)
c
c     Compute the value of phi for macproj 
c
c     (subtract divu_ave twice due to precision problems)
      implicit none

      integer DIMDEC(u0)
      integer DIMDEC(u1)
      integer DIMDEC(div)
      integer DIMDEC(rho)
      integer r_len
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(divuExt)
      integer DIMDEC(rhoExt)
      REAL_T      u0(DIMV(u0))
      REAL_T      u1(DIMV(u1))
      REAL_T divu(DIMV(div))
      REAL_T    rho(DIMV(rho))
      REAL_T  redge(0:r_len-1)
      REAL_T   divuExt(DIMV(divuExt))
      REAL_T   rhoExt(DIMV(rhoExt))
      REAL_T   dx(2)
      integer face
      integer per
      integer zeroIt
      
c     Local variables
      REAL_T small_pert
      REAL_T small_udiff
      parameter ( small_pert = SMALL)
      parameter (small_udiff = SMALL)
      integer i, j, neq,n
      REAL_T divu_ave1,divu_ave2
      REAL_T max_divu, min_divu, max_pert
      REAL_T diff
      REAL_T rc,hx,hy
      integer ics,ice,jcs,jce
      integer ifs,ife,jfs,jfe
      integer if,jf,ic,jc
c     NOTE: Assumes that rho at edge between i, i-1 = half*(rho(i)+rho(i-1))
c             (1) Linear fit of rho between nodes
c             (2) rho, divu on same boxes (box)
c             (3) phi is on box, shifted up one
c             (4) u is edge-based, on surroundingNodes(box)

c     Compute average of divu over outflow bc.  Set trivial solution if average
c     is zero, or if divu is constant
#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3
      ics = ARG_L1(rho)
      ice = ARG_H1(rho)
      jcs = ARG_L2(rho)
      jce = ARG_H2(rho)

      ifs = lo(1)
      ife = hi(1)
      jfs = lo(2)
      jfe = hi(2)

      hx = dx(1)
      hy = dx(2)

      zeroIt = 0

      if (face .eq. XLO) then

         if = ifs
         max_divu = divu(ics,jcs)
         min_divu = max_divu
         do j = jcs, jce
            divuExt(j,if) = divu(ics,j)
            rhoExt(j,if)  = rho(ics,j)
            max_divu = max(max_divu,divuExt(j,if))
            min_divu = min(min_divu,divuExt(j,if))
         end do

c        Here we modify divuExt to include the velocity terms.
         do j = jcs, jce
            divuExt(j,if) = redge(j-jcs)*(divuExt(j,if)*hy*hy - (u1(ics,j+1)-u1(ics,j))*hy)
         end do

         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave2,face)

         max_pert = ABS(divuExt(jcs,if))
         do j = jcs, jce
            max_pert = MAX(max_pert,ABS(divuExt(j,if)))
         end do
      
c        Make sure u_mac is periodic
         if (per .eq. 1) then
           diff = u1(ics,jcs)-u1(ics,jce+1)
           if (ABS(diff) .gt. small_udiff) then
              write(6,*) 'EXTRAPMAC: FACE XLO : umac not periodic'
              write(6,*) 'V AT    TOP: ',u1(ics,jce+1)
              write(6,*) 'V AT BOTTOM: ',u1(ics,jcs  )
              call bl_abort(" ")
           endif
         endif

      else if (face .eq. YLO) then

         jf = jfs
         max_divu = divu(ics,jcs)
         min_divu = max_divu
         do i = ics, ice
            divuExt(i,jf) = divu(i,jcs)
            rhoExt(i,jf)  = rho(i,jcs)
            max_divu = max(max_divu,divuExt(i,jf))
            min_divu = min(min_divu,divuExt(i,jf))
         end do

c        Here we modify divuExt to include the velocity terms.
         do i = ics, ice
            rc = half*(redge(i+1-ics)+redge(i-ics))
            divuExt(i,jf) = rc*divuExt(i,jf)*hx*hx - 
     $                      (redge(i+1-ics)*u0(i+1,jcs)-redge(i-ics)*u0(i,jcs))*hx
         end do

         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(ics,jf))
         do i = ics, ice
            max_pert = MAX(max_pert,ABS(divuExt(i,jf)))
         end do
      
c        Make sure u_mac is periodic
         if (per .eq. 1) then
           diff = u0(ics,jcs)-u0(ice+1,jcs)
           if (ABS(diff) .gt. small_udiff) then
              write(6,*) 'EXTRAPMAC: FACE YLO : umac not periodic'
              write(6,*) 'U AT LEFT: ',u0(ics  ,jcs)
              write(6,*) 'U AT RGHT: ',u0(ice+1,jcs)
              call bl_abort(" ")
           endif
         endif

      else if (face .eq. XHI) then

         if = ife
         max_divu = divu(ice,jcs)
         min_divu = max_divu
         do j = jcs, jce
            divuExt(j,if) = divu(ice,j)
            rhoExt(j,if)  = rho(ice,j)
            max_divu = max(max_divu,divuExt(j,if))
            min_divu = min(min_divu,divuExt(j,if))
         end do

c        Here we modify divuExt to include the velocity terms.
         do j = jcs, jce
            divuExt(j,if) = redge(j-jcs)*(divuExt(j,if)*hy*hy - (u1(ice,j+1)-u1(ice,j))*hy)
         end do

         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave2,face)

         max_pert = ABS(divuExt(jcs,if))
         do j = jcs, jce
            max_pert = MAX(max_pert,ABS(divuExt(j,if)))
         end do
      
c        Make sure u_mac is periodic
         if (per .eq. 1) then
           diff = u1(ice,jcs)-u1(ice,jce+1)
           if (ABS(diff) .gt. small_udiff) then
              write(6,*) 'EXTRAPMAC: FACE XHI : umac not periodic'
              write(6,*) 'V AT    TOP: ',u1(ice,jce+1)
              write(6,*) 'V AT BOTTOM: ',u1(ice,jcs  )
              call bl_abort(" ")
           endif
         endif

      else if (face .eq. YHI) then

         jf = jfe
         max_divu = divu(ics,jce)
         min_divu = max_divu
         do i = ics, ice
            divuExt(i,jf) = divu(i,jce)
            rhoExt(i,jf)  = rho(i,jce)
            max_divu = max(max_divu,divuExt(i,jf))
            min_divu = min(min_divu,divuExt(i,jf))
         end do

c        Here we modify divuExt to include the velocity terms.
         do i = ics, ice
            rc = half*(redge(i+1-ics)+redge(i-ics))
            divuExt(i,jf) = rc*divuExt(i,jf)*hx*hx - 
     $                      (redge(i+1-ics)*u0(i+1,jce)-redge(i-ics)*u0(i,jce))*hx
         end do

         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave2,face)

         max_pert = ABS(divuExt(ics,jf))
         do i = ics, ice
            max_pert = MAX(max_pert,ABS(divuExt(i,jf)))
         end do
      
c        Make sure u_mac is periodic
         if (per .eq. 1) then
           diff = u0(ics,jce)-u0(ice+1,jce)
           if (ABS(diff) .gt. small_udiff) then
              write(6,*) 'EXTRAPMAC: FACE YHI : umac not periodic'
              write(6,*) 'U AT LEFT: ',u0(ics  ,jce)
              write(6,*) 'U AT RGHT: ',u0(ice+1,jce)
              call bl_abort(" ")
           endif
         endif

      endif
      
c  check to see if we should zero phi
         max_pert = max_pert/(ABS(divu_ave1+divu_ave2)+small_pert)
      if ((max_divu.eq.zero.and.min_divu.eq.zero)
     &     .or.(max_pert.le.small_pert)) then
         zeroIt = 1
      end if
      end

c *************************************************************************
c ** SUBTRACTAVG
c *************************************************************************

      subroutine subtractavg(DIMS(divu),divu,redge,r_len,lo,hi,divu_ave,face)
      implicit none
      integer DIMDEC(divu)
      integer r_len
      integer lo(SDIM),hi(SDIM)
      REAL_T  redge(0:r_len-1)
      REAL_T divu(DIMV(divu))
      REAL_T divu_ave
      integer face

      integer i,j
      REAL_T rcen
      REAL_T vtot

      divu_ave = zero
      vtot = zero

      if (face .eq. XLO .or. face .eq. XHI) then
         i = lo(1)
         do j=lo(2),hi(2)
            vtot = vtot+one
            divu_ave = divu_ave+divu(j,i)
         enddo
         divu_ave = divu_ave/vtot
         do j=lo(2),hi(2)
            divu(j,i) = divu(j,i) - divu_ave
         enddo
      elseif (face .eq. YLO .or. face .eq. YHI) then
         j = lo(2)
         do i=lo(1),hi(1)
            rcen = half*(redge(i)+redge(i+1))
            vtot = vtot+rcen
            divu_ave = divu_ave+rcen*divu(i,j)
         enddo
         divu_ave = divu_ave/vtot
         do i=lo(1),hi(1)
            divu(i,j) = divu(i,j) - divu_ave
         enddo
      else 
         print*, "bad value of face in subtractavg"
      endif

      end
#undef XLO
#undef YLO
#undef XHI
#undef YHI


c *************************************************************************
c ** MACRELAX
c *************************************************************************

#define DGX (beta(i)*phi(i-1) - (beta(i)+beta(i+1))*phi(i) \
            +beta(i+1)*phi(i+1))*(hxsqinv)

      subroutine FORT_MACRELAX(DIMS(rhs),rhs,DIMS(beta),beta,DIMS(phi),phi,
     &                        lo,hi,h,isPeriodic,niter)
      implicit none
      integer DIMDEC(beta)
      integer DIMDEC(rhs)
      integer DIMDEC(phi)
      REAL_T beta(DIM1(beta))
      REAL_T rhs(DIM1(rhs))
      REAL_T phi(DIM1(phi))
      REAL_T h(SDIM)
      integer lo(SDIM),hi(SDIM)
      integer isPeriodic(SDIM)

c Local variables
      integer redblack
      integer ics,ice
      integer i,iter
      REAL_T lam, dg
      REAL_T phitot,vtot
      logical setSingularPoint
      REAL_T maxnorm
      REAL_T hxsqinv
      integer niter

      hxsqinv = 1.0D0/(h(1)*h(1))
      ics = lo(1)
      ice = hi(1)

      setSingularPoint = .false.
      call setmacbc(DIMS(phi),phi,lo,hi,isPeriodic,setSingularPoint)

      do iter = 1,niter
         do redblack = 0,1
            do i=ics+redblack,ice,2
               dg = DGX
               lam = -one/(hxsqinv*(beta(i)+beta(i+1)))
               phi(i) = phi(i) + lam*(rhs(i)-dg)
            enddo
           call setmacbc(DIMS(phi),phi,lo,hi,isPeriodic,setSingularPoint)
         end do
      end do

      end

c *************************************************************************
c ** MACSUBTRACTAVGPHI
c *************************************************************************

      subroutine FORT_MACSUBTRACTAVGPHI(DIMS(phi),phi,r_lo,r_hi,r,lo,hi,
     &                                  isPeriodic)
      implicit none
      integer DIMDEC(phi)
      REAL_T phi(DIM1(phi))
      integer r_lo,r_hi
      REAL_T r(r_lo:r_hi)
      integer lo(SDIM),hi(SDIM)
      integer isPeriodic(SDIM)

      REAL_T phitot,vtot
      integer ics,ice,ine
      integer i
      logical setSingularPoint

      ics = lo(1)
      ice = hi(1)
      setSingularPoint = .false.

      phitot = zero
      vtot   = zero
      do i=ics,ice
         phitot = phitot+phi(i)*r(i)
         vtot = vtot + r(i)
      enddo
      phitot = phitot/vtot
      do i=ics,ice
         phi(i) = phi(i) - phitot
      enddo

      call setmacbc(DIMS(phi),phi,lo,hi,isPeriodic,setSingularPoint)
      
      end

c *************************************************************************
c ** MACRESID
c *************************************************************************

      subroutine FORT_MACRESID(DIMS(rhs),rhs,DIMS(beta),beta,DIMS(phi),phi,
     &                    DIMS(resid),resid,lo,hi,h,isPeriodic,maxnorm)
      implicit none
      integer DIMDEC(beta)
      integer DIMDEC(rhs)
      integer DIMDEC(phi)
      integer DIMDEC(resid)
      REAL_T beta(DIM1(beta))
      REAL_T rhs(DIM1(rhs))
      REAL_T phi(DIM1(phi))
      REAL_T resid(DIM1(resid))
      integer isPeriodic(SDIM)
      REAL_T h(SDIM)
      integer lo(SDIM),hi(SDIM)

      REAL_T maxnorm
      integer i
      REAL_T hxsqinv

      hxsqinv = one/(h(1)*h(1))

      maxnorm = zero
      do i = lo(1),hi(1)
         resid(i) = rhs(i)-(DGX)
         maxnorm = max(maxnorm,ABS(resid(i)))
      enddo

      end

c *************************************************************************
c ** SETMACBC
c *************************************************************************

      subroutine setmacbc(DIMS(phi),phi,lo,hi,isPeriodic,setSingularPoint)
      implicit none
      integer DIMDEC(phi)
      REAL_T phi(DIM1(phi))
      integer lo(SDIM),hi(SDIM)
      integer isPeriodic(SDIM)
      integer ics,ice
      logical setSingularPoint
      ics = lo(1)
      ice = hi(1)

      if (isPeriodic(1) .NE. 1 .and. setSingularPoint) then
         phi(ice)= zero
      endif

      if (isPeriodic(1).eq.1) then
         phi(ics-1) = phi(ice)
         phi(ice+1) = phi(ics)
      else
         phi(ics-1) = phi(ics)
         phi(ice+1) = phi(ice)
      endif

      end

c *************************************************************************
c ** MAC_SHIFT_PHI
c *************************************************************************

      subroutine FORT_MAC_SHIFT_PHI(DIMS(out),out,DIMS(in),in,face)
      implicit none
      integer face
      integer DIMDEC(in)
      integer DIMDEC(out)
      REAL_T in(DIMV(in))
      REAL_T out(DIMV(out))
      integer i,j
#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3

      if (face .eq. XLO .or. face .eq. XHI) then
         do j= ARG_L2(out),ARG_H2(out)
            do i = ARG_L1(out),ARG_H1(out)
               out(i,j) = in(j,i)
            enddo
         enddo
      else if (face .eq. YLO .or. face .eq. YHI) then
         do j= ARG_L2(out),ARG_H2(out)
            do i = ARG_L1(out),ARG_H1(out)
               out(i,j) = in(i,j)
            enddo
         enddo
      endif
#undef XLO
#undef YLO
#undef XHI
#undef YHI

      end

c *************************************************************************
c ** MAC_RESHIFT_PHI
c *************************************************************************


      subroutine FORT_MAC_RESHIFT_PHI(DIMS(out),out,DIMS(in),in,face)
      implicit none
      integer face
      integer DIMDEC(in)
      integer DIMDEC(out)
      REAL_T in(DIMV(in))
      REAL_T out(DIMV(out))
      integer i,j
#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3

      if (face .eq. XLO .or. face .eq. XHI) then
         do j= ARG_L2(out),ARG_H2(out)
            do i = ARG_L1(out),ARG_H1(out)
               out(i,j) = in(j,i)
            enddo
         enddo
      else if (face .eq. YLO .or. face .eq. YHI) then
         do j= ARG_L2(out),ARG_H2(out)
            do i = ARG_L1(out),ARG_H1(out)
               out(i,j) = in(i,j)
            enddo
         enddo
      endif
#undef XLO
#undef YLO
#undef XHI
#undef YHI

      end


c *************************************************************************
c ** SOLVEMAC
c *************************************************************************

      subroutine FORT_SOLVEMAC(p, DIMS(p),dest0, DIMS(dest0),
     &                         source,DIMS(source), sigma, DIMS(sigma),
     &                         cen, DIMS(cen),
     $                         r,DIMS(r), w, DIMS(w),z, DIMS(z),
     &                         x, DIMS(x),lo, hi, h,
     $                         isPeriodic, maxiter, tol, abs_tol,max_jump,norm)

      implicit none

      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(p)
      integer DIMDEC(dest0)
      integer DIMDEC(source)
      integer DIMDEC(sigma)
      integer DIMDEC(cen)
      integer DIMDEC(r)
      integer DIMDEC(w)
      integer DIMDEC(z)
      integer DIMDEC(x)
      REAL_T      p(DIM1(p))
      REAL_T  dest0(DIM1(dest0))
      REAL_T source(DIM1(source))
      REAL_T  sigma(DIM1(sigma))
      REAL_T    cen(DIM1(cen))
      REAL_T      r(DIM1(r))
      REAL_T      w(DIM1(w))
      REAL_T      z(DIM1(z))
      REAL_T      x(DIM1(x))
      REAL_T h(SDIM)
      integer isPeriodic(SDIM)
      REAL_T norm
      REAL_T tol
      integer maxiter
      REAL_T abs_tol,max_jump

c     Local variables
      integer i,iter
      integer istart,iend
      REAL_T alpha, beta, rho, rho_old
      REAL_T goal
      REAL_T norm0
      logical setSingularPoint

      istart = lo(1)
      iend = hi(1)

      setSingularPoint = .false.

      do i = lo(1)-1,hi(1)+1 
         dest0(i) = p(i)
         p(i) = zero
      enddo

      do i=ARG_L1(w),ARG_H1(w)
         w(i) = zero
      enddo

      call setmacbc(DIMS(dest0),dest0,lo,hi,isPeriodic,setSingularPoint)
      call makemacdgphi(dest0,DIMS(dest0),w,DIMS(w),sigma,DIMS(sigma),
     &               lo,hi,h,isPeriodic,setSingularPoint)

      do i = istart, iend 
         r(i) = source(i) - w(i)
      enddo
      
      rho = zero
      norm0 = zero
      do i = istart, iend 
         norm0 = max(norm0, abs(r(i)))
         z(i) = r(i) 
         rho = rho + z(i) * r(i)
      enddo
      norm = norm0
      
      goal = max(tol*norm0,abs_tol)
      
      if (norm0 .le. goal) then
         do i = istart, iend 
            p(i) = dest0(i)
         enddo
         return
      endif

      do i = istart, iend 
         x(i) = zero
         p(i) = z(i)
      enddo
      
      iter = 0
      
 100  continue
      
      do i=ARG_L1(w),ARG_H1(w)
         w(i) = zero
      enddo
      
      call setmacbc(DIMS(p),p,lo,hi,isPeriodic,setSingularPoint)
      call makemacdgphi(p,DIMS(p),w,DIMS(w),sigma,DIMS(sigma),
     &     lo,hi,h,isPeriodic,setSingularPoint)
      
      alpha = zero
      do i = istart, iend 
         alpha = alpha + p(i)*w(i)
      enddo
      
      alpha = rho / alpha
      rho_old = rho
      rho = zero
      norm = zero
      do i = istart, iend 
         x(i) = x(i) + alpha * p(i)
         r(i) = r(i) - alpha * w(i)
         z(i) = r(i) 
         rho = rho + z(i) * r(i)
         norm = max(norm,abs(r(i)))
      enddo
      
      iter = iter+1
c      write(6,*) iter,norm

      if (iter .gt. maxiter .or. norm .gt. max_jump*norm0) then
         
         print *, "cg solve in mac failed to converge"
         do i = istart, iend 
            p(i) = x(i) + dest0(i)
         enddo
         call setmacbc(DIMS(p),p,lo,hi,isPeriodic,setSingularPoint)
         
      else if (norm .lt. goal) then

         do i = istart, iend 
            p(i) = x(i) + dest0(i)
         enddo
         call setmacbc(DIMS(p),p,lo,hi,isPeriodic,setSingularPoint)
         
      else

        beta = rho / rho_old
        do i = istart, iend 
           p(i) = z(i) + beta * p(i)
        enddo

        goto 100

      endif

      return
      end

      subroutine makemacdgphi(phi,DIMS(phi),dgphi,DIMS(dgphi),
     &                     beta,DIMS(beta),
     &                     lo,hi,h,isPeriodic,setSingularPoint)
      implicit none

      integer DIMDEC(phi)
      integer DIMDEC(dgphi)
      integer DIMDEC(beta)
      REAL_T phi(DIM1(phi))
      REAL_T dgphi(DIM1(dgphi))
      REAL_T beta(DIM1(beta))
      integer lo(SDIM),hi(SDIM)
      REAL_T h(SDIM)
      integer isPeriodic(SDIM)
      logical setSingularPoint

      integer i,j
      REAL_T hxsqinv

      hxsqinv = one/(h(1)*h(1))

      do i = lo(1),hi(1) 
         dgphi(i) = DGX
      enddo

      end

c *************************************************************************
c ** COARSIGMA **
c ** Coarsen the edge-based sigma coefficients
c *************************************************************************

      subroutine FORT_COARSIGMA(sigma,DIMS(sigma),sigmac,DIMS(sigmac),
     &                          lo,hi,loc,hic)

      implicit none
      integer lo(SDIM),hi(SDIM)
      integer loc(SDIM),hic(SDIM)
      integer DIMDEC(sigma)
      integer DIMDEC(sigmac)
      REAL_T  sigma(DIM1(sigma))
      REAL_T sigmac(DIM1(sigmac))

c     Local variables
      integer i,twoi

      do i = loc(1),hic(1)+1
        twoi = 2*(i-loc(1))+lo(1)
        sigmac(i) = sigma(twoi)
      enddo

      return
      end


c *************************************************************************
c ** RESTRICT **
c ** Conservatively average the residual
c *************************************************************************

      subroutine FORT_RESTRICT(res,DIMS(res),resc,DIMS(resc),
     &                         lo,hi,loc,hic)

      implicit none
      integer lo(SDIM),hi(SDIM)
      integer loc(SDIM),hic(SDIM)
      integer DIMDEC(res)
      integer DIMDEC(resc)
      REAL_T  res(DIM1(res))
      REAL_T resc(DIM1(resc))

c     Local variables
      integer i,twoi

c ::: NOTE: dont need factor of r here for volume-weighting because
c ::: what were calling the residual is really already r*residual

        do i = loc(1),hic(1) 
          twoi = 2*(i-loc(1))+lo(1)
          resc(i) = half*(res(twoi) + res(twoi+1))
        enddo

      return
      end

c *************************************************************************
c ** INTERPOLATE **
c ** Piecewise constant interpolation
c *************************************************************************

      subroutine FORT_INTERPOLATE(phi,DIMS(phi),deltac,DIMS(deltac),
     &                            lo,hi,loc,hic)

      implicit none
      integer lo(SDIM),hi(SDIM)
      integer loc(SDIM),hic(SDIM)
      integer DIMDEC(phi)
      integer DIMDEC(deltac)
      REAL_T    phi(DIM1(phi))
      REAL_T deltac(DIM1(deltac))

c     Local variables
      integer i,twoi

      do i = loc(1), hic(1) 
         twoi = 2*(i-loc(1))+lo(1)
         phi(twoi  ) = phi(twoi  ) + deltac(i)
         phi(twoi+1) = phi(twoi+1) + deltac(i)
      enddo
      
      return
      end


c *************************************************************************
c ** MACPHIBC **
c *************************************************************************

      subroutine FORT_MACPHIBC(phi,length,divuExt,rhoExt,redge,hx,per)
c
c    Compute the value of phi for macproj to be used at an  outflow face,
c    assuming that the tangential velocity on the edges of the outflow boundary
c    are either zero or periodic.
c
      implicit none
      integer length
      integer r_len
      integer per
      REAL_T     phi(length)
      REAL_T divuExt(0:length-1)
      REAL_T  rhoExt(0:length-1)
      REAL_T   redge(0:length)
      REAL_T hx

c     Local variables
      integer NstripMAX
      parameter (NstripMAX = 2000)
      integer i, j, neq,n
      REAL_T a(NstripMAX), b(NstripMAX), c(NstripMAX), s(NstripMAX)
      REAL_T alpha, beta, sVal
      logical rNormed
      REAL_T vtot
      REAL_T rL, rR, diff
      REAL_T phitot
      REAL_T rcen
      integer ics, ice

c     This description assumes outflow at yhi; however, code works for 
c     outflow at any face.
c     NOTE: Assumes that rho at edge between i, i-1 = half*(rho(i)+rho(i-1))
c             (1) Linear fit of rho between nodes
c             (2) rho, divu on same boxes (box)
c             (3) phi is on box, shifted up one

c     Solve d/dx( 1/rho d/dx( phi ) ) = dU/dx - (S - S_ave) [S = divu] 
c     w/periodic or Neumann BC's, using a tridiagonal solve which detects, 
c     and deals with, the singular equations.  In the Neumann case, 
c     arbitrarily set the upper right corner to zero to pin the solution.  
c     Note that the RHS of this equation satisfies the solvability 
c     constraint that Int[RHS.dV] = 0 by construction.
c     This implies that the normal component takes up the slack:
c     
c                        d/dy( 1/rho d/dy( phi ) ) = dV/dy - S_ave
c     
c     This information should be used to construct the normal gradient of the
c     normal velocity, for the advective/diffusive step, for example.
c     In this implementation, use that d/dy == 0 at top, so y-edge centered 
c     values come directly from cell-centers just inside domain

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3

      ics = 0
      ice = length-1
      neq = length

      if ( neq .gt. NStripMax ) then
         call bl_error( 'MACPHIBC: NStripMax too small' )
      end if

c     Carry out non-trivial solve.  First set interior equations, then do BC's
      do n = 2,neq-1
         i = n + ics - 1
         a(n)=two*redge(i  )/(rhoExt(i)+rhoExt(i-1))
         c(n)=two*redge(i+1)/(rhoExt(i)+rhoExt(i+1))
         b(n)=- a(n) - c(n)
         s(n)= -divuExt(i)
      end do
      
      if (per .eq. 1) then
         
c     Do left-side periodic BC (keep r in there to guarantee correct scaling)
         i = ics
         beta=two*redge(i  )/(rhoExt(i)+rhoExt(ice))
         c(1)=two*redge(i+1)/(rhoExt(i)+rhoExt(i+1))
         b(1)=- beta - c(1)
         s(1)= -divuExt(i)
         
c     Do right-side periodic
         i = ice
         a(neq)=two*redge(i  )/(rhoExt(i)+rhoExt(i-1))
         alpha =two*redge(i+1)/(rhoExt(i)+rhoExt(ics))
         b(neq)=- a(neq) - alpha
         s(neq)= -divuExt(i)
         
c     Solve the equations
         call cyclic(a,b,c,alpha,beta,s,phi,neq)
         
      else
c     Solid walls, Neumann conditions (dphi/dx=u=0)
         i = ics
         c(1) = two*redge(i+1)/(rhoExt(i)+rhoExt(i+1))
         b(1) = - c(1)
         s(1)=  -divuExt(i)
         
         i = ice
         a(neq) = two*redge(i)/(rhoExt(i)+rhoExt(i-1))
         b(neq) = - a(neq)
         s(neq)= -divuExt(i)

c     Solve the equations (we know they're singular, pass the arbitrary value,
c     and a flag that we've already normalized the rhs, in the sense that
c     Int[dU/dx - (S-S_ave)] == 0
         sVal = zero
         rNormed = .true.
         call tridag_sing(a,b,c,s,phi,neq,sVal,rNormed)
         
      end if
      
c     Try normalizing phi to average to zero
      phitot = zero
      vtot = zero
      do n = 1,neq
         rcen = half*(redge(ics+n) + redge(ics+n-1))
         phitot = phitot + phi(n)*rcen
         vtot = vtot + rcen
      end do
      phitot = phitot / vtot
      do n = 1,neq
         phi(n) = phi(n) - phitot
      end do

#undef XLO
#undef YLO
#undef XHI
#undef YHI
      end

c *************************************************************************
c ** MACFILL_ONED **
c *************************************************************************

      subroutine FORT_MACFILL_ONED(lenx,leny,length,faces,numOutFlowFaces,
     $                             cc0,cc1,cc2,cc3,
     $                             r0, r1, r2, r3 , conn, redge_conn)

      integer lenx,leny,length
      integer faces(4)
      integer numOutFlowFaces
      REAL_T cc0(leny,2)
      REAL_T cc1(lenx,2)
      REAL_T cc2(leny,2)
      REAL_T cc3(lenx,2)
      REAL_T conn(length,2)
      REAL_T r0(leny+1)
      REAL_T r1(lenx+1)
      REAL_T r2(leny+1)
      REAL_T r3(lenx+1)
      REAL_T redge_conn(length+1)

      integer xlo_outflow,ylo_outflow
      integer xhi_outflow,yhi_outflow
      integer i,ifinal,n
      REAL_T sum

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3

c     Want to find the single non-outflow face.
      xlo_outflow = 0
      ylo_outflow = 0
      xhi_outflow = 0
      yhi_outflow = 0

      do i = 1, numOutFlowFaces
        if (faces(i) .eq. XLO) xlo_outflow = 1
        if (faces(i) .eq. YLO) ylo_outflow = 1
        if (faces(i) .eq. XHI) xhi_outflow = 1
        if (faces(i) .eq. YHI) yhi_outflow = 1
      enddo

c     Possible combinations of faces to come in here:
c       cc0 cc1 cc2 cc3
c       XLO YLO 
c       XLO         YHI 
c           YLO XHI 
c           YLO     YHI 
c       XLO YLO XHI
c       XLO     XHI YHI
c       XLO YLO     YHI
c           YLO XHI YHI
c       XLO YLO XHI YHI

c     We must remember here that the cc* arrays have already been
c       ordered so that the 2nd dimension is one cell wide.

      ifinal = 0
      if (numOutFlowFaces .eq. 4 .or. 
     $    (xlo_outflow .eq. 1 .and. ylo_outflow .eq. 0) ) then
          do i = 1,leny
            conn(i,1) = cc0(i,1)
            conn(i,2) = cc0(i,2)
          enddo
          do i = 1,leny+1
            redge_conn(i) = r0(i)
          enddo
          ifinal = leny
      endif

      if (yhi_outflow .eq. 1 .and. 
     $    .not. (numOutFlowFaces .eq. 3 .and. xhi_outflow .eq. 0) ) then
          do i = 1,lenx
            conn(ifinal+i,1) = cc3(i,1)
            conn(ifinal+i,2) = cc3(i,2)
          enddo
          do i = 1,lenx+1
            redge_conn(ifinal+i) = r3(i)
          enddo
          ifinal = ifinal + lenx
      endif

      if (xhi_outflow .eq. 1) then
          do i = 1,leny
            conn(ifinal+i,1) = cc2(leny+1-i,1)
            conn(ifinal+i,2) = cc2(leny+1-i,2)
          enddo
          do i = 1,leny+1
            redge_conn(ifinal+i) = r2(leny+2-i)
          enddo
          ifinal = ifinal + leny
      endif

      if (ylo_outflow .eq. 1) then
          do i = 1,lenx
            conn(ifinal+i,1) = cc1(lenx+1-i,1)
            conn(ifinal+i,2) = cc1(lenx+1-i,2)
          enddo
          do i = 1,lenx+1
            redge_conn(ifinal+i) = r1(lenx+2-i)
          enddo
          ifinal = ifinal + lenx
      endif

      if (numOutFlowFaces .lt. 4 .and.
     $    (xlo_outflow .eq. 1 .and. ylo_outflow .eq. 1) ) then
          do i = 1,leny
            conn(ifinal+i,1) = cc0(i,1)
            conn(ifinal+i,2) = cc0(i,2)
          enddo
          do i = 1,leny+1
            redge_conn(ifinal+i) = r0(i)
          enddo
          ifinal = ifinal + leny
      endif

      if (yhi_outflow .eq. 1 .and. 
     $    (numOutFlowFaces .eq. 3 .and. xhi_outflow .eq. 0) ) then
          do i = 1,lenx
            conn(ifinal+i,1) = cc3(i,1)
            conn(ifinal+i,2) = cc3(i,2)
          enddo
          do i = 1,lenx+1
            redge_conn(ifinal+i) = r3(i)
          enddo
          ifinal = ifinal + lenx
      endif

      length = ifinal
      sum = zero
      do i = 1,length
        sum = sum + conn(i,2)
      enddo

#undef XLO
#undef YLO
#undef XHI
#undef YHI

      end

c *************************************************************************
c ** MACPHI_FROM_X **
c *************************************************************************

      subroutine FORT_MACPHI_FROM_X(DIMS(phi),phi,length,x)

      integer DIMDEC(phi)
      integer length
      REAL_T phi(DIMV(phi))
      REAL_T x(length)

      integer lenx, leny
      integer i,j

c     We know that the faces are ordered: XLO,XHI,YLO,YHI
      lenx = ARG_H1(phi)-ARG_L1(phi)
      leny = ARG_H2(phi)-ARG_L2(phi)

      if (lenx .eq. 0) then
        do j = 1,length
          phi(ARG_L1(phi),j-1) = x(j)
        enddo
      elseif (leny .eq. 0) then
        do i = 1,length
          phi(i-1,ARG_L2(phi)) = x(i)
        enddo
      endif

#undef XLO
#undef YLO
#undef XHI
#undef YHI

      end 

c *************************************************************************
c ** MACALLPHI_FROM_X **
c *************************************************************************

      subroutine FORT_MACALLPHI_FROM_X(lenx,leny,length,faces,numOutFlowFaces,
     $                                 phi0,phi1,phi2,phi3,x)

      integer lenx,leny,length
      integer numOutFlowFaces
      integer faces(4)
      REAL_T phi0(0:leny-1)
      REAL_T phi1(0:lenx-1)
      REAL_T phi2(0:leny-1)
      REAL_T phi3(0:lenx-1)
      REAL_T x(0:length-1)

      integer xlo_outflow,ylo_outflow
      integer xhi_outflow,yhi_outflow
      integer i,j,ifinal,n

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3

c     Possible combinations of faces to come in here:
c       phi0 phi1 phi2 phi3
c       XLO  YLO 
c       XLO            YHI 
c            YLO  XHI 
c            YLO       YHI 
c       XLO  YLO  XHI
c       XLO       XHI  YHI
c       XLO  YLO       YHI
c            YLO  XHI  YHI
c       XLO  YLO  XHI  YHI

c     Want to find which are outflow faces.
      xlo_outflow = 0
      ylo_outflow = 0
      xhi_outflow = 0
      yhi_outflow = 0

      do i = 1, numOutFlowFaces
        if (faces(i) .eq. XLO) xlo_outflow = 1
        if (faces(i) .eq. YLO) ylo_outflow = 1
        if (faces(i) .eq. XHI) xhi_outflow = 1
        if (faces(i) .eq. YHI) yhi_outflow = 1
      enddo

c     We know that the faces are ordered: XLO,XHI,YLO,YHI
      
      ifinal = 0

      if (numOutFlowFaces .eq. 4 .or. 
     $    (xlo_outflow .eq. 1 .and. ylo_outflow .eq. 0) ) then
        do j = 0,leny-1
          phi0(j) = x(j)
        enddo
        ifinal = leny
      endif

      if (yhi_outflow .eq. 1 .and. 
     $    .not. (numOutFlowFaces .eq. 3 .and. xhi_outflow .eq. 0) ) then
        do i = 0,lenx-1
          phi3(i) = x(i+ifinal)
        enddo
        ifinal = ifinal+lenx
      endif

      if (xhi_outflow .eq. 1) then
        do j = 0,leny-1
          phi2(leny-1-j) = x(ifinal+j)
        enddo
        ifinal = ifinal+leny
      endif

      if (ylo_outflow .eq. 1) then
        if (numOutFlowFaces .eq. 4) then 
          do i = 0,lenx-1
            phi1(lenx-1-i) = x(ifinal+i)
          enddo
          phi1(0) = x(0)
        else
          do i = 0,lenx-1
            phi1(lenx-1-i) = x(ifinal+i)
          enddo
        endif
        ifinal = ifinal+lenx
      endif

      if (numOutFlowFaces .lt. 4 .and.
     $    (xlo_outflow .eq. 1 .and. ylo_outflow .eq. 1) ) then
        do j = 0,leny-1
          phi0(j) = x(j+ifinal)
        enddo
        ifinal = ifinal+leny
      endif

      if (yhi_outflow .eq. 1 .and. 
     $    (numOutFlowFaces .eq. 3 .and. xhi_outflow .eq. 0) ) then
        do i = 0,lenx-1
          phi3(i) = x(i+ifinal)
        enddo
        ifinal = ifinal+lenx
      endif

      end
