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


      subroutine FORT_EXTRAP_MAC(DIMS(u0),u0,DIMS(u1),u1,DIMS(div),divu,DIMS(rho),rho,
     &                         r_lo,r_hi,r,DIMS(uExt),uExt,DIMS(divuExt),divuExt,
     &                         DIMS(rhoExt),rhoExt,lo,hi,face,zeroIt)
c
c     Compute the value of phi for macproj 
c     assuming that the tangential velocity on the edges of the outflow boundary
c     are either zero or periodic.
c     note that u is edge centered

c subtract divu_ave twice due to precision problems
      implicit none

      integer DIMDEC(u0)
      integer DIMDEC(u1)
      integer DIMDEC(div)
      integer DIMDEC(rho)
      integer r_lo,r_hi
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(uExt)
      integer DIMDEC(divuExt)
      integer DIMDEC(rhoExt)
      REAL_T      u0(DIMV(u0))
      REAL_T      u1(DIMV(u1))
      REAL_T divu(DIMV(div))
      REAL_T    rho(DIMV(rho))
      REAL_T      r(r_lo:r_hi)
      REAL_T   uExt(DIMV(uExt))
      REAL_T   divuExt(DIMV(divuExt))
      REAL_T   rhoExt(DIMV(rhoExt))
      integer zeroIt
      integer face
      
c     Local variables
      REAL_T small_pert
      parameter ( small_pert = SMALL)
      integer i, j, neq,n
      REAL_T divu_ave1,divu_ave2
      REAL_T max_divu, min_divu, max_pert
      integer ics,ice,jcs,jce
      integer ifs,ife,jfs,jfe
      integer if,jf,ic,jc
c     NOTE: Assumes that rho at edge between i, i-1 = half*(rho(i)+rho(i-1))
c             (1) Linear fit of rho between nodes
c             (2) rho, divu on same boxes (box), and r(ARG_L1(rho),ARG_L2(rho))
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

      zeroIt = 0

      if (face .eq. XLO) then
         if = ifs
         do j = jcs, jce+1
            uExt(if,j)  = u1(ics,j)
         enddo
         max_divu = divu(ics,jcs)
         min_divu = max_divu
         do j = jcs, jce
            divuExt(if,j) = divu(ics,j)
            rhoExt(if,j)  = rho(ics,j)
            max_divu = max(max_divu,divuExt(if,j))
            min_divu = min(min_divu,divuExt(if,j))
         end do
         call subtractavg(DIMS(divuExt),divuExt,r,r_lo,r_hi,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,r,r_lo,r_hi,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(if,jcs))
         do j = jcs, jce
            max_pert = MAX(max_pert,ABS(divuExt(if,j)))
         end do
      else if (face .eq. YLO) then
         jf = jfs
         do i = ics, ice+1
            uExt(i,jf)  = u0(i,jcs)
         enddo
         max_divu = divu(ics,jcs)
         min_divu = max_divu
         do i = ics, ice
            divuExt(i,jf) = divu(i,jcs)
            rhoExt(i,jf)  = rho(i,jcs)
            max_divu = max(max_divu,divuExt(i,jf))
            min_divu = min(min_divu,divuExt(i,jf))
         end do
         call subtractavg(DIMS(divuExt),divuExt,r,r_lo,r_hi,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,r,r_lo,r_hi,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(ics,jf))
         do i = ics, ice
            max_pert = MAX(max_pert,ABS(divuExt(i,jf)))
         end do
      else if (face .eq. XHI) then
         if = ife
         do j = jcs, jce+1
            uExt(if,j)  = u1(ice,j)
         enddo
         max_divu = divu(ice,jcs)
         min_divu = max_divu
         do j = jcs, jce
            divuExt(if,j) = divu(ice,j)
            rhoExt(if,j)  = rho(ice,j)
            max_divu = max(max_divu,divuExt(if,j))
            min_divu = min(min_divu,divuExt(if,j))
         end do
         call subtractavg(DIMS(divuExt),divuExt,r,r_lo,r_hi,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,r,r_lo,r_hi,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(if,jcs))
         do j = jcs, jce
            max_pert = MAX(max_pert,ABS(divuExt(if,j)))
         end do
      else if (face .eq. YHI) then
         jf = jfe
         do i = ics, ice+1
            uExt(i,jf)  = u0(i,jce)
         enddo
         max_divu = divu(ics,jce)
         min_divu = max_divu
         do i = ics, ice
            divuExt(i,jf) = divu(i,jce)
            rhoExt(i,jf)  = rho(i,jce)
            max_divu = max(max_divu,divuExt(i,jf))
            min_divu = min(min_divu,divuExt(i,jf))
         end do
         call subtractavg(DIMS(divuExt),divuExt,r,r_lo,r_hi,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,r,r_lo,r_hi,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(ics,jf))
         do i = ics, ice
            max_pert = MAX(max_pert,ABS(divuExt(i,jf)))
         end do
      endif
      
c  check to see if we should zero phi
         max_pert = max_pert/(ABS(divu_ave1+divu_ave2)+small_pert)
      if ((max_divu.eq.zero.and.min_divu.eq.zero)
     &     .or.(max_pert.le.small_pert)) then
         zeroIt = 1
      end if
      end

      subroutine subtractavg(DIMS(divu),divu,r,r_lo,r_hi,lo,hi,divu_ave,face)
      implicit none
      integer DIMDEC(divu)
      integer lo(SDIM),hi(SDIM)
      integer r_lo,r_hi
      REAL_T r(r_lo:r_hi)
      REAL_T divu_ave
      REAL_T divu(DIMV(divu))
      integer face

      integer i,j
      REAL_T vtot

      divu_ave = zero
      vtot = zero

      if (face .eq. XLO .or. face .eq. XHI) then
         i = lo(1)
         do j=lo(2),hi(2)
            vtot = vtot+one
            divu_ave = divu_ave+divu(i,j)
         enddo
         divu_ave = divu_ave/vtot
         do j=lo(2),hi(2)
            divu(i,j) = divu(i,j) - divu_ave
         enddo
      elseif (face .eq. YLO .or. face .eq. YHI) then
         j = lo(2)
         do i=lo(1),hi(1)
            vtot = vtot+r(i)
            divu_ave = divu_ave+r(i)*divu(i,j)
         enddo
         divu_ave = divu_ave/vtot
         do i=lo(1),hi(1)
            divu(i,j) = divu(i,j) - divu_ave
         enddo
      else 
         print*, "bad length"
      endif

      end
#undef XLO
#undef YLO
#undef XHI
#undef YHI

      subroutine FORT_COMPUTE_MACCOEFF(DIMS(rhs),rhs,DIMS(beta),beta,
     &                              DIMS(uExt),uExt,DIMS(divuExt),divuExt,
     &                              DIMS(rhoExt),rhoExt,r_lo,r_hi,r,
     &                              lo,hi,h,isPeriodic)
      implicit none
      integer DIMDEC(rhs)
      integer DIMDEC(beta)
      integer DIMDEC(uExt)
      integer DIMDEC(divuExt)
      integer DIMDEC(rhoExt)
      integer lo(SDIM),hi(SDIM)
      integer r_lo, r_hi
      REAL_T uExt(DIM1(uExt))
      REAL_T divuExt(DIM1(divuExt))
      REAL_T rhoExt(DIM1(rhoExt))
      REAL_T beta(DIM1(beta))
      REAL_T rhs(DIM1(rhs))
      REAL_T r(r_lo:r_hi)
      REAL_T h
      integer isPeriodic

      REAL_T a(0:100),b(0:100),c(0:100),s(0:100)
      integer i
      integer ins,ine,ics,ice
      REAL_T rR, rL
      ins = lo(1)
      ine = hi(1)+1
      ics = lo(1)
      ice = hi(1)

      do i=ics+1,ice-1
         rL = half*(r(i)+r(i-1))
         rR = half*(r(i)+r(i+1))
         rhs(i)=(rR*uExt(i+1)-rL*uExt(i))/h - r(i)*(divuExt(i))
      enddo

      do i=ins+1,ine-1
         rL = half*(r(i)+r(i-1))
         beta(i)=two*rL/(rhoExt(i)+rhoExt(i-1))
      enddo

      if (isPeriodic) then
         rL = half*(r(ics)+r(ice))
         rR = half*(r(ics)+r(ics+1))
         rhs(ics)=(rR*uExt(ics+1)-rL*uExt(ics))/h - r(ics)*(divuExt(ics))
         beta(ins) = two*rL/(rhoExt(ics)+rhoExt(ice)) 

         rL = half*(r(ice)+r(ice-1))
         rR = half*(r(ics)+r(ice))
         rhs(ice)=(rR*uExt(ice+1)-rL*uExt(ice))/h - r(ice)*(divuExt(ice))
         beta(ine) =  two*rR/(rhoExt(ics)+rhoExt(ice))
      else
         rR = half*(r(ics)+r(ics+1))
         rhs(ics)= rR*uExt(ics+1)/h - r(ics)*(divuExt(ics))
         beta(ins) = zero

         rL = half*(r(ice)+r(ice-1))
         rhs(ice)=-rL*uExt(ice)/h - r(ice)*(divuExt(ice))
         beta(ine) = zero
      endif


c      note:
c         a(i) = beta(i)
c         c(i) = beta(i+1)
c         b(i) = -a(i)-c(i)
c         s(i) = rhs(i)*h*h
      end 

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

      hxsqinv = 1.0/(h(1)*h(1))
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

      subroutine setmacbc(DIMS(phi),phi,lo,hi,isPeriodic,setSingularPoint)
      implicit none
      integer DIMDEC(phi)
      REAL_T phi(DIM1(phi))
      integer lo(SDIM),hi(SDIM)
      integer isPeriodic
      integer ics,ice
      logical setSingularPoint
      ics = lo(1)
      ice = hi(1)

      if (isPeriodic .NE. 1 .and. setSingularPoint) then
         phi(ice)= zero
      endif

      if (isPeriodic) then
         phi(ics-1) = phi(ice)
         phi(ice+1) = phi(ics)
      else
         phi(ics-1) = phi(ics)
         phi(ice+1) = phi(ice)
      endif

      end


      subroutine FORT_MACTRANSLATE(DIMS(out),out,DIMS(in),in,face)
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


      subroutine FORT_MACPHIBC(DIMS(uExt),uExt,DIMS(divuExt),divuExt,
     &                         DIMS(rhoExt),rhoExt,
     &                         r_lo,r_hi,r,hx,DIMS(phi),phi,lo,hi,face,per)
c
c     Compute the value of phi for macproj to be used at an  upper, y-outflow face,
c     assuming that the tangential velocity on the edges of the outflow boundary
c     are either zero or periodic.
c
      implicit none
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(phi)
      integer DIMDEC(rhoExt)
      integer DIMDEC(divuExt)
      integer DIMDEC(uExt)
      integer r_lo,r_hi
      integer per
      integer face
      REAL_T    phi(DIMV(phi))
      REAL_T      r(r_lo:r_hi)
      REAL_T divuExt(DIM1(divuExt))
      REAL_T    uExt(DIM1(uExt))
      REAL_T  rhoExt(DIM1(rhoExt))
      REAL_T hx

c     Local variables
      integer NstripMAX
      parameter (NstripMAX = 1000)
      REAL_T small_udiff
      parameter (small_udiff = SMALL)
      integer i, j, neq,n
      REAL_T a(NstripMAX), b(NstripMAX), c(NstripMAX), s(NstripMAX)
      REAL_T x(NstripMAX), alpha, beta, sVal
      logical rNormed
      REAL_T vtot
      REAL_T rL, rR, diff
      REAL_T phitot
      integer ics, ice

c     NOTE: Assumes that rho at edge between i, i-1 = half*(rho(i)+rho(i-1))
c             (1) Linear fit of rho between nodes
c             (2) rho, divu on same boxes (box), and r(ARG_L1(rho),ARG_L2(rho))
c             (3) phi is on box, shifted up one
c             (4) u is edge-based, on surroundingNodes(box)

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

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3
      
      ics = lo(1)
      ice = hi(1)
      neq = ice - ics + 1

c     Carry out non-trivial solve.  First set interior equations, then do BC's
      do n = 2,neq-1
         i = n + ics - 1
         rL = half*(r(i)+r(i-1))
         rR = half*(r(i)+r(i+1))
         a(n)=two*rL/(rhoExt(i)+rhoExt(i-1))
         c(n)=two*rR/(rhoExt(i)+rhoExt(i+1))
         b(n)=- a(n) - c(n)
         s(n)=(rR*uExt(i+1)-rL*uExt(i))*hx - r(i)*(divuExt(i))*hx*hx
      end do
      
      if (per .eq. 1) then
         
c     Make sure u_mac is periodic
         diff = uExt(ics) - uExt(ice+1)
         if (ABS(diff) .gt. small_udiff) then
            write(6,*) 'MACPHIBC: periodic umac screwy.  diff =',diff
            stop
         endif
         
c     Do left-side periodic BC (keep r in there to guarantee correct scaling)
         i = ics
         rL = half*(r(i)+r(ice))
         rR = half*(r(i)+r(i+1))
         beta=two*rL/(rhoExt(i)+rhoExt(ice))
         c(1)=two*rR/(rhoExt(i)+rhoExt(i+1))
         b(1)=- beta - c(1)
         s(1)=(rR*uExt(i+1)-rL*uExt(i))*hx - r(i)*(divuExt(i))*hx*hx
         
c     Do right-side periodic
         i = ice
         rL = half*(r(i)+r(i-1))
         rR = half*(r(ics)+r(i))
         a(neq)=two*rL/(rhoExt(i)+rhoExt(i-1))
         alpha=two*rR/(rhoExt(i)+rhoExt(ics))
         b(neq)=- a(neq) - alpha
         s(neq)=(rR*uExt(i+1)-rL*uExt(i))*hx - r(i)*(divuExt(i))*hx*hx
         
c     Solve the equations
         call cyclic(a,b,c,alpha,beta,s,x,neq)
         
      else
c     Solid walls, Nuemann conditions (dphi/dx=u=0)
         i = ics
         rR = half*(r(i)+r(i+1))
         c(1) = two*rR/(rhoExt(i)+rhoExt(i+1))
         b(1) = - c(1)
         s(1)= rR*uExt(i+1)*hx - r(i)*(divuExt(i))*hx*hx
         
         i = ice
         rL = half*(r(i)+r(i-1))
         a(neq) = two*rL/(rhoExt(i)+rhoExt(i-1))
         b(neq) = - a(neq)
         s(neq)=-rL*uExt(i)*hx - r(i)*(divuExt(i))*hx*hx
         
c     Solve the equations (we know they're singular, pass the arbitrary value, and
c     a flag that we've already normalized the rhs, in the sense that
c     Int[dU/dx - (S-S_ave)] == 0
         sVal = zero
         rNormed = .true.
         call tridag_sing(a,b,c,s,x,neq,sVal,rNormed)
         
      end if
      
c     Try normalizing phi to average to zero
      phitot = zero
      vtot = zero
      do n = 1,neq
         phitot = phitot + x(n)*r(ics+n-1)
         vtot = vtot + r(ics+n-1)
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

