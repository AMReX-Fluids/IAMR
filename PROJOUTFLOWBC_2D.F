#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "PROJOUTFLOWBC_F.H"
#include "ArrayLim.H"

#define SDIM 2

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define SMALL 1.0e-10
#else
#define SMALL 1.0d-10
#endif

      subroutine FORT_EXTRAP_PROJ(DIMS(u),u,DIMS(divu),divu,DIMS(rho),rho,
     &     r_lo,r_hi,r,DIMS(uExt),uExt,DIMS(divuExt),divuExt,
     &     DIMS(rhoExt),rhoExt,lo,hi,face, zeroIt)
      implicit none

c subtract divu_ave twice due to precision problems

      integer DIMDEC(u)
      integer DIMDEC(divu)
      integer DIMDEC(rho)
      integer DIMDEC(uExt)
      integer DIMDEC(divuExt)
      integer DIMDEC(rhoExt)
      integer r_lo,r_hi
      integer lo(SDIM),hi(SDIM)
      integer face
      REAL_T      u(DIMV(u),SDIM)
      REAL_T   divu(DIMV(divu))
      REAL_T    rho(DIMV(rho))
      REAL_T   uExt(DIMV(uExt),SDIM-1)
      REAL_T   divuExt(DIMV(divuExt))
      REAL_T   rhoExt(DIMV(rhoExt))
      REAL_T   r(r_lo:r_hi)
      integer  zeroIt

c local variables
      integer ics,ice,jcs,jce
      integer ifs,ife,jfs,jfe
      integer if,jf
      REAL_T divu_ave1,divu_ave2
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

      ifs = lo(1)
      ife = hi(1)
      jfs = lo(2)
      jfe = hi(2)

      zeroIt = 0

      if (face .eq. XLO) then
         if=ife
         max_divu = half*(three*divu(ice-1,jcs) - divu(ice,jcs))
         min_divu = max_divu
         do j = jcs, jce
            uExt(j,if,1)    = half*(three*u(ice-1,j,2)    - u(ice,j,2))
            divuExt(j,if) = half*(three*divu(ice-1,j) - divu(ice,j))
            rhoExt(j,if)  = half*(three*rho(ice-1,j)    - rho(ice,j))
            max_divu = max(max_divu,divuExt(j,if))
            min_divu = min(min_divu,divuExt(j,if))
         end do
         call subtractavg(DIMS(divuExt),divuExt,r,r_lo,r_hi,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,r,r_lo,r_hi,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(jcs,if))
         do j = jcs, jce
            max_pert = MAX(max_pert,ABS(divuExt(j,if)))
         end do
      else if (face .eq. YLO) then
         jf = jfe
         max_divu = half*(three*divu(ics,jce-1) - divu(ics,jce))
         min_divu = max_divu
         do i = ics, ice
            uExt(i,jf,1)    = half*(three*u(i,jce-1,1)    - u(i,jce,1))
            divuExt(i,jf) = half*(three*divu(i,jce-1) - divu(i,jce))
            rhoExt(i,jf)  = half*(three*rho(i,jce-1)    - rho(i,jce))
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
         max_divu = half*(three*divu(ics+1,jcs) - divu(ics,jcs))
         min_divu = max_divu
         do j = jcs, jce
            uExt(j,if,1)    = half*(three*u(ics+1,j,2)    - u(ics,j,2))
            divuExt(j,if) = half*(three*divu(ics+1,j) - divu(ics,j))
            rhoExt(j,if)  = half*(three*rho(ics+1,j)    - rho(ics,j))
            max_divu = max(max_divu,divuExt(j,if))
            min_divu = min(min_divu,divuExt(j,if))
         end do
         call subtractavg(DIMS(divuExt),divuExt,r,r_lo,r_hi,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,r,r_lo,r_hi,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(jcs,if))
         do j = jcs, jce
            max_pert = MAX(max_pert,ABS(divuExt(j,if)))
         end do
      else if (face .eq. YHI) then
         jf = jfe
         max_divu = half*(three*divu(ics,jcs+1) - divu(ics,jcs))
         min_divu = max_divu
         do i = ics, ice
            uExt(i,jf,1)    = half*(three*u(i,jcs+1,1)    - u(i,jcs,1))
            divuExt(i,jf) = half*(three*divu(i,jcs+1) - divu(i,jcs))
            rhoExt(i,jf)  = half*(three*rho(i,jcs+1)    - rho(i,jcs))
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

c check to see if we should zero phi
      max_pert = max_pert/(ABS(divu_ave1+divu_ave2)+small_pert)
      if ((max_divu.eq.zero.and.min_divu.eq.zero)
     &     .or.(max_pert.le.small_pert)) then
         zeroIt = 1
      end if
#undef XLO
#undef YLO
#undef XHI
#undef YHI

      end

      subroutine FORT_COMPUTE_COEFF(DIMS(rhs),rhs,DIMS(beta),beta,
     &                              DIMS(uExt),uExt,DIMS(divuExt),divuExt,
     &                              DIMS(rhoExt),rhoExt,r_lo,r_hi,r,
     &                              lo,hi,h,isPeriodic)
      implicit none
      integer DIMDEC(rhs)
      integer DIMDEC(beta)
      integer DIMDEC(uExt)
      integer DIMDEC(divuExt)
      integer DIMDEC(rhoExt)
      integer r_lo, r_hi
      REAL_T uExt(DIM1(uExt),SDIM-1)
      REAL_T divuExt(DIM1(divuExt))
      REAL_T rhoExt(DIM1(rhoExt))
      REAL_T beta(DIM1(beta))
      REAL_T rhs(DIM1(rhs))
      REAL_T r(r_lo:r_hi)
      REAL_T h(SDIM)
      integer isPeriodic(SDIM)
      integer lo(SDIM),hi(SDIM)
      
      integer i
      integer ics,ice
      integer ins,ine
      REAL_T rnode,hdr,hxm1

      ics = lo(1)
      ice = hi(1)

      ins = lo(1)
      ine = hi(1)+1

      hxm1 = one/h(1)

      do i=ics,ice
        beta(i) = r(i)/rhoExt(i)
      enddo
      do i=ins+1,ine-1
         rnode = half*(r(i)+r(i-1))
         rhs(i) = (r(i)*uExt(i,1)-r(i-1)*uExt(i-1,1))*hxm1
     &        -  rnode*half*(divuExt(i-1)+divuExt(i))
      enddo

      if (isPeriodic(1).eq.1) then
          rnode = half*(r(ics)+r(ice))

          beta(ics-1) = beta(ice)
          rhs(ins)    =  (r(ics)*uExt(ics,1)-r(ice)*uExt(ice,1))*hxm1
     &        -  rnode*half*(divuExt(ice)+divuExt(ics))

          beta(ice+1) = beta(ics)
          rhs(ine)    =  (r(ics)*uExt(ics,1)-r(ice)*uExt(ice,1))*hxm1
     &        -  rnode*half*(divuExt(ice)+divuExt(ics))
      else
         hdr = half*(r(ics+1)-r(ics)) 
         rnode = r(ics)-hdr
         beta(ics-1) = zero
         rhs(ins) =   (r(ics)*uExt(ics,1))*hxm1
     &        -  rnode*half*(divuExt(ics))

         rnode = r(ice)+hdr
         beta(ice+1) = zero
         rhs(ine) =  (-r(ice)*uExt(ice,1))*hxm1
     &        -  rnode*half*(divuExt(ice))

       endif

c double rhs at edges
       if (isPeriodic(1) .ne. 1) then
         rhs(ins) = rhs(ins) * two
         rhs(ine) = rhs(ine) * two
       endif

c      do i=ins,ine
c         a(i) = beta(i-1)
c         c(i) = beta(i)
c         b(i) = -a(i)-c(i)
c         s(i) = rhs(i)*h(1)*h(1)
c      end do

       end 

#define DGX (beta(i-1)*phi(i-1) - (beta(i)+beta(i-1))*phi(i) \
            +beta(i)*phi(i+1))*(hxsqinv)

      subroutine FORT_HGRELAX(DIMS(rhs),rhs,DIMS(beta),beta,DIMS(phi),phi,
     &                        DIMS(dgphi),dgphi,lo,hi,h,isPeriodic,niter)
      implicit none
      integer DIMDEC(beta)
      integer DIMDEC(rhs)
      integer DIMDEC(phi)
      integer DIMDEC(dgphi)
      REAL_T beta(DIM1(beta))
      REAL_T rhs(DIM1(rhs))
      REAL_T phi(DIM1(phi))
      REAL_T dgphi(DIM1(dgphi))
      integer isPeriodic(SDIM)
      REAL_T h(SDIM)
      integer lo(SDIM),hi(SDIM)
      integer niter

      integer ins, ine
      integer i,iter
      REAL_T lam
      logical setSingularPoint
      REAL_T hxsqinv
      integer redblack

      hxsqinv = 1.0D0/(h(1)*h(1))

      ins = lo(1)
      ine = hi(1)+1
      setSingularPoint = .false.

      do iter = 1,niter
         call setprojbc(DIMS(phi),phi,lo,hi,isPeriodic,setSingularPoint)
         do redblack = 0, 1
            do i=ins+redblack,ine,2
               dgphi(i) = DGX
               lam = hxsqinv*(beta(i)+beta(i-1))
c double dgphi at edges
               if ((i .eq. lo(1) .or. i .eq. hi(1)+1)
     &              .and. isPeriodic(1) .ne. 1) then
                  dgphi(i) = dgphi(i) * two
                  lam = lam* two
               endif
               lam = -one/lam
               phi(i) = phi(i) + lam*(rhs(i)-dgphi(i))
            enddo
         end do
      end do

      call setprojbc(DIMS(phi),phi,lo,hi,isPeriodic,setSingularPoint)

      end

      subroutine FORT_HGSUBTRACTAVGPHI(DIMS(phi),phi,r_lo,r_hi,r,lo,hi,
     &                                 isPeriodic)
      implicit none
      integer DIMDEC(phi)
      REAL_T phi(DIM1(phi))
      integer r_lo,r_hi
      REAL_T r(r_lo:r_hi)
      integer lo(SDIM),hi(SDIM)
      integer isPeriodic(SDIM)

      REAL_T phitot,vtot
      REAL_T hdr,rnode
      integer ics,ins,ine
      integer i
      logical setSingularPoint

      ics = lo(1)
      ins = lo(1)
      ine = hi(1)+1
      hdr = half *(r(ics+1)-r(ics))
      setSingularPoint = .false.

      phitot = zero
      vtot = zero

      do i = ins,ine-ins-2
         rnode = r(i+ins-ics) - hdr
         phitot = phitot + phi(i) * rnode
         vtot   = vtot + rnode
      enddo
      do i = ine-ins-1,ine
         rnode = r(i+ins-ics-1) + hdr
         phitot = phitot + phi(i) * rnode
         vtot   = vtot + rnode
      enddo
      phitot = phitot/vtot
      do i = ins,ine
         phi(i) = phi(i) - phitot
      enddo

      call setprojbc(DIMS(phi),phi,lo,hi,isPeriodic,setSingularPoint)

      end


      subroutine FORT_HGRESID(DIMS(rhs),rhs,DIMS(beta),beta,DIMS(phi),phi,
     &                   DIMS(resid),resid,DIMS(dgphi),dgphi,
     &                   lo,hi,h,isPeriodic,maxnorm)

      integer DIMDEC(beta)
      integer DIMDEC(rhs)
      integer DIMDEC(phi)
      integer DIMDEC(resid)
      integer DIMDEC(dgphi)
      REAL_T beta(DIM1(beta))
      REAL_T rhs(DIM1(rhs))
      REAL_T phi(DIM1(phi))
      REAL_T resid(DIM1(resid))
      REAL_T dgphi(DIM1(dgphi))
      integer isPeriodic(SDIM)
      REAL_T h(SDIM)
      integer lo(SDIM),hi(SDIM)
      REAL_T hxsqinv,hysqinv

      integer i,j
      REAL_T maxnorm
      logical setSingularPoint
 
      hxsqinv = one/(h(1)*h(1))
      setSingularPoint = .false.
      maxnorm = zero

      call makeprojdgphi(phi,DIMS(phi),dgphi,DIMS(dgphi),beta,DIMS(beta),
     &                   lo,hi,h,isPeriodic,setSingularPoint)
      do i=lo(1),hi(1)+1
         resid(i) = rhs(i)-dgphi(i)
         maxnorm = max(maxnorm,ABS(resid(i)))         
      enddo
      end


      subroutine setprojbc(DIMS(phi),phi,lo,hi,isPeriodic,setSingularPoint)
      implicit none
      integer DIMDEC(phi)
      REAL_T phi(DIM1(phi))
      integer lo(SDIM),hi(SDIM)
      integer isPeriodic(SDIM)
      logical setSingularPoint

      integer ins,ine
      
      ins = lo(1)
      ine = hi(1)+1

      if (isPeriodic(1) .NE. 1 .AND. setSingularPoint) then
         phi(ine) = zero
      endif

      if (isPeriodic(1).eq.1) then
         phi(ins-1) = phi(ine-1)
         phi(ine)   = phi(ins)
         phi(ine+1) = phi(ins+1)
      else
         phi(ins-1) = phi(ins+1)
         phi(ine+1) = phi(ine-1)
      endif

      end

      subroutine FORT_HG_SHIFT_PHI(DIMS(out),out,DIMS(in),in,face)
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

      subroutine FORT_HG_RESHIFT_PHI(DIMS(out),out,DIMS(in),in,face)
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

      subroutine FORT_SOLVEHG(p,DIMS(p),dest0,DIMS(dest0),
     &     source,DIMS(source),sigma,DIMS(sigma),
     &     cen,DIMS(cen),r,DIMS(r),w,DIMS(w),
     &     z,DIMS(z),x,DIMS(x),
     $     lo,hi,h,isPeriodic,maxiter,tol,abs_tol,max_jump,norm)
      
      implicit none

      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(p)
      integer DIMDEC(dest0)
      integer DIMDEC(source)
      integer DIMDEC(sigma)
      integer DIMDEC(r)
      integer DIMDEC(w)
      integer DIMDEC(z)
      integer DIMDEC(x)
      integer DIMDEC(cen)
      REAL_T   p(DIM1(p))
      REAL_T  dest0(DIM1(dest0))
      REAL_T source(DIM1(source))
      REAL_T  sigma(DIM1(sigma))
      REAL_T      r(DIM1(r))
      REAL_T      w(DIM1(w))
      REAL_T      z(DIM1(z))
      REAL_T      x(DIM1(x))
      REAL_T      cen(DIM1(cen))
      REAL_T h(SDIM)
      integer isPeriodic(SDIM)
      integer maxiter
      REAL_T norm
      REAL_T tol
      REAL_T abs_tol,max_jump

c     Local variables
      REAL_T factor
      REAL_T  alpha,beta, rho, rho_old
      logical testx
      integer i,iter
      integer istart,iend
      REAL_T norm0,goal
      logical setSingularPoint

      istart = lo(1)
      iend = hi(1)+1
      
      setSingularPoint = .false.

      do i = lo(1)-1,hi(1)+2
         dest0(i) = p(i)
         p(i) = zero
      enddo
      
      do i=ARG_L1(w),ARG_H1(w)
         w(i) = zero
      enddo
      
      call setprojbc(DIMS(dest0),dest0,lo,hi,isPeriodic,setSingularPoint)
      call makeprojdgphi(dest0,DIMS(dest0),w,DIMS(w),sigma,DIMS(sigma),
     &     lo,hi,h,isPeriodic,setSingularPoint)
      
      do i = istart, iend 
         r(i) = source(i) - w(i)
      enddo
      
c note that all of this factor stuff is due to the problem being doubled
c at edges -- both the rhs and the operator.
      rho = zero
      norm0 = zero
      do i = istart, iend
         factor = one
         testx  = (i .eq. lo(1) .or. i .eq. hi(1)+1)
         factor = cvmgt(factor*half,factor,testx)
         z(i) = r(i) 
         rho    = rho + factor * z(i) * r(i) 
         norm0 = max(norm0,abs(r(i)))
      enddo
      norm = norm0
      
      goal = max(tol * norm0,abs_tol)
      
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
      
      iter  = 0
      
 100  continue  
      
      do i=ARG_L1(w),ARG_H1(w)
         w(i) = zero
      enddo
      
      call setprojbc(DIMS(p),p,lo,hi,isPeriodic,setSingularPoint)
      call makeprojdgphi(p,DIMS(p),w,DIMS(w),sigma,DIMS(sigma),
     &     lo,hi,h,isPeriodic,setSingularPoint)
      
      alpha = zero
      do i = istart, iend 
         factor = one
         testx  = (i .eq. lo(1) .or. i .eq. hi(1)+1)
         factor = cvmgt(factor*half,factor,testx)
         alpha  = alpha + factor*p(i)*w(i) 
      enddo
      
      if (alpha .ne. zero) then
         alpha = rho / alpha
      else
         print *, "divide by zero"
         goto 101
      endif

      rho_old  = rho
      rho   = zero
      norm = zero
      do i = istart, iend 
         factor = one
         testx  = (i .eq. lo(1) .or. i .eq. hi(1)+1)
         factor = cvmgt(factor*half,factor,testx)
         x(i) = x(i) + alpha * p(i)
         r(i) = r(i) - alpha * w(i)
         z(i) = r(i) 
         rho    = rho + factor * z(i) * r(i) 
         norm = max(norm,abs(r(i)))
      enddo
      
      iter = iter+1
c      write(6,*) iter,norm

      if (iter .gt. maxiter .or. norm .gt. max_jump*norm0) then
         
         goto 101

      else if (norm .lt. goal) then
         
         do i = istart, iend 
            p(i) = x(i) + dest0(i)
         enddo
         
         return

      else
         
         beta = rho / rho_old
         do i = istart, iend 
            p(i) = z(i) + beta * p(i)
         enddo
         
        goto 100
        
      endif
      
 101  print *, "cg solve in proj failed to coverge"
      do i = istart, iend 
         p(i) = dest0(i)
      enddo

      return
      end

      subroutine makeprojdgphi(phi,DIMS(phi),dgphi,DIMS(dgphi),
     &                     beta,DIMS(beta),
     &                     lo,hi,h,isPeriodic,setSingularPoint)

      implicit none

      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(phi)
      integer DIMDEC(dgphi)
      integer DIMDEC(beta)
      REAL_T    phi(DIM1(phi))
      REAL_T  dgphi(DIM1(dgphi))
      REAL_T  beta(DIM1(beta))
      integer isPeriodic(SDIM)
      REAL_T h(SDIM)
      logical setSingularPoint

c     Local variables
      REAL_T hxsqinv
      integer is,ie
      integer i

      is = lo(1)
      ie = hi(1)

      hxsqinv = one/(h(1)*h(1))
      
      if (isPeriodic(1).eq.1 ) then
          phi(ie+2) = phi(is+1)
          phi(is-1) = phi(ie  )
      endif

      do i = is,ie+1
         dgphi(i) = DGX
      enddo
      
c  double dgphi at edges
      if (isPeriodic(1) .ne. 1) then
        dgphi(is) = dgphi(is) * two
        dgphi(ie+1) = dgphi(ie+1) * two
      endif

      if (setSingularPoint .and. isPeriodic(1) .NE. 1) then
         dgphi(hi(1)+1) = zero
      endif

      return
      end

c *************************************************************************
c ** COARSIG **
c ** Coarsening of the sig coefficients
c *************************************************************************

      subroutine FORT_COARSIG(sigma,DIMS(sigma),sigmac,DIMS(sigmac),
     &                        lo,hi,loc,hic,isPeriodic)

      implicit none

      integer lo(SDIM),hi(SDIM)
      integer loc(SDIM),hic(SDIM)
      integer DIMDEC(sigma)
      integer DIMDEC(sigmac)
      REAL_T  sigma(DIM1(sigma))
      REAL_T sigmac(DIM1(sigmac))
      integer isPeriodic(SDIM)

c     Local variables
      integer i,i2 

      do i = loc(1),hic(1) 
         i2 = 2*(i-loc(1))+lo(1)
         sigmac(i) = half*(sigma(i2) +sigma(i2+1))
      enddo

      if (isPeriodic(1).eq.1) then
         sigmac(loc(1)-1) = sigmac(hic(1))
         sigmac(hic(1)+1) = sigmac(loc(1))
      else
         sigmac(loc(1)-1) = zero
         sigmac(hic(1)+1) = zero
      endif

      return
      end

c *************************************************************************
c ** RESTRICT **
c ** Conservative restriction of the residual
c *************************************************************************

      subroutine FORT_RESTRICT(res,DIMS(res),resc,DIMS(resc),
     &                        lo,hi,loc,hic,isPeriodic)

      implicit none

      integer lo(SDIM),hi(SDIM)
      integer loc(SDIM),hic(SDIM)
      integer DIMDEC(res)
      integer DIMDEC(resc)
      REAL_T   res(DIM1(res))
      REAL_T  resc(DIM1(resc))
      integer isPeriodic(SDIM)

c     Local variables
      integer i,ii
      integer istart,iend

      istart = loc(1)
      iend = hic(1)+1
      
      if (isPeriodic(1).eq.1) then
          res(hi(1)+1) = res(lo(1)  )
          res(hi(1)+2) = res(lo(1)+1)
          res(lo(1)-1) = res(hi(1)  )
      endif

      do i = istart,iend
         ii = 2*(i-loc(1))+lo(1)
         resc(i) = half*res(ii) + fourth*(res(ii+1)+res(ii-1))
      enddo

c  the top version is what we use when we double the problem at edges
c  the bottom version (commented out) is what we would use if we did not 
c      double the problem at edges.
      if (isPeriodic(1) .NE. 1) then
        i = loc(1)
        ii = 2*(i-loc(1))+lo(1)
        
        resc(i) = half*(res(ii) + res(ii+1))
c        resc(i) = half*res(ii) + fourth*res(ii+1)
        
        i = hic(1)+1
        ii = 2*(i-loc(1))+lo(1)
        
        resc(i) = half*(res(ii) + res(ii-1))
c        resc(i) = half*res(ii) + fourth*res(ii-1)

      endif

      return
      end


c *************************************************************************
c ** INTERP **
c ** Simple bilinear interpolation
c *************************************************************************

      subroutine FORT_INTERP(phi,DIMS(phi),temp,DIMS(temp),deltac,DIMS(deltac),
     &                       sigma,DIMS(sigma),lo,hi,loc,hic,isPeriodic)

      implicit none

      integer lo(SDIM),hi(SDIM)
      integer loc(SDIM),hic(SDIM)
      integer DIMDEC(phi)
      integer DIMDEC(deltac)
      integer DIMDEC(sigma)
      integer DIMDEC(temp)
      REAL_T     phi(DIM1(phi))
      REAL_T  deltac(DIM1(deltac))
      REAL_T  sigma(DIM1(sigma))
      REAL_T   temp(DIM1(temp))
      integer isPeriodic(SDIM)

c     Local variables
      integer ii,ic
      integer is,ie,isc,iec
      logical setSingularPoint

      is = lo(1)
      ie = hi(1)

      isc = loc(1)
      iec = hic(1)

      do ic = isc, iec+1
         ii = 2*(ic-isc)+is
         temp(ii) = deltac(ic)
      enddo
      do ic = isc, iec 
         ii = 2*(ic-isc)+is
         temp(ii+1) = half*(deltac(ic) + deltac(ic+1)) 
      enddo
      
      do ii = is,ie+1
         phi(ii) = phi(ii) + temp(ii)
      enddo

      return
      end

      subroutine FORT_HGPHIBC(hx,r,uExt, divuExt,rhoExt,length,
     &                        DIMS(phi),phi,lo,hi,face,per)
c
c    Compute the value of phi for hgproj to be used at an outflow face,
c    assuming that the tangential velocity on the edges of the outflow boundary
c    are either zero or periodic.  
c
      implicit none

      integer DIMDEC(phi)
      integer length
      integer per
      integer face
      integer lo(SDIM),hi(SDIM)
      REAL_T    phi(DIMV(phi))
      REAL_T      r(0:length-1)
      REAL_T      uExt(0:length-1)
      REAL_T      divuExt(0:length-1)
      REAL_T      rhoExt(0:length-1)
      REAL_T hx
      
c     Local variables
      integer NstripMAX
      parameter (NstripMAX = 2000)
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

c     This description assumes outflow at yhi; however, code works for 
c     outflow at any face.      
c     Solve d/dx( 1/rho d/dx( phi ) ) = dU/dx - (S - S_ave) [S = divu if U is 
c     zero, S = d/dt(divu) if U = (ustar - uold)/dt] with periodic or Neumann
c     boundary conditions, using a tridiagonal solve which detects, and deals 
c     with the singular equations.  In the Neumann case, arbitrarily set the 
c     upper right corner to zero to pin the solution.  Note that the RHS of 
c     this equation satisfies the solvability constraint that 
c     Int[RHS.dV] = 0 by construction.
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
         
         hdr = 0.0D0 ! FIXME: This wasn't set before, used below

c     Do left-side periodic BC (since first/last node coincide, use first 
c     node only (retain r-stuff here, just to be sure scaling is not destroyed)
         neq = ine - ins
         if ( neq .gt. NStripMax ) then
            STOP 'HGPHIBC: NstripMax too small'
         end if
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
         hdr = half*(r(ics+1) - r(ics))
         neq = ine - ins + 1
         if ( neq .gt. NStripMax ) then
            STOP 'HGPHIBC: NstripMax too small'
         end if
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
         
c     Solve the equations (we know they're singular, pass the arbitrary value, 
c     and a flag that we've already normalized the rhs, in the sense that
c                          Int[dU/dx - (S-S_ave)] == 0
         sVal = zero
         rNormed = .true.
         call tridag_sing(a,b,c,s,x,neq,sVal,rNormed)
      end if
      
c     Try normalizing phi to average to zero
      phitot = zero
      vtot = zero
C     do i = 1, ine-ins-1
      do i = 2, ine-ins-1
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
      do i = 1,neq
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
      PARAMETER (NMAX=2000)
      PARAMETER (eps=1.D-12)
      INTEGER j
      REAL_T bet,gam(NMAX)
      logical rNormed, singular
      singular = .false.
      if(n.gt.NMAX)then
         STOP 'TRIDAG_SING: NMAX too small in tridag'
      end if
      if((b(1).eq.zero).or.(ABS(b(n)).le.eps))then
         STOP 'TRIDAG_SING: rewrite equations'
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
              STOP 'tridag failed'
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
      PARAMETER (NMAX=2000)
      INTEGER j
      REAL_T bet,gam(NMAX)
      if(n.gt.NMAX)then
         STOP 'TRIDAG: NMAX too small in tridag'
      end if
      if(b(1) .eq. 0.0D0) then
         STOP 'TRIDAG: rewrite equations'
      end if
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet .eq. 0.0D0) STOP 'tridag failed'
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
      PARAMETER (NMAX=2000)
      INTEGER i
      REAL_T fact,gamma,bb(NMAX),u(NMAX),z(NMAX),den
      if(n.le.2)then
         STOP 'CYCLIC: n too small in cyclic'
      end if
      if(n.gt.NMAX)then
         STOP 'CYCLIC: NMAX too small in cyclic'
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
        u(i)=0.0D0
12    continue
      call tridag(a,bb,c,u,z,n)
      den=one+z(1)+beta*z(n)/gamma
      fact=cvmgt(zero,(x(1)+beta*x(n)/gamma)/den,den.eq.zero)
      do 13 i=1,n
        x(i)=x(i)-fact*z(i)
13    continue
      return
      END


