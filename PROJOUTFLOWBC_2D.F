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

c *************************************************************************
c ** EXTRAP_PROJ **
c *************************************************************************

      subroutine FORT_EXTRAP_PROJ(DIMS(u),u,DIMS(divu),divu,DIMS(rho),rho,
     &     r_len,redge,DIMS(uExt),uExt,DIMS(divuExt),divuExt,
     &     DIMS(rhoExt),rhoExt,lo,hi,face,zeroIt,hx)
      implicit none

c subtract divu_ave twice due to precision problems

      integer DIMDEC(u)
      integer DIMDEC(divu)
      integer DIMDEC(rho)
      integer DIMDEC(uExt)
      integer DIMDEC(divuExt)
      integer DIMDEC(rhoExt)
      integer r_len
      integer lo(SDIM),hi(SDIM)
      integer face
      REAL_T      u(DIMV(u),SDIM)
      REAL_T   divu(DIMV(divu))
      REAL_T    rho(DIMV(rho))
      REAL_T      uExt(DIMV(uExt),SDIM-1)
      REAL_T   divuExt(DIMV(divuExt))
      REAL_T   r hoExt(DIMV(rhoExt))
      REAL_T   redge(0:r_len-1)
      REAL_T   hx
      integer  zeroIt

c local variables
      integer ics,ice,jcs,jce
      integer ife,jfe
      integer if,jf
      REAL_T divu_ave1,divu_ave2
      REAL_T max_divu
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

      ife = hi(1)
      jfe = hi(2)

      zeroIt = 0

      if (face .eq. XLO) then
         if=ife
         do j = jcs+1, jce-1
               uExt(j,if,1) = half*(three*   u(ice-1,j,2) -    u(ice,j,2))
            divuExt(j,if  ) = half*(three*divu(ice-1,j  ) - divu(ice,j))
             rhoExt(j,if  ) = half*(three* rho(ice-1,j  ) -  rho(ice,j))
         end do
         max_divu = ABS(divuExt(jcs+1,if))
         do j = jcs+1, jce-1
            max_divu = max(max_divu,ABS(divuExt(j,if)))
         end do
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(jcs+1,if))
         do j = jcs+1, jce-1
            max_pert = MAX(max_pert,ABS(divuExt(j,if)))
         end do
      else if (face .eq. YLO) then
         jf = jfe
         do i = ics+1, ice-1
               uExt(i,jf,1) = half*(three*   u(i,jce-1,1) -    u(i,jce,1))
            divuExt(i,jf  ) = half*(three*divu(i,jce-1  ) - divu(i,jce))
             rhoExt(i,jf  ) = half*(three* rho(i,jce-1  ) -  rho(i,jce))
         end do
         max_divu = ABS(divuExt(ics+1,jf))
         do i = ics+1, ice-1
            max_divu = max(max_divu,ABS(divuExt(i,jf)))
         end do
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(ics+1,jf))
         do i = ics+1, ice-1
            max_pert = MAX(max_pert,ABS(divuExt(i,jf)))
         end do
      else if (face .eq. XHI) then
         if = ife
         do j = jcs+1, jce-1
               uExt(j,if,1) = half*(three*   u(ics+1,j,2) -    u(ics,j,2))
            divuExt(j,if  ) = half*(three*divu(ics+1,j  ) - divu(ics,j))
             rhoExt(j,if  ) = half*(three* rho(ics+1,j  ) -  rho(ics,j))
         end do
         max_divu = ABS(divuExt(jcs+1,if))
         do j = jcs+1, jce-1
            max_divu = max(max_divu,ABS(divuExt(j,if)))
         end do
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(jcs+1,if))
         do j = jcs+1, jce-1
            max_pert = MAX(max_pert,ABS(divuExt(j,if)))
         end do
      else if (face .eq. YHI) then
         jf = jfe
         do i = ics+1, ice-1
               uExt(i,jf,1) = half*(three*   u(i,jcs+1,1) -    u(i,jcs,1))
            divuExt(i,jf  ) = half*(three*divu(i,jcs+1  ) - divu(i,jcs))
             rhoExt(i,jf  ) = half*(three* rho(i,jcs+1  ) -  rho(i,jcs))
         end do
         max_divu = ABS(divuExt(ics+1,jf))
         do i = ics+1, ice-1
            max_divu = max(max_divu,ABS(divuExt(i,jf)))
         end do
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(ics+1,jf))
         do i = ics+1, ice-1
            max_pert = MAX(max_pert,ABS(divuExt(i,jf)))
         end do
      endif

c check to see if we should zero phi
      max_pert = max_pert/(ABS(divu_ave1+divu_ave2)+small_pert)
      if ((max_divu.eq.zero) .or. (max_pert.le.small_pert)) zeroIt = 1
#undef XLO
#undef YLO
#undef XHI
#undef YHI

      end

c *************************************************************************
c ** HGRELAX **
c *************************************************************************

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

c *************************************************************************
c ** HGSUBTRACTAVGPHI **
c *************************************************************************

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


c *************************************************************************
c ** HGRESID **
c *************************************************************************

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


c *************************************************************************
c ** SETPROJBC **
c *************************************************************************

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

c *************************************************************************
c ** HG_SHIFT_PHI **
c *************************************************************************

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

c *************************************************************************
c ** HG_RESHIFT_PHI **
c *************************************************************************

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

c *************************************************************************
c ** HG_SOLVEHG **
c *************************************************************************

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

c *************************************************************************
c ** makeprojdgphi **
c *************************************************************************

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

c *************************************************************************
c ** HGPHIBC **
c ** Solution by back substitution
c *************************************************************************

      subroutine FORT_HGPHIBC(hx,sigExt,s,x,length,per)
c
c    Compute the value of phi for hgproj to be used at an outflow face,
c    assuming that the tangential velocity on the edges of the outflow boundary
c    are either zero or periodic.  
c
      implicit none

      integer length
      integer per
      REAL_T       s(  length+1)
      REAL_T       x(  length+1)
      REAL_T       r( 0:length-1)
      REAL_T  sigExt(0:length-1)
      REAL_T hx
      
c     Local variables
      integer NstripMAX
      parameter (NstripMAX = 2000)
      integer ics, ice
      integer i, j, icL, icR, neq
      REAL_T a(NstripMAX), b(NstripMAX), c(NstripMAX)
      REAL_T alpha, beta, sVal
      logical rNormed
      REAL_T vtot, rnode, hdr
      REAL_T phitot
#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3
      
      ics = 0
      ice = length-1

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

      do i = 2,length
         icR = ics + i - 1
         icL = ics + i - 2
         a(i) = sigExt(icL)
         c(i) = sigExt(icR)
         b(i) = - a(i) - c(i)
      end do

      if (per .eq. 1) then
         
         hdr = 0.0d0

c     Do left-side periodic BC (since first/last node coincide, use first 
c     node only (retain r-stuff here, just to be sure scaling is not destroyed)
         neq = length
         if ( neq .gt. NStripMax ) then
            call bl_abort('HGPHIBC: NstripMax too small')
         end if
         icL = ice
         icR = ics
         beta = sigExt(icL)
         c(1) = sigExt(icR)
         b(1) = - beta - c(1)
         
c     Do right-side periodic on penultimate node
         icL = ice - 1
         icR = ice
         a(neq) = sigExt(icL)
         alpha  = sigExt(icR)
         b(neq) = - a(neq) - alpha
         
         call cyclic(a,b,c,alpha,beta,s,x,neq)
      else

c     Solid walls, Neumann conditions

         hdr = half*(r(ics+1) - r(ics))
         neq = length + 1
         if ( neq .gt. NStripMax ) then
            call bl_abort('HGPHIBC: NstripMax too small')
         end if
         icR = ics
         c(1) = sigExt(icR)
         b(1) = - c(1)

         icL = ice
         a(neq) = sigExt(icL)
         b(neq) = - a(neq)

c     Solve the equations (we know they're singular, pass the arbitrary value, 
c     and a flag that we've already normalized the rhs, in the sense that
c                          Int[dU/dx - (S-S_ave)] == 0
         sVal = zero
         rNormed = .true.
         call tridag_sing(a,b,c,s,x,neq,sVal,rNormed)
      end if
      
#if 1
c     Try normalizing phi to average to zero
      phitot = zero
      vtot = zero
C     do i = 1, length-1
      do i = 2, length-1
c        rnode = r(ics+i-1) - hdr
c        phitot = phitot + x(i)*rnode
c        vtot = vtot + rnode

         phitot = phitot + x(i)
           vtot =   vtot + one

      end do
      do i = length,neq

c        rnode = r(ics+i-2) + hdr
c        phitot = phitot + x(i)*rnode
c        vtot = vtot + rnode

         phitot = phitot + x(i)
         vtot = vtot + one

      end do
      phitot = phitot / vtot
      do i = 1,neq
         x(i) = x(i) - phitot
      end do
#endif
      
#undef XLO
#undef YLO
#undef XHI
#undef YHI
      end

c *************************************************************************
c ** OLDHGPHIBC **
c ** Solution by back substitution
c *************************************************************************

      subroutine FORT_OLDHGPHIBC(hx,r,uExt, divuExt,rhoExt,x,flag,length,per)
c
c    Compute the value of phi for hgproj to be used at an outflow face,
c    assuming that the tangential velocity on the edges of the outflow boundary
c    are either zero or periodic.  
c
      implicit none

      integer length
      integer per
      REAL_T       x(  length+1)
      REAL_T       r( 0:length-1)
      REAL_T    uExt(0:length-1)
      REAL_T divuExt(0:length-1)
      REAL_T  rhoExt(0:length-1)
      REAL_T    flag(0:length-1)
      REAL_T hx
      
c     Local variables
      integer NstripMAX
      parameter (NstripMAX = 2000)
      integer ics, ice
      integer i, j, icL, icR, neq
      REAL_T a(NstripMAX), b(NstripMAX), c(NstripMAX), s(NstripMAX)
      REAL_T alpha, beta, sVal
      logical rNormed
      REAL_T vtot, rnode, hdr
      REAL_T phitot
#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3
      
      ics = 0
      ice = length-1

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

      do i = 2,length
         icR = ics + i - 1
         icL = ics + i - 2
         a(i) = r(icL) / rhoExt(icL)
         c(i) = r(icR) / rhoExt(icR)
         b(i) = - a(i) - c(i)
         rnode = half*(r(icL)+r(icR))
         s(i) = (r(icR)*flag(icR)*uExt(icR)-r(icL)*flag(icL)*uExt(icL))*hx
     &        -  rnode*half*(divuExt(icL)+divuExt(icR))*hx*hx
      end do

      if (per .eq. 1) then
         
         hdr = 0.0D0 ! FIXME: This wasn't set before, used below

c     Do left-side periodic BC (since first/last node coincide, use first 
c     node only (retain r-stuff here, just to be sure scaling is not destroyed)
         neq = length
         if ( neq .gt. NStripMax ) then
            call bl_abort('HGPHIBC: NstripMax too small')
         end if
         icL = ice
         icR = ics
         beta = r(icL) / rhoExt(icL)
         c(1) = r(icR) / rhoExt(icR)
         b(1) = - beta - c(1)
         rnode = half*(r(icL)+r(icR))
         s(1) = (flag(icR)*r(icR)*uExt(icR)-flag(icL)*r(icL)*uExt(icL))*hx
     &        -  rnode*half*(divuExt(icL)+divuExt(icR))*hx*hx
         
c     Do right-side periodic on penultimate node
         icL = ice - 1
         icR = ice
         a(neq) = r(icL) / rhoExt(icL)
         alpha  = r(icR) / rhoExt(icR)
         b(neq) = - a(neq) - alpha
         s(neq) = (flag(icR)*r(icR)*uExt(icR)-flag(icL)*r(icL)*uExt(icL))*hx
     &          -  rnode*half*(divuExt(icL)+divuExt(icR))*hx*hx
         
c     Solve the equations
         call cyclic(a,b,c,alpha,beta,s,x,neq)
      else

c     Solid walls, Neumann conditions
         hdr = half*(r(ics+1) - r(ics))
         neq = length + 1
         if ( neq .gt. NStripMax ) then
            call bl_abort('HGPHIBC: NstripMax too small')
         end if
         icR = ics
         c(1) = r(icR) / rhoExt(icR-ics)
         b(1) = - c(1)
         rnode = r(ics) - hdr
         s(1) =  flag(icR)*r(icR)*uExt(icR)*hx
     &        -  rnode*half*divuExt(icR)*hx*hx

         icL = ice
         a(neq) = r(icL) / rhoExt(icL)
         b(neq) = - a(neq)
         rnode = r(ice) + hdr
         s(neq) = -flag(icL)*r(icL)*uExt(icL)*hx
     &        -  rnode*half*divuExt(icL)*hx*hx
         
c     Solve the equations (we know they're singular, pass the arbitrary value, 
c     and a flag that we've already normalized the rhs, in the sense that
c                          Int[dU/dx - (S-S_ave)] == 0
         sVal = zero
         rNormed = .true.
         call tridag_sing(a,b,c,s,x,neq,sVal,rNormed)
      end if
      
#if 1
c     Try normalizing phi to average to zero
      phitot = zero
      vtot = zero
C     do i = 1, length-1
      do i = 2, length-1
         rnode = r(ics+i-1) - hdr
         phitot = phitot + x(i)*rnode
         vtot = vtot + rnode
      end do
      do i = length,neq
         rnode = r(ics+i-2) + hdr
         phitot = phitot + x(i)*rnode
         vtot = vtot + rnode
      end do
      phitot = phitot / vtot
      do i = 1,neq
         x(i) = x(i) - phitot
      end do
#endif
      
#undef XLO
#undef YLO
#undef XHI
#undef YHI
      end

c *************************************************************************
c ** RHOGBC **
c *************************************************************************

      subroutine FORT_RHOGBC(rho,DIMS(rho),phi,DIMS(phi),
     &                       face,gravity,dx,domlo,domhi,lo_bc,hi_bc)
c
c    Compute the contribution of gravity to the boundary conditions
c      for phi at outflow faces only.
c
      implicit none

      integer DIMDEC(rho)
      integer DIMDEC(phi)
      integer face
      integer domlo(2)
      integer domhi(2)
      integer lo_bc(2)
      integer hi_bc(2)
      REAL_T  rho(DIMV(rho))
      REAL_T  phi(DIMV(phi))
      REAL_T  dx(2)
      REAL_T  gravity
      
c     Local variables
      integer i,j
      REAL_T rhog
      REAL_T rhoExt

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3

      if (face .eq. YLO .or. face .eq. YHI) 
     $   call bl_abort('SHOULDNT BE IN RHOGBC WITH FACE IN Y-DIR')

c     Ok to only use low index of phi because phi is only one
c        node wide in i-direction.
      i = ARG_L1(phi)

      if (face .eq. XLO) then

        rhog = zero
        do j = ARG_H2(phi)-1,ARG_L2(phi),-1
          rhoExt = half*(three*rho(i,j)-rho(i+1,j))
          rhog = rhog - gravity * rhoExt * dx(2)
          phi(i,j) = phi(i,j) + rhog
        end do

      else if (face .eq. XHI) then

        rhog = zero
        do j = ARG_H2(phi)-1,ARG_L2(phi),-1
          rhoExt  = half*(three*rho(i-1,j)-rho(i-2,j))
          rhog = rhog - gravity * rhoExt * dx(2)
          phi(i,j) = phi(i,j) + rhog
        end do

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
         call bl_abort('TRIDAG_SING: NMAX too small in tridag')
      end if
      if((b(1).eq.zero).or.(ABS(b(n)).le.eps))then
         call bl_abort('TRIDAG_SING: rewrite equations')
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
              call bl_abort('tridag failed')
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
         call bl_abort('TRIDAG: NMAX too small in tridag')
      end if
      if(b(1) .eq. 0.0D0) then
         call bl_abort('TRIDAG: rewrite equations')
      end if
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet .eq. 0.0D0) call bl_abort('tridag failed')
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
         call bl_abort('CYCLIC: n too small in cyclic')
      end if
      if(n.gt.NMAX)then
         call bl_abort('CYCLIC: NMAX too small in cyclic')
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


c *************************************************************************
c ** FILL_ONED **
c *************************************************************************

      subroutine FORT_FILL_ONED(lenx,leny,length,faces,numOutFlowFaces,
     $                          cc0,cc1,cc2,cc3,r0,r1,r2,r3,cc_conn,nodal_conn,per,hx,hy)

      integer lenx,leny,length
      integer faces(4)
      integer numOutFlowFaces
      integer per
      REAL_T cc0(0:leny+1,3)
      REAL_T cc1(0:lenx+1,3)
      REAL_T cc2(0:leny+1,3)
      REAL_T cc3(0:lenx+1,3)
      REAL_T  r0(leny)
      REAL_T  r1(lenx)
      REAL_T  r2(leny)
      REAL_T  r3(lenx)
      REAL_T     cc_conn(length)
      REAL_T  nodal_conn(length+1)
      REAL_T  hx,hy

      integer xlo_outflow,ylo_outflow
      integer xhi_outflow,yhi_outflow
      integer i,ifinal,n
      REAL_T  rnode,vol,sum

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
c       XLO 
c           YLO 
c               XHI
c                   YHI
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
c       ordered so that the 2nd spatial dimension is one cell wide.

c     cc*(i,1) = rho
c     cc*(i,2) = divu
c     cc*(i,3) = tangential vel.
 
      vol = hx*hy

      ifinal = 0

      if (numOutFlowFaces .eq. 4 .or. 
     $    (xlo_outflow .eq. 1 .and. ylo_outflow .eq. 0) ) then

          rnode = r0(1)
          do i = 1,leny
            cc_conn(i) = rnode/cc0(i,1)
          enddo

          do i = 2,leny
            nodal_conn(i) = rnode*(      (cc0(i,3) - cc0(i-1,3))*hx
     &                             -half*(cc0(i,2)+cc0(i-1,2))*vol )
          enddo
          i = 1
          nodal_conn(i) = rnode*(     (cc0(i,3) - cc0(i-1,3))*hx
     &                           -half*cc0(i,2)*vol )
          if (per .eq. 1) then
            nodal_conn(i) = nodal_conn(i) - half*rnode*cc0(leny,2)*vol 
          endif

          i = leny+1
          nodal_conn(i) = rnode*(     (cc0(i,3) - cc0(i-1,3))*hx
     &                           -half*cc0(i-1,2)*vol )
          ifinal = leny

      endif

      if (yhi_outflow .eq. 1 .and. 
     $    .not. (numOutFlowFaces .eq. 3 .and. xhi_outflow .eq. 0) ) then

          do i = 1,lenx
            cc_conn(ifinal+i) = r3(i)/cc3(i,1)
          enddo

          do i = 2,lenx
            rnode = half*(r3(i)+r3(i-1))
            nodal_conn(ifinal+i) = rnode*(      (cc3(i,3) - cc3(i-1,3))*hy
     &                                    -half*(cc3(i,2)+cc3(i-1,2))*vol )
          enddo
          i = 1
          rnode = 1.5d0*r3(1) - 0.5d0*r3(2)
c         Note: we get away with using r3(i) in place of r3(i-1), which isnt
c               defined, because if this is r-z, then the vel. at i-1 is zero,
c               and if not, then r is identically 1.
          nodal_conn(ifinal+i) = nodal_conn(ifinal+i) + 
     &                           r3(i)*(cc3(i,3) - cc3(i-1,3))*hy
     &                           -half*rnode*cc3(i,2)*vol

          if (per .eq. 1) then
            nodal_conn(ifinal+i) = nodal_conn(ifinal+i) - half*rnode*cc3(lenx,2)*vol 
          endif

          i = lenx+1
          rnode = 1.5d0*r3(lenx) - 0.5d0*r3(lenx-1)
c         Note: same reasoning as above for using r3(i-1)
          nodal_conn(ifinal+i) = 
     &                           r3(i-1)*(cc3(i,3) - cc3(i-1,3))*hy
     &                           -half*rnode*cc3(i-1,2)*vol

          ifinal = ifinal + lenx

      endif

      if (xhi_outflow .eq. 1) then
          rnode = r2(1)
          do i = 1,leny
            cc_conn(ifinal+i) = rnode/cc2(leny+1-i,1)
          enddo

          do i = 2,leny
            nodal_conn(ifinal+i) = rnode*(      (cc2(leny+2-i,3)-cc2(leny+1-i,3))*hx
     &                                    -half*(cc2(leny+2-i,2)+cc2(leny+1-i,2))*vol )
          enddo
          i = 1
          nodal_conn(ifinal+i) = nodal_conn(ifinal+i) +
     &                        rnode*(     (cc2(leny+2-i,3) - cc2(leny+1-i,3))*hx
     &                               -half*cc2(leny+1-i,2)*vol )
          if (per .eq. 1) then
            nodal_conn(ifinal+i) = nodal_conn(ifinal+i) - half*rnode*cc2(1,2)*vol 
          endif

          i = leny+1
          nodal_conn(ifinal+i) = 
     &                        rnode*(     (cc2(leny+2-i,3) - cc2(leny+1-i,3))*hx
     &                               -half*cc2(leny+2-i,2)*vol )

          do i = 1,leny+1
            nodal_conn(ifinal+i) = -nodal_conn(ifinal+i)
          enddo

          ifinal = ifinal + leny
      endif


      if (ylo_outflow .eq. 1) then

          do i = 1,lenx
            cc_conn(ifinal+i) = r1(lenx+1-i)/cc1(lenx+1-i,1)
          enddo

          do i = 2,lenx
            rnode = half*(r1(lenx+2-i)+r1(lenx+1-i))
            nodal_conn(ifinal+i) = rnode*(      (cc1(lenx+2-i,3)-cc1(lenx+1-i,3))*hy
     &                                    -half*(cc1(lenx+2-i,2)+cc1(lenx+1-i,2))*vol )
          enddo
          i = 1
          rnode = 1.5d0*r1(lenx) - 0.5d0*r1(lenx-1)
c         Note: we get away with using r1(lenx)
c               because if this is r-z, then the outside vel.
c               and if not, then r is identically 1.
          nodal_conn(ifinal+i) = nodal_conn(ifinal+i) + 
     &                           r1(lenx)*(cc1(lenx+2-i,3) - cc1(lenx+1-i,3))*hy
     &                           -half*rnode*cc1(lenx+1-i,2)*vol
          if (per .eq. 1) then
            nodal_conn(ifinal+i) = nodal_conn(ifinal+i) - half*rnode*cc1(1,2)*vol 
          endif

          i = lenx+1
          rnode = 1.5d0*r1(2) - 0.5d0*r1(1)
c         Note: same reasoning as above for using r1(1)
          nodal_conn(ifinal+i) = 
     &                           r1(1)*(cc1(lenx+2-i,3) - cc1(lenx+1-i,3))*hy
     &                           -half*rnode*cc1(lenx+2-i,2)*vol

          do i = 1,lenx+1
            nodal_conn(ifinal+i) = -nodal_conn(ifinal+i)
          enddo
          ifinal = ifinal + lenx
      endif

      if (numOutFlowFaces .lt. 4 .and.
     $    (xlo_outflow .eq. 1 .and. ylo_outflow .eq. 1) ) then

          rnode = r0(1)
          do i = 1,leny
            cc_conn(ifinal+i) = rnode/cc0(i,1)
          enddo
          do i = 2,leny
            nodal_conn(ifinal+i) = rnode*(      (cc0(i,3) - cc0(i-1,3))*hx
     &                                    -half*(cc0(i,2)+cc0(i-1,2))*vol )
          enddo
          i = 1
          nodal_conn(ifinal+i) = nodal_conn(ifinal+i) + 
     &                           rnode*(     (cc0(i,3) - cc0(i-1,3))*hx
     &                                  -half*cc0(i,2)*vol )
          i = leny+1
          nodal_conn(ifinal+i) = rnode*(     (cc0(i,3) - cc0(i-1,3))*hx
     &                           -half*cc0(i-1,2)*vol )

          ifinal = ifinal + leny
      endif

      if (yhi_outflow .eq. 1 .and. 
     $    (numOutFlowFaces .eq. 3 .and. xhi_outflow .eq. 0) ) then
          do i = 1,lenx
            cc_conn(ifinal+i) = r3(i)/cc3(i,1)
          enddo

          do i = 2,lenx
            rnode = half*(r3(i)+r3(i-1))
            nodal_conn(ifinal+i) = rnode*(      (cc3(i,3) - cc3(i-1,3))*hy
     &                                    -half*(cc3(i,2)+cc3(i-1,2))*vol )
          enddo
          i = 1
          rnode = 1.5d0*r3(1) - 0.5d0*r3(2)
c         Note: we get away with using r3(i) in place of r3(i-1), which isnt
c               defined, because if this is r-z, then the vel. at i-1 is zero,
c               and if not, then r is identically 1.
          nodal_conn(ifinal+i) = nodal_conn(ifinal+i) + 
     &                           r3(i)*(cc3(i,3) - cc3(i-1,3))*hy
     &                           -half*rnode*cc3(i,2)*vol
          i = lenx+1
          rnode = 1.5d0*r3(lenx) - 0.5d0*r3(lenx-1)
c         Note: same reasoning as above for using r3(i-1)
          nodal_conn(ifinal+i) = 
     &                           r3(i-1)*(cc3(i,3) - cc3(i-1,3))*hy
     &                           -half*rnode*cc3(i-1,2)*vol

          ifinal = ifinal + lenx
      endif

      length = ifinal

#undef XLO
#undef YLO
#undef XHI
#undef YHI

      end

c *************************************************************************
c ** ALLPHI_FROM_X **
c *************************************************************************

      subroutine FORT_ALLPHI_FROM_X(lenx,leny,length,faces,numOutFlowFaces,
     $                              phi0,phi1,phi2,phi3,x)

      integer lenx,leny,length
      integer numOutFlowFaces
      integer faces(4)
      REAL_T phi0(0:leny)
      REAL_T phi1(0:lenx)
      REAL_T phi2(0:leny)
      REAL_T phi3(0:lenx)
      REAL_T x(0:length)

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
        do j = 0,leny
          phi0(j) = x(j)
        enddo
        ifinal = leny
      endif

      if (yhi_outflow .eq. 1 .and. 
     $    .not. (numOutFlowFaces .eq. 3 .and. xhi_outflow .eq. 0) ) then
        do i = 0,lenx
          phi3(i) = x(i+ifinal)
        enddo
        ifinal = ifinal+lenx
      endif

      if (xhi_outflow .eq. 1) then
        do j = 0,leny
          phi2(leny-j) = x(ifinal+j)
        enddo
        ifinal = ifinal+leny
      endif

      if (ylo_outflow .eq. 1) then
        if (numOutFlowFaces .eq. 4) then 
          do i = 0,lenx-1
            phi1(lenx-i) = x(ifinal+i)
          enddo
          phi1(0) = x(0)
        else
          do i = 0,lenx
            phi1(lenx-i) = x(ifinal+i)
          enddo
        endif
        ifinal = ifinal+lenx
      endif

      if (numOutFlowFaces .lt. 4 .and.
     $    (xlo_outflow .eq. 1 .and. ylo_outflow .eq. 1) ) then
        do j = 0,leny
          phi0(j) = x(j+ifinal)
        enddo
        ifinal = ifinal+leny
      endif

      if (yhi_outflow .eq. 1 .and. 
     $    (numOutFlowFaces .eq. 3 .and. xhi_outflow .eq. 0) ) then
        do i = 0,lenx
          phi3(i) = x(i+ifinal)
        enddo
        ifinal = ifinal+lenx
      endif

      end


c *************************************************************************
c ** PHI_FROM_X **
c *************************************************************************

      subroutine FORT_PHI_FROM_X(DIMS(phi),phi,length,x,per)

      integer DIMDEC(phi)
      integer length
      REAL_T phi(DIMV(phi))
      REAL_T x(length+1)
      integer per

      integer lenx, leny
      integer i,j

c     We know that the faces are ordered: XLO,XHI,YLO,YHI
      lenx = ARG_H1(phi)-ARG_L1(phi)
      leny = ARG_H2(phi)-ARG_L2(phi)

      if (lenx .eq. 0) then
        do j = 1,length
          phi(ARG_L1(phi),j-1) = x(j)
        enddo
        if (per .eq. 1) then
          phi(ARG_L1(phi),ARG_H2(phi)) = phi(ARG_L1(phi),ARG_L2(phi))
        else
          phi(ARG_L1(phi),ARG_H2(phi)) = x(length+1)
        endif
      elseif (leny .eq. 0) then
        do i = 1,length
          phi(i-1,ARG_L2(phi)) = x(i)
        enddo
        if (per .eq. 1) then
          phi(ARG_H1(phi),ARG_L2(phi)) = phi(ARG_L1(phi),ARG_L2(phi))
        else
          phi(ARG_H1(phi),ARG_L2(phi)) = x(length+1)
        endif
      endif

#undef XLO
#undef YLO
#undef XHI
#undef YHI

      end
