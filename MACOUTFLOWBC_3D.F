#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "MACOUTFLOWBC_F.H"
#include "ArrayLim.H"

#define SDIM 3

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define SMALL 1.0e-10
#else
#define SMALL 1.0d-10
#endif

      subroutine FORT_EXTRAP_MAC(DIMS(u0),u0,DIMS(u1),u1,DIMS(u2),u2,DIMS(div),
     &                           divu,DIMS(rho),rho,
     &                           DIMS(divuExt),divuExt,DIMS(rhoExt),rhoExt,
     &                           dx,lo,hi,face,per,zeroIt)
c
c     Compute the value of phi for macproj 
c     assuming that the tangential velocity on the edges of the outflow boundary
c     are either zero or periodic.
c     note that u is edge centered

c    compute divu_ave twice due to precision problems

      implicit none

      integer DIMDEC(u0)
      integer DIMDEC(u1)
      integer DIMDEC(u2)
      integer DIMDEC(div)
      integer DIMDEC(rho)
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(divuExt)
      integer DIMDEC(rhoExt)
      REAL_T      u0(DIMV(u0))
      REAL_T      u1(DIMV(u1))
      REAL_T      u2(DIMV(u2))
      REAL_T divu(DIMV(div))
      REAL_T    rho(DIMV(rho))
      REAL_T   divuExt(DIMV(divuExt))
      REAL_T   rhoExt(DIMV(rhoExt))
      REAL_T   dx(3)
      integer per(2)
      integer zeroIt
      integer face

c     Local variables
      REAL_T small_pert
      REAL_T small_udiff
      parameter ( small_pert = SMALL)
      parameter (small_udiff = SMALL)
      REAL_T max_divu, min_divu, max_pert
      REAL_T divu_ave1,divu_ave2
      REAL_T hx,hy,hz  
      REAL_T diff
      integer i,j,k,neq,n
      integer ics,ice,jcs,jce,kcs,kce
      integer ifs,ife,jfs,jfe,kfs,kfe
      integer if,jf,kf,ic,jc,kc
c     NOTE: Assumes that rho at edge between i, i-1 = half*(rho(i)+rho(i-1))
c             (1) Linear fit of rho between nodes
c             (2) rho, divu on same boxes (box)
c             (3) phi is on box, shifted up one
c             (4) u is edge-based, on surroundingNodes(box)

c     Compute average of divu over outflow bc.  Set trivial solution if average
c     is zero, or if divu is constant
#define XLO 0
#define YLO 1
#define ZLO 2
#define XHI 3
#define YHI 4
#define ZHI 5
      ics = ARG_L1(rho)
      ice = ARG_H1(rho)
      jcs = ARG_L2(rho)
      jce = ARG_H2(rho)
      kcs = ARG_L3(rho)
      kce = ARG_H3(rho)

      ifs = lo(1)
      ife = hi(1)
      jfs = lo(2)
      jfe = hi(2)
      kfs = lo(3)
      kfe = hi(3)

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      zeroIt = 0

      if (face .eq. XLO) then
         if = ifs
         max_divu = divu(ics,jcs,kcs)
         min_divu = max_divu
         do k = kcs, kce
         do j = jcs, jce
            divuExt(j,k,if) = divu(ics,j,k)
            rhoExt(j,k,if)  = rho(ics,j,k)
            max_divu = max(max_divu,divuExt(j,k,if))
            min_divu = min(min_divu,divuExt(j,k,if))
         end do
         end do

c        Here we modify divuExt to include the velocity terms.
         do k = kcs, kce
         do j = jcs, jce
            divuExt(j,k,if) = divuExt(j,k,if) 
     $        - (u1(ics,j+1,k)-u1(ics,j,k))/hy 
     $        - (u2(ics,j,k+1)-u2(ics,j,k))/hz 
         end do
         end do

         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave2,face)

         max_pert = ABS(divuExt(jcs,kcs,if))
         do k = kcs, kce
         do j = jcs, jce
            max_pert = MAX(max_pert,ABS(divuExt(j,k,if)))
         end do
         end do

c        Make sure u_mac is periodic
         if (per(1) .eq. 1) then
           diff = abs(u1(ics,jcs,kcs)-u1(ics,jce+1,kcs))
           do k = kcs, kce
             diff = max(diff,abs(u1(ics,jcs,k)-u1(ics,jce+1,k)))
           enddo
           if (diff .gt. small_udiff) then
             write(6,*) 'EXTRAPMAC: FACE XLO : vmac not periodic'
             call bl_abort(" ")
           endif
         endif
         if (per(2) .eq. 1) then
           diff = abs(u2(ics,jcs,kcs)-u2(ics,jcs,kce+1))
           do j = jcs, jce
             diff = max(diff,abs(u2(ics,j,kcs)-u2(ics,j,kce+1)))
           enddo
           if (diff .gt. small_udiff) then
              write(6,*) 'EXTRAPMAC: FACE XLO : wmac not periodic'
              call bl_abort(" ")
           endif
         endif

      else if (face .eq. YLO) then
         jf = jfs
         max_divu = divu(ics,jcs,kcs)
         min_divu = max_divu
         do k = kcs, kce
         do i = ics, ice
            divuExt(i,k,jf) = divu(i,jcs,k)
            rhoExt(i,k,jf)  = rho(i,jcs,k)
            max_divu = max(max_divu,divuExt(i,k,jf))
            min_divu = min(min_divu,divuExt(i,k,jf))
         end do
         end do

c        Here we modify divuExt to include the velocity terms.
         do k = kcs, kce
         do i = ics, ice
            divuExt(i,k,jf) = divuExt(i,k,jf)
     $        - (u0(i+1,jcs,k)-u0(i,jcs,k))/hx
     $        - (u2(i,jcs,k+1)-u2(i,jcs,k))/hz
         end do
         end do

         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(ics,kcs,jf))
         do k = kcs, kce
         do i = ics, ice
            max_pert = MAX(max_pert,ABS(divuExt(i,k,jf)))
         end do
         end do

c        Make sure u_mac is periodic
         if (per(1) .eq. 1) then
           diff = abs(u0(ics,jcs,kcs)-u0(ice+1,jcs,kcs))
           do k = kcs, kce
             diff = max(diff,abs(u0(ics,jcs,k)-u0(ice+1,jcs,k)))
           enddo
           if (diff .gt. small_udiff) then
             write(6,*) 'EXTRAPMAC: FACE YLO : umac not periodic'
             call bl_abort(" ")
           endif
         endif
         if (per(2) .eq. 1) then
           diff = abs(u2(ics,jcs,kcs)-u2(ics,jcs,kce+1))
           do i = ics, ice
             diff = max(diff,abs(u2(i,jcs,kcs)-u2(i,jcs,kce+1)))
           enddo
           if (diff .gt. small_udiff) then
              write(6,*) 'EXTRAPMAC: FACE YLO : wmac not periodic'
              call bl_abort(" ")
           endif
         endif

      else if (face .eq. ZLO) then
         kf = kfs
         max_divu = divu(ics,jcs,kcs)
         min_divu = max_divu
         do j = jcs, jce
         do i = ics, ice
            divuExt(i,j,kf) = divu(i,j,kcs)
            rhoExt(i,j,kf)  = rho(i,j,kcs)
            max_divu = max(max_divu,divuExt(i,j,kf))
            min_divu = min(min_divu,divuExt(i,j,kf))
         end do
         end do

c        Here we modify divuExt to include the velocity terms.
         do j = jcs, jce
         do i = ics, ice
            divuExt(i,j,kf) = divuExt(i,j,kf)
     $        - (u0(i+1,j,kcs)-u0(i,j,kcs))/hx
     $        - (u1(i,j+1,kcs)-u1(i,j,kcs))/hy
         end do
         end do

         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(ics,jcs,kf))
         do j = jcs, jce
         do i = ics, ice
            max_pert = MAX(max_pert,ABS(divuExt(i,j,kf)))
         end do
         end do

c        Make sure u_mac is periodic
         if (per(1) .eq. 1) then
           diff = abs(u0(ics,jcs,kcs)-u0(ice+1,jcs,kcs))
           do j = jcs, jce
             diff = max(diff,abs(u0(ics,j,kcs)-u0(ice+1,j,kcs)))
           enddo
           if (diff .gt. small_udiff) then
             write(6,*) 'EXTRAPMAC: FACE ZLO : umac not periodic'
             call bl_abort(" ")
           endif
         endif
         if (per(2) .eq. 1) then
           diff = abs(u1(ics,jcs,kcs)-u1(ics,jce+1,kcs))
           do i = ics, ice
             diff = max(diff,abs(u1(i,jcs,kcs)-u1(i,jce+1,kcs)))
           enddo
           if (diff .gt. small_udiff) then
              write(6,*) 'EXTRAPMAC: FACE ZLO : vmac not periodic'
              call bl_abort(" ")
           endif
         endif

      else if (face .eq. XHI) then
         if = ife
         max_divu = divu(ice,jcs,kcs)
         min_divu = max_divu

         do k = kcs, kce
         do j = jcs, jce
            divuExt(j,k,if) = divu(ice,j,k)
            rhoExt(j,k,if)  = rho(ice,j,k)
            max_divu = max(max_divu,divuExt(j,k,if))
            min_divu = min(min_divu,divuExt(j,k,if))
         end do
         end do

         do k = kcs, kce
         do j = jcs, jce
            divuExt(j,k,if) = divuExt(j,k,if) 
     $        - (u1(ice,j+1,k)-u1(ice,j,k))/hy 
     $        - (u2(ice,j,k+1)-u2(ice,j,k))/hz 
         end do
         end do

         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave2,face)

         max_pert = ABS(divuExt(jcs,kcs,if))
         do k = kcs, kce
         do j = jcs, jce
            max_pert = MAX(max_pert,ABS(divuExt(j,k,if)))
         end do
         end do

c        Make sure u_mac is periodic
         if (per(1) .eq. 1) then
           diff = abs(u1(ice,jcs,kcs)-u1(ice,jce+1,kcs))
           do k = kcs, kce
             diff = max(diff,abs(u1(ice,jcs,k)-u1(ice,jce+1,k)))
           enddo
           if (diff .gt. small_udiff) then
             write(6,*) 'EXTRAPMAC: FACE XLO : umac not periodic'
             call bl_abort(" ")
           endif
         endif
         if (per(2) .eq. 1) then
           diff = abs(u2(ice,jcs,kcs)-u2(ice,jcs,kce+1))
           do j = jcs, jce
             diff = max(diff,abs(u2(ice,j,kcs)-u2(ice,j,kce+1)))
           enddo
           if (diff .gt. small_udiff) then
              write(6,*) 'EXTRAPMAC: FACE XLO : umac not periodic'
              call bl_abort(" ")
           endif
         endif

      else if (face .eq. YHI) then
         jf = jfe
         max_divu = divu(ics,jce,kcs)
         min_divu = max_divu

         do k = kcs, kce
         do i = ics, ice
            divuExt(i,k,jf) = divu(i,jce,k)
            rhoExt(i,k,jf)  = rho(i,jce,k)
            max_divu = max(max_divu,divuExt(i,k,jf))
            min_divu = min(min_divu,divuExt(i,k,jf))
         end do
         end do

         do k = kcs, kce
         do i = ics, ice
            divuExt(i,k,jf) = divuExt(i,k,jf)
     $        - (u0(i+1,jce,k)-u0(i,jce,k))/hx
     $        - (u2(i,jce,k+1)-u2(i,jce,k))/hz
         end do
         end do

         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(ics,kcs,jf))
         do k = kcs, kce
         do i = ics, ice
            max_pert = MAX(max_pert,ABS(divuExt(i,k,jf)))
         end do
         end do

c        Make sure u_mac is periodic
         if (per(1) .eq. 1) then
           diff = abs(u0(ics,jce,kcs)-u0(ice+1,jce,kcs))
           do k = kcs, kce
             diff = max(diff,abs(u0(ics,jce,k)-u0(ice+1,jce,k)))
           enddo
           if (diff .gt. small_udiff) then
             write(6,*) 'EXTRAPMAC: FACE YLO : umac not periodic'
             call bl_abort(" ")
           endif
         endif
         if (per(2) .eq. 1) then
           diff = abs(u2(ics,jce,kcs)-u2(ics,jce,kce+1))
           do i = ics, ice
             diff = max(diff,abs(u2(i,jce,kcs)-u2(i,jce,kce+1)))
           enddo
           if (diff .gt. small_udiff) then
              write(6,*) 'EXTRAPMAC: FACE YLO : wmac not periodic'
              call bl_abort(" ")
           endif
         endif

      else if (face .eq. ZHI) then
         kf = kfe
         max_divu = divu(ics,jcs,kce)
         min_divu = max_divu

         do j = jcs, jce
         do i = ics, ice
            divuExt(i,j,kf) = divu(i,j,kce)
            rhoExt(i,j,kf)  = rho(i,j,kce)
            max_divu = max(max_divu,divuExt(i,j,kf))
            min_divu = min(min_divu,divuExt(i,j,kf))
         end do
         end do

c        Here we modify divuExt to include the velocity terms.
         do j = jcs, jce
         do i = ics, ice
            divuExt(i,j,kf) = divuExt(i,j,kf)
     $        - (u0(i+1,j,kce)-u0(i,j,kce))/hx
     $        - (u1(i,j+1,kce)-u1(i,j,kce))/hy
         end do
         end do

         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(ics,jcs,kf))
         do j = jcs, jce
         do i = ics, ice
            max_pert = MAX(max_pert,ABS(divuExt(i,j,kf)))
         end do
         end do

c        Make sure u_mac is periodic
         if (per(1) .eq. 1) then
           diff = abs(u0(ics,jcs,kce)-u0(ice+1,jcs,kce))
           do j = jcs, jce
             diff = max(diff,abs(u0(ics,j,kce)-u0(ice+1,j,kce)))
           enddo
           if (diff .gt. small_udiff) then
             write(6,*) 'EXTRAPMAC: FACE ZHI : umac not periodic'
             call bl_abort(" ")
           endif
         endif
         if (per(2) .eq. 1) then
           diff = abs(u1(ics,jcs,kce)-u1(ics,jce+1,kce))
           do i = ics, ice
             diff = max(diff,abs(u1(i,jcs,kce)-u1(i,jce+1,kce)))
           enddo
           if (diff .gt. small_udiff) then
              write(6,*) 'EXTRAPMAC: FACE ZHI : vmac not periodic'
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


      subroutine subtractavg(DIMS(divu),divu,lo,hi,divu_ave,face)
      implicit none
      integer DIMDEC(divu)
      integer lo(SDIM),hi(SDIM)
      REAL_T divu_ave
      REAL_T divu(DIMV(divu))
      integer face

      integer i,j,k
      REAL_T vtot

      divu_ave = zero
      vtot = zero

      if (face .eq. XLO .or. face .eq. XHI) then
         i = lo(1)
         do k = lo(3),hi(3)
            do j=lo(2),hi(2)
               vtot = vtot+one
               divu_ave = divu_ave+divu(j,k,i)
            enddo
         enddo
         divu_ave = divu_ave/vtot
         do k = lo(3),hi(3)
            do j=lo(2),hi(2)
               divu(j,k,i) = divu(j,k,i) - divu_ave
            enddo
         enddo
      elseif (face .eq. YLO .or. face .eq. YHI) then
         j = lo(2)
         do k = lo(3),hi(3)
            do i=lo(1),hi(1)
               vtot = vtot+one
               divu_ave = divu_ave+divu(i,k,j)
            enddo
         enddo
         divu_ave = divu_ave/vtot
         do k = lo(3),hi(3)
            do i=lo(1),hi(1)
               divu(i,k,j) = divu(i,k,j) - divu_ave
            enddo
         enddo
      elseif(face .eq. ZLO .or. face .eq. ZHI) then
         k = lo(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               vtot = vtot+one
               divu_ave = divu_ave+divu(i,j,k)
            enddo
         enddo
         divu_ave = divu_ave/vtot
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               divu(i,j,k) = divu(i,j,k) - divu_ave
            enddo
         enddo
      else 
         print*, "bad length"
      endif

      end

#undef XLO
#undef YLO
#undef ZLO
#undef XHI
#undef YHI
#undef ZHI

      subroutine FORT_COMPUTE_MACCOEFF(DIMS(rhs),rhs,
     &     DIMS(beta0),beta0,DIMS(beta1),beta1,
     &     DIMS(divuExt),divuExt,
     &     DIMS(rhoExt),rhoExt,lo,hi,h,isPeriodic)
      implicit none
      integer DIMDEC(rhs)
      integer DIMDEC(beta0)
      integer DIMDEC(beta1)
      integer DIMDEC(divuExt)
      integer DIMDEC(rhoExt)
      integer lo(SDIM),hi(SDIM)
      REAL_T divuExt(DIM12(divuExt))
      REAL_T rhoExt(DIM12(rhoExt))
      REAL_T beta0(DIM12(beta0))
      REAL_T beta1(DIM12(beta1))
      REAL_T rhs(DIM12(rhs))
      REAL_T h(SDIM)
      integer isPeriodic(SDIM)
      
      integer i,j
      integer ics,jcs,ice,jce

      ics = lo(1)
      ice = hi(1)
      jcs = lo(2)
      jce = hi(2)

      do j=jcs,jce
        do i=ics,ice
           rhs(i,j) = -divuExt(i,j)
        enddo
      enddo

      do j=jcs,jce
         do i=ics+1,ice
            beta0(i,j) = two/(rhoExt(i,j)+rhoExt(i-1,j))
         enddo
      enddo

      do j=jcs+1,jce
         do i=ics,ice
            beta1(i,j) = two/(rhoExt(i,j)+rhoExt(i,j-1))
         enddo
      enddo
      
      if (isPeriodic(1) .eq. 1) then
         do j=jcs,jce
            beta0(ics  ,j) = two/(rhoExt(ics,j)+rhoExt(ice,j))
            beta0(ice+1,j) = two/(rhoExt(ics,j)+rhoExt(ice,j))
         enddo
      else
         do j=jcs,jce
            beta0(ics  ,j) = zero
            beta0(ice+1,j) = zero
         enddo
      endif
      if (isPeriodic(2) .eq. 1) then
         do i=ics,ice
            beta1(i,jcs  ) = two/(rhoExt(i,jcs)+rhoExt(i,jce))
            beta1(i,jce+1) = two/(rhoExt(i,jcs)+rhoExt(i,jce))
         enddo
      else
         do i=ics,ice
            beta1(i,jcs  ) = zero
            beta1(i,jce+1) = zero
         enddo
      endif

      end


#define DGXY (hxsqinv*(beta0(i+1,j)*(phi(i+1,j)-phi(i  ,j))      \
                      -beta0(i  ,j)*(phi(i  ,j)-phi(i-1,j))) +   \
              hysqinv*(beta1(i,j+1)*(phi(i,j+1)-phi(i,j))        \
                      -beta1(i,j  )*(phi(i,j  )-phi(i,j-1))))

      subroutine FORT_MACRELAX(DIMS(rhs),rhs,DIMS(beta0),beta0,
     &                         DIMS(beta1),beta1,DIMS(phi),phi,
     &                         lo,hi,h,isPeriodic,niter)
      implicit none
      integer DIMDEC(beta0)
      integer DIMDEC(beta1)
      integer DIMDEC(rhs)
      integer DIMDEC(phi)
      REAL_T beta0(DIM12(beta0))
      REAL_T beta1(DIM12(beta1))
      REAL_T rhs(DIM12(rhs))
      REAL_T phi(DIM12(phi))
      integer isPeriodic(SDIM)
      REAL_T h(SDIM)
      integer lo(SDIM),hi(SDIM)

c Local variables
      integer ics,ice,jcs,jce
      logical setSingularPoint
      integer i,j,iter
      integer redblack
      REAL_T lam,dg
      REAL_T phitot, vtot
      REAL_T maxnorm
      REAL_T hxsqinv,hysqinv
      integer niter
      integer iinc

      hxsqinv = one/(h(1)*h(1))
      hysqinv = one/(h(2)*h(2))
      ics = lo(1)
      ice = hi(1)
      jcs = lo(2)
      jce = hi(2)

      if (h(2). gt. 1.5D0*h(1)) then
        call bl_abort("line solve for mac not yet implemented")
      else if (h(1) .gt. 1.5D0*h(2)) then
        call bl_abort("line solve for mac not yet implemented")
      endif

      setSingularPoint = .false.
      call setmacbc(DIMS(phi),phi,lo,hi,isPeriodic,setSingularPoint)

      do iter = 1,niter
         do redblack = 0,1
            do j=jcs,jce
               iinc = mod(j+redblack,2)
               do i=ics+iinc,ice,2
                  dg = DGXY 
                  lam = -one/(hxsqinv*(beta0(i,j)+beta0(i+1,j))
     &                 +hysqinv*(beta1(i,j)+beta1(i,j+1)))
                  phi(i,j) = phi(i,j) + lam*(rhs(i,j)-dg)
               enddo
            end do
            call setmacbc(DIMS(phi),phi,lo,hi,isPeriodic,setSingularPoint)
         enddo
      end do
      
      end

      subroutine FORT_MACSUBTRACTAVGPHI(DIMS(phi),phi,lo,hi,isPeriodic)
      implicit none
      integer DIMDEC(phi)
      REAL_T phi(DIM12(phi))
      integer lo(SDIM),hi(SDIM)
      integer isPeriodic(SDIM)

      REAL_T phitot,vtot
      integer ics,ice,jcs,jce
      integer i,j
      logical setSingularPoint

      ics = lo(1)
      ice = hi(1)
      jcs = lo(2)
      jce = hi(2)
      setSingularPoint = .false.

      phitot = zero
      vtot   = zero

      do j=jcs,jce
      do i=ics,ice
         phitot = phitot+phi(i,j)
         vtot = vtot + one
      enddo
      enddo
      phitot = phitot/vtot

      do j=jcs,jce
      do i=ics,ice
         phi(i,j) = phi(i,j) - phitot
      enddo
      enddo

      call setmacbc(DIMS(phi),phi,lo,hi,isPeriodic,setSingularPoint)

      end

      subroutine FORT_MACRESID(DIMS(rhs),rhs,DIMS(beta0),beta0,
     &                    DIMS(beta1),beta1,DIMS(phi),phi,
     &                    DIMS(resid),resid,lo,hi,h,isPeriodic,maxnorm)
      implicit none
      integer DIMDEC(beta0)
      integer DIMDEC(beta1)
      integer DIMDEC(rhs)
      integer DIMDEC(phi)
      integer DIMDEC(resid)
      REAL_T beta0(DIM12(beta0))
      REAL_T beta1(DIM12(beta1))
      REAL_T rhs(DIM12(rhs))
      REAL_T phi(DIM12(phi))
      REAL_T resid(DIM12(resid))
      integer isPeriodic(SDIM)
      REAL_T h(SDIM)
      integer lo(SDIM),hi(SDIM)

      REAL_T maxnorm
      integer i,j
      REAL_T hxsqinv,hysqinv

      hxsqinv = one/(h(1)*h(1))
      hysqinv = one/(h(2)*h(2))

      maxnorm = zero
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            resid(i,j) = rhs(i,j) - (DGXY)
            maxnorm = max(maxnorm,ABS(resid(i,j)))
         enddo
      enddo

      end


      subroutine setmacbc(DIMS(phi),phi,lo,hi,isPeriodic,setSingularPoint)
      implicit none
      integer DIMDEC(phi)
      REAL_T phi(DIM12(phi))
      integer lo(SDIM),hi(SDIM)
      integer isPeriodic(SDIM)
      logical setSingularPoint
      
      integer ics,ice,jcs,jce
      integer i,j

      ics = lo(1)
      ice = hi(1)
      jcs = lo(2)
      jce = hi(2)

      if (isPeriodic(1) .NE. 1 .AND. isPeriodic(2) .NE. 1 .AND. setSingularPoint) then
         phi(ice,jce) = zero
      endif

      if (isPeriodic(1) .eq. 1) then
         do j=jcs,jce
            phi(ics-1,j) = phi(ice,j)
            phi(ice+1,j) = phi(ics,j)
         enddo
      else
         do j=jcs,jce
            phi(ics-1,j) = phi(ics,j)
            phi(ice+1,j) = phi(ice,j)
         enddo
      endif
      
      if (isPeriodic(2) .eq. 1) then
         do i=ics,ice
            phi(i,jcs-1) = phi(i,jce)
            phi(i,jce+1) = phi(i,jcs)
         enddo
      else
         do i=ics,ice
            phi(i,jcs-1) = phi(i,jcs)
            phi(i,jce+1) = phi(i,jce)
         enddo
      endif

      end

      subroutine FORT_MAC_SHIFT_PHI(DIMS(out),out,DIMS(in),in,face)
      implicit none
      integer face
      integer DIMDEC(in)
      integer DIMDEC(out)
      REAL_T in(DIMV(in))
      REAL_T out(DIMV(out))
      integer i,j,k
#define XLO 0
#define YLO 1
#define ZLO 2
#define XHI 3
#define YHI 4
#define ZHI 5
      if (face .eq. XLO .or. face .eq. XHI) then
         do k = ARG_L3(out), ARG_H3(out)
            do j= ARG_L2(out),ARG_H2(out)
               do i = ARG_L1(out),ARG_H1(out)
                  out(i,j,k) = in(k,i,j)
               enddo
            enddo
         enddo
      else if (face .eq. YLO .or. face .eq. YHI) then
         do k = ARG_L3(out), ARG_H3(out)
            do j= ARG_L2(out),ARG_H2(out)
               do i = ARG_L1(out),ARG_H1(out)
                  out(i,j,k) = in(i,k,j)
               enddo
            enddo
         enddo
      else if (face .eq. ZLO .or. face .eq. ZHI) then
         do k = ARG_L3(out), ARG_H3(out)
            do j= ARG_L2(out),ARG_H2(out)
               do i = ARG_L1(out),ARG_H1(out)
                  out(i,j,k) = in(i,j,k)
               enddo
            enddo
         enddo
      endif

#undef XLO
#undef YLO
#undef ZLO
#undef XHI
#undef YHI
#undef ZHI

      end


      subroutine FORT_MAC_RESHIFT_PHI(DIMS(out),out,DIMS(in),in,face)
      implicit none
      integer face
      integer DIMDEC(in)
      integer DIMDEC(out)
      REAL_T in(DIMV(in))
      REAL_T out(DIMV(out))
      integer i,j,k
#define XLO 0
#define YLO 1
#define ZLO 2
#define XHI 3
#define YHI 4
#define ZHI 5
      if (face .eq. XLO .or. face .eq. XHI) then
         do k = ARG_L3(out), ARG_H3(out)
            do j= ARG_L2(out),ARG_H2(out)
               do i = ARG_L1(out),ARG_H1(out)
                  out(i,j,k) = in(j,k,i)
               enddo
            enddo
         enddo
      else if (face .eq. YLO .or. face .eq. YHI) then
         do k = ARG_L3(out), ARG_H3(out)
            do j= ARG_L2(out),ARG_H2(out)
               do i = ARG_L1(out),ARG_H1(out)
                  out(i,j,k) = in(i,k,j)
               enddo
            enddo
         enddo
      else if (face .eq. ZLO .or. face .eq. ZHI) then
         do k = ARG_L3(out), ARG_H3(out)
            do j= ARG_L2(out),ARG_H2(out)
               do i = ARG_L1(out),ARG_H1(out)
                  out(i,j,k) = in(i,j,k)
               enddo
            enddo
         enddo
      endif

#undef XLO
#undef YLO
#undef ZLO
#undef XHI
#undef YHI
#undef ZHI

      end

      subroutine FORT_SOLVEMAC(p, DIMS(p),dest0, DIMS(dest0),
     &                         source,DIMS(source), sigmax, DIMS(sigmax),
     &                         sigmay,DIMS(sigmay), cen, DIMS(cen),
     $                         r,DIMS(r), w, DIMS(w),z, DIMS(z),
     &                         x, DIMS(x),lo, hi, h,
     $                         isPeriodic, maxiter,tol, abs_tol,max_jump,norm)

      implicit none

      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(p)
      integer DIMDEC(dest0)
      integer DIMDEC(source)
      integer DIMDEC(sigmax)
      integer DIMDEC(sigmay)
      integer DIMDEC(cen)
      integer DIMDEC(r)
      integer DIMDEC(w)
      integer DIMDEC(z)
      integer DIMDEC(x)
      REAL_T      p(DIM12(p))
      REAL_T  dest0(DIM12(dest0))
      REAL_T source(DIM12(source))
      REAL_T sigmax(DIM12(sigmax))
      REAL_T sigmay(DIM12(sigmay))
      REAL_T    cen(DIM12(cen))
      REAL_T      r(DIM12(r))
      REAL_T      w(DIM12(w))
      REAL_T      z(DIM12(z))
      REAL_T      x(DIM12(x))
      REAL_T h(SDIM)
      integer isPeriodic(SDIM)
      REAL_T norm
      REAL_T tol
      integer maxiter
      REAL_T abs_tol,max_jump

c     Local variables
      integer i,j,iter
      integer istart,iend,jstart,jend
      REAL_T alpha, beta, rho, rho_old
      REAL_T  goal
      REAL_T norm0
      logical setSingularPoint

      istart = lo(1)
      iend = hi(1)
      jstart = lo(2)
      jend = hi(2)

      setSingularPoint = .false.

      do j = lo(2)-1,hi(2)+1 
        do i = lo(1)-1,hi(1)+1 
           dest0(i,j) = p(i,j)
           p(i,j) = zero
        enddo
      enddo

      do j=ARG_L2(w),ARG_H2(w)
         do i=ARG_L1(w),ARG_H1(w)
            w(i,j) = zero
         enddo
      enddo

      call setmacbc(DIMS(dest0),dest0,lo,hi,isPeriodic,setSingularPoint)
      call makemacdgphi(dest0,DIMS(dest0),w,DIMS(w),sigmax,DIMS(sigmax),
     &               sigmay,DIMS(sigmay),lo,hi,h,isPeriodic,setSingularPoint)

      do j = jstart, jend 
        do i = istart, iend 
          r(i,j) = source(i,j) - w(i,j)
        enddo
      enddo

      rho = zero
      norm0 = zero
      do j = jstart, jend 
        do i = istart, iend 
          norm0 = max(norm0, abs(r(i,j)))
          z(i,j) = r(i,j) 
          rho = rho + z(i,j) * r(i,j)
        enddo
      enddo
      norm = norm0

      goal = max(tol*norm0,abs_tol)

      if (norm0 .le. goal) then
         do j = jstart, jend
            do i = istart, iend 
               p(i,j) = dest0(i,j)
            enddo
         enddo
         return
      endif

      do j = jstart, jend 
        do i = istart, iend 
          x(i,j) = zero
          p(i,j) = z(i,j)
        enddo
      enddo

      iter = 0

100   continue

      do j=ARG_L2(w),ARG_H2(w)
         do i=ARG_L1(w),ARG_H1(w)
            w(i,j) = zero
         enddo
      enddo

      call setmacbc(DIMS(p),p,lo,hi,isPeriodic,setSingularPoint)
      call makemacdgphi(p,DIMS(p),w,DIMS(w),sigmax,DIMS(sigmax),
     &               sigmay,DIMS(sigmay),lo,hi,h,isPeriodic,setSingularPoint)

      alpha = zero
      do j = jstart, jend 
        do i = istart, iend 
          alpha = alpha + p(i,j)*w(i,j)
        enddo
      enddo

      alpha = rho / alpha
      rho_old = rho
      rho = zero
      norm = zero
      do j = jstart, jend 
        do i = istart, iend 
          x(i,j) = x(i,j) + alpha * p(i,j)
          r(i,j) = r(i,j) - alpha * w(i,j)
          z(i,j) = r(i,j) 
          rho = rho + z(i,j) * r(i,j)
          norm = max(norm,abs(r(i,j)))
        enddo
      enddo

      iter = iter+1
c      write(6,*) iter,norm

      if (iter .gt. maxiter .or. norm .gt. max_jump*norm0) then
         
         print *, "cg solve in mac failed to converge"
         do j = jstart, jend 
            do i = istart, iend 
               p(i,j) = dest0(i,j)
            enddo
         enddo
         call setmacbc(DIMS(p),p,lo,hi,isPeriodic,setSingularPoint)

      else if (norm .lt. goal) then

         do j = jstart, jend 
            do i = istart, iend 
               p(i,j) = x(i,j) + dest0(i,j)
            enddo
         enddo
         call setmacbc(DIMS(p),p,lo,hi,isPeriodic,setSingularPoint)
         
      else

        beta = rho / rho_old
        do j = jstart, jend 
           do i = istart, iend 
              p(i,j) = z(i,j) + beta * p(i,j)
           enddo
        enddo

        goto 100

      endif

      return
      end

      subroutine makemacdgphi(phi,DIMS(phi),dgphi,DIMS(dgphi),
     &                     beta0,DIMS(beta0),beta1,DIMS(beta1),
     &                     lo,hi,h,isPeriodic,setSingularPoint)
      implicit none

      integer DIMDEC(phi)
      integer DIMDEC(dgphi)
      integer DIMDEC(beta0)
      integer DIMDEC(beta1)
      REAL_T phi(DIM12(phi))
      REAL_T dgphi(DIM12(dgphi))
      REAL_T beta0(DIM12(beta0))
      REAL_T beta1(DIM12(beta1))
      integer lo(SDIM),hi(SDIM)
      REAL_T h(SDIM)
      integer isPeriodic(SDIM)
      logical setSingularPoint

      integer i,j
      REAL_T hxsqinv, hysqinv

      hxsqinv = one/(h(1)*h(1))
      hysqinv = one/(h(2)*h(2))

      do j = lo(2),hi(2)
        do i = lo(1),hi(1) 
          dgphi(i,j) = DGXY
        enddo
      enddo

      end



c *************************************************************************
c ** COARSIGMA **
c ** Coarsen the edge-based sigma coefficients
c *************************************************************************

      subroutine FORT_COARSIGMA(sigmax,DIMS(sigmax),sigmay,DIMS(sigmay),
     &                          sigmaxc,DIMS(sigmaxc),sigmayc,DIMS(sigmayc),
     &                          lo,hi,loc,hic)

      implicit none
      integer lo(SDIM),hi(SDIM)
      integer loc(SDIM),hic(SDIM)
      integer DIMDEC(sigmax)
      integer DIMDEC(sigmay)
      integer DIMDEC(sigmaxc)
      integer DIMDEC(sigmayc)
      REAL_T  sigmax(DIM12(sigmax))
      REAL_T  sigmay(DIM12(sigmay))
      REAL_T sigmaxc(DIM12(sigmaxc))
      REAL_T sigmayc(DIM12(sigmayc))

c     Local variables
      integer i,j,twoi,twoj

      do j = loc(2),hic(2) 
        do i = loc(1),hic(1)+1 
          twoi = 2*(i-loc(1))+lo(1)
          twoj = 2*(j-loc(2))+lo(2)
          sigmaxc(i,j) = half*(sigmax(twoi,twoj) + sigmax(twoi,twoj+1))
        enddo
      enddo

      do j = loc(2),hic(2)+1 
        do i = loc(1),hic(1) 
          twoi = 2*(i-loc(1))+lo(1)
          twoj = 2*(j-loc(2))+lo(2)
          sigmayc(i,j) = half*(sigmay(twoi,twoj) + sigmay(twoi+1,twoj))
        enddo
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
      REAL_T  res(DIM12(res))
      REAL_T resc(DIM12(resc))

c     Local variables
      integer i,j,twoi,twoj

c ::: NOTE: dont need factor of r here for volume-weighting because
c ::: what were calling the residual is really already r*residual

      do j = loc(2),hic(2) 
        do i = loc(1),hic(1) 
          twoi = 2*(i-loc(1))+lo(1)
          twoj = 2*(j-loc(2))+lo(2)
          resc(i,j) = ((res(twoi  ,twoj) + res(twoi  ,twoj+1)) + 
     $                 (res(twoi+1,twoj) + res(twoi+1,twoj+1)))*fourth
        enddo
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
      REAL_T    phi(DIM12(phi))
      REAL_T deltac(DIM12(deltac))

c     Local variables
      integer i,j,twoi,twoj

      do j = loc(2), hic(2) 
        do i = loc(1), hic(1) 

          twoi = 2*(i-loc(1))+lo(1)
          twoj = 2*(j-loc(2))+lo(2)
          phi(twoi  ,twoj  ) = phi(twoi  ,twoj  ) + deltac(i,j)
          phi(twoi+1,twoj  ) = phi(twoi+1,twoj  ) + deltac(i,j)
          phi(twoi  ,twoj+1) = phi(twoi  ,twoj+1) + deltac(i,j)
          phi(twoi+1,twoj+1) = phi(twoi+1,twoj+1) + deltac(i,j)

        enddo
      enddo

      return
      end

c *************************************************************************
c ** MACFILL_TWOD **
c *************************************************************************

      subroutine FORT_MACFILL_TWOD(lenx,leny,length,width,faces,numOutFlowFaces,
     $                             cc0,cc1,cc2,cc3,cc4,cc5,conn)

      integer lenx,leny,length,width
      integer faces(4)
      integer numOutFlowFaces
      REAL_T cc0(lenx,leny,2)
      REAL_T cc1(lenx,leny,2)
      REAL_T cc2(lenx,leny,2)
      REAL_T cc3(lenx,leny,2)
      REAL_T cc4(lenx,leny,2)
      REAL_T cc5(lenx,leny,2)
      REAL_T conn(length,2)

      integer xlo_outflow,ylo_outflow,zlo_outflow
      integer xhi_outflow,yhi_outflow,zhi_outflow
      integer i,ifinal,n
      REAL_T sum 

#define XLO 0
#define YLO 1
#define ZLO 2
#define XHI 3
#define YHI 4
#define ZHI 5

c     Want to find any non-outflow faces.
      xlo_outflow = 0
      ylo_outflow = 0
      zlo_outflow = 0
      xhi_outflow = 0
      yhi_outflow = 0
      zhi_outflow = 0

      do i = 1, numOutFlowFaces
        if (faces(i) .eq. XLO) xlo_outflow = 1
        if (faces(i) .eq. YLO) ylo_outflow = 1
        if (faces(i) .eq. ZLO) zlo_outflow = 1
        if (faces(i) .eq. XHI) xhi_outflow = 1
        if (faces(i) .eq. YHI) yhi_outflow = 1
        if (faces(i) .eq. ZHI) zhi_outflow = 1
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

       end

