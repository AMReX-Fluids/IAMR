#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "PROJOUTFLOWBC_F.H"
#include "ArrayLim.H"

#define SDIM 3

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define SMALL 1.0e-10
#define sixteenth  .0625e0
#else
#define SMALL 1.0d-10
#define sixteenth  .0625d0
#endif

      subroutine FORT_EXTRAP_PROJ(DIMS(u),u,DIMS(divu),divu,DIMS(rho),rho,
     &     DIMS(uExt),uExt,DIMS(divuExt),divuExt,
     &     DIMS(rhoExt),rhoExt,lo,hi,face, zeroIt)
      implicit none

c    compute divu_ave twice due to precision problems

      integer DIMDEC(u)
      integer DIMDEC(divu)
      integer DIMDEC(rho)
      integer DIMDEC(uExt)
      integer DIMDEC(divuExt)
      integer DIMDEC(rhoExt)
      integer face
      integer lo(SDIM),hi(SDIM)
      REAL_T      u(DIMV(u),SDIM)
      REAL_T   divu(DIMV(divu))
      REAL_T    rho(DIMV(rho))
      REAL_T   uExt(DIMV(uExt),SDIM-1)
      REAL_T   divuExt(DIMV(divuExt))
      REAL_T   rhoExt(DIMV(rhoExt))
      integer  zeroIt

c local variables
      integer ics,ice,jcs,jce,kcs,kce
      integer ife,jfe,kfe
      integer if,jf,kf
      REAL_T divu_ave1,divu_ave2
      REAL_T max_divu
      REAL_T max_pert, small_pert
      parameter ( small_pert = SMALL)
      integer i,j,k

#define XLO 0
#define YLO 1
#define ZLO 2
#define XHI 3
#define YHI 4
#define ZHI 5

      ics = ARG_L1(u)
      ice = ARG_H1(u)
      jcs = ARG_L2(u)
      jce = ARG_H2(u)
      kcs = ARG_L3(u)
      kce = ARG_H3(u)

      ife = hi(1)
      jfe = hi(2)
      kfe = hi(3)

      zeroIt = 0

      if (face .eq. XLO) then
         if=ife
         do k = kcs+1, kce-1
         do j = jcs+1, jce-1
            uExt(j,k,if,1)  = half*(three*u(ice-1,j,k,2)   -  u(ice,j,k,2))
            uExt(j,k,if,2)  = half*(three*u(ice-1,j,k,3)   -  u(ice,j,k,3))
            divuExt(j,k,if) = half*(three*divu(ice-1,j,k)  - divu(ice,j,k))
            rhoExt(j,k,if)  = half*(three*rho(ice-1,j,k)   -  rho(ice,j,k))
         end do
         end do
         max_divu = ABS(divuExt(jcs+1,kcs+1,if))
         do k = kcs+1, kce-1
         do j = jcs+1, jce-1
            max_divu = max(max_divu,ABS(divuExt(j,k,if)))
         end do
         end do
         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(jcs+1,kcs+1,if))
         do k = kcs+1, kce-1
         do j = jcs+1, jce-1
            max_pert = MAX(max_pert,ABS(divuExt(j,k,if)))
         end do
         end do
      else if (face .eq. YLO) then
         jf = jfe
         do k = kcs+1, kce-1
         do i = ics+1, ice-1
            uExt(i,k,jf,1)    = half*(three*u(i,jce-1,k,1)    - u(i,jce,k,1))
            uExt(i,k,jf,2)    = half*(three*u(i,jce-1,k,3)    - u(i,jce,k,3))
            divuExt(i,k,jf) = half*(three*divu(i,jce-1,k) - divu(i,jce,k))
            rhoExt(i,k,jf)  = half*(three*rho(i,jce-1,k)    - rho(i,jce,k))
         end do
         end do
         max_divu = ABS(divuExt(ics+1,kcs+1,jf))
         do k = kcs+1, kce-1
         do i = ics+1, ice-1
            max_divu = max(max_divu,ABS(divuExt(i,k,jf)))
         end do
         end do
         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(ics+1,kcs+1,jf))
         do k = kcs+1, kce-1
         do i = ics+1, ice-1
            max_pert = MAX(max_pert,ABS(divuExt(i,k,jf)))
         end do
         end do
      else if (face .eq. ZLO) then
         kf = kfe
         do j = jcs+1, jce-1
         do i = ics+1, ice-1
            uExt(i,j,kf,1)    = half*(three*u(i,j,kce-1,2)    - u(i,j,kce,2))
            uExt(i,j,kf,2)    = half*(three*u(i,j,kce-1,3)    - u(i,j,kce,3))
            divuExt(i,j,kf) = half*(three*divu(i,j,kce-1) - divu(i,j,kce))
            rhoExt(i,j,kf)  = half*(three*rho(i,j,kce-1)    - rho(i,j,kce))
         end do
         end do
         max_divu = ABS(divuExt(ics+1,jcs+1,kf))
         do j = jcs+1, jce-1
         do i = ics+1, ice-1
            max_divu = max(max_divu,ABS(divuExt(i,j,kf)))
         end do
         end do
         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(ics+1,jcs+1,kf))
         do j = jcs+1, jce-1
         do i = ics+1, ice-1
            max_pert = MAX(max_pert,ABS(divuExt(i,j,kf)))
         end do
         end do
      else if (face .eq. XHI) then
         if = ife
         do k = kcs+1, kce-1
         do j = jcs+1, jce-1
            uExt(j,k,if,1)    = half*(three*u(ics+1,j,k,2)    - u(ics,j,k,2))
            uExt(j,k,if,2)    = half*(three*u(ics+1,j,k,3)    - u(ics,j,k,3))
            divuExt(j,k,if) = half*(three*divu(ics+1,j,k) - divu(ics,j,k))
            rhoExt(j,k,if)  = half*(three*rho(ics+1,j,k)    - rho(ics,j,k))
         end do
         end do
         max_divu = ABS(divuExt(jcs+1,kcs+1,if))
         do k = kcs+1, kce-1
         do j = jcs+1, jce-1
            max_divu = max(max_divu,ABS(divuExt(j,k,if)))
         end do
         end do
         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(jcs+1,kcs+1,if))
         do k = kcs+1, kce-1
         do j = jcs+1, jce-1
            max_pert = MAX(max_pert,ABS(divuExt(j,k,if)))
         end do
         end do
      else if (face .eq. YHI) then
         jf = jfe
         do k = kcs+1, kce-1
         do i = ics+1, ice-1
            uExt(i,k,jf,1)    = half*(three*u(i,jcs+1,k,1)    - u(i,jcs,k,1))
            uExt(i,k,jf,2)    = half*(three*u(i,jcs+1,k,3)    - u(i,jcs,k,3))
            divuExt(i,k,jf) = half*(three*divu(i,jcs+1,k) - divu(i,jcs,k))
            rhoExt(i,k,jf)  = half*(three*rho(i,jcs+1,k)    - rho(i,jcs,k))
         end do
         end do
         max_divu = ABS(divuExt(ics+1,kcs+1,jf))
         do k = kcs+1, kce-1
         do i = ics+1, ice-1
            max_divu = max(max_divu,ABS(divuExt(i,k,jf)))
         end do
         end do
         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(ics+1,kcs+1,jf))
         do k = kcs+1, kce-1
         do i = ics+1, ice-1
            max_pert = MAX(max_pert,ABS(divuExt(i,k,jf)))
         end do
         end do
      else if (face .eq. ZHI) then
         kf = kfe
         do j = jcs+1, jce-1
         do i = ics+1, ice-1
            uExt(i,j,kf,1)    = half*(three*u(i,j,kcs+1,1)    - u(i,j,kcs,1))
            uExt(i,j,kf,2)    = half*(three*u(i,j,kcs+1,2)    - u(i,j,kcs,2))
            divuExt(i,j,kf) = half*(three*divu(i,j,kcs+1) - divu(i,j,kcs))
            rhoExt(i,j,kf)  = half*(three*rho(i,j,kcs+1)    - rho(i,j,kcs))
         end do
         end do
         max_divu = ABS(divuExt(ics+1,jcs+1,kf))
         do j = jcs+1, jce-1
         do i = ics+1, ice-1
            max_divu = max(max_divu,ABS(divuExt(i,j,kf)))
         end do
         end do
         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(ics+1,jcs+1,kf))
         do j = jcs+1, jce-1
         do i = ics+1, ice-1
            uExt(i,j,kf,1)    = half*(three*u(i,j,kcs+1,1)    - u(i,j,kcs,1))
            max_pert = MAX(max_pert,ABS(divuExt(i,j,kf)))
         end do
         end do
      endif

c  check to see if we should zero phi
      max_pert = max_pert/(ABS(divu_ave1+divu_ave2)+small_pert)
      if ((max_divu.eq.zero) .or. (max_pert.le.small_pert)) zeroIt = 1
#undef XLO
#undef YLO
#undef ZLO
#undef XHI
#undef YHI
#undef ZHI
      end

      subroutine FORT_COMPUTE_COEFF(DIMS(rhs),rhs,DIMS(beta),beta,
     &                              DIMS(uExt),uExt,DIMS(divuExt),divuExt,
     &                              DIMS(rhoExt),rhoExt,
     &                              lo,hi,h,isPeriodic)
      implicit none
      integer DIMDEC(rhs)
      integer DIMDEC(beta)
      integer DIMDEC(uExt)
      integer DIMDEC(divuExt)
      integer DIMDEC(rhoExt)
      integer lo(SDIM),hi(SDIM)
      REAL_T uExt(DIM12(uExt),SDIM-1)
      REAL_T divuExt(DIM12(divuExt))
      REAL_T rhoExt(DIM12(rhoExt))
      REAL_T beta(DIM12(beta))
      REAL_T rhs(DIM12(rhs))
      REAL_T h(SDIM)
      integer isPeriodic(SDIM)

      integer ics,ice,jcs,jce
      integer ins,ine,jns,jne
      integer i,j
      REAL_T hxm1,hym1

      ics = lo(1)
      ice = hi(1)
      jcs = lo(2)
      jce = hi(2)
      ins = ics
      jns = jcs
      ine = ice+1
      jne = jce+1

      hxm1 = one/h(1)
      hym1 = one/h(2)
      
      do j=jcs,jce
         do i=ics,ice
            beta(i,j) = one/rhoExt(i,j)
         enddo
      enddo

      do j=jns+1,jne-1
         do i=ins+1,ine-1
            rhs(i,j) = half*
     @           (hxm1 * (uExt(i,j-1,1) - uExt(i-1,j-1,1) +
     @                    uExt(i,j  ,1) - uExt(i-1,j  ,1)) +
     @            hym1 * (uExt(i-1,j,2) - uExt(i-1,j-1,2) +
     @                    uExt(i  ,j,2) - uExt(i,j-1  ,2))) 
     &           - fourth * ( divuExt(i-1,j  ) + divuExt(i-1,j-1) +
     &                        divuExt(i  ,j-1) + divuExt(i  ,j  ) ) 
         enddo
      enddo

      if (isPeriodic(1) .eq. 1) then
         do j=jcs,jce
            beta(ics-1,j) = beta(ice,j)
            beta(ice+1,j) = beta(ics,j)
         enddo 
         do j = jns+1,jne-1
            rhs(ins,j) = half*
     @           (hxm1 * (uExt(ics,j-1,1) - uExt(ice,j-1,1) +
     @                    uExt(ics,j  ,1) - uExt(ice,j,1  )) +
     @            hym1 * (uExt(ice,j  ,2) - uExt(ice,j-1,2) +
     @                    uExt(ics,j  ,2) - uExt(ics,j-1,2))) 
     &           - fourth * ( divuExt(ice,j  ) + divuExt(ice,j-1) +
     &                        divuExt(ics,j-1) + divuExt(ics  ,j  ) ) 

            rhs(ine,j) = rhs(ins,j)

            
         enddo
      else
         do j=jcs,jce
            beta(ics-1,j) = zero
            beta(ice+1,j) = zero
         enddo
         do j = jns+1,jne-1
            i = ins
            rhs(i,j) = half*
     @           (hxm1 * (uExt(i,j-1,1) + uExt(i,j  ,1) ) +
     @            hym1 * (uExt(i  ,j,2) - uExt(i,j-1,2))) 
     &           - fourth * (divuExt(i,j-1) + divuExt(i,j) ) 

            i = ine
            rhs(i,j) = half*
     @           (hxm1 * (- uExt(i-1,j-1,1) - uExt(i-1,j  ,1)) +
     @            hym1 * (  uExt(i-1,j  ,2) - uExt(i-1,j-1,2))) 
     &        - fourth * ( divuExt(i-1,j  ) + divuExt(i-1,j-1) ) 

         enddo
      endif

      if (isPeriodic(2) .eq. 1) then
         do i=ics,ice
            beta(i,jcs-1) = beta(i,jce)
            beta(i,jce+1) = beta(i,jcs)
         enddo
         do i = ins+1,ine-1
            rhs(i,jns) = half*
     @           (hxm1 * (uExt(i  ,jce,1) - uExt(i-1,jce,1) +
     @                    uExt(i  ,jcs,1) - uExt(i-1,jcs,1))+
     @            hym1 * (uExt(i-1,jcs,2) - uExt(i-1,jce,2) +
     @                    uExt(i  ,jcs,2) - uExt(i  ,jce,2))) 
     &           - fourth * ( divuExt(i-1,jcs) + divuExt(i-1,jce) +
     &                        divuExt(i  ,jce) + divuExt(i  ,jcs) ) 

            rhs(i,jne) = rhs(i,jns)
         enddo
      else
         do i=ics,ice
            beta(i,jcs-1) = zero
            beta(i,jce+1) = zero
         enddo
         do i = ins+1,ine-1
            j = jns
            rhs(i,j) = half*
     @           (hxm1 * (uExt(i,j  ,1) - uExt(i-1,j  ,1)) +
     @            hym1 * (uExt(i-1,j,2) + uExt(i  ,j,2) )) 
     &           - fourth * ( divuExt(i-1,j) + divuExt(i,j) )

            j = jne
            rhs(i,j) = half*
     @           (hxm1 * ( uExt(i  ,j-1,1) - uExt(i-1,j-1,1)) +
     @            hym1 * (-uExt(i-1,j-1,2) - uExt(i  ,j-1,2))) 
     &           - fourth * ( divuExt(i-1,j-1) + divuExt(i  ,j-1)) 
         enddo
      endif

      if (isPeriodic(1) .eq. 1 .AND. isPeriodic(2) .eq. 1) then

         beta(ics-1,jcs-1) = beta(ice,jce)
         beta(ics-1,jce+1) = beta(ice,jcs)
         beta(ice+1,jcs-1) = beta(ics,jce)
         beta(ice+1,jce+1) = beta(ics,jcs)

         rhs(ins,jns) = half*
     @        (hxm1 * (uExt(ics,jce,1) - uExt(ice,jce,1) +
     @                 uExt(ics,jcs,1) - uExt(ice,jcs,1)) +
     @         hym1 * (uExt(ice,jcs,2) - uExt(ice,jce,2) +
     @                 uExt(ics,jcs,2) - uExt(ics,jce,2))) 
     &        - fourth * (divuExt(ice,jcs) + divuExt(ice,jce) +
     &                    divuExt(ics,jce) + divuExt(ics,jcs) ) 

         rhs(ins,jne) = rhs(ins,jns) 
         rhs(ine,jns) = rhs(ins,jns) 
         rhs(ine,jne) = rhs(ins,jns) 

      else if (isPeriodic(1) .eq. 1 .AND. isPeriodic(2) .NE. 1) then

         beta(ics-1,jcs-1) = beta(ice,jcs-1)
         beta(ics-1,jce+1) = beta(ice,jce+1)
         beta(ice+1,jcs-1) = beta(ics,jcs-1)
         beta(ice+1,jce+1) = beta(ics,jce+1)

         rhs(ins,jns) = half*
     @        (hxm1 * (uExt(ics,jcs,1) - uExt(ice,jcs,1  )) +
     @         hym1 * (uExt(ice,jcs,2) + uExt(ics,jcs,2) )) 
     &        - fourth * ( divuExt(ice,jcs) + divuExt(ics,jcs) ) 
         
         rhs(ins,jne) = half*
     @        (hxm1 * (  uExt(ics,jce,1) - uExt(ice,jce,1) ) +
     @         hym1 * (- uExt(ice,jce,2) - uExt(ics,jce,2))) 
     &        - fourth * ( divuExt(ice,jce) + divuExt(ics,jce)) 
         
         rhs(ine,jns) = rhs(ins,jns)
         rhs(ine,jne) = rhs(ins,jne)

      else if (isPeriodic(1) .NE. 1 .AND. isPeriodic(2) .eq. 1) then

         beta(ics-1,jcs-1) = beta(ics-1,jce)
         beta(ics-1,jce+1) = beta(ics-1,jcs)
         beta(ice+1,jcs-1) = beta(ice+1,jce)
         beta(ice+1,jce+1) = beta(ice+1,jcs)

         rhs(ins,jns) = half*
     @        (hxm1 * (uExt(ics,jce,1) + uExt(ics,jcs,1))+
     @         hym1 * (uExt(ics,jcs,2) - uExt(ics,jce,2))) 
     &      - fourth * ( divuExt(ics,jce) + divuExt(ics,jcs) ) 
         
         rhs(ine,jns) = half*
     @        (hxm1 * (uExt(ice,jce,1) - uExt(ice,jcs,1))+
     @         hym1 * (uExt(ice,jcs,2) - uExt(ice,jce,2))) 
     &        - fourth * ( divuExt(ice,jcs) + divuExt(ice,jce)) 
         
         rhs(ins,jne) = rhs(ins,jns)
         rhs(ine,jne) = rhs(ine,jns)

      else

         beta(ics-1,jcs-1) = zero
         beta(ics-1,jce+1) = zero
         beta(ice+1,jcs-1) = zero
         beta(ice+1,jce+1) = zero

         i = ins
         j = jns
         rhs(i,j) = half* (hxm1 * (uExt(i,j,1) ) + hym1 * (uExt(i,j,2) )) 
     &        - fourth * ( divuExt(i,j) )

         
         i = ine
         j = jns
         rhs(i,j) = half* (hxm1 * (-uExt(i-1,j,1))+ hym1 * (uExt(i-1,j,2))) 
     &        - fourth * ( divuExt(i-1,j))
         
         i = ins
         j = jne
         rhs(i,j) = half* (hxm1 * (uExt(i,j-1,1) ) + hym1 * (-uExt(i,j-1,2))) 
     &        - fourth * (divuExt(i  ,j-1)) 

         i = ine
         j = jne
         rhs(i,j) = half*(hxm1 *(-uExt(i-1,j-1,1)) + hym1 *(-uExt(i-1,j-1,2))) 
     &        - fourth * ( divuExt(i-1,j-1)) 

      endif

c  double rhs at edges
      if (isPeriodic(1) .ne. 1) then
         do j=jns,jne
            rhs(ins,j) = rhs(ins,j) * two
            rhs(ine,j) = rhs(ine,j) * two
         enddo
      endif
      if (isPeriodic(2) .ne. 1) then
         do i = ins,ine
            rhs(i,jns) = rhs(i,jns) * two
            rhs(i,jne) = rhs(i,jne) * two
         enddo
      endif
 

      end

#define DGXY_5PT    half * (hxsqinv * \
           ((beta(i-1,j-1) + beta(i-1,j)) * (phi(i-1,j) - phi(i,j)) + \
            (beta(i  ,j-1) + beta(i  ,j)) * (phi(i+1,j) - phi(i,j))) + \
            hysqinv * \
           ((beta(i-1,j-1) + beta(i,j-1)) *(phi(i,j-1) - phi(i,j)) + \
            (beta(i-1,j  ) + beta(i,j  )) *(phi(i,j+1) - phi(i,j))))


      subroutine FORT_HGRELAX(DIMS(rhs),rhs,DIMS(beta),beta,DIMS(phi),phi,
     &                        DIMS(dgphi),dgphi,lo,hi,h,isPeriodic,niter)
      implicit none
      integer DIMDEC(beta)
      integer DIMDEC(rhs)
      integer DIMDEC(phi)
      integer DIMDEC(dgphi)
      REAL_T beta(DIM12(beta))
      REAL_T rhs(DIM12(rhs))
      REAL_T phi(DIM12(phi))
      REAL_T dgphi(DIM12(dgphi))
      integer isPeriodic(SDIM)
      REAL_T h(SDIM)
      integer lo(SDIM),hi(SDIM)
      integer niter

      logical setSingularPoint
      REAL_T lam
      integer i,j
      integer redblack,iter
      integer ins,ine,jns,jne
      REAL_T hxsqinv, hysqinv
      integer iinc

      ins = lo(1)
      ine = hi(1)+1
      jns = lo(2)
      jne = hi(2)+1
      hxsqinv = one/(h(1)*h(1))
      hysqinv = one/(h(2)*h(2))

      if (h(2). gt. 1.5D0*h(1)) then
         call bl_abort("line solve for proj_bc not yet implemented")
      elseif (h(1) .gt. 1.5D0*h(2)) then
         call bl_abort("line solve for proj_bc not yet implemented")
      endif
      
      setSingularPoint = .false.

      do iter = 1,niter
         call setprojbc(DIMS(phi),phi,lo,hi,isPeriodic,setSingularPoint)
         do redblack = 0,1
            do j=jns,jne
               iinc = mod(j+redblack,2)
               do i=ins+iinc,ine,2

                  dgphi(i,j) = DGXY_5PT
                  lam =  (hxsqinv +hysqinv) * (beta(i-1,j-1) + beta(i-1,j) +
     @                              beta(i  ,j-1) + beta(i  ,j))
c double dgphi at edges
                  if ((i .eq. lo(1) .or. i .eq. hi(1)+1) 
     &                 .and. isPeriodic(1) .ne. 1) then
                     dgphi(i,j) = dgphi(i,j)*two
                     lam = lam* two
                  endif
                  if ((j .eq. lo(2) .or. j .eq. hi(2)+1) 
     &                 .and. isPeriodic(2) .ne. 1) then
                     dgphi(i,j) = dgphi(i,j)*two
                     lam = lam*two
                  endif
                  lam =  -two/lam
                  phi(i,j) = phi(i,j) + lam*(rhs(i,j)-dgphi(i,j))
               enddo
            enddo
         end do
      end do

      call setprojbc(DIMS(phi),phi,lo,hi,isPeriodic,setSingularPoint)

      end
      
      subroutine FORT_HGSUBTRACTAVGPHI(DIMS(phi),phi,lo,hi,isPeriodic)
      implicit none
      integer DIMDEC(phi)
      REAL_T phi(DIM12(phi))
      integer lo(SDIM),hi(SDIM)
      integer isPeriodic(SDIM)

      REAL_T phitot,vtot
      integer ins,ine,jns,jne
      integer i,j
      logical setSingularPoint

      phitot = zero
      vtot = zero
      ins = lo(1)
      ine = hi(1)+1
      jns = lo(2)
      jne = hi(2)+1
      setSingularPoint = .false.

      do j = jns,jne
         do i = ins,ine
            phitot = phitot + phi(i,j) 
            vtot   = vtot + one
         enddo
      enddo
      phitot = phitot/vtot
      do j = jns,jne
         do i = ins,ine
            phi(i,j) = phi(i,j) - phitot
         enddo
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
      REAL_T beta(DIM12(beta))
      REAL_T rhs(DIM12(rhs))
      REAL_T phi(DIM12(phi))
      REAL_T resid(DIM12(resid))
      REAL_T dgphi(DIM12(dgphi))
      integer isPeriodic(SDIM)
      REAL_T h(SDIM)
      integer lo(SDIM),hi(SDIM)
      REAL_T hxsqinv,hysqinv
      
      integer i,j
      REAL_T maxnorm
      logical setSingularPoint

      hxsqinv = one/(h(1)*h(1))
      hysqinv = one/(h(2)*h(2))
      setSingularPoint = .false.
      maxnorm = zero

      call makeprojdgphi(phi,DIMS(phi),dgphi,DIMS(dgphi),beta,DIMS(beta),
     &                   lo,hi,h,isPeriodic,setSingularPoint)
      do j=lo(2),hi(2)+1
        do i=lo(1),hi(1)+1
           resid(i,j) = rhs(i,j)-dgphi(i,j)
           maxnorm = max(maxnorm,ABS(resid(i,j)))         
        enddo
      enddo
      end

      subroutine setprojbc(DIMS(phi),phi,lo,hi,isPeriodic,setSingularPoint)
      implicit none
      integer DIMDEC(phi)
      REAL_T phi(DIM12(phi))
      integer isPeriodic(SDIM)
      logical setSingularPoint
      integer lo(SDIM),hi(SDIM)

      integer ins,ine,jns,jne
      integer i,j
    
      ins = lo(1)
      ine = hi(1)+1
      jns = lo(2)
      jne = hi(2)+1


      if (isPeriodic(1) .NE. 1 .AND. isPeriodic(2) .NE. 1 
     &     .AND. setSingularPoint) then
         phi(ine,jne) = zero
      endif

      if (isPeriodic(1) .eq. 1) then
         do j= jns,jne
            phi(ins-1,j) = phi(ine-1,j)
            phi(ine,j)   = phi(ins,j)
            phi(ine+1,j) = phi(ins+1,j)
         enddo
      else
         do j= jns,jne
            phi(ins-1,j) = phi(ins+1,j)
            phi(ine+1,j) = phi(ine-1,j)
         enddo
      endif

      if (isPeriodic(2) .eq. 1) then
         do i= ins,ine
            phi(i,jns-1) = phi(i,jne-1)
            phi(i,jne)   = phi(i,jns)
            phi(i,jne+1) = phi(i,jns+1)
         enddo
      else
         do i= ins,ine
            phi(i,jns-1) = phi(i,jns+1)
            phi(i,jne+1) = phi(i,jne-1)
         enddo
      endif

      end

      subroutine FORT_HG_SHIFT_PHI(DIMS(out),out,DIMS(in),in,face)
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

      subroutine FORT_HG_RESHIFT_PHI(DIMS(out),out,DIMS(in),in,face)
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
               out(ARG_L1(out),j,k) = in(j,k,ARG_L3(in))
            enddo
         enddo
      else if (face .eq. YLO .or. face .eq. YHI) then
         do k = ARG_L3(out), ARG_H3(out)
            do i = ARG_L1(out),ARG_H1(out)
               out(i,ARG_L2(out),k) = in(i,k,ARG_L3(in))
            enddo
         enddo
      else if (face .eq. ZLO .or. face .eq. ZHI) then
         do j= ARG_L2(out),ARG_H2(out)
            do i = ARG_L1(out),ARG_H1(out)
               out(i,j,ARG_L3(out)) = in(i,j,ARG_L3(in))
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
      REAL_T   p(DIM12(p))
      REAL_T  dest0(DIM12(dest0))
      REAL_T source(DIM12(source))
      REAL_T  sigma(DIM12(sigma))
      REAL_T      r(DIM12(r))
      REAL_T      w(DIM12(w))
      REAL_T      z(DIM12(z))
      REAL_T      x(DIM12(x))
      REAL_T      cen(DIM12(cen))
      REAL_T h(SDIM)
      integer isPeriodic(SDIM)
      integer maxiter
      REAL_T norm
      REAL_T tol
      REAL_T abs_tol,max_jump

c     Local variables
      REAL_T factor
      REAL_T  alpha,beta, rho, rho_old
      logical testx,testy
      integer i,j,iter
      integer istart,iend,jstart,jend
      REAL_T norm0,goal
      logical setSingularPoint

      istart = lo(1)
      iend = hi(1)+1
      jstart = lo(2)
      jend = hi(2)+1
      
      setSingularPoint = .false.

      do j = lo(2)-1,hi(2)+2
         do i = lo(1)-1,hi(1)+2
            dest0(i,j) = p(i,j)
            p(i,j) = zero
         enddo
      enddo

      do j=ARG_L2(w),ARG_H2(w)
         do i=ARG_L1(w),ARG_H1(w)
            w(i,j) = zero
         enddo
      enddo

      call setprojbc(DIMS(dest0),dest0,lo,hi,isPeriodic,setSingularPoint)
      call makeprojdgphi(dest0,DIMS(dest0),w,DIMS(w),sigma,DIMS(sigma),
     &        lo,hi,h,isPeriodic,setSingularPoint)

      do j = jstart, jend 
        do i = istart, iend 
          r(i,j) = source(i,j) - w(i,j)
        enddo
      enddo

c note that all of this factor stuff is due to the problem being doubled
c at edges -- both the rhs and the operator.
      rho = zero
      norm0 = zero
      do j = jstart, jend 
         do i = istart, iend
            factor = one
            testx  = (i .eq. lo(1) .or. i .eq. hi(1)+1)
            testy  = (j .eq. lo(2) .or. j .eq. hi(2)+1)
            factor = cvmgt(factor*half,factor,testx)
            factor = cvmgt(factor*half,factor,testy)
            z(i,j) = r(i,j) 
            rho    = rho + factor * z(i,j) * r(i,j) 
            norm0 = max(norm0,abs(r(i,j)))
         enddo
      enddo
      norm = norm0

      goal = max(tol * norm0,abs_tol)

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

      iter  = 0

100   continue  

      do j=ARG_L2(w),ARG_H2(w)
         do i=ARG_L1(w),ARG_H1(w)
            w(i,j) = zero
         enddo
      enddo

      call setprojbc(DIMS(p),p,lo,hi,isPeriodic,setSingularPoint)
      call makeprojdgphi(p,DIMS(p),w,DIMS(w),sigma,DIMS(sigma),
     &      lo,hi,h,isPeriodic,setSingularPoint)

      alpha = zero
      do j = jstart, jend 
        do i = istart, iend 
            factor = one
            testx  = (i .eq. lo(1) .or. i .eq. hi(1)+1)
            testy  = (j .eq. lo(2) .or. j .eq. hi(2)+1)
            factor = cvmgt(factor*half,factor,testx)
            factor = cvmgt(factor*half,factor,testy)
            alpha  = alpha + factor*p(i,j)*w(i,j) 
        enddo
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
      do j = jstart, jend 
         do i = istart, iend 
            factor = one
            testx  = (i .eq. lo(1) .or. i .eq. hi(1)+1)
            testy  = (j .eq. lo(2) .or. j .eq. hi(2)+1)
            factor = cvmgt(factor*half,factor,testx)
            factor = cvmgt(factor*half,factor,testy)
            x(i,j) = x(i,j) + alpha * p(i,j)
            r(i,j) = r(i,j) - alpha * w(i,j)
            z(i,j) = r(i,j) 
            rho    = rho + factor * z(i,j) * r(i,j) 
            norm = max(norm,abs(r(i,j)))
         enddo
      enddo

      iter = iter+1
c      write(6,*) iter,norm

      if (iter .gt. maxiter .or. norm .gt. max_jump*norm0) then

         goto 101

      else if (norm .lt. goal) then

         do j = jstart, jend 
            do i = istart, iend 
               p(i,j) = x(i,j) + dest0(i,j)
            enddo
         enddo
         
         return

      else
         
        beta = rho / rho_old
        do j = jstart, jend 
          do i = istart, iend 
            p(i,j) = z(i,j) + beta * p(i,j)
          enddo
        enddo

        goto 100

      endif

 101  print *, "cg solve in proj failed to coverge"
      do j = jstart, jend 
         do i = istart, iend 
            p(i,j) = dest0(i,j)
         enddo
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
      REAL_T    phi(DIM12(phi))
      REAL_T  dgphi(DIM12(dgphi))
      REAL_T  beta(DIM12(beta))
      integer isPeriodic(SDIM)
      REAL_T h(SDIM)
      logical setSingularPoint

c     Local variables
      REAL_T hxsqinv, hysqinv
      integer is,ie,js,je
      integer i,j

      is = lo(1)
      js = lo(2)
      ie = hi(1)
      je = hi(2)

      hxsqinv = one/(h(1)*h(1))
      hysqinv = one/(h(2)*h(2))
      
      if (isPeriodic(1) .eq. 1) then
        do j = js,je+1 
          phi(ie+2,j) = phi(is+1,j)
          phi(is-1,j) = phi(ie  ,j)
        enddo
      endif

      if (isPeriodic(2) .eq. 1) then
        do i = is,ie+1 
          phi(i,je+2) = phi(i,js+1)
          phi(i,js-1) = phi(i,je  )
        enddo
      endif

      if (isPeriodic(1) .eq. 1 .and. isPeriodic(2) .eq. 1) then
        phi(is-1,js-1) = phi(ie  ,je  )
        phi(is-1,je+2) = phi(ie  ,js+1)
        phi(ie+2,js-1) = phi(is+1,je  )
        phi(ie+2,je+2) = phi(is+1,js+1)
      endif    

      do j = js,je+1
        do i = is,ie+1
           dgphi(i,j) = DGXY_5PT
        enddo
      enddo

c double dgphi at edges
      if (isPeriodic(1) .ne. 1) then
         do j=js,je+1
            dgphi(is,j) = dgphi(is,j) * two
            dgphi(ie+1,j) = dgphi(ie+1,j) * two
         enddo
      endif
      if (isPeriodic(2) .ne. 1) then
         do i = is,ie+1
            dgphi(i,js) = dgphi(i,js) * two
            dgphi(i,je+1) = dgphi(i,je+1) * two
         enddo
      endif

      if (setSingularPoint .and. 
     &     isPeriodic(1) .NE. 1 .and. isPeriodic(2) .NE. 1) then
         dgphi(hi(1)+1,hi(2)+1) = zero
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
      REAL_T  sigma(DIM12(sigma))
      REAL_T sigmac(DIM12(sigmac))
      integer isPeriodic(SDIM)

c     Local variables
      integer i ,j
      integer i2,j2

      do j = loc(2),hic(2) 
        do i = loc(1),hic(1) 
          i2 = 2*(i-loc(1))+lo(1)
          j2 = 2*(j-loc(2))+lo(2)
          sigmac(i,j) = (sigma(i2  ,j2) + sigma(i2  ,j2+1)+ 
     $                   sigma(i2+1,j2) + sigma(i2+1,j2+1))*fourth
        enddo
      enddo

      if (isPeriodic(1) .eq. 1) then
         
         do j = loc(2)-1,hic(2)+1 
            sigmac(loc(1)-1,j) = sigmac(hic(1),j)
            sigmac(hic(1)+1,j) = sigmac(loc(1),j)
         enddo
      else
         do j = loc(2)-1,hic(2)+1
            sigmac(loc(1)-1,j) = zero
            sigmac(hic(1)+1,j) = zero
         enddo
      endif

      if (isPeriodic(2) .eq. 1) then
         
         do i = loc(1)-1,hic(1)+1 
            sigmac(i,loc(2)-1) = sigmac(i,hic(2))
            sigmac(i,hic(2)+1) = sigmac(i,loc(2))
         enddo
      else
         do i = loc(1)-1,hic(1)+1
            sigmac(i,loc(2)-1) = zero
            sigmac(i,hic(2)+1) = zero
         enddo
         
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
      REAL_T   res(DIM12(res))
      REAL_T  resc(DIM12(resc))
      integer isPeriodic(SDIM)

c     Local variables
      integer i,j,ii,jj
      integer istart,iend
      integer jstart,jend

      istart = loc(1)
      iend = hic(1)+1
      jstart = loc(2)
      jend = hic(2)+1
      
      if (isPeriodic(1) .eq. 1) then
        do j = lo(2)-1,hi(2)+2
          res(hi(1)+1,j) = res(lo(1)  ,j)
          res(hi(1)+2,j) = res(lo(1)+1,j)
          res(lo(1)-1,j) = res(hi(1)  ,j)
        enddo
      endif

      if (isPeriodic(2) .eq. 1) then
        do i = lo(1)-1,hi(1)+2
          res(i,hi(2)+1) = res(i,lo(2)  )
          res(i,hi(2)+2) = res(i,lo(2)+1)
          res(i,lo(2)-1) = res(i,hi(2)  )
        enddo
      endif

      do j = jstart,jend
        do i = istart,iend

          ii = 2*(i-loc(1))+lo(1)
          jj = 2*(j-loc(2))+lo(2)

          resc(i,j) = fourth*res(ii  ,jj) + 
     $               eighth*(res(ii+1,jj  ) + res(ii-1,jj  ) + 
     $                       res(ii  ,jj+1) + res(ii  ,jj-1) ) +
     $            sixteenth*(res(ii+1,jj+1) + res(ii+1,jj-1) + 
     $                       res(ii-1,jj+1) + res(ii-1,jj-1) )
        enddo
      enddo

c  the top version is what we use when we double the problem at edges
c  the bottom version (commented out) is what we would use if we did not 
c      double the problem at edges.
      if (isPeriodic(1) .NE. 1) then
        i = loc(1)
        ii = 2*(i-loc(1))+lo(1)

        do j = jstart,jend
          jj = 2*(j-loc(2))+lo(2)
          resc(i,j) = fourth*(res(ii,jj  ) + res(ii+1,jj  )) + 
     $                eighth*(res(ii,jj-1) + res(ii+1,jj-1)+
     $                        res(ii,jj+1) + res(ii+1,jj+1) )
cnd          resc(i,j) = fourth*res(ii  ,jj) + 
cnd     $               eighth*(res(ii+1,jj  ) + 
cnd     $                       res(ii  ,jj+1) + res(ii  ,jj-1) ) +
cnd     $            sixteenth*(res(ii+1,jj+1) + res(ii+1,jj-1))
        enddo

        i = hic(1)+1
        ii = 2*(i-loc(1))+lo(1)

        do j = jstart,jend
          jj = 2*(j-loc(2))+lo(2)
          resc(i,j) = fourth*(res(ii,jj  ) + res(ii-1,jj  )) + 
     $                eighth*(res(ii,jj-1) + res(ii-1,jj-1)+
     $                        res(ii,jj+1) + res(ii-1,jj+1) )
cnd          resc(i,j) = fourth*res(ii  ,jj) + 
cnd     $               eighth*(res(ii-1,jj  ) + 
cnd     $                       res(ii  ,jj+1) + res(ii  ,jj-1) ) +
cnd     $            sixteenth*(res(ii-1,jj+1) + res(ii-1,jj-1) )
        enddo


      endif

      if (isPeriodic(2) .NE. 1) then
        j = loc(2)
        jj = 2*(j-loc(2))+lo(2)

        do i = istart,iend
          ii = 2*(i-loc(1))+lo(1)
          resc(i,j) = fourth*(res(ii  ,jj) + res(ii  ,jj+1)) + 
     $                eighth*(res(ii+1,jj) + res(ii+1,jj+1)+
     $                        res(ii-1,jj) + res(ii-1,jj+1) )
cnd          resc(i,j) = fourth*res(ii  ,jj) + 
cnd     $               eighth*(res(ii+1,jj  ) + res(ii-1,jj  ) + 
cnd     $                       res(ii  ,jj+1)) +
cnd     $            sixteenth*(res(ii+1,jj+1) + res(ii-1,jj+1))
        enddo

        j = hic(2)+1
        jj = 2*(j-loc(2))+lo(2)

        do i = istart,iend
          ii = 2*(i-loc(1))+lo(1)
          resc(i,j) = fourth*(res(ii  ,jj) + res(ii  ,jj-1)) + 
     $                eighth*(res(ii+1,jj) + res(ii+1,jj-1)+
     $                        res(ii-1,jj) + res(ii-1,jj-1) )
cnd          resc(i,j) = fourth*res(ii  ,jj) + 
cnd     $               eighth*(res(ii+1,jj  ) + res(ii-1,jj  ) + 
cnd     $                       res(ii  ,jj-1) ) +
cnd     $            sixteenth*(res(ii+1,jj-1) + res(ii-1,jj-1) )
        enddo
      endif

      if ( isPeriodic(1) .NE. 1 .and. isPeriodic(2) .NE. 1 ) then
        i = loc(1)
        j = loc(2)
        ii = 2*(i-loc(1))+lo(1)
        jj = 2*(j-loc(2))+lo(2)
        resc(i,j) = fourth*(res(ii,jj  ) + res(ii+1,jj  ) +
     $                      res(ii,jj+1) + res(ii+1,jj+1) )
cnd          resc(i,j) = fourth*res(ii  ,jj) + 
cnd     $               eighth*(res(ii+1,jj  ) + res(ii  ,jj+1)) +
cnd     $            sixteenth*(res(ii+1,jj+1))

        i = hic(1)+1
        j = hic(2)+1
        ii = 2*(i-loc(1))+lo(1)
        jj = 2*(j-loc(2))+lo(2)
        resc(i,j) = fourth*(res(ii,jj  ) + res(ii-1,jj  ) +
     $                      res(ii,jj-1) + res(ii-1,jj-1) )
cnd          resc(i,j) = fourth*res(ii  ,jj) + 
cnd     $               eighth*(res(ii-1,jj  ) + res(ii  ,jj-1) ) +
cnd     $            sixteenth*(res(ii-1,jj-1) )

        i = loc(1)
        j = hic(2)+1
        ii = 2*(i-loc(1))+lo(1)
        jj = 2*(j-loc(2))+lo(2)

        resc(i,j) = fourth*(res(ii,jj  ) + res(ii+1,jj  ) +
     $                      res(ii,jj-1) + res(ii+1,jj-1) )
cnd          resc(i,j) = fourth*res(ii  ,jj) + 
cnd     $               eighth*(res(ii+1,jj  ) + res(ii  ,jj-1) ) +
cnd     $            sixteenth*(res(ii+1,jj-1) )

        i = hic(1)+1
        j = loc(2)
        ii = 2*(i-loc(1))+lo(1)
        jj = 2*(j-loc(2))+lo(2)

        resc(i,j) = fourth*(res(ii,jj  ) + res(ii-1,jj  ) +
     $                      res(ii,jj+1) + res(ii-1,jj+1) )
cnd          resc(i,j) = fourth*res(ii  ,jj) + 
cnd     $               eighth*(res(ii-1,jj  ) + res(ii  ,jj+1)) +
cnd     $            sixteenth*(res(ii-1,jj+1))

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
      REAL_T     phi(DIM12(phi))
      REAL_T  deltac(DIM12(deltac))
      REAL_T  sigma(DIM12(sigma))
      REAL_T   temp(DIM12(temp))
      integer isPeriodic(SDIM)

c     Local variables
      integer ii,jj,ic,jc
      integer is,ie,js,je,isc,iec,jsc,jec
      logical setSingularPoint

      is = lo(1)
      js = lo(2)
      ie = hi(1)
      je = hi(2)

      isc = loc(1)
      jsc = loc(2)
      iec = hic(1)
      jec = hic(2)

      do jc = jsc, jec+1 
        do ic = isc, iec+1
          ii = 2*(ic-isc)+is
          jj = 2*(jc-jsc)+js
          temp(ii,jj) = deltac(ic,jc)
        enddo
      enddo
      do jc = jsc, jec
        do ic = isc, iec+1
          ii = 2*(ic-isc)+is
          jj = 2*(jc-jsc)+js
          temp(ii,jj+1) = half*(deltac(ic,jc) + deltac(ic,jc+1)) 
        enddo
      enddo
      do jc = jsc, jec+1 
        do ic = isc, iec
          ii = 2*(ic-isc)+is
          jj = 2*(jc-jsc)+js
          temp(ii+1,jj) = half*(deltac(ic,jc) + deltac(ic+1,jc))
        enddo
      enddo
      do jc = jsc, jec 
        do ic = isc, iec 
          ii = 2*(ic-isc)+is
          jj = 2*(jc-jsc)+js
          temp(ii+1,jj+1) = fourth*(deltac(ic,jc  ) + deltac(ic+1,jc  ) + 
     $                             deltac(ic,jc+1) + deltac(ic+1,jc+1) )
        enddo
      enddo

      do jj=js,je+1
         do ii = is,ie+1
            phi(ii,jj) = phi(ii,jj) + temp(ii,jj)
         enddo
      enddo

      return
      end

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
      integer domlo(3)
      integer domhi(3)
      integer lo_bc(3)
      integer hi_bc(3)
      REAL_T  rho(DIMV(rho))
      REAL_T  phi(DIMV(phi))
      REAL_T  dx(3)
      REAL_T  gravity
      
c     Local variables
      integer i,j,k
      REAL_T rhog
      REAL_T rho_i,rho_ip1,rho_im1
      REAL_T rho_j,rho_jp1,rho_jm1
      REAL_T rhoExt
      
#define XLO 0
#define YLO 1
#define ZLO 2
#define XHI 3
#define YHI 4
#define ZHI 5

      if (face .eq. ZLO .or. face .eq. ZHI) 
     $   call bl_abort('SHOULDNT BE IN RHOGBC WITH FACE IN Z-DIR')

c     Ok to only use low index of phi because phi is only one
c        node wide in direction of face.

      if (face .eq. XLO) then

        i = ARG_L1(phi)

        do j = ARG_L2(phi)+1,ARG_H2(phi)-1
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_i   = half * (rho(i  ,j,k) + rho(i  ,j-1,k))
            rho_ip1 = half * (rho(i+1,j,k) + rho(i+1,j-1,k))
            rhoExt  = half * (three*rho_i - rho_ip1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        end do

        j = ARG_L2(phi)
        if (j .eq. domlo(2) .and. 
     $      (lo_bc(2) .eq. HOEXTRAP .or.
     $       lo_bc(2) .eq. FOEXTRAP .or.
     $       lo_bc(2) .eq. EXT_DIR) ) then
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_i   = rho(i  ,j-1,k)
            rho_ip1 = rho(i+1,j-1,k)
            rhoExt  = half * (three*rho_i - rho_ip1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_i   = half * (rho(i  ,j,k) + rho(i  ,j-1,k))
            rho_ip1 = half * (rho(i+1,j,k) + rho(i+1,j-1,k))
            rhoExt  = half * (three*rho_i - rho_ip1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        end if

        j = ARG_H2(phi)
        if (j.eq.domhi(2)+1 .and. 
     $      (hi_bc(2) .eq. HOEXTRAP .or.
     $       hi_bc(2) .eq. FOEXTRAP .or.
     $       hi_bc(2) .eq. EXT_DIR) ) then
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_i   = rho(i  ,j,k)
            rho_ip1 = rho(i+1,j,k)
            rhoExt  = half * (three*rho_i - rho_ip1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_i   = half * (rho(i  ,j,k) + rho(i  ,j-1,k))
            rho_ip1 = half * (rho(i+1,j,k) + rho(i+1,j-1,k))
            rhoExt  = half * (three*rho_i - rho_ip1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        end if

      else if (face .eq. XHI) then

        i = ARG_L1(phi)

        do j = ARG_L2(phi)+1,ARG_H2(phi)-1
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_i   = half * (rho(i-1,j,k) + rho(i-1,j-1,k))
            rho_im1 = half * (rho(i-2,j,k) + rho(i-2,j-1,k))
            rhoExt  = half * (three*rho_i - rho_im1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        end do

        j = ARG_L2(phi)
        if (j .eq. domlo(2) .and. 
     $      (lo_bc(2) .eq. HOEXTRAP .or.
     $       lo_bc(2) .eq. FOEXTRAP .or.
     $       lo_bc(2) .eq. EXT_DIR) ) then
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_i   = rho(i-1,j-1,k)
            rho_im1 = rho(i-2,j-1,k)
            rhoExt  = half * (three*rho_i - rho_im1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_i   = half * (rho(i-1,j,k) + rho(i-1,j-1,k))
            rho_im1 = half * (rho(i-2,j,k) + rho(i-2,j-1,k))
            rhoExt  = half * (three*rho_i - rho_im1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        end if

        j = ARG_H2(phi)
        if (j.eq.domhi(2)+1 .and. 
     $      (hi_bc(2) .eq. HOEXTRAP .or.
     $       hi_bc(2) .eq. FOEXTRAP .or.
     $       hi_bc(2) .eq. EXT_DIR) ) then
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_i   = rho(i-1,j,k)
            rho_im1 = rho(i-2,j,k)
            rhoExt  = half * (three*rho_i - rho_im1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_i   = half * (rho(i-1,j,k) + rho(i-1,j-1,k))
            rho_im1 = half * (rho(i-2,j,k) + rho(i-2,j-1,k))
            rhoExt  = half * (three*rho_i - rho_im1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        end if

      else if (face .eq. YLO) then

        j = ARG_L2(phi)

        do i = ARG_L1(phi)+1,ARG_H1(phi)-1
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_j   = half * (rho(i,j  ,k) + rho(i-1,j  ,k))
            rho_jp1 = half * (rho(i,j+1,k) + rho(i-1,j+1,k))
            rhoExt  = half * (three*rho_j - rho_jp1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        end do

        i = ARG_L1(phi)
        if (i .eq. domlo(1) .and. 
     $      (lo_bc(1) .eq. HOEXTRAP .or.
     $       lo_bc(1) .eq. FOEXTRAP .or.
     $       lo_bc(1) .eq. EXT_DIR) ) then
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_j   = rho(i-1,j  ,k)
            rho_jp1 = rho(i-1,j+1,k)
            rhoExt  = half * (three*rho_j - rho_jp1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_j   = half * (rho(i,j  ,k) + rho(i-1,j  ,k))
            rho_jp1 = half * (rho(i,j+1,k) + rho(i-1,j+1,k))
            rhoExt  = half * (three*rho_j - rho_jp1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        end if

        i = ARG_H1(phi)
        if (i .eq. domhi(1)+1 .and. 
     $      (hi_bc(1) .eq. HOEXTRAP .or.
     $       hi_bc(1) .eq. FOEXTRAP .or.
     $       hi_bc(1) .eq. EXT_DIR) ) then
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_j   = rho(i,j  ,k)
            rho_jp1 = rho(i,j+1,k)
            rhoExt  = half * (three*rho_j - rho_jp1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_j   = half * (rho(i,j  ,k) + rho(i-1,j  ,k))
            rho_jp1 = half * (rho(i,j+1,k) + rho(i-1,j+1,k))
            rhoExt  = half * (three*rho_j - rho_jp1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        end if

      else if (face .eq. YHI) then

        j = ARG_L2(phi)

        do i = ARG_L1(phi)+1,ARG_H1(phi)-1
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_j   = half * (rho(i,j-1,k) + rho(i-1,j-1,k))
            rho_jm1 = half * (rho(i,j-2,k) + rho(i-1,j-2,k))
            rhoExt  = half * (three*rho_j - rho_jm1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        end do

        i = ARG_L1(phi)
        if (i .eq. domlo(1) .and. 
     $      (lo_bc(1) .eq. HOEXTRAP .or.
     $       lo_bc(1) .eq. FOEXTRAP .or.
     $       lo_bc(1) .eq. EXT_DIR) ) then
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_j   = rho(i-1,j-1,k)
            rho_jm1 = rho(i-1,j-2,k)
            rhoExt  = half * (three*rho_j - rho_jm1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_j   = half * (rho(i,j-1,k) + rho(i-1,j-1,k))
            rho_jm1 = half * (rho(i,j-2,k) + rho(i-1,j-2,k))
            rhoExt  = half * (three*rho_j - rho_jm1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        end if

        i = ARG_H1(phi)
        if (i .eq. domhi(1)+1 .and. 
     $      (hi_bc(1) .eq. HOEXTRAP .or.
     $       hi_bc(1) .eq. FOEXTRAP .or.
     $       hi_bc(1) .eq. EXT_DIR) ) then
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_j   = rho(i,j-1,k)
            rho_jm1 = rho(i,j-2,k)
            rhoExt  = half * (three*rho_j - rho_jm1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_j   = half * (rho(i,j-1,k) + rho(i-1,j-1,k))
            rho_jm1 = half * (rho(i,j-2,k) + rho(i-1,j-2,k))
            rhoExt  = half * (three*rho_j - rho_jm1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        end if

      endif

#undef XLO
#undef YLO
#undef ZLO
#undef XHI
#undef YHI
#undef ZHI

      end

c *************************************************************************
c ** FILL_TWOD **
c *************************************************************************

      subroutine FORT_FILL_TWOD(lenx,leny,lenz,length,width,
     $                          faces,numOutFlowFaces,
     $                          cc0,cc1,cc2,cc3,cc4,cc5,conn)

      integer lenx,leny,lenz
      integer length,width
      integer faces(6)
      integer numOutFlowFaces
      REAL_T cc0(leny,lenz,4)
      REAL_T cc1(lenx,lenz,4)
      REAL_T cc2(leny,lenz,4)
      REAL_T cc3(lenx,lenz,4)
      REAL_T cc4(lenx,leny,4)
      REAL_T cc5(lenx,leny,4)
      REAL_T conn(length,width,4)

      integer xlo_outflow,ylo_outflow,zlo_outflow
      integer xhi_outflow,yhi_outflow,zhi_outflow
      integer i,k,ifinal,n

#define XLO 0
#define YLO 1
#define ZLO 2
#define XHI 3
#define YHI 4
#define ZHI 5

c     Want to find the single non-outflow face.
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
c       cc0 cc1 cc2 cc3 cc4 cc5
c       XLO YLO 
c       XLO             YHI 
c           YLO     XHI 
c           YLO     YHI 
c       XLO YLO     XHI
c       XLO         XHI YHI
c       XLO YLO     YHI
c           YLO     XHI YHI
c       XLO YLO     XHI YHI

c     We must remember here that the cc* arrays have already been
c       ordered so that the 3nd dimension is one cell wide.

      ifinal = 0
      if (numOutFlowFaces .eq. 4 .or. 
     $    (xlo_outflow .eq. 1 .and. ylo_outflow .eq. 0) ) then
          do k = 1,lenz
          do i = 1,leny
            conn(i,k,1) = cc0(i,k,1)
            conn(i,k,2) = cc0(i,k,2)
            conn(i,k,3) = cc0(i,k,3)
            conn(i,k,4) = cc0(i,k,4)
          enddo
          enddo
          ifinal = leny
      endif

      if (yhi_outflow .eq. 1 .and. 
     $    .not. (numOutFlowFaces .eq. 3 .and. xhi_outflow .eq. 0) ) then
          do k = 1,lenz
          do i = 1,lenx
            conn(ifinal+i,k,1) = cc4(i,k,1)
            conn(ifinal+i,k,2) = cc4(i,k,2)
            conn(ifinal+i,k,3) = cc4(i,k,3)
            conn(ifinal+i,k,4) = cc4(i,k,4)
          enddo
          enddo
          ifinal = ifinal + lenx
      endif

      if (xhi_outflow .eq. 1) then
          do k = 1,lenz
          do i = 1,leny
            conn(ifinal+i,k,1) = cc3(leny+1-i,k,1)
            conn(ifinal+i,k,2) = cc3(leny+1-i,k,2)
            conn(ifinal+i,k,3) = cc3(leny+1-i,k,3)
            conn(ifinal+i,k,4) = cc3(leny+1-i,k,4)
          enddo
          enddo
          ifinal = ifinal + leny
      endif

      if (ylo_outflow .eq. 1) then
          do k = 1,lenz
          do i = 1,lenx
            conn(ifinal+i,k,1) = cc1(lenx+1-i,k,1)
            conn(ifinal+i,k,2) = cc1(lenx+1-i,k,2)
            conn(ifinal+i,k,3) = cc1(lenx+1-i,k,3)
            conn(ifinal+i,k,4) = cc1(lenx+1-i,k,4)
          enddo
          enddo
          ifinal = ifinal + lenx
      endif

      if (numOutFlowFaces .lt. 4 .and.
     $    (xlo_outflow .eq. 1 .and. ylo_outflow .eq. 1) ) then
          do k = 1,lenz
          do i = 1,leny
            conn(ifinal+i,k,1) = cc0(i,k,1)
            conn(ifinal+i,k,2) = cc0(i,k,2)
            conn(ifinal+i,k,3) = cc0(i,k,3)
            conn(ifinal+i,k,4) = cc0(i,k,4)
          enddo
          enddo
          ifinal = ifinal + leny
      endif

      if (yhi_outflow .eq. 1 .and. 
     $    (numOutFlowFaces .eq. 3 .and. xhi_outflow .eq. 0) ) then
          do k = 1,lenz
          do i = 1,lenx
            conn(ifinal+i,k,1) = cc4(i,k,1)
            conn(ifinal+i,k,2) = cc4(i,k,2)
            conn(ifinal+i,k,3) = cc4(i,k,3)
            conn(ifinal+i,k,4) = cc4(i,k,4)
          enddo
          enddo
          ifinal = ifinal + lenx
      endif

      length = ifinal

#undef XLO
#undef YLO
#undef ZLO
#undef XHI
#undef YHI
#undef ZHI

      end

c *************************************************************************
c ** ALLPHI_FROM_X **
c *************************************************************************

      subroutine FORT_ALLPHI_FROM_X(lenx,leny,lenz,length,width,
     $                              faces,numOutFlowFaces,
     $                              phi0,phi1,phi2,phi3,phi4,phi5,x,DIMS(x))

      integer DIMDEC(x)
      integer lenx,leny,lenz,length,width
      integer numOutFlowFaces
      integer faces(6)
      REAL_T phi0(0:leny,0:lenz)
      REAL_T phi1(0:lenx,0:lenz)
      REAL_T phi2(0:leny,0:lenz)
      REAL_T phi3(0:lenx,0:lenz)
      REAL_T phi4(0:lenx,0:leny)
      REAL_T phi5(0:lenx,0:leny)
      REAL_T x(DIMV(x))

      integer xlo_outflow,ylo_outflow,zlo_outflow
      integer xhi_outflow,yhi_outflow,zhi_outflow
      integer i,j,k,ifinal,n

#define XLO 0
#define YLO 1
#define ZLO 2
#define XHI 3
#define YHI 4
#define ZHI 5

c     Want to find the single non-outflow face.
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
c       phi0 phi1 phi2 phi3 phi4 phi5
c       XLO  YLO 
c       XLO                 YHI 
c            YLO       XHI 
c            YLO            YHI 
c       XLO  YLO       XHI
c       XLO            XHI  YHI
c       XLO  YLO            YHI
c            YLO       XHI  YHI
c       XLO  YLO       XHI  YHI


c     We know that the faces are ordered: XLO,XHI,YLO,YHI
      
      ifinal = 0

      if (numOutFlowFaces .eq. 4 .or. 
     $    (xlo_outflow .eq. 1 .and. ylo_outflow .eq. 0) ) then
        do k = 0,lenz
        do j = 0,leny
          phi0(j,k) = x(j,k,ARG_L3(x))
        enddo
        enddo
        ifinal = leny
      endif

      if (yhi_outflow .eq. 1 .and. 
     $    .not. (numOutFlowFaces .eq. 3 .and. xhi_outflow .eq. 0) ) then
        do k = 0,lenz
        do i = 0,lenx
          phi4(i,k) = x(i+ifinal,k,ARG_L3(x))
        enddo
        enddo
        ifinal = ifinal+lenx
      endif

      if (xhi_outflow .eq. 1) then
        do k = 0,lenz
        do j = 0,leny
          phi3(leny-j,k) = x(ifinal+j,k,ARG_L3(x))
        enddo
        enddo
        ifinal = ifinal+leny
      endif

      if (ylo_outflow .eq. 1) then
        if (numOutFlowFaces .eq. 4) then 
          do k = 0,lenz
          do i = 0,lenx-1
            phi1(lenx-i,k) = x(ifinal+i,k,ARG_L3(x))
          enddo
          phi1(0,k) = x(0,k,ARG_L3(x))
          enddo
        else
          do k = 0,lenz
          do i = 0,lenx
            phi1(lenx-i,k) = x(ifinal+i,k,ARG_L3(x))
          enddo
          enddo
        endif
        ifinal = ifinal+lenx
      endif

      if (numOutFlowFaces .lt. 4 .and.
     $    (xlo_outflow .eq. 1 .and. ylo_outflow .eq. 1) ) then
        do k = 0,lenz
        do j = 0,leny
          phi0(j,k) = x(j+ifinal,k,ARG_L3(x))
        enddo
        enddo
        ifinal = ifinal+leny
      endif

      if (yhi_outflow .eq. 1 .and. 
     $    (numOutFlowFaces .eq. 3 .and. xhi_outflow .eq. 0) ) then
        do k = 0,lenz
        do i = 0,lenx
          phi4(i,k) = x(i+ifinal,k,ARG_L3(x))
        enddo
        enddo
        ifinal = ifinal+lenx
      endif

#undef XLO
#undef YLO
#undef ZLO
#undef XHI
#undef YHI
#undef ZHI
      end
