#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <MACOUTFLOWBC_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 3

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define SMALL 1.0e-10
#else
#define SMALL 1.0d-10
#endif

module macoutflowbc_3d_module
  
  implicit none

  private 

  public :: extrap_mac
            
contains


      subroutine extrap_mac(DIMS(u0),u0,DIMS(u1),u1,DIMS(u2),u2,DIMS(div), &
                            divu,DIMS(rho),rho, &
                            DIMS(divuExt),divuExt,DIMS(rhoExt),rhoExt, &
                            dx,lo,hi,face,per,zeroIt,small_udiff) &
                            bind(C,name="extrap_mac")
!c
!c     Compute the value of phi for macproj 
!c     assuming that the tangential velocity on the edges of the outflow boundary
!c     are either zero or periodic.
!c     note that u is edge centered

!c    compute divu_ave twice due to precision problems

      use projoutflowbc_3d_module, only : subtractavg

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
      REAL_T small_udiff

!c     Local variables
      REAL_T small_pert
      parameter ( small_pert = SMALL)
      REAL_T max_divu, min_divu, max_pert
      REAL_T divu_ave1,divu_ave2
      REAL_T hx,hy,hz  
      REAL_T diff
      integer i,j,k
      integer ics,ice,jcs,jce,kcs,kce
      integer ifs,ife,jfs,jfe,kfs,kfe
      integer if,jf,kf
!c     NOTE: Assumes that rho at edge between i, i-1 = half*(rho(i)+rho(i-1))
!c             (1) Linear fit of rho between nodes
!c             (2) rho, divu on same boxes (box)
!c             (3) phi is on box, shifted up one
!c             (4) u is edge-based, on surroundingNodes(box)

!c     Compute average of divu over outflow bc.  Set trivial solution if average
!c     is zero, or if divu is constant
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

!c        Here we modify divuExt to include the velocity terms.
         do k = kcs, kce
         do j = jcs, jce
            divuExt(j,k,if) = divuExt(j,k,if)  &
             - (u1(ics,j+1,k)-u1(ics,j,k))/hy  &
             - (u2(ics,j,k+1)-u2(ics,j,k))/hz 
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

!c        Make sure u_mac is periodic
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

!c        Here we modify divuExt to include the velocity terms.
         do k = kcs, kce
         do i = ics, ice
            divuExt(i,k,jf) = divuExt(i,k,jf) &
             - (u0(i+1,jcs,k)-u0(i,jcs,k))/hx &
             - (u2(i,jcs,k+1)-u2(i,jcs,k))/hz
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

!c        Make sure u_mac is periodic
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

!c        Here we modify divuExt to include the velocity terms.
         do j = jcs, jce
         do i = ics, ice
            divuExt(i,j,kf) = divuExt(i,j,kf) &
             - (u0(i+1,j,kcs)-u0(i,j,kcs))/hx &
             - (u1(i,j+1,kcs)-u1(i,j,kcs))/hy
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

!c        Make sure u_mac is periodic
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
            divuExt(j,k,if) = divuExt(j,k,if) &
             - (u1(ice,j+1,k)-u1(ice,j,k))/hy  &
             - (u2(ice,j,k+1)-u2(ice,j,k))/hz 
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

!c        Make sure u_mac is periodic
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
            divuExt(i,k,jf) = divuExt(i,k,jf) &
             - (u0(i+1,jce,k)-u0(i,jce,k))/hx &
             - (u2(i,jce,k+1)-u2(i,jce,k))/hz
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

!c        Make sure u_mac is periodic
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

!c        Here we modify divuExt to include the velocity terms.
         do j = jcs, jce
         do i = ics, ice
            divuExt(i,j,kf) = divuExt(i,j,kf) &
             - (u0(i+1,j,kce)-u0(i,j,kce))/hx &
             - (u1(i,j+1,kce)-u1(i,j,kce))/hy
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

!c        Make sure u_mac is periodic
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

!c  check to see if we should zero phi
         max_pert = max_pert/(ABS(divu_ave1+divu_ave2)+small_pert)
      if ((max_divu.eq.zero.and.min_divu.eq.zero) &
          .or.(max_pert.le.small_pert)) then
         zeroIt = 1
      end if

      end subroutine extrap_mac

  end module macoutflowbc_3d_module
