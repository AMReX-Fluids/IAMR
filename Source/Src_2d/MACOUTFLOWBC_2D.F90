#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <MACOUTFLOWBC_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 2

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define SMALL 1.0e-10
#else
#define SMALL 1.0d-10
#endif


module macoutflowbc_2d_module
  
  implicit none

  private 

  public :: extrap_mac
  
contains

!c *************************************************************************
!c ** EXTRAP_MAC
!c *************************************************************************

      subroutine extrap_mac(DIMS(u0),u0,DIMS(u1),u1,DIMS(div),divu,DIMS(rho),rho,&
                              r_len,redge,DIMS(divuExt),divuExt,&
                              DIMS(rhoExt),rhoExt,dx,lo,hi,face,per,zeroIt,&
                              small_udiff) bind(C,name="extrap_mac")
!c
!c     Compute the value of phi for macproj 
!c
!c     (subtract divu_ave twice due to precision problems)
        use projoutflowbc_2d_module, only : subtractavg
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
      REAL_T small_udiff
      
!c     Local variables
      REAL_T small_pert
      parameter ( small_pert = SMALL)
      integer i, j
      REAL_T divu_ave1,divu_ave2
      REAL_T max_divu, min_divu, max_pert
      REAL_T diff
      REAL_T rc,hx,hy
      integer ics,ice,jcs,jce
      integer ifs,ife,jfs,jfe
      integer if,jf
!c     NOTE: Assumes that rho at edge between i, i-1 = half*(rho(i)+rho(i-1))
!c             (1) Linear fit of rho between nodes
!c             (2) rho, divu on same boxes (box)
!c             (3) phi is on box, shifted up one
!c             (4) u is edge-based, on surroundingNodes(box)

!c     Compute average of divu over outflow bc.  Set trivial solution if average
!c     is zero, or if divu is constant
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

!c        Here we modify divuExt to include the velocity terms.
         do j = jcs, jce
            divuExt(j,if) = redge(j-jcs)*(divuExt(j,if)*hy*hy - (u1(ics,j+1)-u1(ics,j))*hy)
         end do

         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave2,face)

         max_pert = ABS(divuExt(jcs,if))
         do j = jcs, jce
            max_pert = MAX(max_pert,ABS(divuExt(j,if)))
         end do
      
!c        Make sure u_ma!c is periodic
         if (per .eq. 1) then
           diff = u1(ics,jcs)-u1(ics,jce+1)
           if (ABS(diff) .gt. small_udiff) then
              write(6,*) 'EXTRAPMAC: FACE XLO : uma!c not periodic'
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

!c        Here we modify divuExt to include the velocity terms.
         do i = ics, ice
            rc = half*(redge(i+1-ics)+redge(i-ics))
            divuExt(i,jf) = rc*divuExt(i,jf)*hx*hx - &
                           (redge(i+1-ics)*u0(i+1,jcs)-redge(i-ics)*u0(i,jcs))*hx
         end do

         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(ics,jf))
         do i = ics, ice
            max_pert = MAX(max_pert,ABS(divuExt(i,jf)))
         end do
      
!c        Make sure u_ma!c is periodic
         if (per .eq. 1) then
           diff = u0(ics,jcs)-u0(ice+1,jcs)
           if (ABS(diff) .gt. small_udiff) then
              write(6,*) 'EXTRAPMAC: FACE YLO : uma!c not periodic'
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

!c        Here we modify divuExt to include the velocity terms.
         do j = jcs, jce
            divuExt(j,if) = redge(j-jcs)*(divuExt(j,if)*hy*hy - (u1(ice,j+1)-u1(ice,j))*hy)
         end do

         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave2,face)

         max_pert = ABS(divuExt(jcs,if))
         do j = jcs, jce
            max_pert = MAX(max_pert,ABS(divuExt(j,if)))
         end do
      
!c        Make sure u_mac is periodic
         if (per .eq. 1) then
           diff = u1(ice,jcs)-u1(ice,jce+1)
           if (ABS(diff) .gt. small_udiff) then
              write(6,*) 'EXTRAPMAC: FACE XHI : uma!c not periodic'
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

!c        Here we modify divuExt to include the velocity terms.
         do i = ics, ice
            rc = half*(redge(i+1-ics)+redge(i-ics))
            divuExt(i,jf) = rc*divuExt(i,jf)*hx*hx - &
                           (redge(i+1-ics)*u0(i+1,jce)-redge(i-ics)*u0(i,jce))*hx
         end do

         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave2,face)

         max_pert = ABS(divuExt(ics,jf))
         do i = ics, ice
            max_pert = MAX(max_pert,ABS(divuExt(i,jf)))
         end do
      
!c        Make sure u_ma!c is periodic
         if (per .eq. 1) then
           diff = u0(ics,jce)-u0(ice+1,jce)
           if (ABS(diff) .gt. small_udiff) then
              write(6,*) 'EXTRAPMAC: FACE YHI : uma!c not periodic'
              write(6,*) 'U AT LEFT: ',u0(ics  ,jce)
              write(6,*) 'U AT RGHT: ',u0(ice+1,jce)
              call bl_abort(" ")
           endif
         endif

      endif
      
!c  check to see if we should zero phi
         max_pert = max_pert/(ABS(divu_ave1+divu_ave2)+small_pert)
      if ((max_divu.eq.zero.and.min_divu.eq.zero)&
          .or.(max_pert.le.small_pert)) then
         zeroIt = 1
      end if
    end subroutine extrap_mac

  end module macoutflowbc_2d_module
