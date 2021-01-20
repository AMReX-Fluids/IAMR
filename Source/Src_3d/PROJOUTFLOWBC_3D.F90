#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PROJOUTFLOWBC_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 3

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define SMALL 1.0e-10
#define sixteenth  .0625e0
#else
#define SMALL 1.0d-10
#define sixteenth  .0625d0
#endif

module projoutflowbc_3d_module
  
  implicit none

  private 

  public  :: extrap_proj, rhogbc, subtractavg

contains

      subroutine extrap_proj(DIMS(u),u,DIMS(divu),divu,DIMS(rho),rho, &
                             DIMS(uExt),uExt,DIMS(divuExt),divuExt, &
                             DIMS(rhoExt),rhoExt,lo,hi,face, zeroIt) &
                             bind(C,name="extrap_proj")
                             
      implicit none

!c    compute divu_ave twice due to precision problems

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

!c local variables
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

!c  check to see if we should zero phi
      max_pert = max_pert/(ABS(divu_ave1+divu_ave2)+small_pert)
      if ((max_divu.eq.zero) .or. (max_pert.le.small_pert)) zeroIt = 1
#undef XLO
#undef YLO
#undef ZLO
#undef XHI
#undef YHI
#undef ZHI
      end subroutine extrap_proj

      subroutine rhogbc(rho,DIMS(rho),phi,DIMS(phi), &
                        face,gravity,dx,domlo,domhi,lo_bc,hi_bc) &
                        bind(C, name="rhogbc")
!c
!c    Compute the contribution of gravity to the boundary conditions
!c      for phi at outflow faces only.
!c
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
      
!c     Local variables
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

      if (face .eq. ZLO .or. face .eq. ZHI) &
        call bl_abort('SHOULDNT BE IN RHOGBC WITH FACE IN Z-DIR')

!c     Ok to only use low index of phi because phi is only one
!c        node wide in direction of face.

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

        if (j .eq. domlo(2) .and. lo_bc(2) .eq. EXT_DIR) then
          rhog = zero
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_i   = rho(i  ,j-1,k)
            rho_ip1 = rho(i+1,j-1,k)
            rhoExt  = half * (three*rho_i - rho_ip1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
         else if (j .eq. domlo(2) .and. lo_bc(2) .eq. HOEXTRAP ) then 
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_i   = half*(three*rho(i  ,j,k) - rho(i  ,j+1,k))
            rho_ip1 = half*(three*rho(i+1,j,k) - rho(i+1,j+1,k))
            rhoExt  = half * (three*rho_i - rho_ip1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else if (j .eq. domlo(2) .and. lo_bc(2) .eq. FOEXTRAP ) then 
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

        j = ARG_H2(phi)
        rhog = zero
        if (j.eq.domhi(2)+1 .and. hi_bc(2) .eq. EXT_DIR) then
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_i   = rho(i  ,j,k)
            rho_ip1 = rho(i+1,j,k)
            rhoExt  = half * (three*rho_i - rho_ip1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else if (j.eq.domhi(2)+1 .and. hi_bc(2) .eq. HOEXTRAP) then
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_i   = half*(three*rho(i  ,j-1,k) - rho(i  ,j-2,k))
            rho_ip1 = half*(three*rho(i+1,j-1,k) - rho(i+1,j-2,k))
            rhoExt  = half * (three*rho_i - rho_ip1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else if (j.eq.domhi(2)+1 .and. hi_bc(2) .eq. FOEXTRAP) then
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
        rhog = zero
        if (j .eq. domlo(2) .and. lo_bc(2) .eq. EXT_DIR) then
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_i   = rho(i-1,j-1,k)
            rho_im1 = rho(i-2,j-1,k)
            rhoExt  = half * (three*rho_i - rho_im1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else if (j .eq. domlo(2) .and. lo_bc(2) .eq. HOEXTRAP) then
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_i   = half*(three*rho(i-1,j,k) - rho(i-1,j+1,k))
            rho_im1 = half*(three*rho(i-2,j,k) - rho(i-2,j+1,k))
            rhoExt  = half * (three*rho_i - rho_im1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else if (j .eq. domlo(2) .and. lo_bc(2) .eq. FOEXTRAP) then
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

        j = ARG_H2(phi)
        rhog = zero
        if (j.eq.domhi(2)+1 .and. hi_bc(2) .eq. EXT_DIR) then
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_i   = rho(i-1,j,k)
            rho_im1 = rho(i-2,j,k)
            rhoExt  = half * (three*rho_i - rho_im1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else if (j.eq.domhi(2)+1 .and. hi_bc(2) .eq. HOEXTRAP) then
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_i   = half*(three*rho(i-1,j-1,k) - rho(i-1,j-2,k))
            rho_im1 = half*(three*rho(i-2,j-1,k) - rho(i-2,j-2,k))
            rhoExt  = half * (three*rho_i - rho_im1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else if (j.eq.domhi(2)+1 .and. hi_bc(2) .eq. FOEXTRAP) then
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
        rhog = zero
        if (i .eq. domlo(1) .and. lo_bc(1) .eq. EXT_DIR) then
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_j   = rho(i-1,j  ,k)
            rho_jp1 = rho(i-1,j+1,k)
            rhoExt  = half * (three*rho_j - rho_jp1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else if (i .eq. domlo(1) .and. lo_bc(1) .eq. HOEXTRAP) then
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_j   = half*(three*rho(i,j  ,k) - rho(i+1,j  ,k))
            rho_jp1 = half*(three*rho(i,j+1,k) - rho(i+1,j+1,k))
            rhoExt  = half * (three*rho_j - rho_jp1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else if (i .eq. domlo(1) .and. lo_bc(1) .eq. FOEXTRAP) then
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

        i = ARG_H1(phi)
        rhog = zero
        if (i .eq. domhi(1)+1 .and. hi_bc(1) .eq. EXT_DIR) then
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_j   = rho(i,j  ,k)
            rho_jp1 = rho(i,j+1,k)
            rhoExt  = half * (three*rho_j - rho_jp1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else if (i .eq. domhi(1)+1 .and. hi_bc(1) .eq. HOEXTRAP) then
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_j   = half*(three*rho(i-1,j  ,k) - rho(i-2,j  ,k))
            rho_jp1 = half*(three*rho(i-1,j+1,k) - rho(i-2,j+1,k))
            rhoExt  = half * (three*rho_j - rho_jp1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else if (i .eq. domhi(1)+1 .and. hi_bc(1) .eq. FOEXTRAP) then
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
        rhog = zero
        if (i .eq. domlo(1) .and. lo_bc(1) .eq. EXT_DIR) then
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_j   = rho(i-1,j-1,k)
            rho_jm1 = rho(i-1,j-2,k)
            rhoExt  = half * (three*rho_j - rho_jm1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else if (i .eq. domlo(1) .and. lo_bc(1) .eq. HOEXTRAP) then
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_j   = half*(three*rho(i,j-1,k) - rho(i+1,j-1,k))
            rho_jm1 = half*(three*rho(i,j-2,k) - rho(i+1,j-2,k))
            rhoExt  = half * (three*rho_j - rho_jm1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else if (i .eq. domlo(1) .and. lo_bc(1) .eq. FOEXTRAP) then
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

        i = ARG_H1(phi)
        rhog = zero
        if (i .eq. domhi(1)+1 .and. hi_bc(1) .eq. EXT_DIR) then
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_j   = rho(i,j-1,k)
            rho_jm1 = rho(i,j-2,k)
            rhoExt  = half * (three*rho_j - rho_jm1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else if (i .eq. domhi(1)+1 .and. hi_bc(1) .eq. HOEXTRAP) then
          do k = ARG_H3(phi)-1,ARG_L3(phi),-1
            rho_j   = half*(three*rho(i-1,j-1,k) - rho(i-2,j-1,k))
            rho_jm1 = half*(three*rho(i-1,j-2,k) - rho(i-2,j-2,k))
            rhoExt  = half * (three*rho_j - rho_jm1 )
            rhog = rhog - gravity * rhoExt * dx(3)
            phi(i,j,k) = phi(i,j,k) + rhog
          end do
        else if (i .eq. domhi(1)+1 .and. hi_bc(1) .eq. FOEXTRAP) then
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

      endif

#undef XLO
#undef YLO
#undef ZLO
#undef XHI
#undef YHI
#undef ZHI

      end subroutine rhogbc

           subroutine subtractavg(DIMS(divu),divu,lo,hi,divu_ave,face)
      implicit none
      integer DIMDEC(divu)
      integer lo(SDIM),hi(SDIM)
      REAL_T divu_ave
      REAL_T divu(DIMV(divu))
      integer face

      integer i,j,k
      REAL_T vtot

#define XLO 0
#define YLO 1
#define ZLO 2
#define XHI 3
#define YHI 4
#define ZHI 5
      
      
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

#undef XLO
#undef YLO
#undef ZLO
#undef XHI
#undef YHI
#undef ZHI
      
      end subroutine subtractavg

end module projoutflowbc_3d_module
