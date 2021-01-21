
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PROJECTION_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 3

module projection_3d_module

  implicit none

  private

  public ::  anelcoeffmpy, rhogbc

contains

      subroutine anelcoeffmpy (lo,hi,a,DIMS(a),domlo,domhi, &
                               anel_coeff,anel_lo,anel_hi,bogus_value,mult) &
                               bind(C, name="anelcoeffmpy")
!c 
!c     multiply A by anel_coeff
!c
!c 
!c     NOTE: THIS ROUTINE HAS BEEN MODIFIED SO THAT ALL VALUES
!c           OUTSIDE THE DOMAIN ARE SET TO BOGUS VALUE
      !c-- only sets boudnary cells at top and bottom to Bogus Val,
      !c   boundary cells on the sides are left alone

      implicit none
      integer    lo(SDIM),hi(SDIM)
      integer    DIMDEC(a)
      integer    domlo(3), domhi(3)
      REAL_T     a(DIMV(a))
      integer    mult, anel_lo, anel_hi
      REAL_T     anel_coeff(anel_lo:anel_hi)
      REAL_T     bogus_value

      integer i, j, k
      integer klo,khi

      klo = lo(3)
      khi = hi(3)
      
      if (lo(3) .lt. domlo(3)) then
         klo = domlo(3)
         do k = lo(3), domlo(3)-1
         do j = lo(2), hi(2)
         do i = lo(1),hi(1)
           a(i,j,k) = bogus_value
         end do
         end do
         end do
      end if
      if (hi(3) .gt. domhi(3)) then
         khi = domhi(3)
         do k = domhi(3)+1, hi(3)
         do j = lo(2), hi(2)
         do i = lo(1),hi(1)
           a(i,j,k) = bogus_value
         end do
         end do
         end do
      end if

      if (mult .eq. 1) then
         do k = klo,khi
         do j = lo(2),hi(2)
         do i = lo(1),hi(1)
           a(i,j,k) = a(i,j,k) * anel_coeff(k)
         end do
         end do
         end do
      else if (mult .eq. 0) then
         do k = klo,khi
         do j = lo(2),hi(2)
         do i = lo(1),hi(1)
           a(i,j,k) = a(i,j,k) / anel_coeff(k)
         end do
         end do
         end do
      else 
         print *,'BOGUS MULT IN ANELCOEFFMULT ',mult
         stop
      end if

      end subroutine anelcoeffmpy

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

end module projection_3d_module
