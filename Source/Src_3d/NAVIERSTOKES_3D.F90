
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <NAVIERSTOKES_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 3

module navierstokes_3d_module
  
  implicit none

  private 

  public :: sumturb, &
            fort_maxval, cen2edg, FORT_AVERAGE_EDGE_STATES
  
contains

!c ::
!c :: ----------------------------------------------------------
!c :: SUMTURB
!c :: ----------------------------------------------------------
!c ::
      subroutine sumturb(dat,pres,DIMS(dat),DIMS(pres),DIMS(grid),delta, &
                             turb,ksize,turbVars) bind(C, name="sumturb")

      implicit none

      integer DIMDEC(dat)
      integer DIMDEC(pres)
      integer DIMDEC(grid)
      integer ksize, turbVars
      REAL_T  delta(SDIM)
      REAL_T  dat(DIMV(dat),16)
      REAL_T  pres(DIMV(pres),4)
      REAL_T  turb(0:ksize*turbVars-1)
      
      integer i, j, k
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      
      REAL_T area
      REAL_T rho, ux, uy, uz, p
      REAL_T drhodx, drhody, drhodz
      REAL_T duxdx, duxdy, duxdz
      REAL_T duydx, duydy, duydz
      REAL_T duzdx, duzdy, duzdz
      REAL_T dpdx, dpdy, dpdz
      REAL_T dx, dy, dz
      
      ilo = ARG_L1(grid)
      ihi = ARG_H1(grid)
      jlo = ARG_L2(grid)
      jhi = ARG_H2(grid)
      klo = ARG_L3(grid)
      khi = ARG_H3(grid)
      
      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      area = dx*dy
      
      do k = klo, khi
         do i = ilo, ihi
            do j = jlo, jhi

               rho  = dat(i,j,k,1)
               ux   = dat(i,j,k,2)
               uy   = dat(i,j,k,3)
               uz   = dat(i,j,k,4)

!c     Here are the derivatives, can't do it here because of zeroed intersections
               drhodx= dat(i,j,k,5)
               duxdx = dat(i,j,k,6)
               duydx = dat(i,j,k,7)
               duzdx = dat(i,j,k,8)

               drhody= dat(i,j,k,9)
               duxdy = dat(i,j,k,10)
               duydy = dat(i,j,k,11)
               duzdy = dat(i,j,k,12)

               drhodz= dat(i,j,k,13)
               duxdz = dat(i,j,k,14)
               duydz = dat(i,j,k,15)
               duzdz = dat(i,j,k,16)

               p    = pres(i,j,k,1)
               dpdx = pres(i,j,k,2)
               dpdy = pres(i,j,k,3)
               dpdz = pres(i,j,k,4)

               turb(k*turbVars+0)  = turb(k*turbVars+0)  + area*rho
               turb(k*turbVars+1)  = turb(k*turbVars+1)  + area*rho*ux
               turb(k*turbVars+2)  = turb(k*turbVars+2)  + area*rho*uy
               turb(k*turbVars+3)  = turb(k*turbVars+3)  + area*rho*uz
               turb(k*turbVars+4)  = turb(k*turbVars+4)  + area*rho*ux*ux
               turb(k*turbVars+5)  = turb(k*turbVars+5)  + area*rho*ux*uy
               turb(k*turbVars+6)  = turb(k*turbVars+6)  + area*rho*ux*uz
               turb(k*turbVars+7)  = turb(k*turbVars+7)  + area*rho*uy*uy
               turb(k*turbVars+8)  = turb(k*turbVars+8)  + area*rho*uy*uz
               turb(k*turbVars+9)  = turb(k*turbVars+9)  + area*rho*uz*uz
               turb(k*turbVars+10) = turb(k*turbVars+10) + area*rho*ux*ux*ux
               turb(k*turbVars+11) = turb(k*turbVars+11) + area*rho*ux*ux*uy
               turb(k*turbVars+12) = turb(k*turbVars+12) + area*rho*ux*ux*uz
               turb(k*turbVars+13) = turb(k*turbVars+13) + area*rho*ux*uy*uy
               turb(k*turbVars+14) = turb(k*turbVars+14) + area*rho*ux*uy*uz
               turb(k*turbVars+15) = turb(k*turbVars+15) + area*rho*ux*uz*uz
               turb(k*turbVars+16) = turb(k*turbVars+16) + area*rho*uy*uy*uy
               turb(k*turbVars+17) = turb(k*turbVars+17) + area*rho*uy*uy*uz
               turb(k*turbVars+18) = turb(k*turbVars+18) + area*rho*uz*uz*uz

            end do
         end do
      end do
      
      end subroutine sumturb


       subroutine fort_maxval(rho,DIMS(rho),DIMS(grid),mxval)bind(C, name="fort_maxval")

       implicit none

       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  rho(DIMV(rho))
       REAL_T  mxval

       integer i, j, k

       mxval = -Huge(0.0d0)

       do k = ARG_L3(grid), ARG_H3(grid)
          do j = ARG_L2(grid), ARG_H2(grid)
             do i = ARG_L1(grid), ARG_H1(grid)
                mxval = max(mxval, rho(i,j,k))
             end do
          end do
       end do

       end subroutine fort_maxval

!c-----------------------------------------------------------------------
!c     This routine fills an edge-centered fab from a cell-centered
!c     fab using simple linear interpolation.
!c
!c     INPUTS / OUTPUTS:
!c     lo,hi      => index limits of the region of the edge-centered fab
!c                   to be filled
!c     DIMS(cfab) => index limits of the cell-centered fab
!c     cfab       => cell-centered data
!c     DIMS(efab) => index limits of the edge-centered fab
!c     efab       => edge-centered fab to fill
!c     nc         => Number of components in the fab to fill
!c     dir        => direction data needs to be shifted to get to edges
!c-----------------------------------------------------------------------
!c
      subroutine cen2edg(lo, hi, &
          DIMS(cfab), cfab, &
          DIMS(efab), efab, nc, dir, isharm &
          ) bind(C,name="cen2edg")
          
      implicit none

      integer lo(SDIM), hi(SDIM), nc, dir, isharm
      integer DIMDEC(cfab)
      integer DIMDEC(efab)
      REAL_T  cfab(DIMV(cfab), nc)
      REAL_T  efab(DIMV(efab), nc)
      integer i,j,k,n

      if ( isharm .eq. 0 ) then
         if (dir .EQ. 0) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        efab(i,j,k,n) = half*(cfab(i,j,k,n) + cfab(i-1,j,k,n))
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 1) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        efab(i,j,k,n) = half*(cfab(i,j,k,n) + cfab(i,j-1,k,n))
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 2) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        efab(i,j,k,n) = half*(cfab(i,j,k,n) + cfab(i,j,k-1,n))
                     end do
                  end do
               end do
            end do
         end if
      else
         if (dir .EQ. 0) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        if((cfab(i,j,k,n) * cfab(i-1,j,k,n)) .gt.zero) &
                            then
                           efab(i,j,k,n) = &
                               2*(cfab(i,j,k,n) * cfab(i-1,j,k,n))/ &
                               (cfab(i,j,k,n) + cfab(i-1,j,k,n))
                        else
                           efab(i,j,k,n) = zero
                        endif
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 1) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        if((cfab(i,j,k,n) * cfab(i,j-1,k,n)).gt.zero) &
                            then
                           efab(i,j,k,n) = &
                               2*(cfab(i,j,k,n) * cfab(i,j-1,k,n))/ &
                               (cfab(i,j,k,n) + cfab(i,j-1,k,n))
                        else
                           efab(i,j,k,n) = zero
                        endif
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 2) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        if((cfab(i,j,k,n) * cfab(i,j,k-1,n)).gt.zero) &
                            then
                           efab(i,j,k,n) = &
                               2*(cfab(i,j,k,n) * cfab(i,j,k-1,n))/ &
                               (cfab(i,j,k,n) + cfab(i,j,k-1,n))
                        else
                           efab(i,j,k,n) = zero
                        endif
                     end do
                  end do
               end do
            end do
         end if
      end if

      end subroutine cen2edg

!c     
!c     
!c     ::: -----------------------------------------------------------
!c     
!c     This routine averages the mac face velocities for makeforce at half time

   subroutine FORT_AVERAGE_EDGE_STATES( vel, v_lo, v_hi,&
                                        umacx, ux_lo, ux_hi,&
                                        umacy, uy_lo, uy_hi,&
#if ( AMREX_SPACEDIM == 3 )
                                        umacz, uz_lo, uz_hi,&
#endif
                                        getForceVerbose)&
                                        bind(C, name="FORT_AVERAGE_EDGE_STATES")

      implicit none

      integer :: v_lo(3), v_hi(3)
      integer :: ux_lo(3), ux_hi(3)
      integer :: uy_lo(3), uy_hi(3)
#if ( AMREX_SPACEDIM == 3 )
      integer :: uz_lo(3), uz_hi(3)
#endif
      integer :: getForceVerbose
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3), SDIM) :: vel
      REAL_T, dimension(ux_lo(1):ux_hi(1),ux_lo(2):ux_hi(2),ux_lo(3):ux_hi(3)) :: umacx
      REAL_T, dimension(uy_lo(1):uy_hi(1),uy_lo(2):uy_hi(2),uy_lo(3):uy_hi(3)) :: umacy
#if ( AMREX_SPACEDIM == 3 )
      REAL_T, dimension(uz_lo(1):uz_hi(1),uz_lo(2):uz_hi(2),uz_lo(3):uz_hi(3)) :: umacz
#endif

      REAL_T  :: velmin(3)
      REAL_T  :: velmax(3)
      integer :: isioproc

      integer :: i, j, k, n

      do n = 1, SDIM
         velmin(n) = 1.d234
         velmax(n) = -1.d234
      enddo

      do k = v_lo(3), v_hi(3)
         do j = v_lo(2), v_hi(2)
            do i = v_lo(1), v_hi(1)
               vel(i,j,k,1) = half*(umacx(i,j,k)+umacx(i+1,j,k))
               vel(i,j,k,2) = half*(umacy(i,j,k)+umacy(i,j+1,k))
#if ( AMREX_SPACEDIM == 3 )
               vel(i,j,k,3) = half*(umacz(i,j,k)+umacz(i,j,k+1))
#endif
               do n = 1, SDIM
                  velmin(n) = min(velmin(n),vel(i,j,k,n))
                  velmax(n) = max(velmax(n),vel(i,j,k,n))
               enddo
            enddo
         enddo
      enddo

      if (getForceVerbose.gt.0) then
         call bl_pd_is_ioproc(isioproc)
         if (isioproc.eq.1) then
            do n = 1, SDIM
               write (6,*) "mac velmin (",n,") = ",velmin(n)
               write (6,*) "mac velmax (",n,") = ",velmax(n)
            enddo
         endif
      endif

   end subroutine FORT_AVERAGE_EDGE_STATES

end module navierstokes_3d_module
