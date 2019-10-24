
#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <DERIVE_F.H>
#include <PROB_NS_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 3


module derive_3d_module
  
  implicit none

  private

  public :: derpresvars, derturbvars, &
#ifdef SUMJET
            derjetpresvars, derjetvars, &
#endif
            dermodgradrho, derudotlapu, derkeng, derlogs, dermvel, &
            derlgrhodust, derdmag, &
            dervortx, dervorty, dervortz, &
            dergrdp, &
            derradvel, derazivel, derxvelrot, deryvelrot, dermagvelrot, &
            dermagvortrot, &
#if defined(DO_IAMR_FORCE) 
            derforcing, derforcex, derforcey, derforcez, &    
#endif
            dernull

contains

!c     -----------------------------------------------------------
!c     This file contains functions which compute derived quantities.  
!c     All of the argument lists have the same template, shown below
!c     
!c     INPUTS/OUTPUTS:
!c     
!c     e         <= the quantity derived
!c     DIMS(e)   => index extent of e array
!c     nv        => number of components in e array (should be 1)
!c     dat       => data neded to derive e
!c     DIMS(dat) => index limits of dat array
!c     ncomp     => number of components of dat array (3)
!c     lo,hi     => subrange of e array where result is requested
!c     domlo,hi  => index extent of problem domain (cell centered)
!c     delta     => cell spacing
!c     xlo       => physical location of lower left hand
!c 	           corner of e array
!c     time      => problem evolution time
!c     bc        => array of bndry types for component values
!c                  valid only if component touches bndry
!c     -----------------------------------------------------------

      subroutine derradvel (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt, &
                              bc,level, grid_no) &
                              bind(C, name= "derradvel")

      implicit none

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      REAL_T     ux, uy
      REAL_T     x, y, r
      REAL_T     hx, hy

#include <probdata.H>

      hx = delta(1)
      hy = delta(2)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            y = xlo(2) + hy * ( dble( j-lo(2) ) + half )
            do i = lo(1), hi(1)
               x = xlo(1) + hx * ( dble( i-lo(1) ) + half )
               r = sqrt(x*x+y*y)
               ux = dat(i,j,k,1)
               uy = dat(i,j,k,2)
               e(i,j,k,1) = (ux*x+uy*y)/r
            end do
         end do
      end do


      end subroutine derradvel

      subroutine derazivel  (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt, &
                             bc,level, grid_no) &
                             bind(C, name="derazivel")

      implicit none

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      REAL_T     ux, uy
      REAL_T     x, y, r
      REAL_T     hx, hy

#include <probdata.H>

      hx = delta(1)
      hy = delta(2)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            y = xlo(2) + hy * ( dble( j-lo(2) ) + half )
            do i = lo(1), hi(1)
               x = xlo(1) + hx * ( dble( i-lo(1) ) + half )
               r = sqrt(x*x+y*y)
               ux = dat(i,j,k,1)
               uy = dat(i,j,k,2)
               e(i,j,k,1) = (uy*x-ux*y)/r
            end do
         end do
      end do

      end subroutine derazivel

      subroutine derxvelrot (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt, &
                              bc,level, grid_no) &
                              bind(C, name="derxvelrot")

      implicit none

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      REAL_T     u, v, w
      REAL_T     y
      REAL_T     hy

#include <probdata.H>

      hy = delta(2)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            y = xlo(2) + hy * ( dble( j-lo(2) ) + half )
            do i = lo(1), hi(1)
               e(i,j,k,1) = dat(i,j,k,1) - omega*y
            end do
         end do
      end do

      end subroutine derxvelrot

      subroutine deryvelrot  (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt, &
                              bc,level, grid_no) &
                              bind(C, name= "deryvelrot")

      implicit none

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      REAL_T     u, v, w
      REAL_T     x
      REAL_T     hx

#include <probdata.H>

      hx = delta(1)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               x = xlo(1) + hx * ( dble( i-lo(1) ) + half )
               e(i,j,k,1) = dat(i,j,k,2) + omega*x
            end do
         end do
      end do

      end subroutine deryvelrot

      subroutine dermagvelrot (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt,  &
                              bc,level, grid_no) &
                              bind(C, name="dermagvelrot")

      implicit none

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      REAL_T     u, v, w
      REAL_T     x, y
      REAL_T     hx, hy

#include <probdata.H>

      hx = delta(1)
      hy = delta(2)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            y = xlo(2) + hy * ( dble( j-lo(2) ) + half )
            do i = lo(1), hi(1)
               x = xlo(1) + hx * ( dble( i-lo(1) ) + half )
               u = dat(i,j,k,1) - omega*y
               v = dat(i,j,k,2) + omega*x
               w = dat(i,j,k,3)
               e(i,j,k,1) = sqrt(u*u+v*v+w*w)
            end do
         end do
      end do

      end subroutine dermagvelrot

      subroutine dermagvortrot (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                               lo,hi,domlo,domhi,delta,xlo,time,dt, &
                               bc,level, grid_no) &
                               bind(C, name="dermagvortrot")

      implicit none

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      REAL_T     u, v, w
      REAL_T     wx, wy, wz
      REAL_T     hx, hy, hz
      REAL_T     dudy, dudz
      REAL_T     dvdx, dvdz
      REAL_T     dwdx, dwdy
      REAL_T     twoihx, twoihy, twoihz

#include <probdata.H>

      hx     = delta(1)
      hy     = delta(2)
      hz     = delta(3)
      twoihx = 1.0d0/(two*hx)
      twoihy = 1.0d0/(two*hy)
      twoihz = 1.0d0/(two*hz)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               u   = dat(i,j,k,1)
               v   = dat(i,j,k,2)
               w   = dat(i,j,k,3)

               dudy = (dat(i,j+1,k,1)-dat(i,j-1,k,1))*twoihy
               dudz = (dat(i,j,k+1,1)-dat(i,j,k-1,1))*twoihz

               dvdx = (dat(i+1,j,k,2)-dat(i-1,j,k,2))*twoihx
               dvdz = (dat(i,j,k+1,2)-dat(i,j,k-1,2))*twoihz

               dwdx = (dat(i+1,j,k,3)-dat(i-1,j,k,3))*twoihx
               dwdy = (dat(i,j+1,k,3)-dat(i,j-1,k,3))*twoihy

               wx = dwdy - dvdz
               wy = dudz - dwdx
               wz = dvdx - dudy - two*omega

               e(i,j,k,1) = sqrt(wx*wx+wy*wy+wz*wz)

            end do
         end do
      end do


      end subroutine dermagvortrot

#if defined(DO_IAMR_FORCE) 
      subroutine derforcing (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt, &
                             bc,level,grid_no) &
                             bind(C, name="derforcing")
!c
!c     This routine will derive the energy being injected by the
!c     forcing term used for generating turbulence in probtype 14
!c     Requires velocity field, time, and the right parameters
!c     for the forcing term, i.e. probin, *somehow*
!c

      implicit none

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no, isioproc

      integer    i,j,k,n
      REAL_T     twicePi
      REAL_T     x, y, z
      REAL_T     hx, hy, hz
      REAL_T     u, v, w
      REAL_T     f1, f2, f3

      integer    kx, ky, kz, xstep, ystep, zstep
      REAL_T     kxd, kyd, kzd
      REAL_T     xt, yt, zt
      REAL_T     zlo, infl_time

      REAL_T     Lx, Ly, Lz, Lmin, HLx, HLy, HLz
      REAL_T     kappa, kappaMax

#include <probdata.H>
#include <forcedata.H>

      call bl_pd_is_ioproc(isioproc)

      hx = delta(1)
      hy = delta(2)
      hz = delta(3)

      twicePi=two*Pi

      if (probtype.eq.14.or.probtype.eq.15) then
!c     Homogeneous Isotropic Turbulence or Inflow

!c     Adjust z offset for probtype 15
         if (probtype.eq.15.and.infl_time_offset.gt.(-half)) then
            infl_time = time + infl_time_offset
            zlo = xlo(3) - (time*adv_vel)
         else
            if (time_offset.gt.zero) then
               infl_time = time + time_offset
            else
               infl_time = time
            endif
            zlo = xlo(3)
         endif
         
         if (probtype.eq.14) then
            Lx = domnhi(1)-domnlo(1)
            Ly = domnhi(2)-domnlo(2)
            Lz = domnhi(3)-domnlo(3)
         else if (probtype.eq.15) then
            Lx = forcing_xlength
            Ly = forcing_ylength
            Lz = forcing_zlength
         endif

         if (hack_lz.eq.1) then
            Lz = Lz/two
         endif

         Lmin = min(Lx,Ly,Lz)
         kappaMax = dble(nmodes)/Lmin + 1.0d-8
         nxmodes = nmodes*int(0.5d0+Lx/Lmin)
         nymodes = nmodes*int(0.5d0+Ly/Lmin)
         nzmodes = nmodes*int(0.5d0+Lz/Lmin)
         xstep = int(Lx/Lmin+0.5d0)
         ystep = int(Ly/Lmin+0.5d0)
         zstep = int(Lz/Lmin+0.5d0)

         if (forcing_twice_wavelength.eq.1) then
            HLx = Lx/two
            HLy = Ly/two
            HLz = Lz/two
         else
            HLx = Lx
            HLy = Ly
            HLz = Lz
         endif

         do k = lo(3), hi(3)
            z = zlo + hz * ( dble( k-lo(3) ) + half )
            do j = lo(2), hi(2)
               y = xlo(2) + hy * ( dble( j-lo(2) ) + half )
               do i = lo(1), hi(1)
                  x = xlo(1) + hx * ( dble( i-lo(1) ) + half )

                  f1 = zero
                  f2 = zero
                  f3 = zero
                  
                  do kz = mode_start*zstep, nzmodes, zstep
                     kzd = dble(kz)
                     do ky = mode_start*ystep, nymodes, ystep
                        kyd = dble(ky)
                        do kx = mode_start*xstep, nxmodes, xstep
                           kxd = dble(kx)
                           kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                           if (kappa.le.kappaMax) then
                              xT = cos(FTX(kx,ky,kz)*infl_time+TAT(kx,ky,kz))
                              if (div_free_force.eq.1) then
                                 f1 = f1 + xT * ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) &
                                     -           FAY(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )
                                 f2 = f2 + xT * ( FAX(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) &
                                     -           FAZ(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )
                                 f3 = f3 + xT * ( FAY(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) &
                                     -           FAX(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )
                              else
                                 f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                                 f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                                 f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                              endif
                           endif
                        enddo
                     enddo
                  enddo
                  do kz = 1, zstep - 1
                     kzd = dble(kz)
                     do ky = mode_start, nymodes
                        kyd = dble(ky)
                        do kx = mode_start, nxmodes
                           kxd = dble(kx)
                           kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                           if (kappa.le.kappaMax) then
                              xT = cos(FTX(kx,ky,kz)*infl_time+TAT(kx,ky,kz))
                              if (div_free_force.eq.1) then
                                 f1 = f1 + xT * ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) &
                                     -           FAY(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )
                                 f2 = f2 + xT * ( FAX(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) &
                                     -           FAZ(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )
                                 f3 = f3 + xT * ( FAY(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) &
                                     -           FAX(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )
                              else
                                 f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                                 f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                                 f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                              endif
                           endif
                        enddo
                     enddo
                  enddo

                  u   = dat(i,j,k,2)
                  v   = dat(i,j,k,3)
                  w   = dat(i,j,k,4)
                  
                  if (use_rho_in_forcing.eq.1) then 
                     e(i,j,k,1) = dat(i,j,k,1) * ( u*f1 + v*f2 + w*f3 )
                  else
                     e(i,j,k,1) = u*f1 + v*f2 + w*f3
                  endif
                  
               end do
            end do
         end do

      endif
      end subroutine derforcing


      subroutine derforcex (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                            lo,hi,domlo,domhi,delta,xlo,time,dt, &
                            bc,level,grid_no) &
                            bind(C, name="derforcex")
!c     
!c     This routine will derive the energy being injected by the
!c     forcing term used for generating turbulence in probtype 14
!c     Requires velocity field, time, and the right parameters
!c     for the forcing term, i.e. probin, *somehow*
!c     

      implicit none

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no, isioproc

      integer    i,j,k,n
      REAL_T     twicePi
      REAL_T     x, y, z
      REAL_T     hx, hy, hz
      REAL_T     f1

      integer    kx, ky, kz, xstep, ystep, zstep
      REAL_T     kxd, kyd, kzd
      REAL_T     xt, yt, zt
      REAL_T     zlo, infl_time

      REAL_T     Lx, Ly, Lz, Lmin, HLx, HLy, HLz
      REAL_T     kappa, kappaMax

#include <probdata.H>
#include <forcedata.H>

      call bl_pd_is_ioproc(isioproc)

      hx = delta(1)
      hy = delta(2)
      hz = delta(3)

      twicePi=two*Pi

      if (probtype.eq.14.or.probtype.eq.15) then
!c     Homogeneous Isotropic Turbulence or Inflow

!c     Adjust z offset for probtype 15
         if (probtype.eq.15.and.infl_time_offset.gt.(-half)) then
            infl_time = time + infl_time_offset
            zlo = xlo(3) - (time*adv_vel)
         else
            if (time_offset.gt.zero) then
               infl_time = time + time_offset
            else
               infl_time = time
            endif
            zlo = xlo(3)
         endif

         if (probtype.eq.14) then
            Lx = domnhi(1)-domnlo(1)
            Ly = domnhi(2)-domnlo(2)
            Lz = domnhi(3)-domnlo(3)
         else if (probtype.eq.15) then
            Lx = forcing_xlength
            Ly = forcing_ylength
            Lz = forcing_zlength
         endif
         if (hack_lz.eq.1) then
            Lz = Lz/two
         endif
         Lmin = min(Lx,Ly,Lz)
         kappaMax = dble(nmodes)/Lmin + 1.0d-8
         nxmodes = nmodes*int(0.5d0+Lx/Lmin)
         nymodes = nmodes*int(0.5d0+Ly/Lmin)
         nzmodes = nmodes*int(0.5d0+Lz/Lmin)
         xstep = int(Lx/Lmin+0.5d0)
         ystep = int(Ly/Lmin+0.5d0)
         zstep = int(Lz/Lmin+0.5d0)
         if (forcing_twice_wavelength.eq.1) then
            HLx = Lx/two
            HLy = Ly/two
            HLz = Lz/two
         else
            HLx = Lx
            HLy = Ly
            HLz = Lz
         endif

         do k = lo(3), hi(3)
            z = zlo + hz * ( dble( k-lo(3) ) + half )
            do j = lo(2), hi(2)
               y = xlo(2) + hy * ( dble( j-lo(2) ) + half )
               do i = lo(1), hi(1)
                  x = xlo(1) + hx * ( dble( i-lo(1) ) + half )
                  
                  f1 = zero
                  
                  do kz = mode_start*zstep, nzmodes, zstep
                     kzd = dble(kz)
                     do ky = mode_start*ystep, nymodes, ystep
                        kyd = dble(ky)
                        do kx = mode_start*xstep, nxmodes, xstep
                           kxd = dble(kx)
                           kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                           if (kappa.le.kappaMax) then
                              xT = cos(FTX(kx,ky,kz)*infl_time+TAT(kx,ky,kz))
                              if (div_free_force.eq.1) then
                                 f1 = f1 + xT * ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) &
                                     -           FAY(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )
                              else
                                 f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                              endif
                           endif
                        enddo
                     enddo
                  enddo
                  do kz = 1, zstep - 1
                     kzd = dble(kz)
                     do ky = mode_start, nymodes
                        kyd = dble(ky)
                        do kx = mode_start, nxmodes
                           kxd = dble(kx)
                           kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                           if (kappa.le.kappaMax) then
                              xT = cos(FTX(kx,ky,kz)*infl_time+TAT(kx,ky,kz))
                              if (div_free_force.eq.1) then
                                 f1 = f1 + xT * ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) &
                                     -           FAY(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )
                              else
                                 f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                              endif
                           endif
                        enddo
                     enddo
                  enddo
                  
                  if (use_rho_in_forcing.eq.1) then 
                     e(i,j,k,1) = f1 * dat(i,j,k,1)
                  else
                     e(i,j,k,1) = f1
                  endif
               enddo
            enddo
         enddo

      endif

      end subroutine derforcex



      subroutine derforcey (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
          lo,hi,domlo,domhi,delta,xlo,time,dt, &
          bc,level,grid_no) &
          bind(C, name="derforcey")
!c     
!c     This routine will derive the energy being injected by the
!c     forcing term used for generating turbulence in probtype 14
!c     Requires velocity field, time, and the right parameters
!c     for the forcing term, i.e. probin, *somehow*
!c     

      implicit none

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no, isioproc

      integer    i,j,k,n
      REAL_T     twicePi
      REAL_T     x, y, z
      REAL_T     hx, hy, hz
      REAL_T     f2

      integer    kx, ky, kz, xstep, ystep, zstep
      REAL_T     kxd, kyd, kzd
      REAL_T     xt, yt, zt
      REAL_T     zlo, infl_time

      REAL_T     Lx, Ly, Lz, Lmin, HLx, HLy, HLz
      REAL_T     kappa, kappaMax

#include <probdata.H>
#include <forcedata.H>

      call bl_pd_is_ioproc(isioproc)

      hx = delta(1)
      hy = delta(2)
      hz = delta(3)

      twicePi=two*Pi

      if (probtype.eq.14.or.probtype.eq.15) then
!c     Homogeneous Isotropic Turbulence or Inflow

!c     Adjust z offset for probtype 15
         if (probtype.eq.15.and.infl_time_offset.gt.(-half)) then
            infl_time = time + infl_time_offset
            zlo = xlo(3) - (time*adv_vel)
         else
            if (time_offset.gt.zero) then
               infl_time = time + time_offset
            else
               infl_time = time
            endif
            zlo = xlo(3)
         endif

         if (probtype.eq.14) then
            Lx = domnhi(1)-domnlo(1)
            Ly = domnhi(2)-domnlo(2)
            Lz = domnhi(3)-domnlo(3)
         else if (probtype.eq.15) then
            Lx = forcing_xlength
            Ly = forcing_ylength
            Lz = forcing_zlength
         endif
         if (hack_lz.eq.1) then
            Lz = Lz/two
         endif
         Lmin = min(Lx,Ly,Lz)
         kappaMax = dble(nmodes)/Lmin + 1.0d-8
         nxmodes = nmodes*int(0.5d0+Lx/Lmin)
         nymodes = nmodes*int(0.5d0+Ly/Lmin)
         nzmodes = nmodes*int(0.5d0+Lz/Lmin)
         xstep = int(Lx/Lmin+0.5d0)
         ystep = int(Ly/Lmin+0.5d0)
         zstep = int(Lz/Lmin+0.5d0)
         if (forcing_twice_wavelength.eq.1) then
            HLx = Lx/two
            HLy = Ly/two
            HLz = Lz/two
         else
            HLx = Lx
            HLy = Ly
            HLz = Lz
         endif

         do k = lo(3), hi(3)
            z = zlo + hz * ( dble( k-lo(3) ) + half )
            do j = lo(2), hi(2)
               y = xlo(2) + hy * ( dble( j-lo(2) ) + half )
               do i = lo(1), hi(1)
                  x = xlo(1) + hx * ( dble( i-lo(1) ) + half )
                  
                  f2 = zero
                  
                  do kz = mode_start*zstep, nzmodes, zstep
                     kzd = dble(kz)
                     do ky = mode_start*ystep, nymodes, ystep
                        kyd = dble(ky)
                        do kx = mode_start*xstep, nxmodes, xstep
                           kxd = dble(kx)
                           kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                           if (kappa.le.kappaMax) then
                              xT = cos(FTX(kx,ky,kz)*infl_time+TAT(kx,ky,kz))
                              if (div_free_force.eq.1) then
                                 f2 = f2 + xT * ( FAX(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) &
                                     -           FAZ(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )
                              else
                                 f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                              endif
                           endif
                        enddo
                     enddo
                  enddo
                  do kz = 1, zstep - 1
                     kzd = dble(kz)
                     do ky = mode_start, nymodes
                        kyd = dble(ky)
                        do kx = mode_start, nxmodes
                           kxd = dble(kx)
                           kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                           if (kappa.le.kappaMax) then
                              xT = cos(FTX(kx,ky,kz)*infl_time+TAT(kx,ky,kz))
                              if (div_free_force.eq.1) then
                                 f2 = f2 + xT * ( FAX(kx,ky,kz)*twicePi*(kzd/HLz)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) &
                                     -           FAZ(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )
                              else
                                 f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                              endif
                           endif
                        enddo
                     enddo
                  enddo
                  
                  if (use_rho_in_forcing.eq.1) then 
                     e(i,j,k,1) = f2 * dat(i,j,k,1)
                  else
                     e(i,j,k,1) = f2
                  endif
               enddo
            enddo
         enddo

      endif

      end subroutine derforcey



      subroutine derforcez (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                            lo,hi,domlo,domhi,delta,xlo,time,dt, &
                            bc,level,grid_no) &
                            bind(C, name="derforcez")
!c     
!c     This routine will derive the energy being injected by the
!c     forcing term used for generating turbulence in probtype 14
!c     Requires velocity field, time, and the right parameters
!c     for the forcing term, i.e. probin, *somehow*
!c     

      implicit none

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no, isioproc

      integer    i,j,k,n
      REAL_T     twicePi
      REAL_T     x, y, z
      REAL_T     hx, hy, hz
      REAL_T     rho, u, v, w
      REAL_T     f3

      integer    kx, ky, kz, xstep, ystep, zstep
      REAL_T     kxd, kyd, kzd
      REAL_T     xt, yt, zt
      REAL_T     zlo, infl_time

      REAL_T     Lx, Ly, Lz, Lmin, HLx, HLy, HLz
      REAL_T     kappa, kappaMax

#include <probdata.H>
#include <forcedata.H>

      call bl_pd_is_ioproc(isioproc)

      hx = delta(1)
      hy = delta(2)
      hz = delta(3)

      twicePi=two*Pi
      
      if (probtype.eq.14.or.probtype.eq.15) then
!c     Homogeneous Isotropic Turbulence or Inflow
         
!c     Adjust z offset for probtype 15
         if (probtype.eq.15.and.infl_time_offset.gt.(-half)) then
            infl_time = time + infl_time_offset
            zlo = xlo(3) - (time*adv_vel)
         else
            if (time_offset.gt.zero) then
               infl_time = time + time_offset
            else
               infl_time = time
            endif
            zlo = xlo(3)
         endif

         if (probtype.eq.14) then
            Lx = domnhi(1)-domnlo(1)
            Ly = domnhi(2)-domnlo(2)
            Lz = domnhi(3)-domnlo(3)
         else if (probtype.eq.15) then
            Lx = forcing_xlength
            Ly = forcing_ylength
            Lz = forcing_zlength
         endif
         if (hack_lz.eq.1) then
            Lz = Lz/two
         endif
         Lmin = min(Lx,Ly,Lz)
         kappaMax = dble(nmodes)/Lmin + 1.0d-8
         nxmodes = nmodes*int(0.5d0+Lx/Lmin)
         nymodes = nmodes*int(0.5d0+Ly/Lmin)
         nzmodes = nmodes*int(0.5d0+Lz/Lmin)
         xstep = int(Lx/Lmin+0.5d0)
         ystep = int(Ly/Lmin+0.5d0)
         zstep = int(Lz/Lmin+0.5d0)
         if (forcing_twice_wavelength.eq.1) then
            HLx = Lx/two
            HLy = Ly/two
            HLz = Lz/two
         else
            HLx = Lx
            HLy = Ly
            HLz = Lz
         endif

         do k = lo(3), hi(3)
            z = zlo + hz * ( dble( k-lo(3) ) + half )
            do j = lo(2), hi(2)
               y = xlo(2) + hy * ( dble( j-lo(2) ) + half )
               do i = lo(1), hi(1)
                  x = xlo(1) + hx * ( dble( i-lo(1) ) + half )
                  
                  f3 = zero
                  
                  do kz = mode_start*zstep, nzmodes, zstep
                     kzd = dble(kz)
                     do ky = mode_start*ystep, nymodes, ystep
                        kyd = dble(ky)
                        do kx = mode_start*xstep, nxmodes, xstep
                           kxd = dble(kx)
                           kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                           if (kappa.le.kappaMax) then
                              xT = cos(FTX(kx,ky,kz)*infl_time+TAT(kx,ky,kz))
                              if (div_free_force.eq.1) then
                                 f3 = f3 + xT * ( FAY(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) &
                                     -           FAX(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )
                              else
                                 f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                              endif
                           endif
                        enddo
                     enddo
                  enddo
                  do kz = 1, zstep - 1
                     kzd = dble(kz)
                     do ky = mode_start, nymodes
                        kyd = dble(ky)
                        do kx = mode_start, nxmodes
                           kxd = dble(kx)
                           kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                           if (kappa.le.kappaMax) then
                              xT = cos(FTX(kx,ky,kz)*infl_time+TAT(kx,ky,kz))
                              if (div_free_force.eq.1) then
                                 f3 = f3 + xT * ( FAY(kx,ky,kz)*twicePi*(kxd/HLx)*cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) &
                                     -           FAX(kx,ky,kz)*twicePi*(kyd/HLy)*sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )
                              else
                                 f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                              endif
                           endif
                        enddo
                     enddo
                  enddo
                  
                  if (use_rho_in_forcing.eq.1) then 
                     e(i,j,k,1) = f3 * dat(i,j,k,1)
                  else
                     e(i,j,k,1) = f3
                  endif
               enddo
            enddo
         enddo

      endif

      end subroutine derforcez

#endif /*defined(DO_IAMR_FORCE) */

      subroutine derpresvars (e,DIMS(gp),nv,dat,DIMS(dat),ncomp, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt, &
                              bc,level,grid_no) &
                              bind(C, name="derpresvars")
!c
!c     This routine computes cell-centered pressure as average of the eight
!c     surrounding nodal values, along with the three gradients.
!c
      implicit none

      integer DIMDEC(gp)
      integer DIMDEC(dat)
      integer nv, ncomp
      REAL_T  e(DIMV(gp),nv)
      REAL_T  dat(DIMV(dat))
      integer lo(SDIM), hi(SDIM)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM)
      REAL_T  xlo(SDIM)
      REAL_T  time, dt, dx, dy, dz
      integer bc(SDIM,2,ncomp)
      integer level
      integer grid_no

      integer i,j,k
 
      dx = fourth/delta(1)
      dy = fourth/delta(2)
      dz = fourth/delta(3)

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            ! Average pressure
            e(i,j,k,1) = eighth*( &
                 dat(i+1,j,k)     + dat(i,j,k)     + &
                 dat(i+1,j+1,k)   + dat(i,j+1,k)   + &
                 dat(i+1,j,k+1)   + dat(i,j,k+1)   + &
                 dat(i+1,j+1,k+1) + dat(i,j+1,k+1) )
          end do
        end do
      end do

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            ! X Gradient
            e(i,j,k,2) = dx*( &
                 dat(i+1,j,k)     - dat(i,j,k)     + &
                 dat(i+1,j+1,k)   - dat(i,j+1,k)   + &
                 dat(i+1,j,k+1)   - dat(i,j,k+1)   + &
                 dat(i+1,j+1,k+1) - dat(i,j+1,k+1) )
          end do
        end do
      end do

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            ! Y Gradient
            e(i,j,k,3) = dy*( &
                 dat(i,j+1,k)     - dat(i,j,k)     + &
                 dat(i+1,j+1,k)   - dat(i+1,j,k)   + &
                 dat(i,j+1,k+1)   - dat(i,j,k+1)   + &
                 dat(i+1,j+1,k+1) - dat(i+1,j,k+1) )
          end do
        end do
      end do

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            ! Z Gradient
            e(i,j,k,4) = dz*( &
                 dat(i,j,k+1)     - dat(i,j,k)     + &
                 dat(i+1,j,k+1)   - dat(i+1,j,k)   + &
                 dat(i,j+1,k+1)   - dat(i,j+1,k)   + &
                 dat(i+1,j+1,k+1) - dat(i+1,j+1,k) )
          end do
        end do
      end do

      end subroutine derpresvars

      subroutine derturbvars (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt, &
                              bc,level,grid_no) &
                              bind(C, name="derturbvars")

      implicit none

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      REAL_T     tdx, tdy, tdz

      tdx = 1.0d0/(two*delta(1))
      tdy = 1.0d0/(two*delta(2))
      tdz = 1.0d0/(two*delta(3))

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = dat(i,j,k,1)
            enddo
         enddo
      enddo

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,5) = (dat(i+1,j,k,5)-dat(i-1,j,k,5))*tdx
            enddo
         enddo
      enddo

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,9) = (dat(i,j+1,k,9)-dat(i,j-1,k,9))*tdy
            enddo
         enddo
      enddo

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,13) = (dat(i,j,k+1,13)-dat(i,j,k-1,13))*tdz
            enddo
         enddo
      enddo

      end subroutine derturbvars

#ifdef SUMJET

      subroutine derjetpresvars (e,DIMS(gp),nv,dat,DIMS(dat),ncomp, &
                                 lo,hi,domlo,domhi,delta,xlo,time,dt, &
                                 bc,level,grid_no) &
                                 bind(C, name="derjetpresvars")
!c
!c     This routine computes cell-centered pressure as average of the eight
!c     surrounding nodal values, along with the three gradients.  The gradients
!c     of each of these is also calculated for slope reconstruction in the
!c     jet diagnostics routine.
!c
      implicit none

      integer DIMDEC(gp)
      integer DIMDEC(dat)
      integer nv, ncomp
      REAL_T  e(DIMV(gp),nv)
      REAL_T  dat(DIMV(dat))
      integer lo(SDIM), hi(SDIM)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM)
      REAL_T  xlo(SDIM)
      REAL_T  time, dt, dx, dy, dz, tdx, tdy, tdz
      integer bc(SDIM,2,ncomp)
      integer level
      integer grid_no

      integer i,j,k
 
      dx = fourth/delta(1)
      dy = fourth/delta(2)
      dz = fourth/delta(3)

      tdx = two*delta(1)
      tdy = two*delta(2)
      tdz = two*delta(3)

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)
!c     Average pressure
            e(i,j,k,1) = eighth*( &
                 dat(i+1,j,k)     + dat(i,j,k)     + &
                 dat(i+1,j+1,k)   + dat(i,j+1,k)   + &
                 dat(i+1,j,k+1)   + dat(i,j,k+1)   + &
                 dat(i+1,j+1,k+1) + dat(i,j+1,k+1) )
!c     X Gradient
            e(i,j,k,2) = dx*( &
                 dat(i+1,j,k)     - dat(i,j,k)     + &
                 dat(i+1,j+1,k)   - dat(i,j+1,k)   + &
                 dat(i+1,j,k+1)   - dat(i,j,k+1)   + &
                 dat(i+1,j+1,k+1) - dat(i,j+1,k+1) )
!c     Y Gradient
            e(i,j,k,3) = dy*( &
                 dat(i,j+1,k)     - dat(i,j,k)     + &
                 dat(i+1,j+1,k)   - dat(i+1,j,k)   + &
                 dat(i,j+1,k+1)   - dat(i,j,k+1)   + &
                 dat(i+1,j+1,k+1) - dat(i+1,j,k+1) )
!c     Z Gradient
            e(i,j,k,4) = dz*( &
                 dat(i,j,k+1)     - dat(i,j,k)     + &
                 dat(i+1,j,k+1)   - dat(i+1,j,k)   + &
                 dat(i,j+1,k+1)   - dat(i,j+1,k)   + &
                 dat(i+1,j+1,k+1) - dat(i+1,j+1,k) )

!c     Average pressure X Gradient
            e(i,j,k,5) = (eighth*( &
                 dat(i+2,j,k)     + dat(i+1,j,k)   + &
                 dat(i+2,j+1,k)   + dat(i+1,j+1,k) + &
                 dat(i+2,j,k+1)   + dat(i+1,j,k+1) + &
                 dat(i+2,j+1,k+1) + dat(i+1,j+1,k+1) ) &
                -eighth*(  &
                 dat(i,j,k)     + dat(i-1,j,k)     + &
                 dat(i,j+1,k)   + dat(i-1,j+1,k)   + &
                 dat(i,j,k+1)   + dat(i-1,j,k+1)   + &
                 dat(i,j+1,k+1) + dat(i-1,j+1,k+1) ) ) / tdx
!c     X Gradient X Gradient
            e(i,j,k,6) = (dx*( &
                 dat(i+2,j,k)     - dat(i+1,j,k)   + &
                 dat(i+2,j+1,k)   - dat(i+1,j+1,k) + &
                 dat(i+2,j,k+1)   - dat(i+1,j,k+1) + &
                 dat(i+2,j+1,k+1) - dat(i+1,j+1,k+1) ) &
                -dx*( &
                 dat(i,j,k)     - dat(i-1,j,k)     + &
                 dat(i,j+1,k)   - dat(i-1,j+1,k)   + &
                 dat(i,j,k+1)   - dat(i-1,j,k+1)   + &
                 dat(i,j+1,k+1) - dat(i-1,j+1,k+1) ) ) / tdx
!c     Y Gradient X Gradient
            e(i,j,k,7) = (dy*( &
                 dat(i+1,j+1,k)   - dat(i+1,j,k)   + &
                 dat(i+2,j+1,k)   - dat(i+2,j,k)   + &
                 dat(i+1,j+1,k+1) - dat(i+1,j,k+1) + &
                 dat(i+2,j+1,k+1) - dat(i+2,j,k+1) ) &
                -dy*( &
                 dat(i-1,j+1,k)   - dat(i-1,j,k)   + &
                 dat(i,j+1,k)     - dat(i,j,k)     + &
                 dat(i-1,j+1,k+1) - dat(i-1,j,k+1) + &
                 dat(i,j+1,k+1)   - dat(i,j,k+1) ) ) / tdx
!c     Z Gradient X Gradient
            e(i,j,k,8) = (dz*( &
                 dat(i+1,j,k+1)   - dat(i+1,j,k)   + &
                 dat(i+2,j,k+1)   - dat(i+2,j,k)   + &
                 dat(i+1,j+1,k+1) - dat(i+1,j+1,k) + &
                 dat(i+2,j+1,k+1) - dat(i+2,j+1,k) ) &
                -dz*( &
                 dat(i-1,j,k+1)   - dat(i-1,j,k)   + &
                 dat(i,j,k+1)     - dat(i,j,k)     + &
                 dat(i-1,j+1,k+1) - dat(i-1,j+1,k) + &
                 dat(i,j+1,k+1)   - dat(i,j+1,k) ) ) / tdx

!c     Average pressure Y Gradient
            e(i,j,k,9) = ( eighth*( &
                 dat(i+1,j+1,k)   + dat(i,j+1,k)   + &
                 dat(i+1,j+2,k)   + dat(i,j+2,k)   + &
                 dat(i+1,j+1,k+1) + dat(i,j+1,k+1) + &
                 dat(i+1,j+2,k+1) + dat(i,j+2,k+1) ) &
                -eighth*(  &
                 dat(i+1,j-1,k)   + dat(i,j-1,k)   + &
                 dat(i+1,j,k)     + dat(i,j,k)     + &
                 dat(i+1,j-1,k+1) + dat(i,j-1,k+1) + &
                 dat(i+1,j,k+1)   + dat(i,j,k+1) ) ) / tdy
!c     X Gradient Y Gradient
            e(i,j,k,10) = (dx*( &
                 dat(i+1,j+1,k)   - dat(i,j+1,k)   + &
                 dat(i+1,j+2,k)   - dat(i,j+2,k)   + &
                 dat(i+1,j+1,k+1) - dat(i,j+1,k+1) + &
                 dat(i+1,j+2,k+1) - dat(i,j+2,k+1) ) &
                -dx*( &
                 dat(i+1,j-1,k)   - dat(i,j-1,k)   + &
                 dat(i+1,j,k)     - dat(i,j,k)     + &
                 dat(i+1,j-1,k+1) - dat(i,j-1,k+1) + &
                 dat(i+1,j,k+1)   - dat(i,j,k+1) ) ) / tdy
!c     Y Gradient Y Gradient
            e(i,j,k,11) = (dy*( &
                 dat(i,j+2,k)     - dat(i,j+1,k)   + &
                 dat(i+1,j+2,k)   - dat(i+1,j+1,k) + &
                 dat(i,j+2,k+1)   - dat(i,j+1,k+1) + &
                 dat(i+1,j+2,k+1) - dat(i+1,j+1,k+1) ) &
                -dy*( &
                 dat(i,j,k)     - dat(i,j-1,k)     + &
                 dat(i+1,j,k)   - dat(i+1,j-1,k)   + &
                 dat(i,j,k+1)   - dat(i,j-1,k+1)   + &
                 dat(i+1,j,k+1) - dat(i+1,j-1,k+1) ) ) / tdy
!c     Z Gradient Y Gradient
            e(i,j,k,12) = (dz*( &
                 dat(i,j+1,k+1)   - dat(i,j+1,k)   + &
                 dat(i+1,j+1,k+1) - dat(i+1,j+1,k) + &
                 dat(i,j+2,k+1)   - dat(i,j+2,k)   + &
                 dat(i+1,j+2,k+1) - dat(i+1,j+2,k) ) &
                -dz*( &
                 dat(i,j-1,k+1)   - dat(i,j-1,k)   + &
                 dat(i+1,j-1,k+1) - dat(i+1,j-1,k) + &
                 dat(i,j,k+1)     - dat(i,j,k)     + &
                 dat(i+1,j,k+1)   - dat(i+1,j,k) ) ) / tdy

!c     Average pressure Z Gradient
            e(i,j,k,13) = ( eighth*(  &
                 dat(i+1,j,k+1)   + dat(i,j,k+1)   + &
                 dat(i+1,j+1,k+1) + dat(i,j+1,k+1) + &
                 dat(i+1,j,k+2)   + dat(i,j,k+2)   + &
                 dat(i+1,j+1,k+2) + dat(i,j+1,k+2) ) &
                -eighth*(  &
                 dat(i+1,j,k-1)   + dat(i,j,k-1)   + &
                 dat(i+1,j+1,k-1) + dat(i,j+1,k-1) + &
                 dat(i+1,j,k)     + dat(i,j,k)     + &
                 dat(i+1,j+1,k)   + dat(i,j+1,k) ) ) / tdz
!c     X Gradient Z Gradient
            e(i,j,k,14) = (dx*( &
                 dat(i+1,j,k+1)   - dat(i,j,k+1)   + &
                 dat(i+1,j+1,k+1) - dat(i,j+1,k+1) + &
                 dat(i+1,j,k+2)   - dat(i,j,k+2)   + &
                 dat(i+1,j+1,k+2) - dat(i,j+1,k+2) ) &
                -dx*( &
                 dat(i+1,j,k-1)   - dat(i,j,k-1)   + &
                 dat(i+1,j+1,k-1) - dat(i,j+1,k-1) + &
                 dat(i+1,j,k)     - dat(i,j,k)     + &
                 dat(i+1,j+1,k)   - dat(i,j+1,k) ) ) / tdz
!c     Y Gradient Z Gradient
            e(i,j,k,15) = (dy*( &
                 dat(i,j+1,k+1)   - dat(i,j,k+1)   + &
                 dat(i+1,j+1,k+1) - dat(i+1,j,k+1) + &
                 dat(i,j+1,k+2)   - dat(i,j,k+2)   + &
                 dat(i+1,j+1,k+2) - dat(i+1,j,k+2) ) &
                -dy*( &
                 dat(i,j+1,k-1)   - dat(i,j,k-1)   + &
                 dat(i+1,j+1,k-1) - dat(i+1,j,k-1) + &
                 dat(i,j+1,k)     - dat(i,j,k)     + &
                 dat(i+1,j+1,k)   - dat(i+1,j,k) ) ) / tdz
!c     Z Gradient Z Gradient
            e(i,j,k,16) = (dz*( &
                 dat(i,j,k+2)     - dat(i,j,k+1)   + &
                 dat(i+1,j,k+2)   - dat(i+1,j,k+1) + &
                 dat(i,j+1,k+2)   - dat(i,j+1,k+1) + &
                 dat(i+1,j+1,k+2) - dat(i+1,j+1,k+1) ) &
                -dz*( &
                 dat(i,j,k)     - dat(i,j,k-1)     + &
                 dat(i+1,j,k)   - dat(i+1,j,k-1)   + &
                 dat(i,j+1,k)   - dat(i,j+1,k-1)   + &
                 dat(i+1,j+1,k) - dat(i+1,j+1,k-1) ) ) / tdz

          end do
        end do
      end do

      end subroutine derjetpresvars

      subroutine derjetvars (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt, &
                             bc,level,grid_no)&
                             bind(C, name="derjetvars")

      implicit none

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k,n,nvars,do_trac2
      REAL_T     hx, hy, hz
      REAL_T     tdx, tdy, tdz
      REAL_T     dxt, dyt, dzt

#include <probdata.H>

      hx = delta(1)
      hy = delta(2)
      hz = delta(3)

      tdx = 1.0d0/(two*hx)
      tdy = 1.0d0/(two*hy)
      tdz = 1.0d0/(two*hz)
      dxt = 1.0d0/(hx*hx)
      dyt = 1.0d0/(hy*hy)
      dzt = 1.0d0/(hz*hz)

      call bl_ns_dotrac2(do_trac2)

      if (do_trac2.eq.1) then
         nvars=6
      else
         nvars=5
      endif

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               do n = 1, nvars
                  e(i,j,k,n+ 0) = dat(i,j,k,n)
                  e(i,j,k,n+ 6) = (dat(i+1,j,k,n)-dat(i-1,j,k,n))*tdx
                  e(i,j,k,n+12) = (dat(i,j+1,k,n)-dat(i,j-1,k,n))*tdy
                  e(i,j,k,n+18) = (dat(i,j,k+1,n)-dat(i,j,k-1,n))*tdz
               enddo
               ! Need these for D = u.nabla^2 u
               e(i,j,k,25) = (dat(i+1,j,k,2)-two*dat(i,j,k,2)+dat(i-1,j,k,2))*dxt
               e(i,j,k,26) = (dat(i,j+1,k,3)-two*dat(i,j,k,3)+dat(i,j-1,k,3))*dyt
               e(i,j,k,27) = (dat(i,j,k+1,4)-two*dat(i,j,k,4)+dat(i,j,k-1,4))*dzt

               if (do_trac2.eq.0) then
                  n = 6
                  e(i,j,k,n+ 0) = zero
                  e(i,j,k,n+ 6) = zero
                  e(i,j,k,n+12) = zero
                  e(i,j,k,n+18) = zero
               endif

            enddo
         enddo
      enddo

      end subroutine derjetvars

#endif

      subroutine dermodgradrho (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt, &
                              bc,level,grid_no) bind(C,name="dermodgradrho")
                              
      implicit none
!c
!c     This routine will derive the modulus of the density gradients
!c
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      REAL_T     drdx, drdy, drdz, tdx, tdy, tdz

      tdx = 1.0d0/(two*delta(1))
      tdy = 1.0d0/(two*delta(2))
      tdz = 1.0d0/(two*delta(3))

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               drdx = (dat(i+1,j,k,1)-dat(i-1,j,k,1))*tdx
               drdy = (dat(i,j+1,k,1)-dat(i,j-1,k,1))*tdy
               drdz = (dat(i,j,k+1,1)-dat(i,j,k-1,1))*tdz
               e(i,j,k,1) = sqrt( drdx*drdx + drdy*drdy + drdz*drdz )
            end do
         end do
      end do

      end subroutine dermodgradrho

      subroutine derudotlapu (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt, &
                              bc,level,grid_no) &
                              bind(C, name="derudotlapu")
      implicit none
!c
!c     This routine will derive u dot laplacian u
!c     from the velocity field.
!c
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k

      REAL_T     dx2, dy2, dz2

      dx2 = one/(delta(1)*delta(1))
      dy2 = one/(delta(2)*delta(2))
      dz2 = one/(delta(3)*delta(3))

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = &
                   dat(i,j,k,1) * ( &
                   ( dat(i+1,j,k,1)-two*dat(i,j,k,1)+dat(i-1,j,k,1) ) * dx2 + &
                   ( dat(i,j+1,k,1)-two*dat(i,j,k,1)+dat(i,j-1,k,1) ) * dy2 + &
                   ( dat(i,j,k+1,1)-two*dat(i,j,k,1)+dat(i,j,k-1,1) ) * dz2 ) + &
                   dat(i,j,k,2) * ( & 
                   ( dat(i+1,j,k,2)-two*dat(i,j,k,2)+dat(i-1,j,k,2) ) * dx2 + &
                   ( dat(i,j+1,k,2)-two*dat(i,j,k,2)+dat(i,j-1,k,2) ) * dy2 + &
                   ( dat(i,j,k+1,2)-two*dat(i,j,k,2)+dat(i,j,k-1,2) ) * dz2 ) + &
                   dat(i,j,k,3) * ( &
                   ( dat(i+1,j,k,3)-two*dat(i,j,k,3)+dat(i-1,j,k,3) ) * dx2 + &
                   ( dat(i,j+1,k,3)-two*dat(i,j,k,3)+dat(i,j-1,k,3) ) * dy2 + &
                   ( dat(i,j,k+1,3)-two*dat(i,j,k,3)+dat(i,j,k-1,3) ) * dz2 )
            end do
         end do
      end do

      end subroutine derudotlapu

      subroutine derkeng (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt, &
                              bc,level,grid_no) bind(C,name="derkeng")
      implicit none
!c
!c     This routine will derive kinetic energy from density
!c     and the velocity field.
!c
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      REAL_T     rho, u, v, w

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rho = dat(i,j,k,1)
               u   = dat(i,j,k,2)
               v   = dat(i,j,k,3)
               w   = dat(i,j,k,4)
               e(i,j,k,1) = half*rho*(u**2 + v**2 + w**2)
            end do
         end do
      end do

      end subroutine derkeng

      subroutine derlogs (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                          lo,hi,domlo,domhi,delta,xlo,time,dt, &
                          bc,level, grid_no) bind(C,name="derlogs")
      implicit none
!c
!c     This routine will derive log of given scalar quantity
!c
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      REAL_T     rho
      REAL_T     sml

      parameter (sml = 1.0D-10)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rho = max(dat(i,j,k,1),sml)
               e(i,j,k,1) = log10(rho)
            end do
         end do
      end do

      end subroutine derlogs

      subroutine dermvel (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                          lo,hi,domlo,domhi,delta,xlo,time,dt, &
                          bc,level, grid_no) bind(C, name="dermvel")
      implicit none
!c
!c ::: This routine will derive the magnitude of the velocity field
!c ::: from the velocity field
!c
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      REAL_T     u, v, w

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               u   = dat(i,j,k,1)
               v   = dat(i,j,k,2)
               w   = dat(i,j,k,3)
               e(i,j,k,1) = sqrt(u**2 + v**2 + w**2)
            end do
         end do
      end do

      end subroutine dermvel

      subroutine derlgrhodust (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                                   lo,hi,domlo,domhi,delta,xlo,time,dt, &
                                   bc,level,grid_no)bind(C,name="derlgrhodust")
      implicit none
!c
!c ::: This routine will derive log(RHO*C)
!c
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      REAL_T     dust, small

      parameter (small = 1.0D-10)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               dust = max(small,dat(i,j,k,2)*dat(i,j,k,1))
               e(i,j,k,1) = log10(dust)
            end do
         end do
      end do

      end subroutine derlgrhodust

      subroutine derdmag (dmag,DIMS(dmag),nv,dat,DIMS(dat),ncomp, &
                          lo,hi,domlo,domhi,delta,xlo,time,dt, &
                          bc,level,grid_no) bind(C, name="derdmag")
      implicit none
!c
!c ::: John's weird diagnostic routine ...
!c
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(dmag)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     dmag(DIMV(dmag),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer   i,j,k
      REAL_T    ux, uy, uz, vx, vy, vz, wx, wy, wz, dx, dy, dz
      REAL_T    uxcen, uycen, uzcen, uxlo, uxhi, uylo, uyhi, uzlo, uzhi
      REAL_T    vxcen, vycen, vzcen, vxlo, vxhi, vylo, vyhi, vzlo, vzhi
      REAL_T    wxcen, wycen, wzcen, wxlo, wxhi, wylo, wyhi, wzlo, wzhi
  
      REAL_T    xi1,xi2,xi3,s11,s12,s13,s22,s23,s33,xmag

      logical   fixulo_x, fixvlo_x, fixwlo_x, fixuhi_x, fixvhi_x, fixwhi_x
      logical   fixulo_y, fixvlo_y, fixwlo_y, fixuhi_y, fixvhi_y, fixwhi_y
      logical   fixulo_z, fixvlo_z, fixwlo_z, fixuhi_z, fixvhi_z, fixwhi_z
!c
!c     ::::: some useful macro definitions
!c
#     define U(i,j,k) dat(i,j,k,1)
#     define V(i,j,k) dat(i,j,k,2)
#     define W(i,j,k) dat(i,j,k,3)

#     define ULOX bc(1,1,1)
#     define UHIX bc(1,2,1)
#     define ULOY bc(2,1,1)
#     define UHIY bc(2,2,1)
#     define ULOZ bc(3,1,1)
#     define UHIZ bc(3,2,1)

#     define VLOX bc(1,1,2)
#     define VHIX bc(1,2,2)
#     define VLOY bc(2,1,2)
#     define VHIY bc(2,2,2)
#     define VLOZ bc(3,1,2)
#     define VHIZ bc(3,2,2)

#     define WLOX bc(1,1,3)
#     define WHIX bc(1,2,3)
#     define WLOY bc(2,1,3)
#     define WHIY bc(2,2,3)
#     define WLOZ bc(3,1,3)
#     define WHIZ bc(3,2,3)
!c
!c     ::::: statement functions that implement stencil
!c
      uxcen(i,j,k) = half*(U(i+1,j,k)-U(i-1,j,k))/dx
      uxlo(i,j,k)  = (U(i+1,j,k)+three*U(i,j,k)-four*U(i-1,j,k))/(three*dx)
      uxhi(i,j,k)  =-(U(i-1,j,k)+three*U(i,j,k)-four*U(i+1,j,k))/(three*dx)

      uycen(i,j,k) = half*(U(i,j+1,k)-U(i,j-1,k))/dy
      uylo(i,j,k)  = (U(i,j+1,k)+three*U(i,j,k)-four*U(i,j-1,k))/(three*dy)
      uyhi(i,j,k)  =-(U(i,j-1,k)+three*U(i,j,k)-four*U(i,j+1,k))/(three*dy)

      uzcen(i,j,k) = half*(U(i,j,k+1)-U(i,j,k-1))/dz
      uzlo(i,j,k)  = (U(i,j,k+1)+three*U(i,j,k)-four*U(i,j,k-1))/(three*dz)
      uzhi(i,j,k)  =-(U(i,j,k-1)+three*U(i,j,k)-four*U(i,j,k+1))/(three*dz)

      vxcen(i,j,k) = half*(V(i+1,j,k)-V(i-1,j,k))/dx
      vxlo(i,j,k)  = (V(i+1,j,k)+three*V(i,j,k)-four*V(i-1,j,k))/(three*dx)
      vxhi(i,j,k)  =-(V(i-1,j,k)+three*V(i,j,k)-four*V(i+1,j,k))/(three*dx)

      vycen(i,j,k) = half*(V(i,j+1,k)-V(i,j-1,k))/dy
      vylo(i,j,k)  = (V(i,j+1,k)+three*V(i,j,k)-four*V(i,j-1,k))/(three*dy)
      vyhi(i,j,k)  =-(V(i,j-1,k)+three*V(i,j,k)-four*V(i,j+1,k))/(three*dy)

      vzcen(i,j,k) = half*(V(i,j,k+1)-V(i,j,k-1))/dz
      vzlo(i,j,k)  = (V(i,j,k+1)+three*V(i,j,k)-four*V(i,j,k-1))/(three*dz)
      vzhi(i,j,k)  =-(V(i,j,k-1)+three*V(i,j,k)-four*V(i,j,k+1))/(three*dz)

      wxcen(i,j,k) = half*(W(i+1,j,k)-W(i-1,j,k))/dx
      wxlo(i,j,k)  = (W(i+1,j,k)+three*W(i,j,k)-four*W(i-1,j,k))/(three*dx)
      wxhi(i,j,k)  =-(W(i-1,j,k)+three*W(i,j,k)-four*W(i+1,j,k))/(three*dx)

      wycen(i,j,k) = half*(W(i,j+1,k)-W(i,j-1,k))/dy
      wylo(i,j,k)  = (W(i,j+1,k)+three*W(i,j,k)-four*W(i,j-1,k))/(three*dy)
      wyhi(i,j,k)  =-(W(i,j-1,k)+three*W(i,j,k)-four*W(i,j+1,k))/(three*dy)

      wzcen(i,j,k) = half*(W(i,j,k+1)-W(i,j,k-1))/dz
      wzlo(i,j,k)  = (W(i,j,k+1)+three*W(i,j,k)-four*W(i,j,k-1))/(three*dz)
      wzhi(i,j,k)  =-(W(i,j,k-1)+three*W(i,j,k)-four*W(i,j,k+1))/(three*dz)

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      fixulo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (ULOX .eq. EXT_DIR .or. ULOX .eq. HOEXTRAP) )
      fixuhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (UHIX .eq. EXT_DIR .or. UHIX .eq. HOEXTRAP) )
      fixvlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (VLOX .eq. EXT_DIR .or. VLOX .eq. HOEXTRAP) )
      fixvhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (VHIX .eq. EXT_DIR .or. VHIX .eq. HOEXTRAP) )
      fixwlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (WLOX .eq. EXT_DIR .or. WLOX .eq. HOEXTRAP) )
      fixwhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (WHIX .eq. EXT_DIR .or. WHIX .eq. HOEXTRAP) )

      fixulo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (ULOY .eq. EXT_DIR .or. ULOY .eq. HOEXTRAP) )
      fixuhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (UHIY .eq. EXT_DIR .or. UHIY .eq. HOEXTRAP) )
      fixvlo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (VLOY .eq. EXT_DIR .or. VLOY .eq. HOEXTRAP) )
      fixvhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (VHIY .eq. EXT_DIR .or. VHIY .eq. HOEXTRAP) )
      fixwlo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (WLOY .eq. EXT_DIR .or. WLOY .eq. HOEXTRAP) )
      fixwhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (WHIY .eq. EXT_DIR .or. WHIY .eq. HOEXTRAP) )

      fixulo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (ULOZ .eq. EXT_DIR .or. ULOZ .eq. HOEXTRAP) )
      fixuhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (UHIZ .eq. EXT_DIR .or. UHIZ .eq. HOEXTRAP) )
      fixvlo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (VLOZ .eq. EXT_DIR .or. VLOZ .eq. HOEXTRAP) )
      fixvhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (VHIZ .eq. EXT_DIR .or. VHIZ .eq. HOEXTRAP) )
      fixwlo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (WLOZ .eq. EXT_DIR .or. WLOZ .eq. HOEXTRAP) )
      fixwhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (WHIZ .eq. EXT_DIR .or. WHIZ .eq. HOEXTRAP) )
               
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ux = uxcen(i,j,k)
               uy = uycen(i,j,k)
               uz = uzcen(i,j,k)
               vx = vxcen(i,j,k)
               vy = vycen(i,j,k)
               vz = vzcen(i,j,k)
               wx = wxcen(i,j,k)
               wy = wycen(i,j,k)
               wz = wzcen(i,j,k)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
            end do
         end do
      end do


!c
!c     First do all the faces
!c
      if (fixvlo_x .or. fixwlo_x.or.fixulo_x) then
         i = lo(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
               vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
               wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
               uy = uycen(i,j,k)
               vy = vycen(i,j,k)
               wy = wycen(i,j,k)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wz = wzcen(i,j,k)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
            end do
         end do
      end if

      if (fixvhi_x .or. fixwhi_x.or.fixuhi_x) then
         i = hi(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               ux = merge(uxhi(i,j,k),vxcen(i,j,k),fixuhi_x)
               vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
               wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
               uy = uycen(i,j,k)
               vy = vycen(i,j,k)
               wy = wycen(i,j,k)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wz = wzcen(i,j,k)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
            end do
         end do
      end if

      if (fixulo_y .or. fixwlo_y.or.fixvlo_y) then
         j = lo(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               ux = uxcen(i,j,k)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
               vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
               wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wz = wzcen(i,j,k)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
            end do
         end do
      end if

      if (fixuhi_y .or. fixwhi_y.or.fixvhi_y) then
         j = hi(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               ux = uxcen(i,j,k)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
               vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
               wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wz = vzcen(i,j,k)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
            end do
         end do
      end if

      if (fixulo_z .or. fixvlo_z.or.fixwlo_z) then
         k = lo(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               ux = uxcen(i,j,k)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = uycen(i,j,k)
               vy = vycen(i,j,k)
               wy = wycen(i,j,k)
               uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
               vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
               wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
            end do
         end do
      end if

      if (fixuhi_z .or. fixvhi_z .or. fixwhi_z) then
         k = hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               ux = uxcen(i,j,k)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = uycen(i,j,k)
               vy = vycen(i,j,k)
               wy = wycen(i,j,k)
               uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
               vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
               wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
            end do
         end do
      end if
!c
!c     Next do all the edges
!c
      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixvlo_y .or. fixwlo_y)) then
         i = lo(1)
         j = lo(2)
         do k = lo(3),hi(3)
            ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            wz = wzcen(i,j,k)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
         end do
      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixvlo_y .or. fixwlo_y)) then
         i = hi(1)
         j = lo(2)
         do k = lo(3),hi(3)
            ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            wz = wzcen(i,j,k)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
         end do
      end if

      if ((fixulo_x.or.fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or.fixvhi_y .or. fixwhi_y)) then
         i = lo(1)
         j = hi(2)
         do k = lo(3),hi(3)
            ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            wz = wzcen(i,j,k)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
         end do
      end if

      if ((fixuhi_x.or.fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or.fixvhi_y .or. fixwhi_y)) then
         i = hi(1)
         j = hi(2)
         do k = lo(3),hi(3)
            ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            wz = vzcen(i,j,k)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
         end do
      end if

      if ((fixulo_x.or.fixvlo_x .or. fixwlo_x) .and. (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = lo(1)
         k = lo(3)
         do j = lo(2),hi(2)
            ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uy = uycen(i,j,k)
            vy = vycen(i,j,k)
            wy = wycen(i,j,k)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
         end do
      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = hi(1)
         k = lo(3)
         do j = lo(2),hi(2)
            ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uy = uycen(i,j,k)
            vy = vycen(i,j,k)
            wy = wycen(i,j,k)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
         end do
      end if


      if ((fixulo_x.or.fixvlo_x .or. fixwlo_x) .and. (fixuhi_z .or. fixvhi_z.or.fixwhi_z)) then
         i = lo(1)
         k = hi(3)
         do j = lo(2),hi(2)
            ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uy = uycen(i,j,k)
            vy = vycen(i,j,k)
            wy = wycen(i,j,k)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
         end do
      end if

      if ((fixuhi_x.or.fixvhi_x .or. fixwhi_x) .and. (fixuhi_z .or. fixvhi_z.or.fixwhi_z)) then
         i = hi(1)
         k = hi(3)
         do j = lo(2),hi(2)
            ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uy = uycen(i,j,k)
            vy = vycen(i,j,k)
            wy = wycen(i,j,k)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
         end do
      end if

      if ((fixulo_y .or.fixvlo_y.or.fixwlo_y) .and. (fixulo_z .or. fixvlo_z.or.fixwlo_z)) then
         j = lo(2)
         k = lo(3)
         do i = lo(1),hi(1)
            ux = uxcen(i,j,k)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
         end do
      end if

      if ((fixuhi_y .or.fixvhi_y.or. fixwhi_y) .and. (fixulo_z .or. fixvlo_z.or.fixvlo_z)) then
         j = hi(2)
         k = lo(3)
         do i = lo(1),hi(1)
            ux = vxcen(i,j,k)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
         end do
      end if

      if ((fixulo_y .or. fixvlo_y.or.fixwlo_y) .and. (fixuhi_z .or. fixvhi_z.or.fixwhi_z)) then
         j = lo(2)
         k = hi(3)
         do i = lo(1),hi(1)
            ux = vxcen(i,j,k)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
         end do
      end if

      if ((fixuhi_y .or. fixvhi_y.or.fixwhi_y) .and. (fixuhi_z .or. fixvhi_z.or.fixwhi_z)) then
         j = hi(2)
         k = hi(3)
         do i = lo(1),hi(1)
            ux = uxcen(i,j,k)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
         end do
      end if
!c
!c     Finally do all the corners
!c
      if ((fixulo_x.or.fixvlo_x .or. fixwlo_x) .and. (fixulo_y.or.fixvlo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z.or.fixwlo_z)) then
         i = lo(1)
         j = lo(2)
         k = lo(3)
         ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
      end if

      if ((fixuhi_x.or. fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixvlo_y.or.fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z.or.fixwlo_z)) then
         i = hi(1)
         j = lo(2)
         k = lo(3)
         ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
      end if

      if ((fixulo_x.or.fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or.fixvhi_y.or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z.or.fixwlo_z)) then
         i = lo(1)
         j = hi(2)
         k = lo(3)
         ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
      end if

      if ((fixuhi_x.or.fixvhi_x .or. fixwhi_x) .and. (fixuhi_y.or.fixvhi_y .or. fixwhi_y) .and.  &
          (fixulo_z .or. fixvlo_z.or.fixwlo_z)) then
         i = hi(1)
         j = hi(2)
         k = lo(3)
         ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
      end if

      if ((fixulo_x.or.fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or.fixvlo_y.or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z.or.fixwhi_z)) then
         i = lo(1)
         j = lo(2)
         k = hi(3)
         ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
      end if

      if ((fixuhi_z.or.fixvhi_x .or. fixwhi_x) .and. (fixulo_y.or.fixvlo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z.or.fixwhi_z)) then
         i = hi(1)
         j = lo(2)
         k = hi(3)
         ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
      end if

      if ((fixulo_x.or.fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixvhi_y.or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z.or.fixwhi_z)) then
         i = lo(1)
         j = hi(2)
         k = hi(3)
         ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
      end if

      if ((fixuhi_x.or.fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixvhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z.or. fixwhi_z)) then
         i = hi(1)
         j = hi(2)
         k = hi(3)
         ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
               xi1 = wy-vz
               xi2 = uz-wx
               xi3 = vx-uy
               xmag = max(sqrt(xi1**2+xi2**2+xi3**2), 1.d-12)
               xmag = one / xmag
               xi1 = xi1*xmag
               xi2 = xi2*xmag
               xi3 = xi3*xmag
               s11 = ux
               s12 = 0.5d0*(uy+vx)
               s13 = 0.5d0*(uz+wx)
               s22 = vy
               s23 = 0.5d0*(vz+wy)
               s33 = wz
               dmag(i,j,k,1) = xi1*(s11*xi1+s12*xi2+xi3*s13) &
                            + xi2*(s12*xi1+s22*xi2+xi3*s23) &
                            + xi3*(s13*xi1+s23*xi2+xi3*s33)
      end if

#     undef U
#     undef V      
#     undef W
#     undef ULOX
#     undef UHIX
#     undef ULOY
#     undef UHIY
#     undef ULOZ
#     undef UHIZ
#     undef VLOX
#     undef VHIX
#     undef VLOY
#     undef VHIY
#     undef VLOZ
#     undef VHIZ
#     undef WLOX
#     undef WHIX
#     undef WLOY
#     undef WHIY
#     undef WLOZ
#     undef WHIZ

      end subroutine derdmag

      subroutine dervortx (vort,DIMS(vort),nv,dat,DIMS(dat),ncomp, &
                                lo,hi,domlo,domhi,delta,xlo,time,dt, &
                                bc,level,grid_no) bind(C, name="dervortx")
      implicit none

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(vort)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     vort(DIMV(vort),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer   i,j,k
      REAL_T    vz, wy, dx, dy, dz
      REAL_T vzcen, vzlo, vzhi, wycen, wylo, wyhi

      logical   fixvlo_x, fixwlo_x, fixvhi_x, fixwhi_x
      logical   fixulo_y, fixwlo_y, fixuhi_y, fixwhi_y
      logical   fixulo_z, fixvlo_z, fixuhi_z, fixvhi_z
!c
!c     ::::: some useful macro definitions
!c
#     define U(i,j,k) dat(i,j,k,1)
#     define V(i,j,k) dat(i,j,k,2)
#     define W(i,j,k) dat(i,j,k,3)

#     define ULOY bc(2,1,1)
#     define UHIY bc(2,2,1)
#     define ULOZ bc(3,1,1)
#     define UHIZ bc(3,2,1)

#     define VLOX bc(1,1,2)
#     define VHIX bc(1,2,2)
#     define VLOZ bc(3,1,2)
#     define VHIZ bc(3,2,2)

#     define WLOX bc(1,1,3)
#     define WHIX bc(1,2,3)
#     define WLOY bc(2,1,3)
#     define WHIY bc(2,2,3)
!c
!c     ::::: statement functions that implement stencil
!c
      vzcen(i,j,k) = half*(V(i,j,k+1)-V(i,j,k-1))/dz
      vzlo(i,j,k)  = (V(i,j,k+1)+three*V(i,j,k)-four*V(i,j,k-1))/(three*dz)
      vzhi(i,j,k)  =-(V(i,j,k-1)+three*V(i,j,k)-four*V(i,j,k+1))/(three*dz)

      wycen(i,j,k) = half*(W(i,j+1,k)-W(i,j-1,k))/dy
      wylo(i,j,k)  = (W(i,j+1,k)+three*W(i,j,k)-four*W(i,j-1,k))/(three*dy)
      wyhi(i,j,k)  =-(W(i,j-1,k)+three*W(i,j,k)-four*W(i,j+1,k))/(three*dy)

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               vz = vzcen(i,j,k)
               wy = wycen(i,j,k)
               vort(i,j,k,1) = wy-vz
            end do
         end do
      end do

      fixvlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (VLOX .eq. EXT_DIR .or. VLOX .eq. HOEXTRAP) )
      fixvhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (VHIX .eq. EXT_DIR .or. VHIX .eq. HOEXTRAP) )
      fixwlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (WLOX .eq. EXT_DIR .or. WLOX .eq. HOEXTRAP) )
      fixwhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (WHIX .eq. EXT_DIR .or. WHIX .eq. HOEXTRAP) )

      fixulo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (ULOY .eq. EXT_DIR .or. ULOY .eq. HOEXTRAP) )
      fixuhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (UHIY .eq. EXT_DIR .or. UHIY .eq. HOEXTRAP) )
      fixwlo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (WLOY .eq. EXT_DIR .or. WLOY .eq. HOEXTRAP) )
      fixwhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (WHIY .eq. EXT_DIR .or. WHIY .eq. HOEXTRAP) )

      fixulo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (ULOZ .eq. EXT_DIR .or. ULOZ .eq. HOEXTRAP) )
      fixuhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (UHIZ .eq. EXT_DIR .or. UHIZ .eq. HOEXTRAP) )
      fixvlo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (VLOZ .eq. EXT_DIR .or. VLOZ .eq. HOEXTRAP) )
      fixvhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (VHIZ .eq. EXT_DIR .or. VHIZ .eq. HOEXTRAP) )
!c
!c     First do all the faces
!c
      if (fixvlo_x .or. fixwlo_x) then
         i = lo(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               wy = wycen(i,j,k)
               vz = vzcen(i,j,k)
               vort(i,j,k,1) = wy-vz
            end do
         end do
      end if

      if (fixvhi_x .or. fixwhi_x) then
         i = hi(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               wy = wycen(i,j,k)
               vz = vzcen(i,j,k)
               vort(i,j,k,1) = wy-vz
            end do
         end do
      end if

      if (fixulo_y .or. fixwlo_y) then
         j = lo(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
               vz = vzcen(i,j,k)
               vort(i,j,k,1) = wy-vz
            end do
         end do
      end if

      if (fixuhi_y .or. fixwhi_y) then
         j = hi(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
               vz = vzcen(i,j,k)
               vort(i,j,k,1) = wy-vz
            end do
         end do
      end if

      if (fixulo_z .or. fixvlo_z) then
         k = lo(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               wy = wycen(i,j,k)
               vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
               vort(i,j,k,1) = wy-vz
            end do
         end do
      end if

      if (fixuhi_z .or. fixvhi_z) then
         k = hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               wy = wycen(i,j,k)
               vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
               vort(i,j,k,1) = wy-vz
            end do
         end do
      end if
!c
!c     Next do all the edges
!c
      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y)) then
         i = lo(1)
         j = lo(2)
         do k = lo(3),hi(3)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            vz = vzcen(i,j,k)
            vort(i,j,k,1) = wy-vz
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y)) then
         i = hi(1)
         j = lo(2)
         do k = lo(3),hi(3)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            vz = vzcen(i,j,k)
            vort(i,j,k,1) = wy-vz
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y)) then
         i = lo(1)
         j = hi(2)
         do k = lo(3),hi(3)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            vz = vzcen(i,j,k)
            vort(i,j,k,1) = wy-vz
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y)) then
         i = hi(1)
         j = hi(2)
         do k = lo(3),hi(3)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            vz = vzcen(i,j,k)
            vort(i,j,k,1) = wy-vz
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         k = lo(3)
         do j = lo(2),hi(2)
            wy = wycen(i,j,k)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            vort(i,j,k,1) = wy-vz
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         k = lo(3)
         do j = lo(2),hi(2)
            wy = wycen(i,j,k)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            vort(i,j,k,1) = wy-vz
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         k = hi(3)
         do j = lo(2),hi(2)
            wy = wycen(i,j,k)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            vort(i,j,k,1) = wy-vz
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         k = hi(3)
         do j = lo(2),hi(2)
            wy = wycen(i,j,k)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            vort(i,j,k,1) = wy-vz
         end do
      end if

      if ((fixulo_y .or. fixwlo_y) .and. (fixulo_z .or. fixvlo_z)) then
         j = lo(2)
         k = lo(3)
         do i = lo(1),hi(1)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            vort(i,j,k,1) = wy-vz
         end do
      end if

      if ((fixuhi_y .or. fixwhi_y) .and. (fixulo_z .or. fixvlo_z)) then
         j = hi(2)
         k = lo(3)
         do i = lo(1),hi(1)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            vort(i,j,k,1) = wy-vz
         end do
      end if

      if ((fixulo_y .or. fixwlo_y) .and. (fixuhi_z .or. fixvhi_z)) then
         j = lo(2)
         k = hi(3)
         do i = lo(1),hi(1)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            vort(i,j,k,1) = wy-vz
         end do
      end if

      if ((fixuhi_y .or. fixwhi_y) .and. (fixuhi_z .or. fixvhi_z)) then
         j = hi(2)
         k = hi(3)
         do i = lo(1),hi(1)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            vort(i,j,k,1) = wy-vz
         end do
      end if
!c
!c     Finally do all the corners
!c
      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         j = lo(2)
         k = lo(3)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         vort(i,j,k,1) = wy-vz
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         j = lo(2)
         k = lo(3)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         vort(i,j,k,1) = wy-vz
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         j = hi(2)
         k = lo(3)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         vort(i,j,k,1) = wy-vz
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         j = hi(2)
         k = lo(3)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         vort(i,j,k,1) = wy-vz
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         j = lo(2)
         k = hi(3)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         vort(i,j,k,1) = wy-vz
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         j = lo(2)
         k = hi(3)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         vort(i,j,k,1) = wy-vz
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         j = hi(2)
         k = hi(3)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         vort(i,j,k,1) = wy-vz
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         j = hi(2)
         k = hi(3)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         vort(i,j,k,1) = wy-vz
      end if

#     undef U
#     undef V      
#     undef W
#     undef ULOY
#     undef UHIY
#     undef ULOZ
#     undef UHIZ
#     undef VLOX
#     undef VHIX
#     undef VLOZ
#     undef VHIZ
#     undef WLOX
#     undef WHIX
#     undef WLOY
#     undef WHIY

      end  subroutine dervortx


      subroutine dervorty (vort,DIMS(vort),nv,dat,DIMS(dat),ncomp, &
                                lo,hi,domlo,domhi,delta,xlo,time,dt, &
                                bc,level,grid_no)bind(C, name="dervorty")
      implicit none

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(vort)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     vort(DIMV(vort),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer   i,j,k
      REAL_T    uz, wx, dx, dy, dz
      REAL_T    uzcen, uzlo, uzhi, wxcen, wxlo, wxhi

      logical   fixvlo_x, fixwlo_x, fixvhi_x, fixwhi_x
      logical   fixulo_y, fixwlo_y, fixuhi_y, fixwhi_y
      logical   fixulo_z, fixvlo_z, fixuhi_z, fixvhi_z
!c
!c     ::::: some useful macro definitions
!c
#     define U(i,j,k) dat(i,j,k,1)
#     define V(i,j,k) dat(i,j,k,2)
#     define W(i,j,k) dat(i,j,k,3)

#     define ULOY bc(2,1,1)
#     define UHIY bc(2,2,1)
#     define ULOZ bc(3,1,1)
#     define UHIZ bc(3,2,1)

#     define VLOX bc(1,1,2)
#     define VHIX bc(1,2,2)
#     define VLOZ bc(3,1,2)
#     define VHIZ bc(3,2,2)

#     define WLOX bc(1,1,3)
#     define WHIX bc(1,2,3)
#     define WLOY bc(2,1,3)
#     define WHIY bc(2,2,3)
!c
!c     ::::: statement functions that implement stencil
!c
      uzcen(i,j,k) = half*(U(i,j,k+1)-U(i,j,k-1))/dz
      uzlo(i,j,k)  = (U(i,j,k+1)+three*U(i,j,k)-four*U(i,j,k-1))/(three*dz)
      uzhi(i,j,k)  =-(U(i,j,k-1)+three*U(i,j,k)-four*U(i,j,k+1))/(three*dz)

      wxcen(i,j,k) = half*(W(i+1,j,k)-W(i-1,j,k))/dx
      wxlo(i,j,k)  = (W(i+1,j,k)+three*W(i,j,k)-four*W(i-1,j,k))/(three*dx)
      wxhi(i,j,k)  =-(W(i-1,j,k)+three*W(i,j,k)-four*W(i+1,j,k))/(three*dx)

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               uz = uzcen(i,j,k)
               wx = wxcen(i,j,k)
               vort(i,j,k,1) = uz-wx
            end do
         end do
      end do

      fixvlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (VLOX .eq. EXT_DIR .or. VLOX .eq. HOEXTRAP) )
      fixvhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (VHIX .eq. EXT_DIR .or. VHIX .eq. HOEXTRAP) )
      fixwlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (WLOX .eq. EXT_DIR .or. WLOX .eq. HOEXTRAP) )
      fixwhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (WHIX .eq. EXT_DIR .or. WHIX .eq. HOEXTRAP) )

      fixulo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (ULOY .eq. EXT_DIR .or. ULOY .eq. HOEXTRAP) )
      fixuhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (UHIY .eq. EXT_DIR .or. UHIY .eq. HOEXTRAP) )
      fixwlo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (WLOY .eq. EXT_DIR .or. WLOY .eq. HOEXTRAP) )
      fixwhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (WHIY .eq. EXT_DIR .or. WHIY .eq. HOEXTRAP) )

      fixulo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (ULOZ .eq. EXT_DIR .or. ULOZ .eq. HOEXTRAP) )
      fixuhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (UHIZ .eq. EXT_DIR .or. UHIZ .eq. HOEXTRAP) )
      fixvlo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (VLOZ .eq. EXT_DIR .or. VLOZ .eq. HOEXTRAP) )
      fixvhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (VHIZ .eq. EXT_DIR .or. VHIZ .eq. HOEXTRAP) )
!c
!c     First do all the faces
!c
      if (fixvlo_x .or. fixwlo_x) then
         i = lo(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
               uz = uzcen(i,j,k)
               vort(i,j,k,1) = uz-wx
            end do
         end do
      end if

      if (fixvhi_x .or. fixwhi_x) then
         i = hi(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
               uz = uzcen(i,j,k)
               vort(i,j,k,1) = uz-wx
            end do
         end do
      end if

      if (fixulo_y .or. fixwlo_y) then
         j = lo(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               wx = wxcen(i,j,k)
               uz = uzcen(i,j,k)
               vort(i,j,k,1) = uz-wx
            end do
         end do
      end if

      if (fixuhi_y .or. fixwhi_y) then
         j = hi(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               wx = wxcen(i,j,k)
               uz = uzcen(i,j,k)
               vort(i,j,k,1) = uz-wx
            end do
         end do
      end if

      if (fixulo_z .or. fixvlo_z) then
         k = lo(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               wx = wxcen(i,j,k)
               uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
               vort(i,j,k,1) = uz-wx
            end do
         end do
      end if

      if (fixuhi_z .or. fixvhi_z) then
         k = hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               wx = wxcen(i,j,k)
               uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
               vort(i,j,k,1) = uz-wx
            end do
         end do
      end if
!c
!c     Next do all the edges
!c
      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y)) then
         i = lo(1)
         j = lo(2)
         do k = lo(3),hi(3)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uz = uzcen(i,j,k)
            vort(i,j,k,1) = uz-wx
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y)) then
         i = hi(1)
         j = lo(2)
         do k = lo(3),hi(3)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uz = uzcen(i,j,k)
            vort(i,j,k,1) = uz-wx
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y)) then
         i = lo(1)
         j = hi(2)
         do k = lo(3),hi(3)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uz = uzcen(i,j,k)
            vort(i,j,k,1) = uz-wx
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y)) then
         i = hi(1)
         j = hi(2)
         do k = lo(3),hi(3)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uz = uzcen(i,j,k)
            vort(i,j,k,1) = uz-wx
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         k = lo(3)
         do j = lo(2),hi(2)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vort(i,j,k,1) = uz-wx
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         k = lo(3)
         do j = lo(2),hi(2)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vort(i,j,k,1) = uz-wx
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         k = hi(3)
         do j = lo(2),hi(2)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vort(i,j,k,1) = uz-wx
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         k = hi(3)
         do j = lo(2),hi(2)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vort(i,j,k,1) = uz-wx
         end do
      end if

      if ((fixulo_y .or. fixwlo_y) .and. (fixulo_z .or. fixvlo_z)) then
         j = lo(2)
         k = lo(3)
         do i = lo(1),hi(1)
            wx = wxcen(i,j,k)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vort(i,j,k,1) = uz-wx
         end do
      end if

      if ((fixuhi_y .or. fixwhi_y) .and. (fixulo_z .or. fixvlo_z)) then
         j = hi(2)
         k = lo(3)
         do i = lo(1),hi(1)
            wx = wxcen(i,j,k)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vort(i,j,k,1) = uz-wx
         end do
      end if

      if ((fixulo_y .or. fixwlo_y) .and. (fixuhi_z .or. fixvhi_z)) then
         j = lo(2)
         k = hi(3)
         do i = lo(1),hi(1)
            wx = wxcen(i,j,k)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vort(i,j,k,1) = uz-wx
         end do
      end if

      if ((fixuhi_y .or. fixwhi_y) .and. (fixuhi_z .or. fixvhi_z)) then
         j = hi(2)
         k = hi(3)
         do i = lo(1),hi(1)
            wx = wxcen(i,j,k)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vort(i,j,k,1) = uz-wx
         end do
      end if
!c
!c     Finally do all the corners
!c
      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         j = lo(2)
         k = lo(3)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vort(i,j,k,1) = uz-wx
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         j = lo(2)
         k = lo(3)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vort(i,j,k,1) = uz-wx
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         j = hi(2)
         k = lo(3)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vort(i,j,k,1) = uz-wx
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         j = hi(2)
         k = lo(3)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vort(i,j,k,1) = uz-wx
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         j = lo(2)
         k = hi(3)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vort(i,j,k,1) = uz-wx
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         j = lo(2)
         k = hi(3)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vort(i,j,k,1) = uz-wx
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         j = hi(2)
         k = hi(3)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vort(i,j,k,1) = uz-wx
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         j = hi(2)
         k = hi(3)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vort(i,j,k,1) = uz-wx
      end if

#     undef U
#     undef V      
#     undef W
#     undef ULOY
#     undef UHIY
#     undef ULOZ
#     undef UHIZ
#     undef VLOX
#     undef VHIX
#     undef VLOZ
#     undef VHIZ
#     undef WLOX
#     undef WHIX
#     undef WLOY
#     undef WHIY

      end subroutine dervorty


      subroutine dervortz (vort,DIMS(vort),nv,dat,DIMS(dat),ncomp, &
                                lo,hi,domlo,domhi,delta,xlo,time,dt, &
                                bc,level,grid_no) bind(C, name="dervortz")
      implicit none

      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(vort)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     vort(DIMV(vort),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer   i,j,k
      REAL_T    uy, vx, dx, dy, dz
      REAL_T    uycen, uylo, uyhi, vxcen, vxlo, vxhi

      logical   fixvlo_x, fixwlo_x, fixvhi_x, fixwhi_x
      logical   fixulo_y, fixwlo_y, fixuhi_y, fixwhi_y
      logical   fixulo_z, fixvlo_z, fixuhi_z, fixvhi_z
!c
!c     ::::: some useful macro definitions
!c
#     define U(i,j,k) dat(i,j,k,1)
#     define V(i,j,k) dat(i,j,k,2)
#     define W(i,j,k) dat(i,j,k,3)

#     define ULOY bc(2,1,1)
#     define UHIY bc(2,2,1)
#     define ULOZ bc(3,1,1)
#     define UHIZ bc(3,2,1)

#     define VLOX bc(1,1,2)
#     define VHIX bc(1,2,2)
#     define VLOZ bc(3,1,2)
#     define VHIZ bc(3,2,2)

#     define WLOX bc(1,1,3)
#     define WHIX bc(1,2,3)
#     define WLOY bc(2,1,3)
#     define WHIY bc(2,2,3)
!c
!c     ::::: statement functions that implement stencil
!c
      uycen(i,j,k) = half*(U(i,j+1,k)-U(i,j-1,k))/dy
      uylo(i,j,k)  = (U(i,j+1,k)+three*U(i,j,k)-four*U(i,j-1,k))/(three*dy)
      uyhi(i,j,k)  =-(U(i,j-1,k)+three*U(i,j,k)-four*U(i,j+1,k))/(three*dy)

      vxcen(i,j,k) = half*(V(i+1,j,k)-V(i-1,j,k))/dx
      vxlo(i,j,k)  = (V(i+1,j,k)+three*V(i,j,k)-four*V(i-1,j,k))/(three*dx)
      vxhi(i,j,k)  =-(V(i-1,j,k)+three*V(i,j,k)-four*V(i+1,j,k))/(three*dx)

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               uy = uycen(i,j,k)
               vx = vxcen(i,j,k)
               vort(i,j,k,1) = vx-uy
            end do
         end do
      end do

      fixvlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (VLOX .eq. EXT_DIR .or. VLOX .eq. HOEXTRAP) )
      fixvhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (VHIX .eq. EXT_DIR .or. VHIX .eq. HOEXTRAP) )
      fixwlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (WLOX .eq. EXT_DIR .or. WLOX .eq. HOEXTRAP) )
      fixwhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (WHIX .eq. EXT_DIR .or. WHIX .eq. HOEXTRAP) )

      fixulo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (ULOY .eq. EXT_DIR .or. ULOY .eq. HOEXTRAP) )
      fixuhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (UHIY .eq. EXT_DIR .or. UHIY .eq. HOEXTRAP) )
      fixwlo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (WLOY .eq. EXT_DIR .or. WLOY .eq. HOEXTRAP) )
      fixwhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (WHIY .eq. EXT_DIR .or. WHIY .eq. HOEXTRAP) )

      fixulo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (ULOZ .eq. EXT_DIR .or. ULOZ .eq. HOEXTRAP) )
      fixuhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (UHIZ .eq. EXT_DIR .or. UHIZ .eq. HOEXTRAP) )
      fixvlo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (VLOZ .eq. EXT_DIR .or. VLOZ .eq. HOEXTRAP) )
      fixvhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (VHIZ .eq. EXT_DIR .or. VHIZ .eq. HOEXTRAP) )
!c
!c     First do all the faces
!c
      if (fixvlo_x .or. fixwlo_x) then
         i = lo(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
               uy = uycen(i,j,k)
               vort(i,j,k,1) = vx-uy
            end do
         end do
      end if

      if (fixvhi_x .or. fixwhi_x) then
         i = hi(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
               uy = uycen(i,j,k)
               vort(i,j,k,1) = vx-uy
            end do
         end do
      end if

      if (fixulo_y .or. fixwlo_y) then
         j = lo(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
               vort(i,j,k,1) = vx-uy
            end do
         end do
      end if

      if (fixuhi_y .or. fixwhi_y) then
         j = hi(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
               vort(i,j,k,1) = vx-uy
            end do
         end do
      end if

      if (fixulo_z .or. fixvlo_z) then
         k = lo(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               uy = uycen(i,j,k)
               vort(i,j,k,1) = vx-uy
            end do
         end do
      end if

      if (fixuhi_z .or. fixvhi_z) then
         k = hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               uy = uycen(i,j,k)
               vort(i,j,k,1) = vx-uy
            end do
         end do
      end if
!c
!c     Next do all the edges
!c
      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y)) then
         i = lo(1)
         j = lo(2)
         do k = lo(3),hi(3)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vort(i,j,k,1) = vx-uy
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y)) then
         i = hi(1)
         j = lo(2)
         do k = lo(3),hi(3)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vort(i,j,k,1) = vx-uy
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y)) then
         i = lo(1)
         j = hi(2)
         do k = lo(3),hi(3)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vort(i,j,k,1) = vx-uy
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y)) then
         i = hi(1)
         j = hi(2)
         do k = lo(3),hi(3)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vort(i,j,k,1) = vx-uy
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         k = lo(3)
         do j = lo(2),hi(2)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            uy = uycen(i,j,k)
            vort(i,j,k,1) = vx-uy
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         k = lo(3)
         do j = lo(2),hi(2)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            uy = uycen(i,j,k)
            vort(i,j,k,1) = vx-uy
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         k = hi(3)
         do j = lo(2),hi(2)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            uy = uycen(i,j,k)
            vort(i,j,k,1) = vx-uy
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         k = hi(3)
         do j = lo(2),hi(2)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            uy = uycen(i,j,k)
            vort(i,j,k,1) = vx-uy
         end do
      end if

      if ((fixulo_y .or. fixwlo_y) .and. (fixulo_z .or. fixvlo_z)) then
         j = lo(2)
         k = lo(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vort(i,j,k,1) = vx-uy
         end do
      end if

      if ((fixuhi_y .or. fixwhi_y) .and. (fixulo_z .or. fixvlo_z)) then
         j = hi(2)
         k = lo(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vort(i,j,k,1) = vx-uy
         end do
      end if

      if ((fixulo_y .or. fixwlo_y) .and. (fixuhi_z .or. fixvhi_z)) then
         j = lo(2)
         k = hi(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vort(i,j,k,1) = vx-uy
         end do
      end if

      if ((fixuhi_y .or. fixwhi_y) .and. (fixuhi_z .or. fixvhi_z)) then
         j = hi(2)
         k = hi(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vort(i,j,k,1) = vx-uy
         end do
      end if
!c
!c     Finally do all the corners
!c
      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         j = lo(2)
         k = lo(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vort(i,j,k,1) = vx-uy
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         j = lo(2)
         k = lo(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vort(i,j,k,1) = vx-uy
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         j = hi(2)
         k = lo(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vort(i,j,k,1) = vx-uy
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         j = hi(2)
         k = lo(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vort(i,j,k,1) = vx-uy
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         j = lo(2)
         k = hi(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vort(i,j,k,1) = vx-uy
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         j = lo(2)
         k = hi(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vort(i,j,k,1) = vx-uy
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         j = hi(2)
         k = hi(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vort(i,j,k,1) = vx-uy
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         j = hi(2)
         k = hi(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vort(i,j,k,1) = vx-uy
      end if

#     undef U
#     undef V      
#     undef W
#     undef ULOY
#     undef UHIY
#     undef ULOZ
#     undef UHIZ
#     undef VLOX
#     undef VHIX
#     undef VLOZ
#     undef VHIZ
#     undef WLOX
#     undef WHIX
#     undef WLOY
#     undef WHIY

      end subroutine dervortz

!c=========================================================

      subroutine dergrdp (grdp,DIMS(gp),nv,p,DIMS(p),ncomp, &
                              lo,hi,domlo,domhi,dx,xlo,time,dt, &
                              bc,level,grid_no) bind(C, name="dergrdp")
      implicit none
!c
!c     This routine computes the magnitude of pressure gradient 
!c
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(gp)
      integer DIMDEC(p)
      integer domlo(SDIM), domhi(SDIM)
      integer nv, ncomp
      integer bc(SDIM,2,ncomp)
      REAL_T  dx(SDIM), xlo(SDIM), time, dt
      REAL_T  grdp(DIMV(gp),nv)
      REAL_T  p(DIMV(p),ncomp)
      integer level, grid_no

      REAL_T     gpx, gpy, gpz, idx, idy, idz
      integer    i,j,k

      idx = 1.0d0 / dx(1)
      idy = 1.0d0 / dx(2)
      idz = 1.0d0 / dx(3)

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            gpx = fourth * (p(i+1,j,k  ,1)-p(i,j,k  ,1)+p(i+1,j+1,k  ,1)-p(i,j+1,k  ,1)+ &
                           p(i+1,j,k+1,1)-p(i,j,k+1,1)+p(i+1,j+1,k+1,1)-p(i,j+1,k+1,1))*idx
            gpy = fourth * (p(i,j+1,k  ,1)-p(i,j,k  ,1)+p(i+1,j+1,k  ,1)-p(i+1,j,k  ,1)+ &
                           p(i,j+1,k+1,1)-p(i,j,k+1,1)+p(i+1,j+1,k+1,1)-p(i+1,j,k+1,1))*idy
            gpz = fourth * (p(i,  j,k+1,1)-p(i,  j,k,1)+p(i,  j+1,k+1,1)-p(i,  j+1,k,1)+ &
                           p(i+1,j,k+1,1)-p(i+1,j,k,1)+p(i+1,j+1,k+1,1)-p(i+1,j+1,k,1))*idz
            grdp(i,j,k,1) = sqrt(gpx**2 + gpy**2 + gpz**2)
          end do
        end do
      end do

      end subroutine dergrdp

!c=========================================================

      subroutine dernull (e,DIMS(e),nv,dat,DIMS(dat),ncomp, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt,bc, &
                              level,grid_no) bind(C, name="dernull")
      implicit none
      !
      ! This is a null derived routine.
      !
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      end subroutine dernull

end module derive_3d_module
