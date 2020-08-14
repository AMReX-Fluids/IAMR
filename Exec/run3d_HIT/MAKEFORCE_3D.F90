#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PROB_NS_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 3

module MakeForce_3D_module

  implicit none

  private

  public :: FORT_MAKEFORCE

contains

! ::: -----------------------------------------------------------
!
!     This is a dummy routine to add the forcing terms to the momentum equation
!     It should be overwritten in the local run folder.  
!
   subroutine FORT_MAKEFORCE( time, &
                              force, f_lo, f_hi,&
                              vel, v_lo, v_hi,&
                              scal, s_lo, s_hi,&
                              dx, xlo, xhi, gravity, scomp, ncomp, &
                              nscal, getForceVerbose ) &
                              bind(C, name="FORT_MAKEFORCE")

      implicit none

! In/Out
      integer :: f_lo(3), f_hi(3)
      integer :: v_lo(3), v_hi(3)
      integer :: s_lo(3), s_hi(3)
      integer :: scomp, ncomp
      integer :: nscal, getForceVerbose
      REAL_T  :: time, dx(3)
      REAL_T  :: xlo(3), xhi(3)
      REAL_T  :: gravity
      REAL_T, dimension(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),scomp:scomp+ncomp-1) :: force
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),0:SDIM-1) :: vel
      REAL_T, dimension(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),0:nscal-1) :: scal

#include <probdata.H>

#ifdef DO_IAMR_FORCE
#include <forcedata.H>
#endif

!c
!c     ::::: local variables
!c
      integer i, j, k, n
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      integer a2, a3, a4, a5
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  sga, cga
      REAL_T  f1, f2, f3
      REAL_T  twicePi
      REAL_T  infl_time
      REAL_T  kxd, kyd, kzd
      REAL_T  xt, yt, zt
      REAL_T  zlo
      integer kx, ky, kz, mode_count, xstep, ystep, zstep, count
      integer isioproc, do_trac2
      integer nXvel, nYvel, nZvel, nRho, nTrac, nTrac2, nRhoScal, nTracScal, nTrac2Scal

      REAL_T  velmin(0:SDIM-1)
      REAL_T  velmax(0:SDIM-1)
      REAL_T  scalmin(0:nscal-1)
      REAL_T  scalmax(0:nscal-1)
      REAL_T  forcemin(scomp:scomp+ncomp-1)
      REAL_T  forcemax(scomp:scomp+ncomp-1)
      REAL_T  tmpMin, tmpMax

      call bl_ns_dotrac2(do_trac2)

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      ilo = force_l1
      jlo = force_l2
      klo = force_l3
      ihi = force_h1
      jhi = force_h2
      khi = force_h3

!c     Assumes components are in the following order
      nXvel  = 0
      nYvel  = 1
      nZvel  = 2
      nRho   = 3
      nTrac  = 4
      nTrac2 = 5

      nRhoScal   = nRho-SDIM
      nTracScal  = nTrac-SDIM
      nTrac2Scal = nTrac2-SDIM

      if (getForceVerbose.gt.0) then
         call bl_pd_is_ioproc(isioproc)
         if (isioproc.eq.1) then

            write (6,*) "In MAKEFORCE"
            
            write (6,*) "probtype = ",probtype
            write (6,*) "gravity = ",gravity
            write (6,*) "scomp = ",scomp
            write (6,*) "ncomp = ",ncomp
            write (6,*) "nscal = ",nscal
            
            do n = 0, SDIM-1
               velmin(n) = 1.d234
               velmax(n) = -1.d234
            enddo
            do n = 0, nscal-1
               scalmin(n) = 1.d234
               scalmax(n) = -1.d234
            enddo

            count = 0
!c     Get min/max
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
!c     Velocities
                     do n = 0, SDIM-1
                        if (vel(i,j,k,n).gt.velmax(n)) then
                           velmax(n)=vel(i,j,k,n)
                        endif
                        if (vel(i,j,k,n).lt.velmin(n)) then
                           velmin(n)=vel(i,j,k,n)
                        endif
                     enddo
!c     Scalars
                     if (scal(i,j,k,0).lt.0.001) then
                        count=count+1
!c                        write (*,*) i,j,k,scal(i,j,k,n)
                     endif
                     do n = 0, nscal-1
                        if (scal(i,j,k,n).gt.scalmax(n)) then
                           scalmax(n)=scal(i,j,k,n)
                        endif
                        if (scal(i,j,k,n).lt.scalmin(n)) then
                           scalmin(n)=scal(i,j,k,n)
                        endif
                     enddo
                     
                  enddo
               enddo
            enddo

            write(*,*) "DODGY DENSITY COUNT = ",count

            do n = 0, SDIM-1
               write (6,*) "velmin (",n,") = ",velmin(n)
               write (6,*) "velmax (",n,") = ",velmax(n)
            enddo
            do n = 0, nscal-1
               write (6,*) "scalmin(",n,") = ",scalmin(n)
               write (6,*) "scalmax(",n,") = ",scalmax(n)
            enddo
         endif
      endif

!c
!c     Here's where the forcing actually gets done
!c
      
      if (scomp.eq.0) then
!c
!c     Do velocity forcing
!c
         if (probtype.eq.20) then 
            do k = klo, khi
               z = xlo(3) + hz*(float(k-klo) + half)
               do j = jlo, jhi
                  y = xlo(2) + hy*(float(j-jlo) + half)
                  do i = ilo, ihi
                     x = xlo(1) + hx*(float(i-ilo) + half)
                     force(i,j,k,nXvel) = zero
                     force(i,j,k,nYvel) = zero
                     force(i,j,k,nZvel) = zero
                     if (do_trac2.eq.1) then
                        force(i,j,k,nZvel) = force(i,j,k,nZvel) + thermal_expansion*scal(i,j,k,nTrac2Scal)
                     endif
                  enddo
               enddo
            enddo
         else if (probtype.eq.18) then 
!c     Round jet/plume
            do k = klo, khi
               z = xlo(3) + hz*(float(k-klo) + half)
               do j = jlo, jhi
                  y = xlo(2) + hy*(float(j-jlo) + half)
                  do i = ilo, ihi
                     x = xlo(1) + hx*(float(i-ilo) + half)
                     force(i,j,k,nXvel) = zero
                     force(i,j,k,nYvel) = zero
                     force(i,j,k,nZvel) = gravity*scal(i,j,k,nRhoScal)
                     if (do_trac2.eq.1) then
                        force(i,j,k,nZvel) = force(i,j,k,nZvel) + thermal_expansion*scal(i,j,k,nTrac2Scal)
                     endif
                     if (do_jet_sponge.eq.1) then
!c                        call bl_pd_is_ioproc(isioproc)
!c                        if (isioproc.eq.1) then
!c                           write (*,*) "jet_sponge_scale = ",jet_sponge_scale
!c                        endif
                        if (z.gt.jet_sponge_height) then
                           force(i,j,k,nXvel) = force(i,j,k,nXvel) - jet_sponge_scale*vel(i,j,k,0)*scal(i,j,k,nRhoScal)
                           force(i,j,k,nYvel) = force(i,j,k,nYvel) - jet_sponge_scale*vel(i,j,k,1)*scal(i,j,k,nRhoScal)
                           force(i,j,k,nZvel) = force(i,j,k,nZvel) - jet_sponge_scale*vel(i,j,k,2)*scal(i,j,k,nRhoScal)
                        else if (sqrt(x*x+y*y).gt.jet_sponge_radius) then
                           force(i,j,k,nZvel) = force(i,j,k,nZvel) - jet_sponge_scale*vel(i,j,k,2)*scal(i,j,k,nRhoScal)
                        endif
                     endif
                  enddo
               enddo
            enddo
         else if (probtype.eq.14.or.probtype.eq.15) then
#ifdef DO_IAMR_FORCE
!c     Homogeneous Isotropic Turbulence
            twicePi=two*Pi
            
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
            kappaMax = dfloat(nmodes)/Lmin + 1.0d-8
            nxmodes = nmodes*int(0.5+Lx/Lmin)
            nymodes = nmodes*int(0.5+Ly/Lmin)
            nzmodes = nmodes*int(0.5+Lz/Lmin)

            xstep = int(Lx/Lmin+0.5)
            ystep = int(Ly/Lmin+0.5)
            zstep = int(Lz/Lmin+0.5)

            if (forcing_twice_wavelength.eq.1) then
               HLx = Lx/two
               HLy = Ly/two
               HLz = Lz/two
            else
               HLx = Lx
               HLy = Ly
               HLz = Lz
            endif

!!$omp parallel do private(i,j,k,x,y,z,f1,f2,f3,kx,ky,kz,kxd,kyd,kzd)
!!$omp&private(xT,kappa)
            do k = klo, khi
               z = zlo + hz*(float(k-klo) + half)
               do j = jlo, jhi
                  y = xlo(2) + hy*(float(j-jlo) + half)
                  do i = ilo, ihi
                     x = xlo(1) + hx*(float(i-ilo) + half)
                     f1 = zero
                     f2 = zero
                     f3 = zero
                     do kz = mode_start*zstep, nzmodes, zstep
                        kzd = dfloat(kz)
                        do ky = mode_start*ystep, nymodes, ystep
                           kyd = dfloat(ky)
                           do kx = mode_start*xstep, nxmodes, xstep
                              kxd = dfloat(kx)
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
                        kzd = dfloat(kz)
                        do ky = mode_start, nymodes
                           kyd = dfloat(ky)
                           do kx = mode_start, nxmodes
                              kxd = dfloat(kx)
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
                                    f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/Lx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/Ly+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/Lz+FPZ(kx,ky,kz))
                                    f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/Lx+FPX(kx,ky,kz)) * cos(twicePi*kyd*y/Ly+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/Lz+FPZ(kx,ky,kz))
                                    f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/Lx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/Ly+FPY(kx,ky,kz)) * cos(twicePi*kzd*z/Lz+FPZ(kx,ky,kz))
                                 endif
                              endif
                           enddo
                        enddo
                     enddo
                     if (use_rho_in_forcing.eq.1) then
                        force(i,j,k,nXvel) = f1*scal(i,j,k,nRhoScal)
                        force(i,j,k,nYvel) = f2*scal(i,j,k,nRhoScal)
                        force(i,j,k,nZvel) = f3*scal(i,j,k,nRhoScal)
                     else
                        force(i,j,k,nXvel) = f1
                        force(i,j,k,nYvel) = f2
                        force(i,j,k,nZvel) = f3
                     endif
                  enddo
               enddo
            enddo
#else
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,k,nXvel) = zero
                     force(i,j,k,nYvel) = zero
                     force(i,j,k,nXvel) = zero
                  enddo
               enddo
            enddo
#endif
         else if (probtype.eq.19) then
!c     Coriolis
            do k = klo, khi
               z = xlo(3) + hz*(float(k-klo) + half)
               do j = jlo, jhi
                  y = xlo(2) + hy*(float(j-jlo) + half)
                  do i = ilo, ihi
                     x = xlo(1) + hx*(float(i-ilo) + half)
                     force(i,j,k,nXvel) = scal(i,j,k,nRhoScal)*( two*omega*vel(i,j,k,nYvel)+omega*omega*x)
                     force(i,j,k,nYvel) = scal(i,j,k,nRhoScal)*(-two*omega*vel(i,j,k,nXvel)+omega*omega*y)
                     force(i,j,k,nZvel) = gravity*scal(i,j,k,nRhoScal)
                     if (do_trac2.eq.1) then
                        force(i,j,k,nZvel) = force(i,j,k,nZvel) + thermal_expansion*scal(i,j,k,nTrac2Scal)
                     endif
                  enddo
               enddo
            enddo
         else if (probtype.eq.99.and.abs(grav_angle).gt.0.001) then
!c     Angled gravity
            sga =  gravity * sin(Pi*grav_angle/180.)
            cga = -gravity * cos(Pi*grav_angle/180.)
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,k,nXvel) = scal(i,j,k,nRhoScal)*sga
                     force(i,j,k,nYvel) = zero
                     force(i,j,k,nZvel) = scal(i,j,k,nRhoScal)*cga
                  enddo
               enddo
            enddo
!c     Default to gravity...
         elseif (abs(gravity).gt.0.0001) then
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,k,nXvel) = zero
                     force(i,j,k,nYvel) = zero
                     force(i,j,k,nZvel) = gravity*scal(i,j,k,nRhoScal)
                  enddo
               enddo
            enddo
!c     else to zero
         else
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,k,nXvel) = zero
                     force(i,j,k,nYvel) = zero
                     force(i,j,k,nZvel) = zero
                  enddo
               enddo
            enddo
         endif
!c     End of velocity forcing
      endif

      if ((scomp+ncomp).gt.BL_SPACEDIM) then
!c
!c     Scalar forcing
!c
         do n = max(scomp,nRho), scomp+ncomp-1
            if (n.eq.nRho) then
!c
!c     Density
!c
               do k = klo, khi
                  do j = jlo, jhi
                     do i = ilo, ihi
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            else if (n.eq.nTrac) then
!c
!c     Tracer
!c
               do k = klo, khi
                  do j = jlo, jhi
                     do i = ilo, ihi
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            else if (n.eq.nTrac2.and.do_trac2.eq.1) then
!c
!c     Other scalar
!c
               if (probtype.eq.20) then 
!c     Temperature perturbation
                  do k = klo, khi
                     do j = jlo, jhi
                        do i = ilo, ihi
                           force(i,j,k,n) = heating_coeff * scal(i,j,k,nTracScal)
                        enddo
                     enddo
                  enddo
               else  if (probtype.eq.18) then 
!c     Round Jet/Plume (18)
                  do k = klo, khi
                     z = xlo(3) + hz*(float(k-klo) + half)
                     do j = jlo, jhi
                        do i = ilo, ihi
                           if (abs(z-heating_centre).lt.heating_radius) then 
                              force(i,j,k,n) = heating_coeff * scal(i,j,k,nTracScal)
                           else
                              force(i,j,k,n) = zero
                           endif
                        enddo
                     enddo
                  enddo
               else  if (probtype.eq.19) then 
!c     Coriolis evaporation (19)
                  do k = klo, khi
                     z = xlo(3) + hz*(float(k-klo) + half)
                     do j = jlo, jhi
                        do i = ilo, ihi
                           if (abs(z-heating_centre).lt.heating_radius) then 
                              force(i,j,k,n) = heating_coeff
                           else
                              force(i,j,k,n) = zero
                           endif
                        enddo
                     enddo
                  enddo
               else
!c     Some other probtype
                  do k = klo, khi
                     do j = jlo, jhi
                        do i = ilo, ihi
                           force(i,j,k,n) = zero
                        enddo
                     enddo
                  enddo
               endif
            else
!c
!c     Other scalar
!c
               do k = klo, khi
                  do j = jlo, jhi
                     do i = ilo, ihi
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            endif
         enddo
      endif
      
      if (getForceVerbose.gt.0 .and. isioproc.eq.1) then
         do n = scomp,scomp+ncomp-1
            forcemin(n) = 1.d234
            forcemax(n) = -1.d234
         enddo
         do k = klo, khi
            do j = jlo, jhi
               do i = ilo, ihi
                  do n = scomp,ncomp+scomp-1
                     forcemin(n) = min(forcemin(n),force(i,j,k,n))
                     forcemax(n) = max(forcemax(n),force(i,j,k,n))
                  enddo
               enddo
            enddo
         enddo
         do n = scomp,ncomp+scomp-1
            write (6,*) "forcemin (",n,") = ",forcemin(n)
            write (6,*) "forcemax (",n,") = ",forcemax(n)
         enddo
      endif

   end subroutine FORT_MAKEFORCE

end module MakeForce_3D_module
