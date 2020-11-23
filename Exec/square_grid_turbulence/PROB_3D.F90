
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PROB_NS_F.H>
#include <AMReX_ArrayLim.H>

#ifdef BL_DO_FLCT
#include <infl_frc.H>
#endif

#define SDIM 3


      block data rt_common
#include <probdata.H>
      data rt_pertamp/0.0D0/
      data rt_nfreq/0/
      data rt_xfrontw/0.0d0/
      end

module prob_3D_module

  implicit none

  private

  public :: amrex_probinit, FORT_INITDATA,  &
            initflowpastcylinder, &
            FORT_DSDTFILL, &
            FORT_ADVERROR, FORT_ADV2ERROR, FORT_TEMPERROR, FORT_MVERROR, &
            FORT_DENFILL, FORT_ADVFILL, FORT_TEMPFILL, FORT_XVELFILL, &
            FORT_YVELFILL, FORT_ZVELFILL, FORT_PRESFILL, FORT_DIVUFILL

  !
  ! Set problo and probhi as module variables so that they can be
  ! used
  !
  REAL_T :: m_problo(SDIM), m_probhi(SDIM)

contains


!c ::: -----------------------------------------------------------
!c ::: This routine is called at problem initialization time
!c ::: and when restarting from a checkpoint file.
!c ::: The purpose is (1) to specify the initial time value
!c ::: (not all problems start at time=0.0) and (2) to read
!c ::: problem specific data from a namelist or other inputcdm
!c ::: files and possibly store them or derived information
!c ::: in FORTRAN common blocks for later use.
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: init      => TRUE if called at start of problem run
!c :::              FALSE if called from restart
!c ::: name      => name of "probin" file
!c ::: namlen    => length of name
!c ::: strttime <=  start problem with this time variable
!c :::
!c ::: -----------------------------------------------------------

      subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)
      implicit none
      integer init,namlen
      integer name(namlen)
      integer untin, i
      REAL_T  problo(SDIM), probhi(SDIM)


#include <probdata.H>

#ifdef DO_IAMR_FORCE
#include <forcedata.H>
#endif

#ifdef BL_DO_FLCT
#include <INFL_FORCE_F.H>
#endif

      INTEGER dimFile(3)
      integer nCompFile
      REAL_T dxFile(3)

      REAL_T  twicePi, kxd, kyd, kzd
      REAL_T  thetaTmp, phiTmp
      REAL_T  cosThetaTmp, cosPhiTmp
      REAL_T  sinThetaTmp, sinPhiTmp
      REAL_T  px, py, pz, mp2, Ekh
      integer kx, ky, kz, mode_count, reduced_mode_count
      integer modx, mody, modz, xstep, ystep, zstep
      integer iseed

      REAL_T  Lx, Ly, Lz, Lmin
      REAL_T  kappa, kappaMax, freqMin, freqMax, freqDiff, pdk

      namelist /fortin/ denerr, vorterr, adverr, temperr, &
     			denfact, xblob, yblob, zblob, radblob,  &
                       velfact, probtype, randfact, bubgrad, &
     			rhozero, rhograd, tempzero, c_d, r_d,  &
                       adv_dir, adv_vel, slot_vel, axis_dir, radvort, &
                       den1,den2,vel1,vel2,delta0,xlev1,zlev1,amag, &
                       vb_unifdir, blrandseed, turb_scale, injection_time, &
                       override_turb_force
#ifdef BL_DO_FLCT
      namelist /fortin/ forceInflow, numInflPlanesStore, strmwse_dir, &
                       forceLo, forceHi, flct_file, nCompInflow, infl_type, &
                       tstart_turb
#endif
      namelist /fortin/ rt_splitx, rt_xfrontw, rt_den_1, rt_den_2, &
                       rt_pertamp, rt_nfreq, rt_graddenerr

      namelist /fortin/ eul_nfreq, eul_pertamp, iseed

      namelist /fortin/ Vco, Rfu, Rtran, tVco_l, tVco_r, Vco_l, Vco_r

      namelist /fortin/ grav_angle, omega, infl_time_offset, ref_centre, ref_radius, time_offset, &
                       thermal_expansion, heating_coeff, heating_centre, heating_radius

      namelist /fortin/ density_pert, interface_height, wavelength_min, wavelength_max, tracer_height

#ifdef DO_IAMR_FORCE
      namelist /fortin/ nmodes, nxmodes, nymodes, nzmodes, mode_start, hack_lz, &
                      forcing_type, spectrum_type, ref_type, forcing_twice_wavelength, &
                       forcing_xlength, forcing_ylength, forcing_zlength, &
                       forcing_time_scale_min, forcing_time_scale_max,  &
                       force_scale, forcing_epsilon, &
                       use_rho_in_forcing, do_mode_division, div_free_force, moderate_zero_modes, &
                       AXY, BXY, CXY, DXY, PXY, QXY, RXY, &
                       AZX, BZX, CZX, DZX, PZX, QZX, RZX, &
                       AYZ, BYZ, CYZ, DYZ, PYZ, QYZ, RYZ, &
                       FTX, FTY, FTZ, TAT, TAP, &
                       FPX, FPY, FPZ, FAX, FAY, FAZ
#endif

      namelist /fortin/ jet_x, jet_y, jet_width, jet_rho, jet_vel, &
                       coflow_rho, coflow_vel, jet_temp, ref_height, ref_height2, plane_jet, &
                       do_jet_sponge, jet_sponge_scale, jet_sponge_height, jet_sponge_radius

      namelist /fortin/ holeRad, holeBLfac, nHolesX, nHolesY, nHolesZ, holeSp, slotWidth, alpha, beta

      namelist /fortin/ tInflowFact_l, tInflowFact_r, InflowFact_l, InflowFact_r
      namelist /fortin/ do_inlet_ref, inlet_ref_height
      namelist /fortin/ lid_vel
!c
!c      Build "probin" filename -- the name of file containing fortin namelist.
!c
      integer maxlen, isioproc, ierr
      parameter (maxlen=256)

      REAL_T frecon, onep7
      parameter (frecon=.219, onep7=1.7)

      character probin*(maxlen)

      integer j, k, m, n
      integer MAXPHASE
      parameter (MAXPHASE = 1000)
      DOUBLE PRECISION rn, permin, permax, xtmp, ytmp, pert

      REAL_T   f0, fmin, fmax, fdif
      integer  idum

      call bl_pd_is_ioproc(isioproc)

      if (namlen .gt. maxlen) call bl_error('probin file name too long')

      do i = 1, namlen
         probin(i:i) = char(name(i))
      end do

      !
      ! Set module variables
      !
      m_problo = problo
      m_probhi = probhi


#ifdef BL_DO_FLCT
      forceInflow = .FALSE.
      numInflPlanesStore = 8
      forceLo = .TRUE.
      forceHi = .FALSE.
      strmwse_dir = FLCT_ZVEL
      flct_file = ""
      turb_scale = 1
      infl_type = infl_periodic_type
#endif

      tVco_l = -1
      tVco_r = -1
      Vco_l = -1
      Vco_r = -1
      eul_nfreq = 2
      iseed = 111397

      holeRad = .0015d0
      holeSp = .006875d0
      slotWidth = 0.006
      nHolesX = 5
      nHolesY = 5
      nHolesZ = 1
      holeBLfac = 0.5d0
      alpha = 0.d0
      beta = 0.d0

      tInflowFact_l = -1.d0
      tInflowFact_r = -1.d0
      InflowFact_l = 1.d0
      InflowFact_r = 1.d0

      do_inlet_ref = 0
      inlet_ref_height = -1.d0

      untin = 9
      if (namlen .eq. 0) then
         open(untin,file='probin',form='formatted',status='old')
      else
         open(untin,file=probin(1:namlen),form='formatted',status='old')
      end if
      read(untin,fortin)
!c      if (isioproc .eq. 1) write(6,fortin)
      close(unit=untin)

#ifdef BL_DO_FLCT
      if (forceInflow .eqv. .FALSE.) then
         forceLo = .FALSE.
         forceHi = .FALSE.
      else
         if (flct_file .ne. "") then
            ierr=0
            untin=20
            if (isioproc .eq. 1) print*, 'Initializing turbulence ...'
            open(untin, file=trim(flct_file)//'/HDR', form='formatted', iostat=ierr)
            if (ierr .ne. 0) then
               call bl_abort('Problem opening file: ' // trim(flct_file) // '/HDR')
            end if
            call RD_SCL_FLCTHD(untin,nCompFile, dimFile, probSizeFile, dxFile)
            close(untin)
         endif
         if (strmwse_dir .NE. 3) then
            call bl_error('ERROR: turbulent inflow needs strmwse_dir=3')
         endif
         if (isioproc .eq. 1) then
            print *, 'dimFile: ',      (dimFile(i),i=1,3)
            print *, 'probSizeFile: ', (probSizeFile(i),i=1,3)
            print *, 'dxFile: ',       (dxFile(i),i=1,3)
         end if
      endif

!c     Need to set adv_vel for the round jet
      if (probtype.eq.18) then
         adv_vel = jet_vel
      endif
#endif

      domnlo(1) = problo(1)
      domnlo(2) = problo(2)
      domnlo(3) = problo(3)
      domnhi(1) = probhi(1)
      domnhi(2) = probhi(2)
      domnhi(3) = probhi(3)

      if (probtype.eq.8) then
         freq(1) = frecon
         do i=2,10
            freq(i)=freq(1)/float(i)
         enddo
         mag(1) = amag
         mag(2) = .75*amag
         mag(3) = .55*amag
         mag(4) = .44*amag
         do i=5,10
            mag(i)=onep7*mag(1)/float(i)
         enddo
      endif
!c
!c     Initialize the common blocks
!c
      do i=1, SDIM
         f_problo(i) = problo(i)
         f_probhi(i) = probhi(i)
      enddo

      if ( probtype .eq. 21 ) then

         call blutilinitrand(iseed)
         do i = 1, eul_nfreq
            do j = 1, eul_nfreq
               do k = 1, eul_nfreq
                 call blutilrand(rn)
                 eul_ranampl(i,j,k,1) = 2.d0*(rn-0.5d0)
                 call blutilrand(rn)
                 eul_ranampl(i,j,k,2) = 2.d0*(rn-0.5d0)
                 call blutilrand(rn)
                 eul_ranampl(i,j,k,3) = 2.d0*(rn-0.5d0)
                 call blutilrand(rn)
                 eul_ranphase(i,j,k,1) = 2.d0*Pi*rn
                 call blutilrand(rn)
                 eul_ranphase(i,j,k,2) = 2.d0*Pi*rn
                 call blutilrand(rn)
                 eul_ranphase(i,j,k,3) = 2.d0*Pi*rn
              end do
            end do
         end do

         eul_ranampl(1,1,1,1) = 0.d0
         eul_ranampl(1,1,1,2) = 0.d0
         eul_ranampl(1,1,1,3) = 0.d0

      endif

      if ( probtype .eq. 10 ) then
         if ( rt_max_freq .lt. rt_nfreq ) then
            stop 'RT_INIT broken: 1'
         end if
         call blutilinitrand(111397)
         do i = 1, rt_nfreq
            do j = 1, rt_nfreq
               call blutilrand(rn)
               rt_ranampl(i,j) = 2.d0*(rn-0.5d0)
               call blutilrand(rn)
               rt_ranphse(i,j,1) = 2.d0*rt_PI*rn
               call blutilrand(rn)
               rt_ranphse(i,j,2) = 2.d0*rt_PI*rn
            end do
         end do
         if ( rt_nfreq .gt. 1 ) then
            permin =  rt_nfreq**2
            permax = -rt_nfreq**2
            do i = 0, MAXPHASE
               do j = 0, MAXPHASE
                  xtmp = 2.d0*rt_PI*dble(i)/dble(MAXPHASE)
                  ytmp = 2.d0*rt_PI*dble(j)/dble(MAXPHASE)
                  pert = 0.d0
                  do n = 1, rt_nfreq
                     do m = 1, rt_nfreq
                        pert = pert &
                            + sin( &
                            2.0D0*rt_PI*dble(n)*xtmp + rt_ranphse(n,m,1) &
                            ) &
                            * sin( &
                            2.0D0*rt_PI*dble(m)*ytmp + rt_ranphse(n,m,2) &
                            ) &
                            *rt_ranampl(n,m)
                     end do
                  end do
                  permin = min(permin, pert)
                  permax = max(permax, pert)
               end do
            end do
            rt_pertamp = 2.d0*rt_pertamp/(permax - permin)
         end if
      end if

#ifdef DO_IAMR_FORCE
      if ( probtype.eq.14 .or. probtype.eq.15 ) then

         if (isioproc .eq. 1) then
            write (*,*) "Initialising random number generator..."
         endif

         twicePi = two*Pi

         if (blrandseed.gt.0) then
            call blutilinitrand(blrandseed)
            call blutilrand(rn)
            call blutilinitrand(blrandseed)
            if (isioproc .eq. 1) then
               write (*,*) "blrandseed = ",blrandseed
               write (*,*) "first random number = ",rn
            endif
         else
            call blutilinitrand(111397)
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

         if (isioproc .eq. 1) then
            write(*,*) "Lx = ",Lx
            write(*,*) "Ly = ",Ly
            write(*,*) "Lz = ",Lz
         endif

         Lmin = min(Lx,Ly,Lz)
         kappaMax = dfloat(nmodes)/Lmin + 1.0d-8
         nxmodes = nmodes*int(0.5+Lx/Lmin)
         nymodes = nmodes*int(0.5+Ly/Lmin)
         nzmodes = nmodes*int(0.5+Lz/Lmin)
         if (isioproc .eq. 1) then
            write(*,*) "Lmin = ",Lmin
            write(*,*) "kappaMax = ",kappaMax
            write(*,*) "nxmodes = ",nxmodes
            write(*,*) "nymodes = ",nymodes
            write(*,*) "nzmodes = ",nzmodes
         endif

         if (forcing_time_scale_min.eq.zero) then
            forcing_time_scale_min = half
         endif
         if (forcing_time_scale_max.eq.zero) then
            forcing_time_scale_max = one
         endif

         freqMin = one/forcing_time_scale_max
         freqMax = one/forcing_time_scale_min
         freqDiff= freqMax-freqMin

         if (isioproc .eq. 1) then
            write(*,*) "forcing_time_scale_min = ",forcing_time_scale_min
            write(*,*) "forcing_time_scale_max = ",forcing_time_scale_max
            write(*,*) "freqMin = ",freqMin
            write(*,*) "freqMax = ",freqMax
            write(*,*) "freqDiff = ",freqDiff
         endif

         mode_count = 0

         xstep = int(Lx/Lmin+0.5)
         ystep = int(Ly/Lmin+0.5)
         zstep = int(Lz/Lmin+0.5)
         if (isioproc .eq. 1) then
            write (*,*) "Mode step ",xstep, ystep, zstep
         endif

         do kz = mode_start*zstep, nzmodes, zstep
            kzd = dfloat(kz)
            do ky = mode_start*ystep, nymodes, ystep
               kyd = dfloat(ky)
               do kx = mode_start*xstep, nxmodes, xstep
                  kxd = dfloat(kx)

                  kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )

                  if (kappa.le.kappaMax) then
                     call blutilrand(rn)
                     FTX(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
                     call blutilrand(rn)
                     FTY(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
                     call blutilrand(rn)
                     FTZ(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
!c     Translation angles, theta=0..2Pi and phi=0..Pi
                     call blutilrand(rn)
                     TAT(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     TAP(kx,ky,kz) = rn*Pi
!c     Phases
                     call blutilrand(rn)
                     FPX(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     FPY(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     FPZ(kx,ky,kz) = rn*twicePi
                     if (div_free_force.eq.1) then
                        call blutilrand(rn)
                        FPXX(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPYX(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPZX(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPXY(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPYY(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPZY(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPXZ(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPYZ(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPZZ(kx,ky,kz) = rn*twicePi
                     endif
!c     Amplitudes (alpha)
                     call blutilrand(rn)
                     thetaTmp      = rn*twicePi
                     cosThetaTmp   = cos(thetaTmp)
                     sinThetaTmp   = sin(thetaTmp)
                     call blutilrand(rn)
                     phiTmp        = rn*Pi
                     cosPhiTmp     = cos(phiTmp)
                     sinPhiTmp     = sin(phiTmp)

                     px = cosThetaTmp * sinPhiTmp
                     py = sinThetaTmp * sinPhiTmp
                     pz =               cosPhiTmp

                     mp2           = px*px + py*py + pz*pz
                     if (kappa .lt. 0.000001) then
                        if (isioproc .eq. 1) then
                           write(*,*) "ZERO AMPLITUDE MODE ",kx,ky,kz
                        endif
                        FAX(kx,ky,kz) = zero
                        FAY(kx,ky,kz) = zero
                        FAZ(kx,ky,kz) = zero
                     else
!c     Count modes that contribute
                        mode_count = mode_count + 1
!c     Set amplitudes
                        if (spectrum_type.eq.1) then
                           Ekh        = one / kappa
                        else if (spectrum_type.eq.2) then
                           Ekh        = one / (kappa*kappa)
                        else
                           Ekh        = one
                        endif
                        if (div_free_force.eq.1) then
                           Ekh = Ekh / kappa
                        endif
                        if (moderate_zero_modes.eq.1) then
                           if (kx.eq.0) Ekh = Ekh / two
                           if (ky.eq.0) Ekh = Ekh / two
                           if (kz.eq.0) Ekh = Ekh / two
                        endif
                        if (force_scale.gt.zero) then
                           FAX(kx,ky,kz) = force_scale * px * Ekh / mp2
                           FAY(kx,ky,kz) = force_scale * py * Ekh / mp2
                           FAZ(kx,ky,kz) = force_scale * pz * Ekh / mp2
                        else
                           FAX(kx,ky,kz) = px * Ekh / mp2
                           FAY(kx,ky,kz) = py * Ekh / mp2
                           FAZ(kx,ky,kz) = pz * Ekh / mp2
                        endif

                        if (isioproc.eq.1) then
                           write (*,*) "Mode"
                           write (*,*) "kappa = ",kx,ky,kz,kappa, sqrt(FAX(kx,ky,kz)**2+FAY(kx,ky,kz)**2+FAZ(kx,ky,kz)**2)
                           write (*,*) "Amplitudes - A"
                           write (*,*) FAX(kx,ky,kz), FAY(kx,ky,kz), FAZ(kx,ky,kz)
                           write (*,*) "Frequencies"
                           write (*,*) FTX(kx,ky,kz), FTY(kx,ky,kz), FTZ(kx,ky,kz)
                        endif
                     endif
                  endif
               enddo
            enddo
         enddo

!c     Now let's break symmetry, have to assume high aspect ratio in z for now
         reduced_mode_count = 0
         do kz = 1, zstep - 1
            kzd = dfloat(kz)
            do ky = mode_start, nymodes
               kyd = dfloat(ky)
               do kx = mode_start, nxmodes
                  kxd = dfloat(kx)

                  kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )

                  if (kappa.le.kappaMax) then
                     call blutilrand(rn)
                     FTX(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
                     call blutilrand(rn)
                     FTY(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
                     call blutilrand(rn)
                     FTZ(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
!c     Translation angles, theta=0..2Pi and phi=0..Pi
                     call blutilrand(rn)
                     TAT(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     TAP(kx,ky,kz) = rn*Pi
!c     Phases
                     call blutilrand(rn)
                     FPX(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     FPY(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     FPZ(kx,ky,kz) = rn*twicePi
                     if (div_free_force.eq.1) then
                        call blutilrand(rn)
                        FPXX(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPYX(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPZX(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPXY(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPYY(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPZY(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPXZ(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPYZ(kx,ky,kz) = rn*twicePi
                        call blutilrand(rn)
                        FPZZ(kx,ky,kz) = rn*twicePi
                     endif
!c     Amplitudes (alpha)
                     call blutilrand(rn)
                     thetaTmp      = rn*twicePi
                     cosThetaTmp   = cos(thetaTmp)
                     sinThetaTmp   = sin(thetaTmp)
                     call blutilrand(rn)
                     phiTmp        = rn*Pi
                     cosPhiTmp     = cos(phiTmp)
                     sinPhiTmp     = sin(phiTmp)

                     px = cosThetaTmp * sinPhiTmp
                     py = sinThetaTmp * sinPhiTmp
                     pz =               cosPhiTmp

                     mp2           = px*px + py*py + pz*pz
                     if (kappa .lt. 0.000001) then
                        if (isioproc .eq. 1) then
                           write(*,*) "ZERO AMPLITUDE MODE ",kx,ky,kz
                        endif
                        FAX(kx,ky,kz) = zero
                        FAY(kx,ky,kz) = zero
                        FAZ(kx,ky,kz) = zero
                     else
!c     Count modes that contribute
                        reduced_mode_count = reduced_mode_count + 1
!c     Set amplitudes
                        if (spectrum_type.eq.1) then
                           Ekh        = one / kappa
                        else if (spectrum_type.eq.2) then
                           Ekh        = one / (kappa*kappa)
                        else
                           Ekh        = one
                        endif
                        if (div_free_force.eq.1) then
                           Ekh = Ekh / kappa
                        endif
                        if (moderate_zero_modes.eq.1) then
                           if (kx.eq.0) Ekh = Ekh / two
                           if (ky.eq.0) Ekh = Ekh / two
                           if (kz.eq.0) Ekh = Ekh / two
                        endif
                        if (force_scale.gt.zero) then
                           FAX(kx,ky,kz) = forcing_epsilon * force_scale * px * Ekh / mp2
                           FAY(kx,ky,kz) = forcing_epsilon * force_scale * py * Ekh / mp2
                           FAZ(kx,ky,kz) = forcing_epsilon * force_scale * pz * Ekh / mp2
                        else
                           FAX(kx,ky,kz) = forcing_epsilon * px * Ekh / mp2
                           FAY(kx,ky,kz) = forcing_epsilon * py * Ekh / mp2
                           FAZ(kx,ky,kz) = forcing_epsilon * pz * Ekh / mp2
                        endif

                        if (isioproc.eq.1) then
                           write (*,*) "Mode"
                           write (*,*) "kappa = ",kx,ky,kz,kappa, sqrt(FAX(kx,ky,kz)**2+FAY(kx,ky,kz)**2+FAZ(kx,ky,kz)**2)
                           write (*,*) "Amplitudes - B"
                           write (*,*) FAX(kx,ky,kz), FAY(kx,ky,kz), FAZ(kx,ky,kz)
                           write (*,*) "Frequencies"
                           write (*,*) FTX(kx,ky,kz), FTY(kx,ky,kz), FTZ(kx,ky,kz)
                        endif
                     endif
                  endif
               enddo
            enddo
         enddo

         if (isioproc .eq. 1) then
            write(*,*) "mode_count = ",mode_count
            write(*,*) "reduced_mode_count = ",reduced_mode_count
            if (spectrum_type.eq.1) then
               write (*,*) "Spectrum type 1"
            else if (spectrum_type.eq.2) then
               write (*,*) "Spectrum type 2"
            else
               write (*,*) "Spectrum type OTHER"
            endif
         endif

!c Override random numbers
         if (override_turb_force.ge.1) then

            if (override_turb_force.eq.1) then
               do kz=0,32
                  do ky=0,32
                     do kx=0,32
                        FAX(kx,ky,kz) = 0.d0
                        FAY(kx,ky,kz) = 0.d0
                        FAZ(kx,ky,kz) = 0.d0
                        FTX(kx,ky,kz) = 0.d0
                        FTY(kx,ky,kz) = 0.d0
                        FTZ(kx,ky,kz) = 0.d0
                        FPX(kx,ky,kz) = 0.d0
                        FPY(kx,ky,kz) = 0.d0
                        FPZ(kx,ky,kz) = 0.d0
                        TAT(kx,ky,kz) = 0.d0
                     enddo
                  enddo
               enddo
            endif

            FAX(1,1,1) = 4.9159487814194529d0
            FAY(1,1,1) = 8.4733283347910291d0
            FAZ(1,1,1) = -4.2206554333498882d0
            FTX(1,1,1) = 5.1975883412614365d0
            FTY(1,1,1) =  5.9058228872665275d0
            FTZ(1,1,1) = 3.6304842628015055d0
            FPX(1,1,1) = 5.8021873618124946d0
            FPY(1,1,1) = 5.0447187825551207d0
            FPZ(1,1,1) = 1.8627240964262393d0
            TAT(1,1,1) = 5.6478239516237787d0

            FAX(2,1,1) = -2.5666438215760063d0
            FAY(2,1,1) = -1.4239169102661757d0
            FAZ(2,1,1) = -4.4530039939649946d0
            FTX(2,1,1) = 3.7679884103646226d0
            FTY(2,1,1) = 5.3250653798614538d0
            FTZ(2,1,1) = 4.5150594485888291d0
            FPX(2,1,1) = 6.2749868772687361d0
            FPY(2,1,1) = 1.3705682802978747d0
            FPZ(2,1,1) = 1.0304654059116896d0
            TAT(2,1,1) = 2.8651852272695173d0

            FAX(3,1,1) = 1.2864170156656984d0
            FAY(3,1,1) = -2.2893975400193582d0
            FAZ(3,1,1) = 1.2516389586915393d0
            FTX(3,1,1) = 5.9033711974860932d0
            FTY(3,1,1) = 4.4281708270687501d0
            FTZ(3,1,1) = 4.2777444944303191d0
            FPX(3,1,1) = 0.22011703989520284d0
            FPY(3,1,1) = 4.2466220672604873d0
            FPZ(3,1,1) = 4.8764721700720708d0
            TAT(3,1,1) = 5.0553217246012565d0

            FAX(1,2,1) = -0.20089188587586918d0
            FAY(1,2,1) = 1.3763908645879952d0
            FAZ(1,2,1) = 5.1487508273864258d0
            FTX(1,2,1) = 3.9485897668695142d0
            FTY(1,2,1) = 4.4692744553694244d0
            FTZ(1,2,1) = 3.9508792365899366d0
            FPX(1,2,1) = 1.7852441923049107d0
            FPY(1,2,1) = 3.0829934386493503d0
            FPZ(1,2,1) = 1.4672406821539532d0
            TAT(1,2,1) = 1.4755814993444105d0

            FAX(2,2,1) = 0.69396713816146716d0
            FAY(2,2,1) = 2.861361830195257d0
            FAZ(2,2,1) = -1.9932369142918278d0
            FTX(2,2,1) = 5.6956730220436205d0
            FTY(2,2,1) = 5.5588867839730547d0
            FTZ(2,2,1) = 5.7649122186530564d0
            FPX(2,2,1) = 2.5552810738283052d0
            FPY(2,2,1) = 3.388311033772359d0
            FPZ(2,2,1) = 1.3010559678548617d0
            TAT(2,2,1) = 1.0339095056702137d0

            FAX(3,2,1) = 1.9350279711990759d0
            FAY(3,2,1) = -1.1279468985303027d0
            FAZ(3,2,1) = -0.45595212543799785d0
            FTX(3,2,1) = 4.772596551965731d0
            FTY(3,2,1) = 3.8777674032401634d0
            FTZ(3,2,1) = 5.3116096780714743d0
            FPX(3,2,1) = 1.780426633370265d0
            FPY(3,2,1) = 5.475785414516678d0
            FPZ(3,2,1) = 4.4419968980233211d0
            TAT(3,2,1) = 0.22540787623529085d0

            FAX(1,3,1) = 2.864674922606242d0
            FAY(1,3,1) = -0.50544520889614819d0
            FAZ(1,3,1) = 3.11872722472469888d-2
            FTX(1,3,1) = 3.7635706940051188d0
            FTY(1,3,1) = 5.208809484692571d0
            FTZ(1,3,1) = 4.7738469453527825d0
            FPX(1,3,1) = 1.864303985046738d0
            FPY(1,3,1) = 6.2540331417366657d0
            FPZ(1,3,1) = 1.6475935810409588d0
            TAT(1,3,1) = 3.3634507598363745d0

            FAX(2,3,1) = 0.77945529357360954d0
            FAY(2,3,1) = -0.91747387637372091d0
            FAZ(2,3,1) = 1.9429824825278852d0
            FTX(2,3,1) = 5.742158566501101d0
            FTY(2,3,1) = 3.7482952632957014d0
            FTZ(2,3,1) = 5.4874587715272858d0
            FPX(2,3,1) = 2.9531052813760001d0
            FPY(2,3,1) = 0.26215970214491496d0
            FPZ(2,3,1) = 1.7493938062755146d0
            TAT(2,3,1) = 1.91306574518684d0

            FAX(1,1,2) = -5.305003847200731d0
            FAY(1,1,2) = 0.50749773957239552d0
            FAZ(1,1,2) = 0.2093434258791288d0
            FTX(1,1,2) = 5.8647362046774285d0
            FTY(1,1,2) = 4.1987859158702072d0
            FTZ(1,1,2) = 3.9029081485618291d0
            FPX(1,1,2) = 5.1002773317760663d0
            FPY(1,1,2) = 2.1221436871767354d0
            FPZ(1,1,2) = 0.6975038860903926d0
            TAT(1,1,2) = 0.44023114810257413d0

            FAX(2,1,2) = -2.4840532177547301d0
            FAY(2,1,2) = -2.2711993844331797d0
            FAZ(2,1,2) = 1.1459093664660789d0
            FTX(2,1,2) = 4.3877095778001758d0
            FTY(2,1,2) = 3.7581829005749277d0
            FTZ(2,1,2) = 4.8766734961267275d0
            FPX(2,1,2) = 4.7357438644376497d0
            FPY(2,1,2) = 1.4783784275512053d0
            FPZ(2,1,2) = 5.3774890550958716d0
            TAT(2,1,2) = 5.7023832517246111d0

            FAX(3,1,2) = 0.29669381505201736d0
            FAY(3,1,2) = -1.9473540911634448d0
            FAZ(3,1,2) = 1.1594285746251218d0
            FTX(3,1,2) = 5.8542809412855927d0
            FTY(3,1,2) = 3.8167045618392166d0
            FTZ(3,1,2) = 5.4244408688497341d0
            FPX(3,1,2) = 3.1993678345714365d0
            FPY(3,1,2) = 0.88139277647217429d0
            FPZ(3,1,2) = 5.4657479049374924d0
            TAT(3,1,2) = 1.3614786992096499d0

            FAX(1,2,2) = -1.9318740901427485d0
            FAY(1,2,2) = -0.80218340917884778d0
            FAZ(1,2,2) = 2.8751242732298206d0
            FTX(1,2,2) = 5.8773266401625346d0
            FTY(1,2,2) = 4.4329190612414395d0
            FTZ(1,2,2) = 3.2881333624326348d0
            FPX(1,2,2) = 5.003162786337116d0
            FPY(1,2,2) = 1.2164277369010781d0
            FPZ(1,2,2) = 0.14669665298429169d0
            TAT(1,2,2) = 5.5595175028166013d0

            FAX(2,2,2) = 0.4297562850971306d0
            FAY(2,2,2) = -1.9466539669692717d0
            FAZ(2,2,2) = 1.7711462332098464d0
            FTX(2,2,2) = 4.6865618035445218d0
            FTY(2,2,2) = 5.9494078155410257d0
            FTZ(2,2,2) = 4.0205365556673778d0
            FPX(2,2,2) = 4.5491812147738679d0
            FPY(2,2,2) = 1.6761367863523891d0
            FPZ(2,2,2) = 5.2026905021309613d0
            TAT(2,2,2) = 4.6252688992576214d0

            FAX(1,3,2) = 2.06658971082688374d-2
            FAY(1,3,2) = 2.74013268199971953d-2
            FAZ(1,3,2) = -2.2854566029359611d0
            FTX(1,3,2) = 5.8157142214007846d0
            FTY(1,3,2) = 4.9600077380032275d0
            FTZ(1,3,2) = 6.2077903202581997d0
            FPX(1,3,2) = 1.5010646369853131d0
            FPY(1,3,2) = 5.3916183465437202d0
            FPZ(1,3,2) = 4.6437631037683511d0
            TAT(1,3,2) = 3.0680376745582629d0

            FAX(1,1,3) = 2.2949361665546917d0
            FAY(1,1,3) = 1.4268107896343492d0
            FAZ(1,1,3) = 1.0771670619628302d0
            FTX(1,1,3) = 3.4163673645862564d0
            FTY(1,1,3) = 3.4248811396554295d0
            FTZ(1,1,3) = 4.6260711993739712d0
            FPX(1,1,3) = 4.5109620235136738d0
            FPY(1,1,3) = 5.8391370634844844d0
            FPZ(1,1,3) = 3.4206386516235008d0
            TAT(1,1,3) = 5.1764802854269849d0

            FAX(2,1,3) = 0.72905204823447511d0
            FAY(2,1,3) = -0.25947497300324907d0
            FAZ(2,1,3) = 2.1507314209980866d0
            FTX(2,1,3) = 3.8967963050404482d0
            FTY(2,1,3) = 6.1779977397507002d0
            FTZ(2,1,3) = 5.9597973778189886d0
            FPX(2,1,3) = 2.991726230902561d0
            FPY(2,1,3) = 5.4725093396771065d0
            FPZ(2,1,3) = 0.95909300483380011d0
            TAT(2,1,3) = 5.4692810554453581d0

            FAX(1,2,3) = -2.1551885177430785d0
            FAY(1,2,3) = 0.69991730842965261d0
            FAZ(1,2,3) = 0.29961310096080823d0
            FTX(1,2,3) = 4.0690501771057352d0
            FTY(1,2,3) = 3.5887121197057796d0
            FTZ(1,2,3) = 5.0908129869173395d0
            FPX(1,2,3) = 2.2658198780716678d0
            FPY(1,2,3) = 1.3269114298503693d0
            FPZ(1,2,3) = 2.5274702857680205d0
            TAT(1,2,3) = 3.3730006261545773d0

            if (isioproc.eq.1) then
               do kz = mode_start*zstep, nzmodes, zstep
                  kzd = dfloat(kz)
                  do ky = mode_start*ystep, nymodes, ystep
                     kyd = dfloat(ky)
                     do kx = mode_start*xstep, nxmodes, xstep
                        kxd = dfloat(kx)
                        kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                        if (kappa.le.kappaMax) then
                           write (*,*) "Mode"
                           write (*,*) "kappa = ",kx,ky,kz,kappa
                           write (*,*) "Amplitudes - C"
                           write (*,*) FAX(kx,ky,kz), FAY(kx,ky,kz), FAZ(kx,ky,kz)
                           write (*,*) "Frequencies"
                           write (*,*) FTX(kx,ky,kz), FTY(kx,ky,kz), FTZ(kx,ky,kz)
                           write (*,*) "Phases"
                           write (*,*) FPX(kx,ky,kz), FPY(kx,ky,kz), FPZ(kx,ky,kz)
                           write (*,*) "TAT"
                           write (*,*) TAT(kx,ky,kz)
                        endif
                     enddo
                  enddo
               enddo
            endif
!c     override
         endif

!c     probtype 14 or 15
      endif

!c     iamr_force
#endif

      if (probtype.eq.23) then
         Lx = f_probhi(1)-f_problo(1)
         Ly = f_probhi(2)-f_problo(2)
         Lz = f_probhi(3)-f_problo(3)

!c     Initialise random number generator
         idum=-1
         f0 = ran1(idum)

!c     Get frequencies higher than fundamental
         f0 = half*(vel1+vel2)/(five*delta0)
         fmin = f0/five
         fmax = f0
         fdif = fmax-fmin

!c     Get random numbers
         do n = 1, 10
            mag(n)  = amag*ran1(idum)
            freq(n) = ran1(idum)*(fmax-fmin)+fmin
            phi1(n) = ran1(idum)*two*Pi
            phi2(n) = ran1(idum)*two*Pi
         end do
      endif

      end subroutine amrex_probinit

!c ::: -----------------------------------------------------------
!c ::: This routine is called at problem setup time and is used
!c ::: to initialize data on each grid.  The velocity field you
!c ::: provide does not have to be divergence free and the pressure
!c ::: field need not be set.  A subsequent projection iteration
!c ::: will define aa divergence free velocity field along with a
!c ::: consistant pressure.
!c :::
!c ::: NOTE:  all arrays have one cell of ghost zones surrounding
!c :::        the grid interior.  Values in these cells need not
!c :::        be set here.
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: level     => amr level of grid
!c ::: time      => time at which to init data
!c ::: lo,hi     => index limits of grid interior (cell centered)
!c ::: nscal     => number of scalar quantities.  You should know
!c :::		   this already!
!c ::: vel      <=  Velocity array
!c ::: scal     <=  Scalar array
!c ::: press    <=  Pressure array
!c ::: dx     => cell size
!c ::: xlo,xhi   => physical locations of lower left and upper
!c :::              right hand corner of grid.  (does not include
!c :::		   ghost region).
!c ::: -----------------------------------------------------------
      subroutine FORT_INITDATA(level,time,lo,hi,nscal, &
                               vel,scal,DIMS(state),press,DIMS(press), &
                               dx,xlo,xhi) &
                               bind(C, name="FORT_INITDATA")

      implicit none
      integer    level, nscal
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(state)
      integer    DIMDEC(press)
      REAL_T     time, dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     vel(DIMV(state),SDIM)
      REAL_T    scal(DIMV(state),nscal)
      REAL_T   press(DIMV(press))


#include <probdata.H>

!c      print *, 'INITDATA ', time, lo(1),hi(1), lo(2), hi(2), lo(3), hi(3)

         call initflowpastcylinder(level,time,lo,hi,nscal, &
          vel,scal,DIMS(state),press,DIMS(press), &
          dx,xlo,xhi, m_probhi )


      end subroutine FORT_INITDATA

!c----------------------------------------------------------------------
!c     A handy statement function: here blend goes from 0 to 1 in x at
!c     rad, over a width of trn

      DOUBLE PRECISION function zblend1(x,rad,trn)
      DOUBLE PRECISION x, rad, trn
      if ( trn .eq. 0.0d0 ) then
         if ( x .lt. rad ) then
            zblend1 = 0.0
         else
            zblend1 = 1.0
         end if
      else
         zblend1 = 0.5D0*(1.0D0 + TANH((x-rad)/trn))
      end if
    end function zblend1

!c
!c ::: -----------------------------------------------------------
!c
      subroutine initflowpastcylinder ( level,time,lo,hi,nscal, &
       &                                vel,scal,DIMS(state),press,DIMS(press), &
       &                                dx,xlo,xhi,probhi)

      implicit none

      integer    level, nscal
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(state)
      integer    DIMDEC(press)
      REAL_T     time, dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM), probhi(SDIM)
      REAL_T     vel(DIMV(state),SDIM)
      REAL_T     scal(DIMV(state),nscal)
      REAL_T     press(DIMV(press))
!c
!c     ::::: local variables
!c
      integer i, j, k, n
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  dist

#include <probdata.H>

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

         do k = lo(3), hi(3)
            z = xlo(3) + hz*(float(k-lo(3)) + half)
            do j = lo(2), hi(2)
               y = xlo(2) + hy*(float(j-lo(2)) + half)
               do i = lo(1), hi(1)
                  x = xlo(1) + hx*(float(i-lo(1)) + half)

                  vel(i,j,k,1) = adv_vel 
                  vel(i,j,k,2) = zero
                  vel(i,j,k,3) = zero

                  scal(i,j,k,1) = denfact

                  do n = 2,nscal-1
                     scal(i,j,k,n) = one
                  end do

                  scal(i,j,k,nscal) = 1.
               end do
            end do
         end do

   end subroutine initflowpastcylinder


!c ::: -----------------------------------------------------------
!c
      subroutine swirl(level,time,lo,hi,nscal, &
     	 	       vel,scal,DIMS(state),press,DIMS(press), &
                      dx,xlo,xhi)
      implicit none

      integer    level, nscal
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(state)
      integer    DIMDEC(press)
      REAL_T     time, dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     vel(DIMV(state),SDIM)
      REAL_T    scal(DIMV(state),nscal)
      REAL_T   press(DIMV(press))
!c
!c     ::::: local variables
!c
      integer i, j, k, n
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  dist

#include <probdata.H>

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
	    do i = lo(1), hi(1)

               x = xlo(1) + hx*(float(i-lo(1)) + half)

               vel(i,j,k,1) = 0.d0
               vel(i,j,k,2) = 0.d0
               vel(i,j,k,3) = Vco

               scal(i,j,k,1) = denfact
               do n = 2,nscal-1
                  scal(i,j,k,n) = one
               end do

	    end do
         end do
      end do

      end subroutine swirl

!c     ::: -----------------------------------------------------------
!c     ::: This routine will tag high error cells based on the
!c     ::: magnitude or gradient of the density
!c     :::
!c     ::: INPUTS/OUTPUTS:
!c     :::
!c     ::: tag      <=  integer tag array
!c     ::: DIMS(tag) => index extent of tag array
!c     ::: set       => integer value to tag cell for refinement
!c     ::: clear     => integer value to untag cell
!c     ::: rho       => density array
!c     ::: DIMS(rho) => index extent of rho array
!c     ::: nvar      => number of components in rho array (should be 1)
!c     ::: lo,hi     => index extent of grid
!c     ::: domlo,hi  => index extent of problem domain
!c     ::: dx        => cell spacing
!c     ::: xlo       => physical location of lower left hand
!c     :::	           corner of tag array
!c     ::: problo    => phys loc of lower left corner of prob domain
!c     ::: time      => problem evolution time
!c     ::: -----------------------------------------------------------
      subroutine FORT_DENERROR (tag,DIMS(tag),set,clear, &
                                rho,DIMS(rho),lo,hi,nvar, &
                                domlo,domhi,dx,xlo,  &
                                problo,time,level)&
                                bind(C, name="FORT_DENERROR")
      implicit none
      integer   DIMDEC(tag)
      integer   DIMDEC(rho)
      integer   lo(SDIM), hi(SDIM)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    rho(DIMV(rho),nvar)

      integer   i, j, k
      REAL_T    ax, ay, az, aerr

#include <probdata.H>

      if ( probtype .eq. 10 ) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(rho(i+1,j,k,1) - rho(i-1,j,k,1))
                  ay = ABS(rho(i,j+1,k,1) - rho(i,j-1,k,1))
                  az = ABS(rho(i,j,k+1,1) - rho(i,j,k-1,1))
                  aerr = MAX(ax,ay,az)
                  if ( aerr .GE. rt_graddenerr ) tag(i,j,k) = set
               end do
            end do
         end do
      else if ( probtype .eq. 19 ) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k),rho(i,j,k,1).gt.denerr)
               end do
            end do
         end do
      else if ( probtype .eq. 32 ) then
         select case (adv_dir)
         case (1)
            do i = lo(1), hi(1)
               if (i<domhi(1)/2) then
                  tag(i,:,:) = set
               end if
            end do
         case (2)
            do j = lo(2), hi(2)
               if (i<domhi(2)/2) then
                  tag(:,j,:) = set
               end if
            end do
         case (3)
            do k = lo(3), hi(3)
               if (i<domhi(3)/2) then
                  tag(:,:,k) = set
               end if
            end do
         end select
      else
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k),rho(i,j,k,1).lt.denerr)
               end do
            end do
         end do
      endif

      end subroutine FORT_DENERROR

!c ::: -----------------------------------------------------------
!c ::: This routine will tag high error cells based on the
!c ::: magnitude of the tracer
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: tag      <=  integer tag array
!c ::: DIMS(tag) => index extent of tag array
!c ::: set       => integer value to tag cell for refinement
!c ::: clear     => integer value to untag cell
!c ::: adv       => scalar array
!c ::: DIMS(adv) => index extent of adv array
!c ::: nvar      => number of components in rho array (should be 1)
!c ::: lo,hi     => index extent of grid
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of tag array
!c ::: problo    => phys loc of lower left corner of prob domain
!c ::: time      => problem evolution time
!c ::: -----------------------------------------------------------
      subroutine FORT_ADVERROR (tag,DIMS(tag),set,clear, &
                               adv,DIMS(adv),lo,hi,nvar, &
                               domlo,domhi,dx,xlo, &
     			                     problo,time,level)&
                               bind(C, name="FORT_ADVERROR")
      implicit none
      integer   DIMDEC(tag)
      integer   DIMDEC(adv)
      integer   lo(SDIM), hi(SDIM)
      integer   ng, nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    adv(DIMV(adv),nvar)

      REAL_T    x, y, z, ax, ay, az, aerr, rh
      integer   i, j, k

#include <probdata.H>

!c     probtype = SPIN
      if (probtype .eq. 1) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

!c     probtype = BUBBLE
      else if (probtype .eq. 2) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

!c     probtype = VORTEX IN A BOX
      else if (probtype .eq. 3) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

!c     probtype = CHANNEL
      else if (probtype .eq. 4) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

!c     probtype = PERIODIC SHEAR LAYER
      else if (probtype .eq. 5) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

!c     probtype = HOT SPOT
      else if (probtype .eq. 6) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

      else if (probtype .eq. 7) then

        return

      else if (probtype .eq. 8) then

        if (level .eq. 0) then
          do k = lo(3), hi(3)
          z = xlo(3) + dx(3)*(float(k-lo(3)) + half)
             do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                  x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                  tag(i,j,k) = merge(set,clear,abs(z).lt.zlev1.and.x.lt.xlev1)
             enddo
             enddo
          enddo
        else
          do k = lo(3), hi(3)
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
            ax = abs(adv(i+1,j,k,1) - adv(i-1,j,k,1))
            ay = abs(adv(i,j+1,k,1) - adv(i,j-1,k,1))
            az = abs(adv(i,j,k+1,1) - adv(i,j,k-1,1))
            aerr = max(ax,ay,az)
            tag(i,j,k) = merge(set,tag(i,j,k),aerr.gt.adverr.and. x.lt.xlev1)
          enddo
          enddo
          enddo
        endif

!c     probtype = VISCOUS BENCHMARK
      else if (probtype .eq. 9) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do


      else if (probtype .eq. 10.or.probtype.eq.11) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

      else if (probtype .eq. 12) then

      else if (probtype .eq. 13) then

      else if (probtype .eq. 14) then

         do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

      else if (probtype .eq. 15) then

      else if (probtype .eq. 16) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

      else if (probtype .eq. 17) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),abs(adv(i,j,k,1)).gt.adverr)
              end do
           end do
        end do

        if (time.lt.injection_time) then
           if (lo(3).eq.domlo(3)) then
              k=lo(3)
              do j = lo(2), hi(2)
                 y = xlo(2) + dx(2)*(float(j-lo(2)) + half)
                 do i = lo(1), hi(1)
                    x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                    if (sqrt(x*x+y*y).lt.jet_width) then
                       tag(i,j,k) = set
                    endif
                 enddo
              enddo
           endif
           if (hi(3).eq.domhi(3)) then
              k=hi(3)
              do j = lo(2), hi(2)
                 y = xlo(2) + dx(2)*(float(j-lo(2)) + half)
                 do i = lo(1), hi(1)
                    x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                    if (sqrt(x*x+y*y).lt.jet_width) then
                       tag(i,j,k) = set
                    endif
                 enddo
              enddo
           endif
        endif

      else if (probtype .eq. 18) then

        if (ref_height2.gt.0 .and. level.eq.1) then
           rh = ref_height2
        else
           rh = ref_height
        endif

        do k = lo(3), hi(3)
           z = xlo(3) + dx(3)*(float(k-lo(3)) + half)
           if (z.lt.rh) then
              do j = lo(2), hi(2)
                 do i = lo(1), hi(1)
                    tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
                 end do
              end do
           endif
        end do

      else if (probtype .eq. 19) then

        do k = lo(3), hi(3)
           z = xlo(3) + dx(3)*(float(k-lo(3)) + half)
           if (z.lt.ref_height) then
              do j = lo(2), hi(2)
                 do i = lo(1), hi(1)
                    tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
                 end do
              end do
           endif
        end do

      else if (probtype .eq. 20) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

      else if (probtype .eq. 21) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

      else if (probtype .eq. 22) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

      else if (probtype .eq. 23) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

      else if (probtype .eq. 24) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

      else if (probtype .eq. 25) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

      else if (probtype .eq. 31) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

      else
        print *,'DONT KNOW THIS PROBTYPE IN FORT_ADVERROR ',probtype
        stop
      end if

      end subroutine FORT_ADVERROR

      subroutine FORT_ADV2ERROR (tag,DIMS(tag),set,clear, &
                                 adv,DIMS(adv),lo,hi,nvar, &
                                 domlo,domhi,dx,xlo, &
     		 	                       problo,time,level)&
                                 bind(C, name="FORT_ADV2ERROR")
      implicit none
      integer   DIMDEC(tag)
      integer   DIMDEC(adv)
      integer   lo(SDIM), hi(SDIM)
      integer   ng, nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    adv(DIMV(adv),nvar)

      REAL_T    x, y, z, ax, ay, az, aerr
      integer   i, j, k

#include <probdata.H>

      if (probtype .eq. 20) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
               end do
            end do
         end do
      else
         print *,'DONT KNOW THIS PROBTYPE IN FORT_ADVERROR ',probtype
         stop
      end if

      end subroutine FORT_ADV2ERROR

!c ::: -----------------------------------------------------------
!c ::: This routine will tag high error cells based on the
!c ::: magnitude or gradient of temperature
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: tag      <=  integer tag array
!c ::: DIMS(tag) => index extent of tag array
!c ::: set       => integer value to tag cell for refinement
!c ::: clear     => integer value to untag cell
!c ::: temp      => density array
!c ::: DIMS(temp)=> index extent of temp array
!c ::: lo,hi     => index extent of grid
!c ::: nvar      => number of components in rho array (should be 1)
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::              corner of tag array
!c ::: problo    => phys loc of lower left corner of prob domain
!c ::: time      => problem evolution time
!c ::: -----------------------------------------------------------
      subroutine FORT_TEMPERROR (tag,DIMS(tag),set,clear, &
                                 temperature,DIMS(temp),lo,hi,nvar, &
                                 domlo,domhi,dx,xlo, &
                                 problo,time,level)&
                                 bind(C, name="FORT_TEMPERROR")
      implicit none

      integer   DIMDEC(tag)
      integer   DIMDEC(temp)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      integer   lo(SDIM), hi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    temperature(DIMV(temp),nvar)

      REAL_T    x, y, z, ax, ay, az, aerr
      integer   i, j, k

#include <probdata.H>

!c     probtype = SPIN
      if (probtype .eq. 1) then

!c     probtype = BUBBLE
      else if (probtype .eq. 2) then

!c     probtype = VORTEX IN A BOX
      else if (probtype .eq. 3) then

!c     probtype = CHANNEL
      else if (probtype .eq. 4) then

!c     probtype = PERIODIC SHEAR LAYER
      else if (probtype .eq. 5) then

!c     probtype = HOT SPOT
      else if (probtype .eq. 6) then

        if (level .eq. 0) then
!c         ::::: refine around entire hot spot
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   tag(i,j,k) = merge(set,tag(i,j,k),temperature(i,j,k,1).gt.temperr)
                end do
             end do
          end do
        else
!c         ::::: refine where there is temperature gradient
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   ax = abs(temperature(i+1,j,k,1) - temperature(i-1,j,k,1))
                   ay = abs(temperature(i,j+1,k,1) - temperature(i,j-1,k,1))
                   az = abs(temperature(i,j,k+1,1) - temperature(i,j,k-1,1))
                   aerr = max(ax,ay,az)
                   tag(i,j,k) = merge(set,tag(i,j,k),aerr.gt.bubgrad)
                end do
             end do
          end do
       end if

!c     probtype = VISCOUS BENCHMARK
      else if (probtype .eq. 9) then

      else if (probtype .eq. 10.or. probtype .eq. 11 .or.probtype.eq.12) then

      else if (probtype .eq. 13) then

      else if (probtype .eq. 14) then

      else if (probtype .eq. 15) then

      else if (probtype .eq. 16) then

      else if (probtype .eq. 17) then

      else if (probtype .eq. 18) then

      else if (probtype .eq. 19) then

      else if (probtype .eq. 20) then

      else if (probtype .eq. 21) then

      else if (probtype .eq. 22) then

      else if (probtype .eq. 23) then

      else if (probtype .eq. 24) then

      else if (probtype .eq. 25) then

      else
        print *,'DONT KNOW THIS PROBTYPE IN FORT_TEMPERROR ',probtype
        stop
      end if

      end subroutine FORT_TEMPERROR
!c ::: -----------------------------------------------------------
!c ::: This routine will tag high error cells based on the
!c ::: magnitude of vorticity
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: tag      <=  integer tag array
!c ::: DIMS(tag) => index extent of tag array
!c ::: set       => integer value to tag cell for refinement
!c ::: clear     => integer value to untag cell
!c ::: vort      => vorticitiy array
!c ::: DIMS(vort)=> index extent of vort array
!c ::: nvar      => number of components in rho array (should be 1)
!c ::: lo,hi     => index extent of grid
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of tag array
!c ::: problo    => phys loc of lower left corner of prob domain
!c ::: time      => problem evolution time
!c ::: -----------------------------------------------------------
      subroutine FORT_MVERROR (tag,DIMS(tag),set,clear, &
                               vort,DIMS(vort),lo,hi,nvar, &
                               domlo,domhi,dx,xlo, &
                               problo,time,level) &
                               bind(C, name="FORT_MVERROR")
      implicit none

      integer   DIMDEC(tag)
      integer   DIMDEC(vort)
      integer   lo(SDIM), hi(SDIM)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    vort(DIMV(vort),nvar)

      REAL_T    x, y, z, dist
      integer   i, j, k, ztag
      REAL_T    radius, maxvort

#include <probdata.H>

!c      write (*,*) "MVERROR: probtype ",probtype," on level ",level

!c     probtype = SPIN
      if (probtype .eq. 1) then

!c     probtype = BUBBLE
      else if (probtype .eq. 2) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

!c     probtype = VORTEX IN A BOX
      else if (probtype .eq. 3) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

!c     probtype = CHANNEL
      else if (probtype .eq. 4) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

!c     probtype = PERIODIC SHEAR LAYER
      else if (probtype .eq. 5) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

!c     probtype = HOT SPOT
      else if (probtype .eq. 6) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do
!c     probtype = EULER
      else if (probtype .eq. 7) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k), &
                        abs(dx(1)*vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

!c     probtype = BROWN ROSHKO
      else if (probtype .eq. 8) then

        return


!c     probtype = VISCOUS BENCHMARK
      else if (probtype .eq. 9) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

      else if (probtype .eq. 10 .or. probtype .eq. 11 .or.probtype.eq.12) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

      else if (probtype .eq. 13) then

      else if (probtype .eq. 14) then
         if (ref_type.eq.1) then
            do k = lo(3), hi(3)
               z = xlo(3) + dx(3)*(float(k-lo(3)) + half)
               ztag = 0
               if (z.lt.4.5) then
                  if (z.gt.3.5 .or. level.eq.0) then
                     ztag = 1
                  endif
               endif
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tag(i,j,k) = merge(set,tag(i,j,k),ztag.eq.1)
                  end do
               end do
            end do
         else if (ref_type.eq.2) then
            write (*,*) "Refining probtype 14 on level",level
            do k = lo(3), hi(3)
               z = xlo(3) + dx(3)*(float(k-lo(3)) + half)
               if (abs(z-ref_centre).lt.ref_radius) then
                  ztag = 1
               else
                  ztag = 0
               endif
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tag(i,j,k) = merge(set,tag(i,j,k),ztag.eq.1)
                  end do
               end do
            end do
         else
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
                  end do
               end do
            end do
         endif

      else if (probtype .eq. 15) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

         if (do_inlet_ref .eq. 1) then
            do k = lo(3), hi(3)
               z = xlo(3) + dx(3)*(float(k-lo(3)) + half)
               if (z .le. inlet_ref_height) then
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        tag(i,j,k) = set
                     end do
                  end do
               endif
            end do
         endif

      else if (probtype .eq. 16) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

      else if (probtype .eq. 17) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

      else if (probtype .eq. 18) then

         do k = lo(3), hi(3)
            z = xlo(3) + dx(3)*(float(k-lo(3)) + half)
            if (z.lt.ref_height) then
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
                  end do
               end do
            endif
         end do

      else if (probtype .eq. 19) then

         do k = lo(3), hi(3)
            z = xlo(3) + dx(3)*(float(k-lo(3)) + half)
            if (z.lt.ref_height) then
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
                  end do
               end do
            endif
         end do

      else if (probtype .eq. 20) then

         do k = lo(3), hi(3)
            z = xlo(3) + dx(3)*(float(k-lo(3)) + half)
            if (z.lt.ref_height) then
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
                  end do
               end do
            endif
         end do

      else if (probtype .eq. 21) then

         do k = lo(3), hi(3)
            z = xlo(3) + dx(3)*(float(k-lo(3)) + half)
            if (z.lt.ref_height) then
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
                  end do
               end do
            endif
         end do

      else if (probtype .eq. 22) then

         do k = lo(3), hi(3)
            z = xlo(3) + dx(3)*(float(k-lo(3)) + half)
            if (z.lt.ref_height) then
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
                  end do
               end do
            endif
         end do

      else if (probtype .eq. 23) then

         do k = lo(3), hi(3)
            z = xlo(3) + dx(3)*(float(k-lo(3)) + half)
            if (z.lt.ref_height) then
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
                  end do
               end do
            endif
         end do

      else if (probtype .eq. 24) then

         do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),dx(1)*abs(vort(i,j,k,1)).gt.vorterr)
              end do
           end do
         end do

      else if (probtype .eq. 25) then

         do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),dx(1)*abs(vort(i,j,k,1)).gt.vorterr)
              end do
           end do
         end do
 
      else if (probtype .eq. 32) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do


      else
        print *,'DONT KNOW THIS PROBTYPE IN FORT_MVERROR ',probtype
        stop
      end if

      end subroutine FORT_MVERROR

!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.
!c :::
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data.
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: rho      <=  density array
!c ::: DIMS(rho) => index extent of rho array
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of rho array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
!c ::: -----------------------------------------------------------

      subroutine FORT_DENFILL (rho,DIMS(rho),domlo,domhi,dx, &
                              xlo,time,bc ) &
                              bind(C, name="FORT_DENFILL")
      implicit none

      integer    DIMDEC(rho)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     rho(DIMV(rho))
      integer    bc(SDIM,2)
      integer    lo(SDIM), hi(SDIM)

      integer    i, j, k, n
      REAL_T  umid,rmid,lamv,lamr,rfact, Ly, z0
      REAL_T  hx, hy, hz, gpert, z, y, x, r, ypert, magwif, constn, rs, rd

      parameter (constn=.22089323)

#include <probdata.H>

      ! filcc fills bc_types foextrap, hoextrap, reflect_odd and reflect_even
      call filcc(rho,DIMS(rho),domlo,domhi,dx,xlo,bc)

      lo(1) = ARG_L1(rho)
      lo(2) = ARG_L2(rho)
      lo(3) = ARG_L3(rho)
      hi(1) = ARG_H1(rho)
      hi(2) = ARG_H2(rho)
      hi(3) = ARG_H3(rho)

      if(probtype .eq.8 ) then

      umid = (vel1+vel2)*half
      rmid = (den1+den2)*half
      lamv = (vel1-vel2)/(vel1+vel2)
      lamr = (den1-den2)/(den1+den2)
      hx = dx(1)
      hy = dx(2)
      hz = dx(3)
      magwif = one/ten

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            do k = lo(3), hi(3)
            z = xlo(3) + hz*(float(k-lo(3)) + half)
            do j = lo(2), hi(2)
               y = xlo(2)+hy*(float(j-lo(2))+half)
               ypert = magwif*sin(freq(1)*time)*sin(constn*y)
               rho(i,j,k) =  rmid*(one+lamr*tanh(two*(z-ypert)/delta0))
            enddo
            enddo
         enddo
      endif

      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               rho(i,j,k) = denfact
            enddo
            enddo
         enddo
      endif

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
           do j = lo(2), domlo(2)-1
              do k = lo(3), hi(3)
              do i = lo(1), hi(1)
                 rho(i,j,k) = denfact
              enddo
              enddo
           enddo
      endif

      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            do k = lo(3), hi(3)
            do i = lo(1), hi(1)
               rho(i,j,k) = denfact
            enddo
            enddo
         enddo
      endif

      if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
           do k = lo(3), domlo(3)-1
              do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 if ((probtype.eq.26) .or. (probtype.eq.28)) then
                    rho(i,j,k) = 1.d0
                 else
                    rho(i,j,k) = den2
                 endif
              enddo
              enddo
           enddo
      endif

      if (bc(3,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
           do k = domhi(3)+1, hi(3)
              do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                rho(i,j,k) = den1
              enddo
              enddo
         enddo
      endif

      else ! not probtype 8, all other probtypes


      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(rho).lt.domlo(1)) then
         if (probtype.eq.24) then
!c     SHEAR LAYER
            Ly = f_probhi(2)-f_problo(2)
            do k = ARG_L3(rho), ARG_H3(rho)
               z = xlo(3) + dx(3)*(float(k-lo(3)) + half)
               do j = ARG_L2(rho), ARG_H2(rho)
                  y = xlo(2) + dx(2)*(float(j-lo(2)) + half)
                  z0 = interface_height
                  do n = 1, 10
                     z0 = z0 + mag(n)*cos(freq(n)*time+phi1(n))*cos(float(n)*y/Ly+phi2(n))
                  enddo
                  do i = ARG_L1(rho), domlo(1)-1
!c                     x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                     rho(i,j,k) = half*((den1+den2)+(den1-den2)*tanh(two*(z-z0)/delta0))
                  end do
               end do
            end do
         else
            do i = ARG_L1(rho), domlo(1)-1
               do k = ARG_L3(rho), ARG_H3(rho)
                  do j = ARG_L2(rho), ARG_H2(rho)
                     rho(i,j,k) = denfact
                  end do
               end do
            end do
         end if
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(rho).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(rho)
            do k = ARG_L3(rho), ARG_H3(rho)
               do j = ARG_L2(rho), ARG_H2(rho)
                  rho(i,j,k) = denfact
               end do
	    end do
	 end do
      end if

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(rho).lt.domlo(2)) then
         do j = ARG_L2(rho), domlo(2)-1
            do k = ARG_L3(rho), ARG_H3(rho)
               do i = ARG_L1(rho), ARG_H1(rho)
                  rho(i,j,k) = denfact
               end do
            end do
         end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(rho).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(rho)
            do k = ARG_L3(rho), ARG_H3(rho)
               do i = ARG_L1(rho), ARG_H1(rho)
                  rho(i,j,k) = denfact
               end do
	    end do
	 end do
      end if

      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(rho).lt.domlo(3)) then
         if (probtype.eq.15) then
!c           write (6,*) "Filling with unit density"
           do k = lo(3), domlo(3)-1
              do j = lo(2), hi(2)
                 do i = lo(1), hi(1)
                    rho(i,j,k) = one
                 enddo
              enddo
           enddo
        else if (probtype.eq.16) then
           do k = ARG_L3(rho), domlo(3)-1
              do j = ARG_L2(rho), ARG_H2(rho)
                 do i = ARG_L1(rho), ARG_H1(rho)
                    x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                    rho(i,j,k) = merge(jet_rho,coflow_rho,abs(x-jet_x).le.jet_width)
                 enddo
              enddo
           enddo
        else if (probtype.eq.17) then
           do k = ARG_L3(rho), domlo(3)-1
              do j = ARG_L2(rho), ARG_H2(rho)
                 do i = ARG_L1(rho), ARG_H1(rho)
                    rho(i,j,k) = rhozero
                 enddo
              enddo
           enddo
        else if (probtype.eq.18) then
           rs = jet_rho+coflow_rho
           rd = jet_rho-coflow_rho
           do k = ARG_L3(rho), domlo(3)-1
              do j = ARG_L2(rho), ARG_H2(rho)
                 y = xlo(2) + dx(2)*(float(j-lo(2)) + half)
                 do i = ARG_L1(rho), ARG_H1(rho)
                    x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                    if (plane_jet.eq.1) then
                       r = abs(x-jet_x)
                    else
                       r = sqrt((x-jet_x)*(x-jet_x)+(y-jet_y)*(y-jet_y))
                    endif
                    rho(i,j,k) = half*(rs-rd*tanh(two*(r-jet_width)/delta0))
                 enddo
              enddo
           enddo
        else if (probtype.eq.26) then
           do k = ARG_L3(rho), domlo(3)-1
              do j = ARG_L2(rho), ARG_H2(rho)
                 do i = ARG_L1(rho), ARG_H1(rho)
                    rho(i,j,k) = 1.d0
                 end do
              end do
           end do
        else if (probtype.eq.28) then
           do k = ARG_L3(rho), domlo(3)-1
              do j = ARG_L2(rho), ARG_H2(rho)
                 do i = ARG_L1(rho), ARG_H1(rho)
                    rho(i,j,k) = 1.d0
                 end do
              end do
           end do
        else
           do k = ARG_L3(rho), domlo(3)-1
              do j = ARG_L2(rho), ARG_H2(rho)
                 do i = ARG_L1(rho), ARG_H1(rho)
                    rho(i,j,k) = denfact
                 end do
              end do
           end do
        endif
      end if

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(rho).gt.domhi(3)) then
         if (probtype.eq.17) then
            do k = domhi(3)+1, ARG_H3(rho)
               do j = ARG_L2(rho), ARG_H2(rho)
                  do i = ARG_L1(rho), ARG_H1(rho)
                     rho(i,j,k) = rhozero
                  enddo
               enddo
            enddo
         else
            do k = domhi(3)+1, ARG_H3(rho)
               do j = ARG_L2(rho), ARG_H2(rho)
                  do i = ARG_L1(rho), ARG_H1(rho)
                     rho(i,j,k) = denfact
                  end do
               end do
            end do
         endif
      end if

      endif

      end subroutine FORT_DENFILL

!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.
!c :::
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data.
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: adv      <=  advected quantity array
!c ::: DIMS(adv) => index extent of adv array
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of adv array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
!c ::: -----------------------------------------------------------

      subroutine FORT_ADVFILL (adv,DIMS(adv),domlo,domhi,dx,&
                               xlo,time,bc )&
                               bind(C, name="FORT_ADVFILL")
      implicit none

      integer    DIMDEC(adv)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     adv(DIMV(adv))
      integer    bc(SDIM,2)
      integer    lo(SDIM), hi(SDIM)

      integer    i, j, k, n
      REAL_T  hx, hy, hz, gpert, z, y, x, ypert, magwif, constn,r, Ly, z0

      parameter (constn=.22089323)

#include <probdata.H>

      call filcc(adv,DIMS(adv),domlo,domhi,dx,xlo,bc)

      lo(1) = ARG_L1(adv)
      lo(2) = ARG_L2(adv)
      lo(3) = ARG_L3(adv)
      hi(1) = ARG_H1(adv)
      hi(2) = ARG_H2(adv)
      hi(3) = ARG_H3(adv)

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      if( probtype .eq. 8) then
      magwif = one/ten
      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            do k = lo(3), hi(3)
            z = xlo(3) + hz*(float(k-lo(3)) + half)
            do j = lo(2), hi(2)
               y = xlo(2)+hy*(float(j-lo(2))+half)
               ypert = magwif*sin(freq(1)*time)*sin(constn*y)
               adv(i,j,k) = merge(one,zero,(z-ypert).gt.zero)
            enddo
            enddo
         enddo
      endif

      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               adv(i,j,k) = merge(one,zero,z.gt.zero)
            enddo
            enddo
         enddo
      endif

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            do k = lo(3), hi(3)
            z = xlo(3) + hz*(float(k-lo(3)) + half)
            do i = lo(1), hi(1)
               adv(i,j,k) =  merge(one,zero,z.gt.zero)
            enddo
            enddo
         enddo
      endif

      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            do k = lo(3), hi(3)
            z = xlo(3) + hz*(float(k-lo(3)) + half)
            do i = lo(1), hi(1)
               adv(i,j,k) =  merge(one,zero,z.gt.zero)
            enddo
            enddo
         enddo
      endif

      if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
         do k = lo(3), domlo(3)-1
            do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               adv(i,j,k) = zero
            enddo
            enddo
         enddo
      endif

      if (bc(3,2).eq.EXT_DIR.and.hi(3).gt.domhi(3)) then
         do k = domhi(3)+1, hi(3)
            do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               adv(i,j,k) = one
            enddo
            enddo
         enddo
      endif


      else

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(adv).lt.domlo(1)) then
         if (probtype.eq.24) then
!c     SHEAR LAYER
            do k = ARG_L3(adv), ARG_H3(adv)
               z = xlo(3) + dx(3)*(float(k-lo(3)) + half)
               do j = ARG_L2(adv), ARG_H2(adv)
!c                  y = xlo(2) + dx(2)*(float(j-lo(2)) + half)
                  z0 = interface_height
!c                  do n = 1, 10
!c                     z0 = z0 + mag(n)*cos(freq(n)*time+phi1(n))*cos(float(n)*y/Ly+phi2(n))
!c                  enddo
                  do i = ARG_L1(adv), domlo(1)-1
!c                     x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                     adv(i,j,k) = exp(-((z-z0)/delta0)**2)
                  end do
               end do
            end do
         else
            do i = ARG_L1(adv), domlo(1)-1
               do k = ARG_L3(adv), ARG_H3(adv)
                  do j = ARG_L2(adv), ARG_H2(adv)
                     adv(i,j,k) = zero
                  end do
               end do
            end do
         end if
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(adv).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(adv)
            do k = ARG_L3(adv), ARG_H3(adv)
               do j = ARG_L2(adv), ARG_H2(adv)
                  adv(i,j,k) = zero
               end do
	    end do
	 end do
      end if

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(adv).lt.domlo(2)) then
         do j = ARG_L2(adv), domlo(2)-1
            do k = ARG_L3(adv), ARG_H3(adv)
               do i = ARG_L1(adv), ARG_H1(adv)
                  adv(i,j,k) = zero
               end do
	    end do
	 end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(adv).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(adv)
            do k = ARG_L3(adv), ARG_H3(adv)
               do i = ARG_L1(adv), ARG_H1(adv)
                  adv(i,j,k) = zero
               end do
	    end do
	 end do
      end if

      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(adv).lt.domlo(3)) then
         if (probtype.eq.16) then
            do k = ARG_L3(adv), domlo(3)-1
               do j = ARG_L2(adv), ARG_H2(adv)
                  y = xlo(2) + dx(2)*(float(j-lo(2)) + half)
                  do i = ARG_L1(adv), ARG_H1(adv)
                     x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                     r = sqrt((x-jet_x)*(x-jet_x)+(y-jet_y)*(y-jet_y))
                     adv(i,j,k) = merge(one,zero,r.le.jet_width)
                  enddo
               enddo
            enddo
         else if (probtype.eq.17) then
            do k = ARG_L3(adv), domlo(3)-1
               do j = ARG_L2(adv), ARG_H2(adv)
                  y = xlo(2) + dx(2)*(float(j-lo(2)) + half)
                  do i = ARG_L1(adv), ARG_H1(adv)
                     x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                     if (time.lt.injection_time) then
                        r = sqrt( x*x + y*y )
                        adv(i,j,k) = merge(one,zero,r.le.jet_width)
                     else
                        adv(i,j,k) = zero
                     endif
                  enddo
               enddo
            enddo
         else if (probtype.eq.18) then
            do k = ARG_L3(adv), domlo(3)-1
               do j = ARG_L2(adv), ARG_H2(adv)
                  y = xlo(2) + dx(2)*(float(j-lo(2)) + half)
                  do i = ARG_L1(adv), ARG_H1(adv)
                     x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                     if (plane_jet.eq.1) then
                        r = abs(x-jet_x)
                     else
                        r = sqrt((x-jet_x)*(x-jet_x)+(y-jet_y)*(y-jet_y))
                     endif
                     adv(i,j,k) = half*(one-tanh(two*(r-jet_width)/delta0))
                  enddo
               enddo
            enddo
         else
            do k = ARG_L3(adv), domlo(3)-1
               do j = ARG_L2(adv), ARG_H2(adv)
                  do i = ARG_L1(adv), ARG_H1(adv)
                     adv(i,j,k) = zero
                  end do
               end do
            end do
         end if
      endif

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(adv).gt.domhi(3)) then
         do k = domhi(3)+1, ARG_H3(adv)
            do j = ARG_L2(adv), ARG_H2(adv)
               do i = ARG_L1(adv), ARG_H1(adv)
                  adv(i,j,k) = zero
               end do
            end do
         end do
      end if

      endif

      end subroutine FORT_ADVFILL

      subroutine FORT_ADV2FILL (adv,DIMS(adv),domlo,domhi,dx, &
                                xlo,time,bc )&
                                bind(C, name="FORT_ADV2FILL")

      implicit none

      integer    DIMDEC(adv)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     adv(DIMV(adv))
      integer    bc(SDIM,2)
      integer    lo(SDIM), hi(SDIM)

      integer    i, j, k
      REAL_T  hx, hy, hz, gpert, z, y, x, ypert, magwif, constn,r

      parameter (constn=.22089323)

#include <probdata.H>

      call filcc(adv,DIMS(adv),domlo,domhi,dx,xlo,bc)

      lo(1) = ARG_L1(adv)
      lo(2) = ARG_L2(adv)
      lo(3) = ARG_L3(adv)
      hi(1) = ARG_H1(adv)
      hi(2) = ARG_H2(adv)
      hi(3) = ARG_H3(adv)


      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(adv).lt.domlo(1)) then
         do i = ARG_L1(adv), domlo(1)-1
            do k = ARG_L3(adv), ARG_H3(adv)
               do j = ARG_L2(adv), ARG_H2(adv)
                  adv(i,j,k) = zero
               end do
	    end do
	 end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(adv).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(adv)
            do k = ARG_L3(adv), ARG_H3(adv)
               do j = ARG_L2(adv), ARG_H2(adv)
                  adv(i,j,k) = zero
               end do
	    end do
	 end do
      end if

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(adv).lt.domlo(2)) then
         do j = ARG_L2(adv), domlo(2)-1
            do k = ARG_L3(adv), ARG_H3(adv)
               do i = ARG_L1(adv), ARG_H1(adv)
                  adv(i,j,k) = zero
               end do
	    end do
	 end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(adv).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(adv)
            do k = ARG_L3(adv), ARG_H3(adv)
               do i = ARG_L1(adv), ARG_H1(adv)
                  adv(i,j,k) = zero
               end do
	    end do
	 end do
      end if

      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(adv).lt.domlo(3)) then
         if (probtype.eq.18) then
            do k = ARG_L3(adv), domlo(3)-1
               do j = ARG_L2(adv), ARG_H2(adv)
                  y = xlo(2) + dx(2)*(float(j-lo(2)) + half)
                  do i = ARG_L1(adv), ARG_H1(adv)
                     x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                     if (plane_jet.eq.1) then
                        r = abs(x-jet_x)
                     else
                        r = sqrt((x-jet_x)*(x-jet_x)+(y-jet_y)*(y-jet_y))
                     endif
                     adv(i,j,k) = jet_temp*half*(one-tanh(two*(r-jet_width)/delta0))
                  enddo
               enddo
            enddo
         else
            do k = ARG_L3(adv), domlo(3)-1
               do j = ARG_L2(adv), ARG_H2(adv)
                  do i = ARG_L1(adv), ARG_H1(adv)
                     adv(i,j,k) = zero
                  end do
               end do
            end do
         end if
      endif

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(adv).gt.domhi(3)) then
         if (probtype.eq.17) then
            do k = domhi(3)+1, ARG_H3(adv)
               do j = ARG_L2(adv), ARG_H2(adv)
                  y = xlo(2) + dx(2)*(float(j-lo(2)) + half)
                  do i = ARG_L1(adv), ARG_H1(adv)
                     x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                     if (time.lt.injection_time) then
                        r = sqrt( x*x + y*y )
                        adv(i,j,k) = merge(one,zero,r.le.jet_width)
                     else
                        adv(i,j,k) = zero
                     endif
                  enddo
               enddo
            enddo
         else
            do k = domhi(3)+1, ARG_H3(adv)
               do j = ARG_L2(adv), ARG_H2(adv)
                  do i = ARG_L1(adv), ARG_H1(adv)
                     adv(i,j,k) = zero
                  end do
               end do
            end do
         endif
      end if

      end subroutine FORT_ADV2FILL

!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.
!c :::
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data.
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: temp      <= temperature array
!c ::: DIMS(temp)=> index extent of temp array
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of temp array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
!c ::: -----------------------------------------------------------

      subroutine FORT_TEMPFILL (temp,DIMS(temp),domlo,domhi,dx,&
                                xlo,time,bc )&
                                bind(C, name="FORT_TEMPFILL")

      implicit none

      integer    DIMDEC(temp)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     temp(DIMV(temp))
      integer    bc(SDIM,2)
      integer    lo(SDIM), hi(SDIM)

      integer    i, j, k
      REAL_T     x, y, z, r

#include <probdata.H>

      call filcc(temp,DIMS(temp),domlo,domhi,dx,xlo,bc)

      lo(1) = ARG_L1(temp)
      lo(2) = ARG_L2(temp)
      lo(3) = ARG_L3(temp)
      hi(1) = ARG_H1(temp)
      hi(2) = ARG_H2(temp)
      hi(3) = ARG_H3(temp)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(temp).lt.domlo(1)) then
         do i = ARG_L1(temp), domlo(1)-1
            do k = ARG_L3(temp), ARG_H3(temp)
               do j = ARG_L2(temp), ARG_H2(temp)
                  temp(i,j,k) = one
               end do
	    end do
	 end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(temp).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(temp)
            do k = ARG_L3(temp), ARG_H3(temp)
               do j = ARG_L2(temp), ARG_H2(temp)
                  temp(i,j,k) = one
               end do
	    end do
	 end do
      end if

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(temp).lt.domlo(2)) then
         do j = ARG_L2(temp), domlo(2)-1
            do k = ARG_L3(temp), ARG_H3(temp)
               do i = ARG_L1(temp), ARG_H1(temp)
                  temp(i,j,k) = one
               end do
	    end do
	 end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(temp).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(temp)
            do k = ARG_L3(temp), ARG_H3(temp)
               do i = ARG_L1(temp), ARG_H1(temp)
                  temp(i,j,k) = one
               end do
	    end do
	 end do
      end if

      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(temp).lt.domlo(3)) then
         if (probtype.eq.18) then
            do k = ARG_L3(temp), domlo(3)-1
               do j = ARG_L2(temp), ARG_H2(temp)
                  y = xlo(2) + dx(2)*(float(j-lo(2)) + half)
                  do i = ARG_L1(temp), ARG_H1(temp)
                     x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                     if (plane_jet.eq.1) then
                        r = abs(x-jet_x)
                     else
                        r = sqrt((x-jet_x)*(x-jet_x)+(y-jet_y)*(y-jet_y))
                     endif
                     temp(i,j,k) = merge(one,zero,r.le.jet_width)
                  end do
               end do
            end do
         else
            do k = ARG_L3(temp), domlo(3)-1
               do j = ARG_L2(temp), ARG_H2(temp)
                  do i = ARG_L1(temp), ARG_H1(temp)
                     temp(i,j,k) = one
                  end do
               end do
            end do
         endif
      end if

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(temp).gt.domhi(3)) then
         do k = domhi(3)+1, ARG_H3(temp)
            do j = ARG_L2(temp), ARG_H2(temp)
               do i = ARG_L1(temp), ARG_H1(temp)
                  temp(i,j,k) = one
               end do
            end do
         end do
      end if

      end subroutine FORT_TEMPFILL

!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.
!c :::
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data.
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: u        <=  x velocity array
!c ::: DIMS(u)   => index extent of u array
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of rho array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
!c ::: -----------------------------------------------------------

      subroutine FORT_XVELFILL (u,DIMS(u),domlo,domhi,dx,xlo,time,bc)&
                                bind(C, name="FORT_XVELFILL")
      implicit none
      integer    DIMDEC(u)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     u(DIMV(u))
      integer    bc(SDIM,2)
      integer    lo(SDIM),hi(SDIM)
      integer    i, j, k, n
      integer    isioproc
      REAL_T     x_vel, Ly, z0
      REAL_T     y, ul, g_t
      REAL_T     umid,rmid,lamv,lamr,rfact
      REAL_T     hx, hy, hz, gpert, z,  ypert, magwif, constn
      REAL_T     jv, s, t, xc, yc

      REAL_T  factor

#include <probdata.H>

#ifdef BL_DO_FLCT
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      integer DIMDEC(uflct)
      REAL_T  t_flct
      REAL_T, allocatable :: uflct(:,:,:)
#include <INFL_FORCE_F.H>
#endif

      REAL_T x,r,u1,u2,u3,u_inf,eta

      parameter (constn=.22089323)

      lo(1) = ARG_L1(u)
      lo(2) = ARG_L2(u)
      lo(3) = ARG_L3(u)
      hi(1) = ARG_H1(u)
      hi(2) = ARG_H2(u)
      hi(3) = ARG_H3(u)

#ifdef BL_DO_FLCT
      if (forceInflow) then
         do i = 1, 3
            loFlctArray(i) = lo(i)
            hiFlctArray(i) = hi(i)
         end do
         loFlctArray(adv_dir) = 1
         hiFlctArray(adv_dir) = 1
         call SET_ARGS(DIMS(uflct), loFlctArray, hiFlctArray)
         allocate(uflct(DIMV(uflct)))
!c
!c        Note that we are 'scaling time' here to step into the fluct file to the
!c        correct depth.  This requires that time is not further scaled inside the
!c        the INFL_FILL routine.  Just to be sure, we set convVel = 1 here again.
!c
         convVel = one
         t_flct = time

         call INFL_FILL(FLCT_XVEL,DIMS(uflct),uflct,xlo,dx,t_flct,bc,domnlo,domnhi)
      end if
#endif

      xc = 0.5*(domnhi(1) - domnlo(1))
      yc = 0.5*(domnhi(2) - domnlo(2))

      factor = 1.d0
      if ( (time .gt. 0.d0).and.(tInflowFact_l.ge.0.d0) ) then
         if (time .le. tInflowFact_l) then
            factor = InflowFact_l
         else if (time .ge. tInflowFact_r) then
            factor = InflowFact_r
         else
            factor = InflowFact_l &
                +(time-tInflowFact_l)*(InflowFact_r-InflowFact_l) &
                /(tInflowFact_r-tInflowFact_l)
         endif
      endif

      if (adv_dir .eq. 1) then
         x_vel = adv_vel
      else
         x_vel = zero
      end if

      call filcc(u,DIMS(u),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(u).lt.domlo(1)) then
         if (probtype.eq.24) then
!c     SHEAR LAYER
            Ly = f_probhi(2)-f_problo(2)
            do k = ARG_L3(u), ARG_H3(u)
               z = xlo(3) + dx(3)*(float(k-lo(3)) + half)
               do j = ARG_L2(u), ARG_H2(u)
                  y = xlo(2) + dx(2)*(float(j-lo(2)) + half)
                  z0 = interface_height
                  do n = 1, 10
                     z0 = z0 + mag(n)*cos(freq(n)*time+phi1(n))*cos(float(n)*y/Ly+phi2(n))
                  enddo
                  do i = ARG_L1(u), domlo(1)-1
!c                     x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                     u(i,j,k) = half*((vel1+vel2)+(vel1-vel2)*tanh(two*(z-z0)/delta0))
                  end do
               end do
            end do
         else if (probtype.eq.32) then
            do i = ARG_L1(u), domlo(1)-1
               do k = ARG_L3(u), ARG_H3(u)
                  do j = ARG_L2(u), ARG_H2(u)
                     u(i,j,k) =  x_vel 
                  end do
               end do
            end do
         else
            do i = ARG_L1(u), domlo(1)-1
               do k = ARG_L3(u), ARG_H3(u)
                  do j = ARG_L2(u), ARG_H2(u)
                     u(i,j,k) = x_vel
                  end do
               end do
            end do
         end if
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(u).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(u)
            do k = ARG_L3(u), ARG_H3(u)
               do j = ARG_L2(u), ARG_H2(u)
                  u(i,j,k) = x_vel
               end do
	    end do
	 end do
      end if

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(u).lt.domlo(2)) then
         do j = ARG_L2(u), domlo(2)-1
            do k = ARG_L3(u), ARG_H3(u)
               do i = ARG_L1(u), ARG_H1(u)
                  u(i,j,k) = zero
               end do
	    end do
	 end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(u).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(u)
            do k = ARG_L3(u), ARG_H3(u)
               do i = ARG_L1(u), ARG_H1(u)
                  u(i,j,k) = zero
               end do
	    end do
	 end do
      end if

      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(u).lt.domlo(3)) then
         if (probtype.eq.17) then
            do k = ARG_L3(u), domlo(3)-1
               do j = ARG_L2(u), ARG_H2(u)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(u), ARG_H1(u)
                     x = domnlo(1) + (i+0.5)*dx(1)
                     r = sqrt( x*x + y*y )
                     u(i,j,k) = zero
#ifdef BL_DO_FLCT
                     if (r.lt.jet_width) then
                        u(i,j,k) = u(i,j,k) + uflct(i,j,1)*turb_scale
                     endif
#endif
                  enddo
               enddo
            enddo
         else if (probtype.eq.18) then
            do k = ARG_L3(u), domlo(3)-1
               do j = ARG_L2(u), ARG_H2(u)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(u), ARG_H1(u)
                     x = domnlo(1) + (i+0.5)*dx(1)
                     if (plane_jet.eq.1) then
                        r = abs(x-jet_x)
                     else
                        r = sqrt( (x-jet_x)*(x-jet_x) + (y-jet_y)*(y-jet_y) )
                     endif
                     jv= half*jet_vel*(one-tanh(two*(r-jet_width)/delta0))
                     u(i,j,k) = zero
#ifdef BL_DO_FLCT
                     if (r.lt.jet_width) then
                        u(i,j,k) = u(i,j,k) + jv*uflct(i,j,1)*turb_scale
                     endif
#endif
                  enddo
               enddo
            enddo
         else if (probtype.eq.26) then
            do k = ARG_L3(u), domlo(3)-1
               z = domnlo(3) + (k+0.5)*dx(3)
               do j = ARG_L2(u), ARG_H2(u)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(u), ARG_H1(u)
                     x = domnlo(1) + (i+0.5)*dx(1)
                     u(i,j,k) = zero
#ifdef BL_DO_FLCT
!c     u(i,j,k) = u(i,j,k) + uflct(i,j,1)*turb_scale
                     u(i,j,k) = plateVel(x, y, z, 1, uflct(i,j,1))

#endif
                  enddo
               enddo
            enddo
         else if (probtype.eq.28) then

            do k = ARG_L3(u), domlo(3)-1
               z = domnlo(3) + (k+0.5)*dx(3)
               do j = ARG_L2(u), ARG_H2(u)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(u), ARG_H1(u)
                     x = domnlo(1) + (i+0.5)*dx(1)
                     r = SQRT( (x-xc)*(x-xc) + (y-yc)*(y-yc))
                     u(i,j,k) = zero
#ifdef BL_DO_FLCT
                     if (r .le. 0.025) then
                        u(i,j,k) = uflct(i,j,1)
                     endif

#endif
                  enddo
               enddo
            enddo


         else


            if (probtype.eq.16) then
               t = time*jet_vel/jet_width
            endif

            do k = ARG_L3(u), domlo(3)-1
               do j = ARG_L2(u), ARG_H2(u)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(u), ARG_H1(u)
                     x = domnlo(1) + (i+0.5)*dx(1)

                     if (probtype.eq.13) then
#ifdef BL_DO_FLCT
                        r = SQRT( x*x + y*y )
!c     call vswirlXYZ(x,y,u1,u2,u3)
                        if (r.gt.Rfu) then
                           u_inf = 0.d0
                           eta = TANH(2*(r-Rfu)/Rtran)
                           u(i,j,k) = eta*u_inf + (1.d0-eta)*u1
                           uflct(i,j,1) = 0.d0
                        else
                           u(i,j,k) = u1
                        endif
#endif
                     else if (probtype.eq.15) then
                        u(i,j,k) = plateVel(x,y,z,1,0.d0) * factor
                     else
                        u(i,j,k) = zero
                     endif

#ifdef BL_DO_FLCT
                     if (probtype.eq.16) then
                        x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                        if (t.lt.ten) then
                           s = tenth*(one+fourth*(sin(two*Pi*x/jet_width)+sin(two*Pi*y/jet_width)))
                           jv = jet_vel * (1-exp(-s*t*t))
                        else
                           jv = jet_vel
                        endif
                        if (forceLo .and. adv_dir .eq. 3 .and. abs(x-jet_x).le.jet_width) then
                           u(i,j,k) = u(i,j,k) + uflct(i,j,1)*turb_scale*jv
                        end if
                     else
                        if (forceLo .and. adv_dir .eq. 3) then
                           u(i,j,k) = u(i,j,k) + uflct(i,j,1)*turb_scale
                        endif
                     endif
#endif

                  end do
               end do
            end do

         endif

      end if

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(u).gt.domhi(3)) then
         if (probtype.eq.17) then
            do k = domhi(3)+1, ARG_H3(u)
               do j = ARG_L2(u), ARG_H2(u)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(u), ARG_H1(u)
                     x = domnlo(1) + (i+0.5)*dx(1)
                     r = sqrt( x*x + y*y )
                     u(i,j,k) = zero
#ifdef BL_DO_FLCT
                     if (r.lt.jet_width) then
                        u(i,j,k) = u(i,j,k) + uflct(i,j,1)*turb_scale
                     endif
#endif
                  enddo
               enddo
            enddo
         else if (probtype .eq. 29) then
!c ::: Lid-driven cavity test case, constant velocity on top of domain
            do k = domhi(3)+1, ARG_H3(u)
               do j = ARG_L2(u), ARG_H2(u)
                  do i = ARG_L1(u), ARG_H1(u)
                     u(i,j,k) = lid_vel
                  end do
               end do
            end do
         else
            do k = domhi(3)+1, ARG_H3(u)
               do j = ARG_L2(u), ARG_H2(u)
                  do i = ARG_L1(u), ARG_H1(u)
                     u(i,j,k) = zero
                  end do
               end do
            end do
         endif
      end if



#ifdef BL_DO_FLCT
      if (forceInflow) deallocate(uflct)
#endif

      end subroutine FORT_XVELFILL

!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.
!c :::
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data.
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: v        <=  y velocity array
!c ::: DIMS(v)   => index extent of v array
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of rho array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
!c ::: -----------------------------------------------------------

      subroutine FORT_YVELFILL (v,DIMS(v),domlo,domhi,dx,xlo,time,bc)&
                                bind(C, name="FORT_YVELFILL")

      implicit none
      integer    DIMDEC(v)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     v(DIMV(v))
      integer    bc(SDIM,2)
      integer    lo(SDIM),hi(SDIM)

      integer    i, j, k
      REAL_T     y_vel
      REAL_T     rn
      REAL_T     jv, s, t, xc, yc
      REAL_T factor
#include <probdata.H>

#ifdef BL_DO_FLCT
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      integer DIMDEC(vflct)
      REAL_T  t_flct
      REAL_T, allocatable :: vflct(:,:,:)
#include <INFL_FORCE_F.H>
#endif
      REAL_T x,y,z,r,u1,u2,u3,u_inf,eta

      lo(1) = ARG_L1(v)
      lo(2) = ARG_L2(v)
      lo(3) = ARG_L3(v)
      hi(1) = ARG_H1(v)
      hi(2) = ARG_H2(v)
      hi(3) = ARG_H3(v)

#ifdef BL_DO_FLCT
      if (forceInflow) then
         do i = 1, 3
            loFlctArray(i) = lo(i)
            hiFlctArray(i) = hi(i)
         end do
         loFlctArray(adv_dir) = 1
         hiFlctArray(adv_dir) = 1
         call SET_ARGS(DIMS(vflct), loFlctArray, hiFlctArray)
         allocate(vflct(DIMV(vflct)))
         convVel = one
         t_flct = time
         call INFL_FILL(FLCT_YVEL,DIMS(vflct),vflct,xlo,dx,t_flct,bc,domnlo,domnhi)
      end if
#endif

      xc = 0.5*(domnhi(1) - domnlo(1))
      yc = 0.5*(domnhi(2) - domnlo(2))

      factor = 1.d0
      if ( (time .gt. 0.d0).and.(tInflowFact_l.ge.0.d0) ) then
         if (time .le. tInflowFact_l) then
            factor = InflowFact_l
         else if (time .ge. tInflowFact_r) then
            factor = InflowFact_r
         else
            factor = InflowFact_l &
                +(time-tInflowFact_l)*(InflowFact_r-InflowFact_l) &
                /(tInflowFact_r-tInflowFact_l)
         endif
      endif

      if (adv_dir .eq. 2) then
         y_vel = adv_vel
      else
         y_vel = zero
      end if

      call filcc(v,DIMS(v),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(v).lt.domlo(1)) then
         do i = ARG_L1(v), domlo(1)-1
            do k = ARG_L3(v), ARG_H3(v)
               do j = ARG_L2(v), ARG_H2(v)
                  v(i,j,k) = zero
               end do
	    end do
	 end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(v).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(v)
            do k = ARG_L3(v), ARG_H3(v)
               do j = ARG_L2(v), ARG_H2(v)
                  v(i,j,k) = zero
               end do
            end do
	 end do
      end if

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(v).lt.domlo(2)) then
         if (probtype.eq.32) then
            do j = ARG_L2(v), domlo(2)-1
               do k = ARG_L3(v), ARG_H3(v)
                  do i = ARG_L1(v), ARG_H1(v)
                     z = (xlo(3) + dx(3)*(float(k-lo(3)) + half))/m_probhi(3)
                     v(i,j,k) =  6.0d0 * y_vel * z * (1.0d0 - z)
                  end do
               end do
            end do
         else
            do j = ARG_L2(v), domlo(2)-1
               do k = ARG_L3(v), ARG_H3(v)
                  do i = ARG_L1(v), ARG_H1(v)
                     v(i,j,k) = y_vel
                  end do
               end do
            end do
         end if
      end if

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(v).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(v)
            do k = ARG_L3(v), ARG_H3(v)
               do i = ARG_L1(v), ARG_H1(v)
                  v(i,j,k) = y_vel
               end do
            end do
	 end do
      end if

      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(v).lt.domlo(3)) then
         if (probtype.eq.17) then
            do k = ARG_L3(v), domlo(3)-1
               do j = ARG_L2(v), ARG_H2(v)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(v), ARG_H1(v)
                     x = domnlo(1) + (i+0.5)*dx(1)
                     r = sqrt( x*x + y*y )
                     v(i,j,k) = zero
#ifdef BL_DO_FLCT
                     if (r.lt.jet_width) then
                        v(i,j,k) = v(i,j,k) + vflct(i,j,1)*turb_scale
                     endif
#endif
                  enddo
               enddo
            enddo
         else if (probtype.eq.18) then
            do k = ARG_L3(v), domlo(3)-1
               do j = ARG_L2(v), ARG_H2(v)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(v), ARG_H1(v)
                     x = domnlo(1) + (i+0.5)*dx(1)
                     if (plane_jet.eq.1) then
                        r = abs(x-jet_x)
                     else
                        r = sqrt( (x-jet_x)*(x-jet_x) + (y-jet_y)*(y-jet_y) )
                     endif
                     jv= half*jet_vel*(one-tanh(two*(r-jet_width)/delta0))
                     v(i,j,k) = zero
#ifdef BL_DO_FLCT
                     if (r.lt.jet_width) then
                        v(i,j,k) = v(i,j,k) + jv*vflct(i,j,1)*turb_scale
                     endif
#endif
                  enddo
               enddo
            enddo
         else if (probtype.eq.26) then
            do k = ARG_L3(v), domlo(3)-1
               z = domnlo(3) + (k+0.5)*dx(3)
               do j = ARG_L2(v), ARG_H2(v)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(v), ARG_H1(v)
                     x = domnlo(1) + (i+0.5)*dx(1)
                     v(i,j,k) = zero
#ifdef BL_DO_FLCT
!c     v(i,j,k) = v(i,j,k) + vflct(i,j,1)*turb_scale
                     v(i,j,k) = plateVel(x, y, z, 2, vflct(i,j,1))
#endif
                  enddo
               enddo
            enddo
         else if (probtype.eq.28) then
            do k = ARG_L3(v), domlo(3)-1
               z = domnlo(3) + (k+0.5)*dx(3)
               do j = ARG_L2(v), ARG_H2(v)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(v), ARG_H1(v)
                     x = domnlo(1) + (i+0.5)*dx(1)
                     r = SQRT( (x-xc)*(x-xc) + (y-yc)*(y-yc))
                     v(i,j,k) = zero
#ifdef BL_DO_FLCT
                     if (r .le. 0.025) then
                        v(i,j,k) = vflct(i,j,1)
                     endif
#endif
                  enddo
               enddo
            enddo
         else
            if (probtype.eq.16) then
               t = time*jet_vel/jet_width
            endif

            do k = ARG_L3(v), domlo(3)-1
               do j = ARG_L2(v), ARG_H2(v)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(v), ARG_H1(v)
                     x = domnlo(1) + (i+0.5)*dx(1)

                     if (probtype.eq.13) then
#ifdef BL_DO_FLCT
                        r = SQRT( x*x + y*y )
!c     call vswirlXYZ(x,y,u1,u2,u3)
                        if (r.gt.Rfu) then
                           u_inf = 0.d0
                           eta = TANH(2*(r-Rfu)/Rtran)
                           v(i,j,k) = eta*u_inf + (1.d0-eta)*u2
                           vflct(i,j,1) = 0.d0
                        else
                           v(i,j,k) = u2
                        endif
#endif
                     else if (probtype.eq.15) then
                        v(i,j,k) = plateVel(x,y,z,2,0.d0) * factor
                     else
                        v(i,j,k) = zero
                     endif

#ifdef BL_DO_FLCT
                     if (probtype.eq.16) then
                        x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                        if (t.lt.ten) then
                           s = tenth*(one+fourth*(sin(two*Pi*x/jet_width)+sin(two*Pi*y/jet_width)))
                           jv = jet_vel * (1-exp(-s*t*t))
                        else
                           jv = jet_vel
                        endif
                        if (forceLo .and. adv_dir .eq. 3 .and. abs(x-jet_x).le.jet_width) then
                           v(i,j,k) = v(i,j,k) + vflct(i,j,1)*turb_scale*jv
                        end if
                     else
                        if (forceLo .and. adv_dir .eq. 3) then
                           v(i,j,k) = v(i,j,k) + vflct(i,j,1)*turb_scale
                        end if
                     endif
#endif

                  end do
               end do
            end do
         endif
      end if

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(v).gt.domhi(3)) then
         if (probtype.eq.17) then
            do k = domhi(3)+1, ARG_H3(v)
               do j = ARG_L2(v), ARG_H2(v)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(v), ARG_H1(v)
                     x = domnlo(1) + (i+0.5)*dx(1)
                     r = sqrt( x*x + y*y )
                     v(i,j,k) = zero
#ifdef BL_DO_FLCT
                     if (r.lt.jet_width) then
                        v(i,j,k) = v(i,j,k) + vflct(i,j,1)*turb_scale
                     endif
#endif
                  enddo
               enddo
            enddo
         else
            do k = domhi(3)+1, ARG_H3(v)
               do j = ARG_L2(v), ARG_H2(v)
                  do i = ARG_L1(v), ARG_H1(v)
                     v(i,j,k) = zero
                  end do
               end do
            end do
         endif
      end if

#ifdef BL_DO_FLCT
      if (forceInflow) deallocate(vflct)
#endif

      end subroutine FORT_YVELFILL

!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.
!c :::
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data.
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: w        <=  z velocity array
!c ::: DIMS(w)   => index extent of v array
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of rho array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
!c ::: -----------------------------------------------------------

      subroutine FORT_ZVELFILL (w,DIMS(w),domlo,domhi,dx,xlo,time,bc)&
                                bind(C, name="FORT_ZVELFILL")

      implicit none
      integer    DIMDEC(w)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     w(DIMV(w))
      integer    bc(SDIM,2)
      integer    lo(SDIM),hi(SDIM)
      integer    i, j, k
      REAL_T     z_vel
      REAL_T     jv, s, t,xc,yc

#include <probdata.H>

#ifdef BL_DO_FLCT
      REAL_T t_flct
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      integer DIMDEC(wflct)
      REAL_T, allocatable :: wflct(:,:,:)
#include <INFL_FORCE_F.H>
#endif
      REAL_T x,y,z,r,u1,u2,u3,u_inf,eta
      REAL_T Lx, Ly, pert, factor

      lo(1) = ARG_L1(w)
      lo(2) = ARG_L2(w)
      lo(3) = ARG_L3(w)
      hi(1) = ARG_H1(w)
      hi(2) = ARG_H2(w)
      hi(3) = ARG_H3(w)

#ifdef BL_DO_FLCT
      if (forceInflow) then
         do i = 1, 3
            loFlctArray(i) = lo(i)
            hiFlctArray(i) = hi(i)
         end do
         loFlctArray(adv_dir) = 1
         hiFlctArray(adv_dir) = 1
         call SET_ARGS(DIMS(wflct), loFlctArray, hiFlctArray)
         allocate(wflct(DIMV(wflct)))
         convVel = one
         t_flct = time
         call INFL_FILL(FLCT_ZVEL,DIMS(wflct),wflct,xlo,dx,t_flct,bc,domnlo,domnhi)
      end if
#endif

      xc = 0.5*(domnhi(1) - domnlo(1))
      yc = 0.5*(domnhi(2) - domnlo(2))

      if ( (time .gt. zero).and.(tVco_l.ge.zero) ) then
         if (time .le. tVco_l) then
            Vco = Vco_l
         else if (time .ge. tVco_r) then
            Vco = Vco_r
         else
            Vco = Vco_l+(time-tVco_l)*(Vco_r-Vco_l)/(tVco_r-tVco_l)
         endif
      endif

      factor = 1.d0
      if ( (time .gt. 0.d0).and.(tInflowFact_l.ge.0.d0) ) then
         if (time .le. tInflowFact_l) then
            factor = InflowFact_l
         else if (time .ge. tInflowFact_r) then
            factor = InflowFact_r
         else
            factor = InflowFact_l &
                +(time-tInflowFact_l)*(InflowFact_r-InflowFact_l) &
                /(tInflowFact_r-tInflowFact_l)
         endif
      endif

      if (adv_dir .eq. 3) then
         z_vel = adv_vel
      else
         z_vel = zero
      end if

      call filcc(w,DIMS(w),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(w).lt.domlo(1)) then
         do i = ARG_L1(w), domlo(1)-1
            do k = ARG_L3(w), ARG_H3(w)
               do j = ARG_L2(w), ARG_H2(w)
                  w(i,j,k) = zero
               end do
	    end do
	 end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(w).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(w)
            do k = ARG_L3(w), ARG_H3(w)
               do j = ARG_L2(w), ARG_H2(w)
                  w(i,j,k) = zero
               end do
	    end do
	 end do
      end if

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(w).lt.domlo(2)) then
         do j = ARG_L2(w), domlo(2)-1
            do k = ARG_L3(w), ARG_H3(w)
               do i = ARG_L1(w), ARG_H1(w)
                  w(i,j,k) = zero
               end do
            end do
	 end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(w).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(w)
            do k = ARG_L3(w), ARG_H3(w)
               do i = ARG_L1(w), ARG_H1(w)
                  w(i,j,k) = zero
               end do
            end do
	 end do
      end if

      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(w).lt.domlo(3)) then
         if (probtype.eq.17) then
            do k = ARG_L3(w), domlo(3)-1
               do j = ARG_L2(w), ARG_H2(w)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(w), ARG_H1(w)
                     x = domnlo(1) + (i+0.5)*dx(1)
                     r = sqrt( x*x + y*y )
                     if (r.lt.jet_width.and.time.lt.injection_time) then
                        w(i,j,k) = jet_vel
                     else
                        w(i,j,k) = zero
                     endif
#ifdef BL_DO_FLCT
                     if (r.lt.jet_width) then
                        w(i,j,k) = w(i,j,k) + wflct(i,j,1)*turb_scale
                     endif
#endif
                  enddo
               enddo
            enddo
         else if (probtype.eq.32) then
            do k = ARG_L3(w), domlo(3)-1
               do j = ARG_L2(w), ARG_H2(w)
                  do i = ARG_L1(w), ARG_H1(w)
                     x = (xlo(1) + dx(1)*(float(i-lo(1)) + half))/m_probhi(1)
                     w(i,j,k) = 6.0d0 * z_vel * x * (1.0d0 - x)
                  end do
               end do
            end do

         else if (probtype.eq.18) then
            do k = ARG_L3(w), domlo(3)-1
               do j = ARG_L2(w), ARG_H2(w)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(w), ARG_H1(w)
                     x = domnlo(1) + (i+0.5)*dx(1)
                     if (plane_jet.eq.1) then
                        r = abs(x-jet_x)
                     else
                        r = sqrt( (x-jet_x)*(x-jet_x) + (y-jet_y)*(y-jet_y) )
                     endif
                     w(i,j,k) = half*((jet_vel+coflow_vel)-(jet_vel-coflow_vel)*tanh(two*(r-jet_width)/delta0))
#ifdef BL_DO_FLCT
                     if (r.lt.jet_width) then
                        w(i,j,k) = w(i,j,k)*(one+wflct(i,j,1)*turb_scale)
                     endif
#endif
                  enddo
               enddo
            enddo

         else if (probtype.eq.26) then
            do k = ARG_L3(w), domlo(3)-1
               z = domnlo(3) + (k+0.5)*dx(3)
               do j = ARG_L2(w), ARG_H2(w)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(w), ARG_H1(w)
                     x = domnlo(1) + (i+0.5)*dx(1)
                     w(i,j,k) = meanPlateVel(x,y)
#ifdef BL_DO_FLCT
                     w(i,j,k) = plateVel(x, y, z, 3, wflct(i,j,1))
!c                     w(i,j,k) = w(i,j,k)*(one+wflct(i,j,1)*turb_scale)
#endif
                  enddo
               enddo
            enddo
         else if (probtype.eq.28) then
            do k = ARG_L3(w), domlo(3)-1
               z = domnlo(3) + (k+0.5)*dx(3)
               do j = ARG_L2(w), ARG_H2(w)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(w), ARG_H1(w)
                     x = domnlo(1) + (i+0.5)*dx(1)
                     r = SQRT( (x-xc)*(x-xc) + (y-yc)*(y-yc))
                     w(i,j,k) = coflow_vel
#ifdef BL_DO_FLCT
                     if( r.le. 0.025) then
                       	w(i,j,k) = wflct(i,j,1)
                     endif
#endif
                  enddo
               enddo
            enddo



         else
            if (probtype.eq.16) then
               t = time*jet_vel/jet_width
            endif

            do k = ARG_L3(w), domlo(3)-1
               do j = ARG_L2(w), ARG_H2(w)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(w), ARG_H1(w)
                     x = domnlo(1) + (i+0.5)*dx(1)

                     if (probtype.eq.13) then
#ifdef BL_DO_FLCT
                        r = SQRT( x*x + y*y )
!c     call vswirlXYZ(x,y,u1,u2,u3)
                        if (r.gt.Rfu) then
                           u_inf = Vco
                           eta = TANH(2*(r-Rfu)/Rtran)
                           w(i,j,k) = eta*u_inf + (1.d0-eta)*u3
                           wflct(i,j,1) = 0.d0
                        else
                           w(i,j,k) = u3
                        endif
#endif
                     else if (probtype.eq.16) then
                        x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                        if (t.lt.ten) then
                           s = tenth*(one+fourth*(sin(two*Pi*x/jet_width)+sin(two*Pi*y/jet_width)))
                           jv = jet_vel * (1-exp(-s*t*t))
                        else
                           jv = jet_vel
                        endif
                        if (abs(x-jet_x).lt.jet_width) then
                           w(i,j,k) = jv * half * (1 - tanh((abs(x-jet_x)-jet_width)/delta0))
                        else
                           w(i,j,k) = zero
                        endif
                     else if (probtype.eq.15) then
                        w(i,j,k) = plateVel(x,y,z,3,0.d0) * factor

                        Lx = domnhi(1)-domnlo(1)
                        Ly = domnhi(2)-domnlo(2)
                        pert = .03*(sin(2*Pi*4*x/Lx) * sin(2*Pi*3*y/Ly) ) * sin(time/.00001) + &
                            .0438*(sin(2*Pi*5*(x-.17)/Lx) * sin(2*Pi*6*(y-.49)/Ly)) * sin(time/.000012)
                        w(i,j,k) = w(i,j,k)*(1.d0 + pert)
                     else
                        w(i,j,k) = z_vel
                     endif

#ifdef BL_DO_FLCT
                     if (probtype.eq.16) then
                        x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                        if (forceLo .and. adv_dir .eq. 3 .and. abs(x-jet_x).le.jet_width) then
                           w(i,j,k) = w(i,j,k) + wflct(i,j,1)*turb_scale*jv
                        endif
                     else
                        if (forceLo .and. adv_dir .eq. 3) then
                           w(i,j,k) = w(i,j,k) + wflct(i,j,1)*turb_scale
                        endif
                     endif
#endif

                  end do
               end do
            end do
         endif
      end if

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(w).gt.domhi(3)) then
         if (probtype.eq.17) then
            do k = domhi(3)+1, ARG_H3(w)
               do j = ARG_L2(w), ARG_H2(w)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(w), ARG_H1(w)
                     x = domnlo(1) + (i+0.5)*dx(1)
                     r = sqrt( x*x + y*y )
                     if (r.lt.jet_width.and.time.lt.injection_time) then
                        w(i,j,k) = -jet_vel
                     else
                        w(i,j,k) = zero
                     endif
#ifdef BL_DO_FLCT
                     if (r.lt.jet_width) then
                        w(i,j,k) = w(i,j,k) + wflct(i,j,1)*turb_scale
                     endif
#endif
                  enddo
               enddo
            enddo
         else
            do k = domhi(3)+1, ARG_H3(w)
               do j = ARG_L2(w), ARG_H2(w)
                  do i = ARG_L1(w), ARG_H1(w)
                     w(i,j,k) = z_vel
                  end do
               end do
            end do
         endif
      end if

#ifdef BL_DO_FLCT
      if (forceInflow) deallocate(wflct)
#endif

      end subroutine FORT_ZVELFILL

      REAL_T function meanPlateVel(x,y)
      implicit none
#include <probdata.H>
      REAL_T x,y
      REAL_T totHoleA, holeSepX, holeSepY, Lx, Ly
      REAL_T jetVel, xHole, yHole, rHole, rFact
      REAL_T hRad, hBL
      integer iHoleX, iHoleY

      Lx = (domnhi(1) - domnlo(1))
      Ly = (domnhi(2) - domnlo(2))
      holeSepX = Lx / nHolesX
      holeSepY = Ly / nHolesY
      hRad = MIN(holeRad, 0.3*MIN(holeSepX,holeSepY))
      hBL = hRad/holeBLfac

      totHoleA = nHolesX*nHolesY*Pi*hRad**2
      jetVel = adv_vel * Lx * Ly / totHoleA

      iHoleX = INT(x/holeSepX)
      iHoleY = INT(y/holeSepY)

      xHole = x - (iHoleX + 0.5)*holeSepX
      yHole = y - (iHoleY + 0.5)*holeSepY
      rHole = SQRT(xHole**2 + yHole**2)

      rFact = 1.d0 - zblend1(rHole,hRad,hBL)
      meanPlateVel = jetVel * rFact
      end function meanPlateVel

!c *****
!c **
!c *****
      REAL_T function plateVel(x,y,z,cord,flct)
      implicit none

      REAL_T   flct
      integer dir(SDIM), numholes(SDIM), idholes(SDIM)
      REAL_T  holespace(SDIM)
      REAL_T  dist(SDIM), loc(SDIM)
!c
!c     ::::: local variables
!c
      integer cord, n
      REAL_T  twicePi
      REAL_T totHoleA
      REAL_T jetVel, radius, sectarea
      REAL_T x, y, z, xcen,ycen, xloc, yloc, xcord(5), ycord(5)
      REAL_T mean, vholeSp, vane, xshift, yshift, delta,vslot
      integer icell,jcell,npts

      integer NringMAX, NperMAX
      parameter (NringMAX=30,NperMAX=300)
      integer NperRing(NringMAX)
      double precision centHole(NringMAX,NperMAX,2)
      double precision ringRad(NringMAX), cent(2)
      integer i,j, Nring
      double precision theta, ringSp


#include <probdata.H>
      mean = 0.40
      dir(1) = 0
      dir(2) = 0
      dir(3) = 0
      vholeSp = 0.d0
      vane = 0.001
      delta = 0.d0
      vslot = 32.d0


      do  npts = 1, 5
         xcord(npts) = 0
         ycord(npts) = 0
      end do

      if (adv_dir.eq.1) then
         dir(1) = 1
      else if (adv_dir.eq.2) then
         dir(2) = 1
      else
         dir(3) = 1
      endif
!c     perforated plate of square pattern
#if 0

      radius = 0.d0
      jetVel = 0.d0

      twicePi = two*Pi

      numholes(1)= nHolesX
      numholes(2)= nHolesY
      numholes(3)= nHolesZ

      totHoleA = numholes(1)*numholes(2)*numholes(3)*Pi*holeRad**2


      sectarea = dir(1)*(domnhi(2) - domnlo(2))*(domnhi(3) - domnlo(3))+ &
          dir(2)*(domnhi(1) - domnlo(1)) * (domnhi(3) - domnlo(3))+ &
          dir(3)*(domnhi(1) - domnlo(1)) * (domnhi(2) - domnlo(2))

      do n = 1, 3
         holespace(n) = (domnhi(n) - domnlo(n)) / numholes(n)
      end do

      jetVel = adv_vel*sectarea/totHoleA

      loc(3) = z
      loc(2) = y
      loc(1) = x

      do n = 1, 3
         idholes(n) = INT(loc(n)/holespace(n))
         dist(n) = loc(n) - (idholes(n) + half)*holespace(n)
      end do

      radius = SQRT(dir(1)*(dist(2)**2 + dist(3)**2)+ &
          dir(2)*(dist(1)**2 + dist(3)**2) + &
          dir(3)*(dist(1)**2 + dist(2)**2))

      if (radius .le. holeRad) then
         plateVel = jetVel*dir(cord)+flct*turb_scale
      else
         plateVel = 0.d0
      endif
#endif
!c     perforated plate of hexagnoal pattern
#if 0

      jetVel = adv_vel/(1.0 - holeBLfac)
      holeSp = SQRT(two*Pi/(sqrt(3.d0)*(1.d0 - holeBLfac)))*holeRad

      vholeSp = sqrt(3.d0)*holeSp

      plateVel = 0.d0

      icell = int (x/holeSp)
      jcell = int (y/vholeSp)

      xcen = holeSp*(dfloat(icell) + half)
      ycen = vholeSp*(dfloat(jcell) + half)

      xloc = x-xcen
      yloc = y-ycen

      xcord(1) = -half*holeSp
      ycord(1) = -half*vholeSp

      xcord(2) = -half*holeSp
      ycord(2) = half*vholeSp

      xcord(3) = half*holeSp
      ycord(3) = half*vholeSp

      xcord(4) = half*holeSp
      ycord(4) = -half*vholeSp

      xcord(5) = 0.d0
      ycord(5) = 0.d0

      do npts = 1, 5

         radius = SQRT((xloc-xcord(npts))**2 + (yloc-ycord(npts))**2)

         if(radius .le. holeRad) then
            plateVel = jetVel*dir(cord) + flct*turb_scale
            plateVel = plateVel *(1.d0 - 2.d0*zblend1(radius,holeRad,holeRad/8.0))
         endif
      end do


#endif
!c     perforated plate of hexagnoal pattern, no flow on boundaries (to enable periodic ok)
#if 0

      jetVel = adv_vel/(1.0 - holeBLfac)
      holeSp = SQRT(two*Pi/(sqrt(3.d0)*(1.d0 - holeBLfac)))*holeRad

      vholeSp = sqrt(3.d0)*holeSp

      plateVel = 0.d0

      icell = int (x/holeSp)
      jcell = int (y/vholeSp)

      xcen = holeSp*(dfloat(icell) + half)
      ycen = vholeSp*(dfloat(jcell) + half)

      xloc = x-xcen
      yloc = y-ycen

      xcord(1) = -half*holeSp
      ycord(1) = -half*vholeSp

      xcord(2) = -half*holeSp
      ycord(2) = half*vholeSp

      xcord(3) = half*holeSp
      ycord(3) = half*vholeSp

      xcord(4) = half*holeSp
      ycord(4) = -half*vholeSp

      xcord(5) = 0.d0
      ycord(5) = 0.d0

      do npts = 1, 5

         radius = SQRT((xloc-xcord(npts))**2 + (yloc-ycord(npts))**2)

         if(radius .le. holeRad) then
            plateVel = jetVel*dir(cord) + flct*turb_scale
!c            plateVel = plateVel *(1.d0 - 2.d0*zblend1(radius,holeRad,holeRad/8.0))
         endif
      end do

      if ((x.lt.domnlo(1) + .001) .or.  &
          (x.gt.domnhi(1) - .001) .or. &
          (y.lt.domnlo(2) + .001) .or.  &
          (y.gt.domnhi(2) - .001)) then
         plateVel = 0.d0
      endif

#endif

!c     perforated plate with radial pattern, no flow on boundaries (to enable periodic ok)
#if 1
      jetVel = adv_vel/(1.0 - holeBLfac)

! FIXME : from compiler
!PROB_3D.F90:7106:0:
!
!          cent(n) = 0.5d0 * (domnlo(n) + domnhi(n))
!
!Warning: iteration 2 invokes undefined behavior [-Waggressive-loop-optimizations]
!PROB_3D.F90:7105:0:
!
!       do n=1,3
!
!note: within this loop
!
      do n=1,3
         cent(n) = 0.5d0 * (domnlo(n) + domnhi(n))
      enddo
      holeRad = 0.0036d0
      ringSp = .0112d0
      Nring = MIN(NringMAX,1 + INT(0.5d0*sqrt(2.d0)*MAX(domnlo(1)+domnhi(1),domnlo(2)+domnhi(2))/ringSp))

      do i=1,Nring
         ringRad(i) = ringSp * (i-1)
         if (i.eq.1) then
            NperRing(i) = 1
            centHole(i,1,1) = cent(1)
            centHole(i,1,2) = cent(2)
         else
            NperRing(i) = 6 * (i-1)
            do j=1,NperRing(i)
               if (j .le. NperMAX) then
                  theta = j*4*ASIN(1.d0)/DBLE(NperRing(i))
                  centHole(i,j,1) = cent(1) + ringRad(i) * SIN(theta)
                  centHole(i,j,2) = cent(2) + ringRad(i) * COS(theta)
               endif
            enddo
         endif
      enddo

      plateVel = 0.d0
      do i=1,Nring
         do j=1,NperRing(i)
            if (j .le. NperMAX) then
               radius = sqrt( (x-centHole(i,j,1))**2 + (y-centHole(i,j,2))**2 )
               if(radius .le. holeRad) then
                  plateVel = jetVel*dir(cord) + flct*turb_scale
               endif
            endif
         enddo
      enddo

      if ((x.lt.domnlo(1) + .001) .or. &
          (x.gt.domnhi(1) - .001) .or. &
          (y.lt.domnlo(2) + .001) .or.  &
          (y.gt.domnhi(2) - .001)) then
         plateVel = 0.d0
      endif

#endif
!c  this is the original code
#if 0

      jetVel = adv_vel/(1.0 - holeBLfac)
      holeSp = SQRT(two*Pi/(sqrt(3.d0)*(1.d0 - holeBLfac)))*holeRad

      if(y.ge.0.d0)then
         jcell = int(y / (half*sqrt(3.d0)*holeSp))
      else
         jcell = -1-int(-y/(half*sqrt(3.d0)*holeSp))
      endif
      yloc = y-dfloat(jcell)*half*sqrt(3.d0)*holeSp
      x = x-y/sqrt(3.d0)
      if(x.ge.0.d0) then
         icell = int(x / (holeSp))
      else
         icell = -1-int(-x/(holeSp))
      endif

      xloc = x - holeSp * dfloat(icell) + yloc/sqrt(3.d0)
      xcen =.75d0*holeSp
      ycen =(sqrt(3.d0)*.25d0)*holeSp
      radius = SQRT((xloc-xcen)**2 + (yloc-ycen)**2)
      if(radius .le. holeRad) then
         plateVel = jetVel*dir(cord) + flct*turb_scale
!c         plateVel = plateVel *(1.d0 - 2.d0*zblend1(radius,holeRad,holeRad/8.0))
!c      else
!c         plateVel = -jetVel*dir(cord)
!c         plateVel = 0.d0
      endif

#endif
!c     perforated plate of hexagnoal pattern with slot  (2 gaps in swirl region)
#if 0

!c     jetVel = adv_vel/(1.0 - holeBLfac)
!c
      jetVel = 18.d0
      plateVel = 0.d0

      holeSp = SQRT(two*Pi/(sqrt(3.d0)*(1.d0 - holeBLfac)))*holeRad
      vholeSp = sqrt(3.d0)*holeSp

      icell = int (x/holeSp)
      jcell = int (y/vholeSp)

      xcen = holeSp*(dfloat(icell) + half)
      ycen = vholeSp*(dfloat(jcell) + half)

      xloc = x-xcen
      yloc = y-ycen

      xcord(1) = -half*holeSp
      ycord(1) = -half*vholeSp

      xcord(2) = -half*holeSp
      ycord(2) = half*vholeSp

      xcord(3) = half*holeSp
      ycord(3) = half*vholeSp

      xcord(4) = half*holeSp
      ycord(4) = -half*vholeSp

      xcord(5) = 0.d0
      ycord(5) = 0.d0

      do npts = 1, 5

         radius = SQRT((xloc-xcord(npts))**2 + (yloc-ycord(npts))**2)

         if(radius .le. holeRad) then
            plateVel = jetVel*dir(cord) + flct*turb_scale
!c            plateVel = plateVel *(1.d0 - 2.d0*zblend1(radius,holeRad,holeRad/8.0))
         endif
      end do
!c     in the slot region (5 mm width)
      if (y.ge.((domnhi(2)-domnlo(2)) - slotWidth)) then
         if (cord .eq. 1) then
            plateVel = 32.d0*tan(37.d0/180.d0*Pi)

            if (x .ge. 0.25*(domnhi(1)-domnlo(1)) .and. (x .le. (0.25*(domnhi(1)-domnlo(1)) + vane))) then
               plateVel = 0.d0
            else  if (x .ge. 0.75*(domnhi(1)-domnlo(1)) .and. (x .le. (0.75*(domnhi(1)-domnlo(1)) + vane))) then
               plateVel = 0.d0
            end if

         else if (cord .eq. 2) then
            plateVel = 0.d0
         else
            if (x .ge. 0.25*(domnhi(1)-domnlo(1)) .and. (x .le. (0.25*(domnhi(1)-domnlo(1)) + vane))) then
               plateVel = 0.d0
            else  if (x .ge. 0.75*(domnhi(1)-domnlo(1)) .and. (x .le. (0.75*(domnhi(1)-domnlo(1)) + vane))) then
               plateVel = 0.d0
            else
               plateVel = 32.d0
            endif
         endif
      endif
      if ((y.ge.((domnhi(2)-domnlo(2)) - slotWidth - vane)) .and. &
          (y.lt.((domnhi(2)-domnlo(2)) - slotWidth))) then
         plateVel = 0.d0
      endif
!c     in the region of 1 mm between the slot region and the plate

#endif


!c     perforated plate of hexagnoal pattern with slot and shifting in x and y directions (4 gaps in swirl region)
#if 0

!c     jetVel = adv_vel/(1.0 - holeBLfac)
!c

      jetVel = adv_vel
      vslot = slot_vel

      plateVel = 0.d0


      holeSp = SQRT(two*Pi/(sqrt(3.d0)*(1.d0 - holeBLfac)))*holeRad
      vholeSp = sqrt(3.d0)*holeSp

      xshift = alpha*holeSp
      yshift = beta*vholeSp

      icell = int ((x-xshift)/holeSp)
      jcell = int ((y-yshift)/vholeSp)

      if ((x-xshift) .le. 0.d0) then
         icell = - 1 - int ((x-xshift)/holeSp)
      endif
      if ((y-yshift) .le. 0.d0) then
         jcell = - 1 - int ((y-yshift)/vholeSp)
      endif

      xcen = holeSp*(dfloat(icell) + half) + xshift
      ycen = vholeSp*(dfloat(jcell) + half) + yshift

      xloc = x-xcen
      yloc = y-ycen

      xcord(1) = -half*holeSp
      ycord(1) = -half*vholeSp

      xcord(2) = -half*holeSp
      ycord(2) = half*vholeSp

      xcord(3) = half*holeSp
      ycord(3) = half*vholeSp

      xcord(4) = half*holeSp
      ycord(4) = -half*vholeSp

      xcord(5) = 0.d0
      ycord(5) = 0.d0

      do npts = 1, 5

         radius = SQRT((xloc-xcord(npts))**2 + (yloc-ycord(npts))**2)

         if(radius .le. holeRad) then
            plateVel = jetVel*dir(cord) + flct*turb_scale
!c            plateVel = plateVel *(1.d0 - 2.d0*zblend1(radius,holeRad,holeRad/8.0))
         endif
      end do
!c     in the slot region (5 mm width)
      if (y.ge.((domnhi(2)-domnlo(2)) - slotWidth)) then
         delta = (domnhi(1)-domnlo(1))
         if (cord .eq. 1) then
!c            plateVel = 32.d0*tan(37.d0/180.d0*Pi)
            plateVel = vslot*tan(37.d0/180.d0*Pi)

            do n = 1, 7, 2
               if (x .ge. (dfloat(n)/8.d0*delta-half*vane) .and. (x .le. (dfloat(n)/8.d0*delta + half*vane))) then
                  plateVel = 0.d0
               end if
            end do

         else if (cord .eq. 2) then
            plateVel = 0.d0
         else
!c            plateVel = 32.d0
            plateVel = vslot
            do n = 1, 7, 2
               if (x .ge. (dfloat(n)/8.d0*delta-half*vane) .and. (x .le. (dfloat(n)/8.d0*delta + half*vane))) then
                  plateVel = 0.d0
               end if
            end do
         endif
      endif
!c     in the thin 1 mm wide no flow region
      if ((y.ge.((domnhi(2)-domnlo(2)) - slotWidth - vane)) .and. &
          (y.lt.((domnhi(2)-domnlo(2)) - slotWidth))) then
         plateVel = 0.d0
      endif


#endif



      end function PlateVel



!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.
!c :::
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data.
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: u        <=  full velocity array
!c ::: DIMS(u)   => index extent of v array
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of rho array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi,comp)
!c ::: -----------------------------------------------------------

      subroutine FORT_VELFILL (u,DIMS(u),domlo,domhi,dx,xlo,time,bc) &
                               bind(C, name="FORT_VELFILL")
      implicit none
      integer    DIMDEC(u)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     u(DIMV(u),SDIM)
      integer    bc(SDIM,2,SDIM)
      integer    lo(SDIM),hi(SDIM)
      integer    i, j, k
      integer    isioproc
      integer    ii, ji, ip, jp
      REAL_T     x, y, z, r, jv
      REAL_T     xi, yi, xp, yp, rp
      REAL_T     uxp, uyp, up, ax, r2

      integer    ilo, ihi, jlo, jhi

#include <probdata.H>

#ifdef BL_DO_FLCT
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      integer DIMDEC(uflct)
      REAL_T  t_flct
      REAL_T, allocatable :: uflct(:,:,:)
      REAL_T, allocatable :: vflct(:,:,:)
      REAL_T, allocatable :: wflct(:,:,:)
#include <INFL_FORCE_F.H>
#endif

      lo(1) = ARG_L1(u)
      lo(2) = ARG_L2(u)
      lo(3) = ARG_L3(u)
      hi(1) = ARG_H1(u)
      hi(2) = ARG_H2(u)
      hi(3) = ARG_H3(u)

#ifdef BL_DO_FLCT
      if (forceInflow) then
         do i = 1, 3
            loFlctArray(i) = lo(i)
            hiFlctArray(i) = hi(i)
         end do
         loFlctArray(adv_dir) = 1
         hiFlctArray(adv_dir) = 1
         call SET_ARGS(DIMS(uflct), loFlctArray, hiFlctArray)
         allocate(uflct(DIMV(uflct)))
         allocate(vflct(DIMV(uflct)))
         allocate(wflct(DIMV(uflct)))
!c
!c        Note that we are 'scaling time' here to step into the fluct file to the
!c        correct depth.  This requires that time is not further scaled inside the
!c        the INFL_FILL routine.  Just to be sure, we set convVel = 1 here again.
!c
         convVel = one
         t_flct = adv_vel*time
         call INFL_FILL(FLCT_XVEL,DIMS(uflct),uflct,xlo,dx,t_flct,bc(1,1,1),domnlo,domnhi)
         call INFL_FILL(FLCT_XVEL,DIMS(uflct),vflct,xlo,dx,t_flct,bc(1,1,2),domnlo,domnhi)
         call INFL_FILL(FLCT_XVEL,DIMS(uflct),wflct,xlo,dx,t_flct,bc(1,1,3),domnlo,domnhi)
      end if
#endif

      call filcc (u(ARG_L1(u),ARG_L2(u),ARG_L3(u),1), &
                 DIMS(u),domlo,domhi,dx,xlo,bc(1,1,1))
      call filcc (u(ARG_L1(u),ARG_L2(u),ARG_L3(u),2),&
                 DIMS(u),domlo,domhi,dx,xlo,bc(1,1,2))
      call filcc (u(ARG_L1(u),ARG_L2(u),ARG_L3(u),3),&
                 DIMS(u),domlo,domhi,dx,xlo,bc(1,1,3))

!c     This forces radial inflow in a plane satisfying potential flow/incompressiblity

      call bl_pd_is_ioproc(isioproc)

!c      if (isioproc.eq.1) then
!c         write (*,*) "In VELFILL..."
!c      endif

      if (probtype.eq.18) then
!c     This is the Jet probtype
!c     Here, we assume all the velocity components have the same bc type

         if (bc(1,1,1).eq.FOEXTRAP.and.ARG_L1(u).lt.domlo(1)) then
!c     Lo x boundary
!c            write (*,*) "Doing 'outflow' boundary conditions for jet"

!c     i on the interior
            ii = domlo(1)
!c     x position of ii
            xi = domnlo(1) + (ii+0.5)*dx(1) - jet_x

            jlo = ARG_L2(u); if (jlo.lt.domlo(2)) jlo=domlo(2);
            jhi = ARG_H2(u); if (jhi.gt.domlo(2)) jhi=domlo(2);

            do i = ARG_L1(u), domlo(1)-1
               x = domnlo(1) + (i+0.5)*dx(1) - jet_x
               do k = ARG_L3(u), ARG_H3(u)
                  do j = jlo, jhi
                     y = domnlo(2) + (j+0.5)*dx(2) - jet_y
!c     radius squared of ghost cell
                     r2 = x*x + y*y
!c     j one cell closer to the origin
                     jp = merge(j+1,j-1,y.lt.zero)
!c     y position of the radius at ii
                     yp = xi*y/x
!c     radius of yp
                     rp = sqrt( xi*xi + yp*yp )
!c     interpolate factor
                     ax = (y-yp)/dx(2)
!c     interpolate velocities to xp
                     uxp = u(ii,j,k,1) + ax*(u(ii,jp,k,1)-u(ii,j,k,1))
                     uyp = u(ii,j,k,2) + ax*(u(ii,jp,k,2)-u(ii,j,k,2))
!c     radial velocity at yp
!c                     up = (uxp*xi-uyp*yp)/rp
                     up = (uxp*xi-uyp*yp)
!c     radial velocity at bc
!c                     u(i,j,k,1) = x*rp*up/r2
!c                     u(i,j,k,2) = y*rp*up/r2
                     u(i,j,k,1) = x*up/r2
                     u(i,j,k,2) = y*up/r2
                     u(i,j,k,3) = zero
                  end do
               end do
            end do
         endif

         if (bc(1,2,1).eq.FOEXTRAP.and.ARG_H1(u).gt.domhi(1)) then
!c     Hi x boundary
!c            write (*,*) "Doing 'outflow' boundary conditions for jet"

!c     i on the interior
            ii = domhi(1)
!c     x position of ji
            xi = domnlo(1) + (ii+0.5)*dx(1) - jet_x

            jlo = ARG_L2(u); if (jlo.lt.domlo(2)) jlo=domlo(2);
            jhi = ARG_H2(u); if (jhi.gt.domlo(2)) jhi=domlo(2);

            do i = domhi(1)+1, ARG_H1(u)
               x = domnlo(1) + (i+0.5)*dx(1) - jet_x
               do k = ARG_L3(u), ARG_H3(u)
                  do j = jlo, jhi
                     y = domnlo(2) + (j+0.5)*dx(2) - jet_y
!c     radius squared of ghost cell
                     r2 = x*x + y*y
!c     j one cell closer to the origin
                     jp = merge(j+1,j-1,y.lt.zero)
!c     y position of the radius at ii
                     yp = xi*y/x
!c     radius of yp
                     rp = sqrt( xi*xi + yp*yp )
!c     interpolate factor
                     ax = (y-yp)/dx(2)
!c     interpolate velocities to xp
                     uxp = u(ii,j,k,1) + ax*(u(ii,jp,k,1)-u(ii,j,k,1))
                     uyp = u(ii,j,k,2) + ax*(u(ii,jp,k,2)-u(ii,j,k,2))
!c     radial velocity at yp
!c                     up = (uxp*xi-uyp*yp)/rp
                     up = (uxp*xi-uyp*yp)/rp
!c     radial velocity at bc
!c                     u(i,j,k,1) = x*rp*up/r2
!c                     u(i,j,k,2) = y*rp*up/r2
                     u(i,j,k,1) = x*up/r2
                     u(i,j,k,2) = y*up/r2
                     u(i,j,k,3) = zero
                  end do
               end do
            end do
         endif

         if (bc(2,1,1).eq.FOEXTRAP.and.ARG_L2(u).lt.domlo(2)) then
!c     Lo y boundary
!c            write (*,*) "Doing 'outflow' boundary conditions for jet"

!c     j on the interior
            ji = domlo(2)
!c     y position of ji
            yi = domnlo(2) + (ji+0.5)*dx(2) - jet_y

            ilo = ARG_L1(u); if (ilo.lt.domlo(1)) ilo=domlo(1);
            ihi = ARG_H1(u); if (ihi.gt.domlo(1)) ihi=domlo(1);

            do j = ARG_L2(u), domlo(2)-1
               y = domnlo(2) + (j+0.5)*dx(2) - jet_y
               do k = ARG_L3(u), ARG_H3(u)
                  do i = ilo, ihi
                     x = domnlo(1) + (i+0.5)*dx(1) - jet_x
!c     radius squared of ghost cell
                     r2 = x*x + y*y
!c     i one cell closer to the origin
                     ip = merge(i+1,i-1,x.lt.zero)
!c     x position of the radius at ji
                     xp = yi*x/y
!c     radius of xp
                     rp = sqrt( xp*xp + yi*yi )
!c     interpolate factor
                     ax = (x-xp)/dx(1)
!c     interpolate velocities to xp (assumes azimuthal component is zero)
                     uxp = u(i,ji,k,1) + ax*(u(ip,ji,k,1)-u(i,ji,k,1))
                     uyp = u(i,ji,k,2) + ax*(u(ip,ji,k,2)-u(i,ji,k,2))
!c     radial velocity at xp
!c                     up = (uxp*xp-uyp*yi)/rp
                     up = (uxp*xp-uyp*yi)
!c     radial velocity at bc
!c                     u(i,j,k,1) = x*rp*up/r2
!c                     u(i,j,k,2) = y*rp*up/r2
                     u(i,j,k,1) = x*up/r2
                     u(i,j,k,2) = y*up/r2
                     u(i,j,k,3) = zero
                  end do
               end do
            end do
         endif

         if (bc(2,2,1).eq.FOEXTRAP.and.ARG_H2(u).gt.domhi(2)) then
!c     Hi y boundary
!c            write (*,*) "Doing 'outflow' boundary conditions for jet"

!c     j on the interior
            ji = domhi(2)
!c     y position of ji
            yi = domnlo(2) + (ji+0.5)*dx(2) - jet_y

            ilo = ARG_L1(u); if (ilo.lt.domlo(1)) ilo=domlo(1);
            ihi = ARG_H1(u); if (ihi.gt.domlo(1)) ihi=domlo(1);

            do j = domhi(2)+1, ARG_H2(u)
               y = domnlo(2) + (j+0.5)*dx(2) - jet_y
               do k = ARG_L3(u), ARG_H3(u)
                  do i = ilo, ihi
                     x = domnlo(1) + (i+0.5)*dx(1) - jet_x
!c     radius squared of ghost cell
                     r2 = x*x + y*y
!c     i one cell closer to the origin
                     ip = merge(i+1,i-1,x.lt.zero)
!c     x position of the radius at ji
                     xp = yi*x/y
!c     radius of xp
                     rp = sqrt( xp*xp + yi*yi )
!c     interpolate factor
                     ax = (x-xp)/dx(1)
!c     interpolate velocities to xp
                     uxp = u(i,ji,k,1) + ax*(u(ip,ji,k,1)-u(i,ji,k,1))
                     uyp = u(i,ji,k,2) + ax*(u(ip,ji,k,2)-u(i,ji,k,2))
!c     radial velocity at xp
!c                     up = (uxp*xp-uyp*yi)/rp
                     up = (uxp*xp-uyp*yi)
!c     radial velocity at bc
!c                     u(i,j,k,1) = x*rp*up/r2
!c                     u(i,j,k,2) = y*rp*up/r2
                     u(i,j,k,1) = x*up/r2
                     u(i,j,k,2) = y*up/r2
                     u(i,j,k,3) = zero
                  end do
               end do
            end do
         endif

         if (bc(3,1,1).eq.EXT_DIR.and.ARG_L3(u).lt.domlo(3)) then
!c     Lo z boundary
            do k = ARG_L3(u), domlo(3)-1
               do j = ARG_L2(u), ARG_H2(u)
                  y = domnlo(2) + (j+0.5)*dx(2) - jet_y
                  do i = ARG_L1(u), ARG_H1(u)
                     x = domnlo(1) + (i+0.5)*dx(1) - jet_x
                     r = sqrt( x*x + y*y )
                     jv= half*((jet_vel+coflow_vel)-(jet_vel-coflow_vel)*tanh(two*(r-jet_width)/delta0))
                     u(i,j,k,1) = zero
                     u(i,j,k,2) = zero
                     u(i,j,k,3) = jv
#ifdef BL_DO_FLCT
                     if (r.lt.jet_width) then
                        u(i,j,k,1) = u(i,j,k,1) + jv*uflct(i,j,1)*turb_scale
                        u(i,j,k,2) = u(i,j,k,2) + jv*vflct(i,j,1)*turb_scale
                        u(i,j,k,3) = u(i,j,k,3) + jv*wflct(i,j,1)*turb_scale
                     endif
#endif
                  enddo
               enddo
            enddo
         endif

         if (bc(3,2,1).eq.EXT_DIR.and.ARG_H3(u).gt.domhi(3)) then
!c     Hi z boundary
            do k = domhi(3)+1, ARG_H3(u)
               do j = ARG_L2(u), ARG_H2(u)
                  do i = ARG_L1(u), ARG_H1(u)
                     u(i,j,k,1) = zero
                     u(i,j,k,2) = zero
                     u(i,j,k,3) = zero
                  end do
               end do
            end do
         endif
      else
!c     Not probtype 18 - do general fill
         call bl_error('General VELFILL not implemented yet')
      endif

#if defined(BL_DO_FLCT)
      if (forceInflow) then
         deallocate(uflct)
         deallocate(vflct)
         deallocate(wflct)
      endif
#endif
      end subroutine FORT_VELFILL

!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.
!c :::
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data.
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: divu     <=  divergence of velocity array
!c ::: DIMS(divu)=> index extent of divu array
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of rho array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
!c ::: -----------------------------------------------------------

      subroutine FORT_DIVUFILL (divu,DIMS(divu),domlo,domhi,dx, &
                               xlo,time,bc)&
                               bind(C, name="FORT_DIVUFILL")
      implicit none

      integer    DIMDEC(divu)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     divu(DIMV(divu))
      integer    bc(SDIM,2)

      integer    i, j, k
      REAL_T     z_vel

#include <probdata.H>

      if (adv_dir .eq. 3) then
         z_vel = adv_vel
      else
         z_vel = zero
      end if

      call filcc(divu,DIMS(divu),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(divu).lt.domlo(1)) then
         do i = ARG_L1(divu), domlo(1)-1
            do k = ARG_L3(divu), ARG_H3(divu)
               do j = ARG_L2(divu), ARG_H2(divu)
                  divu(i,j,k) = z_vel
               end do
	    end do
	 end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(divu).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(divu)
            do k = ARG_L3(divu), ARG_H3(divu)
               do j = ARG_L2(divu), ARG_H2(divu)
                  divu(i,j,k) = z_vel
               end do
	    end do
	 end do
      end if

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(divu).lt.domlo(2)) then
         do j = ARG_L2(divu), domlo(2)-1
            do k = ARG_L3(divu), ARG_H3(divu)
               do i = ARG_L1(divu), ARG_H1(divu)
                  divu(i,j,k) = z_vel
               end do
            end do
	 end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(divu).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(divu)
            do k = ARG_L3(divu), ARG_H3(divu)
               do i = ARG_L1(divu), ARG_H1(divu)
                  divu(i,j,k) = z_vel
               end do
            end do
	 end do
      end if

      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(divu).lt.domlo(3)) then
         do k = ARG_L3(divu), domlo(3)-1
            do j = ARG_L2(divu), ARG_H2(divu)
               do i = ARG_L1(divu), ARG_H1(divu)
                  divu(i,j,k) = z_vel
               end do
            end do
         end do
      end if

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(divu).gt.domhi(3)) then
         do k = domhi(3)+1, ARG_H3(divu)
            do j = ARG_L2(divu), ARG_H2(divu)
               do i = ARG_L1(divu), ARG_H1(divu)
                  divu(i,j,k) = z_vel
               end do
            end do
         end do
      end if

      end subroutine FORT_DIVUFILL

!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.
!c :::
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data.
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: dsdt     <=  dsdt array
!c ::: DIMS(dsdt)=> index extent of dsdt array
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of rho array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
!c ::: -----------------------------------------------------------

      subroutine FORT_DSDTFILL (dsdt,DIMS(dsdt),domlo,domhi,dx, &
                               xlo,time,bc) &
                               bind(C, name="FORT_DSDTFILL")
      implicit none

      integer    DIMDEC(dsdt)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     dsdt(DIMV(dsdt))
      integer    bc(SDIM,2)

      integer    i, j, k

      call filcc(dsdt,DIMS(dsdt),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(dsdt).lt.domlo(1)) then
         do i = ARG_L1(dsdt), domlo(1)-1
            do k = ARG_L3(dsdt), ARG_H3(dsdt)
               do j = ARG_L2(dsdt), ARG_H2(dsdt)
                  dsdt(i,j,k) = zero
               end do
	    end do
	 end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(dsdt).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(dsdt)
            do k = ARG_L3(dsdt), ARG_H3(dsdt)
               do j = ARG_L2(dsdt), ARG_H2(dsdt)
                  dsdt(i,j,k) = zero
               end do
	    end do
	 end do
      end if

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(dsdt).lt.domlo(2)) then
         do j = ARG_L2(dsdt), domlo(2)-1
            do k = ARG_L3(dsdt), ARG_H3(dsdt)
               do i = ARG_L1(dsdt), ARG_H1(dsdt)
                  dsdt(i,j,k) = zero
               end do
            end do
	 end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(dsdt).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(dsdt)
            do k = ARG_L3(dsdt), ARG_H3(dsdt)
               do i = ARG_L1(dsdt), ARG_H1(dsdt)
                  dsdt(i,j,k) = zero
               end do
            end do
	 end do
      end if

      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(dsdt).lt.domlo(3)) then
         do k = ARG_L3(dsdt), domlo(3)-1
            do j = ARG_L2(dsdt), ARG_H2(dsdt)
               do i = ARG_L1(dsdt), ARG_H1(dsdt)
                  dsdt(i,j,k) = zero
               end do
            end do
         end do
      end if

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(dsdt).gt.domhi(3)) then
         do k = domhi(3)+1, ARG_H3(dsdt)
            do j = ARG_L2(dsdt), ARG_H2(dsdt)
               do i = ARG_L1(dsdt), ARG_H1(dsdt)
                  dsdt(i,j,k) = zero
               end do
            end do
         end do
      end if

      end subroutine FORT_DSDTFILL

!c ::: -----------------------------------------------------------
!c ::: This routine is called during a filpatch operation when
!c ::: the patch to be filled falls outside the interior
!c ::: of the problem domain.  You are requested to supply the
!c ::: data outside the problem interior in such a way that the
!c ::: data is consistant with the types of the boundary conditions
!c ::: you specified in the C++ code.
!c :::
!c ::: NOTE:  you can assume all interior cells have been filled
!c :::        with valid data.
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: p        <=  pressure array
!c ::: lo,hi     => index extent of p array
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of rho array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
!c ::: -----------------------------------------------------------

      subroutine FORT_PRESFILL (p,DIMS(p),domlo,domhi,dx, &
                                xlo,time,bc) &
                                bind(C, name="FORT_PRESFILL")
      implicit none

      integer    DIMDEC(p)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     p(DIMV(p))
      integer    bc(SDIM,2)

      integer    i, j, k
      integer    jlo, jhi, ilo, ihi, klo, khi
      logical    fix_xlo, fix_ylo, fix_zlo
      logical    fix_xhi, fix_yhi, fix_zhi
      logical    per_xlo, per_ylo, per_zlo
      logical    per_xhi, per_yhi, per_zhi

#include <probdata.H>

      fix_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .ne. INT_DIR)
      per_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .eq. INT_DIR)
      fix_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .ne. INT_DIR)
      per_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .eq. INT_DIR)
      fix_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .ne. INT_DIR)
      per_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .eq. INT_DIR)
      fix_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .ne. INT_DIR)
      per_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .eq. INT_DIR)
      fix_zlo = (ARG_L3(p) .lt. domlo(3)) .and. (bc(3,1) .ne. INT_DIR)
      per_zlo = (ARG_L3(p) .lt. domlo(3)) .and. (bc(3,1) .eq. INT_DIR)
      fix_zhi = (ARG_H3(p) .gt. domhi(3)) .and. (bc(3,2) .ne. INT_DIR)
      per_zhi = (ARG_H3(p) .gt. domhi(3)) .and. (bc(3,2) .eq. INT_DIR)

      ilo = max(ARG_L1(p),domlo(1))
      ihi = min(ARG_H1(p),domhi(1))
      jlo = max(ARG_L2(p),domlo(2))
      jhi = min(ARG_H2(p),domhi(2))
      Klo = max(ARG_L3(p),domlo(3))
      khi = min(ARG_H3(p),domhi(3))

!c*****************************************************************************
!c SETTING XLO
!c*****************************************************************************

      if (fix_xlo) then
         do i = ARG_L1(p), domlo(1)-1
            do k = klo, khi
               do j = jlo,jhi
                  p(i,j,k) = p(ilo,j,k)
               end do
            end do
	 end do

	 if (fix_ylo) then
	    do i = ARG_L1(p), domlo(1)-1
               do j = ARG_L2(p), domlo(2)-1
                  do k = klo, khi
                     p(i,j,k) = p(ilo,jlo,k)
                  end do
               end do
	    end do

	    if (fix_zlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ilo,jlo,klo)
                     end do
                  end do
               end do
	    else if (per_zlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ilo,jlo,k)
                     end do
                  end do
               end do
	    end if
	    if (fix_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,jlo,khi)
                     end do
                  end do
               end do
	    else if (per_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,jlo,k)
                     end do
                  end do
               end do
	    end if
	 end if

	 if (fix_yhi) then
	    do i = ARG_L1(p), domlo(1)-1
               do j = domhi(2)+1, ARG_H2(p)
                  do k = klo, khi
                     p(i,j,k) = p(ilo,jhi,k)
                  end do
               end do
	    end do
	    if (fix_zlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ilo,jhi,klo)
                     end do
                  end do
               end do
	    else if (per_zlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ilo,jhi,k)
                     end do
                  end do
               end do
	    end if
	    if (fix_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,jhi,khi)
                     end do
                  end do
               end do
	    else if (per_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,jhi,k)
                     end do
                  end do
               end do
	    end if
	 end if

	 if (fix_zlo) then
	    do i = ARG_L1(p), domlo(1)-1
               do j = jlo, jhi
                  do k = ARG_L3(p), domlo(3)-1
                     p(i,j,k) = p(ilo,j,klo)
                  end do
               end do
	    end do
            if (per_ylo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ilo,j,klo)
                     end do
                  end do
               end do
            end if
            if (per_yhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ilo,j,klo)
                     end do
                  end do
               end do
            end if

	 end if

	 if (fix_zhi) then
	    do i = ARG_L1(p), domlo(1)-1
               do j = jlo, jhi
                  do k = domhi(3)+1, ARG_H3(p)
                     p(i,j,k) = p(ilo,j,khi)
                  end do
               end do
	    end do
            if (per_ylo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,j,khi)
                     end do
                  end do
               end do
            end if
            if (per_yhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,j,khi)
                     end do
                  end do
               end do
            end if
	 end if

         if (per_ylo) then
               do i = ARG_L1(p), domlo(1)-1
                  do k = klo,khi
                     do j = ARG_L2(p), domlo(2)-1
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if
         if (per_yhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do k = klo,khi
                     do j = domhi(2)+1, ARG_H2(p)
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if

         if (per_zlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = jlo,jhi
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if
         if (per_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = jlo,jhi
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if

         if (per_ylo .and. per_zlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
	 end if

         if (per_ylo .and. per_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
	 end if

         if (per_yhi .and. per_zlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
	 end if

         if (per_yhi .and. per_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
	 end if

      end if

!c*****************************************************************************
!c SETTING XHI
!c*****************************************************************************

      if (fix_xhi) then
         do i = domhi(1)+1, ARG_H1(p)
            do k = klo, khi
               do j = jlo,jhi
                  p(i,j,k) = p(ihi,j,k)
               end do
            end do
	 end do

	 if (fix_ylo) then
	    do i = domhi(1)+1, ARG_H1(p)
               do j = ARG_L2(p), domlo(2)-1
                  do k = klo, khi
                     p(i,j,k) = p(ihi,jlo,k)
                  end do
               end do
	    end do

	    if (fix_zlo) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ihi,jlo,klo)
                     end do
                  end do
               end do
	    else if (per_zlo) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ihi,jlo,k)
                     end do
                  end do
               end do
	    end if
	    if (fix_zhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,jlo,khi)
                     end do
                  end do
               end do
	    else if (per_zhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,jlo,k)
                     end do
                  end do
               end do
	    end if
	 end if
	 if (fix_yhi) then
	    do i = domhi(1)+1, ARG_H1(p)
               do j = domhi(2)+1, ARG_H2(p)
                  do k = klo, khi
                     p(i,j,k) = p(ihi,jhi,k)
                  end do
               end do
	    end do
	    if (fix_zlo) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ihi,jhi,klo)
                     end do
                  end do
               end do
	    else if (per_zlo) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ihi,jhi,k)
                     end do
                  end do
               end do
	    end if
	    if (fix_zhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,jhi,khi)
                     end do
                  end do
               end do
	    else if (per_zhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,jhi,k)
                     end do
                  end do
               end do
	    end if
	 end if

	 if (fix_zlo) then
	    do i = domhi(1)+1, ARG_H1(p)
               do j = jlo, jhi
                  do k = ARG_L3(p), domlo(3)-1
                     p(i,j,k) = p(ihi,j,klo)
                  end do
               end do
	    end do
            if (per_ylo) then
	       do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ihi,j,klo)
                     end do
                  end do
               end do
            end if
            if (per_yhi) then
	       do i = domhi(1)+1, ARG_H1(p)
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ihi,j,klo)
                     end do
                  end do
               end do
            end if

	 end if

	 if (fix_zhi) then
	    do i = domhi(1)+1, ARG_H1(p)
               do j = jlo, jhi
                  do k = domhi(3)+1, ARG_H3(p)
                     p(i,j,k) = p(ihi,j,khi)
                  end do
               end do
	    end do
            if (per_ylo) then
	       do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,j,khi)
                     end do
                  end do
               end do
            end if
            if (per_yhi) then
	       do i = domhi(1)+1, ARG_H1(p)
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,j,khi)
                     end do
                  end do
               end do
            end if
	 end if

         if (per_ylo) then
	       do i = domhi(1)+1, ARG_H1(p)
                  do k = klo,khi
                     do j = ARG_L2(p), domlo(2)-1
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if
         if (per_yhi) then
	       do i = domhi(1)+1, ARG_H1(p)
                  do k = klo,khi
                     do j = domhi(2)+1, ARG_H2(p)
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if

         if (per_zlo) then
	       do i = domhi(1)+1, ARG_H1(p)
                  do j = jlo,jhi
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if
         if (per_zhi) then
	       do i = domhi(1)+1, ARG_H1(p)
                  do j = jlo,jhi
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if


         if (per_ylo .and. per_zlo) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if

         if (per_ylo .and. per_zhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if

         if (per_yhi .and. per_zlo) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if

         if (per_yhi .and. per_zhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if

      end if

!c*****************************************************************************
!c SETTING YLO
!c*****************************************************************************

      if (fix_ylo) then
         do j = ARG_L2(p), domlo(2)-1
            do k = klo, khi
               do i = ilo, ihi
                  p(i,j,k) = p(i,jlo,k)
               end do
            end do
	 end do

	 if (fix_zlo) then
	    do j = ARG_L2(p), domlo(2)-1
               do k = ARG_L3(p), domlo(3)-1
                  do i = ilo, ihi
                     p(i,j,k) = p(i,jlo,klo)
                  end do
               end do
	    end do
            if (per_xlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jlo,klo)
                     end do
                  end do
               end do
            end if
            if (per_xhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jlo,klo)
                     end do
                  end do
               end do
            end if
	 end if

	 if (fix_zhi) then
	    do j = ARG_L2(p), domlo(2)-1
               do k = domhi(3)+1, ARG_H3(p)
                  do i = ilo, ihi
                     p(i,j,k) = p(i,jlo,khi)
                  end do
               end do
	    end do
            if (per_xlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jlo,khi)
                     end do
                  end do
               end do
            end if
            if (per_xhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jlo,khi)
                     end do
                  end do
               end do
            end if
	 end if

         if (per_xlo) then
               do j = ARG_L2(p), domlo(2)-1
                  do k = klo,khi
                     do i = ARG_L1(p), domlo(1)-1
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if
         if (per_xhi) then
               do j = ARG_L2(p), domlo(2)-1
                  do k = klo,khi
                     do i = domhi(1)+1, ARG_H1(p)
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if

         if (per_zlo) then
               do j = ARG_L2(p), domlo(2)-1
                  do i = ilo,ihi
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if
         if (per_zhi) then
               do j = ARG_L2(p), domlo(2)-1
                  do i = ilo,ihi
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if


         if (per_xlo .and. per_zlo) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_zlo) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_zhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if

      end if

!c*****************************************************************************
!c SETTING YHI
!c*****************************************************************************

      if (fix_yhi) then
         do j = domhi(2)+1, ARG_H2(p)
            do k = klo, khi
               do i = ilo, ihi
                  p(i,j,k) = p(i,jhi,k)
               end do
            end do
	 end do

	 if (fix_zlo) then
	    do j = domhi(2)+1, ARG_H2(p)
               do k = ARG_L3(p), domlo(3)-1
                  do i = ilo, ihi
                     p(i,j,k) = p(i,jhi,klo)
                  end do
               end do
	    end do
            if (per_xlo) then
               do i = ARG_L1(p), domlo(1)-1
	          do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jhi,klo)
                     end do
                  end do
               end do
            end if
            if (per_xhi) then
               do i = domhi(1)+1, ARG_H1(p)
	          do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jhi,klo)
                     end do
                  end do
               end do
            end if
	 end if

	 if (fix_zhi) then
	    do j = domhi(2)+1, ARG_H2(p)
               do k = domhi(3)+1, ARG_H3(p)
                  do i = ilo, ihi
                     p(i,j,k) = p(i,jhi,khi)
                  end do
               end do
	    end do
            if (per_xlo) then
               do i = ARG_L1(p), domlo(1)-1
	          do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jhi,khi)
                     end do
                  end do
               end do
            end if
            if (per_xhi) then
               do i = domhi(1)+1, ARG_H1(p)
	          do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jhi,khi)
                     end do
                  end do
               end do
            end if
	 end if

         if (per_xlo) then
               do j = domhi(2)+1, ARG_H2(p)
                  do k = klo,khi
                     do i = ARG_L1(p), domlo(1)-1
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if
         if (per_xhi) then
               do j = domhi(2)+1, ARG_H2(p)
                  do k = klo,khi
                     do i = domhi(1)+1, ARG_H1(p)
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

         if (per_zlo) then
               do j = domhi(2)+1, ARG_H2(p)
                  do i = ilo,ihi
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if
         if (per_zhi) then
               do j = domhi(2)+1, ARG_H2(p)
                  do i = ilo,ihi
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_zlo) then
               do i = ARG_L1(p), domlo(1)-1
	          do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_zhi) then
               do i = ARG_L1(p), domlo(1)-1
	          do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_zlo) then
               do i = domhi(1)+1, ARG_H1(p)
	          do j = domhi(2)+1, ARG_H2(p)
                     do k = ARG_L3(p), domlo(3)-1
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_zhi) then
               do i = domhi(1)+1, ARG_H1(p)
	          do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

      end if

!c*****************************************************************************
!c SETTING ZLO
!c*****************************************************************************

      if (fix_zlo) then
         do k = ARG_L3(p), domlo(3)-1
            do j = jlo, jhi
               do i = ilo, ihi
                  p(i,j,k) = p(i,j,klo)
               end do
            end do
	 end do

         if (per_xlo) then
               do k = ARG_L3(p), domlo(3)-1
                  do j = jlo,jhi
                     do i = ARG_L1(p), domlo(1)-1
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if
         if (per_xhi) then
               do k = ARG_L3(p), domlo(3)-1
                  do j = jlo,jhi
                     do i = domhi(1)+1, ARG_H1(p)
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_ylo) then
               do k = ARG_L3(p), domlo(3)-1
                  do i = ilo,ihi
                     do j = ARG_L2(p), domlo(2)-1
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if
         if (per_yhi) then
               do k = ARG_L3(p), domlo(3)-1
                  do i = ilo,ihi
                     do j = domhi(2)+1, ARG_H2(p)
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_ylo) then
               do k = ARG_L3(p), domlo(3)-1
                  do i = ARG_L1(p), domlo(1)-1
                     do j = ARG_L2(p), domlo(2)-1
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_yhi) then
               do k = ARG_L3(p), domlo(3)-1
                  do i = ARG_L1(p), domlo(1)-1
                     do j = domhi(2)+1, ARG_H2(p)
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_ylo) then
               do k = ARG_L3(p), domlo(3)-1
                  do i = domhi(1)+1, ARG_H1(p)
                     do j = ARG_L2(p), domlo(2)-1
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_yhi) then
               do k = ARG_L3(p), domlo(3)-1
                  do i = domhi(1)+1, ARG_H1(p)
                     do j = domhi(2)+1, ARG_H2(p)
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

      end if

!c*****************************************************************************
!c SETTING ZHI
!c*****************************************************************************

      if (fix_zhi) then
         do k = domhi(3)+1, ARG_H3(p)
            do j = jlo, jhi
               do i = ilo, ihi
                  p(i,j,k) = p(i,j,khi)
               end do
            end do
	 end do

         if (per_xlo) then
               do k = domhi(3)+1, ARG_H3(p)
                  do j = jlo,jhi
                     do i = ARG_L1(p), domlo(1)-1
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if
         if (per_xhi) then
               do k = domhi(3)+1, ARG_H3(p)
                  do j = jlo,jhi
                     do i = domhi(1)+1, ARG_H1(p)
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

         if (per_ylo) then
               do k = domhi(3)+1, ARG_H3(p)
                  do i = ilo,ihi
                     do j = ARG_L2(p), domlo(2)-1
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if
         if (per_yhi) then
               do k = domhi(3)+1, ARG_H3(p)
                  do i = ilo,ihi
                     do j = domhi(2)+1, ARG_H2(p)
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if


         if (per_xlo .and. per_ylo) then
               do k = domhi(3)+1, ARG_H3(p)
                  do i = ARG_L1(p), domlo(1)-1
                     do j = ARG_L2(p), domlo(2)-1
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_yhi) then
               do k = domhi(3)+1, ARG_H3(p)
                  do i = ARG_L1(p), domlo(1)-1
                     do j = domhi(2)+1, ARG_H2(p)
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_ylo) then
               do k = domhi(3)+1, ARG_H3(p)
                  do i = domhi(1)+1, ARG_H1(p)
                     do j = ARG_L2(p), domlo(2)-1
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_yhi) then
               do k = domhi(3)+1, ARG_H3(p)
                  do i = domhi(1)+1, ARG_H1(p)
                     do j = domhi(2)+1, ARG_H2(p)
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

      end if

!c*****************************************************************************

      end subroutine FORT_PRESFILL

!***************************************************************
!*    "Minimal" random number generator of Park and Miller with
!*    Bays-Durham shuffle and added safeguards.  Returns a uniform random
!*    deviate between 0.0 and 1.0 (exclusive of the endpoint values).
!*    Call with IDUM a negative integer to initialize; thereafter, do not
!*    alter IDUM between successive deviates in a sequence.  RNMX should
!*    approximate the largest floating value that is less than 1.

      function ran1(idum)
      integer idum,ia,im,iq,ir,ntab,ndiv
      REAL_T ran1,am,eps,rnmx
      parameter (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836, &
          ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-7,rnmx=1.-eps)
      integer j,k,iv(ntab),iy
      save iv,iy
      data iv /ntab*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
         idum=max(-idum,1)
         do 11 j=ntab+8,1,-1
            k=idum/iq
            idum=ia*(idum-k*iq)-ir*k
            if (idum.lt.0) idum=idum+im
            if (j.le.ntab) iv(j)=idum
 11      continue
         iy=iv(1)
      endif
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      ran1=min(am*iy,rnmx)
      return
      end function ran1
!*
!*     After Press et al., Numerical Recipes for Fortran


end module prob_3D_module
