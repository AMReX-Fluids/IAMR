
c
c $Id: PROB_3D.F,v 1.43 2007-12-21 18:38:17 aaspden Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "PROB_AMR_F.H"
#include "PROB_NS_F.H"
#include "ArrayLim.H"

#ifdef BL_DO_FLCT
#include "infl_frc.H"
#endif

#define SDIM 3

      block data rt_common
#include "probdata.H"
      data rt_pertamp/0.0D0/
      data rt_nfreq/0/
      data rt_xfrontw/0.0d0/
      end
c ::: -----------------------------------------------------------
c ::: This routine is called at problem initialization time
c ::: and when restarting from a checkpoint file.
c ::: The purpose is (1) to specify the initial time value
c ::: (not all problems start at time=0.0) and (2) to read
c ::: problem specific data from a namelist or other input
c ::: files and possibly store them or derived information
c ::: in FORTRAN common blocks for later use.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: init      => TRUE if called at start of problem run
c :::              FALSE if called from restart
c ::: name      => name of "probin" file
c ::: namlen    => length of name
c ::: strttime <=  start problem with this time variable
c ::: 
c ::: -----------------------------------------------------------
      subroutine FORT_PROBINIT (init,name,namlen,problo,probhi)
      implicit none
      integer init,namlen
      integer name(namlen)
      integer untin, i
      REAL_T  problo(SDIM), probhi(SDIM)

#include "probdata.H"

#ifdef DO_IAMR_FORCE
#include "forcedata.H"
#endif

#ifdef BL_DO_FLCT
#include "INFL_FORCE_F.H"
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

      REAL_T  Lx, Ly, Lz, Lmin
      REAL_T  kappa, kappaMax, freqMin, freqMax, freqDiff, pdk

      namelist /fortin/ denerr, vorterr, adverr, temperr,
     &			denfact, xblob, yblob, zblob, radblob, 
     &                  velfact, probtype, randfact, bubgrad,
     &			rhozero, tempzero, c_d, r_d, 
     &                  adv_dir, adv_vel, axis_dir, radvort,
     &                  den1,den2,vel1,vel2,delta0,xlev1,zlev1,amag,
     &                  vb_unifdir, blrandseed, turb_scale, injection_time
#ifdef BL_DO_FLCT
      namelist /fortin/ forceInflow, numInflPlanesStore, strmwse_dir,
     &                  forceLo, forceHi, flct_file, nCompInflow
#endif
      namelist /fortin/ rt_splitx, rt_xfrontw, rt_den_1, rt_den_2,
     &                  rt_pertamp, rt_nfreq, rt_graddenerr

      namelist /fortin/ Vco, Rfu, Rtran, tVco_l, tVco_r, Vco_l, Vco_r

      namelist /fortin/ grav_angle, omega, infl_time_offset, ref_centre, ref_radius,
     &                  thermal_expansion, heating_coeff, heating_centre, heating_radius


#ifdef DO_IAMR_FORCE
      namelist /fortin/ nmodes, nxmodes, nymodes, nzmodes, mode_start,
     &                  forcing_type, spectrum_type, ref_type, forcing_twice_wavelength,
     &                  forcing_xlength, forcing_ylength, forcing_zlength,
     &                  forcing_time_scale_min, forcing_time_scale_max, 
     &                  force_scale, forcing_epsilon,
     &                  use_rho_in_forcing, do_mode_division,
     &                  AXY, BXY, CXY, DXY, PXY, QXY, RXY,
     &                  AZX, BZX, CZX, DZX, PZX, QZX, RZX,
     &                  AYZ, BYZ, CYZ, DYZ, PYZ, QYZ, RYZ,
     &                  FTX, FTY, FTZ, TAT, TAP,
     &                  FPX, FPY, FPZ, FAX, FAY, FAZ
#endif

      namelist /fortin/ jet_x, jet_y, jet_width, jet_rho, jet_vel, coflow_rho, coflow_vel, jet_temp, ref_height
c
c      Build "probin" filename -- the name of file containing fortin namelist.
c
      integer maxlen, isioproc
      parameter (maxlen=256)

      REAL_T frecon, onep7
      parameter (frecon=.219, onep7=1.7)

      character probin*(maxlen)

      integer j, m, n
      integer MAXPHASE
      parameter (MAXPHASE = 1000)
      DOUBLE PRECISION rn, permin, permax, xtmp, ytmp, pert

      call bl_pd_is_ioproc(isioproc)

      if (namlen .gt. maxlen) call bl_error('probin file name too long')

      do i = 1, namlen
         probin(i:i) = char(name(i))
      end do

#ifdef BL_DO_FLCT
      forceInflow = .FALSE.
      numInflPlanesStore = 16
      forceLo = .TRUE.
      forceHi = .FALSE.
      strmwse_dir = FLCT_ZVEL
      flct_file = ""
      turb_scale = 1
#endif

      tVco_l = -1
      tVco_r = -1
      Vco_l = -1
      Vco_r = -1

      untin = 9
      if (namlen .eq. 0) then
         open(untin,file='probin',form='formatted',status='old')
      else
         open(untin,file=probin(1:namlen),form='formatted',status='old')
      end if
      read(untin,fortin)
c      if (isioproc .eq. 1) write(6,fortin)
      close(unit=untin)

#ifdef BL_DO_FLCT
      if (forceInflow .eqv. .FALSE.) then
         forceLo = .FALSE.
         forceHi = .FALSE.
      else
         if (flct_file .ne. "") then
            untin=20
            if (isioproc .eq. 1) print*, 'Initializing turbulence ...'
c     New file format: change as of Mike's email
c            open(untin, file=flct_file, form='unformatted')
            open(untin, file=trim(flct_file)//'/HDR', form='formatted')
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
c
c     Initialize the common blocks
c
      do i=1, SDIM
         f_problo(i) = problo(i)
         f_probhi(i) = probhi(i)
      enddo

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
                        pert = pert
     &                       + sin(
     &                       2.0D0*rt_PI*dble(n)*xtmp + rt_ranphse(n,m,1)
     &                       )
     &                       * sin(
     &                       2.0D0*rt_PI*dble(m)*ytmp + rt_ranphse(n,m,2)
     &                       )
     &                       *rt_ranampl(n,m)
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
c     Translation angles, theta=0..2Pi and phi=0..Pi
                     call blutilrand(rn)
                     TAT(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     TAP(kx,ky,kz) = rn*Pi
c     Phases
                     call blutilrand(rn)
                     FPX(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     FPY(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     FPZ(kx,ky,kz) = rn*twicePi
c     Amplitudes (alpha)
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
                     if (mp2 .lt. 0.000001) then
                        write(*,*) "ZERO AMPLITUDE MODE ",kx,ky,kz
                        FAX(kx,ky,kz) = zero
                        FAY(kx,ky,kz) = zero
                        FAZ(kx,ky,kz) = zero
                     else
c     Count modes that contribute
                        mode_count = mode_count + 1
c     Set amplitudes
                        if (spectrum_type.eq.1) then
                           Ekh        = one / kappa
                        else if (spectrum_type.eq.2) then
                           Ekh        = one / (kappa*kappa)
                        else
                           Ekh        = one
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
                           write (*,*) "kappa = ",kx,ky,kz,kappa
                           write (*,*) "Amplitudes"
                           write (*,*) FAX(kx,ky,kz), FAY(kx,ky,kz), FAZ(kx,ky,kz)
                           write (*,*) "Frequencies"
                           write (*,*) FTX(kx,ky,kz), FTY(kx,ky,kz), FTZ(kx,ky,kz)
                        endif
                     endif
                  endif   
               enddo
            enddo
         enddo

c     Now let's break symmetry, have to assume high aspect ratio in z for now
         reduced_mode_count = 0
         do kz = mode_start, zstep - 1
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
c     Translation angles, theta=0..2Pi and phi=0..Pi
                     call blutilrand(rn)
                     TAT(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     TAP(kx,ky,kz) = rn*Pi
c     Phases
                     call blutilrand(rn)
                     FPX(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     FPY(kx,ky,kz) = rn*twicePi
                     call blutilrand(rn)
                     FPZ(kx,ky,kz) = rn*twicePi
c     Amplitudes (alpha)
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
                     if (mp2 .lt. 0.000001) then
                        write(*,*) "ZERO AMPLITUDE MODE ",kx,ky,kz
                        FAX(kx,ky,kz) = zero
                        FAY(kx,ky,kz) = zero
                        FAZ(kx,ky,kz) = zero
                     else
c     Count modes that contribute
                        reduced_mode_count = reduced_mode_count + 1
c     Set amplitudes
                        if (spectrum_type.eq.1) then
                           Ekh        = one / kappa
                        else if (spectrum_type.eq.2) then
                           Ekh        = one / (kappa*kappa)
                        else
                           Ekh        = one
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
                           write (*,*) "kappa = ",kx,ky,kz,kappa
                           write (*,*) "Amplitudes"
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
      endif
#endif

      end

c ::: -----------------------------------------------------------
c ::: This routine is called at problem setup time and is used
c ::: to initialize data on each grid.  The velocity field you
c ::: provide does not have to be divergence free and the pressure
c ::: field need not be set.  A subsequent projection iteration
c ::: will define aa divergence free velocity field along with a
c ::: consistant pressure.
c ::: 
c ::: NOTE:  all arrays have one cell of ghost zones surrounding
c :::        the grid interior.  Values in these cells need not
c :::        be set here.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: level     => amr level of grid
c ::: time      => time at which to init data             
c ::: lo,hi     => index limits of grid interior (cell centered)
c ::: nscal     => number of scalar quantities.  You should know
c :::		   this already!
c ::: vel      <=  Velocity array
c ::: scal     <=  Scalar array
c ::: press    <=  Pressure array
c ::: dx     => cell size
c ::: xlo,xhi   => physical locations of lower left and upper
c :::              right hand corner of grid.  (does not include
c :::		   ghost region).
c ::: -----------------------------------------------------------
      subroutine FORT_INITDATA(level,time,lo,hi,nscal,
     &	 	               vel,scal,DIMS(state),press,DIMS(press),
     &                         dx,xlo,xhi)
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

#include "probdata.H"

      if (probtype .eq. 1) then
         call initspin(level,time,lo,hi,nscal,
     &     	       vel,scal,DIMS(state),press,DIMS(press),
     &                 dx,xlo,xhi)

      else if (probtype .eq. 2) then
         call initbubble(level,time,lo,hi,nscal,
     &     	         vel,scal,DIMS(state),press,DIMS(press),
     &                   dx,xlo,xhi)

      else if (probtype .eq. 3) then
         call initvort(level,time,lo,hi,nscal,
     &     	       vel,scal,DIMS(state),press,DIMS(press),
     &                 dx,xlo,xhi)

      else if (probtype .eq. 4) then
         call initchannel(level,time,lo,hi,nscal,
     &     	          vel,scal,DIMS(state),press,DIMS(press),
     &                    dx,xlo,xhi)

      else if (probtype .eq. 5) then
         call initpervort(level,time,lo,hi,nscal,
     &     	          vel,scal,DIMS(state),press,DIMS(press),
     &                    dx,xlo,xhi)

      else if (probtype .eq. 6) then
         call inithotspot(level,time,lo,hi,nscal,
     &     	          vel,scal,DIMS(state),press,DIMS(press),
     &                    dx,xlo,xhi)

      else if (probtype .eq. 7) then
         call initeuler(level,time,lo,hi,nscal,
     &                    vel,scal,DIMS(state),press,DIMS(press),
     &                    dx,xlo,xhi)

      else if (probtype .eq. 8) then
         call initbrnrsh(level,time,lo,hi,nscal,
     &                    vel,scal,DIMS(state),press,DIMS(press),
     &                    dx,xlo,xhi)

      else if (probtype .eq. 9) then
         call initviscbench(level,time,lo,hi,nscal,
     &                      vel,scal,DIMS(state),press,DIMS(press),
     &                      dx,xlo,xhi)


      else if (probtype .eq. 10) then
         call initrt(level,time,lo,hi,nscal,
     &        vel,scal,DIMS(state),press,DIMS(press),
     &        dx,xlo,xhi)


      else if (probtype .eq. 11) then
         call inflowt(level,time,lo,hi,nscal,
     &        vel,scal,DIMS(state),press,DIMS(press),
     &        dx,xlo,xhi)

      else if (probtype .eq. 12) then
         call taylorgreen(level,time,lo,hi,nscal,
     &        vel,scal,DIMS(state),press,DIMS(press),
     &        dx,xlo,xhi)

      else if (probtype .eq. 13) then
         call swirl(level,time,lo,hi,nscal,
     &        vel,scal,DIMS(state),press,DIMS(press),
     &        dx,xlo,xhi)

      else if (probtype .eq. 14) then
         call inithit(level,time,lo,hi,nscal,
     &        vel,scal,DIMS(state),press,DIMS(press),
     &        dx,xlo,xhi)

      else if (probtype .eq. 15) then
         call initinflow(level,time,lo,hi,nscal,
     &        vel,scal,DIMS(state),press,DIMS(press),
     &        dx,xlo,xhi)

      else if (probtype .eq. 16) then
         call initjet(level,time,lo,hi,nscal,
     &        vel,scal,DIMS(state),press,DIMS(press),
     &        dx,xlo,xhi)

      else if (probtype .eq. 17) then
         call initinjection(level,time,lo,hi,nscal,
     &        vel,scal,DIMS(state),press,DIMS(press),
     &        dx,xlo,xhi)

      else if (probtype .eq. 18) then
         call initroundjet(level,time,lo,hi,nscal,
     &        vel,scal,DIMS(state),press,DIMS(press),
     &        dx,xlo,xhi)

      else if (probtype .eq. 19) then
         call initcoriolis(level,time,lo,hi,nscal,
     &        vel,scal,DIMS(state),press,DIMS(press),
     &        dx,xlo,xhi)

      else if (probtype .eq. 20) then
         call initheating(level,time,lo,hi,nscal,
     &        vel,scal,DIMS(state),press,DIMS(press),
     &        dx,xlo,xhi)

      else
         write(6,*) "INITDATA: bad probtype = ",probtype
      end if
      end

c----------------------------------------------------------------------
c     A handy statement function: here blend goes from 0 to 1 in x at
c     rad, over a width of trn

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
      end
c
c ::: -----------------------------------------------------------





c
      subroutine initrt(level,time,lo,hi,nscal,
     &	 	            vel,scal,DIMS(state),press,DIMS(press),
     &                      dx,xlo,xhi)
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
      REAL_T  zblend1
      EXTERNAL zblend1
c
c     ::::: local variables
c
      integer i, j, k, n, m
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  pert,ztemp

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      do j = lo(2), hi(2)
         y = xlo(2) + hy*(dble(j-lo(2)) + 0.5d0)
         do i = lo(1), hi(1)
            x = xlo(1) + hx*(dble(i-lo(1)) + 0.5d0)
            pert = 0.d0
            do n = 1, rt_nfreq
               do m = 1, rt_nfreq
                  pert = pert
     &                 + sin(
     &                 2.0D0*rt_PI*dble(n)*x/f_probhi(1) + rt_ranphse(n,m,1)
     &                 )
     &                 * sin(
     &                 2.0D0*rt_PI*dble(m)*y/f_probhi(2) + rt_ranphse(n,m,2)
     &                 )
     &                 *rt_ranampl(n,m)
               end do
            end do
            ztemp = rt_splitx - pert*rt_pertamp
            do k = lo(3), hi(3)
               z = xlo(3) + hz*(dble(k-lo(3)) + 0.5d0)
               scal(i,j,k,1) = 
     &              rt_den_1 
     &              + (rt_den_2-rt_den_1)*zblend1(z,ztemp,rt_xfrontw)
               vel(i,j,k,1) = 0
               vel(i,j,k,2) = 0
               vel(i,j,k,3) = 0
               do n = 2,nscal-1
                  scal(i,j,k,n) = one
               end do                  
               scal(i,j,k,nscal) = 0
   	    end do
         end do
      end do

      end
c
c ::: -----------------------------------------------------------
c
      subroutine initbubble(level,time,lo,hi,nscal,
     &	 	            vel,scal,DIMS(state),press,DIMS(press),
     &                      dx,xlo,xhi)
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
c
c     ::::: local variables
c
      integer i, j, k, n
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  dist
      REAL_T  x_vel, y_vel, z_vel

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      if (adv_dir .eq. 1) then
         x_vel = adv_vel
         y_vel = zero
         z_vel = zero
      else if (adv_dir .eq. 2) then
         x_vel = zero
         y_vel = adv_vel
         z_vel = zero
c     AJA - This next line originally said 2!
      else if (adv_dir .eq. 3) then
         x_vel = zero
         y_vel = zero
         z_vel = adv_vel
      else 
         write(6,*) "initbubble: adv_dir = ",adv_dir
         stop
      end if

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
            do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half)
               dist = sqrt((x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2)
               vel(i,j,k,1) = x_vel
               vel(i,j,k,2) = y_vel
               vel(i,j,k,3) = z_vel
               scal(i,j,k,1) = cvmgt(denfact,one,dist.lt.radblob)
               do n = 2,nscal-1
                  scal(i,j,k,n) = one
               end do                  
               scal(i,j,k,nscal) = cvmgt(one,zero,dist.lt.radblob)
   	    end do
         end do
      end do

      end
c
c ::: -----------------------------------------------------------
c
      subroutine initpervort(level,time,lo,hi,nscal,
     &                       vel,scal,DIMS(state),press,DIMS(press),
     &                       dx,xlo,xhi)
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
c
c     ::::: local variables
c
      integer i, j, k, n
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  dist

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
            do i = lo(1), hi(1)

               x = xlo(1) + hx*(float(i-lo(1)) + half)

               vel(i,j,k,1) = tanh(30.*(.25-abs(y-.5)))
               vel(i,j,k,2) = .05*sin(two*Pi*x)
               vel(i,j,k,3) = zero

               scal(i,j,k,1) = one
               do n = 2,nscal-1
                  scal(i,j,k,n) = one
               end do

               dist = sqrt((x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2)
               scal(i,j,k,nscal) = cvmgt(one,zero,dist.lt.radblob)
            end do
         end do
      end do

      end
c

c
c ::: -----------------------------------------------------------
c

      subroutine initspin(level,time,lo,hi,nscal,
     &	 	          vel,scal,DIMS(state),press,DIMS(press),
     &                    dx,xlo,xhi)
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
c
c     ::::: local variables
c
      integer i, j, k, n
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  dist
      REAL_T  x_vel, y_vel, z_vel
      REAL_T  spx, spy, cpx, cpy

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)
      hz = dx(2)

      if (adv_dir .eq. 1) then
         x_vel = adv_vel
         y_vel = zero
         z_vel = zero
      else if (adv_dir .eq. 2) then
         x_vel = zero
         y_vel = adv_vel
         z_vel = zero
c     AJA - This next line also originally said 2!
      else if (adv_dir .eq. 3) then
         x_vel = zero
         y_vel = zero
         z_vel = adv_vel
      else
         write(6,*) "INITSPIN: adv_dir = ",adv_dir
         stop
      end if

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
	    spy = sin(Pi*y)
	    cpy = cos(Pi*y)
	    do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half)

 	       spx = sin(Pi*x)
               cpx = cos(Pi*x)

               vel(i,j,k,1) = x_vel - velfact*two*spy*cpy*spx**2
               vel(i,j,k,2) = y_vel + velfact*two*spx*cpx*spy**2
               vel(i,j,k,3) = z_vel

  	       dist = sqrt((x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2)

               scal(i,j,k,1) = one + (denfact-one) * tanh(10.*(dist-radblob))
               do n = 2,nscal-1
                  scal(i,j,k,n) = one
               end do                  
	       scal(i,j,k,nscal) = cvmgt(one,zero,dist.lt.radblob)

	    end do
         end do
      end do

      end
c
c ::: -----------------------------------------------------------
c ::: This case is an unsteady  viscous benchmark for which the
c ::: exact solution is,
c :::     u(x,y,t) = - Cos(Pi x) Sin(Pi y) Exp(-2 Pi^2 Nu t)
c :::     v(x,y,t) =   Sin(Pi x) Cos(Pi y) Exp(-2 Pi^2 Nu t)
c :::     p(x,y,t) = - {Cos(2 Pi x) + Cos(2 Pi y)} Exp(-4 Pi^2 Nu t) / 4
c ::: In the utilities, iamrlib/BenchMarks, there is a
c ::: tool ViscBench2d.cpp that reads a plot file and compares the
c ::: solution against this exact solution.  This benchmark was
c ::: originally derived by G.I. Taylor (Phil. Mag., Vol. 46, No. 274,
c ::: pp. 671-674, 1923) and Ethier and Steinman
c ::: (Intl. J. Num. Meth. Fluids, Vol. 19, pp. 369-375, 1994) give
c ::: the pressure field.
c
      subroutine initviscbench(level,time,lo,hi,nscal,
     &                         vel,scal,DIMS(state),press,DIMS(press),
     &                         dx,xlo,xhi)
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
c
c     ::::: local variables
c
      integer i, j, k, n
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  spx, spy, spz, cpx, cpy, cpz

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      do k = lo(3), hi(3)
        z = xlo(3) + hz*(float(k-lo(3)) + half)
        spz = sin(Pi*z)
        cpz = cos(Pi*z)

        do j = lo(2), hi(2)
          y = xlo(2) + hy*(float(j-lo(2)) + half)
          spy = sin(Pi*y)
          cpy = cos(Pi*y)

          do i = lo(1), hi(1)
            x = xlo(1) + hx*(float(i-lo(1)) + half)

            spx = sin(Pi*x)
            cpx = cos(Pi*x)

c
c           Uniform in the X-direction
c
            if (vb_unifdir .eq. 0) then
              vel(i,j,k,1) =   zero
              vel(i,j,k,2) =   spz*cpy
              vel(i,j,k,3) = - cpz*spy
              scal(i,j,k,1) =   one
              do n = 2, nscal
                scal(i,j,k,n) =   cpz*cpy
              enddo

c
c           Uniform in the Y-direction
c
            elseif (vb_unifdir .eq. 1) then
              vel(i,j,k,1) = - cpx*spz
              vel(i,j,k,2) =   zero
              vel(i,j,k,3) =   spx*cpz
              scal(i,j,k,1) =   one
              do n = 2, nscal
                scal(i,j,k,n) =   cpx*cpz
              enddo

c
c           Uniform in the Z-direction
c
            elseif (vb_unifdir .eq. 2) then
              vel(i,j,k,1) = - cpx*spy
              vel(i,j,k,2) =   spx*cpy
              vel(i,j,k,3) =   zero
              scal(i,j,k,1) =   one
              do n = 2, nscal
                scal(i,j,k,n) =   cpx*cpy
              enddo
            endif
            end do
         end do
      end do
      end
c
c ::: -----------------------------------------------------------
c
      subroutine initvort(level,time,lo,hi,nscal,
     &	 	          vel,scal,DIMS(state),press,DIMS(press),
     &                    dx,xlo,xhi)
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
c
c     ::::: local variables
c
      integer i, j, k, n
      REAL_T  x, y, z, r
      REAL_T  hx, hy, hz
      REAL_T  c, ux, uy, uz
      REAL_T  umagin, umagout, absu, sinth, costh
      REAL_T  small, a, b, r0

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)
      small = 1.0e-10

      r0 = two/three * radvort
      a = one / ((radvort - r0)*(two*radvort - r0))
      b = a * radvort**2 * (radvort - r0)

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half) - zblob
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half) - yblob
	    do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half) - xblob
  	       r = sqrt(x**2 + y**2)
c              umagin = .5*r - 4*r**3
c              umagout = radvort*(.5*radvort - 4*radvort**3)/max(radvort,r)
               umagin = velfact * (one - a*(r - r0)**2)
               umagout = velfact * b/max(radvort,r)
               absu = cvmgp(umagout,umagin,r - radvort)
               sinth = y/max(r,small*radvort)
               costh = x/max(r,small*radvort)
               vel(i,j,k,1) = -absu*sinth
               vel(i,j,k,2) = absu*costh
               vel(i,j,k,3) = zero
               scal(i,j,k,1) = cvmgt(denfact,one,r.lt.radblob)
               do n = 2,nscal-1
                  scal(i,j,k,n) = one
               end do                  
               scal(i,j,k,nscal) = cvmgt(one,zero,r.lt.radblob)
	    end do
         end do
      end do

      end
c
c ::: -----------------------------------------------------------
c
      subroutine initchannel(level,time,lo,hi,nscal,
     &	 	             vel,scal,DIMS(state),press,DIMS(press),
     &                       dx,xlo,xhi)
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
c
c     ::::: local variables
c
      integer i, j, k, n
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  dist

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
            do i = lo(1), hi(1)
               vel(i,j,k,1) = adv_vel
               vel(i,j,k,2) = zero
               vel(i,j,k,3) = zero
               scal(i,j,k,1) = one

               do n = 2,nscal-1
                  scal(i,j,k,n) = one
               end do                  

               x = xlo(1) + hx*(float(i-lo(1)) + half)
  	       dist = sqrt((x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2)
	       scal(i,j,k,nscal) = cvmgt(one,zero,dist.lt.radblob)
            end do
         end do
      end do

      end
c ::: -----------------------------------------------------------
c
      subroutine inithotspot(level,time,lo,hi,nscal,
     &	 	             vel,scal,DIMS(state),press,DIMS(press),
     &                       dx,xlo,xhi)
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


c     ::::: local variables
      integer i, j, k, n
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  x_vel, y_vel, z_vel
      REAL_T  dist

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      if (adv_dir .eq. 1) then
         x_vel = adv_vel
         y_vel = zero
         z_vel = zero
      else if (adv_dir .eq. 2) then
         x_vel = zero
         y_vel = adv_vel
         z_vel = zero
c     AJA - Strangely, this one did not
      else if (adv_dir .eq. 3) then
         x_vel = zero
         y_vel = zero
         z_vel = adv_vel
      else 
         write(6,*) "inithotspot: adv_dir = ",adv_dir
         stop
      end if

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
            do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half)
               dist = sqrt((x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2)
               vel(i,j,k,1) = x_vel
               vel(i,j,k,2) = y_vel
               vel(i,j,k,3) = z_vel
               scal(i,j,k,1) = one/denfact + (one - one/denfact)
     &              *half*(one + tanh(40.*(dist - radblob)))
               scal(i,j,k,2) = cvmgt(one,zero,dist.lt.radblob)
               do n = 3,nscal-1
                  scal(i,j,k,n) = one
               end do
               scal(i,j,k,nscal) = one / scal(i,j,k,1)
            end do
         end do
      end do
      
      end
c
c ::: -----------------------------------------------------------
c
      subroutine initeuler(level,time,lo,hi,nscal,
     &                       vel,scal,DIMS(state),press,DIMS(press),
     &                       dx,xlo,xhi)
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
c
c     ::::: local variables
c
      integer i, j, k, n
      REAL_T  x, y, z
      REAL_T  hx, hy, hz, r_yz

      REAL_T  eps_input, beta_input, rho_input
      REAL_T  delta_input, kappa_input

      parameter (eps_input=0.05, rho_input=0.15)
      parameter (beta_input=15., delta_input=0.0333)
      parameter (kappa_input=500.)

#include "probdata.H"
      
      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half) -half
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half) -half
            r_yz = sqrt(y*y+z*z)
            do i = lo(1), hi(1)

               x = xlo(1) + hx*(float(i-lo(1)) + half) -half
               vel(i,j,k,1) = tanh( (rho_input - r_yz) / delta_input)
               vel(i,j,k,2) = zero
               vel(i,j,k,3) = eps_input * exp(-beta_input * (x*x + y*y) )

               scal(i,j,k,1) = one
               scal(i,j,k,2) = exp( -kappa_input * (rho_input - r_yz)**2 )

               do n = 3,nscal
                  scal(i,j,k,n) = one
               end do

            end do
         end do
      end do

      end
c
c ::: -----------------------------------------------------------
c
      subroutine initbrnrsh(level,time,lo,hi,nscal,
     &                      vel,scal,DIMS(state),press,DIMS(press),
     &                      dx,xlo,xhi)
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
c
c     ::::: local variables
c
      integer i, j, k, n
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  dist
      REAL_T  x_vel, y_vel, z_vel
      REAL_T  umid,rmid,lamv,lamr,rfact


#include "probdata.H"
      hx = dx(1)
      hy = dx(2)
      hz = dx(3)


      umid = (vel1+vel2)*half
      rmid = (den1+den2)*half
      lamv = (vel1-vel2)/(vel1+vel2)
      lamr = (den1-den2)/(den1+den2)

         do k = lo(3), hi(3)
            z = xlo(3) + hz*(float(k-lo(3)) + half)
            do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
            do i = lo(1), hi(1)
               vel(i,j,k,1) = umid*(one+lamv*tanh(two*z/delta0))
               vel(i,j,k,2) = zero
               vel(i,j,k,3) = zero
               scal(i,j,k,1) = rmid*(one+lamr*tanh(two*z/delta0))
               scal(i,j,k,2) = cvmgt(one,zero,z.gt.zero)
            enddo
         enddo
         enddo


      end
c
c ::: -----------------------------------------------------------
c
      subroutine inflowt(level,time,lo,hi,nscal,
     &                       vel,scal,DIMS(state),press,DIMS(press),
     &                       dx,xlo,xhi)
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
c
c     ::::: local variables
c
      integer i, j, k, n
      REAL_T  x, y, z
      REAL_T  hx, hy, hz, r_yz

      REAL_T  eps_input, beta_input, rho_input
      REAL_T  delta_input, kappa_input


#include "probdata.H"
      
      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               vel(i,j,k,1) = zero
               vel(i,j,k,2) = zero
               vel(i,j,k,3) = adv_vel

               scal(i,j,k,1) = denfact
               scal(i,j,k,2) = zero

               do n = 3,nscal
                  scal(i,j,k,n) = one
               end do

            end do
         end do
      end do

      end
c
c ::: -----------------------------------------------------------
c
      subroutine taylorgreen(level,time,lo,hi,nscal,
     &	 	             vel,scal,DIMS(state),press,DIMS(press),
     &                       dx,xlo,xhi)
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
c
c     ::::: local variables
c
      integer i, j, k, n
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  dist, pi

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)
      pi = 4.d0*atan2(1.d0,1.d0)

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
	    do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half)

               vel(i,j,k,1) =  velfact*sin(x)*cos(y)*cos(z)
               vel(i,j,k,2) = -velfact*cos(x)*sin(y)*cos(z)

c               vel(i,j,k,1) = sin(x-pi*0.5d0)*cos(y-pi)*cos(z-pi*0.5d0)
c               vel(i,j,k,2) = -cos(x-pi*0.5d0)*sin(y-pi)*cos(z-pi*0.5d0)
               vel(i,j,k,3) = zero

               scal(i,j,k,1) = denfact
c     This is the theoretical pressure perturbation from p_0
               scal(i,j,k,2) = (denfact*velfact*velfact/16.d0)*(two+cos(two*z))*(cos(two*x)+cos(two*y))
               do n = 2,nscal-1
                  scal(i,j,k,n) = one
               end do
                  
c 	       dist = sqrt((x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2)
c	       scal(i,j,k,nscal) = cvmgt(one,zero,dist.lt.radblob)
	    end do
         end do
      end do

      end
c
c ::: -----------------------------------------------------------
c
      subroutine swirl(level,time,lo,hi,nscal,
     &	 	       vel,scal,DIMS(state),press,DIMS(press),
     &                 dx,xlo,xhi)
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
c
c     ::::: local variables
c
      integer i, j, k, n
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  dist

#include "probdata.H"

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

      end
c
c ::: -----------------------------------------------------------
c
      subroutine inithit(level,time,lo,hi,nscal,
     &	 	       vel,scal,DIMS(state),press,DIMS(press),
     &                 dx,xlo,xhi)
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
c
c     ::::: local variables
c
      integer i, j, k, n
      REAL_T  x, y, z
      REAL_T  f_xlo, f_ylo, f_zlo
      REAL_T  f_xhi, f_yhi, f_zhi
      REAL_T  Lmin
      REAL_T  hx, hy, hz
      REAL_T  dist
      REAL_T  twicePi, rn

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      twicePi = two*Pi

      f_xlo = f_problo(1)
      f_ylo = f_problo(2)
      f_zlo = f_problo(3)
      f_xhi = f_probhi(1)
      f_yhi = f_probhi(2)
      f_zhi = f_probhi(3)
      Lmin= min( f_xhi-f_xlo, f_yhi-f_ylo, f_zhi-f_zlo)

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         z = (z-f_zlo)/Lmin
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
            y = (y-f_ylo)/Lmin
	    do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half)
               x = (x-f_xlo)/Lmin

               if (turb_scale.gt.zero) then
                  vel(i,j,k,1) = turb_scale * cos(twicePi*y)*cos(twicePi*z)
                  vel(i,j,k,2) = turb_scale * cos(twicePi*x)*cos(twicePi*z)
                  vel(i,j,k,3) = turb_scale * cos(twicePi*x)*cos(twicePi*y)
               else
                  vel(i,j,k,1) = cos(twicePi*y)*cos(twicePi*z)
                  vel(i,j,k,2) = cos(twicePi*x)*cos(twicePi*z)
                  vel(i,j,k,3) = cos(twicePi*x)*cos(twicePi*y)
               endif

c               vel(i,j,k,1) = zero
c               vel(i,j,k,2) = zero
c               vel(i,j,k,3) = zero
               
               if (rhozero.gt.zero) then 
                  scal(i,j,k,1) = rhozero
               else
                  scal(i,j,k,1) = one
               endif
               do n = 2,nscal
                  if (abs(z-ref_centre).lt.ref_radius) then
                     scal(i,j,k,n) = one
                  else
                     scal(i,j,k,n) = zero
                  endif
               end do
                  
	    end do
         end do
      end do

      end
c
c ::: -----------------------------------------------------------
c
      subroutine initinflow(level,time,lo,hi,nscal,
     &	 	       vel,scal,DIMS(state),press,DIMS(press),
     &                 dx,xlo,xhi)
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
c
c     ::::: local variables
c
      integer i, j, k, n
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  dist
      REAL_T  x_vel, y_vel, z_vel
      REAL_T  twicePi, rn

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      twicePi = two*Pi

      x_vel = zero
      y_vel = zero
      z_vel = zero
      if (adv_dir.eq.1) then
         x_vel = adv_vel
      else if (adv_dir.eq.2) then
         y_vel = adv_vel
      else if (adv_dir.eq.3) then
         z_vel = adv_vel
      endif

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
	    do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half)

               vel(i,j,k,1) = x_vel
               vel(i,j,k,2) = y_vel
               vel(i,j,k,3) = z_vel

               scal(i,j,k,1) = one
               do n = 2,nscal
                  scal(i,j,k,n) = zero
               end do
                  
	    end do
         end do
      end do

      end
c
c ::: -----------------------------------------------------------
c
      subroutine initjet(level,time,lo,hi,nscal,
     &	 	       vel,scal,DIMS(state),press,DIMS(press),
     &                 dx,xlo,xhi)
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
c
c     ::::: local variables
c
      integer i, j, k, n
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  dist
      REAL_T  x_vel, y_vel, z_vel
      REAL_T  twicePi, rn

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
	    do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half)

               vel(i,j,k,1) = zero
               vel(i,j,k,2) = zero
               vel(i,j,k,3) = zero

               scal(i,j,k,1) = coflow_rho
               scal(i,j,k,2) = cvmgt(one,zero,abs(x-jet_x).le.jet_width .and. z.le.jet_width)

               if (nscal.ge.3) then 
                  do n = 3,nscal
                     scal(i,j,k,n) = zero
                  end do
               endif

	    end do
         end do
      end do

      end
c
c ::: -----------------------------------------------------------
c
      subroutine initinjection(level,time,lo,hi,nscal,
     &	 	       vel,scal,DIMS(state),press,DIMS(press),
     &                 dx,xlo,xhi)
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
c
c     ::::: local variables
c
      integer i, j, k, n
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  dist
      REAL_T  x_vel, y_vel, z_vel
      REAL_T  twicePi, rn

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
	    do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half)

               vel(i,j,k,1) = zero
               vel(i,j,k,2) = zero
               vel(i,j,k,3) = zero

               scal(i,j,k,1) = rhozero
               scal(i,j,k,2) = zero

               if (nscal.ge.3) then 
                  do n = 3,nscal
                     scal(i,j,k,n) = zero
                  end do
               endif

	    end do
         end do
      end do

      end
c
c ::: -----------------------------------------------------------
c
      subroutine initroundjet(level,time,lo,hi,nscal,
     &	 	       vel,scal,DIMS(state),press,DIMS(press),
     &                 dx,xlo,xhi)
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
c
c     ::::: local variables
c
      integer i, j, k, n
      REAL_T  x, y, z, r, r2
      REAL_T  hx, hy, hz
      REAL_T  dist
      REAL_T  x_vel, y_vel, z_vel
      REAL_T  twicePi, rn

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
	    do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half)

               r2= (x-jet_x)*(x-jet_x)+(y-jet_y)*(y-jet_y)
               r = sqrt(r2)

               vel(i,j,k,1) = zero
               vel(i,j,k,2) = zero
               vel(i,j,k,3) = zero

               scal(i,j,k,1) = coflow_rho
               scal(i,j,k,2) = cvmgt(one,zero,sqrt(r2+z*z).le.jet_width+hx)

               if (nscal.ge.3) then 
                  do n = 3,nscal
                     scal(i,j,k,n) = zero
                  end do
               endif

	    end do
         end do
      end do

      end
c
c ::: -----------------------------------------------------------
c
      subroutine initcoriolis(level,time,lo,hi,nscal,
     &	 	              vel,scal,DIMS(state),press,DIMS(press),
     &                        dx,xlo,xhi)
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
c
c     ::::: local variables
c
      integer i, j, k, n, m
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  Lx, Ly, Lz

      integer nx, mx
c     nx is the number of modes to use for velocity perturbation
      parameter (nx=10)
c     mx is the maximum wavenumber in each direction to choose random modes from
      parameter (mx=4)
c     a,b,c are random amplitudes in the range (0,1)
      REAL_T  a(nx,SDIM),b(nx,SDIM),c(nx,SDIM),rn
c     d,e,f are random modes in the range (1,mx)
      integer d(nx,SDIM),e(nx,SDIM),f(nx,SDIM)

      integer isioproc

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      Lx = f_probhi(1)-f_problo(1)
      Ly = f_probhi(1)-f_problo(1)
      Lz = min(ref_height-f_problo(1),f_probhi(1)-f_problo(1))

      call bl_pd_is_ioproc(isioproc)

      call blutilinitrand(blrandseed)
      do n = 1, nx
         do m = 1,3
            call blutilrand(rn)
            a(n,m) = rn
            call blutilrand(rn)
            b(n,m) = rn
            call blutilrand(rn)
            c(n,m) = rn
            call blutilrand(rn)
            d(n,m) = int(rn*mx+one)
            call blutilrand(rn)
            e(n,m) = int(rn*mx+one)
            call blutilrand(rn)
            f(n,m) = int(rn*mx+one)
         enddo
      enddo
      if (lo(1).eq.0.and.lo(2).eq.0.and.lo(3).eq.0.and.isioproc.eq.1) then
         write (6,*) "a",a
         write (6,*) "b",b
         write (6,*) "c",c
         write (6,*) "d",d
         write (6,*) "e",e
         write (6,*) "f",f
      endif

      do k = lo(3), hi(3)
         z = ( xlo(3) + hz*(float(k-lo(3)) + half) )
         do j = lo(2), hi(2)
            y = ( xlo(2) + hy*(float(j-lo(2)) + half) )
	    do i = lo(1), hi(1)
               x = ( xlo(1) + hx*(float(i-lo(1)) + half) )

               vel(i,j,k,1)  =  omega*y
               vel(i,j,k,2)  = -omega*x
               vel(i,j,k,3)  = zero

               if (ref_height.gt.zero.and.z.gt.ref_height) then
c     Maybe try a real free surface? Let's soften the blow with a hyperbolic tangent though
c     To run without the free surface, make sure ref_height is greater than the top of the box
                  scal(i,j,k,1) = half*rhozero*((one+denfact)-(one-denfact)*tanh((z-ref_height)/(two*delta0)))
                  scal(i,j,k,2) = zero
                  do n = 2,nscal-1
                     scal(i,j,k,n) = zero
                  enddo
               else
                  do n = 1, nx
                     do m = 1, SDIM
                        vel(i,j,k,m) = vel(i,j,k,m) + turb_scale*
     &                       a(n,m)*sin(two*Pi*d(n,m)*(x-f_problo(1))/Lx)*
     &                       b(n,m)*sin(two*Pi*e(n,m)*(y-f_problo(2))/Ly)*
     &                       c(n,m)*sin(two*Pi*f(n,m)*(z-f_problo(3))/Lz)
                     enddo
                  enddo
                  scal(i,j,k,1) = rhozero
                  scal(i,j,k,2) = zero
                  do n = 2,nscal-1
                     scal(i,j,k,n) = zero
                  end do
               endif
                  
	    end do
         end do
      end do

      end
c
c ::: -----------------------------------------------------------
c
      subroutine initheating(level,time,lo,hi,nscal,
     &	 	             vel,scal,DIMS(state),press,DIMS(press),
     &                       dx,xlo,xhi)
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
c
c     ::::: local variables
c
      integer i, j, k, n
      integer ii,jj,kk,nn
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  dist, pi

      REAL_T  xx, yy, zz
      REAL_T  hxx, hyy, hzz, ds

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      nn = 8
      ds = one/float(nn*nn*nn)
      hxx = hx/float(nn)
      hyy = hy/float(nn)
      hzz = hz/float(nn)

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)))
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)))
	    do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)))

               vel(i,j,k,1) = zero
               vel(i,j,k,2) = zero
               vel(i,j,k,3) = zero

               scal(i,j,k,1) = rhozero
               scal(i,j,k,2) = zero
               scal(i,j,k,3) = zero

               do kk = 1, nn
                  zz = z + hzz*(float(kk) - half)
                  do jj = 1, nn
                     yy = y + hyy*(float(jj) - half)
                     do ii = 1, nn
                        xx = x + hxx*(float(ii) - half)
                        if (sqrt((xx-0.25)*(xx-0.25)+yy*yy+(zz+0.25)*(zz+0.25)).lt.0.1) then
#ifdef MOREGENGETFORCE
                           scal(i,j,k,2) = scal(i,j,k,2) + ds
#else
                           scal(i,j,k,1) = scal(i,j,k,1) - ds
#endif
                        endif
                        if (sqrt((xx+0.25)*(xx+0.25)+yy*yy+(zz+0.25)*(zz+0.25)).lt.0.1) then
                           scal(i,j,k,3) = scal(i,j,k,3) + ds
                        endif                        
                     enddo
                  enddo
               enddo
	    end do
         end do
      end do

      end
c
c
c ::: -----------------------------------------------------------
c
c     This routine averages the mac face velocities for makeforce at half time
c
      subroutine FORT_AVERAGE_EDGE_STATES(vel,umacx,umacy,umacz,
     &                                    DIMS(vel),DIMS(umacx),DIMS(umacy),DIMS(umacz),
     &                                    getForceVerbose)

      implicit none

      integer    DIMDEC(vel)
      integer    DIMDEC(umacx)
      integer    DIMDEC(umacy)
      integer    DIMDEC(umacz)
      integer    getForceVerbose
      REAL_T     vel  (DIMV(vel),SDIM)
      REAL_T     umacx(DIMV(umacx))
      REAL_T     umacy(DIMV(umacy))
      REAL_T     umacz(DIMV(umacz))

      integer i,j,k,n
      integer ilo,jlo,klo
      integer ihi,jhi,khi

      integer isioproc

      REAL_T  velmin(3)
      REAL_T  velmax(3)

      do n = 1, 3
         velmin(n) = 1.d234
         velmax(n) = -1.d234
      enddo

      ilo = vel_l1
      jlo = vel_l2
      klo = vel_l3
      ihi = vel_h1
      jhi = vel_h2
      khi = vel_h3

      do k = klo, khi
         do j = jlo, jhi
            do i = ilo, ihi
               vel(i,j,k,1) = half*(umacx(i,j,k)+umacx(i+1,j,k))
               vel(i,j,k,2) = half*(umacy(i,j,k)+umacy(i,j+1,k))
               vel(i,j,k,3) = half*(umacz(i,j,k)+umacz(i,j,k+1))
               do n=1, 3
                  velmin(n)=min(velmin(n),vel(i,j,k,n))
                  velmax(n)=max(velmax(n),vel(i,j,k,n))
               enddo
            enddo
         enddo
      enddo
      
      if (getForceVerbose.gt.0) then
         call bl_pd_is_ioproc(isioproc)
         if (isioproc) then
            do n = 1, 3
               write (6,*) "mac velmin (",n,") = ",velmin(n)
               write (6,*) "mac velmax (",n,") = ",velmax(n)
            enddo
         endif
      endif

      end
c
c
c ::: -----------------------------------------------------------
c
c     This routine add the forcing terms to the momentum equation
c
      subroutine FORT_MAKEFORCE(time,force,
#ifdef MOREGENGETFORCE
     &                          vel,
#endif
     &                          scal,
     &                          DIMS(force),
#ifdef MOREGENGETFORCE
     &                          DIMS(vel),
#endif
     &                          DIMS(scal),
     &                          dx,xlo,xhi,gravity,scomp,ncomp
#ifdef MOREGENGETFORCE
     &,                         nscal,getForceVerbose
#endif
     &)

      implicit none

      integer    DIMDEC(force)
      integer    DIMDEC(scal)
      integer    scomp, ncomp
      REAL_T     time, dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     force  (DIMV(force),scomp:scomp+ncomp-1)
      REAL_T     gravity
#ifdef MOREGENGETFORCE
      integer    DIMDEC(vel)
      integer    nscal, getForceVerbose
      REAL_T     vel    (DIMV(vel),0:SDIM-1)
      REAL_T     scal   (DIMV(scal),0:nscal-1)
#else
      REAL_T     scal   (DIMV(scal),0:0)
#endif

      REAL_T     Lx, Ly, Lz, Lmin, HLx, HLy, HLz
      REAL_T     kappa, kappaMax

#include "probdata.H"

#ifdef DO_IAMR_FORCE
#include "forcedata.H"
#endif

c
c     ::::: local variables
c
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
      integer kx, ky, kz, mode_count, xstep, ystep, zstep
      integer isioproc, do_trac2
      integer nXvel, nYvel, nZvel, nRho, nTrac, nTrac2, nRhoScal, nTracScal, nTrac2Scal

#ifdef MOREGENGETFORCE
      REAL_T  velmin(0:SDIM-1)
      REAL_T  velmax(0:SDIM-1)
      REAL_T  scalmin(0:nscal-1)
      REAL_T  scalmax(0:nscal-1)
      REAL_T  forcemin(scomp:scomp+ncomp-1)
      REAL_T  forcemax(scomp:scomp+ncomp-1)
#endif

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

c     Assumes components are in the following order
      nXvel  = 0
      nYvel  = 1
      nZvel  = 2
      nRho   = 3
      nTrac  = 4
      nTrac2 = 5

      nRhoScal   = nRho-SDIM
      nTracScal  = nTrac-SDIM
      nTrac2Scal = nTrac2-SDIM

#ifdef MOREGENGETFORCE
      if (getForceVerbose.gt.0) then
         call bl_pd_is_ioproc(isioproc)
         if (isioproc) then

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
            
c     Get min/max
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
c     Velocities
                     do n = 0, SDIM-1
                        if (vel(i,j,k,n).gt.velmax(n)) then
                           velmax(n)=vel(i,j,k,n)
                        endif
                        if (vel(i,j,k,n).lt.velmin(n)) then
                           velmin(n)=vel(i,j,k,n)
                        endif
                     enddo
c     Scalars
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
#endif

c
c     Here's where the forcing actually gets done
c
      
      if (scomp.eq.0) then
c
c     Do velocity forcing
c
         if (probtype.eq.20) then 
            do k = klo, khi
               z = zlo + hz*(float(k-klo) + half)
               do j = jlo, jhi
                  y = xlo(2) + hy*(float(j-jlo) + half)
                  do i = ilo, ihi
                     x = xlo(1) + hx*(float(i-ilo) + half)
                     force(i,j,k,nXvel) = zero
                     force(i,j,k,nYvel) = zero
                     force(i,j,k,nZvel) = zero
#ifdef MOREGENGETFORCE
                     if (do_trac2) then
                        force(i,j,k,nZvel) = force(i,j,k,nZvel) + thermal_expansion*scal(i,j,k,nTrac2Scal)
                     endif
#endif
                  enddo
               enddo
            enddo
         else if (probtype.eq.18) then 
c     Round jet/plume
            do k = klo, khi
               z = zlo + hz*(float(k-klo) + half)
               do j = jlo, jhi
                  y = xlo(2) + hy*(float(j-jlo) + half)
                  do i = ilo, ihi
                     x = xlo(1) + hx*(float(i-ilo) + half)
                     force(i,j,k,nXvel) = zero
                     force(i,j,k,nYvel) = zero
                     force(i,j,k,nZvel) = gravity*scal(i,j,k,nRhoScal)
#ifdef MOREGENGETFORCE
                     if (do_trac2.eq.1) then
                        force(i,j,k,nZvel) = force(i,j,k,nZvel) + thermal_expansion*scal(i,j,k,nTrac2Scal)
                     endif
#endif
                  enddo
               enddo
            enddo
         else if (probtype.eq.14.or.probtype.eq.15) then
#ifdef DO_IAMR_FORCE
c     Homogeneous Isotropic Turbulence
            twicePi=two*Pi
            
c     Adjust z offset for probtype 15
            if (probtype.eq.15.and.infl_time_offset.gt.(-half)) then
               infl_time = time + infl_time_offset
               zlo = xlo(3) - (time*adv_vel)
            else
               infl_time = time
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
                                 xT = cos(FTX(kx,ky,kz)*infl_time+TAT(kx,kx,kz))
                                 f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                                 f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                                 f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) * cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                              endif
                           enddo
                        enddo
                     enddo
                     do kz = mode_start, zstep - 1
                        kzd = dfloat(kz)
                        do ky = mode_start, nymodes
                           kyd = dfloat(ky)
                           do kx = mode_start, nxmodes
                              kxd = dfloat(kx)
                              kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                              if (kappa.le.kappaMax) then
                                 xT = cos(FTX(kx,ky,kz)*infl_time+TAT(kx,kx,kz))
                                 f1 = f1 + xT*FAX(kx,ky,kz)*cos(twicePi*kxd*x/Lx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/Ly+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/Lz+FPZ(kx,ky,kz))
                                 f2 = f2 + xT*FAY(kx,ky,kz)*sin(twicePi*kxd*x/Lx+FPX(kx,ky,kz)) * cos(twicePi*kyd*y/Ly+FPY(kx,ky,kz)) * sin(twicePi*kzd*z/Lz+FPZ(kx,ky,kz))
                                 f3 = f3 + xT*FAZ(kx,ky,kz)*sin(twicePi*kxd*x/Lx+FPX(kx,ky,kz)) * sin(twicePi*kyd*y/Ly+FPY(kx,ky,kz)) * cos(twicePi*kzd*z/Lz+FPZ(kx,ky,kz))
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
#ifdef MOREGENGETFORCE
c     Coriolis
            do k = klo, khi
               z = zlo + hz*(float(k-klo) + half)
               do j = jlo, jhi
                  y = xlo(2) + hy*(float(j-jlo) + half)
                  do i = ilo, ihi
                     x = xlo(1) + hx*(float(i-ilo) + half)
                     force(i,j,k,nXvel) = scal(i,j,k,nRhoScal)*( two*omega*vel(i,j,k,nYvel)+omega*omega*x)
                     force(i,j,k,nYvel) = scal(i,j,k,nRhoScal)*(-two*omega*vel(i,j,k,nXvel)+omega*omega*y)
                     force(i,j,k,nZvel) = gravity*scal(i,j,k,nRhoScal)
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
         else if (probtype.eq.99.and.abs(grav_angle).gt.0.001) then
c     Angled gravity
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
c     Default to gravity...
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
c     else to zero
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
c     End of velocity forcing
      endif

      if ((scomp+ncomp).gt.BL_SPACEDIM) then
c
c     Scalar forcing
c
         do n = max(scomp,nRho), scomp+ncomp-1
            if (n.eq.nRho) then
c
c     Density
c
               do k = klo, khi
                  do j = jlo, jhi
                     do i = ilo, ihi
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            else if (n.eq.nTrac) then
c
c     Tracer
c
               do k = klo, khi
                  do j = jlo, jhi
                     do i = ilo, ihi
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            else if (n.eq.nTrac2.and.do_trac2.eq.1) then
c
c     Other scalar
c
               if (probtype.eq.20) then 
c     Temperature perturbation
                  do k = klo, khi
                     do j = jlo, jhi
                        do i = ilo, ihi
#ifdef MOREGENGETFORCE
                           force(i,j,k,n) = heating_coeff * scal(i,j,k,nTracScal)
#else
                           force(i,j,k,n) = zero
#endif
                        enddo
                     enddo
                  enddo
               else  if (probtype.eq.18) then 
c     Round Jet/Plume
                  do k = klo, khi
                     z = zlo + hz*(float(k-klo) + half)
                     do j = jlo, jhi
                        do i = ilo, ihi
#ifdef MOREGENGETFORCE
                           if (abs(z-heating_centre).lt.heating_radius) then 
                              force(i,j,k,n) = heating_coeff * scal(i,j,k,nTracScal)
                           else
                              force(i,j,k,n) = zero
                           endif
#else
                           force(i,j,k,n) = zero
#endif
                        enddo
                     enddo
                  enddo
               else
c     Some other probtype
                  do k = klo, khi
                     do j = jlo, jhi
                        do i = ilo, ihi
                           force(i,j,k,n) = zero
                        enddo
                     enddo
                  enddo
               endif
            else
c
c     Other scalar
c
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
      
#ifdef MOREGENGETFORCE
      if (getForceVerbose.gt.0 .and. isioproc) then
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
#endif

      end

c     ::: -----------------------------------------------------------
c     ::: This routine will tag high error cells based on the 
c     ::: density gradient
c     ::: 
c     ::: INPUTS/OUTPUTS:
c     ::: 
c     ::: tag      <=  integer tag array
c     ::: DIMS(tag) => index extent of tag array
c     ::: set       => integer value to tag cell for refinement
c     ::: clear     => integer value to untag cell
c     ::: rho       => density array
c     ::: DIMS(rho) => index extent of rho array
c     ::: nvar      => number of components in rho array (should be 1)
c     ::: lo,hi     => index extent of grid
c     ::: domlo,hi  => index extent of problem domain
c     ::: dx        => cell spacing
c     ::: xlo       => physical location of lower left hand
c     :::	           corner of tag array
c     ::: problo    => phys loc of lower left corner of prob domain
c     ::: time      => problem evolution time
c     ::: -----------------------------------------------------------
      subroutine FORT_DENERROR (tag,DIMS(tag),set,clear,
     &     rho,DIMS(rho),lo,hi,nvar,
     &     domlo,domhi,dx,xlo,
     &     problo,time,level)
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

#include "probdata.H"

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
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),rho(i,j,k,1).gt.denerr)
               end do
            end do
         end do
      else
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),rho(i,j,k,1).lt.denerr)
               end do
            end do
         end do
      endif

      end

c ::: -----------------------------------------------------------
c ::: This routine will tag high error cells based on the 
c ::: density gradient
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: tag      <=  integer tag array
c ::: DIMS(tag) => index extent of tag array
c ::: set       => integer value to tag cell for refinement
c ::: clear     => integer value to untag cell
c ::: adv       => scalar array
c ::: DIMS(adv) => index extent of adv array
c ::: nvar      => number of components in rho array (should be 1)
c ::: lo,hi     => index extent of grid
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of tag array
c ::: problo    => phys loc of lower left corner of prob domain
c ::: time      => problem evolution time
c ::: -----------------------------------------------------------
      subroutine FORT_ADVERROR (tag,DIMS(tag),set,clear,
     &                          adv,DIMS(adv),lo,hi,nvar,
     &                          domlo,domhi,dx,xlo,
     &			        problo,time,level)
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

#include "probdata.H"

c     probtype = SPIN
      if (probtype .eq. 1) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = cvmgt(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

c     probtype = BUBBLE
      else if (probtype .eq. 2) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = cvmgt(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

c     probtype = VORTEX IN A BOX
      else if (probtype .eq. 3) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = cvmgt(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

c     probtype = CHANNEL
      else if (probtype .eq. 4) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = cvmgt(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

c     probtype = PERIODIC SHEAR LAYER
      else if (probtype .eq. 5) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = cvmgt(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

c     probtype = HOT SPOT
      else if (probtype .eq. 6) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = cvmgt(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
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
                  tag(i,j,k) = cvmgt(set,clear,abs(z).lt.zlev1.and.x.lt.xlev1)
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
            tag(i,j,k) = cvmgt(set,tag(i,j,k),aerr.gt.adverr.and. x.lt.xlev1)
          enddo
          enddo
          enddo
        endif

c     probtype = VISCOUS BENCHMARK
      else if (probtype .eq. 9) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = cvmgt(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do


      else if (probtype .eq. 10.or.probtype.eq.11) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = cvmgt(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

      else if (probtype .eq. 12) then

      else if (probtype .eq. 13) then

      else if (probtype .eq. 14) then

         do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = cvmgt(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

      else if (probtype .eq. 15) then

      else if (probtype .eq. 16) then
 
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = cvmgt(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

      else if (probtype .eq. 17) then
 
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = cvmgt(set,tag(i,j,k),abs(adv(i,j,k,1)).gt.adverr)
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
 
        do k = lo(3), hi(3)
           z = xlo(3) + dx(3)*(float(k-lo(3)) + half)
           if (z.lt.ref_height) then
              do j = lo(2), hi(2)
                 do i = lo(1), hi(1)
                    tag(i,j,k) = cvmgt(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
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
                    tag(i,j,k) = cvmgt(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
                 end do
              end do
           endif
        end do

      else if (probtype .eq. 20) then
 
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = cvmgt(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do
        
      else
        print *,'DONT KNOW THIS PROBTYPE IN FORT_ADVERROR ',probtype
        stop
      end if
 
      end

      subroutine FORT_ADV2ERROR (tag,DIMS(tag),set,clear,
     &                           adv,DIMS(adv),lo,hi,nvar,
     &                           domlo,domhi,dx,xlo,
     &			         problo,time,level)
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

#include "probdata.H"

      if (probtype .eq. 20) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
               end do
            end do
         end do
      else
         print *,'DONT KNOW THIS PROBTYPE IN FORT_ADVERROR ',probtype
         stop
      end if
      
      end

c ::: -----------------------------------------------------------
c ::: This routine will tag high error cells based on the
c ::: temperature gradient
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: tag      <=  integer tag array
c ::: DIMS(tag) => index extent of tag array
c ::: set       => integer value to tag cell for refinement
c ::: clear     => integer value to untag cell
c ::: temp      => density array
c ::: DIMS(temp)=> index extent of temp array
c ::: lo,hi     => index extent of grid
c ::: nvar      => number of components in rho array (should be 1)
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::              corner of tag array
c ::: problo    => phys loc of lower left corner of prob domain
c ::: time      => problem evolution time
c ::: -----------------------------------------------------------
      subroutine FORT_TEMPERROR (tag,DIMS(tag),set,clear,
     &                           temperature,DIMS(temp),lo,hi,nvar,
     &                           domlo,domhi,dx,xlo,
     &                           problo,time,level)
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

#include "probdata.H"

c     probtype = SPIN
      if (probtype .eq. 1) then

c     probtype = BUBBLE
      else if (probtype .eq. 2) then

c     probtype = VORTEX IN A BOX
      else if (probtype .eq. 3) then

c     probtype = CHANNEL
      else if (probtype .eq. 4) then

c     probtype = PERIODIC SHEAR LAYER
      else if (probtype .eq. 5) then

c     probtype = HOT SPOT
      else if (probtype .eq. 6) then

        if (level .eq. 0) then
c         ::::: refine around entire hot spot
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   tag(i,j,k) = cvmgt(set,tag(i,j,k),temperature(i,j,k,1).gt.temperr)
                end do
             end do
          end do
        else
c         ::::: refine where there is temperature gradient
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   ax = abs(temperature(i+1,j,k,1) - temperature(i-1,j,k,1))
                   ay = abs(temperature(i,j+1,k,1) - temperature(i,j-1,k,1))
                   az = abs(temperature(i,j,k+1,1) - temperature(i,j,k-1,1))
                   aerr = max(ax,ay,az)
                   tag(i,j,k) = cvmgt(set,tag(i,j,k),aerr.gt.bubgrad)
                end do
             end do
          end do
        end if

c     probtype = VISCOUS BENCHMARK
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

      else
        print *,'DONT KNOW THIS PROBTYPE IN FORT_TEMPERROR ',probtype
        stop
      end if

      end
c ::: -----------------------------------------------------------
c ::: This routine will tag high error cells based on the 
c ::: magnitude of vorticity
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: tag      <=  integer tag array
c ::: DIMS(tag) => index extent of tag array
c ::: set       => integer value to tag cell for refinement
c ::: clear     => integer value to untag cell
c ::: vort      => vorticitiy array
c ::: DIMS(vort)=> index extent of vort array
c ::: nvar      => number of components in rho array (should be 1)
c ::: lo,hi     => index extent of grid
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of tag array
c ::: problo    => phys loc of lower left corner of prob domain
c ::: time      => problem evolution time
c ::: -----------------------------------------------------------
      subroutine FORT_MVERROR (tag,DIMS(tag),set,clear,
     &                         vort,DIMS(vort),lo,hi,nvar,
     &                         domlo,domhi,dx,xlo,
     &			       problo,time,level)
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

#include "probdata.H"

c     probtype = SPIN
      if (probtype .eq. 1) then

c     probtype = BUBBLE
      else if (probtype .eq. 2) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

c     probtype = VORTEX IN A BOX
      else if (probtype .eq. 3) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

c     probtype = CHANNEL
      else if (probtype .eq. 4) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

c     probtype = PERIODIC SHEAR LAYER
      else if (probtype .eq. 5) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

c     probtype = HOT SPOT
      else if (probtype .eq. 6) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do
c     probtype = EULER
      else if (probtype .eq. 7) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),
     &                   abs(dx(1)*vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

c     probtype = BROWN ROSHKO
      else if (probtype .eq. 8) then

        return


c     probtype = VISCOUS BENCHMARK
      else if (probtype .eq. 9) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

      else if (probtype .eq. 10 .or. probtype .eq. 11 .or.probtype.eq.12) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
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
                     tag(i,j,k) = cvmgt(set,tag(i,j,k),ztag.eq.1)
                  end do
               end do
            end do
         else if (ref_type.eq.2) then
            do k = lo(3), hi(3)
               z = xlo(3) + dx(3)*(float(k-lo(3)) + half)
               if (abs(z-ref_centre).lt.ref_radius) then
                  ztag = 1
               else
                  ztag = 0
               endif
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tag(i,j,k) = cvmgt(set,tag(i,j,k),ztag.eq.1)
                  end do
               end do
            end do
         else
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tag(i,j,k) = cvmgt(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
                  end do
               end do
            end do
         endif

      else if (probtype .eq. 15) then

      else if (probtype .eq. 16) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

      else if (probtype .eq. 17) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

      else if (probtype .eq. 18) then

         do k = lo(3), hi(3)
            z = xlo(3) + dx(3)*(float(k-lo(3)) + half)
            if (z.lt.ref_height) then
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tag(i,j,k) = cvmgt(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
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
                     tag(i,j,k) = cvmgt(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
                  end do
               end do
            endif
         end do

      else
        print *,'DONT KNOW THIS PROBTYPE IN FORT_MVERROR ',probtype
        stop
      end if

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: rho      <=  density array
c ::: DIMS(rho) => index extent of rho array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_DENFILL (rho,DIMS(rho),domlo,domhi,dx,
     &                         xlo,time,bc )
      implicit none

      integer    DIMDEC(rho)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     rho(DIMV(rho))
      integer    bc(SDIM,2)
      integer    lo(SDIM), hi(SDIM)

      integer    i, j, k
      REAL_T  umid,rmid,lamv,lamr,rfact
      REAL_T  hx, hy, hz, gpert, z, y, x, r, ypert, magwif, constn, rs, rd

      parameter (constn=.22089323)

#include "probdata.H"

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
                 rho(i,j,k) = den2
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

      else


      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(rho).lt.domlo(1)) then
         do i = ARG_L1(rho), domlo(1)-1
            do k = ARG_L3(rho), ARG_H3(rho)
               do j = ARG_L2(rho), ARG_H2(rho)
                  rho(i,j,k) = denfact
               end do
	    end do
	 end do
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
c           write (6,*) "Filling with unit density"
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
                    rho(i,j,k) = cvmgt(jet_rho,coflow_rho,abs(x-jet_x).le.jet_width)
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
                    r = sqrt((x-jet_x)*(x-jet_x)+(y-jet_y)*(y-jet_y))
                    rho(i,j,k) = half*(rs-rd*tanh(two*(r-jet_width)/delta0))
                 enddo
              enddo
           enddo
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

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: adv      <=  advected quantity array
c ::: DIMS(adv) => index extent of adv array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of adv array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_ADVFILL (adv,DIMS(adv),domlo,domhi,dx,
     &                         xlo,time,bc )
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

#include "probdata.H"

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
               adv(i,j,k) = cvmgt(one,zero,(z-ypert).gt.zero)
            enddo
            enddo
         enddo
      endif

      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               adv(i,j,k) = cvmgt(one,zero,z.gt.zero)
            enddo
            enddo
         enddo
      endif

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            do k = lo(3), hi(3)
            z = xlo(3) + hz*(float(k-lo(3)) + half)
            do i = lo(1), hi(1)
               adv(i,j,k) =  cvmgt(one,zero,z.gt.zero)
            enddo
            enddo
         enddo
      endif

      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            do k = lo(3), hi(3)
            z = xlo(3) + hz*(float(k-lo(3)) + half)
            do i = lo(1), hi(1)
               adv(i,j,k) =  cvmgt(one,zero,z.gt.zero)
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
         if (probtype.eq.16) then
            do k = ARG_L3(adv), domlo(3)-1
               do j = ARG_L2(adv), ARG_H2(adv)
                  y = xlo(2) + dx(2)*(float(j-lo(2)) + half)
                  do i = ARG_L1(adv), ARG_H1(adv)
                     x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                     r = sqrt((x-jet_x)*(x-jet_x)+(y-jet_y)*(y-jet_y))
                     adv(i,j,k) = cvmgt(one,zero,r.le.jet_width)
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
                        adv(i,j,k) = cvmgt(one,zero,r.le.jet_width)
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
                     r = sqrt((x-jet_x)*(x-jet_x)+(y-jet_y)*(y-jet_y))
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

      end

      subroutine FORT_ADV2FILL (adv,DIMS(adv),domlo,domhi,dx,
     &                         xlo,time,bc )
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

#include "probdata.H"

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
                     r = sqrt((x-jet_x)*(x-jet_x)+(y-jet_y)*(y-jet_y))
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
                        adv(i,j,k) = cvmgt(one,zero,r.le.jet_width)
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

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: temp      <= temperature array
c ::: DIMS(temp)=> index extent of temp array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of temp array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_TEMPFILL (temp,DIMS(temp),domlo,domhi,dx,
     &                          xlo,time,bc )
      implicit none

      integer    DIMDEC(temp)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     temp(DIMV(temp))
      integer    bc(SDIM,2)
      integer    lo(SDIM), hi(SDIM)

      integer    i, j, k
      REAL_T     x, y, z, r

#include "probdata.H"

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
                     r = sqrt((x-jet_x)*(x-jet_x)+(y-jet_y)*(y-jet_y))
                     temp(i,j,k) = cvmgt(one,zero,r.le.jet_width)
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

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: u        <=  x velocity array
c ::: DIMS(u)   => index extent of u array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_XVELFILL (u,DIMS(u),domlo,domhi,dx,xlo,time,bc)
      implicit none
      integer    DIMDEC(u)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     u(DIMV(u))
      integer    bc(SDIM,2)
      integer    lo(SDIM),hi(SDIM)
      integer    i, j, k
      REAL_T     x_vel
      REAL_T     y, ul, g_t
      REAL_T     umid,rmid,lamv,lamr,rfact
      REAL_T     hx, hy, hz, gpert, z,  ypert, magwif, constn
      REAL_T     jv, s, t

#include "probdata.H"

#ifdef BL_DO_FLCT
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      integer DIMDEC(uflct)
      REAL_T  t_flct
      REAL_T, allocatable :: uflct(:,:,:)
#include "INFL_FORCE_F.H"
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
c
c        Note that we are 'scaling time' here to step into the fluct file to the
c        correct depth.  This requires that time is not further scaled inside the
c        the INFL_FILL routine.  Just to be sure, we set convVel = 1 here again.
c
         convVel = one
         t_flct = adv_vel*time
         call INFL_FILL(FLCT_XVEL,DIMS(uflct),uflct,xlo,dx,t_flct,bc,domnlo,domnhi)
      end if
#endif

      if (adv_dir .eq. 1) then
         x_vel = adv_vel
      else  
         x_vel = zero
      end if

      call filcc(u,DIMS(u),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(u).lt.domlo(1)) then
         do i = ARG_L1(u), domlo(1)-1
            do k = ARG_L3(u), ARG_H3(u)
               do j = ARG_L2(u), ARG_H2(u)
                  u(i,j,k) = x_vel
               end do
	    end do
	 end do
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
                     r = sqrt( (x-jet_x)*(x-jet_x) + (y-jet_y)*(y-jet_y) )
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
         else 


            if (probtype.eq.16) then
               t = time*jet_vel/jet_width
            endif

            do k = ARG_L3(u), domlo(3)-1
               do j = ARG_L2(u), ARG_H2(u)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(u), ARG_H1(u)

                     if (probtype.eq.13) then
#ifdef BL_DO_FLCT
                        x = domnlo(1) + (i+0.5)*dx(1)
                        r = SQRT( x*x + y*y )
c     call vswirlXYZ(x,y,u1,u2,u3)
                        if (r.gt.Rfu) then
                           u_inf = 0.d0
                           eta = TANH(2*(r-Rfu)/Rtran)
                           u(i,j,k) = eta*u_inf + (1.d0-eta)*u1
                           uflct(i,j,1) = 0.d0
                        else
                           u(i,j,k) = u1
                        endif
#endif
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

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: v        <=  y velocity array
c ::: DIMS(v)   => index extent of v array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_YVELFILL (v,DIMS(v),domlo,domhi,dx,xlo,time,bc)
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
      REAL_T     jv, s, t

#include "probdata.H"

#ifdef BL_DO_FLCT
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      integer DIMDEC(vflct)
      REAL_T  t_flct
      REAL_T, allocatable :: vflct(:,:,:)
#include "INFL_FORCE_F.H"
#endif
      REAL_T x,y,r,u1,u2,u3,u_inf,eta

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
         t_flct = adv_vel*time
         call INFL_FILL(FLCT_YVEL,DIMS(vflct),vflct,xlo,dx,t_flct,bc,domnlo,domnhi)
      end if
#endif

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
         do j = ARG_L2(v), domlo(2)-1
            do k = ARG_L3(v), ARG_H3(v)
               do i = ARG_L1(v), ARG_H1(v)
                  v(i,j,k) = y_vel
               end do
            end do
	 end do
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
                     r = sqrt( (x-jet_x)*(x-jet_x) + (y-jet_y)*(y-jet_y) )
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
         else 
            if (probtype.eq.16) then
               t = time*jet_vel/jet_width
            endif

            do k = ARG_L3(v), domlo(3)-1
               do j = ARG_L2(v), ARG_H2(v)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(v), ARG_H1(v)

                     if (probtype.eq.13) then
#ifdef BL_DO_FLCT
                        x = domnlo(1) + (i+0.5)*dx(1)
                        r = SQRT( x*x + y*y )
c     call vswirlXYZ(x,y,u1,u2,u3)
                        if (r.gt.Rfu) then
                           u_inf = 0.d0
                           eta = TANH(2*(r-Rfu)/Rtran)
                           v(i,j,k) = eta*u_inf + (1.d0-eta)*u2
                           vflct(i,j,1) = 0.d0
                        else
                           v(i,j,k) = u2
                        endif
#endif
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

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: w        <=  z velocity array
c ::: DIMS(w)   => index extent of v array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_ZVELFILL (w,DIMS(w),domlo,domhi,dx,xlo,time,bc)
      implicit none
      integer    DIMDEC(w)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     w(DIMV(w))
      integer    bc(SDIM,2)
      integer    lo(SDIM),hi(SDIM)
      integer    i, j, k
      REAL_T     z_vel
      REAL_T     jv, s, t

#include "probdata.H"

#ifdef BL_DO_FLCT
      REAL_T t_flct
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      integer DIMDEC(wflct)
      REAL_T, allocatable :: wflct(:,:,:)
#include "INFL_FORCE_F.H"
#endif
      REAL_T x,y,r,u1,u2,u3,u_inf,eta

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
         t_flct = adv_vel*time
         call INFL_FILL(FLCT_ZVEL,DIMS(wflct),wflct,xlo,dx,t_flct,bc,domnlo,domnhi)
      end if
#endif

      if ( (time .gt. zero).and.(tVco_l.ge.zero) ) then
         if (time .le. tVco_l) then
            Vco = Vco_l
         else if (time .ge. tVco_r) then
            Vco = Vco_r
         else
            Vco = Vco_l+(time-tVco_l)*(Vco_r-Vco_l)/(tVco_r-tVco_l)
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

         else if (probtype.eq.18) then
            do k = ARG_L3(w), domlo(3)-1
               do j = ARG_L2(w), ARG_H2(w)
                  y = domnlo(2) + (j+0.5)*dx(2)
                  do i = ARG_L1(w), ARG_H1(w)
                     x = domnlo(1) + (i+0.5)*dx(1)
                     r = sqrt( (x-jet_x)*(x-jet_x) + (y-jet_y)*(y-jet_y) )
                     w(i,j,k) = half*jet_vel*(one-tanh(two*(r-jet_width)/delta0))
#ifdef BL_DO_FLCT
                     if (r.lt.jet_width) then
                        w(i,j,k) = w(i,j,k)*(one+wflct(i,j,1)*turb_scale)
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

                     if (probtype.eq.13) then
#ifdef BL_DO_FLCT
                        x = domnlo(1) + (i+0.5)*dx(1)
                        r = SQRT( x*x + y*y )
c     call vswirlXYZ(x,y,u1,u2,u3)
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

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: divu     <=  divergence of velocity array
c ::: DIMS(divu)=> index extent of divu array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_DIVUFILL (divu,DIMS(divu),domlo,domhi,dx,
     &                          xlo,time,bc)
      implicit none

      integer    DIMDEC(divu)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     divu(DIMV(divu))
      integer    bc(SDIM,2)

      integer    i, j, k
      REAL_T     z_vel

#include "probdata.H"

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

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: dsdt     <=  dsdt array
c ::: DIMS(dsdt)=> index extent of dsdt array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_DSDTFILL (dsdt,DIMS(dsdt),domlo,domhi,dx,
     &                          xlo,time,bc)
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

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: p        <=  pressure array
c ::: lo,hi     => index extent of p array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_PRESFILL (p,DIMS(p),domlo,domhi,dx,
     &                          xlo,time,bc)
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

#include "probdata.H"

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

c*****************************************************************************
c SETTING XLO 
c*****************************************************************************

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

c*****************************************************************************
c SETTING XHI
c*****************************************************************************

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

c*****************************************************************************
c SETTING YLO
c*****************************************************************************

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
 
c*****************************************************************************
c SETTING YHI
c*****************************************************************************

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

c*****************************************************************************
c SETTING ZLO
c*****************************************************************************

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

c*****************************************************************************
c SETTING ZHI
c*****************************************************************************

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

c*****************************************************************************

      end
