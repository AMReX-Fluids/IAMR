
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


      block data rt_common
#include <probdata.H>
      data rt_pertamp/0.0D0/
      data rt_nfreq/0/
      data rt_xfrontw/0.0d0/
      end  
    
module prob_3D_module

  implicit none

  private

  public :: amrex_probinit, FORT_INITDATA, initfromrest, &
            initrt, initpervort, &
            initviscbench, initchannel, inithotspot, &
            initeuler, taylorgreen, &
            initrayleightaylor, initraneul, initdenadvect, &
            initshearlayer, &
            FORT_DSDTFILL, &
            FORT_ADVERROR, FORT_TEMPERROR, FORT_MVERROR, &
            FORT_DENFILL, FORT_ADVFILL, FORT_TEMPFILL, FORT_XVELFILL, &
            FORT_YVELFILL, FORT_ZVELFILL, FORT_PRESFILL, FORT_DIVUFILL

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

      INTEGER dimFile(3)
      integer nCompFile
      REAL_T dxFile(3)

      REAL_T  twicePi, kxd, kyd, kzd
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

      namelist /fortin/ rt_splitx, rt_xfrontw, rt_den_1, rt_den_2, &
                       rt_pertamp, rt_nfreq, rt_graddenerr

      namelist /fortin/ eul_nfreq, eul_pertamp, iseed

      namelist /fortin/ Vco, Rfu, Rtran, tVco_l, tVco_r, Vco_l, Vco_r

      namelist /fortin/ grav_angle, omega, infl_time_offset, ref_centre, ref_radius, time_offset, &
                       thermal_expansion, heating_coeff, heating_centre, heating_radius

      namelist /fortin/ density_pert, interface_height, wavelength_min, wavelength_max, tracer_height

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


      domnlo(1) = problo(1)
      domnlo(2) = problo(2)
      domnlo(3) = problo(3) 
      domnhi(1) = probhi(1)
      domnhi(2) = probhi(2)
      domnhi(3) = probhi(3)

!c
!c     Initialize the common blocks
!c
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

      if (probtype .eq. 4) then
         call initchannel(level,time,lo,hi,nscal, &
          	          vel,scal,DIMS(state),press,DIMS(press), &
                         dx,xlo,xhi) 

      else if (probtype .eq. 5) then
         call initpervort(level,time,lo,hi,nscal, &
          	          vel,scal,DIMS(state),press,DIMS(press), &
                         dx,xlo,xhi)

      else if (probtype .eq. 6) then
         call inithotspot(level,time,lo,hi,nscal, &
          	          vel,scal,DIMS(state),press,DIMS(press), &
                         dx,xlo,xhi)

      else if (probtype .eq. 7) then
         call initeuler(level,time,lo,hi,nscal, &
                         vel,scal,DIMS(state),press,DIMS(press), &
                         dx,xlo,xhi)

      else if (probtype .eq. 9) then
         call initviscbench(level,time,lo,hi,nscal, &
                           vel,scal,DIMS(state),press,DIMS(press), &
                           dx,xlo,xhi)

      else if (probtype .eq. 10) then
         call initrt(level,time,lo,hi,nscal, &
             vel,scal,DIMS(state),press,DIMS(press), &
             dx,xlo,xhi)

      else if (probtype .eq. 12) then
         call taylorgreen(level,time,lo,hi,nscal, &
             vel,scal,DIMS(state),press,DIMS(press), &
             dx,xlo,xhi)

      else if (probtype .eq. 22) then
         call initdenadvect(level,time,lo,hi,nscal, &
             vel,scal,DIMS(state),press,DIMS(press), &
             dx,xlo,xhi)

      else if (probtype .eq. 29) then
         call initfromrest(level,time,lo,hi,nscal, &
             vel,scal,DIMS(state),press,DIMS(press), &
             dx,xlo,xhi)
      else
         write(6,*) "INITDATA: bad probtype = ",probtype
      end if
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
      end
!c
!c ::: -----------------------------------------------------------
!c ::: Initialise system from rest. Introduced for the lid-driven cavity
!c ::: test case. 
!c
      subroutine initfromrest(level,time,lo,hi,nscal, &
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

#include <probdata.H>

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            do k = lo(3), hi(3)
               scal(i,j,k,1) = one
               vel(i,j,k,1) = zero
               vel(i,j,k,2) = zero
               vel(i,j,k,3) = zero
               do n = 2,nscal-1
                  scal(i,j,k,n) = one
               end do                  
               scal(i,j,k,nscal) = zero
            end do
         end do
      end do

      end subroutine initfromrest
!c
!c ::: -----------------------------------------------------------
!c
      subroutine initrt(level,time,lo,hi,nscal, &
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
      integer i, j, k, n, m
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  pert,ztemp

#include <probdata.H>

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
                  pert = pert &
                      + sin( &
                      2.0D0*rt_PI*dble(n)*x/f_probhi(1) + rt_ranphse(n,m,1) &
                      ) &
                      * sin( &
                      2.0D0*rt_PI*dble(m)*y/f_probhi(2) + rt_ranphse(n,m,2) &
                      ) &
                      *rt_ranampl(n,m)
               end do
            end do
            ztemp = rt_splitx - pert*rt_pertamp
            do k = lo(3), hi(3)
               z = xlo(3) + hz*(dble(k-lo(3)) + 0.5d0)
               scal(i,j,k,1) = &
                   rt_den_1  &
                   + (rt_den_2-rt_den_1)*zblend1(z,ztemp,rt_xfrontw)
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

      end subroutine initrt
!c
!c ::: -----------------------------------------------------------
!c
      subroutine initpervort(level,time,lo,hi,nscal, &
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

               vel(i,j,k,1) = tanh(30.*(.25-abs(y-.5)))
               vel(i,j,k,2) = .05*sin(two*Pi*x)
               vel(i,j,k,3) = zero

               scal(i,j,k,1) = one
               do n = 2,nscal-1
                  scal(i,j,k,n) = one
               end do

               dist = sqrt((x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2)
               scal(i,j,k,nscal) = merge(one,zero,dist.lt.radblob)
            end do
         end do
      end do

      end subroutine initpervort
!c
!c
!c ::: -----------------------------------------------------------
!c ::: This case is an unsteady  viscous benchmark for which the
!c ::: exact solution is,
!c :::     u(x,y,t) = - Cos(Pi x) Sin(Pi y) Exp(-2 Pi^2 Nu t)
!c :::     v(x,y,t) =   Sin(Pi x) Cos(Pi y) Exp(-2 Pi^2 Nu t)
!c :::     p(x,y,t) = - {Cos(2 Pi x) + Cos(2 Pi y)} Exp(-4 Pi^2 Nu t) / 4
!c ::: In the utilities, iamrlib/BenchMarks, there is a
!c ::: tool ViscBench2d.cpp that reads a plot file and compares the
!c ::: solution against this exact solution.  This benchmark was
!c ::: originally derived by G.I. Taylor (Phil. Mag., Vol. 46, No. 274,
!c ::: pp. 671-674, 1923) and Ethier and Steinman
!c ::: (Intl. J. Num. Meth. Fluids, Vol. 19, pp. 369-375, 1994) give
!c ::: the pressure field.
!c
      subroutine initviscbench(level,time,lo,hi,nscal, &
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
      REAL_T  spx, spy, spz, cpx, cpy, cpz

#include <probdata.H>

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

!c
!c           Uniform in the X-direction
!c
            if (vb_unifdir .eq. 0) then
              vel(i,j,k,1) =   zero
              vel(i,j,k,2) =   spz*cpy
              vel(i,j,k,3) = - cpz*spy
              scal(i,j,k,1) =   one
              do n = 2, nscal
                scal(i,j,k,n) =   cpz*cpy
              enddo

!c
!c           Uniform in the Y-direction
!c
            elseif (vb_unifdir .eq. 1) then
              vel(i,j,k,1) = - cpx*spz
              vel(i,j,k,2) =   zero
              vel(i,j,k,3) =   spx*cpz
              scal(i,j,k,1) =   one
              do n = 2, nscal
                scal(i,j,k,n) =   cpx*cpz
              enddo

!c
!c           Uniform in the Z-direction
!c
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
      end subroutine initviscbench
!c
!c ::: -----------------------------------------------------------
!c
      subroutine initchannel(level,time,lo,hi,nscal, &
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

      if ( adv_dir .ne. 1 ) then
         write(6,*) "initchannel requires adv_dir=1, currently adv_dir=",adv_dir
         stop
      end if

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
            do i = lo(1), hi(1)
               vel(i,j,k,1) = adv_vel
               vel(i,j,k,2) = zero
               vel(i,j,k,3) = zero
               scal(i,j,k,1) = denfact

               do n = 2,nscal-1
                  scal(i,j,k,n) = one
               end do                  

               x = xlo(1) + hx*(float(i-lo(1)) + half)
  	       dist = sqrt((x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2)
	       scal(i,j,k,nscal) = merge(one,zero,dist.lt.radblob)
            end do
         end do
      end do

      end subroutine initchannel
!c ::: -----------------------------------------------------------
!c
      subroutine inithotspot(level,time,lo,hi,nscal, &
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


!c     ::::: local variables
      integer i, j, k, n
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  x_vel, y_vel, z_vel
      REAL_T  dist

#include <probdata.H>

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
!c     AJA - Strangely, this one did not
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
               scal(i,j,k,1) = one/denfact + (one - one/denfact) &
                   *half*(one + tanh(40.*(dist - radblob)))
               scal(i,j,k,2) = merge(one,zero,dist.lt.radblob)
               do n = 3,nscal-1
                  scal(i,j,k,n) = one
               end do
               scal(i,j,k,nscal) = one / scal(i,j,k,1)
            end do
         end do
      end do
      
      end subroutine inithotspot
!c
!c ::: -----------------------------------------------------------
!c
      subroutine initeuler(level,time,lo,hi,nscal, &
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
      REAL_T  hx, hy, hz, r_yz

      REAL_T  eps_input, beta_input, rho_input
      REAL_T  delta_input, kappa_input

      parameter (eps_input=0.05d0, rho_input=0.15d0)
      parameter (beta_input=15.d0, delta_input=0.0333d0)
      parameter (kappa_input=500.d0)

#include <probdata.H>
      
      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(DBLE(k-lo(3)) + half) -half
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(DBLE(j-lo(2)) + half) -half
            r_yz = sqrt(y*y+z*z)
            do i = lo(1), hi(1)

               x = xlo(1) + hx*(DBLE(i-lo(1)) + half) -half
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

      end subroutine initeuler
!c
!c ::: -----------------------------------------------------------
!c
      subroutine taylorgreen(level,time,lo,hi,nscal, &
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
      REAL_T  dist, tpi

#include <probdata.H>

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      tpi = 8.d0*atan(1.d0)

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
	    do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half)

               vel(i,j,k,1) =  velfact*sin(tpi * x)*cos(tpi * y)*cos(tpi * z)
               vel(i,j,k,2) = -velfact*cos(tpi * x)*sin(tpi * y)*cos(tpi * z)
               vel(i,j,k,3) = zero

               scal(i,j,k,1) = denfact

               ! This is the theoretical pressure perturbation from p_0
               scal(i,j,k,2) = (denfact*velfact*velfact/16.d0)*(two+cos(two*z))*(cos(two*x)+cos(two*y))
               do n = 2,nscal-1
                  scal(i,j,k,n) = one
               end do
                  
!  	       dist = sqrt((x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2)
! 	       scal(i,j,k,nscal) = merge(one,zero,dist.lt.radblob)
	    end do
         end do
      end do

      end subroutine taylorgreen
!c
!c ::: -----------------------------------------------------------
!c
      subroutine initrayleightaylor(level,time,lo,hi,nscal, &
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
      integer i, j, k
      REAL_T  x, y, z
      REAL_T  hx, hy, hz, hz0

      integer  idum

      integer x_modes, y_modes, amps, l, m
      REAL_T  pert, Lx, Ly, Lz, lambda_x, lambda_y, lambda, random_number, cx, cy
      REAL_T  pert_min, pert_max, amp_rms, frac
      REAL_T  amp(-32:32,-32:32)
      REAL_T  phs(-32:32,-32:32,2)
      REAL_T  eta(state_l1:state_h1,state_l2:state_h2)

#include <probdata.H>

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      Lx = f_probhi(1)-f_problo(1)
      Ly = f_probhi(2)-f_problo(2)
      Lz = f_probhi(3)-f_problo(3)

!c     Use hz0 for base resolution for consistency (only used in defining perturbation)
      hz0 = Lz/512.d0

      idum=-1
      random_number = ran1(idum)
!
!c      write (*,*) "wavelength_min = ",wavelength_min
!c      write (*,*) "wavelength_max = ",wavelength_max
      
      x_modes = int( 0.5 + Lx / wavelength_min )
      y_modes = int( 0.5 + Ly / wavelength_min )
      
      if (x_modes.gt.32) x_modes=32
      if (y_modes.gt.32) y_modes=32
      
      write (*,*) "x_modes = ",x_modes
      write (*,*) "y_modes = ",y_modes
      
      amps = 0
      amp_rms = zero
      
      do l = -x_modes, x_modes
         if (x_modes.eq.0) then
            lambda_x=zero
         else
            lambda_x = Lx / float(l)
         endif
            
         do m = -y_modes, y_modes
            if (y_modes.eq.0) then
               lambda_y=zero
            else
               lambda_y = Ly / float(m)
            endif
            
            lambda = sqrt( lambda_x*lambda_x + lambda_y*lambda_y )
            
            if (lambda .ge. wavelength_min .and. lambda .le. wavelength_max) then
               amp(l,m) = two*(ran1(idum)-half)*lambda*lambda
               phs(l,m,1) = two*Pi*ran1(idum)
               phs(l,m,2) = two*Pi*ran1(idum)
               amps = amps + 1
            else
               amp(l,m) = zero
            end if

         end do
      end do
      
      write (*,*) "amps = ",amps
         
!c     Need to integrate to normalise the perturbation
      amp_rms = zero
      do j = 1,int(0.5+Ly/hy)
         y = hy*(float(j) + half)/Ly
         do i = 1,int(0.5+Lx/hx)
            x = hx*(float(i) + half)/Lx
            pert = zero
            do l = -x_modes, x_modes
               do m = -y_modes, y_modes
                  if (abs(amp(l,m)).gt.zero) then
                     cx = cos(2*Pi*float(l)*x/Lx+phs(l,m,1))
                     cy = cos(2*Pi*float(m)*y/Ly+phs(l,m,2))
                     pert = pert + amp(l,m)*cx*cy
                  endif
               end do
            end do
            amp_rms = amp_rms + pert*pert
         end do
      end do
      
      amp_rms = sqrt(amp_rms*hx*hy/(Lx*Ly))

      write (*,*) "amp_rms = ",amp_rms

      do l = -x_modes, x_modes
         do m = -y_modes, y_modes
            amp(l,m) = density_pert*hz0*amp(l,m)/amp_rms
         enddo
      enddo
      
      pert_min = 1.d10
      pert_max = -pert_min

      do j = lo(2), hi(2)
         y = xlo(2) + hy*(float(j-lo(2)) + half)
         do i = lo(1), hi(1)
            x = xlo(1) + hx*(float(i-lo(1)) + half)
            pert = zero
            do l = -x_modes, x_modes
               do m = -y_modes, y_modes
                  if (abs(amp(l,m)).gt.zero) then
                     cx = cos(2*Pi*float(l)*x/Lx+phs(l,m,1))
                     cy = cos(2*Pi*float(m)*y/Ly+phs(l,m,2))
                     pert = pert + amp(l,m)*cx*cy
                  endif
               end do
            end do
            eta(i,j) = pert
            if (pert.gt.pert_max) pert_max=pert
            if (pert.lt.pert_min) pert_min=pert
         end do
      end do
      
      write (*,*) "pert_min/hz = ",pert_min/hz
      write (*,*) "pert_max/hz = ",pert_max/hz

!c     Go to work...
      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half) - interface_height
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               
!c     Set velocities to zero if required
               vel(i,j,k,1) = zero
               vel(i,j,k,2) = zero
               vel(i,j,k,3) = zero
               
               if (abs(z-eta(i,j)).lt.half*hz) then
                  frac = half+(eta(i,j)-z)/hz
                  scal(i,j,k,1) = den1*(one-frac) + den2*frac
               else
                  if (z.gt.zero) then
                     scal(i,j,k,1) = den1
                  else
                     scal(i,j,k,1) = den2
                  endif
               endif

!c     Set a tracer on the interface
               if (abs(z).lt.half*tracer_height) then
                  scal(i,j,k,2) = one
               else
                  scal(i,j,k,2) = zero
               endif

            enddo
         enddo
      enddo

      end subroutine initrayleightaylor
!c
!c ::: -----------------------------------------------------------
!c
      subroutine initraneul(level,time,lo,hi,nscal, &
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
      integer i, j, k, n, ifx, jfx, kfx
      integer ift,jft,kft
      REAL_T  x, y, z
      REAL_T  sinx, siny, sinz
      REAL_T  cosx, cosy, cosz
      REAL_T  keng
      REAL_T  hx, hy, hz
      REAL_T  piloc

#include <probdata.H>

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)
      piloc = 4.d0*atan2(1.d0,1.d0)


      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
	    do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half)

               vel(i,j,k,1) = zero
               vel(i,j,k,2) = zero
               vel(i,j,k,3) = zero

!c              do kfx = 1, eul_nfreq
!c              do jfx = 1, eul_nfreq
!c              do ifx = 1, eul_nfreq
               do kft = 1, eul_nfreq
               do jft = 1, eul_nfreq
               do ift = 1, eul_nfreq
                  ifx = ift-1
                  jfx = jft-1
                  kfx = kft-1

                  sinx = sin(2.d0*piloc*dfloat(ifx)*x + eul_ranphase(ifx,jfx,kfx,1))
                  siny = sin(2.d0*piloc*dfloat(jfx)*y + eul_ranphase(ifx,jfx,kfx,2))
                  sinz = sin(2.d0*piloc*dfloat(kfx)*z + eul_ranphase(ifx,jfx,kfx,3))
                  cosx = dfloat(ifx)* cos(2.d0*piloc*dfloat(ifx)*x + eul_ranphase(ifx,jfx,kfx,1))
                  cosy = dfloat(jfx)* cos(2.d0*piloc*dfloat(jfx)*y + eul_ranphase(ifx,jfx,kfx,2))
                  cosz = dfloat(kfx)* cos(2.d0*piloc*dfloat(kfx)*z + eul_ranphase(ifx,jfx,kfx,3))

                  vel(i,j,k,1) = vel(i,j,k,1) &
                     + eul_ranampl(ifx,jfx,kfx,3)*sinx*cosy*sinz &
                     - eul_ranampl(ifx,jfx,kfx,2)*sinx*siny*cosz
                  vel(i,j,k,2) = vel(i,j,k,2)  &
                     + eul_ranampl(ifx,jfx,kfx,1)*sinx*siny*cosz &
                     - eul_ranampl(ifx,jfx,kfx,3)*cosx*siny*sinz
                  vel(i,j,k,3) = vel(i,j,k,3)  &
                     + eul_ranampl(ifx,jfx,kfx,2)*cosx*siny*sinz &
                     - eul_ranampl(ifx,jfx,kfx,1)*sinx*cosy*sinz
               
               enddo
               enddo
               enddo

               vel(i,j,k,1) = eul_pertamp*vel(i,j,k,1)
               vel(i,j,k,2) = eul_pertamp*vel(i,j,k,2)
               vel(i,j,k,3) = eul_pertamp*vel(i,j,k,3)
               
               scal(i,j,k,1) = one
               do n = 2,nscal
                     scal(i,j,k,n) = one
               end do
                  
	    end do
         end do
      end do

      end subroutine initraneul
!c
!c ::: -----------------------------------------------------------
!c
      subroutine initdenadvect(level,time,lo,hi,nscal, &
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

               vel(i,j,k,1) = one
               vel(i,j,k,2) = one
               vel(i,j,k,3) = one

               scal(i,j,k,1) = one
               scal(i,j,k,2) = zero

               do n = 2,nscal-1
                  scal(i,j,k,n) = zero
               end do

	    end do
         end do
      end do

!c     initialize a block of tracer
      scal(7,7,7,2) = 1.0d0
      scal(7,7,8,2) = 1.0d0
      scal(7,8,7,2) = 1.0d0
      scal(7,8,8,2) = 1.0d0
      scal(8,7,7,2) = 1.0d0
      scal(8,7,8,2) = 1.0d0
      scal(8,8,7,2) = 1.0d0
      scal(8,8,8,2) = 1.0d0

      end subroutine initdenadvect
!c
!c
!c ::: -----------------------------------------------------------
!c
      subroutine initshearlayer(level,time,lo,hi,nscal, &
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
      REAL_T  x, y, z, t
      REAL_T  hx, hy, hz
      REAL_T  Lx, Ly, Lz
      REAL_T  z0

#include <probdata.H>

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      Lx = f_probhi(1)-f_problo(1)
      Ly = f_probhi(2)-f_problo(2)
      Lz = f_probhi(3)-f_problo(3)

!c     Go to work...
      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
            do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half)
               
               t = time-half*x*(vel1+vel2)
               z0 = interface_height
               do n = 1, 10
                  z0 = z0 + mag(n)*cos(freq(n)*t+phi1(n))*cos(float(n)*y/Ly+phi2(n))
               enddo

               vel(i,j,k,1) = half*((vel1+vel2)+(vel1-vel2)*tanh(two*(z-z0)/delta0))
               vel(i,j,k,2) = zero
               vel(i,j,k,3) = zero

               scal(i,j,k,1) = half*((den1+den2)+(den1-den2)*tanh(two*(z-z0)/delta0))
               scal(i,j,k,2) = exp(-((z-interface_height)/delta0)**2)

            enddo
         enddo
      enddo

      end subroutine initshearlayer

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

!c     probtype = CHANNEL
      if (probtype .eq. 4) then

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

!c     probtype = VISCOUS BENCHMARK
      else if (probtype .eq. 9) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do


      else if (probtype .eq. 10) then

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
              end do
           end do
        end do

      else if (probtype .eq. 12) then
        
      else if (probtype .eq. 22) then
 
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
      print *,'DONT KNOW THIS PROBTYPE IN FORT_ADV2ERROR ',probtype
      stop
      
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

!c     probtype = CHANNEL
      if (probtype .eq. 4) then

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

      else if (probtype .eq. 10 .or. probtype.eq.12) then

      else if (probtype .eq. 22) then

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


!c     probtype = CHANNEL
      if (probtype .eq. 4) then

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

!c     probtype = VISCOUS BENCHMARK
      else if (probtype .eq. 9) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
         end do

      else if (probtype .eq. 10 .or. probtype.eq.12) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr)
               end do
            end do
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
           do k = ARG_L3(rho), domlo(3)-1
              do j = ARG_L2(rho), ARG_H2(rho)
                 do i = ARG_L1(rho), ARG_H1(rho)
                    rho(i,j,k) = denfact
                 end do
              end do
           end do
      end if
      
      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(rho).gt.domhi(3)) then
            do k = domhi(3)+1, ARG_H3(rho)
               do j = ARG_L2(rho), ARG_H2(rho)
                  do i = ARG_L1(rho), ARG_H1(rho)
                     rho(i,j,k) = denfact
                  end do
               end do
            end do
      end if

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
            do k = ARG_L3(adv), domlo(3)-1
               do j = ARG_L2(adv), ARG_H2(adv)
                  do i = ARG_L1(adv), ARG_H1(adv)
                     adv(i,j,k) = zero
                  end do
               end do
            end do
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
            do k = ARG_L3(adv), domlo(3)-1
               do j = ARG_L2(adv), ARG_H2(adv)
                  do i = ARG_L1(adv), ARG_H1(adv)
                     adv(i,j,k) = zero
                  end do
               end do
            end do
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
            do k = ARG_L3(temp), domlo(3)-1
               do j = ARG_L2(temp), ARG_H2(temp)
                  do i = ARG_L1(temp), ARG_H1(temp)
                     temp(i,j,k) = one
                  end do
               end do
            end do
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

      REAL_T x,r,u1,u2,u3,u_inf,eta

      parameter (constn=.22089323)

      lo(1) = ARG_L1(u)
      lo(2) = ARG_L2(u)
      lo(3) = ARG_L3(u)
      hi(1) = ARG_H1(u)
      hi(2) = ARG_H2(u)
      hi(3) = ARG_H3(u)

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
            do k = ARG_L3(u), domlo(3)-1
               do j = ARG_L2(u), ARG_H2(u)
                  do i = ARG_L1(u), ARG_H1(u)
                     u(i,j,k) = zero
                  end do
               end do
            end do
      end if

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(u).gt.domhi(3)) then
         if (probtype .eq. 29) then
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

      REAL_T x,y,z,r,u1,u2,u3,u_inf,eta

      lo(1) = ARG_L1(v)
      lo(2) = ARG_L2(v)
      lo(3) = ARG_L3(v)
      hi(1) = ARG_H1(v)
      hi(2) = ARG_H2(v)
      hi(3) = ARG_H3(v)
      
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
            do k = ARG_L3(v), domlo(3)-1
               do j = ARG_L2(v), ARG_H2(v)
                  do i = ARG_L1(v), ARG_H1(v)
                        v(i,j,k) = zero
                  end do
               end do
            end do
      end if        

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(v).gt.domhi(3)) then
            do k = domhi(3)+1, ARG_H3(v)
               do j = ARG_L2(v), ARG_H2(v)
                  do i = ARG_L1(v), ARG_H1(v)
                     v(i,j,k) = zero
                  end do
               end do
            end do
      end if        

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

      REAL_T x,y,z,r,u1,u2,u3,u_inf,eta
      REAL_T Lx, Ly, pert, factor

      lo(1) = ARG_L1(w)
      lo(2) = ARG_L2(w)
      lo(3) = ARG_L3(w)
      hi(1) = ARG_H1(w)
      hi(2) = ARG_H2(w)
      hi(3) = ARG_H3(w)

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
            do k = ARG_L3(w), domlo(3)-1
               do j = ARG_L2(w), ARG_H2(w)
                  do i = ARG_L1(w), ARG_H1(w)
                     w(i,j,k) = zero
                  end do
               end do
            end do
      end if        

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(w).gt.domhi(3)) then
            do k = domhi(3)+1, ARG_H3(w)
               do j = ARG_L2(w), ARG_H2(w)
                  do i = ARG_L1(w), ARG_H1(w)
                     w(i,j,k) = z_vel
                  end do
               end do
            end do
      end if        

      end subroutine FORT_ZVELFILL
      
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
