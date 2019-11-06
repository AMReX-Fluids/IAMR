
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PROB_NS_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 2

module prob_2D_module

  implicit none

  private

  public :: amrex_probinit, FORT_INITDATA, initbubble, initspin, &
            initviscbench, initvort, initchannel, initpervort, &
            inithotspot, initrt, inittraceradvect, initfromrest, &
            init_taylorgreen, &
            FORT_DENERROR, FORT_AVERAGE_EDGE_STATES, FORT_MAKEFORCE, &
            FORT_ADVERROR, FORT_ADV2ERROR, FORT_TEMPERROR, FORT_MVERROR, &
            FORT_DENFILL, FORT_ADVFILL, FORT_TEMPFILL, FORT_XVELFILL, &
            FORT_YVELFILL, FORT_PRESFILL, FORT_DIVUFILL, FORT_DSDTFILL

contains


!c ::: -----------------------------------------------------------
!c ::: This routine is called at problem initialization time
!c ::: and when restarting from a checkpoint file.
!c ::: The purpose is (1) to specify the initial time value
!c ::: (not all problems start at time=0.0) and (2) to read
!c ::: problem specific data from a namelist or other input
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

#ifdef BL_DO_FLCT
#include <INFL_FORCE_F.H>
#include <infl_frc.H>
#endif

!c
!c Dimensions of the Inflow file.
!c
      INTEGER dimFile(3)
      integer nCompFile
      parameter (nCompFile = 2)
      REAL_T dxFile(3)

      namelist /fortin/ denerr, vorterr, adverr, temperr, &
     			denfact, xblob, yblob, zblob, radblob, &
                       velfact, probtype, randfact, bubgrad, &
     			rhozero, tempzero, c_d, r_d, grav_angle, &
                       adv_dir, adv_vel, axis_dir, radvort, &
                       lid_vel 
#ifdef BL_DO_FLCT
                       ,forceInflow, numInflPlanesStore, strmwse_dir, &
                       forceLo, forceHi, flct_file, turb_scale
#endif

!c
!c      Build "probin" filename -- the name of file containing fortin namelist.
!c
      integer maxlen, isioproc
      parameter (maxlen=256)

      character probin*(maxlen)

      call bl_pd_is_ioproc(isioproc)

      if (namlen .gt. maxlen) call bl_error('probin file name too long')

      do i = 1, namlen
         probin(i:i) = char(name(i))
      end do

      untin = 9
      if (namlen .eq. 0) then
         open(untin,file='probin',form='formatted',status='old')
      else
         open(untin,file=probin(1:namlen),form='formatted',status='old')
      end if

#ifdef BL_DO_FLCT
      forceInflow = .FALSE.
      numInflPlanesStore = 16
      forceLo = .TRUE.
      forceHi = .FALSE.
      strmwse_dir = FLCT_YVEL
      flct_file = ""
      turb_scale = 1
#endif

      read(untin,fortin)
      if (isioproc .eq. 1) write(6,fortin)
      close(unit=untin)

#ifdef BL_DO_FLCT
      if (forceInflow .eqv. .FALSE.) then
         forceLo = .FALSE.
         forceHi = .FALSE.
      else
         if (flct_file .ne. "") then
            if (isioproc .eq. 1) print*, 'Initializing turbulence ...'
            open(20, file=flct_file, form='unformatted')
            call RD_SCL_FLCTHD(20,nCompFile,dimFile,probSizeFile,dxFile)
            close(20)
         endif
         if (strmwse_dir .NE. 2) then
            call bl_error('turbulent inflow needs strmwse_dir=2')
         end if
         if (isioproc .eq. 1) then
            print *, 'dimFile: ',      (dimFile(i),i=1,3)
            print *, 'probSizeFile: ', (probSizeFile(i),i=1,3)
            print *, 'dxFile: ',       (dxFile(i),i=1,3)
         end if
      endif
#endif

      do i=1, SDIM
        f_problo(i) = problo(i)
        f_probhi(i) = probhi(i)
      enddo

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
!c ::: dx       => cell size
!c ::: xlo,xhi   => physical locations of lower left and upper
!c :::              right hand corner of grid.  (does not include
!c :::		   ghost region).
!c ::: -----------------------------------------------------------
      subroutine FORT_INITDATA(level,time,lo,hi,nscal, &
     	 	               vel,scal,DIMS(state),press,DIMS(press), &
                              dx,xlo,xhi) &
                              bind(C, name="FORT_INITDATA")
      integer    level, nscal
      integer    lo(SDIM),hi(SDIM)
      integer    DIMDEC(state)
      integer    DIMDEC(press)
      REAL_T     time, dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     vel(DIMV(state),SDIM)
      REAL_T    scal(DIMV(state),nscal)
      REAL_T   press(DIMV(press))

#include <probdata.H>

      if (probtype .eq. 1) then
         call initspin(level,time,lo,hi,nscal, &
          	       vel,scal,DIMS(state),press,DIMS(press), &
                      dx,xlo,xhi)

      else if (probtype .eq. 2) then
         call initbubble(level,time,lo,hi,nscal, &
          	         vel,scal,DIMS(state),press,DIMS(press), &
                        dx,xlo,xhi)

      else if (probtype .eq. 3) then
         call initvort(level,time,lo,hi,nscal, &
          	       vel,scal,DIMS(state),press,DIMS(press), &
                      dx,xlo,xhi)

      else if (probtype .eq. 4) then
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
         call initviscbench(level,time,lo,hi,nscal, &
          	            vel,scal,DIMS(state),press,DIMS(press),  &
                           dx,xlo,xhi)

      else if (probtype .eq. 8) then
         call initrt(level,time,lo,hi,nscal, &
                vel,scal,DIMS(state),press,DIMS(press), &
                    dx,xlo,xhi)

      else if (probtype .eq. 9) then
         call inittraceradvect(lo,hi,nscal, &
                          vel,scal,DIMS(state), &
                              dx,xlo,xhi)

      else if (probtype .eq. 10 .or. probtype .eq. 11 .or. &
              probtype .eq. 12) then
         call initfromrest(lo,hi,nscal, &
                       vel,scal,DIMS(state), &
                           dx,xlo,xhi)

      else if (probtype .eq. 13) then
         call init_taylorgreen(level,time,lo,hi,nscal, &
                vel,scal,DIMS(state),press,DIMS(press), &
                    dx,xlo,xhi)
      
      else
         write(6,*) "INITDATA: bad probtype = ",probtype
         stop
      end if

      end subroutine FORT_INITDATA
!c
!c ::: -----------------------------------------------------------
!c
      subroutine initbubble(level,time,lo,hi,nscal, &
     	 	            vel,scal,DIMS(state),press,DIMS(press), &
                           dx,xlo,xhi) &
                           bind(C, name="initbubble")

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
      integer i, j, n
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  dist
      REAL_T  x_vel, y_vel

#include <probdata.H>

      hx = dx(1)
      hy = dx(2)

      if (adv_dir .eq. 1) then
         x_vel = adv_vel
         y_vel = zero
      else if (adv_dir .eq. 2) then
         x_vel = zero
         y_vel = adv_vel
      else 
         write(6,*) "initbubble: adv_dir = ",adv_dir
         stop
      end if

      do j = lo(2), hi(2)
         y = xlo(2) + hy*(float(j-lo(2)) + half)
         do i = lo(1), hi(1)
            x = xlo(1) + hx*(float(i-lo(1)) + half)
             dist = sqrt((x-xblob)**2 + (y-yblob)**2)
            vel(i,j,1) = x_vel
            vel(i,j,2) = y_vel
            scal(i,j,1) = one + half*(denfact-one)*(one-tanh(30.*(dist-radblob)))
!c           scal(i,j,1) = merge(denfact,one,dist.lt.radblob)
            do n = 2,nscal-1
               scal(i,j,n) = one
            end do                  
            scal(i,j,nscal) = merge(one,zero,dist.lt.radblob)
         end do
      end do

      end subroutine initbubble
!c
!c ::: -----------------------------------------------------------
!c
      subroutine initspin(level,time,lo,hi,nscal, &
     	 	          vel,scal,DIMS(state),press,DIMS(press), &
                         dx,xlo,xhi) &
                         bind(C, name="initspin")

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
      integer i, j, n
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  dist
      REAL_T  x_vel, y_vel
      REAL_T  spx, spy, cpx, cpy

#include <probdata.H>

      hx = dx(1)
      hy = dx(2)

      if (adv_dir .eq. 1) then
         x_vel = adv_vel
         y_vel = zero
      else if (adv_dir .eq. 2) then
         x_vel = zero
         y_vel = adv_vel
      else 
         write(6,*) "INITSPIN: adv_dir = ",adv_dir
         stop
      end if

         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
            spy = sin(Pi*y)
            cpy = cos(Pi*y)
            do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half)

               spx = sin(Pi*x)
               cpx = cos(Pi*x)

               vel(i,j,1) = x_vel - velfact*two*spy*cpy*spx**2
               vel(i,j,2) = y_vel + velfact*two*spx*cpx*spy**2

               dist = sqrt((x-xblob)**2 + (y-yblob)**2)

               scal(i,j,1) = one + (denfact-one) * tanh(10.*(dist-radblob))
               do n = 2,nscal-1
                  scal(i,j,n) = one
               end do                  
               scal(i,j,nscal) = merge(one,zero,dist.lt.radblob)

            end do
         end do

      end subroutine initspin
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
                              dx,xlo,xhi) &
                              bind(C, name="initviscbench")

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
      integer i, j, n
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  spx, spy, cpx, cpy

#include <probdata.H>

      hx = dx(1)
      hy = dx(2)

      do j = lo(2), hi(2)
         y = xlo(2) + hy*(float(j-lo(2)) + half)
         spy = sin(Pi*y)
         cpy = cos(Pi*y)

         do i = lo(1), hi(1)
            x = xlo(1) + hx*(float(i-lo(1)) + half)

            spx = sin(Pi*x)
            cpx = cos(Pi*x)

            vel(i,j,1) = - cpx*spy
            vel(i,j,2) =   spx*cpy

            scal(i,j,1) = one
            do n = 2,nscal
               scal(i,j,n) = cpx*cpy
            end do                  

         end do
      end do

      end subroutine initviscbench
!c
!c ::: -----------------------------------------------------------
!c
      subroutine initvort(level,time,lo,hi,nscal, &
     	 	          vel,scal,DIMS(state),press,DIMS(press), &
                         dx,xlo,xhi) &
                         bind(C, name="initvort")
                         
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
      integer i, j, n
      REAL_T  x, y, r
      REAL_T  hx, hy
      REAL_T  c, ux, uy
      REAL_T  umagin, umagout, absu, sinth, costh
      REAL_T  small, a, b, r0

#include <probdata.H>

      hx = dx(1)
      hy = dx(2)
      small = 1.0e-10

      r0 = two/three * radvort
      a = one / ((radvort - r0)*(two*radvort - r0))
      b = a * radvort**2 * (radvort - r0)

         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half) - yblob
            do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half) - xblob
               r = sqrt(x**2 + y**2)
!c              umagin = .5*r - 4*r**3
!c              umagout = radvort*(.5*radvort - 4*radvort**3)/max(radvort,r)
               umagin = velfact * (one - a*(r - r0)**2)
               umagout = velfact * b/max(radvort,r)
               absu = merge(umagout,umagin,(r - radvort) .ge. 0.0d0)
               sinth = y/max(r,small*radvort)
               costh = x/max(r,small*radvort)
               vel(i,j,1) = -absu*sinth
               vel(i,j,2) = absu*costh
               scal(i,j,1) = merge(denfact,one,r.lt.radblob)
               do n = 2,nscal-1
                  scal(i,j,n) = one
               end do                  
               scal(i,j,nscal) = merge(one,zero,r.lt.radblob)
            end do
         end do

      end subroutine initvort
!c
!c ::: -----------------------------------------------------------
!c
      subroutine initchannel(level,time,lo,hi,nscal, &
     	 	             vel,scal,DIMS(state),press,DIMS(press), &
                            dx,xlo,xhi)&
                            bind(C, name="initchannel")

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
      integer i, j, n
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  dist

#include <probdata.H>

      hx = dx(1)
      hy = dx(2)

      if ( adv_dir .ne. 1 ) then
         write(6,*) "initchannel requires adv_dir=1, currently adv_dir=",adv_dir
         stop
      end if

      do j = lo(2), hi(2)
         y = xlo(2) + hy*(float(j-lo(2)) + half)
         do i = lo(1), hi(1)
            vel(i,j,1) = adv_vel
            vel(i,j,2) = zero
            scal(i,j,1) = denfact

            do n = 2,nscal-1
               scal(i,j,n) = one
            end do                  

            x = xlo(1) + hx*(float(i-lo(1)) + half)
            dist = sqrt((x-xblob)**2 + (y-yblob)**2)
            scal(i,j,nscal) = merge(one,zero,dist.lt.radblob)

         end do
      end do

      end subroutine initchannel
!c
!c ::: -----------------------------------------------------------
!c
      subroutine initpervort(level,time,lo,hi,nscal, &
     	 	             vel,scal,DIMS(state),press,DIMS(press), &
                            dx,xlo,xhi) &
                            bind(C, name="initpervort")

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
      integer i, j, n
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  dist

#include <probdata.H>

      hx = dx(1)
      hy = dx(2)

         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
            do i = lo(1), hi(1)

               x = xlo(1) + hx*(float(i-lo(1)) + half)

               vel(i,j,1) = tanh(30.*(.25-abs(y-.5)))
               vel(i,j,2) = .05*sin(two*Pi*x)

               scal(i,j,1) = one
               do n = 2,nscal-1
                  scal(i,j,n) = one
               end do
                  
               dist = sqrt((x-xblob)**2 + (y-yblob)**2)
               scal(i,j,nscal) = merge(one,zero,dist.lt.radblob)
            end do
         end do

      end subroutine initpervort
!c
!c ::: -----------------------------------------------------------
!c
      subroutine inithotspot(level,time,lo,hi,nscal, &
     	 	             vel,scal,DIMS(state),press,DIMS(press), &
                            dx,xlo,xhi) &
                            bind(C, name="inithotspot")

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
      integer i, j, n
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  x_vel, y_vel
      REAL_T  dist

#include <probdata.H>

      hx = dx(1)
      hy = dx(2)

      if (adv_dir .eq. 1) then
         x_vel = adv_vel
         y_vel = zero
      else if (adv_dir .eq. 2) then
         x_vel = zero
         y_vel = adv_vel
      else 
         write(6,*) "inithotspot: adv_dir = ",adv_dir
         stop
      end if

      do j = lo(2), hi(2)
         y = xlo(2) + hy*(float(j-lo(2)) + half)
         do i = lo(1), hi(1)
            x = xlo(1) + hx*(float(i-lo(1)) + half)
            dist = sqrt((x-xblob)**2 + (y-yblob)**2)
            vel(i,j,1) = x_vel
            vel(i,j,2) = y_vel
            scal(i,j,1) = one/denfact + (one - one/denfact) &
                *half*(one + tanh(40.*(dist - radblob)))
            scal(i,j,2) = merge(one,zero,dist.lt.radblob)
            do n = 3,nscal-1
               scal(i,j,n) = one
            end do
            scal(i,j,nscal) = one / scal(i,j,1)
         end do
      end do
      
      end subroutine inithotspot
      
!c
!c ::: -----------------------------------------------------------
!c
      subroutine initrt(level,time,lo,hi,nscal, &
     	 	         vel,scal,DIMS(state),press,DIMS(press),  &
                        dx,xlo,xhi) &
                        bind(C, name="initrt")

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
      integer i, j, n
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  x_vel, y_vel
      REAL_T  dist, pertheight, L_x, rho_1, rho_2

#include <probdata.H>

      hx = dx(1)
      hy = dx(2)

      L_x = 0.5d0

      rho_1 = 1.d0
      rho_2 = 2.d0

      do j = lo(2), hi(2)
         y = hy*(float(j) + half)
         do i = lo(1), hi(1)
            x = hx*(float(i) + half)

            pertheight = 0.5d0 + 0.005d0* &
                (cos(2.d0*Pi*x/L_x)+cos(2.d0*Pi*(L_x-x)/L_x))

            scal(i,j,1) = rho_1 +  &
                ((rho_2-rho_1)/2.d0)*(1+tanh((y-pertheight)/0.005d0))

            vel(i,j,1) = zero
            vel(i,j,2) = zero
            scal(i,j,2) = zero

         end do
      end do
      
      end subroutine initrt

!c
!c ::: -----------------------------------------------------------
!c
      subroutine inittraceradvect(lo,hi,nscal, &
                             vel,scal,DIMS(state), &
                                 dx,xlo,xhi) &
                                 bind(C, name="inittraceradvect")

      integer    nscal
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(state)
      REAL_T     dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     vel(DIMV(state),SDIM)
      REAL_T    scal(DIMV(state),nscal)


!c     ::::: local variables
      integer i, j
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  dist

#include <probdata.H>

      hx = dx(1)
      hy = dx(2)

      do j=lo(2),hi(2)
         y = xlo(2) + hy*(float(j-lo(2)) + half)
         do i=lo(1),hi(1)
            x = xlo(1) + hx*(float(i-lo(1)) + half)

            dist = sqrt((0.5d0-x)**2 + (0.5d0-y)**2)

            vel(i,j,1) = 1.0d0
            vel(i,j,2) = 2.0d0
            scal(i,j,1) = 1.0d0
            scal(i,j,2) = 1.0d0*exp(-(10.0d0*dist)**2)

         end do
      end do
      
      end subroutine inittraceradvect
      
!c
!c ::: -----------------------------------------------------------
!c ::: Initialise system from rest. Introduced for the lid-driven cavity
!c ::: test case, also used for Poiseuille flow in square duct. 
!c
      subroutine initfromrest(lo,hi,nscal, &
                         vel,scal,DIMS(state), &
                             dx,xlo,xhi) &
                             bind(C, name="initfromrest")

      integer    nscal
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(state)
      REAL_T     dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     vel(DIMV(state),SDIM)
      REAL_T    scal(DIMV(state),nscal)


!c     ::::: local variables
      integer i, j

#include <probdata.H>

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)

            vel(i,j,1) = 0.0d0
            vel(i,j,2) = 0.0d0
            scal(i,j,1) = 1.0d0
            scal(i,j,2) = 0.0d0

         end do
      end do
      
      end subroutine initfromrest

!c
!c ::: -----------------------------------------------------------
!c
      subroutine init_taylorgreen(level,time,lo,hi,nscal, &
                                  vel,scal,DIMS(state),press,DIMS(press), &
                                  dx,xlo,xhi) &
                                  bind(C, name="init_taylorgreen")
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
      integer i, j, n
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  dist, tpi

#include <probdata.H>

      hx = dx(1)
      hy = dx(2)

      tpi = 8.d0*atan(1.d0)

      do j = lo(2), hi(2)
        y = xlo(2) + hy*(float(j-lo(2)) + half)
        do i = lo(1), hi(1)
          x = xlo(1) + hx*(float(i-lo(1)) + half)

          vel(i,j,1) =  velfact*sin(x)*cos(y)
          vel(i,j,2) = -velfact*cos(x)*sin(y)

          scal(i,j,1) = denfact

          ! This is the theoretical pressure perturbation from p_0
          scal(i,j,2) = (denfact*velfact*velfact/4.d0)*(cos(two*x)+cos(two*y))
          do n = 2,nscal-1
            scal(i,j,n) = one
           end do

!              dist = sqrt((x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2)
!              scal(i,j,k,nscal) = merge(one,zero,dist.lt.radblob)
        end do
      end do

      end subroutine init_taylorgreen

      
!c ::: -----------------------------------------------------------
!c ::: This routine will tag high error cells based on the 
!c ::: magnitude of the density
!c ::: 
!c ::: INPUTS/OUTPUTS:
!c ::: 
!c ::: tag      <=  integer tag array
!c ::: DIMS(tag) => index extent of tag array
!c ::: set       => integer value to tag cell for refinement
!c ::: clear     => integer value to untag cell
!c ::: rho       => density array
!c ::: DIMS(rho) => index extent of rho array
!c ::: lo,hi     => index extent of grid
!c ::: nvar      => number of components in rho array (should be 1)
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of tag array
!c ::: problo    => phys loc of lower left corner of prob domain
!c ::: time      => problem evolution time
!c ::: -----------------------------------------------------------
      subroutine FORT_DENERROR (tag,DIMS(tag),set,clear, &
                               rho,DIMS(rho),lo,hi,nvar, &
                               domlo,domhi,dx,xlo, &
     			        problo,time,level) &
                  bind(C, name="FORT_DENERROR")

      integer   DIMDEC(rho)
      integer   DIMDEC(tag)
      integer   lo(SDIM), hi(SDIM)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    rho(DIMV(rho), nvar)

      integer   i, j

#include <probdata.H>

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            tag(i,j) = merge(set,tag(i,j),rho(i,j,1).lt.denerr)
         end do
      end do

      end subroutine FORT_DENERROR

!c
!c
!c ::: -----------------------------------------------------------
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

!c
!c
!c
!c ::: -----------------------------------------------------------
!c
!c     This routine add the forcing terms to the momentum equation
!c
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

! Local
      REAL_T  :: velmin(0:SDIM-1)
      REAL_T  :: velmax(0:SDIM-1)
      REAL_T  :: scalmin(0:nscal-1)
      REAL_T  :: scalmax(0:nscal-1)
      REAL_T  :: forcemin(scomp:scomp+ncomp-1)
      REAL_T  :: forcemax(scomp:scomp+ncomp-1)
      REAL_T  :: x, y, z
      REAL_T  :: hx, hy, hz
      REAL_T  :: sga, cga
      integer :: kx, ky, mode_count, xstep, ystep
      integer :: isioproc
      integer :: nXvel, nYvel, nZvel, nRho, nTrac, nTrac2
      integer :: nRhoScal, nTracScal, nTrac2Scal
      integer :: i, j, k, n

      hx = dx(1)
      hy = dx(2)
#if ( AMREX_SPACEDIM == 3 )
      hz = dx(3)
#endif

!     Assumes components are in the following order
      nXvel = 0
      nYvel = 1
      nZvel = SDIM-1
      nRho  = SDIM
      nTrac = SDIM+1
      nTrac2= SDIM+2

      nRhoScal   = nRho-SDIM
      nTracScal  = nTrac-SDIM
      nTrac2Scal = nTrac2-SDIM

      if ( getForceVerbose > 0 ) then
         call bl_pd_is_ioproc(isioproc)
         if ( isioproc == 1 ) then

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
            
            ! Get min/max
            do k = f_lo(3), f_hi(3)
               do j = f_lo(2), f_hi(2)
                  do i = f_lo(1), f_hi(1)

                     ! Velocities
                     do n = 0, SDIM-1
                        if (vel(i,j,k,n)>velmax(n)) then
                           velmax(n)=vel(i,j,k,n)
                        endif
                        if (vel(i,j,k,n)<velmin(n)) then
                           velmin(n)=vel(i,j,k,n)
                        endif
                     enddo
                     ! Scalars
                     do n = 0, nscal-1
                        if (scal(i,j,k,n)>scalmax(n)) then
                           scalmax(n)=scal(i,j,k,n)
                        endif
                        if (scal(i,j,k,n)<scalmin(n)) then
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

!     
!     Here's where the forcing actually gets done
!
      
      if ( scomp == 0 ) then
!
!     Do velocity forcing
!
         if ( probtype == 99 .and. abs(grav_angle)>0.001) then
!     Angled gravity
            sga =  gravity * sin(Pi*grav_angle/180.)
            cga = -gravity * cos(Pi*grav_angle/180.)
            do k = f_lo(3), f_hi(3)
               do j = f_lo(2), f_hi(2)
                  do i = f_lo(1), f_hi(1)
                     force(i,j,k,nXvel) = scal(i,j,k,nRhoScal)*sga
#if ( AMREX_SPACEDIM == 2 )
                     force(i,j,k,nYvel) = scal(i,j,k,nRhoScal)*cga
#elif ( AMREX_SPACEDIM == 3 )
                     force(i,j,k,nYvel) = zero
                     force(i,j,k,nZvel) = scal(i,j,k,nRhoScal)*cga
#endif
                  enddo
               enddo
            enddo
!     Default to gravity...
         else if ( abs(gravity) > 0.0001 ) then
            do k = f_lo(3), f_hi(3)
               do j = f_lo(2), f_hi(2)
                  do i = f_lo(1), f_hi(1)
                     force(i,j,k,nXvel) = zero
#if ( AMREX_SPACEDIM == 2 )
                     force(i,j,k,nYvel) = gravity*scal(i,j,k,nRhoScal)
#elif ( AMREX_SPACEDIM == 3 )
                     force(i,j,k,nYvel) = zero
                     force(i,j,k,nZvel) = gravity*scal(i,j,k,nRhoScal)
#endif
                  enddo
               enddo
            enddo
!     else to zero
         else
            do k = f_lo(3), f_hi(3)
               do j = f_lo(2), f_hi(2)
                  do i = f_lo(1), f_hi(1)
                     force(i,j,k,nXvel) = zero
                     force(i,j,k,nYvel) = zero
#if ( AMREX_SPACEDIM == 3 )
                     force(i,j,k,nZvel) = zero
#endif
                  enddo
               enddo
            enddo
         endif
!     End of velocity forcing
      endif

      if ( (scomp+ncomp) > AMREX_SPACEDIM) then
         ! Scalar forcing
         do n = max(scomp,nRho), scomp+ncomp-1
            if ( n == nRho) then
               ! Density
               do k = f_lo(3), f_hi(3)
                  do j = f_lo(2), f_hi(2)
                     do i = f_lo(1), f_hi(1)
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            else if ( n == nTrac ) then
               ! Tracer
               do k = f_lo(3), f_hi(3)
                  do j = f_lo(2), f_hi(2)
                     do i = f_lo(1), f_hi(1)
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            else
               ! Other scalar
               do k = f_lo(3), f_hi(3)
                  do j = f_lo(2), f_hi(2)
                     do i = f_lo(1), f_hi(1)
                        force(i,j,k,n) = zero
                     enddo
                  enddo
               enddo
            endif
         enddo
      endif

      if ( getForceVerbose>0 .and. isioproc==1) then
         do n = scomp,scomp+ncomp-1
            forcemin(n) = 1.d234
            forcemax(n) = -1.d234
         enddo
         do k = f_lo(3), f_hi(3)
            do j = f_lo(2), f_hi(2)
               do i = f_lo(1), f_hi(1)
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
!c ::: DIMS(adv) => index extent of scalar array
!c ::: lo,hi     => index extent of grid
!c ::: nvar      => number of components in rho array (should be 1)
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of tag array
!c ::: problo    => phys loc of lower left corner of prob domain
!c ::: time      => problem evolution time
!c ::: -----------------------------------------------------------
      subroutine FORT_ADVERROR (tag,DIMS(tag),set,clear, &
                               adv,DIMS(adv),lo,hi,nvar, &
                               domlo,domhi,dx,xlo,  &
     			        problo,time,level)&
                  bind(C, name="FORT_ADVERROR")

      integer   DIMDEC(tag)
      integer   DIMDEC(adv)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      integer   lo(SDIM), hi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    adv(DIMV(adv),nvar)

      REAL_T    x, y, ax, ay, aerr, dy
      integer   i, j

#include <probdata.H>

!c     probtype = SPIN
      if (probtype .eq. 1) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do

!c     probtype = BUBBLE
      else if (probtype .eq. 2) then

        if (level .eq. 0) then
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                tag(i,j) = merge(set,tag(i,j),adv(i,j,1).gt.adverr)
             end do
          end do
        end if

!c     probtype = VORTEX IN A BOX
      else if (probtype .eq. 3) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do

!c     probtype = CHANNEL
      else if (probtype .eq. 4) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do

!c     probtype = PERIODIC SHEAR LAYER
      else if (probtype .eq. 5) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do

!c     probtype = HOT SPOT
      else if (probtype .eq. 6) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do


!c     probtype = VISCOUS BENCHMARK
      else if (probtype .eq. 7) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do

!c     probtype = RT
      else if (probtype .eq. 8) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do

      else
        print *,'DONT KNOW THIS PROBTYPE IN FORT_ADVERROR ',probtype
        stop
      end if
 
      end subroutine FORT_ADVERROR
      
      
      subroutine FORT_ADV2ERROR (tag,DIMS(tag),set,clear, &
                               adv,DIMS(adv),lo,hi,nvar, &
                               domlo,domhi,dx,xlo,  &
     			        problo,time,level) &
                  bind(C, name="FORT_ADV2ERROR")

      integer   DIMDEC(tag)
      integer   DIMDEC(adv)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      integer   lo(SDIM), hi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    adv(DIMV(adv),nvar)

      REAL_T    x, y, ax, ay, aerr, dy
      integer   i, j

#include <probdata.H>

!c     probtype = SPIN
      if (probtype .eq. 1) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do

!c     probtype = BUBBLE
      else if (probtype .eq. 2) then

        if (level .eq. 0) then
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                tag(i,j) = merge(set,tag(i,j),adv(i,j,1).gt.adverr)
             end do
          end do
        end if

!c     probtype = VORTEX IN A BOX
      else if (probtype .eq. 3) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do

!c     probtype = CHANNEL
      else if (probtype .eq. 4) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do

!c     probtype = PERIODIC SHEAR LAYER
      else if (probtype .eq. 5) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do

!c     probtype = HOT SPOT
      else if (probtype .eq. 6) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do


!c     probtype = VISCOUS BENCHMARK
      else if (probtype .eq. 7) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do

!c     probtype = RT
      else if (probtype .eq. 8) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do

      else
        print *,'DONT KNOW THIS PROBTYPE IN FORT_ADVERROR ',probtype
        stop
      end if
 
      end subroutine FORT_ADV2ERROR 

!c ::: -----------------------------------------------------------
!c ::: This routine will tag high error cells based on the
!c ::: magnitude or gradient of the temperature
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
                               problo,time,level) &
                               bind(C, name="FORT_TEMPERROR")

      integer   DIMDEC(tag)
      integer   DIMDEC(temp)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      integer   lo(SDIM), hi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    temperature(DIMV(temp),nvar)

      REAL_T    x, y, ax, ay, aerr
      integer   i, j

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
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                tag(i,j) = merge(set,tag(i,j),temperature(i,j,1).gt.temperr)
             end do
          end do
        else
!c         ::::: refine where there is temperature gradient
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = abs(temperature(i+1,j,1) - temperature(i-1,j,1))
                ay = abs(temperature(i,j+1,1) - temperature(i,j-1,1))
                aerr = max(ax,ay)
                tag(i,j) = merge(set,tag(i,j),aerr.gt.bubgrad)
             end do
          end do
        end if


!c     probtype = VISCOUS BENCHMARK
      else if (probtype .eq. 7) then

!c     probtype = VISCOUS BENCHMARK
      else if (probtype .eq. 8) then

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
!c ::: vort      => array of vorticity values
!c ::: DIMS(vor) => index extent of vort array
!c ::: nvar      => number of components in vort array (should be 1)
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

      integer   DIMDEC(tag)
      integer   DIMDEC(vort)
      integer   nvar, set, clear, level
      integer   lo(SDIM), hi(SDIM)
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    vort(DIMV(vort),nvar)

      REAL_T    x, y
      integer   i, j

#include <probdata.H>

!c     probtype = SPIN
      if (probtype .eq. 1) then

!c     probtype = BUBBLE
      else if (probtype .eq. 2) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),abs(vort(i,j,1)).gt.vorterr*2.d0**level)
           end do
        end do

!c     probtype = VORTEX IN A BOX
      else if (probtype .eq. 3) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),abs(vort(i,j,1)).gt.vorterr)
           end do
        end do

!c     probtype = CHANNEL
      else if (probtype .eq. 4) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),abs(vort(i,j,1)).gt.vorterr)
           end do
        end do

!c     probtype = PERIODIC SHEAR LAYER
      else if (probtype .eq. 5) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),abs(vort(i,j,1)).gt.vorterr)
           end do
        end do

!c     probtype = HOT SPOT
      else if (probtype .eq. 6) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),abs(vort(i,j,1)).gt.vorterr)
           end do
        end do

!c     probtype = VISCOUS BENCHMARK
      else if (probtype .eq. 7) then
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),abs(vort(i,j,1)).gt.vorterr)
           end do
        end do

!c     probtype = RT
      else if (probtype .eq. 8) then
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = merge(set,tag(i,j),abs(vort(i,j,1)).gt.vorterr)
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
!c :::        with valid data and that all non-interior cells have
!c ::         have been filled with a large real number.
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

      integer    DIMDEC(rho)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     rho(DIMV(rho))
      integer    bc(SDIM,2)

      integer    i, j

#include <probdata.H>

      call filcc(rho,DIMS(rho),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(rho).lt.domlo(1)) then
         do i = ARG_L1(rho), domlo(1)-1
            do j = ARG_L2(rho), ARG_H2(rho)
               rho(i,j) = denfact
            end do
         end do
      end if            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(rho).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(rho)
            do j = ARG_L2(rho), ARG_H2(rho)
               rho(i,j) = denfact
            end do
         end do
      end if            


      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(rho).lt.domlo(2)) then
         do j = ARG_L2(rho), domlo(2)-1
            do i = ARG_L1(rho), ARG_H1(rho)
               rho(i,j) = denfact
            end do
         end do
      end if            

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(rho).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(rho)
            do i = ARG_L1(rho), ARG_H1(rho)
               rho(i,j) = denfact
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
!c :::        with valid data and that all non-interior cells have
!c ::         have been filled with a large real number.
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

      subroutine FORT_ADVFILL (adv,DIMS(adv),domlo,domhi,dx,xlo,time,bc)&
                               bind(C, name="FORT_ADVFILL")

      integer    DIMDEC(adv)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     adv(DIMV(adv))
      integer    bc(SDIM,2)

      integer    i, j

#include <probdata.H>

      call filcc(adv,DIMS(adv),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(adv).lt.domlo(1)) then
         do i = ARG_L1(adv), domlo(1)-1
            do j = ARG_L2(adv), ARG_H2(adv)
               adv(i,j) = zero
            end do
         end do
      end if            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(adv).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(adv)
            do j = ARG_L2(adv), ARG_H2(adv)
               adv(i,j) = zero
            end do
         end do
      end if            

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(adv).lt.domlo(2)) then

         if (probtype .eq. 12) then

            do j = ARG_L2(adv), domlo(2)-1
               do i = ARG_L1(adv), ARG_H1(adv)
                  adv(i,j) = 1.d0
               end do
            end do

         else

            do j = ARG_L2(adv), domlo(2)-1
               do i = ARG_L1(adv), ARG_H1(adv)
                  adv(i,j) = zero
               end do
            end do

         end if

      end if            

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(adv).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(adv)
            do i = ARG_L1(adv), ARG_H1(adv)
               adv(i,j) = zero
            end do
         end do
      end if            

      end subroutine FORT_ADVFILL

      subroutine FORT_ADV2FILL (adv,DIMS(adv),domlo,domhi,dx,xlo,time,bc)&
                                bind(C, name="FORT_ADV2FILL")

      integer    DIMDEC(adv)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     adv(DIMV(adv))
      integer    bc(SDIM,2)

      integer    i, j

#include <probdata.H>

      call filcc(adv,DIMS(adv),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(adv).lt.domlo(1)) then
         do i = ARG_L1(adv), domlo(1)-1
            do j = ARG_L2(adv), ARG_H2(adv)
               adv(i,j) = zero
            end do
         end do
      end if            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(adv).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(adv)
            do j = ARG_L2(adv), ARG_H2(adv)
               adv(i,j) = zero
            end do
         end do
      end if            

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(adv).lt.domlo(2)) then

         do j = ARG_L2(adv), domlo(2)-1
            do i = ARG_L1(adv), ARG_H1(adv)
               adv(i,j) = zero
            end do
         end do

      end if            

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(adv).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(adv)
            do i = ARG_L1(adv), ARG_H1(adv)
               adv(i,j) = zero
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
!c :::        with valid data and that all non-interior cells have
!c ::         have been filled with a large real number.
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: temperature <=  temperature array
!c ::: DIMS(temp)   => index extent of adv array
!c ::: domlo,hi     => index extent of problem domain
!c ::: dx           => cell spacing
!c ::: xlo          => physical location of lower left hand
!c :::                 corner of temperature array
!c ::: time         => problem evolution time
!c ::: bc           => array of boundary flags bc(BL_SPACEDIM,lo:hi)
!c ::: -----------------------------------------------------------

      subroutine FORT_TEMPFILL (temperature,DIMS(temp),domlo,domhi,dx, &
                               xlo,time,bc )&
                               bind(C, name="FORT_TEMPFILL")

      integer    DIMDEC(temp)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     temperature(DIMV(temp))
      integer    bc(SDIM,2)

      integer    i, j

#include <probdata.H>

      call filcc(temperature,DIMS(temp),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(temp).lt.domlo(1)) then
         do i = ARG_L1(temp), domlo(1)-1
           do j = ARG_L2(temp), ARG_H2(temp)
               temperature(i,j) = one
           end do
         end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(temp).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(temp)
           do j = ARG_L2(temp), ARG_H2(temp)
               temperature(i,j) = one
           end do
         end do
      end if    

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(temp).lt.domlo(2)) then
         do j = ARG_L2(temp), domlo(2)-1
           do i = ARG_L1(temp), ARG_H1(temp)
               temperature(i,j) = one
          end do
       end do
      end if    

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(temp).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(temp)
           do i = ARG_L1(temp), ARG_H1(temp)
               temperature(i,j) = one
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
!c :::        with valid data and that all non-interior cells have
!c ::         have been filled with a large real number.
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
      REAL_T     u(DIMV(u)), x_vel
      integer    lo(SDIM),hi(SDIM), bc(SDIM,2), i, j

#ifdef BL_DO_FLCT
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      integer DIMDEC(uflct)
      REAL_T  t_flct
      REAL_T, allocatable :: uflct(:,:)
#endif

#include <probdata.H>

#ifdef BL_DO_FLCT
#include <INFL_FORCE_F.H>
#endif

      lo(1) = ARG_L1(u)
      lo(2) = ARG_L2(u)
      hi(1) = ARG_H1(u)
      hi(2) = ARG_H2(u)

#ifdef BL_DO_FLCT
      if (forceInflow) then
         do i = 1, SDIM
            loFlctArray(i) = lo(i)
            hiFlctArray(i) = hi(i)
         end do
         loFlctArray(adv_dir) = 1
         hiFlctArray(adv_dir) = 1
         call SET_ARGS(DIMS(uflct), loFlctArray, hiFlctArray)
         allocate(uflct(DIMV(uflct)))
!c
!c     Note that we are 'scaling time' here to step into the fluct file to the
!c     correct depth.  This requires that time is not further scaled inside the
!c     the INFL_FILL routine.  Just to be sure, we set convVel = 1 here again.
!c
         convVel = one
         t_flct = adv_vel*time
         call INFL_FILL(FLCT_XVEL,DIMS(uflct),uflct,xlo,dx,t_flct,bc,f_problo,f_probhi)
      end if
#endif

      if (adv_dir .eq. 1)then
         x_vel = adv_vel
      else  
         x_vel = zero
      end if

      call filcc(u,DIMS(u),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(u).lt.domlo(1)) then
         do i = ARG_L1(u), domlo(1)-1
            do j = ARG_L2(u), ARG_H2(u)
               u(i,j) = x_vel
            end do
         end do
      end if            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(u).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(u)
            do j = ARG_L2(u), ARG_H2(u)
               u(i,j) = x_vel
            end do
         end do
      end if            

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(u).lt.domlo(2)) then
         do j = ARG_L2(u), domlo(2)-1
            do i = ARG_L1(u), ARG_H1(u)
#ifdef BL_DO_FLCT
               if (forceLo .and. adv_dir .eq. 2) then
                  u(i,j) = zero + uflct(i,1)*turb_scale
               else
                  u(i,j) = zero
               end if
#else
               u(i,j) = zero
#endif
            end do
         end do
      end if            

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(u).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(u)
            do i = ARG_L1(u), ARG_H1(u)
               if (probtype .eq. 10) then
!c ::: Lid-driven cavity test case, constant velocity on top of domain
                  u(i,j) = lid_vel
               else 
                  u(i,j) = zero
               end if
            end do
         end do
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
!c :::        with valid data and that all non-interior cells have
!c ::         have been filled with a large real number.
!c ::: 
!c ::: INPUTS/OUTPUTS:
!c ::: 
!c ::: v        <=  y velocity array
!c ::: DIMS(v)  => index extent of v array
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
      integer    bc(SDIM,2), i, j
      REAL_T     y_vel
      integer    lo(SDIM),hi(SDIM)

      REAL_T  x

#ifdef BL_DO_FLCT
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      integer DIMDEC(vflct)
      REAL_T, allocatable :: vflct(:,:)
      REAL_T  t_flct
#endif

#include <probdata.H>

#ifdef BL_DO_FLCT
#include <INFL_FORCE_F.H>
#endif

      lo(1) = ARG_L1(v)
      lo(2) = ARG_L2(v)
      hi(1) = ARG_H1(v)
      hi(2) = ARG_H2(v)

#ifdef BL_DO_FLCT
      if (forceInflow) then
         do i = 1, SDIM
            loFlctArray(i) = lo(i)
            hiFlctArray(i) = hi(i)
         end do
         loFlctArray(adv_dir) = 1
         hiFlctArray(adv_dir) = 1
         call SET_ARGS(DIMS(vflct), loFlctArray, hiFlctArray)
         allocate(vflct(DIMV(vflct)))
         convVel = one
         t_flct = adv_vel*time
         call INFL_FILL(FLCT_YVEL,DIMS(vflct),vflct,xlo,dx,t_flct,bc,f_problo,f_probhi)
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
           do j = ARG_L2(v),ARG_H2(v)
             v(i,j) = zero
           end do
         end do
      end if            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(v).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(v)
           do j = ARG_L2(v), ARG_H2(v)
             v(i,j) = zero
           end do
         end do
      end if   

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(v).lt.domlo(2)) then

         if (probtype .eq. 12) then

            do j = ARG_L2(v), domlo(2)-1
               do i = ARG_L1(v), ARG_H1(v)
              
                  x = xlo(1) + dx(1)*(float(i-lo(1)) + half)
                  if (x .gt. 0.45d0 .and.  &
                     x .lt. 0.55d0 .and. &
                     time .lt. 1.d0) then
                     v(i,j) = adv_vel
                  else
                     v(i,j) = 0.d0
                  end if
               end do
            end do

         else

            do j = ARG_L2(v), domlo(2)-1
               do i = ARG_L1(v), ARG_H1(v)
#ifdef BL_DO_FLCT
                  if (forceLo .and. adv_dir .eq. 2) then
                     v(i,j) = y_vel + vflct(i,1)*turb_scale
                  else
                     v(i,j) = y_vel
                  end if
#else
                  v(i,j) = y_vel
#endif
               end do
            end do

         end if            

      end if

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(v).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(v)
           do i = ARG_L1(v), ARG_H1(v)
             v(i,j) = y_vel
           end do
         end do
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
!c ::: p        <=  pressure array
!c ::: DIMS(p)   => index extent of p array
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of rho array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi) 
!c ::: -----------------------------------------------------------

      subroutine FORT_PRESFILL (p,DIMS(p),domlo,domhi,dx,xlo,time,bc)&
                                bind(C, name="FORT_PRESFILL")

      integer    DIMDEC(p)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     p(DIMV(p))
      integer    bc(SDIM,2)

      integer    i, j
      integer    ilo, ihi, jlo, jhi
      logical    fix_xlo, fix_xhi, fix_ylo, fix_yhi
      logical    per_xlo, per_xhi, per_ylo, per_yhi

      fix_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .ne. INT_DIR)
      per_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .eq. INT_DIR)
      fix_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .ne. INT_DIR)
      per_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .eq. INT_DIR)
      fix_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .ne. INT_DIR)
      per_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .eq. INT_DIR)
      fix_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .ne. INT_DIR)
      per_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .eq. INT_DIR)

      ilo = max(ARG_L1(p),domlo(1))
      ihi = min(ARG_H1(p),domhi(1))
      jlo = max(ARG_L2(p),domlo(2))
      jhi = min(ARG_H2(p),domhi(2))
!c
!c     ::::: left side
!c
      if (fix_xlo) then
         do i = ARG_L1(p), domlo(1)-1
            do j = jlo,jhi
               p(i,j) = p(ilo,j)
            end do
         end do
         if (fix_ylo) then
            do i = ARG_L1(p), domlo(1)-1
               do j = ARG_L2(p), domlo(2)-1
                  p(i,j) = p(ilo,jlo)
               end do
            end do
         else if (per_ylo) then
            do i = ARG_L1(p), domlo(1)-1
               do j = ARG_L2(p), domlo(2)-1
                  p(i,j) = p(ilo,j)
               end do
            end do
         end if
         if (fix_yhi) then
            do i = ARG_L1(p), domlo(1)-1
               do j = domhi(2)+1, ARG_H2(p)
                  p(i,j) = p(ilo,jhi)
               end do
            end do
         else if (per_yhi) then
            do i = ARG_L1(p), domlo(1)-1
               do j = domhi(2)+1, ARG_H2(p)
                  p(i,j) = p(ilo,j)
               end do
            end do
         end if
      end if
!c
!c     ::::: right side
!c
      if (fix_xhi) then
         do i = domhi(1)+1, ARG_H1(p)
            do j = jlo,jhi
               p(i,j) = p(ihi,j)
            end do
         end do
         if (fix_ylo) then
            do i = domhi(1)+1, ARG_H1(p)
               do j = ARG_L2(p), domlo(2)-1
                  p(i,j) = p(ihi,jlo)
               end do
            end do
         else if (per_ylo) then
            do i = domhi(1)+1, ARG_H1(p)
               do j = ARG_L2(p), domlo(2)-1
                  p(i,j) = p(ihi,j)
               end do
            end do
         end if
         if (fix_yhi) then
            do i = domhi(1)+1, ARG_H1(p)
               do j = domhi(2)+1, ARG_H2(p)
                  p(i,j) = p(ihi,jhi)
               end do
            end do
         else if (per_yhi) then
            do i = domhi(1)+1, ARG_H1(p)
               do j = domhi(2)+1, ARG_H2(p)
                  p(i,j) = p(ihi,j)
               end do
            end do
         end if
      end if
      
!c
!c     ::::: bottom 
!c
      if (fix_ylo) then
         do j = ARG_L2(p), domlo(2)-1
            do i = ilo, ihi
               p(i,j) = p(i,jlo)
            end do
         end do
         if (per_xlo) then
            do j = ARG_L2(p), domlo(2)-1
               do i = ARG_L1(p), domlo(1)-1
                  p(i,j) = p(i,jlo)
               end do
            end do
         end if
         if (per_xhi) then
            do j = ARG_L2(p), domlo(2)-1
               do i = domhi(1)+1, ARG_H1(p)
                  p(i,j) = p(i,jlo)
               end do
            end do
         end if
      end if

!c
!c      ::::: top
!c
      if (fix_yhi) then
         do j = domhi(2)+1, ARG_H2(p)
            do i = ilo, ihi
               p(i,j) = p(i,jhi)
            end do
         end do
         if (per_xlo) then
            do j = domhi(2)+1, ARG_H2(p)
               do i = ARG_L1(p), domlo(1)-1
                  p(i,j) = p(i,jhi)
               end do
            end do
         end if
         if (per_xhi) then
            do j = domhi(2)+1, ARG_H2(p)
               do i = domhi(1)+1, ARG_H1(p)
                  p(i,j) = p(i,jhi)
               end do
            end do
         end if
      end if

      end subroutine FORT_PRESFILL

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
!c ::: divu     <=  divu array
!c ::: DIMS(divu)=> index extent of p array
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of rho array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi) 
!c ::: -----------------------------------------------------------

      subroutine FORT_DIVUFILL (divu,DIMS(divu),domlo,domhi,delta, &
                               xlo,time,bc ) &
                               bind(C, name="FORT_DIVUFILL")

      integer    DIMDEC(divu)
      integer    bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     divu(DIMV(divu))

      integer    i, j
      integer    ilo, ihi, jlo, jhi
      REAL_T     y

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(divu)
      hi(1) = ARG_H1(divu)
      lo(2) = ARG_L2(divu)
      hi(2) = ARG_H2(divu)

      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))

      call filcc (divu,DIMS(divu),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
           if(jlo.le.jhi)then
             do j = jlo, jhi
               do i = lo(1), domlo(1)-1
                 divu(i,j) = divu(domlo(1),j)
               end do
             end do
           end if
           if (lo(2).lt.domlo(2)) then
             do j = lo(2), domlo(2)-1
               do i = lo(1), domlo(1)-1
                 divu(i,j) = divu(domlo(1),domlo(2))
               end do
             end do
           end if
           if(hi(2).gt.domhi(2))then
             do j = domhi(2)+1, hi(2)
               do i = lo(1), domlo(1)-1
                 divu(i,j) = divu(domlo(1),domhi(2))
               end do
             end do
           end if

      end if            

      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
           if(jlo.le.jhi)then
             do j = jlo,jhi
               do i = domhi(1)+1,hi(1)
                 divu(i,j) = divu(domhi(1),j)
               end do
             end do
           end if
           if (lo(2).lt.domlo(2)) then
             do j = lo(2), domlo(2)-1
               do i = domhi(1)+1,hi(1)
                 divu(i,j) = divu(domhi(1),domlo(2))
               end do
             end do
           end if
           if(hi(2).gt.domhi(2))then
             do j = domhi(2)+1, hi(2)
               do i = domhi(1)+1,hi(1)
                 divu(i,j) = divu(domhi(1),domhi(2))
               end do
             end do
           end if
      end if            

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
           if(ilo.le.ihi)then
             do j = lo(2), domlo(2)-1
               do i = ilo,ihi
                 divu(i,j) = divu(i,domlo(2))
               end do
             end do
           end if
           if (lo(1).lt.domlo(1)) then
             do j = lo(2), domlo(2)-1
               do i = lo(1), domlo(1)-1
                 divu(i,j) = divu(domlo(1),domlo(2))
               end do
             end do
           end if
           if(hi(1).gt.domhi(1))then
             do j = lo(2), domlo(2)-1
               do i = domhi(1)+1, hi(1)
                 divu(i,j) = divu(domhi(1),domlo(2))
               end do
             end do
           end if

      end if            

      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
           if(ilo.le.ihi)then
             do j = domhi(2)+1, hi(2)
               do i = ilo,ihi
                 divu(i,j) = divu(i,domhi(2))
               end do
             end do
           end if
           if (lo(1).lt.domlo(1)) then
             do j = domhi(2)+1, hi(2)
               do i = lo(1), domlo(1)-1
                 divu(i,j) = divu(domlo(1),domhi(2))
               end do
             end do
           end if
           if(hi(1).gt.domhi(1))then
             do j = domhi(2)+1, hi(2)
               do i = domhi(1)+1, hi(1)
                 divu(i,j) = divu(domhi(1),domhi(2))
               end do
             end do
           end if
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
!c ::: DIMS(dsdt)=> index extent of p array
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of rho array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi) 
!c ::: -----------------------------------------------------------


      subroutine FORT_DSDTFILL (dsdt,DIMS(dsdt),domlo,domhi,delta, &
                               xlo,time,bc ) &
                               bind(C, name="FORT_DSDTFILL")

      integer    DIMDEC(dsdt)
      integer    bc(SDIM,2)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     dsdt(DIMV(dsdt))

      integer    i, j
      integer    ilo, ihi, jlo, jhi
      REAL_T     y

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(dsdt)
      hi(1) = ARG_H1(dsdt)
      lo(2) = ARG_L2(dsdt)
      hi(2) = ARG_H2(dsdt)

      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))

      call filcc (dsdt,DIMS(dsdt),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
           do i = lo(1), domlo(1)-1
             do j = lo(2), hi(2)
               dsdt(i,j) = zero
             end do
           end do
      end if            

      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
           do i = domhi(1)+1, hi(1)
             do j = lo(2), hi(2)
               dsdt(i,j) = zero
             end do
           end do
      end if            

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
           do j = lo(2), domlo(2)-1
              do i = lo(1), hi(1)
                 dsdt(i,j) = zero
              end do
           end do
      end if            

      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
           do j = domhi(2)+1, hi(2)
              do i = lo(1), hi(1)
                 dsdt(i,j) = zero
              end do
           end do
      end if            

      end subroutine FORT_DSDTFILL

end module prob_2D_module
