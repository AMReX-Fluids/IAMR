c
c $Id: PROB_2D.F,v 1.24 2005-12-15 21:52:16 lijewski Exp $
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

#define SDIM 2

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
#include "INFL_FORCE_F.H"
#include "infl_frc.H"

c
c Dimensions of the Inflow file.
c
      INTEGER dimFile(3)
      integer nCompFile
      parameter (nCompFile = 2)
      REAL_T dxFile(3)

      namelist /fortin/ denerr, vorterr, adverr, temperr,
     &			denfact, xblob, yblob, zblob, radblob, 
     &                  velfact, probtype, randfact, bubgrad,
     &			rhozero, tempzero, c_d, r_d, 
     &                  adv_dir, adv_vel, axis_dir, radvort,
     &                  forceInflow, numInflPlanesStore, strmwse_dir,
     $                  forceLo, forceHi, flct_file, turb_scale
c
c      Build "probin" filename -- the name of file containing fortin namelist.
c
      integer maxlen
      parameter (maxlen=256)

      character probin*(maxlen)

      if (namlen .gt. maxlen) then
         write(6,*) 'probin file name too long'
         stop
      end if

      do i = 1, namlen
         probin(i:i) = char(name(i))
      end do

      untin = 9
      if (namlen .eq. 0) then
         open(untin,file='probin',form='formatted',status='old')
      else
         open(untin,file=probin(1:namlen),form='formatted',status='old')
      end if

      forceInflow = .FALSE.
      numInflPlanesStore = 16
      forceLo = .TRUE.
      forceHi = .FALSE.
      strmwse_dir = FLCT_YVEL
      flct_file = ""
      turb_scale = 1

      read(untin,fortin)
c      write(6,fortin)
      close(unit=untin)

      if (forceInflow .eqv. .FALSE.) then
         forceLo = .FALSE.
         forceHi = .FALSE.
      else
         if (flct_file .ne. "") then
#define FF_UNIT 20
            print*, '...initializing turbulence, reading header info'
            open(FF_UNIT, file=flct_file, form='unformatted')
            call RD_SCL_FLCTHD(FF_UNIT,nCompFile,dimFile,probSizeFile,dxFile)
            close(FF_UNIT)
#undef  FF_UNIT
         endif
         if (strmwse_dir .NE. 2) then
            call bl_error('ERROR: turbulent inflow needs strmwse_dir=3')
         endif
c        print *,(dimFile(i),i=1,3)
c        print *,(probSizeFile(i),i=1,3)
c        print *,(dxFile(i),i=1,3)
      endif

      do i=1, SDIM
        f_problo(i) = problo(i)
        f_probhi(i) = probhi(i)
      enddo

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
c ::: dx       => cell size
c ::: xlo,xhi   => physical locations of lower left and upper
c :::              right hand corner of grid.  (does not include
c :::		   ghost region).
c ::: -----------------------------------------------------------
      subroutine FORT_INITDATA(level,time,lo,hi,nscal,
     &	 	               vel,scal,DIMS(state),press,DIMS(press),
     &                         dx,xlo,xhi)
      integer    level, nscal
      integer    lo(SDIM),hi(SDIM)
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
         call initviscbench(level,time,lo,hi,nscal,
     &     	            vel,scal,DIMS(state),press,DIMS(press),
     &                      dx,xlo,xhi)

      else
         write(6,*) "INITDATA: bad probtype = ",probtype
      end if

      end
c
c ::: -----------------------------------------------------------
c
      subroutine initbubble(level,time,lo,hi,nscal,
     &	 	            vel,scal,DIMS(state),press,DIMS(press),
     &                      dx,xlo,xhi)

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
      integer i, j, n
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  dist
      REAL_T  x_vel, y_vel

#include "probdata.H"

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
c           scal(i,j,1) = cvmgt(denfact,one,dist.lt.radblob)
            do n = 2,nscal-1
               scal(i,j,n) = one
            end do                  
            scal(i,j,nscal) = cvmgt(one,zero,dist.lt.radblob)
	 end do
      end do

      end
c
c ::: -----------------------------------------------------------
c
      subroutine initspin(level,time,lo,hi,nscal,
     &	 	          vel,scal,DIMS(state),press,DIMS(press),
     &                    dx,xlo,xhi)

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
      integer i, j, n
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  dist
      REAL_T  x_vel, y_vel
      REAL_T  spx, spy, cpx, cpy

#include "probdata.H"

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
	       scal(i,j,nscal) = cvmgt(one,zero,dist.lt.radblob)

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
     &	 	               vel,scal,DIMS(state),press,DIMS(press),
     &                         dx,xlo,xhi)

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
      integer i, j, n
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  spx, spy, cpx, cpy

#include "probdata.H"

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

      end
c
c ::: -----------------------------------------------------------
c
      subroutine initvort(level,time,lo,hi,nscal,
     &	 	          vel,scal,DIMS(state),press,DIMS(press),
     &                    dx,xlo,xhi)
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
      integer i, j, n
      REAL_T  x, y, r
      REAL_T  hx, hy
      REAL_T  c, ux, uy
      REAL_T  umagin, umagout, absu, sinth, costh
      REAL_T  small, a, b, r0

#include "probdata.H"

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
c              umagin = .5*r - 4*r**3
c              umagout = radvort*(.5*radvort - 4*radvort**3)/max(radvort,r)
               umagin = velfact * (one - a*(r - r0)**2)
               umagout = velfact * b/max(radvort,r)
               absu = cvmgp(umagout,umagin,r - radvort)
               sinth = y/max(r,small*radvort)
               costh = x/max(r,small*radvort)
               vel(i,j,1) = -absu*sinth
               vel(i,j,2) = absu*costh
               scal(i,j,1) = cvmgt(denfact,one,r.lt.radblob)
               do n = 2,nscal-1
                  scal(i,j,n) = one
               end do                  
               scal(i,j,nscal) = cvmgt(one,zero,r.lt.radblob)
	    end do
         end do

      end
c
c ::: -----------------------------------------------------------
c
      subroutine initchannel(level,time,lo,hi,nscal,
     &	 	             vel,scal,DIMS(state),press,DIMS(press),
     &                       dx,xlo,xhi)

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
      integer i, j, n
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  dist

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)

      do j = lo(2), hi(2)
         y = xlo(2) + hy*(float(j-lo(2)) + half)
         do i = lo(1), hi(1)
            vel(i,j,1) = adv_vel
            vel(i,j,2) = zero
            scal(i,j,1) = one

            do n = 2,nscal-1
               scal(i,j,n) = one
            end do                  

            x = xlo(1) + hx*(float(i-lo(1)) + half)
  	    dist = sqrt((x-xblob)**2 + (y-yblob)**2)
	    scal(i,j,nscal) = cvmgt(one,zero,dist.lt.radblob)

         end do
      end do

      end
c
c ::: -----------------------------------------------------------
c
      subroutine initpervort(level,time,lo,hi,nscal,
     &	 	             vel,scal,DIMS(state),press,DIMS(press),
     &                       dx,xlo,xhi)

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
      integer i, j, n
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  dist

#include "probdata.H"

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
	       scal(i,j,nscal) = cvmgt(one,zero,dist.lt.radblob)
	    end do
         end do

      end
c
c ::: -----------------------------------------------------------
c
      subroutine inithotspot(level,time,lo,hi,nscal,
     &	 	             vel,scal,DIMS(state),press,DIMS(press),
     &                       dx,xlo,xhi)

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
      integer i, j, n
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  x_vel, y_vel
      REAL_T  dist

#include "probdata.H"

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
            scal(i,j,1) = one/denfact + (one - one/denfact)
     &           *half*(one + tanh(40.*(dist - radblob)))
            scal(i,j,2) = cvmgt(one,zero,dist.lt.radblob)
            do n = 3,nscal-1
               scal(i,j,n) = one
            end do
            scal(i,j,nscal) = one / scal(i,j,1)
         end do
      end do
      
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
c ::: rho       => density array
c ::: DIMS(rho) => index extent of rho array
c ::: lo,hi     => index extent of grid
c ::: nvar      => number of components in rho array (should be 1)
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of tag array
c ::: problo    => phys loc of lower left corner of prob domain
c ::: time      => problem evolution time
c ::: -----------------------------------------------------------
      subroutine FORT_DENERROR (tag,DIMS(tag),set,clear,
     &                          rho,DIMS(rho),lo,hi,nvar,
     &                          domlo,domhi,dx,xlo,
     &			        problo,time,level)

      integer   DIMDEC(rho)
      integer   DIMDEC(tag)
      integer   lo(SDIM), hi(SDIM)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    rho(DIMV(rho), nvar)

      integer   i, j

#include "probdata.H"

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            tag(i,j) = cvmgt(set,tag(i,j),rho(i,j,1).lt.denerr)
	 end do
      end do

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
c ::: DIMS(adv) => index extent of scalar array
c ::: lo,hi     => index extent of grid
c ::: nvar      => number of components in rho array (should be 1)
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

#include "probdata.H"

c     probtype = SPIN
      if (probtype .eq. 1) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = cvmgt(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do

c     probtype = BUBBLE
      else if (probtype .eq. 2) then

        if (level .eq. 0) then
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                tag(i,j) = cvmgt(set,tag(i,j),adv(i,j,1).gt.adverr)
             end do
          end do
        end if

c     probtype = VORTEX IN A BOX
      else if (probtype .eq. 3) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = cvmgt(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do

c     probtype = CHANNEL
      else if (probtype .eq. 4) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = cvmgt(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do

c     probtype = PERIODIC SHEAR LAYER
      else if (probtype .eq. 5) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = cvmgt(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do

c     probtype = HOT SPOT
      else if (probtype .eq. 6) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = cvmgt(set,tag(i,j),adv(i,j,1).gt.adverr)
           end do
        end do


c     probtype = VISCOUS BENCHMARK
      else if (probtype .eq. 7) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = cvmgt(set,tag(i,j),adv(i,j,1).gt.adverr)
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
     &                          temperature,DIMS(temp),lo,hi,nvar,
     &                          domlo,domhi,dx,xlo,
     &                          problo,time,level)

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
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                tag(i,j) = cvmgt(set,tag(i,j),temperature(i,j,1).gt.temperr)
             end do
          end do
        else
c         ::::: refine where there is temperature gradient
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = abs(temperature(i+1,j,1) - temperature(i-1,j,1))
                ay = abs(temperature(i,j+1,1) - temperature(i,j-1,1))
                aerr = max(ax,ay)
                tag(i,j) = cvmgt(set,tag(i,j),aerr.gt.bubgrad)
             end do
          end do
        end if


c     probtype = VISCOUS BENCHMARK
      else if (probtype .eq. 7) then

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
c ::: vort      => array of vorticity values
c ::: DIMS(vor) => index extent of vort array
c ::: nvar      => number of components in vort array (should be 1)
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

#include "probdata.H"

c     probtype = SPIN
      if (probtype .eq. 1) then

c     probtype = BUBBLE
      else if (probtype .eq. 2) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = cvmgt(set,tag(i,j),abs(vort(i,j,1)).gt.vorterr)
           end do
        end do

c     probtype = VORTEX IN A BOX
      else if (probtype .eq. 3) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = cvmgt(set,tag(i,j),abs(vort(i,j,1)).gt.vorterr)
           end do
        end do

c     probtype = CHANNEL
      else if (probtype .eq. 4) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = cvmgt(set,tag(i,j),abs(vort(i,j,1)).gt.vorterr)
           end do
        end do

c     probtype = PERIODIC SHEAR LAYER
      else if (probtype .eq. 5) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = cvmgt(set,tag(i,j),abs(vort(i,j,1)).gt.vorterr)
           end do
        end do

c     probtype = HOT SPOT
      else if (probtype .eq. 6) then

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = cvmgt(set,tag(i,j),abs(vort(i,j,1)).gt.vorterr)
           end do
        end do

c     probtype = VISCOUS BENCHMARK
      else if (probtype .eq. 7) then
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              tag(i,j) = cvmgt(set,tag(i,j),abs(vort(i,j,1)).gt.vorterr)
           end do
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
c :::        with valid data and that all non-interior cells have
c ::         have been filled with a large real number.
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

      integer    DIMDEC(rho)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     rho(DIMV(rho))
      integer    bc(SDIM,2)

      integer    i, j

#include "probdata.H"

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
c :::        with valid data and that all non-interior cells have
c ::         have been filled with a large real number.
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

      subroutine FORT_ADVFILL (adv,DIMS(adv),domlo,domhi,dx,xlo,time,bc)

      integer    DIMDEC(adv)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     adv(DIMV(adv))
      integer    bc(SDIM,2)

      integer    i, j

#include "probdata.H"

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
c :::        with valid data and that all non-interior cells have
c ::         have been filled with a large real number.
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: temperature <=  temperature array
c ::: DIMS(temp)   => index extent of adv array
c ::: domlo,hi     => index extent of problem domain
c ::: dx           => cell spacing
c ::: xlo          => physical location of lower left hand
c :::                 corner of temperature array
c ::: time         => problem evolution time
c ::: bc           => array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_TEMPFILL (temperature,DIMS(temp),domlo,domhi,dx,
     &                          xlo,time,bc )

      integer    DIMDEC(temp)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     temperature(DIMV(temp))
      integer    bc(SDIM,2)

      integer    i, j

#include "probdata.H"

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
c :::        with valid data and that all non-interior cells have
c ::         have been filled with a large real number.
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
      REAL_T     u(DIMV(u)), x_vel
      integer    lo(SDIM),hi(SDIM), bc(SDIM,2), i, j

      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      integer DIMDEC(uflct)
      REAL_T  t_flct
      REAL_T, allocatable :: uflct(:,:)

#include "probdata.H"
#include "INFL_FORCE_F.H"

      lo(1) = ARG_L1(u)
      lo(2) = ARG_L2(u)
      hi(1) = ARG_H1(u)
      hi(2) = ARG_H2(u)

      if (forceInflow) then
         do i = 1, SDIM
            loFlctArray(i) = lo(i)
            hiFlctArray(i) = hi(i)
         end do
         loFlctArray(adv_dir) = 1
         hiFlctArray(adv_dir) = 1
         call SET_ARGS(DIMS(uflct), loFlctArray, hiFlctArray)
         allocate(uflct(DIMV(uflct)))
c
c     Note that we are 'scaling time' here to step into the fluct file to the
c     correct depth.  This requires that time is not further scaled inside the
c     the INFL_FILL routine.  Just to be sure, we set convVel = 1 here again.
c
         convVel = one
         t_flct = adv_vel*time
         call INFL_FILL(FLCT_XVEL,DIMS(uflct),uflct,xlo,dx,t_flct,bc,domnlo,domnhi)
      end if

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
               if (forceLo .and. adv_dir .eq. 2) then
                  u(i,j) = zero + uflct(i,1)*turb_scale
               else
                  u(i,j) = zero
               end if
	    end do
	 end do
      end if            

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(u).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(u)
            do i = ARG_L1(u), ARG_H1(u)
	       u(i,j) = zero
	    end do
	 end do
      end if    
        
      if (forceInflow) deallocate(uflct)

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
c :::        with valid data and that all non-interior cells have
c ::         have been filled with a large real number.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: v        <=  y velocity array
c ::: DIMS(v)  => index extent of v array
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
      integer    bc(SDIM,2), i, j
      REAL_T     y_vel

      integer lo(SDIM),hi(SDIM)
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      integer DIMDEC(vflct)
      REAL_T, allocatable :: vflct(:,:)
      REAL_T  t_flct

#include "probdata.H"
#include "INFL_FORCE_F.H"

      lo(1) = ARG_L1(v)
      lo(2) = ARG_L2(v)
      hi(1) = ARG_H1(v)
      hi(2) = ARG_H2(v)

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
         call INFL_FILL(FLCT_YVEL,DIMS(vflct),vflct,xlo,dx,t_flct,bc,domnlo,domnhi)
      end if

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
         do j = ARG_L2(v), domlo(2)-1
           do i = ARG_L1(v), ARG_H1(v)
              if (forceLo .and. adv_dir .eq. 2) then
                 v(i,j) = y_vel + vflct(i,1)*turb_scale
c                 if (j ==  domlo(2)-1) print*, i, j, vflct(i,1)
              else
                 v(i,j) = y_vel
              end if
	   end do
	 end do
      end if            

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(v).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(v)
           do i = ARG_L1(v), ARG_H1(v)
	     v(i,j) = y_vel
	   end do
	 end do
      end if            

      if (forceInflow) deallocate(vflct)

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
c ::: DIMS(p)   => index extent of p array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi) 
c ::: -----------------------------------------------------------

      subroutine FORT_PRESFILL (p,DIMS(p),domlo,domhi,dx,xlo,time,bc)

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
c
c     ::::: left side
c
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
c
c     ::::: right side
c
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
c ::: divu     <=  divu array
c ::: DIMS(divu)=> index extent of p array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi) 
c ::: -----------------------------------------------------------

      subroutine FORT_DIVUFILL (divu,DIMS(divu),domlo,domhi,delta,
     &                          xlo,time,bc )

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
c ::: DIMS(dsdt)=> index extent of p array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi) 
c ::: -----------------------------------------------------------


      subroutine FORT_DSDTFILL (dsdt,DIMS(dsdt),domlo,domhi,delta,
     &                          xlo,time,bc )

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

      end

