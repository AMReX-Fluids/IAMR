c
c $Id: PROB_2D.F,v 1.6 1998-03-31 17:26:39 lijewski Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "PROB_AMR_F.H"
#include "PROB_F.H"
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
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: init      => TRUE if called at start of problem run
c :::              FALSE if called from restart
c ::: strttime <=  start problem with this time variable
c ::: 
c ::: -----------------------------------------------------------
      subroutine FORT_PROBINIT (init)
      integer    init
      integer    untin, i

#include "probdata.H"

      namelist /fortin/ denerr, vorterr, adverr, temperr,
     &			denfact, xblob, yblob, zblob, radblob, 
     &                  xoffset, velfact, probtype, randfact, bubgrad,
     &			rhozero, tempzero, c_d, r_d, 
     &                  adv_dir, adv_vel, axis_dir, radvort

      untin = 9
      open(unit=untin, file='probin' ,form='formatted',status='old')
      read(untin,fortin)
      write(6,fortin)
      close(unit=untin)

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
c
c     ::::: local variables
c
      integer i, j
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  x_vel, y_vel, dist

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
         call initpipe(level,time,lo,hi,nscal,
     &     	       vel,scal,DIMS(state),press,DIMS(press),
     &                 dx,xlo,xhi)

      else if (probtype .eq. 5) then
         call initpervort(level,time,lo,hi,nscal,
     &     	          vel,scal,DIMS(state),press,DIMS(press),
     &                    dx,xlo,xhi)

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
      integer i, j
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  rho, dist, c, ux, uy, rhobub
      REAL_T  x_vel, y_vel
      REAL_T  rho1, rho2, rhoavg, rhodiff
      REAL_T  rn

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)

      rho1 = 999.2
      rho2 = 1.225

      rhoavg  = half * (rho1 + rho2)
      rhodiff = rho1 - rho2

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
	       c   = cvmgt(one,zero,dist.lt.radblob)
               ux = x_vel
               uy = y_vel
               call blutilrand(rn)
               rhobub = denfact*(one + randfact*rn)
               rho = cvmgt(rhobub,one,dist.lt.radblob)
c              rho = rhoavg + half * rhodiff *  tanh(5000.*(dist-radblob))
	       c   = cvmgt(one,zero,dist.lt.radblob)
	       vel(i,j,1) = ux
	       vel(i,j,2) = uy
	       scal(i,j,1) = rho
               scal(i,j,2) = c
	    end do
         end do

      end

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
      integer i, j
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  rho, dist, c, ux, uy
      REAL_T  x_vel, y_vel
      REAL_T  spx, spy, cpx, cpy
      REAL_T  xlen, ylen

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

      xlen = one
      ylen = one

         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
	    spy = sin(Pi*y / ylen)
	    cpy = cos(Pi*y / ylen)
	    do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half)
 	       spx = sin(Pi*x / xlen + Pi*xoffset) 
               cpx = cos(Pi*x / xlen + Pi*xoffset)
  	       dist = sqrt((x-xblob)**2 + (y-yblob)**2)
               rho = cvmgt(denfact,one,dist.lt.radblob)
c              rho = one + denfact * tanh(10.*(dist-radblob))
	       c   = cvmgt(one,zero,dist.lt.radblob)
               ux = x_vel - sign(1.d0,x)*velfact*two*spy*cpy*spx**2
               uy = y_vel + sign(1.d0,x)*velfact*two*spx*cpx*spy**2

               vel(i,j,1) = ux
               vel(i,j,2) = uy
               scal(i,j,1) = rho
	       scal(i,j,2) = c
	    end do
         end do

      end

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
      integer i, j
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
               scal(i,j,2) = cvmgt(one,zero,r.lt.radblob)
	    end do
         end do

      end

      subroutine initpipe(level,time,lo,hi,nscal,
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

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)

         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
              vel(i,j,1) = zero
              vel(i,j,2) = adv_vel
              scal(i,j,1) = one
            end do
         end do

         do n = 2,nscal
          do j = lo(2), hi(2)
            do i = lo(1), hi(1)
              scal(i,j,n) = one
            end do
          end do
         end do

      end

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
      integer i, j
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

  	       dist = sqrt((x-xblob)**2 + (y-yblob)**2)
	       scal(i,j,2) = cvmgt(one,zero,dist.lt.radblob)
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

c     if (level .eq. 0) then
c       ::::: refine around entire bubble
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          tag(i,j) = cvmgt(set,tag(i,j),adv(i,j,1).gt.adverr)
        end do
        end do
#if 0
      else
c       ::::: refine where there is bubble gradient
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          ax = abs(adv(i+1,j,1) - adv(i-1,j,1))
          ay = abs(adv(i,j+1,1) - adv(i,j-1,1))
          aerr = max(ax,ay)
          tag(i,j) = cvmgt(set,tag(i,j),aerr.gt.bubgrad)
        end do
        end do
      end if
#endif
 
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

      if (level .eq. 0) then
c
c       ::::: refine around entire bubble
c
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j) = cvmgt(set,tag(i,j),temperature(i,j,1).gt.temperr)
            end do
         end do
      else
c
c       ::::: refine where there is bubble gradient
c
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ax = abs(temperature(i+1,j,1) - temperature(i-1,j,1))
               ay = abs(temperature(i,j+1,1) - temperature(i,j-1,1))
               aerr = max(ax,ay)
               tag(i,j) = cvmgt(set,tag(i,j),aerr.gt.bubgrad)
            end do
         end do
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

      REAL_T    x, y, dist
      integer   i, j
      REAL_T    radius

#include "probdata.H"

      do j = lo(2), hi(2)
c         y = xlo(2) + dx(2)*(j - lo(2) + 0.5)
         do i = lo(1), hi(1)
c            x = xlo(1) + dx(1)*(i - lo(1) + 0.5)
            tag(i,j) = cvmgt(set,clear,abs(vort(i,j,1)).gt.vorterr)
	 end do
      end do

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
	       adv(i,j) = one
	    end do
	 end do

      end if            

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(adv).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(adv)
            do i = ARG_L1(adv), ARG_H1(adv)
	       adv(i,j) = one
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
               temperature(i,j) = zero
           end do
         end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(temp).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(temp)
           do j = ARG_L2(temp), ARG_H2(temp)
               temperature(i,j) = zero
           end do
         end do
      end if    

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(temp).lt.domlo(2)) then
         do j = ARG_L2(temp), domlo(2)-1
           do i = ARG_L1(temp), ARG_H1(temp)
               temperature(i,j) = zero
          end do
       end do
      end if    

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(temp).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(temp)
           do i = ARG_L1(temp), ARG_H1(temp)
               temperature(i,j) = zero
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

      integer    DIMDEC(u)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     u(DIMV(u))
      integer    bc(SDIM,2)

      integer    i, j
      REAL_T     x_vel

#include "probdata.H"

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
	       u(i,j) = zero
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

      integer    DIMDEC(v)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     v(DIMV(v))
      integer    bc(SDIM,2)

      integer    i, j
      REAL_T     y_vel

#include "probdata.H"

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
             v(i,j) = y_vel
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

      fix_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .ne. INT_DIR)
      fix_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .ne. INT_DIR)
      fix_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .ne. INT_DIR)
      fix_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .ne. INT_DIR)

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
         end if
         if (fix_yhi) then
            do i = ARG_L1(p), domlo(1)-1
               do j = domhi(2)+1, ARG_H2(p)
                  p(i,j) = p(ilo,jhi)
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
         end if
	 if (fix_yhi) then
	    do i = domhi(1)+1, ARG_H1(p)
               do j = domhi(2)+1, ARG_H2(p)
                  p(i,j) = p(ihi,jhi)
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
      end if

      if (fix_yhi) then
         do j = domhi(2)+1, ARG_H2(p)
            do i = ilo, ihi
               p(i,j) = p(i,jhi)
            end do
	 end do
      end if

      end

