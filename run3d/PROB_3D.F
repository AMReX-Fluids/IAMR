c
c $Id: PROB_3D.F,v 1.7 1998-03-31 17:26:39 lijewski Exp $
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

#define SDIM 3

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
     &                  velfact, probtype, randfact, bubgrad,
     &			rhozero, tempzero, c_d, r_d, 
     &                  adv_dir, adv_vel, axis_dir, radvort

       untin = 9
       open(unit=untin, file='probin' ,form='formatted')
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
c ::: dx     => cell size
c ::: xlo,xhi   => physical locations of lower left and upper
c :::              right hand corner of grid.  (does not include
c :::		   ghost region).
c ::: -----------------------------------------------------------
      subroutine FORT_INITDATA(level,time,lo,hi,nscal,
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
      integer i, j, k
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  x_vel, y_vel, z_vel, dist

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
      else
         write(6,*) "INITDATA: bad probtype = ",probtype
      end if

      end

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
      integer i, j, k
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  rho, dist, c, ux, uy, uz, rhobub
      REAL_T  x_vel, y_vel, z_vel
      REAL_T  spx, spy, cpx, cpy, spz, cpz
      REAL_T  rn

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
      else if (adv_dir .eq. 3) then
         x_vel = zero
         y_vel = zero
         z_vel = adv_vel
      else 
         write(6,*) "initbubble: adv_dir = ",adv_dir
         stop
      end if

#if 0
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
	    spy = sin(Pi*y)
	    cpy = cos(Pi*y)
	    do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half)
 	       spx = sin(Pi*x)
	       cpx = cos(Pi*x)
  	       dist = sqrt((x-xblob)**2 + (y-yblob)**2)
	       c   = cvmgt(one,zero,dist.lt.radblob)
               ux = x_vel - velfact*two*spy*cpy*spx**2
               uy = y_vel + velfact*two*spx*cpx*spy**2
               uz = z_vel
               call blutilrand(rn)
               rhobub = denfact*(one + randfact*rn)
               rho = cvmgt(rhobub,one,dist.lt.radblob)
	       c   = cvmgt(one,zero,dist.lt.radblob)
	       vel(i,j,k,1) = ux
	       vel(i,j,k,2) = uy
	       vel(i,j,k,3) = uz
	       scal(i,j,k,1) = rho
               scal(i,j,k,2) = c
	    end do
         end do
      end do
#endif

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         do j = lo(2), hi(2)
            if (axis_dir .eq. 1) then
               y = xlo(2) + hy*(float(j-lo(2)) + half)
               dist = sqrt((z-zblob)**2 + (y-yblob)**2)
            end if
            do i = lo(1), hi(1)
               if (axis_dir .eq. 2) then
                  x = xlo(1) + hx*(float(i-lo(1)) + half)
                  dist = sqrt((z-zblob)**2 + (x-xblob)**2)
               end if
               ux = x_vel
               uy = y_vel
               uz = z_vel
               call blutilrand(rn)
               rhobub = denfact*(one + randfact*rn)
               rho = cvmgt(rhobub,one,dist.lt.radblob)
               c   = cvmgt(one,zero,dist.lt.radblob)
               vel(i,j,k,1) = ux
               vel(i,j,k,2) = uy
               vel(i,j,k,3) = uz
               scal(i,j,k,1) = rho
               scal(i,j,k,2) = c
            end do
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
      integer i, j, k
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  rho, dist, c, ux, uy, uz
      REAL_T  x_vel, y_vel, z_vel
      REAL_T  spx, spy, cpx, cpy, spz, cpz
      REAL_T  xlen, ylen, zlen

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
      else if (adv_dir .eq. 3) then
         x_vel = zero
         y_vel = zero
         z_vel = adv_vel
      else
         write(6,*) "initspin: adv_dir = ",adv_dir
         stop
      end if

      xlen = 1.d0
      ylen = 1.d0
      zlen = 1.d0

      if (axis_dir .eq. 3) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               y = xlo(2) + hy*(float(j-lo(2)) + half)
               spy = sin(Pi*y / ylen)
               cpy = cos(Pi*y / ylen)
               do i = lo(1), hi(1)
                  x = xlo(1) + hx*(float(i-lo(1)) + half)
                  spx = sin(Pi*x / xlen)
                  cpx = cos(Pi*x / xlen)
                  dist = sqrt((x-xblob)**2 + (y-yblob)**2)
                  rho = cvmgt(denfact,one,dist.lt.radblob)
                  c   = cvmgt(one,zero,dist.lt.radblob)
                  ux = x_vel - velfact*two*spy*cpy*spx**2
                  uy = y_vel + velfact*two*spx*cpx*spy**2
                  uz = z_vel
                  vel(i,j,k,1) = ux
                  vel(i,j,k,2) = uy
                  vel(i,j,k,3) = uz
                  scal(i,j,k,1) = rho
                  scal(i,j,k,2) = c
               end do
            end do
         end do
      else if (axis_dir .eq. 2) then
         do j = lo(2), hi(2)
            do k = lo(3), hi(3)
               z = xlo(3) + hz*(float(k-lo(3)) + half)
               spz = sin(Pi*z / zlen)
               cpz = cos(Pi*z / zlen)
               do i = lo(1), hi(1)
                  x = xlo(1) + hx*(float(i-lo(1)) + half)
                  spx = sin(Pi*x / xlen)
                  cpx = cos(Pi*x / xlen)
                  dist = sqrt((x-xblob)**2 + (z-zblob)**2)
                  rho = cvmgt(denfact,one,dist.lt.radblob)
                  c   = cvmgt(one,zero,dist.lt.radblob)
                  ux = x_vel - velfact*two*spz*cpz*spx**2
                  uy = y_vel
                  uz = z_vel + velfact*two*spx*cpx*spz**2
                  vel(i,j,k,1) = ux
                  vel(i,j,k,2) = uy
                  vel(i,j,k,3) = uz
                  scal(i,j,k,1) = rho
                  scal(i,j,k,2) = c
               end do
            end do
         end do
      else if (axis_dir .eq. 1) then
         do i = lo(1), hi(1)
            do k = lo(3), hi(3)
               z = xlo(3) + hz*(float(k-lo(3)) + half)
               spz = sin(Pi*z / zlen)
               cpz = cos(Pi*z / zlen)
               do j = lo(2), hi(2)
                  y = xlo(2) + hy*(float(j-lo(2)) + half)
                  spy = sin(Pi*y / ylen)
                  cpy = cos(Pi*y / ylen)
                  dist = sqrt((z-zblob)**2 + (y-yblob)**2)
                  rho = cvmgt(denfact,one,dist.lt.radblob)
                  c   = cvmgt(one,zero,dist.lt.radblob)
                  ux = x_vel
                  uy = y_vel + velfact*two*spz*cpz*spy**2
                  uz = z_vel - velfact*two*spy*cpy*spz**2
                  vel(i,j,k,1) = ux
                  vel(i,j,k,2) = uy
                  vel(i,j,k,3) = uz
                  scal(i,j,k,1) = rho
                  scal(i,j,k,2) = c
               end do
            end do
         end do
      else 
         print *,'BAD AXIS_DIR IN INITSPIN ',axis_dir
         stop
      end if

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
      integer i, j, k
      REAL_T  x, y, z, r
      REAL_T  hx, hy, hz
      REAL_T  c, ux, uy, uz
      REAL_T  umagin, umagout, absu, sinth, costh
      REAL_T  small

#include "probdata.H"

      hx = dx(1)
      hy = dx(2)
      hz = dx(2)
      small = 1.0e-10

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half) - yblob
	    do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half) - xblob
  	       r = sqrt(x**2 + y**2)
               umagin = .5*r - 4*r**3
               umagout = radvort*(.5*radvort - 4*radvort**3)/max(radvort,r)
               absu = cvmgp(umagout,umagin,r - radvort)
               sinth = y/max(r,small*radvort)
               costh = x/max(r,small*radvort)
               vel(i,j,k,1) = -absu*sinth
               vel(i,j,k,2) = absu*costh
               vel(i,j,k,3) = zero
               scal(i,j,k,1) = cvmgt(denfact,one,r.lt.radblob)
               scal(i,j,k,2) = cvmgt(one,zero,r.lt.radblob)
	    end do
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
c ::: nvar      => number of components in rho array (should be 1)
c ::: lo,hi     => index extent of grid
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

      integer   DIMDEC(tag)
      integer   DIMDEC(rho)
      integer   lo(SDIM), hi(SDIM)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    rho(DIMV(rho),nvar)

      integer   i, j, k

#include "probdata.H"

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j,k) = cvmgt(set,tag(i,j,k),rho(i,j,k,1).lt.denerr)
            end do
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

      if (level .eq. 0) then
c
c       ::::: refine around entire bubble
c
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),adv(i,j,k,1).gt.adverr)
               end do
            end do
         end do
      else
c
c       ::::: refine where there is bubble gradient
c
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = abs(adv(i+1,j,k,1) - adv(i-1,j,k,1))
                  ay = abs(adv(i,j+1,k,1) - adv(i,j-1,k,1))
                  az = abs(adv(i,j,k+1,1) - adv(i,j,k-1,1))
                  aerr = max(ax,ay,az)
                  tag(i,j,k) = cvmgt(set,tag(i,j,k),aerr.gt.bubgrad)
               end do
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

      integer   DIMDEC(tag)
      integer   DIMDEC(vort)
      integer   lo(SDIM), hi(SDIM)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    vort(DIMV(vort),nvar)

      REAL_T    x, y, dist
      integer   i, j, k
      REAL_T    radius

#include "probdata.H"

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j,k) = cvmgt(set,clear,abs(vort(i,j,k,1)).gt.vorterr)
            end do
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

      integer    DIMDEC(rho)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     rho(DIMV(rho))
      integer    bc(SDIM,2)

      integer    i, j, k

#include "probdata.H"

      call filcc(rho,DIMS(rho),domlo,domhi,dx,xlo,bc)

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

      integer    DIMDEC(adv)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     adv(DIMV(adv))
      integer    bc(SDIM,2)

      integer    i, j, k
      REAL_T     z

#include "probdata.H"

      call filcc(adv,DIMS(adv),domlo,domhi,dx,xlo,bc)

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
      end if 

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(adv).gt.domhi(3)) then
         do k = domhi(3)+1, ARG_H3(adv)
            do j = ARG_L2(adv), ARG_H2(adv)
               do i = ARG_L1(adv), ARG_H1(adv)
                  adv(i,j,k) = zero
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

      subroutine FORT_XVELFILL (u,DIMS(u),domlo,domhi,dx,
     &                          xlo,time,bc)

      integer    DIMDEC(u)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     u(DIMV(u))
      integer    bc(SDIM,2)

      integer    i, j, k
      REAL_T     x_vel

#include "probdata.H"

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
                  u(i,j,k) = x_vel
               end do
	    end do
	 end do
      end if            

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(u).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(u)
            do k = ARG_L3(u), ARG_H3(u)
               do i = ARG_L1(u), ARG_H1(u)
                  u(i,j,k) = x_vel
               end do
	    end do
	 end do
      end if            

      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(u).lt.domlo(3)) then
         do k = ARG_L3(u), domlo(3)-1
            do j = ARG_L2(u), ARG_H2(u)
               do i = ARG_L1(u), ARG_H1(u)
                  u(i,j,k) = x_vel
               end do
            end do
         end do
      end if

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(u).gt.domhi(3)) then
         do k = domhi(3)+1, ARG_H3(u)
            do j = ARG_L2(u), ARG_H2(u)
               do i = ARG_L1(u), ARG_H1(u)
                  u(i,j,k) = x_vel
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
c ::: v        <=  y velocity array
c ::: DIMS(v)   => index extent of v array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_YVELFILL (v,DIMS(v),domlo,domhi,dx,
     &                          xlo,time,bc)

      integer    DIMDEC(v)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     v(DIMV(v))
      integer    bc(SDIM,2)

      integer    i, j, k
      REAL_T     y_vel
      REAL_T     rn

#include "probdata.H"

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
                  v(i,j,k) = y_vel
               end do
	    end do
	 end do
      end if            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(v).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(v)
            do k = ARG_L3(v), ARG_H3(v)
               do j = ARG_L2(v), ARG_H2(v)
                  v(i,j,k) = y_vel
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
                  v(i,j,k) = y_vel
               end do
            end do
         end do
      end if        

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(v).gt.domhi(3)) then
         do k = domhi(3)+1, ARG_H3(v)
            do j = ARG_L2(v), ARG_H2(v)
               do i = ARG_L1(v), ARG_H1(v)
                  v(i,j,k) = y_vel
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
c ::: w        <=  z velocity array
c ::: DIMS(w)   => index extent of v array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_ZVELFILL (w,DIMS(w),domlo,domhi,dx,
     &                          xlo,time,bc)

      integer    DIMDEC(w)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     w(DIMV(w))
      integer    bc(SDIM,2)

      integer    i, j, k
      REAL_T     z_vel

#include "probdata.H"

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
                  w(i,j,k) = z_vel
               end do
	    end do
	 end do
      end if            

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(w).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(w)
            do k = ARG_L3(w), ARG_H3(w)
               do j = ARG_L2(w), ARG_H2(w)
                  w(i,j,k) = z_vel
               end do
	    end do
	 end do
      end if            

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(w).lt.domlo(2)) then
         do j = ARG_L2(w), domlo(2)-1
            do k = ARG_L3(w), ARG_H3(w)
               do i = ARG_L1(w), ARG_H1(w)
                  w(i,j,k) = z_vel
               end do
            end do
	 end do
      end if            

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(w).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(w)
            do k = ARG_L3(w), ARG_H3(w)
               do i = ARG_L1(w), ARG_H1(w)
                  w(i,j,k) = z_vel
               end do
            end do
	 end do
      end if            

      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(w).lt.domlo(3)) then
         do k = ARG_L3(w), domlo(3)-1
            do j = ARG_L2(w), ARG_H2(w)
               do i = ARG_L1(w), ARG_H1(w)
                  w(i,j,k) = z_vel
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

      integer    DIMDEC(p)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     p(DIMV(p))
      integer    bc(SDIM,2)

      integer    i, j, k
      integer    jlo, jhi, ilo, ihi, klo, khi
      logical    fix_xlo, fix_ylo, fix_zlo
      logical    fix_xhi, fix_yhi, fix_zhi

#include "probdata.H"

      fix_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .ne. INT_DIR)
      fix_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .ne. INT_DIR)
      fix_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .ne. INT_DIR)
      fix_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .ne. INT_DIR)
      fix_zlo = (ARG_L3(p) .lt. domlo(3)) .and. (bc(3,1) .ne. INT_DIR)
      fix_zhi = (ARG_H3(p) .gt. domhi(3)) .and. (bc(3,2) .ne. INT_DIR)

      ilo = max(ARG_L1(p),domlo(1))
      ihi = min(ARG_H1(p),domhi(1))
      jlo = max(ARG_L2(p),domlo(2))
      jhi = min(ARG_H2(p),domhi(2))
      klo = max(ARG_L3(p),domlo(3))
      khi = min(ARG_H3(p),domhi(3))

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
	    end if
	    if (fix_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,jlo,khi)
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
	    end if
	    if (fix_zhi) then
               do i = ARG_L1(p), domlo(1)-1
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ilo,jhi,khi)
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
	 end if
	 if (fix_zhi) then
	    do i = ARG_L1(p), domlo(1)-1
               do j = jlo, jhi
                  do k = domhi(3)+1, ARG_H3(p)
                     p(i,j,k) = p(ilo,j,khi)
                  end do
               end do
	    end do
	 end if
      end if            

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
	    end if
	    if (fix_zhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = ARG_L2(p), domlo(2)-1
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,jlo,khi)
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
	    end if
	    if (fix_zhi) then
               do i = domhi(1)+1, ARG_H1(p)
                  do j = domhi(2)+1, ARG_H2(p)
                     do k = domhi(3)+1, ARG_H3(p)
                        p(i,j,k) = p(ihi,jhi,khi)
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
	 end if
	 if (fix_zhi) then
	    do i = domhi(1)+1, ARG_H1(p)
               do j = jlo, jhi
                  do k = domhi(3)+1, ARG_H3(p)
                     p(i,j,k) = p(ihi,j,khi)
                  end do
               end do
	    end do
	 end if
      end if            

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
	 end if
	 if (fix_zhi) then
	    do j = ARG_L2(p), domlo(2)-1
               do k = domhi(3)+1, ARG_H3(p)
                  do i = ilo, ihi
                     p(i,j,k) = p(i,jlo,khi)
                  end do
               end do
	    end do
	 end if
      end if            

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
	 end if
	 if (fix_zhi) then
	    do j = domhi(2)+1, ARG_H2(p)
               do k = domhi(3)+1, ARG_H3(p)
                  do i = ilo, ihi
                     p(i,j,k) = p(i,jhi,khi)
                  end do
               end do
	    end do
	 end if
      end if            

      if (fix_zlo) then
         do k = ARG_L3(p), domlo(3)-1
            do j = jlo, jhi
               do i = ilo, ihi
                  p(i,j,k) = p(i,j,klo)
               end do
            end do
	 end do
      end if            

      if (fix_zhi) then
         do k = domhi(3)+1, ARG_H3(p)
            do j = jlo, jhi
               do i = ilo, ihi
                  p(i,j,k) = p(i,j,khi)
               end do
            end do
	 end do
      end if            

      end


      
