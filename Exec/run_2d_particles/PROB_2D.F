
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>
#include <CONSTANTS.H>
#include <BC_TYPES.H>
#include <PROB_AMR_F.H>
#include <PROB_NS_F.H>
#include <ArrayLim.H>

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

#include <probdata.H>

c
c Dimensions of the Inflow file.
c
      integer nCompFile
      parameter (nCompFile = 2)

      namelist /fortin/ vorterr,denfact,velfact
c
c      Build "probin" filename -- the name of file containing fortin namelist.
c
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

      read(untin,fortin)
      if (isioproc .eq. 1) write(6,fortin)
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

      implicit none

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

c      call initvort(level,time,lo,hi,nscal,
c     &              vel,scal,DIMS(state),press,DIMS(press),
c     &              dx,xlo,xhi)

      call initpervort(level,time,lo,hi,nscal,
     &              vel,scal,DIMS(state),press,DIMS(press),
     &              dx,xlo,xhi)

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
      integer i, j, n
      REAL_T  x, y, r
      REAL_T  hx, hy

#include <probdata.H>

      hx = dx(1)
      hy = dx(2)

        do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half) - 0.5d0
            do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half) - 0.5d0
               r = sqrt(x**2 + y**2)

               vel(i,j,1) = 1.d0 
               vel(i,j,2) = 0.d0

               scal(i,j,1) = 1.d0

               do n = 2,nscal-1
                  scal(i,j,n) = one
               end do                  
               scal(i,j,nscal) = merge(one,zero,r.lt.0.25d0)
            end do
         end do

      end
c
c ::: -----------------------------------------------------------
c
      subroutine initpervort(level,time,lo,hi,nscal,
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
      integer i, j, n
      REAL_T  x, y, r
      REAL_T  hx, hy

#include <probdata.H>

      hx = dx(1)
      hy = dx(2)

      do j = lo(2), hi(2)
      do i = lo(1), hi(1)

            y = xlo(2) + hy*(float(j-lo(2)) + half)
            x = xlo(1) + hx*(float(i-lo(1)) + half)
            r = sqrt(x**2 + y**2)

            vel(i,j,1) = tanh(30.*(.25-abs(y-.5)))
            vel(i,j,2) = .05*sin(two*Pi*x)

            scal(i,j,1) = 1.d0

            do n = 2,nscal-1
               scal(i,j,n) = one
            end do                  
            scal(i,j,nscal) = merge(one,zero,r.lt.0.25d0)
      end do
      end do

      end
c
c
c ::: -----------------------------------------------------------
c
c     This routine averages the mac face velocities for makeforce at half time
c
      subroutine FORT_AVERAGE_EDGE_STATES(vel,umacx,umacy,
     &                                    DIMS(vel),DIMS(umacx),DIMS(umacy),
     &                                    getForceVerbose)

      implicit none

      integer    DIMDEC(vel)
      integer    DIMDEC(umacx)
      integer    DIMDEC(umacy)
      integer    getForceVerbose
      REAL_T     vel  (DIMV(vel),SDIM)
      REAL_T     umacx(DIMV(umacx))
      REAL_T     umacy(DIMV(umacy))

      integer i,j,n
      integer ilo,jlo
      integer ihi,jhi

      integer isioproc

      REAL_T  velmin(SDIM)
      REAL_T  velmax(SDIM)

      do n = 1, SDIM
         velmin(n) = 1.d234
         velmax(n) = -1.d234
      enddo

      ilo = vel_l1
      jlo = vel_l2
      ihi = vel_h1
      jhi = vel_h2

      do j = jlo, jhi
         do i = ilo, ihi
            vel(i,j,1) = half*(umacx(i,j)+umacx(i+1,j))
            vel(i,j,2) = half*(umacy(i,j)+umacy(i,j+1))
            do n=1, SDIM
               velmin(n)=min(velmin(n),vel(i,j,n))
               velmax(n)=max(velmax(n),vel(i,j,n))
            enddo
         enddo
      enddo
      
      if (getForceVerbose.gt.0) then
         call bl_pd_is_ioproc(isioproc)
         if (isioproc .eq. 1) then
            do n = 1, SDIM
               write (6,*) "mac velmin (",n,") = ",velmin(n)
               write (6,*) "mac velmax (",n,") = ",velmax(n)
            enddo
         endif
      endif

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

      integer   i, j

#include <probdata.H>

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            tag(i,j) = merge(set,tag(i,j),abs(vort(i,j,1)).gt.vorterr*2.d0**level)
         end do
      end do

      end

