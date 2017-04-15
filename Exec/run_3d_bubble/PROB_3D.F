
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

c ::: -----------------------------------------------------------
c ::: This routine is called at problem initialization time
c ::: and when restarting from a checkpoint file.
c ::: The purpose is (1) to specify the initial time value
c ::: (not all problems start at time=0.0) and (2) to read
c ::: problem specific data from a namelist or other inputcdm
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
      REAL_T  thetaTmp, phiTmp
      REAL_T  cosThetaTmp, cosPhiTmp
      REAL_T  sinThetaTmp, sinPhiTmp
      REAL_T  px, py, pz, mp2, Ekh
      integer kx, ky, kz, mode_count, reduced_mode_count
      integer modx, mody, modz, xstep, ystep, zstep
      integer iseed

      REAL_T  Lx, Ly, Lz, Lmin
      REAL_T  kappa, kappaMax, freqMin, freqMax, freqDiff, pdk

      namelist /fortin/ adverr, xblob, yblob, zblob

      integer maxlen, isioproc, ierr
      parameter (maxlen=256)

      character probin*(maxlen)

      integer j, k, m, n

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
               dist = sqrt((x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2)
               vel(i,j,k,1) = 0.d0
               vel(i,j,k,2) = 0.d0
               vel(i,j,k,3) = 0.d0
  
c              This is density  
               scal(i,j,k,1) = 1.d0

c              This is the scalar that defines the bubble 
               do n = 2,nscal
                  scal(i,j,k,2) = 0.5d0 * (1.d0 - tanh((dist-1.d0)/0.05d0))
               end do                  
   	    end do
         end do
      end do

      end
