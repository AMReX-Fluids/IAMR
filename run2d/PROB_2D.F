
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "PROB_AMR_F.H"
#include "PROB_PM_F.H"
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

#ifdef BL_DO_FLCT
#include "INFL_FORCE_F.H"
#include "infl_frc.H"
#endif

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
     &                  adv_dir, adv_vel, axis_dir, radvort
#ifdef BL_DO_FLCT
     &                  ,forceInflow, numInflPlanesStore, strmwse_dir,
     $                  forceLo, forceHi, flct_file, turb_scale
#endif
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
     &	 	               scal,DIMS(state),press,DIMS(press),
     &                         dx,xlo,xhi)
      integer    level, nscal
      integer    lo(SDIM),hi(SDIM)
      integer    DIMDEC(state)
      integer    DIMDEC(press)
      REAL_T     time, dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T    scal(DIMV(state),nscal)
      REAL_T   press(DIMV(press))

#include "probdata.H"

      if (probtype .eq. 1) then
         call initpm(level,time,lo,hi,nscal,
     &     	     scal,DIMS(state),press,DIMS(press),
     &               dx,xlo,xhi)
      else
         write(6,*) "INITDATA: bad probtype = ",probtype
         stop
      end if

      end
c
c ::: -----------------------------------------------------------
c
      subroutine initpm(level,time,lo,hi,nscal,
     &	 	        scal,DIMS(state),press,DIMS(press),
     &                  dx,xlo,xhi)

      integer    level, nscal
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(state)
      integer    DIMDEC(press)
      REAL_T     time, dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM)
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

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        scal(i,j,1) = 0.75d0
        scal(i,j,2) = 1.0d0 - scal(i,j,1)
        scal(i,j,3) = 1.0d0
        scal(i,j,4) = 1.0d0
        press(i,j) = 1.0d0
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

#if 0
      tag(0,13) = set
      tag(1,13) = set
      tag(1,12) = set
      tag(2,12) = set
      tag(3,12) = set
      tag(3,11) = set
      tag(4,11) = set
      tag(5,11) = set
      tag(5,10) = set
      tag(6,10) = set
      tag(6, 9) = set
      tag(7,10) = set
      tag(7, 9) = set
      tag(7, 8) = set
      tag(8, 8) = set
      tag(8, 7) = set
      tag(8, 6) = set
      tag(9, 6) = set
      tag(9, 5) = set
      tag(9, 4) = set
      tag(10, 4) = set
      tag(10, 3) = set
      tag(10, 2) = set
      tag(11, 2) = set
      tag(11, 1) = set
      tag(11, 0) = set
      tag(12, 0) = set

      return
#endif

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
                ax = abs(adv(i+1,j,1) - adv(i-1,j,1))
                ay = abs(adv(i,j+1,1) - adv(i,j-1,1))
                aerr = max(ax,ay)
                tag(i,j) = cvmgt(set,tag(i,j),aerr.gt.adverr)
!               tag(i,j) = cvmgt(set,tag(i,j),adv(i,j,1).gt.adverr)
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
              tag(i,j) = cvmgt(set,tag(i,j),abs(vort(i,j,1)).gt.vorterr*2.d0**level)
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

      subroutine FORT_N_FILL (rho,DIMS(rho),domlo,domhi,dx,
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

      subroutine FORT_ENTHFILL (enthalpy,DIMS(temp),domlo,domhi,dx,
     &                          xlo,time,bc )

      integer    DIMDEC(temp)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     enthalpy(DIMV(temp))
      integer    bc(SDIM,2)

      integer    i, j

#include "probdata.H"

      call filcc(enthalpy,DIMS(temp),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(temp).lt.domlo(1)) then
         do i = ARG_L1(temp), domlo(1)-1
           do j = ARG_L2(temp), ARG_H2(temp)
               enthalpy(i,j) = one
           end do
         end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(temp).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(temp)
           do j = ARG_L2(temp), ARG_H2(temp)
               enthalpy(i,j) = one
           end do
         end do
      end if    

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(temp).lt.domlo(2)) then
         do j = ARG_L2(temp), domlo(2)-1
           do i = ARG_L1(temp), ARG_H1(temp)
               enthalpy(i,j) = one
          end do
       end do
      end if    

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(temp).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(temp)
           do i = ARG_L1(temp), ARG_H1(temp)
               enthalpy(i,j) = one
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

