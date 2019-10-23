
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

  use amrex_fort_module, only : dim=>amrex_spacedim

  implicit none

  private

  public :: amrex_probinit, FORT_INITDATA, initvort, initpervort, &
            FORT_AVERAGE_EDGE_STATES, FORT_MVERROR, &
            FORT_DENFILL, FORT_ADVFILL, FORT_ADV2FILL, FORT_TEMPFILL, &
            FORT_XVELFILL, FORT_YVELFILL, FORT_PRESFILL, FORT_DIVUFILL, &
            FORT_DSDTFILL
  
contains
  
! c ::: -----------------------------------------------------------
! c ::: This routine is called at problem initialization time
! c ::: and when restarting from a checkpoint file.
! c ::: The purpose is (1) to specify the initial time value
! c ::: (not all problems start at time=0.0) and (2) to read
! c ::: problem specific data from a namelist or other input
! c ::: files and possibly store them or derived information
! c ::: in FORTRAN common blocks for later use.
! c ::: 
! c ::: INPUTS/OUTPUTS:
! c ::: 
! c ::: init      => TRUE if called at start of problem run
! c :::              FALSE if called from restart
! c ::: name      => name of "probin" file
! c ::: namlen    => length of name
! c ::: strttime <=  start problem with this time variable
! c ::: 
! c ::: -----------------------------------------------------------
      subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

      implicit none

      integer init,namlen
      integer name(namlen)
      integer untin, i
      REAL_T  problo(SDIM), probhi(SDIM)

#include <probdata.H>

!c
! Dimensions of the Inflow file.
!c
      integer nCompFile
      parameter (nCompFile = 2)

      namelist /fortin/ vorterr,denfact,velfact
!c
!     Build "probin" filename -- the name of file containing fortin namelist.
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

      read(untin,fortin)
      if (isioproc .eq. 1) write(6,fortin)
      close(unit=untin)

    end subroutine amrex_probinit

!::: -----------------------------------------------------------
!::: This routine is called at problem setup time and is used
!::: to initialize data on each grid.  The velocity field you
!::: provide does not have to be divergence free and the pressure
!::: field need not be set.  A subsequent projection iteration
!::: will define aa divergence free velocity field along with a
!::: consistant pressure.
!::: 
!::: NOTE:  all arrays have one cell of ghost zones surrounding
!:::        the grid interior.  Values in these cells need not
!:::        be set here.
!::: 
!::: INPUTS/OUTPUTS:
!::: 
!::: level     => amr level of grid
!::: time      => time at which to init data             
!::: lo,hi     => index limits of grid interior (cell centered)
!::: nscal     => number of scalar quantities.  You should know
!:::		   this already!
!::: vel      <=  Velocity array
!::: scal     <=  Scalar array
!::: press    <=  Pressure array
!::: dx       => cell size
!::: xlo,xhi   => physical locations of lower left and upper
!:::              right hand corner of grid.  (does not include
!:::		   ghost region).
!::: -----------------------------------------------------------
      subroutine FORT_INITDATA(level,time,lo,hi,nscal, &
     	 	               vel,scal,DIMS(state),press,DIMS(press),&
                              dx,xlo,xhi) bind(c, name="FORT_INITDATA")

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

!     call initvort(level,time,lo,hi,nscal, &
!                  vel,scal,DIMS(state),press,DIMS(press), &
!                  dx,xlo,xhi)

      call initpervort(level,time,lo,hi,nscal,&
                   vel,scal,DIMS(state),press,DIMS(press),&
                   dx,xlo,xhi)

    end subroutine FORT_INITDATA
!c
!::: -----------------------------------------------------------
!c
      subroutine initvort(level,time,lo,hi,nscal,&
     	 	          vel,scal,DIMS(state),press,DIMS(press),&
                         dx,xlo,xhi) bind(c, name="initvort")

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
!    ::::: local variables
!c
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

       end subroutine initvort
!c
!::: -----------------------------------------------------------
!c
      subroutine initpervort(level,time,lo,hi,nscal,&
     	 	             vel,scal,DIMS(state),press,DIMS(press),&
                            dx,xlo,xhi) bind(c,name="initpervort")

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
!    ::::: local variables
!c
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

    end subroutine initpervort

!=========================================================
!  This routine averages the mac face velocities for makeforce at half time
!=========================================================

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
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3), dim) :: vel
      REAL_T, dimension(ux_lo(1):ux_hi(1),ux_lo(2):ux_hi(2),ux_lo(3):ux_hi(3)) :: umacx
      REAL_T, dimension(uy_lo(1):uy_hi(1),uy_lo(2):uy_hi(2),uy_lo(3):uy_hi(3)) :: umacy
#if ( AMREX_SPACEDIM == 3 )
      REAL_T, dimension(uz_lo(1):uz_hi(1),uz_lo(2):uz_hi(2),uz_lo(3):uz_hi(3)) :: umacz
#endif

      REAL_T  :: velmin(3)
      REAL_T  :: velmax(3)
      integer :: isioproc

      integer :: i, j, k, n

      do n = 1, dim
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
               do n = 1, dim
                  velmin(n) = min(velmin(n),vel(i,j,k,n))
                  velmax(n) = max(velmax(n),vel(i,j,k,n))
               enddo
            enddo
         enddo
      enddo

      if (getForceVerbose.gt.0) then
         call bl_pd_is_ioproc(isioproc)
         if (isioproc.eq.1) then
            do n = 1, dim
               write (6,*) "mac velmin (",n,") = ",velmin(n)
               write (6,*) "mac velmax (",n,") = ",velmax(n)
            enddo
         endif
      endif

   end subroutine FORT_AVERAGE_EDGE_STATES

!::: -----------------------------------------------------------
!::: This routine will tag high error cells based on the 
!::: magnitude of vorticity
!::: 
!::: INPUTS/OUTPUTS:
!::: 
!::: tag      <=  integer tag array
!::: DIMS(tag) => index extent of tag array
!::: set       => integer value to tag cell for refinement
!::: clear     => integer value to untag cell
!::: vort      => array of vorticity values
!::: DIMS(vor) => index extent of vort array
!::: nvar      => number of components in vort array (should be 1)
!::: lo,hi     => index extent of grid
!::: domlo,hi  => index extent of problem domain
!::: dx        => cell spacing
!::: xlo       => physical location of lower left hand
!:::	           corner of tag array
!::: problo    => phys lo!of lower left corner of prob domain
!::: time      => problem evolution time
!::: -----------------------------------------------------------
      subroutine FORT_MVERROR (tag,DIMS(tag),set,clear,&
                              vort,DIMS(vort),lo,hi,nvar,&
                              domlo,domhi,dx,xlo,&
     			       problo,time,level) bind(c,name="FORT_MVERROR")

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

    end subroutine FORT_MVERROR

    
!::: -----------------------------------------------------------
!::: This routine is called during a filpatch operation when
!::: the patch to be filled falls outside the interior
!::: of the problem domain.  You are requested to supply the
!::: data outside the problem interior in such a way that the
!::: data is consistant with the types of the boundary conditions
!::: you specified in the C++ code.  
!::: 
!::: NOTE:  you can assume all interior cells have been filled
!:::        with valid data and that all non-interior cells have
!::         have been filled with a large real number.
!::: 
!::: INPUTS/OUTPUTS:
!::: 
!::: rho      <=  density array
!::: DIMS(rho) => index extent of rho array
!::: domlo,hi  => index extent of problem domain
!::: dx        => cell spacing
!::: xlo       => physical location of lower left hand
!:::	           corner of rho array
!::: time      => problem evolution time
!::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
!::: -----------------------------------------------------------

      subroutine FORT_DENFILL (rho,DIMS(rho),domlo,domhi,dx,&
                              xlo,time,bc) bind(c, name="FORT_DENFILL")

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

!::: -----------------------------------------------------------
!::: This routine is called during a filpatch operation when
!::: the patch to be filled falls outside the interior
!::: of the problem domain.  You are requested to supply the
!::: data outside the problem interior in such a way that the
!::: data is consistant with the types of the boundary conditions
!::: you specified in the C++ code.  
!::: 
!::: NOTE:  you can assume all interior cells have been filled
!:::        with valid data and that all non-interior cells have
!::         have been filled with a large real number.
!::: 
!::: INPUTS/OUTPUTS:
!::: 
!::: adv      <=  advected quantity array
!::: DIMS(adv) => index extent of adv array
!::: domlo,hi  => index extent of problem domain
!::: dx        => cell spacing
!::: xlo       => physical location of lower left hand
!:::	           corner of adv array
!::: time      => problem evolution time
!::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
!::: -----------------------------------------------------------

      subroutine FORT_ADVFILL (adv,DIMS(adv),domlo,domhi,dx,xlo,&
                               time,bc) bind(c,name="FORT_ADVFILL")

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

    end subroutine FORT_ADVFILL

      subroutine FORT_ADV2FILL (adv,DIMS(adv),domlo,domhi,dx,xlo,&
                                time,bc) bind(c,name="FORT_ADV2FILL")

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

!::: -----------------------------------------------------------
!::: This routine is called during a filpatch operation when
!::: the patch to be filled falls outside the interior
!::: of the problem domain.  You are requested to supply the
!::: data outside the problem interior in such a way that the
!::: data is consistant with the types of the boundary conditions
!::: you specified in the C++ code.
!:::
!::: NOTE:  you can assume all interior cells have been filled
!:::        with valid data and that all non-interior cells have
!::         have been filled with a large real number.
!:::
!::: INPUTS/OUTPUTS:
!:::
!::: temperature <=  temperature array
!::: DIMS(temp)   => index extent of adv array
!::: domlo,hi     => index extent of problem domain
!::: dx           => cell spacing
!::: xlo          => physical location of lower left hand
!:::                 corner of temperature array
!::: time         => problem evolution time
!::: bc           => array of boundary flags bc(BL_SPACEDIM,lo:hi)
!::: -----------------------------------------------------------

      subroutine FORT_TEMPFILL (temperature,DIMS(temp),domlo,domhi,dx,&
                               xlo,time,bc) bind(c,name="FORT_TEMPFILL")

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

!::: -----------------------------------------------------------
!::: This routine is called during a filpatch operation when
!::: the patch to be filled falls outside the interior
!::: of the problem domain.  You are requested to supply the
!::: data outside the problem interior in such a way that the
!::: data is consistant with the types of the boundary conditions
!::: you specified in the C++ code.  
!::: 
!::: NOTE:  you can assume all interior cells have been filled
!:::        with valid data and that all non-interior cells have
!::         have been filled with a large real number.
!::: 
!::: INPUTS/OUTPUTS:
!::: 
!::: u        <=  x velocity array
!::: DIMS(u)   => index extent of u array
!::: domlo,hi  => index extent of problem domain
!::: dx        => cell spacing
!::: xlo       => physical location of lower left hand
!:::	           corner of rho array
!::: time      => problem evolution time
!::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
!::: -----------------------------------------------------------

      subroutine FORT_XVELFILL (u,DIMS(u),domlo,domhi,dx,xlo,&
                                time,bc) bind(c,name="FORT_XVELFILL")

      implicit none

      integer    DIMDEC(u)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     u(DIMV(u))
      integer    lo(SDIM),hi(SDIM), bc(SDIM,2)

#include <probdata.H>

      lo(1) = ARG_L1(u)
      lo(2) = ARG_L2(u)
      hi(1) = ARG_H1(u)
      hi(2) = ARG_H2(u)

      call filcc(u,DIMS(u),domlo,domhi,dx,xlo,bc)

    end subroutine FORT_XVELFILL

!::: -----------------------------------------------------------
!::: This routine is called during a filpatch operation when
!::: the patch to be filled falls outside the interior
!::: of the problem domain.  You are requested to supply the
!::: data outside the problem interior in such a way that the
!::: data is consistant with the types of the boundary conditions
!::: you specified in the C++ code.  
!::: 
!::: NOTE:  you can assume all interior cells have been filled
!:::        with valid data and that all non-interior cells have
!::         have been filled with a large real number.
!::: 
!::: INPUTS/OUTPUTS:
!::: 
!::: v        <=  y velocity array
!::: DIMS(v)  => index extent of v array
!::: domlo,hi  => index extent of problem domain
!::: dx        => cell spacing
!::: xlo       => physical location of lower left hand
!:::	           corner of rho array
!::: time      => problem evolution time
!::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
!::: -----------------------------------------------------------

      subroutine FORT_YVELFILL (v,DIMS(v),domlo,domhi,dx,xlo,&
                                time,bc) bind(c, name="FORT_YVELFILL")

      implicit none

      integer    DIMDEC(v)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     v(DIMV(v))
      integer    bc(SDIM,2)

#include <probdata.H>

      call filcc(v,DIMS(v),domlo,domhi,dx,xlo,bc)

    end subroutine FORT_YVELFILL

!::: -----------------------------------------------------------
!::: This routine is called during a filpatch operation when
!::: the patch to be filled falls outside the interior
!::: of the problem domain.  You are requested to supply the
!::: data outside the problem interior in such a way that the
!::: data is consistant with the types of the boundary conditions
!::: you specified in the C++ code.  
!::: 
!::: NOTE:  you can assume all interior cells have been filled
!:::        with valid data.
!::: 
!::: INPUTS/OUTPUTS:
!::: 
!::: p        <=  pressure array
!::: DIMS(p)   => index extent of p array
!::: domlo,hi  => index extent of problem domain
!::: dx        => cell spacing
!::: xlo       => physical location of lower left hand
!:::	           corner of rho array
!::: time      => problem evolution time
!::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi) 
!::: -----------------------------------------------------------

      subroutine FORT_PRESFILL (p,DIMS(p),domlo,domhi,dx,xlo,&
                                time,bc) bind(c, name="FORT_PRESFILL")

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
!    ::::: left side
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
!    ::::: right side
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

    end subroutine FORT_PRESFILL

!::: -----------------------------------------------------------
!::: This routine is called during a filpatch operation when
!::: the patch to be filled falls outside the interior
!::: of the problem domain.  You are requested to supply the
!::: data outside the problem interior in such a way that the
!::: data is consistant with the types of the boundary conditions
!::: you specified in the C++ code.  
!::: 
!::: NOTE:  you can assume all interior cells have been filled
!:::        with valid data.
!::: 
!::: INPUTS/OUTPUTS:
!::: 
!::: p        <=  pressure array
!::: DIMS(p)   => index extent of p array
!::: domlo,hi  => index extent of problem domain
!::: dx        => cell spacing
!::: xlo       => physical location of lower left hand
!:::	           corner of rho array
!::: time      => problem evolution time
!::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi) 
!::: -----------------------------------------------------------

      subroutine FORT_DIVUFILL (divu,DIMS(divu),domlo,domhi,dx,xlo,&
                                time,bc) bind(c, name="FORT_DIVUFILL")

      integer    DIMDEC(divu)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     divu(DIMV(divu))
      integer    bc(SDIM,2)

#include <probdata.H>

      call filcc(divu,DIMS(divu),domlo,domhi,dx,xlo,bc)

    end subroutine FORT_DIVUFILL


!::: -----------------------------------------------------------
!::: This routine is called during a filpatch operation when
!::: the patch to be filled falls outside the interior
!::: of the problem domain.  You are requested to supply the
!::: data outside the problem interior in such a way that the
!::: data is consistant with the types of the boundary conditions
!::: you specified in the C++ code.  
!::: 
!::: NOTE:  you can assume all interior cells have been filled
!:::        with valid data.
!::: 
!::: INPUTS/OUTPUTS:
!::: 
!::: p        <=  pressure array
!::: DIMS(p)   => index extent of p array
!::: domlo,hi  => index extent of problem domain
!::: dx        => cell spacing
!::: xlo       => physical location of lower left hand
!:::	           corner of rho array
!::: time      => problem evolution time
!::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi) 
!::: -----------------------------------------------------------

      subroutine FORT_DSDTFILL (dsdt,DIMS(dsdt),domlo,domhi,dx,xlo,&
                                time,bc) bind(c, name="FORT_DSDTFILL")

      integer    DIMDEC(dsdt)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     dsdt(DIMV(dsdt))
      integer    bc(SDIM,2)

#include <probdata.H>

      call filcc(dsdt,DIMS(dsdt),domlo,domhi,dx,xlo,bc)

    end subroutine FORT_DSDTFILL

  end module prob_2D_module
