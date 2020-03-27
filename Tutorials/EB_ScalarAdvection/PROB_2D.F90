
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

   public :: amrex_probinit, FORT_INITDATA,  &
   FORT_DENERROR, FORT_AVERAGE_EDGE_STATES, FORT_MAKEFORCE, &
   FORT_ADVERROR, FORT_ADV2ERROR, FORT_TEMPERROR, FORT_MVERROR, &
   FORT_DENFILL, FORT_ADVFILL, FORT_TEMPFILL, FORT_XVELFILL, &
   FORT_YVELFILL, FORT_PRESFILL, FORT_DIVUFILL, FORT_DSDTFILL

   !
   ! Define some parameters
   ! These could be loaded via a probin file instead of being hardcoded
   ! like we do here
   !
   REAL_T, parameter :: U0       = BL_REAL(4.0)
   REAL_T, parameter :: density  = BL_REAL(1.1798291685390001)
   REAL_T, parameter :: ADV0     = BL_REAL(4.0)

   !
   ! Set problo and probhi as module variables so that they can be
   ! used everywehere in this module
   !
   REAL_T :: m_problo(SDIM), m_probhi(SDIM)

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
      REAL_T  problo(SDIM), probhi(SDIM)

      !
      ! No need to read any probin file for this test
      ! We just Set module variables
      !
      m_problo = problo
      m_probhi = probhi

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

      integer i, j
      REAL_T  x, y, hx, hy

      hx = dx(1)
      hy = dx(2)

       do j = lo(2), hi(2)
          y = xlo(2) + hy*(float(j-lo(2)) + half)
          do i = lo(1), hi(1)
             x = xlo(1) + hx*(float(i-lo(1)) + half)

             vel(i,j,1) = U0
             vel(i,j,2) = zero

             ! Density
             scal(i,j,1) = density

             ! All other scalars
             ! We still initialize the tracers arrays to  zero
             scal(i,j,2:nscal) = zero

          end do
       end do


   end subroutine FORT_INITDATA

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

      integer   i

      ! For this case, refine the first half of the domain
      do i = lo(1), hi(1)
         if (i<domhi(1)/2) then
            tag(i,:) = set
         end if
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
   getForceVerbose) &
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
               do n=1, SDIM
                  velmin(n) = min(velmin(n),vel(i,j,k,n))
                  velmax(n) = max(velmax(n),vel(i,j,k,n))
               enddo
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

   end subroutine FORT_AVERAGE_EDGE_STATES

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
   dx,xlo,xhi,gravity,scomp,ncomp, &
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


      integer :: isioproc, n

      !
      ! No forcing for this problem
      !
      force = zero

      if (getForceVerbose > 0) then
         call bl_pd_is_ioproc(isioproc)
         if (isioproc.eq.1) then
            do n = 1, SDIM
               write (6,*) "No forcing applied"
            enddo
         endif
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

      integer :: isioproc

      call bl_pd_is_ioproc(isioproc)
      if (isioproc.eq.1) then
         stop "FORT_ADVERROR not implemented for this case"
      else
         stop
      endif

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

      integer :: isioproc

      call bl_pd_is_ioproc(isioproc)
      if (isioproc.eq.1) then
         stop "FORT_ADV2ERROR not implemented for this case"
      else
         stop
      endif

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

      integer :: isioproc

      call bl_pd_is_ioproc(isioproc)
      if (isioproc.eq.1) then
         stop "FORT_TEMPERROR not implemented for this case"
      else
         stop
      endif

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

      integer :: isioproc

      call bl_pd_is_ioproc(isioproc)
      if (isioproc.eq.1) then
         stop "FORT_TEMPERROR not implemented for this case"
      else
         stop
      endif

   end subroutine FORT_MVERROR

   !c ::: -----------------------------------------------------------
   !c ::: This routine is called during a fillpatch operation when
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

      ! filcc fills bc_types foextrap, hoextrap, reflect_odd and reflect_even
      call filcc(rho,DIMS(rho),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(rho).lt.domlo(1)) then
         do i = ARG_L1(rho), domlo(1)-1
            do j = ARG_L2(rho), ARG_H2(rho)
               rho(i,j) = density
            end do
         end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(rho).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(rho)
            do j = ARG_L2(rho), ARG_H2(rho)
               rho(i,j) = density
            end do
         end do
      end if

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(rho).lt.domlo(2)) then
         do j = ARG_L2(rho), domlo(2)-1
            do i = ARG_L1(rho), ARG_H1(rho)
               rho(i,j) = density
            end do
         end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(rho).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(rho)
            do i = ARG_L1(rho), ARG_H1(rho)
               rho(i,j) = density
            end do
         end do
      end if

   end subroutine FORT_DENFILL

   !c ::: -----------------------------------------------------------
   !c ::: This routine is called during a fillpatch operation when
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
   subroutine FORT_ADVFILL (adv,DIMS(adv),domlo,domhi,dx,&
   xlo,time,bc )&
   bind(C, name="FORT_ADVFILL")
      implicit none

      integer    DIMDEC(adv)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     adv(DIMV(adv))
      integer    bc(SDIM,2)
      integer    i,j

      call filcc(adv,DIMS(adv),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(adv).lt.domlo(1)) then
         do i = ARG_L1(adv), domlo(1)-1
            do j = ARG_L2(adv), ARG_H2(adv)
               adv(i,j) = ADV0
            end do
         end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(adv).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(adv)
            do j = ARG_L2(adv), ARG_H2(adv)
               adv(i,j) = ADV0
            end do
         end do
      end if

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(adv).lt.domlo(2)) then
         do j = ARG_L2(adv), domlo(2)-1
            do i = ARG_L1(adv), ARG_H1(adv)
               adv(i,j) = ADV0
            end do
         end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(adv).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(adv)
            do i = ARG_L1(adv), ARG_H1(adv)
               adv(i,j) = ADV0
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
      integer    i,j

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
   !c ::: This routine is called during a fillpatch operation when
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
   subroutine FORT_TEMPFILL (temp,DIMS(temp),domlo,domhi,dx,&
   xlo,time,bc )&
   bind(C, name="FORT_TEMPFILL")

      implicit none

      integer    DIMDEC(temp)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     temp(DIMV(temp))
      integer    bc(SDIM,2)

      ! Since we are not advecting any scalar besides density
      ! we do not need this subroutine to do anything

   end subroutine FORT_TEMPFILL

   !c ::: -----------------------------------------------------------
   !c ::: This routine is called during a fillpatch operation when
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
   subroutine FORT_XVELFILL (u,DIMS(u),domlo,domhi,dx,xlo,time,bc) &
   &          bind(C, name="FORT_XVELFILL")

      implicit none

      integer    DIMDEC(u)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     u(DIMV(u))
      integer    bc(SDIM,2)
      integer    lo(SDIM),hi(SDIM)
      integer    i, j

      lo(1) = ARG_L1(u)
      lo(2) = ARG_L2(u)
      hi(1) = ARG_H1(u)
      hi(2) = ARG_H2(u)

      call filcc(u,DIMS(u),domlo,domhi,dx,xlo,bc)

      !
      ! At the inlet we enforce a inflow condition
      ! At the outlet we do not need any condition on U since
      ! we have an outflow condition
      !
      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            do j = lo(2), hi(2)
               u(i,j) = U0
            end do
         end do
      end if

      ! At No-slip walls (= EXT_DIR for tangential component of velocity),
      ! set tangential velocity to zero
      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            do i = lo(1), hi(1)
               u(i,j) = zero
            end do
         end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            do i = lo(1), hi(1)
               u(i,j) = zero
            end do
         end do
      end if

      ! No need to fill along third direction because periodic conditions
      ! apply


   end subroutine FORT_XVELFILL

   !c ::: -----------------------------------------------------------
   !c ::: This routine is called during a fillpatch operation when
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
   subroutine FORT_YVELFILL (v,DIMS(v),domlo,domhi,dx,xlo,time,bc) &
   &          bind(C, name="FORT_YVELFILL")

      implicit none
      integer    DIMDEC(v)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     v(DIMV(v))
      integer    bc(SDIM,2)
      integer    lo(SDIM),hi(SDIM)
      integer    i, j, k

      lo(1) = ARG_L1(v)
      lo(2) = ARG_L2(v)
      hi(1) = ARG_H1(v)
      hi(2) = ARG_H2(v)

      call filcc(v,DIMS(v),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            do j = lo(2), hi(2)
               v(i,j) = zero
            end do
         end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            do j = lo(2), hi(2)
               v(i,j) = zero
            end do
         end do
      end if

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            do i = lo(1), hi(1)
               v(i,j) = zero
            end do
         end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            do i = lo(1), hi(1)
               v(i,j) = zero
            end do
         end do
      end if

      ! No need to fill along third direction because periodic conditions
      ! apply

   end subroutine FORT_YVELFILL

   !c ::: -----------------------------------------------------------
   !c ::: This routine is called during a fillpatch operation when
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

      ! we do not need this subroutine to do anything

   end subroutine FORT_PRESFILL

   !c ::: -----------------------------------------------------------
   !c ::: This routine is called during a fillpatch operation when
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

      ! Since we are not advecting any scalar besides density
      ! we do not need this subroutine to do anything

   end subroutine FORT_DIVUFILL

   !c ::: -----------------------------------------------------------
   !c ::: This routine is called during a fillpatch operation when
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

      ! Since we are not advecting any scalar besides density
      ! we do not need this subroutine to do anything

   end subroutine FORT_DSDTFILL

end module prob_2D_module
