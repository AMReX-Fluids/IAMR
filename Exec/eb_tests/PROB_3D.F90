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

module prob_3D_module

   implicit none

   private

   public :: amrex_probinit, FORT_INITDATA, FORT_DSDTFILL, &   
   &         FORT_ADVERROR, FORT_ADV2ERROR, FORT_TEMPERROR, FORT_MVERROR, &
   &         FORT_DENFILL, FORT_ADVFILL, FORT_TEMPFILL, FORT_XVELFILL, &
   &         FORT_YVELFILL, FORT_ZVELFILL, FORT_PRESFILL, FORT_DIVUFILL, &
   &         FORT_RGERROR


   ! Define some parameters
   ! These could be loaded via a probin file instead of being hardcoded
   REAL_T, parameter :: U0       = 1.0d0
   REAL_T, parameter :: density  = 1.0d0
   REAL_T, parameter :: temp_bot  = 293.0d0 !373.15
   REAL_T, parameter :: temp_top  = 293.0d0

   ! Set problo and probhi as module variables so that they can be
   ! used everywehere in this module
   REAL_T :: m_problo(SDIM), m_probhi(SDIM)

   
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

!   subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)
!      implicit none
!      integer init,namlen
!      integer name(namlen)
!      REAL_T  problo(SDIM), probhi(SDIM)

      !
      ! No need to read any probin file for this test
      ! We just Set module variables
      !
!      m_problo = problo
!      m_probhi = probhi

!    end subroutine amrex_probinit


!-----------------------------------------------------------------
subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)
!-----------------------------------------------------------------
      implicit none
      integer init,namlen
      integer name(namlen)
      integer untin, i
      REAL_T  problo(SDIM), probhi(SDIM)

#include <probdata.H>


      ! Dimensions of the Inflow file.
      integer nCompFile
      parameter (nCompFile = 2)
      namelist /fortin/ rgerr

      
      ! Build "probin" filename -- the name of file containing fortin namelist.
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


!      call initinflow(level,time,lo,hi,nscal, &
!           vel,scal,DIMS(state),press,DIMS(press), &
!           dx,xlo,xhi)
         

    end subroutine amrex_probinit



   !c ::: -----------------------------------------------------------
   !c ::: This routine is called at problem setup time and is used
   !c ::: to initialize data on each grid.  The velocity field you
   !c ::: provide does not have to be divergence free and the pressure
   !c ::: field need not be set.  A subsequent projection iteration
   !c ::: will define a divergence free velocity field along with a
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
      REAL_T   ATG, LTG

      integer i, j, k
      REAL_T  x, y, z, yn, hx, hy, hz

      
      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      do k = lo(3), hi(3)
         z = xlo(3) + hz*(float(k-lo(3)) + half)
         do j = lo(2), hi(2)
            y = xlo(2) + hy*(float(j-lo(2)) + half)
            do i = lo(1), hi(1)
               x = xlo(1) + hx*(float(i-lo(1)) + half)

!               yn = y / m_probhi(2)
!               vel(i,j,k,1) = 6.0d0 * U0 * yn * (1.0 - yn)
               yn = y - 1.0d0 ! channel goes from 0 to 2
               vel(i,j,k,1) = zero !1.1d0 * U0 * (1.0d0 - yn**4)
               vel(i,j,k,2) = zero
               vel(i,j,k,3) = zero

               ! add perturbations
!               ATG = U0/10.0d0
!               LTG = 1.0d0*3.14159265359
!               vel(i,j,k,1) = vel(i,j,k,1) + ATG * 1.0d0 * cos(LTG*x)*sin(LTG*yn)*sin(LTG*z)
!               vel(i,j,k,2) = vel(i,j,k,2) - ATG * 3.0d0 * sin(LTG*x)*cos(LTG*yn)*sin(LTG*z)
!               vel(i,j,k,3) = vel(i,j,k,3) + ATG * 2.0d0 * sin(LTG*x)*sin(LTG*yn)*cos(LTG*z)

               ! Density
               scal(i,j,k,1) = density

               ! All other scalars -- even if this case does not
               ! trace anything, we still initialize the tracers
               ! arrays to  zero
               scal(i,j,k,2:nscal) = zero

               ! Temp, not sure how this gets a 3...
               scal(i,j,k,3) = 293.0d0 !temp_bot + 0.5d0*y*(temp_top-temp_bot)

            end do
         end do
      end do
      

      vel(:,:,:,1) = 0.0d0
      vel(:,:,:,2) = 0.0d0
      vel(:,:,:,3) = 0.0d0      
      scal(:,:,:,1) = 1.0d0
      scal(:,:,:,2) = 0.0d0      
      scal(:,:,:,3) = 293.0d0      

   end subroutine FORT_INITDATA

!------------------------------------------------------------

   
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

      integer   i
      integer :: isioproc

      ! For this case, refine the first half of the domain
!      do i = lo(1), hi(1)
!         if (i<domhi(1)/2) then
!            tag(i,:,:) = set
!         end if
!      end do

      call bl_pd_is_ioproc(isioproc)
      if (isioproc.eq.1) then
         stop "FORT_DENERROR not implemented for this case"
      else
         stop
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

      integer   i, j, k
      integer   DIMDEC(tag)
      integer   DIMDEC(vort)
      integer   nvar, set, clear, level
      integer   lo(SDIM), hi(SDIM)
      integer   domlo(SDIM), domhi(SDIM)
      integer   tag(DIMV(tag))
      REAL_T    vort(DIMV(vort),nvar)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time


#include <probdata.H>

!      do j = lo(2), hi(2)
!         do i = lo(1), hi(1)
!            tag(i,j) = merge(set,tag(i,j),abs(vort(i,j,1)).gt.vorterr*2.d0**level)
!         end do
!      end do

!      print*, " HERE!", tag(i,j,k), vorterr
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j,k) = clear
               tag(i,j,k) = merge(set,tag(i,j,k),abs(vort(i,j,k,1)).gt.vorterr*2.0d0**level)
!               if (abs(vort(i,j,k,1)) .GT. vorterr*2.0d0**level) print*, vort(i,j,k,1), vorterr
            end do
         end do
      end do

    end subroutine FORT_MVERROR


!c ::: -----------------------------------------------------------
      subroutine FORT_RGERROR (tag,DIMS(tag),set,clear,&
                               re_g,DIMS(re_g),lo,hi,nvar,&
                               domlo,domhi,dx,xlo,&
     			       problo,time,level) bind(c,name="FORT_RGERROR")

      integer   i, j, k
      integer   DIMDEC(tag)
      integer   DIMDEC(re_g)
      integer   nvar, set, clear, level
      integer   lo(SDIM), hi(SDIM)
      integer   domlo(SDIM), domhi(SDIM)
      integer   tag(DIMV(tag))
      REAL_T    re_g(DIMV(re_g),nvar)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time


#include <probdata.H>

!      print*, " RGERROR DUMP:"
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
!               tag(i,j,k) = clear

               tag(i,j,k) = merge(set,tag(i,j,k),re_g(i,j,k,1) .GT. rgerr )
!               tag(i,j,k) = merge(set,tag(i,j,k),re_g(i,j,k,1) .GT. rgerr / (10.0d0**level))

!               if (re_g(i,j,k,1) .GT. rgerr) print*, "TEST:", level, re_g(i,j,k,1), rgerr, tag(i,j,k), clear, set
!               if (re_g(i,j,k,1) .GT. rgerr) then 
!                  print*, "LEVEL++:", level, "Re_G:", re_g(i,j,k,1), "TAG:", tag(i,j,k)
!               else
!                  print*, "LEVEL:", level, "Re_G:", re_g(i,j,k,1), "TAG:", tag(i,j,k)
!               endif
!               if (level .eq. 2) print*, "***LEVEL 2 ACTIVE***", re_g(i,j,k,1)
            enddo
         enddo
      enddo

    end subroutine FORT_RGERROR


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

      integer    i, j, k, n

      ! filcc fills bc_types foextrap, hoextrap, reflect_odd and reflect_even
      call filcc(rho,DIMS(rho),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(rho).lt.domlo(1)) then
         do i = ARG_L1(rho), domlo(1)-1
            do k = ARG_L3(rho), ARG_H3(rho)
               do j = ARG_L2(rho), ARG_H2(rho)
                  rho(i,j,k) = density
               end do
            end do
         end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(rho).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(rho)
            do k = ARG_L3(rho), ARG_H3(rho)
               do j = ARG_L2(rho), ARG_H2(rho)
                  rho(i,j,k) = density
               end do
            end do
         end do
      end if

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(rho).lt.domlo(2)) then
         do j = ARG_L2(rho), domlo(2)-1
            do k = ARG_L3(rho), ARG_H3(rho)
               do i = ARG_L1(rho), ARG_H1(rho)
                  rho(i,j,k) = density
               end do
            end do
         end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(rho).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(rho)
            do k = ARG_L3(rho), ARG_H3(rho)
               do i = ARG_L1(rho), ARG_H1(rho)
                  rho(i,j,k) = density
               end do
            end do
         end do
      end if

      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(rho).lt.domlo(3)) then
         do k = ARG_L3(rho), domlo(3)-1
            do j = ARG_L2(rho), ARG_H2(rho)
               do i = ARG_L1(rho), ARG_H1(rho)
                  rho(i,j,k) = density
               end do
            end do
         end do
      endif

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(rho).gt.domhi(3)) then
         do k = domhi(3)+1, ARG_H3(rho)
            do j = ARG_L2(rho), ARG_H2(rho)
               do i = ARG_L1(rho), ARG_H1(rho)
                  rho(i,j,k) = density
               end do
            end do
         end do
      end if

   end subroutine FORT_DENFILL


!-----------------------------------------------------------------
! Everything after here specific for imp-eff problem   
!-----------------------------------------------------------------   
   

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

   ! ADV is for tracer 1 and ADV2 is for tracer 2, i.e. idgaf
   subroutine FORT_ADVFILL (adv,DIMS(adv),domlo,domhi,dx,&
   xlo,time,bc )&
   bind(C, name="FORT_ADVFILL")
      implicit none

      integer    DIMDEC(adv)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     adv(DIMV(adv))
      integer    bc(SDIM,2)

      ! Since we are not advecting any scalar besides density
      ! we do not need this subroutine to do anything

      adv = 0.0d0

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

      ! Since we are not advecting any scalar besides density
      ! we do not need this subroutine to do anything

      adv = 0.0d0

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
                  temp(i,j,k) = temp_bot
               end do
	    end do
	 end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(temp).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(temp)
            do k = ARG_L3(temp), ARG_H3(temp)
               do j = ARG_L2(temp), ARG_H2(temp)
                  temp(i,j,k) = temp_top
               end do
	    end do
	 end do
      end if

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(temp).lt.domlo(2)) then
!      if (bc(2,1).eq.FOEXTRAP.and.ARG_L2(temp).lt.domlo(2)) then            
         do j = ARG_L2(temp), domlo(2)-1
            do k = ARG_L3(temp), ARG_H3(temp)
               do i = ARG_L1(temp), ARG_H1(temp)
                  temp(i,j,k) = temp_bot
               end do
	    end do
	 end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(temp).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(temp)
            do k = ARG_L3(temp), ARG_H3(temp)
               do i = ARG_L1(temp), ARG_H1(temp)
                  temp(i,j,k) = temp_top
               end do
	    end do
	 end do
      end if

      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(temp).lt.domlo(3)) then
         do k = ARG_L3(temp), domlo(3)-1
            do j = ARG_L2(temp), ARG_H2(temp)
               do i = ARG_L1(temp), ARG_H1(temp)
                  temp(i,j,k) = temp_bot
               end do
            end do
         end do
      end if

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(temp).gt.domhi(3)) then
         do k = domhi(3)+1, ARG_H3(temp)
            do j = ARG_L2(temp), ARG_H2(temp)
               do i = ARG_L1(temp), ARG_H1(temp)
                  temp(i,j,k) = temp_top
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

      REAL_T factor, twidth, toffset, ui, u_face
      REAL_T x_top, z_top, r_top, x_bot, z_bot, r_bot, zc

      REAL_T pert, unsteady, pi, Lf
      REAL_T x,r,u1,u2,u3,u_inf,eta      


#include <probdata.H>

!#ifdef BL_DO_FLCT
!      integer loFlctArray(SDIM), hiFlctArray(SDIM)
!      integer DIMDEC(uflct)
!      REAL_T  t_flct
!      REAL_T, allocatable :: uflct(:,:,:)
!#include <INFL_FORCE_F.H>
!#endif

!      parameter (constn=.22089323)

      lo(1) = ARG_L1(u)
      lo(2) = ARG_L2(u)
      lo(3) = ARG_L3(u)
      hi(1) = ARG_H1(u)
      hi(2) = ARG_H2(u)
      hi(3) = ARG_H3(u)

!#ifdef BL_DO_FLCT
!      if (forceInflow) then
!         do i = 1, 3
!            loFlctArray(i) = lo(i)
!            hiFlctArray(i) = hi(i)
!         end do
!         loFlctArray(adv_dir) = 1
!         hiFlctArray(adv_dir) = 1
!         call SET_ARGS(DIMS(uflct), loFlctArray, hiFlctArray)
!         allocate(uflct(DIMV(uflct)))
!!
!!        Note that we are 'scaling time' here to step into the fluct file to the
!!        correct depth.  This requires that time is not further scaled inside the
!!        the INFL_FILL routine.  Just to be sure, we set convVel = 1 here again.
!!
!         convVel = one
!         t_flct = time
!
!         call INFL_FILL(FLCT_XVEL,DIMS(uflct),uflct,xlo,dx,t_flct,bc,domnlo,domnhi)
!      end if
!#endif

!      xc = 0.5*(domnhi(1) - domnlo(1))
!      yc = 0.5*(domnhi(2) - domnlo(2))
!      zc = 0.5*(domnhi(3) - domnlo(3))      


      ! geom descriptions, should make this connected to EB2 file
      x_top = 2.5d0 ! center of top cylinder
      z_top = 1.0d0
      r_top = 0.5d0
      
      x_bot = 1.5d0 ! center of bottom cylinder
      z_bot = 1.0d0
      r_bot = 0.5d0      

      pi = 3.14159265359D0
      Lf = pi/r_top
      twidth = 1.0d0
      toffset = 2.5d0
      factor = 0.5d0 * (tanh(time/twidth-toffset)-1.0d0)

      unsteady = (sin((time/twidth)*pi) + 10.0d0)/10.0d0
      factor = factor * unsteady

      ui = 1.0d0 ! holder for fluctuations
      ui = factor * ui

      call filcc(u,DIMS(u),domlo,domhi,dx,xlo,bc)


      ! x-face low      
      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(u).lt.domlo(1)) then
         do i = ARG_L1(u), domlo(1)-1
            do k = ARG_L3(u), ARG_H3(u)
               do j = ARG_L2(u), ARG_H2(u)
                  !                  u(i,j,k) = 0.0d0
                  u(i,j,k) = u(i+1,j,k)                  
               end do
	    end do
	 end do
      end if

      ! x-face high      
      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(u).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(u)
            do k = ARG_L3(u), ARG_H3(u)
               do j = ARG_L2(u), ARG_H2(u)
                  !                  u(i,j,k) = 0.0d0
                  u(i,j,k) = u(i-1,j,k)                  
               end do
	    end do
	 end do
      endif

      ! y-face low OUTFLOW
!      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(u).lt.domlo(2)) then 
      if (bc(2,1).eq.FOEXTRAP.and.ARG_L2(u).lt.domlo(2)) then
         do j = ARG_L2(u), domlo(2)-1
            do k = ARG_L3(u), ARG_H3(u)
               z = domnlo(3) + (k+0.5)*dx(3)
               zc = z - z_bot            
               do i = ARG_L1(u), ARG_H1(u)
                  x = domnlo(1) + (i+0.5)*dx(1)
                  xc = x - x_bot
                  r = sqrt( xc*xc + zc*zc )
                  u(i,j,k) = 0.0d0
!                  if (r.lt.r_top) then
                     u_face = u(i,domlo(2),k)
!                     u(i,j,k) = u_face ! careful, this may not work => just simple zero gradient
!                  endif
!                     u(i,j,k) = 0.0d0                     
               end do
	    end do
	 end do
      end if

      ! y-face high INFLOW
      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(u).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(u)
            do k = ARG_L3(u), ARG_H3(u)
               z = domnlo(3) + (dble(k)+0.5)*dx(3)
               zc = z - z_top                              
               do i = ARG_L1(u), ARG_H1(u)
                  x = domnlo(1) + (dble(i)+0.5)*dx(1)
                  xc = x - x_top                  
                  r = sqrt( xc*xc + zc*zc )
                  u(i,j,k) = 0.0d0
                  if (r.lt.r_top) then

                    pert = -1.0d0 * cos(Lf*xc) * sin(Lf*time) * sin(Lf*zc)
                    pert = 0.06d0 * pert * factor
                    u(i,j,k) = pert * (1.0d0-(r/r_top)**2)
                    
                  endif
                  
               end do
	    end do
	 end do
      end if



      ! z-face low        
      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(u).lt.domlo(3)) then
         do k = ARG_L3(u), domlo(3)-1
            do j = ARG_L2(u), ARG_H2(u)
               do i = ARG_L1(u), ARG_H1(u)
                  !                  u(i,j,k) = 0.0d0
                  u(i,j,k) = u(i,j,k+1)                  
               end do
            end do
         end do
      end if

      ! z-face high
      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(u).gt.domhi(3)) then
         do k = domhi(3)+1, ARG_H3(u)
            do j = ARG_L2(u), ARG_H2(u)
               do i = ARG_L1(u), ARG_H1(u)
                  !                  u(i,j,k) = 0.0d0
                  u(i,j,k) = u(i,j,k-1)           
               end do
            end do
         end do
      end if



!#ifdef BL_DO_FLCT
!      if (forceInflow) deallocate(uflct)
!#endif

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
      REAL_T factor, twidth, toffset, vi, v_face
      REAL_T x_top, z_top, r_top, x_bot, z_bot, r_bot, zc
      REAL_T uA, gvar, pi, pert, unsteady, Lf
      REAL_T x,y,z,r,u1,u2,u3,u_inf,eta
      REAL_T Qin, Aout      
      
#include <probdata.H>

!#ifdef BL_DO_FLCT
!      integer loFlctArray(SDIM), hiFlctArray(SDIM)
!      integer DIMDEC(vflct)
!      REAL_T  t_flct
!      REAL_T, allocatable :: vflct(:,:,:)
!#include <INFL_FORCE_F.H>
!#endif
      

      lo(1) = ARG_L1(v)
      lo(2) = ARG_L2(v)
      lo(3) = ARG_L3(v)
      hi(1) = ARG_H1(v)
      hi(2) = ARG_H2(v)
      hi(3) = ARG_H3(v)

!#ifdef BL_DO_FLCT
!      if (forceInflow) then
!         do i = 1, 3
!            loFlctArray(i) = lo(i)
!            hiFlctArray(i) = hi(i)
!         end do
!         loFlctArray(adv_dir) = 1
!         hiFlctArray(adv_dir) = 1
!         call SET_ARGS(DIMS(vflct), loFlctArray, hiFlctArray)
!         allocate(vflct(DIMV(vflct)))
!         convVel = one
!         t_flct = time
!         call INFL_FILL(FLCT_YVEL,DIMS(vflct),vflct,xlo,dx,t_flct,bc,domnlo,domnhi)
!      end if
!#endif


!      xc = 0.5*(domnhi(1) - domnlo(1))
!      yc = 0.5*(domnhi(2) - domnlo(2))
!      zc = 0.5*(domnhi(3) - domnlo(3))      

! GET THIS FROM INPUT
      
      ! geom descriptions, should make this connected to EB2 file
      x_top = 0.75d0 ! center of top cylinder 
      z_top = 1.25d0
      r_top = 0.5d0
      
      x_bot = 2.0d0 ! center of bottom cylinder
      z_bot = 2.0d0
      r_bot = 0.5d0

      gvar = r_top/10.0d0

      pi = 3.14159265359D0
      Lf = pi/r_top      
      twidth = 4.0d0 !1.0d0
      toffset = 2.5d0 !2.5d0
      factor = 0.5d0 * (tanh(time/twidth-toffset)+1.0d0)

      unsteady = (sin((time/twidth)*pi) + 10.0d0)/10.0d0
      factor = factor * unsteady

      vi = -1.0d0
      vi = factor * vi


      ! base profile:  v(i,j,k) = vi * (1.0d0-(r/r_top)**4)
      ! (UA)_in = [ r - (r^5/5) *(1/r_top**4) ] * 2pi
      Qin = vi * ( r_top - r_top/5.0d0 ) * 2.0d0 * pi
      Aout = 6.0d0 ! fix these hardcodes
      
      
      call filcc(v,DIMS(v),domlo,domhi,dx,xlo,bc)


      ! x-face low      
      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(v).lt.domlo(1)) then
         do i = ARG_L1(v), domlo(1)-1
            do k = ARG_L3(v), ARG_H3(v)
               do j = ARG_L2(v), ARG_H2(v)
                  !                  v(i,j,k) = 0.0d0
                  v(i,j,k) = v(i+1,j,k) ! fix this
               end do
	    end do
	 end do
      end if

      ! x-face high      
      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(v).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(v)
            do k = ARG_L3(v), ARG_H3(v)
               do j = ARG_L2(v), ARG_H2(v)
                  !                  v(i,j,k) = 0.0d0
                  v(i,j,k) = v(i-1,j,k)                  
               end do
	    end do
	 end do
      endif

      ! y-face low OUTFLOW, apparently this doesnt matterd
!      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(v).lt.domlo(2)) then
       if (bc(2,1).eq.FOEXTRAP .AND. ARG_L2(v).lt.domlo(2)) then
         !         print*, ">>> Setting outflow ghost"
         do j = ARG_L2(v), domlo(2)-1
            do k = ARG_L3(v), ARG_H3(v)
               z = domnlo(3) + (dble(k)+0.5)*dx(3)
               zc = z - z_bot            
               do i = ARG_L1(v), ARG_H1(v)
                  x = domnlo(1) + (dble(i)+0.5)*dx(1)
                  xc = x - x_bot
                  r = sqrt( xc*xc + zc*zc )
!                  v(i,j,k) = 0.0d0
                  v_face = v(i,domlo(2),k)
!                  if (r.lt.r_bot) v_face = v(i,domlo(2),k)                  
!                     print*, ">>> Setting outflow ghost", v_face                     
!                     v_face = vi * (r_top/r_bot)**2 * (1.0d0-r/r_bot)**2
!                     v(i,j,k) = v_face ! careful, this may not work => just simple zero gradient
!                     !                     v(i,domlo(2),k) = 0.5d0 * ( v_face + v(i,domlo(2),k) )
!                  endif
!                  v(i,j,k) = min(v_face,0.0d0)
!                  v(i,j+1,k) = min(v(i,j+1,k), 0.0d0)
!!!                  v(i,j,k) = v_face
!                  v(i,j,k) = Qin/Aout
               end do
	    end do
	 end do
      end if

      ! y-face high INFLOW
      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(v).gt.domhi(2)) then
!         print*, ">>> Setting inflow ghost"         
         do j = domhi(2)+1, ARG_H2(v)
            do k = ARG_L3(v), ARG_H3(v)
               z = domnlo(3) + (dble(k)+0.5)*dx(3)
               zc = z - z_top                              
               do i = ARG_L1(v), ARG_H1(v)
                  x = domnlo(1) + (dble(i)+0.5)*dx(1)
                  xc = x - x_top                  
                  r = sqrt( xc*xc + zc*zc )
                  v(i,j,k) = 0.0d0
                  if (r.lt.r_top) then

                     v(i,j,k) = vi * (1.0d0-(r/r_top)**4) ! for now no fluc

                     pert = 3.0d0 * sin(Lf*xc) * cos(Lf*time) * sin(Lf*zc)
                     pert = 0.06d0 * pert * factor
                     v(i,j,k) = v(i,j,k) + pert * (1.0d0-(r/r_top)**2)
                     
                  endif
               end do
	    end do
	 end do
      end if

      ! z-face low        
      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(v).lt.domlo(3)) then
         do k = ARG_L3(v), domlo(3)-1
            do j = ARG_L2(v), ARG_H2(v)
               do i = ARG_L1(v), ARG_H1(v)
                  !                  v(i,j,k) = 0.0d0
                  v(i,j,k) = v(i,j,k+1)                  
               end do
            end do
         end do
      end if

      ! z-face high
      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(v).gt.domhi(3)) then
         do k = domhi(3)+1, ARG_H3(v)
            do j = ARG_L2(v), ARG_H2(v)
               do i = ARG_L1(v), ARG_H1(v)
                  !                  v(i,j,k) = 0.0d0
                  v(i,j,k) = v(i,j,k-1)                  
               end do
            end do
         end do
      end if


!#ifdef BL_DO_FLCT
!      if (forceInflow) deallocate(uflct)
!#endif

      
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

      REAL_T x,y,z,r,u1,u2,u3,u_inf,eta
      REAL_T Lx, Ly
      REAL_T factor, twidth, toffset, w_face
      REAL_T x_top, z_top, r_top, x_bot, z_bot, r_bot, zc
      REAL_T pert, unsteady, pi, Lf
      

#include <probdata.H>

!#ifdef BL_DO_FLCT
!      REAL_T t_flct
!      integer loFlctArray(SDIM), hiFlctArray(SDIM)
!      integer DIMDEC(wflct)
!      REAL_T, allocatable :: wflct(:,:,:)
!#include <INFL_FORCE_F.H>
!#endif

      lo(1) = ARG_L1(w)
      lo(2) = ARG_L2(w)
      lo(3) = ARG_L3(w)
      hi(1) = ARG_H1(w)
      hi(2) = ARG_H2(w)
      hi(3) = ARG_H3(w)

!#ifdef BL_DO_FLCT
!      if (forceInflow) then
!         do i = 1, 3
!            loFlctArray(i) = lo(i)
!            hiFlctArray(i) = hi(i)
!         end do
!         loFlctArray(adv_dir) = 1
!         hiFlctArray(adv_dir) = 1
!         call SET_ARGS(DIMS(wflct), loFlctArray, hiFlctArray)
!         allocate(wflct(DIMV(wflct)))
!         convVel = one
!         t_flct = time
!         call INFL_FILL(FLCT_ZVEL,DIMS(wflct),wflct,xlo,dx,t_flct,bc,domnlo,domnhi)
!      end if
!#endif


!      xc = 0.5*(domnhi(1) - domnlo(1))
!      yc = 0.5*(domnhi(2) - domnlo(2))
!      zc = 0.5*(domnhi(3) - domnlo(3))      


      ! geom descriptions, should make this connected to EB2 file
      x_top = 2.5d0 ! center of top cylinder
      z_top = 1.0d0
      r_top = 0.5d0
      
      x_bot = 1.5d0 ! center of bottom cylinder
      z_bot = 1.0d0
      r_bot = 0.5d0

      pi = 3.14159265359D0
      Lf = pi/r_top      
      twidth = 1.0d0
      toffset = 2.5d0
      factor = 0.5d0 * (tanh(time/twidth-toffset)+1.0d0)

      unsteady = (sin((time/twidth)*pi) + 10.0d0)/10.0d0
      factor = factor * unsteady

      !wi = -1.0d0
      !wi = factor * wi
      


      
      call filcc(w,DIMS(w),domlo,domhi,dx,xlo,bc)


      ! x-face low      
      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(w).lt.domlo(1)) then
         do i = ARG_L1(w), domlo(1)-1
            do k = ARG_L3(w), ARG_H3(w)
               do j = ARG_L2(w), ARG_H2(w)
                  !                  w(i,j,k) = 0.0d0
                  w(i,j,k) = w(i+1,j,k)                  
               end do
	    end do
	 end do
      end if

      ! x-face high      
      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(w).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(w)
            do k = ARG_L3(w), ARG_H3(w)
               do j = ARG_L2(w), ARG_H2(w)
                  !                  w(i,j,k) = 0.0d0
                  w(i,j,k) = w(i-1,j,k)                  
               end do
	    end do
	 end do
      endif

      ! y-face low OUTFLOW
!      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(w).lt.domlo(2)) then
      if (bc(2,1).eq.FOEXTRAP.and.ARG_L2(w).lt.domlo(2)) then            
         do j = ARG_L2(w), domlo(2)-1
            do k = ARG_L3(w), ARG_H3(w)
               z = domnlo(3) + (k+0.5)*dx(3)
               zc = z - z_bot            
               do i = ARG_L1(w), ARG_H1(w)
                  x = domnlo(1) + (i+0.5)*dx(1)
                  xc = x - x_bot
                  r = sqrt( xc*xc + zc*zc )
                  w(i,j,k) = 0.0d0
!                  if (r.lt.r_top) then
                     w_face = w(i,domlo(2),k)
!                     w(i,j,k) = w_face ! careful, this may not work => just simple zero gradient
!                  endif
!                     w(i,j,k) = 0.0d0                     
               end do
	    end do
	 end do
      end if

      ! y-face high INFLOW
      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(w).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(w)
            do k = ARG_L3(w), ARG_H3(w)
               z = domnlo(3) + (k+0.5)*dx(3)
               zc = z - z_top                              
               do i = ARG_L1(w), ARG_H1(w)
                  x = domnlo(1) + (i+0.5)*dx(1)
                  xc = x - x_top                  
                  r = sqrt( xc*xc + zc*zc )
                  w(i,j,k) = 0.0d0
                  if (r.lt.r_top) then

                    pert = -2.0d0 * sin(Lf*xc) * sin(Lf*time) * cos(Lf*xc)
                    pert = 0.06d0 * pert * factor                     
                    w(i,j,k) = w(i,j,k) + pert * (1.0d0-(r/r_top)**2)
                    
                  endif
               end do
	    end do
	 end do
      end if

      ! z-face low        
      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(w).lt.domlo(3)) then
         do k = ARG_L3(w), domlo(3)-1
            do j = ARG_L2(w), ARG_H2(w)
               do i = ARG_L1(w), ARG_H1(w)
                  !                  w(i,j,k) = 0.0d0
                  w(i,j,k) = w(i,j,k+1)                  
               end do
            end do
         end do
      end if

      ! z-face high
      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(w).gt.domhi(3)) then
         do k = domhi(3)+1, ARG_H3(w)
            do j = ARG_L2(w), ARG_H2(w)
               do i = ARG_L1(w), ARG_H1(w)
                  !                  w(i,j,k) = 0.0d0
                  w(i,j,k) = w(i,j,k-1)                  
               end do
            end do
         end do
      end if

!#ifdef BL_DO_FLCT
!      if (forceInflow) deallocate(uflct)
!#endif

      end subroutine FORT_ZVELFILL

      
!---------------------------------------------------
!      REAL_T function meanPlateVel(x,y)
!        
!      implicit none
!#include <probdata.H>
!      REAL_T x,y
!      REAL_T totHoleA, holeSepX, holeSepY, Lx, Ly
!      REAL_T jetVel, xHole, yHole, rHole, rFact
!      REAL_T hRad, hBL
!      integer iHoleX, iHoleY
!
!      Lx = (domnhi(1) - domnlo(1))
!      Ly = (domnhi(2) - domnlo(2))
!      holeSepX = Lx / nHolesX
!      holeSepY = Ly / nHolesY
!      hRad = MIN(holeRad, 0.3*MIN(holeSepX,holeSepY))
!      hBL = hRad/holeBLfac
!
!      totHoleA = nHolesX*nHolesY*Pi*hRad**2
!      jetVel = adv_vel * Lx * Ly / totHoleA
!
!      iHoleX = INT(x/holeSepX)
!      iHoleY = INT(y/holeSepY)
!
!      xHole = x - (iHoleX + 0.5)*holeSepX
!      yHole = y - (iHoleY + 0.5)*holeSepY
!      rHole = SQRT(xHole**2 + yHole**2)
!
!      rFact = 1.d0 - zblend1(rHole,hRad,hBL)
!      meanPlateVel = jetVel * rFact
!      end function meanPlateVel
!
!
!---------------------------------------------------      
      REAL_T function plateVel(x,y,z,cord,flct)
      implicit none

      REAL_T   flct
      integer dir(SDIM), numholes(SDIM), idholes(SDIM)
      REAL_T  holespace(SDIM)
      REAL_T  dist(SDIM), loc(SDIM)
!c
!c     ::::: local variables
!c
      integer cord, n
      REAL_T  twicePi
      REAL_T totHoleA
      REAL_T jetVel, radius, sectarea
      REAL_T x, y, z, xcen,ycen, xloc, yloc, xcord(5), ycord(5)
      REAL_T mean, vholeSp, vane, xshift, yshift, delta,vslot
      integer icell,jcell,npts

      integer NringMAX, NperMAX
      parameter (NringMAX=30,NperMAX=300)
      integer NperRing(NringMAX)
      double precision centHole(NringMAX,NperMAX,2)
      double precision ringRad(NringMAX), cent(2)
      integer i,j, Nring
      double precision theta, ringSp


#include <probdata.H>
      mean = 0.40
      dir(1) = 0
      dir(2) = 0
      dir(3) = 0
      vholeSp = 0.d0
      vane = 0.001
      delta = 0.d0
      vslot = 32.d0


      do  npts = 1, 5
         xcord(npts) = 0
         ycord(npts) = 0
      end do

      if (adv_dir.eq.1) then
         dir(1) = 1
      else if (adv_dir.eq.2) then
         dir(2) = 1
      else
         dir(3) = 1
      endif
!c     perforated plate of square pattern
#if 0

      radius = 0.d0
      jetVel = 0.d0

      twicePi = two*Pi

      numholes(1)= nHolesX
      numholes(2)= nHolesY
      numholes(3)= nHolesZ

      totHoleA = numholes(1)*numholes(2)*numholes(3)*Pi*holeRad**2


      sectarea = dir(1)*(domnhi(2) - domnlo(2))*(domnhi(3) - domnlo(3))+ &
          dir(2)*(domnhi(1) - domnlo(1)) * (domnhi(3) - domnlo(3))+ &
          dir(3)*(domnhi(1) - domnlo(1)) * (domnhi(2) - domnlo(2))

      do n = 1, 3
         holespace(n) = (domnhi(n) - domnlo(n)) / numholes(n)
      end do

      jetVel = adv_vel*sectarea/totHoleA

      loc(3) = z
      loc(2) = y
      loc(1) = x

      do n = 1, 3
         idholes(n) = INT(loc(n)/holespace(n))
         dist(n) = loc(n) - (idholes(n) + half)*holespace(n)
      end do

      radius = SQRT(dir(1)*(dist(2)**2 + dist(3)**2)+ &
          dir(2)*(dist(1)**2 + dist(3)**2) + &
          dir(3)*(dist(1)**2 + dist(2)**2))

      if (radius .le. holeRad) then
         plateVel = jetVel*dir(cord)+flct*turb_scale
      else
         plateVel = 0.d0
      endif
#endif
!c     perforated plate of hexagnoal pattern
#if 0

      jetVel = adv_vel/(1.0 - holeBLfac)
      holeSp = SQRT(two*Pi/(sqrt(3.d0)*(1.d0 - holeBLfac)))*holeRad

      vholeSp = sqrt(3.d0)*holeSp

      plateVel = 0.d0

      icell = int (x/holeSp)
      jcell = int (y/vholeSp)

      xcen = holeSp*(dfloat(icell) + half)
      ycen = vholeSp*(dfloat(jcell) + half)

      xloc = x-xcen
      yloc = y-ycen

      xcord(1) = -half*holeSp
      ycord(1) = -half*vholeSp

      xcord(2) = -half*holeSp
      ycord(2) = half*vholeSp

      xcord(3) = half*holeSp
      ycord(3) = half*vholeSp

      xcord(4) = half*holeSp
      ycord(4) = -half*vholeSp

      xcord(5) = 0.d0
      ycord(5) = 0.d0

      do npts = 1, 5

         radius = SQRT((xloc-xcord(npts))**2 + (yloc-ycord(npts))**2)

         if(radius .le. holeRad) then
            plateVel = jetVel*dir(cord) + flct*turb_scale
            plateVel = plateVel *(1.d0 - 2.d0*zblend1(radius,holeRad,holeRad/8.0))
         endif
      end do


#endif
!c     perforated plate of hexagnoal pattern, no flow on boundaries (to enable periodic ok)
#if 0

      jetVel = adv_vel/(1.0 - holeBLfac)
      holeSp = SQRT(two*Pi/(sqrt(3.d0)*(1.d0 - holeBLfac)))*holeRad

      vholeSp = sqrt(3.d0)*holeSp

      plateVel = 0.d0

      icell = int (x/holeSp)
      jcell = int (y/vholeSp)

      xcen = holeSp*(dfloat(icell) + half)
      ycen = vholeSp*(dfloat(jcell) + half)

      xloc = x-xcen
      yloc = y-ycen

      xcord(1) = -half*holeSp
      ycord(1) = -half*vholeSp

      xcord(2) = -half*holeSp
      ycord(2) = half*vholeSp

      xcord(3) = half*holeSp
      ycord(3) = half*vholeSp

      xcord(4) = half*holeSp
      ycord(4) = -half*vholeSp

      xcord(5) = 0.d0
      ycord(5) = 0.d0

      do npts = 1, 5

         radius = SQRT((xloc-xcord(npts))**2 + (yloc-ycord(npts))**2)

         if(radius .le. holeRad) then
            plateVel = jetVel*dir(cord) + flct*turb_scale
!c            plateVel = plateVel *(1.d0 - 2.d0*zblend1(radius,holeRad,holeRad/8.0))
         endif
      end do

      if ((x.lt.domnlo(1) + .001) .or.  &
          (x.gt.domnhi(1) - .001) .or. &
          (y.lt.domnlo(2) + .001) .or.  &
          (y.gt.domnhi(2) - .001)) then
         plateVel = 0.d0
      endif

#endif

!c     perforated plate with radial pattern, no flow on boundaries (to enable periodic ok)
#if 1
      jetVel = adv_vel/(1.0 - holeBLfac)

! FIXME : from compiler
!PROB_3D.F90:7106:0:
!
!          cent(n) = 0.5d0 * (domnlo(n) + domnhi(n))
!
!Warning: iteration 2 invokes undefined behavior [-Waggressive-loop-optimizations]
!PROB_3D.F90:7105:0:
!
!       do n=1,3
!
!note: within this loop
!
      do n=1,3
         cent(n) = 0.5d0 * (domnlo(n) + domnhi(n))
      enddo
      holeRad = 0.0036d0
      ringSp = .0112d0
      Nring = MIN(NringMAX,1 + INT(0.5d0*sqrt(2.d0)*MAX(domnlo(1)+domnhi(1),domnlo(2)+domnhi(2))/ringSp))

      do i=1,Nring
         ringRad(i) = ringSp * (i-1)
         if (i.eq.1) then
            NperRing(i) = 1
            centHole(i,1,1) = cent(1)
            centHole(i,1,2) = cent(2)
         else
            NperRing(i) = 6 * (i-1)
            do j=1,NperRing(i)
               if (j .le. NperMAX) then
                  theta = j*4*ASIN(1.d0)/DBLE(NperRing(i))
                  centHole(i,j,1) = cent(1) + ringRad(i) * SIN(theta)
                  centHole(i,j,2) = cent(2) + ringRad(i) * COS(theta)
               endif
            enddo
         endif
      enddo

      plateVel = 0.d0
      do i=1,Nring
         do j=1,NperRing(i)
            if (j .le. NperMAX) then
               radius = sqrt( (x-centHole(i,j,1))**2 + (y-centHole(i,j,2))**2 )
               if(radius .le. holeRad) then
                  plateVel = jetVel*dir(cord) + flct*turb_scale
               endif
            endif
         enddo
      enddo

      if ((x.lt.domnlo(1) + .001) .or. &
          (x.gt.domnhi(1) - .001) .or. &
          (y.lt.domnlo(2) + .001) .or.  &
          (y.gt.domnhi(2) - .001)) then
         plateVel = 0.d0
      endif

#endif
!c  this is the original code
#if 0

      jetVel = adv_vel/(1.0 - holeBLfac)
      holeSp = SQRT(two*Pi/(sqrt(3.d0)*(1.d0 - holeBLfac)))*holeRad

      if(y.ge.0.d0)then
         jcell = int(y / (half*sqrt(3.d0)*holeSp))
      else
         jcell = -1-int(-y/(half*sqrt(3.d0)*holeSp))
      endif
      yloc = y-dfloat(jcell)*half*sqrt(3.d0)*holeSp
      x = x-y/sqrt(3.d0)
      if(x.ge.0.d0) then
         icell = int(x / (holeSp))
      else
         icell = -1-int(-x/(holeSp))
      endif

      xloc = x - holeSp * dfloat(icell) + yloc/sqrt(3.d0)
      xcen =.75d0*holeSp
      ycen =(sqrt(3.d0)*.25d0)*holeSp
      radius = SQRT((xloc-xcen)**2 + (yloc-ycen)**2)
      if(radius .le. holeRad) then
         plateVel = jetVel*dir(cord) + flct*turb_scale
!c         plateVel = plateVel *(1.d0 - 2.d0*zblend1(radius,holeRad,holeRad/8.0))
!c      else
!c         plateVel = -jetVel*dir(cord)
!c         plateVel = 0.d0
      endif

#endif
!c     perforated plate of hexagnoal pattern with slot  (2 gaps in swirl region)
#if 0

!c     jetVel = adv_vel/(1.0 - holeBLfac)
!c
      jetVel = 18.d0
      plateVel = 0.d0

      holeSp = SQRT(two*Pi/(sqrt(3.d0)*(1.d0 - holeBLfac)))*holeRad
      vholeSp = sqrt(3.d0)*holeSp

      icell = int (x/holeSp)
      jcell = int (y/vholeSp)

      xcen = holeSp*(dfloat(icell) + half)
      ycen = vholeSp*(dfloat(jcell) + half)

      xloc = x-xcen
      yloc = y-ycen

      xcord(1) = -half*holeSp
      ycord(1) = -half*vholeSp

      xcord(2) = -half*holeSp
      ycord(2) = half*vholeSp

      xcord(3) = half*holeSp
      ycord(3) = half*vholeSp

      xcord(4) = half*holeSp
      ycord(4) = -half*vholeSp

      xcord(5) = 0.d0
      ycord(5) = 0.d0

      do npts = 1, 5

         radius = SQRT((xloc-xcord(npts))**2 + (yloc-ycord(npts))**2)

         if(radius .le. holeRad) then
            plateVel = jetVel*dir(cord) + flct*turb_scale
!c            plateVel = plateVel *(1.d0 - 2.d0*zblend1(radius,holeRad,holeRad/8.0))
         endif
      end do
!c     in the slot region (5 mm width)
      if (y.ge.((domnhi(2)-domnlo(2)) - slotWidth)) then
         if (cord .eq. 1) then
            plateVel = 32.d0*tan(37.d0/180.d0*Pi)

            if (x .ge. 0.25*(domnhi(1)-domnlo(1)) .and. (x .le. (0.25*(domnhi(1)-domnlo(1)) + vane))) then
               plateVel = 0.d0
            else  if (x .ge. 0.75*(domnhi(1)-domnlo(1)) .and. (x .le. (0.75*(domnhi(1)-domnlo(1)) + vane))) then
               plateVel = 0.d0
            end if

         else if (cord .eq. 2) then
            plateVel = 0.d0
         else
            if (x .ge. 0.25*(domnhi(1)-domnlo(1)) .and. (x .le. (0.25*(domnhi(1)-domnlo(1)) + vane))) then
               plateVel = 0.d0
            else  if (x .ge. 0.75*(domnhi(1)-domnlo(1)) .and. (x .le. (0.75*(domnhi(1)-domnlo(1)) + vane))) then
               plateVel = 0.d0
            else
               plateVel = 32.d0
            endif
         endif
      endif
      if ((y.ge.((domnhi(2)-domnlo(2)) - slotWidth - vane)) .and. &
          (y.lt.((domnhi(2)-domnlo(2)) - slotWidth))) then
         plateVel = 0.d0
      endif
!c     in the region of 1 mm between the slot region and the plate

#endif


!c     perforated plate of hexagnoal pattern with slot and shifting in x and y directions (4 gaps in swirl region)
#if 0

!c     jetVel = adv_vel/(1.0 - holeBLfac)
!c

      jetVel = adv_vel
      vslot = slot_vel

      plateVel = 0.d0


      holeSp = SQRT(two*Pi/(sqrt(3.d0)*(1.d0 - holeBLfac)))*holeRad
      vholeSp = sqrt(3.d0)*holeSp

      xshift = alpha*holeSp
      yshift = beta*vholeSp

      icell = int ((x-xshift)/holeSp)
      jcell = int ((y-yshift)/vholeSp)

      if ((x-xshift) .le. 0.d0) then
         icell = - 1 - int ((x-xshift)/holeSp)
      endif
      if ((y-yshift) .le. 0.d0) then
         jcell = - 1 - int ((y-yshift)/vholeSp)
      endif

      xcen = holeSp*(dfloat(icell) + half) + xshift
      ycen = vholeSp*(dfloat(jcell) + half) + yshift

      xloc = x-xcen
      yloc = y-ycen

      xcord(1) = -half*holeSp
      ycord(1) = -half*vholeSp

      xcord(2) = -half*holeSp
      ycord(2) = half*vholeSp

      xcord(3) = half*holeSp
      ycord(3) = half*vholeSp

      xcord(4) = half*holeSp
      ycord(4) = -half*vholeSp

      xcord(5) = 0.d0
      ycord(5) = 0.d0

      do npts = 1, 5

         radius = SQRT((xloc-xcord(npts))**2 + (yloc-ycord(npts))**2)

         if(radius .le. holeRad) then
            plateVel = jetVel*dir(cord) + flct*turb_scale
!c            plateVel = plateVel *(1.d0 - 2.d0*zblend1(radius,holeRad,holeRad/8.0))
         endif
      end do
!c     in the slot region (5 mm width)
      if (y.ge.((domnhi(2)-domnlo(2)) - slotWidth)) then
         delta = (domnhi(1)-domnlo(1))
         if (cord .eq. 1) then
!c            plateVel = 32.d0*tan(37.d0/180.d0*Pi)
            plateVel = vslot*tan(37.d0/180.d0*Pi)

            do n = 1, 7, 2
               if (x .ge. (dfloat(n)/8.d0*delta-half*vane) .and. (x .le. (dfloat(n)/8.d0*delta + half*vane))) then
                  plateVel = 0.d0
               end if
            end do

         else if (cord .eq. 2) then
            plateVel = 0.d0
         else
!c            plateVel = 32.d0
            plateVel = vslot
            do n = 1, 7, 2
               if (x .ge. (dfloat(n)/8.d0*delta-half*vane) .and. (x .le. (dfloat(n)/8.d0*delta + half*vane))) then
                  plateVel = 0.d0
               end if
            end do
         endif
      endif
!c     in the thin 1 mm wide no flow region
      if ((y.ge.((domnhi(2)-domnlo(2)) - slotWidth - vane)) .and. &
          (y.lt.((domnhi(2)-domnlo(2)) - slotWidth))) then
         plateVel = 0.d0
      endif


#endif



      end function PlateVel



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
!c ::: u        <=  full velocity array
!c ::: DIMS(u)   => index extent of v array
!c ::: domlo,hi  => index extent of problem domain
!c ::: dx        => cell spacing
!c ::: xlo       => physical location of lower left hand
!c :::	           corner of rho array
!c ::: time      => problem evolution time
!c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi,comp)
!c ::: -----------------------------------------------------------

  subroutine FORT_VELFILL (u,DIMS(u),domlo,domhi,dx,xlo,time,bc) &
                               bind(C, name="FORT_VELFILL")
      implicit none
      integer    DIMDEC(u)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     u(DIMV(u),SDIM)
      integer    bc(SDIM,2,SDIM)
      integer    lo(SDIM),hi(SDIM)
      integer    i, j, k
      integer    isioproc
      integer    ii, ji, ip, jp
      REAL_T     x, y, z, r, jv
      REAL_T     xi, yi, xp, yp, rp
      REAL_T     uxp, uyp, up, ax, r2

      integer    ilo, ihi, jlo, jhi

#include <probdata.H>

#ifdef BL_DO_FLCT
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      integer DIMDEC(uflct)
      REAL_T  t_flct
      REAL_T, allocatable :: uflct(:,:,:)
      REAL_T, allocatable :: vflct(:,:,:)
      REAL_T, allocatable :: wflct(:,:,:)
#include <INFL_FORCE_F.H>
#endif

      lo(1) = ARG_L1(u)
      lo(2) = ARG_L2(u)
      lo(3) = ARG_L3(u)
      hi(1) = ARG_H1(u)
      hi(2) = ARG_H2(u)
      hi(3) = ARG_H3(u)


      !shouldnt be in here
      print*, " "
      print*, " >>> STOP STOP STOP STOP STOP STOP. <<<" 
      print*, " >>> WARNING: IN VELFILL, NOT FIXED <<<"
      print*, " >>> STOP STOP STOP STOP STOP STOP. <<<"
      print*, " "      

    
  end subroutine FORT_VELFILL

  
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

      
      ! wtf? this is incompressible, divu should be zero...

!      if (adv_dir .eq. 3) then
!         z_vel = adv_vel
!      else
!         z_vel = zero
!      end if

      call filcc(divu,DIMS(divu),domlo,domhi,dx,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(divu).lt.domlo(1)) then
         do i = ARG_L1(divu), domlo(1)-1
            do k = ARG_L3(divu), ARG_H3(divu)
               do j = ARG_L2(divu), ARG_H2(divu)
                  divu(i,j,k) = 0.0d0
               end do
	    end do
	 end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(divu).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(divu)
            do k = ARG_L3(divu), ARG_H3(divu)
               do j = ARG_L2(divu), ARG_H2(divu)
                  divu(i,j,k) = 0.0d0
               end do
	    end do
	 end do
      end if

      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(divu).lt.domlo(2)) then
         do j = ARG_L2(divu), domlo(2)-1
            do k = ARG_L3(divu), ARG_H3(divu)
               do i = ARG_L1(divu), ARG_H1(divu)
                  divu(i,j,k) = 0.0d0
               end do
            end do
	 end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(divu).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(divu)
            do k = ARG_L3(divu), ARG_H3(divu)
               do i = ARG_L1(divu), ARG_H1(divu)
                  divu(i,j,k) = 0.0d0
               end do
            end do
	 end do
      end if

      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(divu).lt.domlo(3)) then
         do k = ARG_L3(divu), domlo(3)-1
            do j = ARG_L2(divu), ARG_H2(divu)
               do i = ARG_L1(divu), ARG_H1(divu)
                  divu(i,j,k) = 0.0d0
               end do
            end do
         end do
      end if

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(divu).gt.domhi(3)) then
         do k = domhi(3)+1, ARG_H3(divu)
            do j = ARG_L2(divu), ARG_H2(divu)
               do i = ARG_L1(divu), ARG_H1(divu)
                  divu(i,j,k) = 0.0d0
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


      ! what is dsdt?

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

!      fix_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .ne. INT_DIR)
!      per_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .eq. INT_DIR)
!      fix_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .ne. INT_DIR)
!      per_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .eq. INT_DIR)
!      fix_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .ne. INT_DIR)
!      per_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .eq. INT_DIR)
!      fix_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .ne. INT_DIR)
!      per_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .eq. INT_DIR)
!      fix_zlo = (ARG_L3(p) .lt. domlo(3)) .and. (bc(3,1) .ne. INT_DIR)
!      per_zlo = (ARG_L3(p) .lt. domlo(3)) .and. (bc(3,1) .eq. INT_DIR)
!      fix_zhi = (ARG_H3(p) .gt. domhi(3)) .and. (bc(3,2) .ne. INT_DIR)
!      per_zhi = (ARG_H3(p) .gt. domhi(3)) .and. (bc(3,2) .eq. INT_DIR)

      ilo = max(ARG_L1(p),domlo(1))
      ihi = min(ARG_H1(p),domhi(1))
      jlo = max(ARG_L2(p),domlo(2))
      jhi = min(ARG_H2(p),domhi(2))
      Klo = max(ARG_L3(p),domlo(3))
      khi = min(ARG_H3(p),domhi(3))


      ! just assume dpdn = zero at in/out flow
      if (bc(1,1).eq.EXT_DIR.and.ARG_L1(p).lt.domlo(1)) then
         do i = ARG_L1(p), domlo(1)-1
            do k = ARG_L3(p), ARG_H3(p)
               do j = ARG_L2(p), ARG_H2(p)
                  p(i,j,k) = p(domlo(1),j,k)
               end do
	    end do
	 end do
      end if

      if (bc(1,2).eq.EXT_DIR.and.ARG_H1(p).gt.domhi(1)) then
         do i = domhi(1)+1, ARG_H1(p)
            do k = ARG_L3(p), ARG_H3(p)
               do j = ARG_L2(p), ARG_H2(p)
                  p(i,j,k) = p(domhi(1),j,k)
               end do
	    end do
	 end do
      end if

      
      if (bc(2,1).eq.EXT_DIR.and.ARG_L2(p).lt.domlo(2)) then
         do j = ARG_L2(p), domlo(2)-1
            do k = ARG_L3(p), ARG_H3(p)
               do i = ARG_L1(p), ARG_H1(p)
                  p(i,j,k) = p(i,domlo(2),k)
               end do
            end do
	 end do
      end if

      if (bc(2,2).eq.EXT_DIR.and.ARG_H2(p).gt.domhi(2)) then
         do j = domhi(2)+1, ARG_H2(p)
            do k = ARG_L3(p), ARG_H3(p)
               do i = ARG_L1(p), ARG_H1(p)
                  p(i,j,k) = p(i,domhi(2),k)
               end do
            end do
	 end do
      end if

      
      if (bc(3,1).eq.EXT_DIR.and.ARG_L3(p).lt.domlo(3)) then
         do k = ARG_L3(p), domlo(3)-1
            do j = ARG_L2(p), ARG_H2(p)
               do i = ARG_L1(p), ARG_H1(p)
                  p(i,j,k) = p(i,j,domlo(3))
               end do
            end do
         end do
      end if

      if (bc(3,2).eq.EXT_DIR.and.ARG_H3(p).gt.domhi(3)) then
         do k = domhi(3)+1, ARG_H3(p)
            do j = ARG_L2(p), ARG_H2(p)
               do i = ARG_L1(p), ARG_H1(p)
                  p(i,j,k) = p(i,j,domhi(3))
               end do
            end do
         end do
      end if

    end subroutine FORT_PRESFILL
    

!***************************************************************
!    "Minimal" random number generator of Park and Miller with
!     Bays-Durham shuffle and added safeguards.  Returns a
!     uniform random deviate between 0.0 and 1.0 (exclusive of
!     the endpoint values).  Call with IDUM a negative integer
!     to initialize; thereafter, do not alter IDUM between
!     successive deviates in a sequence.  RNMX should approximate
!       the largest floating value that is less than 1.
!***************************************************************      
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
    
!***************************************************************
!*     After Press et al., Numerical Recipes for Fortran
!***************************************************************
   

end module prob_3D_module
