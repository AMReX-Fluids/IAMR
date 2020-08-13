
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

module MakeForce_2D_module

  implicit none

  private

  public :: FORT_MAKEFORCE

contains

! ::: -----------------------------------------------------------
!
!     This is a dummy routine to add the forcing terms to the momentum equation
!     It should be overwritten in the local run folder.  
!
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

end module MakeForce_2D_module
