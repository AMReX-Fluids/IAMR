#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <DERIVE_F.H>
#include <PROB_NS_F.H>
#include <AMReX_ArrayLim.H>

module derive_nd_module
  
  implicit none

  private

  public :: derdvrho, dermprho

contains

!===============================================================
! This file contains functions which compute derived quantities.  
! All of the argument lists have the same template, shown below
! 
! INPUTS/OUTPUTS:
! 
! e         <= the quantity derived
! DIMS(e)   => index extent of e array
! nv        => number of components in e array (should be 1)
! dat       => data neded to derive e
! DIMS(dat) => index limits of dat array
! ncomp     => number of components of dat array (3)
! lo,hi     => subrange of e array where result is requested
! domlo,hi  => index extent of problem domain (cell centered)
! delta     => cell spacing
! xlo       => physical location of lower left hand
!         corner of e array
! time      => problem evolution time
! bc        => array of bndry types for component values
!              valid only if component touches bndry
!===============================================================

!=========================================================
!  Compute C / rho
!=========================================================

   subroutine derdvrho (e,   e_lo, e_hi, nv, &
                        dat, d_lo, d_hi, ncomp, &
                        lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                        level, grid_no) &
                        bind(C, name="derdvrho")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      integer :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = dat(i,j,k,2)/dat(i,j,k,1)
            end do
         end do
      end do

   end subroutine derdvrho

!=========================================================
!  Compute rho * C 
!=========================================================

   subroutine dermprho (e,   e_lo, e_hi, nv, &
                        dat, d_lo, d_hi, ncomp, &
                        lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                        level, grid_no) &
                        bind(C, name="dermprho")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      integer :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = dat(i,j,k,2)*dat(i,j,k,1)
            end do
         end do
      end do

   end subroutine dermprho

!=========================================================
!  Compute cell-centered pressure as average of the 
!  surrounding nodal values 
!=========================================================

   subroutine deravgpres (e,   e_lo, e_hi, nv, &
                          dat, d_lo, d_hi, ncomp, &
                          lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                          level, grid_no) &
                          bind(C, name="deravgpres")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(in), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      REAL_T  :: factor
      REAL_T  :: zfac
      integer :: i, j, k

      zfac = 0.0d0
      factor = 0.25d0
      if ( Delta(3) > 0.0d0 ) then
         factor = 0.125d0
         zfac = 1.0d0
      end if

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) =  factor * (  dat(i+1,j,k,1)     + dat(i,j,k,1)  &
                                       + dat(i+1,j+1,k,1)   + dat(i,j+1,k,1) &
                                       + zfac * (  dat(i+1,j,k+1,1)   + dat(i,j,k+1,1)  &
                                                 + dat(i+1,j+1,k+1,1) + dat(i,j+1,k+1,1) ) )
            end do
         end do
      end do

   end subroutine deravgpres
end module derive_nd_module
