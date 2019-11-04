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

  use amrex_error_module, only : amrex_abort

  implicit none

  private

  public :: derdvrho, dermprho, deravgpres, &
            dermgvort, dermgdivu, &
            dergrdpx, dergrdpy, dergrdpz

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
!  Compute the amagnitude of the vorticity from the 
!  velocity field
!=========================================================

   subroutine dermgvort (e,   e_lo, e_hi, nv, &
                         dat, d_lo, d_hi, ncomp, &
                         lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                         level, grid_no) &
                         bind(C, name="dermgvort")

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(inout), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      REAL_T :: uy, uz, vx, vz, wx, wy, dx, dy, dz
      REAL_T :: uycen, uzcen, uylo, uyhi, uzlo, uzhi
      REAL_T :: vxcen, vzcen, vxlo, vxhi, vzlo, vzhi
      REAL_T :: wxcen, wycen, wxlo, wxhi, wylo, wyhi
      REAL_T :: vorfun

      logical :: fixvlo_x, fixwlo_x, fixvhi_x, fixwhi_x
      logical :: fixulo_y, fixwlo_y, fixuhi_y, fixwhi_y
      logical :: fixulo_z, fixvlo_z, fixuhi_z, fixvhi_z

      integer :: i, j, k

!
!     ::::: some useful macro definitions
!
#     define U(i,j,k) dat(i,j,k,1)
#     define V(i,j,k) dat(i,j,k,2)
#     define W(i,j,k) dat(i,j,k,3)

#     define ULOY bc(2,1,1)
#     define UHIY bc(2,2,1)
#     define ULOZ bc(3,1,1)
#     define UHIZ bc(3,2,1)

#     define VLOX bc(1,1,2)
#     define VHIX bc(1,2,2)
#     define VLOZ bc(3,1,2)
#     define VHIZ bc(3,2,2)

#     define WLOX bc(1,1,3)
#     define WHIX bc(1,2,3)
#     define WLOY bc(2,1,3)
#     define WHIY bc(2,2,3)

!
!     ::::: statement functions that implement stencil
!
      uycen(i,j,k) = half*(U(i,j+1,k)-U(i,j-1,k))/dy
      uylo(i,j,k)  = (U(i,j+1,k)+three*U(i,j,k)-four*U(i,j-1,k))/(three*dy)
      uyhi(i,j,k)  =-(U(i,j-1,k)+three*U(i,j,k)-four*U(i,j+1,k))/(three*dy)

      vxcen(i,j,k) = half*(V(i+1,j,k)-V(i-1,j,k))/dx
      vxlo(i,j,k)  = (V(i+1,j,k)+three*V(i,j,k)-four*V(i-1,j,k))/(three*dx)
      vxhi(i,j,k)  =-(V(i-1,j,k)+three*V(i,j,k)-four*V(i+1,j,k))/(three*dx)

#if ( AMREX_SPACEDIM == 3 )
      uzcen(i,j,k) = half*(U(i,j,k+1)-U(i,j,k-1))/dz
      uzlo(i,j,k)  = (U(i,j,k+1)+three*U(i,j,k)-four*U(i,j,k-1))/(three*dz)
      uzhi(i,j,k)  =-(U(i,j,k-1)+three*U(i,j,k)-four*U(i,j,k+1))/(three*dz)

      vzcen(i,j,k) = half*(V(i,j,k+1)-V(i,j,k-1))/dz
      vzlo(i,j,k)  = (V(i,j,k+1)+three*V(i,j,k)-four*V(i,j,k-1))/(three*dz)
      vzhi(i,j,k)  =-(V(i,j,k-1)+three*V(i,j,k)-four*V(i,j,k+1))/(three*dz)

      wxcen(i,j,k) = half*(W(i+1,j,k)-W(i-1,j,k))/dx
      wxlo(i,j,k)  = (W(i+1,j,k)+three*W(i,j,k)-four*W(i-1,j,k))/(three*dx)
      wxhi(i,j,k)  =-(W(i-1,j,k)+three*W(i,j,k)-four*W(i+1,j,k))/(three*dx)

      wycen(i,j,k) = half*(W(i,j+1,k)-W(i,j-1,k))/dy
      wylo(i,j,k)  = (W(i,j+1,k)+three*W(i,j,k)-four*W(i,j-1,k))/(three*dy)
      wyhi(i,j,k)  =-(W(i,j-1,k)+three*W(i,j,k)-four*W(i,j+1,k))/(three*dy)
#endif

#if ( AMREX_SPACEDIM == 2 )
      vorfun(uy,uz,vx,vz,wx,wy) = vx - uy
#elif ( AMREX_SPACEDIM == 3 )
      vorfun(uy,uz,vx,vz,wx,wy) = sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)
#endif

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      ! Init all logical tests on BC to false
      fixvlo_x = .FALSE. ; fixwlo_x = .FALSE. ; fixvhi_x = .FALSE. ; fixwhi_x = .FALSE.
      fixulo_y = .FALSE. ; fixwlo_y = .FALSE. ; fixuhi_y = .FALSE. ; fixwhi_y = .FALSE.
      fixulo_z = .FALSE. ; fixvlo_z = .FALSE. ; fixuhi_z = .FALSE. ; fixvhi_z = .FALSE.

      ! Init all vorticity comp. In 2d uz, vz, wx, wy will alway be zero
      uy = 0.0d0
      uz = 0.0d0
      vx = 0.0d0
      vz = 0.0d0
      wx = 0.0d0
      wy = 0.0d0

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               uy = uycen(i,j,k)
               vx = vxcen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wx = wxcen(i,j,k)
               wy = wycen(i,j,k)
#endif
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end do

      fixvlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (VLOX .eq. EXT_DIR .or. VLOX .eq. HOEXTRAP) )
      fixvhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (VHIX .eq. EXT_DIR .or. VHIX .eq. HOEXTRAP) )
#if ( AMREX_SPACEDIM == 3 )
      fixwlo_x = ( (lo(1) .eq. domlo(1)) .and. &
                  (WLOX .eq. EXT_DIR .or. WLOX .eq. HOEXTRAP) )
      fixwhi_x = ( (hi(1) .eq. domhi(1)) .and. &
                  (WHIX .eq. EXT_DIR .or. WHIX .eq. HOEXTRAP) )
#endif

      fixulo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (ULOY .eq. EXT_DIR .or. ULOY .eq. HOEXTRAP) )
      fixuhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (UHIY .eq. EXT_DIR .or. UHIY .eq. HOEXTRAP) )
#if ( AMREX_SPACEDIM == 3 )
      fixwlo_y = ( (lo(2) .eq. domlo(2)) .and. &
                  (WLOY .eq. EXT_DIR .or. WLOY .eq. HOEXTRAP) )
      fixwhi_y = ( (hi(2) .eq. domhi(2)) .and. &
                  (WHIY .eq. EXT_DIR .or. WHIY .eq. HOEXTRAP) )

      fixulo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (ULOZ .eq. EXT_DIR .or. ULOZ .eq. HOEXTRAP) )
      fixuhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (UHIZ .eq. EXT_DIR .or. UHIZ .eq. HOEXTRAP) )
      fixvlo_z = ( (lo(3) .eq. domlo(3)) .and. &
                  (VLOZ .eq. EXT_DIR .or. VLOZ .eq. HOEXTRAP) )
      fixvhi_z = ( (hi(3) .eq. domhi(3)) .and. &
                  (VHIZ .eq. EXT_DIR .or. VHIZ .eq. HOEXTRAP) )
#endif

!
!     First do all the faces
!
      if (fixvlo_x .or. fixwlo_x) then
         i = lo(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
               uy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
               wy = wycen(i,j,k)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
#endif
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixvhi_x .or. fixwhi_x) then
         i = hi(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
               uy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
               wy = wycen(i,j,k)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
#endif
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixulo_y .or. fixwlo_y) then
         j = lo(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
               wx = wxcen(i,j,k)
               wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
#endif
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixuhi_y .or. fixwhi_y) then
         j = hi(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
               wx = wxcen(i,j,k)
               wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
#endif
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

#if ( AMREX_SPACEDIM == 3 )
      if (fixulo_z .or. fixvlo_z) then
         k = lo(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = uycen(i,j,k)
               wy = wycen(i,j,k)
               uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
               vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixuhi_z .or. fixvhi_z) then
         k = hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = uycen(i,j,k)
               wy = wycen(i,j,k)
               uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
               vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
               e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if
#endif

!
!     Next do all the edges
!
      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y)) then
         i = lo(1)
         j = lo(2)
         do k = lo(3),hi(3)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y)) then
         i = hi(1)
         j = lo(2)
         do k = lo(3),hi(3)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y)) then
         i = lo(1)
         j = hi(2)
         do k = lo(3),hi(3)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y)) then
         i = hi(1)
         j = hi(2)
         do k = lo(3),hi(3)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         k = lo(3)
         do j = lo(2),hi(2)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            uy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = wycen(i,j,k)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         k = lo(3)
         do j = lo(2),hi(2)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            uy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = wycen(i,j,k)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         k = hi(3)
         do j = lo(2),hi(2)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            uy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            wy = wycen(i,j,k)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         k = hi(3)
         do j = lo(2),hi(2)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            uy = uycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            wy = wycen(i,j,k)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixulo_y .or. fixwlo_y) .and. (fixulo_z .or. fixvlo_z)) then
         j = lo(2)
         k = lo(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = wxcen(i,j,k)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixuhi_y .or. fixwhi_y) .and. (fixulo_z .or. fixvlo_z)) then
         j = hi(2)
         k = lo(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = wxcen(i,j,k)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixulo_y .or. fixwlo_y) .and. (fixuhi_z .or. fixvhi_z)) then
         j = lo(2)
         k = hi(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = wxcen(i,j,k)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixuhi_y .or. fixwhi_y) .and. (fixuhi_z .or. fixvhi_z)) then
         j = hi(2)
         k = hi(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
            wx = wxcen(i,j,k)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
            e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

!
!     Finally do all the corners
!
      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         j = lo(2)
         k = lo(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         j = lo(2)
         k = lo(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         j = hi(2)
         k = lo(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         j = hi(2)
         k = lo(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         j = lo(2)
         k = hi(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         j = lo(2)
         k = hi(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         j = hi(2)
         k = hi(3)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y) .and. &
          (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         j = hi(2)
         k = hi(3)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
#if ( AMREX_SPACEDIM == 3 )
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
#endif
         e(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

#     undef U
#     undef V
#     undef W
#     undef ULOY
#     undef UHIY
#     undef ULOZ
#     undef UHIZ
#     undef VLOX
#     undef VHIX
#     undef VLOZ
#     undef VHIZ
#     undef WLOX
#     undef WHIX
#     undef WLOY
#     undef WHIY

   end subroutine dermgvort

!=========================================================
!  Compute the magnitude of the velocity divergence
!=========================================================

   subroutine dermgdivu (e,   e_lo, e_hi, nv, &
                         dat, d_lo, d_hi, ncomp, &
                         lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                         level, grid_no) &
                         bind(C, name="dermgdivu")

#if ( AMREX_SPACEDIM == 2 )                      
      use prob_2D_module, only : FORT_XVELFILL, FORT_YVELFILL
#elif ( AMREX_SPACEDIM == 3 )
      use prob_3D_module, only : FORT_XVELFILL, FORT_YVELFILL, FORT_ZVELFILL
#endif

      implicit none

!  In/Out
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: e_lo(3), e_hi(3), nv
      integer, intent(in) :: d_lo(3), d_hi(3), ncomp
      integer, intent(in) :: domlo(3), domhi(3)
      integer, intent(in) :: bc(3,2,ncomp)
      REAL_T, intent(in)  :: delta(3), xlo(3), time, dt
      REAL_T, intent(out),dimension(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nv) :: e
      REAL_T, intent(inout), dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp) :: dat
      integer, intent(in) :: level, grid_no

!  Local
      REAL_T  :: ux, vy, wz, dx, dy, dz
      REAL_T  :: uxcen, uxlo, uxhi
      REAL_T  :: vycen, vylo, vyhi
      REAL_T  :: wzcen, wzlo, wzhi
      integer :: bc_dumm(2,2,1)
      integer :: i, j, k

!
!     ::::: some useful macro definitions
!
#     define U(i,j,k) dat(i,j,k,1)
#     define V(i,j,k) dat(i,j,k,2)
#     define W(i,j,k) dat(i,j,k,3)

#     define ULOX bc(1,1,1)
#     define UHIX bc(1,2,1)
#     define VLOY bc(2,1,2)
#     define VHIY bc(2,2,2)
#     define WLOZ bc(3,1,3)
#     define WHIZ bc(3,2,3)

!
!     ::::: statement functions that implement stencil
!
      uxcen(i,j,k) = half*(U(i+1,j,k)-U(i-1,j,k))/dx
      uxlo(i,j,k) = (eight*U(i,j,k)-six*U(i+1,j,k)+U(i+2,j,k))/(three*dx)
      uxhi(i,j,k) = (eight*U(i,j,k)-six*U(i-1,j,k)+U(i-2,j,k))/(three*dx)

#if ( AMREX_SPACEDIM >= 2 )
      vycen(i,j,k) = half*(V(i,j+1,k)-V(i,j-1,k))/dy
      vylo(i,j,k) = (eight*V(i,j,k)-six*V(i,j+1,k)+V(i,j+2,k))/(three*dy)
      vyhi(i,j,k) = (eight*V(i,j,k)-six*V(i,j-1,k)+V(i,j-2,k))/(three*dy)

#if ( AMREX_SPACEDIM == 3 )
      wzcen(i,j,k) = half*(W(i,j,k+1)-W(i,j,k-1))/dz
      wzlo(i,j,k) = (eight*W(i,j,k)-six*W(i,j,k+1)+W(i,j,k+2))/(three*dz)
      wzhi(i,j,k) = (eight*W(i,j,k)-six*W(i,j,k-1)+W(i,j,k-2))/(three*dz)
#endif
#endif

#if ( AMREX_SPACEDIM == 2 )
      bc_dumm(1:2,1:2,1) = bc(1:2,1:2,1)
      call FORT_XVELFILL (dat(:,:,d_lo(3),1), d_lo(1), d_lo(2), d_hi(1), d_hi(2), & 
                          domlo, domhi, delta, xlo, time, bc_dumm(1,1,1))
      bc_dumm(1:2,1:2,1) = bc(1:2,1:2,2)
      call FORT_YVELFILL (dat(:,:,d_lo(3),2), d_lo(1), d_lo(2), d_hi(1), d_hi(2), &
                          domlo, domhi, delta, xlo, time, bc_dumm(1,1,1))
#elif ( AMREX_SPACEDIM == 3 )
      call FORT_XVELFILL (dat(:,:,:,1), d_lo(1), d_lo(2), d_lo(3), d_hi(1), d_hi(2), d_hi(3), & 
                          domlo, domhi, delta, xlo, time, bc(1,1,1))
      call FORT_YVELFILL (dat(:,:,:,2), d_lo(1), d_lo(2), d_lo(3), d_hi(1), d_hi(2), d_hi(3), &
                          domlo, domhi, delta, xlo, time, bc(1,1,2))
      call FORT_ZVELFILL (dat(:,:,:,3), d_lo(1), d_lo(2), d_lo(3), d_hi(1), d_hi(2), d_hi(3), &
                          domlo, domhi, delta, xlo, time, bc(1,1,3))
#endif

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

!
!     :: at physical bndries where an edge value is prescribed,
!     :: set the value in the outside cell so that a central
!     :: difference formula is equivalent to the higher order
!     :: one sided formula
!
      if (lo(1) == domlo(1)) then
         i = lo(1)
         if (ULOX==EXT_DIR) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i-1,j,k) = two*U(i-1,j,k) - U(i,j,k)
               end do
            end do
         else if (ULOX==HOEXTRAP) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i-1,j,k) = uxlo(i,j,k)
               end do
            end do
         end if
      end if
      if (hi(1) == domhi(1)) then
         i = hi(1)
         if (UHIX==EXT_DIR) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i+1,j,k) = two*U(i+1,j,k) - U(i,j,k)
               end do
            end do
         else if (UHIX==HOEXTRAP) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i+1,j,k) = uxhi(i,j,k)
               end do
            end do
         end if
      end if

#if ( AMREX_SPACEDIM >= 2 )
      if (lo(2) == domlo(2)) then
         j = lo(2)
         if (VLOY==EXT_DIR) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j-1,k) = two*V(i,j-1,k) - V(i,j,k)
               end do
            end do
         else if (VLOY==HOEXTRAP) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j-1,k) = vylo(i,j,k)
               end do
            end do
         end if
      end if
      if (hi(2) == domhi(2)) then
         j = hi(2)
         if (VHIY==EXT_DIR) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j+1,k) = two*V(i,j+1,k) - V(i,j,k)
               end do
            end do
         else if (VHIY==HOEXTRAP) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j+1,k) = vyhi(i,j,k)
               end do
            end do
         end if
      end if

#if ( AMREX_SPACEDIM == 3 )
      if (lo(3) == domlo(3)) then
         k = lo(3)
         if (WLOZ==EXT_DIR) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k-1) = two*W(i,j,k-1) - W(i,j,k)
               end do
            end do
         else if (WLOZ==HOEXTRAP) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k-1) = wzlo(i,j,k)
               end do
            end do
         end if
      end if
      if (hi(3) == domhi(3)) then
         k = hi(3)
         if (WHIZ==EXT_DIR) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k+1) = two*W(i,j,k+1) - W(i,j,k)
               end do
            end do
         else if (WHIZ==HOEXTRAP) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k+1) = wzhi(i,j,k)
               end do
            end do
         end if
      end if
#endif
#endif

      ux = 0.0d0
      vy = 0.0d0
      wz = 0.0d0
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ux = uxcen(i,j,k)
#if ( AMREX_SPACEDIM >= 2 )
               vy = vycen(i,j,k)
#if ( AMREX_SPACEDIM == 3 )
               wz = wzcen(i,j,k)
#endif
#endif
               e(i,j,k,1) = ux + vy + wz
            end do
         end do
      end do

!
! we overwrote the ghost cells above, so set them back below
!
#if ( AMREX_SPACEDIM == 2 )
      bc_dumm(1:2,1:2,1) = bc(1:2,1:2,1)
      call FORT_XVELFILL (dat(:,:,d_lo(3),1), d_lo(1), d_lo(2), d_hi(1), d_hi(2), & 
                          domlo, domhi, delta, xlo, time, bc_dumm(1,1,1))
      bc_dumm(1:2,1:2,1) = bc(1:2,1:2,2)
      call FORT_YVELFILL (dat(:,:,d_lo(3),2), d_lo(1), d_lo(2), d_hi(1), d_hi(2), &
                          domlo, domhi, delta, xlo, time, bc_dumm(1,1,1))
#elif ( AMREX_SPACEDIM == 3 )
      call FORT_XVELFILL (dat(:,:,:,1), d_lo(1), d_lo(2), d_lo(3), d_hi(1), d_hi(2), d_hi(3), & 
                          domlo, domhi, delta, xlo, time, bc(1,1,1))
      call FORT_YVELFILL (dat(:,:,:,2), d_lo(1), d_lo(2), d_lo(3), d_hi(1), d_hi(2), d_hi(3), &
                          domlo, domhi, delta, xlo, time, bc(1,1,2))
      call FORT_ZVELFILL (dat(:,:,:,3), d_lo(1), d_lo(2), d_lo(3), d_hi(1), d_hi(2), d_hi(3), &
                          domlo, domhi, delta, xlo, time, bc(1,1,3))
#endif

#     undef U
#     undef V      
#     undef W
#     undef ULOX
#     undef UHIX
#     undef VLOY
#     undef VHIY
#     undef WLOZ
#     undef WHIZ

   end subroutine dermgdivu

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
      integer :: i, j, k

      factor = 0.5d0
#if (AMREX_SPACEDIM >= 2 )
      factor = 0.25d0
#if (AMREX_SPACEDIM == 3 )
      factor = 0.125d0
#endif
#endif

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) =  factor * (  dat(i+1,j,k,1)     + dat(i,j,k,1)  &
#if (AMREX_SPACEDIM >= 2 )
                                       + dat(i+1,j+1,k,1)   + dat(i,j+1,k,1) &
#if (AMREX_SPACEDIM == 3 )
                                       + dat(i+1,j,k+1,1)   + dat(i,j,k+1,1)  &
                                       + dat(i+1,j+1,k+1,1) + dat(i,j+1,k+1,1) &
#endif
#endif
                                      )    
            end do
         end do
      end do

   end subroutine deravgpres

!=========================================================
!  Compute node centered pressure gradient in direction dir
!=========================================================

   subroutine gradp_dir (p, p_lo, p_hi, &
                         gp, g_lo, g_hi, &
                         lo, hi, dir, dx) &
                         bind(C, name="gradp_dir")

      implicit none

! In/Out
      integer :: lo(3),  hi(3)
      integer :: p_lo(3), p_hi(3)
      integer :: g_lo(3), g_hi(3)
      integer :: dir
      REAL_T  :: dx
      REAL_T, dimension(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3)) :: p
      REAL_T, dimension(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3)) :: gp

! Local
      integer :: i, j, k
      REAL_T  :: d

      d = one/dx
#if (AMREX_SPACEDIM >= 2 )
      d = half/dx
#if (AMREX_SPACEDIM == 3 )
      d = fourth/dx
#endif
#endif

!
!     ::::: compute gradient on interior
!
      if (dir == 0) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
#if (AMREX_SPACEDIM == 1 )
                  gp(i,j,k) = d*(p(i+1,j,k)-p(i,j,k))
#elif (AMREX_SPACEDIM == 2 )
                  gp(i,j,k) = d*(p(i+1,j,k)-p(i,j,k)+p(i+1,j+1,k)-p(i,j+1,k))
#elif (AMREX_SPACEDIM == 3 )
                  gp(i,j,k) = d*( p(i+1,j,k  )-p(i,j,k  )+p(i+1,j+1,k  )-p(i,j+1,k  ) + &
                                  p(i+1,j,k+1)-p(i,j,k+1)+p(i+1,j+1,k+1)-p(i,j+1,k+1))
#endif 

               end do
            end do
         end do
      else if (dir == 1) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
#if (AMREX_SPACEDIM == 2 )
                  gp(i,j,k) = d*(p(i,j+1,k)-p(i,j,k)+p(i+1,j+1,k)-p(i+1,j,k))
#elif (AMREX_SPACEDIM == 3 )
                  gp(i,j,k) = d*( p(i,j+1,k  )-p(i,j,k  )+p(i+1,j+1,k  )-p(i+1,j,k  ) + &
                                  p(i,j+1,k+1)-p(i,j,k+1)+p(i+1,j+1,k+1)-p(i+1,j,k+1))
#endif
               end do
            end do
         end do
      else if (dir == 2) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  gp(i,j,k) = d*( p(i,  j,k+1)-p(i,  j,k)+p(i,  j+1,k+1)-p(i,  j+1,k) + &
                                  p(i+1,j,k+1)-p(i+1,j,k)+p(i+1,j+1,k+1)-p(i+1,j+1,k))
               end do
            end do
         end do
      else
         call amrex_abort("gradp_dir: invalid dir = ")
      end if

   end subroutine gradp_dir

!=========================================================
!  Compute node centered pressure gradient in X-dir
!=========================================================

   subroutine dergrdpx (e,   e_lo, e_hi, nv, &
                        dat, d_lo, d_hi, ncomp, &
                        lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                        level, grid_no) &
                        bind(C, name="dergrdpx")

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
      integer :: i, j, k

      call gradp_dir ( dat, d_lo, d_hi, &
                       e,   e_lo, e_hi, &
                       lo, hi, 0, delta(1))

   end subroutine dergrdpx

!=========================================================
!  Compute node centered pressure gradient in Y-dir
!=========================================================

   subroutine dergrdpy (e,   e_lo, e_hi, nv, &
                        dat, d_lo, d_hi, ncomp, &
                        lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                        level, grid_no) &
                        bind(C, name="dergrdpy")

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
      integer :: i, j, k

#if (AMREX_SPACEDIM < 2 )
      call amrex_abort("dergrdpy called but AMREX_SPACEDIM<2 !")
#endif

      call gradp_dir ( dat, d_lo, d_hi, &
                       e,   e_lo, e_hi, &
                       lo, hi, 1, delta(2))

   end subroutine dergrdpy

!=========================================================
!  Compute node centered pressure gradient in Z-dir
!=========================================================

   subroutine dergrdpz (e,   e_lo, e_hi, nv, &
                        dat, d_lo, d_hi, ncomp, &
                        lo, hi, domlo, domhi, delta, xlo, time, dt, bc, &
                        level, grid_no) &
                        bind(C, name="dergrdpz")

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
      integer :: i, j, k

#if (AMREX_SPACEDIM < 3 )
      call amrex_abort("dergrdpz called but AMREX_SPACEDIM<3 !")
#endif

      call gradp_dir ( dat, d_lo, d_hi, &
                       e,   e_lo, e_hi, &
                       lo, hi, 2, delta(3))

   end subroutine dergrdpz

end module derive_nd_module
