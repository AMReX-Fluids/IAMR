c
c $Id: DERIVE_3D.F,v 1.7 1998-03-31 17:26:34 lijewski Exp $
c
#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "DERIVE_F.H"
#include "PROB_F.H"
#include "ArrayLim.H"

#define SDIM 3

c     -----------------------------------------------------------
c     This file contains functions which compute derived quantities.  
c     All of the argument lists have the same template, shown below
c     
c     INPUTS/OUTPUTS:
c     
c     e         <= the quantity derived
c     DIMS(e)   => index extent of e array
c     nv        => number of components in e array (should be 1)
c     dat       => data neded to derive e
c     DIMS(dat) => index limits of dat array
c     ncomp     => number of components of dat array (3)
c     lo,hi     => subrange of e array where result is requested
c     domlo,hi  => index extent of problem domain (cell centered)
c     delta     => cell spacing
c     xlo       => physical location of lower left hand
c 	           corner of e array
c     time      => problem evolution time
c     bc        => array of bndry types for component values
c                  valid only if component touches bndry
c     -----------------------------------------------------------

      subroutine FORT_DERKENG (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                         lo,hi,domlo,domhi,delta,xlo,time,dt,
     &                         bc,level,grid_no)
c
c     This routine will derive kinetic energy from density
c     and the velocity field.
c
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      integer    nxlo, nxhi, nylo, nyhi, nzlo, nzhi
      REAL_T     rho, u, v, w
c      
c     ::::: lets punt if not in domain interior
c
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
         write(6,*) "FORT_DERKENG: outside domain"
	 stop
      end if

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rho = dat(i,j,k,1)
               u   = dat(i,j,k,2)
               v   = dat(i,j,k,3)
               w   = dat(i,j,k,4)
               e(i,j,k,1) = half*rho*(u**2 + v**2 + w**2)
            end do
         end do
      end do

      end

      subroutine FORT_DERLOGS (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                         lo,hi,domlo,domhi,delta,xlo,time,dt,
     &                         bc,level, grid_no)
c
c     This routine will derive log of given scalar quantity
c
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      integer    nxlo, nxhi, nylo, nyhi, nzlo, nzhi
      REAL_T     rho
      REAL_T     sml

      data sml / 1.0e-10 /
c      
c     ::::: lets punt if not in domain interior
c
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
         write(6,*) "FORT_DERLOGS: outside domain"
	 stop
      end if

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rho = max(dat(i,j,k,1),sml)
               e(i,j,k,1) = log10(rho)
            end do
         end do
      end do

      end

      subroutine FORT_DERMVEL (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                         lo,hi,domlo,domhi,delta,xlo,time,dt,
     &                         bc,level, grid_no)
c
c ::: This routine will derive the magnitude of the velocity field
c ::: from the velocity field
c
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      integer    nxlo, nxhi, nylo, nyhi, nzlo, nzhi
      REAL_T     u, v, w
c      
c     ::::: lets punt if not in domain interior
c
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
         write(6,*) "FORT_DERMVEL: outside domain"
	 stop
      end if

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               u   = dat(i,j,k,1)
               v   = dat(i,j,k,2)
               w   = dat(i,j,k,3)
               e(i,j,k,1) = sqrt(u**2 + v**2 + w**2)
            end do
         end do
      end do

      end

      subroutine FORT_DERDVRHO (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                          lo,hi,domlo,domhi,delta,xlo,time,dt,
     &                          bc,level, grid_no)
c
c ::: This routine will derive C/RHO
c
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      integer    nxlo, nxhi, nylo, nyhi, nzlo, nzhi
      
c     ::::: lets punt if not in domain interior
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
         write(6,*) "FORT_DERDVRHO: outside domain"
	 stop
      end if

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = dat(i,j,k,2)/dat(i,j,k,1)
            end do
         end do
      end do

      end

      subroutine FORT_DERMPRHO (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                          lo,hi,domlo,domhi,delta,xlo,time,dt,
     &                          bc,level, grid_no)
c
c ::: This routine will derive RHO*C
c
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      integer    nxlo, nxhi, nylo, nyhi, nzlo, nzhi
c      
c     ::::: lets punt if not in domain interior
c
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
         write(6,*) "FORT_DERMPRHO: outside domain"
	 stop
      end if

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               e(i,j,k,1) = dat(i,j,k,2)*dat(i,j,k,1)
            end do
         end do
      end do

      end

      subroutine FORT_DERLGRHODUST (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                              lo,hi,domlo,domhi,delta,xlo,time,dt,
     &                              bc,level,grid_no)
c
c ::: This routine will derive log(RHO*C)
c
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j,k
      integer    nxlo, nxhi, nylo, nyhi, nzlo, nzhi
      REAL_T     dust, small
c      
c     ::::: lets punt if not in domain interior
c
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
         write(6,*) "FORT_DERLGRHODUST: outside domain"
	 stop
      end if

      small = 1.0e-10

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               dust = max(small,dat(i,j,k,2)*dat(i,j,k,1))
               e(i,j,k,1) = log10(dust)
            end do
         end do
      end do

      end

      subroutine FORT_DERMGVORT (vort,DIMS(vort),nv,dat,DIMS(dat),ncomp,
     &                           lo,hi,domlo,domhi,delta,xlo,time,dt,
     &                           bc,level,grid_no)
c
c ::: This routine will derive magnitude of vorticity from
c ::: the velocity field
c
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(vort)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     vort(DIMV(vort),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer   i,j,k
      integer   nxlo, nxhi, nylo, nyhi, nzlo, nzhi
      REAL_T    uy, uz, vx, vz, wx, wy, dx, dy, dz
      REAL_T    uycen, uzcen, uylo, uyhi, uzlo, uzhi
      REAL_T    vxcen, vzcen, vxlo, vxhi, vzlo, vzhi
      REAL_T    wxcen, wycen, wxlo, wxhi, wylo, wyhi
      REAL_T    vorfun

      logical   fixvlo_x, fixwlo_x, fixvhi_x, fixwhi_x
      logical   fixulo_y, fixwlo_y, fixuhi_y, fixwhi_y
      logical   fixulo_z, fixvlo_z, fixuhi_z, fixvhi_z
c
c     ::::: some useful macro definitions
c
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
c
c     ::::: statement functions that implement stencil
c
      uycen(i,j,k) = half*(U(i,j+1,k)-U(i,j-1,k))/dy
      uylo(i,j,k)  = (U(i,j+1,k)+three*U(i,j,k)-four*U(i,j-1,k))/(three*dy)
      uyhi(i,j,k)  =-(U(i,j-1,k)+three*U(i,j,k)-four*U(i,j+1,k))/(three*dy)

      uzcen(i,j,k) = half*(U(i,j,k+1)-U(i,j,k-1))/dz
      uzlo(i,j,k)  = (U(i,j,k+1)+three*U(i,j,k)-four*U(i,j,k-1))/(three*dz)
      uzhi(i,j,k)  =-(U(i,j,k-1)+three*U(i,j,k)-four*U(i,j,k+1))/(three*dz)

      vxcen(i,j,k) = half*(V(i+1,j,k)-V(i-1,j,k))/dx
      vxlo(i,j,k)  = (V(i+1,j,k)+three*V(i,j,k)-four*V(i-1,j,k))/(three*dx)
      vxhi(i,j,k)  =-(V(i-1,j,k)+three*V(i,j,k)-four*V(i+1,j,k))/(three*dx)

      vzcen(i,j,k) = half*(V(i,j,k+1)-V(i,j,k-1))/dz
      vzlo(i,j,k)  = (V(i,j,k+1)+three*V(i,j,k)-four*V(i,j,k-1))/(three*dz)
      vzhi(i,j,k)  =-(V(i,j,k-1)+three*V(i,j,k)-four*V(i,j,k+1))/(three*dz)

      wxcen(i,j,k) = half*(W(i+1,j,k)-W(i-1,j,k))/dx
      wxlo(i,j,k)  = (W(i+1,j,k)+three*W(i,j,k)-four*W(i-1,j,k))/(three*dx)
      wxhi(i,j,k)  =-(W(i-1,j,k)+three*W(i,j,k)-four*W(i+1,j,k))/(three*dx)

      wycen(i,j,k) = half*(W(i,j+1,k)-W(i,j-1,k))/dy
      wylo(i,j,k)  = (W(i,j+1,k)+three*W(i,j,k)-four*W(i,j-1,k))/(three*dy)
      wyhi(i,j,k)  =-(W(i,j-1,k)+three*W(i,j,k)-four*W(i,j+1,k))/(three*dy)

      vorfun(uy,uz,vx,vz,wx,wy) = sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)
c      
c     ::::: lets punt if not in domain interior
c
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
         write(6,*) "FORT_DERMGVORT: outside domain"
	 stop
      end if

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               uy = uycen(i,j,k)
               uz = uzcen(i,j,k)
               vx = vxcen(i,j,k)
               vz = vzcen(i,j,k)
               wx = wxcen(i,j,k)
               wy = wycen(i,j,k)
               vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end do

      fixvlo_x = ( (lo(1) .eq. domlo(1)) .and.
     &             (VLOX .eq. EXT_DIR .or. VLOX .eq. HOEXTRAP) )
      fixvhi_x = ( (hi(1) .eq. domhi(1)) .and.
     &             (VHIX .eq. EXT_DIR .or. VHIX .eq. HOEXTRAP) )
      fixwlo_x = ( (lo(1) .eq. domlo(1)) .and.
     &             (WLOX .eq. EXT_DIR .or. WLOX .eq. HOEXTRAP) )
      fixwhi_x = ( (hi(1) .eq. domhi(1)) .and.
     &             (WHIX .eq. EXT_DIR .or. WHIX .eq. HOEXTRAP) )

      fixulo_y = ( (lo(2) .eq. domlo(2)) .and.
     &             (ULOY .eq. EXT_DIR .or. ULOY .eq. HOEXTRAP) )
      fixuhi_y = ( (hi(2) .eq. domhi(2)) .and.
     &             (UHIY .eq. EXT_DIR .or. UHIY .eq. HOEXTRAP) )
      fixwlo_y = ( (lo(2) .eq. domlo(2)) .and.
     &             (WLOY .eq. EXT_DIR .or. WLOY .eq. HOEXTRAP) )
      fixwhi_y = ( (hi(2) .eq. domhi(2)) .and.
     &             (WHIY .eq. EXT_DIR .or. WHIY .eq. HOEXTRAP) )

      fixulo_z = ( (lo(3) .eq. domlo(3)) .and.
     &             (ULOZ .eq. EXT_DIR .or. ULOZ .eq. HOEXTRAP) )
      fixuhi_z = ( (hi(3) .eq. domhi(3)) .and.
     &             (UHIZ .eq. EXT_DIR .or. UHIZ .eq. HOEXTRAP) )
      fixvlo_z = ( (lo(3) .eq. domlo(3)) .and.
     &             (VLOZ .eq. EXT_DIR .or. VLOZ .eq. HOEXTRAP) )
      fixvhi_z = ( (hi(3) .eq. domhi(3)) .and.
     &             (VHIZ .eq. EXT_DIR .or. VHIZ .eq. HOEXTRAP) )
c
c     First do all the faces
c
      if (fixvlo_x .or. fixwlo_x) then
         i = lo(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               vx = cvmgt(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
               wx = cvmgt(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
               uy = uycen(i,j,k)
               wy = wycen(i,j,k)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixvhi_x .or. fixwhi_x) then
         i = hi(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               vx = cvmgt(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
               wx = cvmgt(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
               uy = uycen(i,j,k)
               wy = wycen(i,j,k)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixulo_y .or. fixwlo_y) then
         j = lo(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = cvmgt(uylo(i,j,k),uycen(i,j,k),fixulo_y)
               wy = cvmgt(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixuhi_y .or. fixwhi_y) then
         j = hi(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = cvmgt(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
               wy = cvmgt(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if

      if (fixulo_z .or. fixvlo_z) then
         k = lo(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = uycen(i,j,k)
               wy = wycen(i,j,k)
               uz = cvmgt(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
               vz = cvmgt(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
               vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
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
               uz = cvmgt(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
               vz = cvmgt(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
               vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
            end do
         end do
      end if
c
c     Next do all the edges
c
      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y)) then
         i = lo(1)
         j = lo(2)
         do k = lo(3),hi(3)
            vx = cvmgt(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            wx = cvmgt(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uy = cvmgt(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            wy = cvmgt(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y)) then
         i = hi(1)
         j = lo(2)
         do k = lo(3),hi(3)
            vx = cvmgt(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            wx = cvmgt(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uy = cvmgt(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            wy = cvmgt(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y)) then
         i = lo(1)
         j = hi(2)
         do k = lo(3),hi(3)
            vx = cvmgt(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            wx = cvmgt(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uy = cvmgt(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            wy = cvmgt(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y)) then
         i = hi(1)
         j = hi(2)
         do k = lo(3),hi(3)
            vx = cvmgt(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            wx = cvmgt(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uy = cvmgt(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            wy = cvmgt(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         k = lo(3)
         do j = lo(2),hi(2)
            vx = cvmgt(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            wx = cvmgt(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uy = uycen(i,j,k)
            wy = wycen(i,j,k)
            uz = cvmgt(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = cvmgt(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         k = lo(3)
         do j = lo(2),hi(2)
            vx = cvmgt(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            wx = cvmgt(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uy = uycen(i,j,k)
            wy = wycen(i,j,k)
            uz = cvmgt(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = cvmgt(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         k = hi(3)
         do j = lo(2),hi(2)
            vx = cvmgt(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            wx = cvmgt(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uy = uycen(i,j,k)
            wy = wycen(i,j,k)
            uz = cvmgt(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = cvmgt(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         k = hi(3)
         do j = lo(2),hi(2)
            vx = cvmgt(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            wx = cvmgt(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uy = uycen(i,j,k)
            wy = wycen(i,j,k)
            uz = cvmgt(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = cvmgt(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixulo_y .or. fixwlo_y) .and. (fixulo_z .or. fixvlo_z)) then
         j = lo(2)
         k = lo(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = cvmgt(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            wy = cvmgt(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = cvmgt(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = cvmgt(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixuhi_y .or. fixwhi_y) .and. (fixulo_z .or. fixvlo_z)) then
         j = hi(2)
         k = lo(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = cvmgt(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            wy = cvmgt(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = cvmgt(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = cvmgt(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixulo_y .or. fixwlo_y) .and. (fixuhi_z .or. fixvhi_z)) then
         j = lo(2)
         k = hi(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = cvmgt(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            wy = cvmgt(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = cvmgt(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = cvmgt(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if

      if ((fixuhi_y .or. fixwhi_y) .and. (fixuhi_z .or. fixvhi_z)) then
         j = hi(2)
         k = hi(3)
         do i = lo(1),hi(1)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = cvmgt(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            wy = cvmgt(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = cvmgt(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = cvmgt(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
         end do
      end if
c
c     Finally do all the corners
c
      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y) .and. 
     $     (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         j = lo(2)
         k = lo(3)
         vx = cvmgt(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         wx = cvmgt(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uy = cvmgt(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         wy = cvmgt(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = cvmgt(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = cvmgt(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y) .and. 
     $     (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         j = lo(2)
         k = lo(3)
         vx = cvmgt(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         wx = cvmgt(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uy = cvmgt(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         wy = cvmgt(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = cvmgt(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = cvmgt(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y) .and. 
     $     (fixulo_z .or. fixvlo_z)) then
         i = lo(1)
         j = hi(2)
         k = lo(3)
         vx = cvmgt(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         wx = cvmgt(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uy = cvmgt(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         wy = cvmgt(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = cvmgt(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = cvmgt(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y) .and. 
     $     (fixulo_z .or. fixvlo_z)) then
         i = hi(1)
         j = hi(2)
         k = lo(3)
         vx = cvmgt(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         wx = cvmgt(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uy = cvmgt(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         wy = cvmgt(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = cvmgt(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = cvmgt(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixulo_y .or. fixwlo_y) .and. 
     $     (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         j = lo(2)
         k = hi(3)
         vx = cvmgt(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         wx = cvmgt(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uy = cvmgt(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         wy = cvmgt(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = cvmgt(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = cvmgt(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixulo_y .or. fixwlo_y) .and. 
     $     (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         j = lo(2)
         k = hi(3)
         vx = cvmgt(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         wx = cvmgt(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uy = cvmgt(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         wy = cvmgt(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = cvmgt(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = cvmgt(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. (fixuhi_y .or. fixwhi_y) .and. 
     $     (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         j = hi(2)
         k = hi(3)
         vx = cvmgt(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         wx = cvmgt(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uy = cvmgt(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         wy = cvmgt(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = cvmgt(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = cvmgt(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
      end if

      if ((fixvhi_x .or. fixwhi_x) .and. (fixuhi_y .or. fixwhi_y) .and. 
     $     (fixuhi_z .or. fixvhi_z)) then
         i = hi(1)
         j = hi(2)
         k = hi(3)
         vx = cvmgt(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         wx = cvmgt(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uy = cvmgt(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         wy = cvmgt(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = cvmgt(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = cvmgt(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         vort(i,j,k,1) = vorfun(uy,uz,vx,vz,wx,wy)
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

      end

      subroutine FORT_DERMGDIVU (divu,DIMS(divu),nv,dat,DIMS(dat),ncomp,
     &                           lo,hi,domlo,domhi,delta,xlo,time,dt,
     &                           bc,level,grid_no)
c
c ::: This routine will derive magnitude of the divergence of velocity
c
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(divu)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     divu(DIMV(divu),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer   i,j,k
      integer   nxlo, nxhi, nylo, nyhi, nzlo, nzhi
      REAL_T    ux, vy, wz, dx, dy, dz
      REAL_T    uxcen, uxlo, uxhi
      REAL_T    vycen, vylo, vyhi
      REAL_T    wzcen, wzlo, wzhi
c
c     ::::: some useful macro definitions
c
#     define U(i,j,k) dat(i,j,k,1)
#     define V(i,j,k) dat(i,j,k,2)
#     define W(i,j,k) dat(i,j,k,3)

#     define ULOX bc(1,1,1)
#     define UHIX bc(1,2,1)
#     define VLOY bc(2,1,2)
#     define VHIY bc(2,2,2)
#     define WLOZ bc(3,1,2)
#     define WHIZ bc(3,2,2)
c
c     ::::: statement functions that implement stencil
c
      uxcen(i,j,k) = half*(U(i+1,j,k)-U(i-1,j,k))/dx
      uxlo(i,j,k) = (eight*U(i,j,k)-six*U(i+1,j,k)+U(i+2,j,k))/(three*dx)
      uxhi(i,j,k) = (eight*U(i,j,k)-six*U(i-1,j,k)+U(i-2,j,k))/(three*dx)

      vycen(i,j,k) = half*(V(i,j+1,k)-V(i,j-1,k))/dy
      vylo(i,j,k) = (eight*V(i,j,k)-six*V(i,j+1,k)+V(i,j+2,k))/(three*dy)
      vyhi(i,j,k) = (eight*V(i,j,k)-six*V(i,j-1,k)+V(i,j-2,k))/(three*dy)

      wzcen(i,j,k) = half*(W(i,j,k+1)-W(i,j,k-1))/dz
      wzlo(i,j,k) = (eight*W(i,j,k)-six*W(i,j,k+1)+W(i,j,k+2))/(three*dz)
      wzhi(i,j,k) = (eight*W(i,j,k)-six*W(i,j,k-1)+W(i,j,k-2))/(three*dz)

      call FORT_XVELFILL(dat(ARG_L1(dat),ARG_L2(dat),ARG_L3(dat),1),DIMS(dat),
     $                   domlo,domhi,delta,xlo,time,bc(1,1,1))
      call FORT_YVELFILL(dat(ARG_L1(dat),ARG_L2(dat),ARG_L3(dat),2),DIMS(dat),
     $                   domlo,domhi,delta,xlo,time,bc(1,1,2))
      call FORT_ZVELFILL(dat(ARG_L1(dat),ARG_L2(dat),ARG_L3(dat),3),DIMS(dat),
     $                   domlo,domhi,delta,xlo,time,bc(1,1,3))

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)
c      
c     ::::: lets punt if not in domain interior
c
      nxlo = max(0,domlo(1)-lo(1))
      nxhi = max(0,hi(1)-domhi(1))
      nylo = max(0,domlo(2)-lo(2))
      nyhi = max(0,hi(2)-domhi(2))
      nzlo = max(0,domlo(3)-lo(3))
      nzhi = max(0,hi(3)-domhi(3))

      if (nxlo+nxhi+nylo+nyhi+nzlo+nzhi .gt. 0) then
         write(6,*) "FORT_DERMGDIVU: outside domain"
	 stop
      end if
c
c     :: at physical bndries where an edge value is prescribed,
c     :: set the value in the outside cell so that a central
c     :: difference formula is equivalent to the higher order
c     :: one sided formula
c
      if (lo(1) .eq. domlo(1)) then
         i = lo(1)
         if (ULOX.eq.EXT_DIR) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i-1,j,k) = two*U(i-1,j,k) - U(i,j,k)
               end do
            end do
         else if (ULOX.eq.HOEXTRAP) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i-1,j,k) = uxlo(i,j,k)
               end do
            end do
	 end if
      end if
      if (hi(1) .eq. domhi(1)) then
         i = hi(1)
         if (UHIX.eq.EXT_DIR) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i+1,j,k) = two*U(i+1,j,k) - U(i,j,k)
               end do
            end do
         else if (UHIX.eq.HOEXTRAP) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  U(i+1,j,k) = uxhi(i,j,k)
               end do
            end do
	 end if
      end if
      if (lo(2) .eq. domlo(2)) then
         j = lo(2)
	 if (VLOY.eq.EXT_DIR) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j-1,k) = two*V(i,j-1,k) - V(i,j,k)
               end do
            end do
         else if (VLOY.eq.HOEXTRAP) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j-1,k) = vylo(i,j,k)
               end do
            end do
	 end if
      end if
      if (hi(2) .eq. domhi(2)) then
         j = hi(2)
	 if (VHIY.eq.EXT_DIR) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j+1,k) = two*V(i,j+1,k) - V(i,j,k)
               end do
            end do
	 else if (VHIY.eq.HOEXTRAP) then
            do k = lo(3), hi(3)
               do i = lo(1), hi(1)
                  V(i,j+1,k) = vyhi(i,j,k)
               end do
            end do
	 end if
      end if
      if (lo(3) .eq. domlo(3)) then
         k = lo(3)
	 if (WLOZ.eq.EXT_DIR) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k-1) = two*W(i,j,k-1) - W(i,j,k)
               end do
            end do
	 else if (WLOZ.eq.HOEXTRAP) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k-1) = wzlo(i,j,k)
               end do
            end do
	 end if
      end if
      if (hi(3) .eq. domhi(3)) then
         k = hi(3)
	 if (WHIZ.eq.EXT_DIR) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k+1) = two*W(i,j,k+1) - W(i,j,k)
               end do
            end do
	 else if (WHIZ.eq.HOEXTRAP) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  W(i,j,k+1) = wzhi(i,j,k)
               end do
            end do
	 end if
      end if

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ux = uxcen(i,j,k)
               vy = vycen(i,j,k)
               wz = wzcen(i,j,k)
               divu(i,j,k,1) = ux + vy + wz
            end do
         end do
      end do
c
c we overwrote the ghost cells above, so set them back below
c
      call FORT_XVELFILL(dat(ARG_L1(dat),ARG_L2(dat),ARG_L3(dat),1),DIMS(dat),
     $                   domlo,domhi,delta,xlo,time,bc(1,1,1))
      call FORT_YVELFILL(dat(ARG_L1(dat),ARG_L2(dat),ARG_L3(dat),1),DIMS(dat),
     $                   domlo,domhi,delta,xlo,time,bc(1,1,2))
      call FORT_ZVELFILL(dat(ARG_L1(dat),ARG_L2(dat),ARG_L3(dat),1),DIMS(dat),
     $                   domlo,domhi,delta,xlo,time,bc(1,1,3))

#     undef U
#     undef V      
#     undef W
#     undef ULOX
#     undef UHIX
#     undef VLOY
#     undef VHIY
#     undef WLOZ
#     undef WHIZ

      end

      subroutine FORT_GRADP_DIR (
     &     p,DIMS(p),
     &     gp,DIMS(gp),
     &     lo,hi,dir,dx)
c
c     compute a node centered pressure gradient in direction (dir)
c
      integer    DIMDEC(p)
      integer    DIMDEC(gp)
      integer     lo(SDIM),  hi(SDIM)
      integer    dir
      REAL_T     dx
      REAL_T   p(DIMV(p))
      REAL_T  gp(DIMV(gp))

      logical    invalid
      integer    i,j,k
      REAL_T     d
      integer    ilo, ihi, jlo, jhi, klo, khi

      d = fourth/dx
c
c     ::::: compute gradient on interior
c
      if (dir .eq. 0) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  gp(i,j,k) = d*(
     &                 p(i+1,j,k  )-p(i,j,k  )+p(i+1,j+1,k  )-p(i,j+1,k  )+
     $                 p(i+1,j,k+1)-p(i,j,k+1)+p(i+1,j+1,k+1)-p(i,j+1,k+1))
               end do
            end do
         end do
      else if (dir .eq. 1) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  gp(i,j,k) = d*(
     &                 p(i,j+1,k  )-p(i,j,k  )+p(i+1,j+1,k  )-p(i+1,j,k  )+
     $                 p(i,j+1,k+1)-p(i,j,k+1)+p(i+1,j+1,k+1)-p(i+1,j,k+1))
               end do
            end do
         end do
      else if (dir .eq. 2) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  gp(i,j,k) = d*(
     &                 p(i,  j,k+1)-p(i,  j,k)+p(i,  j+1,k+1)-p(i,  j+1,k)+
     $                 p(i+1,j,k+1)-p(i+1,j,k)+p(i+1,j+1,k+1)-p(i+1,j+1,k))
               end do
            end do
         end do
      else
         write(6,*) "FORT_GRADP_DIR: invalid dir = ",dir
	 stop
      end if

      end

      subroutine FORT_DERGRDPX (grdpx,DIMS(gp),nv,dat,DIMS(dat),ncomp,
     &                          lo,hi,domlo,domhi,delta,xlo,time,dt,
     &                          bc,level,grid_no)
c
c     This routine computes pressure gradient in x direciton
c
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(gp)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     grdpx(DIMV(gp),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no
c
      call FORT_GRADP_DIR (
     &     dat,DIMS(dat),grdpx,DIMS(gp),
     &     lo,hi,0,delta)

      end

      subroutine FORT_DERGRDPY (grdpy,DIMS(gp),nv,dat,DIMS(dat),ncomp,
     &                          lo,hi,domlo,domhi,delta,xlo,time,dt,
     &                          bc,level,grid_no)
c
c     This routine computes pressure gradient in Y direciton
c
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(gp)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     grdpy(DIMV(gp),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no
c
      call FORT_GRADP_DIR (
     &     dat,DIMS(dat),grdpy,DIMS(gp),
     &     lo,hi,1,delta)

      end

      subroutine FORT_DERGRDPZ (grdpz,DIMS(gp),nv,dat,DIMS(dat),ncomp,
     &                          lo,hi,domlo,domhi,delta,xlo,time,dt,
     &                          bc,level,grid_no)
c
c     This routine computes pressure gradient in Z direciton
c
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(gp)
      integer    DIMDEC(dat)
      integer    domlo(SDIM), domhi(SDIM)
      integer    nv, ncomp
      integer    bc(SDIM,2,ncomp)
      REAL_T     delta(SDIM), xlo(SDIM)
      REAL_T     time, dt
      REAL_T     grdpz(DIMV(gp),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no
c
      call FORT_GRADP_DIR (
     &     dat,DIMS(dat),grdpz,DIMS(gp),
     &     lo,hi,2,delta)

      end
