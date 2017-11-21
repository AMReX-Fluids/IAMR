
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

#define SDIM 2

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define SMALL 1.0e-10
#else
#define SMALL 1.0d-10
#endif

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

      subroutine FORT_DERMODGRADRHO (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                         lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                         bc,level,grid_no)
      implicit none
c
c     This routine will derive the modulus of the density gradients
c
      integer    lo(2), hi(2)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(2), domhi(2)
      integer    nv, ncomp
      integer    bc(2,2,ncomp)
      REAL_T     delta(2), xlo(2), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j
      REAL_T     drdx, drdy, tdx, tdy

      tdx = two*delta(1)
      tdy = two*delta(2)

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            drdx = (dat(i+1,j,1)-dat(i-1,j,1))/tdx
            drdy = (dat(i,j+1,1)-dat(i,j-1,1))/tdy
            e(i,j,1) = sqrt( drdx*drdx + drdy*drdy )
         end do
      end do

      end

      subroutine FORT_DERKENG (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                         lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                         bc,level,grid_no)
      implicit none
c
c     This routine will derive kinetic energy from density
c     and the velocity field.
c
      integer    lo(2), hi(2)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(2), domhi(2)
      integer    nv, ncomp
      integer    bc(2,2,ncomp)
      REAL_T     delta(2), xlo(2), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j
      REAL_T     rho, u, v

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
	    rho = dat(i,j,1)
	    u   = dat(i,j,2)
	    v   = dat(i,j,3)
	    e(i,j,1) = half*rho*(u**2 + v**2)
	 end do
      end do

      end

      subroutine FORT_DERLOGS (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                         lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                         bc,level,grid_no)
      implicit none
c
c     This routine will derive log of given scalar quantity
c
      integer    lo(2), hi(2)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(2), domhi(2)
      integer    nv, ncomp
      integer    bc(2,2,ncomp)
      REAL_T     delta(2), xlo(2), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j
      REAL_T     rho

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
	    rho = max(dat(i,j,1),SMALL)
	    e(i,j,1) = log10(rho)
	 end do
      end do

      end

      subroutine FORT_DERMVEL (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                         lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                         bc,level,grid_no)
      implicit none
c
c ::: This routine will derive the magnitude of the velocity field
c ::: from the velocity field
c

      integer    lo(2), hi(2)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(2), domhi(2)
      integer    nv, ncomp
      integer    bc(2,2,ncomp)
      REAL_T     delta(2), xlo(2), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j
      REAL_T     u, v

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
	    u   = dat(i,j,1)
	    v   = dat(i,j,2)
	    e(i,j,1) = sqrt(u**2 + v**2)
	 end do
      end do

      end

      subroutine FORT_DERDVRHO (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                         lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                         bc,level,grid_no)
      implicit none
c
c ::: This routine will derive C/RHO
c
      integer    lo(2), hi(2)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(2), domhi(2)
      integer    nv, ncomp
      integer    bc(2,2,ncomp)
      REAL_T     delta(2), xlo(2), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
	    e(i,j,1) = dat(i,j,2)/dat(i,j,1)
	 end do
      end do

      end

      subroutine FORT_DERMPRHO (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                         lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                         bc,level,grid_no)
      implicit none
c
c ::: This routine will derive RHO*C
c
      integer    lo(2), hi(2)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(2), domhi(2)
      integer    nv, ncomp
      integer    bc(2,2,ncomp)
      REAL_T     delta(2), xlo(2), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
	    e(i,j,1) = dat(i,j,2)*dat(i,j,1)
	 end do
      end do

      end

      subroutine FORT_DERLGRHODUST (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                              lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                              bc,level,grid_no)
      implicit none
c
c ::: This routine will derive log(RHO*C)
c
      integer    lo(2), hi(2)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(2), domhi(2)
      integer    nv, ncomp
      integer    bc(2,2,ncomp)
      REAL_T     delta(2), xlo(2), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer    i,j
      REAL_T     dust

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
	    dust = max(SMALL,dat(i,j,2)*dat(i,j,1))
	    e(i,j,1) = log10(dust)
	 end do
      end do

      end

      subroutine FORT_DERMGVORT (vort,DIMS(vort),nv,dat,DIMS(dat),ncomp,
     &                           lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                           bc,level,grid_no)
      implicit none
c
c ::: This routine will derive magnitude of vorticity from
c ::: the velocity field
c
      integer    lo(2), hi(2)
      integer    DIMDEC(vort)
      integer    DIMDEC(dat)
      integer    domlo(2), domhi(2)
      integer    nv, ncomp
      integer    bc(2,2,ncomp)
      REAL_T     delta(2), xlo(2), time, dt
      REAL_T     vort(DIMV(vort),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer   i,j
      REAL_T    vx, uy, dx, dy
      logical   fixlft, fixrgt, fixbot, fixtop
      REAL_T    vxcen, vxlft, vxrgt, uycen, uybot, uytop, vorfun
c
c     ::::: some useful macro definitions
c
#     define U(i,j) dat(i,j,1)
#     define V(i,j) dat(i,j,2)
#     define VLOX bc(1,1,2)
#     define VHIX bc(1,2,2)
#     define ULOY bc(2,1,1)
#     define UHIY bc(2,2,1)
c
c     ::::: statement functions that implement stencil
c
      vxcen(i,j) = half*(V(i+1,j)-V(i-1,j))/dx
      vxlft(i,j) = (V(i+1,j)+three*V(i,j)-four*V(i-1,j))/(three*dx)
      vxrgt(i,j) = (four*V(i+1,j)-three*V(i,j)-V(i-1,j))/(three*dx)
      uycen(i,j) = half*(U(i,j+1)-U(i,j-1))/dy
      uybot(i,j) = (U(i,j+1)+three*U(i,j)-four*U(i,j-1))/(three*dy)
      uytop(i,j) = (four*U(i,j+1)-three*U(i,j)-U(i,j-1))/(three*dy)
      vorfun(vx,uy) = vx - uy

      dx = delta(1)
      dy = delta(2)

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
	    vx  = vxcen(i,j)
	    uy  = uycen(i,j)
	    vort(i,j,1) = vorfun(vx,uy)
	 end do
      end do

      fixlft = ( (lo(1) .eq. domlo(1)) .and. 
     &           (VLOX .eq. EXT_DIR .or. VLOX .eq. HOEXTRAP) )
      fixrgt = ( (hi(1) .eq. domhi(1)) .and. 
     &           (VHIX .eq. EXT_DIR .or. VHIX .eq. HOEXTRAP) )
      fixbot = ( (lo(2) .eq. domlo(2)) .and. 
     &           (ULOY .eq. EXT_DIR .or. ULOY .eq. HOEXTRAP) )
      fixtop = ( (hi(2) .eq. domhi(2)) .and. 
     &           (UHIY .eq. EXT_DIR .or. UHIY .eq. HOEXTRAP) )
c
c     ::::: handle special bndry conditions at the left edge
c
      if (fixlft) then
         i = lo(1)
         do j = lo(2), hi(2)
	    vx  = vxlft(i,j)
	    uy  = uycen(i,j)
    	    vort(i,j,1) = vorfun(vx,uy)
	 end do
      end if
c
c     ::::: handle special bndry conditions on the right
c
      if (fixrgt) then
         i = hi(1)
         do j = lo(2), hi(2)
	    vx  = vxrgt(i,j)
	    uy  = uycen(i,j)
    	    vort(i,j,1) = vorfun(vx,uy)
	 end do
      end if
c
c     ::::: handle special bndry conditions on bottom
c
      if (fixbot) then
         j = lo(2)
         do i = lo(1), hi(1)
	    vx  = vxcen(i,j)
	    uy  = uybot(i,j)
    	    vort(i,j,1) = vorfun(vx,uy)
	 end do
      end if
c
c     ::::: handle special bndry conditions on top
c
      if (fixtop) then
         j = hi(2)
         do i = lo(1), hi(1)
	    vx  = vxcen(i,j)
	    uy  = uytop(i,j)
    	    vort(i,j,1) = vorfun(vx,uy)
	 end do
      end if
c
c     ::::: check corners
c
      if (fixlft .and. fixbot) then
         i = lo(1)
	 j = lo(2)
         vx = vxlft(i,j)
	 uy = uybot(i,j)
	 vort(i,j,1) = vorfun(vx,uy)
      end if
      if (fixlft .and. fixtop) then
         i = lo(1)
	 j = hi(2)
         vx = vxlft(i,j)
         uy = uytop(i,j)
	 vort(i,j,1) = vorfun(vx,uy)
      end if
      if (fixrgt .and. fixtop) then
         i = hi(1)
	 j = hi(2)
	 vx = vxrgt(i,j)
         uy = uytop(i,j)
	 vort(i,j,1) = vorfun(vx,uy)
      end if
      if (fixrgt .and. fixbot) then
         i = hi(1)
	 j = lo(2)
	 vx = vxrgt(i,j)
	 uy = uybot(i,j)
	 vort(i,j,1) = vorfun(vx,uy)
      end if

#     undef U
#     undef V      
#     undef VLOX
#     undef VHIX
#     undef ULOY
#     undef UHIY

      end

      subroutine FORT_DERMGVORT2 (vort,DIMS(vort),nv,dat,DIMS(dat),ncomp,
     &                            lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                            bc,level,grid_no)
      implicit none
c
c ::: This routine will derive magnitude of vorticity from
c ::: the velocity field
c
      integer    lo(2), hi(2)
      integer    DIMDEC(vort)
      integer    DIMDEC(dat)
      integer    domlo(2), domhi(2)
      integer    nv, ncomp
      integer    bc(2,2,ncomp)
      REAL_T     delta(2), xlo(2), time, dt
      REAL_T     vort(DIMV(vort),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer   i,j
      REAL_T    uy, vx, dx, dy
      REAL_T    uycen, uylo, uyhi
      REAL_T    vxcen, vxlo, vxhi
      REAL_T    vorfun
c
c     ::::: some useful macro definitions
c
#     define U(i,j) dat(i,j,1)
#     define V(i,j) dat(i,j,2)

#     define ULOY bc(2,1,1)
#     define UHIY bc(2,2,1)

#     define VLOX bc(1,1,2)
#     define VHIX bc(1,2,2)
c
c     ::::: statement functions that implement stencil
c
      uycen(i,j) = half*(U(i,j+1)-U(i,j-1))/dy
      uylo(i,j) = (eight*U(i,j)-six*U(i,j+1)+U(i,j+2))/(three*dy)
      uyhi(i,j) = (eight*U(i,j)-six*U(i,j-1)+U(i,j-2))/(three*dy)

      vxcen(i,j) = half*(V(i+1,j)-V(i-1,j))/dx
      vxlo(i,j) = (eight*V(i,j)-six*V(i+1,j)+V(i+2,j))/(three*dx)
      vxhi(i,j) = (eight*V(i,j)-six*V(i-1,j)+V(i-2,j))/(three*dx)

      vorfun(vx,uy) = vx - uy

      dx = delta(1)
      dy = delta(2)
c
c     :: at physical bndries where an edge value is prescribed,
c     :: set the value in the outside cell so that a central
c     :: difference formula is equivalent to the higher order
c     :: one sided formula
c
      if (lo(1) .eq. domlo(1)) then
         i = lo(1)
         if (VLOX.eq.EXT_DIR .or. VLOX.eq.HOEXTRAP) then
	   do j = lo(2), hi(2)
	      V(i-1,j) = vxlo(i,j)
	   end do
	 end if
      end if
      if (hi(1) .eq. domhi(1)) then
         i = hi(1)
         if (VHIX.eq.EXT_DIR .or. VHIX.eq.HOEXTRAP) then
	   do j = lo(2), hi(2)
	      V(i+1,j) = vxhi(i,j)
	   end do
	 end if
      end if
      if (lo(2) .eq. domlo(2)) then
         j = lo(2)
	 if (ULOY.eq.EXT_DIR .or. ULOY.eq.HOEXTRAP) then
	   do i = lo(1), hi(1)
	      U(i,j-1) = uylo(i,j)
	   end do
	 end if
      end if
      if (hi(2) .eq. domhi(2)) then
         j = hi(2)
	 if (UHIY.eq.EXT_DIR .or. UHIY.eq.HOEXTRAP) then
	   do i = lo(1), hi(1)
	      U(i,j+1) = uyhi(i,j)
	   end do
	 end if
      end if

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            uy = uycen(i,j)
            vx = vxcen(i,j)
            vort(i,j,1) = vorfun(vx,uy)
         end do
      end do

#     undef U
#     undef V      
#     undef ULOY
#     undef UHIY
#     undef VLOX
#     undef VHIX

      end

      subroutine FORT_DERMGSTRNRT (strnrt,DIMS(strnrt),nv,dat,DIMS(dat),ncomp,
     &                             lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                             bc,level,grid_no)
      implicit none
c
c ::: This routine will derive magnitude of rate-of-strain tensor from
c ::: the velocity field
c
      integer    lo(2), hi(2)
      integer    DIMDEC(strnrt)
      integer    DIMDEC(dat)
      integer    domlo(2), domhi(2)
      integer    nv, ncomp
      integer    bc(2,2,ncomp)
      REAL_T     delta(2), xlo(2), time, dt
      REAL_T     strnrt(DIMV(strnrt),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer   i,j,k
      REAL_T    ux, vx, uy, vy, dx, dy
      logical   fixlft, fixrgt, fixbot, fixtop
      REAL_T    uxcen, uxlft, uxrgt, vxcen, vxlft, vxrgt
      REAL_T    uycen, uybot, uytop, vycen, vybot, vytop
      REAL_T    strnrt_fun
c
c     ::::: some useful macro definitions
c
#     define U(i,j) dat(i,j,1)
#     define V(i,j) dat(i,j,2)

#     define ULOX bc(1,1,1)
#     define UHIX bc(1,2,1)

#     define VLOX bc(1,1,2)
#     define VHIX bc(1,2,2)

#     define ULOY bc(2,1,1)
#     define UHIY bc(2,2,1)

#     define VLOY bc(2,1,2)
#     define VHIY bc(2,2,2)

c
c     ::::: statement functions that implement stencil
c
      uxcen(i,j) = half*(U(i+1,j)-U(i-1,j))/dx
      uxlft(i,j) = (U(i+1,j)+three*U(i,j)-four*U(i-1,j))/(three*dx)
      uxrgt(i,j) = (four*U(i+1,j)-three*U(i,j)-U(i-1,j))/(three*dx)
      
      vxcen(i,j) = half*(V(i+1,j)-V(i-1,j))/dx
      vxlft(i,j) = (V(i+1,j)+three*V(i,j)-four*V(i-1,j))/(three*dx)
      vxrgt(i,j) = (four*V(i+1,j)-three*V(i,j)-V(i-1,j))/(three*dx)

      uycen(i,j) = half*(U(i,j+1)-U(i,j-1))/dy
      uybot(i,j) = (U(i,j+1)+three*U(i,j)-four*U(i,j-1))/(three*dy)
      uytop(i,j) = (four*U(i,j+1)-three*U(i,j)-U(i,j-1))/(three*dy)

      vycen(i,j) = half*(V(i,j+1)-V(i,j-1))/dy
      vybot(i,j) = (V(i,j+1)+three*V(i,j)-four*V(i,j-1))/(three*dy)
      vytop(i,j) = (four*V(i,j+1)-three*V(i,j)-V(i,j-1))/(three*dy)
c
c     ::::: function for computing magnitude of rate of strain tensor
c
      strnrt_fun(ux,vx,uy,vy) = sqrt(ux**2+half*(vx+uy)**2+vy**2)

      dx = delta(1)
      dy = delta(2)

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            ux  = uxcen(i,j)
            vx  = vxcen(i,j)
            uy  = uycen(i,j)
            vy  = vycen(i,j)
            strnrt(i,j,1) = strnrt_fun(ux,vx,uy,vy)
         end do
      end do

      fixlft = ( (lo(1) .eq. domlo(1)) .and. 
     &           ( (ULOX .eq. EXT_DIR .or. ULOX .eq. HOEXTRAP) .or. 
     &             (VLOX .eq. EXT_DIR .or. VLOX .eq. HOEXTRAP) ) )
      fixrgt = ( (hi(1) .eq. domhi(1)) .and. 
     &           ( (UHIX .eq. EXT_DIR .or. UHIX .eq. HOEXTRAP) .or.
     &             (VHIX .eq. EXT_DIR .or. VHIX .eq. HOEXTRAP) ) )
      fixbot = ( (lo(2) .eq. domlo(2)) .and. 
     &           ( (ULOY .eq. EXT_DIR .or. ULOY .eq. HOEXTRAP) .or.
     &           (VLOY .eq. EXT_DIR .or. VLOY .eq. HOEXTRAP) ) )
      fixtop = ( (hi(2) .eq. domhi(2)) .and. 
     &           ( (UHIY .eq. EXT_DIR .or. UHIY .eq. HOEXTRAP) .or.
     &           (VHIY .eq. EXT_DIR .or. VHIY .eq. HOEXTRAP) ) )
c
c     ::::: handle special bndry conditions at the left edge
c
      if (fixlft) then
         i = lo(1)
         do j = lo(2), hi(2)
            ux  = uxlft(i,j)
            vx  = vxlft(i,j)
            uy  = uycen(i,j)
            vy  = vycen(i,j)
            strnrt(i,j,1) = strnrt_fun(ux,vx,uy,vy)
         end do
      end if
c
c     ::::: handle special bndry conditions on the right
c
      if (fixrgt) then
         i = hi(1)
         do j = lo(2), hi(2)
            ux  = uxrgt(i,j)
            vx  = vxrgt(i,j)
            uy  = uycen(i,j)
            vy  = vycen(i,j)
            strnrt(i,j,1) = strnrt_fun(ux,vx,uy,vy)
         end do
      end if
c
c     ::::: handle special bndry conditions on bottom
c
      if (fixbot) then
         j = lo(2)
         do i = lo(1), hi(1)
            ux  = uxcen(i,j)
            vx  = vxcen(i,j)
            uy  = uybot(i,j)
            vy  = vybot(i,j)
            strnrt(i,j,1) = strnrt_fun(ux,vx,uy,vy)
         end do
      end if
c
c     ::::: handle special bndry conditions on top
c
      if (fixtop) then
         j = hi(2)
         do i = lo(1), hi(1)
            ux  = uxcen(i,j)
            vx  = vxcen(i,j)
            uy  = uytop(i,j)
            vy  = vytop(i,j)
            strnrt(i,j,1) = strnrt_fun(ux,vx,uy,vy)
         end do
      end if
c
c     ::::: check corners
c
      if (fixlft .and. fixbot) then
         i = lo(1)
         j = lo(2)
         ux = uxlft(i,j)
         vx = vxlft(i,j)
         uy = uybot(i,j)
         vy = vybot(i,j)
         strnrt(i,j,1) = strnrt_fun(ux,vx,uy,vy)
      end if
      if (fixlft .and. fixtop) then
         i = lo(1)
         j = hi(2)
         ux = uxlft(i,j)
         vx = vxlft(i,j)
         uy = uytop(i,j)
         vy = vytop(i,j)
         strnrt(i,j,1) = strnrt_fun(ux,vx,uy,vy)
      end if
      if (fixrgt .and. fixtop) then
         i = hi(1)
         j = hi(2)
         ux = uxrgt(i,j)
         vx = vxrgt(i,j)
         uy = uytop(i,j)
         vy = vytop(i,j)
         strnrt(i,j,1) = strnrt_fun(ux,vx,uy,vy)
      end if
      if (fixrgt .and. fixbot) then
         i = hi(1)
         j = lo(2)
         ux = uxrgt(i,j)
         vx = vxrgt(i,j)
         uy = uybot(i,j)
         vy = vybot(i,j)
         strnrt(i,j,1) = strnrt_fun(ux,vx,uy,vy)
      end if

#     undef U
#     undef V      
#     undef VLOX
#     undef VHIX
#     undef ULOY
#     undef UHIY

      end

      subroutine FORT_DERVISC (visc,DIMS(visc),nv,dat,DIMS(dat),ncomp,
     &                         lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                         bc,level,grid_no)
      implicit none

#include <NSCOMM_F.H>
c
c ::: This routine will compute the 
c ::: viscosity based on the velocity field
c
      integer    lo(2), hi(2)
      integer    DIMDEC(visc)
      integer    DIMDEC(dat)
      integer    domlo(2), domhi(2)
      integer    nv, ncomp
      integer    bc(2,2,ncomp)
      REAL_T     delta(2), xlo(2), time, dt
      REAL_T     visc(DIMV(visc),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer   i,j,k
      REAL_T    ux, vx, uy, vy, dx, dy
      logical   fixlft, fixrgt, fixbot, fixtop
      REAL_T    uxcen, uxlft, uxrgt, vxcen, vxlft, vxrgt
      REAL_T    uycen, uybot, uytop, vycen, vybot, vytop
      REAL_T    strnrt, strnrt_fun, relax, relax_exp, relax_tlr
      REAL_T    visc_pseudoplastic, visc_bingham, visc_dilatant
c
c     ::::: some useful macro definitions
c
#     define U(i,j) dat(i,j,1)
#     define V(i,j) dat(i,j,2)

#     define ULOX bc(1,1,1)
#     define UHIX bc(1,2,1)

#     define VLOX bc(1,1,2)
#     define VHIX bc(1,2,2)

#     define ULOY bc(2,1,1)
#     define UHIY bc(2,2,1)

#     define VLOY bc(2,1,2)
#     define VHIY bc(2,2,2)

c
c     ::::: statement functions that implement stencil
c
      uxcen(i,j) = half*(U(i+1,j)-U(i-1,j))/dx
      uxlft(i,j) = (U(i+1,j)+three*U(i,j)-four*U(i-1,j))/(three*dx)
      uxrgt(i,j) = (four*U(i+1,j)-three*U(i,j)-U(i-1,j))/(three*dx)
      
      vxcen(i,j) = half*(V(i+1,j)-V(i-1,j))/dx
      vxlft(i,j) = (V(i+1,j)+three*V(i,j)-four*V(i-1,j))/(three*dx)
      vxrgt(i,j) = (four*V(i+1,j)-three*V(i,j)-V(i-1,j))/(three*dx)

      uycen(i,j) = half*(U(i,j+1)-U(i,j-1))/dy
      uybot(i,j) = (U(i,j+1)+three*U(i,j)-four*U(i,j-1))/(three*dy)
      uytop(i,j) = (four*U(i,j+1)-three*U(i,j)-U(i,j-1))/(three*dy)

      vycen(i,j) = half*(V(i,j+1)-V(i,j-1))/dy
      vybot(i,j) = (V(i,j+1)+three*V(i,j)-four*V(i,j-1))/(three*dy)
      vytop(i,j) = (four*V(i,j+1)-three*V(i,j)-V(i,j-1))/(three*dy)
c
c     ::::: function for computing magnitude of rate of strain tensor
c
      strnrt_fun(ux,vx,uy,vy) = sqrt(ux**2+half*(vx+uy)**2+vy**2)
c
c     ::::: Functions for computing the regularised term in the viscosity
c     ::::: When the strain rate is small, Taylor expand the exponential
c     ::::: to avoid floating point errors.
c
      relax_exp(strnrt) = (one-exp(-strnrt/eps))/strnrt
      relax_tlr(strnrt) = one/eps*
     &     (one
     &      -half*strnrt/eps
     &      +sixth*(strnrt/eps)**2
     &      -sixth*fourth*(strnrt/eps)**3
     &      +sixth*fourth*fifth*(strnrt/eps)**4)
c
c     ::::: Functions for computing the regularised stress for different
c     ::::: power law behaviours 
c
      ! Pseudoplastic (shear-thinning)
      ! Regularise power-law term and yield stress term
      visc_pseudoplastic(strnrt,relax) = 
     &      (mu*two**((n-one)/two)*strnrt**n+half*tau)*relax
      ! Bingham fluid
      ! Regularise yield stress term only
      visc_bingham(strnrt,relax) = 
     &      mu+half*tau*relax
      ! Dilatant (shear-thickening)
      ! Regularise yield stress term only
      visc_dilatant(strnrt,relax) = 
     &      mu*two**((n-one)/two)*strnrt**(n-one)+half*tau*relax

      if (varvisc .eq. 0) then
c     ::::: If non-variable, viscosity coefficient equals mu everywhere

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            visc(i,j,1) = mu
         end do
      end do

      else 
c     ::::: Compute Herschel-Bulkley viscosity

      dx = delta(1)
      dy = delta(2)

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            ux  = uxcen(i,j)
            vx  = vxcen(i,j)
            uy  = uycen(i,j)
            vy  = vycen(i,j)
            strnrt = strnrt_fun(ux,vx,uy,vy)
            ! Compute relaxed term based on strain rate magnitude
            if (strnrt > 1.0e-8) then
               relax = relax_exp(strnrt)
            else
               relax = relax_tlr(strnrt)
            end if
            ! Calculate viscosity based on flow index value
            if (n .lt. one) then
               visc(i,j,1) = visc_pseudoplastic(strnrt,relax)
            else if (n .eq. one) then
               visc(i,j,1) = visc_bingham(strnrt,relax)
            else 
               visc(i,j,1) = visc_dilatant(strnrt,relax)
            end if
         end do
      end do

      fixlft = ( (lo(1) .eq. domlo(1)) .and. 
     &           ( (ULOX .eq. EXT_DIR .or. ULOX .eq. HOEXTRAP .or. ULOX.eq.FOEXTRAP) .or. 
     &             (VLOX .eq. EXT_DIR .or. VLOX .eq. HOEXTRAP .or. ULOX.eq.FOEXTRAP) ) )
      fixrgt = ( (hi(1) .eq. domhi(1)) .and. 
     &           ( (UHIX .eq. EXT_DIR .or. UHIX .eq. HOEXTRAP .or. UHIX.eq.FOEXTRAP) .or.
     &             (VHIX .eq. EXT_DIR .or. VHIX .eq. HOEXTRAP .or. VHIX.eq.FOEXTRAP) ) )
      fixbot = ( (lo(2) .eq. domlo(2)) .and. 
     &           ( (ULOY .eq. EXT_DIR .or. ULOY .eq. HOEXTRAP .or. ULOY.eq.FOEXTRAP) .or.
     &             (VLOY .eq. EXT_DIR .or. VLOY .eq. HOEXTRAP .or. VLOY.eq.FOEXTRAP) ) )
      fixtop = ( (hi(2) .eq. domhi(2)) .and. 
     &           ( (UHIY .eq. EXT_DIR .or. UHIY .eq. HOEXTRAP .or. UHIY.eq.FOEXTRAP) .or.
     &             (VHIY .eq. EXT_DIR .or. VHIY .eq. HOEXTRAP .or. VHIY.eq.FOEXTRAP) ) )

c
c     ::::: handle special bndry conditions at the left edge
c
      if (fixlft) then
         i = lo(1)
         do j = lo(2), hi(2)
            ux  = uxlft(i,j)
            vx  = vxlft(i,j)
            uy  = uycen(i,j)
            vy  = vycen(i,j)
            strnrt = strnrt_fun(ux,vx,uy,vy)
            if (strnrt > 1.0e-8) then
               relax = relax_exp(strnrt)
            else
               relax = relax_tlr(strnrt)
            end if
            if (n .lt. one) then
               visc(i,j,1) = visc_pseudoplastic(strnrt,relax)
            else if (n .eq. one) then
               visc(i,j,1) = visc_bingham(strnrt,relax)
            else 
               visc(i,j,1) = visc_dilatant(strnrt,relax)
            end if
         end do
      end if
c
c     ::::: handle special bndry conditions on the right
c
      if (fixrgt) then
         i = hi(1)
         do j = lo(2), hi(2)
            ux  = uxrgt(i,j)
            vx  = vxrgt(i,j)
            uy  = uycen(i,j)
            vy  = vycen(i,j)
            strnrt = strnrt_fun(ux,vx,uy,vy)
            if (strnrt > 1.0e-8) then
               relax = relax_exp(strnrt)
            else
               relax = relax_tlr(strnrt)
            end if
            if (n .lt. one) then
               visc(i,j,1) = visc_pseudoplastic(strnrt,relax)
            else if (n .eq. one) then
               visc(i,j,1) = visc_bingham(strnrt,relax)
            else 
               visc(i,j,1) = visc_dilatant(strnrt,relax)
            end if
         end do
      end if
c
c     ::::: handle special bndry conditions on bottom
c
      if (fixbot) then
         j = lo(2)
         do i = lo(1), hi(1)
            ux  = uxcen(i,j)
            vx  = vxcen(i,j)
            uy  = uybot(i,j)
            vy  = vybot(i,j)
            strnrt = strnrt_fun(ux,vx,uy,vy)
            if (strnrt > 1.0e-8) then
               relax = relax_exp(strnrt)
            else
               relax = relax_tlr(strnrt)
            end if
            if (n .lt. one) then
               visc(i,j,1) = visc_pseudoplastic(strnrt,relax)
            else if (n .eq. one) then
               visc(i,j,1) = visc_bingham(strnrt,relax)
            else 
               visc(i,j,1) = visc_dilatant(strnrt,relax)
            end if
         end do
      end if
c
c     ::::: handle special bndry conditions on top
c
      if (fixtop) then
         j = hi(2)
         do i = lo(1), hi(1)
            ux  = uxcen(i,j)
            vx  = vxcen(i,j)
            uy  = uytop(i,j)
            vy  = vytop(i,j)
            strnrt = strnrt_fun(ux,vx,uy,vy)
            if (strnrt > 1.0e-8) then
               relax = relax_exp(strnrt)
            else
               relax = relax_tlr(strnrt)
            end if
            if (n .lt. one) then
               visc(i,j,1) = visc_pseudoplastic(strnrt,relax)
            else if (n .eq. one) then
               visc(i,j,1) = visc_bingham(strnrt,relax)
            else 
               visc(i,j,1) = visc_dilatant(strnrt,relax)
            end if
         end do
      end if
c
c     ::::: check corners
c
      if (fixlft .and. fixbot) then
         i = lo(1)
         j = lo(2)
         ux = uxlft(i,j)
         vx = vxlft(i,j)
         uy = uybot(i,j)
         vy = vybot(i,j)
         strnrt = strnrt_fun(ux,vx,uy,vy)
         if (strnrt > 1.0e-8) then
             relax = relax_exp(strnrt)
         else
             relax = relax_tlr(strnrt)
         end if
         if (n .lt. one) then
             visc(i,j,1) = visc_pseudoplastic(strnrt,relax)
         else if (n .eq. one) then
             visc(i,j,1) = visc_bingham(strnrt,relax)
         else 
             visc(i,j,1) = visc_dilatant(strnrt,relax)
         end if
      end if
      if (fixlft .and. fixtop) then
         i = lo(1)
         j = hi(2)
         ux = uxlft(i,j)
         vx = vxlft(i,j)
         uy = uytop(i,j)
         vy = vytop(i,j)
         strnrt = strnrt_fun(ux,vx,uy,vy)
         if (strnrt > 1.0e-8) then
             relax = relax_exp(strnrt)
         else
             relax = relax_tlr(strnrt)
         end if
         if (n .lt. one) then
             visc(i,j,1) = visc_pseudoplastic(strnrt,relax)
         else if (n .eq. one) then
             visc(i,j,1) = visc_bingham(strnrt,relax)
         else 
             visc(i,j,1) = visc_dilatant(strnrt,relax)
         end if
      end if
      if (fixrgt .and. fixtop) then
         i = hi(1)
         j = hi(2)
         ux = uxrgt(i,j)
         vx = vxrgt(i,j)
         uy = uytop(i,j)
         vy = vytop(i,j)
         strnrt = strnrt_fun(ux,vx,uy,vy)
         if (strnrt > 1.0e-8) then
             relax = relax_exp(strnrt)
         else
             relax = relax_tlr(strnrt)
         end if
         if (n .lt. one) then
             visc(i,j,1) = visc_pseudoplastic(strnrt,relax)
         else if (n .eq. one) then
             visc(i,j,1) = visc_bingham(strnrt,relax)
         else 
             visc(i,j,1) = visc_dilatant(strnrt,relax)
         end if
      end if
      if (fixrgt .and. fixbot) then
         i = hi(1)
         j = lo(2)
         ux = uxrgt(i,j)
         vx = vxrgt(i,j)
         uy = uybot(i,j)
         vy = vybot(i,j)
         strnrt = strnrt_fun(ux,vx,uy,vy)
         if (strnrt > 1.0e-8) then
             relax = relax_exp(strnrt)
         else
             relax = relax_tlr(strnrt)
         end if
         if (n .lt. one) then
             visc(i,j,1) = visc_pseudoplastic(strnrt,relax)
         else if (n .eq. one) then
             visc(i,j,1) = visc_bingham(strnrt,relax)
         else 
             visc(i,j,1) = visc_dilatant(strnrt,relax)
         end if
      end if

      end if

#     undef U
#     undef V      
#     undef VLOX
#     undef VHIX
#     undef ULOY
#     undef UHIY

      end

      subroutine FORT_DERSTRESS (stress,DIMS(stress),nv,dat,DIMS(dat),ncomp,
     &                           lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                           bc,level,grid_no)
      implicit none

#include <NSCOMM_F.H>

c
c ::: This routine will compute the 
c ::: stressosity based on the velocity field
c
      integer    lo(2), hi(2)
      integer    DIMDEC(stress)
      integer    DIMDEC(dat)
      integer    domlo(2), domhi(2)
      integer    nv, ncomp
      integer    bc(2,2,ncomp)
      REAL_T     delta(2), xlo(2), time, dt
      REAL_T     stress(DIMV(stress),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer   i,j,k
      REAL_T    ux, vx, uy, vy, dx, dy
      logical   fixlft, fixrgt, fixbot, fixtop
      REAL_T    uxcen, uxlft, uxrgt, vxcen, vxlft, vxrgt
      REAL_T    uycen, uybot, uytop, vycen, vybot, vytop
      REAL_T    strnrt, strnrt_fun, relax
      REAL_T    stress_pseudoplastic, stress_bingham, stress_dilatant
c
c     ::::: some useful macro definitions
c
#     define U(i,j) dat(i,j,1)
#     define V(i,j) dat(i,j,2)

#     define ULOX bc(1,1,1)
#     define UHIX bc(1,2,1)

#     define VLOX bc(1,1,2)
#     define VHIX bc(1,2,2)

#     define ULOY bc(2,1,1)
#     define UHIY bc(2,2,1)

#     define VLOY bc(2,1,2)
#     define VHIY bc(2,2,2)

c
c     ::::: statement functions that implement stencil
c
      uxcen(i,j) = half*(U(i+1,j)-U(i-1,j))/dx
      uxlft(i,j) = (U(i+1,j)+three*U(i,j)-four*U(i-1,j))/(three*dx)
      uxrgt(i,j) = (four*U(i+1,j)-three*U(i,j)-U(i-1,j))/(three*dx)
      
      vxcen(i,j) = half*(V(i+1,j)-V(i-1,j))/dx
      vxlft(i,j) = (V(i+1,j)+three*V(i,j)-four*V(i-1,j))/(three*dx)
      vxrgt(i,j) = (four*V(i+1,j)-three*V(i,j)-V(i-1,j))/(three*dx)

      uycen(i,j) = half*(U(i,j+1)-U(i,j-1))/dy
      uybot(i,j) = (U(i,j+1)+three*U(i,j)-four*U(i,j-1))/(three*dy)
      uytop(i,j) = (four*U(i,j+1)-three*U(i,j)-U(i,j-1))/(three*dy)

      vycen(i,j) = half*(V(i,j+1)-V(i,j-1))/dy
      vybot(i,j) = (V(i,j+1)+three*V(i,j)-four*V(i,j-1))/(three*dy)
      vytop(i,j) = (four*V(i,j+1)-three*V(i,j)-V(i,j-1))/(three*dy)

c
c     ::::: function for computing magnitude of rate of strain tensor
c
      strnrt_fun(ux,vx,uy,vy) = sqrt(ux**2+half*(vx+uy)**2+vy**2)

c
c     ::::: Functions for computing the regularised stress for different
c     ::::: power law behaviours 
c
      ! Pseudoplastic (shear-thinning)
      ! Regularise power-law term and yield stress term
      stress_pseudoplastic(strnrt) = 
     &      (mu*two**((n+one)/two)*strnrt**n+tau)*(one-exp(-strnrt/eps))
      ! Bingham fluid
      ! Regularise yield stress term only
      stress_bingham(strnrt) = 
     &      mu*two*strnrt+tau*(one-exp(-strnrt/eps))
      ! Dilatant (shear-thickening)
      ! Regularise yield stress term only
      stress_dilatant(strnrt) = 
     &      mu*two**((n+one)/two)*strnrt**n+tau*(one-exp(-strnrt/eps))


      dx = delta(1)
      dy = delta(2)

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            ux  = uxcen(i,j)
            vx  = vxcen(i,j)
            uy  = uycen(i,j)
            vy  = vycen(i,j)
            strnrt = strnrt_fun(ux,vx,uy,vy)
            if (varvisc .eq. 0) then
               ! If non-variable, stress equals 2 mu strnrt everywhere
               stress(i,j,1) = 2*mu*strnrt
            else
               ! Calculate viscosity based on flow index value
               if (n .lt. one) then
                  stress(i,j,1) = stress_pseudoplastic(strnrt)
               else if (n .eq. one) then
                  stress(i,j,1) = stress_bingham(strnrt)
               else 
                  stress(i,j,1) = stress_dilatant(strnrt)
               end if
            end if
         end do
      end do

      fixlft = ( (lo(1) .eq. domlo(1)) .and. 
     &           ( (ULOX .eq. EXT_DIR .or. ULOX .eq. HOEXTRAP) .or. 
     &             (VLOX .eq. EXT_DIR .or. VLOX .eq. HOEXTRAP) ) )
      fixrgt = ( (hi(1) .eq. domhi(1)) .and. 
     &           ( (UHIX .eq. EXT_DIR .or. UHIX .eq. HOEXTRAP) .or.
     &             (VHIX .eq. EXT_DIR .or. VHIX .eq. HOEXTRAP) ) )
      fixbot = ( (lo(2) .eq. domlo(2)) .and. 
     &           ( (ULOY .eq. EXT_DIR .or. ULOY .eq. HOEXTRAP) .or.
     &           (VLOY .eq. EXT_DIR .or. VLOY .eq. HOEXTRAP) ) )
      fixtop = ( (hi(2) .eq. domhi(2)) .and. 
     &           ( (UHIY .eq. EXT_DIR .or. UHIY .eq. HOEXTRAP) .or.
     &           (VHIY .eq. EXT_DIR .or. VHIY .eq. HOEXTRAP) ) )
c
c     ::::: handle special bndry conditions at the left edge
c
      if (fixlft) then
         i = lo(1)
         do j = lo(2), hi(2)
            ux  = uxlft(i,j)
            vx  = vxlft(i,j)
            uy  = uycen(i,j)
            vy  = vycen(i,j)
            strnrt = strnrt_fun(ux,vx,uy,vy)
            if (varvisc .eq. 0) then
               stress(i,j,1) = 2*mu*strnrt
            else
               if (n .lt. one) then
                  stress(i,j,1) = stress_pseudoplastic(strnrt)
               else if (n .eq. one) then
                  stress(i,j,1) = stress_bingham(strnrt)
               else 
                  stress(i,j,1) = stress_dilatant(strnrt)
               end if
            end if
         end do
      end if
c
c     ::::: handle special bndry conditions on the right
c
      if (fixrgt) then
         i = hi(1)
         do j = lo(2), hi(2)
            ux  = uxrgt(i,j)
            vx  = vxrgt(i,j)
            uy  = uycen(i,j)
            vy  = vycen(i,j)
            strnrt = strnrt_fun(ux,vx,uy,vy)
            if (varvisc .eq. 0) then
               stress(i,j,1) = 2*mu*strnrt
            else
               if (n .lt. one) then
                  stress(i,j,1) = stress_pseudoplastic(strnrt)
               else if (n .eq. one) then
                  stress(i,j,1) = stress_bingham(strnrt)
               else 
                  stress(i,j,1) = stress_dilatant(strnrt)
               end if
            end if
         end do
      end if
c
c     ::::: handle special bndry conditions on bottom
c
      if (fixbot) then
         j = lo(2)
         do i = lo(1), hi(1)
            ux  = uxcen(i,j)
            vx  = vxcen(i,j)
            uy  = uybot(i,j)
            vy  = vybot(i,j)
            strnrt = strnrt_fun(ux,vx,uy,vy)
            if (varvisc .eq. 0) then
               stress(i,j,1) = 2*mu*strnrt
            else
               if (n .lt. one) then
                  stress(i,j,1) = stress_pseudoplastic(strnrt)
               else if (n .eq. one) then
                  stress(i,j,1) = stress_bingham(strnrt)
               else 
                  stress(i,j,1) = stress_dilatant(strnrt)
               end if
            end if
         end do
      end if
c
c     ::::: handle special bndry conditions on top
c
      if (fixtop) then
         j = hi(2)
         do i = lo(1), hi(1)
            ux  = uxcen(i,j)
            vx  = vxcen(i,j)
            uy  = uytop(i,j)
            vy  = vytop(i,j)
            strnrt = strnrt_fun(ux,vx,uy,vy)
            if (varvisc .eq. 0) then
               stress(i,j,1) = 2*mu*strnrt
            else
               if (n .lt. one) then
                  stress(i,j,1) = stress_pseudoplastic(strnrt)
               else if (n .eq. one) then
                  stress(i,j,1) = stress_bingham(strnrt)
               else 
                  stress(i,j,1) = stress_dilatant(strnrt)
               end if
            end if
         end do
      end if
c
c     ::::: check corners
c
      if (fixlft .and. fixbot) then
         i = lo(1)
         j = lo(2)
         ux = uxlft(i,j)
         vx = vxlft(i,j)
         uy = uybot(i,j)
         vy = vybot(i,j)
         strnrt = strnrt_fun(ux,vx,uy,vy)
         if (varvisc .eq. 0) then
            stress(i,j,1) = 2*mu*strnrt
         else
            if (n .lt. one) then
               stress(i,j,1) = stress_pseudoplastic(strnrt)
            else if (n .eq. one) then
               stress(i,j,1) = stress_bingham(strnrt)
            else 
               stress(i,j,1) = stress_dilatant(strnrt)
            end if
         end if
      end if
      if (fixlft .and. fixtop) then
         i = lo(1)
         j = hi(2)
         ux = uxlft(i,j)
         vx = vxlft(i,j)
         uy = uytop(i,j)
         vy = vytop(i,j)
         strnrt = strnrt_fun(ux,vx,uy,vy)
         if (varvisc .eq. 0) then
            stress(i,j,1) = 2*mu*strnrt
         else
            if (n .lt. one) then
               stress(i,j,1) = stress_pseudoplastic(strnrt)
            else if (n .eq. one) then
               stress(i,j,1) = stress_bingham(strnrt)
            else 
               stress(i,j,1) = stress_dilatant(strnrt)
            end if
         end if
      end if
      if (fixrgt .and. fixtop) then
         i = hi(1)
         j = hi(2)
         ux = uxrgt(i,j)
         vx = vxrgt(i,j)
         uy = uytop(i,j)
         vy = vytop(i,j)
         strnrt = strnrt_fun(ux,vx,uy,vy)
         if (varvisc .eq. 0) then
            stress(i,j,1) = 2*mu*strnrt
         else
            if (n .lt. one) then
               stress(i,j,1) = stress_pseudoplastic(strnrt)
            else if (n .eq. one) then
               stress(i,j,1) = stress_bingham(strnrt)
            else 
               stress(i,j,1) = stress_dilatant(strnrt)
            end if
         end if
      end if
      if (fixrgt .and. fixbot) then
         i = hi(1)
         j = lo(2)
         ux = uxrgt(i,j)
         vx = vxrgt(i,j)
         uy = uybot(i,j)
         vy = vybot(i,j)
         strnrt = strnrt_fun(ux,vx,uy,vy)
         if (varvisc .eq. 0) then
            stress(i,j,1) = 2*mu*strnrt
         else
            if (n .lt. one) then
               stress(i,j,1) = stress_pseudoplastic(strnrt)
            else if (n .eq. one) then
               stress(i,j,1) = stress_bingham(strnrt)
            else 
               stress(i,j,1) = stress_dilatant(strnrt)
            end if
         end if
      end if

#     undef U
#     undef V      
#     undef VLOX
#     undef VHIX
#     undef ULOY
#     undef UHIY

      end

      subroutine FORT_DERMGDIVU (divu,DIMS(divu),nv,dat,DIMS(dat),ncomp,
     &                           lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                           bc,level,grid_no)
      implicit none
c
c ::: This routine will derive magnitude of the divergence of velocity
c
      integer    lo(2), hi(2)
      integer    DIMDEC(divu)
      integer    DIMDEC(dat)
      integer    domlo(2), domhi(2)
      integer    nv, ncomp
      integer    bc(2,2,ncomp)
      REAL_T     delta(2), xlo(2), time, dt
      REAL_T     divu(DIMV(divu),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      integer   i,j
      REAL_T    ux, vy, dx, dy
      REAL_T    uxcen, uxlo, uxhi
      REAL_T    vycen, vylo, vyhi
c
c     ::::: some useful macro definitions
c
#     define U(i,j) dat(i,j,1)
#     define V(i,j) dat(i,j,2)

#     define ULOX bc(1,1,1)
#     define UHIX bc(1,2,1)

#     define VLOY bc(2,1,2)
#     define VHIY bc(2,2,2)
c
c     ::::: statement functions that implement stencil
c
      uxcen(i,j) = half*(U(i+1,j)-U(i-1,j))/dx
      uxlo(i,j) = (eight*U(i,j)-six*U(i+1,j)+U(i+2,j))/(three*dx)
      uxhi(i,j) = (eight*U(i,j)-six*U(i-1,j)+U(i-2,j))/(three*dx)

      vycen(i,j) = half*(V(i,j+1)-V(i,j-1))/dy
      vylo(i,j) = (eight*V(i,j)-six*V(i,j+1)+V(i,j+2))/(three*dy)
      vyhi(i,j) = (eight*V(i,j)-six*V(i,j-1)+V(i,j-2))/(three*dy)

      call FORT_XVELFILL(dat(ARG_L1(dat),ARG_L2(dat),1),DIMS(dat),
     $                   domlo,domhi,delta,xlo,time,bc(1,1,1))
      call FORT_YVELFILL(dat(ARG_L1(dat),ARG_L2(dat),2),DIMS(dat),
     $                   domlo,domhi,delta,xlo,time,bc(1,1,2))

      dx = delta(1)
      dy = delta(2)
c
c     :: at physical bndries where an edge value is prescribed,
c     :: set the value in the outside cell so that a central
c     :: difference formula is equivalent to the higher order
c     :: one sided formula
c
c     boundaries handled by fill routines
c
      if (lo(1) .eq. domlo(1)) then
         i = lo(1)
         if (ULOX.eq.EXT_DIR) then
           do j = lo(2), hi(2)
              U(i-1,j) = two*U(i-1,j)-U(i,j)
           end do
         else if (ULOX.eq.HOEXTRAP) then
           do j = lo(2), hi(2)
              U(i-1,j) = uxlo(i,j)
           end do
         end if
      end if
      if (hi(1) .eq. domhi(1)) then
         i = hi(1)
         if (UHIX.eq.EXT_DIR) then
           do j = lo(2), hi(2)
              U(i+1,j) = two*U(i+1,j)-U(i,j)
           end do
         else if (UHIX.eq.HOEXTRAP) then
           do j = lo(2), hi(2)
              U(i+1,j) = uxhi(i,j)
           end do
         end if
      end if
      if (lo(2) .eq. domlo(2)) then
         j = lo(2)
         if (VLOY.eq.EXT_DIR) then
           do i = lo(1), hi(1)
              V(i,j-1) = two*V(i,j-1)-V(i,j)
           end do
         else if (VLOY.eq.HOEXTRAP) then
           do i = lo(1), hi(1)
              V(i,j-1) = vylo(i,j)
           end do
         end if
      end if
      if (hi(2) .eq. domhi(2)) then
         j = hi(2)
         if (VHIY.eq.EXT_DIR) then
           do i = lo(1), hi(1)
              V(i,j+1) = two*V(i,j+1)-V(i,j)
           end do
         else if (VHIY.eq.HOEXTRAP) then
           do i = lo(1), hi(1)
              V(i,j+1) = vyhi(i,j)
           end do
         end if
      end if

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            ux = uxcen(i,j)
            vy = vycen(i,j)
            divu(i,j,1) = ux + vy
         end do
      end do
c
c we overwrote the ghost cells above, so set them back below
c
      call FORT_XVELFILL(dat(ARG_L1(dat),ARG_L2(dat),1),DIMS(dat),
     $                   domlo,domhi,delta,xlo,time,bc(1,1,1))
      call FORT_YVELFILL(dat(ARG_L1(dat),ARG_L2(dat),2),DIMS(dat),
     $                   domlo,domhi,delta,xlo,time,bc(1,1,2))

#     undef U
#     undef V      
#     undef ULOX
#     undef UHIX
#     undef VLOY
#     undef VHIY

      end

      subroutine FORT_GRADP_DIR (
     &     p,DIMS(p),
     &     gp,DIMS(gp),
     &     lo,hi,dir,dx)
      implicit none
c
c     compute a node centered pressure gradient in direction (dir)
c
      integer    DIMDEC(p)
      integer    DIMDEC(gp)
      integer    lo(SDIM),  hi(SDIM)
      integer    dir
      REAL_T     dx
      REAL_T     p(DIMV(p))
      REAL_T     gp(DIMV(gp))
      integer    i,j
      REAL_T     d

#if (BL_PRVERSION == 9)
      d = half
#else
      d = half/dx
#endif

      if (dir .eq. 0) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               gp(i,j) = d*(p(i+1,j)-p(i,j)+p(i+1,j+1)-p(i,j+1))
            end do
         end do
      else if (dir .eq. 1) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               gp(i,j) = d*(p(i,j+1)-p(i,j)+p(i+1,j+1)-p(i+1,j))
            end do
         end do
      else
         call bl_abort("FORT_GRADP_DIR: invalid dir")
      end if

      end

      subroutine FORT_DERAVGPRES (avgpres,DIMS(gp),nv,dat,DIMS(dat),ncomp,
     &                            lo,hi,domlo,domhi,delta,xlo,time,dt,
     &                            bc,level,grid_no)
      implicit none
c
c     This routine computes cell-centered pressure as average of the four
c       surrounding nodal values.
c
      integer DIMDEC(gp)
      integer DIMDEC(dat)
      REAL_T  avgpres(DIMV(gp))
      REAL_T  dat(DIMV(dat))
      integer nv, ncomp
      integer lo(2), hi(2)
      integer domlo(2), domhi(2)
      REAL_T  delta(2)
      REAL_T  xlo(2)
      REAL_T  time, dt
      integer bc(2,2,ncomp)
      integer level
      integer grid_no

      integer i,j

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            avgpres(i,j) = fourth*( dat(i+1,j  )+dat(i,j  ) +
     $                              dat(i+1,j+1)+dat(i,j+1) )
         end do
      end do

      end

c=========================================================


      subroutine FORT_DERGRDPX (grdpx,DIMS(gp),nv,dat,DIMS(dat),ncomp,
     &                          lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                          bc,level,grid_no)
      implicit none
c
c     This routine computes pressure gradient in X direction
c
      integer lo(2), hi(2)
      integer DIMDEC(gp)
      integer DIMDEC(dat)
      integer domlo(2), domhi(2)
      integer nv, ncomp
      integer bc(2,2,ncomp)
      REAL_T  delta(2), xlo(2), time, dt
      REAL_T  grdpx(DIMV(gp),nv)
      REAL_T  dat(DIMV(dat),ncomp)
      integer level, grid_no
c
      call FORT_GRADP_DIR (
     &     dat,DIMS(dat),grdpx,DIMS(gp),
     &     lo,hi,0,delta(1))

      end

      subroutine FORT_DERGRDPY (grdpy,DIMS(gp),nv,dat,DIMS(dat),ncomp,
     &                          lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                          bc,level,grid_no)
      implicit none
c
c     This routine computes pressure gradient in Y direction
c
      integer lo(2), hi(2)
      integer DIMDEC(gp)
      integer DIMDEC(dat)
      integer domlo(2), domhi(2)
      integer nv, ncomp
      integer bc(2,2,ncomp)
      REAL_T  delta(2), xlo(2), time, dt
      REAL_T  grdpy(DIMV(gp),nv)
      REAL_T  dat(DIMV(dat),ncomp)
      integer level, grid_no
c
      call FORT_GRADP_DIR (
     &     dat,DIMS(dat),grdpy,DIMS(gp),
     &     lo,hi,1,delta(2))

      end

c=========================================================

      subroutine FORT_DERGRDP (grdp,DIMS(gp),nv,p,DIMS(p),ncomp,
     &                         lo,hi,domlo,domhi,dx,xlo,time,dt,
     $                         bc,level,grid_no)
      implicit none
c
c     This routine computes the magnitude of pressure gradient 
c
      integer lo(2), hi(2)
      integer DIMDEC(gp)
      integer DIMDEC(p)
      integer domlo(2), domhi(2)
      integer nv, ncomp
      integer bc(2,2,ncomp)
      REAL_T  dx(2), xlo(2), time, dt
      REAL_T  grdp(DIMV(gp),nv)
      REAL_T  p(DIMV(p),ncomp)
      integer level, grid_no

      REAL_T     gpx, gpy
      integer    i,j
c
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            gpx = half * (p(i+1,j,1)-p(i,j,1)+p(i+1,j+1,1)-p(i,j+1,1)) / dx(1)
            gpy = half * (p(i,j+1,1)-p(i,j,1)+p(i+1,j+1,1)-p(i+1,j,1)) / dx(2)
            grdp(i,j,1) = sqrt(gpx**2 + gpy**2)      
         end do
      end do

      end

c=========================================================

      subroutine FORT_DERNULL (e,DIMS(e),nv,dat,DIMS(dat),ncomp,
     &                         lo,hi,domlo,domhi,delta,xlo,time,dt,bc,
     $                         level,grid_no)
      implicit none
      !
      ! This is a null derived routine.
      !
      integer    lo(2), hi(2)
      integer    DIMDEC(e)
      integer    DIMDEC(dat)
      integer    domlo(2), domhi(2)
      integer    nv, ncomp
      integer    bc(2,2,ncomp)
      REAL_T     delta(2), xlo(2), time, dt
      REAL_T     e(DIMV(e),nv)
      REAL_T     dat(DIMV(dat),ncomp)
      integer    level, grid_no

      end

