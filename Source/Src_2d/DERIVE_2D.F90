
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


module derive_2d_module
  
  implicit none

  private

  public dermodgradrho, derkeng,derlogs, dermvel, &
       derlgrhodust, &
       dermgvort2, &
       dergrdp, dernull

contains

! c     -----------------------------------------------------------
! c     This file contains functions which compute derived quantities.  
! c     All of the argument lists have the same template, shown below
! c     
! c     INPUTS/OUTPUTS:
! c     
! c     e         <= the quantity derived
! c     DIMS(e)   => index extent of e array
! c     nv        => number of components in e array (should be 1)
! c     dat       => data neded to derive e
! c     DIMS(dat) => index limits of dat array
! c     ncomp     => number of components of dat array (3)
! c     lo,hi     => subrange of e array where result is requested
! c     domlo,hi  => index extent of problem domain (cell centered)
! c     delta     => cell spacing
! c     xlo       => physical location of lower left hand
! c 	           corner of e array
! c     time      => problem evolution time
! c     bc        => array of bndry types for component values
! c                  valid only if component touches bndry
! c     -----------------------------------------------------------

      subroutine dermodgradrho (e,DIMS(e),nv,dat,DIMS(dat),ncomp,&
                              lo,hi,domlo,domhi,delta,xlo,time,dt,&
                              bc,level,grid_no) bind(C,name="dermodgradrho")
      implicit none
!c
!c     This routine will derive the modulus of the density gradients
!c
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

    end subroutine dermodgradrho

      subroutine derkeng (e,DIMS(e),nv,dat,DIMS(dat),ncomp,&
                              lo,hi,domlo,domhi,delta,xlo,time,dt,&
                              bc,level,grid_no) bind(C,name="derkeng")
      implicit none
!c
!c     This routine will derive kinetic energy from density
!c     and the velocity field.
!c
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

    end subroutine derkeng

      subroutine derlogs (e,DIMS(e),nv,dat,DIMS(dat),ncomp,&
                              lo,hi,domlo,domhi,delta,xlo,time,dt,&
                              bc,level,grid_no) bind(C,name="derlogs")
      implicit none
!c
!c     This routine will derive log of given scalar quantity
!c
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

    end subroutine derlogs

      subroutine dermvel (e,DIMS(e),nv,dat,DIMS(dat),ncomp,&
                              lo,hi,domlo,domhi,delta,xlo,time,dt,&
                              bc,level,grid_no) bind(C,name="dermvel")
      implicit none
!c
!c ::: This routine will derive the magnitude of the velocity field
!c ::: from the velocity field
!c

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

    end subroutine dermvel

      subroutine derlgrhodust (e,DIMS(e),nv,dat,DIMS(dat),ncomp,&
                                   lo,hi,domlo,domhi,delta,xlo,time,dt,&
                                   bc,level,grid_no) bind(C,name="derlgrhodust")
      implicit none
!c
!c ::: This routine will derive log(RHO*C)
!c
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

    end subroutine derlgrhodust

      subroutine dermgvort2 (vort,DIMS(vort),nv,dat,DIMS(dat),ncomp,&
                                 lo,hi,domlo,domhi,delta,xlo,time,dt,&
                                 bc,level,grid_no) bind(C,name="dermgvort2")
      implicit none
!c
!c ::: This routine will derive magnitude of vorticity from
!c ::: the velocity field
!c
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
!c
!c     ::::: some useful macro definitions
!c
#     define U(i,j) dat(i,j,1)
#     define V(i,j) dat(i,j,2)

#     define ULOY bc(2,1,1)
#     define UHIY bc(2,2,1)

#     define VLOX bc(1,1,2)
#     define VHIX bc(1,2,2)
!c
!c     ::::: statement functions that implement stencil
!c
      uycen(i,j) = half*(U(i,j+1)-U(i,j-1))/dy
      uylo(i,j) = (eight*U(i,j)-six*U(i,j+1)+U(i,j+2))/(three*dy)
      uyhi(i,j) = (eight*U(i,j)-six*U(i,j-1)+U(i,j-2))/(three*dy)

      vxcen(i,j) = half*(V(i+1,j)-V(i-1,j))/dx
      vxlo(i,j) = (eight*V(i,j)-six*V(i+1,j)+V(i+2,j))/(three*dx)
      vxhi(i,j) = (eight*V(i,j)-six*V(i-1,j)+V(i-2,j))/(three*dx)

      vorfun(vx,uy) = vx - uy

      dx = delta(1)
      dy = delta(2)
!c
!c     :: at physical bndries where an edge value is prescribed,
!c     :: set the value in the outside cell so that a central
!c     :: difference formula is equivalent to the higher order
!c     :: one sided formula
!c
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

    end subroutine dermgvort2

!c=========================================================

      subroutine dergrdp (grdp,DIMS(gp),nv,p,DIMS(p),ncomp,&
                              lo,hi,domlo,domhi,dx,xlo,time,dt,&
                              bc,level,grid_no) bind(C,name="dergrdp")
      implicit none

!c     This routine computes the magnitude of pressure gradient 

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
!c
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            gpx = half * (p(i+1,j,1)-p(i,j,1)+p(i+1,j+1,1)-p(i,j+1,1)) / dx(1)
            gpy = half * (p(i,j+1,1)-p(i,j,1)+p(i+1,j+1,1)-p(i+1,j,1)) / dx(2)
            grdp(i,j,1) = sqrt(gpx**2 + gpy**2)      
         end do
      end do

    end subroutine dergrdp

!c=========================================================

      subroutine dernull (e,DIMS(e),nv,dat,DIMS(dat),ncomp,&
                              lo,hi,domlo,domhi,delta,xlo,time,dt,bc,&
                              level,grid_no) bind(C,name="dernull")
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

    end subroutine dernull

  end module derive_2d_module
