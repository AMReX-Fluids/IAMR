
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
       derlgrhodust, dergrdp, dernull

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
