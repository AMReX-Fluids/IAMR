#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H" 
#include "CONSTANTS.H" 
#include "BC_TYPES.H" 
#include "DERIVE_F.H" 
#include "PROB_PM_F.H"
#include "ArrayLim.H"

#define SDIM 2

      subroutine FORT_DERGRDPX (grdpx,DIMS(gp),nv,dat,DIMS(dat),ncomp,
     $                          lo,hi,domlo,domhi,delta,xlo,time,dt,
     $                          bc,level,grid_no)
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

      subroutine FORT_GRADP_DIR (
     &     p,DIMS(p),
     &     gp,DIMS(gp),
     &     lo,hi,dir,dx)
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
