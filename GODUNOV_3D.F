c
c $Id: GODUNOV_3D.F,v 1.50 2008-12-08 23:39:32 ajnonaka Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "GODUNOV_F.H"
#include "ArrayLim.H"

#define SDIM 3
#define XVEL 1
#define YVEL 2
#define ZVEL 3
#define RHO  4

#define ALL  999

      subroutine FORT_ESTDT (
     &     vel,DIMS(vel),
     &     tforces,DIMS(tf),
     &     rho,DIMS(rho),
     &     lo,hi,dt,dx,cfl,u_max)
c 
c     ----------------------------------------------------------
c     estimate the timestep for this grid and scale by CFL number
c     This routine sets dt as dt = dt_est*cfl where
c     dt_est is estimated from the actual velocities and their 
c     total forcing
c     ----------------------------------------------------------
c
      implicit none
      integer i, j, k
      REAL_T  u, v, w
      REAL_T  small
      REAL_T  dt_start
      REAL_T  tforce1,tforce2,tforce3
      integer lo(SDIM), hi(SDIM)
      REAL_T  dt,dx(SDIM),cfl,u_max(SDIM)

      integer DIMDEC(vel)
      integer DIMDEC(rho)
      integer DIMDEC(tf)

      REAL_T  vel(DIMV(vel),SDIM)
      REAL_T  rho(DIMV(rho))
      REAL_T  tforces(DIMV(tf),SDIM)

      small   = 1.0D-8
      u       = zero
      v       = zero
      w       = zero
      tforce1 = zero
      tforce2 = zero
      tforce3 = zero

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               u = max(u,abs(vel(i,j,k,1)))	
               v = max(v,abs(vel(i,j,k,2)))
               w = max(w,abs(vel(i,j,k,3)))
               tforce1 = max(tforce1,abs(tforces(i,j,k,1)/rho(i,j,k)))
               tforce2 = max(tforce2,abs(tforces(i,j,k,2)/rho(i,j,k)))
               tforce3 = max(tforce3,abs(tforces(i,j,k,3)/rho(i,j,k)))
            end do
         end do
      end do

      u_max(1) = u
      u_max(2) = v
      u_max(3) = w

      dt_start = 1.0D+20
      dt = dt_start

      if (u .gt. small) dt = min(dt,dx(1)/u)
      if (v .gt. small) dt = min(dt,dx(2)/v)
      if (w .gt. small) dt = min(dt,dx(3)/w)

      if (tforce1 .gt. small) then
         dt  = min(dt,sqrt(two*dx(1)/tforce1))
      end if

      if (tforce2 .gt. small) then
         dt  = min(dt,sqrt(two*dx(2)/tforce2))
      end if

      if (tforce3 .gt. small) then
         dt  = min(dt,sqrt(two*dx(3)/tforce3))
      end if

      if (dt .eq. dt_start) dt = min(dx(1),dx(2),dx(3))

      dt = dt*cfl

      end

      subroutine FORT_TEST_U_RHO(
     &     u,DIMS(u),
     &     v,DIMS(v),
     &     w,DIMS(w),
     &     rho,DIMS(rho),
     &     lo,hi,dt,dx,cflmax,u_max,verbose)
c
c     This subroutine computes the extrema of the density
c     and velocities at cell centers
c
      implicit none
      integer DIMDEC(u)
      integer DIMDEC(v)
      integer DIMDEC(w)
      integer DIMDEC(rho)
      integer imin, imax, jmin, jmax, kmin, kmax
      integer i, j, k
      integer lo(SDIM),hi(SDIM)
      REAL_T  dx(SDIM), u_max(SDIM),cflmax, dt
      REAL_T  hx, hy, hz
      REAL_T  umax, vmax, wmax, rhomax
      REAL_T  umin, vmin, wmin, rhomin
      REAL_T  u(DIMV(u))
      REAL_T  v(DIMV(v))
      REAL_T  w(DIMV(w))
      REAL_T  rho(DIMV(rho))
      integer verbose

      hx   = dx(1)
      hy   = dx(2)
      hz   = dx(3)
      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
      kmin = lo(3)
      kmax = hi(3)
      umax = -1.d200
      vmax = -1.d200
      wmax = -1.d200
      umin =  1.d200
      vmin =  1.d200
      wmin =  1.d200
      rhomax = -1.d200
      rhomin =  1.d200

      do k = kmin, kmax
         do j = jmin, jmax
            do i = imin, imax
               umax = max(umax,u(i,j,k))
               umin = min(umin,u(i,j,k))
               vmax = max(vmax,v(i,j,k))
               vmin = min(vmin,v(i,j,k))
               wmax = max(wmax,w(i,j,k))
               wmin = min(wmin,w(i,j,k))
               rhomax = max(rhomax,rho(i,j,k))
               rhomin = min(rhomin,rho(i,j,k))
            end do
         end do
      end do

      u_max(1) = max(abs(umax), abs(umin))
      u_max(2) = max(abs(vmax), abs(vmin))
      u_max(3) = max(abs(wmax), abs(wmin))
      cflmax   = dt*max(u_max(1)/hx,u_max(2)/hy,u_max(3)/hz)

      if(verbose.eq.1)then
         write(6,1000) umax,umin,u_max(1)
         write(6,1001) vmax,vmin,u_max(2)
         write(6,1002) wmax,wmin,u_max(3)
         write(6,1003) rhomax,rhomin
#ifndef	BL_NO_FORT_FLUSH
c         call flush(6)
#endif
      end if

 1000 format(' U  MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1001 format(' V  MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1002 format(' W  MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1003 format('RHO MAX/MIN      ',e21.14,2x,e21.14)

      end

      subroutine FORT_TEST_UMAC_RHO(
     &     umac,DIMS(umac),
     &     vmac,DIMS(vmac),
     &     wmac,DIMS(wmac),
     &     rho,DIMS(rho),
     &     lo,hi,dt,dx,cflmax,u_max)
c
c     This subroutine computes the extrema of the density
c     and mac edge velocities
c
      implicit none
      integer lo(SDIM),hi(SDIM)
      REAL_T  dt, dx(SDIM), u_max(SDIM), cflmax
      integer imin, imax, jmin, jmax, kmin, kmax
      integer i, j, k
      REAL_T  hx, hy, hz
      REAL_T  umax, vmax, wmax, rhomax
      REAL_T  umin, vmin, wmin, rhomin

      integer DIMDEC(umac)
      integer DIMDEC(vmac)
      integer DIMDEC(wmac)
      integer DIMDEC(rho)

      REAL_T  umac(DIMV(umac))
      REAL_T  vmac(DIMV(vmac))
      REAL_T  wmac(DIMV(wmac))
      REAL_T  rho(DIMV(rho))

      hx   = dx(1)
      hy   = dx(2)
      hz   = dx(3)
      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
      kmin = lo(3)
      kmax = hi(3)
      umax = -1.d200
      vmax = -1.d200
      wmax = -1.d200
      umin =  1.d200
      vmin =  1.d200
      wmin =  1.d200
      rhomax = -1.d200
      rhomin =  1.d200

      do k = kmin, kmax
         do j = jmin, jmax
            do i = imin, imax+1
               umax = max(umax,umac(i,j,k))
               umin = min(umin,umac(i,j,k))
            end do
         end do
      end do

      do k = kmin, kmax
         do j = jmin, jmax+1
            do i = imin, imax
               vmax = max(vmax,vmac(i,j,k))
               vmin = min(vmin,vmac(i,j,k))
            end do
         end do
      end do

      do k = kmin, kmax+1
         do j = jmin, jmax
            do i = imin, imax
               wmax = max(wmax,wmac(i,j,k))
               wmin = min(wmin,wmac(i,j,k))
            end do
         end do
      end do

      do k = kmin, kmax
         do j = jmin, jmax
            do i = imin, imax
               rhomax = max(rhomax,rho(i,j,k))
               rhomin = min(rhomin,rho(i,j,k))
            end do
         end do
      end do

      u_max(1) = max(abs(umax), abs(umin))
      u_max(2) = max(abs(vmax), abs(vmin))
      u_max(3) = max(abs(wmax), abs(wmin))
      cflmax   = dt*max(u_max(1)/hx,u_max(2)/hy,u_max(3)/hz)

      write(6,1000) umax,umin,u_max(1)
      write(6,1001) vmax,vmin,u_max(2)
      write(6,1002) wmax,wmin,u_max(3)
      write(6,1003) rhomax,rhomin

 1000 format('UMAC MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1001 format('VMAC MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1002 format('WMAC MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1003 format('RHO  MAX/MIN      ',e21.14,2x,e21.14)

#ifndef	BL_NO_FORT_FLUSH
c      call flush(6)
#endif

      end

      subroutine FORT_TRANSVEL(
     &     u, ulo, uhi, sx, ubc, slxscr,
     &     v, vlo, vhi, sy, vbc, slyscr,
     &     w, wlo, whi, sz, wbc, slzscr,
     &     DIMS(s), DIMS(work),
     &     lo,hi,dt,dx,use_minion,tforces)
c
c     This subroutine computes the advective velocities used in
c     the transverse derivatives of the Godunov box
c
      implicit none
      integer i,j,k
      integer ubc(SDIM,2),vbc(SDIM,2),wbc(SDIM,2)
      integer lo(SDIM),hi(SDIM)
      integer imin,jmin,kmin,imax,jmax,kmax
      REAL_T hx, hy, hz, dth, dthx, dthy, dthz
      REAL_T dt, dx(SDIM)
      REAL_T eps, eps_for_bc, val, tst
      logical ltm
      parameter (eps        = 1.d-6 )
      parameter (eps_for_bc = 1.d-10)

      integer DIMDEC(s)
      integer DIMDEC(work)
      REAL_T  u(DIMV(s))
      REAL_T  v(DIMV(s))
      REAL_T  w(DIMV(s))
      REAL_T ulo(DIMV(work)),uhi(DIMV(work))
      REAL_T vlo(DIMV(work)),vhi(DIMV(work))
      REAL_T wlo(DIMV(work)),whi(DIMV(work))
      REAL_T sx(DIMV(work))
      REAL_T sy(DIMV(work))
      REAL_T sz(DIMV(work))
      REAL_T slxscr(DIM1(s), 4)
      REAL_T slyscr(DIM2(s), 4)
      REAL_T slzscr(DIM3(s), 4)

      integer use_minion
      REAL_T tforces(DIMV(work),SDIM)

      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)
      dthz = half*dt / dx(3)
      hx   = dx(1)
      hy   = dx(2)
      hz   = dx(3)
      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
c
c     =============== THE SCREWY ORDER is to maximize comparability
c     with the old fortran
c     --------------------------------------------------------------
c     compute the x transverse velocities
c     --------------------------------------------------------------
c     --------------------------------------------------------------
c     compute the y transverse velocities
c     --------------------------------------------------------------
c     --------------------------------------------------------------
c     compute the z transverse velocities
c     --------------------------------------------------------------
c
      call FORT_SLOPES( XVEL,
     &     u,DIMS(s),
     &     sx,sy,sz,DIMS(work),
     &     lo,hi,slxscr,slyscr,slzscr,ubc)

      do k = kmin-1,kmax+1
         do j = jmin-1,jmax+1
            do i = imin,  imax+1
               ulo(i,j,k) = u(i-1,j,k) + (half  - dthx*u(i-1,j,k))*sx(i-1,j,k)
               uhi(i,j,k) = u(i,  j,k) + (-half - dthx*u(i,  j,k))*sx(i,  j,k)
            end do
         end do
      end do

      if (use_minion .eq. 1 )then
         do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin,  imax+1
                  ulo(i,j,k) = ulo(i,j,k) + dth*tforces(i-1,j,k,1)
                  uhi(i,j,k) = uhi(i,j,k) + dth*tforces(i,  j,k,1)
               end do
            end do
         end do
      end if

      call FORT_SLOPES(YVEL,
     &     v,DIMS(s),
     &     sx,sy,sz,DIMS(work),
     &     lo,hi,slxscr,slyscr,slzscr,vbc)

      do k = kmin-1,kmax+1
         do j = jmin,  jmax+1
            do i = imin-1,imax+1
               vlo(i,j,k) = v(i,j-1,k) + (half  - dthy*v(i,j-1,k))*sy(i,j-1,k)
               vhi(i,j,k) = v(i,j,  k) + (-half - dthy*v(i,j,  k))*sy(i,j,  k)
            end do
         end do
      end do

      if (use_minion .eq. 1 )then
         do k = kmin-1,kmax+1
            do j = jmin,    jmax+1
               do i = imin-1,  imax+1
                  vlo(i,j,k) = vlo(i,j,k) + dth*tforces(i,j-1,k,2)
                  vhi(i,j,k) = vhi(i,j,k) + dth*tforces(i,j,  k,2)
               end do
            end do
         end do
      end if

      call FORT_SLOPES(ZVEL,
     &     w,DIMS(s),
     &     sx,sy,sz,DIMS(work),
     &     lo,hi,slxscr,slyscr,slzscr,wbc)

      do k = kmin,kmax+1
         do j = jmin-1,jmax+1
            do i = imin-1,imax+1
               wlo(i,j,k) = w(i,j,k-1) + (half  - dthz*w(i,j,k-1))*sz(i,j,k-1)
               whi(i,j,k) = w(i,j,k  ) + (-half - dthz*w(i,j,k  ))*sz(i,j,k  )
            end do
         end do
      end do

      if (use_minion .eq. 1 )then
         do k = kmin,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,  imax+1
                  wlo(i,j,k) = wlo(i,j,k) + dth*tforces(i,j,k-1,3)
                  whi(i,j,k) = whi(i,j,k) + dth*tforces(i,j,k,  3)
               end do
            end do
         end do
      end if

      call trans_xbc(
     &     u,DIMS(s),
     &     ulo,uhi,DIMS(work),ulo,DIMS(work),
     &     lo,hi,XVEL,ubc,eps_for_bc,.false.,.false.)

      call trans_ybc(
     &     v,DIMS(s),
     &     vlo,vhi,DIMS(work),vlo,DIMS(work),
     &     lo,hi,YVEL,vbc,eps_for_bc,.false.,.false.)

      call trans_zbc(
     &     w,DIMS(s),
     &     wlo,whi,DIMS(work),wlo,DIMS(work),
     &     lo,hi,ZVEL,wbc,eps_for_bc,.false.,.false.)

      do k = kmin-1,kmax+1
         do j = jmin-1,jmax+1
            do i = imin,  imax+1
               tst = ulo(i,j,k)+uhi(i,j,k)
               val = cvmgp(ulo(i,j,k),uhi(i,j,k),tst)
               ltm =
     &              ( (ulo(i,j,k) .le. zero) .and.
     &              (uhi(i,j,k) .ge. zero) ) .or.
     &              (abs(tst)   .lt. eps )
               ulo(i,j,k) = cvmgt(zero,val,ltm)
            end do
         end do
      end do

      do k = kmin-1,kmax+1
         do j = jmin,  jmax+1
            do i = imin-1,imax+1
               tst = vlo(i,j,k)+vhi(i,j,k)
               val = cvmgp(vlo(i,j,k),vhi(i,j,k),tst)
               ltm =
     &              ( (vlo(i,j,k) .le. zero) .and.
     &              (vhi(i,j,k) .ge. zero) ) .or.
     &              (abs(tst)   .lt. eps )
               vlo(i,j,k) = cvmgt(zero,val,ltm)
            end do
         end do
      end do
      
      do k = kmin,kmax+1
         do j = jmin-1,jmax+1
            do i = imin-1,imax+1
               tst = wlo(i,j,k)+whi(i,j,k)
               val = cvmgp(wlo(i,j,k),whi(i,j,k),tst)
               ltm =
     &              ( (wlo(i,j,k) .le. zero) .and.
     &              (whi(i,j,k) .ge. zero) ) .or.
     &              (abs(tst)   .lt. eps )
               wlo(i,j,k) = cvmgt(zero,val,ltm)
            end do
         end do
      end do

      end

      subroutine FORT_ESTATE(s, tforces, DIMS(s),

     &     u, xlo, xhi, sx, uad, slxscr, stxlo, stxhi,
     &     uedge, DIMS(uedge), xstate, DIMS(xstate),

     &     v, ylo, yhi, sy, vad, slyscr, stylo, styhi,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate),

     &     w, zlo, zhi, sz, wad, slzscr, stzlo, stzhi,
     &     wedge, DIMS(wedge), zstate, DIMS(zstate),

     &     xedge, yedge, zedge,
     &     xylo, xzlo, yxlo, yzlo, zxlo, zylo,
     &     xyhi, xzhi, yxhi, yzhi, zxhi, zyhi,

     &     DIMS(work),
     &     bc,lo,hi,dt,dx,n,velpred, use_minion)
c
c     This subroutine computes edges states, right now it uses
c     a lot of memory, but there becomes a trade off between
c     simplicity-efficiency in the new way of computing states
c     and complexity in the old way.  By eliminating loops over
c     state components though, the new way uses much less memory.
c
      implicit none
      integer i,j,k,n,velpred
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer imin,jmin,kmin,imax,jmax,kmax,inc
      REAL_T place_to_break
      REAL_T hx, hy, hz, dt, dth, dthx, dthy, dthz
      REAL_T dt3, dt3x, dt3y, dt3z, dt4, dt4x, dt4y, dt4z
      REAL_T dt6, dt6x, dt6y, dt6z
      REAL_T tr1,tr2,ubar,vbar,wbar,stx,sty,stz,fu,fv,fw,dx(SDIM)
      REAL_T eps,eps_for_bc
      logical ltx,lty,ltz
      parameter (eps        = 1.d-6 )
      parameter (eps_for_bc = 1.d-10)

      integer DIMDEC(s)
      integer DIMDEC(work)
      integer DIMDEC(uedge)
      integer DIMDEC(vedge)
      integer DIMDEC(wedge)
      integer DIMDEC(xstate)
      integer DIMDEC(ystate)
      integer DIMDEC(zstate)

      REAL_T s(DIMV(s))
      REAL_T u(DIMV(s))
      REAL_T v(DIMV(s))
      REAL_T w(DIMV(s))
      REAL_T stxlo(DIM1(s)),stxhi(DIM1(s)),slxscr(DIM1(s),4)
      REAL_T stylo(DIM2(s)),styhi(DIM2(s)),slyscr(DIM2(s),4)
      REAL_T stzlo(DIM3(s)),stzhi(DIM3(s)),slzscr(DIM3(s),4) 
      REAL_T uedge(DIMV(uedge)), xstate(DIMV(xstate))
      REAL_T vedge(DIMV(vedge)), ystate(DIMV(ystate))
      REAL_T wedge(DIMV(wedge)), zstate(DIMV(zstate))

      REAL_T xlo(DIMV(work)), xhi(DIMV(work)), xedge(DIMV(work))
      REAL_T ylo(DIMV(work)), yhi(DIMV(work)), yedge(DIMV(work))
      REAL_T zlo(DIMV(work)), zhi(DIMV(work)), zedge(DIMV(work))
      REAL_T  sx(DIMV(work)), uad(DIMV(work))
      REAL_T  sy(DIMV(work)), vad(DIMV(work))
      REAL_T  sz(DIMV(work)), wad(DIMV(work))
      REAL_T tforces(DIMV(work))

      logical corner_couple

c     used in 3d corner coupling
      REAL_T xylo(DIMV(work)), xyhi(DIMV(work))
      REAL_T xzlo(DIMV(work)), xzhi(DIMV(work))
      REAL_T yxlo(DIMV(work)), yxhi(DIMV(work))
      REAL_T yzlo(DIMV(work)), yzhi(DIMV(work))
      REAL_T zxlo(DIMV(work)), zxhi(DIMV(work))
      REAL_T zylo(DIMV(work)), zyhi(DIMV(work))

      integer use_minion

      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)
      dthz = half*dt / dx(3)
      dt3  = dt / 3.0d0
      dt3x = dt3 / dx(1)
      dt3y = dt3 / dx(2)
      dt3z = dt3 / dx(3)
      dt4  = dt / 4.0d0
      dt4x = dt4 / dx(1)
      dt4y = dt4 / dx(2)
      dt4z = dt4 / dx(3)
      dt6  = dt / 6.0d0
      dt6x = dt6 / dx(1)
      dt6y = dt6 / dx(2)
      dt6z = dt6 / dx(3)
      hx   = dx(1)
      hy   = dx(2)
      hz   = dx(3)
      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

      corner_couple = .false.

c
c     compute the slopes
c
      call FORT_SLOPES(ALL,
     &     s,DIMS(s),
     &     sx,sy,sz,DIMS(work),
     &     lo,hi,slxscr,slyscr,slzscr,bc)
c
c     trace the state to the cell edges
c
      do k = kmin-1,kmax+1
         do j = jmin-1,jmax+1
            do i = imin,  imax+1
               xlo(i,j,k) = s(i-1,j,k) + (half  - dthx*u(i-1,j,k))*sx(i-1,j,k)
               xhi(i,j,k) = s(i,  j,k) + (-half - dthx*u(i,  j,k))*sx(i,  j,k)
            end do
         end do
      end do

      if(use_minion.eq.1)then
         do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin,  imax+1
                  xlo(i,j,k) = xlo(i,j,k) + dth*tforces(i-1,j,k)
                  xhi(i,j,k) = xhi(i,j,k) + dth*tforces(i,  j,k)
               end do
            end do
         end do
      end if

      call trans_xbc(
     &     s,DIMS(s),
     &     xlo,xhi,DIMS(work),uad,DIMS(work),
     &     lo,hi,n,bc,eps_for_bc,.false.,.false.)

      do k = kmin-1,kmax+1
         do j = jmin-1,jmax+1
            do i = imin,  imax+1
               fu  = cvmgt(zero,one,abs(uad(i,j,k)).lt.eps)
               stx = cvmgp(xlo(i,j,k),xhi(i,j,k),uad(i,j,k))
               xedge(i,j,k) = fu*stx + (one - fu)*half*(xhi(i,j,k)+xlo(i,j,k))
            end do
         end do
      end do

      do k = kmin-1,kmax+1
         do j = jmin,  jmax+1
            do i = imin-1,imax+1
               ylo(i,j,k) = s(i,j-1,k) + (half  - dthy*v(i,j-1,k))*sy(i,j-1,k)
               yhi(i,j,k) = s(i,j, k)  + (-half - dthy*v(i,j,  k))*sy(i,j, k)
            end do
         end do
      end do

      if (use_minion.eq.1)then
         do k = kmin-1,kmax+1
            do j = jmin, jmax+1
               do i = imin-1,  imax+1
                  ylo(i,j,k) = ylo(i,j,k) + dth*tforces(i,j-1,k)
                  yhi(i,j,k) = yhi(i,j,k) + dth*tforces(i,j,  k)
               end do
            end do
         end do
      end if

      call trans_ybc(
     &     s,DIMS(s),
     &     ylo,yhi,DIMS(work),vad,DIMS(work),
     &     lo,hi,n,bc,eps_for_bc,.false.,.false.)

      do k = kmin-1,kmax+1
         do j = jmin,  jmax+1
            do i = imin-1,imax+1
               fv  = cvmgt(zero,one,abs(vad(i,j,k)).lt.eps)
               sty = cvmgp(ylo(i,j,k),yhi(i,j,k),vad(i,j,k))
               yedge(i,j,k) = fv*sty + (one - fv)*half*(yhi(i,j,k)+ylo(i,j,k))
            end do
         end do
      end do

      do k = kmin,kmax+1
         do j = jmin-1,jmax+1
            do i = imin-1,imax+1
               zlo(i,j,k) = s(i,j,k-1) + (half  - dthz*w(i,j,k-1))*sz(i,j,k-1)
               zhi(i,j,k) = s(i,j,k  ) + (-half - dthz*w(i,j,k  ))*sz(i,j,k  )
            end do
         end do
      end do

      if (use_minion.eq.1)then
         do k = kmin,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,  imax+1
                  zlo(i,j,k) = zlo(i,j,k) + dth*tforces(i,j,k-1)
                  zhi(i,j,k) = zhi(i,j,k) + dth*tforces(i,j,k)
               end do
            end do
         end do
      end if

      call trans_zbc(
     &     s,DIMS(s),
     &     zlo,zhi,DIMS(work),wad,DIMS(work),
     &     lo,hi,n,bc,eps_for_bc,.false.,.false.)

      do k = kmin,kmax+1
         do j = jmin-1,jmax+1
            do i = imin-1,imax+1
               fw  = cvmgt(zero,one,abs(wad(i,j,k)).lt.eps)
               stz = cvmgp(zlo(i,j,k),zhi(i,j,k),wad(i,j,k))
               zedge(i,j,k) = fw*stz + (one-fw)*half*(zhi(i,j,k)+zlo(i,j,k))
            end do
         end do
      end do

      if (corner_couple) then

c
c     NEW CORNER COUPLING CODE
c

c
c     compute the corner-coupled terms:
c     xylo/hi, xzlo/hi, yxlo/hi, yzlo/hi, zxlo/hi, zylo/hi
c

c     loop over appropriate xy faces
      do k=kmin-1,kmax+1
         do j=jmin,jmax
            do i=imin,imax+1
               xylo(i,j,k) = xlo(i,j,k)
     &              - dt6y*(vad(i-1,j+1,k)+vad(i-1,j,k))
     &              *(yedge(i-1,j+1,k)-yedge(i-1,j,k))
               xyhi(i,j,k) = xhi(i,j,k)
     &              - dt6y*(vad(i  ,j+1,k)+vad(i  ,j,k))
     &              *(yedge(i  ,j+1,k)-yedge(i  ,j,k))
            end do
         end do
      end do

c     boundary conditions
      call trans_xbc(
     &     s,DIMS(s),
     &     xylo,xyhi,DIMS(work),uad,DIMS(work),
     &     lo,hi,n,bc,eps_for_bc,.true.,.false.)

c     upwind
      do k=kmin-1,kmax+1
         do j=jmin,jmax
            do i=imin,imax+1
               fu  = cvmgt(zero,one,abs(uad(i,j,k)).lt.eps)
               stx = cvmgp(xylo(i,j,k),xyhi(i,j,k),uad(i,j,k))
               xylo(i,j,k) = fu*stx + (one - fu)*half*(xyhi(i,j,k)+xylo(i,j,k))
            end do
         end do
      end do

c     loop over appropriate xz faces
      do k=kmin,kmax
         do j=jmin-1,jmax+1
            do i=imin,imax+1
               xzlo(i,j,k) = xlo(i,j,k)
     &              - dt6z*(wad(i-1,j,k+1)+wad(i-1,j,k))
     &              *(zedge(i-1,j,k+1)-zedge(i-1,j,k))
               xzhi(i,j,k) = xhi(i,j,k)
     &              - dt6z*(wad(i  ,j,k+1)+wad(i  ,j,k))
     &              *(zedge(i  ,j,k+1)-zedge(i  ,j,k))
            end do
         end do
      end do

c     boundary conditions
      call trans_xbc(
     &     s,DIMS(s),
     &     xzlo,xzhi,DIMS(work),uad,DIMS(work),
     &     lo,hi,n,bc,eps_for_bc,.false.,.true.)

c     upwind
      do k=kmin,kmax
         do j=jmin-1,jmax+1
            do i=imin,imax+1
               fu  = cvmgt(zero,one,abs(uad(i,j,k)).lt.eps)
               stx = cvmgp(xzlo(i,j,k),xzhi(i,j,k),uad(i,j,k))
               xzlo(i,j,k) = fu*stx + (one - fu)*half*(xzhi(i,j,k)+xzlo(i,j,k))
            end do
         end do
      end do

c     loop over appropriate yx faces
      do k=kmin-1,kmax+1
         do j=jmin,jmax+1
            do i=imin,imax
               yxlo(i,j,k) = ylo(i,j,k)
     &              - dt6x*(uad(i+1,j-1,k)+uad(i,j-1,k))
     &              *(xedge(i+1,j-1,k)-xedge(i,j-1,k))
               yxhi(i,j,k) = yhi(i,j,k)
     &              - dt6x*(uad(i+1,j  ,k)+uad(i,j  ,k))
     &              *(xedge(i+1,j  ,k)-xedge(i,j  ,k))
            end do
         end do
      end do

c     boundary conditions
      call trans_ybc(
     &     s,DIMS(s),
     &     yxlo,yxhi,DIMS(work),vad,DIMS(work),
     &     lo,hi,n,bc,eps_for_bc,.true.,.false.)

c     upwind
      do k=kmin-1,kmax+1
         do j=jmin,jmax+1
            do i=imin,imax
               fv  = cvmgt(zero,one,abs(vad(i,j,k)).lt.eps)
               sty = cvmgp(yxlo(i,j,k),yxhi(i,j,k),vad(i,j,k))
               yxlo(i,j,k) = fv*sty + (one - fv)*half*(yxhi(i,j,k)+yxlo(i,j,k))
            end do
         end do
      end do

c     loop over appropriate yz faces
      do k=kmin,kmax
         do j=jmin,jmax+1
            do i=imin-1,imax+1
               yzlo(i,j,k) = ylo(i,j,k)
     &              - dt6z*(wad(i,j-1,k+1)+wad(i,j-1,k))
     &              *(zedge(i,j-1,k+1)-zedge(i,j-1,k))
               yzhi(i,j,k) = yhi(i,j,k)
     &              - dt6z*(wad(i,j  ,k+1)+wad(i,j  ,k))
     &              *(zedge(i,j  ,k+1)-zedge(i,j  ,k))
            end do
         end do
      end do

c     boundary conditions
      call trans_ybc(
     &     s,DIMS(s),
     &     yzlo,yzhi,DIMS(work),vad,DIMS(work),
     &     lo,hi,n,bc,eps_for_bc,.false.,.true.)

c     upwind
      do k=kmin,kmax
         do j=jmin,jmax+1
            do i=imin-1,imax+1
               fv  = cvmgt(zero,one,abs(vad(i,j,k)).lt.eps)
               sty = cvmgp(yzlo(i,j,k),yzhi(i,j,k),vad(i,j,k))
               yzlo(i,j,k) = fv*sty + (one - fv)*half*(yzhi(i,j,k)+yzlo(i,j,k))
            end do
         end do
      end do

c     loop over appropriate zx faces
      do k=kmin,kmax+1
         do j=jmin-1,jmax+1
            do i=imin,imax
               zxlo(i,j,k) = zlo(i,j,k)
     &              - dt6x*(uad(i+1,j,k-1)+uad(i,j,k-1))
     &              *(xedge(i+1,j,k-1)-xedge(i,j,k-1))
               zxhi(i,j,k) = zhi(i,j,k)
     &              - dt6x*(uad(i+1,j,k  )+uad(i,j,k  ))
     &              *(xedge(i+1,j,k  )-xedge(i,j,k  ))
            end do
         end do
      end do

c     boundary conditions
      call trans_zbc(
     &     s,DIMS(s),
     &     zxlo,zxhi,DIMS(work),wad,DIMS(work),
     &     lo,hi,n,bc,eps_for_bc,.true.,.false.)

c     upwind
      do k=kmin,kmax+1
         do j=jmin-1,jmax+1
            do i=imin,imax
               fw  = cvmgt(zero,one,abs(wad(i,j,k)).lt.eps)
               stz = cvmgp(zxlo(i,j,k),zxhi(i,j,k),wad(i,j,k))
               zxlo(i,j,k) = fw*stz + (one-fw)*half*(zxhi(i,j,k)+zxlo(i,j,k))
            end do
         end do
      end do

c     loop over appropriate zy faces
      do k=kmin,kmax+1
         do j=jmin,jmax
            do i=imin-1,imax+1
               zylo(i,j,k) = zlo(i,j,k)
     &              - dt6y*(vad(i,j+1,k-1)+vad(i,j,k-1))
     &              *(yedge(i,j+1,k-1)-yedge(i,j,k-1))
               zyhi(i,j,k) = zhi(i,j,k)
     &              - dt6y*(vad(i,j+1,k  )+vad(i,j,k  ))
     &              *(yedge(i,j+1,k  )-yedge(i,j,k  ))
            end do
         end do
      end do

c     boundary conditions
      call trans_zbc(
     &     s,DIMS(s),
     &     zylo,zyhi,DIMS(work),wad,DIMS(work),
     &     lo,hi,n,bc,eps_for_bc,.false.,.true.)

c     upwind
      do k=kmin,kmax+1
         do j=jmin,jmax
            do i=imin-1,imax+1
               fw  = cvmgt(zero,one,abs(wad(i,j,k)).lt.eps)
               stz = cvmgp(zylo(i,j,k),zyhi(i,j,k),wad(i,j,k))
               zylo(i,j,k) = fw*stz + (one-fw)*half*(zyhi(i,j,k)+zylo(i,j,k))
            end do
         end do
      end do

c
c     compute the xedge states
c
      if ((velpred.ne.1) .or. (n.eq.XVEL)) then
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax+1

                  stxlo(i) = xlo(i,j,k)
     &                 - dt4y*(vad(i-1,j+1,k  )+vad(i-1,j,k))*
     &                 (yzlo(i-1,j+1,k  )-yzlo(i-1,j,k))
     &                 - dt4z*(wad(i-1,j  ,k+1)+wad(i-1,j,k))*
     &                 (zylo(i-1,j  ,k+1)-zylo(i-1,j,k))
                  stxhi(i) = xhi(i,j,k)
     &                 - (dt4/hy)*(vad(i  ,j+1,k  )+vad(i  ,j,k))*
     &                 (yzlo(i  ,j+1,k  )-yzlo(i  ,j,k))
     &                 - (dt4/hz)*(wad(i  ,j  ,k+1)+wad(i  ,j,k))*
     &                 (zylo(i  ,j  ,k+1)-zylo(i  ,j,k))

                  if (use_minion.eq.0) then
                     stxlo(i) = stxlo(i) + dth*tforces(i-1,j,k)
                     stxhi(i) = stxhi(i) + dth*tforces(i,  j,k)
                  end if

               end do

               if (bc(1,1).eq.EXT_DIR .and. velpred.eq.1) then
                  stxhi(imin) = s(imin-1,j,k)
                  stxlo(imin) = s(imin-1,j,k)
               else if (bc(1,1).eq.EXT_DIR .and. uad(imin,j,k).ge.zero) then
                  stxhi(imin) = s(imin-1,j,k)
                  stxlo(imin) = s(imin-1,j,k)
               else if (bc(1,1).eq.EXT_DIR .and. uad(imin,j,k).lt.zero) then
                  stxlo(imin) = stxhi(imin)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
               else if (bc(1,1).eq.FOEXTRAP.or.bc(1,1).eq.HOEXTRAP) then
                  if (n.eq.XVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWXINFLOW
                        stxhi(imin) = MIN(stxhi(imin),zero)
#endif
                        stxlo(imin) = stxhi(imin)
                     else
                        if (uad(imin,j,k).ge.zero) then
#ifndef ALLOWXINFLOW
                           stxhi(imin) = MIN(stxhi(imin),zero)
#endif
                           stxlo(imin) = stxhi(imin)
                        endif
                     endif
                  else
                     stxlo(imin) = stxhi(imin)
#ifdef NOVERTICALINFLOW
c     Hack for no vertical inflow velocity
                     if (n.eq.ZVEL)then
                        stxlo(imin) = zero
                     endif
#endif
                  endif
               else if (bc(1,1).eq.REFLECT_EVEN) then
                  stxlo(imin) = stxhi(imin)
c               else if (bc(1,1).eq.FOEXTRAP.or.bc(1,1).eq.HOEXTRAP
c     &                 .or.bc(1,1).eq.REFLECT_EVEN) then
c                  stxlo(imin) = stxhi(imin)
               else if (bc(1,1).eq.REFLECT_ODD) then
                  stxhi(imin) = zero
                  stxlo(imin) = zero
               end if
               if (bc(1,2).eq.EXT_DIR .and. velpred.eq.1) then
                  stxlo(imax+1) = s(imax+1,j,k)
                  stxhi(imax+1) = s(imax+1,j,k)
               else if (bc(1,2).eq.EXT_DIR .and. uad(imax+1,j,k).le.zero) then
                  stxlo(imax+1) = s(imax+1,j,k)
                  stxhi(imax+1) = s(imax+1,j,k)
               else if (bc(1,2).eq.EXT_DIR .and. uad(imax+1,j,k).gt.zero) then
                  stxhi(imax+1) = stxlo(imax+1)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
               else if (bc(1,2).eq.FOEXTRAP.or.bc(1,2).eq.HOEXTRAP) then
                  if (n.eq.XVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWXINFLOW
                        stxlo(imax+1) = MAX(stxlo(imax+1),zero)
#endif
                        stxhi(imax+1) = stxlo(imax+1)
                     else
                        if (uad(imax+1,j,k).le.zero) then
#ifndef ALLOWXINFLOW
                           stxlo(imax+1) = MAX(stxlo(imax+1),zero)
#endif
                           stxhi(imax+1) = stxlo(imax+1)
                        endif
                     endif
                  else
                     stxhi(imax+1) = stxlo(imax+1)
#ifdef NOVERTICALINFLOW
c     Hack for no vertical inflow velocity
                     if (n.eq.ZVEL)then
                        stxhi(imax+1) = zero
                     endif
#endif
                  endif
               else if (bc(1,2).eq.REFLECT_EVEN) then
                  stxhi(imax+1) = stxlo(imax+1)
c               else if (bc(1,2).eq.FOEXTRAP.or.bc(1,2).eq.HOEXTRAP
c     &                 .or.bc(1,2).eq.REFLECT_EVEN) then
c                  stxhi(imax+1) = stxlo(imax+1)
               else if (bc(1,2).eq.REFLECT_ODD) then
                  stxlo(imax+1) = zero
                  stxhi(imax+1) = zero
               end if

               if ( velpred .eq. 1 ) then
                  do i = imin, imax+1
                     ltx = stxlo(i) .le. zero  .and.  stxhi(i) .ge. zero
                     ltx = ltx .or. (abs(stxlo(i)+stxhi(i)) .lt. eps)
                     stx = cvmgp(stxlo(i),stxhi(i),stxlo(i)+stxhi(i))
                     xstate(i,j,k) = cvmgt(zero,stx,ltx)
                  end do
               else
                  do i = imin, imax+1
                     xstate(i,j,k) = cvmgp(stxlo(i),stxhi(i),uedge(i,j,k))
                     xstate(i,j,k) = cvmgt(half*(stxlo(i)+stxhi(i)),xstate(i,j,k)
     &                    ,abs(uedge(i,j,k)).lt.eps)
                  end do
               end if
               place_to_break = 1
            end do
         end do
      end if
c
c     compute the yedge states
c
      if ((velpred.ne.1) .or. (n.eq.YVEL)) then
         do k = kmin,kmax
            do i = imin,imax
               do j = jmin,jmax+1

                  stylo(j) = ylo(i,j,k)
     &                 - dt4x*(uad(i+1,j-1,k  )+uad(i,j-1,k))*
     &                 (xzlo(i+1,j-1,k  )-xzlo(i,j-1,k))
     &                 - dt4z*(wad(i  ,j-1,k+1)+wad(i,j-1,k))*
     &                 (zxlo(i  ,j-1,k+1)-zxlo(i,j-1,k))
                  styhi(j) = yhi(i,j,k)
     &                 - dt4x*(uad(i+1,j  ,k  )+uad(i,j  ,k))*
     &                 (xzlo(i+1,j  ,k  )-xzlo(i,j  ,k))
     &                 - dt4z*(wad(i  ,j  ,k+1)+wad(i,j  ,k))*
     &                 (zxlo(i  ,j  ,k+1)-zxlo(i,j  ,k))

                  if (use_minion.eq.0) then
                     stylo(j) = stylo(j) + dth*tforces(i,j-1,k)
                     styhi(j) = styhi(j) + dth*tforces(i,j,  k)
                  end if

               end do

               if (bc(2,1).eq.EXT_DIR .and. velpred.eq.1) then
                  styhi(jmin) = s(i,jmin-1,k)
                  stylo(jmin) = s(i,jmin-1,k)
               else if (bc(2,1).eq.EXT_DIR .and. vad(i,jmin,k).ge.zero) then
                  styhi(jmin) = s(i,jmin-1,k)
                  stylo(jmin) = s(i,jmin-1,k)
               else if (bc(2,1).eq.EXT_DIR .and. vad(i,jmin,k).lt.zero) then
                  stylo(jmin) = styhi(jmin)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
               else if (bc(2,1).eq.FOEXTRAP.or.bc(2,1).eq.HOEXTRAP) then
                  if (n.eq.YVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWYINFLOW
                        styhi(jmin) = MIN(styhi(jmin),zero)
#endif
                        stylo(jmin) = styhi(jmin)
                     else
                        if (vad(i,jmin,k).ge.zero) then
#ifndef ALLOWYINFLOW
                           styhi(jmin) = MIN(styhi(jmin),zero)
#endif
                           stylo(jmin) = styhi(jmin)
                        endif
                     endif
                  else
                     stylo(jmin) = styhi(jmin)
#ifdef NOVERTICALINFLOW
c     Hack for no vertical inflow velocity
                     if (n.eq.ZVEL)then
                        stylo(jmin) = zero
                     endif
#endif
                  endif
               else if (bc(2,1).eq.REFLECT_EVEN) then
                  stylo(jmin) = styhi(jmin)
c               else if (bc(2,1).eq.FOEXTRAP.or.bc(2,1).eq.HOEXTRAP
c     &                 .or.bc(2,1).eq.REFLECT_EVEN) then
c                  stylo(jmin) = styhi(jmin)
               else if (bc(2,1).eq.REFLECT_ODD) then
                  styhi(jmin) = zero
                  stylo(jmin) = zero
               end if
               
               if (bc(2,2).eq.EXT_DIR .and. velpred.eq.1) then
                  stylo(jmax+1) = s(i,jmax+1,k)
                  styhi(jmax+1) = s(i,jmax+1,k)
               else if (bc(2,2).eq.EXT_DIR .and. vad(i,jmax+1,k).le.zero) then
                  stylo(jmax+1) = s(i,jmax+1,k)
                  styhi(jmax+1) = s(i,jmax+1,k)
               else if (bc(2,2).eq.EXT_DIR .and. vad(i,jmax+1,k).gt.zero) then
                  styhi(jmax+1) = stylo(jmax+1)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
               else if (bc(2,2).eq.FOEXTRAP.or.bc(2,2).eq.HOEXTRAP) then
                  if (n.eq.YVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWYINFLOW
                        stylo(jmax+1) = MAX(stylo(jmax+1),zero)
#endif
                        styhi(jmax+1) = stylo(jmax+1)
                     else
                        if (vad(i,jmax+1,k).le.zero) then
#ifndef ALLOWYINFLOW
                           stylo(jmax+1) = MAX(stylo(jmax+1),zero)
#endif
                           styhi(jmax+1) = stylo(jmax+1)
                        endif
                     endif
                  else
                     styhi(jmax+1) = stylo(jmax+1)
#ifdef NOVERTICALINFLOW
c     Hack for no vertical inflow velocity
                     if (n.eq.ZVEL)then
                        styhi(jmax+1) = zero
                     endif
#endif
                  endif
               else if (bc(2,2).eq.REFLECT_EVEN) then
                  styhi(jmax+1) = stylo(jmax+1)
c               else if (bc(2,2).eq.FOEXTRAP.or.bc(2,2).eq.HOEXTRAP
c     &                 .or.bc(2,2).eq.REFLECT_EVEN) then
c                  styhi(jmax+1) = stylo(jmax+1)
               else if (bc(2,2).eq.REFLECT_ODD) then
                  stylo(jmax+1) = zero
                  styhi(jmax+1) = zero
               end if

               if ( velpred .eq. 1 ) then
                  do j = jmin, jmax+1
                     lty = stylo(j) .le. zero  .and.  styhi(j) .ge. zero
                     lty = lty .or. (abs(stylo(j)+styhi(j)) .lt. eps)
                     sty = cvmgp(stylo(j),styhi(j),stylo(j)+styhi(j))
                     ystate(i,j,k) = cvmgt(zero,sty,lty)
                  end do
               else
                  do j=jmin,jmax+1
                     ystate(i,j,k) = cvmgp(stylo(j),styhi(j),vedge(i,j,k))
                     ystate(i,j,k) = cvmgt(half*(stylo(j)+styhi(j)),ystate(i,j,k),
     &                    abs(vedge(i,j,k)).lt.eps)
                  end do
               end if
               place_to_break = 1
            end do
         end do
      end if
c
c     compute the zedge states
c
      if ((velpred.ne.1) .or. (n.eq.ZVEL)) then
         do j = jmin,jmax
            do i = imin,imax
               do k = kmin,kmax+1

                  stzlo(k) = zlo(i,j,k)
     &                 - dt4x*(uad(i+1,j  ,k-1)+uad(i,j,k-1))
     &                 *(xylo(i+1,j  ,k-1)-xylo(i,j,k-1))
     &                 - dt4y*(vad(i  ,j+1,k-1)+vad(i,j,k-1))
     &                 *(yxlo(i  ,j+1,k-1)-yxlo(i,j,k-1))
                  
                  stzhi(k) = zhi(i,j,k)
     &                 - dt4x*(uad(i+1,j  ,k  )+uad(i,j,k  ))
     &                 *(xylo(i+1,j  ,k  )-xylo(i,j,k  ))
     &                 - dt4y*(vad(i  ,j+1,k  )+vad(i,j,k  ))
     &                 *(yxlo(i  ,j+1,k  )-yxlo(i,j,k  ))
                  
                  if (use_minion.eq.0) then
                     stzlo(k) = stzlo(k) + dth*tforces(i,j,k-1)
                     stzhi(k) = stzhi(k) + dth*tforces(i,j,k)
                  end if

               end do

               if (bc(3,1).eq.EXT_DIR .and. velpred.eq.1) then
                  stzlo(kmin) = s(i,j,kmin-1)
                  stzhi(kmin) = s(i,j,kmin-1)
               else if (bc(3,1).eq.EXT_DIR .and. wad(i,j,kmin).ge.zero) then
                  stzlo(kmin) = s(i,j,kmin-1)
                  stzhi(kmin) = s(i,j,kmin-1)
               else if (bc(3,1).eq.EXT_DIR .and. wad(i,j,kmin).lt.zero) then
                  stzlo(kmin) = stzhi(kmin)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
               else if (bc(3,1).eq.FOEXTRAP.or.bc(3,1).eq.HOEXTRAP) then
                  if (n.eq.ZVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWZINFLOW
                        stzhi(kmin) = MIN(stzhi(kmin),zero)
#endif
                        stzlo(kmin) = stzhi(kmin)
                     else
                        if (wad(i,j,kmin).ge.zero) then
#ifndef ALLOWZINFLOW
                           stzhi(kmin) = MIN(stzhi(kmin),zero)
#endif
                           stzlo(kmin) = stzhi(kmin)
                        endif
                     endif
                  else
                     stzlo(kmin) = stzhi(kmin)
                  endif
               else if (bc(3,1).eq.REFLECT_EVEN) then
                  stzlo(kmin) = stzhi(kmin)
c               else if (bc(3,1).eq.FOEXTRAP.or.bc(3,1).eq.HOEXTRAP
c     &                 .or.bc(3,1).eq.REFLECT_EVEN) then
c                  stzlo(kmin) = stzhi(kmin)
               else if (bc(3,1).eq.REFLECT_ODD) then
                  stzlo(kmin) = zero
                  stzhi(kmin) = zero
               end if
               if (bc(3,2).eq.EXT_DIR .and. velpred.eq.1) then
                  stzlo(kmax+1) = s(i,j,kmax+1)
                  stzhi(kmax+1) = s(i,j,kmax+1)
               else if (bc(3,2).eq.EXT_DIR .and. wad(i,j,kmax+1).le.zero) then
                  stzlo(kmax+1) = s(i,j,kmax+1)
                  stzhi(kmax+1) = s(i,j,kmax+1)
               else if (bc(3,2).eq.EXT_DIR .and. wad(i,j,kmax+1).gt.zero) then
                  stzhi(kmax+1) = stzlo(kmax+1)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
               else if (bc(3,2).eq.FOEXTRAP.or.bc(3,2).eq.HOEXTRAP) then
                  if (n.eq.ZVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWZINFLOW
                        stzlo(kmax+1) = MAX(stzlo(kmax+1),zero)
#endif
                        stzhi(kmax+1) = stzlo(kmax+1)
                     else
                        if (wad(i,j,kmax+1).le.zero) then
#ifndef ALLOWZINFLOW
                           stzlo(kmax+1) = MAX(stzlo(kmax+1),zero)
#endif
                           stzhi(kmax+1) = stzlo(kmax+1)
                        endif
                     endif
                  else
                     stzhi(kmax+1) = stzlo(kmax+1)
                  endif
               else if (bc(3,2).eq.REFLECT_EVEN) then
                  stzhi(kmax+1) = stzlo(kmax+1)
c               else if (bc(3,2).eq.FOEXTRAP.or.bc(3,2).eq.HOEXTRAP
c     &                 .or.bc(3,2).eq.REFLECT_EVEN) then
c                  stzhi(kmax+1) = stzlo(kmax+1)
               else if (bc(3,2).eq.REFLECT_ODD) then
                  stzlo(kmax+1) = zero
                  stzhi(kmax+1) = zero
               end if

               if ( velpred .eq. 1 ) then
                  do k = kmin,kmax+1
                     ltz = stzlo(k) .le. zero  .and.  stzhi(k) .ge. zero
                     ltz = ltz .or. (abs(stzlo(k)+stzhi(k)) .lt. eps)
                     stz = cvmgp(stzlo(k),stzhi(k),stzlo(k)+stzhi(k))
                     zstate(i,j,k) = cvmgt(zero,stz,ltz)
                  end do
               else
                  do k = kmin,kmax+1
                     zstate(i,j,k) = cvmgp(stzlo(k),stzhi(k),wedge(i,j,k))
                     zstate(i,j,k) = cvmgt(half*(stzlo(k)+stzhi(k)),zstate(i,j,k),
     &                    abs(wedge(i,j,k)).lt.eps)
                  end do
               end if
               place_to_break = 1
            end do
         end do
      end if

      else

c    
c     ORIGINAL NON-CORNER COUPLING CODE
c

c
c     compute the xedge states
c
      if ((velpred.ne.1) .or. (n.eq.XVEL)) then
         do k = kmin,kmax
            do j = jmin,jmax

               do i = imin-1,imax+1
                  if (vad(i,j,k)*vad(i,j+1,k).lt.0.d0) then
                      vbar = 0.5d0*(vad(i,j,k)+vad(i,j+1,k))
                      if (vbar.lt.0.d0) then
                          inc = 1
                      else
                          inc = 0
                      endif
                      tr1 = vbar*(s(i,j+inc,k)-s(i,j+inc-1,k))/hy
                  else
                      tr1 = half*
     &                 (vad(i,j+1,k) + vad(i,j,k)) *
     &                 (yedge(i,j+1,k) - yedge(i,j,k)) / hy
                  endif

                  if (wad(i,j,k)*wad(i,j,k+1).lt.0.d0) then
                      wbar = 0.5d0*(wad(i,j,k)+wad(i,j,k+1))
                      if (wbar.lt.0.d0) then
                          inc = 1
                      else
                          inc = 0
                      endif
                      tr2 = wbar*(s(i,j,k+inc)-s(i,j,k+inc-1))/hz
                  else
                      tr2 = half*
     &                 (wad(i,j,k+1) + wad(i,j,k)) *
     &                 (zedge(i,j,k+1) - zedge(i,j,k)) / hz
                  endif

                  stxlo(i+1)= s(i,j,k) + (half-dthx*u(i,j,k))*sx(i,j,k)
     &                 - dth*tr1 - dth*tr2
     &                 + dth*tforces(i,j,k)
                  stxhi(i  )= s(i,j,k) - (half+dthx*u(i,j,k))*sx(i,j,k)
     &                 - dth*tr1 - dth*tr2
     &                 + dth*tforces(i,j,k)
               end do

               if (bc(1,1).eq.EXT_DIR .and. velpred.eq.1) then
                  stxhi(imin) = s(imin-1,j,k)
                  stxlo(imin) = s(imin-1,j,k)
               else if (bc(1,1).eq.EXT_DIR .and. uad(imin,j,k).ge.zero) then
                  stxhi(imin) = s(imin-1,j,k)
                  stxlo(imin) = s(imin-1,j,k)
               else if (bc(1,1).eq.EXT_DIR .and. uad(imin,j,k).lt.zero) then
                  stxlo(imin) = stxhi(imin)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
               else if (bc(1,1).eq.FOEXTRAP.or.bc(1,1).eq.HOEXTRAP) then
                  if (n.eq.XVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWXINFLOW
                        stxhi(imin) = MIN(stxhi(imin),zero)
#endif
                        stxlo(imin) = stxhi(imin)
                     else
                        if (uad(imin,j,k).ge.zero) then
#ifndef ALLOWXINFLOW
                           stxhi(imin) = MIN(stxhi(imin),zero)
#endif
                           stxlo(imin) = stxhi(imin)
                        endif
                     endif
                  else
                     stxlo(imin) = stxhi(imin)
#ifdef NOVERTICALINFLOW
c     Hack for no vertical inflow velocity
                     if (n.eq.ZVEL)then
                        stxlo(imin) = zero
                     endif
#endif
                  endif
               else if (bc(1,1).eq.REFLECT_EVEN) then
                  stxlo(imin) = stxhi(imin)
c               else if (bc(1,1).eq.FOEXTRAP.or.bc(1,1).eq.HOEXTRAP
c     &                 .or.bc(1,1).eq.REFLECT_EVEN) then
c                  stxlo(imin) = stxhi(imin)
               else if (bc(1,1).eq.REFLECT_ODD) then
                  stxhi(imin) = zero
                  stxlo(imin) = zero
               end if
               if (bc(1,2).eq.EXT_DIR .and. velpred.eq.1) then
                  stxlo(imax+1) = s(imax+1,j,k)
                  stxhi(imax+1) = s(imax+1,j,k)
               else if (bc(1,2).eq.EXT_DIR .and. uad(imax+1,j,k).le.zero) then
                  stxlo(imax+1) = s(imax+1,j,k)
                  stxhi(imax+1) = s(imax+1,j,k)
               else if (bc(1,2).eq.EXT_DIR .and. uad(imax+1,j,k).gt.zero) then
                  stxhi(imax+1) = stxlo(imax+1)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
               else if (bc(1,2).eq.FOEXTRAP.or.bc(1,2).eq.HOEXTRAP) then
                  if (n.eq.XVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWXINFLOW
                        stxlo(imax+1) = MAX(stxlo(imax+1),zero)
#endif
                        stxhi(imax+1) = stxlo(imax+1)
                     else
                        if (uad(imax+1,j,k).le.zero) then
#ifndef ALLOWXINFLOW
                           stxlo(imax+1) = MAX(stxlo(imax+1),zero)
#endif
                           stxhi(imax+1) = stxlo(imax+1)
                        endif
                     endif
                  else
                     stxhi(imax+1) = stxlo(imax+1)
#ifdef NOVERTICALINFLOW
c     Hack for no vertical inflow velocity
                     if (n.eq.ZVEL)then
                        stxhi(imax+1) = zero
                     endif
#endif
                  endif
               else if (bc(1,2).eq.REFLECT_EVEN) then
                  stxhi(imax+1) = stxlo(imax+1)
c               else if (bc(1,2).eq.FOEXTRAP.or.bc(1,2).eq.HOEXTRAP
c     &                 .or.bc(1,2).eq.REFLECT_EVEN) then
c                  stxhi(imax+1) = stxlo(imax+1)
               else if (bc(1,2).eq.REFLECT_ODD) then
                  stxlo(imax+1) = zero
                  stxhi(imax+1) = zero
               end if

               if ( velpred .eq. 1 ) then
                  do i = imin, imax+1
                     ltx = stxlo(i) .le. zero  .and.  stxhi(i) .ge. zero
                     ltx = ltx .or. (abs(stxlo(i)+stxhi(i)) .lt. eps)
                     stx = cvmgp(stxlo(i),stxhi(i),stxlo(i)+stxhi(i))
                     xstate(i,j,k) = cvmgt(zero,stx,ltx)
                  end do
               else
                  do i = imin, imax+1
                     xstate(i,j,k) = cvmgp(stxlo(i),stxhi(i),uedge(i,j,k))
                     xstate(i,j,k) = cvmgt(half*(stxlo(i)+stxhi(i)),xstate(i,j,k)
     &                    ,abs(uedge(i,j,k)).lt.eps)
                  end do
               end if
               place_to_break = 1
            end do
         end do
      end if
c
c     compute the yedge states
c
      if ((velpred.ne.1) .or. (n.eq.YVEL)) then
         do k = kmin,kmax
            do i = imin,imax
               
               do j = jmin-1,jmax+1
                  if (uad(i,j,k)*uad(i+1,j,k).lt.0.d0) then
                      ubar = 0.5d0*(uad(i,j,k)+uad(i+1,j,k))
                      if (ubar.lt.0.d0) then
                          inc = 1
                      else
                          inc = 0
                      endif
                      tr1 = ubar*(s(i+inc,j,k)-s(i+inc-1,j,k))/hx
                  else
                      tr1 = half*
     &                 (uad(i+1,j,k) + uad(i,j,k)) *
     &                 (xedge(i+1,j,k) - xedge(i,j,k)) / hx
                  endif

                  if (wad(i,j,k)*wad(i,j,k+1).lt.0.d0) then
                      wbar = 0.5d0*(wad(i,j,k)+wad(i,j,k+1))
                      if (wbar.lt.0.d0) then
                          inc = 1
                      else
                          inc = 0
                      endif
                      tr2 = wbar*(s(i,j,k+inc)-s(i,j,k+inc-1))/hz
                  else
                      tr2 = half*
     &                 (wad(i,j,k+1) + wad(i,j,k)) *
     &                 (zedge(i,j,k+1) - zedge(i,j,k)) / hz
                  endif

                  stylo(j+1)= s(i,j,k) + (half-dthy*v(i,j,k))*sy(i,j,k)
     &                 - dth*tr1 - dth*tr2
     &                 + dth*tforces(i,j,k)
                  styhi(j)  = s(i,j,k) - (half+dthy*v(i,j,k))*sy(i,j,k)
     &                 - dth*tr1 - dth*tr2
     &                 + dth*tforces(i,j,k)
               end do

               if (bc(2,1).eq.EXT_DIR .and. velpred.eq.1) then
                  styhi(jmin) = s(i,jmin-1,k)
                  stylo(jmin) = s(i,jmin-1,k)
               else if (bc(2,1).eq.EXT_DIR .and. vad(i,jmin,k).ge.zero) then
                  styhi(jmin) = s(i,jmin-1,k)
                  stylo(jmin) = s(i,jmin-1,k)
               else if (bc(2,1).eq.EXT_DIR .and. vad(i,jmin,k).lt.zero) then
                  stylo(jmin) = styhi(jmin)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
               else if (bc(2,1).eq.FOEXTRAP.or.bc(2,1).eq.HOEXTRAP) then
                  if (n.eq.YVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWYINFLOW
                        styhi(jmin) = MIN(styhi(jmin),zero)
#endif
                        stylo(jmin) = styhi(jmin)
                     else
                        if (vad(i,jmin,k).ge.zero) then
#ifndef ALLOWYINFLOW
                           styhi(jmin) = MIN(styhi(jmin),zero)
#endif
                           stylo(jmin) = styhi(jmin)
                        endif
                     endif
                  else
                     stylo(jmin) = styhi(jmin)
#ifdef NOVERTICALINFLOW
c     Hack for no vertical inflow velocity
                     if (n.eq.ZVEL)then
                        stylo(jmin) = zero
                     endif
#endif
                  endif
               else if (bc(2,1).eq.REFLECT_EVEN) then
                  stylo(jmin) = styhi(jmin)
c               else if (bc(2,1).eq.FOEXTRAP.or.bc(2,1).eq.HOEXTRAP
c     &                 .or.bc(2,1).eq.REFLECT_EVEN) then
c                  stylo(jmin) = styhi(jmin)
               else if (bc(2,1).eq.REFLECT_ODD) then
                  styhi(jmin) = zero
                  stylo(jmin) = zero
               end if
               
               if (bc(2,2).eq.EXT_DIR .and. velpred.eq.1) then
                  stylo(jmax+1) = s(i,jmax+1,k)
                  styhi(jmax+1) = s(i,jmax+1,k)
               else if (bc(2,2).eq.EXT_DIR .and. vad(i,jmax+1,k).le.zero) then
                  stylo(jmax+1) = s(i,jmax+1,k)
                  styhi(jmax+1) = s(i,jmax+1,k)
               else if (bc(2,2).eq.EXT_DIR .and. vad(i,jmax+1,k).gt.zero) then
                  styhi(jmax+1) = stylo(jmax+1)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
               else if (bc(2,2).eq.FOEXTRAP.or.bc(2,2).eq.HOEXTRAP) then
                  if (n.eq.YVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWYINFLOW
                        stylo(jmax+1) = MAX(stylo(jmax+1),zero)
#endif
                        styhi(jmax+1) = stylo(jmax+1)
                     else
                        if (vad(i,jmax+1,k).le.zero) then
#ifndef ALLOWYINFLOW
                           stylo(jmax+1) = MAX(stylo(jmax+1),zero)
#endif
                           styhi(jmax+1) = stylo(jmax+1)
                        endif
                     endif
                  else
                     styhi(jmax+1) = stylo(jmax+1)
#ifdef NOVERTICALINFLOW
c     Hack for no vertical inflow velocity
                     if (n.eq.ZVEL)then
                        styhi(jmax+1) = zero
                     endif
#endif
                  endif
               else if (bc(2,2).eq.REFLECT_EVEN) then
                  styhi(jmax+1) = stylo(jmax+1)
c               else if (bc(2,2).eq.FOEXTRAP.or.bc(2,2).eq.HOEXTRAP
c     &                 .or.bc(2,2).eq.REFLECT_EVEN) then
c                  styhi(jmax+1) = stylo(jmax+1)
               else if (bc(2,2).eq.REFLECT_ODD) then
                  stylo(jmax+1) = zero
                  styhi(jmax+1) = zero
               end if

               if ( velpred .eq. 1 ) then
                  do j = jmin, jmax+1
                     lty = stylo(j) .le. zero  .and.  styhi(j) .ge. zero
                     lty = lty .or. (abs(stylo(j)+styhi(j)) .lt. eps)
                     sty = cvmgp(stylo(j),styhi(j),stylo(j)+styhi(j))
                     ystate(i,j,k) = cvmgt(zero,sty,lty)
                  end do
               else
                  do j=jmin,jmax+1
                     ystate(i,j,k) = cvmgp(stylo(j),styhi(j),vedge(i,j,k))
                     ystate(i,j,k) = cvmgt(half*(stylo(j)+styhi(j)),ystate(i,j,k),
     &                    abs(vedge(i,j,k)).lt.eps)
                  end do
               end if
               place_to_break = 1
            end do
         end do
      end if
c
c     compute the zedge states
c
      if ((velpred.ne.1) .or. (n.eq.ZVEL)) then
         do j = jmin,jmax
            do i = imin,imax
               
               do k = kmin-1,kmax+1
                  if (uad(i,j,k)*uad(i+1,j,k).lt.0.d0) then
                      ubar = 0.5d0*(uad(i,j,k)+uad(i+1,j,k))
                      if (ubar.lt.0.d0) then
                          inc = 1
                      else
                          inc = 0
                      endif
                      tr1 = ubar*(s(i+inc,j,k)-s(i+inc-1,j,k))/hx
                  else
                      tr1 = half*
     &                 (uad(i+1,j,k) + uad(i,j,k)) *
     &                 (xedge(i+1,j,k) - xedge(i,j,k)) / hx
                  endif

                  if (vad(i,j,k)*vad(i,j+1,k).lt.0.d0) then
                      vbar = 0.5d0*(vad(i,j,k)+vad(i,j+1,k))
                      if (vbar.lt.0.d0) then
                          inc = 1
                      else
                          inc = 0
                      endif
                      tr2 = vbar*(s(i,j+inc,k)-s(i,j+inc-1,k))/hy
                  else
                      tr2 = half*
     &                 (vad(i,j+1,k) + vad(i,j,k)) *
     &                 (yedge(i,j+1,k) - yedge(i,j,k)) / hy
                  endif

                  stzlo(k+1)= s(i,j,k) + (half-dthz*w(i,j,k))*sz(i,j,k)
     &                 - dth*tr1 - dth*tr2
     &                 + dth*tforces(i,j,k)
                  stzhi(k)  = s(i,j,k) - (half+dthz*w(i,j,k))*sz(i,j,k)
     &                 - dth*tr1 - dth*tr2
     &                 + dth*tforces(i,j,k)
               end do

               if (bc(3,1).eq.EXT_DIR .and. velpred.eq.1) then
                  stzlo(kmin) = s(i,j,kmin-1)
                  stzhi(kmin) = s(i,j,kmin-1)
               else if (bc(3,1).eq.EXT_DIR .and. wad(i,j,kmin).ge.zero) then
                  stzlo(kmin) = s(i,j,kmin-1)
                  stzhi(kmin) = s(i,j,kmin-1)
               else if (bc(3,1).eq.EXT_DIR .and. wad(i,j,kmin).lt.zero) then
                  stzlo(kmin) = stzhi(kmin)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
               else if (bc(3,1).eq.FOEXTRAP.or.bc(3,1).eq.HOEXTRAP) then
                  if (n.eq.ZVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWZINFLOW
                        stzhi(kmin) = MIN(stzhi(kmin),zero)
#endif
                        stzlo(kmin) = stzhi(kmin)
                     else
                        if (wad(i,j,kmin).ge.zero) then
#ifndef ALLOWZINFLOW
                           stzhi(kmin) = MIN(stzhi(kmin),zero)
#endif
                           stzlo(kmin) = stzhi(kmin)
                        endif
                     endif
                  else
                     stzlo(kmin) = stzhi(kmin)
                  endif
               else if (bc(3,1).eq.REFLECT_EVEN) then
                  stzlo(kmin) = stzhi(kmin)
c               else if (bc(3,1).eq.FOEXTRAP.or.bc(3,1).eq.HOEXTRAP
c     &                 .or.bc(3,1).eq.REFLECT_EVEN) then
c                  stzlo(kmin) = stzhi(kmin)
               else if (bc(3,1).eq.REFLECT_ODD) then
                  stzlo(kmin) = zero
                  stzhi(kmin) = zero
               end if
               if (bc(3,2).eq.EXT_DIR .and. velpred.eq.1) then
                  stzlo(kmax+1) = s(i,j,kmax+1)
                  stzhi(kmax+1) = s(i,j,kmax+1)
               else if (bc(3,2).eq.EXT_DIR .and. wad(i,j,kmax+1).le.zero) then
                  stzlo(kmax+1) = s(i,j,kmax+1)
                  stzhi(kmax+1) = s(i,j,kmax+1)
               else if (bc(3,2).eq.EXT_DIR .and. wad(i,j,kmax+1).gt.zero) then
                  stzhi(kmax+1) = stzlo(kmax+1)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
               else if (bc(3,2).eq.FOEXTRAP.or.bc(3,2).eq.HOEXTRAP) then
                  if (n.eq.ZVEL) then
                     if (velpred.eq.1) then
#ifndef ALLOWZINFLOW
                        stzlo(kmax+1) = MAX(stzlo(kmax+1),zero)
#endif
                        stzhi(kmax+1) = stzlo(kmax+1)
                     else
                        if (wad(i,j,kmax+1).le.zero) then
#ifndef ALLOWZINFLOW
                           stzlo(kmax+1) = MAX(stzlo(kmax+1),zero)
#endif
                           stzhi(kmax+1) = stzlo(kmax+1)
                        endif
                     endif
                  else
                     stzhi(kmax+1) = stzlo(kmax+1)
                  endif
               else if (bc(3,2).eq.REFLECT_EVEN) then
                  stzhi(kmax+1) = stzlo(kmax+1)
c               else if (bc(3,2).eq.FOEXTRAP.or.bc(3,2).eq.HOEXTRAP
c     &                 .or.bc(3,2).eq.REFLECT_EVEN) then
c                  stzhi(kmax+1) = stzlo(kmax+1)
               else if (bc(3,2).eq.REFLECT_ODD) then
                  stzlo(kmax+1) = zero
                  stzhi(kmax+1) = zero
               end if

               if ( velpred .eq. 1 ) then
                  do k = kmin,kmax+1
                     ltz = stzlo(k) .le. zero  .and.  stzhi(k) .ge. zero
                     ltz = ltz .or. (abs(stzlo(k)+stzhi(k)) .lt. eps)
                     stz = cvmgp(stzlo(k),stzhi(k),stzlo(k)+stzhi(k))
                     zstate(i,j,k) = cvmgt(zero,stz,ltz)
                  end do
               else
                  do k = kmin,kmax+1
                     zstate(i,j,k) = cvmgp(stzlo(k),stzhi(k),wedge(i,j,k))
                     zstate(i,j,k) = cvmgt(half*(stzlo(k)+stzhi(k)),zstate(i,j,k),
     &                    abs(wedge(i,j,k)).lt.eps)
                  end do
               end if
               place_to_break = 1
            end do
         end do
      end if

      end if

      end

      subroutine FORT_ESTATE_FPU(s, tforces, divu, DIMS(s),
     &     xlo, xhi, sx, slxscr, stxlo, stxhi,
     &     uedge, DIMS(uedge), xstate, DIMS(xstate),

     &     ylo, yhi, sy, slyscr, stylo, styhi,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate),

     &     zlo, zhi, sz, slzscr, stzlo, stzhi,
     &     wedge, DIMS(wedge), zstate, DIMS(zstate),

     &     xedge, yedge, zedge,
     &     xylo, xzlo, yxlo, yzlo, zxlo, zylo,
     &     xyhi, xzhi, yxhi, yzhi, zxhi, zyhi,

     &     DIMS(work),
     &     bc,lo,hi,dt,dx,n,use_minion,iconserv)
c
c     This subroutine computes edges states, right now it uses
c     a lot of memory, but there becomes a trade off between
c     simplicity-efficiency in the new way of computing states
c     and complexity in the old way.  By eliminating loops over
c     state components though, the new way uses much less memory.
c
c     This routine differs from the default ESTATE function above in that
c     it assumes that the edge velocities are valid in a grow cell outside
c     the box, and no *ad (unprojected) velocities are used.  This routine
c     will fail if the UMAC coming in hasn't been "fillpatched"
c
      implicit none
      integer i,j,k,n,inc
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer imin,jmin,kmin,imax,jmax,kmax
      REAL_T place_to_break
      REAL_T hx, hy, hz, dt, dth, dthx, dthy, dthz
      REAL_T dt3, dt3x, dt3y, dt3z, dt4, dt4x, dt4y, dt4z
      REAL_T dt6, dt6x, dt6y, dt6z
      REAL_T tr,tr1,tr2,ubar,vbar,wbar,stx,sty,stz,fu,fv,fw,dx(SDIM)
      REAL_T eps,eps_for_bc
      logical ltx,lty,ltz
      parameter (eps        = 1.d-6 )
      parameter (eps_for_bc = 1.d-10)

      integer DIMDEC(s)
      integer DIMDEC(work)
      integer DIMDEC(uedge)
      integer DIMDEC(xstate)
      integer DIMDEC(vedge)
      integer DIMDEC(ystate)
      integer DIMDEC(wedge)
      integer DIMDEC(zstate)

      REAL_T s(DIMV(s))
      REAL_T stxlo(DIM1(s)),stxhi(DIM1(s)),slxscr(DIM1(s),4)
      REAL_T stylo(DIM2(s)),styhi(DIM2(s)),slyscr(DIM2(s),4)
      REAL_T stzlo(DIM3(s)),stzhi(DIM3(s)),slzscr(DIM3(s),4)

      REAL_T uedge(DIMV(uedge)), xstate(DIMV(xstate))
      REAL_T vedge(DIMV(vedge)), ystate(DIMV(ystate))
      REAL_T wedge(DIMV(wedge)), zstate(DIMV(zstate))

      REAL_T xlo(DIMV(work)), xhi(DIMV(work)), xedge(DIMV(work))
      REAL_T ylo(DIMV(work)), yhi(DIMV(work)), yedge(DIMV(work))
      REAL_T zlo(DIMV(work)), zhi(DIMV(work)), zedge(DIMV(work))
      REAL_T  sx(DIMV(work))
      REAL_T  sy(DIMV(work))
      REAL_T  sz(DIMV(work))
      REAL_T tforces(DIMV(work))
      REAL_T    divu(DIMV(work))

      logical corner_couple

c     used in 3d corner coupling
      REAL_T xylo(DIMV(work)), xyhi(DIMV(work))
      REAL_T xzlo(DIMV(work)), xzhi(DIMV(work))
      REAL_T yxlo(DIMV(work)), yxhi(DIMV(work))
      REAL_T yzlo(DIMV(work)), yzhi(DIMV(work))
      REAL_T zxlo(DIMV(work)), zxhi(DIMV(work))
      REAL_T zylo(DIMV(work)), zyhi(DIMV(work))

      integer use_minion, iconserv
      REAL_T spx,smx,spy,smy,spz,smz,st,denom

      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)
      dthz = half*dt / dx(3)
      dt3  = dt / 3.0d0
      dt3x = dt3 / dx(1)
      dt3y = dt3 / dx(2)
      dt3z = dt3 / dx(3)
      dt4  = dt / 4.0d0
      dt4x = dt4 / dx(1)
      dt4y = dt4 / dx(2)
      dt4z = dt4 / dx(3)
      dt6  = dt / 6.0d0
      dt6x = dt6 / dx(1)
      dt6y = dt6 / dx(2)
      dt6z = dt6 / dx(3)
      hx   = dx(1)
      hy   = dx(2)
      hz   = dx(3)
      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

      corner_couple = .false.

c
c     compute the slopes
c
      call FORT_SLOPES(ALL,
     &     s,DIMS(s),
     &     sx,sy,sz,DIMS(work),
     &     lo,hi,slxscr,slyscr,slzscr,bc)
c
c     trace the state to the cell edges
c
      do k = kmin-1,kmax+1
         do j = jmin-1,jmax+1
            do i = imin,  imax+1
               xlo(i,j,k) = s(i-1,j,k) + (half  - dthx*uedge(i,j,k))*sx(i-1,j,k)
               xhi(i,j,k) = s(i,  j,k) + (-half - dthx*uedge(i,j,k))*sx(i,  j,k)
            end do
         end do
      end do

      if(use_minion.eq.1)then
         do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin,  imax+1
                  xlo(i,j,k) = xlo(i,j,k) + dth*tforces(i-1,j,k)
                  xhi(i,j,k) = xhi(i,j,k) + dth*tforces(i,  j,k)
               end do
            end do
         end do
         if (iconserv .eq. 1) then
           do k = kmin-1,kmax+1
             do j = jmin-1,jmax+1
               do i = imin,  imax+1
                  xlo(i,j,k) = xlo(i,j,k) - dth*s(i-1,j,k)*divu(i-1,j,k)
                  xhi(i,j,k) = xhi(i,j,k) - dth*s(i  ,j,k)*divu(i,  j,k)
               end do
             end do
           end do
         end if

      end if

      call trans_xbc(
     &     s,DIMS(s),
     &     xlo,xhi,DIMS(work),uedge,DIMS(uedge),
     &     lo,hi,n,bc,eps_for_bc,.false.,.false.)

      do k = kmin-1,kmax+1
         do j = jmin-1,jmax+1
            do i = imin,  imax+1
               fu  = cvmgt(zero,one,abs(uedge(i,j,k)).lt.eps)
               stx = cvmgp(xlo(i,j,k),xhi(i,j,k),uedge(i,j,k))
               xedge(i,j,k) = fu*stx + (one - fu)*half*(xhi(i,j,k)+xlo(i,j,k))
            end do
         end do
      end do

      do k = kmin-1,kmax+1
         do j = jmin,  jmax+1
            do i = imin-1,imax+1
               ylo(i,j,k) = s(i,j-1,k) + (half  - dthy*vedge(i,j,k))*sy(i,j-1,k)
               yhi(i,j,k) = s(i,j, k)  + (-half - dthy*vedge(i,j,k))*sy(i,j, k)
            end do
         end do
      end do

      if (use_minion.eq.1)then
         do k = kmin-1,kmax+1
            do j = jmin, jmax+1
               do i = imin-1,  imax+1
                  ylo(i,j,k) = ylo(i,j,k) + dth*tforces(i,j-1,k)
                  yhi(i,j,k) = yhi(i,j,k) + dth*tforces(i,j,  k)
               end do
            end do
         end do
         if (iconserv .eq. 1) then
           do k = kmin-1,kmax+1
             do j = jmin, jmax+1
               do i = imin-1,  imax+1
                  ylo(i,j,k) = ylo(i,j,k) - dth*s(i,j-1,k)*divu(i,j-1,k)
                  yhi(i,j,k) = yhi(i,j,k) - dth*s(i,j  ,k)*divu(i,j,  k)
               end do
             end do
           end do
         end if
      end if

      call trans_ybc(
     &     s,DIMS(s),
     &     ylo,yhi,DIMS(work),vedge,DIMS(vedge),
     &     lo,hi,n,bc,eps_for_bc,.false.,.false.)

      do k = kmin-1,kmax+1
         do j = jmin,  jmax+1
            do i = imin-1,imax+1
               fv  = cvmgt(zero,one,abs(vedge(i,j,k)).lt.eps)
               sty = cvmgp(ylo(i,j,k),yhi(i,j,k),vedge(i,j,k))
               yedge(i,j,k) = fv*sty + (one - fv)*half*(yhi(i,j,k)+ylo(i,j,k))
            end do
         end do
      end do

      do k = kmin,kmax+1
         do j = jmin-1,jmax+1
            do i = imin-1,imax+1
               zlo(i,j,k) = s(i,j,k-1) + (half  - dthz*wedge(i,j,k))*sz(i,j,k-1)
               zhi(i,j,k) = s(i,j,k  ) + (-half - dthz*wedge(i,j,k))*sz(i,j,k  )
            end do
         end do
      end do

      if (use_minion.eq.1)then
         do k = kmin,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,  imax+1
                  zlo(i,j,k) = zlo(i,j,k) + dth*tforces(i,j,k-1)
                  zhi(i,j,k) = zhi(i,j,k) + dth*tforces(i,j,k)
               end do
            end do
         end do
         if (iconserv .eq. 1) then
           do k = kmin,kmax+1
             do j = jmin-1,jmax+1
               do i = imin-1,  imax+1
                  zlo(i,j,k) = zlo(i,j,k) - dth*s(i,j,k-1)*divu(i,j,k-1)
                  zhi(i,j,k) = zhi(i,j,k) - dth*s(i,j,k  )*divu(i,j,k  )
               end do
             end do
           end do
         end if
      end if

      call trans_zbc(
     &     s,DIMS(s),
     &     zlo,zhi,DIMS(work),wedge,DIMS(wedge),
     &     lo,hi,n,bc,eps_for_bc,.false.,.false.)

      do k = kmin,kmax+1
         do j = jmin-1,jmax+1
            do i = imin-1,imax+1
               fw  = cvmgt(zero,one,abs(wedge(i,j,k)).lt.eps)
               stz = cvmgp(zlo(i,j,k),zhi(i,j,k),wedge(i,j,k))
               zedge(i,j,k) = fw*stz + (one-fw)*half*(zhi(i,j,k)+zlo(i,j,k))
            end do
         end do
      end do

      if (corner_couple) then

c
c     NEW CORNER COUPLING CODE
c

c
c     compute the corner-coupled terms:
c     xylo/hi, xzlo/hi, yxlo/hi, yzlo/hi, zxlo/hi, zylo/hi
c

c     loop over appropriate xy faces
      if (iconserv.eq.1) then
         do k=kmin-1,kmax+1
            do j=jmin,jmax
               do i=imin,imax+1
                  xylo(i,j,k) = xlo(i,j,k)
     &                 - dt3y*(yedge(i-1,j+1,k)*vedge(i-1,j+1,k)
     &                 - yedge(i-1,j,k)*vedge(i-1,j,k))
                  xyhi(i,j,k) = xhi(i,j,k)
     &                 - dt3y*(yedge(i  ,j+1,k)*vedge(i  ,j+1,k)
     &                 - yedge(i  ,j,k)*vedge(i  ,j,k))  
               end do
            end do
         end do
      else
         do k=kmin-1,kmax+1
            do j=jmin,jmax
               do i=imin,imax+1
                  xylo(i,j,k) = xlo(i,j,k)
     &                 - dt6y*(vedge(i-1,j+1,k)+vedge(i-1,j,k))
     &                 *(yedge(i-1,j+1,k)-yedge(i-1,j,k))
                  xyhi(i,j,k) = xhi(i,j,k)
     &                 - dt6y*(vedge(i  ,j+1,k)+vedge(i  ,j,k))
     &                 *(yedge(i  ,j+1,k)-yedge(i  ,j,k))
               end do
            end do
         end do
      end if

c     boundary conditions
      call trans_xbc(
     &     s,DIMS(s),
     &     xylo,xyhi,DIMS(work),uedge,DIMS(uedge),
     &     lo,hi,n,bc,eps_for_bc,.true.,.false.)

c     upwind
      do k=kmin-1,kmax+1
         do j=jmin,jmax
            do i=imin,imax+1
               fu  = cvmgt(zero,one,abs(uedge(i,j,k)).lt.eps)
               stx = cvmgp(xylo(i,j,k),xyhi(i,j,k),uedge(i,j,k))
               xylo(i,j,k) = fu*stx + (one - fu)*half*(xyhi(i,j,k)+xylo(i,j,k))
            end do
         end do
      end do

c     loop over appropriate xz faces
      if (iconserv.eq.1) then
         do k=kmin,kmax
            do j=jmin-1,jmax+1
               do i=imin,imax+1
                  xzlo(i,j,k) = xlo(i,j,k)
     &                 - dt3z*(zedge(i-1,j,k+1)*wedge(i-1,j,k+1)
     &                 - zedge(i-1,j,k)*wedge(i-1,j,k))
                  xzhi(i,j,k) = xhi(i,j,k)
     &                 - dt3z*(zedge(i  ,j,k+1)*wedge(i  ,j,k+1)
     &                 - zedge(i  ,j,k)*wedge(i  ,j,k))
               end do
            end do
         end do
      else
         do k=kmin,kmax
            do j=jmin-1,jmax+1
               do i=imin,imax+1
                  xzlo(i,j,k) = xlo(i,j,k)
     &                 - dt6z*(wedge(i-1,j,k+1)+wedge(i-1,j,k))
     &                 *(zedge(i-1,j,k+1)-zedge(i-1,j,k))
                  xzhi(i,j,k) = xhi(i,j,k)
     &                 - dt6z*(wedge(i  ,j,k+1)+wedge(i  ,j,k))
     &                 *(zedge(i  ,j,k+1)-zedge(i  ,j,k))
               end do
            end do
         end do
      end if

c     boundary conditions
      call trans_xbc(
     &     s,DIMS(s),
     &     xzlo,xzhi,DIMS(work),uedge,DIMS(uedge),
     &     lo,hi,n,bc,eps_for_bc,.false.,.true.)

c     upwind
      do k=kmin,kmax
         do j=jmin-1,jmax+1
            do i=imin,imax+1
               fu  = cvmgt(zero,one,abs(uedge(i,j,k)).lt.eps)
               stx = cvmgp(xzlo(i,j,k),xzhi(i,j,k),uedge(i,j,k))
               xzlo(i,j,k) = fu*stx + (one - fu)*half*(xzhi(i,j,k)+xzlo(i,j,k))
            end do
         end do
      end do

c     loop over appropriate yx faces
      if (iconserv.eq.1) then
         do k=kmin-1,kmax+1
            do j=jmin,jmax+1
               do i=imin,imax
                  yxlo(i,j,k) = ylo(i,j,k)
     &                 - dt3x*(xedge(i+1,j-1,k)*uedge(i+1,j-1,k)
     &                 - xedge(i,j-1,k)*uedge(i,j-1,k))
                  yxhi(i,j,k) = yhi(i,j,k)
     &                 - dt3x*(xedge(i+1,j  ,k)*uedge(i+1,j  ,k)
     &                 - xedge(i,j  ,k)*uedge(i,j  ,k))
               end do
            end do
         end do
      else
         do k=kmin-1,kmax+1
            do j=jmin,jmax+1
               do i=imin,imax
                  yxlo(i,j,k) = ylo(i,j,k)
     &                 - dt6x*(uedge(i+1,j-1,k)+uedge(i,j-1,k))
     &                 *(xedge(i+1,j-1,k)-xedge(i,j-1,k))
                  yxhi(i,j,k) = yhi(i,j,k)
     &                 - dt6x*(uedge(i+1,j  ,k)+uedge(i,j  ,k))
     &                 *(xedge(i+1,j  ,k)-xedge(i,j  ,k))
                  
               end do
            end do
         end do
      end if

c     boundary conditions
      call trans_ybc(
     &     s,DIMS(s),
     &     yxlo,yxhi,DIMS(work),vedge,DIMS(vedge),
     &     lo,hi,n,bc,eps_for_bc,.true.,.false.)

c     upwind
      do k=kmin-1,kmax+1
         do j=jmin,jmax+1
            do i=imin,imax
               fv  = cvmgt(zero,one,abs(vedge(i,j,k)).lt.eps)
               sty = cvmgp(yxlo(i,j,k),yxhi(i,j,k),vedge(i,j,k))
               yxlo(i,j,k) = fv*sty + (one - fv)*half*(yxhi(i,j,k)+yxlo(i,j,k))
            end do
         end do
      end do

c     loop over appropriate yz faces
      if (iconserv.eq.1) then
         do k=kmin,kmax
            do j=jmin,jmax+1
               do i=imin-1,imax+1
                  yzlo(i,j,k) = ylo(i,j,k)
     &                 - dt3z*(zedge(i,j-1,k+1)*wedge(i,j-1,k+1)
     &                 - zedge(i,j-1,k)*wedge(i,j-1,k))
                  yzhi(i,j,k) = yhi(i,j,k)
     &                 - dt3z*(zedge(i,j  ,k+1)*wedge(i,j  ,k+1)
     &                 - zedge(i,j  ,k)*wedge(i,j  ,k))
               end do
            end do
         end do
      else
         do k=kmin,kmax
            do j=jmin,jmax+1
               do i=imin-1,imax+1
                  yzlo(i,j,k) = ylo(i,j,k)
     &                 - dt6z*(wedge(i,j-1,k+1)+wedge(i,j-1,k))
     &                 *(zedge(i,j-1,k+1)-zedge(i,j-1,k))
                  yzhi(i,j,k) = yhi(i,j,k)
     &                 - dt6z*(wedge(i,j  ,k+1)+wedge(i,j  ,k))
     &                 *(zedge(i,j  ,k+1)-zedge(i,j  ,k))
               end do
            end do
         end do
      end if

c     boundary conditions
      call trans_ybc(
     &     s,DIMS(s),
     &     yzlo,yzhi,DIMS(work),vedge,DIMS(vedge),
     &     lo,hi,n,bc,eps_for_bc,.false.,.true.)

c     upwind
      do k=kmin,kmax
         do j=jmin,jmax+1
            do i=imin-1,imax+1
               fv  = cvmgt(zero,one,abs(vedge(i,j,k)).lt.eps)
               sty = cvmgp(yzlo(i,j,k),yzhi(i,j,k),vedge(i,j,k))
               yzlo(i,j,k) = fv*sty + (one - fv)*half*(yzhi(i,j,k)+yzlo(i,j,k))
            end do
         end do
      end do

c     loop over appropriate zx faces
      if (iconserv.eq.1) then
         do k=kmin,kmax+1
            do j=jmin-1,jmax+1
               do i=imin,imax
                  zxlo(i,j,k) = zlo(i,j,k)
     &                 - dt3x*(xedge(i+1,j,k-1)*uedge(i+1,j,k-1)
     &                 - xedge(i,j,k-1)*uedge(i,j,k-1))
                  zxhi(i,j,k) = zhi(i,j,k)
     &                 - dt3x*(xedge(i+1,j,k  )*uedge(i+1,j,k  )
     &                 - xedge(i,j,k  )*uedge(i,j,k  ))
               end do
            end do
         end do
      else
         do k=kmin,kmax+1
            do j=jmin-1,jmax+1
               do i=imin,imax
                  zxlo(i,j,k) = zlo(i,j,k)
     &                 - dt6x*(uedge(i+1,j,k-1)+uedge(i,j,k-1))
     &                 *(xedge(i+1,j,k-1)-xedge(i,j,k-1))
                  zxhi(i,j,k) = zhi(i,j,k)
     &                 - dt6x*(uedge(i+1,j,k  )+uedge(i,j,k  ))
     &                 *(xedge(i+1,j,k  )-xedge(i,j,k  ))
               end do
            end do
         end do
      end if

c     boundary conditions
      call trans_zbc(
     &     s,DIMS(s),
     &     zxlo,zxhi,DIMS(work),wedge,DIMS(wedge),
     &     lo,hi,n,bc,eps_for_bc,.true.,.false.)

c     upwind
      do k=kmin,kmax+1
         do j=jmin-1,jmax+1
            do i=imin,imax
               fw  = cvmgt(zero,one,abs(wedge(i,j,k)).lt.eps)
               stz = cvmgp(zxlo(i,j,k),zxhi(i,j,k),wedge(i,j,k))
               zxlo(i,j,k) = fw*stz + (one-fw)*half*(zxhi(i,j,k)+zxlo(i,j,k))
            end do
         end do
      end do

c     loop over appropriate zy faces
      if (iconserv.eq.1) then
         do k=kmin,kmax+1
            do j=jmin,jmax
               do i=imin-1,imax+1
                  zylo(i,j,k) = zlo(i,j,k)
     &                 - dt3y*(yedge(i,j+1,k-1)*vedge(i,j+1,k-1)
     &                 - yedge(i,j,k-1)*vedge(i,j,k-1))
                  zyhi(i,j,k) = zhi(i,j,k)
     &                 - dt3y*(yedge(i,j+1,k  )*vedge(i,j+1,k  )
     &                 - yedge(i,j,k  )*vedge(i,j,k  ))
               end do
            end do
         end do
      else
         do k=kmin,kmax+1
            do j=jmin,jmax
               do i=imin-1,imax+1
                  zylo(i,j,k) = zlo(i,j,k)
     &                 - dt6y*(vedge(i,j+1,k-1)+vedge(i,j,k-1))
     &                 *(yedge(i,j+1,k-1)-yedge(i,j,k-1))
                  zyhi(i,j,k) = zhi(i,j,k)
     &                 - dt6y*(vedge(i,j+1,k  )+vedge(i,j,k  ))
     &                 *(yedge(i,j+1,k  )-yedge(i,j,k  ))
               end do
            end do
         end do
      end if
         
c     boundary conditions
      call trans_zbc(
     &     s,DIMS(s),
     &     zylo,zyhi,DIMS(work),wedge,DIMS(wedge),
     &     lo,hi,n,bc,eps_for_bc,.false.,.true.)

c     upwind
      do k=kmin,kmax+1
         do j=jmin,jmax
            do i=imin-1,imax+1
               fw  = cvmgt(zero,one,abs(wedge(i,j,k)).lt.eps)
               stz = cvmgp(zylo(i,j,k),zyhi(i,j,k),wedge(i,j,k))
               zylo(i,j,k) = fw*stz + (one-fw)*half*(zyhi(i,j,k)+zylo(i,j,k))
            end do
         end do
      end do
      
c
c     compute the xedge states
c
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax+1
                  
               if (iconserv.eq.1) then

                  stxlo(i) = xlo(i,j,k)
     &                 - dthy*(yzlo(i-1,j+1,k  )*vedge(i-1,j+1,k  )
     &                 - yzlo(i-1,j,k)*vedge(i-1,j,k))
     &                 - dthz*(zylo(i-1,j  ,k+1)*wedge(i-1,j  ,k+1)
     &                 - zylo(i-1,j,k)*wedge(i-1,j,k))
     &                 + dthy*s(i-1,j,k)*(vedge(i-1,j+1,k)-vedge(i-1,j,k))
     &                 + dthz*s(i-1,j,k)*(wedge(i-1,j,k+1)-wedge(i-1,j,k))
                  stxhi(i) = xhi(i,j,k)
     &                 - dthy*(yzlo(i  ,j+1,k  )*vedge(i  ,j+1,  k)
     &                 - yzlo(i  ,j,k)*vedge(i  ,j,k))
     &                 - dthz*(zylo(i  ,j  ,k+1)*wedge(i  ,j  ,k+1)
     &                 - zylo(i  ,j,k)*wedge(i  ,j,k))
     &                 + dthy*s(i  ,j,k)*(vedge(i,j+1,k)-vedge(i,j,k))
     &                 + dthz*s(i  ,j,k)*(wedge(i,j,k+1)-wedge(i,j,k))
                
                  if (use_minion.eq.0) then
                     stxlo(i) = stxlo(i) - dth*s(i-1,j,k)*divu(i-1,j,k)
                     stxhi(i) = stxhi(i) - dth*s(i  ,j,k)*divu(i,  j,k)
                  end if

               else

                  stxlo(i) = xlo(i,j,k)
     &                 - dt4y*(vedge(i-1,j+1,k  )+vedge(i-1,j,k))*
     &                 (yzlo(i-1,j+1,k  )-yzlo(i-1,j,k))
     &                 - dt4z*(wedge(i-1,j  ,k+1)+wedge(i-1,j,k))*
     &                 (zylo(i-1,j  ,k+1)-zylo(i-1,j,k))
                  stxhi(i) = xhi(i,j,k)
     &                 - (dt4/hy)*(vedge(i  ,j+1,k  )+vedge(i  ,j,k))*
     &                 (yzlo(i  ,j+1,k  )-yzlo(i  ,j,k))
     &                 - (dt4/hz)*(wedge(i  ,j  ,k+1)+wedge(i  ,j,k))*
     &                 (zylo(i  ,j  ,k+1)-zylo(i  ,j,k))

               endif

               if (use_minion.eq.0) then
                  stxlo(i) = stxlo(i) + dth*tforces(i-1,j,k)
                  stxhi(i) = stxhi(i) + dth*tforces(i,  j,k)
               end if
               
            end do
            
            if (bc(1,1).eq.EXT_DIR .and. uedge(imin,j,k).ge.zero) then
               stxhi(imin) = s(imin-1,j,k)
               stxlo(imin) = s(imin-1,j,k)
            else if (bc(1,1).eq.EXT_DIR .and. uedge(imin,j,k).lt.zero) then
               stxlo(imin) = stxhi(imin)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
            else if (bc(1,1).eq.FOEXTRAP.or.bc(1,1).eq.HOEXTRAP) then
               if (n.eq.XVEL) then
                  if (uedge(imin,j,k).ge.zero) then
                     stxhi(imin) = MIN(stxhi(imin),zero)
                     stxlo(imin) = stxhi(imin)
                  endif
               else
                  stxlo(imin) = stxhi(imin)
               endif
            else if (bc(1,1).eq.REFLECT_EVEN) then
               stxlo(imin) = stxhi(imin)
c            else if (bc(1,1).eq.FOEXTRAP.or.bc(1,1).eq.HOEXTRAP
c     &              .or.bc(1,1).eq.REFLECT_EVEN) then
c               stxlo(imin) = stxhi(imin)
            else if (bc(1,1).eq.REFLECT_ODD) then
               stxhi(imin) = zero
               stxlo(imin) = zero
            end if
            if (bc(1,2).eq.EXT_DIR .and. uedge(imax+1,j,k).le.zero) then
               stxlo(imax+1) = s(imax+1,j,k)
               stxhi(imax+1) = s(imax+1,j,k)
            else if (bc(1,2).eq.EXT_DIR .and. uedge(imax+1,j,k).gt.zero) then
               stxhi(imax+1) = stxlo(imax+1)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
            else if (bc(1,2).eq.FOEXTRAP.or.bc(1,2).eq.HOEXTRAP) then
               if (n.eq.XVEL) then
                  if (uedge(imax+1,j,k).le.zero) then
                     stxlo(imax+1) = MAX(stxlo(imax+1),zero)
                     stxhi(imax+1) = stxlo(imax+1)
                  endif
               else
                  stxhi(imax+1) = stxlo(imax+1)
               endif
            else if (bc(1,2).eq.REFLECT_EVEN) then
               stxhi(imax+1) = stxlo(imax+1)
c            else if (bc(1,2).eq.FOEXTRAP.or.bc(1,2).eq.HOEXTRAP
c     &              .or.bc(1,2).eq.REFLECT_EVEN) then
c               stxhi(imax+1) = stxlo(imax+1)
            else if (bc(1,2).eq.REFLECT_ODD) then
               stxlo(imax+1) = zero
               stxhi(imax+1) = zero
            end if
            
            do i = imin, imax+1
               xstate(i,j,k) = cvmgp(stxlo(i),stxhi(i),uedge(i,j,k))
               xstate(i,j,k) = cvmgt(half*(stxlo(i)+stxhi(i)),xstate(i,j,k)
     &              ,abs(uedge(i,j,k)).lt.eps)
            end do
            place_to_break = 1
         end do
      end do
c
c     compute the yedge states
c
      do k = kmin,kmax
         do i = imin,imax
            do j = jmin,jmax+1

               if (iconserv.eq.1) then

                  stylo(j) = ylo(i,j,k)
     &                 - dthx*(xzlo(i+1,j-1,k  )*uedge(i+1,j-1,k  )
     &                 - xzlo(i,j-1,k)*uedge(i,j-1,k))
     &                 - dthz*(zxlo(i  ,j-1,k+1)*wedge(i  ,j-1,k+1)
     &                 - zxlo(i,j-1,k)*wedge(i,j-1,k))
     &                 + dthx*s(i,j-1,k)*(uedge(i+1,j-1,k)-uedge(i,j-1,k))
     &                 + dthz*s(i,j-1,k)*(wedge(i,j-1,k+1)-wedge(i,j-1,k))
                  styhi(j) = yhi(i,j,k)
     &                 - dthx*(xzlo(i+1,j  ,k  )*uedge(i+1,j  ,k  )
     &                 - xzlo(i,j  ,k)*uedge(i,j  ,k))
     &                 - dthz*(zxlo(i  ,j  ,k+1)*wedge(i  ,j  ,k+1)
     &                 - zxlo(i,j  ,k)*wedge(i,j  ,k))
     &                 + dthx*s(i,j  ,k)*(uedge(i+1,j,k)-uedge(i,j,k))
     &                 + dthz*s(i,j  ,k)*(wedge(i,j,k+1)-wedge(i,j,k))
                  
                  if (use_minion .eq. 0) then
                     stylo(j) = stylo(j) - dth*s(i,j-1,k)*divu(i,j-1,k)
                     styhi(j) = styhi(j) - dth*s(i,j  ,k)*divu(i,j,  k)
                  end if

               else
                  
                  stylo(j) = ylo(i,j,k)
     &                 - dt4x*(uedge(i+1,j-1,k  )+uedge(i,j-1,k))*
     &                 (xzlo(i+1,j-1,k  )-xzlo(i,j-1,k))
     &                 - dt4z*(wedge(i  ,j-1,k+1)+wedge(i,j-1,k))*
     &                 (zxlo(i  ,j-1,k+1)-zxlo(i,j-1,k))
                  styhi(j) = yhi(i,j,k)
     &                 - dt4x*(uedge(i+1,j  ,k  )+uedge(i,j  ,k))*
     &                 (xzlo(i+1,j  ,k  )-xzlo(i,j  ,k))
     &                 - dt4z*(wedge(i  ,j  ,k+1)+wedge(i,j  ,k))*
     &                 (zxlo(i  ,j  ,k+1)-zxlo(i,j  ,k))

               endif

               if (use_minion.eq.0) then
                  stylo(j) = stylo(j) + dth*tforces(i,j-1,k)
                  styhi(j) = styhi(j) + dth*tforces(i,j,  k)
               end if
               
            end do

            if (bc(2,1).eq.EXT_DIR .and. vedge(i,jmin,k).ge.zero) then
               styhi(jmin) = s(i,jmin-1,k)
               stylo(jmin) = s(i,jmin-1,k)
            else if (bc(2,1).eq.EXT_DIR .and. vedge(i,jmin,k).lt.zero) then
               stylo(jmin) = styhi(jmin)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
            else if (bc(2,1).eq.FOEXTRAP.or.bc(2,1).eq.HOEXTRAP) then
               if (n.eq.YVEL) then
                  if (vedge(i,jmin,k).ge.zero) then
                     styhi(jmin) = MIN(styhi(jmin),zero)
                     stylo(jmin) = styhi(jmin)
                  endif
               else
                  stylo(jmin) = styhi(jmin)
               endif
            else if (bc(2,1).eq.REFLECT_EVEN) then
               stylo(jmin) = styhi(jmin)
c            else if (bc(2,1).eq.FOEXTRAP.or.bc(2,1).eq.HOEXTRAP
c     &              .or.bc(2,1).eq.REFLECT_EVEN) then
c               stylo(jmin) = styhi(jmin)
            else if (bc(2,1).eq.REFLECT_ODD) then
               styhi(jmin) = zero
               stylo(jmin) = zero
            end if
            
            if (bc(2,2).eq.EXT_DIR .and. vedge(i,jmax+1,k).le.zero) then
               stylo(jmax+1) = s(i,jmax+1,k)
               styhi(jmax+1) = s(i,jmax+1,k)
            else if (bc(2,2).eq.EXT_DIR .and. vedge(i,jmax+1,k).le.zero) then
               styhi(jmax+1) = stylo(jmax+1)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
            else if (bc(2,2).eq.FOEXTRAP.or.bc(2,2).eq.HOEXTRAP) then
               if (n.eq.YVEL) then
                  if (vedge(i,jmax+1,k).le.zero) then
                     stylo(jmax+1) = MAX(stylo(jmax+1),zero)
                     styhi(jmax+1) = stylo(jmax+1)
                  endif
               else
                  styhi(jmax+1) = stylo(jmax+1)
               endif
            else if (bc(2,2).eq.REFLECT_EVEN) then
               styhi(jmax+1) = stylo(jmax+1)
c            else if (bc(2,2).eq.FOEXTRAP.or.bc(2,2).eq.HOEXTRAP
c     &              .or.bc(2,2).eq.REFLECT_EVEN) then
c               styhi(jmax+1) = stylo(jmax+1)
            else if (bc(2,2).eq.REFLECT_ODD) then
               stylo(jmax+1) = zero
               styhi(jmax+1) = zero
            end if
            
            do j=jmin,jmax+1
               ystate(i,j,k) = cvmgp(stylo(j),styhi(j),vedge(i,j,k))
               ystate(i,j,k) = cvmgt(half*(stylo(j)+styhi(j)),ystate(i,j,k),
     &              abs(vedge(i,j,k)).lt.eps)
            end do
            place_to_break = 1
         end do
      end do
c     
c     compute the zedge states
c
      do j = jmin,jmax
         do i = imin,imax
            do k = kmin,kmax+1
                  
               if (iconserv.eq.1) then
                 
                  stzlo(k) = zlo(i,j,k)
     &                 - dthx*(xylo(i+1,j  ,k-1)*uedge(i+1,j  ,k-1)
     &                 - xylo(i,j,k-1)*uedge(i,j,k-1))
     &                 - dthy*(yxlo(i  ,j+1,k-1)*vedge(i  ,j+1,k-1)
     &                 - yxlo(i,j,k-1)*vedge(i,j,k-1))
     &                 + dthx*s(i,j,k-1)*(uedge(i+1,j,k-1)-uedge(i,j,k-1))
     &                 + dthy*s(i,j,k-1)*(vedge(i,j+1,k-1)-vedge(i,j,k-1))
                  stzhi(k) = zhi(i,j,k)
     &                 - dthx*(xylo(i+1,j  ,k  )*uedge(i+1,j  ,k  )
     &                 - xylo(i,j,k  )*uedge(i,j,k  ))
     &                 - dthy*(yxlo(i  ,j+1,k  )*vedge(i  ,j+1,k  )
     &                 - yxlo(i,j,k  )*vedge(i,j,k  ))
     &                 + dthx*s(i,j,k  )*(uedge(i+1,j,k)-uedge(i,j,k))
     &                 + dthy*s(i,j,k  )*(vedge(i,j+1,k)-vedge(i,j,k))

                  if (use_minion.eq.0) then
                     stzlo(k) = stzlo(k) - dth*s(i,j,k-1)*divu(i,j,k-1)
                     stzhi(k) = stzhi(k) - dth*s(i,j,k  )*divu(i,j,k  )
                  end if

               else

                  stzlo(k) = zlo(i,j,k)
     &                 - dt4x*(uedge(i+1,j  ,k-1)+uedge(i,j,k-1))
     &                 *(xylo(i+1,j  ,k-1)-xylo(i,j,k-1))
     &                 - dt4y*(vedge(i  ,j+1,k-1)+vedge(i,j,k-1))
     &                 *(yxlo(i  ,j+1,k-1)-yxlo(i,j,k-1))
                  
                  stzhi(k) = zhi(i,j,k)
     &                 - dt4x*(uedge(i+1,j  ,k  )+uedge(i,j,k  ))
     &                 *(xylo(i+1,j  ,k  )-xylo(i,j,k  ))
     &                 - dt4y*(vedge(i  ,j+1,k  )+vedge(i,j,k  ))
     &                 *(yxlo(i  ,j+1,k  )-yxlo(i,j,k  ))

               endif

               if (use_minion.eq.0) then
                  stzlo(k) = stzlo(k) + dth*tforces(i,j,k-1)
                  stzhi(k) = stzhi(k) + dth*tforces(i,j,k)
               end if

            end do

            if (bc(3,1).eq.EXT_DIR .and. wedge(i,j,kmin).ge.zero) then
               stzlo(kmin) = s(i,j,kmin-1)
               stzhi(kmin) = s(i,j,kmin-1)
            else if (bc(3,1).eq.EXT_DIR .and. wedge(i,j,kmin).lt.zero) then
               stzlo(kmin) = stzhi(kmin)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
            else if (bc(3,1).eq.FOEXTRAP.or.bc(3,1).eq.HOEXTRAP) then
               if (n.eq.ZVEL) then
                  if (wedge(i,j,kmin).ge.zero) then
                     stzhi(kmin) = MIN(stzhi(kmin),zero)
                     stzlo(kmin) = stzhi(kmin)
                  endif
               else
                  stzlo(kmin) = stzhi(kmin)
               endif
            else if (bc(3,1).eq.REFLECT_EVEN) then
               stzlo(kmin) = stzhi(kmin)
c            else if (bc(3,1).eq.FOEXTRAP.or.bc(3,1).eq.HOEXTRAP
c     &              .or.bc(3,1).eq.REFLECT_EVEN) then
c               stzlo(kmin) = stzhi(kmin)
            else if (bc(3,1).eq.REFLECT_ODD) then
               stzlo(kmin) = zero
               stzhi(kmin) = zero
            end if
            if (bc(3,2).eq.EXT_DIR .and. wedge(i,j,kmax+1).le.zero) then
               stzlo(kmax+1) = s(i,j,kmax+1)
               stzhi(kmax+1) = s(i,j,kmax+1)
            else if (bc(3,2).eq.EXT_DIR .and. wedge(i,j,kmax+1).gt.zero) then
               stzhi(kmax+1) = stzlo(kmax+1)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
            else if (bc(3,2).eq.FOEXTRAP.or.bc(3,2).eq.HOEXTRAP) then
               if (n.eq.ZVEL) then
                  if (wedge(i,j,kmax+1).le.zero) then
                     stzlo(kmax+1) = MAX(stzlo(kmax+1),zero)
                     stzhi(kmax+1) = stzlo(kmax+1)
                  endif
               else
                  stzhi(kmax+1) = stzlo(kmax+1)
               endif
            else if (bc(3,2).eq.REFLECT_EVEN) then
               stzhi(kmax+1) = stzlo(kmax+1)
c            else if (bc(3,2).eq.FOEXTRAP.or.bc(3,2).eq.HOEXTRAP
c     &              .or.bc(3,2).eq.REFLECT_EVEN) then
c               stzhi(kmax+1) = stzlo(kmax+1)
            else if (bc(3,2).eq.REFLECT_ODD) then
               stzlo(kmax+1) = zero
               stzhi(kmax+1) = zero
            end if
               
            do k = kmin,kmax+1
               zstate(i,j,k) = cvmgp(stzlo(k),stzhi(k),wedge(i,j,k))
               zstate(i,j,k) = cvmgt(half*(stzlo(k)+stzhi(k)),zstate(i,j,k),
     &              abs(wedge(i,j,k)).lt.eps)
            end do
            place_to_break = 1
         end do
      end do
      
      else

c    
c     ORIGINAL NON-CORNER COUPLING CODE
c

c
c     compute the xedge states
c
      do k = kmin,kmax
            do j = jmin,jmax
               do i = imin-1,imax+1

                  spx = s(i,j,k) + half*sx(i,j,k)
                  smx = s(i,j,k) - half*sx(i,j,k)

                  if (iconserv.eq.1) then

                     tr =
     &                    (vedge(i,j+1,k)*yedge(i,j+1,k) - vedge(i,j,k)*yedge(i,j,k))/hy +   
     &                    (wedge(i,j,k+1)*zedge(i,j,k+1) - wedge(i,j,k)*zedge(i,j,k))/hz   

                     st = -dth*tr + dth*(tforces(i,j,k) - s(i,j,k)*divu(i,j,k))
     &                    + dth*s(i,j,k)*(vedge(i,j+1,k)-vedge(i,j,k))/hy
     &                    + dth*s(i,j,k)*(wedge(i,j,k+1)-wedge(i,j,k))/hz
c    &                    - dth*s(i,j,k)*(uedge(i+1,j,k)-uedge(i,j,k))/hx

                  else

                     if (vedge(i,j,k)*vedge(i,j+1,k).le.0.d0) then
                        vbar = 0.5d0*(vedge(i,j,k)+vedge(i,j+1,k))
                        if (vbar.lt.0.d0) then
                           inc = 1
                        else
                           inc = 0
                        endif
                        tr1 = vbar*(s(i,j+inc,k)-s(i,j+inc-1,k))/hy
                     else
                        tr1 = half*(vedge(i,j+1,k) + vedge(i,j,k)) *
     &                               (yedge(i,j+1,k) -   yedge(i,j,k)  ) / hy
                     endif
                     if (wedge(i,j,k)*wedge(i,j,k+1).lt.0.d0) then
                        wbar = 0.5d0*(wedge(i,j,k)+wedge(i,j,k+1))
                        if (wbar.lt.0.d0) then
                           inc = 1
                        else
                           inc = 0
                        endif
                        tr2 = wbar*(s(i,j,k+inc)-s(i,j,k+inc-1))/hz
                     else
                        tr2 = half*(wedge(i,j,k+1) + wedge(i,j,k)) *
     &                               (zedge(i,j,k+1) -   zedge(i,j,k)  ) / hz
                     endif

                     st = -dth*(tr1 + tr2) + dth*tforces(i,j,k)
                  endif

                  stxlo(i+1)= spx - dthx*uedge(i+1,j,k)*sx(i,j,k) + st
                  stxhi(i  )= smx - dthx*uedge(i  ,j,k)*sx(i,j,k) + st

               end do

               if (bc(1,1).eq.EXT_DIR .and. uedge(imin,j,k).ge.zero) then
                  stxhi(imin) = s(imin-1,j,k)
                  stxlo(imin) = s(imin-1,j,k)
               else if (bc(1,1).eq.EXT_DIR .and. uedge(imin,j,k).lt.zero) then
                  stxlo(imin) = stxhi(imin)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
               else if (bc(1,1).eq.FOEXTRAP.or.bc(1,1).eq.HOEXTRAP) then
                  if (n.eq.XVEL) then
                     if (uedge(imin,j,k).ge.zero) then
                        stxhi(imin) = MIN(stxhi(imin),zero)
                        stxlo(imin) = stxhi(imin)
                     endif
                  else
                     stxlo(imin) = stxhi(imin)
                  endif
               else if (bc(1,1).eq.REFLECT_EVEN) then
                  stxlo(imin) = stxhi(imin)
c               else if (bc(1,1).eq.FOEXTRAP.or.bc(1,1).eq.HOEXTRAP
c     &                 .or.bc(1,1).eq.REFLECT_EVEN) then
c                  stxlo(imin) = stxhi(imin)
               else if (bc(1,1).eq.REFLECT_ODD) then
                  stxhi(imin) = zero
                  stxlo(imin) = zero
               end if
               if (bc(1,2).eq.EXT_DIR .and. uedge(imax+1,j,k).le.zero) then
                  stxlo(imax+1) = s(imax+1,j,k)
                  stxhi(imax+1) = s(imax+1,j,k)
               else if (bc(1,2).eq.EXT_DIR .and. uedge(imax+1,j,k).gt.zero) then
                  stxhi(imax+1) = stxlo(imax+1)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
               else if (bc(1,2).eq.FOEXTRAP.or.bc(1,2).eq.HOEXTRAP) then
                  if (n.eq.XVEL) then
                     if (uedge(imax+1,j,k).le.zero) then
                        stxlo(imax+1) = MAX(stxlo(imax+1),zero)
                        stxhi(imax+1) = stxlo(imax+1)
                     endif
                  else
                     stxhi(imax+1) = stxlo(imax+1)
                  endif
               else if (bc(1,2).eq.REFLECT_EVEN) then
                  stxhi(imax+1) = stxlo(imax+1)
c               else if (bc(1,2).eq.FOEXTRAP.or.bc(1,2).eq.HOEXTRAP
c     &                 .or.bc(1,2).eq.REFLECT_EVEN) then
c                  stxhi(imax+1) = stxlo(imax+1)
               else if (bc(1,2).eq.REFLECT_ODD) then
                  stxlo(imax+1) = zero
                  stxhi(imax+1) = zero
               end if

               do i = imin, imax+1
                  xstate(i,j,k) = cvmgp(stxlo(i),stxhi(i),uedge(i,j,k))
                  xstate(i,j,k) = cvmgt(half*(stxlo(i)+stxhi(i)),xstate(i,j,k)
     &                 ,abs(uedge(i,j,k)).lt.eps)
               end do
               place_to_break = 1
            end do
      end do
c
c     compute the yedge states
c
      do k = kmin,kmax
            do i = imin,imax
               
               do j = jmin-1,jmax+1
                  spy = s(i,j,k) + half*sy(i,j,k)
                  smy = s(i,j,k) - half*sy(i,j,k)
                  if (iconserv.eq.1) then

                     tr =
     &                    (uedge(i+1,j,k)*xedge(i+1,j,k) - uedge(i,j,k)*xedge(i,j,k))/hx +   
     &                    (wedge(i,j,k+1)*zedge(i,j,k+1) - wedge(i,j,k)*zedge(i,j,k))/hz   

                     st = -dth*tr + dth*(tforces(i,j,k) - s(i,j,k)*divu(i,j,k))
     &                    + dth*s(i,j,k)*(uedge(i+1,j,k)-uedge(i,j,k))/hx
     &                    + dth*s(i,j,k)*(wedge(i,j,k+1)-wedge(i,j,k))/hz
c    &                    - dth*s(i,j,k)*(vedge(i,j+1,k)-vedge(i,j,k))/hy


                  else

                     if (uedge(i,j,k)*uedge(i+1,j,k).le.0.d0) then
                        ubar = 0.5d0*(uedge(i,j,k)+uedge(i+1,j,k))
                        if (ubar.lt.0.d0) then
                           inc = 1
                        else
                           inc = 0
                        endif
                        tr1 = ubar*(s(i+inc,j,k)-s(i+inc-1,j,k))/hx
                     else
                        tr1 = half*(uedge(i+1,j,k) + uedge(i,j,k)) *
     &                               (xedge(i+1,j,k) -   xedge(i,j,k)  ) / hx
                     endif
                     if (wedge(i,j,k)*wedge(i,j,k+1).lt.0.d0) then
                        wbar = 0.5d0*(wedge(i,j,k)+wedge(i,j,k+1))
                        if (wbar.lt.0.d0) then
                           inc = 1
                        else
                           inc = 0
                        endif
                        tr2 = wbar*(s(i,j,k+inc)-s(i,j,k+inc-1))/hz
                     else
                        tr2 = half*(wedge(i,j,k+1) + wedge(i,j,k)) *
     &                               (zedge(i,j,k+1) -   zedge(i,j,k)  ) / hz
                     endif

                     st = -dth*(tr1 + tr2) + dth*tforces(i,j,k)

                  endif

                  stylo(j+1)= spy - dthy*vedge(i,j+1,k)*sy(i,j,k) + st
                  styhi(j  )= smy - dthy*vedge(i,j  ,k)*sy(i,j,k) + st
               end do

               if (bc(2,1).eq.EXT_DIR .and. vedge(i,jmin,k).ge.zero) then
                  styhi(jmin) = s(i,jmin-1,k)
                  stylo(jmin) = s(i,jmin-1,k)
               else if (bc(2,1).eq.EXT_DIR .and. vedge(i,jmin,k).lt.zero) then
                  stylo(jmin) = styhi(jmin)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
               else if (bc(2,1).eq.FOEXTRAP.or.bc(2,1).eq.HOEXTRAP) then
                  if (n.eq.YVEL) then
                     if (vedge(i,jmin,k).ge.zero) then
                        styhi(jmin) = MIN(styhi(jmin),zero)
                        stylo(jmin) = styhi(jmin)
                     endif
                  else
                     stylo(jmin) = styhi(jmin)
                  endif
               else if (bc(2,1).eq.REFLECT_EVEN) then
                  stylo(jmin) = styhi(jmin)
c               else if (bc(2,1).eq.FOEXTRAP.or.bc(2,1).eq.HOEXTRAP
c     &                 .or.bc(2,1).eq.REFLECT_EVEN) then
c                  stylo(jmin) = styhi(jmin)
               else if (bc(2,1).eq.REFLECT_ODD) then
                  styhi(jmin) = zero
                  stylo(jmin) = zero
               end if
               
               if (bc(2,2).eq.EXT_DIR .and. vedge(i,jmax+1,k).le.zero) then
                  stylo(jmax+1) = s(i,jmax+1,k)
                  styhi(jmax+1) = s(i,jmax+1,k)
               else if (bc(2,2).eq.EXT_DIR .and. vedge(i,jmax+1,k).le.zero) then
                  styhi(jmax+1) = stylo(jmax+1)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
               else if (bc(2,2).eq.FOEXTRAP.or.bc(2,2).eq.HOEXTRAP) then
                  if (n.eq.YVEL) then
                     if (vedge(i,jmax+1,k).le.zero) then
                        stylo(jmax+1) = MAX(stylo(jmax+1),zero)
                        styhi(jmax+1) = stylo(jmax+1)
                     endif
                  else
                     styhi(jmax+1) = stylo(jmax+1)
                  endif
               else if (bc(2,2).eq.REFLECT_EVEN) then
                  styhi(jmax+1) = stylo(jmax+1)
c               else if (bc(2,2).eq.FOEXTRAP.or.bc(2,2).eq.HOEXTRAP
c     &                 .or.bc(2,2).eq.REFLECT_EVEN) then
c                  styhi(jmax+1) = stylo(jmax+1)
               else if (bc(2,2).eq.REFLECT_ODD) then
                  stylo(jmax+1) = zero
                  styhi(jmax+1) = zero
               end if

               do j=jmin,jmax+1
                  ystate(i,j,k) = cvmgp(stylo(j),styhi(j),vedge(i,j,k))
                  ystate(i,j,k) = cvmgt(half*(stylo(j)+styhi(j)),ystate(i,j,k),
     &                 abs(vedge(i,j,k)).lt.eps)
               end do
               place_to_break = 1
            end do
      end do
c
c     compute the zedge states
c
      do j = jmin,jmax
            do i = imin,imax
               
               do k = kmin-1,kmax+1
                  
                  spz = s(i,j,k) + half*sz(i,j,k)
                  smz = s(i,j,k) - half*sz(i,j,k)

                  if (iconserv.eq.1) then
                     tr =
     &                    (uedge(i+1,j,k)*xedge(i+1,j,k) - uedge(i,j,k)*xedge(i,j,k))/hx +   
     &                    (vedge(i,j+1,k)*yedge(i,j+1,k) - vedge(i,j,k)*yedge(i,j,k))/hy   
                     
                     st = -dth*tr + dth*(tforces(i,j,k) - s(i,j,k)*divu(i,j,k))
     &                    + dth*s(i,j,k)*(uedge(i+1,j,k)-uedge(i,j,k))/hx
     &                    + dth*s(i,j,k)*(vedge(i,j+1,k)-vedge(i,j,k))/hy
c    &                    - dth*s(i,j,k)*(wedge(i,j,k+1)-wedge(i,j,k))/hz

                  else
                     if (uedge(i,j,k)*uedge(i+1,j,k).le.0.d0) then
                        ubar = 0.5d0*(uedge(i,j,k)+uedge(i+1,j,k))
                        if (ubar.lt.0.d0) then
                           inc = 1
                        else
                           inc = 0
                        endif
                        tr1 = ubar*(s(i+inc,j,k)-s(i+inc-1,j,k))/hx
                     else
                        tr1 = half*(uedge(i+1,j,k) + uedge(i,j,k)) *
     &                       (xedge(i+1,j,k) - xedge(i,j,k)  ) / hx
                     endif
                     if (vedge(i,j,k)*vedge(i,j+1,k).lt.0.d0) then
                        vbar = 0.5d0*(vedge(i,j,k)+vedge(i,j+1,k))
                        if (vbar.lt.0.d0) then
                           inc = 1
                        else
                           inc = 0
                        endif
                        tr2 = vbar*(s(i,j+inc,k)-s(i,j+inc-1,k))/hy
                     else
                        tr2 = half*(vedge(i,j+1,k) + vedge(i,j,k)) *
     &                       (yedge(i,j+1,k) - yedge(i,j,k)  ) / hy
                     endif

                     st = -dth*(tr1 + tr2) + dth*tforces(i,j,k)
                  endif

                  stzlo(k+1)= spz - dthz*wedge(i,j,k+1)*sz(i,j,k) + st
                  stzhi(k  )= smz - dthz*wedge(i,j,k  )*sz(i,j,k) + st

               end do

               if (bc(3,1).eq.EXT_DIR .and. wedge(i,j,kmin).ge.zero) then
                  stzlo(kmin) = s(i,j,kmin-1)
                  stzhi(kmin) = s(i,j,kmin-1)
               else if (bc(3,1).eq.EXT_DIR .and. wedge(i,j,kmin).lt.zero) then
                  stzlo(kmin) = stzhi(kmin)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
               else if (bc(3,1).eq.FOEXTRAP.or.bc(3,1).eq.HOEXTRAP) then
                  if (n.eq.ZVEL) then
                     if (wedge(i,j,kmin).ge.zero) then
                        stzhi(kmin) = MIN(stzhi(kmin),zero)
                        stzlo(kmin) = stzhi(kmin)
                     endif
                  else
                     stzlo(kmin) = stzhi(kmin)
                  endif
               else if (bc(3,1).eq.REFLECT_EVEN) then
                  stzlo(kmin) = stzhi(kmin)
c               else if (bc(3,1).eq.FOEXTRAP.or.bc(3,1).eq.HOEXTRAP
c     &                 .or.bc(3,1).eq.REFLECT_EVEN) then
c                  stzlo(kmin) = stzhi(kmin)
               else if (bc(3,1).eq.REFLECT_ODD) then
                  stzlo(kmin) = zero
                  stzhi(kmin) = zero
               end if
               if (bc(3,2).eq.EXT_DIR .and. wedge(i,j,kmax+1).le.zero) then
                  stzlo(kmax+1) = s(i,j,kmax+1)
                  stzhi(kmax+1) = s(i,j,kmax+1)
               else if (bc(3,2).eq.EXT_DIR .and. wedge(i,j,kmax+1).gt.zero) then
                  stzhi(kmax+1) = stzlo(kmax+1)
c     HACK down to commented code below, to attempt to prevent backflow on outflow face
               else if (bc(3,2).eq.FOEXTRAP.or.bc(3,2).eq.HOEXTRAP) then
                  if (n.eq.ZVEL) then
                     if (wedge(i,j,kmax+1).le.zero) then
                        stzlo(kmax+1) = MAX(stzlo(kmax+1),zero)
                        stzhi(kmax+1) = stzlo(kmax+1)
                     endif
                  else
                     stzhi(kmax+1) = stzlo(kmax+1)
                  endif
               else if (bc(3,2).eq.REFLECT_EVEN) then
                  stzhi(kmax+1) = stzlo(kmax+1)
c               else if (bc(3,2).eq.FOEXTRAP.or.bc(3,2).eq.HOEXTRAP
c     &                 .or.bc(3,2).eq.REFLECT_EVEN) then
c                  stzhi(kmax+1) = stzlo(kmax+1)
               else if (bc(3,2).eq.REFLECT_ODD) then
                  stzlo(kmax+1) = zero
                  stzhi(kmax+1) = zero
               end if

               do k = kmin,kmax+1
                  zstate(i,j,k) = cvmgp(stzlo(k),stzhi(k),wedge(i,j,k))
                  zstate(i,j,k) = cvmgt(half*(stzlo(k)+stzhi(k)),zstate(i,j,k),
     &                 abs(wedge(i,j,k)).lt.eps)
               end do
               place_to_break = 1
            end do
      end do

      end if

      end

      subroutine FORT_ADV_FORCING(
     &     aofs,DIMS(aofs),
     &     xflux,DIMS(xflux),
     &     uedge,DIMS(uedge),
     &     areax,DIMS(ax),
     &     yflux,DIMS(yflux),
     &     vedge,DIMS(vedge),
     &     areay,DIMS(ay),
     &     zflux,DIMS(zflux),
     &     wedge,DIMS(wedge),
     &     areaz,DIMS(az),
     &     vol,DIMS(vol),
     &     lo,hi,iconserv )
c
c     This subroutine uses scalar edge states to compute
c     an advective tendency
c
      implicit none
      integer i,j,k
      integer iconserv
      REAL_T divux,divuy,divuz,rincr
      integer imin,jmin,kmin,imax,jmax,kmax
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(aofs)
      integer DIMDEC(vol)
      integer DIMDEC(uedge)
      integer DIMDEC(vedge)
      integer DIMDEC(wedge)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(zflux)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      integer DIMDEC(az)
      REAL_T aofs(DIMV(aofs))
      REAL_T vol(DIMV(vol))
      REAL_T uedge(DIMV(uedge))
      REAL_T vedge(DIMV(vedge))
      REAL_T wedge(DIMV(wedge))
      REAL_T xflux(DIMV(xflux))
      REAL_T yflux(DIMV(yflux))
      REAL_T zflux(DIMV(zflux))
      REAL_T areax(DIMV(ax))
      REAL_T areay(DIMV(ay))
      REAL_T areaz(DIMV(az))

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
c
c     if nonconservative initialize the advective tendency as -U*grad(S)
c
      if ( iconserv .ne. 1 ) then
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  divux = (
     &                 areax(i+1,j,k)*uedge(i+1,j,k)-
     &                 areax(i,  j,k)*uedge(i,  j,k))/vol(i,j,k)
                  divuy = (
     &                 areay(i,j+1,k)*vedge(i,j+1,k)-
     &                 areay(i,j,  k)*vedge(i,j,  k))/vol(i,j,k)
                  divuz = (
     &                 areaz(i,j,k+1)*wedge(i,j,k+1)-
     &                 areaz(i,j,k  )*wedge(i,j,k  ))/vol(i,j,k)
                  aofs(i,j,k) =
     &                 - divux*half*(xflux(i+1,j,k)+xflux(i,j,k))
     &                 - divuy*half*(yflux(i,j+1,k)+yflux(i,j,k))
     &                 - divuz*half*(zflux(i,j,k+1)+zflux(i,j,k))
              
               end do
            end do
         end do
      end if
c
c     convert edge states to fluxes
c
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax+1
               xflux(i,j,k) = xflux(i,j,k)*uedge(i,j,k)*areax(i,j,k)
            end do
         end do
      end do
      do k = kmin,kmax
         do j = jmin,jmax+1
            do i = imin,imax
               yflux(i,j,k) = yflux(i,j,k)*vedge(i,j,k)*areay(i,j,k)
            end do
         end do
      end do
      do k = kmin,kmax+1
         do j = jmin,jmax
            do i = imin,imax
               zflux(i,j,k) = zflux(i,j,k)*wedge(i,j,k)*areaz(i,j,k)
            end do
         end do
      end do
c
c     compute the part of the advective tendency 
c     that depends on the flux convergence
c
      if ( iconserv .ne. 1 ) then
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  aofs(i,j,k) = aofs(i,j,k) + (
     &                 xflux(i+1,j,k) - xflux(i,j,k) +
     &                 yflux(i,j+1,k) - yflux(i,j,k) +
     &                 zflux(i,j,k+1) - zflux(i,j,k))/vol(i,j,k)
               end do
            end do
         end do
      else
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  aofs(i,j,k) = (
     &                 xflux(i+1,j,k) - xflux(i,j,k) +
     &                 yflux(i,j+1,k) - yflux(i,j,k) +
     &                 zflux(i,j,k+1) - zflux(i,j,k))/vol(i,j,k)
               end do
            end do
         end do
      end if

      end

      subroutine FORT_SYNC_ADV_FORCING(
     &     sync ,DIMS(sync),
     &     xflux,DIMS(xflux),
     &     ucor ,DIMS(ucor),
     &     areax,DIMS(ax),
     &     yflux,DIMS(yflux),
     &     vcor ,DIMS(vcor),
     &     areay,DIMS(ay),
     &     zflux,DIMS(zflux),
     &     wcor ,DIMS(wcor),
     &     areaz,DIMS(az),
     &     vol ,DIMS(vol),
     &     lo,hi,iconserv )
c
c     This subroutine computes the sync advective tendency
c     for a state variable
c
      implicit none
      integer i,j,k
      integer iconserv
      REAL_T divux,divuy,divuz
      integer imin,jmin,kmin,imax,jmax,kmax
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(sync)
      integer DIMDEC(vol)
      integer DIMDEC(ucor)
      integer DIMDEC(vcor)
      integer DIMDEC(wcor)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(zflux)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      integer DIMDEC(az)
      REAL_T sync(DIMV(sync))
      REAL_T vol(DIMV(vol))
      REAL_T ucor(DIMV(ucor))
      REAL_T vcor(DIMV(vcor))
      REAL_T wcor(DIMV(wcor))
      REAL_T xflux(DIMV(xflux))
      REAL_T yflux(DIMV(yflux))
      REAL_T zflux(DIMV(zflux))
      REAL_T areax(DIMV(ax))
      REAL_T areay(DIMV(ay))
      REAL_T areaz(DIMV(az))

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
c
c     compute corrective fluxes from edge states 
c     and perform conservative update
c
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax+1
               xflux(i,j,k) = xflux(i,j,k)*ucor(i,j,k)*areax(i,j,k)
            end do
         end do
      end do
      do k = kmin,kmax
         do j = jmin,jmax+1
            do i = imin,imax
               yflux(i,j,k) = yflux(i,j,k)*vcor(i,j,k)*areay(i,j,k)
            end do
         end do
      end do
      do k = kmin,kmax+1
         do j = jmin,jmax
            do i = imin,imax
               zflux(i,j,k) = zflux(i,j,k)*wcor(i,j,k)*areaz(i,j,k)
            end do
         end do
      end do

      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               sync(i,j,k) = sync(i,j,k) + (
     &              xflux(i+1,j,k)-xflux(i,j,k) +
     &              yflux(i,j+1,k)-yflux(i,j,k) +
     &              zflux(i,j,k+1)-zflux(i,j,k) )/vol(i,j,k)
            end do
         end do
      end do

      end

      subroutine trans_xbc(
     &     s,DIMS(s),
     &     xlo,xhi,DIMS(xx),uad,DIMS(uad),
     &     lo,hi,n,xbc,eps,ycouple,zcouple)
c
c     This subroutine processes boundary conditions on information
c     traced to cell faces in the x direction.  This is used for
c     computing velocities and edge states used in calculating
c     transverse derivatives
c
      implicit none
      integer DIMDEC(s)
      REAL_T s(DIMV(s))
      integer DIMDEC(xx)
      integer DIMDEC(uad)
      REAL_T xlo(DIMV(xx))
      REAL_T xhi(DIMV(xx))
      REAL_T uad(DIMV(uad))
      REAL_T eps
      logical ycouple
      logical zcouple
      integer lo(SDIM), hi(SDIM)
      integer n
      integer xbc(SDIM,2)

      REAL_T stx
      logical ltest
      integer j,k
      integer imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

c     if we are applying bc's to intermediate corner-copuled terms
c     the bounds are slightly different on the valid data
      if (ycouple) then
         jmin = jmin + 1
         jmax = jmax - 1
      end if

      if (zcouple) then
         kmin = kmin + 1
         kmax = kmax - 1
      end if
c
c     -------------- the lower x boundary
c
      if (xbc(1,1).eq.EXT_DIR) then
         if ( n .eq. XVEL ) then
            do j = jmin-1,jmax+1
             do k = kmin-1,kmax+1
              if (uad(imin,j,k) .ge. zero) then
                  xlo(imin,j,k) = s(imin-1,j,k)
                  xhi(imin,j,k) = s(imin-1,j,k)
              else
                  xlo(imin,j,k) = xhi(imin,j,k)
              endif
             end do
            end do
         else
            do j = jmin-1,jmax+1
               do k = kmin-1,kmax+1
                  ltest = uad(imin,j,k).le.eps
                  stx   = cvmgt(xhi(imin,j,k),s(imin-1,j,k),ltest)
                  xlo(imin,j,k) = stx
                  xhi(imin,j,k) = stx
               end do
            end do
         end if
      else if (xbc(1,1).eq.FOEXTRAP.or.xbc(1,1).eq.HOEXTRAP
     &        .or.xbc(1,1).eq.REFLECT_EVEN) then
         do j = jmin-1,jmax+1
            do k = kmin-1,kmax+1
               xlo(imin,j,k) = xhi(imin,j,k)
            end do
         end do
      else if (xbc(1,1).eq.REFLECT_ODD) then
         do j = jmin-1,jmax+1
            do k = kmin-1,kmax+1
               xhi(imin,j,k) = zero
               xlo(imin,j,k) = zero
            end do
         end do
      end if
c
c     -------------- the upper x boundary
c
      if (xbc(1,2).eq.EXT_DIR) then
         if ( n .eq. XVEL ) then
            do j = jmin-1,jmax+1
             do k = kmin-1,kmax+1
               if (uad(imax+1,j,k) .le. zero) then
                  xlo(imax+1,j,k) = s(imax+1,j,k)
                  xhi(imax+1,j,k) = s(imax+1,j,k)
               else
                  xhi(imax+1,j,k) = xlo(imax+1,j,k)
               endif
             end do
            end do
         else
            do j = jmin-1,jmax+1
               do k = kmin-1,kmax+1
                  ltest = uad(imax+1,j,k).ge.-eps
                  stx   = cvmgt(xlo(imax+1,j,k),s(imax+1,j,k),ltest)
                  xlo(imax+1,j,k) = stx
                  xhi(imax+1,j,k) = stx
               end do
            end do
         end if
      else if (xbc(1,2).eq.FOEXTRAP.or.xbc(1,2).eq.HOEXTRAP
     &        .or.xbc(1,2).eq.REFLECT_EVEN) then
         do j = jmin-1,jmax+1
            do k = kmin-1,kmax+1
               xhi(imax+1,j,k) = xlo(imax+1,j,k)
            end do
         end do
      else if (xbc(1,2).eq.REFLECT_ODD) then
         do j = jmin-1,jmax+1
            do k = kmin-1,kmax+1
               xhi(imax+1,j,k) = zero
               xlo(imax+1,j,k) = zero
            end do
         end do
      end if

      end

      subroutine trans_ybc(
     &     s,DIMS(s),
     &     ylo,yhi,DIMS(yy),vad,DIMS(vad),
     &     lo,hi,n,ybc,eps,xcouple,zcouple)
c
c     This subroutine processes boundary conditions on information
c     traced to cell faces in the y direction.  This is used for
c     computing velocities and edge states used in calculating
c     transverse derivatives
c
      implicit none
      integer DIMDEC(s)
      REAL_T s(DIMV(s))
      integer DIMDEC(yy)
      integer DIMDEC(vad)
      REAL_T ylo(DIMV(yy))
      REAL_T yhi(DIMV(yy))
      REAL_T vad(DIMV(vad))
      REAL_T eps
      logical xcouple
      logical zcouple
      integer lo(SDIM), hi(SDIM)
      integer n
      integer ybc(SDIM,2)

      REAL_T sty
      logical ltest
      integer i,k
      integer imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

c     if we are applying bc's to intermediate corner-copuled terms
c     the bounds are slightly different on the valid data
      if (xcouple) then
         imin = imin + 1
         imax = imax - 1
      end if

      if (zcouple) then
         kmin = kmin + 1
         kmax = kmax - 1
      end if

c
c     -------------- the lower y boundary
c
      if (ybc(2,1).eq.EXT_DIR) then
         if ( n .eq. YVEL ) then
            do i = imin-1,imax+1
             do k = kmin-1,kmax+1
              if (vad(i,jmin,k) .ge. zero) then
                  ylo(i,jmin,k) = s(i,jmin-1,k)
                  yhi(i,jmin,k) = s(i,jmin-1,k)
              else
                  ylo(i,jmin,k) = yhi(i,jmin,k)
              endif
             end do
            end do
         else
            do i = imin-1,imax+1
               do k = kmin-1,kmax+1
                  ltest = vad(i,jmin,k).le.eps
                  sty   = cvmgt(yhi(i,jmin,k),s(i,jmin-1,k),ltest)
                  ylo(i,jmin,k) = sty
                  yhi(i,jmin,k) = sty
               end do
            end do
         end if
      else if (ybc(2,1).eq.FOEXTRAP.or.ybc(2,1).eq.HOEXTRAP
     &        .or.ybc(2,1).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
            do k = kmin-1,kmax+1
               ylo(i,jmin,k) = yhi(i,jmin,k)
            end do
         end do
      else if (ybc(2,1).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
            do k = kmin-1,kmax+1
               yhi(i,jmin,k) = zero
               ylo(i,jmin,k) = zero
            end do
         end do
      end if
c
c     -------------- the upper y boundary
c
      if (ybc(2,2).eq.EXT_DIR) then
         if ( n .eq. YVEL ) then
            do i = imin-1,imax+1
             do k = kmin-1,kmax+1
               if (vad(i,jmax+1,k) .le. zero) then
                  ylo(i,jmax+1,k) = s(i,jmax+1,k)
                  yhi(i,jmax+1,k) = s(i,jmax+1,k)
               else
                  yhi(i,jmax+1,k) = ylo(i,jmax+1,k)
               endif
             end do
            end do
         else
            do i = imin-1,imax+1
               do k = kmin-1,kmax+1
                  ltest = vad(i,jmax+1,k).ge.-eps
                  sty   = cvmgt(ylo(i,jmax+1,k),s(i,jmax+1,k),ltest)
                  ylo(i,jmax+1,k) = sty
                  yhi(i,jmax+1,k) = sty
               end do
            end do
         end if
      else if (ybc(2,2).eq.FOEXTRAP.or.ybc(2,2).eq.HOEXTRAP
     &        .or.ybc(2,2).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
            do k = kmin-1,kmax+1
               yhi(i,jmax+1,k) = ylo(i,jmax+1,k)
            end do
         end do
      else if (ybc(2,2).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
            do k = kmin-1,kmax+1
               ylo(i,jmax+1,k) = zero
               yhi(i,jmax+1,k) = zero
            end do
         end do
      end if

      end

      subroutine trans_zbc(
     &     s,DIMS(s),
     &     zlo,zhi,DIMS(zz),wad,DIMS(wad),
     &     lo,hi,n,zbc,eps,xcouple,ycouple)
c
c     This subroutine processes boundary conditions on information
c     traced to cell faces in the z direction.  This is used for
c     computing velocities and edge states used in calculating
c     transverse derivatives
c
      implicit none
      integer DIMDEC(s)
      REAL_T s(DIMV(s))
      integer DIMDEC(zz)
      integer DIMDEC(wad)
      REAL_T zlo(DIMV(zz))
      REAL_T zhi(DIMV(zz))
      REAL_T wad(DIMV(wad))
      REAL_T eps
      logical xcouple
      logical ycouple
      integer lo(SDIM), hi(SDIM)
      integer n
      integer zbc(SDIM,2)

      REAL_T stz
      logical ltest
      integer i,j
      integer imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

c     if we are applying bc's to intermediate corner-copuled terms
c     the bounds are slightly different on the valid data
      if (xcouple) then
         imin = imin + 1
         imax = imax - 1
      end if

      if (ycouple) then
         jmin = jmin + 1
         jmax = jmax - 1
      end if

c
c     -------------- the lower z boundary
c
      if (zbc(3,1).eq.EXT_DIR) then
         if ( n .eq. ZVEL ) then
            do i = imin-1,imax+1
             do j = jmin-1,jmax+1
               if (wad(i,j,kmin) .ge. zero) then
                  zhi(i,j,kmin) = s(i,j,kmin-1)
                  zlo(i,j,kmin) = s(i,j,kmin-1)
               else
                  zlo(i,j,kmin) = zhi(i,j,kmin)
               endif
             end do
            end do
         else
            do i = imin-1,imax+1
               do j = jmin-1,jmax+1
                  ltest = wad(i,j,kmin).le.eps
                  stz   = cvmgt(zhi(i,j,kmin),s(i,j,kmin-1),ltest)
                  zhi(i,j,kmin) = stz
                  zlo(i,j,kmin) = stz
               end do
            end do
         end if
      else if (zbc(3,1).eq.FOEXTRAP.or.zbc(3,1).eq.HOEXTRAP
     &        .or.zbc(3,1).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
            do j = jmin-1,jmax+1
               zlo(i,j,kmin) = zhi(i,j,kmin)
            end do
         end do
      else if (zbc(3,1).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
            do j = jmin-1,jmax+1
               zhi(i,j,kmin) = zero
               zlo(i,j,kmin) = zero
            end do
         end do
      end if
c
c     -------------- the upper z boundary
c
      if (zbc(3,2).eq.EXT_DIR) then
         if ( n .eq. ZVEL ) then
            do i = imin-1,imax+1
             do j = jmin-1,jmax+1
               if (wad(i,j,kmax+1) .le. zero) then
                  zlo(i,j,kmax+1) = s(i,j,kmax+1)
                  zhi(i,j,kmax+1) = s(i,j,kmax+1)
               else
                  zhi(i,j,kmax+1) = zlo(i,j,kmax+1)
               endif
             end do
            end do
         else
            do i = imin-1,imax+1
               do j = jmin-1,jmax+1
                  ltest = wad(i,j,kmax+1).ge.-eps
                  stz   = cvmgt(zlo(i,j,kmax+1),s(i,j,kmax+1),ltest)
                  zhi(i,j,kmax+1) = stz
                  zlo(i,j,kmax+1) = stz
               end do
            end do
         end if
      else if (zbc(3,2).eq.FOEXTRAP.or.zbc(3,2).eq.HOEXTRAP
     &        .or.zbc(3,2).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
            do j = jmin-1,jmax+1
               zhi(i,j,kmax+1) = zlo(i,j,kmax+1)
            end do
         end do
      else if (zbc(3,2).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
            do j = jmin-1,jmax+1
               zlo(i,j,kmax+1) = zero
               zhi(i,j,kmax+1) = zero
            end do
         end do
      end if

      end

      subroutine FORT_SLOPES( dir,
     &     s,DIMS(s),
     &     slx,sly,slz,DIMS(sl),
     &     lo,hi,slxscr,slyscr,slzscr,bc)
c 
c     this subroutine computes first or forth order slopes of
c     a 3D scalar field.
c
c     (dir) is used to eliminate calculating extra slopes in transvel
c
c     Boundary conditions on interior slopes are handled automatically
c     by the ghost cells
c
c     Boundary conditions on EXT_DIR and HOEXTRAP slopes are implemented
c     by setting them to zero outside of the domain and using a
c     one-sided derivative from the interior
c
      implicit none

#include "GODCOMM_F.H"

      integer dir
      integer DIMDEC(s)
      REAL_T     s(DIMV(s))
      integer DIMDEC(sl)
      REAL_T   slx(DIMV(sl))
      REAL_T   sly(DIMV(sl))
      REAL_T   slz(DIMV(sl))
      integer lo(SDIM), hi(SDIM)
      REAL_T slxscr(DIM1(s), 4)
      REAL_T slyscr(DIM2(s), 4)
      REAL_T slzscr(DIM3(s), 4)
      integer bc(SDIM,2)

      integer imin,jmin,kmin,imax,jmax,kmax,i,j,k
      integer ng
      REAL_T dpls,dmin,ds
      REAL_T del,slim,sflg
      integer cen,lim,flag,fromm
      parameter( cen = 1   )
      parameter( lim = 2   )
      parameter( flag = 3  )
      parameter( fromm = 4 )

      ng = lo(1) - ARG_L1(s)
      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
c
c     Added to prevent underflow for small s values.
c
      do k = lo(3)-ng, hi(3)+ng
         do j = lo(2)-ng, hi(2)+ng
            do i = lo(1)-ng, hi(1)+ng
               s(i,j,k) = cvmgt(s(i,j,k), zero, abs(s(i,j,k)).gt.1.0D-20)
            end do
         end do
      end do
c
c     COMPUTE 0TH order slopes
c
      if (slope_order.eq.1) then
        if (ng .lt. 1) then
	   call bl_abort("FORT_SLOPES: not enough bndry cells for 1st order")
        end if
        do k = kmin-1, kmax+1
           do j = jmin-1, jmax+1 
              do i = imin-1, imax+1
                 slx(i,j,k) = zero
                 sly(i,j,k) = zero
                 slz(i,j,k) = zero
              end do
           end do
        end do
        return
      end if

c
c     COMPUTE 2ND order slopes
c
      if (slope_order.eq.2) then
        if (ng .lt. 2) then
	   call bl_abort("FORT_SLOPES: not enough bndry cells for 2nd order")
        end if
c
c     ------------------------ x slopes
c
        if ( (dir.eq.XVEL) .or. (dir.eq.ALL) ) then
         if (use_unlimited_slopes) then
            do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  slx(i,j,k)= half*(s(i+1,j,k) - s(i-1,j,k))
               end do
            end do
            end do
            if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
               do k = kmin-1, kmax+1
               do j = jmin-1, jmax+1
                  slx(imin-1,j,k) = zero
                  slx(imin  ,j,k) = (s(imin+1,j,k)+three*s(imin,j,k)-four*s(imin-1,j,k))/three
               end do
               end do
            end if
            if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
               do k = kmin-1, kmax+1
               do j = jmin-1, jmax+1
                  slx(imax+1,j,k) = zero
                  slx(imax  ,j,k) = -(s(imax-1,j,k)+three*s(imax,j,k)-four*s(imax+1,j,k))/three
               end do
               end do
            end if
         else
            do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  del  = half*(s(i+1,j,k) - s(i-1,j,k))
                  dpls =  two*(s(i+1,j,k) - s(i  ,j,k))
                  dmin =  two*(s(i  ,j,k) - s(i-1,j,k))
                  slim = min(abs(dpls), abs(dmin))
                  slim = cvmgp(slim, zero, dpls*dmin)
                  sflg = sign(one,del)
                  slx(i,j,k)= sflg*min(slim,abs(del))
               end do
            end do
            end do
            
            if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
               do k = kmin-1, kmax+1
               do j = jmin-1, jmax+1
                  slx(imin-1,j,k) = zero
                  del  = (s(imin+1,j,k)+three*s(imin,j,k)-four*s(imin-1,j,k))/three
                  dpls = two*(s(imin+1,j,k) - s(imin  ,j,k))
                  dmin = two*(s(imin  ,j,k) - s(imin-1,j,k))
                  slim = min(abs(dpls), abs(dmin))
                  slim = cvmgp(slim, zero, dpls*dmin)
                  sflg = sign(one,del)
                  slx(imin,j,k)= sflg*min(slim,abs(del))
               end do
               end do
            end if
            if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
               do k = kmin-1, kmax+1
               do j = jmin-1, jmax+1
                  slx(imax+1,j,k) = zero
                  del  = -(s(imax-1,j,k)+three*s(imax,j,k)-four*s(imax+1,j,k))/three
                  dpls = two*(s(imax+1,j,k) - s(imax  ,j,k))
                  dmin = two*(s(imax  ,j,k) - s(imax-1,j,k))
                  slim = min(abs(dpls), abs(dmin))
                  slim = cvmgp(slim, zero, dpls*dmin)
                  sflg = sign(one,del)
                  slx(imax,j,k)= sflg*min(slim,abs(del))
               end do
               end do
            end if
          end if
         end if
c
c     ------------------------ y slopes
c
        if ( (dir.eq.YVEL) .or. (dir.eq.ALL) ) then
         if (use_unlimited_slopes) then
            do k = kmin-1,kmax+1
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  sly(i,j,k) = half*(s(i,j+1,k)-s(i,j-1,k))
               end do
            end do
            end do
            if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
               do k = kmin-1, kmax+1
               do i = imin-1, imax+1
                  sly(i,jmin-1,k) = zero
                  sly(i,jmin  ,k) = (s(i,jmin+1,k)+three*s(i,jmin,k)-four*s(i,jmin-1,k))/three
               end do
               end do
            end if
            if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
               do k = kmin-1, kmax+1
               do i = imin-1, imax+1
                  sly(i,jmax+1,k) = zero
                  sly(i,jmax  ,k) = -(s(i,jmax-1,k)+three*s(i,jmax,k)-four*s(i,jmax+1,k))/three
               end do
               end do
            end if
         else
            do j = jmin-1,jmax+1
               do k = kmin-1,kmax+1
               do i = imin-1,imax+1
                  del  = half*(s(i,j+1,k) - s(i,j-1,k))
                  dpls =  two*(s(i,j+1,k) - s(i,j  ,k))
                  dmin =  two*(s(i,j  ,k) - s(i,j-1,k))
                  slim = min(abs(dpls),abs(dmin))
                  slim = cvmgp(slim, zero, dpls*dmin)
                  sflg = sign(one,del)
                  sly(i,j,k)= sflg*min(slim,abs(del))
               end do
               end do
            end do

            if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
               do k = kmin-1, kmax+1
               do i = imin-1, imax+1
                  sly(i,jmin-1,k) = zero
                  del  = (s(i,jmin+1,k)+three*s(i,jmin,k)-four*s(i,jmin-1,k))/three
                  dpls = two*(s(i,jmin+1,k) - s(i,jmin  ,k))
                  dmin = two*(s(i,jmin  ,k) - s(i,jmin-1,k))
                  slim = min(abs(dpls), abs(dmin))
                  slim = cvmgp(slim, zero, dpls*dmin)
                  sflg = sign(one,del)
                  sly(i,jmin,k)= sflg*min(slim,abs(del))
               end do
               end do
            end if
            if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
               do k = kmin-1, kmax+1
               do i = imin-1, imax+1
                  sly(i,jmax+1,k) = zero
                  del  = -(s(i,jmax-1,k)+three*s(i,jmax,k)-four*s(i,jmax+1,k))/three
                  dpls = two*(s(i,jmax+1,k) - s(i,jmax  ,k))
                  dmin = two*(s(i,jmax  ,k) - s(i,jmax-1,k))
                  slim = min(abs(dpls), abs(dmin))
                  slim = cvmgp(slim, zero, dpls*dmin)
                  sflg = sign(one,del)
                  sly(i,jmax,k)= sflg*min(slim,abs(del))
               end do
               end do
            end if
         end if
        end if

c
c     ------------------------ z slopes
c
        if ( (dir.eq.ZVEL) .or. (dir.eq.ALL) ) then
         if (use_unlimited_slopes) then
            do j = jmin-1,jmax+1
            do i = imin-1,imax+1
               do k = kmin-1,kmax+1
                  slz(i,j,k) = half*(s(i,j,k+1)-s(i,j,k-1))
               end do
            end do
            end do
            if (bc(3,1) .eq. EXT_DIR .or. bc(3,1) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
               do i = imin-1, imax+1
                  slz(i,j,kmin-1) = zero
                  slz(i,j,kmin  ) = (s(i,j,kmin+1)+three*s(i,j,kmin)-four*s(i,j,kmin-1))/three
               end do
               end do
            end if
            if (bc(3,2) .eq. EXT_DIR .or. bc(3,2) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
               do i = imin-1, imax+1
                  slz(i,j,kmax+1) = zero
                  slz(i,j,kmax  ) = -(s(i,j,kmax-1)+three*s(i,j,kmax)-four*s(i,j,kmax+1))/three
               end do
               end do
            end if
         else
            do k = kmin-1,kmax+1
               do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  del  = half*(s(i,j,k+1) - s(i,j,k-1))
                  dpls =  two*(s(i,j,k+1) - s(i,j,k  ))
                  dmin =  two*(s(i,j,k  ) - s(i,j,k-1))
                  slim = min(abs(dpls),abs(dmin))
                  slim = cvmgp(slim, zero, dpls*dmin)
                  sflg = sign(one,del)
                  slz(i,j,k)= sflg*min(slim,abs(del))
               end do
               end do
            end do

            if (bc(3,1) .eq. EXT_DIR .or. bc(3,1) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
               do i = imin-1, imax+1
                  slz(i,j,kmin-1) = zero
                  del  = (s(i,j,kmin+1)+three*s(i,j,kmin)-four*s(i,j,kmin-1))/three
                  dpls = two*(s(i,j,kmin+1) - s(i,j,kmin  ))
                  dmin = two*(s(i,j,kmin  ) - s(i,j,kmin-1))
                  slim = min(abs(dpls), abs(dmin))
                  slim = cvmgp(slim, zero, dpls*dmin)
                  sflg = sign(one,del)
                  slz(i,j,kmin)= sflg*min(slim,abs(del))
               end do
               end do
            end if
            if (bc(3,2) .eq. EXT_DIR .or. bc(3,2) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
               do i = imin-1, imax+1
                  slz(i,j,kmax+1) = zero
                  del  = -(s(i,j,kmax-1)+three*s(i,j,kmax)-four*s(i,j,kmax+1))/three
                  dpls = two*(s(i,j,kmax+1) - s(i,j,kmax  ))
                  dmin = two*(s(i,j,kmax  ) - s(i,j,kmax-1))
                  slim = min(abs(dpls), abs(dmin))
                  slim = cvmgp(slim, zero, dpls*dmin)
                  sflg = sign(one,del)
                  slz(i,j,kmax)= sflg*min(slim,abs(del))
               end do
               end do
            end if
         end if
        end if
c
c ... end, if slope_order .eq. 2
c
      end if
c
c     COMPUTE 4TH order slopes
c
      if (slope_order.eq.4)then
        if (ng .lt. 3) then
	   call bl_abort("SLOPE_3D: not enough bndry cells for 4th order")
        end if
c
c     ------------------------ x slopes
c
        if ( (dir.eq.XVEL) .or. (dir.eq.ALL) ) then
           if (use_unlimited_slopes) then
              do k = kmin-1,kmax+1
                 do j = jmin-1,jmax+1
                    do i = imin-2,imax+2
                       slxscr(i,cen)  = half*(s(i+1,j,k)-s(i-1,j,k))
                    end do
                    do i = imin-1,imax+1
                       slx(i,j,k) = two * two3rd * slxscr(i,cen) -
     &                      sixth * (slxscr(i+1,cen) + slxscr(i-1,cen))
                    end do
                 end do
              end do
              
              if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
                 do k = kmin-1,kmax+1
                    do j = jmin-1, jmax+1
                       slx(imin,j,k) = -sixteen/fifteen*s(imin-1,j,k) + half*s(imin,j,k) + 
     &                      two3rd*s(imin+1,j,k) - tenth*s(imin+2,j,k)
                       slx(imin-1,j,k) = zero
                    end do
                 end do
              end if
              if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
                 do k = kmin-1,kmax+1
                    do j = jmin-1, jmax+1
                       slx(imax,j,k) = -( -sixteen/fifteen*s(imax+1,j,k) + half*s(imax,j,k) + 
     &                      two3rd*s(imax-1,j,k) - tenth*s(imax-2,j,k) )
                       slx(imax+1,j,k) = zero
                    end do
                 end do
              end if
           else
              do k = kmin-1,kmax+1
                 do j = jmin-1,jmax+1 
                    do i = imin-2,imax+2
                       dmin           =  two*(s(i,  j,k)-s(i-1,j,k))
                       dpls           =  two*(s(i+1,j,k)-s(i,  j,k))
                       slxscr(i,cen)  = half*(s(i+1,j,k)-s(i-1,j,k))
                       slxscr(i,lim)  = min(abs(dmin),abs(dpls))
                       slxscr(i,lim)  = cvmgp(slxscr(i,lim),zero,dpls*dmin)
                       slxscr(i,flag) = sign(one,slxscr(i,cen))
                       slxscr(i,fromm)= slxscr(i,flag)*
     &                      min(slxscr(i,lim),abs(slxscr(i,cen)))
                    end do
                    do i = imin-1,imax+1
                       ds = two * two3rd * slxscr(i,cen) - 
     &                      sixth * (slxscr(i+1,fromm) + slxscr(i-1,fromm))
                       slx(i,j,k) = slxscr(i,flag)*min(abs(ds),slxscr(i,lim))
                    end do

                    if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
                       del  = -sixteen/fifteen*s(imin-1,j,k) + half*s(imin,j,k) + 
     &                      two3rd*s(imin+1,j,k) -  tenth*s(imin+2,j,k)
                       dmin = two*(s(imin  ,j,k)-s(imin-1,j,k))
                       dpls = two*(s(imin+1,j,k)-s(imin  ,j,k))
                       slim = min(abs(dpls), abs(dmin))
                       slim = cvmgp(slim, zero, dpls*dmin)
                       sflg = sign(one,del)
                       slx(imin-1,j,k) = zero
                       slx(imin,  j,k) = sflg*min(slim,abs(del))

c                      Recalculate the slope at imin+1 using the revised slxscr(imin,fromm)
                       slxscr(imin,fromm) = slx(imin,j,k)
                       ds = two * two3rd * slxscr(imin+1,cen) -
     $                    sixth * (slxscr(imin+2,fromm) + slxscr(imin,fromm))
                       slx(imin+1,j,k) = slxscr(imin+1,flag)*min(abs(ds),slxscr(imin+1,lim))
                    end if

                    if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
                       del  = -( -sixteen/fifteen*s(imax+1,j,k) + half*s(imax,j,k) + 
     &                      two3rd*s(imax-1,j,k) - tenth*s(imax-2,j,k) )
                       dmin = two*(s(imax  ,j,k)-s(imax-1,j,k))
                       dpls = two*(s(imax+1,j,k)-s(imax  ,j,k))
                       slim = min(abs(dpls), abs(dmin))
                       slim = cvmgp(slim, zero, dpls*dmin)
                       sflg = sign(one,del)
                       slx(imax,  j,k) = sflg*min(slim,abs(del))
                       slx(imax+1,j,k) = zero

c                      Recalculate the slope at imax-1 using the revised slxscr(imax,fromm)
                       slxscr(imax,fromm) = slx(imax,j,k)
                       ds = two * two3rd * slxscr(imax-1,cen) -
     $                    sixth * (slxscr(imax-2,fromm) + slxscr(imax,fromm))
                       slx(imax-1,j,k) = slxscr(imax-1,flag)*min(abs(ds),slxscr(imax-1,lim))
                    end if
                 end do
              end do
           end if
        end if
c
c     ------------------------ y slopes
c
        if ( (dir.eq.YVEL) .or. (dir.eq.ALL) ) then
           if (use_unlimited_slopes) then
              do k = kmin-1,kmax+1
                 do i = imin-1,imax+1
                    do j = jmin-2,jmax+2
                       slyscr(j,cen)  = half*(s(i,j+1,k)-s(i,j-1,k))
                    end do
                    do j = jmin-1,jmax+1
                       sly(i,j,k) = two * two3rd * slyscr(j,cen) -
     &                      sixth * (slyscr(j+1,cen) + slyscr(j-1,cen))
                    end do
                 end do
              end do
              
              if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
                 do k = kmin-1,kmax+1
                    do i = imin-1, imax+1
                       sly(i,jmin-1,k) = zero
                       sly(i,jmin,k) = -sixteen/fifteen*s(i,jmin-1,k) + half*s(i,jmin,k) + 
     &                      two3rd*s(i,jmin+1,k) - tenth*s(i,jmin+2,k)
                    end do
                 end do
              end if
              if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
                 do k = kmin-1,kmax+1
                    do i = imin-1, imax+1
                       sly(i,jmax,k) = -( -sixteen/fifteen*s(i,jmax+1,k) + half*s(i,jmax,k) + 
     &                      two3rd*s(i,jmax-1,k) - tenth*s(i,jmax-2,k) )
                       sly(i,jmax+1,k) = zero
                    end do
                 end do
              end if
           else
              do k = kmin-1,kmax+1
                 do i = imin-1,imax+1 
                    do j = jmin-2,jmax+2
                       dmin           =  two*(s(i,j,  k)-s(i,j-1,k))
                       dpls           =  two*(s(i,j+1,k)-s(i,j,  k))
                       slyscr(j,cen)  = half*(s(i,j+1,k)-s(i,j-1,k))
                       slyscr(j,lim)  = min(abs(dmin),abs(dpls))
                       slyscr(j,lim)  = cvmgp(slyscr(j,lim),zero,dpls*dmin)
                       slyscr(j,flag) = sign(one,slyscr(j,cen))
                       slyscr(j,fromm)= slyscr(j,flag)*
     &                      min(slyscr(j,lim),abs(slyscr(j,cen)))
                    end do
                    do j = jmin-1,jmax+1
                       ds = two * two3rd * slyscr(j,cen) - 
     &                      sixth * (slyscr(j+1,fromm) + slyscr(j-1,fromm))
                       sly(i,j,k) = slyscr(j,flag)*min(abs(ds),slyscr(j,lim))
                    end do
c
                    if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
                       del  = -sixteen/fifteen*s(i,jmin-1,k) + half*s(i,jmin,k) + 
     &                      two3rd*s(i,jmin+1,k) - tenth*s(i,jmin+2,k)
                       dmin = two*(s(i,jmin  ,k)-s(i,jmin-1,k))
                       dpls = two*(s(i,jmin+1,k)-s(i,jmin  ,k))
                       slim = min(abs(dpls), abs(dmin))
                       slim = cvmgp(slim, zero, dpls*dmin)
                       sflg = sign(one,del)
                       sly(i,jmin-1,k) = zero
                       sly(i,jmin,  k) = sflg*min(slim,abs(del))

c                      Recalculate the slope at jmin+1 using the revised slyscr(jmin,fromm)
                       slyscr(jmin,fromm) = sly(i,jmin,k)
                       ds = two * two3rd * slyscr(jmin+1,cen) -
     $                    sixth * (slyscr(jmin+2,fromm) + slyscr(jmin,fromm))
                       sly(i,jmin+1,k) = slyscr(jmin+1,flag)*min(abs(ds),slyscr(jmin+1,lim))
                    end if
                    if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
                       del  = -( -sixteen/fifteen*s(i,jmax+1,k) + half*s(i,jmax,k) +
     &                      two3rd*s(i,jmax-1,k) - tenth*s(i,jmax-2,k) )
                       dmin = two*(s(i,jmax  ,k)-s(i,jmax-1,k))
                       dpls = two*(s(i,jmax+1,k)-s(i,jmax  ,k))
                       slim = min(abs(dpls), abs(dmin))
                       slim = cvmgp(slim, zero, dpls*dmin)
                       sflg = sign(one,del)
                       sly(i,jmax, k)  = sflg*min(slim,abs(del))
                       sly(i,jmax+1,k) = zero

c                      Recalculate the slope at jmax-1 using the revised slyscr(jmax,fromm)
                       slyscr(jmax,fromm) = sly(i,jmax,k)
                       ds = two * two3rd * slyscr(jmax-1,cen) -
     $                    sixth * (slyscr(jmax-2,fromm) + slyscr(jmax,fromm))
                       sly(i,jmax-1,k) = slyscr(jmax-1,flag)*min(abs(ds),slyscr(jmax-1,lim))
                    end if
                 end do
              end do
           end if
        end if
c
c     ------------------------ z slopes
c
        if ( (dir.eq.ZVEL) .or. (dir.eq.ALL) ) then
           if (use_unlimited_slopes) then
              do j = jmin-1,jmax+1
                 do i = imin-1,imax+1
                    do k = kmin-2,kmax+2
                       slzscr(k,cen)  = half*(s(i,j,k+1)-s(i,j,k-1))
                    end do
                    do k = kmin-1,kmax+1
                       slz(i,j,k) = two * two3rd * slzscr(k,cen) -
     &                      sixth * (slzscr(k+1,cen) + slzscr(k-1,cen))
                    end do
                 end do
              end do
              
              if (bc(3,1) .eq. EXT_DIR .or. bc(3,1) .eq. HOEXTRAP) then
                 do i = imin-1, imax+1
                    do j = jmin-1, jmax+1
                       slz(i,j,kmin-1) = zero
                       slz(i,j,kmin) = -sixteen/fifteen*s(i,j,kmin-1) + half*s(i,j,kmin) + 
     &                      two3rd*s(i,j,kmin+1) - tenth*s(i,j,kmin+2)
                    end do
                 end do
              end if
              if (bc(3,2) .eq. EXT_DIR .or. bc(3,2) .eq. HOEXTRAP) then
                 do j = jmin-1, jmax+1
                    do i = imin-1, imax+1
                       slz(i,j,kmax) = -( -sixteen/fifteen*s(i,j,kmax+1) + half*s(i,j,kmax) + 
     &                      two3rd*s(i,j,kmax-1) - tenth*s(i,j,kmax-2) )
                       slz(i,j,kmax+1) = zero
                    end do
                 end do
              end if
           else
              do j = jmin-1,jmax+1
                 do i = imin-1,imax+1
                    do k = kmin-2,kmax+2
                       dmin           =  two*(s(i,j,k  )-s(i,j,k-1))
                       dpls           =  two*(s(i,j,k+1)-s(i,j,k  ))
                       slzscr(k,cen)  = half*(s(i,j,k+1)-s(i,j,k-1))
                       slzscr(k,lim)  = min(abs(dmin),abs(dpls))
                       slzscr(k,lim)  = cvmgp(slzscr(k,lim),zero,dpls*dmin)
                       slzscr(k,flag) = sign(one,slzscr(k,cen))
                       slzscr(k,fromm)= slzscr(k,flag)*
     &                      min(slzscr(k,lim),abs(slzscr(k,cen)))
                    end do
                    do k = kmin-1,kmax+1
                       ds = two * two3rd * slzscr(k,cen) -
     &                      sixth * (slzscr(k+1,fromm) + slzscr(k-1,fromm))
                       slz(i,j,k) = slzscr(k,flag)*min(abs(ds),slzscr(k,lim))
                    end do
c
                    if (bc(3,1) .eq. EXT_DIR .or. bc(3,1) .eq. HOEXTRAP) then
                       del  = -sixteen/fifteen*s(i,j,kmin-1) + half*s(i,j,kmin) +
     &                      two3rd*s(i,j,kmin+1) - tenth*s(i,j,kmin+2)
                       dmin = two*(s(i,j,kmin  )-s(i,j,kmin-1))
                       dpls = two*(s(i,j,kmin+1)-s(i,j,kmin  ))
                       slim = min(abs(dpls), abs(dmin))
                       slim = cvmgp(slim, zero, dpls*dmin)
                       sflg = sign(one,del)
                       slz(i,j,kmin-1) = zero
                       slz(i,j,kmin  ) = sflg*min(slim,abs(del))

c                      Recalculate the slope at jmin+1 using the revised slzscr(kmin,fromm)
                       slzscr(kmin,fromm) = slz(i,j,kmin)
                       ds = two * two3rd * slzscr(kmin+1,cen) -
     $                    sixth * (slzscr(kmin+2,fromm) + slzscr(kmin,fromm))
                       slz(i,j,kmin+1) = slzscr(kmin+1,flag)*min(abs(ds),slzscr(kmin+1,lim))
                    end if
                    if (bc(3,2) .eq. EXT_DIR .or. bc(3,2) .eq. HOEXTRAP) then
                       del  = sixteen/fifteen*s(i,j,kmax+1) - half*s(i,j,kmax)
     &                      - two3rd*s(i,j,kmax-1) + tenth*s(i,j,kmax-2)
                       dmin = two*(s(i,j,kmax  )-s(i,j,kmax-1))
                       dpls = two*(s(i,j,kmax+1)-s(i,j,kmax  ))
                       slim = min(abs(dpls), abs(dmin))
                       slim = cvmgp(slim, zero, dpls*dmin)
                       sflg = sign(one,del)
                       slz(i,j,kmax  ) = sflg*min(slim,abs(del))
                       slz(i,j,kmax+1) = zero

c                      Recalculate the slope at jmax-1 using the revised slzscr(kmax,fromm)
                       slzscr(kmax,fromm) = slz(i,j,kmax)
                       ds = two * two3rd * slzscr(kmax-1,cen) -
     $                    sixth * (slzscr(kmax-2,fromm) + slzscr(kmax,fromm))
                       slz(i,j,kmax-1) = slzscr(kmax-1,flag)*min(abs(ds),slzscr(kmax-1,lim))
                    end if
                 end do
              end do
           end if
        end if
c
c ... end, if slope_order .eq. 4
c
      end if

      end

      subroutine FORT_CONVSCALMINMAX(s,DIMS(s),sn,DIMS(sn),
     &                               smin,smax,DIMS(smin),lo,hi,bc)
c
c     correct an advected field for under/over shoots
c
      implicit none
      integer  i, j, k, imin, imax, jmin, jmax, kmin, kmax
      integer  DIMDEC(s)
      integer  DIMDEC(sn)
      integer  DIMDEC(smin)
      integer  lo(SDIM), hi(SDIM)
      integer  bc(SDIM,2), do_minmax
      REAL_T   s(DIMV(s))
      REAL_T   sn(DIMV(sn))
      integer  km, kk, kp
      REAL_T   smn, smx
      REAL_T   smin(DIM12(smin),0:2)
      REAL_T   smax(DIM12(smin),0:2)

      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
      kmin = lo(3)
      kmax = hi(3)
c
c     correct the 8 corners
c
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
         s(imin-1,jmin-1,kmin-1) = s(imin,jmin,kmin)
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
         s(imax+1,jmin-1,kmin-1) = s(imax,jmin,kmin)
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
         s(imin-1,jmax+1,kmin-1) = s(imin,jmax,kmin)
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
         s(imax+1,jmax+1,kmin-1) = s(imax,jmax,kmin)
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
         s(imin-1,jmin-1,kmax+1) = s(imin,jmin,kmax)
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
         s(imax+1,jmin-1,kmax+1) = s(imax,jmin,kmax)
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
         s(imin-1,jmax+1,kmax+1) = s(imin,jmax,kmax)
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
         s(imax+1,jmax+1,kmax+1) = s(imax,jmax,kmax)
      end if
c
c     correct the 12 edges
c
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
        do k = kmin,kmax
          s(imin-1,jmin-1,k) = s(imin,jmin,k)
        end do
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
        do k = kmin,kmax
          s(imax+1,jmin-1,k) = s(imax,jmin,k)
        end do
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
        do k = kmin,kmax
          s(imin-1,jmax+1,k) = s(imin,jmax,k)
        end do
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
        do k = kmin,kmax
          s(imax+1,jmax+1,k) = s(imax,jmax,k)
        end do
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do j = jmin,jmax
          s(imin-1,j,kmin-1) = s(imin,j,kmin)
        end do
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do j = jmin,jmax
          s(imax+1,j,kmin-1) = s(imax,j,kmin)
        end do
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do j = jmin,jmax
          s(imin-1,j,kmax+1) = s(imin,j,kmax)
        end do
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do j = jmin,jmax
          s(imax+1,j,kmax+1) = s(imax,j,kmax)
        end do
      end if
      if (bc(2,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do i = imin,imax
          s(i,jmin-1,kmin-1) = s(i,jmin,kmin)
        end do
      end if
      if (bc(2,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do i = imin,imax
          s(i,jmax+1,kmin-1) = s(i,jmax,kmin)
        end do
      end if
      if (bc(2,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do i = imin,imax
          s(i,jmin-1,kmax+1) = s(i,jmin,kmax)
        end do
      end if
      if (bc(2,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do i = imin,imax
          s(i,jmax+1,kmax+1) = s(i,jmax,kmax)
        end do
      end if
c
c     ::::: compute min/max a slab at a time
c     ::::: compute min and max of neighbors on kmin-1 slab
c
      km = 0
      kk = 1
      kp = 2

      k = kmin-1
      do j = jmin, jmax         
         do i = imin, imax
            smin(i,j,km) = min(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &           s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &           s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
            smax(i,j,km) = max(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &           s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &           s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
         end do         
      end do
c
c     ::::: compute min and max of neighbors on kmin slab
c
      k = kmin
      do j = jmin, jmax         
         do i = imin, imax
            smin(i,j,kk) = min(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &           s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &           s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
            smax(i,j,kk) = max(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &           s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &           s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
         end do         
      end do

      do k = kmin, kmax
c
c        ::::: compute min and max of neighbors on k+1 slab
c
         do j = jmin, jmax     
            do i = imin, imax   
               smin(i,j,kp) = min(s(i-1,j-1,k+1),s(i,j-1,k+1),s(i+1,j-1,k+1),
     &              s(i-1,j  ,k+1),s(i,j  ,k+1),s(i+1,j  ,k+1),
     &              s(i-1,j+1,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))
               smax(i,j,kp) = max(s(i-1,j-1,k+1),s(i,j-1,k+1),s(i+1,j-1,k+1),
     &              s(i-1,j  ,k+1),s(i,j  ,k+1),s(i+1,j  ,k+1),
     &              s(i-1,j+1,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))
c
c        ::::: compute min/max of cell
c
               smn = min(smin(i,j,km),smin(i,j,kk),smin(i,j,kp))
               smx = max(smax(i,j,km),smax(i,j,kk),smax(i,j,kp))
               sn(i,j,k) = max(sn(i,j,k),smn)
               sn(i,j,k) = min(sn(i,j,k),smx)
               
            end do
         end do
c
c        ::::: roll indices for next slab
c
         km = mod(km+1,3)
         kk = mod(kk+1,3)
         kp = mod(kp+1,3)
      end do

      end

      subroutine FORT_CONSSCALMINMAX(s,sn,rho,rhon,DIMS(s),
     &                               smin,smax,DIMS(smin),lo,hi,bc)
c
c     correct an conservatively-advected field for under/over shoots
c
      implicit none
      integer  i, j, k, imin, imax, jmin, jmax, kmin, kmax
      integer  DIMDEC(s)
      integer  DIMDEC(smin)
      integer  lo(SDIM), hi(SDIM)
      integer  bc(SDIM,2), do_minmax
      REAL_T      s(DIMV(s))
      REAL_T     sn(DIMV(s))
      REAL_T    rho(DIMV(s))
      REAL_T   rhon(DIMV(s))
      integer  km, kk, kp
      REAL_T   smn, smx
      REAL_T   smin(DIM12(smin),0:2)
      REAL_T   smax(DIM12(smin),0:2)

      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
      kmin = lo(3)
      kmax = hi(3)
c
c     correct the 8 corners
c
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
           s(imin-1,jmin-1,kmin-1) =   s(imin,jmin,kmin)
         rho(imin-1,jmin-1,kmin-1) = rho(imin,jmin,kmin)
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
           s(imax+1,jmin-1,kmin-1) =   s(imax,jmin,kmin)
         rho(imax+1,jmin-1,kmin-1) = rho(imax,jmin,kmin)
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
           s(imin-1,jmax+1,kmin-1) =   s(imin,jmax,kmin)
         rho(imin-1,jmax+1,kmin-1) = rho(imin,jmax,kmin)
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
           s(imax+1,jmax+1,kmin-1) =   s(imax,jmax,kmin)
         rho(imax+1,jmax+1,kmin-1) = rho(imax,jmax,kmin)
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
           s(imin-1,jmin-1,kmax+1) =   s(imin,jmin,kmax)
         rho(imin-1,jmin-1,kmax+1) = rho(imin,jmin,kmax)
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
           s(imax+1,jmin-1,kmax+1) =   s(imax,jmin,kmax)
         rho(imax+1,jmin-1,kmax+1) = rho(imax,jmin,kmax)
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
           s(imin-1,jmax+1,kmax+1) =   s(imin,jmax,kmax)
         rho(imin-1,jmax+1,kmax+1) = rho(imin,jmax,kmax)
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
           s(imax+1,jmax+1,kmax+1) =   s(imax,jmax,kmax)
         rho(imax+1,jmax+1,kmax+1) = rho(imax,jmax,kmax)
      end if
c
c     correct the 12 edges
c
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
        do k = kmin,kmax
            s(imin-1,jmin-1,k) =   s(imin,jmin,k)
          rho(imin-1,jmin-1,k) = rho(imin,jmin,k)
        end do
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
        do k = kmin,kmax
            s(imax+1,jmin-1,k) =   s(imax,jmin,k)
          rho(imax+1,jmin-1,k) = rho(imax,jmin,k)
        end do
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
        do k = kmin,kmax
            s(imin-1,jmax+1,k) =   s(imin,jmax,k)
          rho(imin-1,jmax+1,k) = rho(imin,jmax,k)
        end do
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
        do k = kmin,kmax
            s(imax+1,jmax+1,k) =   s(imax,jmax,k)
          rho(imax+1,jmax+1,k) = rho(imax,jmax,k)
        end do
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do j = jmin,jmax
            s(imin-1,j,kmin-1) =   s(imin,j,kmin)
          rho(imin-1,j,kmin-1) = rho(imin,j,kmin)
        end do
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do j = jmin,jmax
            s(imax+1,j,kmin-1) =   s(imax,j,kmin)
          rho(imax+1,j,kmin-1) = rho(imax,j,kmin)
        end do
      end if
      if (bc(1,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do j = jmin,jmax
            s(imin-1,j,kmax+1) =   s(imin,j,kmax)
          rho(imin-1,j,kmax+1) = rho(imin,j,kmax)
        end do
      end if
      if (bc(1,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do j = jmin,jmax
            s(imax+1,j,kmax+1) =   s(imax,j,kmax)
          rho(imax+1,j,kmax+1) = rho(imax,j,kmax)
        end do
      end if
      if (bc(2,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do i = imin,imax
            s(i,jmin-1,kmin-1) =   s(i,jmin,kmin)
          rho(i,jmin-1,kmin-1) = rho(i,jmin,kmin)
        end do
      end if
      if (bc(2,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do i = imin,imax
            s(i,jmax+1,kmin-1) =   s(i,jmax,kmin)
          rho(i,jmax+1,kmin-1) = rho(i,jmax,kmin)
        end do
      end if
      if (bc(2,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do i = imin,imax
            s(i,jmin-1,kmax+1) =   s(i,jmin,kmax)
          rho(i,jmin-1,kmax+1) = rho(i,jmin,kmax)
        end do
      end if
      if (bc(2,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do i = imin,imax
            s(i,jmax+1,kmax+1) =   s(i,jmax,kmax)
          rho(i,jmax+1,kmax+1) = rho(i,jmax,kmax)
        end do
      end if

      do k = kmin-1, kmax+1
      do j = jmin-1, jmax+1
         do i = imin-1, imax+1
            s(i,j,k) = s(i,j,k) / rho(i,j,k)
         end do
      end do
      end do
c
c     ::::: compute min/max a slab at a time
c     ::::: compute min and max of neighbors on kmin-1 slab
c
      km = 0
      kk = 1
      kp = 2

      k = kmin-1
      do j = jmin, jmax         
         do i = imin, imax
            smin(i,j,km) = min(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &           s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &           s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
            smax(i,j,km) = max(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &           s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &           s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k)) 
         end do         
      end do
c
c     ::::: compute min and max of neighbors on kmin slab
c
      k = kmin
      do j = jmin, jmax         
         do i = imin, imax
            smin(i,j,kk) = min(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &           s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &           s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
            smax(i,j,kk) = max(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &           s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &           s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
         end do         
      end do

      do k = kmin, kmax
c
c        ::::: compute min and max of neighbors on k+1 slab
c
         do j = jmin, jmax     
            do i = imin, imax   
               smin(i,j,kp) = min(s(i-1,j-1,k+1),s(i,j-1,k+1),s(i+1,j-1,k+1),
     &              s(i-1,j  ,k+1),s(i,j  ,k+1),s(i+1,j  ,k+1),
     &              s(i-1,j+1,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))
               smax(i,j,kp) = max(s(i-1,j-1,k+1),s(i,j-1,k+1),s(i+1,j-1,k+1),
     &              s(i-1,j  ,k+1),s(i,j  ,k+1),s(i+1,j  ,k+1),
     &              s(i-1,j+1,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))
c
c        ::::: compute min/max of cell
c
               smn = min(smin(i,j,km),smin(i,j,kk),smin(i,j,kp))
               smx = max(smax(i,j,km),smax(i,j,kk),smax(i,j,kp))
               sn(i,j,k) = max(sn(i,j,k)/rhon(i,j,k),smn) * rhon(i,j,k)
               sn(i,j,k) = min(sn(i,j,k)/rhon(i,j,k),smx) * rhon(i,j,k)
               
            end do
         end do
c
c        ::::: roll indices for next slab
c
         km = mod(km+1,3)
         kk = mod(kk+1,3)
         kp = mod(kp+1,3)
      end do

      do k = kmin-1, kmax+1
      do j = jmin-1, jmax+1
         do i = imin-1, imax+1
            s(i,j,k) = s(i,j,k) * rho(i,j,k)
         end do
      end do
      end do

      end

      subroutine FORT_SUM_TF_GP(
     &     tforces,DIMS(tf),
     &     gp,DIMS(gp),
     &     rho,DIMS(rho),
     &     lo,hi )
c
c     sum pressure forcing into tforces
c
      implicit none
      integer i, j, k, n
      integer DIMDEC(tf)
      integer DIMDEC(visc)
      integer DIMDEC(gp)
      integer DIMDEC(rho)
      integer lo(SDIM), hi(SDIM)
      REAL_T tforces(DIMV(tf),SDIM)
      REAL_T gp(DIMV(gp),SDIM)
      REAL_T rho(DIMV(rho))

      do n = 1, SDIM
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,k,n) = (
     &            tforces(i,j,k,n) - gp(i,j,k,n)
     &                 )/rho(i,j,k)
               end do
            end do
         end do
      end do

      end

      subroutine FORT_SUM_TF_GP_VISC(
     &     tforces,DIMS(tf),
     &     visc,DIMS(visc),
     &     gp,DIMS(gp),
     &     rho,DIMS(rho),
     &     lo,hi )
c
c     sum pressure forcing and viscous forcing into
c     tforces
c
      implicit none
      integer i, j, k, n
      integer DIMDEC(tf)
      integer DIMDEC(visc)
      integer DIMDEC(gp)
      integer DIMDEC(rho)
      integer lo(SDIM), hi(SDIM)
      REAL_T tforces(DIMV(tf),SDIM)
      REAL_T visc(DIMV(visc),SDIM)
      REAL_T gp(DIMV(gp),SDIM)
      REAL_T rho(DIMV(rho))

      do n = 1, SDIM
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,k,n) = (
     &                 tforces(i,j,k,n)
     &                 +  visc(i,j,k,n)
     &                 -    gp(i,j,k,n) )/rho(i,j,k)
               end do
            end do
         end do
      end do

      end

      subroutine FORT_SUM_TF_DIVU(
     &     S,DIMS(S),
     &     tforces,DIMS(tf),
     &     divu,DIMS(divu),
     &     rho,DIMS(rho),
     &     lo,hi,nvar,iconserv )
c
c     sum tforces, viscous forcing and divU*S into tforces
c     depending on the value of iconserv
c
      implicit none
      integer nvar, iconserv
      integer lo(SDIM), hi(SDIM)
      integer i, j, k, n

      integer DIMDEC(S)
      integer DIMDEC(tf)
      integer DIMDEC(divu)
      integer DIMDEC(rho)

      REAL_T S(DIMV(S),nvar)
      REAL_T tforces(DIMV(tf),nvar)
      REAL_T divu(DIMV(divu))
      REAL_T rho(DIMV(rho))

      if ( iconserv .eq. 1 ) then
         do n = 1, nvar
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tforces(i,j,k,n) = 
     &               tforces(i,j,k,n) - S(i,j,k,n)*divu(i,j,k)
                  end do
               end do
            end do
         end do
      else
         do n = 1, nvar
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tforces(i,j,k,n) = 
     &               tforces(i,j,k,n)/rho(i,j,k)
                  end do
               end do
            end do
         end do
      end if

      end

      subroutine FORT_SUM_TF_DIVU_VISC(
     &     S,DIMS(S),
     &     tforces,DIMS(tf),
     &     divu,DIMS(divu),
     &     visc,DIMS(visc),
     &     rho,DIMS(rho),
     &     lo,hi,nvar,iconserv )
c
c     sum tforces, viscous forcing and divU*S into tforces
c     depending on the value of iconserv
c
      implicit none
      integer nvar, iconserv
      integer lo(SDIM), hi(SDIM)
      integer i, j, k, n

      integer DIMDEC(S)
      integer DIMDEC(tf)
      integer DIMDEC(divu)
      integer DIMDEC(visc)
      integer DIMDEC(rho)

      REAL_T S(DIMV(S),nvar)
      REAL_T tforces(DIMV(tf),nvar)
      REAL_T divu(DIMV(divu))
      REAL_T visc(DIMV(visc),nvar)
      REAL_T rho(DIMV(rho))

      if ( iconserv .eq. 1 ) then
         do n = 1, nvar
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tforces(i,j,k,n) = 
     &                    tforces(i,j,k,n)
     &                    +  visc(i,j,k,n)
     &                    -     S(i,j,k,n)*divu(i,j,k)
                  end do
               end do
            end do
         end do
      else
         do n = 1, nvar
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tforces(i,j,k,n) = (
     &                    tforces(i,j,k,n)
     &                    +  visc(i,j,k,n) )/rho(i,j,k)
                  end do
               end do
            end do
         end do
      end if

      end

      subroutine FORT_UPDATE_TF(
     &     s,       DIMS(s),
     &     sn,      DIMS(sn),
     &     tforces, DIMS(tf),
     &     lo,hi,dt,nvar)
c
c     update a field with a forcing term
c
      implicit none
      integer i, j, k, n, nvar
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(tf)
      integer lo(SDIM), hi(SDIM)
      REAL_T dt
      REAL_T s(DIMV(s),nvar)
      REAL_T sn(DIMV(sn),nvar)
      REAL_T tforces(DIMV(tf),nvar)

      do n = 1,nvar
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  sn(i,j,k,n) = s(i,j,k,n)
     &                 + dt*tforces(i,j,k,n)
               end do
            end do
         end do
      end do

      end

      subroutine FORT_CORRECT_TF(
     &     ss,  sp,  DIMS(ss),
     &     tfs, tfn, DIMS(tfs),
     &     lo,hi,dt,nvar)
c
c     correct 1st order rk to second-order
c
      implicit none
      integer i, j, k, n, nvar
      integer lo(SDIM), hi(SDIM)
      REAL_T dt,hdt

      integer DIMDEC(ss)
      integer DIMDEC(tfs)
      REAL_T  ss(DIMV(ss),nvar)
      REAL_T  sp(DIMV(ss),nvar)
      REAL_T tfs(DIMV(tfs),nvar)
      REAL_T tfn(DIMV(tfs),nvar)

      hdt = half*dt
      do n = 1,nvar
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  sp(i,j,k,n) = ss(i,j,k,n)
     &                 + hdt*(tfs(i,j,k,n)-tfn(i,j,k,n))
               end do
            end do
         end do
      end do

      end

      subroutine FORT_UPDATE_AOFS_TF(
     &     s,       DIMS(s),
     &     sn,      DIMS(sn),
     &     aofs,    DIMS(aofs),
     &     tforces, DIMS(tf),
     &     lo,hi,dt,nvar)
c
c     update a field with an advective tendency
c     and a forcing term
c
      implicit none
      integer i, j, k, n, nvar
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(aofs)
      integer DIMDEC(tf)
      integer lo(SDIM), hi(SDIM)
      REAL_T dt
      REAL_T s(DIMV(s),nvar)
      REAL_T sn(DIMV(sn),nvar)
      REAL_T aofs(DIMV(aofs),nvar)
      REAL_T tforces(DIMV(tf),nvar)

      do n = 1,nvar
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  sn(i,j,k,n) = s(i,j,k,n)
     &                 - dt*aofs(i,j,k,n)
     &                 + dt*tforces(i,j,k,n)
               end do
            end do
         end do
      end do

      end

      subroutine FORT_UPDATE_AOFS_TF_GP(
     &     u,       DIMS(u),
     &     un,      DIMS(un),
     &     aofs,    DIMS(aofs),
     &     tforces, DIMS(tf),
     &     gp,      DIMS(gp),
     &     rho,     DIMS(rho),
     &     lo, hi, dt)

c
c     update the velocities
c
      implicit none
      integer i, j, k, n
      integer DIMDEC(u)
      integer DIMDEC(un)
      integer DIMDEC(aofs)
      integer DIMDEC(rho)
      integer DIMDEC(gp)
      integer DIMDEC(tf)
      integer lo(SDIM), hi(SDIM)
      REAL_T u(DIMV(u),SDIM)
      REAL_T un(DIMV(un),SDIM)
      REAL_T aofs(DIMV(aofs),SDIM)
      REAL_T rho(DIMV(rho))
      REAL_T gp(DIMV(gp),SDIM)
      REAL_T tforces(DIMV(tf),SDIM)
      REAL_T dt

      do n = 1, SDIM
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  un(i,j,k,n) = u(i,j,k,n)
     &                 - dt*   aofs(i,j,k,n)
     &                 + dt*tforces(i,j,k,n)/rho(i,j,k)
     &                 - dt*     gp(i,j,k,n)/rho(i,j,k)
               end do
            end do
         end do
      end do

      end
