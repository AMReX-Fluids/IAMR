#undef BL_LANG_CC
#define BL_LANG_FORT

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "GODUNOV_F.H"
#include "DIMS.H"

#define SDIM 3
#define XVEL 1
#define YVEL 2
#define ZVEL 3
#define RHO  4

c=========================================================



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
      integer i, j, k
      REAL_T  u, v, w
      REAL_T  small
      REAL_T  dt_start
      REAL_T  tforce1,tforce2,tforce3
      integer lo(SDIM), hi(SDIM)
      REAL_T  dt,dx(SDIM),cfl,u_max(SDIM)
c
      integer DIMDEC(vel)
      integer DIMDEC(rho)
      integer DIMDEC(tf)
c
      REAL_T  vel(DIMV(vel),SDIM)
      REAL_T  rho(DIMV(rho))
      REAL_T  tforces(DIMV(tf),SDIM)
c
      small   = 1.0e-8
      u       = zero
      v       = zero
      w       = zero
      tforce1 = zero
      tforce2 = zero
      tforce3 = zero
c
      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
         u = max(u,abs(vel(i,j,k,1)))	
	 v = max(v,abs(vel(i,j,k,2)))
	 w = max(w,abs(vel(i,j,k,3)))
         tforce1 = max(tforce1,abs(tforces(i,j,k,1)/rho(i,j,k)))
         tforce2 = max(tforce2,abs(tforces(i,j,k,2)/rho(i,j,k)))
         tforce3 = max(tforce3,abs(tforces(i,j,k,3)/rho(i,j,k)))
      enddo
      enddo
      enddo

      u_max(1) = u
      u_max(2) = v
      u_max(3) = w

      dt_start = 1.0e+20
      dt = dt_start

      if (u .gt. small) dt = min(dt,dx(1)/u)
      if (v .gt. small) dt = min(dt,dx(2)/v)
      if (w .gt. small) dt = min(dt,dx(3)/w)

      if (tforce1 .gt. small) then
        dt  = min(dt,sqrt(two*dx(1)/tforce1))
      endif

      if (tforce2 .gt. small) then
        dt  = min(dt,sqrt(two*dx(2)/tforce2))
      endif

      if (tforce3 .gt. small) then
        dt  = min(dt,sqrt(two*dx(3)/tforce3))
      endif

      if (dt .eq. dt_start) dt = min(dx(1),dx(2),dx(3))

      dt = dt*cfl
c
c     end of routine
c
      return
      end


c=========================================================

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
c
      hx   = dx(1)
      hy   = dx(2)
      hz   = dx(3)
      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
      kmin = lo(3)
      kmax = hi(3)
      umax = -1.e30
      vmax = -1.e30
      wmax = -1.e30
      umin =  1.e30
      vmin =  1.e30
      wmin =  1.e30
      rhomax = -1.e30
      rhomin =  1.e30
c
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
      enddo
      enddo
      enddo

      u_max(1) = max(abs(umax), abs(umin))
      u_max(2) = max(abs(vmax), abs(vmin))
      u_max(3) = max(abs(wmax), abs(wmin))
      cflmax   = dt*max(u_max(1)/hx,u_max(2)/hy,u_max(3)/hz)

      if(verbose.eq.1)then
        write(6,1000) umax,umin,u_max(1)
        write(6,1001) vmax,vmin,u_max(2)
        write(6,1002) wmax,wmin,u_max(3)
        write(6,1003) rhomax,rhomin
        call flush(6)
      endif

 1000 format(' U  MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1001 format(' V  MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1002 format(' W  MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1003 format('RHO MAX/MIN      ',e21.14,2x,e21.14)

c
c     end of routine
c
      return
      end



c=========================================================

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
      integer lo(SDIM),hi(SDIM)
      REAL_T  dt, dx(SDIM), u_max(SDIM), cflmax
      integer imin, imax, jmin, jmax, kmin, kmax
      integer i, j, k
      REAL_T  hx, hy, hz
      REAL_T  umax, vmax, wmax, rhomax
      REAL_T  umin, vmin, wmin, rhomin
c
      integer DIMDEC(umac)
      integer DIMDEC(vmac)
      integer DIMDEC(wmac)
      integer DIMDEC(rho)
c
      REAL_T  umac(DIMV(umac))
      REAL_T  vmac(DIMV(vmac))
      REAL_T  wmac(DIMV(wmac))
      REAL_T  rho(DIMV(rho))
c
      hx   = dx(1)
      hy   = dx(2)
      hz   = dx(3)
      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
      kmin = lo(3)
      kmax = hi(3)
      umax = -1.e30
      vmax = -1.e30
      wmax = -1.e30
      umin =  1.e30
      vmin =  1.e30
      wmin =  1.e30
      rhomax = -1.e30
      rhomin =  1.e30
c
      do k = kmin, kmax
      do j = jmin, jmax
      do i = imin, imax+1
         umax = max(umax,umac(i,j,k))
         umin = min(umin,umac(i,j,k))
      enddo
      enddo
      enddo
c
      do k = kmin, kmax
      do j = jmin, jmax+1
      do i = imin, imax
         vmax = max(vmax,vmac(i,j,k))
         vmin = min(vmin,vmac(i,j,k))
      enddo
      enddo
      enddo
c
      do k = kmin, kmax+1
      do j = jmin, jmax
      do i = imin, imax
         wmax = max(wmax,wmac(i,j,k))
         wmin = min(wmin,wmac(i,j,k))
      enddo
      enddo
      enddo
c
      do k = kmin, kmax
      do j = jmin, jmax
      do i = imin, imax
         rhomax = max(rhomax,rho(i,j,k))
         rhomin = min(rhomin,rho(i,j,k))
      enddo
      enddo
      enddo

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

      call flush(6)
c
c     end of routine
c
      return
      end


c=========================================================

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
      REAL_T eps, val, tst
      logical ltm
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
c
      eps  = 1.e-6
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

      call FORT_SLOPES( XVEL,
     &     u,DIMS(s),
     &     sx,sy,sz,DIMS(work),
     &     lo,hi,slxscr,slyscr,slzscr,ubc)

      do k = kmin-1,kmax+1
         do j = jmin-1,jmax+1
         do i = imin,  imax+1
         ulo(i,j,k) = u(i-1,j,k) + (half  - dthx*u(i-1,j,k))*sx(i-1,j,k)
         uhi(i,j,k) = u(i,  j,k) + (-half - dthx*u(i,  j,k))*sx(i,  j,k)
         enddo
         enddo
      enddo

      if(use_minion .eq. 1 )then
        do k = kmin-1,kmax+1
           do j = jmin-1,jmax+1
           do i = imin,  imax+1
           ulo(i,j,k) = ulo(i,j,k) + dth*tforces(i-1,j,k,1)
           uhi(i,j,k) = uhi(i,j,k) + dth*tforces(i,  j,k,1)
           enddo
           enddo
        enddo
      endif

      call FORT_SLOPES(YVEL,
     &     v,DIMS(s),
     &     sx,sy,sz,DIMS(work),
     &     lo,hi,slxscr,slyscr,slzscr,vbc)

      do k = kmin-1,kmax+1
         do j = jmin,  jmax+1
         do i = imin-1,imax+1
         vlo(i,j,k) = v(i,j-1,k) + (half  - dthy*v(i,j-1,k))*sy(i,j-1,k)
         vhi(i,j,k) = v(i,j,  k) + (-half - dthy*v(i,j,  k))*sy(i,j,  k)
         enddo
         enddo
      enddo

      if(use_minion .eq. 1 )then
        do k = kmin-1,kmax+1
           do j = jmin,    jmax+1
           do i = imin-1,  imax+1
           vlo(i,j,k) = vlo(i,j,k) + dth*tforces(i,j-1,k,2)
           vhi(i,j,k) = vhi(i,j,k) + dth*tforces(i,j,  k,2)
           enddo
           enddo
        enddo
      endif

      call FORT_SLOPES(ZVEL,
     &     w,DIMS(s),
     &     sx,sy,sz,DIMS(work),
     &     lo,hi,slxscr,slyscr,slzscr,wbc)

      do k = kmin,kmax+1
         do j = jmin-1,jmax+1
         do i = imin-1,imax+1
         wlo(i,j,k) = w(i,j,k-1) + (half  - dthz*w(i,j,k-1))*sz(i,j,k-1)
         whi(i,j,k) = w(i,j,k  ) + (-half - dthz*w(i,j,k  ))*sz(i,j,k  )
         enddo
         enddo
      enddo

      if(use_minion .eq. 1 )then
        do k = kmin,kmax+1
           do j = jmin-1,jmax+1
           do i = imin-1,  imax+1
           wlo(i,j,k) = wlo(i,j,k) + dth*tforces(i,j,k-1,3)
           whi(i,j,k) = whi(i,j,k) + dth*tforces(i,j,k,  3)
           enddo
           enddo
        enddo
      endif

      call trans_xbc(
     &     u,DIMS(s),
     &     ulo,uhi,ulo,DIMS(work),
     &     lo,hi,XVEL,ubc)

      call trans_ybc(
     &     v,DIMS(s),
     &     vlo,vhi,vlo,DIMS(work),
     &     lo,hi,YVEL,vbc)

      call trans_zbc(
     &     w,DIMS(s),
     &     wlo,whi,wlo,DIMS(work),
     &     lo,hi,ZVEL,wbc)

      do k = kmin-1,kmax+1
         do j = jmin-1,jmax+1
         do i = imin,  imax+1
            tst = ulo(i,j,k)+uhi(i,j,k)
            val = cvmgp(ulo(i,j,k),uhi(i,j,k),tst)
            ltm =
     &         ( (ulo(i,j,k) .le. zero) .and.
     &           (uhi(i,j,k) .ge. zero) ) .or.
     &           (abs(tst)   .lt. eps )
            ulo(i,j,k) = cvmgt(zero,val,ltm)
         enddo
         enddo
      enddo

      do k = kmin-1,kmax+1
         do j = jmin,  jmax+1
         do i = imin-1,imax+1
         tst = vlo(i,j,k)+vhi(i,j,k)
         val = cvmgp(vlo(i,j,k),vhi(i,j,k),tst)
         ltm =
     &      ( (vlo(i,j,k) .le. zero) .and.
     &        (vhi(i,j,k) .ge. zero) ) .or.
     &        (abs(tst)   .lt. eps )
         vlo(i,j,k) = cvmgt(zero,val,ltm)
         enddo
         enddo
      enddo
      
      do k = kmin,kmax+1
         do j = jmin-1,jmax+1
         do i = imin-1,imax+1
         tst = wlo(i,j,k)+whi(i,j,k)
         val = cvmgp(wlo(i,j,k),whi(i,j,k),tst)
         ltm =
     &      ( (wlo(i,j,k) .le. zero) .and.
     &        (whi(i,j,k) .ge. zero) ) .or.
     &        (abs(tst)   .lt. eps )
         wlo(i,j,k) = cvmgt(zero,val,ltm)
         enddo
         enddo
      enddo
c     
c     end of routine
c
      return
      end



      
c=========================================================

      subroutine FORT_ESTATE(s, tforces, DIMS(s),

     &     u, xlo, xhi, sx, uad, slxscr, stxlo, stxhi,
     &     uedge, xstate, DIMS(uedge),

     &     v, ylo, yhi, sy, vad, slyscr, stylo, styhi,
     &     vedge, ystate, DIMS(vedge),

     &     w, zlo, zhi, sz, wad, slzscr, stzlo, stzhi,
     &     wedge, zstate, DIMS(wedge),

     &     DIMS(work),
     &     bc,lo,hi,dt,dx,n,velpred, use_minion)
c
c     This subroutine computes edges states, right now it uses
c     a lot of memory, but there becomes a trade off between
c     simplicity-efficiency in the new way of computing states
c     and complexity in the old way.  By eliminating loops over
c     state components though, the new way uses much less memory.
c
      integer i,j,k,n,velpred
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer imin,jmin,kmin,imax,jmax,kmax
      REAL_T place_to_break
      REAL_T hx, hy, hz, dt, dth, dthx, dthy, dthz
      REAL_T tr,stx,sty,stz,fu,fv,fw,dx(SDIM)
      REAL_T eps
      logical ltx,lty,ltz
      integer DIMDEC(s)
      integer DIMDEC(work)
      integer DIMDEC(uedge)
      integer DIMDEC(vedge)
      integer DIMDEC(wedge)
c
      REAL_T s(DIMV(s))
      REAL_T u(DIMV(s))
      REAL_T v(DIMV(s))
      REAL_T w(DIMV(s))
      REAL_T stxlo(DIM1(s)),stxhi(DIM1(s)),slxscr(DIM1(s),4)
      REAL_T stylo(DIM2(s)),styhi(DIM2(s)),slyscr(DIM2(s),4)
      REAL_T stzlo(DIM3(s)),stzhi(DIM3(s)),slzscr(DIM3(s),4)
c
      REAL_T uedge(DIMV(uedge)), xstate(DIMV(uedge))
      REAL_T vedge(DIMV(vedge)), ystate(DIMV(vedge))
      REAL_T wedge(DIMV(wedge)), zstate(DIMV(wedge))
c
      REAL_T xlo(DIMV(work)), xhi(DIMV(work))
      REAL_T ylo(DIMV(work)), yhi(DIMV(work))
      REAL_T zlo(DIMV(work)), zhi(DIMV(work))
      REAL_T  sx(DIMV(work)), uad(DIMV(work))
      REAL_T  sy(DIMV(work)), vad(DIMV(work))
      REAL_T  sz(DIMV(work)), wad(DIMV(work))
      REAL_T tforces(DIMV(work))

      integer use_minion
c
      eps  = 1.0e-6
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

c     --------------------------------------------------------------
c     compute the slopes
c     --------------------------------------------------------------
      call FORT_SLOPES(10,
     &     s,DIMS(s),
     &     sx,sy,sz,DIMS(work),
     &     lo,hi,slxscr,slyscr,slzscr,bc)

c     --------------------------------------------------------------
c     trace the state to the cell edges
c     --------------------------------------------------------------
      do k = kmin-1,kmax+1
         do j = jmin-1,jmax+1
         do i = imin,  imax+1
         xlo(i,j,k) = s(i-1,j,k) + (half  - dthx*u(i-1,j,k))*sx(i-1,j,k)
         xhi(i,j,k) = s(i,  j,k) + (-half - dthx*u(i,  j,k))*sx(i,  j,k)
         enddo
         enddo
      enddo

      if(use_minion.eq.1)then
        do k = kmin-1,kmax+1
           do j = jmin-1,jmax+1
           do i = imin,  imax+1
           xlo(i,j,k) = xlo(i,j,k) + dth*tforces(i-1,j,k)
           xhi(i,j,k) = xhi(i,j,k) + dth*tforces(i,  j,k)
           enddo
           enddo
        enddo
      endif

      call trans_xbc(
     &     s,DIMS(s),
     &     xlo,xhi,uad,DIMS(work),
     &     lo,hi,n,bc)

      do k = kmin-1,kmax+1
         do j = jmin-1,jmax+1
         do i = imin,  imax+1
         fu  = cvmgt(zero,one,abs(uad(i,j,k)).lt.eps)
         stx = cvmgp(xlo(i,j,k),xhi(i,j,k),uad(i,j,k))
         xlo(i,j,k) = fu*stx + (one - fu)*half*(xhi(i,j,k)+xlo(i,j,k))
         enddo
         enddo
      enddo
c     --------------------------------------------------------------
      do k = kmin-1,kmax+1
         do j = jmin,  jmax+1
         do i = imin-1,imax+1
         ylo(i,j,k) = s(i,j-1,k) + (half  - dthy*v(i,j-1,k))*sy(i,j-1,k)
         yhi(i,j,k) = s(i,j, k)  + (-half - dthy*v(i,j,  k))*sy(i,j, k)
         enddo
         enddo
      enddo

      if(use_minion.eq.1)then
        do k = kmin-1,kmax+1
           do j = jmin, jmax+1
           do i = imin-1,  imax+1
           ylo(i,j,k) = ylo(i,j,k) + dth*tforces(i,j-1,k)
           yhi(i,j,k) = yhi(i,j,k) + dth*tforces(i,j,  k)
           enddo
           enddo
        enddo
      endif

      call trans_ybc(
     &     s,DIMS(s),
     &     ylo,yhi,vad,DIMS(work),
     &     lo,hi,n,bc)

      do k = kmin-1,kmax+1
         do j = jmin,  jmax+1
         do i = imin-1,imax+1
         fv  = cvmgt(zero,one,abs(vad(i,j,k)).lt.eps)
         sty = cvmgp(ylo(i,j,k),yhi(i,j,k),vad(i,j,k))
         ylo(i,j,k) = fv*sty + (one - fv)*half*(yhi(i,j,k)+ylo(i,j,k))
         enddo
         enddo
      enddo
c     ------------------------------------------------------------
      do k = kmin,kmax+1
         do j = jmin-1,jmax+1
         do i = imin-1,imax+1
         zlo(i,j,k) = s(i,j,k-1) + (half  - dthz*w(i,j,k-1))*sz(i,j,k-1)
         zhi(i,j,k) = s(i,j,k  ) + (-half - dthz*w(i,j,k  ))*sz(i,j,k  )
         enddo
         enddo
      enddo

      if(use_minion.eq.1)then
        do k = kmin,kmax+1
           do j = jmin-1,jmax+1
           do i = imin-1,  imax+1
           zlo(i,j,k) = zlo(i,j,k) + dth*tforces(i,j,k-1)
           zhi(i,j,k) = zhi(i,j,k) + dth*tforces(i,j,k)
           enddo
           enddo
        enddo
      endif

      call trans_zbc(
     &     s,DIMS(s),
     &     zlo,zhi,wad,DIMS(work),
     &     lo,hi,n,bc)

      do k = kmin,kmax+1
         do j = jmin-1,jmax+1
         do i = imin-1,imax+1
         fw  = cvmgt(zero,one,abs(wad(i,j,k)).lt.eps)
         stz = cvmgp(zlo(i,j,k),zhi(i,j,k),wad(i,j,k))
         zlo(i,j,k) = fw*stz + (one-fw)*half*(zhi(i,j,k)+zlo(i,j,k))
         enddo
         enddo
      enddo
            

c     --------------------------------------------------------------
c     compute the xedge states
c     --------------------------------------------------------------
      if ((velpred.ne.1) .or. (n.eq.XVEL)) then
      do 110 k = kmin,kmax
      do 100 j = jmin,jmax

      do i = imin-1,imax+1
         tr =
     &        half*
     &        (vad(i,j+1,k) + vad(i,j,k)) *
     &        (ylo(i,j+1,k) - ylo(i,j,k)) / hy +
     &        half*
     &        (wad(i,j,k+1) + wad(i,j,k)) *
     &        (zlo(i,j,k+1) - zlo(i,j,k)) / hz
         
         stxlo(i+1)= s(i,j,k) + (half-dthx*u(i,j,k))*sx(i,j,k)
     &        - dth*tr
     &        + dth*tforces(i,j,k)
         stxhi(i  )= s(i,j,k) - (half+dthx*u(i,j,k))*sx(i,j,k)
     &        - dth*tr
     &        + dth*tforces(i,j,k)
      enddo

      if (bc(1,1).eq.EXT_DIR) then
         stxhi(imin) = s(imin-1,j,k)
         stxlo(imin) = s(imin-1,j,k)
      elseif (bc(1,1).eq.EXTRAP.or.bc(1,1).eq.HOEXTRAP
     &        .or.bc(1,1).eq.REFLECT_EVEN) then
         stxlo(imin) = stxhi(imin)
      elseif (bc(1,1).eq.REFLECT_ODD) then
         stxhi(imin) = zero
         stxlo(imin) = zero
      endif
      if (bc(1,2).eq.EXT_DIR) then
        stxlo(imax+1) = s(imax+1,j,k)
        stxhi(imax+1) = s(imax+1,j,k)
      elseif (bc(1,2).eq.EXTRAP.or.bc(1,2).eq.HOEXTRAP
     &    .or.bc(1,2).eq.REFLECT_EVEN) then
        stxhi(imax+1) = stxlo(imax+1)
      elseif (bc(1,2).eq.REFLECT_ODD) then
        stxlo(imax+1) = zero
        stxhi(imax+1) = zero
      endif

      if ( velpred .eq. 1 ) then
         do i = imin, imax+1
            ltx = stxlo(i) .le. zero  .and.  stxhi(i) .ge. zero
            ltx = ltx .or. (abs(stxlo(i)+stxhi(i)) .lt. eps)
            stx = cvmgp(stxlo(i),stxhi(i),stxlo(i)+stxhi(i))
            xstate(i,j,k) = cvmgt(zero,stx,ltx)
         enddo
      else
         do i = imin, imax+1
            xstate(i,j,k) = cvmgp(stxlo(i),stxhi(i),uedge(i,j,k))
            xstate(i,j,k) = cvmgt(half*(stxlo(i)+stxhi(i)),xstate(i,j,k)
     &           ,abs(uedge(i,j,k)).lt.eps)
         enddo
      endif
      place_to_break = 1
  100 continue
  110 continue
      endif

c     --------------------------------------------------------------
c     compute the yedge states
c     --------------------------------------------------------------
      if ((velpred.ne.1) .or. (n.eq.YVEL)) then
      do 210 k = kmin,kmax
      do 200 i = imin,imax
      
      do j = jmin-1,jmax+1
        tr =
     &       half*
     &       (uad(i+1,j,k) + uad(i,j,k)) *
     &       (xlo(i+1,j,k) - xlo(i,j,k)) / hx +
     &       half*
     &       (wad(i,j,k+1) + wad(i,j,k)) *
     &       (zlo(i,j,k+1) - zlo(i,j,k)) / hz

        stylo(j+1)= s(i,j,k) + (half-dthy*v(i,j,k))*sy(i,j,k)
     &       - dth*tr
     &       + dth*tforces(i,j,k)
        styhi(j)  = s(i,j,k) - (half+dthy*v(i,j,k))*sy(i,j,k)
     &       - dth*tr
     &       + dth*tforces(i,j,k)
      enddo

      if (bc(2,1).eq.EXT_DIR) then
         styhi(jmin) = s(i,jmin-1,k)
         stylo(jmin) = s(i,jmin-1,k)
      elseif (bc(2,1).eq.EXTRAP.or.bc(2,1).eq.HOEXTRAP
     &        .or.bc(2,1).eq.REFLECT_EVEN) then
         stylo(jmin) = styhi(jmin)
      elseif (bc(2,1).eq.REFLECT_ODD) then
         styhi(jmin) = zero
         stylo(jmin) = zero
      endif
      
      if (bc(2,2).eq.EXT_DIR) then
         stylo(jmax+1) = s(i,jmax+1,k)
         styhi(jmax+1) = s(i,jmax+1,k)
      elseif (bc(2,2).eq.EXTRAP.or.bc(2,2).eq.HOEXTRAP
     &        .or.bc(2,2).eq.REFLECT_EVEN) then
        styhi(jmax+1) = stylo(jmax+1)
      elseif (bc(2,2).eq.REFLECT_ODD) then
        stylo(jmax+1) = zero
        styhi(jmax+1) = zero
      endif

      if ( velpred .eq. 1 ) then
         do j = jmin, jmax+1
         lty = stylo(j) .le. zero  .and.  styhi(j) .ge. zero
	 lty = lty .or. (abs(stylo(j)+styhi(j)) .lt. eps)
         sty = cvmgp(stylo(j),styhi(j),stylo(j)+styhi(j))
         ystate(i,j,k) = cvmgt(zero,sty,lty)
         enddo
      else
         do j=jmin,jmax+1
            ystate(i,j,k) = cvmgp(stylo(j),styhi(j),vedge(i,j,k))
            ystate(i,j,k) = cvmgt(half*(stylo(j)+styhi(j)),ystate(i,j,k),
     &                           abs(vedge(i,j,k)).lt.eps)
         enddo
      endif
      place_to_break = 1
  200 continue
  210 continue
      endif

c     --------------------------------------------------------------
c     compute the zedge states
c     --------------------------------------------------------------
      if ((velpred.ne.1) .or. (n.eq.ZVEL)) then
      do 310 j = jmin,jmax
      do 300 i = imin,imax
      
      do k = kmin-1,kmax+1
         tr = 
     &        half*
     &        (uad(i+1,j,k) + uad(i,j,k)) *
     &        (xlo(i+1,j,k) - xlo(i,j,k)) / hx +
     &        half*
     &        (vad(i,j+1,k) + vad(i,j,k)) *
     &        (ylo(i,j+1,k) - ylo(i,j,k)) / hy

         stzlo(k+1)= s(i,j,k) + (half-dthz*w(i,j,k))*sz(i,j,k)
     &        - dth*tr
     &        + dth*tforces(i,j,k)
         stzhi(k)  = s(i,j,k) - (half+dthz*w(i,j,k))*sz(i,j,k)
     &        - dth*tr
     &        + dth*tforces(i,j,k)
      enddo

      if (bc(3,1).eq.EXT_DIR) then
         stzlo(kmin) = s(i,j,kmin-1)
         stzhi(kmin) = s(i,j,kmin-1)
      elseif (bc(3,1).eq.EXTRAP.or.bc(3,1).eq.HOEXTRAP
     &        .or.bc(3,1).eq.REFLECT_EVEN) then
         stzlo(kmin) = stzhi(kmin)
      elseif (bc(3,1).eq.REFLECT_ODD) then
         stzlo(kmin) = zero
         stzhi(kmin) = zero
      endif
      if (bc(3,2).eq.EXT_DIR) then
         stzlo(kmax+1) = s(i,j,kmax+1)
         stzhi(kmax+1) = s(i,j,kmax+1)
      elseif (bc(3,2).eq.EXTRAP.or.bc(3,2).eq.HOEXTRAP
     &        .or.bc(3,2).eq.REFLECT_EVEN) then
         stzhi(kmax+1) = stzlo(kmax+1)
      elseif (bc(3,2).eq.REFLECT_ODD) then
         stzlo(kmax+1) = zero
         stzhi(kmax+1) = zero
      endif

      if ( velpred .eq. 1 ) then
         do k = kmin,kmax+1
            ltz = stzlo(k) .le. zero  .and.  stzhi(k) .ge. zero
            ltz = ltz .or. (abs(stzlo(k)+stzhi(k)) .lt. eps)
            stz = cvmgp(stzlo(k),stzhi(k),stzlo(k)+stzhi(k))
            zstate(i,j,k) = cvmgt(zero,stz,ltz)
         enddo
      else
         do k = kmin,kmax+1
            zstate(i,j,k) = cvmgp(stzlo(k),stzhi(k),wedge(i,j,k))
            zstate(i,j,k) = cvmgt(half*(stzlo(k)+stzhi(k)),zstate(i,j,k),
     &                           abs(wedge(i,j,k)).lt.eps)
         enddo
      endif
      place_to_break = 1
  300 continue
  310 continue
      endif
c
c     end of routine
c
      return
      end
      

c=========================================================

      subroutine FORT_ADV_FORCING(
     &     aofs,DIMS(aofs),
     &     xflux,uedge,DIMS(uedge),
     &     areax,DIMS(ax),
     &     yflux,vedge,DIMS(vedge),
     &     areay,DIMS(ay),
     &     zflux,wedge,DIMS(wedge),
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
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      integer DIMDEC(az)
      REAL_T aofs(DIMV(aofs))
      REAL_T vol(DIMV(vol))
      REAL_T uedge(DIMV(uedge))
      REAL_T vedge(DIMV(vedge))
      REAL_T wedge(DIMV(wedge))
      REAL_T xflux(DIMV(uedge))
      REAL_T yflux(DIMV(vedge))
      REAL_T zflux(DIMV(wedge))
      REAL_T areax(DIMV(ax))
      REAL_T areay(DIMV(ay))
      REAL_T areaz(DIMV(az))
c
      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
c
c     if nonconservative initialize the advective tendency 
c     as -U*grad(S)
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
              
               enddo
            enddo
         enddo
      endif
c
c     convert edge states to fluxes
c
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax+1
               xflux(i,j,k) = xflux(i,j,k)*uedge(i,j,k)*areax(i,j,k)
            enddo
         enddo
      enddo
      do k = kmin,kmax
         do j = jmin,jmax+1
            do i = imin,imax
               yflux(i,j,k) = yflux(i,j,k)*vedge(i,j,k)*areay(i,j,k)
            enddo
         enddo
      enddo
      do k = kmin,kmax+1
         do j = jmin,jmax
            do i = imin,imax
               zflux(i,j,k) = zflux(i,j,k)*wedge(i,j,k)*areaz(i,j,k)
            enddo
         enddo
      enddo
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
               enddo
            enddo
         enddo
      else
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  aofs(i,j,k) = (
     &                 xflux(i+1,j,k) - xflux(i,j,k) +
     &                 yflux(i,j+1,k) - yflux(i,j,k) +
     &                 zflux(i,j,k+1) - zflux(i,j,k))/vol(i,j,k)
               enddo
            enddo
         enddo
      endif
c     
c     end of routine
c
      return
      end
      

c=========================================================

      subroutine FORT_SYNC_ADV_FORCING(
     &     sync,DIMS(sync),
     &     xflux,ucor,DIMS(ucor),
     &     areax,DIMS(ax),
     &     yflux,vcor,DIMS(vcor),
     &     areay,DIMS(ay),
     &     zflux,wcor,DIMS(wcor),
     &     areaz,DIMS(az),
     &     vol,DIMS(vol),
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
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      integer DIMDEC(az)
      REAL_T sync(DIMV(sync))
      REAL_T vol(DIMV(vol))
      REAL_T ucor(DIMV(ucor))
      REAL_T vcor(DIMV(vcor))
      REAL_T wcor(DIMV(wcor))
      REAL_T xflux(DIMV(ucor))
      REAL_T yflux(DIMV(vcor))
      REAL_T zflux(DIMV(wcor))
      REAL_T areax(DIMV(ax))
      REAL_T areay(DIMV(ay))
      REAL_T areaz(DIMV(az))
c
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
            enddo
         enddo
      enddo
      do k = kmin,kmax
         do j = jmin,jmax+1
            do i = imin,imax
               yflux(i,j,k) = yflux(i,j,k)*vcor(i,j,k)*areay(i,j,k)
            enddo
         enddo
      enddo
      do k = kmin,kmax+1
         do j = jmin,jmax
            do i = imin,imax
               zflux(i,j,k) = zflux(i,j,k)*wcor(i,j,k)*areaz(i,j,k)
            enddo
         enddo
      enddo
c
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               sync(i,j,k) = sync(i,j,k) + (
     &              xflux(i+1,j,k)-xflux(i,j,k) +
     &              yflux(i,j+1,k)-yflux(i,j,k) +
     &              zflux(i,j,k+1)-zflux(i,j,k) )/vol(i,j,k)
            enddo
         enddo
      enddo
c
c     end of routine
c
      return
      end


      
c=========================================================


      subroutine trans_xbc(
     &     s,DIMS(s),
     &     xlo,xhi,uad,DIMS(xx),
     &     lo,hi,n,xbc)
c
c     This subroutine processes boundary conditions on information
c     traced to cell faces in the x direction.  This is used for
c     computing velocities and edge states used in calculating
c     transverse derivatives
c
      integer DIMDEC(s)
      REAL_T s(DIMV(s))
      integer DIMDEC(xx)
      REAL_T xlo(DIMV(xx))
      REAL_T xhi(DIMV(xx))
      REAL_T uad(DIMV(xx))
      integer lo(SDIM), hi(SDIM)
      integer n
      integer xbc(SDIM,2)
c
      REAL_T stx
      logical ltest
      integer j,k
      integer imin,jmin,kmin,imax,jmax,kmax
c
      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
c
c     -------------- the lower x boundary
c
      if (xbc(1,1).eq.EXT_DIR) then
         if ( n .eq. XVEL ) then
            do j = jmin-1,jmax+1
            do k = kmin-1,kmax+1
               xlo(imin,j,k) = s(imin-1,j,k)
               xhi(imin,j,k) = s(imin-1,j,k)
            enddo
            enddo
         else
            do j = jmin-1,jmax+1
            do k = kmin-1,kmax+1
               ltest = (abs(uad(imin,j,k)).le.1.0e-10)
               stx   = cvmgt(xhi(imin,j,k),s(imin-1,j,k),ltest)
               xlo(imin,j,k) = stx
               xhi(imin,j,k) = stx
            enddo
            enddo
         endif
      elseif (xbc(1,1).eq.EXTRAP.or.xbc(1,1).eq.HOEXTRAP
     &        .or.xbc(1,1).eq.REFLECT_EVEN) then
         do j = jmin-1,jmax+1
         do k = kmin-1,kmax+1
            xlo(imin,j,k) = xhi(imin,j,k)
         enddo
         enddo
      elseif (xbc(1,1).eq.REFLECT_ODD) then
         do j = jmin-1,jmax+1
         do k = kmin-1,kmax+1
            xhi(imin,j,k) = zero
            xlo(imin,j,k) = zero
         enddo
         enddo
      endif
c
c     -------------- the upper x boundary
c
      if (xbc(1,2).eq.EXT_DIR) then
         if ( n .eq. XVEL ) then
            do j = jmin-1,jmax+1
            do k = kmin-1,kmax+1
               xlo(imax+1,j,k) = s(imax+1,j,k)
               xhi(imax+1,j,k) = s(imax+1,j,k)
            enddo
            enddo
         else
            do j = jmin-1,jmax+1
            do k = kmin-1,kmax+1
               ltest = (abs(uad(imax+1,j,k)).le.1.0e-10)
               stx   = cvmgt(xlo(imax+1,j,k),s(imax+1,j,k),ltest)
               xlo(imax+1,j,k) = stx
               xhi(imax+1,j,k) = stx
            enddo
            enddo
         endif
      elseif (xbc(1,2).eq.EXTRAP.or.xbc(1,2).eq.HOEXTRAP
     &        .or.xbc(1,2).eq.REFLECT_EVEN) then
         do j = jmin-1,jmax+1
         do k = kmin-1,kmax+1
            xhi(imax+1,j,k) = xlo(imax+1,j,k)
         enddo
         enddo
      elseif (xbc(1,2).eq.REFLECT_ODD) then
         do j = jmin-1,jmax+1
         do k = kmin-1,kmax+1
            xhi(imax+1,j,k) = zero
            xlo(imax+1,j,k) = zero
         enddo
         enddo
      endif
c
c     end of routine
c
      return
      end



c=========================================================


      subroutine trans_ybc(
     &     s,DIMS(s),
     &     ylo,yhi,vad,DIMS(yy),
     &     lo,hi,n,ybc)
c
c     This subroutine processes boundary conditions on information
c     traced to cell faces in the y direction.  This is used for
c     computing velocities and edge states used in calculating
c     transverse derivatives
c
      integer DIMDEC(s)
      REAL_T s(DIMV(s))
      integer DIMDEC(yy)
      REAL_T ylo(DIMV(yy))
      REAL_T yhi(DIMV(yy))
      REAL_T vad(DIMV(yy))
      integer lo(SDIM), hi(SDIM)
      integer n
      integer ybc(SDIM,2)
c
      REAL_T sty
      logical ltest
      integer i,k
      integer imin,jmin,kmin,imax,jmax,kmax
c
      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
c
c     -------------- the lower y boundary
c
      if (ybc(2,1).eq.EXT_DIR) then
         if ( n .eq. YVEL ) then
            do i = imin-1,imax+1
            do k = kmin-1,kmax+1
               ylo(i,jmin,k) = s(i,jmin-1,k)
               yhi(i,jmin,k) = s(i,jmin-1,k)
            enddo
            enddo
         else
            do i = imin-1,imax+1
            do k = kmin-1,kmax+1
               ltest = (abs(vad(i,jmin,k)).le.1.0e-10)
               sty   = cvmgt(yhi(i,jmin,k),s(i,jmin-1,k),ltest)
               ylo(i,jmin,k) = sty
               yhi(i,jmin,k) = sty
            enddo
            enddo
         endif
      elseif (ybc(2,1).eq.EXTRAP.or.ybc(2,1).eq.HOEXTRAP
     &        .or.ybc(2,1).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
         do k = kmin-1,kmax+1
            ylo(i,jmin,k) = yhi(i,jmin,k)
         enddo
         enddo
      elseif (ybc(2,1).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
         do k = kmin-1,kmax+1
            yhi(i,jmin,k) = zero
            ylo(i,jmin,k) = zero
         enddo
         enddo
      endif
c
c     -------------- the upper y boundary
c
      if (ybc(2,2).eq.EXT_DIR) then
         if ( n .eq. YVEL ) then
            do i = imin-1,imax+1
            do k = kmin-1,kmax+1
               ylo(i,jmax+1,k) = s(i,jmax+1,k)
               yhi(i,jmax+1,k) = s(i,jmax+1,k)
            enddo
            enddo
         else
            do i = imin-1,imax+1
            do k = kmin-1,kmax+1
               ltest = (abs(vad(i,jmax+1,k)).le.1.0e-10)
               sty   = cvmgt(ylo(i,jmax+1,k),s(i,jmax+1,k),ltest)
               ylo(i,jmax+1,k) = sty
               yhi(i,jmax+1,k) = sty
            enddo
            enddo
         endif
      elseif (ybc(2,2).eq.EXTRAP.or.ybc(2,2).eq.HOEXTRAP
     &        .or.ybc(2,2).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
         do k = kmin-1,kmax+1
            yhi(i,jmax+1,k) = ylo(i,jmax+1,k)
         enddo
         enddo
      elseif (ybc(2,2).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
         do k = kmin-1,kmax+1
            ylo(i,jmax+1,k) = zero
            yhi(i,jmax+1,k) = zero
         enddo
         enddo
      endif
c
c     end of routine
c
      return
      end

c=========================================================


      subroutine trans_zbc(
     &     s,DIMS(s),
     &     zlo,zhi,wad,DIMS(zz),
     &     lo,hi,n,zbc)
c
c     This subroutine processes boundary conditions on information
c     traced to cell faces in the z direction.  This is used for
c     computing velocities and edge states used in calculating
c     transverse derivatives
c
      integer DIMDEC(s)
      REAL_T s(DIMV(s))
      integer DIMDEC(zz)
      REAL_T zlo(DIMV(zz))
      REAL_T zhi(DIMV(zz))
      REAL_T wad(DIMV(zz))
      integer lo(SDIM), hi(SDIM)
      integer n
      integer zbc(SDIM,2)
c
      REAL_T stz
      logical ltest
      integer i,j
      integer imin,jmin,kmin,imax,jmax,kmax
c
      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
c
c     -------------- the lower z boundary
c
      if (zbc(3,1).eq.EXT_DIR) then
         if ( n .eq. ZVEL ) then
            do i = imin-1,imax+1
            do j = jmin-1,jmax+1
               zhi(i,j,kmin) = s(i,j,kmin-1)
               zlo(i,j,kmin) = s(i,j,kmin-1)
            enddo
            enddo
         else
            do i = imin-1,imax+1
            do j = jmin-1,jmax+1
               ltest = (abs(wad(i,j,kmin)).le.1.0e-10)
               stz   = cvmgt(zhi(i,j,kmin),s(i,j,kmin-1),ltest)
               zhi(i,j,kmin) = stz
               zlo(i,j,kmin) = stz
            enddo
            enddo
         endif
      elseif (zbc(3,1).eq.EXTRAP.or.zbc(3,1).eq.HOEXTRAP
     &        .or.zbc(3,1).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
         do j = jmin-1,jmax+1
            zlo(i,j,kmin) = zhi(i,j,kmin)
         enddo
         enddo
      elseif (zbc(3,1).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
         do j = jmin-1,jmax+1
            zhi(i,j,kmin) = zero
            zlo(i,j,kmin) = zero
         enddo
         enddo
      endif
c
c     -------------- the upper z boundary
c
      if (zbc(3,2).eq.EXT_DIR) then
         if ( n .eq. ZVEL ) then
            do i = imin-1,imax+1
            do j = jmin-1,jmax+1
               zlo(i,j,kmax+1) = s(i,j,kmax+1)
               zhi(i,j,kmax+1) = s(i,j,kmax+1)
            enddo
            enddo
         else
            do i = imin-1,imax+1
            do j = jmin-1,jmax+1
               ltest = (abs(wad(i,j,kmax+1)).le.1.0e-10)
               stz   = cvmgt(zlo(i,j,kmax+1),s(i,j,kmax+1),ltest)
               zhi(i,j,kmax+1) = stz
               zlo(i,j,kmax+1) = stz
            enddo
            enddo
         endif
      elseif (zbc(3,2).eq.EXTRAP.or.zbc(3,2).eq.HOEXTRAP
     &        .or.zbc(3,2).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
         do j = jmin-1,jmax+1
            zhi(i,j,kmax+1) = zlo(i,j,kmax+1)
         enddo
         enddo
      elseif (zbc(3,2).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
         do j = jmin-1,jmax+1
            zlo(i,j,kmax+1) = zero
            zhi(i,j,kmax+1) = zero
         enddo
         enddo
      endif
c
c     end of routine
c
      return
      end




c=========================================================


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
#include "GODUNOV.FH"
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
c
      integer imin,jmin,kmin,imax,jmax,kmax,i,j,k
      integer ng
      REAL_T dpls,dmin,ds
      REAL_T del,slim,sflg
      integer cen,lim,flag,fromm
      parameter( cen = 1   )
      parameter( lim = 2   )
      parameter( flag = 3  )
      parameter( fromm = 4 )
c
      ng = lo(1) - ARG_L1(s)
      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

c     -------------------------------------------------
c     COMPUTE 0TH order slopes
c     -------------------------------------------------
      if (slope_order.eq.1) then
        if (ng .lt. 1) then
           write(6,*) "SLOPE_3D: not enough bndry cells for 1st order"
	   stop
        endif
        do k = kmin-1, kmax+1
           do j = jmin-1, jmax+1 
              do i = imin-1, imax+1
                 slx(i,j,k) = zero
                 sly(i,j,k) = zero
                 slz(i,j,k) = zero
              enddo
           enddo
        enddo
        return
      endif

c     -------------------------------------------------
c     COMPUTE 4TH order slopes
c     -------------------------------------------------
      if(slope_order.eq.4)then
        if (ng .lt. 3) then
           write(6,*) "SLOPE_3D: not enough bndry cells for 4th order"
	   stop
        endif
c
c     ------------------------ x slopes
c
        if ( (dir.eq.XVEL) .or. (dir .gt. 5 ) ) then
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
     &           min(slxscr(i,lim),abs(slxscr(i,cen)))
           enddo
           do i = imin-1,imax+1
              ds = two * two3rd * slxscr(i,cen) - 
     &           sixth * (slxscr(i+1,fromm) + slxscr(i-1,fromm))
              slx(i,j,k) = slxscr(i,flag)*min(abs(ds),slxscr(i,lim))
           enddo
        enddo
        enddo
c
        if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
           do k = kmin-1,kmax+1
           do j = jmin-1,jmax+1
              del  = -sixteen/fifteen*s(imin-1,j,k) + half*s(imin,j,k) + 
     &           two3rd*s(imin+1,j,k) -  tenth*s(imin+2,j,k)
              dmin = two*(s(imin  ,j,k)-s(imin-1,j,k))
              dpls = two*(s(imin+1,j,k)-s(imin  ,j,k))
              slim = min(abs(dpls), abs(dmin))
              slim = cvmgp(slim, zero, dpls*dmin)
              sflg = sign(one,del)
	      slx(imin-1,j,k) = zero
              slx(imin,  j,k) = sflg*min(slim,abs(del))
           enddo
           enddo
        endif
        if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
           do k = kmin-1,kmax+1
           do j = jmin-1,jmax+1
              del  = -( -sixteen/fifteen*s(imax+1,j,k) + half*s(imax,j,k) + 
     &           two3rd*s(imax-1,j,k) - tenth*s(imax-2,j,k) )
              dmin = two*(s(imax  ,j,k)-s(imax-1,j,k))
              dpls = two*(s(imax+1,j,k)-s(imax  ,j,k))
              slim = min(abs(dpls), abs(dmin))
              slim = cvmgp(slim, zero, dpls*dmin)
              sflg = sign(one,del)
              slx(imax,  j,k) = sflg*min(slim,abs(del))
	      slx(imax+1,j,k) = zero
           enddo
           enddo
        endif
        endif
c
c     ------------------------ y slopes
c
        if ( (dir.eq.YVEL) .or. (dir .gt. 5 ) ) then
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
     &             min(slyscr(j,lim),abs(slyscr(j,cen)))
           enddo
           do j = jmin-1,jmax+1
              ds = two * two3rd * slyscr(j,cen) - 
     &             sixth * (slyscr(j+1,fromm) + slyscr(j-1,fromm))
              sly(i,j,k) = slyscr(j,flag)*min(abs(ds),slyscr(j,lim))
           enddo
        enddo
        enddo
c
        if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
           do k = kmin-1,kmax+1
           do i = imin-1,imax+1
              del  = -sixteen/fifteen*s(i,jmin-1,k) + half*s(i,jmin,k) + 
     &             two3rd*s(i,jmin+1,k) - tenth*s(i,jmin+2,k)
              dmin = two*(s(i,jmin  ,k)-s(i,jmin-1,k))
              dpls = two*(s(i,jmin+1,k)-s(i,jmin  ,k))
              slim = min(abs(dpls), abs(dmin))
              slim = cvmgp(slim, zero, dpls*dmin)
              sflg = sign(one,del)
	      sly(i,jmin-1,k) = zero
              sly(i,jmin,  k) = sflg*min(slim,abs(del))
	   enddo
           enddo
        endif
        if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
           do k = kmin-1,kmax+1
           do i = imin-1,imax+1
              del  = -( -sixteen/fifteen*s(i,jmax+1,k) + half*s(i,jmax,k) +
     &             two3rd*s(i,jmax-1,k) - tenth*s(i,jmax-2,k) )
              dmin = two*(s(i,jmax  ,k)-s(i,jmax-1,k))
              dpls = two*(s(i,jmax+1,k)-s(i,jmax  ,k))
              slim = min(abs(dpls), abs(dmin))
              slim = cvmgp(slim, zero, dpls*dmin)
              sflg = sign(one,del)
              sly(i,jmax, k)  = sflg*min(slim,abs(del))
	      sly(i,jmax+1,k) = zero
           enddo
	   enddo
        endif
        endif
c
c     ------------------------ z slopes
c
        if ( (dir.eq.ZVEL) .or. (dir .gt. 5 ) ) then
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
     &            min(slzscr(k,lim),abs(slzscr(k,cen)))
          enddo
          do k = kmin-1,kmax+1
             ds = two * two3rd * slzscr(k,cen) -
     &            sixth * (slzscr(k+1,fromm) + slzscr(k-1,fromm))
             slz(i,j,k) = slzscr(k,flag)*min(abs(ds),slzscr(k,lim))
          enddo
        enddo
        enddo
c
        if (bc(3,1) .eq. EXT_DIR .or. bc(3,1) .eq. HOEXTRAP) then
           do j = jmin-1,jmax+1
	   do i = imin-1,imax+1
              del  = -sixteen/fifteen*s(i,j,kmin-1) + half*s(i,j,kmin) +
     &             two3rd*s(i,j,kmin+1) - tenth*s(i,j,kmin+2)
              dmin = two*(s(i,j,kmin  )-s(i,j,kmin-1))
              dpls = two*(s(i,j,kmin+1)-s(i,j,kmin  ))
              slim = min(abs(dpls), abs(dmin))
              slim = cvmgp(slim, zero, dpls*dmin)
              sflg = sign(one,del)
              slz(i,j,kmin-1) = zero
              slz(i,j,kmin  ) = sflg*min(slim,abs(del))
           enddo
           enddo
        endif
        if (bc(3,2) .eq. EXT_DIR .or. bc(3,2) .eq. HOEXTRAP) then
           do j = jmin-1,jmax+1
	   do i = imin-1,imax+1
              del  = sixteen/fifteen*s(i,j,kmax+1) - half*s(i,j,kmax)
     &           - two3rd*s(i,j,kmax-1) + tenth*s(i,j,kmax-2)
              dmin = two*(s(i,j,kmax  )-s(i,j,kmax-1))
              dpls = two*(s(i,j,kmax+1)-s(i,j,kmax  ))
              slim = min(abs(dpls), abs(dmin))
              slim = cvmgp(slim, zero, dpls*dmin)
              sflg = sign(one,del)
              slz(i,j,kmax  ) = sflg*min(slim,abs(del))
              slz(i,j,kmax+1) = zero
           enddo
           enddo
        endif
        endif
      endif
c
c     end of routine
c   
      return
      end


c=========================================================

      subroutine FORT_SCALMINMAX(s,DIMS(s),sn,DIMS(sn),
     &     smin,smax,DIMS(smin),
     &     lo,hi,bc)
c
c     correct an advected field for under/over shoots
c
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
c
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
      endif
      if (bc(1,2).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
         s(imax+1,jmin-1,kmin-1) = s(imax,jmin,kmin)
      endif
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
         s(imin-1,jmax+1,kmin-1) = s(imin,jmax,kmin)
      endif
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
         s(imax+1,jmax+1,kmin-1) = s(imax,jmax,kmin)
      endif
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
         s(imin-1,jmin-1,kmax+1) = s(imin,jmin,kmax)
      endif
      if (bc(1,2).ne.INT_DIR .or. bc(2,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
         s(imax+1,jmin-1,kmax+1) = s(imax,jmin,kmax)
      endif
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
         s(imin-1,jmax+1,kmax+1) = s(imin,jmax,kmax)
      endif
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
         s(imax+1,jmax+1,kmax+1) = s(imax,jmax,kmax)
      endif
c
c     correct the 12 edges
c
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
        do k = kmin,kmax
          s(imin-1,jmin-1,k) = s(imin,jmin,k)
        enddo
      endif
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
        do k = kmin,kmax
          s(imax+1,jmin-1,k) = s(imax,jmin,k)
        enddo
      endif
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
        do k = kmin,kmax
          s(imin-1,jmax+1,k) = s(imin,jmax,k)
        enddo
      endif
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
        do k = kmin,kmax
          s(imax+1,jmax+1,k) = s(imax,jmax,k)
        enddo
      endif
      if (bc(1,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do j = jmin,jmax
          s(imin-1,j,kmin-1) = s(imin,j,kmin)
        enddo
      endif
      if (bc(1,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do j = jmin,jmax
          s(imax+1,j,kmin-1) = s(imax,j,kmin)
        enddo
      endif
      if (bc(1,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do j = jmin,jmax
          s(imin-1,j,kmax+1) = s(imin,j,kmax)
        enddo
      endif
      if (bc(1,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do j = jmin,jmax
          s(imax+1,j,kmax+1) = s(imax,j,kmax)
        enddo
      endif
      if (bc(2,1).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do i = imin,imax
          s(i,jmin-1,kmin-1) = s(i,jmin,kmin)
        enddo
      endif
      if (bc(2,2).ne.INT_DIR .or. bc(3,1).ne.INT_DIR) then
        do i = imin,imax
          s(i,jmax+1,kmin-1) = s(i,jmax,kmin)
        enddo
      endif
      if (bc(2,1).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do i = imin,imax
          s(i,jmin-1,kmax+1) = s(i,jmin,kmax)
        enddo
      endif
      if (bc(2,2).ne.INT_DIR .or. bc(3,2).ne.INT_DIR) then
        do i = imin,imax
          s(i,jmax+1,kmax+1) = s(i,jmax,kmax)
        enddo
      endif

c     ::::: compute min/max a slab at a time
c     ::::: compute min and max of neighbors on kmin-1 slab
      km = 0
      kk = 1
      kp = 2

      k = kmin-1
      do j = jmin, jmax         
      do i = imin, imax
	 smin(i,j,km) = min(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &                     s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &                     s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
	 smax(i,j,km) = max(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &                     s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &                     s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
      enddo         
      enddo         

c     ::::: compute min and max of neighbors on kmin slab
      k = kmin
      do j = jmin, jmax         
      do i = imin, imax
	 smin(i,j,kk) = min(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &                     s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &                     s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
	 smax(i,j,kk) = max(s(i-1,j-1,k),s(i,j-1,k),s(i+1,j-1,k),
     &                     s(i-1,j  ,k),s(i,j  ,k),s(i+1,j  ,k),
     &                     s(i-1,j+1,k),s(i,j+1,k),s(i+1,j+1,k))
      enddo         
      enddo        

      do k = kmin, kmax         
c        ::::: compute min and max of neighbors on k+1 slab
         do j = jmin, jmax     
         do i = imin, imax   
	 smin(i,j,kp) = min(s(i-1,j-1,k+1),s(i,j-1,k+1),s(i+1,j-1,k+1),
     &                      s(i-1,j  ,k+1),s(i,j  ,k+1),s(i+1,j  ,k+1),
     &                      s(i-1,j+1,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))
	 smax(i,j,kp) = max(s(i-1,j-1,k+1),s(i,j-1,k+1),s(i+1,j-1,k+1),
     &                      s(i-1,j  ,k+1),s(i,j  ,k+1),s(i+1,j  ,k+1),
     &                      s(i-1,j+1,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))

c        ::::: compute min/max of cell
         smn = min(smin(i,j,km),smin(i,j,kk),smin(i,j,kp))
         smx = max(smax(i,j,km),smax(i,j,kk),smax(i,j,kp))
         sn(i,j,k) = max(sn(i,j,k),smn)
         sn(i,j,k) = min(sn(i,j,k),smx)
         
         enddo
         enddo

c        ::::: roll indices for next slab
         km = mod(km+1,3)
         kk = mod(kk+1,3)
         kp = mod(kp+1,3)
      enddo
c
c     end of routine
c	 
      return
      end




c=========================================================

      subroutine FORT_SUM_TF_GP(
     &     tforces,DIMS(tf),
     &     gp,DIMS(gp),
     &     rho,DIMS(rho),
     &     lo,hi )
c
c     sum pressure forcing into tforces
c
      integer i, j, k, n
      integer DIMDEC(tf)
      integer DIMDEC(visc)
      integer DIMDEC(gp)
      integer DIMDEC(rho)
      integer lo(SDIM), hi(SDIM)
      REAL_T tforces(DIMV(tf),SDIM)
      REAL_T gp(DIMV(gp),SDIM)
      REAL_T rho(DIMV(rho))
c
      do n = 1, SDIM
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,k,n) = (
     &            tforces(i,j,k,n) - gp(i,j,k,n)
     &                 )/rho(i,j,k)
               enddo
            enddo
         enddo
      enddo
c     
c     end of routine
c     
      return
      end



c=========================================================

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
c
      do n = 1, SDIM
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,k,n) = (
     &                 tforces(i,j,k,n)
     &                 +  visc(i,j,k,n)
     &                 -    gp(i,j,k,n) )/rho(i,j,k)
               enddo
            enddo
         enddo
      enddo
c     
c     end of routine
c     
      return
      end

c=========================================================

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
      integer nvar, iconserv
      integer lo(SDIM), hi(SDIM)
      integer i, j, k, n
c
      integer DIMDEC(S)
      integer DIMDEC(tf)
      integer DIMDEC(divu)
      integer DIMDEC(rho)
c
      REAL_T S(DIMV(S),nvar)
      REAL_T tforces(DIMV(tf),nvar)
      REAL_T divu(DIMV(divu))
      REAL_T rho(DIMV(rho))
c
      if ( iconserv .eq. 1 ) then

         do n = 1, nvar
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tforces(i,j,k,n) = 
     &               tforces(i,j,k,n) - S(i,j,k,n)*divu(i,j,k)
                  enddo
               enddo
            enddo
         enddo

      else

         do n = 1, nvar
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tforces(i,j,k,n) = 
     &               tforces(i,j,k,n)/rho(i,j,k)
                  enddo
               enddo
            enddo
         enddo

      endif
c     
c     end of routine
c     
      return
      end


c=========================================================

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
      integer nvar, iconserv
      integer lo(SDIM), hi(SDIM)
      integer i, j, k, n
c
      integer DIMDEC(S)
      integer DIMDEC(tf)
      integer DIMDEC(divu)
      integer DIMDEC(visc)
      integer DIMDEC(rho)
c
      REAL_T S(DIMV(S),nvar)
      REAL_T tforces(DIMV(tf),nvar)
      REAL_T divu(DIMV(divu))
      REAL_T visc(DIMV(visc),nvar)
      REAL_T rho(DIMV(rho))
c
      if ( iconserv .eq. 1 ) then

         do n = 1, nvar
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tforces(i,j,k,n) = 
     &                    tforces(i,j,k,n)
     &                    +  visc(i,j,k,n)
     &                    -     S(i,j,k,n)*divu(i,j,k)
                  enddo
               enddo
            enddo
         enddo

      else

         do n = 1, nvar
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     tforces(i,j,k,n) = (
     &                    tforces(i,j,k,n)
     &                    +  visc(i,j,k,n) )/rho(i,j,k)
                  enddo
               enddo
            enddo
         enddo

      endif
c     
c     end of routine
c     
      return
      end


c=========================================================


      subroutine FORT_UPDATE_TF(
     &     s,       DIMS(s),
     &     sn,      DIMS(sn),
     &     tforces, DIMS(tf),
     &     lo,hi,dt,nvar)
c
c     update a field with a forcing term
c
      integer i, j, k, n, nvar
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(tf)
      integer lo(SDIM), hi(SDIM)
      REAL_T dt
      REAL_T s(DIMV(s),nvar)
      REAL_T sn(DIMV(sn),nvar)
      REAL_T tforces(DIMV(tf),nvar)
c
      do n = 1,nvar
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  sn(i,j,k,n) = s(i,j,k,n)
     &                 + dt*tforces(i,j,k,n)
               enddo
            enddo
         enddo
      enddo
c
c     end of routine
c
      return
      end



c=========================================================

      subroutine FORT_CORRECT_TF(
     &     ss,  sp,  DIMS(ss),
     &     tfs, tfn, DIMS(tfs),
     &     lo,hi,dt,nvar)
c
c     correct 1st order rk to second-order
c
      integer i, j, k, n, nvar
      integer lo(SDIM), hi(SDIM)
      REAL_T dt,hdt
c
      integer DIMDEC(ss)
      integer DIMDEC(tfs)
      REAL_T  ss(DIMV(ss),nvar)
      REAL_T  sp(DIMV(ss),nvar)
      REAL_T tfs(DIMV(tfs),nvar)
      REAL_T tfn(DIMV(tfs),nvar)
c
      hdt = half*dt
      do n = 1,nvar
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  sp(i,j,k,n) = ss(i,j,k,n)
     &                 + hdt*(tfs(i,j,k,n)-tfn(i,j,k,n))
               enddo
            enddo
         enddo
      enddo
c
c     end of routine
c
      return
      end



c=========================================================


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
c
      do n = 1,nvar
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  sn(i,j,k,n) = s(i,j,k,n)
     &                 - dt*aofs(i,j,k,n)
     &                 + dt*tforces(i,j,k,n)
               enddo
            enddo
         enddo
      enddo
c
c     end of routine
c
      return
      end

c=========================================================

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
c
      do n = 1, SDIM
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  un(i,j,k,n) = u(i,j,k,n)
     &                 - dt*   aofs(i,j,k,n)
     &                 + dt*tforces(i,j,k,n)/rho(i,j,k)
     &                 - dt*     gp(i,j,k,n)/rho(i,j,k)
               enddo
            enddo
         enddo
      enddo
c     
c     end of routine
c
      return
      end

