c $Id: GODUNOV_2D.F,v 1.3 1997-07-18 20:58:24 car Exp $
#undef BL_LANG_CC
#define BL_LANG_FORT

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "GODUNOV_F.H"
#include "DIMS.H"

#define SDIM 2
#define XVEL 1
#define YVEL 2

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
      REAL_T    u, v
      REAL_T    small
      REAL_T    dt_start
      REAL_T    tforce1
      REAL_T    tforce2
      integer   i, j
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
      small   = 1.0d-8
      u       = zero
      v       = zero
      tforce1 = zero
      tforce2 = zero
c
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
         u = max(u,abs(vel(i,j,1)))
         v = max(v,abs(vel(i,j,2)))
         tforce1 = max(tforce1,abs(tforces(i,j,1)/rho(i,j)))
         tforce2 = max(tforce2,abs(tforces(i,j,2)/rho(i,j)))
      enddo
      enddo

      u_max(1) = u
      u_max(2) = v

      dt_start = 1.0e+20
      dt = dt_start

      if (u .gt. small) dt = min(dt,dx(1)/u)
      if (v .gt. small) dt = min(dt,dx(2)/v)

      if (tforce1 .gt. small) then
        dt  = min(dt,sqrt(two*dx(1)/tforce1))
      endif

      if (tforce2 .gt. small) then
        dt  = min(dt,sqrt(two*dx(2)/tforce2))
      endif

      if (dt .eq. dt_start) dt = min(dx(1),dx(2))

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
     &     rho,DIMS(rho),
     &     lo,hi,dt,dx,cflmax,u_max,verbose )
c
c     This subroutine computes the extrema of the density
c     and velocities at cell centers
c
      integer DIMDEC(u)
      integer DIMDEC(v)
      integer DIMDEC(rho)
      REAL_T  u(DIMV(u))
      REAL_T  v(DIMV(v))
      REAL_T  rho(DIMV(rho))
      integer lo(SDIM)
      integer hi(SDIM)
      REAL_T  dt
      REAL_T  dx(SDIM)
      REAL_T  cflmax
      REAL_T  u_max(SDIM)
      integer verbose
c
      REAL_T  hx, hy
      REAL_T  umax, vmax, rhomax
      REAL_T  umin, vmin, rhomin
      integer imin, imax, jmin, jmax
      integer i, j
c
      hx   = dx(1)
      hy   = dx(2)
      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
      umax = -1.e30
      vmax = -1.e30
      umin =  1.e30
      vmin =  1.e30
      rhomax = -1.e30
      rhomin =  1.e30
c
      do j = jmin, jmax
         do i = imin, imax
            umax = max(umax,u(i,j))
            umin = min(umin,u(i,j))
            vmax = max(vmax,v(i,j))
            vmin = min(vmin,v(i,j))
            rhomax = max(rhomax,rho(i,j))
            rhomin = min(rhomin,rho(i,j))
         enddo
      enddo

      u_max(1) = max(abs(umax), abs(umin))
      u_max(2) = max(abs(vmax), abs(vmin))
      cflmax   = dt*max(u_max(1)/hx,u_max(2)/hy)

      if(verbose.eq.1)then
        write(6,1000) umax,umin,u_max(1) 
        write(6,1001) vmax,vmin,u_max(2) 
        write(6,1002) rhomax,rhomin
#ifndef	BL_NO_FORT_FLUSH
        call flush(6)
#endif
      endif

 1000 format(' U  MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1001 format(' V  MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1002 format('RHO MAX/MIN      ',e21.14,2x,e21.14)

c
c     end of routine
c
      return
      end



c=========================================================

      subroutine FORT_TEST_UMAC_RHO(
     &     umac,DIMS(umac),
     &     vmac,DIMS(vmac),
     &     rho,DIMS(rho),
     &     lo,hi,dt,dx,cflmax,u_max)
c
c     This subroutine computes the extrema of the density
c     and mac edge velocities
c
      integer DIMDEC(umac)
      integer DIMDEC(vmac)
      integer DIMDEC(rho)
      integer imin, imax, jmin, jmax
      integer i, j
      integer lo(SDIM),hi(SDIM)
      REAL_T  dt, dx(SDIM), u_max(SDIM), cflmax
      REAL_T  hx, hy
      REAL_T  umax, vmax, rhomax
      REAL_T  umin, vmin, rhomin
      REAL_T  umac(DIMV(umac))
      REAL_T  vmac(DIMV(vmac))
      REAL_T  rho(DIMV(rho))
c
      hx   = dx(1)
      hy   = dx(2)
      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
      umax = -1.e30
      vmax = -1.e30
      umin =  1.e30
      vmin =  1.e30
      rhomax = -1.e30
      rhomin =  1.e30
c
      do j = jmin, jmax
         do i = imin, imax+1
            umax = max(umax,umac(i,j))
            umin = min(umin,umac(i,j))
         enddo
      enddo
      do j = jmin, jmax+1
         do i = imin, imax
            vmax = max(vmax,vmac(i,j))
            vmin = min(vmin,vmac(i,j))
         enddo
      enddo
      do j = jmin, jmax
         do i = imin, imax
            rhomax = max(rhomax,rho(i,j))
            rhomin = min(rhomin,rho(i,j))
         enddo
      enddo

      u_max(1) = max(abs(umax), abs(umin))
      u_max(2) = max(abs(vmax), abs(vmin))
      cflmax   = dt*max(u_max(1)/hx,u_max(2)/hy)

      write(6,1000) umax,umin,u_max(1) 
      write(6,1001) vmax,vmin,u_max(2) 
      write(6,1002) rhomax,rhomin

 1000 format('UMAC MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1001 format('VMAC MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1002 format('RHO  MAX/MIN      ',e21.14,2x,e21.14)

#ifndef	BL_NO_FORT_FLUSH
      call flush(6)
#endif
c
c     end of routine
c
      return
      end




c=========================================================

      subroutine FORT_TRANSVEL(
     &     u, ulo, uhi, sx, ubc, slxscr,
     &     v, vlo, vhi, sy, vbc, slyscr,
     &     DIMS(s), DIMS(work),
     &     lo,hi,dt,dx,use_minion,tforces)
c
c     This subroutine computes the advective velocities used in
c     the transverse derivatives of the Godunov box
c
      implicit none
      integer i,j
      integer ubc(SDIM,2),vbc(SDIM,2)
      integer lo(SDIM),hi(SDIM)
      integer imin,jmin,imax,jmax
      REAL_T hx, hy, dt, dth, dthx, dthy,dx(SDIM)
      REAL_T uad,vad
      REAL_T eps
      logical ltm
      parameter( eps = 1.0e-6 )
      integer DIMDEC(s)
      integer DIMDEC(work)
      REAL_T  u(DIMV(s))
      REAL_T  v(DIMV(s))
      REAL_T ulo(DIMV(work)),uhi(DIMV(work))
      REAL_T vlo(DIMV(work)),vhi(DIMV(work))
      REAL_T sx(DIMV(work))
      REAL_T sy(DIMV(work))
      REAL_T slxscr(DIM1(s), 4)
      REAL_T slyscr(DIM2(s), 4)

      integer use_minion
      REAL_T tforces(DIMV(work),SDIM)
c
      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)
      hx   = dx(1)
      hy   = dx(2)
      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)
c
c     --------------------------------------------------------------
c     compute the x transverse velocities
c     --------------------------------------------------------------
c
      call FORT_SLOPES(u,DIMS(s),
     &     sx,sy,DIMS(work),lo,hi,slxscr,slyscr,ubc)

      do i = imin, imax+1
         do j = jmin-1,jmax+1
            ulo(i,j) = u(i-1,j) + (half  - dthx*u(i-1,j))*sx(i-1,j)
            uhi(i,j) = u(i,j)   + (-half - dthx*u(i,  j))*sx(i,j)
         enddo
      enddo

      if(use_minion.eq.1)then
        do i = imin, imax+1
          do j = jmin-1,jmax+1
            ulo(i,j) = ulo(i,j) + dth*tforces(i-1,j,1)
            uhi(i,j) = uhi(i,j) + dth*tforces(i,  j,1)
          enddo
        enddo
      endif

      call trans_xbc(u,DIMS(s),
     &     ulo,uhi,ulo,DIMS(work),lo,hi,XVEL,ubc)

      do i = imin,imax+1
         do j = jmin-1,jmax+1
            uad = cvmgp(ulo(i,j),uhi(i,j),ulo(i,j)+uhi(i,j))
            ltm = ulo(i,j) .le. zero  .and.  uhi(i,j) .ge. zero
            ltm = ltm .or. (abs(ulo(i,j)+uhi(i,j)) .lt. eps)
            ulo(i,j) = cvmgt(zero,uad,ltm)
         enddo
      enddo

c     --------------------------------------------------------------
c     compute the y transverse velocities
c     --------------------------------------------------------------
      call FORT_SLOPES(v,DIMS(s),
     &     sx,sy,DIMS(work),lo,hi,slxscr,slyscr,vbc)

      do i = imin-1,imax+1
         do j = jmin,jmax+1
            vlo(i,j) = v(i,j-1) + (half  - dthy*v(i,j-1))*sy(i,j-1)
            vhi(i,j) = v(i,j)   + (-half - dthy*v(i,j  ))*sy(i,j)
         enddo
      enddo

      if(use_minion.eq.1)then
        do i = imin-1, imax+1
          do j = jmin,jmax+1
            vlo(i,j) = vlo(i,j) + dth*tforces(i,j-1,2)
            vhi(i,j) = vhi(i,j) + dth*tforces(i,j,  2)
          enddo
        enddo
      endif

      call trans_ybc(v,DIMS(s),
     &     vlo,vhi,vlo,DIMS(work),lo,hi,YVEL,vbc)

      do i = imin-1,imax+1
         do j = jmin,jmax+1
            vad = cvmgp(vlo(i,j),vhi(i,j),vlo(i,j)+vhi(i,j))
            ltm = vlo(i,j) .le. zero  .and.  vhi(i,j) .ge. zero
            ltm = ltm .or. (abs(vlo(i,j)+vhi(i,j)) .lt. eps)
            vlo(i,j) = cvmgt(zero,vad,ltm)
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

     &     DIMS(work),
     &     bc,lo,hi,dt,dx,n,velpred, use_minion)
c
c     This subroutine computes edges states, right now it uses
c     a lot of memory, but there becomes a trade off between
c     simplicity-efficiency in the new way of computing states
c     and complexity in the old way.  By eliminating loops over
c     state components though, the new way uses much less memory.
c
      integer i,j,n,velpred
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer imin,jmin,imax,jmax
      REAL_T place_to_break
      REAL_T hx, hy, dt, dth, dthx, dthy, dx(SDIM)
      REAL_T tr,stx,sty,fu,fv
      REAL_T eps
      logical ltx,lty
      parameter( eps = 1.0e-6 )
c
      integer DIMDEC(s)
      integer DIMDEC(work)
      integer DIMDEC(uedge)
      integer DIMDEC(vedge)
c
      REAL_T s(DIMV(s))
      REAL_T u(DIMV(s))
      REAL_T v(DIMV(s))
      REAL_T stxlo(DIM1(s)),stxhi(DIM1(s)),slxscr(DIM1(s),4)
      REAL_T stylo(DIM2(s)),styhi(DIM2(s)),slyscr(DIM2(s),4)
c
      REAL_T uedge(DIMV(uedge)), xstate(DIMV(uedge))
      REAL_T vedge(DIMV(vedge)), ystate(DIMV(vedge))
c
      REAL_T xlo(DIMV(work)), xhi(DIMV(work))
      REAL_T ylo(DIMV(work)), yhi(DIMV(work))
      REAL_T  sx(DIMV(work)), uad(DIMV(work)) 
      REAL_T  sy(DIMV(work)), vad(DIMV(work)) 
      REAL_T tforces(DIMV(work))

      integer use_minion
c
      dth  = half*dt
      dthx = half*dt/dx(1)
      dthy = half*dt/dx(2)
      hx = dx(1)
      hy = dx(2)
      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

c     --------------------------------------------------------------
c     compute the slopes
c     --------------------------------------------------------------
      call FORT_SLOPES(s,DIMS(s),
     &     sx,sy,DIMS(work),lo,hi,slxscr,slyscr,bc)

c     --------------------------------------------------------------
c     trace the state to the cell edges
c     --------------------------------------------------------------
      do i = imin, imax+1
      do j = jmin-1,jmax+1
        xlo(i,j) = s(i-1,j) + (half - dthx*u(i-1,j))*sx(i-1,j)
        xhi(i,j) = s(i  ,j) - (half + dthx*u(i  ,j))*sx(i  ,j)
      enddo
      enddo

      if(use_minion.eq.1)then
        do i = imin, imax+1
          do j = jmin-1,jmax+1
            xlo(i,j) = xlo(i,j) + dth*tforces(i-1,j)
            xhi(i,j) = xhi(i,j) + dth*tforces(i,  j)
          enddo
        enddo
      endif

      call trans_xbc(s,DIMS(s),
     &     xlo,xhi,uad,DIMS(work),lo,hi,n,bc)
      
      do i = imin, imax+1
      do j = jmin-1,jmax+1
        fu  = cvmgt(zero,one,abs(uad(i,j)).lt.eps)
        stx = cvmgp(xlo(i,j),xhi(i,j),uad(i,j))
        xlo(i,j) = fu*stx + (one - fu)*half*(xhi(i,j)+xlo(i,j))
      enddo
      enddo
c     --------------------------------------------------------------
      do j = jmin,jmax+1
      do i = imin-1,imax+1
        ylo(i,j) = s(i,j-1) + (half - dthy*v(i,j-1))*sy(i,j-1)
        yhi(i,j) = s(i,j  ) - (half + dthy*v(i,j  ))*sy(i,j  )
      enddo
      enddo

      if(use_minion.eq.1)then
        do i = imin-1, imax+1
          do j = jmin,jmax+1
            ylo(i,j) = ylo(i,j) + dth*tforces(i,j-1)
            yhi(i,j) = yhi(i,j) + dth*tforces(i,j)
          enddo
        enddo
      endif

      call trans_ybc(s,DIMS(s),
     &     ylo,yhi,vad,DIMS(work),lo,hi,n,bc)

      do j = jmin,jmax+1
      do i = imin-1,imax+1
            fv  = cvmgt(zero,one,abs(vad(i,j)).lt.eps)
            sty = cvmgp(ylo(i,j),yhi(i,j),vad(i,j))
            ylo(i,j) = fv*sty + (one - fv)*half*(yhi(i,j)+ylo(i,j))
      enddo
      enddo

c     --------------------------------------------------------------
c     compute the xedge state
c     --------------------------------------------------------------

      if ((velpred.ne.1) .or. (n.eq.XVEL)) then
      do 100 j = jmin,jmax
      do i = imin-1,imax+1
         tr = half*
     &        (vad(i,j+1)+vad(i,j))*
     &        (ylo(i,j+1)-ylo(i,j))/hy
         
         stxlo(i+1) = s(i,j) + (half-dthx*u(i,j))*sx(i,j)
     &        - dth*tr
     &        + dth*tforces(i,j)
         stxhi(i  ) = s(i,j) - (half+dthx*u(i,j))*sx(i,j)
     &        - dth*tr
     &        + dth*tforces(i,j)
      enddo

      if (bc(1,1).eq.EXT_DIR) then
         stxhi(imin) = s(imin-1,j)
         stxlo(imin) = s(imin-1,j)
      elseif (bc(1,1).eq.EXTRAP.or.bc(1,1).eq.HOEXTRAP
     &        .or.bc(1,1).eq.REFLECT_EVEN) then
         stxlo(imin) = stxhi(imin)
      elseif (bc(1,1).eq.REFLECT_ODD) then
         stxhi(imin) = zero
         stxlo(imin) = stxhi(imin)
      endif
      if (bc(1,2).eq.EXT_DIR) then
        stxlo(imax+1) = s(imax+1,j)
        stxhi(imax+1) = s(imax+1,j)
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
            xstate(i,j) = cvmgt(zero,stx,ltx)
         enddo
      else
         do i = imin, imax+1
            xstate(i,j) = cvmgp(stxlo(i),stxhi(i),uedge(i,j))
            xstate(i,j) = cvmgt(half*(stxlo(i)+stxhi(i)),xstate(i,j),
     &                           abs(uedge(i,j)).lt.eps)
         enddo
      endif
      place_to_break = 1
  100 continue
      endif

c     --------------------------------------------------------------
c     compute the yedge states
c     --------------------------------------------------------------
      
      if ((velpred.ne.1) .or. (n.eq.YVEL)) then
      do 200 i = imin, imax
      do j = jmin-1,jmax+1
         tr = half*
     &        (uad(i+1,j)+uad(i,j))*
     &        (xlo(i+1,j)-xlo(i,j))/hx

         stylo(j+1)= s(i,j) + (half-dthy*v(i,j))*sy(i,j)
     &        - dth*tr
     &        + dth*tforces(i,j)
         styhi(j  )= s(i,j) - (half+dthy*v(i,j))*sy(i,j)
     &        - dth*tr
     &        + dth*tforces(i,j)

      enddo

      if (bc(2,1).eq.EXT_DIR) then
         styhi(jmin) = s(i,jmin-1)
         stylo(jmin) = s(i,jmin-1)
      elseif (bc(2,1).eq.EXTRAP.or.bc(2,1).eq.HOEXTRAP
     &        .or.bc(2,1).eq.REFLECT_EVEN) then
         stylo(jmin) = styhi(jmin)
      elseif (bc(2,1).eq.REFLECT_ODD) then
         styhi(jmin) = zero
         stylo(jmin) = zero
      endif
      
      if (bc(2,2).eq.EXT_DIR) then
         stylo(jmax+1) = s(i,jmax+1)
         styhi(jmax+1) = s(i,jmax+1)
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
           ystate(i,j)=cvmgt(zero,sty,lty)
         enddo
      else
         do j=jmin,jmax+1
            ystate(i,j) = cvmgp(stylo(j),styhi(j),vedge(i,j))
            ystate(i,j) = cvmgt(half*(stylo(j)+styhi(j)),ystate(i,j),
     &                           abs(vedge(i,j)).lt.eps)
         enddo
      endif
      place_to_break = 1
  200 continue
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
     &     vol,DIMS(vol),
     &     lo,hi,iconserv )
c
c     This subroutine uses scalar edge states to compute
c     an advective tendency
c
      implicit none
      integer i,j
      integer iconserv
      REAL_T divux,divuy
      integer imin,jmin,imax,jmax
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(aofs)
      integer DIMDEC(vol)
      integer DIMDEC(uedge)
      integer DIMDEC(vedge)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      REAL_T aofs(DIMV(aofs))
      REAL_T vol(DIMV(vol))
      REAL_T uedge(DIMV(uedge))
      REAL_T vedge(DIMV(vedge))
      REAL_T xflux(DIMV(uedge))
      REAL_T yflux(DIMV(vedge))
      REAL_T areax(DIMV(ax))
      REAL_T areay(DIMV(ay))
c
      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)
c
c     if nonconservative initialize the advective tendency 
c     as -U*grad(S)
c
      if ( iconserv .ne. 1 ) then
         do j = jmin,jmax
            do i = imin,imax
               divux = (
     &              areax(i+1,j)*uedge(i+1,j) -
     &              areax(i,  j)*uedge(i,  j) )/vol(i,j)
               divuy = (
     &              areay(i,j+1)*vedge(i,j+1) -
     &              areay(i,j  )*vedge(i,j  ) )/vol(i,j)
               aofs(i,j) =
     &              - divux*half*(xflux(i+1,j) + xflux(i,j))
     &              - divuy*half*(yflux(i,j+1) + yflux(i,j))
            enddo
         enddo
      endif
c
c     convert edge states to fluxes
c
      do j = jmin,jmax
         do i = imin,imax+1
            xflux(i,j) = xflux(i,j)*uedge(i,j)*areax(i,j)
         enddo
      enddo
      do j = jmin,jmax+1
         do i = imin,imax
            yflux(i,j) = yflux(i,j)*vedge(i,j)*areay(i,j)
         enddo
      enddo
c     
c     compute the part of the advective tendency 
c     that depends on the flux convergence
c     
      if ( iconserv .ne. 1 ) then
         do j = jmin,jmax
            do i = imin,imax
               aofs(i,j) = aofs(i,j) + (
     &              xflux(i+1,j) - xflux(i,j) +
     &              yflux(i,j+1) - yflux(i,j))/vol(i,j)
            enddo
         enddo
      else
         do j = jmin,jmax
            do i = imin,imax
               aofs(i,j) = (
     &              xflux(i+1,j) - xflux(i,j) +
     &              yflux(i,j+1) - yflux(i,j))/vol(i,j)
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
     &     vol,DIMS(vol),
     &     lo,hi,iconserv )
c
c     This subroutine computes the sync advective tendency
c     for a state variable
c
      implicit none
      integer i,j
      integer iconserv
      integer imin,jmin,imax,jmax
      integer lo(SDIM),hi(SDIM)
      integer DIMDEC(sync)
      integer DIMDEC(vol)
      integer DIMDEC(ucor)
      integer DIMDEC(vcor)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      REAL_T sync(DIMV(sync))
      REAL_T vol(DIMV(vol))
      REAL_T ucor(DIMV(ucor))
      REAL_T vcor(DIMV(vcor))
      REAL_T xflux(DIMV(ucor))
      REAL_T yflux(DIMV(vcor))
      REAL_T areax(DIMV(ax))
      REAL_T areay(DIMV(ay))
c
      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)
c
c     compute corrective fluxes from edge states 
c     and perform conservative update
c
      do j = jmin,jmax
         do i = imin,imax+1
            xflux(i,j) = xflux(i,j)*ucor(i,j)*areax(i,j)
         enddo
      enddo
      do j = jmin,jmax+1
         do i = imin,imax
            yflux(i,j) = yflux(i,j)*vcor(i,j)*areay(i,j)
         enddo
      enddo
c     
      do j = jmin,jmax
         do i = imin,imax
            sync(i,j) = sync(i,j) + (
     &           xflux(i+1,j) - xflux(i,j) +
     &           yflux(i,j+1) - yflux(i,j) )/vol(i,j)
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
      integer DIMDEC(s), DIMDEC(xx)
      integer n
      integer lo(SDIM), hi(SDIM)
      REAL_T stx
      REAL_T s(DIMV(s))
      REAL_T xlo(DIMV(xx))
      REAL_T xhi(DIMV(xx))
      REAL_T uad(DIMV(xx))
      integer xbc(SDIM,2)
      logical ltest
      integer j
      integer imin,jmin,imax,jmax
c
      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)
c
      if (xbc(1,1).eq.EXT_DIR) then
         if (n .eq. XVEL) then
            do j = jmin-1,jmax+1
               xhi(imin,j) = s(imin-1,j)
               xlo(imin,j) = s(imin-1,j)
            enddo
         else
            do j = jmin-1,jmax+1
               ltest = abs(uad(imin,j)).le.1.0e-10
               stx   = cvmgt(xhi(imin,j),s(imin-1,j),ltest)
               xhi(imin,j) = stx
               xlo(imin,j) = stx
            enddo
         endif
      elseif (xbc(1,1).eq.EXTRAP.or.xbc(1,1).eq.HOEXTRAP
     &        .or.xbc(1,1).eq.REFLECT_EVEN) then
         do j = jmin-1,jmax+1
            xlo(imin,j) = xhi(imin,j)
         enddo
      elseif (xbc(1,1).eq.REFLECT_ODD) then
         do j = jmin-1,jmax+1
            xhi(imin,j) = zero
            xlo(imin,j) = zero
         enddo
      endif
c
      if (xbc(1,2).eq.EXT_DIR) then
         if (n .eq. XVEL) then
            do j = jmin-1,jmax+1
               xhi(imax+1,j) = s(imax+1,j)
               xlo(imax+1,j) = s(imax+1,j)
            enddo
         else
            do j = jmin-1,jmax+1
               ltest = abs(uad(imax+1,j)).le.1.0e-10
               stx   = cvmgt(xlo(imax+1,j),s(imax+1,j),ltest)
               xhi(imax+1,j) = stx
               xlo(imax+1,j) = stx
            enddo
         endif
      elseif (xbc(1,2).eq.EXTRAP.or.xbc(1,2).eq.HOEXTRAP
     &        .or.xbc(1,2).eq.REFLECT_EVEN) then
         do j = jmin-1,jmax+1
            xhi(imax+1,j) = xlo(imax+1,j)
         enddo
      elseif (xbc(1,2).eq.REFLECT_ODD) then
         do j = jmin-1,jmax+1
            xhi(imax+1,j) = zero
            xlo(imax+1,j) = zero
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
      integer DIMDEC(s), DIMDEC(yy)
      integer n
      integer lo(SDIM), hi(SDIM)
      REAL_T sty
      REAL_T s(DIMV(s))
      REAL_T ylo(DIMV(yy))
      REAL_T yhi(DIMV(yy))
      REAL_T vad(DIMV(yy))
      integer ybc(SDIM,2)
      logical ltest
      integer i
      integer imin,jmin,imax,jmax
c
      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)
c
      if (ybc(2,1).eq.EXT_DIR) then
         if (n .eq. YVEL) then
            do i = imin-1,imax+1
               yhi(i,jmin) = s(i,jmin-1)
               ylo(i,jmin) = s(i,jmin-1)
            enddo
         else
            do i = imin-1,imax+1
               ltest = abs(vad(i,jmin)).le.1.0e-10
               sty   = cvmgt(yhi(i,jmin),s(i,jmin-1),ltest)
               yhi(i,jmin) = sty
               ylo(i,jmin) = sty
            enddo
         endif
      elseif (ybc(2,1).eq.EXTRAP.or.ybc(2,1).eq.HOEXTRAP
     &        .or.ybc(2,1).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
            ylo(i,jmin) = yhi(i,jmin)
         enddo
      elseif (ybc(2,1).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
            yhi(i,jmin) = zero
            ylo(i,jmin) = zero
         enddo
      endif
c
      if (ybc(2,2).eq.EXT_DIR) then
         if (n .eq. YVEL) then
            do i = imin-1,imax+1
               ylo(i,jmax+1) = s(i,jmax+1)
               yhi(i,jmax+1) = s(i,jmax+1)
            enddo
         else
            do i = imin-1,imax+1
               ltest = abs(vad(i,jmax)).le.1.0e-10
               sty   = cvmgt(ylo(i,jmax+1),s(i,jmax+1),ltest)
               yhi(i,jmax+1) = sty
               ylo(i,jmax+1) = sty
            enddo
         endif
      elseif (ybc(2,2).eq.EXTRAP.or.ybc(2,2).eq.HOEXTRAP
     &        .or.ybc(2,2).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
            yhi(i,jmax+1) = ylo(i,jmax+1)
         enddo
      elseif (ybc(2,2).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
            ylo(i,jmax+1) = zero
            yhi(i,jmax+1) = zero
         enddo
      endif
c
c     end of routine
c
      return
      end


c=========================================================


      subroutine FORT_SLOPES (s,DIMS(s),slx,sly,DIMS(sl),
     &     lo,hi,slxscr,slyscr,bc)
c 
c     this subroutine computes first, second or forth order slopes of
c     a 2D scalar field.
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
      integer lo(SDIM), hi(SDIM)
      integer bc(SDIM,2)
      integer DIMDEC(s)
      integer DIMDEC(sl)
      REAL_T     s(DIMV(s))
      REAL_T   slx(DIMV(sl))
      REAL_T   sly(DIMV(sl))
      REAL_T slxscr(DIM1(s), 4)
      REAL_T slyscr(DIM2(s), 4)

      integer cen,lim,flag,fromm
      parameter( cen = 1 )
      parameter( lim = 2 )
      parameter( flag = 3 )
      parameter( fromm = 4 )

      integer imin,jmin,imax,jmax,i,j
      integer ng
      REAL_T dpls,dmin,ds
      REAL_T del,slim,sflg

      ng = lo(1) - ARG_L1(s)

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

c ::: ::::: added to prevent underflow for small s values
      do j = lo(2)-ng, hi(2)+ng
        do i = lo(1)-ng, hi(1)+ng
	  s(i,j) = cvmgt(s(i,j), zero, abs(s(i,j)).gt.1.0e-20)
        enddo
      enddo

c     -------------------------------------------------
c     COMPUTE 0TH order slopes
c     -------------------------------------------------
      if (slope_order.eq.1) then
        do j = jmin-1,jmax+1
           do i = imin-1,imax+1
              slx(i,j) = zero
              sly(i,j) = zero
  	   enddo
        enddo
      endif

c     -------------------------------------------------
c     COMPUTE 2nd order slopes
c     -------------------------------------------------
      if (slope_order.eq.2) then
        if (ng .lt. 2) then
           write(6,*) "SLOPE_2D: not enough bndry cells for 2nd order"
  	   stop
        endif
c     ------------------------ x slopes
        do j = jmin-1,jmax+1
          do i = imin-1,imax+1
            del  = half*(s(i+1,j)-s(i-1,j))
            dpls = two*(s(i+1,j) - s(i ,j))
            dmin = two*(s(i ,j) - s(i-1,j))
            slim = min(abs(dpls), abs(dmin))
            slim = cvmgp(slim, zero, dpls*dmin)
            sflg = sign(one,del)
            slx(i,j)= sflg*min(slim,abs(del))
  	enddo
        enddo
  
        if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
           do j = jmin-1, jmax+1
  	    slx(imin-1,j) = zero
              del  = (s(imin+1,j)+three*s(imin,j)-four*s(imin-1,j))/three
              dpls = two*(s(imin+1,j) - s(imin  ,j))
              dmin = two*(s(imin  ,j) - s(imin-1,j))
              slim = min(abs(dpls), abs(dmin))
              slim = cvmgp(slim, zero, dpls*dmin)
              sflg = sign(one,del)
              slx(imin,j)= sflg*min(slim,abs(del))
  	 enddo
        endif
        if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
           do j = jmin-1, jmax+1
  	    slx(imax+1,j) = zero
              del  = -(s(imax-1,j)+three*s(imax,j)-four*s(imax+1,j))/three
              dpls = two*(s(imax+1,j) - s(imax  ,j))
              dmin = two*(s(imax  ,j) - s(imax-1,j))
              slim = min(abs(dpls), abs(dmin))
              slim = cvmgp(slim, zero, dpls*dmin)
              sflg = sign(one,del)
              slx(imax,j)= sflg*min(slim,abs(del))
  	 enddo
        endif

c     ------------------------ y slopes
        do j = jmin-1,jmax+1
          do i = imin-1,imax+1
            del  = half*(s(i,j+1)-s(i,j-1))
            dpls = two*(s(i,j+1) - s(i,j ))
            dmin = two*(s(i,j ) - s(i,j-1))
            slim = min(abs(dpls),abs(dmin))
            slim = cvmgp(slim, zero, dpls*dmin)
            sflg = sign(one,del)
            sly(i,j)= sflg*min(slim,abs(del))
  	  enddo
        enddo

        if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
           do i = imin-1, imax+1
  	    sly(i,jmin-1) = zero
              del  = (s(i,jmin+1)+three*s(i,jmin)-four*s(i,jmin-1))/three
              dpls = two*(s(i,jmin+1) - s(i,jmin  ))
              dmin = two*(s(i,jmin  ) - s(i,jmin-1))
              slim = min(abs(dpls), abs(dmin))
              slim = cvmgp(slim, zero, dpls*dmin)
              sflg = sign(one,del)
              sly(i,jmin)= sflg*min(slim,abs(del))
           enddo
        endif
        if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
           do i = imin-1, imax+1
  	    sly(i,jmax+1) = zero
              del  = -(s(i,jmax-1)+three*s(i,jmax)-four*s(i,jmax+1))/three
              dpls = two*(s(i,jmax+1) - s(i,jmax  ))
              dmin = two*(s(i,jmax  ) - s(i,jmax-1))
              slim = min(abs(dpls), abs(dmin))
              slim = cvmgp(slim, zero, dpls*dmin)
              sflg = sign(one,del)
              sly(i,jmax)= sflg*min(slim,abs(del))
           enddo
        endif
c ... end, if slope_order .eq. 2
      endif

c     -------------------------------------------------
c     COMPUTE 4TH order slopes
c     -------------------------------------------------
      if (slope_order.eq.4) then
        if (ng .lt. 3) then
          write(6,*) "SLOPE_2D: not enough bndry cells for 4th order"
 	  stop
        endif
 
c     ------------------------ x slopes
        do j = jmin-1,jmax+1
          do i = imin-2,imax+2
             dmin           =  two*(s(i,  j)-s(i-1,j))
             dpls           =  two*(s(i+1,j)-s(i  ,j))
             slxscr(i,cen)  = half*(s(i+1,j)-s(i-1,j))
             slxscr(i,lim)  = min(abs(dmin),abs(dpls))
             slxscr(i,lim)  = cvmgp(slxscr(i,lim),zero,dpls*dmin)
             slxscr(i,flag) = sign(one,slxscr(i,cen))
             slxscr(i,fromm)= slxscr(i,flag)*
     &             min(slxscr(i,lim),abs(slxscr(i,cen)))
          enddo
          do i = imin-1,imax+1
             ds = two * two3rd * slxscr(i,cen) -
     &           sixth * (slxscr(i+1,fromm) + slxscr(i-1,fromm))
             slx(i,j) = slxscr(i,flag)*min(abs(ds),slxscr(i,lim))
          enddo
        enddo

        if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
           do j = jmin-1, jmax+1
              del  = -sixteen/fifteen*s(imin-1,j) + half*s(imin,j) + 
     &             two3rd*s(imin+1,j) - tenth*s(imin+2,j)
              dmin = two*(s(imin  ,j)-s(imin-1,j))
              dpls = two*(s(imin+1,j)-s(imin  ,j))
              slim = min(abs(dpls), abs(dmin))
              slim = cvmgp(slim, zero, dpls*dmin)
              sflg = sign(one,del)
	      slx(imin-1,j) = zero
              slx(imin,  j) = sflg*min(slim,abs(del))
           enddo
        endif
        if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
           do j = jmin-1, jmax+1
              del  = -( -sixteen/fifteen*s(imax+1,j) + half*s(imax,j) + 
     &                two3rd*s(imax-1,j) - tenth*s(imax-2,j) )
              dmin = two*(s(imax  ,j)-s(imax-1,j))
              dpls = two*(s(imax+1,j)-s(imax  ,j))
              slim = min(abs(dpls), abs(dmin))
              slim = cvmgp(slim, zero, dpls*dmin)
              sflg = sign(one,del)
              slx(imax,  j) = sflg*min(slim,abs(del))
	      slx(imax+1,j) = zero
           enddo
        endif

c     ------------------------ y slopes
        do i = imin-1,imax+1
           do j = jmin-2,jmax+2
              dmin           =  two*(s(i,j  )-s(i,j-1))
              dpls           =  two*(s(i,j+1)-s(i,j  ))
              slyscr(j,cen)  = half*(s(i,j+1)-s(i,j-1))
              slyscr(j,lim)  = min(abs(dmin),abs(dpls))
              slyscr(j,lim)  = cvmgp(slyscr(j,lim),zero,dpls*dmin)
              slyscr(j,flag) = sign(one,slyscr(j,cen))
              slyscr(j,fromm)= slyscr(j,flag)*
     &             min(slyscr(j,lim),abs(slyscr(j,cen)))
           enddo
           do j = jmin-1,jmax+1
              ds = two * two3rd * slyscr(j,cen) -
     &           sixth * (slyscr(j+1,fromm) + slyscr(j-1,fromm))
              sly(i,j) = slyscr(j,flag)*min(abs(ds),slyscr(j,lim))
           enddo
        enddo

        if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
           do i = imin-1, imax+1
              del  = -sixteen/fifteen*s(i,jmin-1) + half*s(i,jmin) + 
     &           two3rd*s(i,jmin+1) - tenth*s(i,jmin+2)
              dmin = two*(s(i,jmin  )-s(i,jmin-1))
              dpls = two*(s(i,jmin+1)-s(i,jmin  ))
              slim = min(abs(dpls), abs(dmin))
              slim = cvmgp(slim, zero, dpls*dmin)
              sflg = sign(one,del)
	      sly(i,jmin-1) = zero
              sly(i,jmin  ) = sflg*min(slim,abs(del))
	   enddo
        endif
        if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
           do i = imin-1, imax+1
              del  = -( -sixteen/fifteen*s(i,jmax+1) + half*s(i,jmax) + 
     &           two3rd*s(i,jmax-1) - tenth*s(i,jmax-2) )
              dmin = two*(s(i,jmax  )-s(i,jmax-1))
              dpls = two*(s(i,jmax+1)-s(i,jmax  ))
              slim = min(abs(dpls), abs(dmin))
              slim = cvmgp(slim, zero, dpls*dmin)
              sflg = sign(one,del)
              sly(i,jmax  ) = sflg*min(slim,abs(del))
	      sly(i,jmax+1) = zero
	   enddo
        endif
c ... end, if slope_order .eq. 4
      endif
c
c     end of routine
c
      return
      end




c=========================================================

      subroutine FORT_SCALMINMAX(s,DIMS(s),sn,DIMS(sn),
     &     lo,hi,bc)
c
c     correct an advected field for under/over shoots
c
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer lo(SDIM), hi(SDIM)
      integer bc(SDIM,2)
      REAL_T s(DIMV(s))
      REAL_T sn(DIMV(sn))
      integer  i, j, imin, imax, jmin, jmax
      REAL_T   smin, smax
c
      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
c
c     set corner values for s
c
c     :: corner ilo, jlo
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
         s(imin-1,jmin-1) = s(imin,jmin)
      endif
c     :: corner ihi, jlo
      if (bc(1,2).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
         s(imax+1,jmin-1) = s(imax,jmin)
      endif
c     :: corner ilo, jhi
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
         s(imin-1,jmax+1) = s(imin,jmax)
      endif
c     :: corner ihi, jhi
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
         s(imax+1,jmax+1) = s(imax,jmax)
      endif
c
c     compute extrema of s
c
      do j = jmin, jmax         
         do i = imin, imax
            smin = min(
     &           s(i-1,j-1),s(i  ,j-1),s(i+1,j-1),
     &           s(i-1,j  ),s(i  ,j  ),s(i+1,j  ),
     &           s(i-1,j+1),s(i  ,j+1),s(i+1,j+1))
            smax = max(
     &           s(i-1,j-1),s(i  ,j-1),s(i+1,j-1),
     &           s(i-1,j  ),s(i  ,j  ),s(i+1,j  ),
     &           s(i-1,j+1),s(i  ,j+1),s(i+1,j+1))
            sn(i,j) = max(sn(i,j),smin)
            sn(i,j) = min(sn(i,j),smax)
         enddo
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
      integer i, j, n
      integer DIMDEC(tf)
      integer DIMDEC(gp)
      integer DIMDEC(rho)
      integer lo(SDIM), hi(SDIM)
      REAL_T tforces(DIMV(tf),SDIM)
      REAL_T gp(DIMV(gp),SDIM)
      REAL_T rho(DIMV(rho))
c
      do n = 1, SDIM
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tforces(i,j,n) = (
     &              tforces(i,j,n)
     &              -    gp(i,j,n))/rho(i,j)
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
      integer i, j, n
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
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tforces(i,j,n) = (
     &              tforces(i,j,n)
     &              +  visc(i,j,n)
     &              -    gp(i,j,n))/rho(i,j)
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
c     sum divU*S into tforces or divide tforces by rho
c     depending on the value of iconserv
c
      integer nvar, iconserv
      integer lo(SDIM), hi(SDIM)
      integer i, j, n
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
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,n) = 
     &            tforces(i,j,n) - S(i,j,n)*divu(i,j)
               enddo
            enddo
         enddo

      else

         do n = 1, nvar
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,n) =
     &            tforces(i,j,n)/rho(i,j)
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
c     sum tforces, viscous foricing and divU*S into tforces
c     depending on the value of iconserv
c
      integer nvar, iconserv
      integer lo(SDIM), hi(SDIM)
      integer i, j, n
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
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,n) = 
     &                 tforces(i,j,n)
     &                 +  visc(i,j,n)
     &                 -     S(i,j,n)*divu(i,j)
               enddo
            enddo
         enddo

      else

         do n = 1, nvar
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,n) = (
     &                 tforces(i,j,n)
     &                 +  visc(i,j,n) )/rho(i,j)
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
      integer i, j, n, nvar
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
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               sn(i,j,n) = s(i,j,n) + dt*tforces(i,j,n)
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
      integer i, j, n, nvar
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
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               sp(i,j,n) = ss(i,j,n)
     &              + hdt*(tfs(i,j,n)-tfn(i,j,n))
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
      integer i, j, n, nvar
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
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               sn(i,j,n) = s(i,j,n)
     &              - dt*aofs(i,j,n)
     &              + dt*tforces(i,j,n)
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
      integer i, j, n
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
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               un(i,j,n) = u(i,j,n) 
     &              - dt*   aofs(i,j,n)
     &              + dt*tforces(i,j,n)/rho(i,j)
     &              - dt*     gp(i,j,n)/rho(i,j)
            enddo
         enddo
      enddo
c
c     end of routine
c
      return
      end


