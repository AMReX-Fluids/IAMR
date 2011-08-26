
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "GODUNOV_F.H"
#include "ArrayLim.H"

#define SDIM 2
#define XVEL 1
#define YVEL 2

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
      REAL_T    u, v
      REAL_T    small
      REAL_T    dt_start
      REAL_T    tforce1
      REAL_T    tforce2
      integer   i, j
      integer lo(SDIM), hi(SDIM)
      REAL_T  dt,dx(SDIM),cfl,u_max(SDIM)

      integer DIMDEC(vel)
      integer DIMDEC(rho)
      integer DIMDEC(tf)

      REAL_T  vel(DIMV(vel),SDIM)
      REAL_T  rho(DIMV(rho))
      REAL_T  tforces(DIMV(tf),SDIM)

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
      small   = 1.0e-8
#else
      small   = 1.0d-8
#endif
      u       = zero
      v       = zero
      tforce1 = zero
      tforce2 = zero

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            u = max(u,abs(vel(i,j,1)))
            v = max(v,abs(vel(i,j,2)))
             tforce1 = max(tforce1,abs(tforces(i,j,1)/rho(i,j)))
             tforce2 = max(tforce2,abs(tforces(i,j,2)/rho(i,j)))
         end do
      end do

      u_max(1) = u
      u_max(2) = v

      dt_start = 1.0D+20
      dt = dt_start

      if (u .gt. small) dt = min(dt,dx(1)/u)
      if (v .gt. small) dt = min(dt,dx(2)/v)

      if (tforce1 .gt. small) then
        dt  = min(dt,sqrt(two*dx(1)/tforce1))
      end if

      if (tforce2 .gt. small) then
        dt  = min(dt,sqrt(two*dx(2)/tforce2))
      end if

      if (dt .eq. dt_start) dt = min(dx(1),dx(2))

      dt = dt*cfl

      end

      subroutine FORT_TEST_U_RHO(
     &     u,DIMS(u),
     &     v,DIMS(v),
     &     rho,DIMS(rho),
     &     lo,hi,dt,dx,cflmax,u_max,verbose )
c
c     This subroutine computes the extrema of the density
c     and velocities at cell centers
c
      implicit none
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
      umax = -1.d200
      vmax = -1.d200
      umin =  1.d200
      vmin =  1.d200
      rhomax = -1.d200
      rhomin =  1.d200
c
      do j = jmin, jmax
         do i = imin, imax
            umax = max(umax,u(i,j))
            umin = min(umin,u(i,j))
            vmax = max(vmax,v(i,j))
            vmin = min(vmin,v(i,j))
            rhomax = max(rhomax,rho(i,j))
            rhomin = min(rhomin,rho(i,j))
         end do
      end do

      u_max(1) = max(abs(umax), abs(umin))
      u_max(2) = max(abs(vmax), abs(vmin))
      cflmax   = dt*max(u_max(1)/hx,u_max(2)/hy)

      if(verbose.eq.1)then
        write(6,1000) umax,umin,u_max(1) 
        write(6,1001) vmax,vmin,u_max(2) 
        write(6,1002) rhomax,rhomin
#ifndef	BL_NO_FORT_FLUSH
c        call flush(6)
#endif
      end if

 1000 format(' U  MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1001 format(' V  MAX/MIN/AMAX ',e21.14,2x,e21.14,2x,e21.14)
 1002 format('RHO MAX/MIN      ',e21.14,2x,e21.14)

      end

      subroutine FORT_TEST_UMAC_RHO(
     &     umac,DIMS(umac),
     &     vmac,DIMS(vmac),
     &     rho,DIMS(rho),
     &     lo,hi,dt,dx,cflmax,u_max)
c
c     This subroutine computes the extrema of the density
c     and mac edge velocities
c
      implicit none
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

      hx   = dx(1)
      hy   = dx(2)
      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
      umax = -1.d200
      vmax = -1.d200
      umin =  1.d200
      vmin =  1.d200
      rhomax = -1.d200
      rhomin =  1.d200

      do j = jmin, jmax
         do i = imin, imax+1
            umax = max(umax,umac(i,j))
            umin = min(umin,umac(i,j))
         end do
      end do
      do j = jmin, jmax+1
         do i = imin, imax
            vmax = max(vmax,vmac(i,j))
            vmin = min(vmin,vmac(i,j))
         end do
      end do
      do j = jmin, jmax
         do i = imin, imax
            rhomax = max(rhomax,rho(i,j))
            rhomin = min(rhomin,rho(i,j))
         end do
      end do

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
c      call flush(6)
#endif

      end

      subroutine FORT_TRANSVEL(
     &     u, ulo, uhi, sx, ubc, slxscr, Imx, Ipx, sedgex, DIMS(sedgex),
     &     v, vlo, vhi, sy, vbc, slyscr, Imy, Ipy, sedgey, DIMS(sedgey),
     &     DIMS(s), DIMS(work), DIMS(I),
     &     dsvl, DIMS(dsvl), sm, sp, DIMS(smp),
     &     lo,hi,dt,dx,use_minion,tforces,ppm_type)
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
      REAL_T eps,eps_for_bc
      logical ltm
      parameter( eps        = 1.0D-6 )
      parameter( eps_for_bc = 1.0D-10 )
      integer DIMDEC(s)
      integer DIMDEC(work)
      integer DIMDEC(I)
      integer DIMDEC(dsvl)
      integer DIMDEC(smp)
      integer DIMDEC(sedgex)
      integer DIMDEC(sedgey)
      REAL_T  u(DIMV(s))
      REAL_T  v(DIMV(s))
      REAL_T ulo(DIMV(work)),uhi(DIMV(work))
      REAL_T vlo(DIMV(work)),vhi(DIMV(work))
      REAL_T sx(DIMV(work))
      REAL_T sy(DIMV(work))
      REAL_T slxscr(DIM1(s), 4)
      REAL_T slyscr(DIM2(s), 4)
      REAL_T Imx(DIMV(I))
      REAL_T Ipx(DIMV(I))
      REAL_T Imy(DIMV(I))
      REAL_T Ipy(DIMV(I))
      REAL_T sedgex(DIMV(sedgex))
      REAL_T sedgey(DIMV(sedgey))
      REAL_T dsvl(DIMV(dsvl)), sm(DIMV(smp)), sp(DIMV(smp))

      integer use_minion, ppm_type
      REAL_T tforces(DIMV(work),SDIM)

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
      if (ppm_type .gt. 0) then
         call FORT_PPM(u,u,v,DIMS(s),
     &        Ipx,Imx,Ipy,Imy,DIMS(work),sm,sp,DIMS(smp),dsvl,DIMS(dsvl),
     &        sedgex,DIMS(sedgex),sedgey,DIMS(sedgey),lo,hi,dx,dt,ubc,
     &        eps_for_bc,ppm_type)

      else
         call FORT_SLOPES(u,DIMS(s),
     &        sx,sy,DIMS(work),lo,hi,slxscr,slyscr,ubc)
      end if
         
      if (ppm_type .gt. 0) then
         do i = imin, imax+1
            do j = jmin-1,jmax+1
               ulo(i,j) = Ipx(i-1,j)
               uhi(i,j) = Imx(i  ,j)
            end do
         end do
      else
         do i = imin, imax+1
            do j = jmin-1,jmax+1
               ulo(i,j) = u(i-1,j) + (half  - dthx*u(i-1,j))*sx(i-1,j)
               uhi(i,j) = u(i,j)   + (-half - dthx*u(i,  j))*sx(i,j)
            end do
         end do
      end if

      if(use_minion.eq.1)then
        do i = imin, imax+1
          do j = jmin-1,jmax+1
            ulo(i,j) = ulo(i,j) + dth*tforces(i-1,j,1)
            uhi(i,j) = uhi(i,j) + dth*tforces(i,  j,1)
          end do
        end do
      end if

      call trans_xbc(u,DIMS(s),
     &     ulo,uhi,DIMS(work),ulo,DIMS(work),lo,hi,XVEL,ubc,eps_for_bc)

      do i = imin,imax+1
         do j = jmin-1,jmax+1
            uad = merge(ulo(i,j),uhi(i,j),(ulo(i,j)+uhi(i,j)) .ge. 0.0d0)
            ltm = ulo(i,j) .le. zero  .and.  uhi(i,j) .ge. zero
            ltm = ltm .or. (abs(ulo(i,j)+uhi(i,j)) .lt. eps)
            ulo(i,j) = merge(zero,uad,ltm)
         end do
      end do
c
c     compute the y transverse velocities
c
      if (ppm_type .gt. 0) then
         call FORT_PPM(v,u,v,DIMS(s),
     &        Ipx,Imx,Ipy,Imy,DIMS(work),sm,sp,DIMS(smp),dsvl,DIMS(dsvl),
     &        sedgex,DIMS(sedgex),sedgey,DIMS(sedgey),lo,hi,dx,dt,vbc,
     &        eps_for_bc,ppm_type)
      else
         call FORT_SLOPES(v,DIMS(s),
     &        sx,sy,DIMS(work),lo,hi,slxscr,slyscr,vbc)
      end if

      if (ppm_type .gt. 0) then
         do i = imin-1,imax+1
            do j = jmin,jmax+1
               vlo(i,j) = Ipy(i,j-1)
               vhi(i,j) = Imy(i,j  )
            end do
         end do
      else
         do i = imin-1,imax+1
            do j = jmin,jmax+1
               vlo(i,j) = v(i,j-1) + (half  - dthy*v(i,j-1))*sy(i,j-1)
               vhi(i,j) = v(i,j)   + (-half - dthy*v(i,j  ))*sy(i,j)
            end do
         end do
      end if

      if(use_minion.eq.1)then
         do i = imin-1, imax+1
            do j = jmin,jmax+1
               vlo(i,j) = vlo(i,j) + dth*tforces(i,j-1,2)
               vhi(i,j) = vhi(i,j) + dth*tforces(i,j,  2)
            end do
         end do
      end if

      call trans_ybc(v,DIMS(s),
     &     vlo,vhi,DIMS(work),vlo,DIMS(work),lo,hi,YVEL,vbc,eps_for_bc)

      do i = imin-1,imax+1
         do j = jmin,jmax+1
            vad = merge(vlo(i,j),vhi(i,j),(vlo(i,j)+vhi(i,j)) .ge. 0.0d0)
            ltm = vlo(i,j) .le. zero  .and.  vhi(i,j) .ge. zero
            ltm = ltm .or. (abs(vlo(i,j)+vhi(i,j)) .lt. eps)
            vlo(i,j) = merge(zero,vad,ltm)
         end do
      end do

      end

      subroutine FORT_ESTATE(s, tforces, DIMS(s),

     &     u, xlo, xhi, sx, uad, slxscr, stxlo, stxhi,
     &     uedge, DIMS(uedge), xstate, DIMS(xstate), Imx, Ipx, sedgex, DIMS(sedgex),

     &     v, ylo, yhi, sy, vad, slyscr, stylo, styhi,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate), Imy, Ipy, sedgey, DIMS(sedgey),

     &     DIMS(work), DIMS(I), dsvl, DIMS(dsvl), sm, sp, DIMS(smp),
     &     bc,lo,hi,dt,dx,n,velpred, use_minion, ppm_type)
c
      implicit none
      integer i,j,n,velpred
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer imin,jmin,imax,jmax
      REAL_T place_to_break
      REAL_T hx, hy, dt, dth, dthx, dthy, dx(SDIM)
      REAL_T tr,stx,sty,fu,fv
      REAL_T eps,eps_for_bc
      logical ltx,lty
      parameter( eps        = 1.0D-6 )
      parameter( eps_for_bc = 1.0D-10 )

      integer DIMDEC(s)
      integer DIMDEC(work)
      integer DIMDEC(uedge)
      integer DIMDEC(vedge)
      integer DIMDEC(xstate)
      integer DIMDEC(ystate)
      integer DIMDEC(I)
      integer DIMDEC(sedgex)
      integer DIMDEC(sedgey)
      integer DIMDEC(dsvl)
      integer DIMDEC(smp)

      REAL_T s(DIMV(s))
      REAL_T u(DIMV(s))
      REAL_T v(DIMV(s))
      REAL_T stxlo(DIM1(s)),stxhi(DIM1(s)),slxscr(DIM1(s),4)
      REAL_T stylo(DIM2(s)),styhi(DIM2(s)),slyscr(DIM2(s),4)

      REAL_T uedge(DIMV(uedge)), xstate(DIMV(xstate)), Imx(DIMV(I)), Ipx(DIMV(I)), sedgex(DIMV(sedgex))
      REAL_T vedge(DIMV(vedge)), ystate(DIMV(ystate)), Imy(DIMV(I)), Ipy(DIMV(I)), sedgey(DIMV(sedgey))
      REAL_T dsvl(DIMV(dsvl)), sm(DIMV(smp)), sp(DIMV(smp))

      REAL_T xlo(DIMV(work)), xhi(DIMV(work))
      REAL_T ylo(DIMV(work)), yhi(DIMV(work))
      REAL_T  sx(DIMV(work)), uad(DIMV(work)) 
      REAL_T  sy(DIMV(work)), vad(DIMV(work)) 
      REAL_T tforces(DIMV(work))

      integer use_minion, ppm_type

      dth  = half*dt
      dthx = half*dt/dx(1)
      dthy = half*dt/dx(2)
      hx = dx(1)
      hy = dx(2)
      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

c
c     compute the slopes
c
c     trace the state to the cell edges
c
      if (ppm_type .gt. 0) then
         call FORT_PPM(s,u,v,DIMS(s),
     &        Ipx,Imx,Ipy,Imy,DIMS(work),sm,sp,DIMS(smp),dsvl,DIMS(dsvl),
     &        sedgex,DIMS(sedgex),sedgey,DIMS(sedgey),lo,hi,dx,dt,bc,
     &        eps_for_bc,ppm_type)
      else
         call FORT_SLOPES(s,DIMS(s),
     &        sx,sy,DIMS(work),lo,hi,slxscr,slyscr,bc)
      end if
c
c     trace the state to the cell edges
c
      if (ppm_type .gt. 0) then
         do i = imin, imax+1
            do j = jmin-1,jmax+1
               xlo(i,j) = Ipx(i-1,j)
               xhi(i,j) = Imx(i  ,j)
            end do
         end do
      else
         do i = imin, imax+1
            do j = jmin-1,jmax+1
               xlo(i,j) = s(i-1,j) + (half - dthx*u(i-1,j))*sx(i-1,j)
               xhi(i,j) = s(i  ,j) - (half + dthx*u(i  ,j))*sx(i  ,j)
            end do
         end do
      end if

      if(use_minion.eq.1)then
         do i = imin, imax+1
            do j = jmin-1,jmax+1
               xlo(i,j) = xlo(i,j) + dth*tforces(i-1,j)
               xhi(i,j) = xhi(i,j) + dth*tforces(i,  j)
            end do
         end do
      end if

      call trans_xbc(
     &     s,DIMS(s),
     &     xlo,xhi,DIMS(work),uad,DIMS(work),
     &     lo,hi,n,bc,eps_for_bc)
      
      do j = jmin-1,jmax+1
         do i = imin, imax+1
            fu  = merge(zero,one,abs(uad(i,j)).lt.eps)
            stx = merge(xlo(i,j),xhi(i,j),uad(i,j) .ge. 0.0d0)
            xlo(i,j) = fu*stx + (one - fu)*half*(xhi(i,j)+xlo(i,j))
         end do
      end do

      if (ppm_type .gt. 0) then
         do j = jmin,jmax+1
            do i = imin-1,imax+1
               ylo(i,j) = Ipy(i,j-1)
               yhi(i,j) = Ipy(i,j  )
            end do
         end do
      else
         do j = jmin,jmax+1
            do i = imin-1,imax+1
               ylo(i,j) = s(i,j-1) + (half - dthy*v(i,j-1))*sy(i,j-1)
               yhi(i,j) = s(i,j  ) - (half + dthy*v(i,j  ))*sy(i,j  )
            end do
         end do
      end if

      if(use_minion.eq.1)then
         do i = imin-1, imax+1
            do j = jmin,jmax+1
               ylo(i,j) = ylo(i,j) + dth*tforces(i,j-1)
               yhi(i,j) = yhi(i,j) + dth*tforces(i,j)
            end do
         end do
      end if

      call trans_ybc(
     &     s,DIMS(s),
     &     ylo,yhi,DIMS(work),vad,DIMS(work),
     &     lo,hi,n,bc,eps_for_bc)

      do j = jmin,jmax+1
         do i = imin-1,imax+1
            fv  = merge(zero,one,abs(vad(i,j)).lt.eps)
            sty = merge(ylo(i,j),yhi(i,j),vad(i,j) .ge. 0.0d0)
            ylo(i,j) = fv*sty + (one - fv)*half*(yhi(i,j)+ylo(i,j))
         end do
      end do
c
c     compute the xedge state
c
      if ((velpred.ne.1) .or. (n.eq.XVEL)) then
         do 100 j = jmin,jmax
            do i = imin-1,imax+1
               tr = half*
     &              (vad(i,j+1)+vad(i,j))*
     &              (ylo(i,j+1)-ylo(i,j))/hy
               
               if (ppm_type .gt. 0) then
                  stxlo(i+1) = Ipx(i,j)
     &                 - dth*tr
     &                 + dth*tforces(i,j)
                  stxhi(i  ) = Imx(i,j)
     &                 - dth*tr
     &                 + dth*tforces(i,j)
               else
                  stxlo(i+1) = s(i,j) + (half-dthx*u(i,j))*sx(i,j)
     &                 - dth*tr
     &                 + dth*tforces(i,j)
                  stxhi(i  ) = s(i,j) - (half+dthx*u(i,j))*sx(i,j)
     &                 - dth*tr
     &                 + dth*tforces(i,j)
               end if
            end do

            if (bc(1,1).eq.EXT_DIR .and. velpred.eq.1) then
               stxhi(imin) = s(imin-1,j)
               stxlo(imin) = s(imin-1,j)
            else if (bc(1,1).eq.EXT_DIR .and. uad(imin,j).ge.zero) then
               stxhi(imin) = s(imin-1,j)
               stxlo(imin) = s(imin-1,j)
            else if (bc(1,1).eq.EXT_DIR .and. uad(imin,j).lt.zero) then
               stxlo(imin) = stxhi(imin)
            else if (bc(1,1).eq.FOEXTRAP.or.bc(1,1).eq.HOEXTRAP
     &              .or.bc(1,1).eq.REFLECT_EVEN) then
               stxlo(imin) = stxhi(imin)
            else if (bc(1,1).eq.REFLECT_ODD) then
               stxhi(imin) = zero
               stxlo(imin) = stxhi(imin)
            end if

            if (bc(1,2).eq.EXT_DIR .and. velpred.eq.1) then
               stxlo(imax+1) = s(imax+1,j)
               stxhi(imax+1) = s(imax+1,j)
            else if (bc(1,2).eq.EXT_DIR .and. uad(imax+1,j).le.zero) then
               stxlo(imax+1) = s(imax+1,j)
               stxhi(imax+1) = s(imax+1,j)
            else if (bc(1,2).eq.EXT_DIR .and. uad(imax+1,j).gt.zero) then
               stxhi(imax+1) = stxlo(imax+1) 
            else if (bc(1,2).eq.FOEXTRAP.or.bc(1,2).eq.HOEXTRAP
     &              .or.bc(1,2).eq.REFLECT_EVEN) then
               stxhi(imax+1) = stxlo(imax+1)
            else if (bc(1,2).eq.REFLECT_ODD) then
               stxlo(imax+1) = zero
               stxhi(imax+1) = zero
            end if

            if ( velpred .eq. 1 ) then
               do i = imin, imax+1
                  ltx = stxlo(i) .le. zero  .and.  stxhi(i) .ge. zero
                  ltx = ltx .or. (abs(stxlo(i)+stxhi(i)) .lt. eps)
                  stx = merge(stxlo(i),stxhi(i),(stxlo(i)+stxhi(i)) .ge. 0.0d0)
                  xstate(i,j) = merge(zero,stx,ltx)
               end do
            else
               do i = imin, imax+1
                  xstate(i,j) = merge(stxlo(i),stxhi(i),uedge(i,j) .ge. 0.0d0)
                  xstate(i,j) = merge(half*(stxlo(i)+stxhi(i)),xstate(i,j),
     &                 abs(uedge(i,j)).lt.eps)
               end do
            end if
            place_to_break = 1
100      continue
      end if
c
c     compute the yedge states
c      
      if ((velpred.ne.1) .or. (n.eq.YVEL)) then
         do 200 i = imin, imax
            do j = jmin-1,jmax+1
               tr = half*
     &              (uad(i+1,j)+uad(i,j))*
     &              (xlo(i+1,j)-xlo(i,j))/hx

               if (ppm_type .gt. 0) then
                  stylo(j+1)= Ipy(i,j)
     &                 - dth*tr
     &                 + dth*tforces(i,j)
                  styhi(j  )= Imy(i,j)
     &                 - dth*tr
     &                 + dth*tforces(i,j)
               else
                  stylo(j+1)= s(i,j) + (half-dthy*v(i,j))*sy(i,j)
     &                 - dth*tr
     &                 + dth*tforces(i,j)
                  styhi(j  )= s(i,j) - (half+dthy*v(i,j))*sy(i,j)
     &                 - dth*tr
     &                 + dth*tforces(i,j)
               end if
            end do

            if (bc(2,1).eq.EXT_DIR .and. velpred .eq. 1) then
               styhi(jmin) = s(i,jmin-1)
               stylo(jmin) = s(i,jmin-1)
            else if (bc(2,1).eq.EXT_DIR .and. vad(i,jmin).ge.zero) then
               styhi(jmin) = s(i,jmin-1)
               stylo(jmin) = s(i,jmin-1)
            else if (bc(2,1).eq.EXT_DIR .and. vad(i,jmin).lt.zero) then
               stylo(jmin) = styhi(jmin)
            else if (bc(2,1).eq.FOEXTRAP.or.bc(2,1).eq.HOEXTRAP
     &              .or.bc(2,1).eq.REFLECT_EVEN) then
               stylo(jmin) = styhi(jmin)
            else if (bc(2,1).eq.REFLECT_ODD) then
               styhi(jmin) = zero
               stylo(jmin) = zero
            end if
            
            if (bc(2,2).eq.EXT_DIR .and. velpred .eq. 1) then
               stylo(jmax+1) = s(i,jmax+1)
               styhi(jmax+1) = s(i,jmax+1)
            else if (bc(2,2).eq.EXT_DIR .and. vad(i,jmax+1).le.zero) then
               stylo(jmax+1) = s(i,jmax+1)
               styhi(jmax+1) = s(i,jmax+1)
            else if (bc(2,2).eq.EXT_DIR .and. vad(i,jmax+1).gt.zero) then
               styhi(jmax+1) = stylo(jmax+1)
            else if (bc(2,2).eq.FOEXTRAP.or.bc(2,2).eq.HOEXTRAP
     &              .or.bc(2,2).eq.REFLECT_EVEN) then
               styhi(jmax+1) = stylo(jmax+1)
            else if (bc(2,2).eq.REFLECT_ODD) then
               stylo(jmax+1) = zero
               styhi(jmax+1) = zero
            end if

            if ( velpred .eq. 1 ) then
               do j = jmin, jmax+1
                  lty = stylo(j) .le. zero  .and.  styhi(j) .ge. zero
                  lty = lty .or. (abs(stylo(j)+styhi(j)) .lt. eps)
                  sty = merge(stylo(j),styhi(j),(stylo(j)+styhi(j)) .ge. 0.0d0)
                  ystate(i,j)=merge(zero,sty,lty)
               end do
            else
               do j=jmin,jmax+1
                  ystate(i,j) = merge(stylo(j),styhi(j),vedge(i,j) .ge. 0.0d0)
                  ystate(i,j) = merge(half*(stylo(j)+styhi(j)),ystate(i,j),
     &                 abs(vedge(i,j)).lt.eps)
               end do
            end if
            place_to_break = 1
200      continue
      end if
      end

      subroutine FORT_ESTATE_FPU(s, tforces, divu, DIMS(s),
     &     xlo, xhi, sx, slxscr, stxlo, stxhi,
     &     uedge, DIMS(uedge), xstate, DIMS(xstate), Imx, Ipx, sedgex, DIMS(sedgex),

     &     ylo, yhi, sy, slyscr, stylo, styhi,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate), Imy, Ipy, sedgey, DIMS(sedgey),

     &     DIMS(work), DIMS(I), dsvl, DIMS(dsvl), sm, sp, DIMS(smp),
     &     bc,lo,hi,dt,dx,n,use_minion,iconserv,ppm_type)
c
c     This subroutine computes edges states, right now it uses
c     a lot of memory, but there becomes a trade off between
c     simplicity-efficiency in the new way of computing states
c     and complexity in the old way.  By eliminating loops over
c     state components though, the new way uses much less memory.
c
      implicit none
      integer i,j,n
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer imin,jmin,imax,jmax
      REAL_T place_to_break
      REAL_T hx, hy, dt, dth, dthx, dthy
      REAL_T tr,ubar,vbar,stx,sty,fu,fv,dx(SDIM)
      REAL_T eps,eps_for_bc
      logical ltx,lty
      parameter( eps        = 1.0D-6 )
      parameter( eps_for_bc = 1.0D-10 )

      integer DIMDEC(s)
      integer DIMDEC(work)
      integer DIMDEC(uedge)
      integer DIMDEC(xstate)
      integer DIMDEC(vedge)
      integer DIMDEC(ystate)
      integer DIMDEC(I)
      integer DIMDEC(sedgex)
      integer DIMDEC(sedgey)
      integer DIMDEC(dsvl)
      integer DIMDEC(smp)

      REAL_T s(DIMV(s))
      REAL_T stxlo(DIM1(s)),stxhi(DIM1(s)),slxscr(DIM1(s),4)
      REAL_T stylo(DIM2(s)),styhi(DIM2(s)),slyscr(DIM2(s),4)

      REAL_T uedge(DIMV(uedge)), xstate(DIMV(xstate)), Imx(DIMV(I)), Ipx(DIMV(I)), sedgex(DIMV(sedgex))
      REAL_T vedge(DIMV(vedge)), ystate(DIMV(ystate)), Imy(DIMV(I)), Ipy(DIMV(I)), sedgey(DIMV(sedgey))
      REAL_T dsvl(DIMV(dsvl)), sm(DIMV(smp)), sp(DIMV(smp))

      REAL_T xlo(DIMV(work)), xhi(DIMV(work))
      REAL_T ylo(DIMV(work)), yhi(DIMV(work))
      REAL_T  sx(DIMV(work))
      REAL_T  sy(DIMV(work))
      REAL_T tforces(DIMV(work))
      REAL_T divu(DIMV(work))

#if 0
      REAL_T extra_work(DIMV(work),6)
#endif
      
      REAL_T st,denom
      integer ppm_type

      integer use_minion, iconserv, inc

      REAL_T umax

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
c     compute the slopes
c
      if (ppm_type .gt. 0) then
         call FORT_PPM_FPU(s,DIMS(s),uedge,vedge,DIMS(uedge),
     &        Ipx,Imx,Ipy,Imy,DIMS(work),sm,sp,DIMS(smp),dsvl,DIMS(dsvl),
     &        sedgex,DIMS(sedgex),sedgey,DIMS(sedgey),lo,hi,dx,dt,bc, 
     &        eps_for_bc,ppm_type)

      else
         call FORT_SLOPES(
     &        s,DIMS(s),
     &        sx,sy,DIMS(work),
     &        lo,hi,slxscr,slyscr,bc)
      end if

c
c     trace the state to the cell edges
c

      if (ppm_type .gt. 0) then
         do j = jmin-1,jmax+1
            do i = imin,  imax+1
               xlo(i,j) = Ipx(i-1,j)
               xhi(i,j) = Imx(i  ,j)
            end do
         end do
      else
         do j = jmin-1,jmax+1
            do i = imin,  imax+1
               xlo(i,j) = s(i-1,j) + (half - dthx*uedge(i,j))*sx(i-1,j)
               xhi(i,j) = s(i,  j) - (half + dthx*uedge(i,j))*sx(i,  j)
            end do
         end do
      end if

      if(use_minion.eq.1)then
         do j = jmin-1,jmax+1
            do i = imin,  imax+1
               xlo(i,j) = xlo(i,j) + dth*tforces(i-1,j)
               xhi(i,j) = xhi(i,j) + dth*tforces(i,  j)
            end do
         end do
         if (iconserv .eq. 1) then
           do j = jmin-1,jmax+1
            do i = imin,  imax+1
               xlo(i,j) = xlo(i,j) - dth*s(i-1,j)*divu(i-1,j)
               xhi(i,j) = xhi(i,j) - dth*s(i  ,j)*divu(i,  j)
            end do
           end do
         end if
      end if

      call trans_xbc(
     &     s,DIMS(s),
     &     xlo,xhi,DIMS(work),uedge,DIMS(uedge),
     &     lo,hi,n,bc,eps_for_bc)

      do j = jmin-1,jmax+1
         do i = imin,  imax+1
            fu  = merge(zero,one,abs(uedge(i,j)).lt.eps)
            stx = merge(xlo(i,j),xhi(i,j),uedge(i,j) .ge. 0.0d0)
            xlo(i,j) = fu*stx + (one - fu)*half*(xhi(i,j)+xlo(i,j))
         end do
      end do

      if (ppm_type .gt. 0) then
         do j = jmin,  jmax+1
            do i = imin-1,imax+1
               ylo(i,j) = Ipy(i,j-1)
               yhi(i,j) = Imy(i,j  )
            end do
         end do
      else
         do j = jmin,  jmax+1
            do i = imin-1,imax+1
               ylo(i,j) = s(i,j-1) + (half - dthy*vedge(i,j))*sy(i,j-1)
               yhi(i,j) = s(i,j)   - (half + dthy*vedge(i,j))*sy(i,j)
            end do
         end do
      end if

      if (use_minion.eq.1)then
         do j = jmin, jmax+1
            do i = imin-1,  imax+1
               ylo(i,j) = ylo(i,j) + dth*tforces(i,j-1)
               yhi(i,j) = yhi(i,j) + dth*tforces(i,j  )
            end do
         end do
         if (iconserv .eq. 1) then
           do j = jmin-1,jmax+1
            do i = imin,  imax+1
               ylo(i,j) = ylo(i,j) - dth*s(i,j-1)*divu(i,j-1)
               yhi(i,j) = yhi(i,j) - dth*s(i,j  )*divu(i,j  )
            end do
           end do
         end if
      end if

      call trans_ybc(
     &     s,DIMS(s),
     &     ylo,yhi,DIMS(work),vedge,DIMS(vedge),
     &     lo,hi,n,bc,eps_for_bc)

      do j = jmin,  jmax+1
         do i = imin-1,imax+1
            fv  = merge(zero,one,abs(vedge(i,j)).lt.eps)
            sty = merge(ylo(i,j),yhi(i,j),vedge(i,j) .ge. 0.0d0)
            ylo(i,j) = fv*sty + (one - fv)*half*(yhi(i,j)+ylo(i,j))
         end do
      end do

#if 0
      call bdsslope(s,lo(1),lo(2),hi(1),hi(2),
     $              extra_work(lo(1)-1,lo(2)-1,1),
     $              extra_work(lo(1)-1,lo(2)-1,2),
     $              extra_work(lo(1)-1,lo(2)-1,3),
     $              dx)
#endif
c
c     compute the xedge states
c
      do j = jmin,jmax

            do i = imin-1,imax+1

               if (iconserv.eq.1) then

                  st = -(vedge(i,j+1)*ylo(i,j+1) - vedge(i,j)*ylo(i,j))/hy
     &                 + s(i,j)*(vedge(i,j+1)-vedge(i,j))/hy
     &                 - s(i,j)*divu(i,j)

               else

#if 1
                  if (vedge(i,j)*vedge(i,j+1).le.0.d0) then
                     vbar = 0.5d0*(vedge(i,j)+vedge(i,j+1))
                     if (vbar.lt.0.d0) then
                        inc = 1
                     else
                        inc = 0
                     endif
                     tr = vbar*(s(i,j+inc)-s(i,j+inc-1))/hy
                  else
                     tr = half*(vedge(i,j+1) + vedge(i,j)) *
     &                    (  ylo(i,j+1) - ylo(i,j)  ) / hy
                  endif

#else
                  tr = half*(vedge(i,j+1) + vedge(i,j)) * extra_work(i,j,2)
#endif
                  st =  -tr

               endif

               if (ppm_type .gt. 0) then
                  stxlo(i+1)= Ipx(i,j)
     &                 + dth*(st + tforces(i,j))
                  stxhi(i  )= Imx(i,j)
     &                 + dth*(st + tforces(i,j))
               else
                  stxlo(i+1)= s(i,j) + (half-dthx*uedge(i+1,j))*sx(i,j)
     &                 + dth*(st + tforces(i,j))
                  stxhi(i  )= s(i,j) - (half+dthx*uedge(i  ,j))*sx(i,j)
     &                 + dth*(st + tforces(i,j))
               end if

            end do

            if (bc(1,1).eq.EXT_DIR .and. uedge(imin,j).ge.zero) then
               stxhi(imin) = s(imin-1,j)
               stxlo(imin) = s(imin-1,j)
            else if (bc(1,1).eq.EXT_DIR .and. uedge(imin,j).lt.zero) then
               stxlo(imin) = stxhi(imin)
            else if (bc(1,1).eq.FOEXTRAP.or.bc(1,1).eq.HOEXTRAP
     &              .or.bc(1,1).eq.REFLECT_EVEN) then
               stxlo(imin) = stxhi(imin)
            else if (bc(1,1).eq.REFLECT_ODD) then
               stxhi(imin) = zero
               stxlo(imin) = zero
            end if
            if (bc(1,2).eq.EXT_DIR .and. uedge(imax+1,j).le.zero) then
               stxlo(imax+1) = s(imax+1,j)
               stxhi(imax+1) = s(imax+1,j)
            else if (bc(1,2).eq.EXT_DIR .and. uedge(imax+1,j).gt.zero) then
               stxhi(imax+1) = stxlo(imax+1)
            else if (bc(1,2).eq.FOEXTRAP.or.bc(1,2).eq.HOEXTRAP
     &              .or.bc(1,2).eq.REFLECT_EVEN) then
               stxhi(imax+1) = stxlo(imax+1)
            else if (bc(1,2).eq.REFLECT_ODD) then
               stxlo(imax+1) = zero
               stxhi(imax+1) = zero
            end if
            
            do i = imin, imax+1
               xstate(i,j) = merge(stxlo(i),stxhi(i),uedge(i,j) .ge. 0.0d0)
               xstate(i,j) = merge(half*(stxlo(i)+stxhi(i)),xstate(i,j)
     &              ,abs(uedge(i,j)).lt.eps)
            end do
            place_to_break = 1
      end do
c
c     compute the yedge states
c
      do i = imin,imax
            
            do j = jmin-1,jmax+1

               if (iconserv.eq.1) then

                  st = -(uedge(i+1,j)*xlo(i+1,j) - uedge(i,j)*xlo(i,j))/hx
     &                 + s(i,j)*(uedge(i+1,j)-uedge(i,j))/hx
     &                 - s(i,j)*divu(i,j)

               else

#if 1
                  if (uedge(i,j)*uedge(i+1,j).le.0.d0) then
                     ubar = 0.5d0*(uedge(i,j)+uedge(i+1,j))
                     if (ubar.lt.0.d0) then
                        inc = 1
                     else
                        inc = 0
                     endif
                     tr = ubar*(s(i+inc,j)-s(i+inc-1,j))/hx
                  else
                     tr = half*(uedge(i+1,j) + uedge(i,j)) *
     &                           (xlo(i+1,j) -   xlo(i,j)  ) / hx
                  endif
#else
                  tr = half*(uedge(i+1,j) + uedge(i,j)) * extra_work(i,j,1)
#endif
                  st = -tr

               endif

               if (ppm_type .gt. 0) then
                  stylo(j+1)= Ipy(i,j)
     &                 + dth*(st + tforces(i,j))
                  styhi(j  )= Imy(i,j)
     &                 + dth*(st + tforces(i,j))
               else
                  stylo(j+1)= s(i,j) + (half-dthy*vedge(i,j+1))*sy(i,j)
     &                 + dth*(st + tforces(i,j))
                  styhi(j  )= s(i,j) - (half+dthy*vedge(i,j  ))*sy(i,j)
     &                 + dth*(st + tforces(i,j))
               end if

            end do
            
            if (bc(2,1).eq.EXT_DIR .and. vedge(i,jmin).ge.zero) then
               styhi(jmin) = s(i,jmin-1)
               stylo(jmin) = s(i,jmin-1)
            else if (bc(2,1).eq.EXT_DIR .and. vedge(i,jmin).lt.zero) then
               stylo(jmin) = styhi(jmin)
            else if (bc(2,1).eq.FOEXTRAP.or.bc(2,1).eq.HOEXTRAP
     &              .or.bc(2,1).eq.REFLECT_EVEN) then
               stylo(jmin) = styhi(jmin)
            else if (bc(2,1).eq.REFLECT_ODD) then
               styhi(jmin) = zero
               stylo(jmin) = zero
            end if
            
            if (bc(2,2).eq.EXT_DIR .and. vedge(i,jmax+1).le.zero) then
               stylo(jmax+1) = s(i,jmax+1)
               styhi(jmax+1) = s(i,jmax+1)
            else if (bc(2,2).eq.EXT_DIR .and. vedge(i,jmax+1).gt.zero) then
               styhi(jmax+1) = stylo(jmax+1)
            else if (bc(2,2).eq.FOEXTRAP.or.bc(2,2).eq.HOEXTRAP
     &              .or.bc(2,2).eq.REFLECT_EVEN) then
               styhi(jmax+1) = stylo(jmax+1)
            else if (bc(2,2).eq.REFLECT_ODD) then
               stylo(jmax+1) = zero
               styhi(jmax+1) = zero
            end if
            
            do j=jmin,jmax+1
               ystate(i,j) = merge(stylo(j),styhi(j),vedge(i,j) .ge. 0.0d0)
               ystate(i,j) = merge(half*(stylo(j)+styhi(j)),ystate(i,j),
     &              abs(vedge(i,j)).lt.eps)
            end do
            place_to_break = 1
      end do

      end


      subroutine FORT_ADV_FORCING(
     &     aofs,DIMS(aofs),
     &     xflux,DIMS(xflux),
     &     uedge,DIMS(uedge),
     &     areax,DIMS(ax),
     &     yflux,DIMS(yflux),
     &     vedge,DIMS(vedge),
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
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      REAL_T aofs(DIMV(aofs))
      REAL_T vol(DIMV(vol))
      REAL_T uedge(DIMV(uedge))
      REAL_T vedge(DIMV(vedge))
      REAL_T xflux(DIMV(xflux))
      REAL_T yflux(DIMV(yflux))
      REAL_T areax(DIMV(ax))
      REAL_T areay(DIMV(ay))

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)
c
c     if nonconservative initialize the advective tendency as -U*grad(S)
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
            end do
         end do
      end if
c
c     convert edge states to fluxes
c
      do j = jmin,jmax
         do i = imin,imax+1
            xflux(i,j) = xflux(i,j)*uedge(i,j)*areax(i,j)
         end do
      end do
      do j = jmin,jmax+1
         do i = imin,imax
            yflux(i,j) = yflux(i,j)*vedge(i,j)*areay(i,j)
         end do
      end do
c     
c     compute part of advective tendency that depends on the flux convergence
c     
      if ( iconserv .ne. 1 ) then
         do j = jmin,jmax
            do i = imin,imax
               aofs(i,j) = aofs(i,j) + (
     &              xflux(i+1,j) - xflux(i,j) +
     &              yflux(i,j+1) - yflux(i,j))/vol(i,j)
            end do
         end do
      else
         do j = jmin,jmax
            do i = imin,imax
               aofs(i,j) = (
     &              xflux(i+1,j) - xflux(i,j) +
     &              yflux(i,j+1) - yflux(i,j))/vol(i,j)
            end do
         end do
      end if

      end

      subroutine FORT_SYNC_ADV_FORCING(
     &     sync,DIMS(sync),
     &     xflux,DIMS(xflux),
     &     ucor,DIMS(ucor),
     &     areax,DIMS(ax),
     &     yflux,DIMS(yflux),
     &     vcor,DIMS(vcor),
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
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      REAL_T sync(DIMV(sync))
      REAL_T vol(DIMV(vol))
      REAL_T ucor(DIMV(ucor))
      REAL_T vcor(DIMV(vcor))
      REAL_T xflux(DIMV(xflux))
      REAL_T yflux(DIMV(yflux))
      REAL_T areax(DIMV(ax))
      REAL_T areay(DIMV(ay))

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
         end do
      end do
      do j = jmin,jmax+1
         do i = imin,imax
            yflux(i,j) = yflux(i,j)*vcor(i,j)*areay(i,j)
         end do
      end do

      do j = jmin,jmax
         do i = imin,imax
            sync(i,j) = sync(i,j) + (
     &           xflux(i+1,j) - xflux(i,j) +
     &           yflux(i,j+1) - yflux(i,j) )/vol(i,j)
         end do
      end do

      end

      subroutine trans_xbc(
     &     s,DIMS(s),
     &     xlo,xhi,DIMS(xx),uad,DIMS(uad),
     &     lo,hi,n,xbc,eps)
c
c     This subroutine processes boundary conditions on information
c     traced to cell faces in the x direction.  This is used for
c     computing velocities and edge states used in calculating
c     transverse derivatives
c
      implicit none
      integer DIMDEC(s), DIMDEC(xx), DIMDEC(uad)
      integer n
      integer lo(SDIM), hi(SDIM)
      REAL_T stx
      REAL_T s(DIMV(s))
      REAL_T xlo(DIMV(xx))
      REAL_T xhi(DIMV(xx))
      REAL_T uad(DIMV(uad))
      REAL_T eps
      integer xbc(SDIM,2)
      logical ltest
      integer j
      integer imin,jmin,imax,jmax

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

      if (xbc(1,1).eq.EXT_DIR) then
         if (n .eq. XVEL) then
            do j = jmin-1,jmax+1
              if (uad(imin,j) .ge. zero) then
                 xhi(imin,j) = s(imin-1,j)
                 xlo(imin,j) = s(imin-1,j)
              else 
                 xlo(imin,j) = xhi(imin,j)
              end if
            end do
         else
            do j = jmin-1,jmax+1
               ltest = uad(imin,j).le.eps
               stx   = merge(xhi(imin,j),s(imin-1,j),ltest)
               xhi(imin,j) = stx
               xlo(imin,j) = stx
            end do
         end if
      else if (xbc(1,1).eq.FOEXTRAP.or.xbc(1,1).eq.HOEXTRAP
     &        .or.xbc(1,1).eq.REFLECT_EVEN) then
         do j = jmin-1,jmax+1
            xlo(imin,j) = xhi(imin,j)
         end do
      else if (xbc(1,1).eq.REFLECT_ODD) then
         do j = jmin-1,jmax+1
            xhi(imin,j) = zero
            xlo(imin,j) = zero
         end do
      end if
c
      if (xbc(1,2).eq.EXT_DIR) then
         if (n .eq. XVEL) then
            do j = jmin-1,jmax+1
              if (uad(imax+1,j) .le. zero) then
                 xhi(imax+1,j) = s(imax+1,j)
                 xlo(imax+1,j) = s(imax+1,j)
               else 
                 xhi(imax+1,j) = xlo(imax+1,j)
               end if
             end do
         else
            do j = jmin-1,jmax+1
               ltest = uad(imax+1,j).ge.-eps
               stx   = merge(xlo(imax+1,j),s(imax+1,j),ltest)
               xhi(imax+1,j) = stx
               xlo(imax+1,j) = stx
            end do
         end if
      else if (xbc(1,2).eq.FOEXTRAP.or.xbc(1,2).eq.HOEXTRAP
     &        .or.xbc(1,2).eq.REFLECT_EVEN) then
         do j = jmin-1,jmax+1
            xhi(imax+1,j) = xlo(imax+1,j)
         end do
      else if (xbc(1,2).eq.REFLECT_ODD) then
         do j = jmin-1,jmax+1
            xhi(imax+1,j) = zero
            xlo(imax+1,j) = zero
         end do
      end if

      end

      subroutine trans_ybc(
     &     s,DIMS(s),
     &     ylo,yhi,DIMS(yy),vad,DIMS(vad),
     &     lo,hi,n,ybc,eps)
c
c     This subroutine processes boundary conditions on information
c     traced to cell faces in the y direction.  This is used for
c     computing velocities and edge states used in calculating
c     transverse derivatives
c
      implicit none
      integer DIMDEC(s), DIMDEC(yy), DIMDEC(vad)
      integer n
      integer lo(SDIM), hi(SDIM)
      REAL_T sty
      REAL_T s(DIMV(s))
      REAL_T ylo(DIMV(yy))
      REAL_T yhi(DIMV(yy))
      REAL_T vad(DIMV(vad))
      REAL_T eps
      integer ybc(SDIM,2)
      logical ltest
      integer i
      integer imin,jmin,imax,jmax

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

      if (ybc(2,1).eq.EXT_DIR) then
         if (n .eq. YVEL) then
            do i = imin-1,imax+1
              if (vad(i,jmin).ge.zero) then
                 yhi(i,jmin) = s(i,jmin-1)
                 ylo(i,jmin) = s(i,jmin-1)
              else
                 ylo(i,jmin) = yhi(i,jmin)
              end if
            end do
         else
            do i = imin-1,imax+1
               ltest = vad(i,jmin).le.eps
               sty   = merge(yhi(i,jmin),s(i,jmin-1),ltest)
               yhi(i,jmin) = sty
               ylo(i,jmin) = sty
            end do
         end if
      else if (ybc(2,1).eq.FOEXTRAP.or.ybc(2,1).eq.HOEXTRAP
     &        .or.ybc(2,1).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
            ylo(i,jmin) = yhi(i,jmin)
         end do
      else if (ybc(2,1).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
            yhi(i,jmin) = zero
            ylo(i,jmin) = zero
         end do
      end if

      if (ybc(2,2).eq.EXT_DIR) then
         if (n .eq. YVEL) then
            do i = imin-1,imax+1
              if (vad(i,jmax+1).le.zero) then
                 ylo(i,jmax+1) = s(i,jmax+1)
                 yhi(i,jmax+1) = s(i,jmax+1)
              else
                 yhi(i,jmax+1) = ylo(i,jmax+1)
              end if
            end do
         else
            do i = imin-1,imax+1
               ltest = vad(i,jmax+1).ge.-eps
               sty   = merge(ylo(i,jmax+1),s(i,jmax+1),ltest)
               yhi(i,jmax+1) = sty
               ylo(i,jmax+1) = sty
            end do
         end if
      else if (ybc(2,2).eq.FOEXTRAP.or.ybc(2,2).eq.HOEXTRAP
     &        .or.ybc(2,2).eq.REFLECT_EVEN) then
         do i = imin-1,imax+1
            yhi(i,jmax+1) = ylo(i,jmax+1)
         end do
      else if (ybc(2,2).eq.REFLECT_ODD) then
         do i = imin-1,imax+1
            ylo(i,jmax+1) = zero
            yhi(i,jmax+1) = zero
         end do
      end if

      end

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

#include "GODCOMM_F.H"

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
c
c ::: ::::: added to prevent underflow for small s values
c
      do j = lo(2)-ng, hi(2)+ng
        do i = lo(1)-ng, hi(1)+ng
	  s(i,j) = merge(s(i,j), zero, abs(s(i,j)).gt.1.0D-20)
        end do
      end do
c
c     COMPUTE 0TH order slopes
c
      if (slope_order.eq.1) then
        if (ng .lt. 1) then
	   call bl_abort("FORT_SLOPES: not enough bndry cells for 1st order")
        end if
        do j = jmin-1,jmax+1
           do i = imin-1,imax+1
              slx(i,j) = zero
              sly(i,j) = zero
  	   end do
        end do
      end if
c
c     COMPUTE 2nd order slopes
c
      if (slope_order.eq.2) then
         if (ng .lt. 2) then
            call bl_abort("SLOPE_2D: not enough bndry cells for 2nd order")
         end if
c
c     ------------------------ x slopes
c
         if (use_unlimited_slopes) then
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  slx(i,j) = half*(s(i+1,j)-s(i-1,j))
               end do
            end do
            if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
                  slx(imin-1,j) = zero
                  slx(imin,j)   = (s(imin+1,j)+three*s(imin,j)-four*s(imin-1,j))/three
               end do
            end if
            if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
                  slx(imax+1,j) = zero
                  slx(imax,j)   = -(s(imax-1,j)+three*s(imax,j)-four*s(imax+1,j))/three
               end do
            end if
         else
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  del  = half*(s(i+1,j)-s(i-1,j))
                  dpls = two*(s(i+1,j) - s(i ,j))
                  dmin = two*(s(i ,j) - s(i-1,j))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  slx(i,j)= sflg*min(slim,abs(del))
               end do
            end do
            
            if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
                  slx(imin-1,j) = zero
                  del  = (s(imin+1,j)+three*s(imin,j)-four*s(imin-1,j))/three
                  dpls = two*(s(imin+1,j) - s(imin  ,j))
                  dmin = two*(s(imin  ,j) - s(imin-1,j))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  slx(imin,j)= sflg*min(slim,abs(del))
               end do
            end if
            if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
                  slx(imax+1,j) = zero
                  del  = -(s(imax-1,j)+three*s(imax,j)-four*s(imax+1,j))/three
                  dpls = two*(s(imax+1,j) - s(imax  ,j))
                  dmin = two*(s(imax  ,j) - s(imax-1,j))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  slx(imax,j)= sflg*min(slim,abs(del))
               end do
            end if
         end if
c
c     ------------------------ y slopes
c
         if (use_unlimited_slopes) then
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  sly(i,j) = half*(s(i,j+1)-s(i,j-1))
               end do
            end do
            if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
               do i = imin-1, imax+1
                  sly(i,jmin-1) = zero
                  sly(i,jmin) = (s(i,jmin+1)+three*s(i,jmin)-four*s(i,jmin-1))/three
               end do
            end if
            if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
               do i = imin-1, imax+1
                  sly(i,jmax+1) = zero
                  sly(i,jmax) = -(s(i,jmax-1)+three*s(i,jmax)-four*s(i,jmax+1))/three
               end do
            end if
         else
            do j = jmin-1,jmax+1
               do i = imin-1,imax+1
                  del  = half*(s(i,j+1)-s(i,j-1))
                  dpls = two*(s(i,j+1) - s(i,j ))
                  dmin = two*(s(i,j ) - s(i,j-1))
                  slim = min(abs(dpls),abs(dmin))
                  slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  sly(i,j)= sflg*min(slim,abs(del))
               end do
            end do

            if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
               do i = imin-1, imax+1
                  sly(i,jmin-1) = zero
                  del  = (s(i,jmin+1)+three*s(i,jmin)-four*s(i,jmin-1))/three
                  dpls = two*(s(i,jmin+1) - s(i,jmin  ))
                  dmin = two*(s(i,jmin  ) - s(i,jmin-1))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  sly(i,jmin)= sflg*min(slim,abs(del))
               end do
            end if
            if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
               do i = imin-1, imax+1
                  sly(i,jmax+1) = zero
                  del  = -(s(i,jmax-1)+three*s(i,jmax)-four*s(i,jmax+1))/three
                  dpls = two*(s(i,jmax+1) - s(i,jmax  ))
                  dmin = two*(s(i,jmax  ) - s(i,jmax-1))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  sly(i,jmax)= sflg*min(slim,abs(del))
               end do
            end if
         end if
c
c ... end, if slope_order .eq. 2
c
      end if
c
c     COMPUTE 4TH order slopes
c
      if (slope_order.eq.4) then
         if (ng .lt. 3) then
            call bl_abort("SLOPE_2D: not enough bndry cells for 4th order")
         end if
c         
c     ------------------------ x slopes
c
         if (use_unlimited_slopes) then
            do j = jmin-1,jmax+1
               do i = imin-2,imax+2
                  slxscr(i,cen)  = half*(s(i+1,j)-s(i-1,j))
               end do
               do i = imin-1,imax+1
                  slx(i,j) = two * two3rd * slxscr(i,cen) -
     &                 sixth * (slxscr(i+1,cen) + slxscr(i-1,cen))
               end do
            end do
            
            if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
                  slx(imin,j) = -sixteen/fifteen*s(imin-1,j) + half*s(imin,j) + 
     &                 two3rd*s(imin+1,j) - tenth*s(imin+2,j)
                  slx(imin-1,j) = zero
               end do
            end if
            if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
               do j = jmin-1, jmax+1
                  slx(imax,j) = -( -sixteen/fifteen*s(imax+1,j) + half*s(imax,j) + 
     &                 two3rd*s(imax-1,j) - tenth*s(imax-2,j) )
                  slx(imax+1,j) = zero
               end do
            end if
         else
            do j = jmin-1,jmax+1
               do i = imin-2,imax+2
                  dmin           =  two*(s(i,  j)-s(i-1,j))
                  dpls           =  two*(s(i+1,j)-s(i  ,j))
                  slxscr(i,cen)  = half*(s(i+1,j)-s(i-1,j))
                  slxscr(i,lim)  = min(abs(dmin),abs(dpls))
                  slxscr(i,lim)  = merge(slxscr(i,lim),zero,(dpls*dmin) .ge. 0.0d0)
                  slxscr(i,flag) = sign(one,slxscr(i,cen))
                  slxscr(i,fromm)= slxscr(i,flag)*
     &                 min(slxscr(i,lim),abs(slxscr(i,cen)))
               end do

               do i = imin-1,imax+1
                  ds = two * two3rd * slxscr(i,cen) -
     &                 sixth * (slxscr(i+1,fromm) + slxscr(i-1,fromm))
                  slx(i,j) = slxscr(i,flag)*min(abs(ds),slxscr(i,lim))
               end do

               if (bc(1,1) .eq. EXT_DIR .or. bc(1,1) .eq. HOEXTRAP) then
                  del  = -sixteen/fifteen*s(imin-1,j) + half*s(imin,j) + 
     &                 two3rd*s(imin+1,j) - tenth*s(imin+2,j)
                  dmin = two*(s(imin  ,j)-s(imin-1,j))
                  dpls = two*(s(imin+1,j)-s(imin  ,j))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  slx(imin-1,j) = zero
                  slx(imin,  j) = sflg*min(slim,abs(del))

c                 Recalculate the slope at imin+1 using the revised slxscr(imin,fromm)
                  slxscr(imin,fromm) = slx(imin,j)
                  ds = two * two3rd * slxscr(imin+1,cen) -
     $               sixth * (slxscr(imin+2,fromm) + slxscr(imin,fromm))
                  slx(imin+1,j) = slxscr(imin+1,flag)*min(abs(ds),slxscr(imin+1,lim))
               end if

               if (bc(1,2) .eq. EXT_DIR .or. bc(1,2) .eq. HOEXTRAP) then
                  del  = -( -sixteen/fifteen*s(imax+1,j) + half*s(imax,j) + 
     &                 two3rd*s(imax-1,j) - tenth*s(imax-2,j) )
                  dmin = two*(s(imax  ,j)-s(imax-1,j))
                  dpls = two*(s(imax+1,j)-s(imax  ,j))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  slx(imax,  j) = sflg*min(slim,abs(del))
                  slx(imax+1,j) = zero

c                 Recalculate the slope at imax-1 using the revised slxscr(imax,fromm)
                  slxscr(imax,fromm) = slx(imax,j)
                  ds = two * two3rd * slxscr(imax-1,cen) -
     $               sixth * (slxscr(imax-2,fromm) + slxscr(imax,fromm))
                  slx(imax-1,j) = slxscr(imax-1,flag)*min(abs(ds),slxscr(imax-1,lim))
               end if
            end do
         end if
c
c     ------------------------ y slopes
c
         if (use_unlimited_slopes) then
            do i = imin-1,imax+1
               do j = jmin-2,jmax+2
                  slyscr(j,cen)  = half*(s(i,j+1)-s(i,j-1))
               end do
               do j = jmin-1,jmax+1
                  sly(i,j) = two * two3rd * slyscr(j,cen) -
     &                 sixth * (slyscr(j+1,cen) + slyscr(j-1,cen))
               end do
            end do
            
            if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
               do i = imin-1, imax+1
                  sly(i,jmin-1) = zero
                  sly(i,jmin) = -sixteen/fifteen*s(i,jmin-1) + half*s(i,jmin) + 
     &                 two3rd*s(i,jmin+1) - tenth*s(i,jmin+2)
               end do
            end if
            if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
               do i = imin-1, imax+1
                  sly(i,jmax) = -( -sixteen/fifteen*s(i,jmax+1) + half*s(i,jmax) + 
     &                 two3rd*s(i,jmax-1) - tenth*s(i,jmax-2) )
                  sly(i,jmax+1) = zero
               end do
            end if
         else
            do i = imin-1,imax+1
              do j = jmin-2,jmax+2
                  dmin           =  two*(s(i,j  )-s(i,j-1))
                  dpls           =  two*(s(i,j+1)-s(i,j  ))
                  slyscr(j,cen)  = half*(s(i,j+1)-s(i,j-1))
                  slyscr(j,lim)  = min(abs(dmin),abs(dpls))
                  slyscr(j,lim)  = merge(slyscr(j,lim),zero,(dpls*dmin) .ge. 0.0d0)
                  slyscr(j,flag) = sign(one,slyscr(j,cen))
                  slyscr(j,fromm)= slyscr(j,flag)*
     &                 min(slyscr(j,lim),abs(slyscr(j,cen)))
               end do
               do j = jmin-1,jmax+1
                  ds = two * two3rd * slyscr(j,cen) -
     &                 sixth * (slyscr(j+1,fromm) + slyscr(j-1,fromm))
                  sly(i,j) = slyscr(j,flag)*min(abs(ds),slyscr(j,lim))
               end do

               if (bc(2,1) .eq. EXT_DIR .or. bc(2,1) .eq. HOEXTRAP) then
                  del  = -sixteen/fifteen*s(i,jmin-1) + half*s(i,jmin) + 
     &                 two3rd*s(i,jmin+1) - tenth*s(i,jmin+2)
                  dmin = two*(s(i,jmin  )-s(i,jmin-1))
                  dpls = two*(s(i,jmin+1)-s(i,jmin  ))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  sly(i,jmin-1) = zero
                  sly(i,jmin  ) = sflg*min(slim,abs(del))

c                 Recalculate the slope at jmin+1 using the revised slyscr(jmin,fromm)
                  slyscr(jmin,fromm) = sly(i,jmin)
                  ds = two * two3rd * slyscr(jmin+1,cen) -
     $               sixth * (slyscr(jmin+2,fromm) + slyscr(jmin,fromm))
                  sly(i,jmin+1) = slyscr(jmin+1,flag)*min(abs(ds),slyscr(jmin+1,lim))
               end if

               if (bc(2,2) .eq. EXT_DIR .or. bc(2,2) .eq. HOEXTRAP) then
                  del  = -( -sixteen/fifteen*s(i,jmax+1) + half*s(i,jmax) + 
     &                 two3rd*s(i,jmax-1) - tenth*s(i,jmax-2) )
                  dmin = two*(s(i,jmax  )-s(i,jmax-1))
                  dpls = two*(s(i,jmax+1)-s(i,jmax  ))
                  slim = min(abs(dpls), abs(dmin))
                  slim = merge(slim, zero, (dpls*dmin) .ge. 0.0d0)
                  sflg = sign(one,del)
                  sly(i,jmax  ) = sflg*min(slim,abs(del))
                  sly(i,jmax+1) = zero

c                 Recalculate the slope at jmax-1 using the revised slyscr(jmax,fromm)
                  slyscr(jmax,fromm) = sly(i,jmax)
                  ds = two * two3rd * slyscr(jmax-1,cen) -
     $               sixth * (slyscr(jmax-2,fromm) + slyscr(jmax,fromm))
                  sly(i,jmax-1) = slyscr(jmax-1,flag)*min(abs(ds),slyscr(jmax-1,lim))
               end if
           end do
         end if
c
c ... end, if slope_order .eq. 4
c
      end if

      end

      subroutine FORT_PPM(s,u,v,DIMS(s),
     &     Ipx,Imx,Ipy,Imy,DIMS(work),sm,sp,DIMS(smp),dsvl,DIMS(dsvl),
     &     sedgex,DIMS(sedgex),sedgey,DIMS(sedgey),lo,hi,dx,dt,bc,eps,
     &     ppm_type)

      implicit none

      integer lo(SDIM), hi(SDIM)
      integer bc(SDIM,2)
      integer DIMDEC(s)
      integer DIMDEC(work)
      integer DIMDEC(smp)
      integer DIMDEC(dsvl)
      integer DIMDEC(sedgex)
      integer DIMDEC(sedgey)
      REAL_T  s(DIMV(s))
      REAL_T  u(DIMV(s))
      REAL_T  v(DIMV(s))
      REAL_T  Ipx(DIMV(work))
      REAL_T  Imx(DIMV(work))
      REAL_T  Ipy(DIMV(work))
      REAL_T  Imy(DIMV(work))
      REAL_T   sm(DIMV(smp))
      REAL_T   sp(DIMV(smp))
      REAL_T dsvl(DIMV(dsvl))
      REAL_T sedgex(DIMV(sedgex))
      REAL_T sedgey(DIMV(sedgey))
      REAL_T eps
      REAL_T dx(SDIM)
      REAL_T dt
      integer ppm_type

c     local
      integer i,j

      logical extremum, bigp, bigm

      REAL_T dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, C, alphap, alpham
      REAL_T sgn, sigma, s6, w0cc, amax, delam, delap
      REAL_T dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin, dachkm, dachkp

c     constant used in Colella 2008
      C = 1.25d0

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     x-direction
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c     compute s at x-edges
      if (ppm_type .eq. 1) then

c     compute van Leer slopes in x-direction
         dsvl = 0.d0
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-2,hi(1)+2
               dsc = 0.5d0 * (s(i+1,j) - s(i-1,j))
               dsl = 2.d0  * (s(i  ,j) - s(i-1,j))
               dsr = 2.d0  * (s(i+1,j) - s(i  ,j))
               if (dsl*dsr .gt. 0.d0)
     &              dsvl(i,j) = sign(1.d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
            end do
         end do
         
c     interpolate s to x-edges
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+2
               sedgex(i,j) = 0.5d0*(s(i,j)+s(i-1,j)) - (1.d0/6.d0)*(dsvl(i,j)-dsvl(i-1,j))
c     make sure sedgex lies in between adjacent cell-centered values
               sedgex(i,j) = max(sedgex(i,j),min(s(i,j),s(i-1,j)))
               sedgex(i,j) = min(sedgex(i,j),max(s(i,j),s(i-1,j)))
            end do
         end do

c     copy sedgex into sp and sm
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               sp(i,j) = sedgex(i+1,j)
               sm(i,j) = sedgex(i  ,j)
            end do
         end do

c     modify using quadratic limiters
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end do
         
c     different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
         if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
c     the value in the first cc ghost cell represents the edge value
            sm(lo(1),lo(2)-1:hi(2)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1)

c     use a modified stencil to get sedgex on the first interior edge
            sedgex(lo(1)+1,lo(2)-1:hi(2)+1) = 
     &           -(1.d0/5.d0)  *s(lo(1)-1,lo(2)-1:hi(2)+1)
     &           +(3.d0/4.d0)  *s(lo(1)  ,lo(2)-1:hi(2)+1)
     &           +0.5d0        *s(lo(1)+1,lo(2)-1:hi(2)+1)
     &           -(1.d0/20.0d0)*s(lo(1)+2,lo(2)-1:hi(2)+1)

c     make sure sedgex lies in between adjacent cell-centered values
            do j=lo(2)-1,hi(2)+1
               sedgex(lo(1)+1,j) = max(sedgex(lo(1)+1,j),min(s(lo(1)+1,j),s(lo(1),j)))
               sedgex(lo(1)+1,j) = min(sedgex(lo(1)+1,j),max(s(lo(1)+1,j),s(lo(1),j)))
            end do

c     copy sedgex into sp and sm
            do j=lo(2)-1,hi(2)+1
               sp(lo(1)  ,j) = sedgex(lo(1)+1,j)
               sm(lo(1)+1,j) = sedgex(lo(1)+1,j)
            end do

c     reset sp on second interior edge
            do j=lo(2)-1,hi(2)+1
               sp(lo(1)+1,j) = sedgex(lo(1)+2,j)
            end do

c     modify using quadratic limiters
            do j=lo(2)-1,hi(2)+1
               i = lo(1)+1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end if
         
         if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
c     the value in the first cc ghost cell represents the edge value
            sp(hi(1),lo(2)-1:hi(2)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1)

c     use a modified stencil to get sedgex on the first interior edge
            sedgex(hi(1),lo(2)-1:hi(2)+1) =
     &           -(1.d0/5.d0)  *s(hi(1)+1,lo(2)-1:hi(2)+1)
     &           +(3.d0/4.d0)  *s(hi(1)  ,lo(2)-1:hi(2)+1)
     &           +0.5d0        *s(hi(1)-1,lo(2)-1:hi(2)+1)
     &           -(1.d0/20.0d0)*s(hi(1)-2,lo(2)-1:hi(2)+1)

c     make sure sedgex lies in between adjacent cell-centered values
            do j=lo(2)-1,hi(2)+1
               sedgex(hi(1),j) = max(sedgex(hi(1),j),min(s(hi(1)-1,j),s(hi(1),j)))
               sedgex(hi(1),j) = min(sedgex(hi(1),j),max(s(hi(1)-1,j),s(hi(1),j)))
            end do

c     copy sedgex into sp and sm
            do j=lo(2)-1,hi(2)+1
               sp(hi(1)-1,j) = sedgex(hi(1),j)
               sm(hi(1)  ,j) = sedgex(hi(1),j)
            end do

c     reset sm on second interior edge
            do j=lo(2)-1,hi(2)+1
               sm(hi(1)-1,j) = sedgex(hi(1)-1,j)
            end do

c     modify using quadratic limiters
            do j=lo(2)-1,hi(2)+1
               i = hi(1)-1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end if

      else if (ppm_type .eq. 2) then

c     interpolate s to x-edges
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-2,hi(1)+3
               sedgex(i,j) = (7.d0/12.d0)*(s(i-1,j)+s(i,j)) - (1.d0/12.d0)*(s(i-2,j)+s(i+1,j))
c     limit sedgex
               if ((sedgex(i,j)-s(i-1,j))*(s(i,j)-sedgex(i,j)) .lt. 0.d0) then
                  D2  = 3.d0*(s(i-1,j)-2.d0*sedgex(i,j)+s(i,j))
                  D2L = s(i-2,j)-2.d0*s(i-1,j)+s(i,j)
                  D2R = s(i-1,j)-2.d0*s(i,j)+s(i+1,j)
                  sgn = sign(1.d0,D2)
                  D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
                  sedgex(i,j) = 0.5d0*(s(i-1,j)+s(i,j)) - (1.d0/6.d0)*D2LIM
               end if
            end do
         end do

c     use Colella 2008 limiters
c     This is a new version of the algorithm 
c     to eliminate sensitivity to roundoff.
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1

               alphap = sedgex(i+1,j)-s(i,j)
               alpham = sedgex(i  ,j)-s(i,j)
               bigp = abs(alphap).gt.2.d0*abs(alpham)
               bigm = abs(alpham).gt.2.d0*abs(alphap)
               extremum = .false.

               if (alpham*alphap .ge. 0.d0) then
                  extremum = .true.
               else if (bigp .or. bigm) then
c     Possible extremum. We look at cell centered values and face
c     centered values for a change in sign in the differences adjacent to
c     the cell. We use the pair of differences whose minimum magnitude is the
c     largest, and thus least susceptible to sensitivity to roundoff.
                  dafacem = sedgex(i,j) - sedgex(i-1,j)
                  dafacep = sedgex(i+2,j) - sedgex(i+1,j)
                  dabarm = s(i,j) - s(i-1,j)
                  dabarp = s(i+1,j) - s(i,j)
                  dafacemin = min(abs(dafacem),abs(dafacep))
                  dabarmin= min(abs(dabarm),abs(dabarp))
                  if (dafacemin.ge.dabarmin) then
                     dachkm = dafacem
                     dachkp = dafacep
                  else
                     dachkm = dabarm
                     dachkp = dabarp
                  endif
                  extremum = (dachkm*dachkp .le. 0.d0)
               end if

               if (extremum) then
                  D2  = 6.d0*(alpham + alphap)
                  D2L = s(i-2,j)-2.d0*s(i-1,j)+s(i,j)
                  D2R = s(i,j)-2.d0*s(i+1,j)+s(i+2,j)
                  D2C = s(i-1,j)-2.d0*s(i,j)+s(i+1,j)
                  sgn = sign(1.d0,D2)
                  D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                  alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                  alphap = alphap*D2LIM/max(abs(D2),1.d-10)
               else
                  if (bigp) then
                     sgn = sign(1.d0,alpham)
                     amax = -alphap**2 / (4*(alpham + alphap))
                     delam = s(i-1,j) - s(i,j)
                     if (sgn*amax .ge. sgn*delam) then
                        if (sgn*(delam - alpham).ge.1.d-10) then
                           alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                        else 
                           alphap = -2.d0*alpham
                        endif
                     endif
                  end if
                  if (bigm) then
                     sgn = sign(1.d0,alphap)
                     amax = -alpham**2 / (4*(alpham + alphap))
                     delap = s(i+1,j) - s(i,j)
                     if (sgn*amax .ge. sgn*delap) then
                        if (sgn*(delap - alphap).ge.1.d-10) then
                           alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                        else
                           alpham = -2.d0*alphap
                        endif
                     endif
                  end if
               end if

               sm(i,j) = s(i,j) + alpham
               sp(i,j) = s(i,j) + alphap

            end do
         end do

c     different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
         if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
c     the value in the first cc ghost cell represents the edge value
            sm(lo(1),lo(2)-1:hi(2)+1)    = s(lo(1)-1,lo(2)-1:hi(2)+1)
            sedgex(lo(1),lo(2)-1:hi(2)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1)

c     use a modified stencil to get sedgex on the first interior edge
            sedgex(lo(1)+1,lo(2)-1:hi(2)+1) =
     &           -(1.d0/5.d0)  *s(lo(1)-1,lo(2)-1:hi(2)+1)
     &           +(3.d0/4.d0)  *s(lo(1)  ,lo(2)-1:hi(2)+1)
     &           +0.5d0        *s(lo(1)+1,lo(2)-1:hi(2)+1)
     &           -(1.d0/20.0d0)*s(lo(1)+2,lo(2)-1:hi(2)+1)

c     make sure sedgex lies in between adjacent cell-centered values
            do j=lo(2)-1,hi(2)+1
               sedgex(lo(1)+1,j) = max(sedgex(lo(1)+1,j),min(s(lo(1)+1,j),s(lo(1),j)))
               sedgex(lo(1)+1,j) = min(sedgex(lo(1)+1,j),max(s(lo(1)+1,j),s(lo(1),j)))
            end do

c     copy sedgex into sp
            do j=lo(2)-1,hi(2)+1
               sp(lo(1)  ,j) = sedgex(lo(1)+1,j)
            end do

c     apply Colella 2008 limiters to compute sm and sp in the second
c     and third inner cells
            do j=lo(2)-1,hi(2)+1
               do i=lo(1)+1,lo(1)+2

                  alphap = sedgex(i+1,j)-s(i,j)
                  alpham = sedgex(i  ,j)-s(i,j)
                  bigp = abs(alphap).gt.2.d0*abs(alpham)
                  bigm = abs(alpham).gt.2.d0*abs(alphap)
                  extremum = .false.

                  if (alpham*alphap .ge. 0.d0) then
                     extremum = .true.
                  else if (bigp .or. bigm) then
c     Possible extremum. We look at cell centered values and face
c     centered values for a change in sign in the differences adjacent to
c     the cell. We use the pair of differences whose minimum magnitude is the
c     largest, and thus least susceptible to sensitivity to roundoff.
                     dafacem = sedgex(i,j) - sedgex(i-1,j)
                     dafacep = sedgex(i+2,j) - sedgex(i+1,j)
                     dabarm = s(i,j) - s(i-1,j)
                     dabarp = s(i+1,j) - s(i,j)
                     dafacemin = min(abs(dafacem),abs(dafacep))
                     dabarmin= min(abs(dabarm),abs(dabarp))
                     if (dafacemin.ge.dabarmin) then
                        dachkm = dafacem
                        dachkp = dafacep
                     else
                        dachkm = dabarm
                        dachkp = dabarp
                     endif
                     extremum = (dachkm*dachkp .le. 0.d0)
                  end if

                  if (extremum) then
                     D2  = 6.d0*(alpham + alphap)
                     D2L = s(i-2,j)-2.d0*s(i-1,j)+s(i,j)
                     D2R = s(i,j)-2.d0*s(i+1,j)+s(i+2,j)
                     D2C = s(i-1,j)-2.d0*s(i,j)+s(i+1,j)
                     sgn = sign(1.d0,D2)
                     D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                     alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                     alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                  else
                     if (bigp) then
                        sgn = sign(1.d0,alpham)
                        amax = -alphap**2 / (4*(alpham + alphap))
                        delam = s(i-1,j) - s(i,j)
                        if (sgn*amax .ge. sgn*delam) then
                           if (sgn*(delam - alpham).ge.1.d-10) then
                              alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                           else 
                              alphap = -2.d0*alpham
                           endif
                        endif
                     end if
                     if (bigm) then
                        sgn = sign(1.d0,alphap)
                        amax = -alpham**2 / (4*(alpham + alphap))
                        delap = s(i+1,j) - s(i,j)
                        if (sgn*amax .ge. sgn*delap) then
                           if (sgn*(delap - alphap).ge.1.d-10) then
                              alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                           else
                              alpham = -2.d0*alphap
                           endif
                        endif
                     end if
                  end if

                  sm(i,j) = s(i,j) + alpham
                  sp(i,j) = s(i,j) + alphap

               end do
            end do
         end if

         if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
c     the value in the first cc ghost cell represents the edge value
            sp(hi(1),lo(2)-1:hi(2)+1)      = s(hi(1)+1,lo(2)-1:hi(2)+1)
            sedgex(hi(1)+1,lo(2)-1:hi(2)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1)

c     use a modified stencil to get sedgex on the first interior edge
            sedgex(hi(1),lo(2)-1:hi(2)+1) =
     &           -(1.d0/5.d0)  *s(hi(1)+1,lo(2)-1:hi(2)+1)
     &           +(3.d0/4.d0)  *s(hi(1)  ,lo(2)-1:hi(2)+1)
     &           +0.5d0        *s(hi(1)-1,lo(2)-1:hi(2)+1)
     &           -(1.d0/20.0d0)*s(hi(1)-2,lo(2)-1:hi(2)+1)

c     make sure sedgex lies in between adjacent cell-centered values
            do j=lo(2)-1,hi(2)+1
               sedgex(hi(1),j) = max(sedgex(hi(1),j),min(s(hi(1)-1,j),s(hi(1),j)))
               sedgex(hi(1),j) = min(sedgex(hi(1),j),max(s(hi(1)-1,j),s(hi(1),j)))
            end do

c     copy sedgex into sm
            do j=lo(2)-1,hi(2)+1
               sm(hi(1)  ,j) = sedgex(hi(1),j)
            end do

c     reset sm on second interior edge
            do j=lo(2)-1,hi(2)+1
               sm(hi(1)-1,j) = sedgex(hi(1)-1,j)
            end do

c     apply Colella 2008 limiters to compute sm and sp in the second
c     and third inner cells
            do j=lo(2)-1,hi(2)+1
               do i=hi(1)-2,hi(1)-1

                  alphap = sedgex(i+1,j)-s(i,j)
                  alpham = sedgex(i  ,j)-s(i,j)
                  bigp = abs(alphap).gt.2.d0*abs(alpham)
                  bigm = abs(alpham).gt.2.d0*abs(alphap)
                  extremum = .false.

                  if (alpham*alphap .ge. 0.d0) then
                     extremum = .true.
                  else if (bigp .or. bigm) then
c     Possible extremum. We look at cell centered values and face
c     centered values for a change in sign in the differences adjacent to
c     the cell. We use the pair of differences whose minimum magnitude is the
c     largest, and thus least susceptible to sensitivity to roundoff.
                     dafacem = sedgex(i,j) - sedgex(i-1,j)
                     dafacep = sedgex(i+2,j) - sedgex(i+1,j)
                     dabarm = s(i,j) - s(i-1,j)
                     dabarp = s(i+1,j) - s(i,j)
                     dafacemin = min(abs(dafacem),abs(dafacep))
                     dabarmin= min(abs(dabarm),abs(dabarp))
                     if (dafacemin.ge.dabarmin) then
                        dachkm = dafacem
                        dachkp = dafacep
                     else
                        dachkm = dabarm
                        dachkp = dabarp
                     endif
                     extremum = (dachkm*dachkp .le. 0.d0)
                  end if

                  if (extremum) then
                     D2  = 6.d0*(alpham + alphap)
                     D2L = s(i-2,j)-2.d0*s(i-1,j)+s(i,j)
                     D2R = s(i,j)-2.d0*s(i+1,j)+s(i+2,j)
                     D2C = s(i-1,j)-2.d0*s(i,j)+s(i+1,j)
                     sgn = sign(1.d0,D2)
                     D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                     alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                     alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                  else
                     if (bigp) then
                        sgn = sign(1.d0,alpham)
                        amax = -alphap**2 / (4*(alpham + alphap))
                        delam = s(i-1,j) - s(i,j)
                        if (sgn*amax .ge. sgn*delam) then
                           if (sgn*(delam - alpham).ge.1.d-10) then
                              alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                           else 
                              alphap = -2.d0*alpham
                           endif
                        endif
                     end if
                     if (bigm) then
                        sgn = sign(1.d0,alphap)
                        amax = -alpham**2 / (4*(alpham + alphap))
                        delap = s(i+1,j) - s(i,j)
                        if (sgn*amax .ge. sgn*delap) then
                           if (sgn*(delap - alphap).ge.1.d-10) then
                              alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                           else
                              alpham = -2.d0*alphap
                           endif
                        endif
                     end if
                  end if

                  sm(i,j) = s(i,j) + alpham
                  sp(i,j) = s(i,j) + alphap

               end do
            end do
         end if

      end if

c     compute x-component of Ip and Im
      do j=lo(2)-1,hi(2)+1
         do i=lo(1)-1,hi(1)+1
            sigma = abs(u(i,j))*dt/dx(1)
            s6 = 6.0d0*s(i,j) - 3.0d0*(sm(i,j)+sp(i,j))
            if (u(i,j) .gt. eps) then
               Ipx(i,j) = sp(i,j) - (sigma/2.0d0)*
     &              (sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
               Imx(i,j) = s(i,j)
            else if (u(i,j) .lt. -eps) then
               Ipx(i,j) = s(i,j)
               Imx(i,j) = sm(i,j) + (sigma/2.0d0)*
     &              (sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
            else
               Ipx(i,j) = s(i,j)
               Imx(i,j) = s(i,j)
            end if
         end do
      end do

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     y-direction
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c     compute s at y-edges
      if (ppm_type .eq. 1) then

c     compute van Leer slopes in y-direction
         dsvl = 0.d0
         do j=lo(2)-2,hi(2)+2
            do i=lo(1)-1,hi(1)+1
               dsc = 0.5d0 * (s(i,j+1) - s(i,j-1))
               dsl = 2.d0  * (s(i,j  ) - s(i,j-1))
               dsr = 2.d0  * (s(i,j+1) - s(i,j  ))
               if (dsl*dsr .gt. 0.d0)
     &              dsvl(i,j) = sign(1.d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
            end do
         end do
         
c     interpolate s to y-edges
         do j=lo(2)-1,hi(2)+2
            do i=lo(1)-1,hi(1)+1
               sedgey(i,j) = 0.5d0*(s(i,j)+s(i,j-1)) - (1.d0/6.d0)*(dsvl(i,j)-dsvl(i,j-1))
c     make sure sedgey lies in between adjacent cell-centered values
               sedgey(i,j) = max(sedgey(i,j),min(s(i,j),s(i,j-1)))
               sedgey(i,j) = min(sedgey(i,j),max(s(i,j),s(i,j-1)))
            end do
         end do

c     copy sedgey into sp and sm
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               sp(i,j) = sedgey(i,j+1)
               sm(i,j) = sedgey(i,j  )
            end do
         end do

c     modify using quadratic limiters
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end do
         
c     different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
         if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
c     the value in the first cc ghost cell represents the edge value
            sm(lo(1)-1:hi(1)+1,lo(2)) = s(lo(1)-1:hi(1)+1,lo(2)-1)

c     use a modified stencil to get sedgey on the first interior edge
            sedgey(lo(1)-1:hi(1)+1,lo(2)+1) =
     &           -(1.d0/5.d0)  *s(lo(1)-1:hi(1)+1,lo(2)-1)
     &           +(3.d0/4.d0)  *s(lo(1)-1:hi(1)+1,lo(2)  )
     &           +0.5d0        *s(lo(1)-1:hi(1)+1,lo(2)+1)
     &           -(1.d0/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)+2)

c     make sure sedgey lies in between adjacent cell-centered values
            do i=lo(1)-1,hi(1)+1
               sedgey(i,lo(2)+1) = max(sedgey(i,lo(2)+1),min(s(i,lo(2)+1),s(i,lo(2))))
               sedgey(i,lo(2)+1) = min(sedgey(i,lo(2)+1),max(s(i,lo(2)+1),s(i,lo(2))))
            end do

c     copy sedgey into sp and sm
            do i=lo(1)-1,hi(1)+1
               sp(i,lo(2)  ) = sedgey(i,lo(2)+1)
               sm(i,lo(2)+1) = sedgey(i,lo(2)+1)
            end do

c     reset sp on second interior edge
            do i=lo(1)-1,hi(1)+1
               sp(i,lo(2)+1) = sedgey(i,lo(2)+2)
            end do

c     modify using quadratic limiters
            do i=lo(1)-1,hi(1)+1
               j = lo(2)+1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end if

         if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
c     the value in the first cc ghost cell represents the edge value
            sp(lo(1)-1:hi(1)+1,hi(2)) = s(lo(1)-1:hi(1)+1,hi(2)+1)

c     use a modified stencil to get sedgey on the first interior edge
            sedgey(lo(1)-1:hi(1)+1,hi(2)) =
     &           -(1.d0/5.d0)  *s(lo(1)-1:hi(1)+1,hi(2)+1)
     &           +(3.d0/4.d0)  *s(lo(1)-1:hi(1)+1,hi(2)  )
     &           +0.5d0        *s(lo(1)-1:hi(1)+1,hi(2)-1)
     &           -(1.d0/20.0d0)*s(lo(1)-1:hi(1)+1,hi(2)-2)

c     make sure sedgey lies in between adjacent cell-centered values
            do i=lo(1)-1,hi(1)+1
               sedgey(i,hi(2)) = max(sedgey(i,hi(2)),min(s(i,hi(2)-1),s(i,hi(2))))
               sedgey(i,hi(2)) = min(sedgey(i,hi(2)),max(s(i,hi(2)-1),s(i,hi(2))))
            end do

c     copy sedgey into sp and sm
            do i=lo(1)-1,hi(1)+1
               sp(i,hi(2)-1) = sedgey(i,hi(2))
               sm(i,hi(2)  ) = sedgey(i,hi(2))
            end do

c     reset sm on second interior edge
            do i=lo(1)-1,hi(1)+1
               sm(i,hi(2)-1) = sedgey(i,hi(2)-1)
            end do

c     modify using quadratic limiters
            do i=lo(1)-1,hi(1)+1
               j = hi(2)-1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end if

      else if (ppm_type .eq. 2) then
         
c     interpolate s to y-edges
         do j=lo(2)-2,hi(2)+3
            do i=lo(1)-1,hi(1)+1
               sedgey(i,j) = (7.d0/12.d0)*(s(i,j-1)+s(i,j)) - (1.d0/12.d0)*(s(i,j-2)+s(i,j+1))
c     limit sedgey
               if ((sedgey(i,j)-s(i,j-1))*(s(i,j)-sedgey(i,j)) .lt. 0.d0) then
                  D2  = 3.d0*(s(i,j-1)-2.d0*sedgey(i,j)+s(i,j))
                  D2L = s(i,j-2)-2.d0*s(i,j-1)+s(i,j)
                  D2R = s(i,j-1)-2.d0*s(i,j)+s(i,j+1)
                  sgn = sign(1.d0,D2)
                  D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
                  sedgey(i,j) = 0.5d0*(s(i,j-1)+s(i,j)) - (1.d0/6.d0)*D2LIM
               end if
            end do
         end do

c     use Colella 2008 limiters
c     This is a new version of the algorithm 
c     to eliminate sensitivity to roundoff.
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1

               alphap = sedgey(i,j+1)-s(i,j)
               alpham = sedgey(i,j  )-s(i,j)
               bigp = abs(alphap).gt.2.d0*abs(alpham)
               bigm = abs(alpham).gt.2.d0*abs(alphap)
               extremum = .false.

               if (alpham*alphap .ge. 0.d0) then
                  extremum = .true.
               else if (bigp .or. bigm) then
c     Possible extremum. We look at cell centered values and face
c     centered values for a change in sign in the differences adjacent to
c     the cell. We use the pair of differences whose minimum magnitude is the
c     largest, and thus least susceptible to sensitivity to roundoff.
                  dafacem = sedgey(i,j) - sedgey(i,j-1)
                  dafacep = sedgey(i,j+2) - sedgey(i,j+1)
                  dabarm = s(i,j) - s(i,j-1)
                  dabarp = s(i,j+1) - s(i,j)
                  dafacemin = min(abs(dafacem),abs(dafacep))
                  dabarmin= min(abs(dabarm),abs(dabarp))
                  if (dafacemin.ge.dabarmin) then
                     dachkm = dafacem
                     dachkp = dafacep
                  else
                     dachkm = dabarm
                     dachkp = dabarp
                  endif
                  extremum = (dachkm*dachkp .le. 0.d0)
               end if

               if (extremum) then
                  D2  = 6.d0*(alpham + alphap)
                  D2L = s(i,j-2)-2.d0*s(i,j-1)+s(i,j)
                  D2R = s(i,j)-2.d0*s(i,j+1)+s(i,j+2)
                  D2C = s(i,j-1)-2.d0*s(i,j)+s(i,j+1)
                  sgn = sign(1.d0,D2)
                  D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                  alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                  alphap = alphap*D2LIM/max(abs(D2),1.d-10)
               else
                  if (bigp) then
                     sgn = sign(1.d0,alpham)
                     amax = -alphap**2 / (4*(alpham + alphap))
                     delam = s(i,j-1) - s(i,j)
                     if (sgn*amax .ge. sgn*delam) then
                        if (sgn*(delam - alpham).ge.1.d-10) then
                           alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                        else 
                           alphap = -2.d0*alpham
                        endif
                     endif
                  end if
                  if (bigm) then
                     sgn = sign(1.d0,alphap)
                     amax = -alpham**2 / (4*(alpham + alphap))
                     delap = s(i,j+1) - s(i,j)
                     if (sgn*amax .ge. sgn*delap) then
                        if (sgn*(delap - alphap).ge.1.d-10) then
                           alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                        else
                           alpham = -2.d0*alphap
                        endif
                     endif
                  end if
               end if

               sm(i,j) = s(i,j) + alpham
               sp(i,j) = s(i,j) + alphap

            end do
         end do

c     different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
         if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
c     the value in the first cc ghost cell represents the edge value
            sm(lo(1)-1:hi(1)+1,lo(2))    = s(lo(1)-1:hi(1)+1,lo(2)-1)
            sedgey(lo(1)-1:hi(1)+1,lo(2)) = s(lo(1)-1:hi(1)+1,lo(2)-1)

c     use a modified stencil to get sedgey on the first interior edge
            sedgey(lo(1)-1:hi(1)+1,lo(2)+1) =
     &           -(1.d0/5.d0)  *s(lo(1)-1:hi(1)+1,lo(2)-1)
     &           +(3.d0/4.d0)  *s(lo(1)-1:hi(1)+1,lo(2)  )
     &           +0.5d0        *s(lo(1)-1:hi(1)+1,lo(2)+1)
     &           -(1.d0/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)+2)

c     make sure sedgey lies in between adjacent cell-centered values
            do i=lo(1)-1,hi(1)+1
               sedgey(i,lo(2)+1) = max(sedgey(i,lo(2)+1),min(s(i,lo(2)+1),s(i,lo(2))))
               sedgey(i,lo(2)+1) = min(sedgey(i,lo(2)+1),max(s(i,lo(2)+1),s(i,lo(2))))
            end do

c     copy sedgey into sp
            do i=lo(1)-1,hi(1)+1
               sp(i,lo(2)  ) = sedgey(i,lo(2)+1)
            end do

c     apply Colella 2008 limiters to compute sm and sp in the second
c     and third inner cells
            do j=lo(2)+1,lo(2)+1
               do i=lo(1)-1,hi(1)+1

                  alphap = sedgey(i,j+1)-s(i,j)
                  alpham = sedgey(i,j  )-s(i,j)
                  bigp = abs(alphap).gt.2.d0*abs(alpham)
                  bigm = abs(alpham).gt.2.d0*abs(alphap)
                  extremum = .false.

                  if (alpham*alphap .ge. 0.d0) then
                     extremum = .true.
                  else if (bigp .or. bigm) then
c     Possible extremum. We look at cell centered values and face
c     centered values for a change in sign in the differences adjacent to
c     the cell. We use the pair of differences whose minimum magnitude is the
c     largest, and thus least susceptible to sensitivity to roundoff.
                     dafacem = sedgey(i,j) - sedgey(i,j-1)
                     dafacep = sedgey(i,j+2) - sedgey(i,j+1)
                     dabarm = s(i,j) - s(i,j-1)
                     dabarp = s(i,j+1) - s(i,j)
                     dafacemin = min(abs(dafacem),abs(dafacep))
                     dabarmin= min(abs(dabarm),abs(dabarp))
                     if (dafacemin.ge.dabarmin) then
                        dachkm = dafacem
                        dachkp = dafacep
                     else
                        dachkm = dabarm
                        dachkp = dabarp
                     endif
                     extremum = (dachkm*dachkp .le. 0.d0)
                  end if

                  if (extremum) then
                     D2  = 6.d0*(alpham + alphap)
                     D2L = s(i,j-2)-2.d0*s(i,j-1)+s(i,j)
                     D2R = s(i,j)-2.d0*s(i,j+1)+s(i,j+2)
                     D2C = s(i,j-1)-2.d0*s(i,j)+s(i,j+1)
                     sgn = sign(1.d0,D2)
                     D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                     alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                     alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                  else
                     if (bigp) then
                        sgn = sign(1.d0,alpham)
                        amax = -alphap**2 / (4*(alpham + alphap))
                        delam = s(i,j-1) - s(i,j)
                        if (sgn*amax .ge. sgn*delam) then
                           if (sgn*(delam - alpham).ge.1.d-10) then
                              alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                           else 
                              alphap = -2.d0*alpham
                           endif
                        endif
                     end if
                     if (bigm) then
                        sgn = sign(1.d0,alphap)
                        amax = -alpham**2 / (4*(alpham + alphap))
                        delap = s(i,j+1) - s(i,j)
                        if (sgn*amax .ge. sgn*delap) then
                           if (sgn*(delap - alphap).ge.1.d-10) then
                              alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                           else
                              alpham = -2.d0*alphap
                           endif
                        endif
                     end if
                  end if

                  sm(i,j) = s(i,j) + alpham
                  sp(i,j) = s(i,j) + alphap

               end do
            end do
         end if

         if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
c     the value in the first cc ghost cell represents the edge value
            sp(lo(1)-1:hi(1)+1,hi(2))      = s(lo(1)-1:hi(1)+1,hi(2)+1)
            sedgey(lo(1)-1:hi(1)+1,hi(2)+1) = s(lo(1)-1:hi(1)+1,hi(2)+1)

c     use a modified stencil to get sedgey on the first interior edge
            sedgey(lo(1)-1:hi(1)+1,hi(2)) =
     &           -(1.d0/5.d0)  *s(lo(1)-1:hi(1)+1,hi(2)+1)
     &           +(3.d0/4.d0)  *s(lo(1)-1:hi(1)+1,hi(2)  )
     &           +0.5d0        *s(lo(1)-1:hi(1)+1,hi(2)-1)
     &           -(1.d0/20.0d0)*s(lo(1)-1:hi(1)+1,hi(2)-2)

c     make sure sedgey lies in between adjacent cell-centered values
            do i=lo(1)-1,hi(1)+1
               sedgey(i,hi(2)) = max(sedgey(i,hi(2)),min(s(i,hi(2)-1),s(i,hi(2))))
               sedgey(i,hi(2)) = min(sedgey(i,hi(2)),max(s(i,hi(2)-1),s(i,hi(2))))
            end do

c     copy sedgey into sm
            do i=lo(1)-1,hi(1)+1
               sm(i,hi(2)  ) = sedgey(i,hi(2))
            end do

c     apply Colella 2008 limiters to compute sm and sp in the second
c     and third inner cells
            do j=hi(2)-2,hi(2)-1
               do i=lo(1)-1,hi(1)+1

                  alphap = sedgey(i,j+1)-s(i,j)
                  alpham = sedgey(i,j  )-s(i,j)
                  bigp = abs(alphap).gt.2.d0*abs(alpham)
                  bigm = abs(alpham).gt.2.d0*abs(alphap)
                  extremum = .false.

                  if (alpham*alphap .ge. 0.d0) then
                     extremum = .true.
                  else if (bigp .or. bigm) then
c     Possible extremum. We look at cell centered values and face
c     centered values for a change in sign in the differences adjacent to
c     the cell. We use the pair of differences whose minimum magnitude is the
c     largest, and thus least susceptible to sensitivity to roundoff.
                     dafacem = sedgey(i,j) - sedgey(i,j-1)
                     dafacep = sedgey(i,j+2) - sedgey(i,j+1)
                     dabarm = s(i,j) - s(i,j-1)
                     dabarp = s(i,j+1) - s(i,j)
                     dafacemin = min(abs(dafacem),abs(dafacep))
                     dabarmin= min(abs(dabarm),abs(dabarp))
                     if (dafacemin.ge.dabarmin) then
                        dachkm = dafacem
                        dachkp = dafacep
                     else
                        dachkm = dabarm
                        dachkp = dabarp
                     endif
                     extremum = (dachkm*dachkp .le. 0.d0)
                  end if

                  if (extremum) then
                     D2  = 6.d0*(alpham + alphap)
                     D2L = s(i,j-2)-2.d0*s(i,j-1)+s(i,j)
                     D2R = s(i,j)-2.d0*s(i,j+1)+s(i,j+2)
                     D2C = s(i,j-1)-2.d0*s(i,j)+s(i,j+1)
                     sgn = sign(1.d0,D2)
                     D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                     alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                     alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                  else
                     if (bigp) then
                        sgn = sign(1.d0,alpham)
                        amax = -alphap**2 / (4*(alpham + alphap))
                        delam = s(i,j-1) - s(i,j)
                        if (sgn*amax .ge. sgn*delam) then
                           if (sgn*(delam - alpham).ge.1.d-10) then
                              alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                           else 
                              alphap = -2.d0*alpham
                           endif
                        endif
                     end if
                     if (bigm) then
                        sgn = sign(1.d0,alphap)
                        amax = -alpham**2 / (4*(alpham + alphap))
                        delap = s(i,j+1) - s(i,j)
                        if (sgn*amax .ge. sgn*delap) then
                           if (sgn*(delap - alphap).ge.1.d-10) then
                              alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                           else
                              alpham = -2.d0*alphap
                           endif
                        endif
                     end if
                  end if

                  sm(i,j) = s(i,j) + alpham
                  sp(i,j) = s(i,j) + alphap

               end do
            end do
         end if

      end if

c     compute y-component of Ip and Im
      do j=lo(2)-1,hi(2)+1
         do i=lo(1)-1,hi(1)+1
            sigma = abs(v(i,j))*dt/dx(2)
            s6 = 6.0d0*s(i,j) - 3.0d0*(sm(i,j)+sp(i,j))
            if (v(i,j) .gt. eps) then
               Ipy(i,j) = sp(i,j) - (sigma/2.0d0)*
     &              (sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
               Imy(i,j) = s(i,j)
            else if (v(i,j) .lt. -eps) then
               Ipy(i,j) = s(i,j)
               Imy(i,j) = sm(i,j) + (sigma/2.0d0)*
     &              (sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
            else
               Ipy(i,j) = s(i,j)
               Imy(i,j) = s(i,j)
            end if
         end do
      end do

      end

      subroutine FORT_PPM_FPU(s,DIMS(s),uedge,vedge,DIMS(uedge),
     &     Ipx,Imx,Ipy,Imy,DIMS(work),sm,sp,DIMS(smp),dsvl,DIMS(dsvl),
     &     sedgex,DIMS(sedgex),sedgey,DIMS(sedgey),lo,hi,dx,dt,bc,eps,
     &     ppm_type)

      implicit none

      integer lo(SDIM), hi(SDIM)
      integer bc(SDIM,2)
      integer DIMDEC(s)
      integer DIMDEC(uedge)
      integer DIMDEC(work)
      integer DIMDEC(smp)
      integer DIMDEC(dsvl)
      integer DIMDEC(sedgex)
      integer DIMDEC(sedgey)
      REAL_T  s(DIMV(s))
      REAL_T  uedge(DIMV(uedge))
      REAL_T  vedge(DIMV(uedge))
      REAL_T  Ipx(DIMV(work))
      REAL_T  Imx(DIMV(work))
      REAL_T  Ipy(DIMV(work))
      REAL_T  Imy(DIMV(work))
      REAL_T   sm(DIMV(smp))
      REAL_T   sp(DIMV(smp))
      REAL_T dsvl(DIMV(dsvl))
      REAL_T sedgex(DIMV(sedgex))
      REAL_T sedgey(DIMV(sedgey))
      REAL_T eps
      REAL_T dx(SDIM)
      REAL_T dt
      integer ppm_type

c     local
      integer i,j

      logical extremum, bigp, bigm

      REAL_T dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, C, alphap, alpham
      REAL_T sgn, sigmam, sigmap, s6, w0cc, amax, delam, delap
      REAL_T dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin, dachkm, dachkp

c     constant used in Colella 2008
      C = 1.25d0

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     x-direction
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c     compute s at x-edges
      if (ppm_type .eq. 1) then

c     compute van Leer slopes in x-direction
         dsvl = 0.d0
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-2,hi(1)+2
               dsc = 0.5d0 * (s(i+1,j) - s(i-1,j))
               dsl = 2.d0  * (s(i  ,j) - s(i-1,j))
               dsr = 2.d0  * (s(i+1,j) - s(i  ,j))
               if (dsl*dsr .gt. 0.d0)
     &              dsvl(i,j) = sign(1.d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
            end do
         end do
         
c     interpolate s to x-edges
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+2
               sedgex(i,j) = 0.5d0*(s(i,j)+s(i-1,j)) - (1.d0/6.d0)*(dsvl(i,j)-dsvl(i-1,j))
c     make sure sedgex lies in between adjacent cell-centered values
               sedgex(i,j) = max(sedgex(i,j),min(s(i,j),s(i-1,j)))
               sedgex(i,j) = min(sedgex(i,j),max(s(i,j),s(i-1,j)))
            end do
         end do

c     copy sedgex into sp and sm
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               sp(i,j) = sedgex(i+1,j)
               sm(i,j) = sedgex(i  ,j)
            end do
         end do

c     modify using quadratic limiters
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end do
         
c     different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
         if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
c     the value in the first cc ghost cell represents the edge value
            sm(lo(1),lo(2)-1:hi(2)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1)

c     use a modified stencil to get sedgex on the first interior edge
            sedgex(lo(1)+1,lo(2)-1:hi(2)+1) = 
     &           -(1.d0/5.d0)  *s(lo(1)-1,lo(2)-1:hi(2)+1)
     &           +(3.d0/4.d0)  *s(lo(1)  ,lo(2)-1:hi(2)+1)
     &           +0.5d0        *s(lo(1)+1,lo(2)-1:hi(2)+1)
     &           -(1.d0/20.0d0)*s(lo(1)+2,lo(2)-1:hi(2)+1)

c     make sure sedgex lies in between adjacent cell-centered values
            do j=lo(2)-1,hi(2)+1
               sedgex(lo(1)+1,j) = max(sedgex(lo(1)+1,j),min(s(lo(1)+1,j),s(lo(1),j)))
               sedgex(lo(1)+1,j) = min(sedgex(lo(1)+1,j),max(s(lo(1)+1,j),s(lo(1),j)))
            end do

c     copy sedgex into sp and sm
            do j=lo(2)-1,hi(2)+1
               sp(lo(1)  ,j) = sedgex(lo(1)+1,j)
               sm(lo(1)+1,j) = sedgex(lo(1)+1,j)
            end do

c     reset sp on second interior edge
            do j=lo(2)-1,hi(2)+1
               sp(lo(1)+1,j) = sedgex(lo(1)+2,j)
            end do

c     modify using quadratic limiters
            do j=lo(2)-1,hi(2)+1
               i = lo(1)+1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end if
         
         if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
c     the value in the first cc ghost cell represents the edge value
            sp(hi(1),lo(2)-1:hi(2)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1)

c     use a modified stencil to get sedgex on the first interior edge
            sedgex(hi(1),lo(2)-1:hi(2)+1) =
     &           -(1.d0/5.d0)  *s(hi(1)+1,lo(2)-1:hi(2)+1)
     &           +(3.d0/4.d0)  *s(hi(1)  ,lo(2)-1:hi(2)+1)
     &           +0.5d0        *s(hi(1)-1,lo(2)-1:hi(2)+1)
     &           -(1.d0/20.0d0)*s(hi(1)-2,lo(2)-1:hi(2)+1)

c     make sure sedgex lies in between adjacent cell-centered values
            do j=lo(2)-1,hi(2)+1
               sedgex(hi(1),j) = max(sedgex(hi(1),j),min(s(hi(1)-1,j),s(hi(1),j)))
               sedgex(hi(1),j) = min(sedgex(hi(1),j),max(s(hi(1)-1,j),s(hi(1),j)))
            end do

c     copy sedgex into sp and sm
            do j=lo(2)-1,hi(2)+1
               sp(hi(1)-1,j) = sedgex(hi(1),j)
               sm(hi(1)  ,j) = sedgex(hi(1),j)
            end do

c     reset sm on second interior edge
            do j=lo(2)-1,hi(2)+1
               sm(hi(1)-1,j) = sedgex(hi(1)-1,j)
            end do

c     modify using quadratic limiters
            do j=lo(2)-1,hi(2)+1
               i = hi(1)-1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end if

      else if (ppm_type .eq. 2) then

c     interpolate s to x-edges
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-2,hi(1)+3
               sedgex(i,j) = (7.d0/12.d0)*(s(i-1,j)+s(i,j)) - (1.d0/12.d0)*(s(i-2,j)+s(i+1,j))
c     limit sedgex
               if ((sedgex(i,j)-s(i-1,j))*(s(i,j)-sedgex(i,j)) .lt. 0.d0) then
                  D2  = 3.d0*(s(i-1,j)-2.d0*sedgex(i,j)+s(i,j))
                  D2L = s(i-2,j)-2.d0*s(i-1,j)+s(i,j)
                  D2R = s(i-1,j)-2.d0*s(i,j)+s(i+1,j)
                  sgn = sign(1.d0,D2)
                  D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
                  sedgex(i,j) = 0.5d0*(s(i-1,j)+s(i,j)) - (1.d0/6.d0)*D2LIM
               end if
            end do
         end do

c     use Colella 2008 limiters
c     This is a new version of the algorithm 
c     to eliminate sensitivity to roundoff.
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1

               alphap = sedgex(i+1,j)-s(i,j)
               alpham = sedgex(i  ,j)-s(i,j)
               bigp = abs(alphap).gt.2.d0*abs(alpham)
               bigm = abs(alpham).gt.2.d0*abs(alphap)
               extremum = .false.

               if (alpham*alphap .ge. 0.d0) then
                  extremum = .true.
               else if (bigp .or. bigm) then
c     Possible extremum. We look at cell centered values and face
c     centered values for a change in sign in the differences adjacent to
c     the cell. We use the pair of differences whose minimum magnitude is the
c     largest, and thus least susceptible to sensitivity to roundoff.
                  dafacem = sedgex(i,j) - sedgex(i-1,j)
                  dafacep = sedgex(i+2,j) - sedgex(i+1,j)
                  dabarm = s(i,j) - s(i-1,j)
                  dabarp = s(i+1,j) - s(i,j)
                  dafacemin = min(abs(dafacem),abs(dafacep))
                  dabarmin= min(abs(dabarm),abs(dabarp))
                  if (dafacemin.ge.dabarmin) then
                     dachkm = dafacem
                     dachkp = dafacep
                  else
                     dachkm = dabarm
                     dachkp = dabarp
                  endif
                  extremum = (dachkm*dachkp .le. 0.d0)
               end if

               if (extremum) then
                  D2  = 6.d0*(alpham + alphap)
                  D2L = s(i-2,j)-2.d0*s(i-1,j)+s(i,j)
                  D2R = s(i,j)-2.d0*s(i+1,j)+s(i+2,j)
                  D2C = s(i-1,j)-2.d0*s(i,j)+s(i+1,j)
                  sgn = sign(1.d0,D2)
                  D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                  alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                  alphap = alphap*D2LIM/max(abs(D2),1.d-10)
               else
                  if (bigp) then
                     sgn = sign(1.d0,alpham)
                     amax = -alphap**2 / (4*(alpham + alphap))
                     delam = s(i-1,j) - s(i,j)
                     if (sgn*amax .ge. sgn*delam) then
                        if (sgn*(delam - alpham).ge.1.d-10) then
                           alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                        else 
                           alphap = -2.d0*alpham
                        endif
                     endif
                  end if
                  if (bigm) then
                     sgn = sign(1.d0,alphap)
                     amax = -alpham**2 / (4*(alpham + alphap))
                     delap = s(i+1,j) - s(i,j)
                     if (sgn*amax .ge. sgn*delap) then
                        if (sgn*(delap - alphap).ge.1.d-10) then
                           alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                        else
                           alpham = -2.d0*alphap
                        endif
                     endif
                  end if
               end if

               sm(i,j) = s(i,j) + alpham
               sp(i,j) = s(i,j) + alphap

            end do
         end do

c     different stencil needed for x-component of EXT_DIR and HOEXTRAP bc's
         if (bc(1,1) .eq. EXT_DIR  .or. bc(1,1) .eq. HOEXTRAP) then
c     the value in the first cc ghost cell represents the edge value
            sm(lo(1),lo(2)-1:hi(2)+1)    = s(lo(1)-1,lo(2)-1:hi(2)+1)
            sedgex(lo(1),lo(2)-1:hi(2)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1)

c     use a modified stencil to get sedgex on the first interior edge
            sedgex(lo(1)+1,lo(2)-1:hi(2)+1) =
     &           -(1.d0/5.d0)  *s(lo(1)-1,lo(2)-1:hi(2)+1)
     &           +(3.d0/4.d0)  *s(lo(1)  ,lo(2)-1:hi(2)+1)
     &           +0.5d0        *s(lo(1)+1,lo(2)-1:hi(2)+1)
     &           -(1.d0/20.0d0)*s(lo(1)+2,lo(2)-1:hi(2)+1)

c     make sure sedgex lies in between adjacent cell-centered values
            do j=lo(2)-1,hi(2)+1
               sedgex(lo(1)+1,j) = max(sedgex(lo(1)+1,j),min(s(lo(1)+1,j),s(lo(1),j)))
               sedgex(lo(1)+1,j) = min(sedgex(lo(1)+1,j),max(s(lo(1)+1,j),s(lo(1),j)))
            end do

c     copy sedgex into sp
            do j=lo(2)-1,hi(2)+1
               sp(lo(1)  ,j) = sedgex(lo(1)+1,j)
            end do

c     apply Colella 2008 limiters to compute sm and sp in the second
c     and third inner cells
            do j=lo(2)-1,hi(2)+1
               do i=lo(1)+1,lo(1)+2

                  alphap = sedgex(i+1,j)-s(i,j)
                  alpham = sedgex(i  ,j)-s(i,j)
                  bigp = abs(alphap).gt.2.d0*abs(alpham)
                  bigm = abs(alpham).gt.2.d0*abs(alphap)
                  extremum = .false.

                  if (alpham*alphap .ge. 0.d0) then
                     extremum = .true.
                  else if (bigp .or. bigm) then
c     Possible extremum. We look at cell centered values and face
c     centered values for a change in sign in the differences adjacent to
c     the cell. We use the pair of differences whose minimum magnitude is the
c     largest, and thus least susceptible to sensitivity to roundoff.
                     dafacem = sedgex(i,j) - sedgex(i-1,j)
                     dafacep = sedgex(i+2,j) - sedgex(i+1,j)
                     dabarm = s(i,j) - s(i-1,j)
                     dabarp = s(i+1,j) - s(i,j)
                     dafacemin = min(abs(dafacem),abs(dafacep))
                     dabarmin= min(abs(dabarm),abs(dabarp))
                     if (dafacemin.ge.dabarmin) then
                        dachkm = dafacem
                        dachkp = dafacep
                     else
                        dachkm = dabarm
                        dachkp = dabarp
                     endif
                     extremum = (dachkm*dachkp .le. 0.d0)
                  end if

                  if (extremum) then
                     D2  = 6.d0*(alpham + alphap)
                     D2L = s(i-2,j)-2.d0*s(i-1,j)+s(i,j)
                     D2R = s(i,j)-2.d0*s(i+1,j)+s(i+2,j)
                     D2C = s(i-1,j)-2.d0*s(i,j)+s(i+1,j)
                     sgn = sign(1.d0,D2)
                     D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                     alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                     alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                  else
                     if (bigp) then
                        sgn = sign(1.d0,alpham)
                        amax = -alphap**2 / (4*(alpham + alphap))
                        delam = s(i-1,j) - s(i,j)
                        if (sgn*amax .ge. sgn*delam) then
                           if (sgn*(delam - alpham).ge.1.d-10) then
                              alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                           else 
                              alphap = -2.d0*alpham
                           endif
                        endif
                     end if
                     if (bigm) then
                        sgn = sign(1.d0,alphap)
                        amax = -alpham**2 / (4*(alpham + alphap))
                        delap = s(i+1,j) - s(i,j)
                        if (sgn*amax .ge. sgn*delap) then
                           if (sgn*(delap - alphap).ge.1.d-10) then
                              alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                           else
                              alpham = -2.d0*alphap
                           endif
                        endif
                     end if
                  end if

                  sm(i,j) = s(i,j) + alpham
                  sp(i,j) = s(i,j) + alphap

               end do
            end do
         end if

         if (bc(1,2) .eq. EXT_DIR  .or. bc(1,2) .eq. HOEXTRAP) then
c     the value in the first cc ghost cell represents the edge value
            sp(hi(1),lo(2)-1:hi(2)+1)      = s(hi(1)+1,lo(2)-1:hi(2)+1)
            sedgex(hi(1)+1,lo(2)-1:hi(2)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1)

c     use a modified stencil to get sedgex on the first interior edge
            sedgex(hi(1),lo(2)-1:hi(2)+1) =
     &           -(1.d0/5.d0)  *s(hi(1)+1,lo(2)-1:hi(2)+1)
     &           +(3.d0/4.d0)  *s(hi(1)  ,lo(2)-1:hi(2)+1)
     &           +0.5d0        *s(hi(1)-1,lo(2)-1:hi(2)+1)
     &           -(1.d0/20.0d0)*s(hi(1)-2,lo(2)-1:hi(2)+1)

c     make sure sedgex lies in between adjacent cell-centered values
            do j=lo(2)-1,hi(2)+1
               sedgex(hi(1),j) = max(sedgex(hi(1),j),min(s(hi(1)-1,j),s(hi(1),j)))
               sedgex(hi(1),j) = min(sedgex(hi(1),j),max(s(hi(1)-1,j),s(hi(1),j)))
            end do

c     copy sedgex into sm
            do j=lo(2)-1,hi(2)+1
               sm(hi(1)  ,j) = sedgex(hi(1),j)
            end do

c     reset sm on second interior edge
            do j=lo(2)-1,hi(2)+1
               sm(hi(1)-1,j) = sedgex(hi(1)-1,j)
            end do

c     apply Colella 2008 limiters to compute sm and sp in the second
c     and third inner cells
            do j=lo(2)-1,hi(2)+1
               do i=hi(1)-2,hi(1)-1

                  alphap = sedgex(i+1,j)-s(i,j)
                  alpham = sedgex(i  ,j)-s(i,j)
                  bigp = abs(alphap).gt.2.d0*abs(alpham)
                  bigm = abs(alpham).gt.2.d0*abs(alphap)
                  extremum = .false.

                  if (alpham*alphap .ge. 0.d0) then
                     extremum = .true.
                  else if (bigp .or. bigm) then
c     Possible extremum. We look at cell centered values and face
c     centered values for a change in sign in the differences adjacent to
c     the cell. We use the pair of differences whose minimum magnitude is the
c     largest, and thus least susceptible to sensitivity to roundoff.
                     dafacem = sedgex(i,j) - sedgex(i-1,j)
                     dafacep = sedgex(i+2,j) - sedgex(i+1,j)
                     dabarm = s(i,j) - s(i-1,j)
                     dabarp = s(i+1,j) - s(i,j)
                     dafacemin = min(abs(dafacem),abs(dafacep))
                     dabarmin= min(abs(dabarm),abs(dabarp))
                     if (dafacemin.ge.dabarmin) then
                        dachkm = dafacem
                        dachkp = dafacep
                     else
                        dachkm = dabarm
                        dachkp = dabarp
                     endif
                     extremum = (dachkm*dachkp .le. 0.d0)
                  end if

                  if (extremum) then
                     D2  = 6.d0*(alpham + alphap)
                     D2L = s(i-2,j)-2.d0*s(i-1,j)+s(i,j)
                     D2R = s(i,j)-2.d0*s(i+1,j)+s(i+2,j)
                     D2C = s(i-1,j)-2.d0*s(i,j)+s(i+1,j)
                     sgn = sign(1.d0,D2)
                     D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                     alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                     alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                  else
                     if (bigp) then
                        sgn = sign(1.d0,alpham)
                        amax = -alphap**2 / (4*(alpham + alphap))
                        delam = s(i-1,j) - s(i,j)
                        if (sgn*amax .ge. sgn*delam) then
                           if (sgn*(delam - alpham).ge.1.d-10) then
                              alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                           else 
                              alphap = -2.d0*alpham
                           endif
                        endif
                     end if
                     if (bigm) then
                        sgn = sign(1.d0,alphap)
                        amax = -alpham**2 / (4*(alpham + alphap))
                        delap = s(i+1,j) - s(i,j)
                        if (sgn*amax .ge. sgn*delap) then
                           if (sgn*(delap - alphap).ge.1.d-10) then
                              alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                           else
                              alpham = -2.d0*alphap
                           endif
                        endif
                     end if
                  end if

                  sm(i,j) = s(i,j) + alpham
                  sp(i,j) = s(i,j) + alphap

               end do
            end do
         end if

      end if

c     compute x-component of Ip and Im
      do j=lo(2)-1,hi(2)+1
         do i=lo(1)-1,hi(1)+1
            sigmap = abs(uedge(i+1,j))*dt/dx(1)
            sigmam = abs(uedge(i,j))*dt/dx(1)
            s6 = 6.0d0*s(i,j) - 3.0d0*(sm(i,j)+sp(i,j))
            if (uedge(i+1,j) .gt. eps) then
               Ipx(i,j) = sp(i,j) - (sigmap/2.0d0)*
     &              (sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigmap)*s6)
            else
               Ipx(i,j) = s(i,j)
            end if
            if (uedge(i,j) .lt. -eps) then
               Imx(i,j) = sm(i,j) + (sigmam/2.0d0)*
     &              (sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigmam)*s6)
            else
               Imx(i,j) = s(i,j)
            end if
         end do
      end do

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     y-direction
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c     compute s at y-edges
      if (ppm_type .eq. 1) then

c     compute van Leer slopes in y-direction
         dsvl = 0.d0
         do j=lo(2)-2,hi(2)+2
            do i=lo(1)-1,hi(1)+1
               dsc = 0.5d0 * (s(i,j+1) - s(i,j-1))
               dsl = 2.d0  * (s(i,j  ) - s(i,j-1))
               dsr = 2.d0  * (s(i,j+1) - s(i,j  ))
               if (dsl*dsr .gt. 0.d0)
     &              dsvl(i,j) = sign(1.d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
            end do
         end do
         
c     interpolate s to y-edges
         do j=lo(2)-1,hi(2)+2
            do i=lo(1)-1,hi(1)+1
               sedgey(i,j) = 0.5d0*(s(i,j)+s(i,j-1)) - (1.d0/6.d0)*(dsvl(i,j)-dsvl(i,j-1))
c     make sure sedgey lies in between adjacent cell-centered values
               sedgey(i,j) = max(sedgey(i,j),min(s(i,j),s(i,j-1)))
               sedgey(i,j) = min(sedgey(i,j),max(s(i,j),s(i,j-1)))
            end do
         end do

c     copy sedgey into sp and sm
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               sp(i,j) = sedgey(i,j+1)
               sm(i,j) = sedgey(i,j  )
            end do
         end do

c     modify using quadratic limiters
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end do
         
c     different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
         if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
c     the value in the first cc ghost cell represents the edge value
            sm(lo(1)-1:hi(1)+1,lo(2)) = s(lo(1)-1:hi(1)+1,lo(2)-1)

c     use a modified stencil to get sedgey on the first interior edge
            sedgey(lo(1)-1:hi(1)+1,lo(2)+1) =
     &           -(1.d0/5.d0)  *s(lo(1)-1:hi(1)+1,lo(2)-1)
     &           +(3.d0/4.d0)  *s(lo(1)-1:hi(1)+1,lo(2)  )
     &           +0.5d0        *s(lo(1)-1:hi(1)+1,lo(2)+1)
     &           -(1.d0/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)+2)

c     make sure sedgey lies in between adjacent cell-centered values
            do i=lo(1)-1,hi(1)+1
               sedgey(i,lo(2)+1) = max(sedgey(i,lo(2)+1),min(s(i,lo(2)+1),s(i,lo(2))))
               sedgey(i,lo(2)+1) = min(sedgey(i,lo(2)+1),max(s(i,lo(2)+1),s(i,lo(2))))
            end do

c     copy sedgey into sp and sm
            do i=lo(1)-1,hi(1)+1
               sp(i,lo(2)  ) = sedgey(i,lo(2)+1)
               sm(i,lo(2)+1) = sedgey(i,lo(2)+1)
            end do

c     reset sp on second interior edge
            do i=lo(1)-1,hi(1)+1
               sp(i,lo(2)+1) = sedgey(i,lo(2)+2)
            end do

c     modify using quadratic limiters
            do i=lo(1)-1,hi(1)+1
               j = lo(2)+1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end if

         if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
c     the value in the first cc ghost cell represents the edge value
            sp(lo(1)-1:hi(1)+1,hi(2)) = s(lo(1)-1:hi(1)+1,hi(2)+1)

c     use a modified stencil to get sedgey on the first interior edge
            sedgey(lo(1)-1:hi(1)+1,hi(2)) =
     &           -(1.d0/5.d0)  *s(lo(1)-1:hi(1)+1,hi(2)+1)
     &           +(3.d0/4.d0)  *s(lo(1)-1:hi(1)+1,hi(2)  )
     &           +0.5d0        *s(lo(1)-1:hi(1)+1,hi(2)-1)
     &           -(1.d0/20.0d0)*s(lo(1)-1:hi(1)+1,hi(2)-2)

c     make sure sedgey lies in between adjacent cell-centered values
            do i=lo(1)-1,hi(1)+1
               sedgey(i,hi(2)) = max(sedgey(i,hi(2)),min(s(i,hi(2)-1),s(i,hi(2))))
               sedgey(i,hi(2)) = min(sedgey(i,hi(2)),max(s(i,hi(2)-1),s(i,hi(2))))
            end do

c     copy sedgey into sp and sm
            do i=lo(1)-1,hi(1)+1
               sp(i,hi(2)-1) = sedgey(i,hi(2))
               sm(i,hi(2)  ) = sedgey(i,hi(2))
            end do

c     reset sm on second interior edge
            do i=lo(1)-1,hi(1)+1
               sm(i,hi(2)-1) = sedgey(i,hi(2)-1)
            end do

c     modify using quadratic limiters
            do i=lo(1)-1,hi(1)+1
               j = hi(2)-1
               if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. 0.d0) then
                  sp(i,j) = s(i,j)
                  sm(i,j) = s(i,j)
               else if (abs(sp(i,j)-s(i,j)) .ge. 2.d0*abs(sm(i,j)-s(i,j))) then
                  sp(i,j) = 3.d0*s(i,j) - 2.d0*sm(i,j)
               else if (abs(sm(i,j)-s(i,j)) .ge. 2.d0*abs(sp(i,j)-s(i,j))) then
                  sm(i,j) = 3.d0*s(i,j) - 2.d0*sp(i,j)
               end if
            end do
         end if

      else if (ppm_type .eq. 2) then
         
c     interpolate s to y-edges
         do j=lo(2)-2,hi(2)+3
            do i=lo(1)-1,hi(1)+1
               sedgey(i,j) = (7.d0/12.d0)*(s(i,j-1)+s(i,j)) - (1.d0/12.d0)*(s(i,j-2)+s(i,j+1))
c     limit sedgey
               if ((sedgey(i,j)-s(i,j-1))*(s(i,j)-sedgey(i,j)) .lt. 0.d0) then
                  D2  = 3.d0*(s(i,j-1)-2.d0*sedgey(i,j)+s(i,j))
                  D2L = s(i,j-2)-2.d0*s(i,j-1)+s(i,j)
                  D2R = s(i,j-1)-2.d0*s(i,j)+s(i,j+1)
                  sgn = sign(1.d0,D2)
                  D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
                  sedgey(i,j) = 0.5d0*(s(i,j-1)+s(i,j)) - (1.d0/6.d0)*D2LIM
               end if
            end do
         end do

c     use Colella 2008 limiters
c     This is a new version of the algorithm 
c     to eliminate sensitivity to roundoff.
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1

               alphap = sedgey(i,j+1)-s(i,j)
               alpham = sedgey(i,j  )-s(i,j)
               bigp = abs(alphap).gt.2.d0*abs(alpham)
               bigm = abs(alpham).gt.2.d0*abs(alphap)
               extremum = .false.

               if (alpham*alphap .ge. 0.d0) then
                  extremum = .true.
               else if (bigp .or. bigm) then
c     Possible extremum. We look at cell centered values and face
c     centered values for a change in sign in the differences adjacent to
c     the cell. We use the pair of differences whose minimum magnitude is the
c     largest, and thus least susceptible to sensitivity to roundoff.
                  dafacem = sedgey(i,j) - sedgey(i,j-1)
                  dafacep = sedgey(i,j+2) - sedgey(i,j+1)
                  dabarm = s(i,j) - s(i,j-1)
                  dabarp = s(i,j+1) - s(i,j)
                  dafacemin = min(abs(dafacem),abs(dafacep))
                  dabarmin= min(abs(dabarm),abs(dabarp))
                  if (dafacemin.ge.dabarmin) then
                     dachkm = dafacem
                     dachkp = dafacep
                  else
                     dachkm = dabarm
                     dachkp = dabarp
                  endif
                  extremum = (dachkm*dachkp .le. 0.d0)
               end if

               if (extremum) then
                  D2  = 6.d0*(alpham + alphap)
                  D2L = s(i,j-2)-2.d0*s(i,j-1)+s(i,j)
                  D2R = s(i,j)-2.d0*s(i,j+1)+s(i,j+2)
                  D2C = s(i,j-1)-2.d0*s(i,j)+s(i,j+1)
                  sgn = sign(1.d0,D2)
                  D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                  alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                  alphap = alphap*D2LIM/max(abs(D2),1.d-10)
               else
                  if (bigp) then
                     sgn = sign(1.d0,alpham)
                     amax = -alphap**2 / (4*(alpham + alphap))
                     delam = s(i,j-1) - s(i,j)
                     if (sgn*amax .ge. sgn*delam) then
                        if (sgn*(delam - alpham).ge.1.d-10) then
                           alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                        else 
                           alphap = -2.d0*alpham
                        endif
                     endif
                  end if
                  if (bigm) then
                     sgn = sign(1.d0,alphap)
                     amax = -alpham**2 / (4*(alpham + alphap))
                     delap = s(i,j+1) - s(i,j)
                     if (sgn*amax .ge. sgn*delap) then
                        if (sgn*(delap - alphap).ge.1.d-10) then
                           alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                        else
                           alpham = -2.d0*alphap
                        endif
                     endif
                  end if
               end if

               sm(i,j) = s(i,j) + alpham
               sp(i,j) = s(i,j) + alphap

            end do
         end do

c     different stencil needed for y-component of EXT_DIR and HOEXTRAP bc's
         if (bc(2,1) .eq. EXT_DIR  .or. bc(2,1) .eq. HOEXTRAP) then
c     the value in the first cc ghost cell represents the edge value
            sm(lo(1)-1:hi(1)+1,lo(2))    = s(lo(1)-1:hi(1)+1,lo(2)-1)
            sedgey(lo(1)-1:hi(1)+1,lo(2)) = s(lo(1)-1:hi(1)+1,lo(2)-1)

c     use a modified stencil to get sedgey on the first interior edge
            sedgey(lo(1)-1:hi(1)+1,lo(2)+1) =
     &           -(1.d0/5.d0)  *s(lo(1)-1:hi(1)+1,lo(2)-1)
     &           +(3.d0/4.d0)  *s(lo(1)-1:hi(1)+1,lo(2)  )
     &           +0.5d0        *s(lo(1)-1:hi(1)+1,lo(2)+1)
     &           -(1.d0/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)+2)

c     make sure sedgey lies in between adjacent cell-centered values
            do i=lo(1)-1,hi(1)+1
               sedgey(i,lo(2)+1) = max(sedgey(i,lo(2)+1),min(s(i,lo(2)+1),s(i,lo(2))))
               sedgey(i,lo(2)+1) = min(sedgey(i,lo(2)+1),max(s(i,lo(2)+1),s(i,lo(2))))
            end do

c     copy sedgey into sp
            do i=lo(1)-1,hi(1)+1
               sp(i,lo(2)  ) = sedgey(i,lo(2)+1)
            end do

c     apply Colella 2008 limiters to compute sm and sp in the second
c     and third inner cells
            do j=lo(2)+1,lo(2)+1
               do i=lo(1)-1,hi(1)+1

                  alphap = sedgey(i,j+1)-s(i,j)
                  alpham = sedgey(i,j  )-s(i,j)
                  bigp = abs(alphap).gt.2.d0*abs(alpham)
                  bigm = abs(alpham).gt.2.d0*abs(alphap)
                  extremum = .false.

                  if (alpham*alphap .ge. 0.d0) then
                     extremum = .true.
                  else if (bigp .or. bigm) then
c     Possible extremum. We look at cell centered values and face
c     centered values for a change in sign in the differences adjacent to
c     the cell. We use the pair of differences whose minimum magnitude is the
c     largest, and thus least susceptible to sensitivity to roundoff.
                     dafacem = sedgey(i,j) - sedgey(i,j-1)
                     dafacep = sedgey(i,j+2) - sedgey(i,j+1)
                     dabarm = s(i,j) - s(i,j-1)
                     dabarp = s(i,j+1) - s(i,j)
                     dafacemin = min(abs(dafacem),abs(dafacep))
                     dabarmin= min(abs(dabarm),abs(dabarp))
                     if (dafacemin.ge.dabarmin) then
                        dachkm = dafacem
                        dachkp = dafacep
                     else
                        dachkm = dabarm
                        dachkp = dabarp
                     endif
                     extremum = (dachkm*dachkp .le. 0.d0)
                  end if

                  if (extremum) then
                     D2  = 6.d0*(alpham + alphap)
                     D2L = s(i,j-2)-2.d0*s(i,j-1)+s(i,j)
                     D2R = s(i,j)-2.d0*s(i,j+1)+s(i,j+2)
                     D2C = s(i,j-1)-2.d0*s(i,j)+s(i,j+1)
                     sgn = sign(1.d0,D2)
                     D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                     alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                     alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                  else
                     if (bigp) then
                        sgn = sign(1.d0,alpham)
                        amax = -alphap**2 / (4*(alpham + alphap))
                        delam = s(i,j-1) - s(i,j)
                        if (sgn*amax .ge. sgn*delam) then
                           if (sgn*(delam - alpham).ge.1.d-10) then
                              alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                           else 
                              alphap = -2.d0*alpham
                           endif
                        endif
                     end if
                     if (bigm) then
                        sgn = sign(1.d0,alphap)
                        amax = -alpham**2 / (4*(alpham + alphap))
                        delap = s(i,j+1) - s(i,j)
                        if (sgn*amax .ge. sgn*delap) then
                           if (sgn*(delap - alphap).ge.1.d-10) then
                              alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                           else
                              alpham = -2.d0*alphap
                           endif
                        endif
                     end if
                  end if

                  sm(i,j) = s(i,j) + alpham
                  sp(i,j) = s(i,j) + alphap

               end do
            end do
         end if

         if (bc(2,2) .eq. EXT_DIR  .or. bc(2,2) .eq. HOEXTRAP) then
c     the value in the first cc ghost cell represents the edge value
            sp(lo(1)-1:hi(1)+1,hi(2))      = s(lo(1)-1:hi(1)+1,hi(2)+1)
            sedgey(lo(1)-1:hi(1)+1,hi(2)+1) = s(lo(1)-1:hi(1)+1,hi(2)+1)

c     use a modified stencil to get sedgey on the first interior edge
            sedgey(lo(1)-1:hi(1)+1,hi(2)) =
     &           -(1.d0/5.d0)  *s(lo(1)-1:hi(1)+1,hi(2)+1)
     &           +(3.d0/4.d0)  *s(lo(1)-1:hi(1)+1,hi(2)  )
     &           +0.5d0        *s(lo(1)-1:hi(1)+1,hi(2)-1)
     &           -(1.d0/20.0d0)*s(lo(1)-1:hi(1)+1,hi(2)-2)

c     make sure sedgey lies in between adjacent cell-centered values
            do i=lo(1)-1,hi(1)+1
               sedgey(i,hi(2)) = max(sedgey(i,hi(2)),min(s(i,hi(2)-1),s(i,hi(2))))
               sedgey(i,hi(2)) = min(sedgey(i,hi(2)),max(s(i,hi(2)-1),s(i,hi(2))))
            end do

c     copy sedgey into sm
            do i=lo(1)-1,hi(1)+1
               sm(i,hi(2)  ) = sedgey(i,hi(2))
            end do

c     apply Colella 2008 limiters to compute sm and sp in the second
c     and third inner cells
            do j=hi(2)-2,hi(2)-1
               do i=lo(1)-1,hi(1)+1

                  alphap = sedgey(i,j+1)-s(i,j)
                  alpham = sedgey(i,j  )-s(i,j)
                  bigp = abs(alphap).gt.2.d0*abs(alpham)
                  bigm = abs(alpham).gt.2.d0*abs(alphap)
                  extremum = .false.

                  if (alpham*alphap .ge. 0.d0) then
                     extremum = .true.
                  else if (bigp .or. bigm) then
c     Possible extremum. We look at cell centered values and face
c     centered values for a change in sign in the differences adjacent to
c     the cell. We use the pair of differences whose minimum magnitude is the
c     largest, and thus least susceptible to sensitivity to roundoff.
                     dafacem = sedgey(i,j) - sedgey(i,j-1)
                     dafacep = sedgey(i,j+2) - sedgey(i,j+1)
                     dabarm = s(i,j) - s(i,j-1)
                     dabarp = s(i,j+1) - s(i,j)
                     dafacemin = min(abs(dafacem),abs(dafacep))
                     dabarmin= min(abs(dabarm),abs(dabarp))
                     if (dafacemin.ge.dabarmin) then
                        dachkm = dafacem
                        dachkp = dafacep
                     else
                        dachkm = dabarm
                        dachkp = dabarp
                     endif
                     extremum = (dachkm*dachkp .le. 0.d0)
                  end if

                  if (extremum) then
                     D2  = 6.d0*(alpham + alphap)
                     D2L = s(i,j-2)-2.d0*s(i,j-1)+s(i,j)
                     D2R = s(i,j)-2.d0*s(i,j+1)+s(i,j+2)
                     D2C = s(i,j-1)-2.d0*s(i,j)+s(i,j+1)
                     sgn = sign(1.d0,D2)
                     D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
                     alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                     alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                  else
                     if (bigp) then
                        sgn = sign(1.d0,alpham)
                        amax = -alphap**2 / (4*(alpham + alphap))
                        delam = s(i,j-1) - s(i,j)
                        if (sgn*amax .ge. sgn*delam) then
                           if (sgn*(delam - alpham).ge.1.d-10) then
                              alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                           else 
                              alphap = -2.d0*alpham
                           endif
                        endif
                     end if
                     if (bigm) then
                        sgn = sign(1.d0,alphap)
                        amax = -alpham**2 / (4*(alpham + alphap))
                        delap = s(i,j+1) - s(i,j)
                        if (sgn*amax .ge. sgn*delap) then
                           if (sgn*(delap - alphap).ge.1.d-10) then
                              alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                           else
                              alpham = -2.d0*alphap
                           endif
                        endif
                     end if
                  end if

                  sm(i,j) = s(i,j) + alpham
                  sp(i,j) = s(i,j) + alphap

               end do
            end do
         end if

      end if

c     compute y-component of Ip and Im
      do j=lo(2)-1,hi(2)+1
         do i=lo(1)-1,hi(1)+1
            sigmap = abs(vedge(i,j+1))*dt/dx(2)
            sigmam = abs(vedge(i,j))*dt/dx(2)
            s6 = 6.0d0*s(i,j) - 3.0d0*(sm(i,j)+sp(i,j))
            if (vedge(i,j+1) .gt. eps) then
               Ipy(i,j) = sp(i,j) - (sigmap/2.0d0)*
     &              (sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigmap)*s6)
            else
               Ipy(i,j) = s(i,j)
            end if
            if (vedge(i,j) .lt. -eps) then
               Imy(i,j) = sm(i,j) + (sigmam/2.0d0)*
     &              (sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigmam)*s6)
            else
               Imy(i,j) = s(i,j)
            end if
         end do
      end do

      end

      subroutine FORT_CONVSCALMINMAX(s,DIMS(s),sn,DIMS(sn),
     &                               lo,hi,bc)
c
c     correct an convectively-advected field for under/over shoots
c
      implicit none
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer lo(SDIM), hi(SDIM)
      integer bc(SDIM,2)
      REAL_T s(DIMV(s))
      REAL_T sn(DIMV(sn))
      integer  i, j, imin, imax, jmin, jmax
      REAL_T   smin, smax

      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
c
c     set corner values for s
c
c     :: corner ilo, jlo
c
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
         s(imin-1,jmin-1) = s(imin,jmin)
      end if
c
c     :: corner ihi, jlo
c
      if (bc(1,2).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
         s(imax+1,jmin-1) = s(imax,jmin)
      end if
c
c     :: corner ilo, jhi
c
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
         s(imin-1,jmax+1) = s(imin,jmax)
      end if
c
c     :: corner ihi, jhi
c
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
         s(imax+1,jmax+1) = s(imax,jmax)
      end if
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
         end do
      end do

      end

      subroutine FORT_CONSSCALMINMAX(s,sn,rho,rhon,DIMS(s),lo,hi,bc)
c
c     correct an conservatively-advected field for under/over shoots
c
      implicit none
      integer DIMDEC(s)
      integer lo(SDIM), hi(SDIM)
      integer bc(SDIM,2)
      REAL_T    s(DIMV(s))
      REAL_T   sn(DIMV(s))
      REAL_T  rho(DIMV(s))
      REAL_T rhon(DIMV(s))
      integer  i, j, imin, imax, jmin, jmax
      REAL_T   slo_min,smd_min,shi_min
      REAL_T   slo_max,smd_max,shi_max
      REAL_T   smin, smax

      imin = lo(1)
      imax = hi(1)
      jmin = lo(2)
      jmax = hi(2)
c
c     set corner values for s
c
c     :: corner ilo, jlo
c
      if (bc(1,1).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
           s(imin-1,jmin-1) =   s(imin,jmin)
         rho(imin-1,jmin-1) = rho(imin,jmin)
      end if
c
c     :: corner ihi, jlo
c
      if (bc(1,2).ne.INT_DIR .or. bc(2,1).ne.INT_DIR) then
           s(imax+1,jmin-1) =   s(imax,jmin)
         rho(imax+1,jmin-1) = rho(imax,jmin)
      end if
c
c     :: corner ilo, jhi
c
      if (bc(1,1).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
           s(imin-1,jmax+1) =   s(imin,jmax)
         rho(imin-1,jmax+1) = rho(imin,jmax)
      end if
c
c     :: corner ihi, jhi
c
      if (bc(1,2).ne.INT_DIR .or. bc(2,2).ne.INT_DIR) then
           s(imax+1,jmax+1) =   s(imax,jmax)
         rho(imax+1,jmax+1) = rho(imax,jmax)
      end if
c
c     compute extrema of s
c
      do j = jmin-1, jmax+1         
         do i = imin-1, imax+1
            s(i,j) = s(i,j) / rho(i,j)
         end do
      end do

      do j = jmin, jmax         
         do i = imin, imax
            slo_min = min(s(i-1,j-1),s(i  ,j-1),s(i+1,j-1))
            smd_min = min(s(i-1,j  ),s(i  ,j  ),s(i+1,j  ))
            shi_min = min(s(i-1,j+1),s(i  ,j+1),s(i+1,j+1))
            smin = min(slo_min,smd_min,shi_min)

            slo_max = max(s(i-1,j-1),s(i  ,j-1),s(i+1,j-1))
            smd_max = max(s(i-1,j  ),s(i  ,j  ),s(i+1,j  ))
            shi_max = max(s(i-1,j+1),s(i  ,j+1),s(i+1,j+1))
            smax = max(slo_max,smd_max,shi_max)

            sn(i,j) = max(sn(i,j)/rhon(i,j),smin) * rhon(i,j)
            sn(i,j) = min(sn(i,j)/rhon(i,j),smax) * rhon(i,j)
         end do
      end do

      do j = jmin-1, jmax+1         
         do i = imin-1, imax+1
            s(i,j) = s(i,j) * rho(i,j)
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
      integer i, j, n
      integer DIMDEC(tf)
      integer DIMDEC(gp)
      integer DIMDEC(rho)
      integer lo(SDIM), hi(SDIM)
      REAL_T tforces(DIMV(tf),SDIM)
      REAL_T gp(DIMV(gp),SDIM)
      REAL_T rho(DIMV(rho))

      do n = 1, SDIM
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tforces(i,j,n) = (
     &              tforces(i,j,n)
     &              -    gp(i,j,n))/rho(i,j)
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

      do n = 1, SDIM
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tforces(i,j,n) = (
     &              tforces(i,j,n)
     &              +  visc(i,j,n)
     &              -    gp(i,j,n))/rho(i,j)
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
c     sum divU*S into tforces or divide tforces by rho
c     depending on the value of iconserv
c
      implicit none
      integer nvar, iconserv
      integer lo(SDIM), hi(SDIM)
      integer i, j, n

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
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,n) = 
     &            tforces(i,j,n) - S(i,j,n)*divu(i,j)
               end do
            end do
         end do
      else
         do n = 1, nvar
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,n) =
     &            tforces(i,j,n)/rho(i,j)
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
c     sum tforces, viscous foricing and divU*S into tforces
c     depending on the value of iconserv
c
      implicit none
      integer nvar, iconserv
      integer lo(SDIM), hi(SDIM)
      integer i, j, n

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
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,n) = 
     &                 tforces(i,j,n)
     &                 +  visc(i,j,n)
     &                 -     S(i,j,n)*divu(i,j)
               end do
            end do
         end do
      else
         do n = 1, nvar
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  tforces(i,j,n) = (
     &                 tforces(i,j,n)
     &                 +  visc(i,j,n) )/rho(i,j)
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
      integer i, j, n, nvar
      integer DIMDEC(s)
      integer DIMDEC(sn)
      integer DIMDEC(tf)
      integer lo(SDIM), hi(SDIM)
      REAL_T dt
      REAL_T s(DIMV(s),nvar)
      REAL_T sn(DIMV(sn),nvar)
      REAL_T tforces(DIMV(tf),nvar)

      do n = 1,nvar
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               sn(i,j,n) = s(i,j,n) + dt*tforces(i,j,n)
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
      integer i, j, n, nvar
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
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               sp(i,j,n) = ss(i,j,n)
     &              + hdt*(tfs(i,j,n)-tfn(i,j,n))
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

      do n = 1,nvar
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               sn(i,j,n) = s(i,j,n)
     &              - dt*aofs(i,j,n)
     &              + dt*tforces(i,j,n)
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

      do n = 1, SDIM
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               un(i,j,n) = u(i,j,n) 
     &              - dt*   aofs(i,j,n)
     &              + dt*tforces(i,j,n)/rho(i,j)
     &              - dt*     gp(i,j,n)/rho(i,j)
            end do
         end do
      end do

      end


      subroutine bdsslope(s,lo_1,lo_2,hi_1,hi_2,slx,sly,sc,dx)

      implicit none
      integer lo_1,lo_2,hi_1,hi_2
      REAL_T      s(lo_1-3:hi_1+3,lo_2-3:hi_2+3)
      REAL_T    slx(lo_1-1:hi_1+1,lo_2-1:hi_2+1)
      REAL_T    sly(lo_1-1:hi_1+1,lo_2-1:hi_2+1)
      REAL_T   slxy(lo_1-1:hi_1+1,lo_2-1:hi_2+1)
      REAL_T   sint(lo_1-2:hi_1+2,lo_2-2:hi_2+2)
      REAL_T     sc(lo_1-1:hi_1+1,lo_2-1:hi_2+1,4)
      REAL_T dx(2)     

      REAL_T   diff(lo_1-1:hi_1+1,4)
      REAL_T   smin(lo_1-1:hi_1+1,4)
      REAL_T   smax(lo_1-1:hi_1+1,4)
      REAL_T sumdif(lo_1-1:hi_1+1)
      REAL_T sgndif(lo_1-1:hi_1+1)
      integer   kdp(lo_1-1:hi_1+1)

      REAL_T hx,hy,sumloc,redfac,redmax,div
      REAL_T choice1,choice2
      REAL_T eps
      integer inc1, inc2, inc3, inc4
      integer i,j,k,ll,is,ie,js,je

      hx = dx(1)
      hy = dx(2)
      is = lo_1
      ie = hi_1
      js = lo_2
      je = hi_2

      eps = 1.d-8

      do i = is-2,ie+1
        do j = js-2,je+1
          sint(i,j) = (
     $              s(i-1,j-1) + s(i-1,j+2) + s(i+2,j-1) + s(i+2,j+2)
     $     - seven*(s(i-1,j  ) + s(i-1,j+1) + s(i  ,j-1) + s(i+1,j-1) + 
     $              s(i  ,j+2) + s(i+1,j+2) + s(i+2,j  ) + s(i+2,j+1)) + 
     $       49.d0*(s(i  ,j  ) + s(i+1,j  ) + s(i  ,j+1) + s(i+1,j+1)) ) / 144.d0
        enddo
      enddo

      do j = js-1,je+1
        do i = is-1,ie+1 

          slx(i,j) = half*(sint(i  ,j) + sint(i  ,j-1) - 
     $                     sint(i-1,j) - sint(i-1,j-1) ) / hx
          sly(i,j) = half*(sint(i  ,j) - sint(i  ,j-1) + 
     $                     sint(i-1,j) - sint(i-1,j-1) ) / hy
          slxy(i,j) = (sint(i,j  ) - sint(i  ,j-1) - 
     $                 sint(i-1,j) + sint(i-1,j-1) ) / (hx*hy)
        enddo
      enddo

      do j = js-1,je+1

        do i = is-1,ie+1
          smin(i,4) = min(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))
          smax(i,4) = max(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))
          smin(i,3) = min(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))
          smax(i,3) = max(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))
          smin(i,2) = min(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))
          smax(i,2) = max(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))
          smin(i,1) = min(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))
          smax(i,1) = max(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))

          sc(i,j,4) = s(i,j) + half*(hx*slx(i,j) + hy*sly(i,j)) 
     $                + fourth*hx*hy*slxy(i,j)
          sc(i,j,3) = s(i,j) + half*(hx*slx(i,j) - hy*sly(i,j)) 
     $                - fourth*hx*hy*slxy(i,j)
          sc(i,j,2) = s(i,j) - half*(hx*slx(i,j) - hy*sly(i,j))
     $                - fourth*hx*hy*slxy(i,j)
          sc(i,j,1) = s(i,j) - half*(hx*slx(i,j) + hy*sly(i,j))
     $                + fourth*hx*hy*slxy(i,j)

          sc(i,j,4) = max(min(sc(i,j,4), smax(i,4)), smin(i,4))
          sc(i,j,3) = max(min(sc(i,j,3), smax(i,3)), smin(i,3))
          sc(i,j,2) = max(min(sc(i,j,2), smax(i,2)), smin(i,2))
          sc(i,j,1) = max(min(sc(i,j,1), smax(i,1)), smin(i,1))
        enddo

        do ll = 1,3 
          do i = is-1,ie+1 
            sumloc = fourth*(sc(i,j,4) + sc(i,j,3) + 
     $                       sc(i,j,2) + sc(i,j,1))
            sumdif(i) = (sumloc - s(i,j))*4.
            sgndif(i) = sign(one,sumdif(i))

            diff(i,4) = (sc(i,j,4) - s(i,j))*sgndif(i)
            diff(i,3) = (sc(i,j,3) - s(i,j))*sgndif(i)
            diff(i,2) = (sc(i,j,2) - s(i,j))*sgndif(i)
            diff(i,1) = (sc(i,j,1) - s(i,j))*sgndif(i)

            inc1 = merge(1,0,(diff(i,1) - eps) .ge. 0.0d0)
            inc2 = merge(1,0,(diff(i,2) - eps) .ge. 0.0d0)
            inc3 = merge(1,0,(diff(i,3) - eps) .ge. 0.0d0)
            inc4 = merge(1,0,(diff(i,4) - eps) .ge. 0.0d0)
            kdp(i) = inc1 + inc2 + inc3 + inc4
          enddo

          do k = 1,4 
            do i = is-1,ie+1 
              div = merge(1,kdp(i),kdp(i) .lt. 1)
	      choice1 = sumdif(i)*sgndif(i)/div
              redfac = merge(choice1,zero,diff(i,k).gt.eps)
              kdp(i) = merge(kdp(i) - 1,kdp(i),diff(i,k) .gt. eps)
              choice1 = sc(i,j,k) - smin(i,k)
              choice2 = smax(i,k) - sc(i,j,k)
              redmax = merge(choice1,choice2,sgndif(i) .ge. 0.0d0)
              redfac = min(redfac,redmax)
              sumdif(i) = sumdif(i) - redfac*sgndif(i)
              sc(i,j,k) = sc(i,j,k) - redfac*sgndif(i)
            enddo
          enddo
        enddo

        do i = is-1,ie+1 
          slx(i,j) = half*( sc(i,j,4) + sc(i,j,3)
     $                     -sc(i,j,1) - sc(i,j,2))/hx
          sly(i,j) = half*( sc(i,j,4) + sc(i,j,2)
     $                     -sc(i,j,1) - sc(i,j,3))/hy
          slxy(i,j) = ( sc(i,j,4) + sc(i,j,1)
     $                 -sc(i,j,2) - sc(i,j,3) ) / (hx*hy)
        enddo
      enddo

      return
      end

      subroutine FORT_ESTATE_BDS(s, tforces, divu, DIMS(s),
     &     xlo, xhi, sx, slxscr, stxlo, stxhi,
     &     uedge, DIMS(uedge), xstate, DIMS(xstate),

     &     ylo, yhi, sy, slyscr, stylo, styhi,
     &     vedge, DIMS(vedge), ystate, DIMS(ystate),

     &     DIMS(work),
     &     bc,lo,hi,dt,dx,n,use_minion,iconserv)
c
c     This subroutine computes edges states, right now it uses
c     a lot of memory, but there becomes a trade off between
c     simplicity-efficiency in the new way of computing states
c     and complexity in the old way.  By eliminating loops over
c     state components though, the new way uses much less memory.
c
      implicit none
      integer i,j,n
      integer lo(SDIM),hi(SDIM),bc(SDIM,2)
      integer is,js,ie,je
      REAL_T place_to_break
      REAL_T hx, hy, dt, dth, dthx, dthy
      REAL_T tr,ubar,vbar,stx,sty,fu,fv,dx(SDIM)
      REAL_T eps,eps_for_bc
      logical ltx,lty
      parameter( eps        = 1.0D-6 )
      parameter( eps_for_bc = 1.0D-10 )

      integer DIMDEC(s)
      integer DIMDEC(work)
      integer DIMDEC(uedge)
      integer DIMDEC(xstate)
      integer DIMDEC(vedge)
      integer DIMDEC(ystate)

      REAL_T s(DIMV(s))
      REAL_T stxlo(DIM1(s)),stxhi(DIM1(s)),slxscr(DIM1(s),4)
      REAL_T stylo(DIM2(s)),styhi(DIM2(s)),slyscr(DIM2(s),4)

      REAL_T uedge(DIMV(uedge)), xstate(DIMV(xstate))
      REAL_T vedge(DIMV(vedge)), ystate(DIMV(ystate))

      REAL_T xlo(DIMV(work)), xhi(DIMV(work))
      REAL_T ylo(DIMV(work)), yhi(DIMV(work))
      REAL_T  sx(DIMV(work))
      REAL_T  sy(DIMV(work))
      REAL_T tforces(DIMV(work))
      REAL_T    divu(DIMV(work))
      
      REAL_T spx,smx,spy,smy,st,denom
      REAL_T smin,smax

      REAL_T   sc(DIMV(work),4)
      REAL_T gamp(DIM1(work))
      REAL_T gamm(DIM1(work))
      REAL_T   xm(DIM1(work))
      REAL_T   ym(DIM1(work))
      REAL_T    c(DIM1(work),4)
      REAL_T dt3rd

      integer choose_hi(DIM1(work))

      integer use_minion, iconserv
      integer inc
      integer iup,isign,jup,jsign
      REAL_T vtrans,vmult,vdif,stem,vaddif

      dth  = half*dt
      dthx = half*dt / dx(1)
      dthy = half*dt / dx(2)
      hx   = dx(1)
      hy   = dx(2)
      is = lo(1)
      js = lo(2)
      ie = hi(1)
      je = hi(2)

      dt3rd = dt * third

      smax = s(is,js)
      smin = s(is,js)
      do j = js,je
        do i = is,ie
          smax = max(smax,s(i,j))
          smin = min(smin,s(i,j))
        enddo
      enddo

      call bdsslope(s,lo(1),lo(2),hi(1),hi(2),sx,sy,sc,dx)

      do j = js,je 

        do i = is-1,ie 
          choose_hi(i+1) = merge(0,1,uedge(i+1,j) .ge. 0.0d0)
        end do 
        if (bc(1,1).eq.FOEXTRAP.or.bc(1,1).eq.HOEXTRAP
     &          .or.bc(1,1).eq.REFLECT_EVEN) then
           choose_hi(is) = 1
        endif
        if (bc(1,2).eq.FOEXTRAP.or.bc(1,2).eq.HOEXTRAP
     &          .or.bc(1,2).eq.REFLECT_EVEN) then
           choose_hi(ie+1) = 0
        endif

        do i = is-1,ie 
          iup   = merge(i,i+1,choose_hi(i+1).eq.0)
          isign = merge(1,-1, choose_hi(i+1).eq.0)
          vtrans = vedge(iup,j+1)
          jup   = merge(j,j+1,vtrans .ge. 0.0d0)
          jsign = merge(1,-1,vtrans .ge. 0.0d0)
          vmult = merge(uedge(i+1,j),half*(uedge(i+1,j)+uedge(i+1,j+1)),vtrans .ge. 0.0d0)
          xm(i) = isign*half*hx - two3rd*dt*vmult
          xm(i) = merge(xm(i), isign*min(isign*xm(i),hx*half), vtrans .ge. 0.0d0)
          ym(i) = jsign*half*hy - dt3rd*vedge(iup,j+1)

          c(i,1) = sc(iup,jup,1)
          c(i,2) = sc(iup,jup,2)
          c(i,3) = sc(iup,jup,3)
          c(i,4) = sc(iup,jup,4)

        enddo

        call bilin(gamp,c,xm,ym,hx,hy,is-1,ie,lo(1),hi(1))

c       Impose BCs on gamp.
        if (j .eq. je) then
          do i = is-1, ie  
            if (bc(2,2).eq.EXT_DIR .and.(vedge(i,je+1).le.zero) ) then
               gamp(i) = s(i,je+1)
            else if (bc(2,2).eq.REFLECT_ODD) then
               gamp(i) = zero
            end if
          end do
        end if

        do i = is-1,ie 
          iup   = merge(i,i+1,choose_hi(i+1).eq.0)
          isign = merge(1,-1, choose_hi(i+1).eq.0)
          vtrans = vedge(iup,j)
          jup   = merge(j-1,j,vtrans .ge. 0.0d0)
          jsign = merge(1,-1,vtrans .ge. 0.0d0)
          vmult = merge(half*(uedge(i+1,j)+uedge(i+1,j-1)),uedge(i+1,j),vtrans .ge. 0.0d0)
          xm(i) = isign*half*hx - two3rd*dt*vmult
          xm(i) = merge(isign*min(isign*xm(i),hx*half), xm(i), vtrans .ge. 0.0d0)
          ym(i) = jsign*half*hy - dt3rd*vedge(iup,j)

          c(i,1) = sc(iup,jup,1)
          c(i,2) = sc(iup,jup,2)
          c(i,3) = sc(iup,jup,3)
          c(i,4) = sc(iup,jup,4)

        enddo

        call bilin(gamm,c,xm,ym,hx,hy,is-1,ie,lo(1),hi(1))

c       Impose BCs on gamm.
        if (j .eq. js) then
          do i = is-1, ie 
            if (bc(2,1).eq.EXT_DIR .and.(vedge(i,js).ge.zero) ) then
               gamm(i) = s(i,js-1)
            else if (bc(2,1).eq.REFLECT_ODD) then
               gamm(i) = zero
            end if
          end do
        end if

        do i = is-1, ie 
          iup   = merge(i,i+1,choose_hi(i+1).eq.0)
          isign = merge(1,-1, choose_hi(i+1).eq.0)
          vdif   = dth*(vedge(iup,j+1)*gamp(i) - vedge(iup,j  )*gamm(i) ) / hy
          stem   = s(iup,j) + (isign*hx - uedge(i+1,j)*dt)*half*sx(iup,j)

c         vaddif = dth*stem*(uedge(iup+1,j) - uedge(iup,j))/hx
c         vaddif = dth*stem*(divu(iup,j) - (vedge(iup,j+1)-vedge(iup,j))/hy)
c         Change per JBB - dont use stem here!
          vaddif = dth*s(iup,j)*(divu(iup,j) - (vedge(iup,j+1)-vedge(iup,j))/hy)

          xstate(i+1,j) = stem - vdif - vaddif + dth*tforces(iup,j)

          if (iconserv .eq. 0) 
     $      xstate(i+1,j) = xstate(i+1,j) + dth*stem*divu(iup,j)

        enddo
      enddo

      do j = js-1,je 

        do i = is,ie 
          choose_hi(i) = merge(0,1,vedge(i,j+1) .ge. 0.0d0)
          if ( (j .eq. js-1) .and. 
     &       (bc(2,1).eq.FOEXTRAP.or.bc(2,1).eq.HOEXTRAP
     &                           .or.bc(2,1).eq.REFLECT_EVEN) )
     &     choose_hi(i) = 1
          if ( (j .eq. je) .and. 
     &       (bc(2,2).eq.FOEXTRAP.or.bc(2,2).eq.HOEXTRAP
     &                           .or.bc(2,2).eq.REFLECT_EVEN) )
     &     choose_hi(i) = 0
        end do 

        do i = is,ie 
          jup   = merge(j,j+1,choose_hi(i).eq.0)
          jsign = merge(1, -1,choose_hi(i).eq.0)
          vtrans = uedge(i+1,jup)
          iup   = merge(i,i+1,vtrans .ge. 0.0d0)
          isign = merge(1,-1,vtrans .ge. 0.0d0)
          vmult = merge(vedge(i,j+1),half*(vedge(i,j+1)+vedge(i+1,j+1)),vtrans .ge. 0.0d0)
          xm(i) = isign*half*hx - dt3rd*uedge(i+1,jup)
          ym(i) = jsign*half*hy - two3rd*dt*vmult
          ym(i) = merge(ym(i),jsign*min(jsign*ym(i),hy*half),vtrans .ge. 0.0d0)

          c(i,1) = sc(iup,jup,1)
          c(i,2) = sc(iup,jup,2)
          c(i,3) = sc(iup,jup,3)
          c(i,4) = sc(iup,jup,4)
        enddo

        call bilin(gamp,c,xm,ym,hx,hy,is,ie,lo(1),hi(1))

c       Impose BCs on gamp.
        if (bc(1,2).eq.EXT_DIR .and.(uedge(ie+1,j).le.zero) ) then
           gamp(ie) = s(ie+1,j)
        else if (bc(1,2).eq.REFLECT_ODD) then
           gamp(ie) = zero
        end if

        do i = is,ie 
          jup   = merge(j,j+1,choose_hi(i).eq.0)
          jsign = merge(1, -1,choose_hi(i).eq.0)
          vtrans = uedge(i,jup)
          iup   = merge(i-1,i,vtrans .ge. 0.0d0)
          isign = merge(1,-1,vtrans .ge. 0.0d0)
          vmult = merge(half*(vedge(i,j+1)+vedge(i-1,j+1)),vedge(i,j+1),vtrans .ge. 0.0d0)
          xm(i) = isign*half*hx - dt3rd*uedge(i,jup)
          ym(i) = jsign*half*hy - two3rd*dt*vmult
          ym(i) = merge(jsign*min(jsign*ym(i),hy*half),ym(i),vtrans .ge. 0.0d0)

          c(i,1) = sc(iup,jup,1)
          c(i,2) = sc(iup,jup,2)
          c(i,3) = sc(iup,jup,3)
          c(i,4) = sc(iup,jup,4)
        enddo

        call bilin(gamm,c,xm,ym,hx,hy,is,ie,lo(1),hi(1))

c       Impose BCs on gamm.
        if (bc(1,1).eq.EXT_DIR .and.(uedge(is,j).ge.zero) ) then
           gamm(is) = s(is,j)
        else if (bc(1,1).eq.REFLECT_ODD) then
           gamm(is) = zero
        end if

        do i = is,ie 
          jup   = merge(j,j+1,choose_hi(i).eq.0)
          jsign = merge(1, -1,choose_hi(i).eq.0)
          vdif   = dth*(uedge(i+1,jup)*gamp(i)-uedge(i,jup)*gamm(i))/hx
          stem   = s(i,jup) + (jsign*hy - vedge(i,j+1)*dt)*half*sy(i,jup)

c         vaddif = dth*stem*(vedge(i,jup+1) - vedge(i,jup))/hy
c         vaddif = dth*stem*(divu(i,jup) - (uedge(i+1,jup)-uedge(i,jup))/hx)
c         Change per JBB - dont use stem here!
          vaddif = dth*s(i,jup)*(divu(i,jup) - (uedge(i+1,jup)-uedge(i,jup))/hx)

          ystate(i,j+1) = stem - vdif - vaddif + dth*tforces(i,jup)

          if (iconserv .eq. 0) 
     $      ystate(i,j+1) = ystate(i,j+1) + dth*stem*divu(i,jup)

        enddo
      enddo

c     
c     Impose BCs on xstate.
c     
      do j = js,je 
        if (bc(1,1).eq.EXT_DIR .and. uedge(is,j).ge.zero) then
           xstate(is,j) = s(is-1,j)
        else if (bc(1,1).eq.REFLECT_ODD) then
           xstate(is,j) = zero
        end if

        if (bc(1,2).eq.EXT_DIR .and. uedge(ie+1,j).le.zero) then
           xstate(ie+1,j) = s(ie+1,j)
        else if (bc(1,2).eq.REFLECT_ODD) then
           xstate(ie+1,j) = zero
        end if
      enddo

c     
c     Impose BCs on ystate.
c     
      do i = is,ie 
        if (bc(2,1).eq.EXT_DIR .and.(vedge(i,js).ge.zero) ) then
           ystate(i,js) = s(i,js-1)
        else if (bc(2,1).eq.REFLECT_ODD) then
            ystate(i,js) = zero
        end if

        if (bc(2,2).eq.EXT_DIR .and. (vedge(i,je+1).le.zero) ) then
           ystate(i,je+1) = s(i,je+1)
        else if (bc(2,2).eq.REFLECT_ODD) then
            ystate(i,je+1) = zero
        end if
      enddo

      end

      subroutine bilin(poly,c,x,y,hx,hy,istart,iend,lo_1,hi_1)

      implicit none
      integer lo_1,hi_1
      integer istart,iend
      REAL_T poly(lo_1-1:hi_1+1)
      REAL_T    x(lo_1-1:hi_1+1)
      REAL_T    y(lo_1-1:hi_1+1)
      REAL_T    c(lo_1-1:hi_1+1,4)
      REAL_T hx, hy
      REAL_T centx, centy
      REAL_T l1,l2,l3,l4

      integer i

      centx = hx/2
      centy = hy/2

      do i = istart,iend

        l1 = (-centx - x(i))*(-centy - y(i))
        l2 = (-centx - x(i))*( centy - y(i))
        l3 = ( centx - x(i))*(-centy - y(i))
        l4 = ( centx - x(i))*( centy - y(i))

        poly(i) = l2*c(i,3) + l3*c(i,2) - l1*c(i,4) - l4*c(i,1)
        poly(i) = -poly(i)/(hx*hy)
      enddo

      end

