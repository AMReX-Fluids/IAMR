c
c $Id: DIFFUSION_2D.F,v 1.7 1998-11-04 23:33:11 lijewski Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "DIFFUSION_F.H"
#include "ArrayLim.H"

#define SDIM 2

c ::
c :: ----------------------------------------------------------
c :: compute the time averaged viscous flux at the given edge
c ::  for variable transport coefficients
c ::
c :: NOTE:
c ::   These are only computed for a reflux operation, and require
c ::   that the ghost cells have been filled according to the 
c ::   diffusion operator (maybe this isn't such a good idea...)
c :: ----------------------------------------------------------
c ::

      subroutine FORT_VISCFLUX_VC (s_o,s_n,DIMS(s),lo,hi,
     $                             flux, DIMS(flux), area, DIMS(area),
     $                             bn, bnp1, DIMS(b),
     $                             dx,mult,theta,dir)

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(s)
      integer DIMDEC(flux)
      integer DIMDEC(area)
      integer DIMDEC(b)
      REAL_T  s_o(DIMV(s))
      REAL_T  s_n(DIMV(s))
      REAL_T  flux(DIMV(flux))
      REAL_T  area(DIMV(area))
      REAL_T  bn(DIMV(b))
      REAL_T  bnp1(DIMV(b))
      REAL_T  dx
      REAL_T  mult,theta
      integer dir

      integer i, j
      REAL_T  sx1, sx2, sy1, sy2
      REAL_T  onemintheta

      onemintheta = one-theta
      if (dir .EQ. 0) then
c     ::::: compute X fluxes
         do j = lo(2), hi(2)
         do i = lo(1), hi(1)+1
            sx1 = (s_o(i,j) - s_o(i-1,j))*bn(i,j)
            sx2 = (s_n(i,j) - s_n(i-1,j))*bnp1(i,j)
            flux(i,j) = mult*(onemintheta*sx1 + theta*sx2)*area(i,j)/dx
         end do
         end do
      else
c     ::::: compute Y fluxes
         do j = lo(2), hi(2)+1
         do i = lo(1), hi(1)
            sy1 = (s_o(i,j) - s_o(i,j-1))*bn(i,j)
            sy2 = (s_n(i,j) - s_n(i,j-1))*bnp1(i,j)
            flux(i,j) = mult*(onemintheta*sy1 + theta*sy2)*area(i,j)/dx
         end do
         end do
      end if
      end

c ::
c :: ----------------------------------------------------------
c :: compute the time averaged viscous flux at the given edge
c ::  for constant transport coefficients
c ::
c :: NOTE:
c ::   These are only computed for a reflux operation, and require
c ::   that the ghost cells have been filled according to the 
c ::   diffusion operator (maybe this isn't such a good idea...)
c :: ----------------------------------------------------------
c ::

      subroutine FORT_VISCFLUX_CC (s_o,s_n,DIMS(s),lo,hi,
     $                             flux, DIMS(flux), area, DIMS(area),
     $                             dx,mult,theta,dir)

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(s)
      integer DIMDEC(flux)
      integer DIMDEC(area)
      REAL_T  s_o(DIMV(s))
      REAL_T  s_n(DIMV(s))
      REAL_T  flux(DIMV(flux))
      REAL_T  area(DIMV(area))
      REAL_T  dx
      REAL_T  mult,theta
      integer dir

      integer i, j
      REAL_T  sx1, sx2, sy1, sy2
      REAL_T  onemintheta

      onemintheta = one-theta
      if (dir .EQ. 0) then
c     ::::: compute X fluxes
         do j = lo(2), hi(2)
         do i = lo(1), hi(1)+1
            sx1 = s_o(i,j) - s_o(i-1,j)
            sx2 = s_n(i,j) - s_n(i-1,j)
            flux(i,j) = mult*(onemintheta*sx1 + theta*sx2)*area(i,j)/dx
         end do
         end do
      else
c     ::::: compute Y fluxes
         do j = lo(2), hi(2)+1
         do i = lo(1), hi(1)
            sy1 = s_o(i,j) - s_o(i,j-1)
            sy2 = s_n(i,j) - s_n(i,j-1)
            flux(i,j) = mult*(onemintheta*sy1 + theta*sy2)*area(i,j)/dx
         end do
         end do
      end if
      end

      subroutine FORT_VISCSYNCFLUX (ssync,DIMS(ssync),lo,hi,
     $                              xflux,DIMS(xf),yflux,DIMS(yf),
     $                              xarea,DIMS(ax),yarea,DIMS(ay),dx,mult)

      integer lo(2), hi(2)
      integer DIMDEC(ssync)
      integer DIMDEC(xf)
      integer DIMDEC(yf)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      REAL_T  ssync(DIMV(ssync))
      REAL_T  xflux(DIMV(xf))
      REAL_T  yflux(DIMV(yf))
      REAL_T  xarea(DIMV(ax))
      REAL_T  yarea(DIMV(ay))
      REAL_T  dx(2)
      REAL_T  mult

      REAL_T  sx
      REAL_T  sy
      integer i, j
c
c     ::::: compute X fluxes
c
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)+1
	       sx = ssync(i,j) - ssync(i-1,j)
               xflux(i,j) = mult*sx*xarea(i,j)/dx(1)
            end do
         end do
c
c     ::::: compute Y fluxes
c
         do j = lo(2), hi(2)+1
            do i = lo(1), hi(1)
	       sy = ssync(i,j) - ssync(i,j-1)
               yflux(i,j) = mult*sy*yarea(i,j)/dx(2)
	    end do
         end do

      end

c :: ----------------------------------------------------------
c :: HOOPSRC
c ::             fab(i,j) = fab(i,j) - mu*u/(r(i)^2)
c ::
c :: INPUTS / OUTPUTS:
c ::  fab       <=  array to be modified
c ::  DIMS(fab)  => index limits of fab
c ::  mu         => viscous coefficient
c :: ----------------------------------------------------------
c ::
       subroutine FORT_HOOPSRC (DIMS(grid), fab, DIMS(fab), u, DIMS(u), r, mu)

       integer DIMDEC(grid)
       integer DIMDEC(fab)
       integer DIMDEC(u)
       REAL_T  fab(DIMV(fab))
       REAL_T  u(DIMV(u))
       REAL_T  r(DIM1(grid))
       REAL_T  mu

       integer i, j

c      if (ARG_L1(u) .lt. ARG_L1(fab) .or. ARG_H1(u) .gt. ARG_H1(fab)) then
c         write(6,*) "FORT_HOOPSRC: bad index limits"
c         stop
c      end if

       do j = ARG_L2(grid), ARG_H2(grid)
          do i = ARG_L1(grid), ARG_H1(grid)
             fab(i,j) = fab(i,j) - mu*u(i,j)/(r(i)*r(i))
          end do
       end do

       end

c :: ----------------------------------------------------------
c :: HOOPRHS
c ::             rhs(i,j) = rhs(i,j) - (one-theta)*dt*u*mu*vol/(r(i)^2)
c ::
c :: INPUTS / OUTPUTS:
c ::  fab       <=  array to be modified
c ::  DIMS(fab)  => index limits of fab
c ::  u          => array to be modified
c ::  DIMS(u)    => index limits of u
c ::  r          => 1-D r array (in first coordinate direction)
c ::  mu         => scalar viscosity
c ::  dt         => time step
c ::  vol        => volume array
c ::  DIMS(vol)  => index limits of vol
c ::  b          => (one-theta)*dt
c :: ----------------------------------------------------------
c ::
       subroutine FORT_HOOPRHS (fab, DIMS(fab), u, DIMS(u), r, b,
     &                          vol, DIMS(vol))
       integer DIMDEC(fab)
       integer DIMDEC(u)
       integer DIMDEC(vol)
       REAL_T  fab(DIMV(fab))
       REAL_T  u(DIMV(u))
       REAL_T  vol(DIMV(vol))
       REAL_T  r(DIM1(fab))
       REAL_T  b

       integer i, j

       do j = ARG_L2(fab), ARG_H2(fab)
          do i = ARG_L1(fab), ARG_H1(fab)
             fab(i,j) = fab(i,j) - b*vol(i,j)*u(i,j)/(r(i)*r(i))
          end do
       end do

       end

c :: ----------------------------------------------------------
c :: TENSOR_HOOPRHS
c ::             rhs(i,j) = rhs(i,j) - (1-theta)*dt*u*two*mu_cen*vol/(r(i)^2)
c ::                                                  ^^^ yes, that is correct
c ::                                                      for variable mu
c ::
c :: INPUTS / OUTPUTS:
c ::  fab       <=  array to be modified
c ::  DIMS(fab)  => index limits of fab
c ::  u          => array to be modified
c ::  DIMS(u)    => index limits of u
c ::  r          => 1-D r array (in first coordinate direction)
c ::  mu         => scalar viscosity
c ::  dt         => time step
c ::  vol        => volume array
c ::  DIMS(vol)  => index limits of vol
c ::  b          => (1-theta)*dt
c :: ----------------------------------------------------------
c ::
       subroutine FORT_TENSOR_HOOPRHS (xvelcomp, fab, DIMS(fab), u, DIMS(u), 
     &                          r, b,
     &                          vol, DIMS(vol), betax, DIMS(betax),
     &                          betay, DIMS(betay))
       integer xvelcomp
       integer DIMDEC(fab)
       integer DIMDEC(u)
       integer DIMDEC(vol)
       integer DIMDEC(betax)
       integer DIMDEC(betay)
       REAL_T  fab(DIMV(fab),2)
       REAL_T  u(DIMV(u),2)
       REAL_T  vol(DIMV(vol))
       REAL_T  r(DIM1(fab))
       REAL_T  betax(DIMV(betax))
       REAL_T  betay(DIMV(betay))
       REAL_T  b

       REAL_T  betacen
       integer i, j

       do j = ARG_L2(fab), ARG_H2(fab)
          do i = ARG_L1(fab), ARG_H1(fab)
             betacen = fourth*(betax(i,j)+betax(i+1,j)+
     &                           betay(i,j)+betay(i,j+1))
             fab(i,j,xvelcomp) = fab(i,j,xvelcomp) - 
     &           b*two*betacen*vol(i,j)*u(i,j,xvelcomp)/(r(i)*r(i))
          end do
       end do

       end

c :: ----------------------------------------------------------
c :: TENSOR_HOOPSRC
c ::             fab(i,j) = fab(i,j) - two*mu*u/(r(i)^2)
c ::                                   ^^^ yes, that is correct
c ::                                       for variable mu
c ::
c :: INPUTS / OUTPUTS:
c ::  fab       <=  array to be modified
c ::  DIMS(fab)  => index limits of fab
c ::  mu         => viscous coefficient
c :: ----------------------------------------------------------
c ::
       subroutine FORT_TENSOR_HOOPSRC (comp, DIMS(grid), fab, DIMS(fab), u, 
     &      DIMS(u), r, betax, DIMS(betax), betay, DIMS(betay))

       integer comp
       integer DIMDEC(grid)
       integer DIMDEC(fab)
       integer DIMDEC(u)
       integer DIMDEC(betax)
       integer DIMDEC(betay)
       REAL_T  fab(DIMV(fab),2)
       REAL_T  u(DIMV(u),2)
       REAL_T  r(DIM1(grid))
       REAL_T  betax(DIMV(betax))
       REAL_T  betay(DIMV(betay))

       integer i, j
       REAL_T  betacen

       do j = ARG_L2(grid), ARG_H2(grid)
          do i = ARG_L1(grid), ARG_H1(grid)
             betacen  = fourth*(betax(i,j)+betax(i+1,j)+
     &                           betay(i,j)+betay(i,j+1))
             fab(i,j,comp) = fab(i,j,comp) - two*betacen*u(i,j,comp)/(r(i)*r(i))
          end do
       end do

       end

c :: ----------------------------------------------------------
c :: SETALPHA
c ::             alpha(i,j) = vol*(1+b/(r(i)^2)) / density
c ::
c :: INPUTS / OUTPUTS:
c ::  fab       <=  array to be modified
c ::  DIMS(fab) => index limits of fab
c ::  lo,hi     => index limits of box
c ::  r         =>  1-d array of radius
c ::  b         =>  either theta*dt*mu or -(1-theta)*dt*mu
c ::  vol       =>  volume array
c ::  DIMS(vol) => index limits of fab
c ::  denfab    => array of density at time n+1/2
c ::  DIMS(den) => index limits of fab
c ::  usehoop   => do we add hoop stress?   (only if x-vel component)
c ::  useden    => do we divide by density? (only if velocity component)
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SETALPHA (fab, DIMS(fab), lo, hi, r, DIMS(r),
     $                           b, vol, DIMS(vol),
     &                           denfab,DIMS(den),usehoop,useden)

       integer DIMDEC(fab)
       integer lo(SDIM), hi(SDIM)
       integer DIMDEC(vol)
       integer DIMDEC(den)
       integer DIMDEC(r)
       REAL_T  fab(DIMV(fab))
       REAL_T  vol(DIMV(vol))
       REAL_T  denfab(DIMV(den))
       REAL_T  r(DIM1(r))
       REAL_T  b
       integer usehoop,useden

       integer i, j

       if (usehoop .eq. 0) then
          if (useden .eq. 0) then
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   fab(i,j) = vol(i,j)
                end do
             end do
          else 
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   fab(i,j) = vol(i,j) * denfab(i,j)
                end do
             end do
          end if
       else
          if (useden .eq. 0) then
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   fab(i,j) = vol(i,j) * (one + (b / (r(i)*r(i))))
                end do
             end do
          else
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   fab(i,j) = vol(i,j) * (denfab(i,j) + (b / (r(i)*r(i))))
                end do
             end do
          end if
       end if

       end

c :: ----------------------------------------------------------
c :: SET_TENSOR_ALPHA
c ::             alpha(i,j) = vol*density+b*dr*dz*two*mu_cen/r(i)
c ::                        = vol*(density+b*two*mu_cen*r(i)**2)
c ::                                         ^^^ yes, that is correct
c ::                                             for variable mu
c ::
c :: INPUTS / OUTPUTS:
c ::  fab       <=  array to be modified
c ::  DIMS(fab) => index limits of fab
c ::  lo,hi     => index limits of box
c ::  r         =>  1-d array of radius
c ::  b         =>  theta*dt or -(1-theta)*dt
c ::  vol       =>  volume array
c ::  DIMS(vol) => index limits of fab
c ::  denfab    => array of density at time n+1/2
c ::  DIMS(den) => index limits of fab
c ::  usehoop   => do we add hoop stress?   (only if x-vel component)
c ::  useden    => do we divide by density? (only if velocity component)
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SET_TENSOR_ALPHA (alpha, DIMS(alpha), lo, hi, r, DIMS(r),
     $                           b, vol, DIMS(vol),
     &                           denfab,DIMS(den),betax,DIMS(betax),
     &                           betay,DIMS(betay),isrz)

       integer DIMDEC(alpha)
       integer lo(SDIM), hi(SDIM)
       integer DIMDEC(vol)
       integer DIMDEC(den)
       integer DIMDEC(betax)
       integer DIMDEC(betay)
       integer DIMDEC(r)
       REAL_T  alpha(DIMV(alpha),2)
       REAL_T  vol(DIMV(vol))
       REAL_T  denfab(DIMV(den))
       REAL_T  betax(DIMV(betax))
       REAL_T  betay(DIMV(betay))
       REAL_T  r(DIM1(r))
       REAL_T  b, betacen
       integer isrz

       integer i, j

       if (isrz .eq. 0) then
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                alpha(i,j,1) = vol(i,j) * denfab(i,j)
                alpha(i,j,2) = vol(i,j) * denfab(i,j)
             end do
          end do
       else
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                alpha(i,j,2) = vol(i,j) * denfab(i,j)
                betacen = fourth*(betax(i,j)+betax(i+1,j)+
     &               betay(i,j)+betay(i,j+1))
                alpha(i,j,1) = vol(i,j) * (denfab(i,j)+
     &               b*two*betacen/(r(i)**2))
             end do
          end do
       end if

       end

      subroutine FORT_DIV_MU_SI(lo, hi, dx, mu, DIMS(divu), divu,
     &     DIMS(divmusi), divmusi)

      implicit none
c
c ... inputs
c
      integer lo(SDIM), hi(SDIM)
      REAL_T  dx(SDIM)
      integer DIMDEC(divu)
      REAL_T  divu(DIMV(divu))      
      REAL_T  mu
c
c ... outputs
c
      integer DIMDEC(divmusi)
      REAL_T  divmusi(DIMV(divmusi),SDIM)
c
c ... local 
c
      integer i,j
      REAL_T sleft, sright, stop, sbot
c
c ... Note: the following IS correct for r-z. Terms from the hoop stress
c           cancel with terms from tau_rr to eliminate all r dependence.
c
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            sleft = half*(divu(i-1,j)+divu(i,j))
            sright = half*(divu(i+1,j)+divu(i,j))

            divmusi(i,j,1) = mu*(sright-sleft)/dx(1)

            stop = half*(divu(i,j)+divu(i,j+1))
            sbot = half*(divu(i,j-1)+divu(i,j))

            divmusi(i,j,2) = mu*(stop-sbot)/dx(2)

         end do
      end do

      end

      subroutine FORT_DIV_VARMU_SI(lo, hi, dx, DIMS(divu), divu,
     &     DIMS(betax), betax, DIMS(betay), betay, DIMS(divmusi), divmusi)

      implicit none
c
c ... inputs
c
      integer lo(SDIM), hi(SDIM)
      REAL_T  dx(SDIM)
      integer DIMDEC(divu)
      REAL_T  divu(DIMV(divu))      
      integer DIMDEC(betax)
      REAL_T  betax(DIMV(betax))
      integer DIMDEC(betay)
      REAL_T  betay(DIMV(betay))
c
c ... outputs
c
      integer DIMDEC(divmusi)
      REAL_T  divmusi(DIMV(divmusi),SDIM)
c
c ... local 
c
      integer i,j
      REAL_T sleft, sright, stop, sbot
c
c ... Note: the following IS correct for r-z. Terms from the hoop stress
c           cancel with terms from tau_rr to eliminate all r dependence.
c
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            sleft = half*(divu(i-1,j)+divu(i,j))
            sright = half*(divu(i+1,j)+divu(i,j))

            divmusi(i,j,1) = (betax(i+1,j)*sright-
     &           betax(i,j)*sleft)/dx(1)

            stop = half*(divu(i,j)+divu(i,j+1))
            sbot = half*(divu(i,j-1)+divu(i,j))

            divmusi(i,j,2) = (betay(i,j+1)*stop-
     &           betay(i,j)*sbot)/dx(2)

         end do
      end do

      end

      subroutine FORT_CCTOEDGE(lo, hi, DIMS(beta), beta,
     &   DIMS(betax), betax, DIMS(betay), betay)

      implicit none
c
c ... inputs
c
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(beta)
      REAL_T  beta(DIMV(beta))
      integer DIMDEC(betax)
      integer DIMDEC(betay)
c
c ... outputs
c
      REAL_T  betax(DIMV(betax))
      REAL_T  betay(DIMV(betay))
c
c ... local
c
      integer i,j

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            betax(i,j) = .5*(beta(i-1,j)+beta(i,j))
         end do
      end do

      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            betay(i,j) = .5*(beta(i,j-1)+beta(i,j))
         end do
      end do

      end

