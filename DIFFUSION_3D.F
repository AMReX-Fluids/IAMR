c
c $Id: DIFFUSION_3D.F,v 1.10 2002-08-29 22:14:52 car Exp $
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

#define SDIM 3

c ::
c :: ----------------------------------------------------------
c :: compute the time averaged viscous flux at the given edge
c ::  for variable transport coefficients
c ::
c :: NOTE:
c ::   These are only computed for a reflux operation, and require
c ::   that the ghost cells have been filled according to the 
c ::   diffusion operator (maybe this isnt such a good idea...)
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

      integer i, j, k
      REAL_T  sx1, sx2, sy1, sy2, sz1, sz2
      REAL_T  onemintheta

      onemintheta = one-theta
      if (dir .EQ. 0) then
c     ::::: compute X fluxes
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
         do i = lo(1), hi(1)+1
	    sx1 = (s_o(i,j,k) - s_o(i-1,j,k))*bn(i,j,k)
	    sx2 = (s_n(i,j,k) - s_n(i-1,j,k))*bnp1(i,j,k)
	    flux(i,j,k) = mult*(onemintheta*sx1 + theta*sx2)*area(i,j,k)/dx
         end do
         end do
         end do
      else if (dir .EQ. 1) then
c     ::::: compute Y fluxes
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)+1
         do i = lo(1), hi(1)
            sy1 = (s_o(i,j,k) - s_o(i,j-1,k))*bn(i,j,k)
            sy2 = (s_n(i,j,k) - s_n(i,j-1,k))*bnp1(i,j,k)
            flux(i,j,k) = mult*(onemintheta*sy1 + theta*sy2)*area(i,j,k)/dx
	 end do
         end do
         end do
      else if (dir .EQ. 2) then
c     ::::: compute Z fluxes
         do k = lo(3), hi(3)+1
         do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            sz1 = (s_o(i,j,k) - s_o(i,j,k-1))*bn(i,j,k)
            sz2 = (s_n(i,j,k) - s_n(i,j,k-1))*bnp1(i,j,k)
            flux(i,j,k) = mult*(onemintheta*sz1 + theta*sz2)*area(i,j,k)/dx
	 end do
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
c ::   diffusion operator (maybe this isnt such a good idea...)
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

      integer i, j, k
      REAL_T  sx1, sx2, sy1, sy2, sz1, sz2
      REAL_T  onemintheta

      onemintheta = one-theta
      if (dir .EQ. 0) then
c     ::::: compute X fluxes
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
         do i = lo(1), hi(1)+1
            sx1 = s_o(i,j,k) - s_o(i-1,j,k)
            sx2 = s_n(i,j,k) - s_n(i-1,j,k)
            flux(i,j,k) = mult*(onemintheta*sx1 + theta*sx2)*area(i,j,k)/dx
         end do
         end do
         end do
      else if (dir .EQ. 1) then
c     ::::: compute Y fluxes
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)+1
         do i = lo(1), hi(1)
            sy1 = s_o(i,j,k) - s_o(i,j-1,k)
            sy2 = s_n(i,j,k) - s_n(i,j-1,k)
            flux(i,j,k) = mult*(onemintheta*sy1 + theta*sy2)*area(i,j,k)/dx
         end do
         end do
         end do
      else
c     ::::: compute Z fluxes
         do k = lo(3), hi(3)+1
         do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            sz1 = s_o(i,j,k) - s_o(i,j,k-1)
            sz2 = s_n(i,j,k) - s_n(i,j,k-1)
            flux(i,j,k) = mult*(onemintheta*sz1 + theta*sz2)*area(i,j,k)/dx
	 end do
         end do
         end do
      end if
      end

      subroutine FORT_VISCSYNCFLUX (ssync,DIMS(ssync),lo,hi,
     $                              xflux,DIMS(xf),yflux,DIMS(yf),
     $                              zflux,DIMS(zf),
     $                              xarea,DIMS(ax),yarea,DIMS(ay),
     $                              zarea,DIMS(az),dx,mult)

      integer lo(3), hi(3)
      integer DIMDEC(ssync)
      integer DIMDEC(xf)
      integer DIMDEC(yf)
      integer DIMDEC(zf)
      integer DIMDEC(ax)
      integer DIMDEC(ay)
      integer DIMDEC(az)
      REAL_T  ssync(DIMV(ssync))
      REAL_T  xflux(DIMV(xf))
      REAL_T  yflux(DIMV(yf))
      REAL_T  zflux(DIMV(zf))
      REAL_T  xarea(DIMV(ax))
      REAL_T  yarea(DIMV(ay))
      REAL_T  zarea(DIMV(az))
      REAL_T  dx(3)
      REAL_T  mult

      integer i, j, k, n
      REAL_T  sx
      REAL_T  sy
      REAL_T  sz
c
c     ::::: compute X fluxes
c
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)+1
               sx = ssync(i,j,k) - ssync(i-1,j,k)
               xflux(i,j,k) = mult*sx*xarea(i,j,k)/dx(1)
            end do
         end do
      end do
c
c     ::::: compute Y fluxes
c
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)+1
            do i = lo(1), hi(1)
               sy = ssync(i,j,k) - ssync(i,j-1,k)
               yflux(i,j,k) = mult*sy*yarea(i,j,k)/dx(2)
            end do
         end do
      end do
c
c     ::::: compute Z fluxes
c
      do k = lo(3), hi(3)+1
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               sz = ssync(i,j,k) - ssync(i,j,k-1)
               zflux(i,j,k) = half*mult*sz*zarea(i,j,k)/dx(3)
            end do
         end do
      end do
      
      end

c :: ----------------------------------------------------------
c :: SETALPHA
c ::             alpha(i,j,k) = vol*(1+b/(r(i)^2)) / density
c ::
c :: INPUTS / OUTPUTS:
c ::  fab       <=  array to be modified
c ::  DIMS(fab) => index limits of fab
c ::  lo,hi     => index limits of box
c ::  r         =>  1-d array of radius
c ::  DIMS(r)   => index limits of r
c ::  b         =>  viscous coefficient
c ::  vol       =>  volume array
c ::  DIMS(vol) => index limits of fab
c ::  denfab    => array of density at time n+1/2
c ::  DIMS(den) => index limits of fab
c ::  usehoop   => do we add hoop stress?   NOT IN 3-D
c ::  useden    => do we divide by density? (only if velocity component)
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SETALPHA (fab, DIMS(fab), lo, hi, r, DIMS(r), b, 
     $                           vol, DIMS(vol), denfab, DIMS(den),
     &                           usehoop,useden)

       integer DIMDEC(fab)
       integer DIMDEC(r)
       integer DIMDEC(vol)
       integer DIMDEC(den)
       integer lo(SDIM), hi(SDIM)
       integer usehoop,useden
       REAL_T  fab(DIMV(fab))
       REAL_T  vol(DIMV(vol))
       REAL_T  denfab(DIMV(den))
       REAL_T  r(DIM1(r))
       REAL_T  b

       integer i, j, k

       if (useden .eq. 0) then

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   fab(i,j,k) = vol(i,j,k)
                end do
             end do
          end do

       else 

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   fab(i,j,k) = vol(i,j,k) * denfab(i,j,k)
                end do
             end do
          end do

       end if

       end

c :: ----------------------------------------------------------
c :: SET_TENSOR_ALPHA
c ::             alpha(i,j) = vol*density
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
     &                           betay,DIMS(betay),betaz,DIMS(betaz),isrz)

       integer DIMDEC(alpha)
       integer lo(SDIM), hi(SDIM)
       integer DIMDEC(vol)
       integer DIMDEC(den)
       integer DIMDEC(betax)
       integer DIMDEC(betay)
       integer DIMDEC(betaz)
       integer DIMDEC(r)
       REAL_T  alpha(DIMV(alpha),1)
       REAL_T  vol(DIMV(vol))
       REAL_T  denfab(DIMV(den))
       REAL_T  betax(DIMV(betax))
       REAL_T  betay(DIMV(betay))
       REAL_T  betaz(DIMV(betaz))
       REAL_T  r(DIM1(r))
       REAL_T  b, betacen
       integer isrz

       integer i, j, k

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                alpha(i,j,k,1) = vol(i,j,k) * denfab(i,j,k)
             end do
          end do
       end do

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
      integer i,j,k
      REAL_T sleft, sright, stop, sbot, sfront, sback
c
c ... Note: the following IS correct for r-z. Terms from the hoop stress
c           cancel with terms from tau_rr to eliminate all r dependence.
c
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               sleft = half*(divu(i-1,j,k)+divu(i,j,k))
               sright = half*(divu(i+1,j,k)+divu(i,j,k))

               divmusi(i,j,k,1) = mu*(sright-sleft)/dx(1)

               stop = half*(divu(i,j,k)+divu(i,j+1,k))
               sbot = half*(divu(i,j-1,k)+divu(i,j,k))

               divmusi(i,j,k,2) = mu*(stop-sbot)/dx(2)

               sfront = half*(divu(i,j,k)+divu(i,j,k+1))
               sback  = half*(divu(i,j,k-1)+divu(i,j,k))

               divmusi(i,j,k,3) = mu*(sfront-sback)/dx(2)
            end do
         end do
      end do

      end

      subroutine FORT_DIV_VARMU_SI(lo, hi, dx, DIMS(divu), divu,
     &     DIMS(betax), betax, DIMS(betay), betay,  DIMS(betaz), 
     &     betaz, DIMS(divmusi), divmusi)

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
      integer DIMDEC(betaz)
      REAL_T  betaz(DIMV(betaz))
c
c ... outputs
c
      integer DIMDEC(divmusi)
      REAL_T  divmusi(DIMV(divmusi),SDIM)
c
c ... local 
c
      integer i,j,k
      REAL_T sleft, sright, stop, sbot, sfront, sback

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               sleft = half*(divu(i-1,j,k)+divu(i,j,k))
               sright = half*(divu(i+1,j,k)+divu(i,j,k))

               divmusi(i,j,k,1) = (betax(i+1,j,k)*sright-
     &              betax(i,j,k)*sleft)/dx(1)

               stop = half*(divu(i,j,k)+divu(i,j+1,k))
               sbot = half*(divu(i,j-1,k)+divu(i,j,k))

               divmusi(i,j,k,2) = (betay(i,j+1,k)*stop-
     &              betay(i,j,k)*sbot)/dx(2)

               sfront = half*(divu(i,j,k)+divu(i,j,k+1))
               sback  = half*(divu(i,j,k-1)+divu(i,j,k))

               divmusi(i,j,k,3) = (betaz(i,j,k+1)*sfront-
     &              betaz(i,j,k)*sback)/dx(3)

            end do
         end do
      end do

      end

      subroutine FORT_CCTOEDGE(lo, hi, DIMS(beta), beta,
     &   DIMS(betax), betax, DIMS(betay), betay, DIMS(betaz), betaz)

      implicit none
c
c ... inputs
c
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(beta)
      REAL_T  beta(DIMV(beta))
      integer DIMDEC(betax)
      integer DIMDEC(betay)
      integer DIMDEC(betaz)
c
c ... outputs
c
      REAL_T  betax(DIMV(betax))
      REAL_T  betay(DIMV(betay))
      REAL_T  betaz(DIMV(betaz))
c
c ... local
c
      integer i,j,k

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
               betax(i,j,k) = .5D0*(beta(i-1,j,k)+beta(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
               betay(i,j,k) = .5D0*(beta(i,j-1,k)+beta(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               betaz(i,j,k) = .5D0*(beta(i,j,k-1)+beta(i,j,k))
            end do
         end do
      end do

      end


