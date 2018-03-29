
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <DIFFUSION_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 3

module diffusion_3d_module
  
  implicit none

  private

  public :: viscsyncflux, fort_setalpha, set_tensor_alpha, &
            div_mu_si, div_varmu_si

contains

      subroutine viscsyncflux (ssync,DIMS(ssync), &
                                   xlo,xhi,ylo,yhi,zlo,zhi, &
                                   xflux,DIMS(xf),yflux,DIMS(yf), &
                                   zflux,DIMS(zf), &
                                   xarea,DIMS(ax),yarea,DIMS(ay),&
                                   zarea,DIMS(az),dx,mult)&
                                   bind(C,name="viscsyncflux")

      implicit none
      integer xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3)
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

      integer i, j, k
      REAL_T  sx, sy, sz, idx, idy, idz

      idx = 1.0d0 / dx(1)
      idy = 1.0d0 / dx(2)
      idz = 1.0d0 / dx(3)

!c
!c     ::::: compute X fluxes
!c
      do       k = xlo(3), xhi(3)
         do    j = xlo(2), xhi(2)
            do i = xlo(1), xhi(1)
               sx = ssync(i,j,k) - ssync(i-1,j,k)
               xflux(i,j,k) = mult*sx*xarea(i,j,k)*idx
            end do
         end do
      end do
!c
!c     ::::: compute Y fluxes
!c
      do       k = ylo(3), yhi(3)
         do    j = ylo(2), yhi(2)
            do i = ylo(1), yhi(1)
               sy = ssync(i,j,k) - ssync(i,j-1,k)
               yflux(i,j,k) = mult*sy*yarea(i,j,k)*idy
            end do
         end do
      end do
!c
!c     ::::: compute Z fluxes
!c
      do       k = zlo(3), zhi(3)
         do    j = zlo(2), zhi(2)
            do i = zlo(1), zhi(1)
               sz = ssync(i,j,k) - ssync(i,j,k-1)
               zflux(i,j,k) = half*mult*sz*zarea(i,j,k)*idz
            end do
         end do
      end do
      
      end subroutine viscsyncflux

!c :: ----------------------------------------------------------
!c :: SETALPHA
!c ::             alpha(i,j,k) = vol*(1+b/(r(i)^2)) / density
!c ::
!c :: INPUTS / OUTPUTS:
!c ::  fab       <=  array to be modified
!c ::  DIMS(fab) => index limits of fab
!c ::  lo,hi     => index limits of box
!c ::  r         =>  1-d array of radius
!c ::  DIMS(r)   => index limits of r
!c ::  b         =>  viscous coefficient
!c ::  vol       =>  volume array
!c ::  DIMS(vol) => index limits of fab
!c ::  denfab    => array of density at time n+1/2
!c ::  DIMS(den) => index limits of fab
!c ::  usehoop   => do we add hoop stress?   NOT IN 3-D
!c ::  useden    => do we divide by density? (only if velocity component)
!c :: ----------------------------------------------------------
!c ::
       subroutine fort_setalpha (fab, DIMS(fab), lo, hi, r, DIMS(r), b, &
                                 vol, DIMS(vol), denfab, DIMS(den), &
                                 usehoop,useden) &
                                 bind(C,name="fort_setalpha")

       implicit none
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

       end subroutine fort_setalpha

!c :: ----------------------------------------------------------
!c :: SET_TENSOR_ALPHA
!c ::             alpha(i,j) = vol*density
!c ::
!c :: INPUTS / OUTPUTS:
!c ::  fab       <=  array to be modified
!c ::  DIMS(fab) => index limits of fab
!c ::  lo,hi     => index limits of box
!c ::  r         =>  1-d array of radius
!c ::  b         =>  theta*dt or -(1-theta)*dt
!c ::  vol       =>  volume array
!c ::  DIMS(vol) => index limits of fab
!c ::  denfab    => array of density at time n+1/2
!c ::  DIMS(den) => index limits of fab
!c ::  usehoop   => do we add hoop stress?   (only if x-vel component)
!c ::  useden    => do we divide by density? (only if velocity component)
!c :: ----------------------------------------------------------
!c ::
       subroutine set_tensor_alpha (alpha, DIMS(alpha), lo, hi, r, DIMS(r), &
                                b, vol, DIMS(vol), &
                                denfab,DIMS(den),betax,DIMS(betax), &
                                betay,DIMS(betay),betaz,DIMS(betaz),isrz) &
                                bind(C,name="set_tensor_alpha")

       implicit none
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
       REAL_T  b
       integer isrz

       integer i, j, k

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                alpha(i,j,k,1) = vol(i,j,k) * denfab(i,j,k)
             end do
          end do
       end do
 
       end subroutine set_tensor_alpha

      subroutine div_mu_si (lo, hi, dx, mu, DIMS(divu), divu, &
          DIMS(divmusi), divmusi) bind(C,name="div_mu_si")

      implicit none

      integer lo(SDIM), hi(SDIM)
      REAL_T  dx(SDIM)
      integer DIMDEC(divu)
      REAL_T  divu(DIMV(divu))      
      REAL_T  mu

      integer DIMDEC(divmusi)
      REAL_T  divmusi(DIMV(divmusi),SDIM)

      integer i,j,k
      REAL_T sleft, sright, stp, sbot, sfront, sback, idx, idy, idz

      idx = 1.0d0 / dx(1)
      idy = 1.0d0 / dx(2)
      idz = 1.0d0 / dx(3)
!c
!c ... Note: the following IS correct for r-z. Terms from the hoop stress
!c           cancel with terms from tau_rr to eliminate all r dependence.
!c
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               sleft  = (divu(i-1,j,k)+divu(i,j,k))
               sright = (divu(i+1,j,k)+divu(i,j,k))

               divmusi(i,j,k,1) = mu*(sright-sleft)*idx*half

               stp  = (divu(i,j,k)+divu(i,j+1,k))
               sbot = (divu(i,j-1,k)+divu(i,j,k))

               divmusi(i,j,k,2) = mu*(stp-sbot)*idy*half

               sfront = (divu(i,j,k)+divu(i,j,k+1))
               sback  = (divu(i,j,k-1)+divu(i,j,k))

               divmusi(i,j,k,3) = mu*(sfront-sback)*idz*half
            end do
         end do
      end do

      end subroutine div_mu_si

      subroutine div_varmu_si(lo, hi, dx, DIMS(divu), divu, &
          DIMS(betax), betax, DIMS(betay), betay,  DIMS(betaz),  &
          betaz, DIMS(divmusi), divmusi)&
           bind(C,name="div_varmu_si")

      implicit none

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

      integer DIMDEC(divmusi)
      REAL_T  divmusi(DIMV(divmusi),SDIM)

      integer i,j,k
      REAL_T sleft, sright, stp, sbot, sfront, sback, idx, idy, idz

      idx = 1.0d0 / dx(1)
      idy = 1.0d0 / dx(2)
      idz = 1.0d0 / dx(3)

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               sleft  = (divu(i-1,j,k)+divu(i,j,k))
               sright = (divu(i+1,j,k)+divu(i,j,k))

               divmusi(i,j,k,1) = (betax(i+1,j,k)*sright- &
                   betax(i,j,k)*sleft)*idx*half

               stp  = (divu(i,j,k)+divu(i,j+1,k))
               sbot = (divu(i,j-1,k)+divu(i,j,k))

               divmusi(i,j,k,2) = (betay(i,j+1,k)*stp- &
                   betay(i,j,k)*sbot)*idy*half

               sfront = (divu(i,j,k)+divu(i,j,k+1))
               sback  = (divu(i,j,k-1)+divu(i,j,k))

               divmusi(i,j,k,3) = (betaz(i,j,k+1)*sfront- &
                   betaz(i,j,k)*sback)*idz*half

            end do
         end do
      end do

      end subroutine div_varmu_si

  end module diffusion_3d_module
