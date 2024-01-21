
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_ArrayLim.H>
#include <EXACT_F.H>

#ifdef BL_USE_FLOAT
#define REAL_T  real
#define zero    0.0e0
#define one     1.0e0
#define two     2.0e0
#define four    4.0e0
#define half    0.5e0
#define Pi      3.1415926535897932e0
#else
#define REAL_T  double precision
#define zero    0.0d0
#define one     1.0d0
#define two     2.0d0
#define four    4.0d0
#define half    0.5d0
#define Pi      3.1415926535897932d0
#endif

#define SDIM 3
c
c ::: -----------------------------------------------------------
c ::: This case is an unsteady  viscous benchmark for which the
c ::: exact solution is,
c :::     u(x,y,t) = - Cos(Pi x) Sin(Pi y) Exp(-2 Pi^2 Nu t)
c :::     v(x,y,t) =   Sin(Pi x) Cos(Pi y) Exp(-2 Pi^2 Nu t)
c :::     p(x,y,t) = - {Cos(2 Pi x) + Cos(2 Pi y)} Exp(-4 Pi^2 Nu t) / 4
c ::: In this dircetory, ViscBench3d.cpp, reads a plot file and compares
c ::: the solution against this exact solution.  This benchmark was
c ::: originally derived by G.I. Taylor (Phil. Mag., Vol. 46, No. 274,
c ::: pp. 671-674, 1923) and Ethier and Steinman
c ::: (Intl. J. Num. Meth. Fluids, Vol. 19, pp. 369-375, 1994) give
c ::: the pressure field.
c
      subroutine FORT_VISCBENCH(time, nu, unifdir, lo, hi,
     &                          ncomp, state,
     &                          state_l1,state_l2,state_l3,
     &                          state_h1,state_h2,state_h3,
     &                          dx, xlo, xhi)

      implicit none

      integer    ncomp, unifdir
      integer    lo(SDIM), hi(SDIM)
      integer    state_l1,state_l2,state_l3
      integer    state_h1,state_h2,state_h3
      REAL_T     time, nu, dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     state(state_l1:state_h1,
     &                 state_l2:state_h2,
     &                 state_l3:state_h3,ncomp)
c
c     ::::: local variables
c
      integer i, j, k
      REAL_T  x, y, z
      REAL_T  hx, hy, hz
      REAL_T  spx, spy, spz, cpx, cpy, cpz
      REAL_T expterm

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      expterm = exp(-two*four*Pi**2*nu*time)

      do k = lo(3), hi(3)
        z = xlo(3) + hz*(float(k-lo(3)) + half)
        spz = sin(two*Pi*z)
        cpz = cos(two*Pi*z)

        do j = lo(2), hi(2)
          y = xlo(2) + hy*(float(j-lo(2)) + half)
          spy = sin(two*Pi*y)
          cpy = cos(two*Pi*y)

          do i = lo(1), hi(1)
            x = xlo(1) + hx*(float(i-lo(1)) + half)

            spx = sin(two*Pi*x)
            cpx = cos(two*Pi*x)

c
c           Uniform in the X-direction
c
            if (unifdir .eq. 0) then
              state(i,j,k,1) =   zero
              state(i,j,k,2) =   spz*cpy * expterm
              state(i,j,k,3) = - cpz*spy * expterm
              state(i,j,k,4) =   one

c
c           Uniform in the Y-direction
c
            elseif (unifdir .eq. 1) then
              state(i,j,k,1) = - cpx*spz * expterm
              state(i,j,k,2) =   zero
              state(i,j,k,3) =   spx*cpz * expterm
              state(i,j,k,4) =   one

c
c           Uniform in the Z-direction
c
            elseif (unifdir .eq. 2) then
              state(i,j,k,1) =   spx*cpy * expterm
              state(i,j,k,2) = - cpx*spy * expterm
              state(i,j,k,3) =   zero
              state(i,j,k,4) =   one
            endif
          end do
        end do
      end do

      end

