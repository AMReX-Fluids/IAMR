
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_ArrayLim.H>
#include <EXACT_F.H>

#define SDIM 2
c
c ::: -----------------------------------------------------------
c ::: This case is an unsteady  viscous benchmark for which the 
c ::: exact solution is,
c :::     u(x,y,t) = - Cos(Pi x) Sin(Pi y) Exp(-2 Pi^2 Nu t)
c :::     v(x,y,t) =   Sin(Pi x) Cos(Pi y) Exp(-2 Pi^2 Nu t)
c :::     p(x,y,t) = - {Cos(2 Pi x) + Cos(2 Pi y)} Exp(-4 Pi^2 Nu t) / 4
c ::: In this dircetory, ViscBench2d.cpp, reads a plot file and compares
c ::: the solution against this exact solution.  This benchmark was
c ::: originally derived by G.I. Taylor (Phil. Mag., Vol. 46, No. 274, 
c ::: pp. 671-674, 1923) and Ethier and Steinman
c ::: (Intl. J. Num. Meth. Fluids, Vol. 19, pp. 369-375, 1994) give
c ::: the pressure field.
c
      subroutine FORT_VISCBENCH(time, nu, lo, hi, 
     &                          ncomp, state, DIMS(state),
     &                          dx, xlo, xhi)

      implicit none

      integer    ncomp
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(state)
      REAL_T     time, nu, dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     state(DIMV(state),ncomp)
c
c     ::::: local variables
c
      integer i, j
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  spx, spy, cpx, cpy
      REAL_T  expterm

      hx = dx(1)
      hy = dx(2)

      expterm = exp(-two*four*Pi**2*nu*time)

      do j = lo(2), hi(2)
         y = xlo(2) + hy*(float(j-lo(2)) + half)
         spy = sin(two*Pi*y)
         cpy = cos(two*Pi*y)

         do i = lo(1), hi(1)
            x = xlo(1) + hx*(float(i-lo(1)) + half)

            spx = sin(two*Pi*x)
            cpx = cos(two*Pi*x)

            state(i,j,1) =   spx*cpy * expterm
            state(i,j,2) = - cpx*spy * expterm
            state(i,j,3) =   one
         end do
      end do

      end

