c
c $Id: EXACT_2D.F,v 1.1 1999-02-18 00:56:35 sstanley Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "PROB_AMR_F.H"
#include "ArrayLim.H"
#include "EXACT_F.H"

#define SDIM 2
c
c ::: -----------------------------------------------------------
c ::: This case is an unsteady  viscous benchmark for which the 
c ::: exact solution is,
c :::     u(x,y,t) = - Cos(Pi x) Sin(Pi y) Exp(-2 Pi^2 Nu t)
c :::     v(x,y,t) =   Sin(Pi x) Cos(Pi y) Exp(-2 Pi^2 Nu t)
c :::     p(x,y,t) = - {Cos(2 Pi x) + Cos(2 Pi y)} Exp(-4 Pi^2 Nu t) / 4
c ::: In the utilities, Parallel/utils/Convergence. there is a 
c ::: tool Inc2dBench.cpp that reads a plot file and compares the
c ::: solution against this exact solution.  This benchmark was
c ::: originally derived by G.I. Taylor (Phil. Mag., Vol. 46, No. 274, 
c ::: pp. 671-674, 1923) and Ethier and Steinman
c ::: (Intl. J. Num. Meth. Fluids, Vol. 19, pp. 369-375, 1994) give
c ::: the pressure field.
c
      subroutine FORT_VISCBENCH(time, nu, lo, hi, ncomp, vel, DIMS(vel),
     &                          dx, xlo, xhi)

      implicit none

      integer    ncomp
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(vel)
      REAL_T     time, nu, dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     vel(DIMV(vel),ncomp)
c
c     ::::: local variables
c
      integer i, j
      REAL_T  x, y
      REAL_T  hx, hy
      REAL_T  spx, spy, cpx, cpy
      REAL_T expterm

      hx = dx(1)
      hy = dx(2)

      expterm = exp(-two*Pi**2*nu*time)

      do j = lo(2), hi(2)
         y = xlo(2) + hy*(float(j-lo(2)) + half)
         spy = sin(Pi*y)
         cpy = cos(Pi*y)

         do i = lo(1), hi(1)
            x = xlo(1) + hx*(float(i-lo(1)) + half)

            spx = sin(Pi*x)
            cpx = cos(Pi*x)

            vel(i,j,1) = - cpx*spy * expterm
            vel(i,j,2) =   spx*cpy * expterm
         end do
      end do

      end

