c
c $Id: VISCOPERATOR_2D.F,v 1.3 1997-12-12 22:44:12 lijewski Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "VISCOPERATOR_F.H"
#include "ArrayLim.H"

#define SDIM 2

c :: ----------------------------------------------------------
c :: VISCEXTRAP
c ::             extrapolate viscous terms to ghost zones
c ::
c :: INPUTS / OUTPUTS:
c ::  visc        <=> array containing visc terms
c ::  DIMS(visc)   => index limits for visc
c ::  lo,hi        => index limits for grid interior
c ::  ncomp        => number of components
c :: ----------------------------------------------------------
c ::
       subroutine FORT_VISCEXTRAP (visc,DIMS(visc),lo,hi,ncomp)

       integer DIMDEC(visc)
       integer lo(SDIM), hi(SDIM)
       integer ncomp
       REAL_T  visc(DIMV(visc),ncomp)

       integer i, j, n

       do n = 1, ncomp
c         ::::: left side
          i = lo(1)-1
          do j = lo(2), hi(2)
c            visc(i,j,n) = two*visc(i+1,j,n) - visc(i+2,j,n)
             visc(i,j,n) = visc(i+1,j,n)
          enddo

c         ::::: right side
          i = hi(1)+1
          do j = lo(2), hi(2)
c            visc(i,j,n) = two*visc(i-1,j,n) - visc(i-2,j,n)
             visc(i,j,n) = visc(i-1,j,n)
          enddo

c         ::::: bottom side
          j = lo(2)-1
          do i = lo(1), hi(1)
c            visc(i,j,n) = two*visc(i,j+1,n) - visc(i,j+2,n)
             visc(i,j,n) = visc(i,j+1,n)
          enddo

c         ::::: top side
          j = hi(2)+1
          do i = lo(1), hi(1)
c            visc(i,j,n) = two*visc(i,j-1,n) - visc(i,j-2,n)
             visc(i,j,n) = visc(i,j-1,n)
          enddo

c         ::::: corners
          i = lo(1)-1
          j = lo(2)-1
c         visc(i,j,n) = two*visc(i+1,j+2,n) - visc(i+2,j+2,n)
          visc(i,j,n) = visc(i+1,j+1,n)

          j = hi(2)+1
c         visc(i,j,n) = two*visc(i+1,j-1,n) - visc(i+2,j-2,n)
          visc(i,j,n) = visc(i+1,j-1,n)

          i = hi(1)+1
c         visc(i,j,n) = two*visc(i-1,j-1,n) - visc(i-2,j-2,n) 
          visc(i,j,n) = visc(i-1,j-1,n)

          j = lo(2)-1
c         visc(i,j,n) = two*visc(i-1,j+1,n) - visc(i-2,j+2,n)
          visc(i,j,n) = visc(i-1,j+1,n)
       enddo

       return
       end

