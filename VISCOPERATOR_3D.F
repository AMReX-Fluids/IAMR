c $Id: VISCOPERATOR_3D.F,v 1.1 1997-07-08 23:08:15 vince Exp $
#undef BL_LANG_CC
#define LANG_FORT

#include "REAL.H"
#include "CONSTANTS.H"
#include "VISCOPERATOR_F.H"
#include "DIMS.H"

#define SDIM 3

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

       integer i, j, k, n

       do n = 1, ncomp
c         ::::: left side
          i = lo(1)-1
          do k = lo(3), hi(3)
          do j = lo(2), hi(2)
c            visc(i,j,k,n) = two*visc(i+1,j,k,n) - visc(i+2,j,k,n)
             visc(i,j,k,n) = visc(i+1,j,k,n)
          enddo
          enddo

c         ::::: right side
          i = hi(1)+1
          do k = lo(3), hi(3)
          do j = lo(2), hi(2)
c            visc(i,j,k,n) = two*visc(i-1,j,k,n) - visc(i-2,j,k,n)
             visc(i,j,k,n) = visc(i-1,j,k,n)
          enddo
          enddo

c         ::::: bottom side
          j = lo(2)-1
          do k = lo(3), hi(3)
          do i = lo(1), hi(1)
c            visc(i,j,k,n) = two*visc(i,j+1,k,n) - visc(i,j+2,k,n)
             visc(i,j,k,n) = visc(i,j+1,k,n)
          enddo
          enddo

c         ::::: top side
          j = hi(2)+1
          do k = lo(3), hi(3)
          do i = lo(1), hi(1)
c            visc(i,j,k,n) = two*visc(i,j-1,k,n) - visc(i,j-2,k,n)
             visc(i,j,k,n) = visc(i,j-1,k,n)
          enddo
          enddo

c         ::::: down side
          k = lo(3)-1
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
c            visc(i,j,k,n) = two*visc(i,j,k+1,n) - visc(i,j,k+2,n)
             visc(i,j,k,n) = visc(i,j,k+1,n)
          enddo
          enddo

c         ::::: up side
          k = hi(3)+1
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
c            visc(i,j,k,n) = two*visc(i,j,k-1,n) - visc(i,j,k-2,n)
             visc(i,j,k,n) = visc(i,j,k-1,n)
          enddo
          enddo

c         ::::: k-edges
          i = lo(1)-1
          j = lo(2)-1
          do k = lo(3), hi(3)
            visc(i,j,k,n) = visc(i+1,j+1,k,n)
          enddo
 
          i = lo(1)-1
          j = hi(2)+1
          do k = lo(3), hi(3)
            visc(i,j,k,n) = visc(i+1,j-1,k,n)
          enddo

          i = hi(1)+1
          j = lo(2)-1
          do k = lo(3), hi(3)
            visc(i,j,k,n) = visc(i-1,j+1,k,n)
          enddo
 
          i = hi(1)+1
          j = hi(2)+1
          do k = lo(3), hi(3)
            visc(i,j,k,n) = visc(i-1,j-1,k,n)
          enddo

c         ::::: j-edges
          i = lo(1)-1
          k = lo(3)-1
          do j = lo(2), hi(2)
            visc(i,j,k,n) = visc(i+1,j,k+1,n)
          enddo
 
          i = lo(1)-1
          k = hi(3)+1
          do j = lo(2), hi(2)
            visc(i,j,k,n) = visc(i+1,j,k-1,n)
          enddo

          i = hi(1)+1
          k = lo(3)-1
          do j = lo(2), hi(2)
            visc(i,j,k,n) = visc(i-1,j,k+1,n)
          enddo
 
          i = hi(1)+1
          k = hi(3)+1
          do j = lo(2), hi(2)
            visc(i,j,k,n) = visc(i-1,j,k-1,n)
          enddo

c         ::::: i-edges
          j = lo(2)-1
          k = lo(3)-1
          do i = lo(1), hi(1)
            visc(i,j,k,n) = visc(i,j+1,k+1,n)
          enddo
 
          j = lo(2)-1
          k = hi(3)+1
          do i = lo(1), hi(1)
            visc(i,j,k,n) = visc(i,j+1,k-1,n)
          enddo

          j = hi(2)+1
          k = lo(3)-1
          do i = lo(1), hi(1)
            visc(i,j,k,n) = visc(i,j-1,k+1,n)
          enddo
 
          j = hi(2)+1
          k = hi(3)+1
          do i = lo(1), hi(1)
            visc(i,j,k,n) = visc(i,j-1,k-1,n)
          enddo

c         ::::: corners
          i = lo(1)-1
          j = lo(2)-1
          k = lo(3)-1
          visc(i,j,k,n) = visc(i+1,j+1,k+1,n)

          i = lo(1)-1
          j = hi(2)+1
          k = lo(3)-1
          visc(i,j,k,n) = visc(i+1,j-1,k+1,n)

          i = hi(1)+1
          j = hi(2)+1
          k = lo(3)-1
          visc(i,j,k,n) = visc(i-1,j-1,k+1,n)

          i = hi(1)+1
          j = lo(2)-1
          k = lo(3)-1
          visc(i,j,k,n) = visc(i-1,j+1,k+1,n)

          i = lo(1)-1
          j = lo(2)-1
          k = hi(3)+1
          visc(i,j,k,n) = visc(i+1,j+1,k-1,n)

          i = lo(1)-1
          j = hi(2)+1
          k = hi(3)+1
          visc(i,j,k,n) = visc(i+1,j-1,k-1,n)

          i = hi(1)+1
          j = lo(2)-1
          k = hi(3)+1
          visc(i,j,k,n) = visc(i-1,j+1,k-1,n)

          i = hi(1)+1
          j = hi(2)+1
          k = hi(3)+1
          visc(i,j,k,n) = visc(i-1,j-1,k-1,n)

       enddo

       return
       end

