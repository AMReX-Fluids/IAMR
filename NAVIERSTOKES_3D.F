c
c $Id: NAVIERSTOKES_3D.F,v 1.3 1997-12-12 22:44:07 lijewski Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "NAVIERSTOKES_F.H"
#include "ArrayLim.H"

#define SDIM 3

c=========================================================


      subroutine FORT_TEST_UMAC_PERIODIC(
     &     lo,hi,
     &     u,DIMS(u),xper,udiff,
     &     v,DIMS(v),yper,vdiff,
     &     w,DIMS(w),zper,wdiff,
     &     level)
c     
c     This function tests u_mac for periodicity
c
      implicit none
      integer i,j,k
      integer imin,jmin,kmin,imax,jmax,kmax
c
      REAL_T udiff,vdiff,wdiff
      integer xper,yper,zper
      integer lo(SDIM),hi(SDIM),level
      integer DIMDEC(u),DIMDEC(v),DIMDEC(w)
      REAL_T u(DIMV(u)),v(DIMV(v)),w(DIMV(w))
c
      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)
      udiff = zero
      vdiff = zero
      wdiff = zero
c
      if ( xper .eq. 1 ) then
         do k = kmin,kmax
            do j = jmin,jmax
               udiff = MAX( udiff, abs(u(imax+1,j,k)-u(imin,j,k)) )
            enddo
         enddo
         if ( udiff .gt. 1.e-10 ) then
            udiff = udiff
         endif
      endif
c     
      if ( yper .eq. 1 ) then
         do k = kmin,kmax
            do i = imin,imax
               vdiff = MAX( vdiff, abs(v(i,jmax+1,k)-v(i,jmin,k)) )
            enddo
         enddo
         if ( vdiff .gt. 1.e-10 ) then
            vdiff = vdiff
         endif
      endif
c     
      if ( zper .eq. 1 ) then
         do j = jmin,jmax
            do i = imin,imax
               wdiff = MAX( wdiff, abs(w(i,j,kmax+1)-w(i,j,kmin)) )
            enddo
         enddo
         if ( vdiff .gt. 1.e-10 ) then
            vdiff = vdiff
         endif
      endif
c     
c     end of routine
c
         return
         end



c=========================================================

      subroutine FORT_GRADP(
     &     p,DIMS(p),
     &     gp,DIMS(gp),
     &     lo,hi,dx)
c ::
c :: ----------------------------------------------------------
c :: Compute a cell centered gradient from a node
c :: centered field.  Returns all components of GRADP
c :: ----------------------------------------------------------
c ::
      integer DIMDEC(p)  
      integer DIMDEC(gp)  
      integer    lo(SDIM),  hi(SDIM)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  dx(SDIM)
      REAL_T  p(DIMV(p))
      REAL_T  gp(DIMV(gp),SDIM)
      logical    invalid
      integer    i,j,k
      REAL_T     ddx, ddy, ddz
      integer    ilo, ihi, jlo, jhi, klo, khi
c
      ddx = fourth/dx(1)
      ddy = fourth/dx(2)
      ddz = fourth/dx(3)
c
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               gp(i,j,k,1) = ddx*(
     &              p(i+1,j,k  )-p(i,j,k  )+p(i+1,j+1,k  )-p(i,j+1,k  )+
     &              p(i+1,j,k+1)-p(i,j,k+1)+p(i+1,j+1,k+1)-p(i,j+1,k+1))

               gp(i,j,k,2) = ddy*(
     &              p(i,j+1,k  )-p(i,j,k  )+p(i+1,j+1,k  )-p(i+1,j,k  )+
     &              p(i,j+1,k+1)-p(i,j,k+1)+p(i+1,j+1,k+1)-p(i+1,j,k+1))

               gp(i,j,k,3) = ddz*(
     &              p(i,  j,k+1)-p(i,  j,k)+p(i,  j+1,k+1)-p(i,  j+1,k)+
     &              p(i+1,j,k+1)-p(i+1,j,k)+p(i+1,j+1,k+1)-p(i+1,j+1,k))
            enddo
         enddo
      enddo
c
c     end of routine
c     
      return
      end



c=========================================================

      subroutine FORT_AVGDOWN (
     &     crse,DIMS(crse),nvar,
     &     fine,DIMS(fine),
     &     cv,DIMS(cv),
     &     fv,DIMS(fv),
     &     lo,hi,ratios)
c     ----------------------------------------------------------
c     Volume-weight average the fine grid data onto the coarse
c     grid.  Overlap is given in coarse grid coordinates.
c     
c     crse      =  coarse grid data
c     nvar	= number of components in arrays
c     fine      = fine grid data
c     cv        = coarse grid volume array
c     fv        = fine grid volume array
c     lo,hi     = index limits of overlap (crse grid)
c     ratios    = IntVect refinement ratio
c     ----------------------------------------------------------
      integer  DIMDEC(crse)
      integer  DIMDEC(cv)
      integer  DIMDEC(fine)
      integer  DIMDEC(fv)
      integer  lo(SDIM), hi(SDIM)
      integer  nvar
      integer  ratios(SDIM)
      REAL_T   crse(DIMV(crse),nvar)
      REAL_T     cv(DIMV(cv))
      REAL_T   fine(DIMV(fine),nvar)
      REAL_T     fv(DIMV(fv))

      integer  i, j, k, n, ic, jc, kc, ioff, joff, koff
      integer  lratx, lraty, lratz

      lratx = ratios(1)
      lraty = ratios(2)
      lratz = ratios(3)
c
      do 100 n = 1, nvar

c     set coarse grid to zero on overlap
      do kc = lo(3), hi(3)
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,kc,n) = zero
            enddo
         enddo
      enddo
      
c     sum fine data
      do koff = 0, lratz-1
         do kc = lo(3),hi(3)
	    k = kc*lratz + koff
	    do joff = 0, lraty-1
               do jc = lo(2), hi(2)
                  j = jc*lraty + joff
                  do ioff = 0, lratx-1
                     do ic = lo(1), hi(1)
                        i = ic*lratx + ioff
                        crse(ic,jc,kc,n) = crse(ic,jc,kc,n) +
     &                       fv(i,j,k)*fine(i,j,k,n)
                     enddo
                  enddo
               enddo
	    enddo
         enddo
      enddo

c     divide out by volume weight
      do kc = lo(3), hi(3)
         do ic = lo(1), hi(1)
            do jc = lo(2), hi(2)
               crse(ic,jc,kc,n) = crse(ic,jc,kc,n)/cv(ic,jc,kc)
            enddo
         enddo
      enddo
  100 continue
c
c     end of function
c
      return
      end



c===========================================================


      subroutine FORT_TESTINJECT(
     &     crse,DIMS(crse),
     &     fine,DIMS(fine),lo,hi,ratios)
c
c     test for consistency of pressure nodes
c
      integer  DIMDEC(crse)
      integer  DIMDEC(fine)
      integer  lo(3), hi(3)
      integer  ratios(3)
      REAL_T   crse(DIMV(crse))
      REAL_T   fine(DIMV(fine))
c
      REAL_T diff
      integer  ic, jc, kc
      integer  lratx, lraty, lratz
c
      lratx = ratios(1)
      lraty = ratios(2)
      lratz = ratios(3)
c
      diff = zero
      do kc = lo(3), hi(3)
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               diff = max( diff,
     &              abs(crse(ic,jc,kc)-fine(lratx*ic,lraty*jc,lratz*kc)))
            enddo
         enddo
      enddo
c
      if ( diff .gt. 1.0e-10 ) then
         print *, 'mismatch between fine and coarse nodes'
         stop
      endif
c
c     end of routine
c
      return
      end


c :: ----------------------------------------------------------
c :: Replace coarse grid pressure data with corresponding
c :: fine grid pressure data.
c ::
c :: INPUTS / OUTPUTS:
c ::  crse      <=  coarse grid data
c ::  DIMS(crse) => index limits of crse
c ::  fine       => fine grid data
c ::  DIMS(fine) => index limits of fine
c ::  lo,hi      => index limits of overlap (crse grid)
c ::  ratios     => refinement ratio
c ::
c :: NOTE:
c ::  Assumes pressure fields are node based
c :: ----------------------------------------------------------
c ::
      subroutine FORT_PUTDOWN (crse,DIMS(crse),
     &			       fine,DIMS(fine),lo,hi,ratios)

      integer  DIMDEC(crse)
      integer  DIMDEC(fine)
      integer  lo(SDIM), hi(SDIM)
      integer  ratios(SDIM)
      REAL_T   crse(DIMV(crse))
      REAL_T   fine(DIMV(fine))

      integer  ic, jc, kc
      integer  lenx, leny, lenz, mxlen
      integer  lratx, lraty, lratz

      lratx = ratios(1)
      lraty = ratios(2)
      lratz = ratios(3)

      lenx = hi(1)-lo(1)+1
      leny = hi(2)-lo(2)+1
      lenz = hi(3)-lo(3)+1
      mxlen = max(lenx,leny,lenz)

      if (lenx .eq. mxlen) then

          do kc = lo(3), hi(3)
          do jc = lo(2), hi(2)
	  do ic = lo(1), hi(1)
	    crse(ic,jc,kc) = fine(lratx*ic,lraty*jc,lratz*kc)
	  enddo
	  enddo
	  enddo

      elseif (leny .eq. mxlen) then

          do kc = lo(3), hi(3)
	  do ic = lo(1), hi(1)
          do jc = lo(2), hi(2)
	    crse(ic,jc,kc) = fine(lratx*ic,lraty*jc,lratz*kc)
	  enddo
	  enddo
	  enddo

      else

	  do ic = lo(1), hi(1)
          do jc = lo(2), hi(2)
          do kc = lo(3), hi(3)
	    crse(ic,jc,kc) = fine(lratx*ic,lraty*jc,lratz*kc)
	  enddo
	  enddo
	  enddo

      endif


      return
      end

c :: ----------------------------------------------------------
c :: UTILITY ROUTINE: compute:
c ::             A += alpha*B on subrange
c ::
c :: INPUTS / OUTPUTS:
c ::  a         <=  output array
c ::  b          => input array
c ::  alo,ahi    => index limits of a array
c ::  blo,bhi    => index limits of a array
c ::  lo,hi      => index limits of update region
c ::  alpha      => multiplicative factor
c :: ----------------------------------------------------------
c ::
       subroutine FORT_INCRMULT(a,DIMS(a),b,DIMS(b),lo,hi,alpha)

       integer    DIMDEC(a)
       integer    DIMDEC(b)
       integer    lo(SDIM), hi(SDIM)
       REAL_T     alpha
       REAL_T     a(DIMV(a))
       REAL_T     b(DIMV(b))

       integer i, j, k
       integer lenx, leny, lenz, mxlen

       lenx = hi(1)-lo(1)+1
       leny = hi(2)-lo(2)+1
       lenz = hi(3)-lo(3)+1
       mxlen = max(lenx,leny,lenz)

       if (lenx .eq. mxlen) then
	  do i = lo(1), hi(1)
          do j = lo(2), hi(2)
          do k = lo(3), hi(3)
	     a(i,j,k) = a(i,j,k) + alpha*b(i,j,k)
	  enddo
	  enddo
	  enddo
       elseif (leny .eq. mxlen) then
	  do i = lo(1), hi(1)
          do j = lo(2), hi(2)
          do k = lo(3), hi(3)
	     a(i,j,k) = a(i,j,k) + alpha*b(i,j,k)
	  enddo
	  enddo
	  enddo
       else
	  do i = lo(1), hi(1)
          do j = lo(2), hi(2)
          do k = lo(3), hi(3)
	     a(i,j,k) = a(i,j,k) + alpha*b(i,j,k)
	  enddo
	  enddo
	  enddo
       endif

       return
       end


c :: ----------------------------------------------------------
c :: SUMMASS
c ::             MASS = sum{ vol(i,j)*rho(i,j) }
c ::
c :: INPUTS / OUTPUTS:
c ::  rho        => density field
c ::  DIMS(rho)  => index limits of rho aray
c ::  lo,hi      => index limits of grid interior
c ::  delta	 => cell size
c ::  mass      <=  total mass
c ::  r		 => radius at cell center
c ::  tmp        => temp column array
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SUMMASS(rho,DIMS(rho),DIMS(grid),delta,mass,tmp)

       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  mass, delta(SDIM)
       REAL_T  rho(DIMV(rho))
       REAL_T  tmp(DIM2(grid))

       integer i, j, k
       REAL_T  vol

       vol = delta(1)*delta(2)*delta(3)

       do j = ARG_L2(grid), ARG_H2(grid)
          tmp(j) = zero
       enddo

       do k = ARG_L3(grid), ARG_H3(grid)
       do i = ARG_L1(grid), ARG_H1(grid)
       do j = ARG_L2(grid), ARG_H2(grid)
          tmp(j) = tmp(j) + vol*rho(i,j,k)
       enddo
       enddo
       enddo

       mass = zero
       do j = ARG_L2(grid), ARG_H2(grid)
          mass = mass + tmp(j)
       enddo

       return
       end


c :: ----------------------------------------------------------
c :: SUMMASSSQ
c ::             MASS = sum{ vol(i,j)*rho(i,j)*rho(i,j) }
c ::
c :: INPUTS / OUTPUTS:
c ::  rho        => density field
c ::  DIMS(rho)  => index limits of rho aray
c ::  lo,hi      => index limits of grid interior
c ::  delta	 => cell size
c ::  mass      <=  total mass
c ::  r		 => radius at cell center
c ::  tmp        => temp column array
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SUMMASSSQ(rho,DIMS(rho),DIMS(grid),delta,mass,tmp)

       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  mass, delta(SDIM)
       REAL_T  rho(DIMV(rho))
       REAL_T  tmp(DIM2(grid))

       integer i, j, k
       REAL_T  vol

       vol = delta(1)*delta(2)*delta(3)

       do j = ARG_L2(grid), ARG_H2(grid)
          tmp(j) = zero
       enddo

       do k = ARG_L3(grid), ARG_H3(grid)
       do i = ARG_L1(grid), ARG_H1(grid)
       do j = ARG_L2(grid), ARG_H2(grid)
          tmp(j) = tmp(j) + vol*rho(i,j,k)*rho(i,j,k)
       enddo
       enddo
       enddo

       mass = zero
       do j = ARG_L2(grid), ARG_H2(grid)
          mass = mass + tmp(j)
       enddo

       return
       end
