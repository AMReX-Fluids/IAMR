c
c $Id: NAVIERSTOKES_3D.F,v 1.11 2002-09-19 21:03:59 car Exp $
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

      ddx = fourth/dx(1)
      ddy = fourth/dx(2)
      ddz = fourth/dx(3)

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
            end do
         end do
      end do

      end

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

      do n = 1, nvar
c
c     set coarse grid to zero on overlap
c
         do kc = lo(3), hi(3)
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,kc,n) = zero
               end do
            end do
         end do
c         
c     sum fine data
c
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
     &                          fv(i,j,k)*fine(i,j,k,n)
                        end do
                     end do
                  end do
               end do
            end do
         end do
c
c     divide out by volume weight
c
         do kc = lo(3), hi(3)
            do ic = lo(1), hi(1)
               do jc = lo(2), hi(2)
                  crse(ic,jc,kc,n) = crse(ic,jc,kc,n)/cv(ic,jc,kc)
               end do
            end do
         end do
      end do

      end

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

      REAL_T diff
      integer  ic, jc, kc
      integer  lratx, lraty, lratz

      lratx = ratios(1)
      lraty = ratios(2)
      lratz = ratios(3)

      diff = zero
      do kc = lo(3), hi(3)
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               diff = max( diff,
     &              abs(crse(ic,jc,kc)-fine(lratx*ic,lraty*jc,lratz*kc)))
            end do
         end do
      end do

      if ( diff .gt. 1.0D-10 ) then
         print *, 'mismatch between fine and coarse nodes'
         stop
      end if

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
      integer  lratx, lraty, lratz

      lratx = ratios(1)
      lraty = ratios(2)
      lratz = ratios(3)

      do kc = lo(3), hi(3)
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,kc) = fine(lratx*ic,lraty*jc,lratz*kc)
            end do
         end do
      end do

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

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                a(i,j,k) = a(i,j,k) + alpha*b(i,j,k)
             end do
          end do
       end do

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
       end do

       do k = ARG_L3(grid), ARG_H3(grid)
          do i = ARG_L1(grid), ARG_H1(grid)
             do j = ARG_L2(grid), ARG_H2(grid)
                tmp(j) = tmp(j) + vol*rho(i,j,k)
             end do
          end do
       end do

       mass = zero
       do j = ARG_L2(grid), ARG_H2(grid)
          mass = mass + tmp(j)
       end do

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
       end do

       do k = ARG_L3(grid), ARG_H3(grid)
          do i = ARG_L1(grid), ARG_H1(grid)
             do j = ARG_L2(grid), ARG_H2(grid)
                tmp(j) = tmp(j) + vol*rho(i,j,k)*rho(i,j,k)
             end do
          end do
       end do

       mass = zero
       do j = ARG_L2(grid), ARG_H2(grid)
          mass = mass + tmp(j)
       end do

       end

c :: ----------------------------------------------------------
c :: SUMPRESS
c ::             SUM = sum{ press(i,j,k) * surroundingVolume(i,j,k) }
c ::
c :: INPUTS / OUTPUTS:
c ::  press      => pressure field on nodes
c ::  DIMS(press)=> index limits of press aray
c ::  vol        => volume field (zero-ed out under finer grids)
c ::  DIMS(vol)  => index limits of vol aray
c ::  lo,hi      => index limits of grid interior in nodal space
c ::  total_press<= sum of volume-weighted nodal pressures
c ::  total_vol  <= sum of volumes
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SUMPRESS(press,DIMS(press),vol,DIMS(vol),lo,hi,
     $                          total_press,total_vol)

       integer lo(3),hi(3)
       integer DIMDEC(press)
       integer DIMDEC(vol)
       integer DIMDEC(grid)
       REAL_T  press(DIMV(press))
       REAL_T  vol(DIMV(vol))
       REAL_T  total_press,total_vol

       integer i, j, k
       REAL_T  surroundingVolume

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                surroundingVolume = vol(i,j  ,k  ) + vol(i-1,j  ,k  ) + 
     $               vol(i,j-1,k  ) + vol(i-1,j-1,k  ) +
     $               vol(i,j  ,k-1) + vol(i-1,j  ,k-1) + 
     $               vol(i,j-1,k-1) + vol(i-1,j-1,k-1)
                total_press = total_press + press(i,j,k) * surroundingVolume
                total_vol   = total_vol   +                surroundingVolume

             end do
          end do
       end do

       end


c ::
c :: ----------------------------------------------------------
c :: This routine fills an edge-centered fab from a cell-centered
c :: fab using simple linear interpolation.
c ::
c :: INPUTS / OUTPUTS:
c ::  lo,hi      => index limits of the region of the edge-centered fab
c ::                to be filled
c ::  DIMS(cfab) => index limits of the cell-centered fab
c ::  cfab       => cell-centered data
c ::  DIMS(efab) => index limits of the edge-centered fab
c ::  efab       => edge-centered fab to fill
c ::  nc         => Number of components in the fab to fill
c ::  dir        => direction data needs to be shifted to get to edges
c :: ----------------------------------------------------------
c ::
      subroutine FORT_CEN2EDG(lo, hi, 
     &     DIMS(cfab), cfab,
     &     DIMS(efab), efab, nc, dir, isharm
     &     )
      integer lo(SDIM), hi(SDIM), nc, dir, isharm
      integer DIMDEC(cfab)
      integer DIMDEC(efab)
      REAL_T  cfab(DIMV(cfab), nc)
      REAL_T  efab(DIMV(efab), nc)

      integer i,j,k,n

      if ( isharm .ne. 0 ) STOP 'ISHARM /= 0 not supported in 3D'
      
      if (dir .EQ. 0) then
         do n = 1,nc
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1)+1, hi(1)
                     efab(i,j,k,n) = half*(cfab(i,j,k,n) + cfab(i-1,j,k,n))
                  end do
               end do
            end do
         end do
      else if (dir .EQ. 1) then
         do n = 1,nc
            do k = lo(3), hi(3)
               do j = lo(2)+1, hi(2)
                  do i = lo(1), hi(1)
                     efab(i,j,k,n) = half*(cfab(i,j,k,n) + cfab(i,j-1,k,n))
                  end do
               end do
            end do
         end do
      else if (dir .EQ. 2) then
         do n = 1,nc
            do k = lo(3)+1, hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     efab(i,j,k,n) = half*(cfab(i,j,k,n) + cfab(i,j,k-1,n))
                  end do
               end do
            end do
         end do
      end if
      end
