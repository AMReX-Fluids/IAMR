c
c $Id: NAVIERSTOKES_2D.F,v 1.14 2003-01-30 07:16:29 almgren Exp $
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

#define SDIM 2

      subroutine FORT_GRADP (
     &     p,DIMS(p),
     &     gp,DIMS(gp),
     &     lo,hi,dx,is_full)
c 
c     Compute a cell centered gradient from a node
c     centered field.  Returns all components of GRADP
c     
      integer    DIMDEC(p)
      integer    DIMDEC(gp)
      integer    lo(SDIM), hi(SDIM)
      REAL_T     dx(SDIM)
      REAL_T     p(DIMV(p))
      REAL_T     gp(DIMV(gp),SDIM)
      integer    is_full
      integer    i,j
      REAL_T     ddx, ddy

      ddx = half/dx(1)
      ddy = half/dx(2)

      if (is_full .eq. 0) then
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          gp(i,j,1) = ddx*(p(i+1,j)-p(i,j)+p(i+1,j+1)-p(i,j+1))
          gp(i,j,2) = ddy*(p(i,j+1)-p(i,j)+p(i+1,j+1)-p(i+1,j))
        end do
        end do
      else
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          gp(i,j,1) = (p(i+1,j)-p(i,j)+p(i+1,j+1)-p(i,j+1))
          gp(i,j,2) = (p(i,j+1)-p(i,j)+p(i+1,j+1)-p(i+1,j))
        end do
        end do
      endif

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

      integer  i, j, n, ic, jc, ioff, joff
      integer  lratx,lraty

      lratx = ratios(1)
      lraty = ratios(2)

      do n = 1, nvar
c
c     set coarse grid to zero on overlap
c
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,n) = zero
            end do
         end do
c
c     sum fine data
c
         do joff = 0, lraty-1
            do jc = lo(2), hi(2)
               j = jc*lraty + joff
               do ioff = 0, lratx-1
                  do ic = lo(1), hi(1)
                     i = ic*lratx + ioff
                     crse(ic,jc,n) = crse(ic,jc,n) +
     &                    fv(i,j)*fine(i,j,n)
                  end do
               end do
            end do
         end do
c
c     divide out by volume weight
c
         do ic = lo(1), hi(1)
            do jc = lo(2), hi(2)
               crse(ic,jc,n) = crse(ic,jc,n)/cv(ic,jc)
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
      integer  lo(2), hi(2)
      integer  ratios(2)
      REAL_T   crse(DIMV(crse))
      REAL_T   fine(DIMV(fine))

      REAL_T diff
      integer ic, jc
      integer lratx, lraty

      lratx = ratios(1)
      lraty = ratios(2)

      diff = zero
      do jc = lo(2), hi(2)
         do ic = lo(1), hi(1)
            diff = max(diff,abs(crse(ic,jc)-fine(lratx*ic,lraty*jc)))
         end do
      end do

      if ( diff .gt. 1.0D-10 ) then
         call bl_abort('mismatch between fine and coarse nodes')
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
c ::  ratios     => IntVect refinement ratio
c ::
c :: NOTE:
c ::  Assumes pressure fields are node based
c :: ----------------------------------------------------------
c ::
      subroutine FORT_PUTDOWN (crse,DIMS(crse),
     &			       fine,DIMS(fine),lo,hi,ratios)
      integer  DIMDEC(crse)
      integer  DIMDEC(fine)
      integer  lo(2), hi(2)
      integer  ratios(2)
      REAL_T   crse(DIMV(crse))
      REAL_T   fine(DIMV(fine))

      integer  ic, jc
      integer  lratx, lraty

      lratx = ratios(1)
      lraty = ratios(2)

      do jc = lo(2), hi(2)
         do ic = lo(1), hi(1)
            crse(ic,jc) = fine(lratx*ic,lraty*jc)
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
c ::  DIMS(a)    => index limits of a array
c ::  DIMS(b)    => index limits of a array
c ::  lo,hi      => index limits of update region
c ::  alpha      => multiplicative factor
c :: ----------------------------------------------------------
c ::
       subroutine FORT_INCRMULT(a,DIMS(a),b,DIMS(b),lo,hi,alpha)

       integer    DIMDEC(a)
       integer    DIMDEC(b)
       integer    lo(2), hi(2)
       REAL_T     alpha
       REAL_T     a(DIMV(a))
       REAL_T     b(DIMV(b))

       integer i, j

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             a(i,j) = a(i,j) + alpha*b(i,j)
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
c ::  dx	 => cell size
c ::  mass      <=  total mass
c ::  r		 => radius at cell center
c ::  irlo,hi    => index limits of r array
c ::  rz_flag    => == 1 if R_Z coords
c ::  tmp        => temp column array
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SUMMASS(rho,DIMS(rho),DIMS(grid),dx,mass,
     &                         r,irlo,irhi,rz_flag,tmp)

       integer irlo, irhi, rz_flag
       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  mass, dx(2)
       REAL_T  rho(DIMV(rho))
       REAL_T  r(irlo:irhi)
       REAL_T  tmp(DIM2(grid))

       integer i, j
       REAL_T  dr, dz, vol

       dr = dx(1)
       dz = dx(2)

       do j = ARG_L2(grid), ARG_H2(grid)
          tmp(j) = zero
       end do

       do i = ARG_L1(grid), ARG_H1(grid)
          vol = dr*dz
	  if (rz_flag .eq. 1) vol = vol*two*Pi*r(i)
          do j = ARG_L2(grid), ARG_H2(grid)
	     tmp(j) = tmp(j) + vol*rho(i,j)
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
c ::  ds   	 => cell size
c ::  mass      <=  total mass
c ::  r		 => radius at cell center
c ::  irlo,hi    => index limits of r array
c ::  rz_flag    => == 1 if R_Z coords
c ::  tmp        => temp column array
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SUMMASSSQ(rho,DIMS(rho),DIMS(grid),dx,mass,
     &                           r,irlo,irhi,rz_flag,tmp)

       integer irlo, irhi, rz_flag
       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  mass, dx(SDIM)
       REAL_T  rho(DIMV(rho))
       REAL_T  r(irlo:irhi)
       REAL_T  tmp(DIM2(grid))

       integer i, j
       REAL_T  dr, dz, vol

       dr = dx(1)
       dz = dx(2)

       do j = ARG_L2(grid), ARG_H2(grid)
          tmp(j) = zero
       end do

       do i = ARG_L1(grid), ARG_H1(grid)
          vol = dr*dz
	  if (rz_flag .eq. 1) vol = vol*two*Pi*r(i)
          do j = ARG_L2(grid), ARG_H2(grid)
             tmp(j) = tmp(j) + vol*rho(i,j)*rho(i,j)
          end do
       end do

       mass = zero
       do j = ARG_L2(grid), ARG_H2(grid)
          mass = mass + tmp(j)
       end do

       end

c :: ----------------------------------------------------------
c :: SUMPRESS
c ::             SUM = sum{ press(i,j) * surroundingVolume(i,j) }
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

       integer lo(2),hi(2)
       integer DIMDEC(press)
       integer DIMDEC(vol)
       REAL_T  press(DIMV(press))
       REAL_T  vol(DIMV(vol))
       REAL_T  total_press,total_vol

       integer i, j
       REAL_T  surroundingVolume

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             surroundingVolume = vol(i,j) + vol(i-1,j) + vol(i,j-1) + vol(i-1,j-1)
             total_press = total_press + press(i,j) * surroundingVolume
             total_vol   = total_vol   +              surroundingVolume

          end do
       end do

       end


c ::
c :: ----------------------------------------------------------
c :: This routine fills an edge-centered fab from a cell-centered
c :: fab using simple linear interpolation.
c ::
c :: INPUTS / OUTPUTS:
c ::  lo,hi      => index limits of the of the cell-centered fab
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
     &     DIMS(efab), efab, nc, dir,
     &     isharm)
      integer lo(SDIM), hi(SDIM), nc, dir, isharm
      integer DIMDEC(cfab)
      integer DIMDEC(efab)
      REAL_T  cfab(DIMV(cfab), nc)
      REAL_T  efab(DIMV(efab), nc)

      integer i,j,n

      if ( isharm .eq. 0 ) then
         if (dir .EQ. 0) then
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1)+1, hi(1)
                     efab(i,j,n) = half*(cfab(i,j,n) + cfab(i-1,j,n))
                  end do
               end do
            end do
         else
            do n = 1,nc
               do j = lo(2)+1, hi(2)
                  do i = lo(1), hi(1)
                     efab(i,j,n) = half*(cfab(i,j,n) + cfab(i,j-1,n))
                  end do
               end do
            end do
         end if
      else
         if (dir .EQ. 0) then
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1)+1, hi(1)
                     if((cfab(i,j,n) * cfab(i-1,j,n)).gt.0.d0)then
                        efab(i,j,n)
     &                       = 2.d0*(cfab(i,j,n) * cfab(i-1,j,n))/
     &                       (cfab(i,j,n) + cfab(i-1,j,n))
                     else
                        efab(i,j,n)=0.d0
                     endif
                  end do
               end do
            end do
         else
            do n = 1,nc
               do j = lo(2)+1, hi(2)
                  do i = lo(1), hi(1)
                     if((cfab(i,j,n) * cfab(i,j-1,n)).gt.0.d0)then
                        efab(i,j,n)
     &                       = 2.d0*(cfab(i,j,n) * cfab(i,j-1,n))/
     &                       (cfab(i,j,n) + cfab(i,j-1,n))
                     else
                        efab(i,j,n)=0.d0
                     endif
                  end do
               end do
            end do
         end if
      end if
      end

c-----------------------------------------------------------------------
      subroutine EDGE_INTERP(flo, fhi, nc, ratio, dir,
     &     fine, fine_l0, fine_l1, fine_h0, fine_h1)
      implicit none
      integer flo(0:2-1), fhi(0:2-1), nc, ratio(0:2-1), dir
      integer fine_l0, fine_l1, fine_h0, fine_h1
      DOUBLE PRECISION fine(fine_l0:fine_h0,fine_l1:fine_h1,nc)
      integer i,j,ii,jj,n,P,M,clo(0:2-1),chi(0:2-1)
      DOUBLE PRECISION val, df

c     Do linear in dir, pc transverse to dir, leave alone the fine values
c     lining up with coarse edges--assume these have been set to hold the 
c     values you want to interpolate to the rest.
      if (dir.eq.0) then
         do n=1,nc
            do j=flo(1),fhi(1),ratio(1)
               do i=flo(0),fhi(0)-ratio(dir),ratio(0)
                  df = fine(i+ratio(dir),j,n)-fine(i,j,n)
                  do M=1,ratio(dir)-1
                     val = fine(i,j,n) + df*dble(M)/dble(ratio(dir))
                     do P=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                        fine(i+M,P,n) = val
                     enddo
                  enddo                     
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do j=flo(1),fhi(1)-ratio(dir),ratio(1)
               do i=flo(0),fhi(0)
                  df = fine(i,j+ratio(dir),n)-fine(i,j,n)
                  do M=1,ratio(dir)-1
                     val = fine(i,j,n) + df*dble(M)/dble(ratio(dir))
                     do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                        fine(P,j+M,n) = val
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif

      end
c-----------------------------------------------------------------------
      subroutine PC_EDGE_INTERP(lo, hi, nc, ratio, dir,
     &     crse, crse_l0, crse_l1, crse_h0, crse_h1,
     &     fine, fine_l0, fine_l1, fine_h0, fine_h1)
      implicit none
      integer lo(2),hi(2), nc, ratio(0:2-1), dir
      integer crse_l0, crse_l1, crse_h0, crse_h1
      integer fine_l0, fine_l1, fine_h0, fine_h1
      DOUBLE PRECISION crse(crse_l0:crse_h0,crse_l1:crse_h1,nc)
      DOUBLE PRECISION fine(fine_l0:fine_h0,fine_l1:fine_h1,nc)
      integer i,j,ii,jj,n,L
      DOUBLE PRECISION val, dc

c     For edge-based data, fill fine values with piecewise-constant interp of coarse data.
c     Operate only on faces that overlap--ie, only fill the fine faces that make up each
c     coarse face, leave the in-between faces alone.
      if (dir.eq.0) then
         do n=1,nc
            do j=lo(2),hi(2)
               jj = ratio(1)*j
               do i=lo(1),hi(1)
                  ii = ratio(0)*i
                  do L=0,ratio(1)-1
                     fine(ii,jj+L,n) = crse(i,j,n)
                  enddo
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do j=lo(2),hi(2)
               jj = ratio(1)*j
               do i=lo(1),hi(1)
                  ii = ratio(0)*i
                  do L=0,ratio(0)-1
                     fine(ii+L,jj,n) = crse(i,j,n)
                  enddo
               enddo
            enddo
         enddo
      endif

      end

