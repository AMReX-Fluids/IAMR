c
c $Id: NAVIERSTOKES_3D.F,v 1.18 2007-08-20 19:36:29 almgren Exp $
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
      implicit none

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
      implicit none

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
      implicit none

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
c ::  ratios     => refinement ratio
c ::
c :: NOTE:
c ::  Assumes pressure fields are node based
c :: ----------------------------------------------------------
c ::
      subroutine FORT_PUTDOWN (crse,DIMS(crse),
     &			       fine,DIMS(fine),lo,hi,ratios)

      implicit none

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

       implicit none
      
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

c ::
c :: ----------------------------------------------------------
c :: SUMTURB
c :: ----------------------------------------------------------
c ::
      subroutine FORT_SUMTURB(dat,pres,DIMS(dat),DIMS(grid),delta,
     &                        turb,ksize,turbVars)

      implicit none

      integer DIMDEC(dat)
      integer DIMDEC(grid)
      integer ksize, turbVars
      REAL_T  delta(SDIM)
      REAL_T  dat(DIMV(dat),4)
      REAL_T  pres(DIMV(dat),4)
      REAL_T  turb(0:ksize*turbVars-1)
      
      integer i, j, k, v
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      
      REAL_T area
      REAL_T rho, ux, uy, uz, p, XC, XMg
      REAL_T e, T, f, wC, wMg
      REAL_T duxdx, duxdy, duxdz
      REAL_T duydx, duydy, duydz
      REAL_T duzdx, duzdy, duzdz
      REAL_T dpdx, dpdy, dpdz
      REAL_T dx, dy, dz
      REAL_T tdx, tdy, tdz
      
      ilo = ARG_L1(grid)
      ihi = ARG_H1(grid)
      jlo = ARG_L2(grid)
      jhi = ARG_H2(grid)
      klo = ARG_L3(grid)
      khi = ARG_H3(grid)
      
      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      tdx = two*dx
      tdy = two*dy
      tdz = two*dz
      
      area = dx*dy
      
      do k = klo, khi
         do i = ilo, ihi
            do j = jlo, jhi

               rho  = dat(i,j,k,1)
               ux   = dat(i,j,k,2)
               uy   = dat(i,j,k,3)
               uz   = dat(i,j,k,4)
               p    = pres(i,j,k,1)
               dpdx = pres(i,j,k,2)
               dpdy = pres(i,j,k,3)
               dpdz = pres(i,j,k,4)
               
               duxdx = ( dat(i+1,j,k,2) - dat(i-1,j,k,2) ) / (tdx)
               duydx = ( dat(i+1,j,k,3) - dat(i-1,j,k,2) ) / (tdx)
               duzdx = ( dat(i+1,j,k,4) - dat(i-1,j,k,2) ) / (tdx)
               
               duxdy = ( dat(i,j+1,k,2) - dat(i,j-1,k,2) ) / (tdy)
               duydy = ( dat(i,j+1,k,3) - dat(i,j-1,k,2) ) / (tdy)
               duzdy = ( dat(i,j+1,k,4) - dat(i,j-1,k,2) ) / (tdy)
               
               duxdz = ( dat(i,j,k+1,2) - dat(i,j,k-1,2) ) / (tdz)
               duydz = ( dat(i,j,k+1,3) - dat(i,j,k-1,2) ) / (tdz)
               duzdz = ( dat(i,j,k+1,4) - dat(i,j,k-1,2) ) / (tdz)
               
               turb(k*turbVars+0)  = turb(k*turbVars+0)  + area*rho
               turb(k*turbVars+1)  = turb(k*turbVars+1)  + area*rho*ux
               turb(k*turbVars+2)  = turb(k*turbVars+2)  + area*rho*uy
               turb(k*turbVars+3)  = turb(k*turbVars+3)  + area*rho*uz
               turb(k*turbVars+4)  = turb(k*turbVars+4)  + area*rho*ux*ux
               turb(k*turbVars+5)  = turb(k*turbVars+5)  + area*rho*ux*uy
               turb(k*turbVars+6)  = turb(k*turbVars+6)  + area*rho*ux*uz
               turb(k*turbVars+7)  = turb(k*turbVars+7)  + area*rho*uy*uy
               turb(k*turbVars+8)  = turb(k*turbVars+8)  + area*rho*uy*uz
               turb(k*turbVars+9)  = turb(k*turbVars+9)  + area*rho*uz*uz
               turb(k*turbVars+10) = turb(k*turbVars+10) + area*rho*ux*ux*ux
               turb(k*turbVars+11) = turb(k*turbVars+11) + area*rho*ux*ux*uy
               turb(k*turbVars+12) = turb(k*turbVars+12) + area*rho*ux*ux*uz
               turb(k*turbVars+13) = turb(k*turbVars+13) + area*rho*ux*uy*uy
               turb(k*turbVars+14) = turb(k*turbVars+14) + area*rho*ux*uy*uz
               turb(k*turbVars+15) = turb(k*turbVars+15) + area*rho*ux*uz*uz
               turb(k*turbVars+16) = turb(k*turbVars+16) + area*rho*uy*uy*uy
               turb(k*turbVars+17) = turb(k*turbVars+17) + area*rho*uy*uy*uz
               turb(k*turbVars+18) = turb(k*turbVars+18) + area*rho*uz*uz*uz

            end do
         end do
      end do
      
      end

       subroutine FORT_MAXVAL(rho,DIMS(rho),DIMS(grid),tmp,mass)

       implicit none

       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  rho(DIMV(rho))
       REAL_T  tmp(DIM3(grid))
       REAL_T  mass

       integer i, j, k

       do k = ARG_L3(grid), ARG_H3(grid)
          tmp(k) = zero
       end do

       do k = ARG_L3(grid), ARG_H3(grid)
          do i = ARG_L1(grid), ARG_H1(grid)
             do j = ARG_L2(grid), ARG_H2(grid)
                tmp(k) = max( tmp(k) , rho(i,j,k) )
             end do
          end do
       end do

       mass = zero
       do j = ARG_L3(grid), ARG_H3(grid)
          mass = max(mass , tmp(k))
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

       implicit none

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

       implicit none

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

       implicit none

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
c-----------------------------------------------------------------------
c     This routine fills an edge-centered fab from a cell-centered
c     fab using simple linear interpolation.
c
c     INPUTS / OUTPUTS:
c     lo,hi      => index limits of the region of the edge-centered fab
c                   to be filled
c     DIMS(cfab) => index limits of the cell-centered fab
c     cfab       => cell-centered data
c     DIMS(efab) => index limits of the edge-centered fab
c     efab       => edge-centered fab to fill
c     nc         => Number of components in the fab to fill
c     dir        => direction data needs to be shifted to get to edges
c-----------------------------------------------------------------------
c
      subroutine FORT_CEN2EDG(lo, hi, 
     &     DIMS(cfab), cfab,
     &     DIMS(efab), efab, nc, dir, isharm
     &     )
      implicit none

      integer lo(SDIM), hi(SDIM), nc, dir, isharm
      integer DIMDEC(cfab)
      integer DIMDEC(efab)
      REAL_T  cfab(DIMV(cfab), nc)
      REAL_T  efab(DIMV(efab), nc)
      integer i,j,k,n

      if ( isharm .eq. 0 ) then
         if (dir .EQ. 0) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1)+1, hi(1)
                        efab(i,j,k,n) =
     &                       half*(cfab(i,j,k,n) + cfab(i-1,j,k,n))
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 1) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2)+1, hi(2)
                     do i = lo(1), hi(1)
                        efab(i,j,k,n) =
     &                       half*(cfab(i,j,k,n) + cfab(i,j-1,k,n))
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 2) then
            do n = 1,nc
               do k = lo(3)+1, hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        efab(i,j,k,n) =
     &                       half*(cfab(i,j,k,n) + cfab(i,j,k-1,n))
                     end do
                  end do
               end do
            end do
         end if
      else
         if (dir .EQ. 0) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1)+1, hi(1)
                        if((cfab(i,j,k,n) * cfab(i-1,j,k,n)) .gt.0.d0)
     &                       then
                           efab(i,j,k,n) =
     &                          2.d0*(cfab(i,j,k,n) * cfab(i-1,j,k,n))/
     &                          (cfab(i,j,k,n) + cfab(i-1,j,k,n))
                        else
                           efab(i,j,k,n)=0.d0
                        endif
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 1) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2)+1, hi(2)
                     do i = lo(1), hi(1)
                        if((cfab(i,j,k,n) * cfab(i,j-1,k,n)).gt.0.d0)
     &                       then
                           efab(i,j,k,n) =
     &                          2.d0*(cfab(i,j,k,n) * cfab(i,j-1,k,n))/
     &                          (cfab(i,j,k,n) + cfab(i,j-1,k,n))
                        else
                           efab(i,j,k,n)=0.d0
                        endif
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 2) then
            do n = 1,nc
               do k = lo(3)+1, hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        if((cfab(i,j,k,n) * cfab(i,j,k-1,n)).gt.0.d0)
     &                       then
                           efab(i,j,k,n) =
     &                          2.d0*(cfab(i,j,k,n) * cfab(i,j,k-1,n))/
     &                          (cfab(i,j,k,n) + cfab(i,j,k-1,n))
                        else
                           efab(i,j,k,n)=0.d0
                        endif
                     end do
                  end do
               end do
            end do
         end if
      end if
      end
c-----------------------------------------------------------------------
      subroutine EDGE_INTERP(flo, fhi, nc, ratio, dir,
     &     fine, fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2)
      implicit none
      integer flo(0:3-1), fhi(0:3-1), nc, ratio(0:3-1), dir
      integer fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2
      DOUBLE PRECISION
     &     fine(fine_l0:fine_h0,fine_l1:fine_h1,fine_l2:fine_h2,nc)
      integer i,j,k,ii,jj,kk,n,P,M,L,clo(0:3-1),chi(0:3-1)
      DOUBLE PRECISION val, df

c     Do linear in dir, pc transverse to dir, leave alone the fine values
c     lining up with coarse edges--assume these have been set to hold the 
c     values you want to interpolate to the rest.
      if (dir.eq.0) then
         do n=1,nc
            do k=flo(2),fhi(2),ratio(2)
               do j=flo(1),fhi(1),ratio(1)
                  do i=flo(0),fhi(0)-ratio(dir),ratio(0)
                     df = fine(i+ratio(dir),j,k,n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n)
     &                       + df*dble(M)/dble(ratio(dir))
                        do P=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                           do L=MAX(k,flo(2)),MIN(k+ratio(2)-1,fhi(2))
                              fine(i+M,P,L,n) = val
                           enddo
                        enddo
                     enddo                     
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do n=1,nc
            do k=flo(2),fhi(2),ratio(2)
               do j=flo(1),fhi(1)-ratio(dir),ratio(1)
                  do i=flo(0),fhi(0)
                     df = fine(i,j+ratio(dir),k,n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n)
     &                       + df*dble(M)/dble(ratio(dir))
                        do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                           do L=MAX(k,flo(2)),MIN(k+ratio(2)-1,fhi(2))
                              fine(P,j+M,L,n) = val
                           enddo
                        enddo
                     enddo                     
                  enddo
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do k=flo(2),fhi(2)-ratio(dir),ratio(2)
               do j=flo(1),fhi(1),ratio(1)
                  do i=flo(0),fhi(0),ratio(0)
                     df = fine(i,j,k+ratio(dir),n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n)
     &                       + df*dble(M)/dble(ratio(dir))
                        do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                           do L=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                              fine(P,L,k+M,n) = val
                           enddo
                        enddo
                     enddo                     
                  enddo
               enddo
            enddo
         enddo
      endif
      end
c-----------------------------------------------------------------------
      subroutine PC_EDGE_INTERP(lo, hi, nc, ratio, dir,
     &     crse, crse_l0, crse_l1, crse_l2, crse_h0, crse_h1, crse_h2,
     &     fine, fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2)
      implicit none
      integer lo(3),hi(3), nc, ratio(0:3-1), dir
      integer crse_l0, crse_l1, crse_l2, crse_h0, crse_h1, crse_h2
      integer fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2
      DOUBLE PRECISION
     &     crse(crse_l0:crse_h0,crse_l1:crse_h1,crse_l2:crse_h2,nc)
      DOUBLE PRECISION
     &     fine(fine_l0:fine_h0,fine_l1:fine_h1,fine_l2:fine_h2,nc)
      integer i,j,k,ii,jj,kk,n,L, P
      DOUBLE PRECISION val, dc

c     For edge-based data, fill fine values with piecewise-constant interp of coarse data.
c     Operate only on faces that overlap--ie, only fill the fine faces that make up each
c     coarse face, leave the in-between faces alone.
      if (dir.eq.0) then
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     do P=0,ratio(2)-1
                        do L=0,ratio(1)-1
                           fine(ii,jj+L,kk+P,n) = crse(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     do P=0,ratio(2)-1
                        do L=0,ratio(0)-1
                           fine(ii+L,jj,kk+P,n) = crse(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     do P=0,ratio(1)-1
                        do L=0,ratio(0)-1
                           fine(ii+L,jj+P,kk,n) = crse(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
      end
