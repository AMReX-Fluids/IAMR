c
c $Id: PROJECTION_3D.F,v 1.5 1998-03-31 17:26:40 lijewski Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "PROJECTION_F.H"
#include "ArrayLim.H"

#define SDIM 3

       subroutine FORT_FILTERP (
     &     pnew,scratch,DIMS(p),
     &     lo,hi,dx,fac)
c
c     ----------------------------------------------------------
c     filter the pressure coming out of the level projection
c     
c     pnew := pnew - fac * lap_7pt(pnew)
c     
c     INPUTS / OUTPUTS:
c     pnew       => 
c     pnew      <=> 
c     p_lo,p_hi  => index limits of pressure array
c     lo,hi      => index limits of cell centered interior
c     fac        => fudge factor for filtering
c     
c     ----------------------------------------------------------
c 
       integer    DIMDEC(p)
       REAL_T     pnew(DIMV(p))
       REAL_T     scratch(DIMV(p))
       integer    lo(3), hi(3)
       REAL_T     dx(3)
       REAL_T     fac

       integer i, j, k
       REAL_T     twelfth

       twelfth = one / 12.d0

       do k = lo(3)+1, hi(3)
          do j = lo(2)+1, hi(2)
             do i = lo(1)+1, hi(1)
                
                scratch(i,j,k) = six * pnew(i,j,k) - 
     $               (pnew(i-1,j,k) + pnew(i+1,j,k) + 
     $               pnew(i,j-1,k) + pnew(i,j+1,k) +
     $               pnew(i,j,k-1) + pnew(i,j,k+1) )
                scratch(i,j,k) = - twelfth * scratch(i,j,k) * fac

             end do
          end do
       end do

       do k = lo(3)+1, hi(3)
          do j = lo(2)+1, hi(2)
             do i = lo(1)+1, hi(1)

                pnew(i,j,k) = pnew(i,j,k) + scratch(i,j,k)

             end do
          end do
       end do
 
       end

       subroutine FORT_ACCEL_TO_VEL( lo, hi, 
     &     uold,DIMS(uold),
     &     dt,
     &     unew,DIMS(unew) )
c
c     This function converts unew into a velocity via
c     Unew = Uold + alpha*Unew
c
       integer    lo(SDIM), hi(SDIM)
       REAL_T     dt
       integer    DIMDEC(uold),DIMDEC(unew)
       REAL_T     uold(DIMV(uold),SDIM)
       REAL_T     unew(DIMV(unew),SDIM)

       integer i, j, k, n

       do n = 1, SDIM
          do k = lo(3)-1, hi(3)+1
             do j = lo(2)-1, hi(2)+1
                do i = lo(1)-1, hi(1)+1
                   unew(i,j,k,n) = uold(i,j,k,n) + dt*unew(i,j,k,n)
                end do
             end do
          end do
       end do

       end

      subroutine FORT_VEL_TO_ACCEL( lo, hi, 
     &     unew,DIMS(unew),
     &     uold,DIMS(uold),
     &     dt )
c     
c     This function converts unew into an acceleration
c     
      integer    lo(SDIM), hi(SDIM)
      REAL_T     dt
      integer    DIMDEC(unew),DIMDEC(uold)
      REAL_T     uold(DIMV(uold),SDIM)
      REAL_T     unew(DIMV(unew),SDIM)

      integer i, j, k, n

      do n = 1, SDIM
         do k = lo(3)-1, hi(3)+1
            do j = lo(2)-1, hi(2)+1
               do i = lo(1)-1, hi(1)+1
                  unew(i,j,k,n) = (unew(i,j,k,n)-uold(i,j,k,n))/dt
               end do
            end do
         end do
      end do

      end

      subroutine FORT_PROJ_UPDATE(
     &     boxlo, boxhi, nvar, ngrow,
     &     un, DIMS(un),
     &     alpha,
     &     uo, DIMS(uo) )
c     
c     This function updates un via un = un + alpha*uo
c     The loop bounds are determined in the C++
c     
      integer    boxlo(SDIM), boxhi(SDIM), nvar, ngrow
      REAL_T     alpha
      integer    DIMDEC(un),DIMDEC(uo)
      REAL_T     un(DIMV(un),nvar)
      REAL_T     uo(DIMV(uo),nvar)

      integer i, j, k, n

      do n = 1, nvar
         do k = boxlo(3), boxhi(3)
            do j = boxlo(2), boxhi(2)
               do i = boxlo(1), boxhi(1)
                  un(i,j,k,n) = un(i,j,k,n) + alpha*uo(i,j,k,n)
               end do
            end do
         end do
      end do

       end

       subroutine FORT_RADMPY(a,DIMS(grid),ng,r,nr,n)
c 
c     multiply A by a Radius r
c
       integer    ng,nr
       integer    n(SDIM)
       integer    DIMDEC(grid)
       REAL_T     a(ARG_L1(grid)-ng:ARG_H1(grid)+ng,
     $              ARG_L2(grid)-ng:ARG_H2(grid)+ng,
     $              ARG_L3(grid)-ng:ARG_H3(grid)+ng)
       REAL_T     r(ARG_L1(grid)-nr:ARG_H1(grid)+nr)

       integer i, j, k

       do k = ARG_L3(grid)-n(3), ARG_H3(grid)+n(3)
          do j = ARG_L2(grid)-n(2), ARG_H2(grid)+n(2)
             do i = ARG_L1(grid)-n(1), ARG_H1(grid)+n(1)
                a(i,j,k) = r(i)*a(i,j,k)
             end do
          end do
       end do

       end

       subroutine FORT_RADDIV(a,DIMS(grid),ng,r,nr,n)
c 
c     divide A by a Radius r
c
       integer    ng,nr
       integer    n(SDIM)
       integer    DIMDEC(grid)
       REAL_T     a(ARG_L1(grid)-ng:ARG_H1(grid)+ng,
     $              ARG_L2(grid)-ng:ARG_H2(grid)+ng,
     $              ARG_L3(grid)-ng:ARG_H3(grid)+ng)
       REAL_T     r(ARG_L1(grid)-nr:ARG_H1(grid)+nr)

       integer i, j, k

       do k = ARG_L3(grid)-n(3), ARG_H3(grid)+n(3)
          do j = ARG_L2(grid)-n(2), ARG_H2(grid)+n(2)
             do i = ARG_L1(grid)-n(1), ARG_H1(grid)+n(1)
                a(i,j,k) = a(i,j,k)/r(i)
             end do
          end do
       end do

       end
