
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <MACOPERATOR_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 3

module macoperator_3d_module
  
  implicit none

  private 

  public :: maccoef, macrhs, macupdate, macsyncrhs
  
contains

!c :: ----------------------------------------------------------
!c :: MACCOEF
!c ::             Compute the coefficents for MAC solve
!c ::
!c :: INPUTS / OUTPUTS:
!c ::  cx,cy,cz    <=  edge coef arrays
!c ::  DIMS(cx)     => index limits for cx
!c ::  DIMS(cy)     => index limits for cy
!c ::  DIMS(cz)     => index limits for cz
!c ::  ax,ay,az     => edge based area arrays
!c ::  DIMS(ax)     => index limits for ax
!c ::  DIMS(ay)     => index limits for ay
!c ::  DIMS(az)     => index limits for az
!c ::  rho          => cell centered density array
!c ::  DIMS(rho)    => index limits of rho array
!c ::  lo,hi        => index limits for rhs
!c ::  dx           => cell size
!c :: ----------------------------------------------------------
!c ::

       subroutine maccoef (cx,DIMS(cx),cy,DIMS(cy),cz,DIMS(cz), &
                           ax,DIMS(ax),ay,DIMS(ay),az,DIMS(az), &
                           rho,DIMS(rho),lo,hi,vbxhi,dx) &
                           bind(C,name="maccoef")
       implicit none
       integer DIMDEC(cx)
       integer DIMDEC(cy)
       integer DIMDEC(cz)
       integer DIMDEC(ax)
       integer DIMDEC(ay)
       integer DIMDEC(az)
       integer DIMDEC(rho)
       integer lo(SDIM), hi(SDIM), vbxhi(SDIM)
       REAL_T  dx(SDIM)
       REAL_T  cx(DIMV(cx))
       REAL_T  cy(DIMV(cy))
       REAL_T  cz(DIMV(cz))
       REAL_T  ax(DIMV(ax))
       REAL_T  ay(DIMV(ay))
       REAL_T  az(DIMV(az))
       REAL_T  rho(DIMV(rho))

       integer i, j, k, hi1, hi2, hi3
       REAL_T  rhoavg

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (rho(i,j,k) .lt. zero) then
                   print *,' '
                   print *,'TESTING in MACCOEF '
                   print *,'RHO HAS GONE NEGATIVE AT ',i,j,k,rho(i,j,k)
                   call flush(6)
                   call bl_abort(" ")
                end if
             end do
          end do
       end do

       ! check to see if we're at a cc box boundary
       ! if so, we need to include 1 more point at high end because
       ! c is nodal in one dim
       if (hi(1) .eq. vbxhi(1)) then
          hi1 = hi(1)+1
       else
          hi1=hi(1)
       endif
       if (hi(2) .eq. vbxhi(2)) then
          hi2 = hi(2)+1
       else
          hi2=hi(2)
       endif
       if (hi(3) .eq. vbxhi(3)) then
          hi3 = hi(3)+1
       else
          hi3=hi(3)
       endif

!c
!c      ::::: finish coef in X direction (part 2)
!c

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi1
                rhoavg = half * ( rho(i,j,k) + rho(i-1,j,k) )
                cx(i,j,k) = dx(1) * ax(i,j,k) / rhoavg 
             end do
          end do
       end do

!c
!c      ::::: finish coef in Y direction (part 2)
!c

       do k = lo(3), hi(3)
          do j = lo(2), hi2
             do i = lo(1), hi(1)
                rhoavg = half * ( rho(i,j,k) + rho(i,j-1,k) )
                cy(i,j,k) = dx(2) * ay(i,j,k) / rhoavg 
             end do
          end do
       end do

!c
!c      ::::: finish coef in Z direction (part 2)
!c

       do k = lo(3), hi3
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rhoavg = half * ( rho(i,j,k) + rho(i,j,k-1) )
                cz(i,j,k) = dx(3) * az(i,j,k) / rhoavg 
             end do
          end do
       end do

       end subroutine maccoef

!c :: ----------------------------------------------------------
!c :: MACRHS
!c ::             Compute the RHS for MAC solve
!c ::
!c :: INPUTS / OUTPUTS:
!c ::  ux,uy,uz    <=  edge velocity arrays
!c ::  DIMS(ux)     => index limits for ux
!c ::  DIMS(uy)     => index limits for uy
!c ::  DIMS(uz)     => index limits for uz
!c ::  ax,ay,az     => edge based area arrays
!c ::  DIMS(ax)     => index limits for ax
!c ::  DIMS(ay)     => index limits for ay
!c ::  DIMS(az)     => index limits for az
!c ::  vol          => cell centered volume array
!c ::  DIMS(vol)    => index limits of vol array
!c ::  rhs         <=> cell centered rhs array
!c ::  DIMS(rhs)    => index limits of rhs array
!c ::  lo,hi        => index limits for rhs
!c ::  scale        => scale factor
!c :: ----------------------------------------------------------
!c ::
       subroutine macrhs (ux,DIMS(ux),uy,DIMS(uy),uz,DIMS(uz), &
                          ax,DIMS(ax),ay,DIMS(ay),az,DIMS(az), &
                          vol,DIMS(vol),rhs,DIMS(rhs),lo,hi,scale) &
                          bind(C,name="macrhs")
      
       implicit none
       integer DIMDEC(ux)
       integer DIMDEC(uy)
       integer DIMDEC(uz)
       integer DIMDEC(ax)
       integer DIMDEC(ay)
       integer DIMDEC(az)
       integer DIMDEC(vol)
       integer DIMDEC(rhs)
       integer lo(SDIM), hi(SDIM)
       REAL_T  scale
       REAL_T  ux(DIMV(ux))
       REAL_T  uy(DIMV(uy))
       REAL_T  uz(DIMV(uz))
       REAL_T  ax(DIMV(ax))
       REAL_T  ay(DIMV(ay))
       REAL_T  az(DIMV(az))
       REAL_T  vol(DIMV(vol))
       REAL_T  rhs(DIMV(rhs))

       integer i, j, k
       REAL_T  divu
!c
!c      ::::: rhs holds the divergence condition (possibly zero)
!c

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                divu = ax(i+1,j,k)*ux(i+1,j,k) - ax(i,j,k)*ux(i,j,k) &
                    + ay(i,j+1,k)*uy(i,j+1,k) - ay(i,j,k)*uy(i,j,k) &
                    + az(i,j,k+1)*uz(i,j,k+1) - az(i,j,k)*uz(i,j,k)
                rhs(i,j,k) = scale*(divu - vol(i,j,k)*rhs(i,j,k))
             end do
          end do
       end do

       end subroutine macrhs

!c :: ----------------------------------------------------------
!c :: MACUPDATE
!c ::             Compute the update to velocity field to
!c ::             make it divergence free
!c ::
!c :: INPUTS / OUTPUTS:
!c ::  ux,uy,uz    <=  edge based velocity arrays
!c ::  DIMS(ux)     => index limits for ux
!c ::  DIMS(uy)     => index limits for uy
!c ::  DIMS(uz)     => index limits for uz
!c ::  phi          => soln from MAC project
!c ::  DIMS(phi)    => index limits for phi
!c ::  rho          => density at time N
!c ::  DIMS(rho)    => index limits for rho
!c ::  dx           => cell size
!c ::  mult         => scalar multiplier
!c :: ----------------------------------------------------------
!c ::
       subroutine macupdate( &
          init, &
          ux,DIMS(ux),uy,DIMS(uy),uz,DIMS(uz), &
          phi,DIMS(phi),rho,DIMS(rho), &
          lo,hi,vbxhi,dx,mult)bind(C,name="macupdate")

       implicit none
       integer DIMDEC(ux)
       integer DIMDEC(uy)
       integer DIMDEC(uz)
       integer DIMDEC(phi)
       integer DIMDEC(rho)
       integer lo(SDIM), hi(SDIM), vbxhi(SDIM)
       REAL_T  dx(SDIM), mult
       REAL_T  ux(DIMV(ux))
       REAL_T  uy(DIMV(uy))
       REAL_T  uz(DIMV(uz))
       REAL_T  phi(DIMV(phi))
       REAL_T  rho(DIMV(rho))
       integer init

       integer i, j, k, hi1, hi2, hi3
       REAL_T  rhoavg, gp, idx, idy, idz
!c
!c     set gradient to zero if initializing
!c
       if ( init .eq. 1 ) then
          do k = ARG_L3(ux), ARG_H3(ux)
             do j = ARG_L2(ux), ARG_H2(ux)
                do i = ARG_L1(ux), ARG_H1(ux)
                   ux(i,j,k) = zero
                end do
             end do
          end do
          do k = ARG_L3(uy), ARG_H3(uy)
             do j = ARG_L2(uy), ARG_H2(uy)
                do i = ARG_L1(uy), ARG_H1(uy)
                   uy(i,j,k) = zero
                end do
             end do
          end do
          do k = ARG_L3(uz), ARG_H3(uz)
             do j = ARG_L2(uz), ARG_H2(uz)
                do i = ARG_L1(uz), ARG_H1(uz)
                   uz(i,j,k) = zero
                end do
             end do
          end do
       end if

       idx = 1.0d0 / dx(1)
       idy = 1.0d0 / dx(2)
       idz = 1.0d0 / dx(3)

       ! check to see if we're at a cc box boundary
       ! if so, we need to include 1 more point at high end because
       ! u is nodal in one dim
       if (hi(1) .eq. vbxhi(1)) then
          hi1 = hi(1)+1
       else
          hi1=hi(1)
       endif
       if (hi(2) .eq. vbxhi(2)) then
          hi2 = hi(2)+1
       else
          hi2=hi(2)
       endif
       if (hi(3) .eq. vbxhi(3)) then
          hi3 = hi(3)+1
       else
          hi3=hi(3)
       endif
!c
!c      compute x MAC gradient
!c
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi1
                rhoavg = half*(rho(i,j,k) + rho(i-1,j,k))
                gp = (phi(i,j,k)-phi(i-1,j,k))*idx
                ux(i,j,k) = ux(i,j,k) + mult * gp / rhoavg
             end do
          end do
       end do
!c     
!c     compute y mac gradient
!c
       do k = lo(3),hi(3)
          do j = lo(2),hi2
             do i = lo(1),hi(1)
                rhoavg = half*(rho(i,j,k) + rho(i,j-1,k))
                gp = (phi(i,j,k)-phi(i,j-1,k))*idy
                uy(i,j,k) = uy(i,j,k) + mult * gp / rhoavg
             end do
          end do
       end do
!c       
!c     compute z mac gradient
!c
       do k = lo(3),hi3
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                rhoavg = half*(rho(i,j,k) + rho(i,j,k-1))
                gp = (phi(i,j,k)-phi(i,j,k-1))*idz
                uz(i,j,k) = uz(i,j,k) + mult * gp / rhoavg
             end do
          end do
       end do

       end subroutine macupdate

!c :: ----------------------------------------------------------
!c :: MACSYNCRHS
!c ::        Modify the RHS for MAC SYNC solve
!c ::
!c :: INPUTS / OUTPUTS:
!c ::  rhs         <=  right hand side array
!c ::  lo,hi        => index limits for rhs
!c ::  vol          => cell centered volume array
!c ::  vlo,vhi      => index limits of vol array
!c ::  rhsscale     => const multiplier to rhs
!c :: ----------------------------------------------------------
!c ::
       subroutine macsyncrhs(rhs,DIMS(rhs),lo,hi, &
                                  vol,DIMS(vol),rhsscale) &
                                  bind(C,name="macsyncrhs")
                                  
       implicit none
       integer DIMDEC(rhs)
       integer DIMDEC(vol)
       integer lo(SDIM), hi(SDIM)
       REAL_T  rhsscale
       REAL_T  rhs(DIMV(rhs))
       REAL_T  vol(DIMV(vol))

       integer i, j, k
!c
!c      ::::: multiply by volume since reflux step (which computed rhs)
!c      ::::: divided by volume.
!c

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rhs(i,j,k) = rhsscale*vol(i,j,k)*rhs(i,j,k)
             end do
          end do
       end do
       
       end subroutine macsyncrhs

end module macoperator_3d_module
