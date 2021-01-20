#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif
  
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PROJOUTFLOWBC_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 2

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define SMALL 1.0e-10
#else
#define SMALL 1.0d-10
#endif


module projoutflowbc_2d_module
  
  implicit none

  private 

  public ::  extrap_proj, rhogbc, subtractavg

contains

!c *************************************************************************
!c ** EXTRAP_PROJ **
!c *************************************************************************

      subroutine extrap_proj(DIMS(u),u,DIMS(divu),divu,DIMS(rho),rho,&
          r_len,redge,DIMS(uExt),uExt,DIMS(divuExt),divuExt,&
          DIMS(rhoExt),rhoExt,lo,hi,face,zeroIt) bind(C,name="extrap_proj")

        implicit none

!c subtract divu_ave twice due to precision problems

      integer DIMDEC(u)
      integer DIMDEC(divu)
      integer DIMDEC(rho)
      integer DIMDEC(uExt)
      integer DIMDEC(divuExt)
      integer DIMDEC(rhoExt)
      integer r_len
      integer lo(SDIM),hi(SDIM)
      integer face
      REAL_T      u(DIMV(u),SDIM)
      REAL_T   divu(DIMV(divu))
      REAL_T    rho(DIMV(rho))
      REAL_T      uExt(DIMV(uExt),SDIM-1)
      REAL_T   divuExt(DIMV(divuExt))
      REAL_T   rhoExt(DIMV(rhoExt))
      REAL_T   redge(0:r_len-1)
      integer  zeroIt

!c local variables
      integer ics,ice,jcs,jce
      integer ife,jfe
      integer if,jf
      REAL_T divu_ave1,divu_ave2
      REAL_T max_divu
      REAL_T max_pert, small_pert
      parameter ( small_pert = SMALL)
      integer i,j

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3

      ics = ARG_L1(u)
      ice = ARG_H1(u)
      jcs = ARG_L2(u)
      jce = ARG_H2(u)

      ife = hi(1)
      jfe = hi(2)

      zeroIt = 0

      if (face .eq. XLO) then
         if=ife
         do j = jcs+1, jce-1
               uExt(j,if,1) = half*(three*   u(ice-1,j,2) -    u(ice,j,2))
            divuExt(j,if  ) = half*(three*divu(ice-1,j  ) - divu(ice,j))
             rhoExt(j,if  ) = half*(three* rho(ice-1,j  ) -  rho(ice,j))
         end do
         max_divu = ABS(divuExt(jcs+1,if))
         do j = jcs+1, jce-1
            max_divu = max(max_divu,ABS(divuExt(j,if)))
         end do
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(jcs+1,if))
         do j = jcs+1, jce-1
            max_pert = MAX(max_pert,ABS(divuExt(j,if)))
         end do
      else if (face .eq. YLO) then
         jf = jfe
         do i = ics+1, ice-1
               uExt(i,jf,1) = half*(three*   u(i,jce-1,1) -    u(i,jce,1))
            divuExt(i,jf  ) = half*(three*divu(i,jce-1  ) - divu(i,jce))
             rhoExt(i,jf  ) = half*(three* rho(i,jce-1  ) -  rho(i,jce))
         end do
         max_divu = ABS(divuExt(ics+1,jf))
         do i = ics+1, ice-1
            max_divu = max(max_divu,ABS(divuExt(i,jf)))
         end do
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(ics+1,jf))
         do i = ics+1, ice-1
            max_pert = MAX(max_pert,ABS(divuExt(i,jf)))
         end do
      else if (face .eq. XHI) then
         if = ife
         do j = jcs+1, jce-1
               uExt(j,if,1) = half*(three*   u(ics+1,j,2) -    u(ics,j,2))
            divuExt(j,if  ) = half*(three*divu(ics+1,j  ) - divu(ics,j))
             rhoExt(j,if  ) = half*(three* rho(ics+1,j  ) -  rho(ics,j))
         end do
         max_divu = ABS(divuExt(jcs+1,if))
         do j = jcs+1, jce-1
            max_divu = max(max_divu,ABS(divuExt(j,if)))
         end do
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(jcs+1,if))
         do j = jcs+1, jce-1
            max_pert = MAX(max_pert,ABS(divuExt(j,if)))
         end do
      else if (face .eq. YHI) then
         jf = jfe
         do i = ics+1, ice-1
               uExt(i,jf,1) = half*(three*   u(i,jcs+1,1) -    u(i,jcs,1))
            divuExt(i,jf  ) = half*(three*divu(i,jcs+1  ) - divu(i,jcs))
             rhoExt(i,jf  ) = half*(three* rho(i,jcs+1  ) -  rho(i,jcs))
         end do
         max_divu = ABS(divuExt(ics+1,jf))
         do i = ics+1, ice-1
            max_divu = max(max_divu,ABS(divuExt(i,jf)))
         end do
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave1,face)
         call subtractavg(DIMS(divuExt),divuExt,redge,r_len,lo,hi,divu_ave2,face)
         max_pert = ABS(divuExt(ics+1,jf))
         do i = ics+1, ice-1
            max_pert = MAX(max_pert,ABS(divuExt(i,jf)))
         end do
      endif

!c check to see if we should zero phi
      max_pert = max_pert/(ABS(divu_ave1+divu_ave2)+small_pert)
      if ((max_divu.eq.zero) .or. (max_pert.le.small_pert)) zeroIt = 1
#undef XLO
#undef YLO
#undef XHI
#undef YHI

    end subroutine extrap_proj
!c *************************************************************************
!c ** RHOGBC **
!c *************************************************************************

      subroutine rhogbc(rho,DIMS(rho),phi,DIMS(phi),&
           face,gravity,dx,domlo,domhi,lo_bc,hi_bc) bind(C,name="rhogbc")
!c
!c    Compute the contribution of gravity to the boundary conditions
!c      for phi at outflow faces only.
!c
      implicit none

      integer DIMDEC(rho)
      integer DIMDEC(phi)
      integer face
      integer domlo(2)
      integer domhi(2)
      integer lo_bc(2)
      integer hi_bc(2)
      REAL_T  rho(DIMV(rho))
      REAL_T  phi(DIMV(phi))
      REAL_T  dx(2)
      REAL_T  gravity
      
!c     Local variables
      integer i,j
      REAL_T rhog
      REAL_T rhoExt

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3

      if (face .eq. YLO .or. face .eq. YHI) &
        call bl_abort('SHOULDNT BE IN RHOGBC WITH FACE IN Y-DIR')

!c     Ok to only use low index of phi because phi is only one
!c        node wide in i-direction.
      i = ARG_L1(phi)

      if (face .eq. XLO) then

        rhog = zero
        do j = ARG_H2(phi)-1,ARG_L2(phi),-1
          rhoExt = half*(three*rho(i,j)-rho(i+1,j))
          rhog = rhog - gravity * rhoExt * dx(2)
          phi(i,j) = phi(i,j) + rhog
        end do

      else if (face .eq. XHI) then

        rhog = zero
        do j = ARG_H2(phi)-1,ARG_L2(phi),-1
          rhoExt  = half*(three*rho(i-1,j)-rho(i-2,j))
          rhog = rhog - gravity * rhoExt * dx(2)
          phi(i,j) = phi(i,j) + rhog
        end do

      endif

#undef XLO
#undef YLO
#undef XHI
#undef YHI

    end subroutine rhogbc

       
!c *************************************************************************
!c ** SUBTRACTAVG
!c *************************************************************************

      subroutine subtractavg(DIMS(divu),divu,redge,r_len,lo,hi,divu_ave,face)


        implicit none
      integer DIMDEC(divu)
      integer r_len
      integer lo(SDIM),hi(SDIM)
      REAL_T  redge(0:r_len-1)
      REAL_T divu(DIMV(divu))
      REAL_T divu_ave
      integer face

      integer i,j
      REAL_T rcen
      REAL_T vtot

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3
      
      divu_ave = zero
      vtot = zero

      if (face .eq. XLO .or. face .eq. XHI) then
         i = lo(1)
         do j=lo(2),hi(2)
            vtot = vtot+one
            divu_ave = divu_ave+divu(j,i)
         enddo
         divu_ave = divu_ave/vtot
         do j=lo(2),hi(2)
            divu(j,i) = divu(j,i) - divu_ave
         enddo
      elseif (face .eq. YLO .or. face .eq. YHI) then
         j = lo(2)
         do i=lo(1),hi(1)
            rcen = half*(redge(i)+redge(i+1))
            vtot = vtot+rcen
            divu_ave = divu_ave+rcen*divu(i,j)
         enddo
         divu_ave = divu_ave/vtot
         do i=lo(1),hi(1)
            divu(i,j) = divu(i,j) - divu_ave
         enddo
      else 
         print*, "bad value of face in subtractavg"
      endif

    end subroutine subtractavg
#undef XLO
#undef YLO
#undef XHI
#undef YHI

  end module projoutflowbc_2d_module
