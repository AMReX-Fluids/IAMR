
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PROJECTION_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 2

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define SMALL 1.0e-10
#else
#define SMALL 1.0d-10
#endif


module projection_2d_module
  
  implicit none

  private 

  public :: radmpyscal, radmpyvel, fort_raddiv, &
            anelcoeffmpy
  
contains

    subroutine radmpyscal(lo,hi,a,DIMS(a),domlo,domhi,r,rlo,rhi) &
         bind(C,name="radmpyscal")
!c 
!c     multiply A by Radius r
!c
!c 
!c     NOTE: THIS ROUTINE HAS BEEN MODIFIED SO THAT ALL VALUES
!c           OUTSIDE THE DOMAIN ARE SET TO ZERO
!c
      implicit none
      integer    rlo,rhi
      integer    DIMDEC(a)
      integer    domlo(2), domhi(2)
      integer    lo(SDIM),hi(SDIM)
      REAL_T     a(DIMV(a))
      REAL_T     r(rlo:rhi)
      REAL_T     bogus_value

      integer i, j

      do j = lo(2),hi(2)
         do i = lo(1),min(domhi(1),hi(1))
           a(i,j) = r(i)*a(i,j)
         end do
      end do

      !c     NOTE: We used to set these to bogus_value to be sure that we
      !c           didn't use them. But now in the divu routine in the F90
      !c           solvers we need to include these values in the stencil
      !c           because they might contain inflow values, for example, 
      !c           and the only test is on the BC for the pressure solve,
      !c           which doesn't differentiate between inflow, reflecting
      !c           and symmetry.

      if (lo(1) .lt. domlo(1)) then
         do j = lo(2),hi(2)
         do i = lo(1), domlo(1)-1
           a(i,j) = 0.d0
         end do
         end do
      end if

      if (hi(1) .gt. domhi(1)) then
         do j = lo(2),hi(2)
         do i = domhi(1)+1, hi(1)
           a(i,j) = 0.d0
         end do
         end do
      end if

    end subroutine radmpyscal

    subroutine radmpyvel(lo,hi,a,DIMS(a),domlo,domhi,r,rlo,rhi,ndim) &
         bind(C,name="radmpyvel")
!c 
!c     multiply A by Radius r
!c
!c 
!c     NOTE: THIS ROUTINE HAS BEEN MODIFIED SO THAT ALL VALUES
!c           OUTSIDE THE DOMAIN ARE SET TO ZERO
!c
      implicit none
      integer    ndim
      integer    rlo,rhi
      integer    DIMDEC(a)
      integer    domlo(2), domhi(2)
      integer    lo(SDIM),hi(SDIM)
      REAL_T     a(DIMV(a))
      REAL_T     r(rlo:rhi)

      integer    i,j 
      REAL_T     dr

      do j = lo(2),hi(2)
         do i = lo(1), min(hi(1),domhi(1))
           a(i,j) = r(i)*a(i,j)
         end do
      end do

!c     NOTE: We used to set these to bogus_value to be sure that we didn't use them.
!c           But now in the divu routine in the F90 solvers we need to include these
!c           values in the stencil because they might contain inflow values, for
!c           example, and the only test is on the BC for the pressure solve, which 
!c           doesn't differentiate between inflow, reflecting and symmetry.

      if (lo(1) .lt. domlo(1)) then
         do j = lo(2),hi(2)
         do i = lo(1), domlo(1)-1
           a(i,j) = 0.d0
         end do
         end do
      end if

!c     Here we only multiply a possibly inflow x-velocity from the hi-r side 
      if (ndim .eq. 0) then
         if (hi(1) .gt. domhi(1)) then
            dr = r(hi(1)) - r(hi(1)-1)
            do j = lo(2),hi(2)
            do i = domhi(1)+1, hi(1)
              a(i,j) = (r(domhi(1)) + (i-domhi(1))*dr - 0.5*dr)*a(i,j)
            end do
            end do
         end if
      else
         if (hi(1) .gt. domhi(1)) then
            do j = lo(2),hi(2)
            do i = domhi(1)+1, hi(1)
              a(i,j) = 0.d0
            end do
            end do
         end if
      end if

    end subroutine radmpyvel

    subroutine fort_raddiv(lo,hi,a,DIMS(a),domlo,domhi,r,rlo,rhi,bogus_value)&
         bind(C,name="fort_raddiv")
!c 
!c     divide A by Radius r
!c
!c 
!c     NOTE: THIS ROUTINE HAS BEEN MODIFIED SO THAT ALL VALUES
!c           OUTSIDE THE DOMAIN ARE SET TO BOGUS VALUE
!c
      implicit none
      integer    rlo,rhi
      integer    DIMDEC(a)
      integer    domlo(2), domhi(2)
      integer    lo(SDIM),hi(SDIM)
      REAL_T     a(DIMV(a))
      REAL_T     r(rlo:rhi)
      REAL_T     bogus_value

      integer i, j

      do j = lo(2), hi(2)
         do i = lo(1),min(domhi(1),hi(1))
           a(i,j) = a(i,j)/r(i)
         end do
      end do

      if (lo(1) .lt. domlo(1)) then
         do j = lo(2),hi(2)
         do i = lo(1), domlo(1)-1
           a(i,j) = bogus_value
         end do
         end do
      end if

      if (hi(1) .gt. domhi(1)) then
         do j = lo(2),hi(2)
         do i = domhi(1)+1, hi(1)
           a(i,j) = bogus_value
         end do
         end do
      end if

    end subroutine fort_raddiv

    subroutine anelcoeffmpy(lo,hi,a,DIMS(a),domlo,domhi,anel_coeff,&
         anel_lo,anel_hi,bogus_value,mult) bind(C,name="anelcoeffmpy")
!c 
!c     multiply A by the anelastic coefficient
!c
!c 
!c     NOTE: THIS ROUTINE HAS BEEN MODIFIED SO THAT ALL VALUES
!c           OUTSIDE THE DOMAIN ARE SET TO BOGUS VALUE
      !c-- only sets boudnary cells at top and bottom to Bogus Val,
      !c   boundary cells on the sides are left alone

      implicit none
      integer    lo(SDIM),hi(SDIM)
      integer    DIMDEC(a)
      integer    domlo(2), domhi(2)
      REAL_T     a(DIMV(a))
      integer    mult,anel_lo,anel_hi
      REAL_T     anel_coeff(anel_lo:anel_hi)
      REAL_T     bogus_value

      integer i, j

      integer jlo,jhi
      jlo = lo(2)
      jhi = hi(2)

      ! set the top and bottom ghost cells to bogus val
      if (lo(2) .lt. domlo(2)) then
         jlo = domlo(2)
         do j = lo(2), domlo(2)-1
         do i = lo(1),hi(1)
           a(i,j) = bogus_value
         end do
         end do
      end if
      if (hi(2) .gt. domhi(2)) then
         jhi = domhi(2)
         do j = domhi(2)+1, hi(2)
         do i = lo(1),hi(1)
           a(i,j) = bogus_value
         end do
         end do
      end if

      ! scale by anelastic coefficient
      if (mult .eq. 1) then
         do j = jlo,jhi
         do i = lo(1),hi(1)
           a(i,j) = a(i,j) * anel_coeff(j)
         end do
         end do
      else if (mult .eq. 0) then
         do j = jlo,jhi
         do i = lo(1),hi(1)
           a(i,j) = a(i,j) / anel_coeff(j)
         end do
         end do
      else 
         print *,'BOGUS MULT IN ANELCOEFFMULT ',mult
         stop
      end if

    end subroutine anelcoeffmpy

      subroutine hgn2c(&
          isrz,lrweighted, DIMS(nodedat), nodedat,&
          DIMS(ccdat), lo, hi, ccdat) bind(C,name="hgn2c")

!c     ----------------------------------------------------------
!c     HGN2C
!c     averages node centered data to cell centers for use in 
!c     holy grail projection

      implicit none
      integer isrz,lrweighted
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(ccdat)
      integer DIMDEC(nodedat)
      REAL_T  nodedat(DIMV(nodedat))
      REAL_T  ccdat(DIMV(ccdat))

      integer i,j

      if (ARG_H1(ccdat)   .lt. lo(1) .or. &
         ARG_L1(ccdat)   .gt. hi(1) .or. &
         ARG_H2(ccdat)   .lt. lo(2) .or. &
         ARG_L2(ccdat)   .gt. hi(2) .or. &
         ARG_H1(nodedat) .lt. lo(1)+1 .or.& 
         ARG_L1(nodedat) .gt. hi(1) .or. &
         ARG_H2(nodedat) .lt. lo(2)+1 .or. &
         ARG_L2(nodedat) .gt. hi(2) ) then 
        call bl_abort("FORT_HG_CELL_TO_NODE: bad index limits")
      end if

      if(isrz.eq.1.and.lrweighted.ne.1)then
        call bl_abort('FORT_HGN2C: isrz=1 and lrweighted!=1 not implmented')
      end if

      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          ccdat(i,j) = fourth*(nodedat(i,j)+nodedat(i+1,j)+&
                              nodedat(i,j+1)+nodedat(i+1,j+1))
        end do
      end do

    end subroutine hgn2c

      subroutine hgc2n(&
          nghost, DIMS(dat), dat, rcen,&
          DIMS(rhs), rhs,&
          domlo, domhi, dr, is_rz, imax)  bind(C,name="hgc2n")
!c
!c     ----------------------------------------------------------
!c     HGC2N
!c     averages cell centered data to nodes for use in 
!c     holy grail projection
!c     
!c     INPUTS / OUTPUTS:
!c     nghost      => indicates buffer of rhs that does not need values
!c     dat         => cell centered array to be averaged
!c     DIMS(dat)   => index limits of dat
!c     rcen        => r-coordinate cell centers if geoem is r-z; 
!c     otherwise, should be 1
!c     rhslo,rhshi => index extents of rhs
!c     rhs         <= node centered array with results
!c     ----------------------------------------------------------
!c 
      implicit none
      integer nghost 
      integer domlo(SDIM), domhi(SDIM)
      integer DIMDEC(dat)
      integer DIMDEC(rhs)
      REAL_T  dr
      REAL_T  rcen(DIM1(dat))
      REAL_T  dat(DIMV(dat))
      REAL_T  rhs(DIMV(rhs))
      integer is_rz, imax

      integer i, j
      REAL_T  rhi, rlo

#if BL_PRVERSION == 9
      REAL_T  factor
#endif

      if (ARG_L1(rhs)+1 .lt. ARG_L1(dat) .or. &
          ARG_H1(rhs)-1 .gt. ARG_H1(dat) .or.&
          ARG_L2(rhs)+1 .lt. ARG_L2(dat) .or. &
          ARG_H2(rhs)-1 .gt. ARG_H2(dat)) then
         call bl_abort("FORT_HG_CELL_TO_NODE: bad index limits")
      end if

      if (is_rz.ne.1) then
         do j=ARG_L2(rhs)+nghost,ARG_H2(rhs)-nghost
            do i=ARG_L1(rhs)+nghost,ARG_H1(rhs)-nghost
               rhs(i,j) = fourth*(dat(i-1,j-1)+dat(i-1,j)+&
                                 dat(i  ,j-1)+dat(i  ,j) )
            end do
         end do

      else

#if BL_PRVERSION == 9
         do j=ARG_L2(rhs)+nghost,ARG_H2(rhs)-nghost
            do i=ARG_L1(rhs)+nghost,ARG_H1(rhs)-nghost
               if (i .eq. imax) then
                  rhi = rcen(i-1)
               else 
                  rhi = rcen(i)
               end if
               if (i .eq. 0) then
                  rlo = rcen(i)
               else
                  rlo = rcen(i-1)
               end if

               rhs(i,j) = fourth*(rlo * (dat(i-1,j-1) + dat(i-1,j)) + &
                                 rhi * (dat(i  ,j-1) + dat(i  ,j)))
            end do
         end do
         factor = dr/24.0D0
         do j=ARG_L2(rhs)+nghost,ARG_H2(rhs)-nghost
            do i=ARG_L1(rhs)+nghost,ARG_H1(rhs)-nghost
               if (i .eq. imax) then
                  rhi = -one
               else 
                  rhi = one
               end if
               if (i .eq. 0) then
                  rlo = -one
               else
                  rlo = one
               end if
               rhs(i,j) = rhs(i,j) + factor *&
                   (rlo * (dat(i-1,j-1) + dat(i-1,j)) -&
                    rhi * (dat(i  ,j-1) + dat(i  ,j)))
            end do
         end do

#else
         do j=ARG_L2(rhs)+nghost,ARG_H2(rhs)-nghost
            do i=ARG_L1(rhs)+nghost,ARG_H1(rhs)-nghost
               if (i .eq. imax) then
                  rhi = rcen(i-1)
               else 
                  rhi = rcen(i)
               end if
               if (i .eq. 0) then
                  rlo = rcen(i)
               else 
                  rlo = rcen(i-1)
               end if
               rhs(i,j) = fourth*(rlo * (dat(i-1,j-1) + dat(i-1,j)) +&
                                 rhi * (dat(i  ,j-1) + dat(i  ,j)))
               if (i .eq. 0) rhs(i,j) = half * rhs(i,j)
            end do
         end do
#endif
      end if

    end subroutine hgc2n

  end module projection_2d_module
