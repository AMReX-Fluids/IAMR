
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif


#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_ArrayLim.H"
#include "SYNCREG_F.H"

#define SDIM 2


module syncreg_2d_module
  
  implicit none

  private 

  public :: srcrsereg, makemask, convertmask

contains
  
!c ::: -----------------------------------------------------------
!c ::: coarsen fine grid node centered data along edge of node
!c ::: centered fine grid array.
!c :::
!c ::: INPUTS/OUTPUTS:
!c ::: crse      <=  node centered coarse data
!c ::: DIMS(crse) => index limits for crse
!c ::: fine       => node centered fine data
!c ::: DIMS(fine) => index limits for fine
!c ::: lo,hi      => node centered subregion of crse to define
!c ::: dir        => index direction of normal (0 based)
!c ::: ratio      => refinement ratio      
!c ::: -----------------------------------------------------------

      subroutine srcrsereg(fine,DIMS(fine),crse,DIMS(crse),lo,hi,dir,&
                           ratios) bind(C,name="srcrsereg")

      implicit none
      integer    DIMDEC(fine)
      integer    DIMDEC(crse)
      integer    lo(SDIM) 
      integer    hi(SDIM)
      integer    dir
      integer    ratios(0:1)
      REAL_T     fine(DIMV(fine))
      REAL_T     crse(DIMV(crse))
      
      integer    i, j, ic, jc
      integer    m
      integer    ratiox, ratioy
      REAL_T     coeff, denom

!c     NOTE: the reason that the coefficients add up to 1/ratio and NOT to 1
!c           is because the divergences and DGphi were computed using the
!c           local dx, so in order to add the fine contribution to the
!c           coarse contribution the fine contribution must be weighted
!c           by dx_fine/dx_coarse = 1/ratio
      
      ratiox = ratios(0)
      ratioy = ratios(1)

!c     NOTE: we halve fine on the corners since the corners are
!c           counted in both directions; later we double them to restore
!c           them to their previous values.

      if (dir .eq. 0) then

!c        ::::: sum in j direction
         ic = lo(1)
         i = ratiox*ic

         do jc = lo(2), hi(2)
           crse(ic,jc) = zero
         end do

         j = ratioy*lo(2)
         fine(i,j) = half*fine(i,j)

         j = ratioy*hi(2)
         fine(i,j) = half*fine(i,j)

         denom = one / (ratiox * ratioy**2)
         do m = 0, (ratioy-1)
           coeff = (ratioy - m) * denom
           if (m .eq. 0) coeff = half * coeff
           do jc = lo(2), hi(2)
             j = ratioy*jc
             crse(ic,jc) = crse(ic,jc) + coeff *&
                 ( fine(i,j+m) + fine(i,j-m) )
           end do
         end do

         j = ratioy*lo(2)
         fine(i,j) = two*fine(i,j)

         j = ratioy*hi(2)
         fine(i,j) = two*fine(i,j)

      else

!c        ::::: sum in i direction
         jc = lo(2)
         j = ratioy*jc

         do ic = lo(1), hi(1)
           crse(ic,jc) = zero
         end do

         i = ratiox*lo(1)
         fine(i,j) = half*fine(i,j)

         i = ratiox*hi(1)
         fine(i,j) = half*fine(i,j)

         denom = one / (ratioy * ratiox**2)
         do m = 0, (ratiox-1)
           coeff = (ratiox - m) * denom
           if (m .eq. 0) coeff = half * coeff
           do ic = lo(1), hi(1)
             i = ratiox*ic
             crse(ic,jc) = crse(ic,jc) + coeff *&
                 ( fine(i+m,j) + fine(i-m,j) )
           end do
         end do

         i = ratiox*lo(1)
         fine(i,j) = two*fine(i,j)

         i = ratiox*hi(1)
         fine(i,j) = two*fine(i,j)

      end if
      
    end subroutine srcrsereg


!c ::: -----------------------------------------------------------
!c ::: create mask at nodes using values at surrounding cells
!c :::
!c ::: INPUTS/OUTPUTS:
!c ::: mask       <=  node-centered mask array
!c ::: DIMS(mask)  => index limits for mask
!c ::: cells       => cell-centered array (each value is 0. or 1.)
!c ::: DIMS(cells) => index limits for cells
!c ::: -----------------------------------------------------------

      subroutine makemask(mask,DIMS(mask),cells,DIMS(cells)) bind(C,name="makemask")

      implicit none
      integer    DIMDEC(mask)
      integer    DIMDEC(cells)
      REAL_T     mask(DIMV(mask))
      REAL_T     cells(DIMV(cells))
      
      integer    i, j

      do j = ARG_L2(mask), ARG_H2(mask)
        do i = ARG_L1(mask), ARG_H1(mask)
 
          mask(i,j) = cells(i,j  ) + cells(i-1,j  ) &
                   + cells(i,j-1) + cells(i-1,j-1)

        end do
      end do
      
    end subroutine makemask

!c ::: -----------------------------------------------------------
!c ::: modify mask from values 0. through 4. to values 0. or 1.
!c :::
!c ::: INPUTS/OUTPUTS:
!c ::: mask       <=  node-centered mask array
!c ::: DIMS(mask)  => index limits for mask
!c ::: -----------------------------------------------------------

      subroutine convertmask(mask,DIMS(mask)) bind(C,name="convertmask")

      implicit none
      integer    DIMDEC(mask)
      REAL_T     mask(DIMV(mask))
      
      integer    i, j

      do j = ARG_L2(mask), ARG_H2(mask)
        do i = ARG_L1(mask), ARG_H1(mask)
      
          if (mask(i,j) .gt. 3.5D0) then
            mask(i,j) = zero
          else
            mask(i,j) = one
          end if

        end do
      end do
      
    end subroutine convertmask
  end module syncreg_2d_module
