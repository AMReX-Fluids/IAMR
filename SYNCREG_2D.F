c
c $Id: SYNCREG_2D.F,v 1.14 1999-08-01 15:12:38 almgren Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "SYNCREG_F.H"
#include "ArrayLim.H"

#define SDIM 2

c ::: -----------------------------------------------------------
c ::: coarsen fine grid node centered data along edge of node
c ::: centered fine grid array.
c :::
c ::: INPUTS/OUTPUTS:
c ::: crse      <=  node centered coarse data
c ::: DIMS(crse) => index limits for crse
c ::: fine       => node centered fine data
c ::: DIMS(fine) => index limits for fine
c ::: lo,hi      => node centered subregion of crse to define
c ::: dir        => index direction of normal (0 based)
c ::: ratio      => refinement ratio      
c ::: -----------------------------------------------------------

      subroutine FORT_SRCRSEREG(fine,DIMS(fine),crse,DIMS(crse),lo,hi,dir,
     &                          ratios,is_per)

      integer    DIMDEC(fine)
      integer    DIMDEC(crse)
      integer    lo(SDIM) 
      integer    hi(SDIM)
      integer    dir
      integer    ratios(0:1)
      REAL_T     fine(DIMV(fine))
      REAL_T     crse(DIMV(crse))
      
      integer    i, j, ic, jc
      integer    ratiox, ratioy
      integer    is_per

c     NOTE: the reason that the coefficients add up to 1/ratio and NOT to 1
c           is because the divergences and DGphi were computed using the
c           local dx, so in order to add the fine contribution to the
c           coarse contribution the fine contribution must be weighted
c           by dx_fine/dx_coarse = 1/ratio
      
      ratiox = ratios(0)
      ratioy = ratios(1)

c     NOTE: we halve fine on the corners since the corners are
c           counted in both directions; later we double them to restore
c           them to their previous values.

      if (dir .eq. 0) then

c        ::::: sum in j direction
         ic = lo(1)
         i = ratiox*ic

         j = ratioy*lo(2)
         fine(i,j) = half*fine(i,j)

         j = ratioy*hi(2)
         fine(i,j) = half*fine(i,j)

         if (ratioy .eq. 1) then
            do jc = lo(2), hi(2)
               crse(ic,jc) = fine(i,jc) / ratiox
            end do

         else if (ratioy .eq. 2) then

            do jc = lo(2), hi(2)
               j = ratioy*jc
               crse(ic,jc) =   half* fine(i,j) 
     &                     + fourth*(fine(i,j+1)+fine(i,j-1))
               crse(ic,jc) = crse(ic,jc) / ratiox
            end do

         else if (ratioy .eq. 4) then

            do jc = lo(2), hi(2)
               j = ratioy*jc
               crse(ic,jc) =  four* fine(i,j) 
     &                     + three*(fine(i,j+1)+fine(i,j-1))
     &                     +   two*(fine(i,j+2)+fine(i,j-2))
     &                     +       (fine(i,j+3)+fine(i,j-3))
               crse(ic,jc) = crse(ic,jc) / (sixteen * ratiox)
            end do
         else 
           print *,'BOGUS RATIOY IN SRCRSEREG ',ratioy
           stop
         end if

         j = ratioy*lo(2)
         fine(i,j) = two*fine(i,j)

         j = ratioy*hi(2)
         fine(i,j) = two*fine(i,j)

         if (is_per .eq. 1) then
            if (ratiox .eq. 2) then
              do jc = lo(2)+1, hi(2)-1
                j = ratioy*jc
                crse(ic,jc) =   crse(ic,jc) + 
     &              (fine(i+1,j)+fine(i-1,j)) / (four*ratioy)
              end do
            else if (ratiox .eq. 4) then
              do jc = lo(2)+1, hi(2)-1
                j = ratioy*jc
                crse(ic,jc) =  crse(ic,jc) + 
     &                      ( three*(fine(i+1,j)+fine(i-1,j))
     &                      +   two*(fine(i+2,j)+fine(i-2,j))
     &                      +       (fine(i+3,j)+fine(i-3,j)) ) / (sixteen*ratioy)
              end do
            endif
         endif

      else

c        ::::: sum in i direction
         jc = lo(2)
         j = ratioy*jc

         i = ratiox*lo(1)
         fine(i,j) = half*fine(i,j)

         i = ratiox*hi(1)
         fine(i,j) = half*fine(i,j)

         if (ratiox .eq. 1) then

            do ic = lo(1), hi(1)
               crse(ic,jc) = fine(ic,j) / ratiox
            end do

         else if (ratiox .eq. 2) then

            do ic = lo(1), hi(1)
               i = ratiox*ic
               crse(ic,jc) =   half* fine(i,j) 
     &                     + fourth*(fine(i+1,j)+fine(i-1,j))
               crse(ic,jc) = crse(ic,jc) / ratioy
            end do

         else if (ratiox .eq. 4) then

            do ic = lo(1), hi(1)
               i = ratiox*ic
               crse(ic,jc) =  four* fine(i,j) 
     &                     + three*(fine(i+1,j)+fine(i-1,j))
     &                     +   two*(fine(i+2,j)+fine(i-2,j))
     &                     +       (fine(i+3,j)+fine(i-3,j))
               crse(ic,jc) = crse(ic,jc) / (sixteen * ratioy)
            end do
         else 
           print *,'BOGUS RATIOX IN SRCRSEREG ',ratiox
           stop
         end if

         i = ratiox*lo(1)
         fine(i,j) = two*fine(i,j)

         i = ratiox*hi(1)
         fine(i,j) = two*fine(i,j)

         if (is_per .eq. 1) then
            if (ratioy .eq. 2) then
              do ic = lo(1)+1, hi(1)-1
                 i = ratiox*ic
                 crse(ic,jc) =   crse(ic,jc) + 
     &               (fine(i,j+1)+fine(i,j-1)) / (four*ratiox)
              end do
            else if (ratioy .eq. 4) then
              do ic = lo(1)+1, hi(1)-1
                 i = ratiox*ic
                 crse(ic,jc) =  crse(ic,jc) + 
     &                       ( three*(fine(i,j+1)+fine(i,j-1))
     &                       +   two*(fine(i,j+2)+fine(i,j-2))
     &                       +       (fine(i,j+3)+fine(i,j-3)) ) / (sixteen*ratiox)
              end do
            endif
         endif

      end if
      
      end


c ::: -----------------------------------------------------------
c ::: create mask at nodes using values at surrounding cells
c :::
c ::: INPUTS/OUTPUTS:
c ::: mask       <=  node-centered mask array
c ::: DIMS(mask)  => index limits for mask
c ::: cells       => cell-centered array (each value is 0. or 1.)
c ::: DIMS(cells) => index limits for cells
c ::: -----------------------------------------------------------

      subroutine FORT_MAKEMASK(mask,DIMS(mask),cells,DIMS(cells))

      integer    DIMDEC(mask)
      integer    DIMDEC(cells)
      REAL_T     mask(DIMV(mask))
      REAL_T     cells(DIMV(cells))
      
      integer    i, j

      do j = ARG_L2(mask), ARG_H2(mask)
        do i = ARG_L1(mask), ARG_H1(mask)
 
          mask(i,j) = cells(i,j  ) + cells(i-1,j  ) 
     $              + cells(i,j-1) + cells(i-1,j-1)

        end do
      end do
      
      end

c ::: -----------------------------------------------------------
c ::: modify mask from values 0. through 4. to values 0. or 1.
c :::
c ::: INPUTS/OUTPUTS:
c ::: mask       <=  node-centered mask array
c ::: DIMS(mask)  => index limits for mask
c ::: -----------------------------------------------------------

      subroutine FORT_CONVERTMASK(mask,DIMS(mask))

      integer    DIMDEC(mask)
      REAL_T     mask(DIMV(mask))
      
      integer    i, j

      do j = ARG_L2(mask), ARG_H2(mask)
        do i = ARG_L1(mask), ARG_H1(mask)
      
          if (mask(i,j) .gt. 3.5) then
            mask(i,j) = zero
          else
            mask(i,j) = one
          end if

        end do
      end do
      
      end
