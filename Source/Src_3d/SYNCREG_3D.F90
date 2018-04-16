
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <SYNCREG_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 3

module syncreg_3d_module

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
!c ::: DIMS(fine)=> index limits for fine
!c ::: lo,hi      => node centered subregion of crse to define
!c ::: dir        => index direction of normal (0 based)
!c ::: ratios     => IntVect refinement ratio      
!c ::: -----------------------------------------------------------

      subroutine srcrsereg(fine,DIMS(fine),crse,DIMS(crse),lo,hi,dir,&
                            ratios)bind(C,name="srcrsereg")

      implicit none
      integer    DIMDEC(fine)
      integer    DIMDEC(crse)
      integer    lo(SDIM)
      integer    hi(SDIM)
      integer    dir
      integer    ratios(0:SDIM-1)
      REAL_T     fine(DIMV(fine))
      REAL_T     crse(DIMV(crse))
      
      integer    i, j, k, ic, jc, kc
      integer    m, n
      integer    ratiox, ratioy, ratioz
      REAL_T     one28, fourthirds, threefourths
      REAL_T     coeff, denom

      parameter(one28        = eight * sixteen)
      parameter(fourthirds   = four  * third  )
      parameter(threefourths = three * fourth )
!c
!c     NOTE: the reason that the coefficients add up to 1/ratio_norm and NOT to 1
!c           is because the divergences and DGphi were computed using the
!c           local dx, so in order to add the fine contribution to the
!c           coarse contribution the fine contribution must be weighted
!c           by dx_fine/dx_coarse = 1/ratio
!c
      ratiox = ratios(0)
      ratioy = ratios(1)
      ratioz = ratios(2)
!c
!c     NOTE: we multiply fine by 1/2 on the edges and 1/3 on the corners
!c           since they are counted in all three directions; later we 
!c           multiply by 2 and 3, respectively, to restore them to their
!c           previous values.
!c     
      if (dir .eq. 0) then
       
!c        ::::: sum in j and k directions
         ic = lo(1)
         i = ratiox*ic

         do kc = lo(3), hi(3)
           do jc = lo(2), hi(2)
             crse(ic,jc,kc) = zero
           end do
         end do

         j = ratioy*lo(2)
         do k = fine_l3, fine_h3
           fine(i,j,k) = half*fine(i,j,k)
         enddo

         j = ratioy*hi(2)
         do k = fine_l3, fine_h3
           fine(i,j,k) = half*fine(i,j,k)
         enddo

         k = ratioz*lo(3)
         do j = fine_l2, fine_h2
           fine(i,j,k) = half*fine(i,j,k)
         enddo

         k = ratioz*hi(3)
         do j = fine_l2, fine_h2
           fine(i,j,k) = half*fine(i,j,k)
         enddo

         j = ratioy*lo(2)
         k = ratioz*lo(3)
         fine(i,j,k) = fourthirds*fine(i,j,k)

         j = ratioy*hi(2)
         k = ratioz*lo(3)
         fine(i,j,k) = fourthirds*fine(i,j,k)

         j = ratioy*lo(2)
         k = ratioz*hi(3)
         fine(i,j,k) = fourthirds*fine(i,j,k)

         j = ratioy*hi(2)
         k = ratioz*hi(3)
         fine(i,j,k) = fourthirds*fine(i,j,k)

         denom = one / (ratiox * ratioy**2 * ratioz**2)
         do n = 0, (ratioz-1)
           do m = 0, (ratioy-1)
             coeff = (ratioy - m) * (ratioz - n) * denom
             if (n .eq. 0) coeff = half * coeff
             if (m .eq. 0) coeff = half * coeff
             do kc = lo(3), hi(3)
               k = ratioz*kc
               do jc = lo(2), hi(2)
                 j = ratioy*jc
                 crse(ic,jc,kc) = crse(ic,jc,kc) + coeff * & 
                    ( fine(i,j+m,k+n) + fine(i,j-m,k+n) &
                     +fine(i,j+m,k-n) + fine(i,j-m,k-n) )
               end do
             end do
           end do
         end do

         j = ratioy*lo(2)
         do k = fine_l3, fine_h3
           fine(i,j,k) = two*fine(i,j,k)
         enddo

         j = ratioy*hi(2)
         do k = fine_l3, fine_h3
           fine(i,j,k) = two*fine(i,j,k)
         enddo

         k = ratioz*lo(3)
         do j = fine_l2, fine_h2
           fine(i,j,k) = two*fine(i,j,k)
         enddo

         k = ratioz*hi(3)
         do j = fine_l2, fine_h2
           fine(i,j,k) = two*fine(i,j,k)
         enddo

         j = ratioy*lo(2)
         k = ratioz*lo(3)
         fine(i,j,k) = threefourths*fine(i,j,k)

         j = ratioy*hi(2)
         k = ratioz*lo(3)
         fine(i,j,k) = threefourths*fine(i,j,k)

         j = ratioy*lo(2)
         k = ratioz*hi(3)
         fine(i,j,k) = threefourths*fine(i,j,k)

         j = ratioy*hi(2)
         k = ratioz*hi(3)
         fine(i,j,k) = threefourths*fine(i,j,k)
      
      else if (dir .eq. 1) then
!c
!c        ::::: sum in i and k directions
!c
         jc = lo(2)
         j = ratioy*jc

         do kc = lo(3), hi(3)
           do ic = lo(1), hi(1)
             crse(ic,jc,kc) = zero
           end do
         end do

         i = ratiox*lo(1)
         do k = fine_l3, fine_h3
           fine(i,j,k) = half*fine(i,j,k)
         enddo

         i = ratiox*hi(1)
         do k = fine_l3, fine_h3
           fine(i,j,k) = half*fine(i,j,k)
         enddo

         k = ratioz*lo(3)
         do i = fine_l1, fine_h1
           fine(i,j,k) = half*fine(i,j,k)
         enddo

         k = ratioz*hi(3)
         do i = fine_l1, fine_h1
           fine(i,j,k) = half*fine(i,j,k)
         enddo

         i = ratiox*lo(1)
         k = ratioz*lo(3)
         fine(i,j,k) = four*third*fine(i,j,k)

         i = ratiox*hi(1)
         k = ratioz*lo(3)
         fine(i,j,k) = four*third*fine(i,j,k)

         i = ratiox*lo(1)
         k = ratioz*hi(3)
         fine(i,j,k) = four*third*fine(i,j,k)

         i = ratiox*hi(1)
         k = ratioz*hi(3)
         fine(i,j,k) = four*third*fine(i,j,k)

         denom = one / (ratioy * ratiox**2 * ratioz**2)
         do n = 0, (ratioz-1)
           do m = 0, (ratiox-1)
             coeff = (ratiox - m) * (ratioz - n) * denom
             if (n .eq. 0) coeff = half * coeff
             if (m .eq. 0) coeff = half * coeff
             do kc = lo(3), hi(3)
               k = ratioz*kc
               do ic = lo(1), hi(1)
                 i = ratiox*ic
                 crse(ic,jc,kc) = crse(ic,jc,kc) + coeff * &
                    ( fine(i+m,j,k+n) + fine(i-m,j,k+n) &
                     +fine(i+m,j,k-n) + fine(i-m,j,k-n) )
               end do
             end do
           end do
         end do

         i = ratiox*lo(1)
         do k = fine_l3, fine_h3
           fine(i,j,k) = two*fine(i,j,k)
         enddo

         i = ratiox*hi(1)
         do k = fine_l3, fine_h3
           fine(i,j,k) = two*fine(i,j,k)
         enddo

         k = ratioz*lo(3)
         do i = fine_l1, fine_h1
           fine(i,j,k) = two*fine(i,j,k)
         enddo

         k = ratioz*hi(3)
         do i = fine_l1, fine_h1
           fine(i,j,k) = two*fine(i,j,k)
         enddo

         i = ratiox*lo(1)
         k = ratioz*lo(3)
         fine(i,j,k) = threefourths*fine(i,j,k)

         i = ratiox*hi(1)
         k = ratioz*lo(3)
         fine(i,j,k) = threefourths*fine(i,j,k)

         i = ratiox*lo(1)
         k = ratioz*hi(3)
         fine(i,j,k) = threefourths*fine(i,j,k)

         i = ratiox*hi(1)
         k = ratioz*hi(3)
         fine(i,j,k) = threefourths*fine(i,j,k)

       else if (dir .eq. 2) then

!c        ::::: sum in i and j directions
         kc = lo(3)
         k = ratioz*kc

         do jc = lo(2), hi(2)
           do ic = lo(1), hi(1)
             crse(ic,jc,kc) = zero
           end do
         end do

         i = ratiox*lo(1)
         do j = fine_l2, fine_h2
           fine(i,j,k) = half*fine(i,j,k)
         enddo

         i = ratiox*hi(1)
         do j = fine_l2, fine_h2
           fine(i,j,k) = half*fine(i,j,k)
         enddo

         j = ratioy*lo(2)
         do i = fine_l1, fine_h1
           fine(i,j,k) = half*fine(i,j,k)
         enddo

         j = ratioy*hi(2)
         do i = fine_l1, fine_h1
           fine(i,j,k) = half*fine(i,j,k)
         enddo

         i = ratiox*lo(1)
         j = ratioy*lo(2)
         fine(i,j,k) = fourthirds*fine(i,j,k)

         i = ratiox*hi(1)
         j = ratioy*lo(2)
         fine(i,j,k) = fourthirds*fine(i,j,k)

         i = ratiox*lo(1)
         j = ratioy*hi(2)
         fine(i,j,k) = fourthirds*fine(i,j,k)

         i = ratiox*hi(1)
         j = ratioy*hi(2)
         fine(i,j,k) = fourthirds*fine(i,j,k)

         denom = one / (ratioz * ratiox**2 * ratioy**2)
         do n = 0, (ratioy-1)
           do m = 0, (ratiox-1)
             coeff = (ratiox - m) * (ratioy - n) * denom
             if (n .eq. 0) coeff = half * coeff
             if (m .eq. 0) coeff = half * coeff
             do jc = lo(2), hi(2)
               j = ratioy*jc
               do ic = lo(1), hi(1)
                 i = ratiox*ic
                 crse(ic,jc,kc) = crse(ic,jc,kc) + coeff * &
                    ( fine(i+m,j+n,k) + fine(i-m,j+n,k) &
                     +fine(i+m,j-n,k) + fine(i-m,j-n,k) )
               end do
             end do
           end do
         end do

         i = ratiox*lo(1)
         do j = fine_l2, fine_h2
           fine(i,j,k) = two*fine(i,j,k)
         enddo

         i = ratiox*hi(1)
         do j = fine_l2, fine_h2
           fine(i,j,k) = two*fine(i,j,k)
         enddo

         j = ratioy*lo(2)
         do i = fine_l1, fine_h1
           fine(i,j,k) = two*fine(i,j,k)
         enddo

         j = ratioy*hi(2)
         do i = fine_l1, fine_h1
           fine(i,j,k) = two*fine(i,j,k)
         enddo

         i = ratiox*lo(1)
         j = ratioy*lo(2)
         fine(i,j,k) = threefourths*fine(i,j,k)

         i = ratiox*hi(1)
         j = ratioy*lo(2)
         fine(i,j,k) = threefourths*fine(i,j,k)

         i = ratiox*lo(1)
         j = ratioy*hi(2)
         fine(i,j,k) = threefourths*fine(i,j,k)

         i = ratiox*hi(1)
         j = ratioy*hi(2)
         fine(i,j,k) = threefourths*fine(i,j,k)

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

      subroutine makemask(mask,DIMS(mask),cells,DIMS(cells))&
                          bind(C,name="makemask")

      implicit none
      integer    DIMDEC(mask)
      integer    DIMDEC(cells)
      REAL_T     mask(DIMV(mask))
      REAL_T     cells(DIMV(cells))

      integer    i, j, k

      do k = ARG_L3(mask), ARG_H3(mask)
       do j = ARG_L2(mask), ARG_H2(mask)
        do i = ARG_L1(mask), ARG_H1(mask)
 
          mask(i,j,k) = cells(i,j  ,k  ) + cells(i-1,j  ,k  ) &
                     + cells(i,j-1,k  ) + cells(i-1,j-1,k  ) &
                     + cells(i,j  ,k-1) + cells(i-1,j  ,k-1) &
                     + cells(i,j-1,k-1) + cells(i-1,j-1,k-1)

        end do
       end do
      end do

      end subroutine makemask

!c ::: -----------------------------------------------------------
!c ::: modify mask from values 0. through 8. to values 0. or 1.
!c :::
!c ::: INPUTS/OUTPUTS:
!c ::: mask       <=  node-centered mask array
!c ::: DIMS(mask)  => index limits for mask
!c ::: -----------------------------------------------------------

      subroutine convertmask(mask,DIMS(mask))bind(C,name="convertmask")

      implicit none
      integer    DIMDEC(mask)
      REAL_T     mask(DIMV(mask))
 
      integer    i, j, k

      do k = ARG_L3(mask), ARG_H3(mask)
         do j = ARG_L2(mask), ARG_H2(mask)
            do i = ARG_L1(mask), ARG_H1(mask)
               
               if (mask(i,j,k) .gt. 7.5D0) then
                  mask(i,j,k) = zero
               else
                  mask(i,j,k) = one
               end if

            end do
         end do
      end do
 
      end subroutine convertmask

end module syncreg_3d_module

