c
c $Id: SYNCREG_3D.F,v 1.6 1998-10-06 22:53:31 lijewski Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "SYNCREG_F.H"
#include "ArrayLim.H"

#define SDIM 3

c ::: -----------------------------------------------------------
c ::: This routine will compute a node centered divergence
c ::: from a cell centered field
c :::
c ::: NOTE: node (i,j,k) is in lower left corner of cell (i,j,k)
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: u          => cell centered velocity field
c ::: DIMS(u)    => index limits for u
c ::: divu      <=  node centered divergence
c ::: DIMS(divu) => index limits for divu
c ::: lo,hi      => node centered subregion of divu to define
c ::: dx         => cell sizes      
c ::: mult       => scalar multiplicative factor      
c ::: -----------------------------------------------------------

      subroutine FORT_SRDIVU(u,DIMS(u),divu,DIMS(divu),lo,hi,dx,mult)

      integer    DIMDEC(u)
      integer    DIMDEC(divu)
      integer    lo(SDIM), hi(SDIM)
      REAL_T     dx(SDIM), mult
      REAL_T     u(DIMV(u),SDIM)
      REAL_T     divu(DIMV(divu))
      
      integer    i, j, k, lenx, leny, lenz
      REAL_T     ux, vy, wz

      integer    ii, jj, kk

c     NOTE: the factor of mult here accounts for the fact that one
c           takes ratio time steps on the fine grid for every one
c           timestep on the coarse grid.  mult = 1/ratio.
      
      lenx = hi(1)-lo(1)+1
      leny = hi(2)-lo(2)+1
      lenz = hi(3)-lo(3)+1

      if (lenx .ge. leny .and. lenx .ge. lenz) then
        do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ux = fourth/dx(1) * (
     $          (u(i,j,k,1) + u(i,j,k-1,1) + u(i,j-1,k,1) + u(i,j-1,k-1,1)) -
     $          (u(i-1,j,k,1) + u(i-1,j,k-1,1) + 
     $           u(i-1,j-1,k,1) + u(i-1,j-1,k-1,1)) )
               vy = fourth/dx(2) * (
     $          (u(i,j,k,2) + u(i,j,k-1,2) + u(i-1,j,k,2) + u(i-1,j,k-1,2)) -
     $          (u(i,j-1,k,2) + u(i,j-1,k-1,2) + 
     $           u(i-1,j-1,k,2) + u(i-1,j-1,k-1,2)) )
               wz = fourth/dx(3) * (
     $          (u(i,j,k,3) + u(i-1,j,k,3) + u(i,j-1,k,3) + u(i-1,j-1,k,3)) -
     $          (u(i,j,k-1,3) + u(i-1,j,k-1,3) + 
     $           u(i,j-1,k-1,3) + u(i-1,j-1,k-1,3)) )
               divu(i,j,k) = mult*(ux+vy+wz)
            end do
           end do
          end do
      else if (leny .ge. lenx .and. leny .ge. lenz) then 
         do k = lo(3), hi(3)
          do i = lo(1), hi(1)
            do j = lo(2), hi(2)
               ux = fourth/dx(1) * (
     $          (u(i,j,k,1) + u(i,j,k-1,1) + u(i,j-1,k,1) + u(i,j-1,k-1,1)) -
     $          (u(i-1,j,k,1) + u(i-1,j,k-1,1) + 
     $           u(i-1,j-1,k,1) + u(i-1,j-1,k-1,1)) )
               vy = fourth/dx(2) * (
     $          (u(i,j,k,2) + u(i,j,k-1,2) + u(i-1,j,k,2) + u(i-1,j,k-1,2)) -
     $          (u(i,j-1,k,2) + u(i,j-1,k-1,2) + 
     $           u(i-1,j-1,k,2) + u(i-1,j-1,k-1,2)) )
               wz = fourth/dx(3) * (
     $          (u(i,j,k,3) + u(i-1,j,k,3) + u(i,j-1,k,3) + u(i-1,j-1,k,3)) -
     $          (u(i,j,k-1,3) + u(i-1,j,k-1,3) + 
     $           u(i,j-1,k-1,3) + u(i-1,j-1,k-1,3)) )
               divu(i,j,k) = mult*(ux+vy+wz)
            end do
           end do
          end do
      else if (lenz .ge. lenx .and. lenz .ge. leny) then 
          do j = lo(2), hi(2)
           do i = lo(1), hi(1)
            do k = lo(3), hi(3)
               ux = fourth/dx(1) * (
     $          (u(i,j,k,1) + u(i,j,k-1,1) + u(i,j-1,k,1) + u(i,j-1,k-1,1)) -
     $          (u(i-1,j,k,1) + u(i-1,j,k-1,1) + 
     $           u(i-1,j-1,k,1) + u(i-1,j-1,k-1,1)) )
               vy = fourth/dx(2) * (
     $          (u(i,j,k,2) + u(i,j,k-1,2) + u(i-1,j,k,2) + u(i-1,j,k-1,2)) -
     $          (u(i,j-1,k,2) + u(i,j-1,k-1,2) + 
     $           u(i-1,j-1,k,2) + u(i-1,j-1,k-1,2)) )
               wz = fourth/dx(3) * (
     $          (u(i,j,k,3) + u(i-1,j,k,3) + u(i,j-1,k,3) + u(i-1,j-1,k,3)) -
     $          (u(i,j,k-1,3) + u(i-1,j,k-1,3) + 
     $           u(i,j-1,k-1,3) + u(i-1,j-1,k-1,3)) )
               divu(i,j,k) = mult*(ux+vy+wz)
            end do
           end do
          end do
      end if
      
      end

c ::: -----------------------------------------------------------
c ::: This routine will compute one version of DIV(sigma*GRAD(phi))
c ::: from phi and sigma.
c ::: This version uses GRAD-PERP.      
c :::
c ::: NOTE: dx == dy == dz is assumed
c :::       phi is node centered
c :::       sigma is cell centered      
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: phi        => node centered data
c ::: DIMS(phi)  => index limits for phi
c ::: sigma      => cell centered coeficient array
c ::: DIMS(sigma)=> index limits for sigma      
c ::: divgp     <=  node centered divergence of gradient
c ::: DIMS(divgp)=> index limits for divgp
c ::: lo,hi      => node centered subregion of divgp to define
c ::: dx         => cell sizes      
c ::: mult       => scalar multiplicative factor      
c ::: -----------------------------------------------------------

      subroutine FORT_SRDGPPHI(phi,DIMS(phi),sigma,DIMS(sigma),
     $                         divgp,DIMS(divgp),lo,hi,dx,mult)

      integer    DIMDEC(phi)
      integer    DIMDEC(sigma)
      integer    DIMDEC(divgp)
      integer    lo(SDIM), hi(SDIM)
      REAL_T     dx(SDIM), mult
      REAL_T     phi(DIMV(phi))
      REAL_T     sigma(DIMV(sigma))
      REAL_T     divgp(DIMV(divgp))
      
      integer    i, j, k, lenx, leny, lenz
      REAL_T     p1, p2, p3, p4, p5, p6, p7, p8, denom

c     NOTE: the factor of mult here accounts for the fact that one
c           takes ratio time steps on the fine grid for every one
c           timestep on the coarse grid.  mult = 1/ratio.
      
      denom = -mult/(sixteen*dx(1)**2)
      lenx = hi(1)-lo(1)+1
      leny = hi(2)-lo(2)+1
      lenz = hi(3)-lo(3)+1
      if (lenx .ge. leny) then
         do k = lo(3), hi(3)
          do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               p1 = sigma(i,j,k)*(
     &          five*(phi(i+1,j,k) + phi(i,j+1,k) + phi(i,j,k+1)) 
     &             - phi(i+1,j+1,k) - phi(i,j+1,k+1) - phi(i+1,j,k+1) 
     &             - three*phi(i+1,j+1,k+1) - nine*phi(i,j,k) )
               p2 = sigma(i-1,j,k)*(
     &          five*(phi(i-1,j,k) + phi(i,j+1,k) + phi(i,j,k+1)) 
     &             - phi(i-1,j+1,k) - phi(i,j+1,k+1) - phi(i-1,j,k+1) 
     &             - three*phi(i-1,j+1,k+1) - nine*phi(i,j,k) )
               p3 = sigma(i,j-1,k)*(
     &          five*(phi(i+1,j,k) + phi(i,j-1,k) + phi(i,j,k+1)) 
     &             - phi(i+1,j-1,k) - phi(i,j-1,k+1) - phi(i+1,j,k+1) 
     &             - three*phi(i+1,j-1,k+1) - nine*phi(i,j,k) )
               p4 = sigma(i-1,j-1,k)*(
     &          five*(phi(i-1,j,k) + phi(i,j-1,k) + phi(i,j,k+1)) 
     &             - phi(i-1,j-1,k) - phi(i,j-1,k+1) - phi(i-1,j,k+1) 
     &             - three*phi(i-1,j-1,k+1) - nine*phi(i,j,k) )
               p5 = sigma(i,j,k-1)*(
     &          five*(phi(i+1,j,k) + phi(i,j+1,k) + phi(i,j,k-1)) 
     &             - phi(i+1,j+1,k) - phi(i,j+1,k-1) - phi(i+1,j,k-1) 
     &             - three*phi(i+1,j+1,k-1) - nine*phi(i,j,k) )
               p6 = sigma(i-1,j,k-1)*(
     &          five*(phi(i-1,j,k) + phi(i,j+1,k) + phi(i,j,k-1)) 
     &             - phi(i-1,j+1,k) - phi(i,j+1,k-1) - phi(i-1,j,k-1) 
     &             - three*phi(i-1,j+1,k-1) - nine*phi(i,j,k) )
               p7 = sigma(i,j-1,k-1)*(
     &          five*(phi(i+1,j,k) + phi(i,j-1,k) + phi(i,j,k-1)) 
     &             - phi(i+1,j-1,k) - phi(i,j-1,k-1) - phi(i+1,j,k-1) 
     &             - three*phi(i+1,j-1,k-1) - nine*phi(i,j,k) )
               p8 = sigma(i-1,j-1,k-1)*(
     &          five*(phi(i-1,j,k) + phi(i,j-1,k) + phi(i,j,k-1)) 
     &             - phi(i-1,j-1,k) - phi(i,j-1,k-1) - phi(i-1,j,k-1) 
     &             - three*phi(i-1,j-1,k-1) - nine*phi(i,j,k) )
               divgp(i,j,k) = denom*(p1+p2+p3+p4+p5+p6+p7+p8)
            end do
          end do
         end do
      else if (leny .ge. lenx .and. leny .ge. lenz) then 
         do k = lo(3), hi(3)
          do i = lo(1), hi(1)
            do j = lo(2), hi(2)
               p1 = sigma(i,j,k)*(
     &          five*(phi(i+1,j,k) + phi(i,j+1,k) + phi(i,j,k+1)) 
     &             - phi(i+1,j+1,k) - phi(i,j+1,k+1) - phi(i+1,j,k+1) 
     &             - three*phi(i+1,j+1,k+1) - nine*phi(i,j,k) )
               p2 = sigma(i-1,j,k)*(
     &          five*(phi(i-1,j,k) + phi(i,j+1,k) + phi(i,j,k+1)) 
     &             - phi(i-1,j+1,k) - phi(i,j+1,k+1) - phi(i-1,j,k+1) 
     &             - three*phi(i-1,j+1,k+1) - nine*phi(i,j,k) )
               p3 = sigma(i,j-1,k)*(
     &          five*(phi(i+1,j,k) + phi(i,j-1,k) + phi(i,j,k+1)) 
     &             - phi(i+1,j-1,k) - phi(i,j-1,k+1) - phi(i+1,j,k+1) 
     &             - three*phi(i+1,j-1,k+1) - nine*phi(i,j,k) )
               p4 = sigma(i-1,j-1,k)*(
     &          five*(phi(i-1,j,k) + phi(i,j-1,k) + phi(i,j,k+1)) 
     &             - phi(i-1,j-1,k) - phi(i,j-1,k+1) - phi(i-1,j,k+1) 
     &             - three*phi(i-1,j-1,k+1) - nine*phi(i,j,k) )
               p5 = sigma(i,j,k-1)*(
     &          five*(phi(i+1,j,k) + phi(i,j+1,k) + phi(i,j,k-1)) 
     &             - phi(i+1,j+1,k) - phi(i,j+1,k-1) - phi(i+1,j,k-1) 
     &             - three*phi(i+1,j+1,k-1) - nine*phi(i,j,k) )
               p6 = sigma(i-1,j,k-1)*(
     &          five*(phi(i-1,j,k) + phi(i,j+1,k) + phi(i,j,k-1)) 
     &             - phi(i-1,j+1,k) - phi(i,j+1,k-1) - phi(i-1,j,k-1) 
     &             - three*phi(i-1,j+1,k-1) - nine*phi(i,j,k) )
               p7 = sigma(i,j-1,k-1)*(
     &          five*(phi(i+1,j,k) + phi(i,j-1,k) + phi(i,j,k-1)) 
     &             - phi(i+1,j-1,k) - phi(i,j-1,k-1) - phi(i+1,j,k-1) 
     &             - three*phi(i+1,j-1,k-1) - nine*phi(i,j,k) )
               p8 = sigma(i-1,j-1,k-1)*(
     &          five*(phi(i-1,j,k) + phi(i,j-1,k) + phi(i,j,k-1)) 
     &             - phi(i-1,j-1,k) - phi(i,j-1,k-1) - phi(i-1,j,k-1) 
     &             - three*phi(i-1,j-1,k-1) - nine*phi(i,j,k) )
               divgp(i,j,k) = denom*(p1+p2+p3+p4+p5+p6+p7+p8)
            end do
          end do
         end do
      else if (lenz .ge. lenx .and. lenz .ge. leny) then 
         do j = lo(2), hi(2)
          do i = lo(1), hi(1)
           do k = lo(3), hi(3)
               p1 = sigma(i,j,k)*(
     &          five*(phi(i+1,j,k) + phi(i,j+1,k) + phi(i,j,k+1)) 
     &             - phi(i+1,j+1,k) - phi(i,j+1,k+1) - phi(i+1,j,k+1) 
     &             - three*phi(i+1,j+1,k+1) - nine*phi(i,j,k) )
               p2 = sigma(i-1,j,k)*(
     &          five*(phi(i-1,j,k) + phi(i,j+1,k) + phi(i,j,k+1)) 
     &             - phi(i-1,j+1,k) - phi(i,j+1,k+1) - phi(i-1,j,k+1) 
     &             - three*phi(i-1,j+1,k+1) - nine*phi(i,j,k) )
               p3 = sigma(i,j-1,k)*(
     &          five*(phi(i+1,j,k) + phi(i,j-1,k) + phi(i,j,k+1)) 
     &             - phi(i+1,j-1,k) - phi(i,j-1,k+1) - phi(i+1,j,k+1) 
     &             - three*phi(i+1,j-1,k+1) - nine*phi(i,j,k) )
               p4 = sigma(i-1,j-1,k)*(
     &          five*(phi(i-1,j,k) + phi(i,j-1,k) + phi(i,j,k+1)) 
     &             - phi(i-1,j-1,k) - phi(i,j-1,k+1) - phi(i-1,j,k+1) 
     &             - three*phi(i-1,j-1,k+1) - nine*phi(i,j,k) )
               p5 = sigma(i,j,k-1)*(
     &          five*(phi(i+1,j,k) + phi(i,j+1,k) + phi(i,j,k-1)) 
     &             - phi(i+1,j+1,k) - phi(i,j+1,k-1) - phi(i+1,j,k-1) 
     &             - three*phi(i+1,j+1,k-1) - nine*phi(i,j,k) )
               p6 = sigma(i-1,j,k-1)*(
     &          five*(phi(i-1,j,k) + phi(i,j+1,k) + phi(i,j,k-1)) 
     &             - phi(i-1,j+1,k) - phi(i,j+1,k-1) - phi(i-1,j,k-1) 
     &             - three*phi(i-1,j+1,k-1) - nine*phi(i,j,k) )
               p7 = sigma(i,j-1,k-1)*(
     &          five*(phi(i+1,j,k) + phi(i,j-1,k) + phi(i,j,k-1)) 
     &             - phi(i+1,j-1,k) - phi(i,j-1,k-1) - phi(i+1,j,k-1) 
     &             - three*phi(i+1,j-1,k-1) - nine*phi(i,j,k) )
               p8 = sigma(i-1,j-1,k-1)*(
     &          five*(phi(i-1,j,k) + phi(i,j-1,k) + phi(i,j,k-1)) 
     &             - phi(i-1,j-1,k) - phi(i,j-1,k-1) - phi(i-1,j,k-1) 
     &             - three*phi(i-1,j-1,k-1) - nine*phi(i,j,k) )
               divgp(i,j,k) = denom*(p1+p2+p3+p4+p5+p6+p7+p8)
            end do
          end do
         end do
      end if
      
      end

c ::: -----------------------------------------------------------
c ::: This routine will compute one version of DIV(sigma*GRAD(phi))
c ::: from phi and sigma.
c ::: This is the non-PERP version.
c :::
c ::: NOTE: dx == dy == dz is assumed
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: phi        => node centered data
c ::: DIMS(phi)  => index limits for phi
c ::: sigma      => cell centered coeficient array
c ::: DIMS(sigma)=> index limits for sigma      
c ::: divgp     <=  node centered divergence of gradient
c ::: DIMS(divgp)=> index limits for divgp
c ::: lo,hi      => node centered subregion of divgp to define
c ::: dx         => cell sizes      
c ::: mult       => scalar multiplicative factor      
c ::: -----------------------------------------------------------

      subroutine FORT_SRDGPHI(phi,DIMS(phi),sigma,DIMS(sigma),
     $                        divgp,DIMS(divgp),lo,hi,dx,mult)

      integer    DIMDEC(phi)
      integer    DIMDEC(sigma)
      integer    DIMDEC(divgp)
      integer    lo(SDIM), hi(SDIM)
      REAL_T     dx(SDIM), mult
      REAL_T     phi(DIMV(phi))
      REAL_T     sigma(DIMV(sigma))
      REAL_T     divgp(DIMV(divgp))
      
      integer    i, j, k
      REAL_T     p1, p2, p3, p4, p5, p6
      REAL_T     denomx, denomy, denomz

c     NOTE: the factor of mult here accounts for the fact that one
c           takes ratio time steps on the fine grid for every one
c           timestep on the coarse grid.  mult = 1/ratio.
      
      denomx = -mult/(four*dx(1)**2)
      denomy = -mult/(four*dx(2)**2)
      denomz = -mult/(four*dx(3)**2)

      do k = lo(3), hi(3)
       do j = lo(2), hi(2)
        do i = lo(1), hi(1)

            p1 = (sigma(i  ,j  ,k  ) + sigma(i  ,j-1,k  ) +
     $            sigma(i  ,j  ,k-1) + sigma(i  ,j-1,k-1) ) * 
     $           (phi(i+1,j  ,k  ) - phi(i,j,k) )

            p2 = (sigma(i-1,j  ,k  ) + sigma(i-1,j-1,k  ) +
     $            sigma(i-1,j  ,k-1) + sigma(i-1,j-1,k-1) ) * 
     $           (phi(i-1,j  ,k  ) - phi(i,j,k) )

            p3 = (sigma(i  ,j  ,k  ) + sigma(i-1,j  ,k  ) +
     $            sigma(i  ,j  ,k-1) + sigma(i-1,j  ,k-1) ) * 
     $           (phi(i  ,j+1,k  ) - phi(i,j,k) )

            p4 = (sigma(i  ,j-1,k  ) + sigma(i-1,j-1,k  ) +
     $            sigma(i  ,j-1,k-1) + sigma(i-1,j-1,k-1) ) * 
     $           (phi(i  ,j-1,k  ) - phi(i,j,k) )

            p5 = (sigma(i  ,j  ,k  ) + sigma(i-1,j  ,k  ) +
     $            sigma(i  ,j-1,k  ) + sigma(i-1,j-1,k  ) ) * 
     $           (phi(i  ,j  ,k+1) - phi(i,j,k) )

            p6 = (sigma(i  ,j  ,k-1) + sigma(i-1,j  ,k-1) +
     $            sigma(i  ,j-1,k-1) + sigma(i-1,j-1,k-1) ) * 
     $           (phi(i  ,j  ,k-1) - phi(i,j,k) )

            divgp(i,j,k) = denomx*(p1+p2) + denomy*(p3+p4) + denomz*(p5+p6)

        end do
       end do
      end do

      end

c ::: -----------------------------------------------------------
c ::: coarsen fine grid node centered data along edge of node
c ::: centered fine grid array.
c :::
c ::: INPUTS/OUTPUTS:
c ::: crse      <=  node centered coarse data
c ::: DIMS(crse) => index limits for crse
c ::: fine       => node centered fine data
c ::: DIMS(fine)=> index limits for fine
c ::: lo,hi      => node centered subregion of crse to define
c ::: dir        => index direction of normal (0 based)
c ::: ratios     => IntVect refinement ratio      
c ::: -----------------------------------------------------------

      subroutine FORT_SRCRSEREG(fine,DIMS(fine),crse,DIMS(crse),lo,hi,dir,
     &                          ratios)

      integer    DIMDEC(fine)
      integer    DIMDEC(crse)
      integer    lo(SDIM)
      integer    hi(SDIM)
      integer    dir
      integer    ratios(0:SDIM-1)
      REAL_T     fine(DIMV(fine))
      REAL_T     crse(DIMV(crse))
      
      integer    i, j, k, ic, jc, kc
      integer    ratiox, ratioy, ratioz
      REAL_T     one28

      one28 = eight * sixteen

c     NOTE: the reason that the coefficients add up to 1/ratio_norm and NOT to 1
c           is because the divergences and DGphi were computed using the
c           local dx, so in order to add the fine contribution to the
c           coarse contribution the fine contribution must be weighted
c           by dx_fine/dx_coarse = 1/ratio

      ratiox = ratios(0)
      ratioy = ratios(1)
      ratioz = ratios(2)
      
      if (dir .eq. 0) then
       
c        ::::: sum in j and k directions
         ic = lo(1)
         i = ratiox*ic

         if (ratioy .eq. 1 .and. ratioz .eq. 1) then

           do kc = lo(3), hi(3)
             do jc = lo(2), hi(2)
               crse(ic,jc,kc) = fine(i,jc,kc) / ratiox
            end do
           end do

         else if (ratioy .eq. 1 .and. ratioz .eq. 2) then

           do kc = lo(3), hi(3)
             k = ratioz*kc
             do jc = lo(2), hi(2)
               j = ratioy*jc
               crse(ic,jc,kc) =   half* fine(i,j,k)
     &                        + fourth*(fine(i,j,k+1)+fine(i,j,k-1))
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratiox
            end do
           end do

         else if (ratioy .eq. 1 .and. ratioz .eq. 4) then

           do kc = lo(3), hi(3)
             k = ratioz*kc
             do jc = lo(2), hi(2)
               j = ratioy*jc
               crse(ic,jc,kc) = (four* fine(i,j,k)
     &                        + three*(fine(i,j,k+1)+fine(i,j,k-1))
     &                        +   two*(fine(i,j,k+2)+fine(i,j,k-2))
     &                        +       (fine(i,j,k+3)+fine(i,j,k-3)))/sixteen
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratiox
            end do
           end do

         else if (ratioy .eq. 2 .and. ratioz .eq. 1) then

           do kc = lo(3), hi(3)
             k = ratioz*kc
             do jc = lo(2), hi(2)
               j = ratioy*jc
               crse(ic,jc,kc) =   half* fine(i,j  ,k)
     &                        + fourth*(fine(i,j+1,k)+fine(i,j-1,k))
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratiox
            end do
           end do

         else if (ratioy .eq. 4 .and. ratioz .eq. 1) then

           do kc = lo(3), hi(3)
             k = ratioz*kc
             do jc = lo(2), hi(2)
               j = ratioy*jc
               crse(ic,jc,kc) = (four* fine(i,j,k)
     &                        + three*(fine(i,j+1,k)+fine(i,j-1,k))
     &                        +   two*(fine(i,j+2,k)+fine(i,j-2,k))
     &                        +       (fine(i,j+3,k)+fine(i,j-3,k)))/sixteen
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratiox
            end do
           end do

         else if (ratioy .eq. 2 .and. ratioz .eq. 2) then

           do kc = lo(3), hi(3)
             k = ratioz*kc
             do jc = lo(2), hi(2)
               j = ratioy*jc
               crse(ic,jc,kc) = half*half*fine(i,j,k) 
     &              + half  *fourth*(fine(i,j+1,k  ) + fine(i,j-1,k  )
     &              +                fine(i,j  ,k+1) + fine(i,j  ,k-1))
     &              + fourth*fourth*(fine(i,j+1,k+1) + fine(i,j-1,k+1)
     &              +                fine(i,j+1,k-1) + fine(i,j-1,k-1)) 
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratiox
            end do
           end do

         else if (ratioy .eq. 2 .and. ratioz .eq. 4) then

           do kc = lo(3), hi(3)
             k = ratioz*kc
             do jc = lo(2), hi(2)
               j = ratioy*jc
               crse(ic,jc,kc) = ( four*four*fine(i,j,k)  +
     &              four*three*(fine(i,j,k+1) + fine(i,j,k-1)) +
     &              four*two  *(fine(i,j,k+2) + fine(i,j,k-2)) +
     &              four      *(fine(i,j,k+3) + fine(i,j,k-3)) +
     &              two *four *(fine(i,j+1,k) + fine(i,j-1,k)) +
     &              two *three*(fine(i,j+1,k+1) + fine(i,j-1,k+1)+
     &                          fine(i,j+1,k-1) + fine(i,j-1,k-1)) +
     &              two *  two*(fine(i,j+1,k+2) + fine(i,j-1,k+2)+
     &                          fine(i,j+1,k-2) + fine(i,j-1,k-2)) +
     &              two       *(fine(i,j+1,k+3) + fine(i,j-1,k+3)+
     &                          fine(i,j+1,k-3) + fine(i,j-1,k-3)) ) / one28
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratiox
             end do
           end do

         else if (ratioy .eq. 4 .and. ratioz .eq. 2) then

           do kc = lo(3), hi(3)
             k = ratioz*kc
             do jc = lo(2), hi(2)
               j = ratioy*jc
               crse(ic,jc,kc) = ( four*four*fine(i,j,k)  +
     &              four*three*(fine(i,j+1,k) + fine(i,j-1,k)) +
     &              four*two  *(fine(i,j+2,k) + fine(i,j-2,k)) +
     &              four      *(fine(i,j+3,k) + fine(i,j-3,k)) +
     &              two *four *(fine(i,j,k+1) + fine(i,j,k-1)) +
     &              two *three*(fine(i,j+1,k+1) + fine(i,j+1,k-1)+
     &                          fine(i,j-1,k+1) + fine(i,j-1,k-1)) +
     &              two *  two*(fine(i,j+2,k+1) + fine(i,j+2,k-1)+
     &                          fine(i,j-2,k+1) + fine(i,j-2,k-1)) +
     &              two       *(fine(i,j+3,k+1) + fine(i,j+3,k-1)+
     &                          fine(i,j-3,k+1) + fine(i,j-3,k-1)) ) / one28
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratiox
             end do
           end do

         else if (ratioy .eq. 4 .and. ratioz .eq. 4) then

           do kc = lo(3), hi(3)
             k = ratioz*kc
             do jc = lo(2), hi(2)
               j = ratioy*jc
               crse(ic,jc,kc) = four*four*fine(i,j,k)  +
     &              four *three*(fine(i,j+1,k  ) + fine(i,j-1,k  )+
     &                           fine(i,j  ,k+1) + fine(i,j  ,k-1)) +
     &              four *two  *(fine(i,j+2,k  ) + fine(i,j-2,k  )+
     &                           fine(i,j  ,k+2) + fine(i,j  ,k-2)) +
     &              four       *(fine(i,j+3,k  ) + fine(i,j-3,k  )+
     &                           fine(i,j  ,k+3) + fine(i,j  ,k-3)) +
     &              three*three*(fine(i,j+1,k+1) + fine(i,j-1,k+1)+
     &                           fine(i,j+1,k-1) + fine(i,j-1,k-1)) +
     &              three*two  *(fine(i,j+1,k+2) + fine(i,j-1,k+2)+
     &                           fine(i,j+1,k-2) + fine(i,j-1,k-2)+
     &                           fine(i,j+2,k+1) + fine(i,j-2,k+1)+
     &                           fine(i,j+2,k-1) + fine(i,j-2,k-1)) +
     &              three      *(fine(i,j+1,k+3) + fine(i,j-1,k+3)+
     &                           fine(i,j+1,k-3) + fine(i,j-1,k-3)+
     &                           fine(i,j+3,k+1) + fine(i,j-3,k+1)+
     &                           fine(i,j+3,k-1) + fine(i,j-3,k-1)) +
     &              two  *two  *(fine(i,j+2,k+2) + fine(i,j-2,k+2)+
     &                           fine(i,j+2,k-2) + fine(i,j-2,k-2)) +
     &              two        *(fine(i,j+2,k+3) + fine(i,j-2,k+3)+
     &                           fine(i,j+2,k-3) + fine(i,j-2,k-3)+
     &                           fine(i,j+3,k+2) + fine(i,j-3,k+2)+
     &                           fine(i,j+3,k-2) + fine(i,j-3,k-2)) +
     &                          (fine(i,j+3,k+3) + fine(i,j-3,k+3)+
     &                           fine(i,j+3,k-3) + fine(i,j-3,k-3)) 
               crse(ic,jc,kc) = crse(ic,jc,kc) / (sixteen*sixteen*ratiox)
             end do
           end do
         else
           print *,'BOGUS RATIOY AND/OR RATIOZ IN SRCRSEREG ',ratioy,ratioz
           stop
         end if
      
      else if (dir .eq. 1) then
c
c        ::::: sum in i and k directions
c
         jc = lo(2)
         j = ratioy*jc

         if (ratiox .eq. 1 .and. ratioz .eq. 1) then

           do kc = lo(3), hi(3)
             do ic = lo(1), hi(1)
               crse(ic,jc,kc) = fine(ic,j,kc) / ratioy
            end do
           end do

         else if (ratiox .eq. 1 .and. ratioz .eq. 2) then

           do kc = lo(3), hi(3)
             k = ratioz*kc
             do ic = lo(1), hi(1)
               i = ratiox*ic
               crse(ic,jc,kc) =   half* fine(i,j,k)
     &                        + fourth*(fine(i,j,k+1)+fine(i,j,k-1))
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratioy
            end do
           end do

         else if (ratiox .eq. 1 .and. ratioz .eq. 4) then

           do kc = lo(3), hi(3)
             k = ratioz*kc
             do ic = lo(1), hi(1)
               i = ratiox*ic
               crse(ic,jc,kc) = (four* fine(i,j,k)
     &                        + three*(fine(i,j,k+1)+fine(i,j,k-1))
     &                        +   two*(fine(i,j,k+2)+fine(i,j,k-2))
     &                        +       (fine(i,j,k+3)+fine(i,j,k-3)))/sixteen
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratioy
            end do
           end do

         else if (ratiox .eq. 2 .and. ratioz .eq. 1) then

           do kc = lo(3), hi(3)
             k = ratioz*kc
             do ic = lo(1), hi(1)
               i = ratiox*ic
               crse(ic,jc,kc) =   half* fine(i  ,j,k)
     &                        + fourth*(fine(i+1,j,k)+fine(i-1,j,k))
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratioy
            end do
           end do

         else if (ratiox .eq. 4 .and. ratioz .eq. 1) then

           do kc = lo(3), hi(3)
             k = ratioz*kc
             do ic = lo(1), hi(1)
               i = ratiox*ic
               crse(ic,jc,kc) = (four* fine(i,j,k)
     &                        + three*(fine(i+1,j,k)+fine(i-1,j,k))
     &                        +   two*(fine(i+2,j,k)+fine(i-2,j,k))
     &                        +       (fine(i+3,j,k)+fine(i-3,j,k)))/sixteen
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratioy
            end do
           end do

         else if (ratiox .eq. 2 .and. ratioz .eq. 2) then
           do kc = lo(3), hi(3)
             k = ratioz*kc
             do ic = lo(1), hi(1)
               i = ratiox*ic
               crse(ic,jc,kc) = half*half*fine(i,j,k) 
     &              + half  *fourth*(fine(i+1,j,k  ) + fine(i-1,j,k  )
     &                              +fine(i  ,j,k+1) + fine(i  ,j,k-1))
     &              + fourth*fourth*(fine(i+1,j,k+1) + fine(i-1,j,k+1)
     &                              +fine(i+1,j,k-1) + fine(i-1,j,k-1)) 
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratioy
            end do
           end do

         else if (ratiox .eq. 2 .and. ratioz .eq. 4) then

           do kc = lo(3), hi(3)
             k = ratioz*kc
             do ic = lo(1), hi(1)
               i = ratiox*ic
               crse(ic,jc,kc) = ( four*four*fine(i,j,k)  +
     &              four*three*(fine(i  ,j,k+1) + fine(i  ,j,k-1)) +
     &              four*two  *(fine(i  ,j,k+2) + fine(i  ,j,k-2)) +
     &              four      *(fine(i  ,j,k+3) + fine(i  ,j,k-3)) +
     &              two *four *(fine(i+1,j,k  ) + fine(i-1,j,k  )) +
     &              two *three*(fine(i+1,j,k+1) + fine(i-1,j,k+1)+
     &                          fine(i+1,j,k-1) + fine(i-1,j,k-1)) +
     &              two *  two*(fine(i+1,j,k+2) + fine(i-1,j,k+2)+
     &                          fine(i+1,j,k-2) + fine(i-1,j,k-2)) +
     &              two       *(fine(i+1,j,k+3) + fine(i-1,j,k+3)+
     &                          fine(i+1,j,k-3) + fine(i-1,j,k-3)) ) / one28
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratioy 
             end do
           end do

         else if (ratiox .eq. 4 .and. ratioz .eq. 2) then

           do kc = lo(3), hi(3)
             k = ratioz*kc
             do ic = lo(1), hi(1)
               i = ratiox*ic
               crse(ic,jc,kc) = ( four*four*fine(i,j,k)  +
     &              four*three*(fine(i+1,j,k  ) + fine(i-1,j,k  )) +
     &              four*two  *(fine(i+2,j,k  ) + fine(i-2,j,k    )) +
     &              four      *(fine(i+3,j,k  ) + fine(i-3,j,k  )) +
     &              two *four *(fine(i  ,j,k+1) + fine(i  ,j,k-1)) +
     &              two *three*(fine(i+1,j,k+1) + fine(i+1,j,k-1)+
     &                          fine(i-1,j,k+1) + fine(i-1,j,k-1)) +
     &              two *  two*(fine(i+2,j,k+1) + fine(i+2,j,k-1)+
     &                          fine(i-2,j,k+1) + fine(i-2,j,k-1)) +
     &              two       *(fine(i+3,j,k+1) + fine(i+3,j,k-1)+
     &                          fine(i-3,j,k+1) + fine(i-3,j,k-1)) ) / one28
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratioy 
             end do
           end do

         else if (ratiox .eq. 4 .and. ratioz .eq. 4) then

           do kc = lo(3), hi(3)
             k = ratioz*kc
             do ic = lo(1), hi(1)
               i = ratiox*ic
               crse(ic,jc,kc) = four*four*fine(i,j,k)  +
     &              four *three*(fine(i+1,j,k  ) + fine(i-1,j,k  )+
     &                           fine(i  ,j,k+1) + fine(i  ,j,k-1)) +
     &              four *two  *(fine(i+2,j,k  ) + fine(i-2,j,k  )+
     &                           fine(i  ,j,k+2) + fine(i  ,j,k-2)) +
     &              four       *(fine(i+3,j,k  ) + fine(i-3,j,k  )+
     &                           fine(i  ,j,k+3) + fine(i  ,j,k-3)) +
     &              three*three*(fine(i+1,j,k+1) + fine(i-1,j,k+1)+
     &                           fine(i+1,j,k-1) + fine(i-1,j,k-1)) +
     &              three*two  *(fine(i+1,j,k+2) + fine(i-1,j,k+2)+
     &                           fine(i+1,j,k-2) + fine(i-1,j,k-2)+
     &                           fine(i+2,j,k+1) + fine(i-2,j,k+1)+
     &                           fine(i+2,j,k-1) + fine(i-2,j,k-1)) +
     &              three      *(fine(i+1,j,k+3) + fine(i-1,j,k+3)+
     &                           fine(i+1,j,k-3) + fine(i-1,j,k-3)+
     &                           fine(i+3,j,k+1) + fine(i-3,j,k+1)+
     &                           fine(i+3,j,k-1) + fine(i-3,j,k-1)) +
     &              two  *two  *(fine(i+2,j,k+2) + fine(i-2,j,k+2)+
     &                           fine(i+2,j,k-2) + fine(i-2,j,k-2)) +
     &              two        *(fine(i+2,j,k+3) + fine(i-2,j,k+3)+
     &                           fine(i+2,j,k-3) + fine(i-2,j,k-3)+
     &                           fine(i+3,j,k+2) + fine(i-3,j,k+2)+
     &                           fine(i+3,j,k-2) + fine(i-3,j,k-2)) +
     &                          (fine(i+3,j,k+3) + fine(i-3,j,k+3)+
     &                           fine(i+3,j,k-3) + fine(i-3,j,k-3)) 
               crse(ic,jc,kc) = crse(ic,jc,kc) / (sixteen*sixteen*ratioy)
             end do
           end do
         else
           print *,'BOGUS RATIOX AND/OR RATIOZ IN SRCRSEREG ',ratiox,ratioz
           stop
         end if
       else if (dir .eq. 2) then

c        ::::: sum in i and j directions
         kc = lo(3)
         k = ratioz*kc

         if (ratiox .eq. 1 .and. ratioy .eq. 1) then

           do jc = lo(2), hi(2)
             do ic = lo(1), hi(1)
               crse(ic,jc,kc) = fine(ic,jc,k) / ratioz
            end do
           end do

         else if (ratiox .eq. 1 .and. ratioy .eq. 2) then

           do jc = lo(2), hi(2)
             j = ratioy*jc
             do ic = lo(1), hi(1)
               i = ratiox*ic
               crse(ic,jc,kc) =   half* fine(i,j,k)
     &                        + fourth*(fine(i,j+1,k)+fine(i,j-1,k))
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratioz
            end do
           end do

         else if (ratiox .eq. 1 .and. ratioy .eq. 4) then

           do jc = lo(2), hi(2)
             j = ratioy*jc
             do ic = lo(1), hi(1)
               i = ratiox*ic
               crse(ic,jc,kc) = (four* fine(i,j,k)
     &                        + three*(fine(i,j+1,k)+fine(i,j-1,k))
     &                        +   two*(fine(i,j+2,k)+fine(i,j-2,k))
     &                        +       (fine(i,j+3,k)+fine(i,j-3,k)))/sixteen
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratioz
            end do
           end do

         else if (ratiox .eq. 2 .and. ratioy .eq. 1) then

           do jc = lo(2), hi(2)
             j = ratioy*jc
             do ic = lo(1), hi(1)
               i = ratiox*ic
               crse(ic,jc,kc) =   half* fine(i  ,j,k)
     &                        + fourth*(fine(i+1,j,k)+fine(i-1,j,k))
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratioz
            end do
           end do

         else if (ratiox .eq. 4 .and. ratioy .eq. 1) then

           do jc = lo(2), hi(2)
             j = ratioy*jc
             do ic = lo(1), hi(1)
               i = ratiox*ic
               crse(ic,jc,kc) = (four* fine(i,j,k)
     &                        + three*(fine(i+1,j,k)+fine(i-1,j,k))
     &                        +   two*(fine(i+2,j,k)+fine(i-2,j,k))
     &                        +       (fine(i+3,j,k)+fine(i-3,j,k)))/sixteen
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratioz
            end do
           end do

         else if (ratiox .eq. 2 .and. ratioy .eq. 2) then
           do jc = lo(2), hi(2)
             j = ratioy*jc
             do ic = lo(1), hi(1)
               i = ratiox*ic
               crse(ic,jc,kc) = half*half*fine(i,j,k) 
     &              + half  *fourth*(fine(i+1,j  ,k) + fine(i-1,j  ,k)
     &                              +fine(i  ,j+1,k) + fine(i  ,j-1,k))
     &              + fourth*fourth*(fine(i+1,j+1,k) + fine(i-1,j+1,k)
     &                              +fine(i+1,j-1,k) + fine(i-1,j-1,k)) 
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratioz
            end do
           end do

         else if (ratiox .eq. 2 .and. ratioy .eq. 4) then

           do jc = lo(2), hi(2)
             j = ratioy*jc
             do ic = lo(1), hi(1)
               i = ratiox*ic
               crse(ic,jc,kc) = ( four*four*fine(i,j,k)  +
     &              four*three*(fine(i  ,j+1,k) + fine(i  ,j-1,k)) +
     &              four*two  *(fine(i  ,j+2,k) + fine(i  ,j-2,k)) +
     &              four      *(fine(i  ,j+3,k) + fine(i  ,j-3,k)) +
     &              two *four *(fine(i+1,j  ,k) + fine(i-1,j  ,k)) +
     &              two *three*(fine(i+1,j+1,k) + fine(i-1,j+1,k)+
     &                          fine(i+1,j-1,k) + fine(i-1,j-1,k)) +
     &              two *  two*(fine(i+1,j+2,k) + fine(i-1,j+2,k)+
     &                          fine(i+1,j-2,k) + fine(i-1,j-2,k)) +
     &              two       *(fine(i+1,j+3,k) + fine(i-1,j+3,k)+
     &                          fine(i+1,j-3,k) + fine(i-1,j-3,k)) ) / one28
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratioz 
             end do
           end do

         else if (ratiox .eq. 4 .and. ratioy .eq. 2) then

           do jc = lo(2), hi(2)
             j = ratioy*jc
             do ic = lo(1), hi(1)
               i = ratiox*ic
               crse(ic,jc,kc) = ( four*four*fine(i,j,k)  +
     &              four*three*(fine(i+1,j  ,k) + fine(i-1,j  ,k)) +
     &              four*two  *(fine(i+2,j  ,k) + fine(i-2,j  ,k)) +
     &              four      *(fine(i+3,j  ,k) + fine(i-3,j  ,k)) +
     &              two *four *(fine(i  ,j+1,k) + fine(i  ,j-1,k)) +
     &              two *three*(fine(i+1,j+1,k) + fine(i+1,j-1,k)+
     &                          fine(i-1,j+1,k) + fine(i-1,j-1,k)) +
     &              two *  two*(fine(i+2,j+1,k) + fine(i+2,j-1,k)+
     &                          fine(i-2,j+1,k) + fine(i-2,j-1,k)) +
     &              two       *(fine(i+3,j+1,k) + fine(i+3,j-1,k)+
     &                          fine(i-3,j+1,k) + fine(i-3,j-1,k)) ) / one28
               crse(ic,jc,kc) = crse(ic,jc,kc) / ratioz 
             end do
           end do

         else if (ratiox .eq. 4 .and. ratioy .eq. 4) then

           do jc = lo(2), hi(2)
             j = ratioy*jc
             do ic = lo(1), hi(1)
               i = ratiox*ic
               crse(ic,jc,kc) = four*four*fine(i,j,k)  +
     &              four *three*(fine(i+1,j  ,k) + fine(i-1,j  ,k)+
     &                           fine(i  ,j+1,k) + fine(i  ,j-1,k)) +
     &              four *two  *(fine(i+2,j  ,k) + fine(i-2,j  ,k)+
     &                           fine(i  ,j+2,k) + fine(i  ,j-2,k)) +
     &              four       *(fine(i+3,j  ,k) + fine(i-3,j  ,k)+
     &                           fine(i  ,j+3,k) + fine(i  ,j-3,k)) +
     &              three*three*(fine(i+1,j+1,k) + fine(i-1,j+1,k)+
     &                           fine(i+1,j-1,k) + fine(i-1,j-1,k)) +
     &              three*two  *(fine(i+1,j+2,k) + fine(i-1,j+2,k)+
     &                           fine(i+1,j-2,k) + fine(i-1,j-2,k)+
     &                           fine(i+2,j+1,k) + fine(i-2,j+1,k)+
     &                           fine(i+2,j-1,k) + fine(i-2,j-1,k)) +
     &              three      *(fine(i+1,j+3,k) + fine(i-1,j+3,k)+
     &                           fine(i+1,j-3,k) + fine(i-1,j-3,k)+
     &                           fine(i+3,j+1,k) + fine(i-3,j+1,k)+
     &                           fine(i+3,j-1,k) + fine(i-3,j-1,k)) +
     &              two  *two  *(fine(i+2,j+2,k) + fine(i-2,j+2,k)+
     &                           fine(i+2,j-2,k) + fine(i-2,j-2,k)) +
     &              two        *(fine(i+2,j+3,k) + fine(i-2,j+3,k)+
     &                           fine(i+2,j-3,k) + fine(i-2,j-3,k)+
     &                           fine(i+3,j+2,k) + fine(i-3,j+2,k)+
     &                           fine(i+3,j-2,k) + fine(i-3,j-2,k)) +
     &                          (fine(i+3,j+3,k) + fine(i-3,j+3,k)+
     &                           fine(i+3,j-3,k) + fine(i-3,j-3,k)) 
               crse(ic,jc,kc) = crse(ic,jc,kc) / (sixteen*sixteen*ratioz)
             end do
           end do
         else
           print *,'BOGUS RATIOX AND/OR RATIOY IN SRCRSEREG ',ratiox,ratioy
           stop
         end if
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

      integer    i, j, k

      do k = ARG_L3(mask), ARG_H3(mask)
       do j = ARG_L2(mask), ARG_H2(mask)
        do i = ARG_L1(mask), ARG_H1(mask)
 
          mask(i,j,k) = cells(i,j  ,k  ) + cells(i-1,j  ,k  )
     $                + cells(i,j-1,k  ) + cells(i-1,j-1,k  )
     $                + cells(i,j  ,k-1) + cells(i-1,j  ,k-1)
     $                + cells(i,j-1,k-1) + cells(i-1,j-1,k-1)

        end do
       end do
      end do

      end

c ::: -----------------------------------------------------------
c ::: modify mask from values 0. through 8. to values 0. or 1.
c :::
c ::: INPUTS/OUTPUTS:
c ::: mask       <=  node-centered mask array
c ::: DIMS(mask)  => index limits for mask
c ::: -----------------------------------------------------------

      subroutine FORT_CONVERTMASK(mask,DIMS(mask))

      integer    DIMDEC(mask)
      REAL_T     mask(DIMV(mask))
 
      integer    i, j, k

      do k = ARG_L3(mask), ARG_H3(mask)
       do j = ARG_L2(mask), ARG_H2(mask)
        do i = ARG_L1(mask), ARG_H1(mask)
 
          if (mask(i,j,k) .gt. 7.5) then
             mask(i,j,k) = zero
          else
             mask(i,j,k) = one
          end if

        end do
       end do
      end do
 
      end

