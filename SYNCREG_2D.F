c
c $Id: SYNCREG_2D.F,v 1.9 1999-06-03 18:17:45 almgren Exp $
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
c ::: This routine will compute a node centered divergence
c ::: from a cell centered field
c :::
c ::: NOTE: node (i,j) is in lower left corner of cell (i,j)      
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

      subroutine FORT_SRDIVU(u,DIMS(u),divu,DIMS(divu),
     $                       lo,hi,dx,mult,is_rz)

      integer    DIMDEC(u)
      integer    DIMDEC(divu)
      integer    lo(SDIM), hi(SDIM)
      REAL_T     dx(SDIM), mult
      REAL_T     u(DIMV(u),SDIM)
      REAL_T     divu(DIMV(divu))
      integer    is_rz
      
      integer    i, j, lenx, leny
      REAL_T     ux, vy
      REAL_T     rp, rm, vext
      REAL_T     rpinv, rminv

c     NOTE: the factor of mult here accounts for the fact that one
c           takes ratio time steps on the fine grid for every one
c           timestep on the coarse grid.  mult = 1/ratio. 
      
      lenx = hi(1)-lo(1)+1
      leny = hi(2)-lo(2)+1

      if (lenx .ge. leny) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ux = half*(u(i,j,1)-u(i-1,j,1)+u(i,j-1,1)-u(i-1,j-1,1))/dx(1)
               vy = half*(u(i,j,2)-u(i,j-1,2)+u(i-1,j,2)-u(i-1,j-1,2))/dx(2)
               divu(i,j) = mult*(ux+vy)
            end do
         end do
      else 
         do i = lo(1), hi(1)
            do j = lo(2), hi(2)
               ux = half*(u(i,j,1)-u(i-1,j,1)+u(i,j-1,1)-u(i-1,j-1,1))/dx(1)
               vy = half*(u(i,j,2)-u(i,j-1,2)+u(i-1,j,2)-u(i-1,j-1,2))/dx(2)
               divu(i,j) = mult*(ux+vy)
            end do
         end do
      end if

#if BL_PRVERSION==9
      if (is_rz .eq. 1) then

       if (lenx .ge. leny) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rm = twelve*(dfloat(i)-half)*dx(1)
               rp = twelve*(dfloat(i)+half)*dx(1)
               rminv = dx(1) / (rm * dx(2))
               rpinv = dx(1) / (rp * dx(2))
               vext = (u(i,j-1,2)-u(i  ,j  ,2))*rpinv + 
     $                (u(i-1,j,2)-u(i-1,j-1,2))*rminv
               divu(i,j) = divu(i,j) + mult*vext
            end do
         end do
       else 
         do i = lo(1), hi(1)
            do j = lo(2), hi(2)
               rm = twelve*(dfloat(i)-half)*dx(1)
               rp = twelve*(dfloat(i)+half)*dx(1)
               rminv = dx(1) / (rm * dx(2))
               rpinv = dx(1) / (rp * dx(2))
               vext = (u(i,j-1,2)-u(i  ,j  ,2))*rpinv + 
     $                (u(i-1,j,2)-u(i-1,j-1,2))*rminv
               divu(i,j) = divu(i,j) + mult*vext
            end do
         end do
       end if

      end if
#endif
      
      end

c ::: -----------------------------------------------------------
c ::: This routine will compute one version of DIV(sigma*GRAD(phi))
c ::: from phi and sigma.
c ::: This version uses GRAD-PERP.      
c :::
c ::: NOTE: phi is node-centered
c :::       sigma is cell-centered      
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: phi        => node centered data
c ::: DIMS(phi)  => index limits for phi
c ::: sigma      => cell-centered coeficient array
c ::: DIMS(sigma)=> index limits for sigma      
c ::: divgp     <=  node-centered divergence of gradient
c ::: DIMS(divgp)=> index limits for divgp
c ::: lo,hi      => node-centered subregion of divgp to define
c ::: dx         => cell sizes      
c ::: mult       => scalar multiplicative factor      
c ::: -----------------------------------------------------------

      subroutine FORT_SRDGPPHI(phi,DIMS(phi),sigma,DIMS(sigma),
     &                         divgp,DIMS(divgp),lo,hi,dx,mult,is_rz)

      integer    DIMDEC(phi)
      integer    DIMDEC(sigma)
      integer    DIMDEC(divgp)
      integer    lo(SDIM), hi(SDIM)
      REAL_T     dx(SDIM), mult
      REAL_T     phi(DIMV(phi))
      REAL_T     sigma(DIMV(sigma))
      REAL_T     divgp(DIMV(divgp))
      integer    is_rz
      
      integer    i, j, lenx, leny
      REAL_T     p1, p2, p3, p4, denom
      REAL_T     rm,rp,rminv,rpinv
      REAL_T     gpy_ul,gpy_ur,gpy_ll,gpy_lr
      REAL_T     pext,vext

#if BL_PRVERSION==9
c     FOR NINE-POINT STENCIL
      denom = one/(twelve*dx(1)**2) + one/(twelve*dx(2)**2)
#elif BL_PRVERSION==5
c     FOR FIVE-POINT STENCIL
      denom = one/(four*dx(1)**2) + one/(four*dx(2)**2)
#else
#error "Unknown PRVERSION setting"
#endif

      lenx = hi(1)-lo(1)+1
      leny = hi(2)-lo(2)+1
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            p1 = sigma(i-1,j-1)*(-phi(i-1,j-1) + phi(i  ,j-1)
     &                           +phi(i-1,j  ) - phi(i  ,j  ))
            p2 = sigma(i-1,j  )*(-phi(i-1,j+1) + phi(i  ,j+1)
     &                           +phi(i-1,j)   - phi(i  ,j  ))
            p3 = sigma(i  ,j-1)*(-phi(i+1,j-1) + phi(i  ,j-1)
     &                           +phi(i+1,j)   - phi(i  ,j  ))
            p4 = sigma(i  ,j  )*(-phi(i+1,j+1) + phi(i  ,j+1)
     &                           +phi(i+1,j)   - phi(i  ,j  ))
            divgp(i,j) = denom*(p1+p2+p3+p4)
         end do
      end do

#if BL_PRVERSION==9
      if (is_rz .eq. 1) then

c      extra term added to DGphi for cylindrical coordinates
       do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            rp = dx(1) / (two * (dfloat(i)+half)*dx(1))
            rm = dx(1) / (two * (dfloat(i)-half)*dx(1))

            pext = ((sigma(i-1,j-1) * (phi(i,j-1) - phi(i,j)) +
     $               sigma(i-1,j  ) * (phi(i,j+1) - phi(i,j))) * rm -
     $              (sigma(i  ,j-1) * (phi(i,j-1) - phi(i,j)) +
     $               sigma(i  ,j  ) * (phi(i,j+1) - phi(i,j))) * rp)

            divgp(i,j) = divgp(i,j) + denom*pext
         end do
       end do

c      extra term added to DGbarphi for cylindrical coordinates --
c        here we subtract it since we want D G^perp phi
       do j = lo(2), hi(2)
       do i = lo(1), hi(1)
           rm = twelve*(dfloat(i)-half)*dx(1)
           rp = twelve*(dfloat(i)+half)*dx(1)
           rminv = dx(1) / (rm * dx(2)**2)
           rpinv = dx(1) / (rp * dx(2)**2)
   
	   gpy_ur = half * sigma(i  ,j  ) * ( phi(i+1,j+1) + phi(i,j+1) 
     $                                       -phi(i+1,j  ) - phi(i,j  ))
	   gpy_ul = half * sigma(i-1,j  ) * ( phi(i-1,j+1) + phi(i,j+1) 
     $                                       -phi(i-1,j  ) - phi(i,j  ))
	   gpy_lr = half * sigma(i  ,j-1) * ( phi(i+1,j  ) + phi(i,j  ) 
     $                                       -phi(i+1,j-1) - phi(i,j-1))
	   gpy_ll = half * sigma(i-1,j-1) * ( phi(i-1,j  ) + phi(i,j  ) 
     $                                       -phi(i-1,j-1) - phi(i,j-1))
           vext =  (gpy_lr - gpy_ur)*rpinv 
     $            +(gpy_ul - gpy_ll)*rminv
           divgp(i,j) = divgp(i,j) - vext
        end do
        end do
      endif
#endif

      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
           divgp(i,j) = divgp(i,j) * mult
      end do
      end do

      end

c ::: -----------------------------------------------------------
c ::: This routine will compute one version of DIV(sigma*GRAD(phi))
c ::: from phi and sigma.
c ::: This is the non-PERP version.
c :::
c ::: INPUTS/OUTPUTS:
c ::: phi        => node-centered data
c ::: DIMS(phi)  => index limits for phi
c ::: sigma      => cell-centered coeficient array
c ::: DIMS(sigma)=> index limits for sigma      
c ::: divgp     <=  node-centered divergence of gradient
c ::: DIMS(dgphi)=> index limits for divgp
c ::: lo,hi      => node-centered subregion of divgp to define
c ::: dx         => cell sizes      
c ::: mult       => scalar multiplicative factor      
c ::: -----------------------------------------------------------

      subroutine FORT_SRDGPHI(phi,DIMS(phi),sigma,DIMS(sigma),
     $                        divgp,DIMS(divgp),lo,hi,dx,mult,is_rz)

      integer    DIMDEC(phi)
      integer    DIMDEC(sigma)
      integer    DIMDEC(divgp)
      integer    lo(SDIM), hi(SDIM)
      REAL_T     dx(SDIM), mult
      REAL_T     phi(DIMV(phi))
      REAL_T     sigma(DIMV(sigma))
      REAL_T     divgp(DIMV(divgp))
      integer    is_rz
      
      integer    i, j
      REAL_T     p1, p2, p3, p4
      REAL_T     denomx, denomy
      REAL_T     rp, rm, pext

c     NOTE: the factor of mult here accounts for the fact that one
c           takes ratio time steps on the fine grid for every one
c           timestep on the coarse grid.  mult = 1/ratio. 
      
#if BL_PRVERSION==9
c
c     FOR NINE-POINT STENCIL
      denomx = -mult/(six*dx(1)**2)
      denomy = -mult/(six*dx(2)**2)

      if (dx(1) .eq. dx(2)) then

         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               p1 = sigma(i-1,j-1)*( two*phi(i-1,j-1) +      phi(i  ,j-1)
     &              +    phi(i-1,j  ) - four*phi(i  ,j  ))
               p2 = sigma(i-1,j  )*(two*phi(i-1,j+1) +      phi(i  ,j+1)
     &              +    phi(i-1,j)   - four*phi(i  ,j  ))
               p3 = sigma(i  ,j-1)*(two*phi(i+1,j-1) +      phi(i  ,j-1)
     &              +    phi(i+1,j)   - four*phi(i  ,j  ))
               p4 = sigma(i  ,j  )*(two*phi(i+1,j+1) +      phi(i  ,j+1)
     &              +    phi(i+1,j)   - four*phi(i  ,j  ))
               divgp(i,j) = denomx*(p1+p2+p3+p4)

            end do
         end do
         
      else

         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               p1 = sigma(i-1,j-1)* (
     $              denomx * (     phi(i-1,j-1)     -     phi(i  ,j-1)
     $              +two*phi(i-1,j  )     - two*phi(i  ,j  ) ) +
     $              denomy * (     phi(i-1,j-1)     -     phi(i-1,j  )
     $              +two*phi(i  ,j-1)     - two*phi(i  ,j  ) ) )

               p2 = sigma(i-1,j  )* (
     $              denomx * (     phi(i-1,j+1)     -     phi(i  ,j+1)
     $              +two*phi(i-1,j  )     - two*phi(i  ,j  ) ) +
     $              denomy * (     phi(i-1,j+1)     -     phi(i-1,j  )
     $              +two*phi(i  ,j+1)     - two*phi(i  ,j  ) ) )

               p3 = sigma(i  ,j-1)* (
     $              denomx * (     phi(i+1,j-1)     -     phi(i  ,j-1)
     $              +two*phi(i+1,j  )     - two*phi(i  ,j  ) ) +
     $              denomy * (     phi(i+1,j-1)     -     phi(i+1,j  )
     $              +two*phi(i  ,j-1)     - two*phi(i  ,j  ) ) )

               p4 = sigma(i  ,j  )* (
     $              denomx * (     phi(i+1,j+1)     -     phi(i  ,j+1)
     $              +two*phi(i+1,j  )     - two*phi(i  ,j  ) ) +
     $              denomy * (     phi(i+1,j+1)     -     phi(i+1,j  )
     $              +two*phi(i  ,j+1)     - two*phi(i  ,j  ) ) )

               divgp(i,j) = p1+p2+p3+p4

            end do
         end do
      end if

c     extra term added for cylindrical coordinates
      if (is_rz .eq. 1) then

         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rp = dx(1) / (two * (dfloat(i)+half)*dx(1))
               rm = dx(1) / (two * (dfloat(i)-half)*dx(1))

               pext = ((sigma(i-1,j-1) * (phi(i,j-1) - phi(i,j)) +
     $              sigma(i-1,j  ) * (phi(i,j+1) - phi(i,j))) * rm -
     $              (sigma(i  ,j-1) * (phi(i,j-1) - phi(i,j)) +
     $              sigma(i  ,j  ) * (phi(i,j+1) - phi(i,j))) * rp)

               divgp(i,j) = divgp(i,j) + denomy*pext
            end do
         end do

      end if
#elif BL_PRVERSION==5
c
c     FOR FIVE-POINT STENCIL
c
      denomx = -mult/(two*dx(1)**2)
      denomy = -mult/(two*dx(2)**2)

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            p1 = (sigma(i-1,j-1) + sigma(i-1,j)) * 
     $           (  phi(i-1,j  ) -   phi(i  ,j))
            p2 = (sigma(i  ,j-1) + sigma(i  ,j)) * 
     $           (  phi(i+1,j  ) -   phi(i  ,j))
            p3 = (sigma(i-1,j-1) + sigma(i,j-1)) * 
     $           (  phi(i  ,j-1) -   phi(i,j  ))
            p4 = (sigma(i-1,j  ) + sigma(i,j  )) * 
     $           (  phi(i  ,j+1) -   phi(i,j  ))
            divgp(i,j) = denomx*(p1+p2)+ denomy*(p3+p4)
         end do
      end do
#else
#error "Unknown PRVERSION setting"
#endif
      
      end

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
     &                          ratios)

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

c     NOTE: the reason that the coefficients add up to 1/ratio and NOT to 1
c           is because the divergences and DGphi were computed using the
c           local dx, so in order to add the fine contribution to the
c           coarse contribution the fine contribution must be weighted
c           by dx_fine/dx_coarse = 1/ratio
      
      ratiox = ratios(0)
      ratioy = ratios(1)
      
      if (dir .eq. 0) then

c        ::::: sum in j direction
         ic = lo(1)
         i = ratiox*ic

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

      else

c        ::::: sum in i direction
         jc = lo(2)
         j = ratioy*jc

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
      end if
      
      end

#if 0
// legacy code that could prove useful for debugging
c :: ----------------------------------------------------------
c :: FIX_SYNC_FINE_CORNER_DSDT
c ::  fixes domain corner nodes for use in sync projection
c ::
c :: INPUTS / OUTPUTS:
c ::  nghost      => number of ghost nodes
c ::  DIMS(dat)   => index limits of dat
c ::  dat        <=> sync register values
c ::  domlo,domhi => index limits of domain
c :: ----------------------------------------------------------
c ::
      subroutine FORT_FIX_SYNC_FINE_DSDT(nghost, DIMS(dat), dat, 
     &                                   domlo, domhi) 

      implicit none
      integer nghost
      integer DIMDEC(dat)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  dat(DIMV(dat))
      integer i,j

c ... fix up sides
c     Note: double counting of corners is intended

      if(ARG_L1(dat)+nghost .eq. domlo(1)) then
        i = ARG_L1(dat)+nghost
        do j = ARG_L2(dat)+nghost, ARG_H2(dat)-nghost
          dat(i,j) = .5*dat(i,j)
        end do
      end if

      if(ARG_H1(dat)-nghost .eq. domhi(1)+1) then
        i = ARG_H1(dat)-nghost
        do j = ARG_L2(dat)+nghost, ARG_H2(dat)-nghost
          dat(i,j) = .5*dat(i,j)
        end do
      end if

      if(ARG_L2(dat)+nghost .eq. domlo(2)) then
        j = ARG_L2(dat)+nghost
        do i = ARG_L1(dat)+nghost, ARG_H1(dat)-nghost
          dat(i,j) = .5*dat(i,j)
        end do
      end if

      if(ARG_H2(dat)-nghost .eq. domhi(2)+1) then
        j = ARG_H2(dat)-nghost
        do i = ARG_L1(dat)+nghost, ARG_H1(dat)-nghost
          dat(i,j) = .5*dat(i,j)
        end do
      end if

      end
#endif

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
