c
c $Id: PROJECTION_2D.F,v 1.34 1999-07-27 00:09:36 propp Exp $
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

#define SDIM 2

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define SMALL 1.0e-10
#else
#define SMALL 1.0d-10
#endif

       subroutine FORT_FILTERP (
     &     pnew,scratch,DIMS(p),
     &     lo,hi,dx,fac,
     &     wrap_around_x,wrap_around_y)
c
c     ----------------------------------------------------------
c     filter the pressure coming out of the level projection
c     
c     pnew := pnew - fac * lap_5pt(pnew)
c     
c     INPUTS / OUTPUTS:
c     pnew       => 
c     pnew      <=> 
c     p_lo,p_hi  => index limits of pressure array
c     lo,hi      => index limits of cell centered interior
c     fac        => fudge factor for filtering
c     wrap_around_x => is this a single periodic box in the x-direction
c     wrap_around_y => is this a single periodic box in the y-direction
c     
c     ----------------------------------------------------------
c 
       integer    DIMDEC(p)
       REAL_T     pnew(DIMV(p))
       REAL_T     scratch(DIMV(p))
       integer    lo(2), hi(2)
       REAL_T     dx(2)
       REAL_T     fac
       integer    wrap_around_x
       integer    wrap_around_y
       integer i, j

       if ( abs(dx(1)-dx(2)) .gt. 1.e-8 .and. fac .gt. 1.e-4 )then
          print *, 'is filter properly formatted for dx != dy ?'
       end if

       do j = lo(2)+1, hi(2)
          do i = lo(1)+1, hi(1)
             scratch(i,j) = four * pnew(i,j) - 
     $         (pnew(i-1,j) + pnew(i+1,j) + pnew(i,j-1) + pnew(i,j+1) )
             scratch(i,j) = - eighth * scratch(i,j) * fac
          end do
       end do
 
       if (wrap_around_x .eq. 1) then
         i = lo(1)
         do j = lo(2)+1, hi(2)
             scratch(i,j) = four * pnew(i,j) - 
     $         (pnew(hi(1),j) + pnew(i+1,j) + pnew(i,j-1) + pnew(i,j+1) )
             scratch(i,j) = - eighth * scratch(i,j) * fac
         end do

         i = hi(1)+1
         do j = lo(2)+1, hi(2)
             scratch(i,j) = four * pnew(i,j) - 
     $         (pnew(i-1,j) + pnew(lo(1)+1,j) + pnew(i,j-1) + pnew(i,j+1) )
             scratch(i,j) = - eighth * scratch(i,j) * fac
         end do
       end if
 
       if (wrap_around_y .eq. 1) then
         j = lo(2)
         do i = lo(1)+1, hi(1)
             scratch(i,j) = four * pnew(i,j) - 
     $         (pnew(i-1,j) + pnew(i+1,j) + pnew(i,hi(2)) + pnew(i,j+1) )
             scratch(i,j) = - eighth * scratch(i,j) * fac
         end do

         j = hi(2)+1
         do i = lo(1)+1, hi(1)
             scratch(i,j) = four * pnew(i,j) - 
     $         (pnew(i-1,j) + pnew(i+1,j) + pnew(i,j-1) + pnew(i,lo(2)+1) )
             scratch(i,j) = - eighth * scratch(i,j) * fac
         end do
       end if
 
       if (wrap_around_x .eq. 1 .and. wrap_around_y .eq. 1) then
         i = lo(1)
         j = lo(2)
           scratch(i,j) = four * pnew(i,j) - 
     $       (pnew(hi(1),j) + pnew(i+1,j) + pnew(i,hi(2)) + pnew(i,j+1) )
           scratch(i,j) = - eighth * scratch(i,j) * fac

         i = hi(1)+1
         j = lo(2)
           scratch(i,j) = four * pnew(i,j) - 
     $       (pnew(i-1,j) + pnew(lo(1)+1,j) + pnew(i,hi(2)) + pnew(i,j+1) )
           scratch(i,j) = - eighth * scratch(i,j) * fac

         i = lo(1)
         j = hi(2)+1
           scratch(i,j) = four * pnew(i,j) - 
     $       (pnew(hi(1),j) + pnew(i+1,j) + pnew(i,j-1) + pnew(i,lo(2)+1) )
           scratch(i,j) = - eighth * scratch(i,j) * fac

         i = hi(1)+1
         j = hi(2)+1
           scratch(i,j) = four * pnew(i,j) - 
     $       (pnew(i-1,j) + pnew(lo(1)+1,j) + pnew(i,j-1) + pnew(i,lo(2)+1) )
           scratch(i,j) = - eighth * scratch(i,j) * fac
       end if

       do j = lo(2)+1, hi(2)
          do i = lo(1)+1, hi(1)
             pnew(i,j) = pnew(i,j) + scratch(i,j)
          end do
       end do

       if (wrap_around_x .eq. 1) then
         do j = lo(2)+1, hi(2)
           pnew(lo(1)  ,j) = pnew(lo(1)  ,j) + scratch(lo(1)  ,j)
           pnew(hi(1)+1,j) = pnew(hi(1)+1,j) + scratch(hi(1)+1,j)
         end do
       end if

       if (wrap_around_y .eq. 1) then
         do i = lo(1)+1, hi(1)
           pnew(i,lo(2)  ) = pnew(i,lo(2)  ) + scratch(i,lo(2)  )
           pnew(i,hi(2)+1) = pnew(i,hi(2)+1) + scratch(i,hi(2)+1)
         end do
       end if

       if (wrap_around_x .eq. 1 .and. wrap_around_y .eq. 1) then
         pnew(lo(1)  ,lo(2)  ) = pnew(lo(1)  ,lo(2)  ) + scratch(lo(1)  ,lo(2)  )
         pnew(hi(1)+1,lo(2)  ) = pnew(hi(1)+1,lo(2)  ) + scratch(hi(1)+1,lo(2)  )
         pnew(lo(1)  ,hi(2)+1) = pnew(lo(1)  ,hi(2)+1) + scratch(lo(1)  ,hi(2)+1)
         pnew(hi(1)+1,hi(2)+1) = pnew(hi(1)+1,hi(2)+1) + scratch(hi(1)+1,hi(2)+1)
       end if

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

       integer i, j, n

       do n = 1, SDIM
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1
                unew(i,j,n) = uold(i,j,n) + dt*unew(i,j,n)
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
     
      integer i, j, n
      REAL_T dt_inv
   
      dt_inv = one/dt

      do n = 1, SDIM
         do j = lo(2)-1, hi(2)+1
            do i = lo(1)-1, hi(1)+1
               unew(i,j,n) = (unew(i,j,n)-uold(i,j,n))*dt_inv
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

      integer i, j, n

      do n = 1, nvar
         do j = boxlo(2), boxhi(2)
            do i = boxlo(1), boxhi(1)
               un(i,j,n) = un(i,j,n) + alpha*uo(i,j,n)
            end do
         end do
      end do

      end

      subroutine FORT_RADMPY(a,DIMS(grid),domlo,domhi,ng,r,nr,bogus_value)
c 
c     multiply A by Radius r
c
c 
c     NOTE: THIS ROUTINE HAS BEEN MODIFIED SO THAT ALL VALUES
c           OUTSIDE THE DOMAIN ARE SET TO BOGUS VALUE
c
      integer    ng,nr
      integer    DIMDEC(grid)
      integer    domlo(2), domhi(2)
      REAL_T     a(ARG_L1(grid)-ng:ARG_H1(grid)+ng, 
     $             ARG_L2(grid)-ng:ARG_H2(grid)+ng)
      REAL_T     r(ARG_L1(grid)-nr:ARG_H1(grid)+nr)
      REAL_T     bogus_value

      integer i, j

      do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do i = ARG_L1(grid)-ng, ARG_H1(grid)+ng
           a(i,j) = r(i)*a(i,j)
         end do
      end do

      if (ARG_L1(grid)-ng .lt. domlo(1)) then
         do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do i = ARG_L1(grid)-ng, domlo(1)-1
           a(i,j) = bogus_value
         end do
         end do
      end if

      if (ARG_H1(grid)+ng .gt. domhi(1)) then
         do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do i = domhi(1)+1, ARG_H1(grid)+ng
           a(i,j) = bogus_value
         end do
         end do
      end if

      end

      subroutine FORT_RADDIV(a,DIMS(grid),domlo,domhi,ng,r,nr,bogus_value)
c 
c     divide A by Radius r
c
c 
c     NOTE: THIS ROUTINE HAS BEEN MODIFIED SO THAT ALL VALUES
c           OUTSIDE THE DOMAIN ARE SET TO BOGUS VALUE
c
      integer    ng,nr
      integer    DIMDEC(grid)
      integer    domlo(2), domhi(2)
      REAL_T     a(ARG_L1(grid)-ng:ARG_H1(grid)+ng, 
     $             ARG_L2(grid)-ng:ARG_H2(grid)+ng)
      REAL_T     r(ARG_L1(grid)-nr:ARG_H1(grid)+nr)
      REAL_T     bogus_value

      integer i, j

      do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do i = ARG_L1(grid)-ng, ARG_H1(grid)+ng
           a(i,j) = a(i,j)/r(i)
         end do
      end do

      if (ARG_L1(grid)-ng .lt. domlo(1)) then
         do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do i = ARG_L1(grid)-ng, domlo(1)-1
           a(i,j) = bogus_value
         end do
         end do
      end if

      if (ARG_H1(grid)+ng .gt. domhi(1)) then
         do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do i = domhi(1)+1, ARG_H1(grid)+ng
           a(i,j) = bogus_value
         end do
         end do
      end if

      end

      subroutine FORT_HGN2C(
     &     isrz,lrweighted, DIMS(nodedat), nodedat,
     &     DIMS(ccdat), lo, hi, ccdat)

c     ----------------------------------------------------------
c     HGN2C
c     averages node centered data to cell centers for use in 
c     holy grail projection

      implicit none
      integer isrz,lrweighted
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(ccdat)
      integer DIMDEC(nodedat)
      REAL_T  nodedat(DIMV(nodedat))
      REAL_T  ccdat(DIMV(ccdat))

      integer i,j

      if (ARG_H1(ccdat)   .lt. lo(1) .or. 
     &    ARG_L1(ccdat)   .gt. hi(1) .or. 
     &    ARG_H2(ccdat)   .lt. lo(2) .or. 
     &    ARG_L2(ccdat)   .gt. hi(2) .or. 
     &    ARG_H1(nodedat) .lt. lo(1)+1 .or. 
     &    ARG_L1(nodedat) .gt. hi(1) .or. 
     &    ARG_H2(nodedat) .lt. lo(2)+1 .or. 
     &    ARG_L2(nodedat) .gt. hi(2) ) then 
        write(6,*) "FORT_HG_CELL_TO_NODE: bad index limits"
        stop
      end if

      if(isrz.eq.1.and.lrweighted.ne.1)then
        write(6,*)
     &    'FORT_HGN2C : isrz=1 and lrweighted!=1 not yet implemented'
#ifndef	BL_NO_FORT_FLUSH
c        call flush(6)
#endif
        stop
      end if

      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          ccdat(i,j) = fourth*(nodedat(i,j)+nodedat(i+1,j)+
     &                         nodedat(i,j+1)+nodedat(i+1,j+1))
        end do
      end do

      end

      subroutine FORT_HGC2N(
     &     nghost, DIMS(dat), dat, rcen,
     &     DIMS(rhs), rhs,
     &     domlo, domhi, dr, is_rz, imax) 
c
c     ----------------------------------------------------------
c     HGC2N
c     averages cell centered data to nodes for use in 
c     holy grail projection
c     
c     INPUTS / OUTPUTS:
c     nghost      => indicates buffer of rhs that does not need values
c     dat         => cell centered array to be averaged
c     DIMS(dat)   => index limits of dat
c     rcen        => r-coordinate cell centers if geoem is r-z; 
c     otherwise, should be 1
c     rhslo,rhshi => index extents of rhs
c     rhs         <= node centered array with results
c     ----------------------------------------------------------
c 
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
      REAL_T  factor, rhi, rlo

      if (ARG_L1(rhs)+1 .lt. ARG_L1(dat) .or. 
     $     ARG_H1(rhs)-1 .gt. ARG_H1(dat) .or.
     &     ARG_L2(rhs)+1 .lt. ARG_L2(dat) .or. 
     $     ARG_H2(rhs)-1 .gt. ARG_H2(dat)) then
         write(6,*) "FORT_HG_CELL_TO_NODE: bad index limits"
         stop
      end if

      if (is_rz.ne.1) then
         do j=ARG_L2(rhs)+nghost,ARG_H2(rhs)-nghost
            do i=ARG_L1(rhs)+nghost,ARG_H1(rhs)-nghost
               rhs(i,j) = fourth*(dat(i-1,j-1)+dat(i-1,j)+
     &                            dat(i  ,j-1)+dat(i  ,j) )
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

               rhs(i,j) = fourth*(rlo * (dat(i-1,j-1) + dat(i-1,j)) +
     &                            rhi * (dat(i  ,j-1) + dat(i  ,j)))
            end do
         end do
         factor = dr/24.0
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
               rhs(i,j) = rhs(i,j) + factor *
     &              (rlo * (dat(i-1,j-1) + dat(i-1,j)) -
     &               rhi * (dat(i  ,j-1) + dat(i  ,j)))
            end do
         end do

#elif BL_PRVERSION == 5
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
               rhs(i,j) = fourth*(rlo * (dat(i-1,j-1) + dat(i-1,j)) +
     &                            rhi * (dat(i  ,j-1) + dat(i  ,j)))
               if (i .eq. 0) rhs(i,j) = half * rhs(i,j)
            end do
         end do
#endif
      end if

      end

c ::: -----------------------------------------------------------
c ::: This routine will compute D (U / dt + sigma G^perp phi)
c :::
c ::: NOTE: phi is node-centered
c :::       sigma is cell-centered      
c :::       U is cell-centered      
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: phi        => node centered data
c ::: DIMS(phi)  => index limits for phi
c ::: sigma      => cell-centered coefficient array
c ::: DIMS(sigma)=> index limits for sigma      
c ::: u          => cell-centered velocity array
c ::: DIMS(u)    => index limits for u      
c ::: divgp     <=  node-centered divergence of gradient
c ::: DIMS(divgp)=> index limits for divgp
c ::: lo,hi      => node-centered subregion of divgp to define
c ::: dx         => cell sizes      
c ::: mult       => scalar multiplicative factor      
c ::: -----------------------------------------------------------

      subroutine FORT_FILTRHS(phi,DIMS(phi),sigma,DIMS(sigma),
     &                        divgp,DIMS(divgp),
     &                        lo,hi,dx,mult,is_rz)
      implicit none

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
      REAL_T     p1, p2, p3, p4, denom
      REAL_T     rminv,rpinv,fac
      REAL_T     gpy_ul,gpy_ur,gpy_ll,gpy_lr
      REAL_T     pext, vext

#if BL_PRVERSION==9
c     FOR NINE-POINT STENCIL
      denom = one/(twelve*dx(1)**2) + one/(twelve*dx(2)**2) 
#elif BL_PRVERSION==5
c     FOR FIVE-POINT STENCIL
      denom = one/(four*dx(1)**2) + one/(four*dx(2)**2)
#else
#error "Unknown PRVERSION setting"
#endif

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
         enddo
      enddo

      if (is_rz .eq. 1) then

c      extra term added to DGphi for cylindrical coordinates
       do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            rpinv = half / (dfloat(i)+half)
            rminv = half / (dfloat(i)-half)

            pext = ((sigma(i-1,j-1) * (phi(i,j-1) - phi(i,j)) +
     $               sigma(i-1,j  ) * (phi(i,j+1) - phi(i,j))) * rminv -
     $              (sigma(i  ,j-1) * (phi(i,j-1) - phi(i,j)) +
     $               sigma(i  ,j  ) * (phi(i,j+1) - phi(i,j))) * rpinv)

            divgp(i,j) = divgp(i,j) + denom*pext
         enddo
       enddo

c      extra term added to DGbarphi for cylindrical coordinates --
c        here we subtract it since we want D G^perp phi
       fac = one / (twelve * dx(2)**2)
       do j = lo(2), hi(2)
       do i = lo(1), hi(1)
           rminv = fac / (dfloat(i)-half)
           rpinv = fac / (dfloat(i)+half)
   
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
        enddo
        enddo

      endif

      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
           divgp(i,j) = divgp(i,j) * mult
      enddo
      enddo

      end

