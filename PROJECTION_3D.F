c
c $Id: PROJECTION_3D.F,v 1.14 2004-09-08 20:58:54 almgren Exp $
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

       twelfth = one / twelve

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
     $             ARG_L2(grid)-ng:ARG_H2(grid)+ng,
     $             ARG_L3(grid)-ng:ARG_H3(grid)+ng)
      REAL_T     r(ARG_L1(grid)-nr:ARG_H1(grid)+nr)

      call bl_abort('SHOULD NEVER BE IN RADMPY IN 3D')

      end

      subroutine FORT_RADDIV(a,DIMS(grid),ng,r,nr,n)
c 
c     divide A by a Radius r
c
      integer    ng,nr
      integer    n(SDIM)
      integer    DIMDEC(grid)
      REAL_T     a(ARG_L1(grid)-ng:ARG_H1(grid)+ng,
     $             ARG_L2(grid)-ng:ARG_H2(grid)+ng,
     $             ARG_L3(grid)-ng:ARG_H3(grid)+ng)
      REAL_T     r(ARG_L1(grid)-nr:ARG_H1(grid)+nr)

      call bl_abort('SHOULD NEVER BE IN RADDIV IN 3D')
      end

      subroutine FORT_ANELCOEFFMPY(a,DIMS(grid),domlo,domhi,ng,rho0,nr,bogus_value,mult)
c 
c     multiply A by rho0
c
c 
c     NOTE: THIS ROUTINE HAS BEEN MODIFIED SO THAT ALL VALUES
c           OUTSIDE THE DOMAIN ARE SET TO BOGUS VALUE
c
      integer    ng,nr
      integer    DIMDEC(grid)
      integer    domlo(3), domhi(3)
      REAL_T     a(ARG_L1(grid)-ng:ARG_H1(grid)+ng, 
     $             ARG_L2(grid)-ng:ARG_H2(grid)+ng,
     $             ARG_L3(grid)-ng:ARG_H3(grid)+ng)
      REAL_T     rho0(ARG_L3(grid)-nr:ARG_H3(grid)+nr)
      REAL_T     bogus_value
      integer    mult

      integer i, j, k

      if (mult .eq. 1) then
         do k = ARG_L3(grid)-ng, ARG_H3(grid)+ng
         do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do i = ARG_L1(grid)-ng, ARG_H1(grid)+ng
           a(i,j,k) = a(i,j,k) * rho0(k)
         end do
         end do
         end do
      else if (mult .eq. 0) then
         do k = ARG_L3(grid)-ng, ARG_H3(grid)+ng
         do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do i = ARG_L1(grid)-ng, ARG_H1(grid)+ng
           a(i,j,k) = a(i,j,k) / rho0(k)
         end do
         end do
         end do
      else 
         print *,'BOGUS MULT IN ANELCOEFFMULT ',mult
         stop
      end if

      if (ARG_L3(grid)-ng .lt. domlo(3)) then
         do i = ARG_L1(grid)-ng, ARG_H1(grid)+ng
         do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do k = ARG_L3(grid)-ng, domlo(3)-1
           a(i,j,k) = bogus_value
         end do
         end do
         end do
      end if

      if (ARG_H3(grid)+ng .gt. domhi(3)) then
         do i = ARG_L1(grid)-ng, ARG_H1(grid)+ng
         do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do k = domhi(3)+1, ARG_H3(grid)+ng
           a(i,j,k) = bogus_value
         end do
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

      integer i,j,k

      if (ARG_H1(ccdat)   .lt. lo(1) .or. 
     &    ARG_L1(ccdat)   .gt. hi(1) .or. 
     &    ARG_H2(ccdat)   .lt. lo(2) .or. 
     &    ARG_L2(ccdat)   .gt. hi(2) .or. 
     &    ARG_H3(ccdat)   .lt. lo(3) .or. 
     &    ARG_L3(ccdat)   .gt. hi(3) .or. 
     &    ARG_H1(nodedat) .lt. lo(1)+1 .or. 
     &    ARG_L1(nodedat) .gt. hi(1) .or. 
     &    ARG_H2(nodedat) .lt. lo(2)+1 .or. 
     &    ARG_L2(nodedat) .gt. hi(2) .or.
     &    ARG_H3(nodedat) .lt. lo(3)+1 .or. 
     &    ARG_L3(nodedat) .gt. hi(3) ) then 
        call bl_abort("FORT_HG_CELL_TO_NODE: bad index limits")
      end if

      if(isrz.eq.1.and.lrweighted.ne.1)then
        call bl_abort('FORT_HGN2C : isrz=1 and lrweighted!=1 not implemented')
      end if

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               ccdat(i,j,k)=eighth*(nodedat(i  ,j,  k  )+nodedat(i+1,j  ,k  )+
     &                              nodedat(i  ,j+1,k  )+nodedat(i+1,j+1,k  )+
     &                              nodedat(i  ,j,  k+1)+nodedat(i+1,j  ,k+1)+
     &                              nodedat(i  ,j+1,k+1)+nodedat(i+1,j+1,k+1))
            end do
         end do
      end do

      end


       subroutine FORT_HGC2N(
     &     nghost, DIMS(dat), dat, rcen,
     &     DIMS(rhs), rhs,
     &     domlo, domhi, lowfix, hifix, dr, is_rz) 
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
c     rhslo,rhshi => index extents of rhs
c     rhs         <= node centered array with results
c     rcen,lowfix,hifix are artifacts of 2D r-z version and are not used
c     ----------------------------------------------------------
c 
      implicit none
      integer nghost
      integer domlo(SDIM), domhi(SDIM)
      integer DIMDEC(dat)
      integer DIMDEC(rhs)
      REAL_T  dat(DIMV(dat))
      REAL_T  rhs(DIMV(rhs))

c not used
      integer is_rz, lowfix, hifix
      REAL_T  rcen(DIM1(dat))
      REAL_T  dr

c local variables
      integer i, j, k

      if (ARG_L1(rhs)+1 .lt. ARG_L1(dat) .or. 
     $    ARG_H1(rhs)-1 .gt. ARG_H1(dat) .or.
     &    ARG_L2(rhs)+1 .lt. ARG_L2(dat) .or. 
     $    ARG_H2(rhs)-1 .gt. ARG_H2(dat) .or.
     &    ARG_L3(rhs)+1 .lt. ARG_L3(dat) .or. 
     $    ARG_H3(rhs)-1 .gt. ARG_H3(dat) ) then
         call bl_abort("FORT_HG_CELL_TO_NODE: bad index limits")
      end if

      do k=ARG_L3(rhs)+nghost,ARG_H3(rhs)-nghost
         do j=ARG_L2(rhs)+nghost,ARG_H2(rhs)-nghost
            do i=ARG_L1(rhs)+nghost,ARG_H1(rhs)-nghost
               rhs(i,j,k)=eighth*(dat(i-1,j-1,k  )+dat(i-1,j,k  )+
     &                            dat(i  ,j-1,k  )+dat(i  ,j,k  )+
     &                            dat(i-1,j-1,k-1)+dat(i-1,j,k-1)+
     &                            dat(i  ,j-1,k-1)+dat(i  ,j,k-1))
            end do
         end do
      end do

      end


c ::: -----------------------------------------------------------
c ::: This routine will compute D (U / dt + sigma G^perp phi)
c :::
c ::: NOTE: phi is node-centered
c :::       sigma is cell-centered      
c :::       U is cell-centered      
c :::
c ::: As of now, this routine assumes dx == dy == dz
c ::: Eventually, this will be changed..
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

      integer    DIMDEC(phi)
      integer    DIMDEC(sigma)
      integer    DIMDEC(divgp)
      integer    lo(SDIM), hi(SDIM)
      REAL_T     dx(SDIM), mult
      REAL_T     phi(DIMV(phi))
      REAL_T     sigma(DIMV(sigma))
      REAL_T     divgp(DIMV(divgp))
      integer    is_rz

      integer    i, j, k
      REAL_T     p1, p2, p3, p4, p5, p6, p7, p8, denom

c     NOTE: the factor of mult here accounts for the fact that one
c           takes ratio time steps on the fine grid for every one
c           timestep on the coarse grid.  mult = 1/ratio.
      
      denom = -mult/(sixteen*dx(1)**2)

c check and make sure dx == dy == dz
      if ( dx(1) .NE. dx(2) .or. dx(1) .NE. dx(3) ) then
         call bl_abort('FILTRHS:  code for dx != dy != dz not implemented')
      endif

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               p1 = sigma(i,j,k)*(
     &              five*(phi(i+1,j,k) + phi(i,j+1,k) + phi(i,j,k+1)) 
     &              - phi(i+1,j+1,k) - phi(i,j+1,k+1) - phi(i+1,j,k+1) 
     &              - three*phi(i+1,j+1,k+1) - nine*phi(i,j,k) )
               p2 = sigma(i-1,j,k)*(
     &              five*(phi(i-1,j,k) + phi(i,j+1,k) + phi(i,j,k+1)) 
     &              - phi(i-1,j+1,k) - phi(i,j+1,k+1) - phi(i-1,j,k+1) 
     &              - three*phi(i-1,j+1,k+1) - nine*phi(i,j,k) )
               p3 = sigma(i,j-1,k)*(
     &              five*(phi(i+1,j,k) + phi(i,j-1,k) + phi(i,j,k+1)) 
     &              - phi(i+1,j-1,k) - phi(i,j-1,k+1) - phi(i+1,j,k+1) 
     &              - three*phi(i+1,j-1,k+1) - nine*phi(i,j,k) )
               p4 = sigma(i-1,j-1,k)*(
     &              five*(phi(i-1,j,k) + phi(i,j-1,k) + phi(i,j,k+1)) 
     &              - phi(i-1,j-1,k) - phi(i,j-1,k+1) - phi(i-1,j,k+1) 
     &              - three*phi(i-1,j-1,k+1) - nine*phi(i,j,k) )
               p5 = sigma(i,j,k-1)*(
     &              five*(phi(i+1,j,k) + phi(i,j+1,k) + phi(i,j,k-1)) 
     &              - phi(i+1,j+1,k) - phi(i,j+1,k-1) - phi(i+1,j,k-1) 
     &              - three*phi(i+1,j+1,k-1) - nine*phi(i,j,k) )
               p6 = sigma(i-1,j,k-1)*(
     &              five*(phi(i-1,j,k) + phi(i,j+1,k) + phi(i,j,k-1)) 
     &              - phi(i-1,j+1,k) - phi(i,j+1,k-1) - phi(i-1,j,k-1) 
     &              - three*phi(i-1,j+1,k-1) - nine*phi(i,j,k) )
               p7 = sigma(i,j-1,k-1)*(
     &              five*(phi(i+1,j,k) + phi(i,j-1,k) + phi(i,j,k-1)) 
     &              - phi(i+1,j-1,k) - phi(i,j-1,k-1) - phi(i+1,j,k-1) 
     &              - three*phi(i+1,j-1,k-1) - nine*phi(i,j,k) )
               p8 = sigma(i-1,j-1,k-1)*(
     &              five*(phi(i-1,j,k) + phi(i,j-1,k) + phi(i,j,k-1)) 
     &              - phi(i-1,j-1,k) - phi(i,j-1,k-1) - phi(i-1,j,k-1) 
     &              - three*phi(i-1,j-1,k-1) - nine*phi(i,j,k) )
               divgp(i,j,k) = denom*(p1+p2+p3+p4+p5+p6+p7+p8)

            end do
         end do
      end do

      end





