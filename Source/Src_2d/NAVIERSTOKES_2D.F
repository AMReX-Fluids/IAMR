
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <NAVIERSTOKES_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 2

      subroutine FORT_GRADP (
     &     p,DIMS(p),
     &     gp,DIMS(gp),
     &     lo,hi,dx,is_full)
c 
c     Compute a cell centered gradient from a node
c     centered field.  Returns all components of GRADP
c
      implicit none
      integer    DIMDEC(p)
      integer    DIMDEC(gp)
      integer    lo(SDIM), hi(SDIM)
      REAL_T     dx(SDIM)
      REAL_T     p(DIMV(p))
      REAL_T     gp(DIMV(gp),SDIM)
      integer    is_full
      integer    i,j
      REAL_T     ddx, ddy

      ddx = half/dx(1)
      ddy = half/dx(2)

      if (is_full .eq. 0) then
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          gp(i,j,1) = ddx*(p(i+1,j)-p(i,j)+p(i+1,j+1)-p(i,j+1))
          gp(i,j,2) = ddy*(p(i,j+1)-p(i,j)+p(i+1,j+1)-p(i+1,j))
        end do
        end do
      else
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          gp(i,j,1) = (p(i+1,j)-p(i,j)+p(i+1,j+1)-p(i,j+1))
          gp(i,j,2) = (p(i,j+1)-p(i,j)+p(i+1,j+1)-p(i+1,j))
        end do
        end do
      endif

      end

c :: ----------------------------------------------------------
c :: Replace coarse grid pressure data with corresponding
c :: fine grid pressure data.
c ::
c :: INPUTS / OUTPUTS:
c ::  crse      <=  coarse grid data
c ::  DIMS(crse) => index limits of crse
c ::  fine       => fine grid data
c ::  DIMS(fine) => index limits of fine 
c ::  lo,hi      => index limits of overlap (crse grid)
c ::  ratios     => IntVect refinement ratio
c ::
c :: NOTE:
c ::  Assumes pressure fields are node based
c :: ----------------------------------------------------------
c ::
      subroutine FORT_PUTDOWN (crse,DIMS(crse),
     &			       fine,DIMS(fine),lo,hi,ratios)
      implicit none
      integer  DIMDEC(crse)
      integer  DIMDEC(fine)
      integer  lo(2), hi(2)
      integer  ratios(2)
      REAL_T   crse(DIMV(crse))
      REAL_T   fine(DIMV(fine))

      integer  ic, jc
      integer  lratx, lraty

      lratx = ratios(1)
      lraty = ratios(2)

      do jc = lo(2), hi(2)
         do ic = lo(1), hi(1)
            crse(ic,jc) = fine(lratx*ic,lraty*jc)
         end do
      end do

      end

c :: ----------------------------------------------------------
c :: UTILITY ROUTINE: compute:
c ::             A += alpha*B on subrange
c ::
c :: INPUTS / OUTPUTS:
c ::  a         <=  output array
c ::  b          => input array
c ::  DIMS(a)    => index limits of a array
c ::  DIMS(b)    => index limits of a array
c ::  lo,hi      => index limits of update region
c ::  alpha      => multiplicative factor
c :: ----------------------------------------------------------
c ::
       subroutine FORT_INCRMULT(a,DIMS(a),b,DIMS(b),lo,hi,alpha)

       implicit none
       integer    DIMDEC(a)
       integer    DIMDEC(b)
       integer    lo(2), hi(2)
       REAL_T     alpha
       REAL_T     a(DIMV(a))
       REAL_T     b(DIMV(b))

       integer i, j

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             a(i,j) = a(i,j) + alpha*b(i,j)
          end do
       end do

       end

c :: ----------------------------------------------------------
c :: MAXVAL
c ::             maxval = max{ rho(i,j) }
c ::
c :: ----------------------------------------------------------
c ::
       subroutine FORT_MAXVAL(rho,DIMS(rho),DIMS(grid),mxval)

       implicit none
       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  rho(DIMV(rho))
       REAL_T  mxval

       integer i,j

       mxval = -Huge(zero)

       do j = ARG_L2(grid), ARG_H2(grid)
          do i = ARG_L1(grid), ARG_H1(grid)
             mxval = max(mxval, rho(i,j))
          end do
       end do

       end

c :: ----------------------------------------------------------
c :: SUMMASS
c ::             MASS = sum{ vol(i,j)*rho(i,j) }
c ::
c :: INPUTS / OUTPUTS:
c ::  rho        => density field
c ::  DIMS(rho)  => index limits of rho aray
c ::  lo,hi      => index limits of grid interior
c ::  dx	 => cell size
c ::  mass      <=  total mass
c ::  r		 => radius at cell center
c ::  irlo,hi    => index limits of r array
c ::  rz_flag    => == 1 if R_Z coords
c ::  tmp        => temp column array
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SUMMASS(rho,DIMS(rho),DIMS(grid),dx,mass,
     &                         r,irlo,irhi,rz_flag)

       implicit none
       integer irlo, irhi, rz_flag
       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  mass, dx(2)
       REAL_T  rho(DIMV(rho))
       REAL_T  r(irlo:irhi)

       integer i, j
       REAL_T  dr, dz, vol

       dr = dx(1)
       dz = dx(2)

       mass = zero

       do i = ARG_L1(grid), ARG_H1(grid)
          vol = dr*dz
	  if (rz_flag .eq. 1) vol = vol*two*Pi*r(i)
          do j = ARG_L2(grid), ARG_H2(grid)
	     mass = mass + vol*rho(i,j)
	  end do
       end do

       end

c :: ----------------------------------------------------------
c :: SUMMASSCYL
c ::    MASS = sum{ vol(i,j)*rho(i,j) } over subregion cylinder
c ::
c :: INPUTS / OUTPUTS:
c ::  rho        => density field
c ::  DIMS(rho)  => index limits of rho aray
c ::  lo,hi      => index limits of grid interior
c ::  dx	 => cell size
c ::  mass      <=  total mass
c ::  r		 => radius at cell center
c ::  irlo,hi    => index limits of r array
c ::  rz_flag    => == 1 if R_Z coords
c ::  plo        => domain lo end
c ::  vws_dz     => height of subregion, from plo
c ::  vws_R      => radius of subregion
c ::  tmp        => temp column array
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SUMMASS_CYL(rho,DIMS(rho),DIMS(grid),dx,mass,
     &                             r,irlo,irhi,rz_flag,plo,vws_dz,vws_R)

       implicit none
       integer irlo, irhi, rz_flag
       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  mass, dx(2)
       REAL_T  rho(DIMV(rho))
       REAL_T  r(irlo:irhi)
       REAL_T  plo(2),vws_dz, vws_R

       integer i, j
       REAL_T  dr, dz, vol, rr

       dr = dx(1)
       dz = dx(2)

       mass = zero

       do i = ARG_L1(grid), ARG_H1(grid)
          if (rz_flag .eq. 1) then
             rr = r(i)
          else
             rr = plo(1) + (i+half)*dx(1)
          endif

          if (rr.le.vws_R) then 
             vol = dr*dz
             if (rz_flag .eq. 1) vol = vol*two*Pi*rr
             do j = ARG_L2(grid), ARG_H2(grid)
                mass = mass + vol*rho(i,j)
             end do
          endif
       end do

       end

c ::
c :: ----------------------------------------------------------
c :: This routine fills an edge-centered fab from a cell-centered
c :: fab using simple linear interpolation.
c ::
c :: INPUTS / OUTPUTS:
c ::  lo,hi      => index limits of the of the cell-centered fab
c ::  DIMS(cfab) => index limits of the cell-centered fab
c ::  cfab       => cell-centered data
c ::  DIMS(efab) => index limits of the edge-centered fab
c ::  efab       => edge-centered fab to fill
c ::  nc         => Number of components in the fab to fill
c ::  dir        => direction data needs to be shifted to get to edges
c :: ----------------------------------------------------------
c ::
      subroutine FORT_CEN2EDG(lo, hi, 
     &     DIMS(cfab), cfab,
     &     DIMS(efab), efab, nc, dir,
     &     isharm)
      implicit none
      integer lo(SDIM), hi(SDIM), nc, dir, isharm
      integer DIMDEC(cfab)
      integer DIMDEC(efab)
      REAL_T  cfab(DIMV(cfab), nc)
      REAL_T  efab(DIMV(efab), nc)

      integer i,j,n

      if ( isharm .eq. 0 ) then
         if (dir .EQ. 0) then
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1)+1, hi(1)
                     efab(i,j,n) = half*(cfab(i,j,n) + cfab(i-1,j,n))
                  end do
               end do
            end do
         else
            do n = 1,nc
               do j = lo(2)+1, hi(2)
                  do i = lo(1), hi(1)
                     efab(i,j,n) = half*(cfab(i,j,n) + cfab(i,j-1,n))
                  end do
               end do
            end do
         end if
      else
         if (dir .EQ. 0) then
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1)+1, hi(1)
                     if((cfab(i,j,n) * cfab(i-1,j,n)).gt.zero)then
                        efab(i,j,n)
     &                       = 2*(cfab(i,j,n) * cfab(i-1,j,n))/
     &                       (cfab(i,j,n) + cfab(i-1,j,n))
                     else
                        efab(i,j,n)=zero
                     endif
                  end do
               end do
            end do
         else
            do n = 1,nc
               do j = lo(2)+1, hi(2)
                  do i = lo(1), hi(1)
                     if((cfab(i,j,n) * cfab(i,j-1,n)).gt.zero)then
                        efab(i,j,n)
     &                       = 2*(cfab(i,j,n) * cfab(i,j-1,n))/
     &                       (cfab(i,j,n) + cfab(i,j-1,n))
                     else
                        efab(i,j,n)=zero
                     endif
                  end do
               end do
            end do
         end if
      end if
      end

c-----------------------------------------------------------------------
      subroutine EDGE_INTERP(flo, fhi, nc, ratio, dir,
     &     fine, fine_l0, fine_l1, fine_h0, fine_h1)
      implicit none
      integer flo(0:2-1), fhi(0:2-1), nc, ratio(0:2-1), dir
      integer fine_l0, fine_l1, fine_h0, fine_h1
      DOUBLE PRECISION fine(fine_l0:fine_h0,fine_l1:fine_h1,nc)
      integer i,j,n,P,M
      DOUBLE PRECISION val, df

c     Do linear in dir, pc transverse to dir, leave alone the fine values
c     lining up with coarse edges--assume these have been set to hold the 
c     values you want to interpolate to the rest.
      if (dir.eq.0) then
         do n=1,nc
            do j=flo(1),fhi(1),ratio(1)
               do i=flo(0),fhi(0)-ratio(dir),ratio(0)
                  df = fine(i+ratio(dir),j,n)-fine(i,j,n)
                  do M=1,ratio(dir)-1
                     val = fine(i,j,n) + df*dble(M)/dble(ratio(dir))
                     do P=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                        fine(i+M,P,n) = val
                     enddo
                  enddo                     
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do j=flo(1),fhi(1)-ratio(dir),ratio(1)
               do i=flo(0),fhi(0)
                  df = fine(i,j+ratio(dir),n)-fine(i,j,n)
                  do M=1,ratio(dir)-1
                     val = fine(i,j,n) + df*dble(M)/dble(ratio(dir))
                     do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                        fine(P,j+M,n) = val
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif

      end
c-----------------------------------------------------------------------
      subroutine PC_EDGE_INTERP(lo, hi, nc, ratio, dir,
     &     crse, crse_l0, crse_l1, crse_h0, crse_h1,
     &     fine, fine_l0, fine_l1, fine_h0, fine_h1)
      implicit none
      integer lo(2),hi(2), nc, ratio(0:2-1), dir
      integer crse_l0, crse_l1, crse_h0, crse_h1
      integer fine_l0, fine_l1, fine_h0, fine_h1
      DOUBLE PRECISION crse(crse_l0:crse_h0,crse_l1:crse_h1,nc)
      DOUBLE PRECISION fine(fine_l0:fine_h0,fine_l1:fine_h1,nc)
      integer i,j,ii,jj,n,L

c     For edge-based data, fill fine values with piecewise-constant interp of coarse data.
c     Operate only on faces that overlap--ie, only fill the fine faces that make up each
c     coarse face, leave the in-between faces alone.
      if (dir.eq.0) then
         do n=1,nc
            do j=lo(2),hi(2)
               jj = ratio(1)*j
               do i=lo(1),hi(1)
                  ii = ratio(0)*i
                  do L=0,ratio(1)-1
                     fine(ii,jj+L,n) = crse(i,j,n)
                  enddo
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do j=lo(2),hi(2)
               jj = ratio(1)*j
               do i=lo(1),hi(1)
                  ii = ratio(0)*i
                  do L=0,ratio(0)-1
                     fine(ii+L,jj,n) = crse(i,j,n)
                  enddo
               enddo
            enddo
         enddo
      endif

      end

      subroutine FORT_BINGHAM(visc, DIMS(visc), vel, DIMS(vel), 
     &     lo, hi, domlo, domhi, delta, bc,
     &     mu, tau, eps)
      implicit none
c
c ::: This routine will compute the Papanastasiou regularised Bingham
c ::: viscosity based on the velocity field
c
      integer    lo(2), hi(2)
      integer    DIMDEC(visc)
      integer    DIMDEC(vel)
      integer    domlo(2), domhi(2)
      integer    bc(2,2,2)
      REAL_T     delta(2)
      REAL_T     visc(DIMV(visc),1)
      REAL_T     vel(DIMV(vel),2)
      REAL_T     mu, tau, eps

      integer   i,j
      REAL_T    ux, vx, uy, vy, dx, dy
      logical   fixlft, fixrgt, fixbot, fixtop
      REAL_T    uxcen, uxlft, uxrgt, vxcen, vxlft, vxrgt
      REAL_T    uycen, uybot, uytop, vycen, vybot, vytop
      REAL_T    strnrt, strnrt_fun, bingham_fun
c
c     ::::: some useful macro definitions
c
#     define U(i,j) dat(i,j,1)
#     define V(i,j) dat(i,j,2)
#     define VLOX bc(1,1,2)
#     define VHIX bc(1,2,2)
#     define ULOY bc(2,1,1)
#     define UHIY bc(2,2,1)
c
c     ::::: statement functions that implement stencil
c
      uxcen(i,j) = half*(U(i+1,j)-U(i-1,j))/dx
      uxlft(i,j) = (U(i+1,j)+three*U(i,j)-four*U(i-1,j))/(three*dx)
      uxrgt(i,j) = (four*U(i+1,j)-three*U(i,j)-U(i-1,j))/(three*dx)
      vxcen(i,j) = half*(V(i+1,j)-V(i-1,j))/dx
      vxlft(i,j) = (V(i+1,j)+three*V(i,j)-four*V(i-1,j))/(three*dx)
      vxrgt(i,j) = (four*V(i+1,j)-three*V(i,j)-V(i-1,j))/(three*dx)
      uycen(i,j) = half*(U(i,j+1)-U(i,j-1))/dy
      uybot(i,j) = (U(i,j+1)+three*U(i,j)-four*U(i,j-1))/(three*dy)
      uytop(i,j) = (four*U(i,j+1)-three*U(i,j)-U(i,j-1))/(three*dy)
      vycen(i,j) = half*(V(i,j+1)-V(i,j-1))/dy
      vybot(i,j) = (V(i,j+1)+three*V(i,j)-four*V(i,j-1))/(three*dy)
      vytop(i,j) = (four*V(i,j+1)-three*V(i,j)-V(i,j-1))/(three*dy)

c
c     ::::: function for computing magnitude of rate of strain tensor
c
      strnrt_fun(ux,vx,uy,vy) = sqrt(ux**2+half*(vx+uy)**2+vy**2)
      bingham_fun(x) = mu+half*tau*(1-exp(-x/eps))/x

      dx = delta(1)
      dy = delta(2)

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            ux  = uxcen(i,j)
            vx  = vxcen(i,j)
            uy  = uycen(i,j)
            vy  = vycen(i,j)
            strnrt = strnrt_fun(ux,vx,uy,vy)
            visc(i,j,1) =
     &          bingham_fun(visc_coef,yield_stress,strnrt,reg_param)
         end do
      end do

      fixlft = ( (lo(1) .eq. domlo(1)) .and. 
     &           (VLOX .eq. EXT_DIR .or. VLOX .eq. HOEXTRAP) )
      fixrgt = ( (hi(1) .eq. domhi(1)) .and. 
     &           (VHIX .eq. EXT_DIR .or. VHIX .eq. HOEXTRAP) )
      fixbot = ( (lo(2) .eq. domlo(2)) .and. 
     &           (ULOY .eq. EXT_DIR .or. ULOY .eq. HOEXTRAP) )
      fixtop = ( (hi(2) .eq. domhi(2)) .and. 
     &           (UHIY .eq. EXT_DIR .or. UHIY .eq. HOEXTRAP) )
c
c     ::::: handle special bndry conditions at the left edge
c
      if (fixlft) then
         i = lo(1)
         do j = lo(2), hi(2)
            ux  = uxlft(i,j)
            vx  = vxlft(i,j)
            uy  = uycen(i,j)
            vy  = vycen(i,j)
            strnrt = strnrt_fun(ux,vx,uy,vy)
            visc(i,j,1) =
     &          bingham_fun(visc_coef,yield_stress,strnrt,reg_param)
         end do
      end if
c
c     ::::: handle special bndry conditions on the right
c
      if (fixrgt) then
         i = hi(1)
         do j = lo(2), hi(2)
            ux  = uxrgt(i,j)
            vx  = vxrgt(i,j)
            uy  = uycen(i,j)
            vy  = vycen(i,j)
            strnrt = strnrt_fun(ux,vx,uy,vy)
            visc(i,j,1) =
     &          bingham_fun(visc_coef,yield_stress,strnrt,reg_param)
         end do
      end if
c
c     ::::: handle special bndry conditions on bottom
c
      if (fixbot) then
         j = lo(2)
         do i = lo(1), hi(1)
            ux  = uxcen(i,j)
            vx  = vxcen(i,j)
            uy  = uybot(i,j)
            vy  = vybot(i,j)
            strnrt = strnrt_fun(ux,vx,uy,vy)
            visc(i,j,1) =
     &          bingham_fun(visc_coef,yield_stress,strnrt,reg_param)
         end do
      end if
c
c     ::::: handle special bndry conditions on top
c
      if (fixtop) then
         j = hi(2)
         do i = lo(1), hi(1)
            ux  = uxcen(i,j)
            vx  = vxcen(i,j)
            uy  = uytop(i,j)
            vy  = vytop(i,j)
            strnrt = strnrt_fun(ux,vx,uy,vy)
            visc(i,j,1) =
     &          bingham_fun(visc_coef,yield_stress,strnrt,reg_param)
         end do
      end if
c
c     ::::: check corners
c
      if (fixlft .and. fixbot) then
         i = lo(1)
         j = lo(2)
         ux = uxlft(i,j)
         vx = vxlft(i,j)
         uy = uybot(i,j)
         vy = vybot(i,j)
         strnrt = strnrt_fun(ux,vx,uy,vy)
         visc(i,j,1) =
     &       bingham_fun(visc_coef,yield_stress,strnrt,reg_param)
      end if
      if (fixlft .and. fixtop) then
         i = lo(1)
         j = hi(2)
         ux = uxlft(i,j)
         vx = vxlft(i,j)
         uy = uytop(i,j)
         vy = vytop(i,j)
         strnrt = strnrt_fun(ux,vx,uy,vy)
         visc(i,j,1) =
     &       bingham_fun(visc_coef,yield_stress,strnrt,reg_param)
      end if
      if (fixrgt .and. fixtop) then
         i = hi(1)
         j = hi(2)
         ux = uxrgt(i,j)
         vx = vxrgt(i,j)
         uy = uytop(i,j)
         vy = vytop(i,j)
         strnrt = strnrt_fun(ux,vx,uy,vy)
         visc(i,j,1) =
     &       bingham_fun(visc_coef,yield_stress,strnrt,reg_param)
      end if
      if (fixrgt .and. fixbot) then
         i = hi(1)
         j = lo(2)
         ux = uxrgt(i,j)
         vx = vxrgt(i,j)
         uy = uybot(i,j)
         vy = vybot(i,j)
         strnrt = strnrt_fun(ux,vx,uy,vy)
         visc(i,j,1) =
     &       bingham_fun(visc_coef,yield_stress,strnrt,reg_param)
      end if

#     undef U
#     undef V      
#     undef VLOX
#     undef VHIX
#     undef ULOY
#     undef UHIY

      end

