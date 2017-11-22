
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <NAVIERSTOKES_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 3

      subroutine FORT_GRADP(
     &     p,DIMS(p),
     &     gp,DIMS(gp),
     &     lo,hi,dx)
c ::
c :: ----------------------------------------------------------
c :: Compute a cell centered gradient from a node
c :: centered field.  Returns all components of GRADP
c :: ----------------------------------------------------------
c ::
      implicit none

      integer DIMDEC(p)  
      integer DIMDEC(gp)  
      integer lo(SDIM),  hi(SDIM)
      REAL_T  dx(SDIM)
      REAL_T  p(DIMV(p))
      REAL_T  gp(DIMV(gp),SDIM)
      integer i,j,k
      REAL_T  ddx, ddy, ddz

      ddx = fourth/dx(1)
      ddy = fourth/dx(2)
      ddz = fourth/dx(3)

!$omp parallel private(i,j,k)
!$omp do
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               gp(i,j,k,1) = ddx*(
     &              p(i+1,j,k  )-p(i,j,k  )+p(i+1,j+1,k  )-p(i,j+1,k  )+
     &              p(i+1,j,k+1)-p(i,j,k+1)+p(i+1,j+1,k+1)-p(i,j+1,k+1))
            end do
         end do
      end do
!$omp end do nowait
!$omp do
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               gp(i,j,k,2) = ddy*(
     &              p(i,j+1,k  )-p(i,j,k  )+p(i+1,j+1,k  )-p(i+1,j,k  )+
     &              p(i,j+1,k+1)-p(i,j,k+1)+p(i+1,j+1,k+1)-p(i+1,j,k+1))
            end do
         end do
      end do
!$omp end do nowait
!$omp do
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               gp(i,j,k,3) = ddz*(
     &              p(i,  j,k+1)-p(i,  j,k)+p(i,  j+1,k+1)-p(i,  j+1,k)+
     &              p(i+1,j,k+1)-p(i+1,j,k)+p(i+1,j+1,k+1)-p(i+1,j+1,k))
            end do
         end do
      end do
!$omp end do nowait
!$omp end parallel

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
c ::  ratios     => refinement ratio
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
      integer  lo(SDIM), hi(SDIM)
      integer  ratios(SDIM)
      REAL_T   crse(DIMV(crse))
      REAL_T   fine(DIMV(fine))

      integer  ic, jc, kc, i, j, k
      integer  lratx, lraty, lratz

      lratx = ratios(1)
      lraty = ratios(2)
      lratz = ratios(3)

      do kc = lo(3), hi(3)
         k = lratz*kc
         do jc = lo(2), hi(2)
            j = lraty*jc
            do ic = lo(1), hi(1)
               i = lratx*ic
               crse(ic,jc,kc) = fine(i,j,k)
            end do
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
c ::  alo,ahi    => index limits of a array
c ::  blo,bhi    => index limits of a array
c ::  lo,hi      => index limits of update region
c ::  alpha      => multiplicative factor
c :: ----------------------------------------------------------
c ::
       subroutine FORT_INCRMULT(a,DIMS(a),b,DIMS(b),lo,hi,alpha)

       implicit none
      
       integer    DIMDEC(a)
       integer    DIMDEC(b)
       integer    lo(SDIM), hi(SDIM)
       REAL_T     alpha
       REAL_T     a(DIMV(a))
       REAL_T     b(DIMV(b))

       integer i, j, k

!$omp parallel do private(i,j,k)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                a(i,j,k) = a(i,j,k) + alpha*b(i,j,k)
             end do
          end do
       end do
!$omp end parallel do

       end

c ::
c :: ----------------------------------------------------------
c :: SUMTURB
c :: ----------------------------------------------------------
c ::
      subroutine FORT_SUMTURB(dat,pres,DIMS(dat),DIMS(pres),DIMS(grid),delta,
     &                        turb,ksize,turbVars)

      implicit none

      integer DIMDEC(dat)
      integer DIMDEC(pres)
      integer DIMDEC(grid)
      integer ksize, turbVars
      REAL_T  delta(SDIM)
      REAL_T  dat(DIMV(dat),16)
      REAL_T  pres(DIMV(pres),4)
      REAL_T  turb(0:ksize*turbVars-1)
      
      integer i, j, k
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      
      REAL_T area
      REAL_T rho, ux, uy, uz, p
      REAL_T drhodx, drhody, drhodz
      REAL_T duxdx, duxdy, duxdz
      REAL_T duydx, duydy, duydz
      REAL_T duzdx, duzdy, duzdz
      REAL_T dpdx, dpdy, dpdz
      REAL_T dx, dy, dz
      
      ilo = ARG_L1(grid)
      ihi = ARG_H1(grid)
      jlo = ARG_L2(grid)
      jhi = ARG_H2(grid)
      klo = ARG_L3(grid)
      khi = ARG_H3(grid)
      
      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      area = dx*dy
      
      do k = klo, khi
         do i = ilo, ihi
            do j = jlo, jhi

               rho  = dat(i,j,k,1)
               ux   = dat(i,j,k,2)
               uy   = dat(i,j,k,3)
               uz   = dat(i,j,k,4)

c     Here are the derivatives, can't do it here because of zeroed intersections
               drhodx= dat(i,j,k,5)
               duxdx = dat(i,j,k,6)
               duydx = dat(i,j,k,7)
               duzdx = dat(i,j,k,8)

               drhody= dat(i,j,k,9)
               duxdy = dat(i,j,k,10)
               duydy = dat(i,j,k,11)
               duzdy = dat(i,j,k,12)

               drhodz= dat(i,j,k,13)
               duxdz = dat(i,j,k,14)
               duydz = dat(i,j,k,15)
               duzdz = dat(i,j,k,16)

               p    = pres(i,j,k,1)
               dpdx = pres(i,j,k,2)
               dpdy = pres(i,j,k,3)
               dpdz = pres(i,j,k,4)

               turb(k*turbVars+0)  = turb(k*turbVars+0)  + area*rho
               turb(k*turbVars+1)  = turb(k*turbVars+1)  + area*rho*ux
               turb(k*turbVars+2)  = turb(k*turbVars+2)  + area*rho*uy
               turb(k*turbVars+3)  = turb(k*turbVars+3)  + area*rho*uz
               turb(k*turbVars+4)  = turb(k*turbVars+4)  + area*rho*ux*ux
               turb(k*turbVars+5)  = turb(k*turbVars+5)  + area*rho*ux*uy
               turb(k*turbVars+6)  = turb(k*turbVars+6)  + area*rho*ux*uz
               turb(k*turbVars+7)  = turb(k*turbVars+7)  + area*rho*uy*uy
               turb(k*turbVars+8)  = turb(k*turbVars+8)  + area*rho*uy*uz
               turb(k*turbVars+9)  = turb(k*turbVars+9)  + area*rho*uz*uz
               turb(k*turbVars+10) = turb(k*turbVars+10) + area*rho*ux*ux*ux
               turb(k*turbVars+11) = turb(k*turbVars+11) + area*rho*ux*ux*uy
               turb(k*turbVars+12) = turb(k*turbVars+12) + area*rho*ux*ux*uz
               turb(k*turbVars+13) = turb(k*turbVars+13) + area*rho*ux*uy*uy
               turb(k*turbVars+14) = turb(k*turbVars+14) + area*rho*ux*uy*uz
               turb(k*turbVars+15) = turb(k*turbVars+15) + area*rho*ux*uz*uz
               turb(k*turbVars+16) = turb(k*turbVars+16) + area*rho*uy*uy*uy
               turb(k*turbVars+17) = turb(k*turbVars+17) + area*rho*uy*uy*uz
               turb(k*turbVars+18) = turb(k*turbVars+18) + area*rho*uz*uz*uz

            end do
         end do
      end do
      
      end

#ifdef SUMJET
c ::
c :: ----------------------------------------------------------
c :: SUMJET
c :: ----------------------------------------------------------
c ::
      subroutine FORT_SUMJET(dat,pres,DIMS(dat),DIMS(pres),DIMS(grid),delta,
     &                      jetData,levRsize,levKsize,rsize,ksize,jetVars,numSplit,
     &                      xlo,xhi);

      implicit none

      integer DIMDEC(dat)
      integer DIMDEC(pres)
      integer DIMDEC(grid)
      integer levRsize, levKsize, rsize, ksize, jetVars
      REAL_T  delta(SDIM)
      REAL_T  xlo(SDIM)
      REAL_T  xhi(SDIM)
      REAL_T  dat(DIMV(dat),27)
      REAL_T  pres(DIMV(pres),4)
      REAL_T  jetData(0:ksize*jetVars-1)
      REAL_T  gridDx
      
      integer i, j, k, v
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      
      REAL_T  rho, ux, uy, uz, p
      REAL_T  dpdx, dpdy, dpdz, dpdr, dpdt
      REAL_T  dx, dy, dz

      integer ii,jj,kk,kn,nn,numSplit,kklo,kkhi,ridx, idx
      REAL_T  dxnn,dynn,dznn
      REAL_T  xlow,ylow,zlow,xctr,yctr,zctr,xx,yy,zz,rr
      REAL_T  rhoctr, uxctr, uyctr, uzctr, tracctr, tempctr
      REAL_T  pctr, dpdxctr, dpdyctr, dpdzctr

      REAL_T  deltax, deltay, deltaz
      REAL_T  ur, ut, trac, temp
      
      REAL_T  mdx, mdy
      REAL_T  drx, dry, dtx, dty
      REAL_T  xxp, yyp, xxm, yym, rrp, rrm
      REAL_T  uxp, uyp, uxm, uym, urp, urm, utp, utm
      REAL_T  durdr, dutdr, durdt, dutdt

      REAL_T  KEK, KED, KEP

      integer idx_rho, idx_ux, idx_uy, idx_uz, idx_trac, idx_temp
      integer idx_p, idx_dpdx, idx_dpdy, idx_dpdz
      integer var
      
      integer isioproc
      
#include <probdata.H>

      call bl_pd_is_ioproc(isioproc)

c     Hack to zero
      isioproc = 0

      ilo = ARG_L1(grid)
      ihi = ARG_H1(grid)
      jlo = ARG_L2(grid)
      jhi = ARG_H2(grid)
      klo = ARG_L3(grid)
      khi = ARG_H3(grid)
      
      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      kn   = ksize/levKsize
      nn   = numSplit*kn
      dxnn = dx/dble(nn)
      dynn = dy/dble(nn)
c     Note this is right, it really should be kn not nn
      dznn = dz/dble(kn)
c     Mod dx
      mdx = sqrt(dxnn*dxnn+dynn*dynn)
      
      gridDx=dx/dble(kn)
      
      if (isioproc.eq.1) then 
         write (*,*) "In SUMJET:"
         write (*,*) "   ilo/ihi:",ilo,ihi
         write (*,*) "   jlo/jhi:",jlo,jhi
         write (*,*) "   klo/khi:",klo,khi
         write (*,*) "   xlo/xhi:",xlo(1),xhi(1)
         write (*,*) "   ylo/yhi:",xlo(2),xhi(2)
         write (*,*) "   zlo/zhi:",xlo(3),xhi(3)
         write (*,*) "   dx:",dx,dy,dz
         write (*,*) "   ksize:",ksize
         write (*,*) "   levKsize:",levKsize
         write (*,*) "   rsize:",rsize
         write (*,*) "   levRsize:",levRsize
         write (*,*) "   kn = ksize / levKsize:",kn
         write (*,*) "   nn = kn * numSplit:",nn
         write (*,*) "   dxnn = dx / nn:"
         write (*,*) "      ",dxnn,dynn,dznn
         write (*,*) "   dat l/h 1 = ", dat_l1, dat_h1
         write (*,*) "   dat l/h 2 = ", dat_l2, dat_h2
         write (*,*) "   dat l/h 3 = ", dat_l3, dat_h3
         write (*,*) "   pres l/h 1 = ", pres_l1, pres_h1
         write (*,*) "   pres l/h 2 = ", pres_l2, pres_h2
         write (*,*) "   pres l/h 3 = ", pres_l3, pres_h3
      endif

      idx_rho  = 1
      idx_ux   = 2
      idx_uy   = 3
      idx_uz   = 4
      idx_trac = 5
      idx_temp = 6

      idx_p    = 1
      idx_dpdx = 2
      idx_dpdy = 3
      idx_dpdz = 4

      do k = klo, khi
c     Bounds for vertical numerical integration
         kklo=k*kn
         kkhi=kklo+kn-1
c     Calculate left hand edge of cell
         zlow = xlo(3) + dz*dble(k-klo)
c     And x centre/er
         zctr = zlow + 0.5d0*dz
         
         do j = jlo, jhi
c     Calculate front of cell
            ylow = xlo(2) + dy*dble(j-jlo) - jet_y
c     And y centre/er
            yctr = ylow + 0.5d0*dy
            
            do i = ilo, ihi
c     Calculate left hand edge of cell
               xlow = xlo(1) + dx*dble(i-ilo) - jet_x
c     And x centre/er
               xctr = xlow + 0.5d0*dx
               
               rhoctr  = dat(i,j,k,idx_rho)
               uxctr   = dat(i,j,k,idx_ux)
               uyctr   = dat(i,j,k,idx_uy)
               uzctr   = dat(i,j,k,idx_uz)
               tracctr = dat(i,j,k,idx_trac)
               tempctr = dat(i,j,k,idx_temp)

               pctr    = pres(i,j,k,idx_p)
               dpdxctr = pres(i,j,k,idx_dpdx)
               dpdyctr = pres(i,j,k,idx_dpdy)
               dpdzctr = pres(i,j,k,idx_dpdz)

               if (rhoctr.gt.zero) then
c     We're not in a zeroed out bit, let's integrate...

c     Integrate numerically
                  do jj=1,nn
                     yy = ylow + dynn*(dble(jj)-half)
                     
                     do ii=1,nn
                        xx = xlow + dxnn*(dble(ii)-half)
                        
                        rr = sqrt(xx*xx+yy*yy)
                        
                        ridx = int(rr/gridDx)
                        
                        if (ridx.lt.rsize) then
                           
                           do kk=kklo,kkhi
                              
                              zz = zlow + dznn*(dble(kk-kklo)+half)
c     How far are we from the centre of the cell?
                              deltax = xx-xctr
                              deltay = yy-yctr
                              deltaz = zz-zctr
c     Do slope reconstruction
                              rho  = rhoctr  + deltax*dat(i,j,k,idx_rho +06) + deltay*dat(i,j,k,idx_rho +12) + deltaz*dat(i,j,k,idx_rho +18)
                              ux   = uxctr   + deltax*dat(i,j,k,idx_ux  +06) + deltay*dat(i,j,k,idx_ux  +12) + deltaz*dat(i,j,k,idx_ux  +18)
                              uy   = uyctr   + deltax*dat(i,j,k,idx_uy  +06) + deltay*dat(i,j,k,idx_uy  +12) + deltaz*dat(i,j,k,idx_uy  +18)
                              uz   = uzctr   + deltax*dat(i,j,k,idx_uz  +06) + deltay*dat(i,j,k,idx_uz  +12) + deltaz*dat(i,j,k,idx_uz  +18)
                              trac = tracctr + deltax*dat(i,j,k,idx_trac+06) + deltay*dat(i,j,k,idx_trac+12) + deltaz*dat(i,j,k,idx_trac+18)
                              temp = tempctr + deltax*dat(i,j,k,idx_temp+06) + deltay*dat(i,j,k,idx_temp+12) + deltaz*dat(i,j,k,idx_temp+18)
c     Radial and azimuthal velocities
                              ur = (ux*xx+uy*yy)/rr
                              ut = (uy*xx-ux*yy)/rr

c     Get radial derivative of radial and azimuthal velocity (more complicated than it sounds)
                              drx = mdx*xx/rr
                              dry = mdx*yy/rr
                              xxp = xx + drx
                              yyp = yy + dry
                              rrp = sqrt(xxp*xxp+yyp*yyp)
                              xxm = xx - drx
                              yym = yy - dry
                              rrm = sqrt(xxm*xxm+yym*yym)
                              uxp = uxctr   + (xxp-xctr)*dat(i,j,k,idx_ux  +06) + (yyp-yctr)*dat(i,j,k,idx_ux  +12) + deltaz*dat(i,j,k,idx_ux  +18)
                              uyp = uyctr   + (xxp-xctr)*dat(i,j,k,idx_uy  +06) + (yyp-yctr)*dat(i,j,k,idx_uy  +12) + deltaz*dat(i,j,k,idx_uy  +18)
                              uxm = uxctr   + (xxm-xctr)*dat(i,j,k,idx_ux  +06) + (yym-yctr)*dat(i,j,k,idx_ux  +12) + deltaz*dat(i,j,k,idx_ux  +18)
                              uym = uyctr   + (xxm-xctr)*dat(i,j,k,idx_uy  +06) + (yym-yctr)*dat(i,j,k,idx_uy  +12) + deltaz*dat(i,j,k,idx_uy  +18)

                              urp = (uxp*xxp+uyp*yyp)/rrp
                              utp = (uyp*xxp-uxp*yyp)/rrp
                              urm = (uxm*xxm+uym*yym)/rrm
                              utm = (uym*xxm-uxm*yym)/rrm

                              durdr = (urp-urm)/(two*mdx)
                              dutdr = (utp-utm)/(two*mdx)
c     Get azimuthal derivative of radial and azimuthal velocity (more complicated than it sounds)
                              dtx =-mdx*yy/rr
                              dty = mdx*xx/rr
                              xxp = xx + dtx
                              yyp = yy + dty
                              rrp = sqrt(xxp*xxp+yyp*yyp)
                              xxm = xx - drx
                              yym = yy - dry
                              rrm = sqrt(xxm*xxm+yym*yym)
                              uxp = uxctr   + (xxp-xctr)*dat(i,j,k,idx_ux  +06) + (yyp-yctr)*dat(i,j,k,idx_ux  +12) + deltaz*dat(i,j,k,idx_ux  +18)
                              uyp = uyctr   + (xxp-xctr)*dat(i,j,k,idx_uy  +06) + (yyp-yctr)*dat(i,j,k,idx_uy  +12) + deltaz*dat(i,j,k,idx_uy  +18)
                              uxm = uxctr   + (xxm-xctr)*dat(i,j,k,idx_ux  +06) + (yym-yctr)*dat(i,j,k,idx_ux  +12) + deltaz*dat(i,j,k,idx_ux  +18)
                              uym = uyctr   + (xxm-xctr)*dat(i,j,k,idx_uy  +06) + (yym-yctr)*dat(i,j,k,idx_uy  +12) + deltaz*dat(i,j,k,idx_uy  +18)
                              urp = (uxp*xxp+uyp*yyp)/rrp
                              utp = (uyp*xxp-uxp*yyp)/rrp
                              urm = (uxm*xxm+uym*yym)/rrm
                              utm = (uym*xxm-uxm*yym)/rrm
                              durdt = (urp-urm)/(two*mdx)
                              dutdt = (utp-utm)/(two*mdx)
c     Pressure and its derivatives
                              p      = pctr    + deltax*pres(i,j,k,idx_p   +04) + deltay*pres(i,j,k,idx_p   +08) + deltaz*pres(i,j,k,idx_p   +12)
                              dpdx   = dpdxctr + deltax*pres(i,j,k,idx_dpdx+04) + deltay*pres(i,j,k,idx_dpdx+08) + deltaz*pres(i,j,k,idx_dpdx+12)
                              dpdy   = dpdyctr + deltax*pres(i,j,k,idx_dpdy+04) + deltay*pres(i,j,k,idx_dpdy+08) + deltaz*pres(i,j,k,idx_dpdy+12)
                              dpdz   = dpdzctr + deltax*pres(i,j,k,idx_dpdz+04) + deltay*pres(i,j,k,idx_dpdz+08) + deltaz*pres(i,j,k,idx_dpdz+12)
c     Radial and azimuthal derivatives of pressure
                              dpdr   = ( xx*dpdx + yy*dpdy)/rr
                              dpdt   = (-yy*dpdx + xx*dpdy)
c     KE terms
                              KEK  = ux*ux   + uy*uy   + uz*uz
                              KEP  = ux*dpdx + uy*dpdy + uz*dpdz
                              KED  = ux*dat(i,j,k,25) + uy*dat(i,j,k,26) + uz*dat(i,j,k,27)
                              
c     Start counter
                              idx = (kk*rsize+ridx)*jetVars+00
c     0 - Let's keep track of how many cells contribute to the integral
                              jetData(idx) = jetData(idx) + one
c     1 - Denstiy
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho
c     2 - Flow vars not density-weighted
                              idx = idx + 1
                              jetData(idx) = jetData(idx) +     ur
                              idx = idx + 1
                              jetData(idx) = jetData(idx) +     ut
                              idx = idx + 1
                              jetData(idx) = jetData(idx) +     uz
                              idx = idx + 1
                              jetData(idx) = jetData(idx) +     trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) +     temp
c     7 - Flow vars density-weighted
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ut
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*uz
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*temp
c     12 - Second order (ur)
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*ur
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*ut
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*uz
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*temp
                              
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ut*ut
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ut*uz
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ut*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ut*temp

                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*uz*uz
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*uz*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*uz*temp

                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*trac*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*trac*temp

                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*temp*temp
c     27 - Third-order (ur)
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*ur*ur
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*ur*ut
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*ur*uz
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*ur*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*ur*temp
                              
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*ut*ut
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*ut*uz
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*ut*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*ut*temp

                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*uz*uz
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*uz*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*uz*temp

                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*trac*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*trac*temp

                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ur*temp*temp
c     42 - Third-order (ut)
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ut*ut*ut
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ut*ut*uz
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ut*ut*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ut*ut*temp

                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ut*uz*uz
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ut*uz*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ut*uz*temp

                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ut*trac*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ut*trac*temp

                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*ut*temp*temp
c     52 - Third-order (uz)
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*uz*uz*uz
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*uz*uz*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*uz*uz*temp

                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*uz*trac*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*uz*trac*temp

                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*uz*temp*temp
c     58 - Third-order (trac)
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*trac*trac*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*trac*trac*temp

                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*trac*temp*temp
c     61 - Third-order (temp)
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*temp*temp*temp
c     62 - Pressure
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*ur
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*ut
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*uz
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*temp
c     68 - Presure gradient (r)
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + dpdr*ur
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + dpdr*ut
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + dpdr*uz
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + dpdr*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + dpdr*temp
c     73 - Presure gradient (t)
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + dpdt*ur
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + dpdt*ut
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + dpdt*uz
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + dpdt*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + dpdt*temp
c     78 - Presure gradient (z)idx_uy
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + dpdz*ur
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + dpdz*ut
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + dpdz*uz
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + dpdz*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + dpdz*temp
c     83 - Pressure . gradient (r)
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*durdr
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*dutdr
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*((xx/rr)*dat(i,j,k,idx_uz  +06)+(yy/rr)*dat(i,j,k,idx_uz  +12))
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*((xx/rr)*dat(i,j,k,idx_trac+06)+(yy/rr)*dat(i,j,k,idx_trac+12))
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*((xx/rr)*dat(i,j,k,idx_temp+06)+(yy/rr)*dat(i,j,k,idx_temp+12))
c     88 - Pressure . gradient (t)
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*durdt
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*dutdt
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*(-yy*dat(i,j,k,idx_uz  +06)+xx*dat(i,j,k,idx_uz  +12))
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*(-yy*dat(i,j,k,idx_trac+06)+xx*dat(i,j,k,idx_trac+12))
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*(-yy*dat(i,j,k,idx_temp+06)+xx*dat(i,j,k,idx_temp+12))
c     93 - Pressure . gradient (z)
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*( xx*dat(i,j,k,idx_ux  +18)+yy*dat(i,j,k,idx_uy  +18))/rr
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*(-yy*dat(i,j,k,idx_ux  +18)+xx*dat(i,j,k,idx_uy  +18))/rr
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*dat(i,j,k,idx_uz  +18)
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*dat(i,j,k,idx_trac+18)
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*dat(i,j,k,idx_temp+18)
c     98 - Density fluctuations
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*rho
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*rho*ur
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*rho*ut
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*rho*uz
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*rho*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*rho*temp
c     103 - Pressue fluctuations
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*p
c     104 - Kinetic Energy Equation (shouldn't need these)
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*KEK
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*KEK*ur
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*KEK*ut
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*KEK*uz
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*KEP
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*KED

c     zz loop
                        enddo
c     r < rsize
                     endif
c     yy
                  enddo
c     xx
               enddo
c     zeroed out
            endif

            end do
         end do
      end do

      end
#endif

       subroutine FORT_MAXVAL(rho,DIMS(rho),DIMS(grid),mxval)

       implicit none

       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  rho(DIMV(rho))
       REAL_T  mxval

       integer i, j, k

       mxval = -Huge(0.0d0)

!$omp parallel do private(i,j,k) reduction(max : mxval)
       do k = ARG_L3(grid), ARG_H3(grid)
          do j = ARG_L2(grid), ARG_H2(grid)
             do i = ARG_L1(grid), ARG_H1(grid)
                mxval = max(mxval, rho(i,j,k))
             end do
          end do
       end do
!$omp end parallel do

       end

c :: ----------------------------------------------------------
c :: SUMMASS
c ::             MASS = sum{ vol(i,j)*rho(i,j) }
c ::
c :: INPUTS / OUTPUTS:
c ::  rho        => density field
c ::  DIMS(rho)  => index limits of rho aray
c ::  lo,hi      => index limits of grid interior
c ::  delta	 => cell size
c ::  mass      <=  total mass
c ::  r		 => radius at cell center
c ::  tmp        => temp column array
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SUMMASS(rho,DIMS(rho),DIMS(grid),delta,mass)

       implicit none

       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  mass, delta(SDIM)
       REAL_T  rho(DIMV(rho))

       integer i, j, k
       REAL_T  vol

       vol = delta(1)*delta(2)*delta(3)

       mass = zero

!$omp parallel do private(i,j,k) reduction(+ : mass)
       do k = ARG_L3(grid), ARG_H3(grid)
          do j = ARG_L2(grid), ARG_H2(grid)
             do i = ARG_L1(grid), ARG_H1(grid)
                mass = mass + rho(i,j,k)
             end do
          end do
       end do
!$omp end parallel do

       mass = vol*mass

       end

c :: ----------------------------------------------------------
c :: SUMMASSCYL
c ::    MASS = sum{ vol(i,j,k)*rho(i,j,k) } over subregion cylinder
c ::
c :: INPUTS / OUTPUTS:
c ::  rho        => density field
c ::  DIMS(rho)  => index limits of rho aray
c ::  lo,hi      => index limits of grid interior
c ::  delta	 => cell size
c ::  mass      <=  total mass
c ::  r		 => radius at cell center
c ::  tmp        => temp column array
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SUMMASS_CYL(rho,DIMS(rho),DIMS(grid),delta,
     &                             plo,vws_dz,vws_Rcyl,mass)

       implicit none

       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  mass, delta(SDIM), plo(SDIM), vws_dz, vws_Rcyl
       REAL_T  rho(DIMV(rho))

       integer i, j, k
       REAL_T  vol, x, y, z, r

       vol = delta(1)*delta(2)*delta(3)

       mass = zero

!$omp parallel do private(i,j,k,x,y,z,r) reduction(+ : mass)
       do k = ARG_L3(grid), ARG_H3(grid)
          z = plo(3) + (k+half)*delta(3)
          if (z-plo(3) .le. vws_dz) then
             do j = ARG_L2(grid), ARG_H2(grid)
                y = plo(2) + (j+half)*delta(2) 
                do i = ARG_L1(grid), ARG_H1(grid)
                   x = plo(1) + (i+half)*delta(1)
                   r = SQRT(x*x + y*y)
                   if (r .le. vws_Rcyl) then
                      mass = mass + rho(i,j,k)
                   end if
                end do
             end do
          end if
       end do
!$omp end parallel do

       mass = vol*mass

       end

c-----------------------------------------------------------------------
c     This routine fills an edge-centered fab from a cell-centered
c     fab using simple linear interpolation.
c
c     INPUTS / OUTPUTS:
c     lo,hi      => index limits of the region of the edge-centered fab
c                   to be filled
c     DIMS(cfab) => index limits of the cell-centered fab
c     cfab       => cell-centered data
c     DIMS(efab) => index limits of the edge-centered fab
c     efab       => edge-centered fab to fill
c     nc         => Number of components in the fab to fill
c     dir        => direction data needs to be shifted to get to edges
c-----------------------------------------------------------------------
c
      subroutine FORT_CEN2EDG(lo, hi, 
     &     DIMS(cfab), cfab,
     &     DIMS(efab), efab, nc, dir, isharm
     &     )
      implicit none

      integer lo(SDIM), hi(SDIM), nc, dir, isharm
      integer DIMDEC(cfab)
      integer DIMDEC(efab)
      REAL_T  cfab(DIMV(cfab), nc)
      REAL_T  efab(DIMV(efab), nc)
      integer i,j,k,n

!$omp parallel private(i,j,k)

      if ( isharm .eq. 0 ) then
         if (dir .EQ. 0) then
            do n = 1,nc
!$omp do
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1)+1, hi(1)
                        efab(i,j,k,n) = half*(cfab(i,j,k,n) + cfab(i-1,j,k,n))
                     end do
                  end do
               end do
!$omp end do nowait
            end do
         else if (dir .EQ. 1) then
            do n = 1,nc
!$omp do
               do k = lo(3), hi(3)
                  do j = lo(2)+1, hi(2)
                     do i = lo(1), hi(1)
                        efab(i,j,k,n) = half*(cfab(i,j,k,n) + cfab(i,j-1,k,n))
                     end do
                  end do
               end do
!$omp end do nowait
            end do
         else if (dir .EQ. 2) then
            do n = 1,nc
!$omp do
               do k = lo(3)+1, hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        efab(i,j,k,n) = half*(cfab(i,j,k,n) + cfab(i,j,k-1,n))
                     end do
                  end do
               end do
!$omp end do nowait
            end do
         end if
      else
         if (dir .EQ. 0) then
            do n = 1,nc
!$omp do
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1)+1, hi(1)
                        if((cfab(i,j,k,n) * cfab(i-1,j,k,n)) .gt.zero)
     &                       then
                           efab(i,j,k,n) =
     &                          2*(cfab(i,j,k,n) * cfab(i-1,j,k,n))/
     &                          (cfab(i,j,k,n) + cfab(i-1,j,k,n))
                        else
                           efab(i,j,k,n) = zero
                        endif
                     end do
                  end do
               end do
!$omp end do nowait
            end do
         else if (dir .EQ. 1) then
            do n = 1,nc
!$omp do
               do k = lo(3), hi(3)
                  do j = lo(2)+1, hi(2)
                     do i = lo(1), hi(1)
                        if((cfab(i,j,k,n) * cfab(i,j-1,k,n)).gt.zero)
     &                       then
                           efab(i,j,k,n) =
     &                          2*(cfab(i,j,k,n) * cfab(i,j-1,k,n))/
     &                          (cfab(i,j,k,n) + cfab(i,j-1,k,n))
                        else
                           efab(i,j,k,n) = zero
                        endif
                     end do
                  end do
               end do
!$omp end do nowait
            end do
         else if (dir .EQ. 2) then
            do n = 1,nc
!$omp do
               do k = lo(3)+1, hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        if((cfab(i,j,k,n) * cfab(i,j,k-1,n)).gt.zero)
     &                       then
                           efab(i,j,k,n) =
     &                          2*(cfab(i,j,k,n) * cfab(i,j,k-1,n))/
     &                          (cfab(i,j,k,n) + cfab(i,j,k-1,n))
                        else
                           efab(i,j,k,n) = zero
                        endif
                     end do
                  end do
               end do
!$omp end do nowait
            end do
         end if
      end if

!$omp end parallel

      end
c-----------------------------------------------------------------------
      subroutine EDGE_INTERP(flo, fhi, nc, ratio, dir,
     &     fine, fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2)
      implicit none
      integer flo(0:3-1), fhi(0:3-1), nc, ratio(0:3-1), dir
      integer fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2
      DOUBLE PRECISION
     &     fine(fine_l0:fine_h0,fine_l1:fine_h1,fine_l2:fine_h2,nc)
      integer i,j,k,n,P,M,L
      DOUBLE PRECISION val, df
c
c     Do linear in dir, pc transverse to dir, leave alone the fine values
c     lining up with coarse edges--assume these have been set to hold the 
c     values you want to interpolate to the rest.
c
      if (dir.eq.0) then
         do n=1,nc
            do k=flo(2),fhi(2),ratio(2)
               do j=flo(1),fhi(1),ratio(1)
                  do i=flo(0),fhi(0)-ratio(dir),ratio(0)
                     df = fine(i+ratio(dir),j,k,n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n)
     &                       + df*dble(M)/dble(ratio(dir))
                        do P=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                           do L=MAX(k,flo(2)),MIN(k+ratio(2)-1,fhi(2))
                              fine(i+M,P,L,n) = val
                           enddo
                        enddo
                     enddo                     
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do n=1,nc
            do k=flo(2),fhi(2),ratio(2)
               do j=flo(1),fhi(1)-ratio(dir),ratio(1)
                  do i=flo(0),fhi(0)
                     df = fine(i,j+ratio(dir),k,n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n)
     &                       + df*dble(M)/dble(ratio(dir))
                        do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                           do L=MAX(k,flo(2)),MIN(k+ratio(2)-1,fhi(2))
                              fine(P,j+M,L,n) = val
                           enddo
                        enddo
                     enddo                     
                  enddo
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do k=flo(2),fhi(2)-ratio(dir),ratio(2)
               do j=flo(1),fhi(1),ratio(1)
                  do i=flo(0),fhi(0),ratio(0)
                     df = fine(i,j,k+ratio(dir),n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n)
     &                       + df*dble(M)/dble(ratio(dir))
                        do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                           do L=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                              fine(P,L,k+M,n) = val
                           enddo
                        enddo
                     enddo                     
                  enddo
               enddo
            enddo
         enddo
      endif
      end
c-----------------------------------------------------------------------
      subroutine PC_EDGE_INTERP(lo, hi, nc, ratio, dir,
     &     crse, crse_l0, crse_l1, crse_l2, crse_h0, crse_h1, crse_h2,
     &     fine, fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2)
      implicit none
      integer lo(3),hi(3), nc, ratio(0:3-1), dir
      integer crse_l0, crse_l1, crse_l2, crse_h0, crse_h1, crse_h2
      integer fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2
      DOUBLE PRECISION
     &     crse(crse_l0:crse_h0,crse_l1:crse_h1,crse_l2:crse_h2,nc)
      DOUBLE PRECISION
     &     fine(fine_l0:fine_h0,fine_l1:fine_h1,fine_l2:fine_h2,nc)
      integer i,j,k,ii,jj,kk,n,L, P
c
c     For edge-based data, fill fine values with piecewise-constant interp of coarse data.
c     Operate only on faces that overlap--ie, only fill the fine faces that make up each
c     coarse face, leave the in-between faces alone.
c
      if (dir.eq.0) then
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     do P=0,ratio(2)-1
                        do L=0,ratio(1)-1
                           fine(ii,jj+L,kk+P,n) = crse(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     do P=0,ratio(2)-1
                        do L=0,ratio(0)-1
                           fine(ii+L,jj,kk+P,n) = crse(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     do P=0,ratio(1)-1
                        do L=0,ratio(0)-1
                           fine(ii+L,jj+P,kk,n) = crse(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
      end



      subroutine FORT_HERSCHEL_BULKLEY(visc, DIMS(visc), vel, DIMS(vel),
     &                                 lo, hi, domlo, domhi, delta, vel_bc)
      use viscoplasticity_module
      implicit none

c
c ::: This routine will compute the Papanastasiou regularised Bingham
c ::: viscosity based on the velocity field
c
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(visc)
      integer  DIMDEC(vel)
      integer  domlo(SDIM), domhi(SDIM)
      integer  vel_bc(SDIM,2,SDIM)
      REAL_T   delta(SDIM)
      REAL_T   visc(DIMV(visc),1)
      REAL_T   vel(DIMV(vel),SDIM)

      integer i,j,k
      REAL_T  ux, uy, uz, vx, vy, vz, wx, wy, wz, dx, dy, dz
      REAL_T  uxcen, uycen, uzcen, uxlo, uylo, uzlo, uxhi, uyhi, uzhi
      REAL_T  vxcen, vycen, vzcen, vxlo, vylo, vzlo, vxhi, vyhi, vzhi
      REAL_T  wxcen, wycen, wzcen, wxlo, wylo, wzlo, wxhi, wyhi, wzhi
      REAL_T  strnrt

      logical fixulo_x, fixvlo_x, fixwlo_x, fixuhi_x, fixvhi_x, fixwhi_x
      logical fixulo_y, fixvlo_y, fixwlo_y, fixuhi_y, fixvhi_y, fixwhi_y
      logical fixulo_z, fixvlo_z, fixwlo_z, fixuhi_z, fixvhi_z, fixwhi_z
c
c     ::::: some useful macro definitions
c
#     define U(i,j,k) vel(i,j,k,1)
#     define V(i,j,k) vel(i,j,k,2)
#     define W(i,j,k) vel(i,j,k,3)

#     define ULOX vel_bc(1,1,1)
#     define UHIX vel_bc(1,2,1)
#     define ULOY vel_bc(2,1,1)
#     define UHIY vel_bc(2,2,1)
#     define ULOZ vel_bc(3,1,1)
#     define UHIZ vel_bc(3,2,1)

#     define VLOX vel_bc(1,1,2)
#     define VHIX vel_bc(1,2,2)
#     define VLOY vel_bc(2,1,2)
#     define VHIY vel_bc(2,2,2)
#     define VLOZ vel_bc(3,1,2)
#     define VHIZ vel_bc(3,2,2)

#     define WLOX vel_bc(1,1,3)
#     define WHIX vel_bc(1,2,3)
#     define WLOY vel_bc(2,1,3)
#     define WHIY vel_bc(2,2,3)
#     define WLOZ vel_bc(3,1,3)
#     define WHIZ vel_bc(3,2,3)
c
c     ::::: statement functions that implement stencil
c
      uxcen(i,j,k) = half*(U(i+1,j,k)-U(i-1,j,k))/dx
      uxlo(i,j,k)  = (U(i+1,j,k)+three*U(i,j,k)-four*U(i-1,j,k))/(three*dx)
      uxhi(i,j,k)  =-(U(i-1,j,k)+three*U(i,j,k)-four*U(i+1,j,k))/(three*dx)

      uycen(i,j,k) = half*(U(i,j+1,k)-U(i,j-1,k))/dy
      uylo(i,j,k)  = (U(i,j+1,k)+three*U(i,j,k)-four*U(i,j-1,k))/(three*dy)
      uyhi(i,j,k)  =-(U(i,j-1,k)+three*U(i,j,k)-four*U(i,j+1,k))/(three*dy)

      uzcen(i,j,k) = half*(U(i,j,k+1)-U(i,j,k-1))/dz
      uzlo(i,j,k)  = (U(i,j,k+1)+three*U(i,j,k)-four*U(i,j,k-1))/(three*dz)
      uzhi(i,j,k)  =-(U(i,j,k-1)+three*U(i,j,k)-four*U(i,j,k+1))/(three*dz)

      vxcen(i,j,k) = half*(V(i+1,j,k)-V(i-1,j,k))/dx
      vxlo(i,j,k)  = (V(i+1,j,k)+three*V(i,j,k)-four*V(i-1,j,k))/(three*dx)
      vxhi(i,j,k)  =-(V(i-1,j,k)+three*V(i,j,k)-four*V(i+1,j,k))/(three*dx)

      vycen(i,j,k) = half*(V(i,j+1,k)-V(i,j-1,k))/dy
      vylo(i,j,k)  = (V(i,j+1,k)+three*V(i,j,k)-four*V(i,j-1,k))/(three*dy)
      vyhi(i,j,k)  =-(V(i,j-1,k)+three*V(i,j,k)-four*V(i,j+1,k))/(three*dy)

      vzcen(i,j,k) = half*(V(i,j,k+1)-V(i,j,k-1))/dz
      vzlo(i,j,k)  = (V(i,j,k+1)+three*V(i,j,k)-four*V(i,j,k-1))/(three*dz)
      vzhi(i,j,k)  =-(V(i,j,k-1)+three*V(i,j,k)-four*V(i,j,k+1))/(three*dz)

      wxcen(i,j,k) = half*(W(i+1,j,k)-W(i-1,j,k))/dx
      wxlo(i,j,k)  = (W(i+1,j,k)+three*W(i,j,k)-four*W(i-1,j,k))/(three*dx)
      wxhi(i,j,k)  =-(W(i-1,j,k)+three*W(i,j,k)-four*W(i+1,j,k))/(three*dx)

      wycen(i,j,k) = half*(W(i,j+1,k)-W(i,j-1,k))/dy
      wylo(i,j,k)  = (W(i,j+1,k)+three*W(i,j,k)-four*W(i,j-1,k))/(three*dy)
      wyhi(i,j,k)  =-(W(i,j-1,k)+three*W(i,j,k)-four*W(i,j+1,k))/(three*dy)

      wzcen(i,j,k) = half*(W(i,j,k+1)-W(i,j,k-1))/dz
      wzlo(i,j,k)  = (W(i,j,k+1)+three*W(i,j,k)-four*W(i,j,k-1))/(three*dz)
      wzhi(i,j,k)  =-(W(i,j,k-1)+three*W(i,j,k)-four*W(i,j,k+1))/(three*dz)


      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ux = uxcen(i,j,k)
               uy = uycen(i,j,k)
               uz = uzcen(i,j,k)
               vx = vxcen(i,j,k)
               vy = vycen(i,j,k)
               vz = vzcen(i,j,k)
               wx = wxcen(i,j,k)
               wy = wycen(i,j,k)
               wz = wzcen(i,j,k)
               strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
               visc(i,j,k,1) = visc_fun(strnrt)
            end do
         end do
      end do

      fixulo_x = ( (lo(1) .eq. domlo(1)) .and.
     &             (ULOX .eq. EXT_DIR .or. ULOX .eq. HOEXTRAP) )
      fixuhi_x = ( (hi(1) .eq. domhi(1)) .and.
     &             (UHIX .eq. EXT_DIR .or. UHIX .eq. HOEXTRAP) )
      fixvlo_x = ( (lo(1) .eq. domlo(1)) .and.
     &             (VLOX .eq. EXT_DIR .or. VLOX .eq. HOEXTRAP) )
      fixvhi_x = ( (hi(1) .eq. domhi(1)) .and.
     &             (VHIX .eq. EXT_DIR .or. VHIX .eq. HOEXTRAP) )
      fixwlo_x = ( (lo(1) .eq. domlo(1)) .and.
     &             (WLOX .eq. EXT_DIR .or. WLOX .eq. HOEXTRAP) )
      fixwhi_x = ( (hi(1) .eq. domhi(1)) .and.
     &             (WHIX .eq. EXT_DIR .or. WHIX .eq. HOEXTRAP) )

      fixulo_y = ( (lo(2) .eq. domlo(2)) .and.
     &             (ULOY .eq. EXT_DIR .or. ULOY .eq. HOEXTRAP) )
      fixuhi_y = ( (hi(2) .eq. domhi(2)) .and.
     &             (UHIY .eq. EXT_DIR .or. UHIY .eq. HOEXTRAP) )
      fixvlo_y = ( (lo(2) .eq. domlo(2)) .and.
     &             (VLOY .eq. EXT_DIR .or. VLOY .eq. HOEXTRAP) )
      fixvhi_y = ( (hi(2) .eq. domhi(2)) .and.
     &             (VHIY .eq. EXT_DIR .or. VHIY .eq. HOEXTRAP) )
      fixwlo_y = ( (lo(2) .eq. domlo(2)) .and.
     &             (WLOY .eq. EXT_DIR .or. WLOY .eq. HOEXTRAP) )
      fixwhi_y = ( (hi(2) .eq. domhi(2)) .and.
     &             (WHIY .eq. EXT_DIR .or. WHIY .eq. HOEXTRAP) )

      fixulo_z = ( (lo(3) .eq. domlo(3)) .and.
     &             (ULOZ .eq. EXT_DIR .or. ULOZ .eq. HOEXTRAP) )
      fixuhi_z = ( (hi(3) .eq. domhi(3)) .and.
     &             (UHIZ .eq. EXT_DIR .or. UHIZ .eq. HOEXTRAP) )
      fixvlo_z = ( (lo(3) .eq. domlo(3)) .and.
     &             (VLOZ .eq. EXT_DIR .or. VLOZ .eq. HOEXTRAP) )
      fixvhi_z = ( (hi(3) .eq. domhi(3)) .and.
     &             (VHIZ .eq. EXT_DIR .or. VHIZ .eq. HOEXTRAP) )
      fixwlo_z = ( (lo(3) .eq. domlo(3)) .and.
     &             (WLOZ .eq. EXT_DIR .or. WLOZ .eq. HOEXTRAP) )
      fixwhi_z = ( (hi(3) .eq. domhi(3)) .and.
     &             (WHIZ .eq. EXT_DIR .or. WHIZ .eq. HOEXTRAP) )
c
c     First do all the faces
c
      if (fixulo_x .or. fixvlo_x .or. fixwlo_x) then
         i = lo(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
               vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
               wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
               uy = uycen(i,j,k)
               vy = vycen(i,j,k)
               wy = wycen(i,j,k)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wz = wzcen(i,j,k)
               strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
               visc(i,j,k,1) = visc_fun(strnrt)
               ! Linear extrapolation 
               visc(i-1,j,k,1) = visc(i,j,k,1)
            end do
         end do
      end if

      if (fixuhi_x .or. fixvhi_x .or. fixwhi_x) then
         i = hi(1)
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
               vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
               wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
               uy = uycen(i,j,k)
               vy = vycen(i,j,k)
               wy = wycen(i,j,k)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wz = wzcen(i,j,k)
               strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
               visc(i,j,k,1) = visc_fun(strnrt)
               ! Linear extrapolation 
               visc(i+1,j,k,1) = visc(i,j,k,1)
            end do
         end do
      end if

      if (fixulo_y .or. fixvlo_y .or. fixwlo_y) then
         j = lo(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               ux = uxcen(i,j,k)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
               vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
               wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wz = wzcen(i,j,k)
               strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
               visc(i,j,k,1) = visc_fun(strnrt)
               ! Linear extrapolation 
               visc(i,j-1,k,1) = visc(i,j,k,1)
            end do
         end do
      end if

      if (fixuhi_y .or. fixvhi_y .or. fixwhi_y) then
         j = hi(2)
         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               ux = uxcen(i,j,k)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
               vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
               wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
               uz = uzcen(i,j,k)
               vz = vzcen(i,j,k)
               wz = wzcen(i,j,k)
               strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
               visc(i,j,k,1) = visc_fun(strnrt)
               ! Linear extrapolation 
               visc(i,j+1,k,1) = visc(i,j,k,1)
            end do
         end do
      end if

      if (fixulo_z .or. fixvlo_z .or. fixwlo_z) then
         k = lo(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               ux = uxcen(i,j,k)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = uycen(i,j,k)
               vy = vycen(i,j,k)
               wy = wycen(i,j,k)
               uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
               vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
               wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
               strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
               visc(i,j,k,1) = visc_fun(strnrt)
               ! Linear extrapolation 
               visc(i,j,k-1,1) = visc(i,j,k,1)
            end do
         end do
      end if

      if (fixuhi_z .or. fixvhi_z .or. fixwhi_z) then
         k = hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               ux = uxcen(i,j,k)
               vx = vxcen(i,j,k)
               wx = wxcen(i,j,k)
               uy = uycen(i,j,k)
               vy = vycen(i,j,k)
               wy = wycen(i,j,k)
               uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
               vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
               wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
               strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
               visc(i,j,k,1) = visc_fun(strnrt)
               ! Linear extrapolation 
               visc(i,j,k+1,1) = visc(i,j,k,1)
            end do
         end do
      end if
c
c     Next do all the edges
c
      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. 
     &    (fixulo_y .or. fixvlo_y .or. fixwlo_y)) then
         i = lo(1)
         j = lo(2)
         do k = lo(3),hi(3)
            ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            wz = wzcen(i,j,k)
            strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
            visc(i,j,k,1) = visc_fun(strnrt)
            ! Linear extrapolation 
            visc(i-1,j-1,k,1) = visc(i,j,k,1)
         end do
      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. 
     &    (fixulo_y .or. fixvlo_y .or. fixwlo_y)) then
         i = hi(1)
         j = lo(2)
         do k = lo(3),hi(3)
            ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            wz = wzcen(i,j,k)
            strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
            visc(i,j,k,1) = visc_fun(strnrt)
            ! Linear extrapolation 
            visc(i+1,j-1,k,1) = visc(i,j,k,1)
         end do
      end if

      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. 
     &    (fixuhi_y .or. fixvhi_y .or. fixwhi_y)) then
         i = lo(1)
         j = hi(2)
         do k = lo(3),hi(3)
            ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            wz = wzcen(i,j,k)
            strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
            visc(i,j,k,1) = visc_fun(strnrt)
            ! Linear extrapolation 
            visc(i-1,j+1,k,1) = visc(i,j,k,1)
         end do
      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. 
     &    (fixuhi_y .or. fixvhi_y .or. fixwhi_y)) then
         i = hi(1)
         j = hi(2)
         do k = lo(3),hi(3)
            ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = uzcen(i,j,k)
            vz = vzcen(i,j,k)
            wz = wzcen(i,j,k)
            strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
            visc(i,j,k,1) = visc_fun(strnrt)
            ! Linear extrapolation 
            visc(i+1,j+1,k,1) = visc(i,j,k,1)
         end do
      end if

      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. 
     &    (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = lo(1)
         k = lo(3)
         do j = lo(2),hi(2)
            ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uy = uycen(i,j,k)
            vy = vycen(i,j,k)
            wy = wycen(i,j,k)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
            strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
            visc(i,j,k,1) = visc_fun(strnrt)
            ! Linear extrapolation 
            visc(i-1,j,k-1,1) = visc(i,j,k,1)
         end do
      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. 
     &    (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = hi(1)
         k = lo(3)
         do j = lo(2),hi(2)
            ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uy = uycen(i,j,k)
            vy = vycen(i,j,k)
            wy = wycen(i,j,k)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
            strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
            visc(i,j,k,1) = visc_fun(strnrt)
            ! Linear extrapolation 
            visc(i+1,j,k-1,1) = visc(i,j,k,1)
         end do
      end if

      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. 
     &    (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         i = lo(1)
         k = hi(3)
         do j = lo(2),hi(2)
            ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
            vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
            wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
            uy = uycen(i,j,k)
            vy = vycen(i,j,k)
            wy = wycen(i,j,k)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
            strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
            visc(i,j,k,1) = visc_fun(strnrt)
            ! Linear extrapolation 
            visc(i-1,j,k+1,1) = visc(i,j,k,1)
         end do
      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. 
     &    (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         i = hi(1)
         k = hi(3)
         do j = lo(2),hi(2)
            ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
            vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
            wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
            uy = uycen(i,j,k)
            vy = vycen(i,j,k)
            wy = wycen(i,j,k)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
            strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
            visc(i,j,k,1) = visc_fun(strnrt)
            ! Linear extrapolation 
            visc(i+1,j,k+1,1) = visc(i,j,k,1)
         end do
      end if

      if ((fixulo_y .or. fixvlo_y .or. fixwlo_y) .and. 
     &    (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         j = lo(2)
         k = lo(3)
         do i = lo(1),hi(1)
            ux = uxcen(i,j,k)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
            strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
            visc(i,j,k,1) = visc_fun(strnrt)
            ! Linear extrapolation 
            visc(i,j-1,k-1,1) = visc(i,j,k,1)
         end do
      end if

      if ((fixuhi_y .or. fixvhi_y .or. fixwhi_y) .and. 
     &    (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         j = hi(2)
         k = lo(3)
         do i = lo(1),hi(1)
            ux = uxcen(i,j,k)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
            vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
            wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
            strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
            visc(i,j,k,1) = visc_fun(strnrt)
            ! Linear extrapolation 
            visc(i,j+1,k-1,1) = visc(i,j,k,1)
         end do
      end if

      if ((fixulo_y .or. fixvlo_y .or. fixwlo_y) .and. 
     &    (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         j = lo(2)
         k = hi(3)
         do i = lo(1),hi(1)
            ux = uxcen(i,j,k)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
            vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
            wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
            strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
            visc(i,j,k,1) = visc_fun(strnrt)
            ! Linear extrapolation 
            visc(i,j-1,k+1,1) = visc(i,j,k,1)
         end do
      end if

      if ((fixuhi_y .or. fixvhi_y .or. fixwhi_y) .and. 
     &    (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         j = hi(2)
         k = hi(3)
         do i = lo(1),hi(1)
            ux = vxcen(i,j,k)
            vx = vxcen(i,j,k)
            wx = wxcen(i,j,k)
            uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
            vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
            wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
            uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
            vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
            wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
            strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
            visc(i,j,k,1) = visc_fun(strnrt)
            ! Linear extrapolation 
            visc(i,j+1,k+1,1) = visc(i,j,k,1)
         end do
      end if
c
c     Finally do all the corners
c
      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. 
     &    (fixulo_y .or. fixvlo_y .or. fixwlo_y) .and. 
     &    (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = lo(1)
         j = lo(2)
         k = lo(3)
         ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
         strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
         visc(i,j,k,1) = visc_fun(strnrt)
         ! Linear extrapolation 
         visc(i-1,j-1,k-1,1) = visc(i,j,k,1)
      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. 
     &    (fixulo_y .or. fixvlo_y .or. fixwlo_y) .and. 
     &    (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = hi(1)
         j = lo(2)
         k = lo(3)
         ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
         strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
         visc(i,j,k,1) = visc_fun(strnrt)
         ! Linear extrapolation 
         visc(i+1,j-1,k-1,1) = visc(i,j,k,1)
      end if

      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. 
     &    (fixuhi_y .or. fixvhi_y .or. fixwhi_y) .and. 
     &    (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = lo(1)
         j = hi(2)
         k = lo(3)
         ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
         strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
         visc(i,j,k,1) = visc_fun(strnrt)
         ! Linear extrapolation 
         visc(i-1,j+1,k-1,1) = visc(i,j,k,1)
      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. 
     &    (fixuhi_y .or. fixvhi_y .or. fixwhi_y) .and. 
     &    (fixulo_z .or. fixvlo_z .or. fixwlo_z)) then
         i = hi(1)
         j = hi(2)
         k = lo(3)
         ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzlo(i,j,k),uzcen(i,j,k),fixulo_z)
         vz = merge(vzlo(i,j,k),vzcen(i,j,k),fixvlo_z)
         wz = merge(wzlo(i,j,k),wzcen(i,j,k),fixwlo_z)
         strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
         visc(i,j,k,1) = visc_fun(strnrt)
         ! Linear extrapolation 
         visc(i+1,j+1,k-1,1) = visc(i,j,k,1)
      end if

      if ((fixulo_x .or. fixvlo_x .or. fixwlo_x) .and. 
     &    (fixulo_y .or. fixvlo_y .or. fixwlo_y) .and. 
     &    (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         i = lo(1)
         j = lo(2)
         k = hi(3)
         ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
         strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
         visc(i,j,k,1) = visc_fun(strnrt)
         ! Linear extrapolation 
         visc(i-1,j-1,k+1,1) = visc(i,j,k,1)
      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. 
     &    (fixulo_y .or. fixvlo_y .or. fixwlo_y) .and. 
     &    (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         i = hi(1)
         j = lo(2)
         k = hi(3)
         ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uy = merge(uylo(i,j,k),uycen(i,j,k),fixulo_y)
         vy = merge(vylo(i,j,k),vycen(i,j,k),fixvlo_y)
         wy = merge(wylo(i,j,k),wycen(i,j,k),fixwlo_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
         strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
         visc(i,j,k,1) = visc_fun(strnrt)
         ! Linear extrapolation 
         visc(i+1,j-1,k+1,1) = visc(i,j,k,1)
      end if

      if ((fixvlo_x .or. fixwlo_x) .and. 
     &    (fixuhi_y .or. fixwhi_y) .and. 
     &    (fixuhi_z .or. fixvhi_z)) then
         i = lo(1)
         j = hi(2)
         k = hi(3)
         ux = merge(uxlo(i,j,k),uxcen(i,j,k),fixulo_x)
         vx = merge(vxlo(i,j,k),vxcen(i,j,k),fixvlo_x)
         wx = merge(wxlo(i,j,k),wxcen(i,j,k),fixwlo_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
         strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
         visc(i,j,k,1) = visc_fun(strnrt)
         ! Linear extrapolation 
         visc(i-1,j+1,k+1,1) = visc(i,j,k,1)
      end if

      if ((fixuhi_x .or. fixvhi_x .or. fixwhi_x) .and. 
     &    (fixuhi_y .or. fixvhi_y .or. fixwhi_y) .and. 
     &    (fixuhi_z .or. fixvhi_z .or. fixwhi_z)) then
         i = hi(1)
         j = hi(2)
         k = hi(3)
         ux = merge(uxhi(i,j,k),uxcen(i,j,k),fixuhi_x)
         vx = merge(vxhi(i,j,k),vxcen(i,j,k),fixvhi_x)
         wx = merge(wxhi(i,j,k),wxcen(i,j,k),fixwhi_x)
         uy = merge(uyhi(i,j,k),uycen(i,j,k),fixuhi_y)
         vy = merge(vyhi(i,j,k),vycen(i,j,k),fixvhi_y)
         wy = merge(wyhi(i,j,k),wycen(i,j,k),fixwhi_y)
         uz = merge(uzhi(i,j,k),uzcen(i,j,k),fixuhi_z)
         vz = merge(vzhi(i,j,k),vzcen(i,j,k),fixvhi_z)
         wz = merge(wzhi(i,j,k),wzcen(i,j,k),fixwhi_z)
         strnrt = strnrt_fun_3d(ux,uy,uz,vx,vy,vz,wx,wy,wz)
         visc(i,j,k,1) = visc_fun(strnrt)
         ! Linear extrapolation 
         visc(i+1,j+1,k+1,1) = visc(i,j,k,1)
      end if

#     undef U
#     undef V      
#     undef W
#     undef ULOX
#     undef UHIX
#     undef ULOY
#     undef UHIY
#     undef ULOZ
#     undef UHIZ
#     undef VLOX
#     undef VHIX
#     undef VLOY
#     undef VHIY
#     undef VLOZ
#     undef VHIZ
#     undef WLOX
#     undef WHIX
#     undef WLOY
#     undef WHIY
#     undef WLOZ
#     undef WHIZ

      end
