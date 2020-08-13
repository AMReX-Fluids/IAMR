
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

module navierstokes_3d_module
  
  implicit none

  private 

  public :: sumturb, &
#ifdef SUMJET
            sum_jet, &
#endif
            fort_maxval, cen2edg, FORT_AVERAGE_EDGE_STATES
  
contains

!c ::
!c :: ----------------------------------------------------------
!c :: SUMTURB
!c :: ----------------------------------------------------------
!c ::
      subroutine sumturb(dat,pres,DIMS(dat),DIMS(pres),DIMS(grid),delta, &
                             turb,ksize,turbVars) bind(C, name="sumturb")

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

!c     Here are the derivatives, can't do it here because of zeroed intersections
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
      
      end subroutine sumturb

#ifdef SUMJET
!c ::
!c :: ----------------------------------------------------------
!c :: SUMJET
!c :: ----------------------------------------------------------
!c ::
      subroutine sum_jet(dat,pres,DIMS(dat),DIMS(pres),DIMS(grid),delta, &
                           jetData,levRsize,levKsize,rsize,ksize,jetVars,numSplit, &
                           xlo,xhi) bind(C, name="sum_jet")

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

!c     Hack to zero
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
!c     Note this is right, it really should be kn not nn
      dznn = dz/dble(kn)
!c     Mod dx
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
!c     Bounds for vertical numerical integration
         kklo=k*kn
         kkhi=kklo+kn-1
!c     Calculate left hand edge of cell
         zlow = xlo(3) + dz*dble(k-klo)
!c     And x centre/er
         zctr = zlow + 0.5d0*dz
         
         do j = jlo, jhi
!c     Calculate front of cell
            ylow = xlo(2) + dy*dble(j-jlo) - jet_y
!c     And y centre/er
            yctr = ylow + 0.5d0*dy
            
            do i = ilo, ihi
!c     Calculate left hand edge of cell
               xlow = xlo(1) + dx*dble(i-ilo) - jet_x
!c     And x centre/er
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
!c     We're not in a zeroed out bit, let's integrate...

!c     Integrate numerically
                  do jj=1,nn
                     yy = ylow + dynn*(dble(jj)-half)
                     
                     do ii=1,nn
                        xx = xlow + dxnn*(dble(ii)-half)
                        
                        rr = sqrt(xx*xx+yy*yy)
                        
                        ridx = int(rr/gridDx)
                        
                        if (ridx.lt.rsize) then
                           
                           do kk=kklo,kkhi
                              
                              zz = zlow + dznn*(dble(kk-kklo)+half)
!c     How far are we from the centre of the cell?
                              deltax = xx-xctr
                              deltay = yy-yctr
                              deltaz = zz-zctr
!c     Do slope reconstruction
                              rho  = rhoctr  + deltax*dat(i,j,k,idx_rho +06) + deltay*dat(i,j,k,idx_rho +12) + deltaz*dat(i,j,k,idx_rho +18)
                              ux   = uxctr   + deltax*dat(i,j,k,idx_ux  +06) + deltay*dat(i,j,k,idx_ux  +12) + deltaz*dat(i,j,k,idx_ux  +18)
                              uy   = uyctr   + deltax*dat(i,j,k,idx_uy  +06) + deltay*dat(i,j,k,idx_uy  +12) + deltaz*dat(i,j,k,idx_uy  +18)
                              uz   = uzctr   + deltax*dat(i,j,k,idx_uz  +06) + deltay*dat(i,j,k,idx_uz  +12) + deltaz*dat(i,j,k,idx_uz  +18)
                              trac = tracctr + deltax*dat(i,j,k,idx_trac+06) + deltay*dat(i,j,k,idx_trac+12) + deltaz*dat(i,j,k,idx_trac+18)
                              temp = tempctr + deltax*dat(i,j,k,idx_temp+06) + deltay*dat(i,j,k,idx_temp+12) + deltaz*dat(i,j,k,idx_temp+18)
!c     Radial and azimuthal velocities
                              ur = (ux*xx+uy*yy)/rr
                              ut = (uy*xx-ux*yy)/rr

!c     Get radial derivative of radial and azimuthal velocity (more complicated than it sounds)
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
!c     Get azimuthal derivative of radial and azimuthal velocity (more complicated than it sounds)
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
!c     Pressure and its derivatives
                              p      = pctr    + deltax*pres(i,j,k,idx_p   +04) + deltay*pres(i,j,k,idx_p   +08) + deltaz*pres(i,j,k,idx_p   +12)
                              dpdx   = dpdxctr + deltax*pres(i,j,k,idx_dpdx+04) + deltay*pres(i,j,k,idx_dpdx+08) + deltaz*pres(i,j,k,idx_dpdx+12)
                              dpdy   = dpdyctr + deltax*pres(i,j,k,idx_dpdy+04) + deltay*pres(i,j,k,idx_dpdy+08) + deltaz*pres(i,j,k,idx_dpdy+12)
                              dpdz   = dpdzctr + deltax*pres(i,j,k,idx_dpdz+04) + deltay*pres(i,j,k,idx_dpdz+08) + deltaz*pres(i,j,k,idx_dpdz+12)
!c     Radial and azimuthal derivatives of pressure
                              dpdr   = ( xx*dpdx + yy*dpdy)/rr
                              dpdt   = (-yy*dpdx + xx*dpdy)
!c     KE terms
                              KEK  = ux*ux   + uy*uy   + uz*uz
                              KEP  = ux*dpdx + uy*dpdy + uz*dpdz
                              KED  = ux*dat(i,j,k,25) + uy*dat(i,j,k,26) + uz*dat(i,j,k,27)
                              
!c     Start counter
                              idx = (kk*rsize+ridx)*jetVars+00
!c     0 - Let's keep track of how many cells contribute to the integral
                              jetData(idx) = jetData(idx) + one
!c     1 - Denstiy
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho
!c     2 - Flow vars not density-weighted
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
!c     7 - Flow vars density-weighted
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
!c     12 - Second order (ur)
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
!c     27 - Third-order (ur)
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
!c     42 - Third-order (ut)
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
!c     52 - Third-order (uz)
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
!c     58 - Third-order (trac)
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*trac*trac*trac
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*trac*trac*temp

                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*trac*temp*temp
!c     61 - Third-order (temp)
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + rho*temp*temp*temp
!c     62 - Pressure
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
!c     68 - Presure gradient (r)
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
!c     73 - Presure gradient (t)
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
!c     78 - Presure gradient (z)idx_uy
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
!c     83 - Pressure . gradient (r)
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
!c     88 - Pressure . gradient (t)
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
!c     93 - Pressure . gradient (z)
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
!c     98 - Density fluctuations
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
!c     103 - Pressue fluctuations
                              idx = idx + 1
                              jetData(idx) = jetData(idx) + p*p
!c     104 - Kinetic Energy Equation (shouldn't need these)
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

!c     zz loop
                        enddo
!c     r < rsize
                     endif
!c     yy
                  enddo
!c     xx
               enddo
!c     zeroed out
            endif

            end do
         end do
      end do

      end subroutine sum_jet
#endif

       subroutine fort_maxval(rho,DIMS(rho),DIMS(grid),mxval)bind(C, name="fort_maxval")

       implicit none

       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  rho(DIMV(rho))
       REAL_T  mxval

       integer i, j, k

       mxval = -Huge(0.0d0)

       do k = ARG_L3(grid), ARG_H3(grid)
          do j = ARG_L2(grid), ARG_H2(grid)
             do i = ARG_L1(grid), ARG_H1(grid)
                mxval = max(mxval, rho(i,j,k))
             end do
          end do
       end do

       end subroutine fort_maxval

!c-----------------------------------------------------------------------
!c     This routine fills an edge-centered fab from a cell-centered
!c     fab using simple linear interpolation.
!c
!c     INPUTS / OUTPUTS:
!c     lo,hi      => index limits of the region of the edge-centered fab
!c                   to be filled
!c     DIMS(cfab) => index limits of the cell-centered fab
!c     cfab       => cell-centered data
!c     DIMS(efab) => index limits of the edge-centered fab
!c     efab       => edge-centered fab to fill
!c     nc         => Number of components in the fab to fill
!c     dir        => direction data needs to be shifted to get to edges
!c-----------------------------------------------------------------------
!c
      subroutine cen2edg(lo, hi, &
          DIMS(cfab), cfab, &
          DIMS(efab), efab, nc, dir, isharm &
          ) bind(C,name="cen2edg")
          
      implicit none

      integer lo(SDIM), hi(SDIM), nc, dir, isharm
      integer DIMDEC(cfab)
      integer DIMDEC(efab)
      REAL_T  cfab(DIMV(cfab), nc)
      REAL_T  efab(DIMV(efab), nc)
      integer i,j,k,n

      if ( isharm .eq. 0 ) then
         if (dir .EQ. 0) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        efab(i,j,k,n) = half*(cfab(i,j,k,n) + cfab(i-1,j,k,n))
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 1) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        efab(i,j,k,n) = half*(cfab(i,j,k,n) + cfab(i,j-1,k,n))
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 2) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        efab(i,j,k,n) = half*(cfab(i,j,k,n) + cfab(i,j,k-1,n))
                     end do
                  end do
               end do
            end do
         end if
      else
         if (dir .EQ. 0) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        if((cfab(i,j,k,n) * cfab(i-1,j,k,n)) .gt.zero) &
                            then
                           efab(i,j,k,n) = &
                               2*(cfab(i,j,k,n) * cfab(i-1,j,k,n))/ &
                               (cfab(i,j,k,n) + cfab(i-1,j,k,n))
                        else
                           efab(i,j,k,n) = zero
                        endif
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 1) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        if((cfab(i,j,k,n) * cfab(i,j-1,k,n)).gt.zero) &
                            then
                           efab(i,j,k,n) = &
                               2*(cfab(i,j,k,n) * cfab(i,j-1,k,n))/ &
                               (cfab(i,j,k,n) + cfab(i,j-1,k,n))
                        else
                           efab(i,j,k,n) = zero
                        endif
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 2) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        if((cfab(i,j,k,n) * cfab(i,j,k-1,n)).gt.zero) &
                            then
                           efab(i,j,k,n) = &
                               2*(cfab(i,j,k,n) * cfab(i,j,k-1,n))/ &
                               (cfab(i,j,k,n) + cfab(i,j,k-1,n))
                        else
                           efab(i,j,k,n) = zero
                        endif
                     end do
                  end do
               end do
            end do
         end if
      end if

      end subroutine cen2edg

!c     
!c     
!c     ::: -----------------------------------------------------------
!c     
!c     This routine averages the mac face velocities for makeforce at half time

   subroutine FORT_AVERAGE_EDGE_STATES( vel, v_lo, v_hi,&
                                        umacx, ux_lo, ux_hi,&
                                        umacy, uy_lo, uy_hi,&
#if ( AMREX_SPACEDIM == 3 )
                                        umacz, uz_lo, uz_hi,&
#endif
                                        getForceVerbose)&
                                        bind(C, name="FORT_AVERAGE_EDGE_STATES")

      implicit none

      integer :: v_lo(3), v_hi(3)
      integer :: ux_lo(3), ux_hi(3)
      integer :: uy_lo(3), uy_hi(3)
#if ( AMREX_SPACEDIM == 3 )
      integer :: uz_lo(3), uz_hi(3)
#endif
      integer :: getForceVerbose
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3), SDIM) :: vel
      REAL_T, dimension(ux_lo(1):ux_hi(1),ux_lo(2):ux_hi(2),ux_lo(3):ux_hi(3)) :: umacx
      REAL_T, dimension(uy_lo(1):uy_hi(1),uy_lo(2):uy_hi(2),uy_lo(3):uy_hi(3)) :: umacy
#if ( AMREX_SPACEDIM == 3 )
      REAL_T, dimension(uz_lo(1):uz_hi(1),uz_lo(2):uz_hi(2),uz_lo(3):uz_hi(3)) :: umacz
#endif

      REAL_T  :: velmin(3)
      REAL_T  :: velmax(3)
      integer :: isioproc

      integer :: i, j, k, n

      do n = 1, SDIM
         velmin(n) = 1.d234
         velmax(n) = -1.d234
      enddo

      do k = v_lo(3), v_hi(3)
         do j = v_lo(2), v_hi(2)
            do i = v_lo(1), v_hi(1)
               vel(i,j,k,1) = half*(umacx(i,j,k)+umacx(i+1,j,k))
               vel(i,j,k,2) = half*(umacy(i,j,k)+umacy(i,j+1,k))
#if ( AMREX_SPACEDIM == 3 )
               vel(i,j,k,3) = half*(umacz(i,j,k)+umacz(i,j,k+1))
#endif
               do n = 1, SDIM
                  velmin(n) = min(velmin(n),vel(i,j,k,n))
                  velmax(n) = max(velmax(n),vel(i,j,k,n))
               enddo
            enddo
         enddo
      enddo

      if (getForceVerbose.gt.0) then
         call bl_pd_is_ioproc(isioproc)
         if (isioproc.eq.1) then
            do n = 1, SDIM
               write (6,*) "mac velmin (",n,") = ",velmin(n)
               write (6,*) "mac velmax (",n,") = ",velmax(n)
            enddo
         endif
      endif

   end subroutine FORT_AVERAGE_EDGE_STATES

end module navierstokes_3d_module
