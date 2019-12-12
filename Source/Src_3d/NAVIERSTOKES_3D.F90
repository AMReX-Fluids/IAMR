
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

  public :: gradp, fort_putdown, sumturb, &
#ifdef SUMJET
            sum_jet, &
#endif
            fort_maxval, summass, summass_eb, summass_cyl, cen2edg, edge_interp, &
            pc_edge_interp, filcc_tile
  
contains

      subroutine gradp ( &
          p,DIMS(p), &
          gp,DIMS(gp), &
          lo,hi,dx)bind(C,name="gradp")
!c ::
!c :: ----------------------------------------------------------
!c :: Compute a cell centered gradient from a node
!c :: centered field.  Returns all components of GRADP
!c :: ----------------------------------------------------------
!c ::
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

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               gp(i,j,k,1) = ddx*( &
                   p(i+1,j,k  )-p(i,j,k  )+p(i+1,j+1,k  )-p(i,j+1,k  )+ &
                   p(i+1,j,k+1)-p(i,j,k+1)+p(i+1,j+1,k+1)-p(i,j+1,k+1))
            end do
         end do
      end do

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               gp(i,j,k,2) = ddy*( &
                   p(i,j+1,k  )-p(i,j,k  )+p(i+1,j+1,k  )-p(i+1,j,k  )+ &
                   p(i,j+1,k+1)-p(i,j,k+1)+p(i+1,j+1,k+1)-p(i+1,j,k+1))
            end do
         end do
      end do

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               gp(i,j,k,3) = ddz*( &
                   p(i,  j,k+1)-p(i,  j,k)+p(i,  j+1,k+1)-p(i,  j+1,k)+ &
                   p(i+1,j,k+1)-p(i+1,j,k)+p(i+1,j+1,k+1)-p(i+1,j+1,k))
            end do
         end do
      end do

      end subroutine gradp 

!c :: ----------------------------------------------------------
!c :: Replace coarse grid pressure data with corresponding
!c :: fine grid pressure data.
!c ::
!c :: INPUTS / OUTPUTS:
!c ::  crse      <=  coarse grid data
!c ::  DIMS(crse) => index limits of crse
!c ::  fine       => fine grid data
!c ::  DIMS(fine) => index limits of fine
!c ::  lo,hi      => index limits of overlap (crse grid)
!c ::  ratios     => refinement ratio
!c ::
!c :: NOTE:
!c ::  Assumes pressure fields are node based
!c :: ----------------------------------------------------------
!c ::
      subroutine fort_putdown (crse,DIMS(crse), &
     			           fine,DIMS(fine),lo,hi,ratios)&
                     bind(C,name="fort_putdown")

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

      end subroutine fort_putdown 

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

!c :: ----------------------------------------------------------
!c :: SUMMASS
!c ::             MASS = sum{ vol(i,j)*rho(i,j) }
!c ::
!c :: INPUTS / OUTPUTS:
!c ::  rho        => density field
!c ::  DIMS(rho)  => index limits of rho aray
!c ::  lo,hi      => index limits of grid interior
!c ::  delta	 => cell size
!c ::  mass      <=  total mass
!c ::  r		 => radius at cell center
!c ::  tmp        => temp column array
!c :: ----------------------------------------------------------
!c ::
       subroutine summass(rho,DIMS(rho),DIMS(grid),delta,mass)&
                          bind(C,name="summass")

       implicit none

       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  mass, delta(SDIM)
       REAL_T  rho(DIMV(rho))

       integer i, j, k
       REAL_T  vol

       vol = delta(1)*delta(2)*delta(3)

       mass = zero

       do k = ARG_L3(grid), ARG_H3(grid)
          do j = ARG_L2(grid), ARG_H2(grid)
             do i = ARG_L1(grid), ARG_H1(grid)
                mass = mass + rho(i,j,k)
             end do
          end do
       end do

       mass = vol*mass

       end subroutine summass

       
!c :: ----------------------------------------------------------
!c :: SUMMASS
!c ::             MASS = sum{ vol(i,j)*rho(i,j) }
!c ::
!c :: INPUTS / OUTPUTS:
!c ::  rho        => density field
!c ::  DIMS(rho)  => index limits of rho aray
!c ::  lo,hi      => index limits of grid interior
!c ::  delta	 => cell size
!c ::  mass      <=  total mass
!c ::  r		 => radius at cell center
!c ::  tmp        => temp column array
!c :: ----------------------------------------------------------
!c ::
       subroutine summass_eb(rho,DIMS(rho),DIMS(grid),vf,DIMS(vf),delta,mass)&
            bind(C,name="summass_eb")

       implicit none

       integer DIMDEC(rho)
       integer DIMDEC(vf)
       integer DIMDEC(grid)
       REAL_T  mass, delta(SDIM)
       REAL_T  rho(DIMV(rho))
       REAL_T  vf(DIMV(vf))

       integer i, j, k
       REAL_T  vol

       vol = delta(1)*delta(2)*delta(3)

       mass = zero

       do k = ARG_L3(grid), ARG_H3(grid)
          do j = ARG_L2(grid), ARG_H2(grid)
             do i = ARG_L1(grid), ARG_H1(grid)
                mass = mass +  vf(i,j,k)*rho(i,j,k)
             end do
          end do
       end do

       mass = vol*mass

     end subroutine summass_eb

!c :: ----------------------------------------------------------
!c :: SUMMASSCYL
!c ::    MASS = sum{ vol(i,j,k)*rho(i,j,k) } over subregion cylinder
!c ::
!c :: INPUTS / OUTPUTS:
!c ::  rho        => density field
!c ::  DIMS(rho)  => index limits of rho aray
!c ::  lo,hi      => index limits of grid interior
!c ::  delta	 => cell size
!c ::  mass      <=  total mass
!c ::  r		 => radius at cell center
!c ::  tmp        => temp column array
!c :: ----------------------------------------------------------
!c ::
       subroutine summass_cyl(rho,DIMS(rho),DIMS(grid),delta, &
                              plo,vws_dz,vws_Rcyl,mass) &
                              bind(C,name="summass_cyl")

       implicit none

       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  mass, delta(SDIM), plo(SDIM), vws_dz, vws_Rcyl
       REAL_T  rho(DIMV(rho))

       integer i, j, k
       REAL_T  vol, x, y, z, r

       vol = delta(1)*delta(2)*delta(3)

       mass = zero

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

       mass = vol*mass

       end subroutine summass_cyl

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
      
!c-----------------------------------------------------------------------

      subroutine edge_interp (flo, fhi, nc, ratio, dir, &
          fine, fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2)&
           bind(C,name="edge_interp")
           
      implicit none
      integer flo(0:3-1), fhi(0:3-1), nc, ratio(0:3-1), dir
      integer fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2
      DOUBLE PRECISION &
          fine(fine_l0:fine_h0,fine_l1:fine_h1,fine_l2:fine_h2,nc)
      integer i,j,k,n,P,M,L
      DOUBLE PRECISION val, df
!c
!c     Do linear in dir, pc transverse to dir, leave alone the fine values
!c     lining up with coarse edges--assume these have been set to hold the 
!c     values you want to interpolate to the rest.
!c
      if (dir.eq.0) then
         do n=1,nc
            do k=flo(2),fhi(2),ratio(2)
               do j=flo(1),fhi(1),ratio(1)
                  do i=flo(0),fhi(0)-ratio(dir),ratio(0)
                     df = fine(i+ratio(dir),j,k,n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n) &
                            + df*dble(M)/dble(ratio(dir))
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
                        val = fine(i,j,k,n) &
                            + df*dble(M)/dble(ratio(dir))
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
                        val = fine(i,j,k,n) &
                            + df*dble(M)/dble(ratio(dir))
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
      end subroutine edge_interp
      
!c-----------------------------------------------------------------------

      subroutine pc_edge_interp(lo, hi, nc, ratio, dir, &
          crse, crse_l0, crse_l1, crse_l2, crse_h0, crse_h1, crse_h2, &
          fine, fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2)&
          bind(C,name="pc_edge_interp")
          
      implicit none
      integer lo(3),hi(3), nc, ratio(0:3-1), dir
      integer crse_l0, crse_l1, crse_l2, crse_h0, crse_h1, crse_h2
      integer fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2
      DOUBLE PRECISION crse(crse_l0:crse_h0,crse_l1:crse_h1,crse_l2:crse_h2,nc)
      DOUBLE PRECISION fine(fine_l0:fine_h0,fine_l1:fine_h1,fine_l2:fine_h2,nc)
      integer i,j,k,ii,jj,kk,n,L, P
!c
!c     For edge-based data, fill fine values with piecewise-constant interp of coarse data.
!c     Operate only on faces that overlap--ie, only fill the fine faces that make up each
!c     coarse face, leave the in-between faces alone.
!c
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
      end subroutine pc_edge_interp
      
      
! ::: -----------------------------------------------------------
! ::: This routine is intended to be a generic fill function
! ::: for cell-centered data.  It knows how to extrapolate
! ::: and reflect data and is used to supplement the problem-specific
! ::: fill functions which call it.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: q           <=  array to fill
! ::: lo,hi        => index extent of loops
! ::: q_l,q_h      => index extent of q array
! ::: domlo,domhi  => index extent of problem domain
! ::: dx           => cell spacing
! ::: xlo          => physical location of lower left hand
! :::	              corner of q array
! ::: bc	   => array of boundary flags bc(SPACEDIM,lo:hi)
! ::: 
! ::: NOTE: all corner as well as edge data is filled if not EXT_DIR
! ::: -----------------------------------------------------------

    subroutine filcc_tile(l1,l2,l3,h1,h2,h3,q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
         domlo,domhi,dx,xlo,bc) bind(C,name="filcc_tile")

      use amrex_filcc_module, only: filccn
      
      implicit none

      integer    l1,l2,l3,h1,h2,h3
      integer    q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
      integer    domlo(SDIM), domhi(SDIM)
      integer    bc(SDIM,2)
      REAL_T     xlo(SDIM), dx(SDIM)
      REAL_T     q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)

      integer :: q_lo(3), q_hi(3)
      integer    lo(3),hi(3)
      
      lo   = [l1, l2, l3]
      hi   = [h1, h2, h3]
      q_lo = [q_l1, q_l2, q_l3]
      q_hi = [q_h1, q_h2, q_h3]

      call filccn(lo, hi, q, q_lo, q_hi, 1, domlo, domhi, dx, xlo, bc)

    end subroutine filcc_tile
    
end module navierstokes_3d_module
