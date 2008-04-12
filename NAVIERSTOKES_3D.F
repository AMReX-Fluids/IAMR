c
c $Id: NAVIERSTOKES_3D.F,v 1.22 2008-04-12 00:27:49 marc Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "NAVIERSTOKES_F.H"
#include "ArrayLim.H"

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
      integer    lo(SDIM),  hi(SDIM)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  dx(SDIM)
      REAL_T  p(DIMV(p))
      REAL_T  gp(DIMV(gp),SDIM)
      logical    invalid
      integer    i,j,k
      REAL_T     ddx, ddy, ddz
      integer    ilo, ihi, jlo, jhi, klo, khi

      ddx = fourth/dx(1)
      ddy = fourth/dx(2)
      ddz = fourth/dx(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               gp(i,j,k,1) = ddx*(
     &              p(i+1,j,k  )-p(i,j,k  )+p(i+1,j+1,k  )-p(i,j+1,k  )+
     &              p(i+1,j,k+1)-p(i,j,k+1)+p(i+1,j+1,k+1)-p(i,j+1,k+1))

               gp(i,j,k,2) = ddy*(
     &              p(i,j+1,k  )-p(i,j,k  )+p(i+1,j+1,k  )-p(i+1,j,k  )+
     &              p(i,j+1,k+1)-p(i,j,k+1)+p(i+1,j+1,k+1)-p(i+1,j,k+1))

               gp(i,j,k,3) = ddz*(
     &              p(i,  j,k+1)-p(i,  j,k)+p(i,  j+1,k+1)-p(i,  j+1,k)+
     &              p(i+1,j,k+1)-p(i+1,j,k)+p(i+1,j+1,k+1)-p(i+1,j+1,k))
            end do
         end do
      end do

      end

      subroutine FORT_AVGDOWN (
     &     crse,DIMS(crse),nvar,
     &     fine,DIMS(fine),
     &     cv,DIMS(cv),
     &     fv,DIMS(fv),
     &     lo,hi,ratios)
c     ----------------------------------------------------------
c     Volume-weight average the fine grid data onto the coarse
c     grid.  Overlap is given in coarse grid coordinates.
c     
c     crse      =  coarse grid data
c     nvar	= number of components in arrays
c     fine      = fine grid data
c     cv        = coarse grid volume array
c     fv        = fine grid volume array
c     lo,hi     = index limits of overlap (crse grid)
c     ratios    = IntVect refinement ratio
c     ----------------------------------------------------------
      implicit none

      integer  DIMDEC(crse)
      integer  DIMDEC(cv)
      integer  DIMDEC(fine)
      integer  DIMDEC(fv)
      integer  lo(SDIM), hi(SDIM)
      integer  nvar
      integer  ratios(SDIM)
      REAL_T   crse(DIMV(crse),nvar)
      REAL_T     cv(DIMV(cv))
      REAL_T   fine(DIMV(fine),nvar)
      REAL_T     fv(DIMV(fv))

      integer  i, j, k, n, ic, jc, kc, ioff, joff, koff
      integer  lratx, lraty, lratz

      lratx = ratios(1)
      lraty = ratios(2)
      lratz = ratios(3)

      do n = 1, nvar
c
c     set coarse grid to zero on overlap
c
         do kc = lo(3), hi(3)
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,kc,n) = zero
               end do
            end do
         end do
c         
c     sum fine data
c
         do koff = 0, lratz-1
            do kc = lo(3),hi(3)
               k = kc*lratz + koff
               do joff = 0, lraty-1
                  do jc = lo(2), hi(2)
                     j = jc*lraty + joff
                     do ioff = 0, lratx-1
                        do ic = lo(1), hi(1)
                           i = ic*lratx + ioff
                           crse(ic,jc,kc,n) = crse(ic,jc,kc,n) +
     &                          fv(i,j,k)*fine(i,j,k,n)
                        end do
                     end do
                  end do
               end do
            end do
         end do
c
c     divide out by volume weight
c
         do kc = lo(3), hi(3)
            do ic = lo(1), hi(1)
               do jc = lo(2), hi(2)
                  crse(ic,jc,kc,n) = crse(ic,jc,kc,n)/cv(ic,jc,kc)
               end do
            end do
         end do
      end do

      end

      subroutine FORT_TESTINJECT(
     &     crse,DIMS(crse),
     &     fine,DIMS(fine),lo,hi,ratios)
c
c     test for consistency of pressure nodes
c
      implicit none

      integer  DIMDEC(crse)
      integer  DIMDEC(fine)
      integer  lo(3), hi(3)
      integer  ratios(3)
      REAL_T   crse(DIMV(crse))
      REAL_T   fine(DIMV(fine))

      REAL_T diff
      integer  ic, jc, kc
      integer  lratx, lraty, lratz

      lratx = ratios(1)
      lraty = ratios(2)
      lratz = ratios(3)

      diff = zero
      do kc = lo(3), hi(3)
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               diff = max( diff,
     &              abs(crse(ic,jc,kc)-fine(lratx*ic,lraty*jc,lratz*kc)))
            end do
         end do
      end do

      if ( diff .gt. 1.0D-10 ) then
         call bl_abort('mismatch between fine and coarse nodes')
      end if

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

      integer  ic, jc, kc
      integer  lratx, lraty, lratz

      lratx = ratios(1)
      lraty = ratios(2)
      lratz = ratios(3)

      do kc = lo(3), hi(3)
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,kc) = fine(lratx*ic,lraty*jc,lratz*kc)
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

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                a(i,j,k) = a(i,j,k) + alpha*b(i,j,k)
             end do
          end do
       end do

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
      REAL_T  dat(DIMV(dat),4)
      REAL_T  pres(DIMV(pres),4)
      REAL_T  turb(0:ksize*turbVars-1)
      
      integer i, j, k, v
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      
      REAL_T area
      REAL_T rho, ux, uy, uz, p, XC, XMg
      REAL_T e, T, f, wC, wMg
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
      REAL_T  dat(DIMV(dat),24)
      REAL_T  pres(DIMV(pres),4)
      REAL_T  jetData(0:ksize*jetVars-1)
      REAL_T  gridDx
      
      integer i, j, k, v
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      
      REAL_T  rho, ux, uy, uz, p, XC, XMg
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

      integer idx_rho, idx_ux, idx_uy, idx_uz, idx_trac, idx_temp
      integer idx_p, idx_dpdx, idx_dpdy, idx_dpdz
      integer var
      
      integer isioproc
      
#include "probdata.H"

      call bl_pd_is_ioproc(isioproc)

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
      dxnn = dx/float(nn)
      dynn = dy/float(nn)
c     Note this is right, it really should be kn not nn
      dznn = dz/float(kn)
c     Mod dx
      mdx = sqrt(dxnn*dxnn+dynn*dynn)
      
      gridDx=dx/dfloat(kn)
      
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
         zlow = xlo(3) + dz*dfloat(k-klo)
c     And x centre/er
         zctr = zlow + 0.5*dz
         
         do j = jlo, jhi
c     Calculate front of cell
            ylow = xlo(2) + dy*dfloat(j-jlo) - jet_y
c     And y centre/er
            yctr = ylow + 0.5*dy
            
            do i = ilo, ihi
c     Calculate left hand edge of cell
               xlow = xlo(1) + dx*dfloat(i-ilo) - jet_x
c     And x centre/er
               xctr = xlow + 0.5*dx
               
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
                     yy = ylow + dynn*(dfloat(jj)-half)
                     
                     do ii=1,nn
                        xx = xlow + dxnn*(dfloat(ii)-half)
                        
                        rr = sqrt(xx*xx+yy*yy)
                        
                        ridx = int(rr/gridDx)
                        
                        if (ridx.lt.rsize) then
                           
                           do kk=kklo,kkhi
                              
                              zz = zlow + dznn*(dfloat(kk-kklo)+half)
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

       subroutine FORT_MAXVAL(rho,DIMS(rho),DIMS(grid),tmp,mass)

       implicit none

       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  rho(DIMV(rho))
       REAL_T  tmp(DIM3(grid))
       REAL_T  mass

       integer i, j, k

       do k = ARG_L3(grid), ARG_H3(grid)
          tmp(k) = zero
       end do

       do k = ARG_L3(grid), ARG_H3(grid)
          do i = ARG_L1(grid), ARG_H1(grid)
             do j = ARG_L2(grid), ARG_H2(grid)
                tmp(k) = max( tmp(k) , rho(i,j,k) )
             end do
          end do
       end do

       mass = zero
       do j = ARG_L3(grid), ARG_H3(grid)
          mass = max(mass , tmp(k))
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
c ::  delta	 => cell size
c ::  mass      <=  total mass
c ::  r		 => radius at cell center
c ::  tmp        => temp column array
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SUMMASS(rho,DIMS(rho),DIMS(grid),delta,mass,tmp)

       implicit none

       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  mass, delta(SDIM)
       REAL_T  rho(DIMV(rho))
       REAL_T  tmp(DIM2(grid))

       integer i, j, k
       REAL_T  vol

       vol = delta(1)*delta(2)*delta(3)

       do j = ARG_L2(grid), ARG_H2(grid)
          tmp(j) = zero
       end do

       do k = ARG_L3(grid), ARG_H3(grid)
          do i = ARG_L1(grid), ARG_H1(grid)
             do j = ARG_L2(grid), ARG_H2(grid)
                tmp(j) = tmp(j) + vol*rho(i,j,k)
             end do
          end do
       end do

       mass = zero
       do j = ARG_L2(grid), ARG_H2(grid)
          mass = mass + tmp(j)
       end do

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
     &                             plo,vws_dz,vws_Rcyl,mass,tmp)

       implicit none

       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  mass, delta(SDIM), plo(SDIM), vws_dz, vws_Rcyl
       REAL_T  rho(DIMV(rho))
       REAL_T  tmp(DIM2(grid))

       integer i, j, k
       REAL_T  vol, x, y, z, r

       vol = delta(1)*delta(2)*delta(3)

       do j = ARG_L2(grid), ARG_H2(grid)
          tmp(j) = zero
       end do

       do k = ARG_L3(grid), ARG_H3(grid)
          z = plo(3) + (k+half)*delta(3)
          if (z-plo(3) .le. vws_dz) then 
             do i = ARG_L1(grid), ARG_H1(grid)
                x = plo(1) + (i+half)*delta(1)
                do j = ARG_L2(grid), ARG_H2(grid)
                   y = plo(2) + (j+half)*delta(2)
                   r = SQRT(x*x + y*y)
                   if (r .le. vws_Rcyl) then
                      tmp(j) = tmp(j) + vol*rho(i,j,k)
                   end if
                end do
             end do
          end if
       end do

       mass = zero
       do j = ARG_L2(grid), ARG_H2(grid)
          mass = mass + tmp(j)
       end do

       end

c :: ----------------------------------------------------------
c :: SUMMASSSQ
c ::             MASS = sum{ vol(i,j)*rho(i,j)*rho(i,j) }
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
       subroutine FORT_SUMMASSSQ(rho,DIMS(rho),DIMS(grid),delta,mass,tmp)

       implicit none

       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  mass, delta(SDIM)
       REAL_T  rho(DIMV(rho))
       REAL_T  tmp(DIM2(grid))

       integer i, j, k
       REAL_T  vol

       vol = delta(1)*delta(2)*delta(3)

       do j = ARG_L2(grid), ARG_H2(grid)
          tmp(j) = zero
       end do

       do k = ARG_L3(grid), ARG_H3(grid)
          do i = ARG_L1(grid), ARG_H1(grid)
             do j = ARG_L2(grid), ARG_H2(grid)
                tmp(j) = tmp(j) + vol*rho(i,j,k)*rho(i,j,k)
             end do
          end do
       end do

       mass = zero
       do j = ARG_L2(grid), ARG_H2(grid)
          mass = mass + tmp(j)
       end do

       end

c :: ----------------------------------------------------------
c :: SUMPRESS
c ::             SUM = sum{ press(i,j,k) * surroundingVolume(i,j,k) }
c ::
c :: INPUTS / OUTPUTS:
c ::  press      => pressure field on nodes
c ::  DIMS(press)=> index limits of press aray
c ::  vol        => volume field (zero-ed out under finer grids)
c ::  DIMS(vol)  => index limits of vol aray
c ::  lo,hi      => index limits of grid interior in nodal space
c ::  total_press<= sum of volume-weighted nodal pressures
c ::  total_vol  <= sum of volumes
c :: ----------------------------------------------------------
c ::
       subroutine FORT_SUMPRESS(press,DIMS(press),vol,DIMS(vol),lo,hi,
     $                          total_press,total_vol)

       implicit none

       integer lo(3),hi(3)
       integer DIMDEC(press)
       integer DIMDEC(vol)
       integer DIMDEC(grid)
       REAL_T  press(DIMV(press))
       REAL_T  vol(DIMV(vol))
       REAL_T  total_press,total_vol

       integer i, j, k
       REAL_T  surroundingVolume

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                surroundingVolume = vol(i,j  ,k  ) + vol(i-1,j  ,k  ) + 
     $               vol(i,j-1,k  ) + vol(i-1,j-1,k  ) +
     $               vol(i,j  ,k-1) + vol(i-1,j  ,k-1) + 
     $               vol(i,j-1,k-1) + vol(i-1,j-1,k-1)
                total_press = total_press + press(i,j,k) * surroundingVolume
                total_vol   = total_vol   +                surroundingVolume

             end do
          end do
       end do

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

      if ( isharm .eq. 0 ) then
         if (dir .EQ. 0) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1)+1, hi(1)
                        efab(i,j,k,n) =
     &                       half*(cfab(i,j,k,n) + cfab(i-1,j,k,n))
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 1) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2)+1, hi(2)
                     do i = lo(1), hi(1)
                        efab(i,j,k,n) =
     &                       half*(cfab(i,j,k,n) + cfab(i,j-1,k,n))
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 2) then
            do n = 1,nc
               do k = lo(3)+1, hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        efab(i,j,k,n) =
     &                       half*(cfab(i,j,k,n) + cfab(i,j,k-1,n))
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
                     do i = lo(1)+1, hi(1)
                        if((cfab(i,j,k,n) * cfab(i-1,j,k,n)) .gt.0.d0)
     &                       then
                           efab(i,j,k,n) =
     &                          2.d0*(cfab(i,j,k,n) * cfab(i-1,j,k,n))/
     &                          (cfab(i,j,k,n) + cfab(i-1,j,k,n))
                        else
                           efab(i,j,k,n)=0.d0
                        endif
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 1) then
            do n = 1,nc
               do k = lo(3), hi(3)
                  do j = lo(2)+1, hi(2)
                     do i = lo(1), hi(1)
                        if((cfab(i,j,k,n) * cfab(i,j-1,k,n)).gt.0.d0)
     &                       then
                           efab(i,j,k,n) =
     &                          2.d0*(cfab(i,j,k,n) * cfab(i,j-1,k,n))/
     &                          (cfab(i,j,k,n) + cfab(i,j-1,k,n))
                        else
                           efab(i,j,k,n)=0.d0
                        endif
                     end do
                  end do
               end do
            end do
         else if (dir .EQ. 2) then
            do n = 1,nc
               do k = lo(3)+1, hi(3)
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        if((cfab(i,j,k,n) * cfab(i,j,k-1,n)).gt.0.d0)
     &                       then
                           efab(i,j,k,n) =
     &                          2.d0*(cfab(i,j,k,n) * cfab(i,j,k-1,n))/
     &                          (cfab(i,j,k,n) + cfab(i,j,k-1,n))
                        else
                           efab(i,j,k,n)=0.d0
                        endif
                     end do
                  end do
               end do
            end do
         end if
      end if
      end
c-----------------------------------------------------------------------
      subroutine EDGE_INTERP(flo, fhi, nc, ratio, dir,
     &     fine, fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2)
      implicit none
      integer flo(0:3-1), fhi(0:3-1), nc, ratio(0:3-1), dir
      integer fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2
      DOUBLE PRECISION
     &     fine(fine_l0:fine_h0,fine_l1:fine_h1,fine_l2:fine_h2,nc)
      integer i,j,k,ii,jj,kk,n,P,M,L,clo(0:3-1),chi(0:3-1)
      DOUBLE PRECISION val, df

c     Do linear in dir, pc transverse to dir, leave alone the fine values
c     lining up with coarse edges--assume these have been set to hold the 
c     values you want to interpolate to the rest.
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
      DOUBLE PRECISION val, dc

c     For edge-based data, fill fine values with piecewise-constant interp of coarse data.
c     Operate only on faces that overlap--ie, only fill the fine faces that make up each
c     coarse face, leave the in-between faces alone.
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
