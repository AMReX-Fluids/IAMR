
#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>
#include <CONSTANTS.H>
#include <BC_TYPES.H>
#include <ArrayLim.H>
#include <infl_frc.H>
#include <FLUCTFILE.H>

#define SDIM 3

c ::: -----------------------------------------------------------
c ::: This routine does the interpolation of the inflow data
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: time              =>  Time at which to fill the data
c ::: xlo               =>  Lower physical location of the inflDat array
c ::: dx                =>  Grid spacing in the inflDat array
c ::: storePnt          =>  Interpolation point in the strmwse_dir direction
c ::: fillComp          =>  Component from the storDat array to use
c ::: timePnt           =>  Physical time corresponding to the storePnt point
c :::                         in the data in storDat
c ::: dtFile            =>  Time step in the array storDat
c ::: nCompStorDat      =>  Number of components in the storDat array
c ::: FF_DIMS(storDat)  =>  Dimensions of storDat
c ::: storDat           =>  Array to fill
c ::: DIMS(inflDat)     =>  Dimensions of the inflDat array
c ::: inflDat          <=   Array to fill by interpolating from the storDat 
c :::                         array
c ::: -----------------------------------------------------------
c :::
c ::: NOTE:  When x, y and z are calculated for each i, j, and k index
c :::        in this routine, they are calculated based on formulaes
c :::        of the form,
c :::            x = probLo(1) + dx(1) * (FLOAT(i) + half)
c :::        rather than on the more common form,
c :::            x = xlo(1) + dx(1) * (FLOAT(i - lo(1)) + half)
c :::        In other words, they are calculated based on the origin
c :::        of the problem domain rather than the origin of the FAB
c :::        we are working on.  This forces the physical location for
c :::        some specified index to be identically the same no matter
c :::        which FAB we are working on.  Seems like a moot point, but
c :::        it isn't when we do an integer shift to force points outside 
c :::        of the domain in periodic directions.
c :::
c :::        For instance if we have a domain from 0->63 in a periodic 
c :::        direction and pass in a FAB to fill with indices from -1 to 5,
c :::        the -1 index is integer shifted over to an index of 63 to do
c :::        the interpolation for the forcing.
c :::
      subroutine INTRP_DATA(time, xlo, dx, storePnt, fillComp, filePnt,
     $     timePnt, dtFile, dxFile, xloFile, xhiFile, nCompStorDat, 
     $     FF_DIMS(storDat), storDat,
     $     DIMS(inflDat), inflDat, bc, probLo, probHi)

      implicit none

      integer storePnt, fillComp, filePnt, nCompStorDat
      integer loStoreDim(3), hiStoreDim(3)
      integer DIMDEC(inflDat)
      integer FF_DIMDEC(storDat)
      integer bc(SDIM,2)

      REAL_T time, timePnt, dtFile
      REAL_T xlo(SDIM), dx(SDIM), dxFile(3), xloFile(3), xhiFile(3)
      REAL_T probLo(SDIM), probHi(SDIM)
      REAL_T inflDat(DIMV(inflDat))
      REAL_T storDat(FF_DIMV(storDat),nCompStorDat)
      REAL_T lengthFile1, lengthFile2

      integer i, j, k, iloc, jloc, kloc, ixcalc, jycalc, kzcalc
      integer lo(SDIM), hi(SDIM), n_lev_cells(SDIM)
      integer loStorDim(3), hiStorDim(3)

      REAL_T tm1, tp1, ctm1, ct0, ctp1, coff,
     $       x, xm1, x0, xp1, cxm1, cx0, cxp1,
     $       y, ym1, y0, yp1, cym1, cy0, cyp1,
     $       z, zm1, z0, zp1, czm1, cz0, czp1,
     $       valm1, val0, valp1

#include <INFL_FORCE_F.H>

      call SET_LOHI(DIMS(inflDat), lo, hi)
      call FF_SET_LOHI(FF_DIMS(storDat), loStorDim, hiStorDim)

      do i = 1, SDIM
        n_lev_cells(i) = (probHi(i) - probLo(i)) / dx(i)
      enddo

      if (infl_type .eq. infl_periodic_type) then
         coff = one
         tm1  = timePnt - dtFile
         tp1  = timePnt + dtFile
      else if (infl_type .eq. infl_swirl_type) then
         coff = half
         tm1  = fluct_times(filePnt - 1)
         tp1  = fluct_times(filePnt + 1)
      else
         call bl_abort('INTRP_DATA: unknown infl_type')
      endif

      ctm1 = (time - timePnt) * (time - tp1) / (tm1 - timePnt) / (tm1 - tp1)
      ct0  = (time - tm1) * (time - tp1) / (timePnt - tm1) / (timePnt - tp1)
      ctp1 = (time - tm1) * (time - timePnt) / (tp1 - tm1) / (tp1 - timePnt)
c
c     ----------------------------------------------
c     :::: Interpolate to fill the inflow array ::::
c     ----------------------------------------------
c
      if (strmwse_dir .eq. FLCT_XVEL) then
c
c     :::: Streamwise X-direction ::::
c
      do k = lo(3), hi(3)
        kzcalc = k
        if (bc(3,1) .eq. INT_DIR .and. k .lt. 0) then
          kzcalc = k + n_lev_cells(3)
        else if (bc(3,2) .eq. INT_DIR .and. k .ge. n_lev_cells(3)) then
          kzcalc = k - n_lev_cells(3)
        endif
        z = probLo(3) + dx(3)*(float(kzcalc) + half)

        if (xloFile(3) .le. z .and. z .le. xhiFile(3)) then
          kloc = (z - xloFile(3) + half*dxFile(3))/dxFile(3) + 1
          kloc = MAX(kloc, loStorDim(3)+1)
          kloc = MIN(kloc, hiStorDim(3)-1)
          z0 = (FLOAT(kloc)-coff) * dxFile(3) + xloFile(3)
          zm1 = z0 - dxFile(3)
          zp1 = z0 + dxFile(3)
          czm1 = (z - z0)  * (z - zp1) / (zm1 - z0)  / (zm1 - zp1)
          cz0  = (z - zm1) * (z - zp1) / (z0  - zm1) / (z0  - zp1)
          czp1 = (z - zm1) * (z - z0)  / (zp1 - zm1) / (zp1 - z0)
        endif

        do j = lo(2), hi(2)
          jycalc = j
          if (bc(2,1) .eq. INT_DIR .and. j .lt. 0) then
            jycalc = j + n_lev_cells(2)
          else if (bc(2,2) .eq. INT_DIR .and. j .ge. n_lev_cells(2)) then
            jycalc = j - n_lev_cells(2)
          endif
          y = probLo(2) + dx(2)*(float(jycalc) + half)

          if ((xloFile(2) .le. y .and. y .le. xhiFile(2)) .and.
     $        (xloFile(3) .le. z .and. z .le. xhiFile(3))) then
            jloc = (y - xloFile(2) + half*dxFile(2))/dxFile(2) + 1
            jloc = MAX(jloc, loStorDim(2)+1)
            jloc = MIN(jloc, hiStorDim(2)-1)
            y0 = (FLOAT(jloc)-coff) * dxFile(2) + xloFile(2)
            ym1 = y0 - dxFile(2)
            yp1 = y0 + dxFile(2)
            cym1 = (y - y0)  * (y - yp1) / (ym1 - y0)  / (ym1 - yp1)
            cy0  = (y - ym1) * (y - yp1) / (y0  - ym1) / (y0  - yp1)
            cyp1 = (y - ym1) * (y - y0)  / (yp1 - ym1) / (yp1 - y0)

            val0 = cym1 * (ctm1 * storDat(storePnt-1,jloc-1,kloc,fillComp)
     $                   +  ct0 * storDat(storePnt,  jloc-1,kloc,fillComp)
     $                   + ctp1 * storDat(storePnt+1,jloc-1,kloc,fillComp))
     $           +  cy0 * (ctm1 * storDat(storePnt-1,jloc,  kloc,fillComp)
     $                   +  ct0 * storDat(storePnt,  jloc,  kloc,fillComp)
     $                   + ctp1 * storDat(storePnt+1,jloc,  kloc,fillComp))
     $           + cyp1 * (ctm1 * storDat(storePnt-1,jloc+1,kloc,fillComp)
     $                   +  ct0 * storDat(storePnt,  jloc+1,kloc,fillComp)
     $                   + ctp1 * storDat(storePnt+1,jloc+1,kloc,fillComp))

            valm1 = cym1 * (ctm1 * storDat(storePnt-1,jloc-1,kloc-1,fillComp)
     $                    +  ct0 * storDat(storePnt,  jloc-1,kloc-1,fillComp)
     $                    + ctp1 * storDat(storePnt+1,jloc-1,kloc-1,fillComp))
     $            +  cy0 * (ctm1 * storDat(storePnt-1,jloc,  kloc-1,fillComp)
     $                    +  ct0 * storDat(storePnt,  jloc,  kloc-1,fillComp)
     $                    + ctp1 * storDat(storePnt+1,jloc,  kloc-1,fillComp))
     $            + cyp1 * (ctm1 * storDat(storePnt-1,jloc+1,kloc-1,fillComp)
     $                    +  ct0 * storDat(storePnt,  jloc+1,kloc-1,fillComp)
     $                    + ctp1 * storDat(storePnt+1,jloc+1,kloc-1,fillComp))

            valp1 = cym1 * (ctm1 * storDat(storePnt-1,jloc-1,kloc+1,fillComp)
     $                    +  ct0 * storDat(storePnt,  jloc-1,kloc+1,fillComp)
     $                    + ctp1 * storDat(storePnt+1,jloc-1,kloc+1,fillComp))
     $            +  cy0 * (ctm1 * storDat(storePnt-1,jloc,  kloc+1,fillComp)
     $                    +  ct0 * storDat(storePnt,  jloc,  kloc+1,fillComp)
     $                    + ctp1 * storDat(storePnt+1,jloc,  kloc+1,fillComp))
     $            + cyp1 * (ctm1 * storDat(storePnt-1,jloc+1,kloc+1,fillComp)
     $                    +  ct0 * storDat(storePnt,  jloc+1,kloc+1,fillComp)
     $                    + ctp1 * storDat(storePnt+1,jloc+1,kloc+1,fillComp))
            val0 = czm1 * valm1 + cz0 * val0 + czp1 * valp1
          else
            val0 = zero
          endif

          do i = lo(1), hi(1)
            inflDat(i,j,k) = val0
          enddo
        enddo
      enddo
      
      elseif (strmwse_dir .eq. FLCT_YVEL) then
c
c     :::: Streamwise Y-direction ::::
c
      do k = lo(3), hi(3)
        kzcalc = k
        if (bc(3,1) .eq. INT_DIR .and. k .lt. 0) then
          kzcalc = k + n_lev_cells(3)
        else if (bc(3,2) .eq. INT_DIR .and. k .ge. n_lev_cells(3)) then
          kzcalc = k - n_lev_cells(3)
        endif
        z = probLo(3) + dx(3)*(float(kzcalc) + half)

        if (xloFile(3) .le. z .and. z .le. xhiFile(3)) then
          kloc = (z - xloFile(3) + half*dxFile(3))/dxFile(3) + 1
          kloc = MAX(kloc, loStorDim(3)+1)
          kloc = MIN(kloc, hiStorDim(3)-1)
          z0 = (FLOAT(kloc)-coff) * dxFile(3) + xloFile(3)
          zm1 = z0 - dxFile(3)
          zp1 = z0 + dxFile(3)
          czm1 = (z - z0)  * (z - zp1) / (zm1 - z0)  / (zm1 - zp1)
          cz0  = (z - zm1) * (z - zp1) / (z0  - zm1) / (z0  - zp1)
          czp1 = (z - zm1) * (z - z0)  / (zp1 - zm1) / (zp1 - z0)
        endif

        do i = lo(1), hi(1)
          ixcalc = i
          if (bc(1,1) .eq. INT_DIR .and. i .lt. 0) then
            ixcalc = i + n_lev_cells(1)
          else if (bc(1,2) .eq. INT_DIR .and. i .ge. n_lev_cells(1)) then
            ixcalc = i - n_lev_cells(1)
          endif
          x = probLo(1) + dx(1)*(float(ixcalc) + half)

          if ((xloFile(1) .le. x .and. x .le. xhiFile(1)) .and.
     $        (xloFile(3) .le. z .and. z .le. xhiFile(3))) then
            iloc = (x - xloFile(1) + half*dxFile(1))/dxFile(1) + 1
            iloc = MAX(iloc, loStorDim(1)+1)
            iloc = MIN(iloc, hiStorDim(1)-1)
            x0 = (FLOAT(iloc)-coff) * dxFile(1) + xloFile(1)
            xm1 = x0 - dxFile(1)
            xp1 = x0 + dxFile(1)
            cxm1 = (x - x0)  * (x - xp1) / (xm1 - x0)  / (xm1 - xp1)
            cx0  = (x - xm1) * (x - xp1) / (x0  - xm1) / (x0  - xp1)
            cxp1 = (x - xm1) * (x - x0)  / (xp1 - xm1) / (xp1 - x0)

            val0 = cxm1 * (ctm1 * storDat(iloc-1,storePnt-1,kloc,fillComp)
     $                   +  ct0 * storDat(iloc-1,storePnt,  kloc,fillComp)
     $                   + ctp1 * storDat(iloc-1,storePnt+1,kloc,fillComp))
     $           +  cx0 * (ctm1 * storDat(iloc,  storePnt-1,kloc,fillComp)
     $                   +  ct0 * storDat(iloc,  storePnt,  kloc,fillComp)
     $                   + ctp1 * storDat(iloc,  storePnt+1,kloc,fillComp))
     $           + cxp1 * (ctm1 * storDat(iloc+1,storePnt-1,kloc,fillComp)
     $                   +  ct0 * storDat(iloc+1,storePnt,  kloc,fillComp)
     $                   + ctp1 * storDat(iloc+1,storePnt+1,kloc,fillComp))

            valm1 = cxm1 * (ctm1 * storDat(iloc-1,storePnt-1,kloc-1,fillComp)
     $                    +  ct0 * storDat(iloc-1,storePnt,  kloc-1,fillComp)
     $                    + ctp1 * storDat(iloc-1,storePnt+1,kloc-1,fillComp))
     $            +  cx0 * (ctm1 * storDat(iloc,  storePnt-1,kloc-1,fillComp)
     $                    +  ct0 * storDat(iloc,  storePnt,  kloc-1,fillComp)
     $                    + ctp1 * storDat(iloc,  storePnt+1,kloc-1,fillComp))
     $            + cxp1 * (ctm1 * storDat(iloc+1,storePnt-1,kloc-1,fillComp)
     $                    +  ct0 * storDat(iloc+1,storePnt,  kloc-1,fillComp)
     $                    + ctp1 * storDat(iloc+1,storePnt+1,kloc-1,fillComp))

            valp1 = cxm1 * (ctm1 * storDat(iloc-1,storePnt-1,kloc+1,fillComp)
     $                    +  ct0 * storDat(iloc-1,storePnt,  kloc+1,fillComp)
     $                    + ctp1 * storDat(iloc-1,storePnt+1,kloc+1,fillComp))
     $            +  cx0 * (ctm1 * storDat(iloc,  storePnt-1,kloc+1,fillComp)
     $                    +  ct0 * storDat(iloc,  storePnt,  kloc+1,fillComp)
     $                    + ctp1 * storDat(iloc,  storePnt+1,kloc+1,fillComp))
     $            + cxp1 * (ctm1 * storDat(iloc+1,storePnt-1,kloc+1,fillComp)
     $                    +  ct0 * storDat(iloc+1,storePnt,  kloc+1,fillComp)
     $                    + ctp1 * storDat(iloc+1,storePnt+1,kloc+1,fillComp))
            val0 = czm1 * valm1 + cz0 * val0 + czp1 * valp1
          else
            val0 = zero
          endif

          do j = lo(2), hi(2)
            inflDat(i,j,k) = val0
          enddo
        enddo
      enddo
      
      elseif (strmwse_dir .eq. FLCT_ZVEL) then
c
c     :::: Streamwise Z-direction ::::
c
c      write (*,*) "dx,lo,hi(1):",dxFile(1),xloFile(1),xhiFile(1)
c      write (*,*) "dx,lo,hi(2):",dxFile(2),xloFile(2),xhiFile(2)
      lengthFile1 = (xhiFile(1)-xloFile(1))-dxFile(1)*two
      lengthFile2 = (xhiFile(2)-xloFile(2))-dxFile(2)*two
      do j = lo(2), hi(2)
#if 1
        jycalc = j
        if (bc(2,1) .eq. INT_DIR .and. j .lt. 0) then
          jycalc = j + n_lev_cells(2)
        else if (bc(2,2) .eq. INT_DIR .and. j .ge. n_lev_cells(2)) then
          jycalc = j - n_lev_cells(2)
        endif
        y = probLo(2) + dx(2)*(float(jycalc) + half)
#else
        y = probLo(2) + dx(2)*(float(j) + half)
        do while (y.gt.xhiFile(2)-dxFile(2))
           y = y - lengthFile2
        end do
        do while (y.lt.xloFile(2)+dxFile(2))
           y = y + lengthFile2
        end do
#endif

        if (xloFile(2) .le. y .and. y .le. xhiFile(2)) then
          jloc = (y - xloFile(2) + half*dxFile(2))/dxFile(2) + 1
          jloc = MAX(jloc, loStorDim(2)+1)
          jloc = MIN(jloc, hiStorDim(2)-1)
          y0 = (FLOAT(jloc)-coff) * dxFile(2) + xloFile(2)
          ym1 = y0 - dxFile(2)
          yp1 = y0 + dxFile(2)
          cym1 = (y - y0)  * (y - yp1) / (ym1 - y0)  / (ym1 - yp1)
          cy0  = (y - ym1) * (y - yp1) / (y0  - ym1) / (y0  - yp1)
          cyp1 = (y - ym1) * (y - y0)  / (yp1 - ym1) / (yp1 - y0)
        endif

        do i = lo(1), hi(1)
#if 1
          ixcalc = i
          if (bc(1,1) .eq. INT_DIR .and. i .lt. 0) then
            ixcalc = i + n_lev_cells(1)
          else if (bc(1,2) .eq. INT_DIR .and. i .ge. n_lev_cells(1)) then
            ixcalc = i - n_lev_cells(1)
          endif
          x = probLo(1) + dx(1)*(float(ixcalc) + half)
#else
          x = probLo(1) + dx(1)*(float(i) + half)
          do while (x.gt.xhiFile(1)-dxFile(1))
             x = x - lengthFile1
          end do
          do while (x.lt.xloFile(1)+dxFile(1))
             x = x + lengthFile1
          end do
#endif

          if ((xloFile(1) .le. x .and. x .le. xhiFile(1)) .and. 
     $        (xloFile(2) .le. y .and. y .le. xhiFile(2))) then
            iloc = (x - xloFile(1) + half*dxFile(1))/dxFile(1) + 1
            iloc = MAX(iloc, loStorDim(1)+1)
            iloc = MIN(iloc, hiStorDim(1)-1)
            x0 = (FLOAT(iloc)-coff) * dxFile(1) + xloFile(1)
            xm1 = x0 - dxFile(1)
            xp1 = x0 + dxFile(1)
            cxm1 = (x - x0)  * (x - xp1) / (xm1 - x0)  / (xm1 - xp1)
            cx0  = (x - xm1) * (x - xp1) / (x0  - xm1) / (x0  - xp1)
            cxp1 = (x - xm1) * (x - x0)  / (xp1 - xm1) / (xp1 - x0)

            val0 = cxm1 * (ctm1 * storDat(iloc-1,jloc,storePnt-1,fillComp)
     $                   +  ct0 * storDat(iloc-1,jloc,storePnt,  fillComp)
     $                   + ctp1 * storDat(iloc-1,jloc,storePnt+1,fillComp))
     $           +  cx0 * (ctm1 * storDat(iloc,  jloc,storePnt-1,fillComp)
     $                   +  ct0 * storDat(iloc,  jloc,storePnt,  fillComp)
     $                   + ctp1 * storDat(iloc,  jloc,storePnt+1,fillComp))
     $           + cxp1 * (ctm1 * storDat(iloc+1,jloc,storePnt-1,fillComp)
     $                   +  ct0 * storDat(iloc+1,jloc,storePnt,  fillComp)
     $                   + ctp1 * storDat(iloc+1,jloc,storePnt+1,fillComp))

            valm1 = cxm1 * (ctm1 * storDat(iloc-1,jloc-1,storePnt-1,fillComp)
     $                    +  ct0 * storDat(iloc-1,jloc-1,storePnt,  fillComp)
     $                    + ctp1 * storDat(iloc-1,jloc-1,storePnt+1,fillComp))
     $            +  cx0 * (ctm1 * storDat(iloc,  jloc-1,storePnt-1,fillComp)
     $                    +  ct0 * storDat(iloc,  jloc-1,storePnt,  fillComp)
     $                    + ctp1 * storDat(iloc,  jloc-1,storePnt+1,fillComp))
     $            + cxp1 * (ctm1 * storDat(iloc+1,jloc-1,storePnt-1,fillComp)
     $                    +  ct0 * storDat(iloc+1,jloc-1,storePnt,  fillComp)
     $                    + ctp1 * storDat(iloc+1,jloc-1,storePnt+1,fillComp))

            valp1 = cxm1 * (ctm1 * storDat(iloc-1,jloc+1,storePnt-1,fillComp)
     $                    +  ct0 * storDat(iloc-1,jloc+1,storePnt,  fillComp)
     $                    + ctp1 * storDat(iloc-1,jloc+1,storePnt+1,fillComp))
     $            +  cx0 * (ctm1 * storDat(iloc,  jloc+1,storePnt-1,fillComp)
     $                    +  ct0 * storDat(iloc,  jloc+1,storePnt,  fillComp)
     $                    + ctp1 * storDat(iloc,  jloc+1,storePnt+1,fillComp))
     $            + cxp1 * (ctm1 * storDat(iloc+1,jloc+1,storePnt-1,fillComp)
     $                    +  ct0 * storDat(iloc+1,jloc+1,storePnt,  fillComp)
     $                    + ctp1 * storDat(iloc+1,jloc+1,storePnt+1,fillComp))
            val0 = cym1 * valm1 + cy0 * val0 + cyp1 * valp1
          else
            val0 = zero
          endif

          do k = lo(3), hi(3)
            inflDat(i,j,k) = val0
          enddo
        enddo
      enddo

      endif

      end
