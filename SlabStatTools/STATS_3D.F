c
c $Id: STATS_3D.F,v 1.1 2000-07-19 21:01:29 sstanley Exp $
c

#include "REAL.H"
#include "CONSTANTS.H"
#include "ArrayLim.H"

#include "StatTypes.H"
#define SDIM 3


c ::: -----------------------------------------------------------
c ::: This routine calculates the velocity statistics using Reynolds
c ::: averaging.
c :::
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: ssDat        =>  (const)  SlabStat data
c ::: nSSComp      =>  (const)  Number of components in ssDat
c ::: DIMS(dat)    =>  (const)  Dimensions of SlabStat data
c ::: nStats       =>  (const)  Number of statistics to calculate
c ::: nStns        =>  (const)  Number of stations in the statistic array
c ::: stats       <=   (modify) Array to hold statistics
c :::
c :::                  Output Values:
c :::                      <U> <V> <uu> <vv> <uv> TKE
c :::                  
c ::: vblo, vbhi   =>  (const)  subregion where statistics are calculated
c ::: axialDir     =>  (const)  streamwise direction.  This must be either 
c :::                           1 or 2.
c ::: -----------------------------------------------------------
c
      subroutine FORT_VEL_RA_RND(ssDat, nSSComp, DIMS(ssdat),
     $                           nStats, nStns, stats, physStn,
     $                           vblo, vbhi, dx, probLo, probHi,
     $                           axialDir, nStations)

      implicit none

c
c     :::: Passed Variables ::::
c
      integer nSSComp, nStats, nStns, axialDir, nStations
      integer vblo(SDIM), vbhi(SDIM)
      integer DIMDEC(ssdat)
      REAL_T ssDat(DIMV(ssdat),nSSComp)
      REAL_T stats(nStns,nStats)
      REAL_T physStn(nStns)
      REAL_T dx(SDIM), probLo(SDIM), probHi(SDIM)

c
c     ::::: local variables
c
      integer i, j, k, n, r, profDir, nSSexpect, nStatsExpect
      integer flo(SDIM), fhi(SDIM)
      integer rho, rho2, rho3, rho4, 
     $        u, u2, u3, u4, rhoU, rhoU2, rhoU3, rhoU4,
     $        v, v2, v3, v4, rhoV, rhoV2, rhoV3, rhoV4,
     $        w, w2, w3, w4, rhoW, rhoW2, rhoW3, rhoW4,
     $        tr, tr2, tr3, tr4, rhoTr, rhoTr2, rhoTr3, rhoTr4, 
     $        uV, uW, vW, rhoUV, rhoUW, rhoVW, 
     $        uTr, vTr, wTr, rhoUTr, rhoVTr, rhoWTr, 
     $        p, p2, p3, p4, uP, vP, wP
      integer stnCnt(:)
      allocatable stnCnt

      REAL_T tol, x, y, z, rad, xoff, yoff, zoff, VaBar, VrBar, VtBar, 
     $       VaVrBar, VaVtBar, VrVtBar, theta, Va2Bar, Vr2Bar, Vt2Bar
      REAL_T cl(SDIM)


c
c     ------------------------------
c     ::: Define Local Constants :::
c     ------------------------------
c
      parameter (nSSexpect = 35, nStatsExpect = 10)

      call SET_LOHI(DIMS(ssdat), flo, fhi)

      do n = 1, SDIM
        cl(n) = half * (probLo(n) + probHi(n))
      enddo

      allocate(stnCnt(nStns))


c
c     -----------------------------
c     ::: Perform Sanity Checks :::
c     -----------------------------
c
      call SANITY_TEST_RND(nStns, nStats, nSSComp, profDir, 
     $                     nSSexpect, nStatsExpect, vblo, vbhi,
     $                     axialDir)


c
c     ------------------------------------------
c     ::: Define SlabStat Variable Locations :::
c     ------------------------------------------
c
      rho    = 1
      u      = 2
      rhoU   = 3
      v      = 4
      rhoV   = 5
      w      = 6
      rhoW   = 7
      tr     = 8
      rhoTr  = 9
      p      = 10
      rho2   = 11
      u2     = 12
      rhoU2  = 13
      v2     = 14
      rhoV2  = 15
      w2     = 16
      rhoW2  = 17
      tr2    = 18
      rhoTr2 = 19
      p2     = 20
      uV     = 21
      uW     = 22
      vW     = 23
      rhoUV  = 24
      rhoUW  = 25
      rhoVW  = 26
      uTr    = 27
      vTr    = 28
      wTr    = 29
      rhoUTr = 30
      rhoVTr = 31
      rhoWTr = 32
      uP     = 33
      vP     = 34
      wP     = 35


c
c     ---------------------------------------
c     ::: Calculate The Physical Location :::
c     ---------------------------------------
c
c     This logic is for round flowfields only.  For planar flowfields the
c     physical locations are easy.  But, for round flowfields we need to 
c     determine which radiuses are defined by the grid spacing of the domain.
c
      if (profDir .ne. axialDir) then
        if (axialDir .eq. 1) then
          tol = ten**(-10) * SQRT(dx(2)**2 + dx(3)**2)
          call CLC_RADIUS(vbhi(2) - vblo(2) + 1, vbhi(3) - vblo(3) + 1,
     $                    dx(2), dx(3), cl(2), cl(3), tol, 
     $                    nStns, nStations, physStn)
        else
          tol = ten**(-10) * SQRT(dx(1)**2 + dx(3)**2)
          call CLC_RADIUS(vbhi(1) - vblo(1) + 1, vbhi(3) - vblo(3) + 1,
     $                    dx(1), dx(3), cl(1), cl(3), tol, 
     $                    nStns, nStations, physStn)
        endif
      else
        tol = ten**(-10) * dx(profDir)
        nStations = vbhi(profDir) - vblo(profDir) + 1
        do n = 0, nStations-1
          physStn(n+1) = (FLOAT(n) + half) * dx(profDir)
        enddo
      endif

      do n = 1, nStations
        do i = 1, nStatsExpect
          stats(n,i) = zero
        enddo
      enddo

c
c     ----------------------------
c     ::: Calculate Statistics :::
c     ----------------------------
c
      do k = vblo(3), vbhi(3)
        z = (half + FLOAT(k)) * dx(3)

        do j = vblo(2), vbhi(2)
          y = (half + FLOAT(j)) * dx(2)

          do i = vblo(1), vbhi(1)
            x = (half + FLOAT(i)) * dx(1)

c
c           ::: Calculate radius and set profile index :::
c
            if (profDir .ne. axialDir) then
              xoff = x - cl(1)
              zoff = z - cl(3)

              rad = SQRT(xoff**2 + zoff**2)
              xoff = MIN(one, xoff/rad)
              xoff = MAX(-one, xoff)
              theta = SIGN(one, -zoff) * ACOS(xoff)

              do r = 1, nStations
                if (physStn(r)-tol .LT. rad .AND. 
     $              rad .LT. physStn(r)+tol)          n = r
              enddo

            else
              xoff = x - cl(1)
              zoff = z - cl(3)

              rad = SQRT(xoff**2 + zoff**2)
              xoff = MIN(one, xoff/rad)
              xoff = MAX(-one, xoff)
              theta = SIGN(one, -zoff) * ACOS(xoff)

              n = j + 1
            endif



c
c           ::: Calculate Statistics :::
c
            stnCnt(n) = stnCnt(n) + 1

            VaBar = ssDat(i,j,k,v)
            VrBar =   ssDat(i,j,k,u)*cos(theta) - ssDat(i,j,k,w)*sin(theta)
            VtBar = - ssDat(i,j,k,w)*cos(theta) - ssDat(i,j,k,u)*sin(theta)
            Va2Bar = ssDat(i,j,k,v2)
            Vr2Bar = ssDat(i,j,k,u2) * cos(theta)**2
     $             - two * ssDat(i,j,k,uW) * cos(theta) * sin(theta)
     $             + ssDat(i,j,k,w2) * sin(theta)**2
            Vt2Bar = ssDat(i,j,k,w2) * cos(theta)**2
     $             + two * ssDat(i,j,k,uW) * cos(theta) * sin(theta)
     $             + ssDat(i,j,k,u2) * sin(theta)**2
            VaVrBar =   ssDat(i,j,k,uV) * cos(theta) 
     $                - ssDat(i,j,k,vW) * sin(theta)
            VaVtBar = - ssDat(i,j,k,vW) * cos(theta) 
     $                - ssDat(i,j,k,uV) * sin(theta)
            VrVtBar = - ssDat(i,j,k,uW) * cos(theta)**2
     $                + (ssDat(i,j,k,w2) - ssDat(i,j,k,u2))
     $                                            * cos(theta) * sin(theta)
     $                - ssDat(i,j,k,uW) * sin(theta)**2

            stats(n,1) = stats(n,1) + VaBar
            stats(n,2) = stats(n,2) + VrBar
            stats(n,3) = stats(n,3) + VtBar
            stats(n,4) = stats(n,4) + Va2Bar - VaBar**2
            stats(n,5) = stats(n,5) + Vr2Bar - VrBar**2
            stats(n,6) = stats(n,6) + Vt2Bar - VtBar**2
            stats(n,7) = stats(n,7) + VaVrBar - VaBar * VrBar
            stats(n,8) = stats(n,8) + VaVtBar - VaBar * VtBar
            stats(n,9) = stats(n,9) + VrVtBar - VrBar * VtBar
            stats(n,10) = stats(n,10) + half * ( Va2Bar - VaBar**2
     $                                         + Vr2Bar - VrBar**2
     $                                         + Vt2Bar - VtBar**2)
          enddo
        enddo
      enddo

      do r = 1, nStations
        do n = 1, nStats
          stats(r,n) = stats(r,n) / FLOAT(stnCnt(r))
        enddo
      enddo


c
c
      return
      end



c ::: ---------------------------------------------------------------------
c ::: This routine does general sanity testing on the inputs for the
c ::: SlabStat statistics routines.  This should be called by each of the 
c ::: actual fortran routines used to calculate statistics.  This routine
c ::: also sets the profile direction, profDir.
c :::
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: nStns         => (const)  Number of stations in the statistic array
c ::: nStats        => (const)  Number of statistics to calculate
c ::: nSSComp       => (const)  Number of components in ssDat
c ::: profDir      <=           Profile direction
c ::: nSSexpect     => (const)  Number of components in ssDat
c ::: nStatsExpect  => (const)  Te number of statistics expected to calculate
c ::: vblo, vbhi    => (const)  subregion where statistics are calculated
c ::: axialDir      => (const)  streamwise direction
c ::: ---------------------------------------------------------------------
c
      subroutine SANITY_TEST_RND(nStns, nStats, nSSComp, profDir,
     $                           nSSexpect, nStatsExpect, vblo, vbhi,
     $                           axialDir)

      implicit none
c
c     ::: Passed Variables :::
c
      integer nStns, nStats, nSSComp, profDir, nSSexpect, nStatsExpect,
     $        axialDir
      integer vblo(SDIM), vbhi(SDIM)

c
c     ::: Local Variables :::
c
      integer n, profDim

c
c     -----------------------------
c     ::: Perform Sanity Checks :::
c     -----------------------------
c
      profDir = -1
      profDim = 1
      do n = 1, SDIM - 1
        if (vbhi(n) - vblo(n) + 1 .gt. 1) then
          if (profDir .ne. -1) then
            write(*,1000) vblo, vbhi
 1000       format("Error: Valid region of SlabStat data is dimensioned",
     $            /"       greater than one in more than one direction.",
     $            /"       vblo = ", SDIM(I4,1x),
     $            /"       vbhi = ", SDIM(I4,1x))
            call BL_PD_ABORT()
          endif

          profDim = vbhi(n) - vblo(n) + 1
          profDir = n
        endif
      enddo
      profDim = profDim * (vbhi(3) - vblo(3) + 1)

      if (nStns .lt. profDim) then
        write(*,1010) nStns, profDir, profDim, vblo, vbhi
 1010   format("Error: The dimensions of the statistics array are insufficient",
     $        /"       to hold the profile.",
     $        /"       nStns = ", I4, 6x, "profDir = ", I1, 
     $                                6x, "profDim = ", I4,
     $        /"       vblo = ", SDIM(I4,1x),
     $        /"       vbhi = ", SDIM(I4,1x))
        call BL_PD_ABORT()
      endif

      if (nSSComp .lt. nSSexpect) then
        write(*,1020) nSSexpect, nSSComp
 1020   format("Error: There are fewer SlabStat components than are expected",
     $         "       from the NavierStokes statistics routines.  The",
     $         "       SlabStat components defined in this routine should be",
     $         "       updated to match those saved by the code.",
     $         "       nSSexpect = ", I4, "     nSSComp = ", I4)
        call BL_PD_ABORT()
      endif

      if (nStats .ne. nStatsExpect) then
        write(*,1030) nStatsExpect, nStats
 1030   format("Error: The number of statistics passed in do not match what",
     $         "       was expected.",
     $         "       nStatsExpect = ", I4, "     nStats = ", I4)
        call BL_PD_ABORT()
      endif

      if (axialDir .ne. 2) then
        write(*,1040) axialDir
 1040   format("Error: The running statistics routines are only implemented",
     $        /"       for axialDir=2.",
     $        /"       axialDir = ", I1)
        call BL_PD_ABORT()
      endif

c
c
      return
      end




      SUBROUTINE CLC_RADIUS(dim1, dim2, dx1, dx2, cl1, cl2, tol,
     $                      nAlloc, nStns, physLoc)

c
c     *****************************
c     *** Variable Declarations ***
c     *****************************
c
      IMPLICIT NONE

c
c     *** Passed Variables ***
c
      INTEGER dim1, dim2, nAlloc, nStns
      REAL_T cl1, cl2, tol, dx1, dx2
      REAL_T physLoc(nAlloc)

c
c     *** Local Variables ***
c
      INTEGER n, m, r
      REAL_T rad, xyz1, xyz2
      LOGICAL found


c
c     **************************************************
c     *** Determine the Values of the Radius Present ***
c     **************************************************
c
c     The tolerence is used here to determine whether a radius has already 
c     been added to the physLoc array.  All of the distinct radiuses are
c     added to the list here and then it is sorted below from lowest to 
c     highest.
c
      nStns = 0
      DO n = 1, dim1
        xyz1 = (half + FLOAT(n)) * dx1
        DO m = 1, dim2
          xyz2 = (half + FLOAT(m)) * dx2
          rad = SQRT( (xyz1 - cl1)**2 + (xyz2 - cl2)**2 )

          found = .FALSE.
          DO r = 1, nStns
            IF (physLoc(r)-tol .LT. rad .AND. rad .LT. physLoc(r)+tol) THEN
              found = .TRUE.
            ENDIF
          ENDDO

          IF (.NOT.found) THEN
            nStns = nStns + 1
            physLoc(nStns) = rad
          ENDIF
        ENDDO
      ENDDO

c
c     *********************************
c     *** Sort the List of Radiuses ***
c     *********************************
c
      DO n = 1, nStns - 1
        DO m = 1, nStns - n
          IF (physLoc(m) .GT. physLoc(m+1)) THEN
            rad = physLoc(m)
            physLoc(m) = physLoc(m+1)
            physLoc(m+1) = rad
          ENDIF
        ENDDO
      ENDDO

c
c
      RETURN
      END
