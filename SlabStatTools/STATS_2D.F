
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_ArrayLim.H>

#include <StatTypes.H>
#define SDIM 2


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
      integer i, j, n, profDir, nSSexpect, nStatsExpect
      integer flo(SDIM), fhi(SDIM)
      integer rho, rho2,
     $        u, u2, rhoU, rhoU2,
     $        v, v2, rhoV, rhoV2,
     $        tr, tr2, rhoTr, rhoTr2,
     $        uV, rhoUV, uTr, vTr, rhoUTr, rhoVTr, p, p2, uP, vP


c
c     ------------------------------
c     ::: Define Local Constants :::
c     ------------------------------
c
      parameter (nSSexpect = 24, nStatsExpect = 6)

      call SET_LOHI(DIMS(ssdat), flo, fhi)


c
c     -----------------------------
c     ::: Perform Sanity Checks :::
c     -----------------------------
c
      call SANITY_TEST(nStns, nStats, nSSComp, profDir, 
     $                 nSSexpect, nStatsExpect, vblo, vbhi)

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
      tr     = 6
      rhoTr  = 7
      p      = 8
      rho2   = 9
      u2     = 10
      rhoU2  = 11
      v2     = 12
      rhoV2  = 13
      tr2    = 14
      rhoTr2 = 15
      p2     = 16
      uV     = 17
      rhoUV  = 18
      uTr    = 19
      vTr    = 20
      rhoUTr = 21
      rhoVTr = 22
      uP     = 23
      vP     = 24


c
c     ----------------------------
c     ::: Calculate Statistics :::
c     ----------------------------
c
      nStations = vbhi(profDir) - vblo(profDir) + 1
      do n = 1, vbhi(profDir) - vblo(profDir) + 1
        physStn(n) = (FLOAT(n) + half) * dx(profDir)

        if (profDir .eq. 1) then
          i = n - 1 + vblo(profDir)
          j = vblo(2)
        else
          i = vblo(1)
          j = n - 1 + vblo(profDir)
        endif

        stats(n,1) = ssDat(i,j,u)
        stats(n,2) = ssDat(i,j,v)
        stats(n,3) = ssDat(i,j,u2) - ssDat(i,j,u)**2
        stats(n,4) = ssDat(i,j,v2) - ssDat(i,j,v)**2
        stats(n,5) = ssDat(i,j,uV) - ssDat(i,j,u) * ssDat(i,j,v)
        stats(n,6) = half * ( stats(n,3) + stats(n,4) )
      enddo


c
c
      return
      end


c ::: -----------------------------------------------------------
c ::: This routine calculates the velocity statistics using Favre
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
c ::: -----------------------------------------------------------
c
      subroutine FORT_VEL_FA_RND(ssDat, nSSComp, DIMS(ssdat),
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
      integer i, j, n, profDir, nSSexpect, nStatsExpect
      integer flo(SDIM), fhi(SDIM)
      integer rho, rho2,
     $        u, u2, rhoU, rhoU2,
     $        v, v2, rhoV, rhoV2,
     $        tr, tr2, rhoTr, rhoTr2,
     $        uV, rhoUV, uTr, vTr, rhoUTr, rhoVTr, p, p2, uP, vP


c
c     ------------------------------
c     ::: Define Local Constants :::
c     ------------------------------
c
      parameter (nSSexpect = 24, nStatsExpect = 6)

      call SET_LOHI(DIMS(ssdat), flo, fhi)


c
c     -----------------------------
c     ::: Perform Sanity Checks :::
c     -----------------------------
c
      call SANITY_TEST(nStns, nStats, nSSComp, profDir, 
     $                 nSSexpect, nStatsExpect, vblo, vbhi)

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
      tr     = 6
      rhoTr  = 7
      p      = 8
      rho2   = 9
      u2     = 10
      rhoU2  = 11
      v2     = 12
      rhoV2  = 13
      tr2    = 14
      rhoTr2 = 15
      p2     = 16
      uV     = 17
      rhoUV  = 18
      uTr    = 19
      vTr    = 20
      rhoUTr = 21
      rhoVTr = 22
      uP     = 23
      vP     = 24


c
c     ----------------------------
c     ::: Calculate Statistics :::
c     ----------------------------
c
      nStations = vbhi(profDir) - vblo(profDir) + 1
      do n = 1, vbhi(profDir) - vblo(profDir) + 1
        physStn(n) = (FLOAT(n) + half) * dx(profDir)

        if (profDir .eq. 1) then
          i = n - 1 + vblo(profDir)
          j = vblo(2)
        else
          i = vblo(1)
          j = n - 1 + vblo(profDir)
        endif

        stats(n,1) = ssDat(i,j,rhoU) / ssDat(i,j,rho)
        stats(n,2) = ssDat(i,j,rhoV) / ssDat(i,j,rho)
        stats(n,3) = ssDat(i,j,rhoU2) / ssDat(i,j,rho)
     $                               - ssDat(i,j,rhoU)**2 / ssDat(i,j,rho)**2
        stats(n,4) = ssDat(i,j,rhoV2) / ssDat(i,j,rho)
     $                               - ssDat(i,j,rhoV)**2 / ssDat(i,j,rho)**2
        stats(n,5) = ssDat(i,j,rhoUV) / ssDat(i,j,rho)
     $                - ssDat(i,j,rhoU) * ssDat(i,j,rhoV) / ssDat(i,j,rho)**2
        stats(n,6) = half * ( stats(n,3) + stats(n,4) )
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
c ::: ---------------------------------------------------------------------
c
      subroutine SANITY_TEST(nStns, nStats, nSSComp, profDir,
     $                       nSSexpect, nStatsExpect, vblo, vbhi)

      implicit none
c
c     ::: Passed Variables :::
c
      integer nStns, nStats, nSSComp, profDir, nSSexpect, nStatsExpect
      integer vblo(SDIM), vbhi(SDIM)

c
c     ::: Local Variables :::
c
      integer n

c
c     -----------------------------
c     ::: Perform Sanity Checks :::
c     -----------------------------
c
      profDir = -1
      do n = 1, SDIM
        if (vbhi(n) - vblo(n) + 1 .gt. 1) then
          if (profDir .ne. -1) then
            write(*,1000) vblo, vbhi
 1000       format("Error: Valid region of SlabStat data is dimensioned",
     $            /"       greater than one in more than one direction.",
     $            /"       vblo = ", SDIM(I4,1x),
     $            /"       vbhi = ", SDIM(I4,1x))
            call BL_PD_ABORT()
          endif

          profDir = n
        endif
      enddo

      if (nStns .lt. vbhi(profDir)-vblo(profDir)) then
        write(*,1010) nStns, profDir, vblo, vbhi
 1010   format("Error: The dimensions of the statistics array are insufficient",
     $        /"       to hold the profile.",
     $        /"       nStns = ", I4, "      profDir = ", I1,
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

c
c
      return
      end
