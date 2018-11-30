
#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_ArrayLim.H>
#include <SLABSTAT_NS_F.H>

#define SDIM 3

module slabstat_ns_3d_module

  implicit none

  private

  public :: ns_basicstats_nctrac, ns_basicstats_ctrac


contains

!c
!c ::: -----------------------------------------------------------
!c ::: This is a general routine to calculate the basic running statistics
!c ::: based on the velocity, density, tracer and pressure.  The tracer
!c ::: passed in to this routine is not a conserved quantity.  The data
!c ::: saved by this routine are sufficient to calculate a core set of
!c ::: statistics on these fields using both Reynolds and Favre averages.
!c ::: The state should be passed in to this routine in the order,
!c :::   Rho, U, V, W, Tr, P
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: runStats      <=  Array to contain the running statistics
!c ::: DIMS(runStats) => Dimensions of the running statistics array
!c ::: nStats         => Number of components in the statistics array
!c ::: state          =) Array containing the state data
!c ::: DIMS(state)    => Dimensions of the state data array
!c ::: nStateComp     => Number of components in the state array
!c ::: dt             => Time step
!c ::: delta          => Grid spacing
!c ::: -----------------------------------------------------------
!c
      subroutine ns_basicstats_nctrac(state, DIMS(state), nStateComp, &
                                      runStats, DIMS(runStats), nStats, &
                                      dt, delta) &
                                      bind(C,name="ns_basicstats_nctrac")

      implicit none

!c
!c     ::::: passed variables
!c
      integer nStats, nStateComp
      integer DIMDEC(state)
      integer DIMDEC(runStats)
      REAL_T dt
      REAL_T delta(SDIM)
      REAL_T runStats(DIMV(runStats),nStats)
      REAL_T state(DIMV(state),nStateComp)

!c
!c     ::::: local variables
!c
      integer i, j, k, n, Rho, xVel, yVel, zVel, Trac, Pres, &
             nStateExpect, nStatsExpect, nGhostExpect, num
      integer loStats(SDIM), hiStats(SDIM)
      integer loState(SDIM), hiState(SDIM)

      PARAMETER (nStateExpect = 6, nStatsExpect = 35, nGhostExpect = 0)

!c
!c     ===================================
!c     ::: Set the Values of LO and HI :::
!c     ===================================
!c
      call SET_LOHI(DIMS(runStats), loStats, hiStats)
      call SET_LOHI(DIMS(state), loState, hiState)

#ifndef NDEBUG
      if (nStats .NE. nStatsExpect .or. nStateComp .ne. nStateExpect) then
        write(*,1000) nStatsExpect, nStateExpect, nStats, nStateComp
 1000   format('Incorrect number of statistic and/or state components', &
             /'FORT_NS_BASICSTATS_NCTRAC.  Should have nStats = ',I5,'and', &
             /'nStateComp = ',I5, &
             /'   nStats = ', I5, 5x, 'nStateComp = ', I5)
        call bl_abort(" ")
      endif

      do n = 1, SDIM
        if (loState(n) .GT. loStats(n) - nGhostExpect .OR. &
           hiState(n) .LT. hiStats(n) + nGhostExpect) then
          write(*,1010) n, nGhostExpect, loStats, hiStats, loState, hiState
 1010     format('Incorrect number of ghost cells in the state date in', &
               /'FORT_NS_BASICSTATS_NCTRAC.', &
               /'   Direction = ', I2, 5x, 'nGhostExpect = ', I2, &
               /'   loStats = ', SDIM (I2,1x), 5x, 'hiStats = ', SDIM (I2,1x), &
               /'   loState = ', SDIM (I2,1x), 5x, 'hiState = ', SDIM (I2,1x))
          call bl_abort(" ")
        endif
      enddo
#endif


!c
!c     =========================
!c     ::: Set State Indices :::
!c     =========================
!c
      Rho  = 1
      xVel = 2
      yVel = 3
      zVel = 4
      Trac = 5
      Pres = 6

!c
!c     ====================================
!c     ::: Calculate Running Statistics :::
!c     ====================================
!c
      do k = loStats(3), hiStats(3)
        do j = loStats(2), hiStats(2)
          do i = loStats(1), hiStats(1)
            num = 1
            do n = 1, 2
              runStats(i,j,k,num) = runStats(i,j,k,num) &
                                                   + dt * state(i,j,k,Rho)**n
              runStats(i,j,k,num+1) = runStats(i,j,k,num+1) &
                                                   + dt * state(i,j,k,xVel)**n
              runStats(i,j,k,num+2) = runStats(i,j,k,num+2) &
                                + dt * state(i,j,k,Rho) * state(i,j,k,xVel)**n
              runStats(i,j,k,num+3) = runStats(i,j,k,num+3) &
                                                   + dt * state(i,j,k,yVel)**n
              runStats(i,j,k,num+4) = runStats(i,j,k,num+4) &
                                + dt * state(i,j,k,Rho) * state(i,j,k,yVel)**n
              runStats(i,j,k,num+5) = runStats(i,j,k,num+5) &
                                                   + dt * state(i,j,k,zVel)**n
              runStats(i,j,k,num+6) = runStats(i,j,k,num+6) &
                                + dt * state(i,j,k,Rho) * state(i,j,k,zVel)**n
              runStats(i,j,k,num+7) = runStats(i,j,k,num+7) &
                                                   + dt * state(i,j,k,Trac)**n
              runStats(i,j,k,num+8) = runStats(i,j,k,num+8) &
                                + dt * state(i,j,k,Rho) * state(i,j,k,Trac)**n
              runStats(i,j,k,num+9) = runStats(i,j,k,num+9) &
                                                   + dt * state(i,j,k,Pres)**n
              num = num + 10
            enddo

            runStats(i,j,k,num) = runStats(i,j,k,num) + dt &
                                       * state(i,j,k,xVel) * state(i,j,k,yVel)
            runStats(i,j,k,num+1) = runStats(i,j,k,num+1) + dt &
                                       * state(i,j,k,xVel) * state(i,j,k,zVel)
            runStats(i,j,k,num+2) = runStats(i,j,k,num+2) + dt &
                                       * state(i,j,k,yVel) * state(i,j,k,zVel)
            num = num + 3

            runStats(i,j,k,num) = runStats(i,j,k,num) + dt * state(i,j,k,Rho) &
                                       * state(i,j,k,xVel) * state(i,j,k,yVel)
            runStats(i,j,k,num+1) = runStats(i,j,k,num+1) + dt *state(i,j,k,Rho) &
                                       * state(i,j,k,xVel) * state(i,j,k,zVel)
            runStats(i,j,k,num+2) = runStats(i,j,k,num+2) + dt *state(i,j,k,Rho) &
                                       * state(i,j,k,yVel) * state(i,j,k,zVel)
            num = num + 3

            runStats(i,j,k,num) = runStats(i,j,k,num) &
                              + dt * state(i,j,k,xVel) * state(i,j,k,Trac)
            runStats(i,j,k,num+1) = runStats(i,j,k,num+1) &
                              + dt * state(i,j,k,yVel) * state(i,j,k,Trac)
            runStats(i,j,k,num+2) = runStats(i,j,k,num+2) &
                              + dt * state(i,j,k,zVel) * state(i,j,k,Trac)
            num = num + 3

            runStats(i,j,k,num) = runStats(i,j,k,num) & 
                              + dt * state(i,j,k,Rho) * state(i,j,k,xVel) &
                                                      * state(i,j,k,Trac)
            runStats(i,j,k,num+1) = runStats(i,j,k,num+1) &
                              + dt * state(i,j,k,Rho) * state(i,j,k,yVel) &
                                                      * state(i,j,k,Trac)
            runStats(i,j,k,num+2) = runStats(i,j,k,num+2) &
                              + dt * state(i,j,k,Rho) * state(i,j,k,zVel) &
                                                      * state(i,j,k,Trac)
            num = num + 3

            runStats(i,j,k,num) = runStats(i,j,k,num) &
                              + dt * state(i,j,k,xVel) * state(i,j,k,Pres)
            runStats(i,j,k,num+1) = runStats(i,j,k,num+1) &
                              + dt * state(i,j,k,yVel) * state(i,j,k,Pres)
            runStats(i,j,k,num+2) = runStats(i,j,k,num+2) &
                              + dt * state(i,j,k,zVel) * state(i,j,k,Pres)
            num = num + 3
          enddo
        enddo
      enddo

!c
!c
      return
      end subroutine ns_basicstats_nctrac



!c
!c ::: -----------------------------------------------------------
!c ::: This is a general routine to calculate the basic running statistics
!c ::: based on the velocity, density, tracer and pressure.  The tracer
!c ::: passed in to this routine is a conserved quantity.  The data
!c ::: saved by this routine are sufficient to calculate a core set of
!c ::: statistics on these fields using both Reynolds and Favre averages.
!c ::: The state should be passed in to this routine in the order,
!c :::   Rho, U, V, W, Tr, P
!c :::
!c ::: INPUTS/OUTPUTS:
!c :::
!c ::: runStats      <=  Array to contain the running statistics
!c ::: DIMS(runStats) => Dimensions of the running statistics array
!c ::: nStats         => Number of components in the statistics array
!c ::: state          =) Array containing the state data
!c ::: DIMS(state)    => Dimensions of the state data array
!c ::: nStateComp     => Number of components in the state array
!c ::: dt             => Time step
!c ::: delta          => Grid spacing
!c ::: -----------------------------------------------------------
!c
      subroutine ns_basicstats_ctrac (state, DIMS(state), nStateComp, &
                  runStats, DIMS(runStats), nStats, &
                  dt, delta) &
                  bind(C,name="ns_basicstats_ctrac")

      implicit none

!c
!c     ::::: passed variables
!c
      integer nStats, nStateComp
      integer DIMDEC(state)
      integer DIMDEC(runStats)
      REAL_T dt
      REAL_T delta(SDIM)
      REAL_T runStats(DIMV(runStats),nStats)
      REAL_T state(DIMV(state),nStateComp)

!c
!c     ::::: local variables
!c
      integer i, j, k, n, Rho, xVel, yVel, zVel, RhoTrac, Pres, &
             nStateExpect, nStatsExpect, nGhostExpect, num
      integer loStats(SDIM), hiStats(SDIM)
      integer loState(SDIM), hiState(SDIM)
      REAL_T  tracer

      PARAMETER (nStateExpect = 6, nStatsExpect = 35, nGhostExpect = 0)

!c
!c     ===================================
!c     ::: Set the Values of LO and HI :::
!c     ===================================
!c
      call SET_LOHI(DIMS(runStats), loStats, hiStats)
      call SET_LOHI(DIMS(state), loState, hiState)

#ifndef NDEBUG
      if (nStats .NE. nStatsExpect .or. nStateComp .ne. nStateExpect) then
        write(*,1000) nStatsExpect, nStateExpect, nStats, nStateComp
 1000   format('Incorrect number of statistic and/or state components', &
             /'FORT_NS_BASICSTATS_NCTRAC.  Should have nStats = ',I5,'and', &
             /'nStateComp = ',I5, &
             /'   nStats = ', I5, 5x, 'nStateComp = ', I5)
        call bl_abort(" ")
      endif

      do n = 1, SDIM
        if (loState(n) .GT. loStats(n) - nGhostExpect .OR. &
           hiState(n) .LT. hiStats(n) + nGhostExpect) then
          write(*,1010) n, nGhostExpect, loStats, hiStats, loState, hiState
 1010     format('Incorrect number of ghost cells in the state date in', &
               /'FORT_NS_BASICSTATS_NCTRAC.', &
               /'   Direction = ', I2, 5x, 'nGhostExpect = ', I2, &
               /'   loStats = ', SDIM (I2,1x), 5x, 'hiStats = ', SDIM (I2,1x), &
               /'   loState = ', SDIM (I2,1x), 5x, 'hiState = ', SDIM (I2,1x))
          call bl_abort(" ")
        endif
      enddo
#endif


!c
!c     =========================
!c     ::: Set State Indices :::
!c     =========================
!c
      Rho  = 1
      xVel = 2
      yVel = 3
      zVel = 4
      RhoTrac = 5
      Pres = 6

!c
!c     ====================================
!c     ::: Calculate Running Statistics :::
!c     ====================================
!c
      do k = loStats(3), hiStats(3)
        do j = loStats(2), hiStats(2)
          do i = loStats(1), hiStats(1)
            num = 1

            tracer = state(i,j,k,RhoTrac) / state(i,j,k,Rho)
            do n = 1, 2
              runStats(i,j,k,num) = runStats(i,j,k,num) + dt &
                                                        * state(i,j,k,Rho)**n
              runStats(i,j,k,num+1) = runStats(i,j,k,num+1) + dt &
                                                        * state(i,j,k,xVel)**n
              runStats(i,j,k,num+2) = runStats(i,j,k,num+2) + dt &
                                     * state(i,j,k,Rho) * state(i,j,k,xVel)**n
              runStats(i,j,k,num+3) = runStats(i,j,k,num+3) + dt &
                                                        * state(i,j,k,yVel)**n
              runStats(i,j,k,num+4) = runStats(i,j,k,num+4) + dt &
                                     * state(i,j,k,Rho) * state(i,j,k,yVel)**n
              runStats(i,j,k,num+5) = runStats(i,j,k,num+5) + dt &
                                                        * state(i,j,k,zVel)**n
              runStats(i,j,k,num+6) = runStats(i,j,k,num+6) + dt &
                                     * state(i,j,k,Rho) * state(i,j,k,zVel)**n
              runStats(i,j,k,num+7) = runStats(i,j,k,num+7) + dt * tracer**n
              runStats(i,j,k,num+8) = runStats(i,j,k,num+8) + dt &
                                                * state(i,j,k,Rho) * tracer**n
              runStats(i,j,k,num+9) = runStats(i,j,k,num+9) + dt &
                                                        * state(i,j,k,Pres)**n
              num = num + 10
            enddo

            runStats(i,j,k,num) = runStats(i,j,k,num) + dt &
                                       * state(i,j,k,xVel) * state(i,j,k,yVel)
            runStats(i,j,k,num+1) = runStats(i,j,k,num+1) + dt &
                                       * state(i,j,k,xVel) * state(i,j,k,zVel)
            runStats(i,j,k,num+2) = runStats(i,j,k,num+2) + dt &
                                       * state(i,j,k,yVel) * state(i,j,k,zVel)
            num = num + 3

            runStats(i,j,k,num) = runStats(i,j,k,num) + dt * state(i,j,k,Rho) &
                                       * state(i,j,k,xVel) * state(i,j,k,yVel)
            runStats(i,j,k,num+1) = runStats(i,j,k,num+1) + dt *state(i,j,k,Rho) &
                                       * state(i,j,k,xVel) * state(i,j,k,zVel)
            runStats(i,j,k,num+2) = runStats(i,j,k,num+2) + dt *state(i,j,k,Rho) &
                                       * state(i,j,k,yVel) * state(i,j,k,zVel)
            num = num + 3

            runStats(i,j,k,num) = runStats(i,j,k,num) + dt &
                                     * state(i,j,k,xVel) * tracer
            runStats(i,j,k,num+1) = runStats(i,j,k,num+1) + dt &
                                     * state(i,j,k,yVel) * tracer
            runStats(i,j,k,num+2) = runStats(i,j,k,num+2) + dt &
                                     * state(i,j,k,zVel) * tracer
            num = num + 3

            runStats(i,j,k,num) = runStats(i,j,k,num) + dt &
                                     * state(i,j,k,xVel) * state(i,j,k,RhoTrac)
            runStats(i,j,k,num+1) = runStats(i,j,k,num+1) + dt &
                                     * state(i,j,k,yVel) * state(i,j,k,RhoTrac)
            runStats(i,j,k,num+2) = runStats(i,j,k,num+2) + dt &
                                     * state(i,j,k,zVel) * state(i,j,k,RhoTrac)
            num = num + 3

            runStats(i,j,k,num) = runStats(i,j,k,num) + dt &
                                     * state(i,j,k,xVel) * state(i,j,k,Pres)
            runStats(i,j,k,num+1) = runStats(i,j,k,num+1) + dt &
                                     * state(i,j,k,yVel) * state(i,j,k,Pres)
            runStats(i,j,k,num+2) = runStats(i,j,k,num+2) + dt &
                                     * state(i,j,k,zVel) * state(i,j,k,Pres)
            num = num + 3
          enddo
        enddo
      enddo

!c
!c
      return
      end subroutine ns_basicstats_ctrac

  end module slabstat_ns_3d_module

