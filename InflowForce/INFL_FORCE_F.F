c
c $Id: INFL_FORCE_F.F,v 1.2 2000-06-23 18:24:29 sstanley Exp $
c

#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "ArrayLim.H"
#include "infl_frc.H"
#include "FLUCTFILE.H"
#include "BC_TYPES.H"

#define SDIM BL_SPACEDIM

#define FF_UNIT       20

#if defined(BL_T3E)
#define OPEN_APPEND position='append'
#elif defined(BL_AIX)
#define OPEN_APPEND position='append'
#else
#define OPEN_APPEND access='append'
#endif

c ::: -----------------------------------------------------------
c ::: This routine is used by INITDATA and the fill routines to 
c ::: extrapolate the perturbations from the flct_file to 
c ::: fill the data required for forcing the inflow.  Mostly this 
c ::: routine manages the reading of the data from the flct_file
c ::: and then passes off to the XTR_DAT routine to actually
c ::: extrapolate the data and fill the arrays.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: fillComp       =>  Component to fill
c ::: DIMS(inflDat)  =>  Dimensions of inflDat
c ::: inflDat       <=   Array to fill
c ::: dx             =>  Grid spacing
c ::: time           =>  Time for the fill
c ::: -----------------------------------------------------------
c
      subroutine INFL_FILL(fillComp, DIMS(inflDat), inflDat, xlo, dx, time,
     $                     bc, probLo, probHi)
      
c
c     :::: Passed Variables ::::
c
      implicit none
      integer fillComp
      integer DIMDEC(inflDat)
      integer bc(SDIM,2)

      REAL_T time
      REAL_T xlo(SDIM), dx(SDIM)
      REAL_T inflDat(DIMV(inflDat))
      REAL_T probLo(SDIM), probHi(SDIM)

c
c     ::::: local variables
c
      integer lo(SDIM), hi(SDIM)
      integer n, npass, filePnt, storePnt, filePntOld, filePntReqd,
     $        maxFilePntChng, npassOld
      REAL_T dtFile, timeMaxFile, timeOffset, timePnt,
     $       dtEstimate
      logical is_periodic

c
c     The arrays below are dimensioned as 3-d no matter what the 
c     BL_SPACEDIM is.  This is to allow for the fact that the inflow forcing
c     data arrays are always 3-d.
c
      integer dimFile(3), loStoreDim(3), hiStoreDim(3)
      integer FF_DIMDEC(storDat)
      REAL_T dxFile(3), probSizeFile(3), xloFile(3), xhiFile(3)
      REAL_T storDat(:,:,:,:)
      allocatable storDat
      save loStoreDim, hiStoreDim, storDat,
     $     dimFile, xloFile, xhiFile, dxFile, dtFile, timeMaxFile, 
     $     storePnt, filePnt, maxFilePntChng, timePnt, timeOffset, npass

c
c     ::::: common blocks
c
#include "INFL_FORCE_F.H"


c
c     --------------------------------------------------------
c     :::: Allocate and Fill the Arrays on the First Call ::::
c     --------------------------------------------------------
c
      if (.NOT. ALLOCATED(storDat)) then
        open(FF_UNIT, file=flct_file, form='unformatted')
        call RD_FLCTHD(FF_UNIT, dimFile, probSizeFile, dxFile)
        close(FF_UNIT)

c 
c  do some size checking
c
        call SET_LOHI(DIMS(inflDat), lo, hi)

        do n= 1,SDIM
          if (n .NE. strmwse_dir) then
            is_periodic = bc(n,1) .eq. INT_DIR .and. 
     $                    xlo(n) + dx(n)*half .lt. probLo(n)
            is_periodic = is_periodic .or. (bc(n,2) .eq. INT_DIR .and.
     $             xlo(n) + dx(n)*(float(hi(n)-lo(n)) + half) .gt. probHi(n))

            if (.NOT.is_periodic .AND.
     $          (probHi(n) - probLo(n) .lt. probSizeFile(n))) then
              print *, "turbFile is too big in dir ",n,probSizeFile(n),
     $                     probHi(n)-probLo(n)
              call BL_PD_ABORT()
            endif

            if (is_periodic .and. 
     $          probHi(n)-probLo(n) - (probSizeFile(n)-two*dxFile(n)) 
     $                                             .gt. 0.001d0*dxFile(n)) then
              print *, "size wrong in periodic dir",n,
     $                      probSizeFile(n)-two*dxFile(n), probHi(n)-probLo(n)
              call BL_PD_ABORT()
            endif
          endif
        enddo

        do n = 1, SDIM
          xloFile(n) = half * (probHi(n) - probLo(n)) - half * probSizeFile(n)
          xhiFile(n) = xloFile(n) + probSizeFile(n)
          if (n .ne. strmwse_dir) then
             write(6,1001) n,probLo(n),probHi(n),xloFile(n),xhiFile(n)
          endif
        enddo
 1001   format('dir: ',i2,' problem: ',f10.5,2x,f10.5,' file:',f10.5,2x,f10.5)
#if (BL_SPACEDIM==2)
        xloFile(3) = zero
        xhiFile(3) = zero
#endif

c
c       :::: Set the Dimensions for the inflow data array and allocate space
c
        do n = 1, 3
          loStoreDim(n) = 1
          hiStoreDim(n) = dimFile(n)
        enddo
        if (numInflPlanesStore .GT. 0 .AND.
     $      numInflPlanesStore .LT. hiStoreDim(strmwse_dir))
     $    hiStoreDim(strmwse_dir) = numInflPlanesStore

        call FF_SET_ARGS(FF_DIMS(storDat), loStoreDim, hiStoreDim)

        ALLOCATE(storDat(FF_DIMV(storDat),COMP_FLCTFILE))
        

c
c       ::::   Convert the streamwise direction lengths to times   ::::
c       :::: and set up the pointers into the data arrays and file ::::
c
        timeMaxFile = probSizeFile(strmwse_dir) / convVel
        dtFile = dxFile(strmwse_dir) / convVel

        npass = INT(time / timeMaxFile)
        timeOffset = FLOAT(npass) * timeMaxFile

        filePnt = INT((time - timeOffset) / dtFile) + 1
        timePnt = timeOffset + FLOAT(filePnt-1) * dtFile
        storePnt = 2
        maxFilePntChng = 0

c
c       ::: Fill the data Arrays :::
c
        call FILL_FRCARRYS(filePnt - 1, 1, dimFile, COMP_FLCTFILE,
     $                     FF_DIMS(storDat), storDat)
      endif

c
c     ------------------------------------------
c     :::: Find Correct Interpolation Point ::::
c     ------------------------------------------
c
      call FF_SET_ARGS(FF_DIMS(storDat), loStoreDim, hiStoreDim)

#ifdef INFL_FRC_DIAGS
      OPEN(12,file='LOG.Inflow_Forcing',OPEN_APPEND)
      if (fillComp .eq. FLCT_XVEL) then
        WRITE(12,2000) 'X',time, storePnt, filePnt, npass, timeOffset,
     $                 dimFile(strmwse_dir), numInflPlanesStore, maxFilePntChng
      else if (fillComp .eq. FLCT_YVEL) then
        WRITE(12,2000) 'Y',time, storePnt, filePnt, npass, timeOffset,
     $                 dimFile(strmwse_dir), numInflPlanesStore, maxFilePntChng
      else if (fillComp .eq. FLCT_ZVEL) then
        WRITE(12,2000) 'Z',time, storePnt, filePnt, npass, timeOffset,
     $                 dimFile(strmwse_dir), numInflPlanesStore, maxFilePntChng
      endif
      close(12)
#endif
      filePntOld = filePnt
      npassOld = npass
      do while ( (time .lt. timePnt - half*dtFile) .or.
     $           (time .gt. timePnt + half*dtFile) )
        if (time .gt. timePnt + half*dtFile) then
          filePnt = filePnt + 1
          storePnt = storePnt + 1

          if (filePnt .gt. dimFile(strmwse_dir) - 1) then
            npass = npass + 1
            timeOffset = FLOAT(npass) * timeMaxFile

            filePnt = filePnt - dimFile(strmwse_dir) + 1
          endif

        else if (time .lt. timePnt - half*dtFile) then
          filePnt = filePnt - 1
          storePnt = storePnt - 1

          if (filePnt .lt. 1) then
            npass = npass - 1
            timeOffset = FLOAT(npass) * timeMaxFile

            filePnt = filePnt + dimFile(strmwse_dir) - 1
          endif
        endif

        timePnt = timeOffset + FLOAT(filePnt-1) * dtFile
      enddo

      if (npass .eq. npassOld) 
     $  maxFilePntChng = MAX(maxFilePntChng, ABS(filePnt-filePntOld))

#ifdef INFL_FRC_DIAGS
      OPEN(12,file='LOG.Inflow_Forcing',OPEN_APPEND)
      if (fillComp .eq. FLCT_XVEL) then
        WRITE(12,2000) 'X',time, storePnt, filePnt, npass, timeOffset,
     $                 dimFile(strmwse_dir), numInflPlanesStore, maxFilePntChng
      else if (fillComp .eq. FLCT_YVEL) then
        WRITE(12,2000) 'Y',time, storePnt, filePnt, npass, timeOffset,
     $                 dimFile(strmwse_dir), numInflPlanesStore, maxFilePntChng
      else if (fillComp .eq. FLCT_ZVEL) then
        WRITE(12,2000) 'Z',time, storePnt, filePnt, npass, timeOffset,
     $                 dimFile(strmwse_dir), numInflPlanesStore, maxFilePntChng
      endif
 2000 FORMAT(A1,1xE11.3,1x,3(I3,1x),E11.3,1x,3(I3,1x))
      close(12)
#endif


c
c     -------------------------------------------
c     ::: Load Data if off end of stored data :::
c     -------------------------------------------
c
c     *** Stepped off the right edge of the stored data ***
c
      if (storePnt + 1 .gt. hiStoreDim(strmwse_dir)) then
        filePntReqd = filePnt - maxFilePntChng

        call FILL_FRCARRYS(filePntReqd - 1, 1, dimFile, COMP_FLCTFILE,
     $                     FF_DIMS(storDat), storDat)
        storePnt = 2 - filePntReqd + filePnt

#ifdef INFL_FRC_DIAGS
        OPEN(12,file='LOG.Inflow_Forcing',OPEN_APPEND)
        WRITE(12,2010) storePnt, filePnt, filePntReqd
        close(12)
 2010 FORMAT(12x,3(I3,1x))
#endif

c
c     *** Stepped off the left edge of the stored data ***
c
c     In this case we flush the data array and load it entirely from scratch.
c     I really should fix the shifting routine to allow the data to be right 
c     shifted as well as left shifted. SAS
c
      else if (storePnt - 1 .lt. 1) then
        open(11,file='WARNING.Inflow_Forcing',OPEN_APPEND)
        write(11,1000) storePnt, filePnt, timePnt, npass, timeOffset,
     $                 dtFile, dimFile(strmwse_dir), time
 1000   format('** Stepped off the left end of the stored data.  Reloading',
     $        /3x,'the entire data array.',
     $        /3x,'storePnt = ',I3,2x,'filePnt = ',I3,2x,'timePnt = ',E12.4,
     $        /3x,'npass = ',I3,2x,'timeOffset = ',E12.4,2x,'dtFile = ',E12.4,
     $        /3x,'dimFile = ',I3,2x,'time = ',E12.4)
        close(11)
        call FILL_FRCARRYS(filePnt - 1, 1, dimFile, COMP_FLCTFILE,
     $                     FF_DIMS(storDat), storDat)
        storePnt = 2

#ifdef INFL_FRC_DIAGS
        OPEN(12,file='LOG.Inflow_Forcing',OPEN_APPEND)
        WRITE(12,2020) storePnt, filePnt
        close(12)
 2020 FORMAT('BACK LOAD',3x,2(I3,1x))
#endif
      endif


c
c     -----------------------------------------
c     :::: Interpolate Data to Fill inflDat :::
c     -----------------------------------------
c
      call INTRP_DATA(time, xlo, dx, storePnt, fillComp, timePnt, dtFile,
     $                dxFile, xloFile, xhiFile, COMP_FLCTFILE, 
     $                FF_DIMS(storDat), storDat, DIMS(inflDat), inflDat,
     $                bc, probLo, probHi)

c
c
      return
      end

c
c ::: -----------------------------------------------------------
c ::: This routine fills the inflow forcing data array from the file.
c ::: A basepoint is specified for the array as well as for the file.
c ::: These are the points in the strmwse_dir at which reading from the
c ::: file is started and at which the array is filled from.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: baseFilePnt   => Basepoint in the file in the strmwse_dir to 
c :::                    start reading from.
c ::: baseArrayPnt  => Basepoint in the array in the strmwse_dir to 
c :::                    start filling from.  The array is filled from
c :::                    this point in strmwse_dir to the end of the 
c :::                    array.
c ::: dimFile       => Dimensions from the header of the fluctuations
c ::: nComp         => Number of components in the array
c ::: FF_DIMS(dat)  => Dimensions of the array dat
c ::: dat          <=  Array to fill
c ::: -----------------------------------------------------------
c
      subroutine FILL_FRCARRYS(baseFilePnt, baseArrayPnt, dimFile, 
     $                         nComp, FF_DIMS(dat), dat)

c
c     *****************************
c     *** Variable Declarations ***
c     *****************************
c
      implicit none

c
c     *** Passed Variables ***
c
      integer baseFilePnt, baseArrayPnt, nComp, dimFile(3)
      integer FF_DIMDEC(dat)
      REAL_T dat(FF_DIMV(dat),nComp)

c
c     *** Local Variables ***
c
      integer intDummy, n
      integer lo(3), hi(3), loRd(3), hiRd(3), filLo(3)
      REAL_T realDummy

c
c     ::: Common Blocks :::
c
#include "INFL_FORCE_F.H"

c
c     *** Define the LO and HI arrays ***
c
      call FF_SET_LOHI(FF_DIMS(dat), lo, hi)

c
c     ------------------------------------------------------
c     *** Wrap baseFilePnt into the box if it is outside ***
c     ------------------------------------------------------
c
c     Note: This assumes the first and last point in the file 
c           are the same.  This is the case for data generated with mkInitFlct.
c
      if (baseFilePnt .lt. 1) then
        baseFilePnt = dimFile(strmwse_dir) - (1 - baseFilePnt)
      endif

      if (baseFilePnt .gt. dimFile(strmwse_dir) - 1) then
        baseFilePnt = baseFilePnt - (dimFile(strmwse_dir) - 1)
      endif

c
c     ----------------------------------------------------
c     *** If the data can be filled in one pass, do so ***
c     ----------------------------------------------------
c
      if (hi(strmwse_dir) - baseArrayPnt + 1 
     $                        .LE. dimFile(strmwse_dir) - baseFilePnt + 1) then

        do n = 1, 3
          filLo(n) = 1
          loRd(n) = lo(n)
          hiRd(n) = hi(n)
        enddo
        loRd(strmwse_dir) = baseArrayPnt
        filLo(strmwse_dir) = baseFilePnt

        open(FF_UNIT, file=flct_file, form='unformatted')
        read(FF_UNIT) intDummy
        read(FF_UNIT) realDummy
        read(FF_UNIT) intDummy
        do n = 1, nComp
          call RD_FLCTREC(FF_UNIT, dimFile, loRd, hiRd, filLo, FF_DIMS(dat),
     $                    dat(lo(1),lo(2),lo(3),n))
        enddo
        close(FF_UNIT)


c
c     --------------------------------------------------------------
c     *** If the data can not be filled in one pass, fill in two ***
c     ***   passes wrapping around in the streamwise direction   ***
c     --------------------------------------------------------------
c
c     Note: In this case, we are guaranteed that the arrays can be filled 
c           in two passes since the array is guaranteed not to have 
c           dimensions larger than the data in the file.
c
      else
c
c       Fill as much data as can be read without reading beyond the end of
c       the file
c
        do n = 1, 3
          filLo(n) = 1
          loRd(n) = lo(n)
          hiRd(n) = hi(n)
        enddo
        filLo(strmwse_dir) = baseFilePnt
        loRd(strmwse_dir) = baseArrayPnt
        hiRd(strmwse_dir) = loRd(strmwse_dir) + dimFile(strmwse_dir) 
     $                                                         - baseFilePnt

        open(FF_UNIT, file=flct_file, form='unformatted')
        read(FF_UNIT) intDummy
        read(FF_UNIT) realDummy
        read(FF_UNIT) intDummy
        do n = 1, nComp
          call RD_FLCTREC(FF_UNIT, dimFile, loRd, hiRd, filLo, FF_DIMS(dat),
     $                    dat(lo(1),lo(2),lo(3),n))
        enddo
        close(FF_UNIT)

c
c       Now fill the rest of the array starting from the beginning of the file
c
c       Note: The first point in the file in the streamwise direction is 
c             skipped because it is identical to the last point in the file.
c
        filLo(strmwse_dir) = 2
        loRd(strmwse_dir) = hiRd(strmwse_dir) + 1
        hiRd(strmwse_dir) = hi(strmwse_dir)

        open(FF_UNIT, file=flct_file, form='unformatted')
        read(FF_UNIT) intDummy
        read(FF_UNIT) realDummy
        read(FF_UNIT) intDummy
        do n = 1, nComp
          call RD_FLCTREC(FF_UNIT, dimFile, loRd, hiRd, filLo, FF_DIMS(dat),
     $                    dat(lo(1),lo(2),lo(3),n))
        enddo
        close(FF_UNIT)
      endif

c
c
      RETURN
      END
