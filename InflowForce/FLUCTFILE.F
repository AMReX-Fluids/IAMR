c
c $Id: FLUCTFILE.F,v 1.2 2000-06-14 21:03:45 sstanley Exp $
c

#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"
#include "FLUCTFILE.H"


c ::: -----------------------------------------------------------
c ::: This routine sets the values for the lo() and hi() arrays
c ::: from the ARG_L1, ARG_H1, ... macros.  This is done since
c ::: it is more convenient to use the lo() and hi() arrays.
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: DIMS(holder)  => index extent of place holder array
c ::: lo(3)        <=  lower index limits
c ::: hi(3)        <=  upper index limits
c ::: -----------------------------------------------------------

      subroutine FF_SET_LOHI(FF_DIMS(holder), lo, hi)

      implicit none

c
c     :::: Passed Variables ::::
c
      integer FF_DIMDEC(holder)
      integer lo(3), hi(3)

c
c     --------------------------------------
c     :::: Set Values for lo() and hi() ::::
c     --------------------------------------
c
      lo(1) = FF_ARG_L1(holder)
      hi(1) = FF_ARG_H1(holder)
      lo(2) = FF_ARG_L2(holder)
      hi(2) = FF_ARG_H2(holder)
      lo(3) = FF_ARG_L3(holder)
      hi(3) = FF_ARG_H3(holder)

c
c
      return
      end



c ::: -----------------------------------------------------------
c ::: This routine sets the values for the ARG_L1, ARG_H1, ... macros
c ::: from the lo() and hi() arrays.  This is done since
c ::: it is more convenient to use the macros to dimension arrays.
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: FF_DIMS(holder) <=  index extent of place holder array
c ::: lo(3)            => lower index limits
c ::: hi(3)            => upper index limits
c ::: -----------------------------------------------------------

      subroutine FF_SET_ARGS(FF_DIMS(holder), lo, hi)

      implicit none

c
c     :::: Passed Variables ::::
c
      integer FF_DIMDEC(holder)
      integer lo(3), hi(3)

c
c     --------------------------------------
c     :::: Set Values for lo() and hi() ::::
c     --------------------------------------
c
      FF_ARG_L1(holder) = lo(1)
      FF_ARG_H1(holder) = hi(1)
      FF_ARG_L2(holder) = lo(2)
      FF_ARG_H2(holder) = hi(2)
      FF_ARG_L3(holder) = lo(3)
      FF_ARG_H3(holder) = hi(3)

c
c
      return
      end



c
c ::: -----------------------------------------------------------
c ::: This routine reads the information from the header of a
c ::: inflow/initial conditions fluctuations file.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: iunit      => Fortran unit for the input fluctuations file
c ::: dimFile    => Dimensions from the fluctuations file
c ::: probSizeFile  => Domain size from the fluctuations file
c ::: dxFile     => Grid spacing from the fluctuations file
c ::: -----------------------------------------------------------
c
      SUBROUTINE RD_FLCTHD(iunit, dimFile, probSizeFile, dxFile)

c
c     *****************************
c     *** Variable Declarations ***
c     *****************************
c
      IMPLICIT NONE

c
c     *** Passed Variables ***
c
      INTEGER iunit, dimFile(3)
      REAL_T probSizeFile(3), dxFile(3)

c
c     *** Local Variables ***
c
      integer i1bc, j1bc, k1bc


c 
c     *************************************
c     *** Read Header from Restart File ***
c     *************************************
c
      READ(iunit) dimFile
      READ(iunit) probSizeFile
      READ(iunit) i1bc, j1bc, k1bc

      dxFile(1) = probSizeFile(1)/FLOAT(dimFile(1)-1)
      dxFile(2) = probSizeFile(2)/FLOAT(dimFile(2)-1)
#if (BL_SPACEDIM==3)
      dxFile(3) = probSizeFile(3)/FLOAT(dimFile(3)-1)
#else
      dxFile(3) = zero
#endif

c
c
      RETURN
      END



c ::: -----------------------------------------------------------
c ::: This routine reads a record of data from an inflow/initial
c ::: conditions fluctuations file.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: iunit         => Fortran unit for the input fluctuations file
c ::: dimFile       => Dimensions from the header of the fluctuations
c ::: arrLo, arrHi  => Range of the array to fill.  This must satisfy
c :::                    rdLgth(1) = arrHi(1) - arrLo(1) + 1
c :::                    dimFile(1) >= fileLo(1) + rdLgth(1) - 1
c ::: fileLo         => Offset into the array in the file indicating the
c :::                    first point to be read.
c ::: FF_DIMS(dat)  => Dimensions of the array dat
c ::: dat          <=  Array to fill
c ::: -----------------------------------------------------------
c
      SUBROUTINE RD_FLCTREC(iunit, dimFile, arrLo, arrHi, fileLo,
     $                      FF_DIMS(dat), dat)

c
c     *****************************
c     *** Variable Declarations ***
c     *****************************
c
      implicit none

c
c     *** Passed Variables ***
c
      integer iunit, dimFile(3)
      integer arrLo(3), arrHi(3), fileLo(3)
      integer FF_DIMDEC(dat)
      REAL_T dat(FF_DIMV(dat))

c
c     *** Local Variables ***
c
      integer i, j, k, ijoff, ipad, kpad
      REAL_T dummy
      REAL_T locDat(:,:,:)
      allocatable locDat
      allocate (locDat(dimFile(1),dimFile(2),dimFile(3)))


c
c     ***********************************
c     *** Read Data From Restart File ***
c     ***********************************
c
c     The logic here seems quite complicated.  This is done so that
c     only the data which is needed to fill the FAB is read from the
c     disk and stored into memory.  This should speed up the I/O as
c     well as keep the memory usage down.
c
c     *** Skip through records to the K-range desired as ***
c     ***   specified by fileLo(3)                        ***
c
      do k = 1, fileLo(3)-1
        READ(iunit) dummy
      enddo

c
c     *** Read the Necessary Data ***
c
      ijoff = (fileLo(2) - 1) * dimFile(1)
      ipad    = dimFile(1) - (fileLo(1)-1 + arrHi(1) - arrLo(1) + 1)

      do k = arrLo(3), arrHi(3)
        READ(iunit) (dummy, i=1, ijoff),
     $             ((dummy, i=1, fileLo(1)-1),
     $              (dat(i,j,k), i=arrLo(1),arrHi(1)),
     $              (dummy, i=1, ipad),                 j=arrLo(2),arrHi(2))
      enddo

c
c     *** Skip through the remaining records for this variable ***
c
      kpad    = dimFile(3) - (fileLo(3)-1 + arrHi(3) - arrLo(3) + 1)
      do k = 1, kpad
        READ(iunit) dummy
      enddo

c
c
      RETURN
      END

