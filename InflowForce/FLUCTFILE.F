c
c $Id: FLUCTFILE.F,v 1.9 2003-12-22 23:28:22 marc Exp $
c

#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"
#include "FLUCTFILE.H"

      block data Inflow_Force_Data

#include "INFL_FORCE_F.H"

      data nCompInflow          / 3/
      data infl_type            /-1/
      data infl_periodic_type   / 1/
      data infl_swirl_type      / 2/
      data tstart_turb          /-1/

      end

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

      integer FF_DIMDEC(holder)
      integer lo(3), hi(3)

      lo(1) = FF_ARG_L1(holder)
      hi(1) = FF_ARG_H1(holder)
      lo(2) = FF_ARG_L2(holder)
      hi(2) = FF_ARG_H2(holder)
      lo(3) = FF_ARG_L3(holder)
      hi(3) = FF_ARG_H3(holder)

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

      integer FF_DIMDEC(holder)
      integer lo(3), hi(3)

      FF_ARG_L1(holder) = lo(1)
      FF_ARG_H1(holder) = hi(1)
      FF_ARG_L2(holder) = lo(2)
      FF_ARG_H2(holder) = hi(2)
      FF_ARG_L3(holder) = lo(3)
      FF_ARG_H3(holder) = hi(3)

      end
c
c ::: -----------------------------------------------------------
c ::: This routine reads the information from the header of a
c ::: inflow/initial conditions fluctuations file.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: iunit      => Fortran unit for the input fluctuations file
c ::: nCmpFile   => Number of components in the file
c ::: dimFile    => Dimensions from the fluctuations file
c ::: probSizeFile  => Domain size from the fluctuations file
c ::: dxFile     => Grid spacing from the fluctuations file
c ::: -----------------------------------------------------------
c
      subroutine RD_SCL_FLCTHD(iunit, nCmpFile, dimFile, probSizeFile, dxFile)

      implicit none

#include "INFL_FORCE_F.H"

      integer iunit, dimFile(3), nCmpFile
      REAL_T probSizeFile(3), dxFile(3)

      integer i1bc, j1bc, k1bc, dummy, ierr, i
      logical first

      first = (infl_type .eq. -1)

      ierr = 0
      read(iunit,iostat=ierr) dimFile, nCmpFile
      if (ierr.ne.0) nCmpFile = BL_SPACEDIM
      if (nCmpFile.ne.nCompInflow) then
         write(6,*) 'nCmpFile=',nCmpFile,'  nCompInflow=',nCompInflow
         call bl_abort('RD_SCL_FLCTHD: nCmpFile != nCompInflow')
      endif

c      if (dimFile(3) .lt. 3) then
c         call bl_abort('RD_SCL_FLCTHD: dimFile(3) must be >= 3')
c      endif

      read(iunit) probSizeFile
      !
      ! Assume infl_swirl_type file until proven otherwise.
      !
      infl_type = infl_swirl_type
      !
      ! infl_periodic_type only has three integers here.
      !
      ierr = 0
      read(iunit,iostat=ierr) i1bc, j1bc, k1bc, dummy
      if (ierr.ne.0) infl_type = infl_periodic_type

      dxFile(1) = probSizeFile(1)/FLOAT(dimFile(1)-1)
      dxFile(2) = probSizeFile(2)/FLOAT(dimFile(2)-1)
      if (dimFile(3) .gt. 1) then
         dxFile(3) = probSizeFile(3)/FLOAT(dimFile(3)-1)
      else
         dxFile(3) = zero
      endif

      if (infl_type .eq. infl_swirl_type) then
         if (.NOT. ASSOCIATED(fluct_times)) then
            allocate(fluct_times(1:dimFile(3)))
         endif
         read(iunit) (fluct_times(i),i=1,dimFile(3))
      endif

      if (first) then
         if (infl_type .eq. infl_periodic_type) then
            write(6,*) 'RD_SCL_FLCTHD: infl_type: infl_periodic_type'
         else if (infl_type .eq. infl_swirl_type) then
            write(6,*) 'RD_SCL_FLCTHD: infl_type: infl_swirl_type'
         else
            call bl_abort('RD_SCL_FLCTHD: unknown infl_type')
         endif
      endif

      end
c
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

      implicit none

      integer iunit, dimFile(3)
      integer arrLo(3), arrHi(3), fileLo(3)
      integer FF_DIMDEC(dat)
      REAL_T dat(FF_DIMV(dat))

      integer i, j, k, ijoff, ipad, kpad
      REAL_T dummy
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
     $              (dummy, i=1, ipad), j=arrLo(2),arrHi(2))
      enddo
c
c     *** Skip through the remaining records for this variable ***
c
      kpad    = dimFile(3) - (fileLo(3)-1 + arrHi(3) - arrLo(3) + 1)
      do k = 1, kpad
        READ(iunit) dummy
      enddo

      END

