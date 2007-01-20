c
c $Id: FLUCTFILE.F,v 1.12 2007-01-20 18:32:45 lijewski Exp $
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

      integer iunit, dimFile(3), nCmpFile, proc
      REAL_T probSizeFile(3), dxFile(3)

      integer i1bc, j1bc, k1bc, dummy, i
      logical first

      call bl_pd_myproc(proc)

      if (proc.eq.0) print*, 'Entered RD_SCL_FLCTHD'

      first = (infl_type .eq. -1)

      read(iunit,*) dimFile

      if (proc.eq.0) print*, 'dimFile: ', dimFile

      nCmpFile = BL_SPACEDIM

      nCompInflow = nCmpFile

      read(iunit,*) probSizeFile

      if (proc.eq.0) print*, 'probSizeFile: ', probSizeFile
      !
      ! infl_periodic_type only has three integers here.
      !
      read(iunit,*) i1bc, j1bc, k1bc
      !
      ! TODO -- fix this
      !
      infl_type = infl_periodic_type

      if (infl_type .eq. infl_swirl_type) then
         !
         ! Assumes cell-centered data.
         !
         dxFile(1) = probSizeFile(1)/FLOAT(dimFile(1))
         dxFile(2) = probSizeFile(2)/FLOAT(dimFile(2))
         dxFile(3) = probSizeFile(3)/FLOAT(dimFile(3))

         if (.NOT. ASSOCIATED(fluct_times)) then
            allocate(fluct_times(1:dimFile(strmwse_dir)))
         endif
         read(iunit,*) (fluct_times(i),i=1,dimFile(strmwse_dir))
      else
         !
         ! Assumes node-centered data.
         !
         dxFile(1) = probSizeFile(1)/FLOAT(dimFile(1)-1)
         dxFile(2) = probSizeFile(2)/FLOAT(dimFile(2)-1)
         if (dimFile(3) .gt. 1) then
            dxFile(3) = probSizeFile(3)/FLOAT(dimFile(3)-1)
         else
            dxFile(3) = zero
         endif
      endif

      if (first .and. proc.eq.0) then
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
      SUBROUTINE RD_FLCTREC(dimFile, lo, hi, fileLo, FF_DIMS(dat), dat, ncomp)

      implicit none

      integer dimFile(3), i, j, k, n, proc, ncomp
      integer lo(3), hi(3), fileLo(3), FF_DIMDEC(dat), iflctfile(300)
      REAL_T dat(FF_DIMV(dat))

#include "INFL_FORCE_F.H"

      call bl_pd_myproc(proc)
c
c     Build integer version of flct_file name passable to C++
c
      n = len(trim(flct_file))

      do i = 1, n
         iflctfile(i) = ichar(flct_file(i:i))
      end do
c
c     Read the Necessary Data
c
      do k = lo(3), hi(3)
        call getplane(iflctfile(1), n, dat(lo(1),lo(2),k), k, ncomp)
      enddo

      END
