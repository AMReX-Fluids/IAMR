#include "REAL.H"
#include "main_F.H"
#include "ArrayLim.H"
#include "BC_TYPES.H"

#define SDIM 3

      subroutine FORT_READ_SWIRL_HEADER (fname,flen,imax,jmax,kmax,
     &     ncomp,scalex,scaley,scalez,times,NN)

      implicit none

      integer imax, jmax, kmax, ncomp, NN, flen, fname(flen)
      double precision scalex, scaley, scalez, times(NN)

      integer maxlen, i, ierr, ibc, jbc, kbc, dummy
      character filet*(flen)

      if (flen.LE.0) then
         call bl_abort('FORT_READ_SWIRL_HEADER: no filename provided')
      endif

      do i = 1, flen
         filet(i:i) = char(fname(i))
      end do

      open(1,status='old',form='unformatted',file=filet(1:flen))

      ierr = 0
      read(1,iostat=ierr) imax, jmax, kmax, ncomp
      if (ierr.ne.0) ncomp = 3

      read(1) scalex, scaley, scalez

      ierr = 0
      read(1,iostat=ierr) ibc, jbc, kbc, dummy
      if (ierr.eq.0) then
         if (kmax.gt.NN) then
            write(6,*) 'NN = ', NN, ', kmax = ', kmax
            call bl_abort('FORT_READ_SWIRL_HEADER: kmax.gt.NN')
         endif
         read(1) (times(i),i=1,kmax)
      else
         call bl_abort('FORT_READ_SWIRL_HEADER: not a swirl-type file')
      endif

      close(1)

      end

      subroutine FORT_FILL(dat,imax,jmax,kmax,ncomp,time,putidx,
     &     getidx,fname,flen,count)

      implicit none

      integer imax,jmax,kmax,ncomp,getidx,putidx,flen,fname(flen),count
      double precision dat(imax,jmax,kmax,ncomp), time

      integer imaxt, jmaxt, kmaxt, ncompt
      integer i, j, k, n, nn
      double precision scalex, scaley, scalez, dummy
      double precision, allocatable :: times(:)
      character filet*(flen)

      if (flen.LE.0) then
         call bl_abort('FORT_FILL: no filename provided')
      endif

      do i = 1, flen
         filet(i:i) = char(fname(i))
      end do

      open(1,status='old',form='unformatted',file=filet(1:flen))

      rewind(1)

      read(1) imaxt, jmaxt, kmaxt, ncompt
      read(1) dummy ! scalex, scaley, scalez
      read(1) dummy ! ibc, jbc, kbc, dummy

      if (imax.ne.imaxt)   call bl_abort('FORT_FILL: imax.ne.imaxt')
      if (jmax.ne.jmaxt)   call bl_abort('FORT_FILL: jmax.ne.jmaxt')
      if (ncomp.ne.ncompt) call bl_abort('FORT_FILL: ncomp.ne.ncompt')

      allocate(times(kmaxt))

      read(1) (times(i),i=1,kmaxt)

      if (times(getidx).ne.time) then
         write(6,*) 'times(getidx): ' , times(getidx), ', time: ', time
         call bl_abort('FORT_FILL: times(getidx).ne.time');
      endif

      deallocate(times)

      do n = 1, ncomp

         if (getidx.gt.1) then
            do k = 1, getidx-1
               read(1) dummy
            enddo
         endif

         do nn = 1,count
            read(1) ((dat(i,j,putidx+nn-1,n), i=1,imax), j=1,jmax)
         enddo

         if ((getidx+count-1).lt.kmaxt) then
            do k = getidx+count, kmaxt
               read(1) dummy
            enddo
         endif

      enddo

100   close(1)

      end

      subroutine FORT_WRITE(dat,imax,jmax,kmax,ncomp,times,
     &     scalex,scaley,scalez,fname,flen)

      implicit none

      integer imax, jmax, kmax, ncomp, flen, fname(flen)
      double precision scalex, scaley, scalez
      double precision dat(imax,jmax,kmax,ncomp), times(kmax)

      integer i, j, k, n, dummy
      character filet*(flen)

      if (flen.LE.0) then
         call bl_abort('FORT_WRITE: no filename provided')
      endif

      do i = 1, flen
         filet(i:i) = char(fname(i))
      end do

      open(1,form='unformatted',file=filet(1:flen))

      rewind(1)

      dummy = 1

      write(1) imax, jmax, kmax, ncomp
      write(1) scalex, scaley, scalez
      write(1) dummy, dummy, dummy, dummy
      write(1) (times(i),i=1,kmax)

      do n=1, ncomp
         do k=1, kmax
            write(1) ((dat(i,j,k,n), i=1,imax), j=1,jmax)
         enddo
      enddo

100   close(1)

      end
