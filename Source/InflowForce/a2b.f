      program ascii2bin

      implicit none

      integer imax, jmax, kmax, ibc, jbc, kbc, n, i, j, k, NN
      parameter (NN=512)
      integer dummy, ierr, ncomp
      double precision scalex, scaley, scalez, dat(NN,NN), times(NN)
      character(80) infile, outfile
      character char
      logical swirl

      write(*,*) 'input file?'
      read(*,1000) infile
 1000 format(A)
      write(*,*) 'output file?'
      read(*,1000) outfile
      write(*,*) 'ncomp?'
      read(*,*) ncomp
      write(*,*) 'swirl-type turbulence file (y|n)?'
      read(*,1000) char

      swirl = .false.
      if ((char .eq. 'y').or.(char .eq. 'Y')) then
         write(6,*) 'Assuming a swirl-type turbulence file ...'
         swirl = .true.
      endif

      open(1,file=outfile,form='unformatted')
      open(2,file=infile)

      if (ncomp .eq. 3) then
         read(2,*) imax, jmax, kmax
      else
         read(2,*) imax, jmax, kmax, dummy
         if (ncomp .ne. dummy) then
            write(6,*) 'ncomp: ', ncomp, ', dummy: ', dummy
            write(6,*) 'ncomp MUST equal dummy!'
            stop
         endif
      endif

      write(6,*) 'imax, jmax, kmax, ncomp:', imax, jmax, kmax, ncomp

      if (NN.lt.max(imax,max(jmax,kmax))) then
         write(6,*) 'NN must be .ge. max(imax,max(jmax,kmax))'
         stop
      endif

      read(2,*) scalex, scaley, scalez

      write(6,*) 'scalex, scaley, scalez:', scalex, scaley, scalez

      if (swirl) then
         read(2,*) ibc, jbc, kbc, dummy
         read(2,*) (times(i),i=1,kmax)
         write(6,*) 'times:'
         write(6,*) (times(i),i=1,kmax)
      else
         read(2,*) ibc, jbc, kbc
      endif

      write(1) imax, jmax, kmax, ncomp
      write(1) scalex, scaley, scalez

      if (swirl) then
         write(1) ibc, jbc, kbc, dummy
         write(1) (times(i),i=1,kmax)
      else
         write(1) ibc, jbc, kbc
      endif

      do n=1, ncomp
         do k=1, kmax
            do j=1, jmax
               do i=1, imax
                  read(2,*) dat(i,j)
               enddo
            enddo
            write(1) ((dat(i,j), i=1,imax), j=1,jmax)
         enddo
      enddo

      close(1)
      close(2)

      end
