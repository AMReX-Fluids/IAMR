      program bin2ascii

      implicit none

      integer imax, jmax, kmax, ibc, jbc, kbc, n, i, j, k
      integer dummy, ierr, ncomp, NN
      parameter (NN=512)
      double precision scalex, scaley, scalez, dat(NN,NN), times(NN)
      character(80) infile, outfile

      write(*,*) 'input file?'
      read(*,1000) infile
 1000 format(A)
      write(*,*) 'output file?'
      read(*,1000) outfile

      open(1,file=infile,form='unformatted')
      open(2,file=outfile)

      ierr = 0
      read(1,iostat=ierr) imax, jmax, kmax, ncomp
      if (ierr.ne.0) ncomp = 3

      read(1) scalex, scaley, scalez

      ierr = 0
      read(1,iostat=ierr) ibc, jbc, kbc, dummy
      if (ierr.eq.0) then
         read(1) (times(i),i=1,kmax)
      endif

      if (ncomp .eq. 3) then
         write(2,*) imax, jmax, kmax
      else
         write(2,*) imax, jmax, kmax, ncomp
      endif

      write(2,*) scalex, scaley, scalez

      write(6,*) 'imax, jmax, kmax, ncomp: ', imax, jmax, kmax, ncomp
      write(6,*) 'scalex, scaley, scalez: ' , scalex, scaley, scalez

      if (NN.lt.max(imax,max(jmax,kmax))) then
         write(6,*) 'NN must be .ge. max(imax,max(jmax,kmax))'
         stop
      endif

      if (ierr.eq.0) then
         write(2,*) ibc, jbc, kbc, dummy
         write(2,*) (times(i),i=1,kmax)
         write(6,*) 'times:'
         write(6,*) (times(i),i=1,kmax)
      else
         write(2,*) ibc, jbc, kbc
      endif

      do n=1, ncomp
         do k=1, kmax
            read(1) ((dat(i,j), i=1,imax), j=1,jmax)
            do j=1, jmax
               do i=1, imax
                  write(2,*) dat(i,j)
               enddo
            enddo
         enddo
      enddo

      close(1)
      close(2)

      end
