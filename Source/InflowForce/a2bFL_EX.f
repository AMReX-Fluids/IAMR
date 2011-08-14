      PROGRAM BIN2ASCII

c
c     *****************************
c     *** Variable Declarations ***
c     *****************************
c
      IMPLICIT NONE

c
c     *** Passed Variables ***
c
      integer imax, jmax, kmax, ibc, jbc, kbc, n, i, j, k
      real*8 scalex, scaley, scalez, dat(1027,1027)
      character*80 infile, outfile


c 
c     *************************************
c     *** Read Header from Restart File ***
c     *************************************
c
      WRITE(*,*) 'input file?'
      READ(*,1000) infile
 1000 FORMAT(A)
      WRITE(*,*) 'output file?'
      READ(*,1000) outfile

      OPEN(1,file=outfile,form='unformatted')
      OPEN(2,file=infile)

      read(2,*) imax, jmax, kmax
      read(2,*) scalex, scaley, scalez
      read(2,*) ibc, jbc, kbc
c     write(1) imax+3, jmax, kmax
      write(1) imax+3, jmax+1, kmax
      write(1) scalex*dfloat(imax+2)/dfloat(imax),
     &    scaley, scalez
      write(1) ibc, jbc, 0

      DO n=1, 2
        DO k=1, kmax
          DO j=1, jmax
            DO i=1, imax
              read(2,*) dat(i,j)
            ENDDO
          ENDDO
          do j=1,jmax
          do i=1,3
          dat(imax+i,j)=dat(i,j) 
          enddo
          enddo
          do j=1,1
          do i=1,imax+3
          dat(i,j+jmax)=dat(i,j) 
          enddo
          enddo
          write(1) ((dat(i,j), i=1,imax+3), j=jmax+1,1,-1)
        ENDDO
      ENDDO

      CLOSE(1)
      CLOSE(2)

c
c
      STOP
      END


