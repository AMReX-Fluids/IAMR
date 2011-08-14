      program dummydata

      implicit none

      integer imax, jmax, kmax
      parameter (imax=49, jmax=25, kmax=1)

      integer ibc, jbc, kbc, i, j, k
      real*8 scalex, scaley, scalez, convVel, dt, t, pi, tMax, dx, x
      real*8 udat(imax,jmax,kmax), vdat(imax,jmax,kmax), 
     $       wdat(imax,jmax,kmax)

c
c
c
      convVel = 1.5d0
      scalex = 0.7500d+01
      scaley = convVel
      scalez = 0.0000d+00
      ibc = 1
      jbc = 1
      kbc = 1

      dt = scaley / FLOAT(jmax-1) / convVel
      dx = scalex / FLOAT(imax-1)
      tMax = scaley / convVel
      pi = ACOS(-1.0d0)

      do k = 1, kmax
        do j = 1, jmax
          do i = 1, imax
            x = FLOAT(i-1) * dx
            t = FLOAT(j-1) * dt

c           udat(i,j,k) = SIN(2.0d0 * pi * t / tMax)
c           vdat(i,j,k) = COS(2.0d0 * pi * t / tMax)
c           wdat(i,j,k) = 0.0d0
            udat(i,j,k) = SIN(2.0d0 * pi * x / scalex)
            vdat(i,j,k) = COS(2.0d0 * pi * x / scalex)
            wdat(i,j,k) = 0.0d0
          enddo
        enddo
      enddo
c
c
      open(1,file='DummyData',form='unformatted')
      write(1) imax, jmax, kmax
      write(1) scalex, scaley, scalez
      write(1) ibc, jbc, kbc

      do k=1, kmax
        WRITE(1) ((udat(i,j,k), i=1,imax), j=1,jmax)
      enddo
      do k=1, kmax
        WRITE(1) ((vdat(i,j,k), i=1,imax), j=1,jmax)
      enddo
      do k=1, kmax
        WRITE(1) ((wdat(i,j,k), i=1,imax), j=1,jmax)
      enddo

      close(1)

c
c
      stop
      end
