#include <CONSTANTS.H>
#include <REAL.H>
#include <DIMS.H>

#include "mainIBD_F.H"

      subroutine FORT_FILLCRSE (
     $     crse, DIMS(crse),
     $     h, nc
     $     )
      integer nc
      integer DIMDEC(crse)
      REAL_T crse(DIMV(crse),nc)
      REAL_T h(BL_SPACEDIM)
c
      integer  i, j, k, n
      REAL_T x, y, z
      

      do n=1,nc
         do k=ARG_L3(crse), ARG_H3(crse)
            z = (k + 0.5)*h(3)
            do j=ARG_L2(crse), ARG_H2(crse)
               y = (j + 0.5)*h(2)
               do i=ARG_L1(crse), ARG_H1(crse)
                  x = (i + 0.5)*h(1)
                  crse(i,j,k,n)=y*y
               enddo
            enddo
         enddo
      enddo
      end
      

      
      subroutine FORT_FILLFINE (
     $     fine, DIMS(fine),
     $     h, nc
     $     )
      integer nc
      integer DIMDEC(fine)
      REAL_T fine(DIMV(fine),nc)
      REAL_T h(BL_SPACEDIM)
c
      integer  i, j, k, n
      REAL_T x, y, z
      

      do n=1,nc
         do k=ARG_L3(fine), ARG_H3(fine)
            z = (k + 0.5)*h(3)
            do j=ARG_L2(fine), ARG_H2(fine)
               y = (j + 0.5)*h(2)
               do i=ARG_L1(fine), ARG_H1(fine)
                  x = (i + 0.5)*h(1)
                  fine(i,j,k,n)=y*y
               enddo
            enddo
         enddo
      enddo
      end
      
      
