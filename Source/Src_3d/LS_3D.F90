#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <LS_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 3
#define BOGUS 1.d30
#define LARGEINT 100000000

module LS_3d_module

  implicit none

  private

  public :: phiupd, switch, lscfl, findintrfce, UPDATEF, FINDDIST, &
         POLYVAL, GRADPVAL, GEPP, narrowband, retypify, fastmarch, &
         UPDATE, EVAL, ADDNODE, UPDATENODE, RMVNODE, fastmarch2, UPDATE2, &
         EVAL2, mine, nbandnumify

contains

      integer function phiupd(phi, DIMS(phi), phin, DIMS(phin), &
                                  uadv, DIMS(uadv), vadv, DIMS(vadv), wadv, DIMS(wadv), &
                                  nband, nbandsize, mine, minesize, &
                                  lo, hi, dt, dx, type,DIMS(type)) bind(C,name="phiupd")

      implicit none
      INTEGER  DIMDEC(phi)
      INTEGER  DIMDEC(phin)
      INTEGER  DIMDEC(uadv)
      INTEGER  DIMDEC(vadv)
      INTEGER  DIMDEC(wadv)
      REAL_T    phi(DIMV(phi))
      REAL_T   phin(DIMV(phin))
      REAL_T   uadv(DIMV(uadv))
      REAL_T   vadv(DIMV(vadv))
      REAL_T   wadv(DIMV(wadv))
      integer  nbandsize, minesize
      integer  nband(nbandsize,SDIM)
      integer  mine(minesize,SDIM)
      INTEGER  lo(SDIM), hi(SDIM)
      REAL_T   dt
      REAL_T   dx(SDIM)
      integer  DIMDEC(type)
      integer  type(DIMV(type))

      
#include <probdata.H>

      integer  i, j,k, p,s,t
      REAL_T   Dxm, Dxp
      REAL_T   Dym, Dyp
      REAL_T   Dzm, Dzp
      REAL_T   phix, phiy,phiz, phixx, phiyy, phizz, phixy, phixz, phiyz, kappa
      REAL_T   uavg, vavg, wavg, Fo, Fkappa, Fadv
      REAL_T   SWITCH
      REAL_T   Dxpp, Dxmm, Dxpm, Dypp, Dymm, Dypm, Dzpp, Dzmm, Dzpm
      REAL_T   Ad, Bd, Cd, Dd, Ed, Fd, Delp, Delm      
      p=1
      do while (nband(p,1) .GT. -LARGEINT)
      
        i = nband(p,1)
        j = nband(p,2)
        k = nband(p,3)
        p = p + 1
        
        if (max(type(i+1,j,k), type(i-1,j,k), type(i,j+1,k), type(i,j-1,k), type(i,j,k+1), type(i,j,k-1), &
            type(i+1,j+1,k), type(i-1,j-1,k), type(i+1,j-1,k), type(i-1,j+1,k), &
             type(i+1,j,k+1), type(i-1,j,k+1), type(i+1,j,k-1), type(i-1,j,k+1),  &
             type(i,j+1,k+1), type(i,j-1,k-1), type(i,j-1,k+1), type(i,j+1,k-1), &
             type(i+1,j+1,k+1), type(i-1,j-1,k-1), type(i+1,j-1,k-1), type(i-1,j+1,k+1), &
             type(i+1,j+1,k-1), type(i-1,j-1,k+1), type(i+1,j-1,k+1), type(i-1,j+1,k-1)) .LE. 1 ) then

           Dxm   = ( phi(i,j,k) - phi(i-1,j,k) ) / dx(1)
           Dxp   = ( phi(i+1,j,k) - phi(i,j,k) ) / dx(1)
           Dym   = ( phi(i,j,k) - phi(i,j-1,k) ) / dx(2)
           Dyp   = ( phi(i,j+1,k) - phi(i,j,k) ) / dx(2)
           Dzm   = ( phi(i,j,k) - phi(i,j,k-1) ) / dx(3)
           Dzp   = ( phi(i,j,k+1) - phi(i,j,k) ) / dx(3)
        
           phix  = ( phi(i+1,j,k) - phi(i-1,j,k) ) / (2*dx(1))
           phiy  = ( phi(i,j+1,k) - phi(i,j-1,k) ) / (2*dx(2))
           phiz  = ( phi(i,j,k+1) - phi(i,j,k-1) ) / (2*dx(3))
        
           phixx = ( phi(i+1,j,k) - 2*phi(i,j,k) + phi(i-1,j,k) ) &
               / (dx(1)**2)
           phiyy = ( phi(i,j+1,k) - 2*phi(i,j,k) + phi(i,j-1,k) ) &
               / (dx(2)**2)
           phizz = ( phi(i,j,k+1) - 2*phi(i,j,k) + phi(i,j,k-1) ) &
               / (dx(3)**2)     
           phixy = ( phi(i+1,j+1,k) + phi(i-1,j-1,k) &
               - phi(i+1,j-1,k) - phi(i-1,j+1,k) ) &
               / (4*dx(1)*dx(2))
           phixz = ( phi(i+1,j,k+1) + phi(i-1,j,k-1) &
               - phi(i+1,j,k-1) - phi(i-1,j,k+1) ) &
               / (4*dx(1)*dx(3))     
           phiyz = ( phi(i,j+1,k+1) + phi(i,j-1,k-1) &
               - phi(i,j-1,k+1) - phi(i,j+1,k-1) ) &
               / (4*dx(1)*dx(3))     
    
        
           if (phix**2 + phiy**2 + phizz**2 .GT. 0 ) then
              Fkappa = -kapb*( (phiyy + phizz)*phixx**2 + (phixx + phizz)*phiyy**2 + (phixx + phiyy)*phizz**2 & 
                      - 2*phix*phiy*phixy - 2*phix*phiz*phixz - 2*phiy*phiz*phiyz) &
                      / (phix**2 + phiy**2 + phiz**2)
      
           else
              Fkappa = 0
           endif
        

           if (LSorder .EQ. 2) then
              
              Dxpp = 0
              if (i+2 .LE. hi(1) +1) then
                 if(type(i+2,j,k) .LT. 2) then
                    Dxpp = (phi(i+2,j,k) - 2*phi(i+1,j,k) + phi(i,j,k))/dx(1)
                 endif
              endif
              
              Dxmm = 0
              if (i-2 .GE. lo(1)) then
                 if (type(i-2,j,k) .LT. 2) then
                    Dxmm = (phi(i,j,k) - 2* phi(i-1,j,k) + phi(i-2,j,k))/dx(1)
                 endif
              endif
       
              Dxpm = (phi(i+1,j,k) - 2* phi(i,j,k) + phi(i-1,j,k))/dx(1)
              
              Dypp = 0
              if(j+2 .LE. hi(2)) then
                 if (type(i,j+2,k) .LT. 2) then
                    Dypp = (phi(i,j+2,k) - 2* phi(i,j+1,k) + phi(i,j,k))/dx(2)
                 endif        
              endif
              
              Dymm = 0
              if(j-2 .GE. lo(2)) then
                 if(type(i,j-2,k) .LT. 2) then
                    Dymm = (phi(i,j,k) - 2* phi(i,j-1,k) + phi(i,j-2,k))/dx(2)
                 endif
              endif
              
              Dypm = (phi(i,j+1,k) - 2* phi(i,j,k) + phi(i,j-1,k))/dx(2)
              
              Dzpp = 0
              if (k+2 .LE. hi(3) +1) then
                 if(type(i,j,k+2) .LT. 2) then
                    Dzpp = (phi(i,j,k+2) - 2*phi(i,j,k+1) + phi(i,j,k))/dx(3)
                 endif
              endif
              
              Dzmm = 0
              if (k-2 .GE. lo(3)) then
                 if (type(i,j,k-2) .LT. 2) then
                    Dzmm = (phi(i,j,k) - 2* phi(i,j,k-1) + phi(i,j,k-2))/dx(3)
                 endif
              endif
       
              Dzpm = (phi(i,j,k+1) - 2* phi(i,j,k) + phi(i,j,k-1))/dx(3)
                            
        

              Ad = Dxm + .5d0*SWITCH(Dxmm,Dxpm)
              Bd = Dxp + .5d0*SWITCH(Dxpp,Dxpm)
              Cd = Dym + .5d0*SWITCH(Dymm,Dypm)
              Dd = Dyp + .5d0*SWITCH(Dypp,Dypm)
              Ed = Dzm + .5d0*SWITCH(Dzmm,Dzpm)
              Fd = Dzp + .5d0*SWITCH(Dzpp,Dzpm)
              
              Delp = (max(Ad,0.d0)**2 + min(Bd,0.d0)**2 + max(Cd,0.d0)**2 + min(Dd,0.d0)**2 + max(Ed,0.d0)**2 + min(Fd,0.d0)**2)**(.5d0)
              
              Delm = (max(Bd,0.d0)**2 + min(Ad,0.d0)**2 + max(Dd,0.d0)**2 + min(Cd,0.d0)**2 + max(Fd,0.d0)**2 + min(Ed,0.d0)**2)**(.5d0)
              
           endif
        
           Fo = kapa
           if (LSorder .EQ. 1) then
              if (Fo .GT. 0) then
                 Fo = Fo*( ( max(Dxm,0.d0) + min(Dxp,0.d0) )**2 & 
                     + ( max(Dym,0.d0) + min(Dyp,0.d0) )**2  &
          + ( max(Dzm,0.d0) + min(Dzp,0.d0) )**2)**(.5d0) 
              else
                 Fo = Fo*( ( min(Dxm,0.d0) + max(Dxp,0.d0) )**2 &
                     + ( min(Dym,0.d0) + max(Dyp,0.d0) )**2   &
           + ( min(Dzm,0.d0) + max(Dzp,0.d0) )**2)**(.5d0) 
              endif
           else if (LSorder .EQ. 2) then
              if (Fo .GT. 0) then
                 Fo = Fo*Delp
              else
                 Fo = Fo*Delm
              endif
           endif
           
           Fadv = 0
           uavg = ( uadv(i,j,k) + uadv(i+1,j,k) ) * .5d0
           vavg = ( vadv(i,j,k) + vadv(i,j+1,k) ) * .5d0
           wavg = ( wadv(i,j,k) + wadv(i,j+1,k) ) * .5d0
           
           if (LSorder .EQ. 1) then

              if (uavg .GT. 0) then
                 Fadv = Fadv + uavg*Dxm
              else 
                 Fadv = Fadv + uavg*Dxp
              endif
              
              if (vavg .GT. 0) then
                 Fadv = Fadv + vavg*Dym
              else 
                 Fadv = Fadv + vavg*Dyp
              endif    
              
              if (wavg .GT. 0) then
                 Fadv = Fadv + wavg*Dzm
              else 
                 Fadv = Fadv + wavg*Dzp
              endif                 
              
           else if (LSorder .EQ. 2) then
              Fadv = uavg*phix + vavg*phiy + wavg*phiz
           endif
           
           phin(i,j,k) = phi(i,j,k) - dt*( Fo + Fkappa + Fadv )
           
        endif
      enddo
      
      PHIUPD = 0
      
      p = 1
      
      do while (mine(p,1) .GT. -LARGEINT)
         i = mine(p,1)
         j = mine(p,2)
         k = mine(p,3)
         p = p + 1
         
         if (sign(1.d0,phi(i,j,k))*sign(1.d0,phin(i,j,k)) .LE. 0) then
            PHIUPD = 1
            exit
         endif 
      enddo
      
      return 
      end function phiupd     
      

      REAL_T function switch(x,y)bind(C, name="switch")
      implicit none      
      REAL_T x,y
      
      if (x*y .GE. 0) then
        if (abs(x) .LE. abs(y)) then
         SWITCH = x
        else 
         SWITCH = y
        endif
      else
        SWITCH = 0
      endif
      end function switch 


      REAL_T function lscfl(phi, DIMS(phi), uadv, DIMS(uadv), vadv, DIMS(vadv), wadv, DIMS(wadv), &
                            nband, nbandsize, mine, minesize, &
                            lo, hi, phit, dx, type, DIMS(type))  &
                            bind(C,name="lscfl")
            
      implicit none
      INTEGER  DIMDEC(phi)
      INTEGER  DIMDEC(uadv)
      INTEGER  DIMDEC(vadv)
      INTEGER  DIMDEC(wadv)
      REAL_T    phi(DIMV(phi))
      REAL_T   uadv(DIMV(uadv))
      REAL_T   vadv(DIMV(vadv))
      REAL_T   wadv(DIMV(wadv))
      integer  nbandsize, minesize
      integer  nband(nbandsize,SDIM)
      integer  mine(minesize,SDIM)
      INTEGER  lo(SDIM), hi(SDIM)
      REAL_T   dt
      REAL_T   dx(SDIM)
      REAL_T   phit
      integer  DIMDEC(type)
      integer  type(DIMV(type))
#include <probdata.H>

      integer i, j, k, p
      REAL_T phix, phiy, phiz, phixx, phiyy, phizz, phixy, phixz, phiyz 
      REAL_T kappa, speed    
      REAL_T phidt

      phidt = phit
      
      p = 1
      do while (nband(p,1) .GT. -LARGEINT)
         
         i = nband(p,1)
         j = nband(p,2)
         k = nband(p,3)
         p = p + 1
         
       if (max(type(i+1,j,k), type(i-1,j,k), type(i,j+1,k), type(i,j-1,k), type(i,j,k+1), type(i,j,k-1), &
          type(i+1,j+1,k), type(i-1,j-1,k), type(i+1,j-1,k), type(i-1,j+1,k),  &
           type(i+1,j,k+1), type(i-1,j,k+1), type(i+1,j,k-1), type(i-1,j,k+1),  &
           type(i,j+1,k+1), type(i,j-1,k-1), type(i,j-1,k+1), type(i,j+1,k-1), &
           type(i+1,j+1,k+1), type(i-1,j-1,k-1), type(i+1,j-1,k-1), type(i-1,j+1,k+1), &
           type(i+1,j+1,k-1), type(i-1,j-1,k+1), type(i+1,j-1,k+1), type(i-1,j+1,k-1)) .LE. 1 ) then

           phix  = ( phi(i+1,j,k) - phi(i-1,j,k) ) / (2*dx(1))
           phiy  = ( phi(i,j+1,k) - phi(i,j-1,k) ) / (2*dx(2))
           phiz  = ( phi(i,j,k+1) - phi(i,j,k-1) ) / (2*dx(3))
        
           phixx = ( phi(i+1,j,k) - 2*phi(i,j,k) + phi(i-1,j,k) ) &
               / (dx(1)**2)
           phiyy = ( phi(i,j+1,k) - 2*phi(i,j,k) + phi(i,j-1,k) ) &
               / (dx(2)**2)
           phizz = ( phi(i,j,k+1) - 2*phi(i,j,k) + phi(i,j,k-1) ) &
               / (dx(3)**2)     
           phixy = ( phi(i+1,j+1,k) + phi(i-1,j-1,k) &
               - phi(i+1,j-1,k) - phi(i-1,j+1,k) ) &
               / (4*dx(1)*dx(2))
           phixz = ( phi(i+1,j,k+1) + phi(i-1,j,k-1) &
               - phi(i+1,j,k-1) - phi(i-1,j,k+1) ) &
               / (4*dx(1)*dx(3))     
           phiyz = ( phi(i,j+1,k+1) + phi(i,j-1,k-1) &
               - phi(i,j-1,k+1) - phi(i,j+1,k-1) ) &
               / (4*dx(1)*dx(3))     
    
        
           if (phix**2 + phiy**2 + phizz**2 .GT. 0 ) then
              kappa = ( (phiyy + phizz)*phixx**2 + (phixx + phizz)*phiyy**2 + (phixx + phiyy)*phizz**2 &
                      - 2*phix*phiy*phixy - 2*phix*phiz*phixz - 2*phiy*phiz*phiyz) &
                      / ((phix**2 + phiy**2 + phiz**2)**(3./2.))
     
           else
              kappa = 0.d0
           endif          
            
            speed = ( (.5d0*( uadv(i,j,k) + uadv(i+1,j,k) ) )**2 + ( .5d0*( vadv(i,j,k) +vadv(i,j+1,k) ) )**2 + ( .5d0*( wadv(i,j,k) +wadv(i,j,k+1) ) )**2)**.5d0 &
                + abs(kapa-kapb*kappa)
             
            phidt = min( phit, &
                        1/(max(1.d0/phidt,speed/(.8d0*min(dx(1),dx(2),dx(3))), 6*abs(kapb)/(.8d0*min(dx(1),dx(2),dx(3))**2) ))  )     

         endif
      enddo
        
      lscfl = phidt
      end function lscfl 

      


      subroutine findintrfce( phi, DIMS(phi), phin, DIMS(phin), type, DIMS(type), &
                              lo, hi, dx, intfacenump, intfacenumn, intfacep,intfacen, &
                              nband, nbandsize, intfacesize) &
                              bind(C, name="findintrfce")
      implicit none
      INTEGER    DIMDEC(phi)
      INTEGER    DIMDEC(phin)
      INTEGER    DIMDEC(type)
      REAL_T      phi(DIMV(phi))
      REAL_T     phin(DIMV(phin))
      integer    type(DIMV(type))
      INTEGER    lo(SDIM), hi(SDIM)    
      REAL_T     dx(SDIM)
      INTEGER    intfacenump, intfacenumn, intfacesize
      INTEGER    intfacep(intfacesize,SDIM), intfacen(intfacesize,SDIM)
      integer    nbandsize
      integer    nband(nbandsize,SDIM)

      
!     Local variables
      INTEGER    i, j, k, r
      
      intfacenump=0
      intfacenumn=0
      
      r = 1
      do while(nband(r,1) .GT. -LARGEINT)
        i = nband(r,1)
        j = nband(r,2)
        k = nband(r,3)
        r = r + 1      
        
        phin(i,j,k) = sign(BOGUS,phi(i,j,k))
      enddo
      
      r = 1
      do while(nband(r,1) .GT. -LARGEINT)
        i = nband(r,1)
        j = nband(r,2)
        k = nband(r,3)
        r = r + 1
        
        call UPDATEF(i, j, k, phi, DIMS(phi), phin, DIMS(phin), type, DIMS(type), &
                    intfacep, intfacen, intfacenump, intfacenumn, &
                    dx, intfacesize, lo, hi) 
           
        if( (i .EQ. lo(1)) .OR. (j .EQ. lo(2)) .OR. k .EQ. lo(3) ) then
           
           if(i .EQ. lo(1)) then
                call UPDATEF(i-1, j, k, phi, DIMS(phi), phin, DIMS(phin), type, DIMS(type), &
                            intfacep, intfacen, intfacenump, intfacenumn,  &
                            dx, intfacesize, lo, hi )
             endif
             
             if(j .EQ. lo(2)) then 
              call UPDATEF(i, j-1, k, phi, DIMS(phi), phin, DIMS(phin), type, DIMS(type), &
                         intfacep, intfacen, intfacenump, intfacenumn, &
                         dx, intfacesize, lo, hi )
             endif
             
             if(k .EQ. lo(3)) then 
              call UPDATEF(i, j, k-1, phi, DIMS(phi), phin, DIMS(phin), type, DIMS(type), &
                         intfacep, intfacen, intfacenump, intfacenumn, &
                         dx, intfacesize, lo, hi )
             endif             
             
             if (i .EQ. lo(1) .AND. j .EQ. lo(2)) then
               call UPDATEF(i-1, j-1, k, phi, DIMS(phi), phin, DIMS(phin), type, DIMS(type), &
                          intfacep, intfacen, intfacenump, intfacenumn, &
                          dx, intfacesize, lo, hi )
             endif

             if (i .EQ. lo(1) .AND. k .EQ. lo(3)) then
               call UPDATEF(i-1, j, k-1, phi, DIMS(phi), phin, DIMS(phin), type, DIMS(type), &
                          intfacep, intfacen, intfacenump, intfacenumn, &
                          dx, intfacesize, lo, hi )
             endif             

             if (j .EQ. lo(2) .AND. k .EQ. lo(3)) then
               call UPDATEF(i, j-1, k-1, phi, DIMS(phi), phin, DIMS(phin), type, DIMS(type), &
                          intfacep, intfacen, intfacenump, intfacenumn, &
                          dx, intfacesize, lo, hi )
             endif   


             if (i .EQ. lo(1) .AND. j .EQ. lo(2) .AND. k .EQ. lo(3)) then
               call UPDATEF(i-1, j-1, k-1, phi, DIMS(phi), phin, DIMS(phin), type, DIMS(type), &
                          intfacep, intfacen, intfacenump, intfacenumn, &
                          dx, intfacesize, lo, hi )
             endif   
             
          endif
      enddo
      
      nband(1,1) = -LARGEINT
      nband(1,2) = -LARGEINT
      nband(1,3) = -LARGEINT
      end subroutine findintrfce
      

      subroutine UPDATEF(i,j,k, phi, DIMS(phi), phin, DIMS(phin), type, DIMS(type), &
                        intfacep, intfacen, intfacenump, intfacenumn,  &
                        dx, intfacesize, lo, hi  )
      
      implicit none      
      INTEGER    DIMDEC(phi)
      INTEGER    DIMDEC(phin)
      INTEGER    DIMDEC(type)
      REAL_T      phi(DIMV(phi))
      REAL_T     phin(DIMV(phin))
      integer    type(DIMV(type))
      INTEGER    intfacenump, intfacenumn, intfacesize
      INTEGER    intfacep(intfacesize,SDIM), intfacen(intfacesize,SDIM)
      REAL_T     dx(SDIM)
      integer    lo(SDIM), hi(SDIM)
      
!     Local variables
      INTEGER    c,d,e,i, j, k, ii, jj, kk, iii, jjj, kkk,l, m,n,p,q,r,s,t
      REAL_T     A(64,64)
      REAL_T     B(64)
      REAL_T     x,y,z
      REAL_T     FINDDIST, distance
      REAL_T     grad(SDIM)
      integer    max0
       
      if(max( abs(phi(i+1,j,k)),abs(phi(i,j+1,k)),abs(phi(i+1,j+1,k)),abs(phi(i,j,k+1)), &
         abs(phi(i+1,j+1,k+1)), abs(phi(i+1,j,k+1)), abs(phi(i,j+1,k+1))  ) .LT. BOGUS) then
         
         
         if ( phi(i,j,k)*phi(i+1,j+1,k)   .LT. 0 &
            .OR. phi(i,j,k)*phi(i+1,j,k) .LT. 0  &
            .OR. phi(i,j,k)*phi(i,j+1,k) .LT. 0  &
            .OR. phi(i,j,k)*phi(i,j,k+1) .LT. 0  &
            .OR. phi(i,j,k)*phi(i,j+1,k+1) .LT. 0  &
            .OR. phi(i,j,k)*phi(i+1,j,k+1) .LT. 0  &
            .OR. phi(i,j,k)*phi(i+1,j+1,k+1) .LT. 0  &     
         ) then

            m=0
            do ii=0,1
               do jj=0,1
                  do kk =0,1
                     x = ii*dx(1)
                     y = jj*dx(2)
                     z = kk*dx(3)
                     do n = 0,15
                       c = n/16
                       d = n/4 - 4*(n/16)
                       e = n   - 4*(n/4)
                       A(m+1,n+1) = x**c * y**d * z**e
                       A(m+2,n+1) = c * x**(max(c-1,0)) * y**d * z**e
                       A(m+3,n+1) = d * x**c * y**(max(d-1,0)) * z**e
                       A(m+4,n+1) = e * x**c * y**d * z**(max(e-1,0))
                       A(m+5,n+1) = c * d * x**(max(c-1,0)) * y**(max(d-1,0)) * z**e                    
                       A(m+6,n+1) = c * e * x**(max(c-1,0)) * y**d * z**(max(e-1,0))  
                       A(m+7,n+1) = d * e * x**c * y**(max(d-1,0)) * z**(max(e-1,0))  
                       A(m+8,n+1) = c * d * e * x**(max(c-1,0)) * y**(max(d-1,0)) * z**(max(e-1,0)) 
                     enddo
                     B(m+1) = phi(i+ii,j+jj,k+kk)
                     B(m+2) = (phi(i+ii+1,j+jj,k+kk) - phi(i+ii-1,j+jj,k+kk))/(2*dx(1))
                     B(m+3) = (phi(i+ii,j+jj+1,k+kk) - phi(i+ii,j+jj-1,k+kk))/(2*dx(2))
                     B(m+4) = (phi(i+ii,j+jj,k+kk+1) - phi(i+ii,j+jj,k+kk-1))/(2*dx(3))
                     B(m+5) = (phi(i+ii+1,j+jj+1,k+kk) - phi(i+ii-1,j+jj+1,k+kk) - phi(i+ii+1,j+jj-1,k+kk) + phi(i+ii-1,j+jj-1,k+kk))/(4*dx(1)*dx(2))
                     B(m+6) = (phi(i+ii+1,j+jj,k+kk+1) - phi(i+ii-1,j+jj,k+kk+1) - phi(i+ii+1,j+jj,k+kk-1) + phi(i+ii-1,j+jj,k+kk-1))/(4*dx(1)*dx(3))
                     B(m+7) = (phi(i+ii,j+jj+1,k+kk+1) - phi(i+ii,j+jj-1,k+kk+1) - phi(i+ii,j+jj+1,k+kk-1) + phi(i+ii,j+jj-1,k+kk-1))/(4*dx(2)*dx(3))    
                     B(m+8) = (phi(i+ii+1,j+jj+1,k+kk+1) - phi(i+ii+1,j+jj+1,k+kk-1) + phi(i+ii+1,j+jj-1,k+kk-1) - phi(i+ii+1,j+jj-1,k+kk+1)  &          
                           -  phi(i+ii-1,j+jj+1,k+kk+1) + phi(i+ii-1,j+jj+1,k+kk-1) - phi(i+ii-1,j+jj-1,k+kk-1) + phi(i+ii-1,j+jj-1,k+kk+1)) &
                              /(8*dx(1)*dx(2)*dx(3))
                     m = m + 8
  enddo
               enddo
            enddo
            
            CALL GEPP(A,B,64)
            
            do ii=0,1
               do jj=0,1
                  do kk=0,1
                  
                     iii = i + ii
                     jjj = j + jj
                     kkk = k + kk
                  
                     distance = FINDDIST(grad,B,sign(1.d0,phi(iii,jjj,kkk)),ii*dx(1),jj*dx(2),kk*dx(3),dx)              
                    
                     if (type(iii,jjj,kkk) .NE. 0 .AND. phi(iii,jjj,kkk) .GE. 0 .AND. distance .GE. 0 &
                         .AND. iii .GE. lo(1) .AND. iii .LE. hi(1) &
                         .AND. jjj .GE. lo(2) .AND. jjj .LE. hi(2) &
                         .AND. kkk .GE. lo(3) .AND. kkk .LE. hi(3)) then
                       
                        intfacenump = intfacenump + 1
                        intfacep(intfacenump,1)=iii
                        intfacep(intfacenump,2)=jjj
                        intfacep(intfacenump,3)=kkk
                     
                     
                     else if (type(iii,jjj,kkk) .NE. 0 .AND. distance .GE. 0 &
                            .AND. iii .GE. lo(1) .AND. iii .LE. hi(1) &
                            .AND. jjj .GE. lo(2) .AND. jjj .LE. hi(2) &
                            .AND. kkk. GE. lo(3) .AND. kkk .LE. hi(3))  then
                     
                        intfacenumn = intfacenumn + 1
                        intfacen(intfacenumn,1)=iii
                        intfacen(intfacenumn,2)=jjj
                        intfacen(intfacenumn,3)=kkk
                     
                     endif
                  
                     if (distance .GE. 0) then
                        type(iii,jjj,kkk) = 0
                        phin(iii,jjj,kkk) = min(abs(phin(iii,jjj,kkk)),distance)*sign(1.d0,phi(iii,jjj,kkk))
                     endif
                  enddo   
               enddo
            enddo
            
         endif
      endif
      return
      end subroutine UPDATEF 

      REAL_T FUNCTION FINDDIST(grad,B,sgn,x0,y0,z0,dx)
      implicit none
      REAL_T grad(SDIM)
      REAL_T B(64)
      REAL_T sgn
      REAL_T x0,y0,z0
      REAL_T dx(SDIM)

!     Local variables
      REAL_T t,tp
      REAL_T POLYVAL
      REAL_T DPOLYVAL
      REAL_T FA,FP
      REAL_T a,d,p
      REAL_T x,y,z
      REAL_T delta1(3), delta2(3)
      INTEGER i    
      integer ITERMAX
      parameter (ITERMAX=30)
      
      x = x0
      y = y0 
      z = z0
      i = 0
      
      do while  ( (delta1(1)**2 + delta1(2)**2 + delta1(3)**2 + delta2(1)**2 + delta2(2)**2 + delta2(3)**2)**(.5d0) &
                  .LT. 10.0**(-3.0)*dx(1)*dx(2)*dx(3) .AND. i .LT. ITERMAX)
        
        CALL GRADPVAL(B,grad,x,y,z)
        
      delta1(1) = -polyval(B,x,y,z) * grad(1)/(grad(1)**2 + grad(2)**2 + grad(3)**2)
      delta1(2) = -polyval(B,x,y,z) * grad(2)/(grad(1)**2 + grad(2)**2 + grad(3)**2)
      delta1(3) = -polyval(B,x,y,z) * grad(3)/(grad(1)**2 + grad(2)**2 + grad(3)**2)
      
      delta2(1) = (x0 - x) - grad(1)*( (x0 - x)*grad(1) + (y0 - y)*grad(2) + (z0-z)*grad(3) ) /(grad(1)**2 + grad(2)**2 + grad(3)**2)
      delta2(2) = (y0 - y) - grad(2)*( (x0 - x)*grad(1) + (y0 - y)*grad(2) + (z0-z)*grad(3) ) /(grad(1)**2 + grad(2)**2 + grad(3)**2) 
      delta2(3) = (z0 - z) - grad(3)*( (x0 - x)*grad(1) + (y0 - y)*grad(2) + (z0-z)*grad(3) ) /(grad(1)**2 + grad(2)**2 + grad(3)**2)       
      
      x = x + delta1(1) + delta2(1)
      y = y + delta1(2) + delta2(2)
      z = z + delta1(3) + delta2(3)
      
      i = i + 1
      
      enddo
      
      if (i .GE. 30) then
      
      FINDDIST = -1
      
      else
      
      FINDDIST = ( (x - x0)**2 + (y-y0)**2 +(z-z0)**2 )**(.5d0)
      
      endif      
      end FUNCTION FINDDIST  


      REAL_T FUNCTION POLYVAL(B,x,y,z)
      REAL_T B(64)
      REAL_T x,y,z
      integer c,d,e,n  
      POLYVAL=0.d0
      do n=0,63
        c = n/16
        d = n/4 - 4*(n/16)
        e = n   - 4*(n/4)
        POLYVAL = POLYVAL + B(n+1) * x**c * y**d * z**e
      enddo      
      end FUNCTION POLYVAL
      
      
      
      subroutine GRADPVAL(B,grad,x,y,z)
      implicit none
      REAL_T B(64)
      REAL_T grad(3)
      REAL_T x,y,z      

!     Local variables
      INTEGER c,d,e,n
      
      grad(1) = 0
      grad(2) = 0
      grad(3) = 0
      
      do n = 0,63
      
        c = n/16
        d = n/4 - 4*(n/16)
        e = n   - 4*(n/4)
      
      grad(1) = grad(1) + B(n+1) * c * x**(max(c-1,0)) * y**d * z**e
      grad(2) = grad(2) + B(n+1) * d * x**c * y**(max(d-1,0)) * z**e
      grad(3) = grad(3) + B(n+1) * e * x**c * y**d * z**(max(e-1,0))
      
      enddo
      return
      end subroutine GRADPVAL   
      
      subroutine GEPP(A,B,N)
      implicit none
      INTEGER N
      REAL_T  A(N,N)
      REAL_T  B(N)
      INTEGER i,j,k,p
      INTEGER ncopy
      INTEGER nrow(N)
      REAL_T  maximum
      REAL_T  X(N)
      
      do i = 1,N
        nrow(i) = i
      enddo
      
      do i = 1,N-1
         maximum = abs(A(nrow(i),i))
         p=i
         do j = i+1, N
            if (maximum < abs(A(nrow(j),i))) then
               maximum = abs(A(nrow(j),i))
               p=j
            endif
         enddo
         
         if (nrow(i) .NE. nrow(p)) then
            ncopy = nrow(i)
            nrow(i) = nrow(p)
            nrow(p) = ncopy
         endif
         
         do j = i+1, N
            do k = 1,N
               A(nrow(j),k) = A(nrow(j),k)- (A(nrow(j),i)/A(nrow(i),i))*A(nrow(i),k)
          enddo
        enddo
      enddo
        
      X(N) = B(nrow(N))/A(nrow(N),N)
      do i= N, 1, -1
        X(i) = B(nrow(i))
        do j = i+1,N
          X(i) = X(i) - A(nrow(i),j)*X(j)
        enddo
        X(i) = X(i)/A(nrow(i),i)
      enddo
      
      do i = 1, N
        B(i) = X(i)
      enddo
      end subroutine GEPP

      subroutine narrowband(type, DIMS(type), &
                           nband, nbandsize, &
                           mine, minesize, &
                           lo, hi)& bind(C,name="narrowband")
      implicit none     
      integer DIMDEC(type)
      integer  type(DIMV(type))
      integer nbandsize, minesize
      integer lo(SDIM), hi(SDIM)
      integer nband(nbandsize,SDIM)
      integer  mine(minesize,SDIM)


      integer numband
      integer nummine
      integer i,j,k
     
      numband = 0
      nummine = 0
     
      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)            
             if(type(i,j,k) .eq. 0  ) then
                numband = numband + 1
                nband(numband,1) = i
                nband(numband,2) = j
                nband(numband,3) = k
             else if (type(i,j,k) .eq. 1) then
                numband = numband + 1
                nband(numband,1) = i
                nband(numband,2) = j          
                nband(numband,3) = k
                
                nummine = nummine + 1
                mine(nummine,1) = i
                mine(nummine,2) = j
                mine(nummine,3) = k
             endif
          enddo
        enddo
      enddo
     
      nband(numband + 1,1)  = -LARGEINT
      nband(numband + 1,2)  = -LARGEINT
      nband(numband + 1,3)  = -LARGEINT      
      
    
      mine(nummine + 1,1) = -LARGEINT
      mine(nummine + 1,2) = -LARGEINT
      mine(nummine + 1,3) = -LARGEINT      
      
      end subroutine narrowband
     
 
 
 
      subroutine retypify(type, DIMS(type), nband, nbandsize) &
                          bind(C,name="retypify")

      implicit none  
      integer DIMDEC(type)
      integer type(DIMV(type))
      integer nbandsize
      integer nband(nbandsize,SDIM)
      integer i, j, k, p


      p = 1
      do while (nband(p,1) .GT. -LARGEINT)
        i = nband(p,1)
        j = nband(p,2)
        k = nband(p,3)
        p = p + 1
        type(i,j,k) = 3
      enddo
      end subroutine retypify
     
     


      subroutine fastmarch(phi, DIMS(phi), type, DIMS(type), &
                               lo, hi, dx, intfacenum, intface,  &
                               nband, nbandsize, nbandnum, mine, minenum, &
                               sgn, intfacesize,heap, heaploc) &
                               bind(C,name="fastmarch")

      implicit none
      integer     DIMS(phi)
      integer     DIMS(type)
      REAL_T      phi(DIMV(phi))
      integer     type(DIMV(type))
      integer     heaploc(DIMV(type))
      integer     lo(SDIM), hi(SDIM)
      REAL_T      dx(SDIM)
      integer     intfacenum, intfacesize
      integer     intface(intfacesize,SDIM)
      integer     nbandsize
      integer     nband(nbandsize,SDIM)
      integer     mine(nbandsize,SDIM)
      integer     nbandnum, minenum
      integer     sgn
      integer heap(nbandsize,SDIM)


#include <probdata.H>

      integer i, j, k, n, p
      integer numtent
      
      numtent = 0
      do n = 1, intfacenum
        i = intface(n,1)
        j = intface(n,2)
        k = intface(n,3)
      
        CALL UPDATE(phi,i,j,k,sgn, type, heap,numtent,DIMS(phi),DIMS(type), &
                    nbandsize, lo, hi, dx, heaploc)

        nbandnum = nbandnum + 1
        nband(nbandnum,1) = i
        nband(nbandnum,2) = j
        nband(nbandnum,3) = k
      enddo
      
      i = intface(1,1)
      j = intface(1,2)
      k = intface(1,3)
      do while (numtent .GT. 0  )

         CALL RMVNODE(heap,i,j,k,numtent, phi, DIMS(phi), &
                      nbandsize,heaploc,DIMS(type))
     
         if (abs(phi(i,j,k)) .LT. nbandwidth) then
            nbandnum = nbandnum + 1
            nband(nbandnum,1) = i
            nband(nbandnum,2) = j
            nband(nbandnum,3) = k
            type(i,j,k) = 0
         else
            type(i,j,k) = 3
            phi(i,j,k) = sign(BOGUS,phi(i,j,k))
            exit
         endif
         
         if (abs(phi(i,j,k)) .GT. mineloc &
             .AND. abs(phi(i,j,k)) .LT. nbandwidth ) then
          
            type(i,j,k) = 1
            
         endif
         
         if (abs(phi(i,j,k)) .LT. nbandwidth) then
            
            CALL UPDATE(phi,i,j,k,sgn, type, heap, numtent, DIMS(phi), DIMS(type), &
                       nbandsize, lo, hi, dx,heaploc)
            
         else
            type(i,j,k) = 3
            phi(i,j,k) = sign(BOGUS,phi(i,j,k))
         endif
         
      enddo
      
      nband(nbandnum+1,1) = -LARGEINT
      nband(nbandnum+1,2) = -LARGEINT      
      nband(nbandnum+1,3) = -LARGEINT
      
      do while (numtent .GT. 0)

         CALL RMVNODE(heap, i, j, k,numtent, phi, DIMS(phi), nbandsize, heaploc, DIMS(type))         
         type(i,j,k) =3
         phi(i,j,k) = sign(BOGUS,phi(i,j,k))

      enddo
      end subroutine fastmarch
      
      subroutine UPDATE(phi, i, j, k, sgn, type, heap, numtent, DIMS(phi), DIMS(type), &
                        nbandsize, lo, hi, dx, heaploc)
      
      implicit none
      integer i,j,k
      integer DIMDEC(phi)
      integer DIMDEC(type)
      REAL_T  phi(DIMV(phi))
      integer type(DIMV(type))
      integer nbandsize
      integer heap(nbandsize,SDIM)
      integer lo(SDIM), hi(SDIM)
      integer numtent
      integer sgn
      integer n,ii,jj,kk
      REAL_T dx(SDIM)
      integer heaploc(DIMV(type))
      
      do n = 1,6
         
         ii = i + (2*n-3)*(1-(n/3)+(n/6))
         jj = j + (2*n-7)*(n/3-(n/5)-(n/6))
         kk = k + (2*n-11)*(n/5)
         
         if ( ii .GE. lo(1) .AND. ii .LE. hi(1) &
              .AND. jj .GE. lo(2) .AND. jj .LE. hi(2) &
              .AND. kk .GE. lo(3) .AND. kk .LE. hi(3)) then
            
            if (type(ii,jj,kk) .GT. 1   .AND.   sgn*phi(ii,jj,kk) .GE. 0) then
               
               CALL EVAL(phi, ii, jj, kk, DIMS(phi), DIMS(type), lo, hi, type, sgn, dx)
               
               if (type(ii,jj,kk) .GT. 2) then
                  
                  type(ii,jj,kk) = 2
                  CALL ADDNODE(heap, ii, jj, kk, numtent, phi, DIMS(phi), lo, hi, &
                      nbandsize,heaploc,DIMS(type))
                  
               else
                  
                  CALL UPDATENODE(heap, ii, jj, kk, numtent, phi, DIMS(phi), lo, hi, & 
                      nbandsize,heaploc,DIMS(type))
                  
               endif
            endif    
         endif
      enddo
      end subroutine UPDATE
            

      subroutine EVAL(phi,i,j,k,DIMS(phi),DIMS(type),lo,hi, type,sgn,dx)
      
      implicit none
      integer i,j,k     
      integer DIMDEC(phi)
      integer DIMDEC(type)
      integer lo(SDIM), hi(SDIM)
      REAL_T phi(DIMV(phi))
      integer  type(DIMV(type))
      REAL_T a,b,c
      REAL_T dx(SDIM)
      integer sgn
      integer  left,right,up,down,front,back
      LOGICAL  lok, rok, uok, dok, fok, bok
      
      a = 0.d0
      b = 0.d0
      c = - 1
      
      left  = i - 1
      right = i + 1
      up    = j + 1
      down  = j - 1
      front = k + 1
      back  = k - 1
      
      lok  = left .GE. lo(1) .AND. sgn*phi(left,j,k) .GE. 0 .AND. type(left,j,k) .LE. 1
      rok  = right. LE. hi(1) .AND. sgn*phi(right,j,k) .GE. 0 . AND. type(right,j,k) .LE. 1 
      uok  = up .LE. hi(2) .AND. sgn*phi(i,up,k) .GE. 0  .AND. type(i,up,k) .LE. 1
      dok  = down . GE. lo(2) .AND. sgn*phi(i,down,k) .GE. 0  .AND. type(i,down,k) .LE. 1
      fok  = front .LE. hi(3) .AND. sgn*phi(i,j,front) .GE. 0  .AND. type(i,j,front) .LE. 1
      bok  = back . GE. lo(3) .AND. sgn*phi(i,j,back) .GE. 0  .AND. type(i,j,back) .LE. 1      
      
!     FIXME: The following had a sign(right,...), mistake?
      if (lok .AND. rok) then
                    
        a = a + 1/dx(1)**2
        b = b + min(sgn*phi(left,j,k),sgn*phi(right,j,k))/(dx(1)**2)
        c = c + min(phi(left,j,k)**2,phi(right,j,k)**2)/(dx(1)**2)
        
      else if (lok) then
      
        a = a + 1/dx(1)**2
        b = b + sgn*phi(left,j,k)/(dx(1)**2)
        c = c + phi(left,j,k)**2/(dx(1)**2)
      
      else if (rok ) then

        a = a + 1/dx(1)**2
        b = b + sgn*phi(right,j,k)/(dx(1)**2)
        c = c + phi(right,j,k)**2/(dx(1)**2)
        
      endif
      
      
      if (dok .AND. uok) then
      
        a = a + 1/dx(2)**2
        b = b + min(sgn*phi(i,down,k),sgn*phi(i,up,k))/(dx(2)**2)
        c = c + min(phi(i,down,k)**2,phi(i,up,k)**2)/(dx(2)**2)
        
      else if (dok) then
      
        a = a + 1/dx(2)**2
        b = b + sgn*phi(i,down,k)/(dx(2)**2)
        c = c + phi(i,down,k)**2/(dx(2)**2)
      
      else if (uok ) then

        a = a + 1/dx(2)**2
        b = b + sgn*phi(i,up,k)/(dx(2)**2)
        c = c + phi(i,up,k)**2/(dx(2)**2)
        
      endif   
      
      if (fok .AND. bok) then
      
        a = a + 1/dx(3)**2
        b = b + min(sgn*phi(i,j,front),sgn*phi(i,j,back))/(dx(3)**2)
        c = c + min(phi(i,j,front)**2,phi(i,j,back)**2)/(dx(3)**2)
        
      else if (fok) then
      
        a = a + 1/dx(2)**2
        b = b + sgn*phi(i,front,k)/(dx(3)**2)
        c = c + phi(i,front,k)**2/(dx(3)**2)
      
      else if (bok ) then

        a = a + 1/dx(2)**2
        b = b + sgn*phi(i,back,k)/(dx(3)**2)
        c = c + phi(i,back,k)**2/(dx(3)**2)
        
      endif       
        
      b = -2*b
      
      if (b**2 - 4*a*c .LT. 0) then
      
         phi(i,j,k) = sgn*(-b)/(2*a)
         return
      
      endif
      
      phi(i,j,k) = sgn*(-b + SQRT(b**2-4*a*c))/(2*a)

      end subroutine EVAL
 

      subroutine ADDNODE(heap,i,j,k,n,phi,DIMS(phi), lo, hi,  &
                        nbandsize, heaploc, DIMS(type))
      implicit none
      integer    DIMDEC(phi)
      REAL_T     phi(DIMV(phi))
      integer    nbandsize
      integer    lo(SDIM), hi(SDIM)
      integer    heap(nbandsize,SDIM)
      integer    i,j,k,n
      integer    DIMDEC(type)
      integer    heaploc(DIMV(type))

      integer index
      integer parent
      
      index = n + 1
      parent = index/2

      if (n .EQ. 0) then
         
         heap(index,1) = i
         heap(index,2) = j
         heap(index,3) = k
         heaploc(i,j,k) = index        
         
         n=n+1
         
         return
      endif
      
      do while ( ABS(phi(heap(parent,1), heap(parent,2), heap(parent,3))) .GT. ABS(phi(i,j,k)))
            
        heap(index,1) = heap(parent,1)
        heap(index,2) = heap(parent,2)
        heap(index,3) = heap(parent,3)
        heaploc(heap(index,1),heap(index,2),heap(index,3)) = index        
        index = parent
        parent = index/2
         
        if (parent .EQ. 0) exit
      enddo
      
      heap(index,1) = i
      heap(index,2) = j
      heap(index,3) = k
      heaploc(i,j,k) = index
      n = n + 1
      end subroutine ADDNODE
      
      subroutine UPDATENODE(heap,i,j,k,n,phi,DIMS(phi), lo, hi, &
                           nbandsize, heaploc,DIMS(type))
      implicit none
      integer    DIMDEC(phi)
      REAL_T     phi(DIMV(phi))
      integer    nbandsize
      integer    lo(SDIM), hi(SDIM)
      integer    heap(nbandsize,SDIM)
      integer    i,j,k, n
      integer    DIMDEC(type)
      integer    heaploc(DIMV(type))

      integer index
      integer parent
      
      index = heaploc(i,j,k)
      parent = index/2

      if (index .EQ. 1) then
        return
      endif
      
      do while ( ABS(phi(heap(parent,1 ), heap(parent,2), heap(parent,3) )) .GT. ABS(phi(i,j,k)))
            
        heap(index,1) = heap(parent,1)
        heap(index,2) = heap(parent,2)
        heap(index,3) = heap(parent,3)
        heaploc(heap(index,1),heap(index,2),heap(index,3)) = index
        index = parent
        parent = index/2
         
        if (parent .EQ. 0) exit         
      enddo
      
      heap(index,1) = i
      heap(index,2) = j
      heap(index,3) = k
      heaploc(i,j,k) = index
      end subroutine UPDATENODE
      

      subroutine RMVNODE(heap, i,j,k,n, phi, DIMS(phi), &
                        nbandsize,heaploc, DIMS(type))
      implicit none
      integer    DIMDEC(phi)
      REAL_T     phi(DIMV(phi))
      integer    nbandsize
      integer    heap(nbandsize,SDIM)
      integer    i,j,k, n
      integer    DIMDEC(type)
      integer    heaploc(DIMV(type))
      integer    index, left, right

      i=heap(1,1)
      j=heap(1,2)
      k=heap(1,3)
      
      heaploc(i,j,k) = -1
      
      index = 1
      left  = 2*index
      right = 2*index + 1
   
      do while (.TRUE.)
         
         if (left. LE. n-1) then
            
            if( ABS(phi(heap(left,1),heap(left,2),heap(left,3))) .LT. ABS(phi(heap(n,1),heap(n,2),heap(n,3)))) then
               
               if (right .LE. n-1) then
                  if (ABS(phi(heap(left,1),heap(left,2),heap(left,3))) .LT. ABS(phi(heap(right,1),heap(right,2),heap(right,3))) ) then
                     heap(index,1) = heap(left,1)
                     heap(index,2) = heap(left,2)
                     heap(index,3) = heap(left,3)
                     heaploc(heap(index,1),heap(index,2),heap(index,3)) = index
                     index = left
                     left=2*index
                     right=2*index+1
                  else
                     heap(index,1) = heap(right,1)
                     heap(index,2) = heap(right,2)
                     heap(index,3) = heap(right,3)
                     heaploc(heap(index,1),heap(index,2),heap(index,3)) = index
                     index = right
                     left=2*index
                     right=2*index+1  
                  endif
               else                  
                  heap(index,1) = heap(left,1)
                  heap(index,2) = heap(left,2)
                  heap(index,3) = heap(left,3)
                  heaploc(heap(index,1),heap(index,2),heap(index,3)) = index
                  
                  index = left
                  left=2*index
                  right=2*index+1  
               endif     
               
            else if (right .LE. n-1) then
         
               if (abs(phi(heap(right,1),heap(right,2),heap(right,3))) .LT. abs( phi(heap(n,1),heap(n,2),heap(n,3)))) then
                  heap(index,1) = heap(right,1)
                  heap(index,2) = heap(right,2)
                  heap(index,3) = heap(right,3)
                  heaploc(heap(index,1),heap(index,2),heap(index,3)) = index 
                  index = right
                  left=2*index
                  right=2*index+1
               else       
                  exit       
               endif
            else    
               exit      
            endif 
         else       
            exit      
         endif
      enddo
      
      heap(index,1) = heap(n,1)
      heap(index,2) = heap(n,2)
      heap(index,3) = heap(n,3)
      if(n .GT. 1) then
        heaploc(heap(index,1),heap(index,2),heap(index,3)) = index
      endif
      n = n - 1
      end subroutine RMVNODE

      INTEGER FUNCTION fastmarch2(phi,DIMS(phi),type,DIMS(type), &
                                  lo, hi, dx, nband, nbandsize, nbandnum, &
                                  sgn, heaploc) bind(C,name="fastmarch2")
      implicit none
      integer     DIMDEC(phi)
      integer     DIMDEC(type)
      REAL_T      phi(DIMV(phi))
      integer     type(DIMV(type))
      integer     lo(SDIM), hi(SDIM)
      REAL_T      dx(SDIM)
      integer     heaploc(DIMV(type))
      integer     nbandsize
      integer     nband(nbandsize,SDIM)
      integer     nbandnum
      integer     sgn


#include <probdata.H>
      
      integer i, j,k, n
      integer heap(nbandsize,SDIM)
      integer numtent
      
      numtent = 0
      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             heaploc(i,j,k) = -1
          enddo
        enddo
      enddo
      
      i = lo(1)-1 
      do k = lo(3), hi(3)
        do j = lo(2), hi(2)      
           if((type(i,j,k) .EQ. 0 .OR. type(i,j,k) .EQ. 1) .AND. sgn*phi(i,j,k) .GE. 0) then
              CALL UPDATE2(phi,i,j,k,sgn, type,heap,numtent,DIMS(phi), &
                         DIMS(type), nbandsize, lo, hi, dx,heaploc)
           endif
        enddo 
      enddo
      
      i = hi(1)+1 
      do k = lo(3), hi(3)
        do j = lo(2), hi(2)      
           if((type(i,j,k) .EQ. 0 .OR. type(i,j,k) .EQ. 1) .AND. sgn*phi(i,j,k) .GE. 0) then
              CALL UPDATE2(phi,i,j,k,sgn, type,heap,numtent,DIMS(phi), &
                         DIMS(type), nbandsize, lo, hi, dx,heaploc)
           endif
        enddo 
      enddo   
      
      
      j = lo(2)-1 
      do k = lo(3), hi(3)
        do i = lo(1), hi(1)      
           if((type(i,j,k) .EQ. 0 .OR. type(i,j,k) .EQ. 1) .AND. sgn*phi(i,j,k) .GE. 0) then
              CALL UPDATE2(phi,i,j,k,sgn, type,heap,numtent,DIMS(phi), &
                         DIMS(type), nbandsize, lo, hi, dx,heaploc)
           endif
        enddo 
      enddo
      
      
      j = hi(2)+1       
      do k = lo(3), hi(3)
        do i = lo(1), hi(1)      
           if((type(i,j,k) .EQ. 0 .OR. type(i,j,k) .EQ. 1) .AND. sgn*phi(i,j,k) .GE. 0) then
              CALL UPDATE2(phi,i,j,k,sgn, type,heap,numtent,DIMS(phi), &
                         DIMS(type), nbandsize, lo, hi, dx,heaploc)
           endif
        enddo 
      enddo
      
      k = lo(3)-1       
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)      
           if((type(i,j,k) .EQ. 0 .OR. type(i,j,k) .EQ. 1) .AND. sgn*phi(i,j,k) .GE. 0) then
              CALL UPDATE2(phi,i,j,k,sgn, type,heap,numtent,DIMS(phi), &
                         DIMS(type), nbandsize, lo, hi, dx,heaploc)
           endif
        enddo 
      enddo      

      k = hi(3) + 1       
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)      
           if((type(i,j,k) .EQ. 0 .OR. type(i,j,k) .EQ. 1) .AND. sgn*phi(i,j,k) .GE. 0) then
              CALL UPDATE2(phi,i,j,k,sgn, type,heap,numtent,DIMS(phi), &
                         DIMS(type), nbandsize, lo, hi, dx,heaploc)
           endif
        enddo 
      enddo      
      
      FASTMARCH2 = 0
      do while (numtent .GT. 0 )
         CALL RMVNODE(heap, i,j,k,numtent, phi, DIMS(phi), &
                      nbandsize,heaploc,DIMS(type))
         if (abs(phi(i,j,k)) .LT. nbandwidth ) then
         
            FASTMARCH2 = 1

            if (type(i,j,k) .GE. 2) then
               nbandnum = nbandnum + 1         
               nband(nbandnum,1) = i
               nband(nbandnum,2) = j            
               nband(nbandnum,3) = k                           
            endif
          
            type(i,j,k) = 0
            CALL UPDATE2(phi,i,j,k,sgn, type, heap,numtent,DIMS(phi), &
                        DIMS(type), nbandsize, lo, hi, dx,heaploc)

         else

            type(i,j,k) = 3
            phi(i,j,k) = sign(BOGUS,phi(i,j,k))
            exit

         endif
         
         if (abs(phi(i,j,k)) .GT. mineloc &
             .AND. abs(phi(i,j,k)) .LT. nbandwidth) then
          
            type(i,j,k) = 1
        endif
      enddo
      
      nband(nbandnum + 1,1) = -LARGEINT
      nband(nbandnum + 1,2) = -LARGEINT
      nband(nbandnum + 1,3) = -LARGEINT
      
      do while (numtent .GT. 0)
            
         CALL RMVNODE(heap, i,j,k,numtent, phi, DIMS(phi), &
                      nbandsize,heaploc,DIMS(type))
         type(i,j,k) = 3
         phi(i,j,k) = sign(BOGUS,phi(i,j,k))    
      enddo
      end FUNCTION fastmarch2
      
      
      
      subroutine UPDATE2(phi,i,j,k,sgn, type, heap,numtent,DIMS(phi), &
                         DIMS(type), nbandsize, lo, hi, dx, heaploc)      
      implicit none
      integer i,j,k
      integer DIMDEC(phi)
      integer DIMDEC(type)
      REAL_T  phi(DIMV(phi))
      integer type(DIMV(type))
      integer nbandsize
      integer heap(nbandsize,SDIM)
      integer lo(SDIM), hi(SDIM)
      integer numtent
      integer sgn
      integer n,ii,jj,kk
      REAL_T dx(SDIM)
      integer heaploc(DIMV(type))

      logical EVAL2
      logical isnew
      integer min

      
      do n = 1,6
         
         ii = i + (2*n-3)*(1-(n/3)+(n/6))
         jj = j + (2*n-7)*(n/3-(n/5)-(n/6))
         kk = k + (2*n-11)*(n/5)
         
              
        if (ii .GE. lo(1) .AND. ii .LE. hi(1) &
            .AND. jj .GE. lo(2) .AND. jj .LE. hi(2) &
            .AND. kk .GE. lo(3) .AND. kk .LE. hi(3)) then
        
           if (sgn*sign(1.d0,phi(ii,jj,kk)) .GE. 0 .AND. abs(phi(ii,jj,kk)) .GT. abs(phi(i,j,k)) &
               .AND. (abs(phi(ii,jj,kk)) .GE. BOGUS .OR. .NOT.(sgn*phi(ii+1,jj,kk) .LE. 0   &
               .OR. sgn*phi(ii-1,jj,kk) .LE. 0 .OR. sgn*phi(ii,jj+1,kk) .LE. 0 .OR. sgn*phi(ii,jj-1,kk) .LE. 0 &
               .OR. sgn*phi(ii+1,jj+1,kk) .LE. 0 .OR. sgn*phi(ii-1,jj-1,kk) .LE. 0 .OR. sgn*phi(ii+1,jj-1,kk) .LE. 0 &
               .OR. sgn*phi(ii-1,jj+1,kk) .LE. 0 .OR. sgn*phi(ii,jj,kk+1) .LE. 0 .OR. sgn*phi(ii,jj,kk-1) .LE. 0 &
               .OR. sgn*phi(ii+1,jj,kk+1) .LE. 0 .OR. sgn*phi(ii+1,jj,kk-1) .LE. 0 .OR. sgn*phi(ii-1,jj,kk+1) .LE. 0 &
               .OR. sgn*phi(ii-1,jj,kk-1) .LE. 0 .OR. sgn*phi(ii,jj+1,kk+1) .LE. 0 .OR. sgn*phi(ii,jj+1,kk-1) .LE. 0 &
               .OR. sgn*phi(ii,jj-1,kk+1) .LE. 0 .OR. sgn*phi(ii,jj-1,kk-1) .LE. 0 .OR. sgn*phi(ii+1,jj+1,kk+1) .LE. 0 &
               .OR. sgn*phi(ii+1,jj+1,kk-1) .LE. 0 .OR. sgn*phi(ii+1,jj-1,kk+1) .LE. 0 .OR. sgn*phi(ii+1,jj-1,kk-1) .LE. 0 &
               .OR. sgn*phi(ii-1,jj+1,kk+1) .LE. 0 .OR. sgn*phi(ii-1,jj+1,kk-1) .LE. 0 .OR. sgn*phi(ii-1,jj-1,kk+1) .LE. 0 &
               .OR. sgn*phi(ii-1,jj-1,kk-1) .LE. 0))) then
              
              isnew = EVAL2(phi,ii,jj,kk,DIMS(phi), DIMS(type), &
                            lo,hi, type,sgn,dx, phi(i,j,k))
           
              isnew = isnew .OR. EVAL2(phi,ii,jj,kk,DIMS(phi), DIMS(type), &
                                      lo,hi, type,sgn,dx, phi(ii,jj,kk))
              
              if (isnew ) then
                 type(ii,jj,kk) = min(2,type(ii,jj,kk))
                 if(heaploc(ii,jj,kk) .EQ. -1) then
                    CALL ADDNODE(heap,ii,jj, kk, numtent,phi,DIMS(phi), &
                                lo, hi, nbandsize, heaploc,DIMS(type))
                 else
                    CALL UPDATENODE(heap,ii,jj, kk, numtent,phi,DIMS(phi), &
                                   lo, hi, nbandsize, heaploc,DIMS(type))   
                 endif
              endif
           endif    
        endif        
      enddo
      end subroutine UPDATE2
      
    
      LOGICAL FUNCTION EVAL2(phi,i,j,k, DIMS(phi), DIMS(type),lo,hi, type,sgn,dx,phisrc)      
      implicit none
      integer i,j,k     
      integer DIMDEC(phi)
      integer DIMDEC(type)
      integer lo(SDIM), hi(SDIM)
      REAL_T   phi(DIMV(phi))
      integer type(DIMV(type))
      REAL_T a,b,c
      REAL_T dx(SDIM)
      integer sgn
      integer  left,right,up,down,front, back
      REAL_T phisrc
      LOGICAL  lok, rok, uok, dok, fok, bok
      
      a = 0.d0
      b = 0.d0
      c = - 1
      
      left  = i - 1
      right = i + 1
      up    = j + 1
      down  = j - 1
      front = k + 1
      back  = k - 1
      
      lok  = sgn*phi(left,j,k) .GE. 0 .AND. type(left,j,k) .LE. 1 .AND. abs(phi(left,j,k)) .LE. abs(phisrc)
      rok  = sgn*phi(right,j,k) .GE. 0 . AND. type(right,j,k) .LE. 1 .AND. abs(phi(right,j,k)) .LE. abs(phisrc)
      uok  = sgn*phi(i,up,k) .GE. 0  .AND. type(i,up,k) .LE. 1 .AND. abs(phi(i,up,k)) .LE. abs(phisrc)
      dok  = sgn*phi(i,down,k) .GE. 0  .AND. type(i,down,k) .LE. 1 .AND. abs(phi(i,down,k)) .LE. abs(phisrc)
      fok  = sgn*phi(i,j,front) .GE. 0  .AND. type(i,j,front) .LE. 1 .AND. abs(phi(i,j,front)) .LE. abs(phisrc)
      bok  = sgn*phi(i,j,back) .GE. 0  .AND. type(i,j,back) .LE. 1 .AND. abs(phi(i,j,back)) .LE. abs(phisrc)         


      if (left .LT. lo(1)-1) then
         
         if (rok )then
        
            a = a + 1/dx(1)**2
            b = b + sgn*phi(right,j,k)/(dx(1)**2)
            c = c + phi(right,j,k)**2/(dx(1)**2)        
            
         endif
         
      else if (right .GT. hi(1) +1) then
      
         if (lok ) then
        
            a = a + 1/dx(1)**2
            b = b + sgn*phi(left,j,k)/(dx(1)**2)
            c = c + phi(left,j,k)**2/(dx(1)**2)        
            
         endif
         
      else
   
         if (lok .AND. rok ) then
      
            a = a + 1/dx(1)**2        
            b = b + min(sgn*phi(left,j,k),sgn*phi(right,j,k))/(dx(1)**2)
            c = c + min(phi(left,j,k)**2,phi(right,j,k)**2)/(dx(1)**2)
            
         else if (lok) then
      
            a = a + 1/dx(1)**2
            b = b + sgn*phi(left,j,k)/(dx(1)**2)
            c = c + phi(left,j,k)**2/(dx(1)**2)
      
         else if (rok ) then
 
            a = a + 1/dx(1)**2
            b = b + sgn*phi(right,j,k)/(dx(1)**2)
            c = c + phi(right,j,k)**2/(dx(1)**2)
        
         endif
         
      endif


      if (down . LT. lo(2)-1) then
         
         if (uok) then
          
            a = a + 1/dx(2)**2
            b = b + sgn*phi(i,up,k)/(dx(2)**2)
            c = c + phi(i,up,k)**2/(dx(2)**2)
            
         endif
         
      else if (up . GT. hi(2) + 1) then
         
         if (dok ) then
            
            a = a + 1/dx(2)**2
            b = b + sgn*phi(i,down,k)/(dx(2)**2)
            c = c + phi(i,down,k)**2/(dx(2)**2) 
            
         endif
        
      else
         
         if ( uok .AND. dok ) then
            
            a = a + 1/dx(2)**2
            b = b + min(sgn*phi(i,down,k),sgn*phi(i,up,k))/(dx(2)**2)
            c = c + min(phi(i,down,k)**2,phi(i,up,k)**2)/(dx(2)**2)
            
        else if (dok) then
           
           a = a + 1/dx(2)**2
           b = b + sgn*phi(i,down,k)/(dx(2)**2)
           c = c + phi(i,down,k)**2/(dx(2)**2)
           
        else if (uok ) then
           
           a = a + 1/dx(2)**2
           b = b + sgn*phi(i,up,k)/(dx(2)**2)
           c = c + phi(i,up,k)**2/(dx(2)**2)
           
        endif      
        
      endif
      
      
      if (back . LT. lo(2)-1) then
         
         if (fok) then
          
            a = a + 1/dx(3)**2
            b = b + sgn*phi(i,j,front)/(dx(3)**2)
            c = c + phi(i,j,front)**2/(dx(3)**2)
            
         endif
         
      else if (front . GT. hi(2) + 1) then
         
         if (bok ) then
            
            a = a + 1/dx(3)**2
            b = b + sgn*phi(i,j,back)/(dx(3)**2)
            c = c + phi(i,j,back)**2/(dx(3)**2) 
            
         endif
        
      else
         
         if ( fok .AND. bok ) then
            
            a = a + 1/dx(3)**2
            b = b + min(sgn*phi(i,j,back),sgn*phi(i,j,front))/(dx(3)**2)
            c = c + min(phi(i,j,back)**2,phi(i,j,front)**2)/(dx(3)**2)
            
        else if (bok) then
           
           a = a + 1/dx(3)**2
           b = b + sgn*phi(i,j,back)/(dx(3)**2)
           c = c + phi(i,j,back)**2/(dx(3)**2)
           
        else if (fok ) then
           
           a = a + 1/dx(3)**2
           b = b + sgn*phi(i,j,front)/(dx(3)**2)
           c = c + phi(i,j,front)**2/(dx(3)**2)
           
        endif      
        
      endif      
      
      b = -2*b
      if (a .EQ. 0.d0) then

        EVAL2 = .FALSE.
        return
       
      endif
      
      EVAL2 = .FALSE.
      
      if (b**2 - 4*a*c .LT. 0.d0) then
         if (ABS(phi(i,j,k)) .GT. (-b)/(2*a) + 1.d-10) then
            phi(i,j,k) = sgn*(-b)/(2*a)
            EVAL2 = .TRUE.
            return
         endif
      endif
      
      if (abs(phi(i,j,k)) .GT. (-b+SQRT(b**2-4*a*c))/(2*a)+1.d-10) then
        phi(i,j,k) = sgn*(-b+SQRT(b**2-4*a*c))/(2*a)
        EVAL2 = .TRUE.
      endif
      
      end FUNCTION EVAL2
      
      subroutine mine(type, DIMS(type), nband, nbandsize, &
                       mine, minesize,lo, hi) bind(C,name="mine") 
      implicit none     
      integer DIMDEC(type)
      integer nbandsize, minesize
      integer lo(SDIM), hi(SDIM)
      integer  type(DIMV(type))
      integer  nband(nbandsize,SDIM)
      integer  mine(minesize,SDIM)


      integer i, j,k, p
      integer nummine
      
      nummine = 0
      p = 1
      do while(nband(p,1) .GT. -LARGEINT)
      
        i = nband(p,1)
        j = nband(p,2)
        k = nband(p,3)
        p = p + 1
        
        if (type(i,j,k) .EQ. 1) then
          nummine = nummine +1
          mine(nummine,1) = i
          mine(nummine,2) = j
          mine(nummine,3) = k
        endif
        
      enddo
      
      mine(nummine+1,1) = -LARGEINT
      mine(nummine+1,2) = -LARGEINT      
      mine(nummine+1,3) = -LARGEINT
      end subroutine mine
           
     

     
      subroutine nbandnumify(nband, nbandsize,nbandnum) &
                 bind(C,name="nbandnumify")
      implicit none      
      integer nbandsize, nbandnum
      integer nband(nbandsize,SDIM)


      integer p
      
      p = 0
      do while(nband(p+1,1) .GT. -LARGEINT)
        p = p + 1
      enddo
      nbandnum = p
      end subroutine nbandnumify
      
end module LS_3d_module
