#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif
  
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <LS_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 2
#define BOGUS 1.d30
#define LARGEINT 100000000


module LS_2d_module
  
  implicit none

  private 

  public :: phiupd, SWITCH, lscfl, findintrfce, UPDATEF, FINDDIST, POLYVAL, &
         GRADPVAL, DPOLYVAL, GEPP, narrowband, retypify, fastmarch, &
         UPDATE, EVAL, ADDNODE, UPDATENODE, RMVNODE, fastmarch2, &
         UPDATE2, EVAL2, mine, nbandnumify 
       
contains

      INTEGER FUNCTION phiupd(phi, DIMS(phi), phin, DIMS(phin),&
                              uadv, DIMS(uadv), vadv, DIMS(vadv),&
                              nband, nbandsize, mine, minesize,&
                              lo, hi, dt, dx, type,DIMS(type)) &
                              bind(C,name="phiupd")

      implicit none

      INTEGER  DIMDEC(phi)
      INTEGER  DIMDEC(phin)
      INTEGER  DIMDEC(uadv)
      INTEGER  DIMDEC(vadv)
      REAL_T    phi(DIMV(phi))
      REAL_T   phin(DIMV(phin))
      REAL_T   uadv(DIMV(uadv))
      REAL_T   vadv(DIMV(vadv))
      integer  nbandsize, minesize
      integer  nband(nbandsize,SDIM)
      integer  mine(minesize,SDIM)
      INTEGER  lo(SDIM), hi(SDIM)
      REAL_T   dt
      REAL_T   dx(SDIM)
      integer  DIMDEC(type)
      integer  type(DIMV(type))
#include <probdata.H>

      integer  i, j, p,s,t
      REAL_T   Dxm, Dxp
      REAL_T   Dym, Dyp
      REAL_T   phix, phiy, phixx, phiyy, phixy, kappa
      REAL_T   uavg, vavg, Fo, Fkappa, Fadv
      REAL_T   SWITCH
      REAL_T   Dxpp, Dxmm, Dxpm, Dypp, Dymm, Dypm
      REAL_T   Ad, Bd, Cd, Dd, Delp, Delm

      p=1
      do while (nband(p,1) .GT. -LARGEINT)
      
        i = nband(p,1)
        j = nband(p,2)
        p = p + 1
        
        if (max(type(i+1,j), type(i-1,j), type(i,j+1), type(i,j-1),&
            type(i+1,j+1), type(i-1,j-1), type(i+1,j-1), type(i-1,j+1)) .LE. 1 ) then

           Dxm   = ( phi(i,j) - phi(i-1,j) ) / dx(1)
           Dxp   = ( phi(i+1,j) - phi(i,j) ) / dx(1)
           Dym   = ( phi(i,j) - phi(i,j-1) ) / dx(2)
           Dyp   = ( phi(i,j+1) - phi(i,j) ) / dx(2)
        
           phix  = ( phi(i+1,j) - phi(i-1,j) ) / (2*dx(1))
           phiy  = ( phi(i,j+1) - phi(i,j-1) ) / (2*dx(2))
        
           phixx = ( phi(i+1,j) - 2*phi(i,j) + phi(i-1,j) ) &
               / (dx(1)**2)
           phiyy = ( phi(i,j+1) - 2*phi(i,j) + phi(i,j-1) )&
               / (dx(2)**2)
           phixy = ( phi(i+1,j+1) + phi(i-1,j-1)&
               - phi(i+1,j-1) - phi(i-1,j+1) )&
               / (4*dx(1)*dx(2))
        
           if (phix**2 + phiy**2 .GT. 0 ) then
              Fkappa = -kapb*( phixx*phiy**2 - 2*phiy*phix*phixy &
                  + phiyy*phix**2 ) / (phix**2 + phiy**2) 
           else
              Fkappa = 0
           endif
        

           if (LSorder .EQ. 2) then
              
              Dxpp = 0
              if (i+2 .LE. hi(1) +1) then
                 if(type(i+2,j) .LT. 2) then
                    Dxpp = (phi(i+2,j) - 2*phi(i+1,j) + phi(i,j))/dx(1)
                 endif
              endif
              
              Dxmm = 0
              if (i-2 .GE. lo(1)) then
                 if (type(i-2,j) .LT. 2) then
                    Dxmm = (phi(i,j) - 2* phi(i-1,j) + phi(i-2,j))/dx(1)
                 endif
              endif
       
              Dxpm = (phi(i+1,j) - 2* phi(i,j) + phi(i-1,j))/dx(1)
              
              Dypp = 0
              if(j+2 .LE. hi(2)) then
                 if (type(i,j+2) .LT. 2) then
                    Dypp = (phi(i,j+2) - 2* phi(i,j+1) + phi(i,j))/dx(2)
                 endif        
              endif
              
              Dymm = 0
              if(j-2 .GE. lo(2)) then
                 if(type(i,j-2) .LT. 2) then
                    Dymm = (phi(i,j) - 2* phi(i,j-1) + phi(i,j-2))/dx(2)
                 endif
              endif
              
              Dypm = (phi(i,j+1) - 2* phi(i,j) + phi(i,j-1))/dx(2)
        

              Ad = Dxm + .5d0*SWITCH(Dxmm,Dxpm)
              Bd = Dxp + .5d0*SWITCH(Dxpp,Dxpm)
              Cd = Dym + .5d0*SWITCH(Dymm,Dypm)
              Dd = Dyp + .5d0*SWITCH(Dypp,Dypm)
              
              Delp = (max(Ad,0.d0)**2 + min(Bd,0.d0)**2 + max(Cd,0.d0)**2 + min(Dd,0.d0)**2)**(.5)
              
              Delm = (max(Bd,0.d0)**2 + min(Ad,0.d0)**2 + max(Dd,0.d0)**2 + min(Cd,0.d0)**2)**(.5)
              
           endif
        
           Fo = kapa
           if (LSorder .EQ. 1) then
              if (Fo .GT. 0) then
                 Fo = Fo*( ( max(Dxm,0.d0) + min(Dxp,0.d0) )**2 &
                     + ( max(Dym,0.d0) + min(Dyp,0.d0) )**2  )**(.5d0) 
              else
                 Fo = Fo*( ( min(Dxm,0.d0) + max(Dxp,0.d0) )**2 &
                     + ( min(Dym,0.d0) + max(Dyp,0.d0) )**2  )**(.5d0) 
              endif
           else if (LSorder .EQ. 2) then
              if (Fo .GT. 0) then
                 Fo = Fo*Delp
              else
                 Fo = Fo*Delm
              endif
           endif
           
           Fadv = 0
           uavg = ( uadv(i,j) + uadv(i+1,j) ) * .5d0
           vavg = ( vadv(i,j) + vadv(i,j+1) ) * .5d0
           
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
              
           else if (LSorder .EQ. 2) then
              Fadv = uavg*phix + vavg*phiy
           endif
           
           phin(i,j) = phi(i,j) - dt*( Fo + Fkappa + Fadv )
           
        endif
      enddo
      
      PHIUPD = 0
      
      p = 1
      
      do while (mine(p,1) .GT. -LARGEINT)
         i = mine(p,1)
         j = mine(p,2)
         p = p + 1
         
         if (sign(1.d0,phi(i,j))*sign(1.d0,phin(i,j)) .LE. 0) then
            PHIUPD = 1
            exit
         endif 
      enddo
      
      return 
    end FUNCTION phiupd



      REAL_T FUNCTION SWITCH(x,y)
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
    end FUNCTION SWITCH



    REAL_T FUNCTION lscfl(phi, DIMS(phi), uadv, DIMS(uadv), vadv,&
                          DIMS(vadv),nband, nbandsize, mine, minesize,&
                          lo, hi, phit, dx, type, DIMS(type)) &
                          bind(C,name="lscfl")

      implicit none
      INTEGER  DIMDEC(phi)
      INTEGER  DIMDEC(uadv)
      INTEGER  DIMDEC(vadv)
      REAL_T    phi(DIMV(phi))
      REAL_T   uadv(DIMV(uadv))
      REAL_T   vadv(DIMV(vadv))
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

      integer i, j, p
      REAL_T phix, phiy, phixx, phiyy, phixy, kappa, speed
      REAL_T phidt

      phidt = phit
      
      p = 1
      do while (nband(p,1) .GT. -LARGEINT)
         
         i = nband(p,1)
         j = nband(p,2)
         p = p + 1
         
         if(max(type(i+1,j),type(i-1,j),type(i,j+1),type(i,j-1),&
             type(i+1,j+1),type(i-1,j-1), type(i+1,j-1), type(i-1,j+1)) .LE. 1 ) then    
          
            phix = (phi(i+1,j)-phi(i-1,j))/(2*dx(1))
            phiy = (phi(i,j+1)-phi(i,j-1))/(2*dx(2))
            phixx = (phi(i+1,j)-2*phi(i,j)+phi(i-1,j))/(dx(1)**2)
            phiyy = ( phi(i,j+1) - 2*phi(i,j) + phi(i,j-1) )/(dx(2)**2)
            phixy = (phi(i+1,j+1) + phi(i-1,j-1) - phi(i+1,j-1) - phi(i-1,j+1))/(4*dx(1)*dx(2))
            if (phix**2 + phiy**2 .GT. 0 ) then
               kappa = (phixx*phiy**2 - 2*phiy*phix*phixy + phiyy*phix**2)/((phix**2+phiy**2)**(3.d0/2.d0))
            else
               kappa = 0.d0
            endif
            
            speed = ( (.5d0*( uadv(i,j) + uadv(i+1,j) ) )**2 + ( .5d0*( vadv(i,j) +vadv(i,j+1) ) )**2)**.5d0&
                + abs(kapa-kapb*kappa)
             
            phidt = min( phit,&
                        1/(max(1.d0/phidt,speed/(.8d0*min(dx(1),dx(2))), 4*abs(kapb)/(.8d0*min(dx(1),dx(2))**2) ))  )     

         endif
      enddo
        
      lscfl = phidt
    end FUNCTION lscfl
      


      subroutine findintrfce( phi, DIMS(phi), phin, DIMS(phin), type, DIMS(type),&
           lo, hi, dx, intfacenump, intfacenumn, intfacep,intfacen,&
           nband, nbandsize, intfacesize) & bind(C,name="findintrfce")

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
      
!c     Local variables
      INTEGER    i, j, r
      
      intfacenump=0
      intfacenumn=0
      
      r = 1
      do while(nband(r,1) .GT. -LARGEINT)
        i = nband(r,1)
        j = nband(r,2)
        r = r + 1      
        
        phin(i,j) = sign(BOGUS,phi(i,j))
      enddo
      
      
      r = 1
      do while(nband(r,1) .GT. -LARGEINT)
        i = nband(r,1)
        j = nband(r,2)
        r = r + 1
        
        call UPDATEF(i, j, phi, DIMS(phi), phin, DIMS(phin), type, DIMS(type),&
                    intfacep, intfacen, intfacenump, intfacenumn,&
                    dx, intfacesize, lo, hi) 
           
        if( (i .EQ. lo(1)) .OR. (j .EQ. lo(2)) ) then
           
           if(i .EQ. lo(1)) then
                call UPDATEF(i-1, j, phi, DIMS(phi), phin, DIMS(phin), type, DIMS(type),&
                            intfacep, intfacen, intfacenump, intfacenumn, &
                            dx, intfacesize, lo, hi )
             endif
             
             if(j. EQ. lo(2)) then 
              call UPDATEF(i, j-1, phi, DIMS(phi), phin, DIMS(phin), type, DIMS(type),&
                         intfacep, intfacen, intfacenump, intfacenumn, &
                         dx, intfacesize, lo, hi )
             endif
             
             if (i .EQ. lo(1) .AND. j .EQ. lo(2)) then
               call UPDATEF(i-1, j-1, phi, DIMS(phi), phin, DIMS(phin), type, DIMS(type),&
                          intfacep, intfacen, intfacenump, intfacenumn, &
                          dx, intfacesize, lo, hi )
             endif
          endif
      enddo
      
      nband(1,1) = -LARGEINT
      nband(1,2) = -LARGEINT
    end subroutine findintrfce


      subroutine UPDATEF(i,j, phi, DIMS(phi), phin, DIMS(phin), type, DIMS(type),&
                        intfacep, intfacen, intfacenump, intfacenumn, &
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
      
!c     Local variables
      INTEGER    c,d,i, j, ii, jj, iii, jjj, k,l, m,n,p,q,r,s,t
      REAL_T     A(16,16)
      REAL_T     B(16)
      REAL_T     x,y
      REAL_T     FINDDIST, distance
      REAL_T     grad(SDIM)
      integer    max0
      integer    nrow(16)
      REAL_T     geppX(16)
      
       
      if(max( abs(phi(i+1,j)),abs(phi(i,j+1)),abs(phi(i+1,j+1))) .LT. BOGUS) then
         
         
         if ( phi(i,j)*phi(i+1,j+1)   .LT. 0&
            .OR. phi(i,j)*phi(i+1,j) .LT. 0 &
            .OR. phi(i,j)*phi(i,j+1) .LT. 0 ) then


            m=0
            do ii=0,1
               do jj=0,1
                  x = ii*dx(1)
                  y = jj*dx(2)
                  do n = 0,15
                    c = n/4
                    d = n-4*(n/4)
                    A(m+1,n+1) = x**(c) * y**(d)
                    A(m+2,n+1) = c * x**(max(c-1,0)) * y**(d)
                    A(m+3,n+1) = d * x**(c) * y**(max(d-1,0))
                    A(m+4,n+1) = c * d * x**(max(c-1,0)) * y**(max(d-1,0))                    
                  enddo
                  B(m+1) = phi(i+ii,j+jj)
                  B(m+2) = (phi(i+ii+1,j+jj) - phi(i+ii-1,j+jj))/(2*dx(1))
                  B(m+3) = (phi(i+ii,j+jj+1) - phi(i+ii,j+jj-1))/(2*dx(2))
                  B(m+4) = (phi(i+ii+1,j+jj+1) - phi(i+ii-1,j+jj+1) - phi(i+ii+1,j+jj-1) + phi(i+ii-1,j+jj-1))/(4*dx(1)*dx(2))
                  m = m + 4
               enddo
            enddo
            
            CALL GEPP(A,B,16,geppX,nrow)
            
            do ii=0,1
               do jj=0,1
                  
                  iii = i + ii
                  jjj = j + jj
                                    
                  distance = FINDDIST(grad,B,sign(1.d0,phi(iii,jjj)),ii*dx(1),jj*dx(2),dx)              
                  
                  if (type(iii,jjj) .NE. 0 .AND. phi(iii,jjj) .GE. 0 .AND. distance .GE. 0&
                      .AND. iii .GE. lo(1) .AND. iii .LE. hi(1)&
                      .AND. jjj .GE. lo(2) .AND. jjj .LE. hi(2)) then
                     
                     intfacenump = intfacenump + 1
                     intfacep(intfacenump,1)=iii
                     intfacep(intfacenump,2)=jjj
                     
                     
                  else if (type(iii,jjj) .NE. 0 .AND. distance .GE. 0&
                         .AND. iii .GE. lo(1) .AND. iii .LE. hi(1)&
                         .AND. jjj .GE. lo(2) .AND. jjj .LE. hi(2))  then
                     
                     intfacenumn = intfacenumn + 1
                     intfacen(intfacenumn,1)=iii
                     intfacen(intfacenumn,2)=jjj
                     
                  endif
                  
                  if (distance .GE. 0) then
                     type(iii,jjj) = 0
                     phin(iii,jjj) = min(abs(phin(iii,jjj)),distance)*sign(1.d0,phi(iii,jjj))
                  endif
                  
               enddo
            enddo
            
         endif
      endif
      return
    end subroutine UPDATEF
      


      REAL_T FUNCTION FINDDIST(grad,B,sgn,x0,y0,dx)
      implicit none
      REAL_T grad(SDIM)
      REAL_T B(16)
      REAL_T sgn
      REAL_T x0,y0
      REAL_T dx(SDIM)

!c     Local variables
      REAL_T t,tp
      REAL_T POLYVAL
      REAL_T DPOLYVAL
      REAL_T FA,FP
      REAL_T a,d,p
      REAL_T x,y
      REAL_T delta1(2), delta2(2)
      INTEGER i    
      integer ITERMAX
      parameter (ITERMAX=30)
      
      x = x0
      y = y0 
      
      i = 0
      
      delta1(1) = BOGUS
      delta1(2) = BOGUS
      delta2(1) = BOGUS
      delta2(2) = BOGUS
      
      do while  ( (delta1(1)**2 + delta1(2)**2 + delta2(1)**2 + delta2(2)**2)**(.5) .GT. 10.0**(-6.0)*dx(1)*dx(2) .AND. i .LT. ITERMAX)
        
        CALL GRADPVAL(B,grad,x,y)
        
      	delta1(1) = -polyval(B,x,y) * grad(1)/(grad(1)**2 + grad(2)**2)
      	delta1(2) = -polyval(B,x,y) * grad(2)/(grad(1)**2 + grad(2)**2)
      	
      	delta2(1) = (x0 - x) - grad(1)*( (x0 - x)*grad(1) + (y0 - y)*grad(2) ) /(grad(1)**2 + grad(2)**2)
      	delta2(2) = (y0 - y) - grad(2)*( (x0 - x)*grad(1) + (y0 - y)*grad(2) ) /(grad(1)**2 + grad(2)**2) 
      	
      	x = x + delta1(1) + delta2(1)
      	y = y + delta1(2) + delta2(2)
      	
      	i = i + 1
      
      enddo
      
      if (i .GE. 30 .OR. x .LT. 0 .OR. x .GT. dx(1) .OR. y .LT. 0 .OR. y .GT. dx(2)) then
      
      	FINDDIST = -1
      	
      else
      
      	FINDDIST = ( (x - x0)**2 + (y-y0)**2 )**(.5)
      
      endif
      

    end FUNCTION FINDDIST
      
      
      REAL_T FUNCTION POLYVAL(B,x,y)
      REAL_T B(16)
      REAL_T x,y
      integer c,d,n  
      POLYVAL=0.d0
      do n=0,15
        c = n/4
        d = n-4*(n/4)
        POLYVAL = POLYVAL + B(n+1)*x**(c)*y**(d)
      enddo
      
    end FUNCTION POLYVAL
      
      
      subroutine GRADPVAL(B,grad,x,y)

      implicit none

      REAL_T B(16)
      REAL_T grad(2)
      REAL_T x,y      

!c     Local variables
      INTEGER c,d,n
      
      grad(1) = 0
      grad(2) = 0
      
      do n = 0,15
      
        c = n/4
        d = n-4*(n/4)
      
        grad(1) = grad(1) + B(n+1) * c * x**(max(c-1,0)) * y**(d)
        grad(2) = grad(2) + B(n+1) * d * x**(c) * y**(max(d-1,0))
      
      enddo

      return
    end subroutine GRADPVAL
      
      REAL_T FUNCTION DPOLYVAL(B,grad,x,y)
      implicit none
      REAL_T B(16)
      REAL_T grad(2)
      REAL_T x,y      
      INTEGER n
      DPOLYVAL= (3*x**2*(B(1)*y**3 + B(2)*y**2 + B(3)*y + B(4)) + 2*x*(B(5)*y**3 + B(6)*y*2* &
                + B(7)*y + B(8))  + (B(9)*y**3 + B(10)*y**2 + B(11)*y + B(12)) )*grad(1)  &
                + (3*y**2*(B(1)*x**3 + B(5)*x**2 + B(9)*x + B(13)) + 2*y*(B(2)*x**3 &
                + B(6)*x**2 + B(10)*x + B(14))  + (B(3)*x**3 + B(7)*x**2 + B(11)*x + B(15)) )*grad(2)
   end FUNCTION DPOLYVAL
      
      
      subroutine GEPP(A,B,N,X,nrow)
      implicit none
      INTEGER N
      REAL_T  A(N,N)
      REAL_T  B(N)
      INTEGER i,j,k,p
      INTEGER ncopy
      INTEGER nrow(N)
      REAL_T  maximum
      REAL_T  X(N)
      REAL_T  m
      
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
            m = A(nrow(j),i)/A(nrow(i),i)
            do k = 1,N
               A(nrow(j),k) = A(nrow(j),k)- m*A(nrow(i),k)        
            enddo
            B(nrow(j)) = B(nrow(j))- m*B(nrow(i))
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
      

    subroutine narrowband(type, DIMS(type),&
                          nband, nbandsize,&
                          mine, minesize,&
                          lo, hi) & bind(C,name="narrowband")
      implicit none     
     
      integer DIMDEC(type)
      integer  type(DIMV(type))
      integer nbandsize, minesize
      integer lo(SDIM), hi(SDIM)
      integer nband(nbandsize,SDIM)
      integer  mine(minesize,SDIM)
     
      integer numband
      integer nummine
      integer i,j
     
      numband = 0
      nummine = 0
 
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            
            if(type(i,j) .eq. 0  ) then
               numband = numband + 1
               nband(numband,1) = i
               nband(numband,2) = j
            else if (type(i,j) .eq. 1) then
               numband = numband + 1
               nband(numband,1) = i
               nband(numband,2) = j          
               
               nummine = nummine + 1
               mine(nummine,1) = i
               mine(nummine,2) = j
            endif
         enddo
      enddo
     
      nband(numband + 1,1)  = -LARGEINT
      nband(numband + 1,2)  = -LARGEINT
    
      mine(nummine + 1,1) = -LARGEINT
      mine(nummine + 1,2) = -LARGEINT
      
    end subroutine narrowband
     
 
 
 
    subroutine retypify(type, DIMS(type), nband, nbandsize) &
         bind(C,name="retypify")
      implicit none  
      integer DIMDEC(type)
      integer type(DIMV(type))
      integer nbandsize
      integer nband(nbandsize,SDIM)
      integer i, j, p

      p = 1
      do while (nband(p,1) .GT. -LARGEINT)
        i = nband(p,1)
        j = nband(p,2)
        p = p + 1
        type(i,j) = 3
      enddo
    end subroutine retypify
     
     


      subroutine fastmarch(phi, DIMS(phi), type, DIMS(type),&
                           lo, hi, dx, intfacenum, intface, &
                           nband, nbandsize, nbandnum, mine,&
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
      integer     nbandnum
      integer     sgn
      integer heap(nbandsize,SDIM)
#include <probdata.H>

      integer i, j, n, p
      integer numtent
      
      numtent = 0
      do n = 1, intfacenum
        i = intface(n,1)
        j = intface(n,2)
      
        CALL UPDATE(phi,i,j,sgn, type, heap,numtent,DIMS(phi),DIMS(type),&
                   nbandsize, lo, hi, dx, heaploc)

        nbandnum = nbandnum + 1
        nband(nbandnum,1) = i
        nband(nbandnum,2) = j
      enddo
      
      do while (numtent .GT. 0  )

         CALL RMVNODE(heap, i,j,numtent, phi, DIMS(phi),&
                     nbandsize,heaploc,DIMS(type))
     
         if (abs(phi(i,j)) .LT. nbandwidth) then
            nbandnum = nbandnum + 1
            nband(nbandnum,1) = i
            nband(nbandnum,2) = j
            type(i,j) = 0
         else
            type(i,j) = 3
            phi(i,j) = sign(BOGUS,phi(i,j))
            exit
         endif
         
         if (abs(phi(i,j)) .GT. mineloc&
             .AND. abs(phi(i,j)) .LT. nbandwidth ) then
          
            type(i,j) = 1
            
         endif
         
         CALL UPDATE(phi,i,j,sgn, type, heap, numtent, DIMS(phi), DIMS(type),&
                    nbandsize, lo, hi, dx,heaploc)
         
      enddo
      
      nband(nbandnum+1,1) = -LARGEINT
      nband(nbandnum+1,2) = -LARGEINT      
      do while (numtent .GT. 0)

         CALL RMVNODE(heap, i, j, numtent, phi, DIMS(phi), nbandsize, heaploc, DIMS(type))         
         type(i,j) =3
         phi(i,j) = sign(BOGUS,phi(i,j))

      enddo
    end subroutine fastmarch



      subroutine UPDATE(phi, i, j, sgn, type, heap, numtent, DIMS(phi), DIMS(type),&
                       nbandsize, lo, hi, dx, heaploc)
      
      implicit none
      integer i,j
      integer DIMDEC(phi)
      integer DIMDEC(type)
      REAL_T  phi(DIMV(phi))
      integer type(DIMV(type))
      integer nbandsize
      integer heap(nbandsize,SDIM)
      integer lo(SDIM), hi(SDIM)
      integer numtent
      integer sgn
      integer n,ii,jj
      REAL_T dx(SDIM)
      integer heaploc(DIMV(type))
      
      do n = 1,4
         
         ii = i - 1 + n/2
         jj = j + 2*(n/3) - n/2
         
         if ( ii .GE. lo(1) .AND. ii .LE. hi(1)&
             .AND. jj .GE. lo(2) .AND. jj .LE. hi(2)) then
            
            if (type(ii,jj) .GT. 1   .AND.   sgn*phi(ii,jj) .GE. 0) then
               
               CALL EVAL(phi, ii, jj, DIMS(phi), DIMS(type), lo, hi, type, sgn, dx)
               
               if (type(ii,jj) .GT. 2) then
                  
                  type(ii,jj) = 2
                  CALL ADDNODE(heap, ii, jj, numtent, phi, DIMS(phi), lo, hi, &
                      nbandsize,heaploc,DIMS(type))
                  
               else
                  
                  CALL UPDATENODE(heap, ii, jj, numtent, phi, DIMS(phi), lo, hi, &
                      nbandsize,heaploc,DIMS(type))
                  
               endif
            endif    
         endif
      enddo
    end subroutine UPDATE
      
      
      
      subroutine EVAL(phi,i,j,DIMS(phi),DIMS(type),lo,hi, type,sgn,dx)
      
      implicit none
      integer i,j     
      integer DIMDEC(phi)
      integer DIMDEC(type)
      integer lo(SDIM), hi(SDIM)
      REAL_T phi(DIMV(phi))
      integer  type(DIMV(type))
      REAL_T a,b,c
      REAL_T dx(SDIM)
      integer sgn
      integer  left,right,up,down
      LOGICAL  lok, rok, uok, dok
      
      a = 0.d0
      b = 0.d0
      c = -dx(1)*dx(2)
      
      left  = i - 1
      right = i + 1
      up    = j + 1
      down  = j - 1

      
      lok  = left .GE. lo(1) .AND. sgn*phi(left,j) .GE. 0 .AND. type(left,j) .LE. 1
      rok  = right. LE. hi(1) .AND. sgn*phi(right,j) .GE. 0 . AND. type(right,j) .LE. 1 
      uok  = up .LE. hi(2) .AND. sgn*phi(i,up) .GE. 0  .AND. type(i,up) .LE. 1
      dok  = down . GE. lo(2) .AND. sgn*phi(i,down) .GE. 0  .AND. type(i,down) .LE. 1
!c      llok = FMMorder .EQ. 2 . AND. ll .GE. lo(1) .AND. sgn*phi(ll,j) .GE. 0 .AND. type(ll,j) .LE. 1
!c      rrok = FMMorder .EQ. 2 . AND. rr. LE. hi(1) .AND. sgn*phi(rr,j) .GE. 0 . AND. type(rr,j) .LE. 1 
!c      uuok = FMMorder .EQ. 2 . AND. uu .LE. hi(2) .AND. sgn*phi(i,uu) .GE. 0  .AND. type(i,uu) .LE. 1
!c      ddok = FMMorder .EQ. 2 . AND. dd . GE. lo(2) .AND. sgn*phi(i,dd) .GE. 0  .AND. type(i,dd) .LE. 1      
            
!c     FIXME: The following had a sign(right,...), mistake?
      if (lok .AND. rok) then
                    
        a = a + dx(2)/dx(1)
        b = b + min(sgn*phi(left,j),sgn*phi(right,j))*(dx(2)/dx(1))
        c = c + min(phi(left,j)**2,phi(right,j)**2)*(dx(2)/dx(1))
        
      else if (lok) then
      
        a = a + dx(2)/dx(1)
        b = b + sgn*phi(left,j)*(dx(2)/dx(1))
        c = c + phi(left,j)**2*(dx(2)/dx(1))
      
      else if (rok ) then

        a = a + dx(2)/dx(1)
        b = b + sgn*phi(right,j)*(dx(2)/dx(1))
        c = c + phi(right,j)**2*(dx(2)/dx(1))
        
      endif
      
      
      if (dok .AND. uok) then
      
        a = a + dx(1)/dx(2)
        b = b + min(sgn*phi(i,down),sgn*phi(i,up))*(dx(1)/dx(2))
        c = c + min(phi(i,down)**2,phi(i,up)**2)*(dx(1)/dx(2))
        
      else if (dok) then
      
        a = a + dx(1)/dx(2)
        b = b + sgn*phi(i,down)*(dx(1)/dx(2))
        c = c + phi(i,down)**2*(dx(1)/dx(2))
      
      else if (uok ) then

        a = a + dx(1)/dx(2)
        b = b + sgn*phi(i,up)*(dx(1)/dx(2))
        c = c + phi(i,up)**2*(dx(1)/dx(2))
        
      endif      
        
      b = -2*b
      
      if (b**2 - 4*a*c .LT. 0) then
      
         phi(i,j) = sgn*(-b)/(2*a)
         return
      
      endif
      
      phi(i,j) = sgn*(-b + SQRT(b**2-4*a*c))/(2*a)

    end subroutine EVAL
 
 
 
      subroutine ADDNODE(heap,i,j,n,phi,DIMS(phi), lo, hi, &
                        nbandsize, heaploc, DIMS(type))
      implicit none
      integer    DIMDEC(phi)
      REAL_T     phi(DIMV(phi))
      integer    nbandsize
      integer    lo(SDIM), hi(SDIM)
      integer    heap(nbandsize,SDIM)
      integer    i,j, n
      integer    DIMDEC(type)
      integer    heaploc(DIMV(type))

      integer index
      integer parent
      
      index = n + 1
      parent = index/2

      if (n .EQ. 0) then
         
         heap(index,1) = i
         heap(index,2) = j
         heaploc(i,j) = index        
         
         n=n+1
         
         return
      endif
      
      do while ( ABS(phi(heap(parent,1), heap(parent,2))) .GT. ABS(phi(i,j)))
            
        heap(index,1) = heap(parent,1)
        heap(index,2) = heap(parent,2)
        heaploc(heap(index,1),heap(index,2)) = index        
        index = parent
        parent = index/2
         
        if (parent .EQ. 0) exit
      enddo
      
      heap(index,1) = i
      heap(index,2) = j
      heaploc(i,j) = index
      n = n + 1
    end subroutine ADDNODE


      
      subroutine UPDATENODE(heap,i,j,n,phi,DIMS(phi), lo, hi,&
                           nbandsize, heaploc,DIMS(type))
      implicit none
      integer    DIMDEC(phi)
      REAL_T     phi(DIMV(phi))
      integer    nbandsize
      integer    lo(SDIM), hi(SDIM)
      integer    heap(nbandsize,SDIM)
      integer    i,j, n
      integer    DIMDEC(type)
      integer    heaploc(DIMV(type))

      integer index
      integer parent
      
      index = heaploc(i,j)
      parent = index/2

      if (index .EQ. 1) then
        return
      endif
      
      do while ( ABS(phi(heap(parent,1 ), heap(parent,2))) .GT. ABS(phi(i,j)))
            
        heap(index,1) = heap(parent,1)
        heap(index,2) = heap(parent,2)
        heaploc(heap(index,1),heap(index,2)) = index
        index = parent
        parent = index/2
         
        if (parent .EQ. 0) exit         
      enddo
      
      heap(index,1) = i
      heap(index,2) = j
      heaploc(i,j) = index
    end subroutine UPDATENODE
            
      
      
      
      subroutine RMVNODE(heap, i,j,n, phi, DIMS(phi),&
                        nbandsize,heaploc, DIMS(type))
      implicit none
      integer    DIMDEC(phi)
      REAL_T     phi(DIMV(phi))
      integer    nbandsize
      integer    heap(nbandsize,SDIM)
      integer    i,j, n
      integer    DIMDEC(type)
      integer    heaploc(DIMV(type))
      integer    index, left, right

      i=heap(1,1)
      j=heap(1,2)
      heaploc(i,j) = -1
      
      index = 1
      left  = 2*index
      right = 2*index + 1
   
      do while (.TRUE.)
         
         if (left. LE. n-1) then
            
            if( ABS(phi(heap(left,1),heap(left,2))) .LT. ABS(phi(heap(n,1),heap(n,2)))) then
               
               if (right .LE. n-1) then
                  if (ABS(phi(heap(left,1),heap(left,2))) .LT. ABS(phi(heap(right,1),heap(right,2))) ) then
                     heap(index,1) = heap(left,1)
                     heap(index,2) = heap(left,2)
                     heaploc(heap(index,1),heap(index,2)) = index
                     index = left
                     left=2*index
                     right=2*index+1
                  else
                     heap(index,1) = heap(right,1)
                     heap(index,2) = heap(right,2)
                     heaploc(heap(index,1),heap(index,2)) = index
                     index = right
                     left=2*index
                     right=2*index+1  
                  endif
               else                  
                  heap(index,1) = heap(left,1)
                  heap(index,2) = heap(left,2)
                  heaploc(heap(index,1),heap(index,2)) = index
                  
                  index = left
                  left=2*index
                  right=2*index+1  
               endif     
               
            else if (right .LE. n-1) then
         
               if (abs(phi(heap(right,1),heap(right,2))) .LT. abs( phi(heap(n,1),heap(n,2)))) then
                  heap(index,1) = heap(right,1)
                  heap(index,2) = heap(right,2)
                  heaploc(heap(index,1),heap(index,2)) = index
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
      if(n .GT. 1) then
        heaploc(heap(index,1),heap(index,2)) = index
      endif
      n = n - 1
    end subroutine RMVNODE
      
      
      
      INTEGER FUNCTION fastmarch2(phi,DIMS(phi),type,DIMS(type),&
                                  lo, hi, dx, nband, nbandsize, nbandnum,&
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
      
      integer i, j, n
      integer heap(nbandsize,SDIM)
      integer numtent
      
      numtent = 0
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            heaploc(i,j) = -1
         enddo
      enddo
      
      i = lo(1)-1 
      do j = lo(2), hi(2)
         if((type(i,j) .EQ. 0 .OR. type(i,j) .EQ. 1) .AND. sgn*phi(i,j) .GE. 0) then
            CALL UPDATE2(phi,i,j,sgn, type,heap,numtent,DIMS(phi),& 
                        DIMS(type), nbandsize, lo, hi, dx,heaploc)
         endif
      enddo
      
      i = hi(1)+1 
      do j = lo(2), hi(2) 
         if((type(i,j) .EQ. 0 .OR. type(i,j) .EQ. 1) .AND. sgn*phi(i,j) .GE. 0) then
            CALL UPDATE2(phi,i,j,sgn, type, heap,numtent,DIMS(phi),&
                        DIMS(type), nbandsize, lo, hi, dx,heaploc)
         endif
      enddo
      
      j = lo(2)-1 
      do i = lo(1) , hi(1)
         if((type(i,j) .EQ. 0 .OR. type(i,j) .EQ. 1) .AND. sgn*phi(i,j) .GE. 0) then
            CALL UPDATE2(phi,i,j,sgn, type, heap,numtent,DIMS(phi),&
                        DIMS(type), nbandsize, lo, hi, dx,heaploc)
         endif
      enddo
      
      j = hi(2)+1 
      do i = lo(1) , hi(1) 
         if((type(i,j) .EQ. 0 .OR. type(i,j) .EQ. 1) .AND. sgn*phi(i,j) .GE. 0) then
            CALL UPDATE2(phi,i,j,sgn, type, heap,numtent,DIMS(phi),&
      	                 DIMS(type), nbandsize, lo, hi, dx,heaploc)
         endif
      enddo
      
      FASTMARCH2 = 0
      do while (numtent .GT. 0 )
         CALL RMVNODE(heap, i,j,numtent, phi, DIMS(phi), &
                     nbandsize,heaploc,DIMS(type))
         if (abs(phi(i,j)) .LT. nbandwidth ) then
         
            FASTMARCH2 = 1

            if (type(i,j) .GE. 2) then
               nbandnum = nbandnum + 1         
               nband(nbandnum,1) = i
               nband(nbandnum,2) = j            
            endif
          
            type(i,j) = 0
            CALL UPDATE2(phi,i,j,sgn, type, heap,numtent,DIMS(phi),&
                        DIMS(type), nbandsize, lo, hi, dx,heaploc)

         else

            type(i,j) = 3
            phi(i,j) = sign(BOGUS,phi(i,j))
            exit

         endif
         
         if (abs(phi(i,j)) .GT. mineloc &
             .AND. abs(phi(i,j)) .LT. nbandwidth) then
          
            type(i,j) = 1
        endif
      enddo
      
      nband(nbandnum + 1,1) = -LARGEINT
      nband(nbandnum + 1,2) = -LARGEINT
      
      do while (numtent .GT. 0)
            
         CALL RMVNODE(heap, i,j,numtent, phi, DIMS(phi),&
                     nbandsize,heaploc,DIMS(type))
         type(i,j) = 3
         phi(i,j) = sign(BOGUS,phi(i,j))    
      enddo

    end FUNCTION fastmarch2



      subroutine UPDATE2(phi,i,j,sgn, type, heap,numtent,DIMS(phi),&
                        DIMS(type), nbandsize, lo, hi, dx, heaploc)      
      implicit none
      integer i,j
      integer DIMDEC(phi)
      integer DIMDEC(type)
      REAL_T  phi(DIMV(phi))
      integer type(DIMV(type))
      integer nbandsize
      integer heap(nbandsize,SDIM)
      integer lo(SDIM), hi(SDIM)
      integer numtent
      integer sgn
      integer n,ii,jj
      REAL_T dx(SDIM)
      integer heaploc(DIMV(type))

      logical EVAL2
      logical isnew
      integer min
      
      do n = 1,4
      
        ii = i -1 +n/2
        jj = j + 2*(n/3) - n/2
              
        if (ii .GE. lo(1) .AND. ii .LE. hi(1)&
            .AND. jj .GE. lo(2) .AND. jj .LE. hi(2)) then
        
           if (sgn*sign(1.d0,phi(ii,jj)) .GE. 0 .AND. abs(phi(ii,jj)) .GT. abs(phi(i,j)) &
               .AND. (abs(phi(ii,jj)) .GE. BOGUS .OR. .NOT.(sgn*phi(ii+1,jj) .LE. 0  &
               .OR. sgn*phi(ii-1,jj) .LE. 0 .OR. sgn*phi(ii,jj+1) .LE. 0 .OR. sgn*phi(ii,jj-1) .LE. 0) ))then 
!c     &          .OR. sgn*phi(ii+1,jj+1) .LE. 0 .OR. sgn*phi(ii-1,jj-1) .LE. 0 .OR. sgn*phi(ii+1,jj-1) .LE. 0 
!c     &          .OR. sgn*phi(ii-1,jj+1) .LE. 0) )) then
              
              isnew = EVAL2(phi,ii,jj,DIMS(phi), DIMS(type),&
                           lo,hi, type,sgn,dx, phi(i,j))
           
              isnew = isnew .OR. EVAL2(phi,ii,jj,DIMS(phi), DIMS(type),&
                                      lo,hi, type,sgn,dx, phi(ii,jj))
              
              if (isnew ) then
                 type(ii,jj) = min(2,type(ii,jj))
                 if(heaploc(ii,jj) .EQ. -1) then
                    CALL ADDNODE(heap,ii,jj, numtent,phi,DIMS(phi), &
                                lo, hi, nbandsize, heaploc,DIMS(type))
                 else
                    CALL UPDATENODE(heap,ii,jj, numtent,phi,DIMS(phi),& 
                                   lo, hi, nbandsize, heaploc,DIMS(type))   
                 endif
              endif
           endif    
        endif        
      enddo
    end subroutine UPDATE2
      
      
      
      LOGICAL FUNCTION EVAL2(phi,i,j,DIMS(phi), DIMS(type),lo,hi, type,sgn,dx,phisrc)      
      implicit none
      integer i,j     
      integer DIMDEC(phi)
      integer DIMDEC(type)
      integer lo(SDIM), hi(SDIM)
      REAL_T   phi(DIMV(phi))
      integer type(DIMV(type))
      REAL_T a,b,c
      REAL_T dx(SDIM)
      integer sgn
      integer  left,right,up,down
      REAL_T phisrc
      
      a = 0.0d0
      b = 0.0d0
      c = -dx(1)*dx(2)
      
      left  = i - 1
      right = i + 1
      up    = j + 1
      down  = j - 1
   
      if (sgn*phi(left,j) .GE. 0 .AND. type(left,j) .LE. 1 .AND. abs(phi(left,j)) .LE. abs(phisrc) &
          .AND. sgn*phi(right,j) .GE. 0 . AND. type(right,j) .LE. 1 .AND. abs(phi(right,j)) .LE. abs(phisrc) ) then
      
         a = a + dx(2)/dx(1)        
         b = b + min(sgn*phi(left,j),sgn*phi(right,j))*(dx(2)/dx(1))
         c = c + min(phi(left,j)**2,phi(right,j)**2)*(dx(2)/dx(1))
         
      else if (sgn*phi(left,j) .GE. 0  .AND. type(left,j) .LE. 1 .AND. abs(phi(left,j)) .LE. abs(phisrc) ) then
         
         a = a + dx(2)/dx(1)
         b = b + sgn*phi(left,j)*(dx(2)/dx(1))
         c = c + phi(left,j)**2*(dx(2)/dx(1))
         
      else if (sgn*phi(right,j) .GE. 0  .AND. type(right,j) .LE. 1 .AND. abs(phi(right,j)) .LE. abs(phisrc) ) then
         
         a = a + dx(2)/dx(1)
         b = b + sgn*phi(right,j)*(dx(2)/dx(1))
         c = c + phi(right,j)**2*(dx(2)/dx(1))
         
      endif

      if ( sgn*phi(i,down) .GE. 0  .AND. type(i,down) .LE. 1 .AND. abs(phi(i,down)) .LE. abs(phisrc) &
          .AND. sgn*phi(i,up) .GE. 0  .AND. type(i,up) .LE. 1 .AND. abs(phi(i,up)) .LE. abs(phisrc) ) then
         
         a = a + dx(1)/dx(2)
         b = b + min(sgn*phi(i,down),sgn*phi(i,up))*(dx(1)/dx(2))
         c = c + min(phi(i,down)**2,phi(i,up)**2)*(dx(1)/dx(2))
         
      else if (sgn*phi(i,down) .GE. 0  .AND. type(i,down) .LE. 1 .AND. abs(phi(i,down)) .LE. abs(phisrc) ) then
         
         a = a + dx(1)/dx(2)
         b = b + sgn*phi(i,down)*(dx(1)/dx(2))
         c = c + phi(i,down)**2*(dx(1)/dx(2))
         
      else if (sgn*phi(i,up) .GE. 0  .AND. type(i,up) .LE. 1 .AND. abs(phi(i,up)) .LE. abs(phisrc) ) then
         
         a = a + dx(1)/dx(2)
         b = b + sgn*phi(i,up)*(dx(1)/dx(2))
         c = c + phi(i,up)**2*(dx(1)/dx(2))
         
      endif      
      
      b = -2*b
      if (a .EQ. 0.d0) then

        EVAL2 = .FALSE.
        return
       
      endif
      
      EVAL2 = .FALSE.
      
      if (b**2 - 4*a*c .LT. 0.d0) then
         if (ABS(phi(i,j)) .GT. (-b)/(2*a) + 1.d-10) then
            phi(i,j) = sgn*(-b)/(2*a)
            EVAL2 = .TRUE.
            return
         endif
      endif
      
      if (abs(phi(i,j)) .GT. (-b+SQRT(b**2-4*a*c))/(2*a)+1.d-10) then
        phi(i,j) = sgn*(-b+SQRT(b**2-4*a*c))/(2*a)
        EVAL2 = .TRUE.
      endif
      
    end FUNCTION EVAL2


     
      subroutine mine(type, DIMS(type), nband, nbandsize,&
                      mine, minesize,lo, hi) bind(C,name="mine")
      implicit none     
      integer DIMDEC(type)
      integer nbandsize, minesize
      integer lo(SDIM), hi(SDIM)
      integer  type(DIMV(type))
      integer  nband(nbandsize,SDIM)
      integer  mine(minesize,SDIM)
      
      integer i, j, p
      integer nummine
      
      nummine = 0
      p = 1
      do while(nband(p,1) .GT. -LARGEINT)
      
        i = nband(p,1)
        j = nband(p,2)        
        p = p + 1
        
        if (type(i,j) .EQ. 1) then
          nummine = nummine +1
          mine(nummine,1) = i
          mine(nummine,2) = j
        endif
        
      enddo
      
      mine(nummine+1,1) = -LARGEINT
      mine(nummine+1,2) = -LARGEINT      
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
  
end module LS_2d_module
