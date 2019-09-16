
#define SDIM 2

module convection_eb_mod

   use amrex_fort_module,       only: ar => amrex_real
   use iso_c_binding ,          only: c_int
   use param,                   only: zero, half, one, my_huge, &
                                      minf_, nsw_, fsw_, psw_, pinf_, pout_
   use amrex_error_module,      only: amrex_abort
   use amrex_ebcellflag_module, only: is_covered_cell, is_single_valued_cell, &
        &                             get_neighbor_cells

   implicit none
   private

contains

   !
   ! EB x-direction MAC velocity
   !
   subroutine compute_velocity_at_x_faces_eb ( lo, hi, u, ulo, uhi,  &
        vel, vello, velhi, &
        slopes, slo, shi, areafrac, alo, ahi,     &
        cent, clo, chi, flags, flo, fhi, &
        ubc, &
        domlo, domhi ) bind(C)

      use convection_mod 

      ! Tile bounds ( x-face centered )
      integer(c_int),  intent(in   ) :: lo(SDIM),  hi(SDIM)

      ! Array Bounds
      integer(c_int),  intent(in   ) :: slo(SDIM), shi(SDIM)
      integer(c_int),  intent(in   ) :: ulo(SDIM), uhi(SDIM)
      integer(c_int),  intent(in   ) :: vello(SDIM), velhi(SDIM)
      integer(c_int),  intent(in   ) :: alo(SDIM), ahi(SDIM)
      integer(c_int),  intent(in   ) :: clo(SDIM), chi(SDIM)
      integer(c_int),  intent(in   ) :: flo(SDIM), fhi(SDIM)

      ! Domain bounds
      integer(c_int),  intent(in   ) :: domlo(SDIM), domhi(SDIM)

      ! Arrays
      real(ar),        intent(in   ) ::                                     &
           & vel(vello(1):velhi(1),vello(2):velhi(2),SDIM),  &
           & slopes(slo(1):shi(1),slo(2):shi(2),SDIM),           &
           & areafrac(alo(1):ahi(1),alo(2):ahi(2)),           &
           & cent(clo(1):chi(1),clo(2):chi(2),2)
      integer(c_int),  intent(in   ) :: &
           & flags(flo(1):fhi(1),flo(2):fhi(2))

      ! Staggered velocity
      real(ar),        intent(inout) ::                    &
           & u(ulo(1):uhi(1),ulo(2):uhi(2))

      ! BC types
      integer(c_int), intent(in   ) ::  ubc(SDIM,2)
      ! integer(c_int), intent(in   ) ::  &
      !      & bc_ilo(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
      !      & bc_ihi(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2)

      ! Local variables
      integer(c_int)                 :: i, j
      integer, parameter             :: bc_list(6) = [MINF_, NSW_, FSW_, PSW_, PINF_, POUT_]
      real(ar)                       :: upls, umns
      
      ! First we compute the face centered MAC velocity
      ! We need 1 layer of ghosts cells in y and z for interpolation (see next step)
      do j = lo(2)-1, hi(2)+1
         do i = lo(1), hi(1)
            if ( areafrac(i,j) > zero ) then
               !
               ! FIXME - only periodic right now. 
               !         for other phys bcs will have to update below to use
               !         bc(SDIM, 2)      
               ! if ( ( i == domlo(1) ) .and. any( bc_ilo(j,k,1) == bc_list) ) then
               !    u(i,j,k) = vel(i-1,j,k,1)
               ! else if ( ( i == domhi(1)+1 ) .and. any( bc_ihi(j,k,1) == bc_list) ) then
               !    u(i,j,k) = vel(i,j,k,1)
               ! else
               upls     = vel(i  ,j,1) - half * slopes(i  ,j,1)
               umns     = vel(i-1,j,1) + half * slopes(i-1,j,1)
               u(i,j) = upwind_normal( umns, upls )
               ! end if
            else
               u(i,j) = my_huge
            end if
         end do
      end do

   end subroutine compute_velocity_at_x_faces_eb

   !
   ! EB y-direction MAC velocity
   !
   subroutine compute_velocity_at_y_faces_eb ( lo, hi, v, vlo, vhi,  &
        vel, vello, velhi, slopes, slo, shi, areafrac, alo, ahi,     &
        cent, clo, chi, flags, flo, fhi, vbc, &
        domlo, domhi ) bind(C)

      use convection_mod 

      ! Tile bounds ( x-face centered )
      integer(c_int),  intent(in   ) :: lo(SDIM),  hi(SDIM)

      ! Array Bounds
      integer(c_int),  intent(in   ) :: slo(SDIM), shi(SDIM)
      integer(c_int),  intent(in   ) :: vlo(SDIM), vhi(SDIM)
      integer(c_int),  intent(in   ) :: vello(SDIM), velhi(SDIM)
      integer(c_int),  intent(in   ) :: alo(SDIM), ahi(SDIM)
      integer(c_int),  intent(in   ) :: clo(SDIM), chi(SDIM)
      integer(c_int),  intent(in   ) :: flo(SDIM), fhi(SDIM)

      ! Domain bounds
      integer(c_int),  intent(in   ) :: domlo(SDIM), domhi(SDIM)

      ! Arrays
      real(ar),        intent(in   ) ::                                     &
           & vel(vello(1):velhi(1),vello(2):velhi(2),SDIM),  &
           & slopes(slo(1):shi(1),slo(2):shi(2),SDIM),           &
           & areafrac(alo(1):ahi(1),alo(2):ahi(2)),           &
           & cent(clo(1):chi(1),clo(2):chi(2),2)
      integer(c_int),  intent(in   ) :: &
           & flags(flo(1):fhi(1),flo(2):fhi(2))

      ! Staggered velocity
      real(ar),        intent(inout) ::                    &
           & v(vlo(1):vhi(1),vlo(2):vhi(2))

      ! BC types
       integer(c_int), intent(in   ) ::  vbc(SDIM,2)
      ! integer(c_int), intent(in   ) ::  &
      !      & bc_jlo(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
      !      & bc_jhi(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2)

      ! Local variables
      integer(c_int)                 :: i, j
      integer, parameter             :: bc_list(6) = [MINF_, NSW_, FSW_, PSW_, PINF_, POUT_]
      real(ar)                       :: vpls, vmns

      do j = lo(2), hi(2)
         do i = lo(1)-1, hi(1)+1
            if ( areafrac(i,j) > zero ) then
               !
               ! FIXME - only periodic right now. 
               !         for other phys bcs will have to update below to use
               !         bc(SDIM, 2)
               
               ! if ( ( j == domlo(2) ) .and. any(bc_jlo(i,k,1) == bc_list) ) then
               !    v(i,j,k) = vel(i,j-1,k,2)
               ! else if ( ( j == domhi(2)+1 ) .and. any(bc_jhi(i,k,1) == bc_list) ) then
               !    v(i,j,k) = vel(i,j,k,2)
               ! else
               vpls     = vel(i,j  ,2) - half * slopes(i,j  ,2)
               vmns     = vel(i,j-1,2) + half * slopes(i,j-1,2)
               v(i,j) = upwind_normal( vmns, vpls )
               ! end if
            else
               v(i,j) = my_huge
            end if
         end do
      end do

   end subroutine compute_velocity_at_y_faces_eb
   
   ! $************************************************
   ! $ WARNING:
   ! $
   ! $ For now the convective term is only div(u_mac . u)
   ! $ The term -u*grad(u_mac) is omitted for velocity
   ! $
   ! $************************************************
   subroutine compute_aofs_eb ( lo, hi, &
                                 aofs, glo, ghi, &
                                 vel, vello, velhi, &
                                 u, ulo, uhi, &
                                 xflx, xflxlo, xflxhi, &
                                 xstate, xstatelo, xstatehi, &
                                 afrac_x, axlo, axhi, &
                                 cent_x,  cxlo, cxhi, &
                                 xslopes, sxlo, sxhi, &                        
                                 v, vlo, vhi, &
                                 yflx, yflxlo, yflxhi, &
                                 ystate, ystatelo, ystatehi, &
                                 afrac_y, aylo, ayhi, &
                                 cent_y,  cylo, cyhi, &
                                 yslopes, sylo, syhi, &
                                 flags,    flo,  fhi, &
                                 vfrac,   vflo, vfhi, &
                                 bcent,    blo,  bhi, &
                                 domlo, domhi, &
                                 dx, nc, ng, known_edgestate ) bind(C)

      use divop_mod, only: compute_divop

      ! Tile bounds
      integer(c_int),  intent(in   ) :: lo(SDIM),  hi(SDIM)

      ! Array Bounds
      integer(c_int),  intent(in   ) :: sxlo(SDIM), sxhi(SDIM)
      integer(c_int),  intent(in   ) :: sylo(SDIM), syhi(SDIM)
      integer(c_int),  intent(in   ) :: glo(SDIM), ghi(SDIM)
      integer(c_int),  intent(in   ) :: vello(SDIM), velhi(SDIM)
      integer(c_int),  intent(in   ) :: ulo(SDIM), uhi(SDIM)
      integer(c_int),  intent(in   ) :: vlo(SDIM), vhi(SDIM)
      integer(c_int),  intent(in   ) :: xflxlo(SDIM), xflxhi(SDIM)
      integer(c_int),  intent(in   ) :: yflxlo(SDIM), yflxhi(SDIM)
      integer(c_int),  intent(in   ) :: xstatelo(SDIM), xstatehi(SDIM)
      integer(c_int),  intent(in   ) :: ystatelo(SDIM), ystatehi(SDIM)
      integer(c_int),  intent(in   ) :: axlo(SDIM), axhi(SDIM)
      integer(c_int),  intent(in   ) :: aylo(SDIM), ayhi(SDIM)
      integer(c_int),  intent(in   ) :: cxlo(SDIM), cxhi(SDIM)
      integer(c_int),  intent(in   ) :: cylo(SDIM), cyhi(SDIM)
      integer(c_int),  intent(in   ) ::  flo(SDIM),  fhi(SDIM)
      integer(c_int),  intent(in   ) :: vflo(SDIM), vfhi(SDIM)
      integer(c_int),  intent(in   ) ::  blo(SDIM),  bhi(SDIM)
      integer(c_int),  intent(in   ) :: domlo(SDIM), domhi(SDIM)
      integer(c_int),  intent(in   ) :: nc, ng, known_edgestate

      ! Grid
      real(ar),        intent(in   ) :: dx(SDIM)

      ! Arrays
      real(ar),        intent(in   ) ::                            &
           & vel(vello(1):velhi(1),vello(2):velhi(2),nc)    , &
           & xslopes(sxlo(1):sxhi(1),sxlo(2):sxhi(2),nc), &
           & yslopes(sylo(1):syhi(1),sylo(2):syhi(2),nc), &
           & u(ulo(1):uhi(1),ulo(2):uhi(2)), &
           & v(vlo(1):vhi(1),vlo(2):vhi(2)), &
           & afrac_x(axlo(1):axhi(1),axlo(2):axhi(2)),       &
           & afrac_y(aylo(1):ayhi(1),aylo(2):ayhi(2)),       &
           & cent_x(cxlo(1):cxhi(1),cxlo(2):cxhi(2),2),&
           & cent_y(cylo(1):cyhi(1),cylo(2):cyhi(2),2),&
           & vfrac(vflo(1):vfhi(1),vflo(2):vfhi(2)), &
           & bcent(blo(1):bhi(1),blo(2):bhi(2),SDIM)

      real(ar),        intent(  out) ::                           &
           & aofs(glo(1):ghi(1),glo(2):ghi(2),nc), &
           & xflx(xflxlo(1):xflxhi(1),xflxlo(2):xflxhi(2),nc), &
           & yflx(yflxlo(1):yflxhi(1),yflxlo(2):yflxhi(2),nc)

      real(ar),        intent(inout) ::                           &
           & xstate(xstatelo(1):xstatehi(1),xstatelo(2):xstatehi(2),nc), &
           & ystate(ystatelo(1):ystatehi(1),ystatelo(2):ystatehi(2),nc)
      
      integer(c_int), intent(in   ) ::  &
           & flags(flo(1):fhi(1),flo(2):fhi(2))
      
      ! BC types
      ! integer(c_int), intent(in   ) ::  &
      integer(c_int)  ::  &
           & bc_ilo(domlo(2)-ng:domhi(2)+ng,2), &
           & bc_ihi(domlo(2)-ng:domhi(2)+ng,2), &
           & bc_jlo(domlo(1)-ng:domhi(1)+ng,2), &
           & bc_jhi(domlo(1)-ng:domhi(1)+ng,2)
           

      ! Temporary array to handle convective fluxes at the cell faces (staggered)
      ! Just reserve space for the tile + 3 ghost layers
      integer, parameter :: nh = 3 ! Number of Halo layers
      real(ar) :: fx(lo(1)-nh:hi(1)+nh+1,lo(2)-nh:hi(2)+nh    ,nc)
      real(ar) :: fy(lo(1)-nh:hi(1)+nh  ,lo(2)-nh:hi(2)+nh+1  ,nc)

      ! Check number of ghost cells
      if (ng < 5) call amrex_abort( "compute_divop(): ng must be >= 5")

      !
      ! First compute the convective fluxes at the face center
      ! Do this on ALL faces on the tile, i.e. INCLUDE as many ghost faces as
      ! possible
      !
      block
         real(ar)               :: u_face, v_face
         real(ar)               :: upls, umns, vpls, vmns
         integer                :: i, j, n
         integer, parameter     :: bc_list(6) = [MINF_, NSW_, FSW_, PSW_, PINF_, POUT_]

      ! FIXME only period for now, so hack bc arrays
      bc_ilo = -1 
      bc_ihi = -1  
      bc_jlo = -1 
      bc_jhi = -1 

      
!      write(*,*) 'DEBUG IN  conv_eb ',lbound(xslopes),ubound(xslopes)
!      write(*,*) 'DEBUG conv eb',lo,hi,nh,nc
      do n = 1, nc

         !
         ! ===================   X   ===================
         !

         do j = lo(2)-nh, hi(2)+nh
            do i = lo(1)-nh, hi(1)+nh+1
            
  !  write(*,*) 'DEBUG conv eb ',i,j,n,xslopes(i  ,j,n)

            
            !write(*,*) 'DEBUG conv eb ',i,j,n,vel(i  ,j,n)
            if (known_edgestate == 0) then
            
               if ( afrac_x(i,j) > zero ) then
                  if ( i <= domlo(1) .and. any(bc_ilo(j,1) == bc_list) ) then
                     u_face =  vel(domlo(1)-1,j,n)
                  else if ( i >= domhi(1)+1 .and. any(bc_ihi(j,1) == bc_list ) ) then
                     u_face =  vel(domhi(1)+1,j,n)
                  else
                     upls  = vel(i  ,j,n) - half * xslopes(i  ,j,n)
                     umns  = vel(i-1,j,n) + half * xslopes(i-1,j,n)
                     
                     u_face = upwind( umns, upls, u(i,j) )
                  end if
               else
                  u_face = my_huge
                  !write(*,*) 'DEBUG conv eb ',i,j,n,u_face,afrac_x(i,j)
               end if
               
               xstate(i,j,n)   = u_face
            
            end if
            !write(*,*) 'DEBUG conv eb ',i,j,n,u(i,j) , xstate(i,j,n)
               fx(i,j,n) = u(i,j) * xstate(i,j,n)
               xflx(i,j,n)     = fx(i,j,n) * dx(1)


            end do
         end do

         !
         ! ===================   Y   ===================
         !
         do j = lo(2)-nh, hi(2)+nh+1
            do i = lo(1)-nh, hi(1)+nh
            
            if (known_edgestate == 0) then
            
               if ( afrac_y(i,j) > zero ) then
                  if ( j <= domlo(2) .and. any(bc_jlo(i,1) == bc_list) ) then
                     v_face =  vel(i,domlo(2)-1,n)
                  else if ( j >= domhi(2)+1 .and. any(bc_jhi(i,1) == bc_list ) ) then
                     v_face =  vel(i,domhi(2)+1,n)
                  else
                     vpls  = vel(i,j  ,n) - half * yslopes(i,j  ,n)
                     vmns  = vel(i,j-1,n) + half * yslopes(i,j-1,n)

                     v_face = upwind( vmns, vpls, v(i,j) )
                  end if
               else
                  v_face = my_huge
               end if
            
               ystate(i,j,n)   = v_face
            
             end if   
            
               fy(i,j,n) = v(i,j) * ystate(i,j,n)
               yflx(i,j,n)     = fy(i,j,n)  * dx(2)

            end do
         end do
         
      end do ! end do n = 1, nc

      end block

      divop: block
         ! Compute div(tau) with EB algorithm
         integer(c_int)  :: fxlo(SDIM), fxhi(SDIM), fylo(SDIM), fyhi(SDIM), fzlo(SDIM), fzhi(SDIM)

         fxlo = lo - nh
         fylo = lo - nh

         fxhi = hi + nh + [1,0]
         fyhi = hi + nh + [0,1]


         call compute_divop(lo, hi, &
                            aofs, glo, ghi, &
                            vel, vlo, vhi, &
                            fx, fxlo, fxhi, &
                            fy, fylo, fyhi, &
                            afrac_x, axlo, axhi, &
                            afrac_y, aylo, ayhi, &
                            cent_x, cxlo, cxhi, &
                            cent_y, cylo, cyhi, &
                            flags, flo, fhi, &
                            vfrac, vflo, vfhi, &
                            bcent, blo, bhi, &
                            domlo, domhi, &
                            dx, ng, nc )
      end block divop


   end subroutine compute_aofs_eb

   
   ! Upwind non-normal velocity
   function upwind ( velmns, velpls, uedge ) result (ev)
     
      ! Small value to protect against tiny velocities used in upwinding
      real(ar),        parameter     :: small_vel = 1.0d-10

      real(ar), intent(in) :: velmns, velpls, uedge
      real(ar)             :: ev

      if ( abs(uedge) .lt. small_vel) then
         ev = half * ( velpls + velmns )
      else
         ev = merge ( velmns, velpls, uedge >= zero )
      end if

   end function upwind

end module convection_eb_mod
