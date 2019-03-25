!
!
!  This module contains the subroutines to compute the three components
!  of the convection term (U.grad)U
!
!

#define SDIM 2

module convection_mod

   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one, &
                               minf_, nsw_, fsw_, psw_, pinf_, pout_

   implicit none
   public upwind, upwind_normal
   private

contains

  subroutine compute_velocity_at_faces ( lo, hi, vel, vello, velhi, &
       u, ulo, uhi, v, vlo, vhi, &
       xslopes, yslopes, slo, shi, ubc, vbc, &
      ! bc_ilo, bc_ihi, bc_jlo, bc_jhi, bc_klo, bc_khi,ng, &
       domlo, domhi ) bind(C)

      ! Tile bounds
      integer(c_int),  intent(in   ) :: lo(SDIM),  hi(SDIM)

      ! Array Bounds
      integer(c_int),  intent(in   ) :: slo(SDIM), shi(SDIM)
      integer(c_int),  intent(in   ) :: ulo(SDIM), uhi(SDIM)
      integer(c_int),  intent(in   ) :: vlo(SDIM), vhi(SDIM)
      integer(c_int),  intent(in   ) :: vello(SDIM), velhi(SDIM)

      ! Domain bounds
      integer(c_int),  intent(in   ) :: domlo(SDIM), domhi(SDIM)

      ! Nghost
!      integer(c_int),  intent(in   ) :: ng

      ! Velocity Array
      real(ar),        intent(in   ) ::                            &
           & vel(vello(1):velhi(1),vello(2):velhi(2),SDIM)    , &
           & xslopes(slo(1):shi(1),slo(2):shi(2),SDIM), &
           & yslopes(slo(1):shi(1),slo(2):shi(2),SDIM)

      ! Staggered velocity
      real(ar),        intent(inout) ::                      &
           & u(ulo(1):uhi(1),ulo(2):uhi(2)), &
           & v(vlo(1):vhi(1),vlo(2):vhi(2))

      ! BC types
       integer(c_int), intent(in   ) ::  ubc(SDIM,2),vbc(SDIM,2)
      ! integer(c_int), intent(in   ) ::  &
      !      & bc_ilo(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
      !      & bc_ihi(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
      !      & bc_jlo(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
      !      & bc_jhi(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
      !      & bc_klo(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
      !      & bc_khi(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

      ! Local variables
      integer(c_int)                 :: i, j
      integer, parameter             :: bc_list(6) = [MINF_, NSW_, FSW_, PSW_, PINF_, POUT_]
      real(ar)                       :: upls, umns, vpls, vmns

      !
      ! FIXME - the bcs here are not tested in any way
      !

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)+1
            ! if ( ( i == domlo(1) ) .and. any(bc_ilo(j,k,1) == bc_list) ) then
            if ( ( i == domlo(1) ) .and. any(ubc(1,1) == bc_list) ) then
               u(i,j) = vel(i-1,j,1)
               ! else if ( ( i == domhi(1)+1 ) .and. any(bc_ihi(j,k,1) == bc_list) ) then
            else if ( ( i == domhi(1)+1 ) .and. any(ubc(1,2) == bc_list) ) then
               u(i,j) = vel(i,j,1)
            else
               upls     = vel(i  ,j,1) - half * xslopes(i  ,j,1)
               umns     = vel(i-1,j,1) + half * xslopes(i-1,j,1)
               u(i,j) = upwind_normal( umns, upls )
            end if
         end do
      end do

      do j = lo(2), hi(2)+1
         do i = lo(1), hi(1)
            ! if ( ( j == domlo(2) ) .and. any(bc_jlo(i,k,1) == bc_list) ) then
            if ( ( j == domlo(2) ) .and. any(vbc(2,1) == bc_list) ) then
               v(i,j) = vel(i,j-1,2)
               ! else if ( ( j == domhi(2)+1 ) .and. any(bc_jhi(i,k,1) == bc_list) ) then
            else if ( ( j == domhi(2)+1 ) .and. any(vbc(2,2) == bc_list) ) then
               v(i,j) = vel(i,j,2)
            else
               vpls     = vel(i,j  ,2) - half * yslopes(i,j  ,2)
               vmns     = vel(i,j-1,2) + half * yslopes(i,j-1,2)
               v(i,j) = upwind_normal( vmns, vpls )
            end if
         end do
      end do

   end subroutine compute_velocity_at_faces

   !#####################################################

   !  MAC VERSION

   !#####################################################
   subroutine compute_ugradu(lo, hi, &
                             ugradu, glo, ghi, &
                             vel, vello, velhi, &
                             u, ulo, uhi, &
                             v, vlo, vhi, &
                             xslopes, yslopes, slo, shi, &
                             domlo, domhi, &
                             ! bc_ilo_type, bc_ihi_type, &
                             ! bc_jlo_type, bc_jhi_type, &
                             dx, ng) bind(C)

      ! Tile bounds
      integer(c_int),  intent(in   ) :: lo(SDIM),  hi(SDIM)

      ! Array Bounds
      integer(c_int),  intent(in   ) :: slo(SDIM), shi(SDIM)
      integer(c_int),  intent(in   ) :: glo(SDIM), ghi(SDIM)
      integer(c_int),  intent(in   ) :: vello(SDIM), velhi(SDIM)
      integer(c_int),  intent(in   ) :: ulo(SDIM), uhi(SDIM)
      integer(c_int),  intent(in   ) :: vlo(SDIM), vhi(SDIM)
      integer(c_int),  intent(in   ) :: domlo(SDIM), domhi(SDIM), ng

      ! Grid
      real(ar),        intent(in   ) :: dx(SDIM)

      ! Velocity Array
      real(ar),        intent(in   ) ::                            &
           & vel(vello(1):velhi(1),vello(2):velhi(2),SDIM)    , &
           & xslopes(slo(1):shi(1),slo(2):shi(2),SDIM), &
           & yslopes(slo(1):shi(1),slo(2):shi(2),SDIM), &
           & u(ulo(1):uhi(1),ulo(2):uhi(2)), &
           & v(vlo(1):vhi(1),vlo(2):vhi(2))

      real(ar),        intent(  out) ::                           &
           & ugradu(glo(1):ghi(1),glo(2):ghi(2),SDIM)

      ! BC types
      !integer(c_int), intent(in   ) ::  &
      integer(c_int) ::  &
           & bc_ilo_type(domlo(2)-ng:domhi(2)+ng,2), &
           & bc_ihi_type(domlo(2)-ng:domhi(2)+ng,2), &
           & bc_jlo_type(domlo(1)-ng:domhi(1)+ng,2), &
           & bc_jhi_type(domlo(1)-ng:domhi(1)+ng,2)
      
      ! Local variables
      integer(c_int)                 :: i, j
      real(ar)                       :: idx, idy
      real(ar)                       :: upls, umns, vpls, vmns
      real(ar)                       :: u_e, u_w, u_s, u_n
      real(ar)                       :: v_e, v_w, v_s, v_n 
      real(ar)                       :: divumac
      integer, parameter             :: bc_list(6) = [MINF_, NSW_, FSW_, PSW_, PINF_, POUT_]


      ! FIXME only period for now, so hack bc arrays
      bc_ilo_type = -1 
      bc_ihi_type = -1  
      bc_jlo_type = -1 
      bc_jhi_type = -1 
      
      idx = one / dx(1)
      idy = one / dx(2)

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)

            ! ****************************************************
            ! West face
            ! ****************************************************
            
            ! In the case of MINF, NSW, FSW, PSW we are using the prescribed Dirichlet value
            ! In the case of PINF, POUT          we are using the upwind value
            if (i.eq.domlo(1) .and. any(bc_ilo_type(j,1) == bc_list ) ) then
               u_w =  vel(i-1,j,1)
               v_w =  vel(i-1,j,2)
            else
               upls  = vel(i  ,j,1) - half * xslopes(i  ,j,1)
               umns  = vel(i-1,j,1) + half * xslopes(i-1,j,1)
               vpls  = vel(i  ,j,2) - half * xslopes(i  ,j,2)
               vmns  = vel(i-1,j,2) + half * xslopes(i-1,j,2)

               u_w   = upwind( umns, upls, u(i,j) )
               v_w   = upwind( vmns, vpls, u(i,j) )
            endif

            ! ****************************************************
            ! East face
            ! ****************************************************
            
            ! In the case of MINF, NSW, FSW, PSW we are using the prescribed Dirichlet value
            ! In the case of PINF, POUT          we are using the upwind value
            if (i.eq.domhi(1) .and. any(bc_ihi_type(j,1) == bc_list ) ) then
               u_e =  vel(i+1,j,1)
               v_e =  vel(i+1,j,2)
            else
               upls  = vel(i+1,j,1) - half * xslopes(i+1,j,1)
               umns  = vel(i  ,j,1) + half * xslopes(i  ,j,1)
               vpls  = vel(i+1,j,2) - half * xslopes(i+1,j,2)
               vmns  = vel(i  ,j,2) + half * xslopes(i  ,j,2)
               
               u_e   = upwind( umns, upls, u(i+1,j) )
               v_e   = upwind( vmns, vpls, u(i+1,j) )
            endif

            ! ****************************************************
            ! South face
            ! ****************************************************
            
            ! In the case of MINF, NSW, FSW, PSW we are using the prescribed Dirichlet value
            ! In the case of PINF, POUT          we are using the upwind value
            if (j.eq.domlo(2) .and. any(bc_jlo_type(i,1) == bc_list ) ) then
               u_s =  vel(i,j-1,1)
               v_s =  vel(i,j-1,2)
            else
               upls  = vel(i,j  ,1) - half * yslopes(i,j  ,1)
               umns  = vel(i,j-1,1) + half * yslopes(i,j-1,1)
               vpls  = vel(i,j  ,2) - half * yslopes(i,j  ,2)
               vmns  = vel(i,j-1,2) + half * yslopes(i,j-1,2)

               v_s   = upwind( vmns, vpls, v(i,j) )
               u_s   = upwind( umns, upls, v(i,j) )
            endif

            ! ****************************************************
            ! North face
            ! ****************************************************
            
            ! In the case of MINF, NSW, FSW, PSW we are using the prescribed Dirichlet value
            ! In the case of PINF, POUT          we are using the upwind value
            if (j.eq.domhi(2) .and.  any(bc_jhi_type(i,1) == bc_list ) ) then
               u_n =  vel(i,j+1,1)
               v_n =  vel(i,j+1,2)
            else
               upls  = vel(i,j+1,1) - half * yslopes(i,j+1,1)
               umns  = vel(i,j  ,1) + half * yslopes(i,j  ,1)
               vpls  = vel(i,j+1,2) - half * yslopes(i,j+1,2)
               vmns  = vel(i,j  ,2) + half * yslopes(i,j  ,2)

               v_n   = upwind( vmns, vpls, v(i,j+1) )
               u_n   = upwind( umns, upls, v(i,j+1) )
            endif


            ! ****************************************************
            ! Define convective terms -- conservatively
            !   ugradu = ( div(u^MAC u^cc) - u^cc div(u^MAC) )
            ! ****************************************************
            
            divumac = (u(i+1,j) - u(i,j)) * idx + &
                      (v(i,j+1) - v(i,j)) * idy 

            ugradu(i,j,1) = (u(i+1,j) * u_e - u(i,j) * u_w) * idx + &
                            (v(i,j+1) * u_n - v(i,j) * u_s) * idy - &
                            vel(i,j,1) * divumac
            ugradu(i,j,2) = (u(i+1,j) * v_e - u(i,j) * v_w) * idx + &
                            (v(i,j+1) * v_n - v(i,j) * v_s) * idy - &
                            vel(i,j,2) * divumac

               ! ****************************************************
               ! Return the negative
               ! ****************************************************

               !ugradu(i,j,k,1) = -ugradu(i,j,k,1)
               !ugradu(i,j,k,2) = -ugradu(i,j,k,2)
               !ugradu(i,j,k,3) = -ugradu(i,j,k,3)

         end do
      end do

   end subroutine compute_ugradu

   !
   ! Compute the conservative advective term for scalar
   ! assumes ghost cells filled
   !
   !FIXME? -- single component
   subroutine compute_divuc(lo, hi, &
                            divuc, glo, ghi, &
                            vel, vello, velhi, &
                            u, ulo, uhi, &
                            v, vlo, vhi, &
                            xslopes, yslopes, slo, shi, &
                            domlo, domhi, &
                            ! bc_ilo_type, bc_ihi_type, &
                            ! bc_jlo_type, bc_jhi_type, &
                            dx, ng) bind(C)

      ! Tile bounds
      integer(c_int),  intent(in   ) :: lo(SDIM),  hi(SDIM)

      ! Array Bounds
      integer(c_int),  intent(in   ) :: slo(SDIM), shi(SDIM)
      integer(c_int),  intent(in   ) :: glo(SDIM), ghi(SDIM)
      integer(c_int),  intent(in   ) :: vello(SDIM), velhi(SDIM)
      integer(c_int),  intent(in   ) :: ulo(SDIM), uhi(SDIM)
      integer(c_int),  intent(in   ) :: vlo(SDIM), vhi(SDIM)
      integer(c_int),  intent(in   ) :: domlo(SDIM), domhi(SDIM), ng

      ! Grid
      real(ar),        intent(in   ) :: dx(SDIM)

      integer(c_int),     parameter  :: nc = 1
      ! Velocity Array
      real(ar),        intent(in   ) ::                            &
           & vel(vello(1):velhi(1),vello(2):velhi(2),nc), &
           & xslopes(slo(1):shi(1),slo(2):shi(2),nc), &
           & yslopes(slo(1):shi(1),slo(2):shi(2),nc), &
           & u(ulo(1):uhi(1),ulo(2):uhi(2)), &
           & v(vlo(1):vhi(1),vlo(2):vhi(2))

      real(ar),        intent(  out) ::                           &
           & divuc(glo(1):ghi(1),glo(2):ghi(2),nc)

      ! BC types
      !integer(c_int), intent(in   ) ::  &
      integer(c_int) ::  &
           & bc_ilo_type(domlo(2)-ng:domhi(2)+ng,2), &
           & bc_ihi_type(domlo(2)-ng:domhi(2)+ng,2), &
           & bc_jlo_type(domlo(1)-ng:domhi(1)+ng,2), &
           & bc_jhi_type(domlo(1)-ng:domhi(1)+ng,2)

      ! Local variables
      integer(c_int)                 :: i, j
      real(ar)                       :: idx, idy
      real(ar)                       :: upls, umns, vpls, vmns
      real(ar)                       :: u_e, u_w, u_s, u_n
      real(ar)                       :: v_e, v_w, v_s, v_n
      real(ar)                       :: divumac
      integer, parameter             :: bc_list(6) = [MINF_, NSW_, FSW_, PSW_, PINF_, POUT_]


      ! FIXME only period for now, so hack bc arrays
      bc_ilo_type = -1 
      bc_ihi_type = -1  
      bc_jlo_type = -1 
      bc_jhi_type = -1 
      
      
      idx = one / dx(1)
      idy = one / dx(2)
      
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)

            ! ****************************************************
            ! West face
            ! ****************************************************
            
            ! In the case of MINF, NSW, FSW, PSW we are using the prescribed Dirichlet value
            ! In the case of PINF, POUT          we are using the upwind value
            if (i.eq.domlo(1) .and. any(bc_ilo_type(j,1) == bc_list ) ) then
               u_w =  vel(i-1,j,1)
            else
               upls  = vel(i  ,j,1) - half * xslopes(i  ,j,1)
               umns  = vel(i-1,j,1) + half * xslopes(i-1,j,1)

               u_w   = upwind( umns, upls, u(i,j) )
            endif

            ! ****************************************************
            ! East face
            ! ****************************************************
            
            ! In the case of MINF, NSW, FSW, PSW we are using the prescribed Dirichlet value
            ! In the case of PINF, POUT          we are using the upwind value
            if (i.eq.domhi(1) .and. any(bc_ihi_type(j,1) == bc_list ) ) then
               u_e =  vel(i+1,j,1)
            else
               upls  = vel(i+1,j,1) - half * xslopes(i+1,j,1)
               umns  = vel(i  ,j,1) + half * xslopes(i  ,j,1)

               u_e   = upwind( umns, upls, u(i+1,j) )
            endif

            ! ****************************************************
            ! South face
            ! ****************************************************
            
            ! In the case of MINF, NSW, FSW, PSW we are using the prescribed Dirichlet value
            ! In the case of PINF, POUT          we are using the upwind value
            if (j.eq.domlo(2) .and. any(bc_jlo_type(i,1) == bc_list ) ) then
               u_s =  vel(i,j-1,1)
            else
               upls  = vel(i,j  ,1) - half * yslopes(i,j  ,1)
               umns  = vel(i,j-1,1) + half * yslopes(i,j-1,1)

               u_s   = upwind( umns, upls, v(i,j) )
            endif

            ! ****************************************************
            ! North face
            ! ****************************************************

            ! In the case of MINF, NSW, FSW, PSW we are using the prescribed Dirichlet value
            ! In the case of PINF, POUT          we are using the upwind value
            if (j.eq.domhi(2) .and.  any(bc_jhi_type(i,1) == bc_list ) ) then
               u_n =  vel(i,j+1,1)
            else
               upls  = vel(i,j+1,1) - half * yslopes(i,j+1,1)
               umns  = vel(i,j  ,1) + half * yslopes(i,j  ,1)

               u_n   = upwind( umns, upls, v(i,j+1) )
            endif

            ! ****************************************************
            ! Define convective terms -- conservatively
            !   divuc =  div(u^MAC c_edge) 
            ! ****************************************************
            

            divuc(i,j,1) = (u(i+1,j) * u_e - u(i,j) * u_w) * idx + &
                           (v(i,j+1) * u_n - v(i,j) * u_s) * idy
                               

            ! ****************************************************
            ! Return the negative
            ! ****************************************************

            !divuc(i,j,1) = -divuc(i,j,1)

         end do
      end do

   end subroutine compute_divuc


   ! Upwind along direction normal to velocity component
   function upwind_normal ( umns, upls ) result (ev)

      ! Small value to protect against tiny velocities used in upwinding
      real(ar),        parameter     :: small_vel = 1.0d-10

      real(ar), intent(in) :: umns, upls
      real(ar)             :: ev, avg

      if ( umns < zero .and. upls > zero ) then
         ev = zero
      else
         avg = half * ( upls + umns )
         if ( abs(avg) .lt. small_vel ) then
            ev = zero
         else
            ev = merge ( umns, upls, avg >= zero )
         end if
      end if

   end function upwind_normal

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

end module convection_mod
