
#define SDIM 2

module divop_mod

   use amrex_fort_module,       only: ar => amrex_real
   use iso_c_binding ,          only: c_int
   use param,                   only: zero, half, one, my_huge
   use amrex_error_module,      only: amrex_abort
   use amrex_ebcellflag_module, only: is_covered_cell, is_single_valued_cell, &
        &                             get_neighbor_cells

   implicit none

contains

   ! Compute the divergence operator for EB geometries
   !
   ! OUTPUTS:
   !
   !   div    the cell-centered divergence
   !
   ! INPUTS:
   !
   !   fx       the fluxes at the x-face CENTER (not centroid!)
   !   fy       the fluxes at the y-face CENTER (not centroid!)
   !   fz       the fluxes at the z-face CENTER (not centroid!)
   !
   ! WARNING: fx, fy, fz HAS to be filled with at least 3 GHOST nodes
   !
   !
   !
   subroutine compute_divop ( lo, hi, &
                             div, dlo, dhi, &
                             vel,vllo,vlhi, &
                             fx, fxlo, fxhi, &
                             fy, fylo, fyhi, &
                             afrac_x, axlo, axhi, &
                             afrac_y, aylo, ayhi, &
                             cent_x,  cxlo, cxhi, &
                             cent_y,  cylo, cyhi, &
                             flags,    flo,  fhi, &
                             vfrac,   vflo, vfhi, &
                             bcent,    blo,  bhi, &
                             domlo, domhi,        &
                             dx, ng, nc, eta,          &
                             do_explicit_diffusion ) bind(C)

     ! use bc ! from incflo/src/boundary_conditions
     ! bc module contains lots of stuff; divop only needs cyclic?
     !use eb_wallflux_mod,      only: compute_diff_wallflux

      ! Tile bounds (cell centered)
      integer(c_int),  intent(in   ) :: lo(SDIM),  hi(SDIM)

      ! Array Bounds
      integer(c_int),  intent(in   ) :: dlo(SDIM), dhi(SDIM)
      integer(c_int),  intent(in   ) :: vllo(SDIM), vlhi(SDIM)
      integer(c_int),  intent(in   ) :: fxlo(SDIM), fxhi(SDIM)
      integer(c_int),  intent(in   ) :: fylo(SDIM), fyhi(SDIM)
      integer(c_int),  intent(in   ) :: axlo(SDIM), axhi(SDIM)
      integer(c_int),  intent(in   ) :: aylo(SDIM), ayhi(SDIM)
      integer(c_int),  intent(in   ) :: cxlo(SDIM), cxhi(SDIM)
      integer(c_int),  intent(in   ) :: cylo(SDIM), cyhi(SDIM)
      integer(c_int),  intent(in   ) ::  flo(SDIM),  fhi(SDIM)
      integer(c_int),  intent(in   ) :: vflo(SDIM), vfhi(SDIM)
      integer(c_int),  intent(in   ) ::  blo(SDIM),  bhi(SDIM)

      integer(c_int),  intent(in   ) :: domlo(SDIM), domhi(SDIM)

      ! Number of ghost cells
      integer(c_int),  intent(in   ) :: nc, ng

      ! Grid
      real(ar),        intent(in   ) :: dx(SDIM)

      ! Arrays
      real(ar),        intent(in   ) ::                       &
           & fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),nc), &
           & fy(fylo(1):fyhi(1),fylo(2):fyhi(2),nc), &
           & afrac_x(axlo(1):axhi(1),axlo(2):axhi(2)), &
           & afrac_y(aylo(1):ayhi(1),aylo(2):ayhi(2)), &
           & cent_x(cxlo(1):cxhi(1),cxlo(2):cxhi(2),2),&
           & cent_y(cylo(1):cyhi(1),cylo(2):cyhi(2),2),&
           & vfrac(vflo(1):vfhi(1),vflo(2):vfhi(2)),   &
           &   vel(vllo(1):vlhi(1),vllo(2):vlhi(2),nc), &
           & bcent(blo(1):bhi(1),blo(2):bhi(2),3)

      ! Optional arrays (only for viscous calculations)
      real(ar),        intent(in   ), optional  ::                &
           &     eta(vflo(1):vfhi(1),vflo(2):vfhi(2))

      real(ar),        intent(inout) ::                           &
           & div(dlo(1):dhi(1),dlo(2):dhi(2),nc)

      ! BC types
      integer(c_int), intent(in   ) ::  &
           & flags(flo(1):fhi(1),flo(2):fhi(2))

      ! If true  then we include all the diffusive terms in this explicit result
      ! If false then we include all only the off-diagonal terms here -- we do this
      !     by computing the full tensor then subtracting the diagonal terms
      integer(c_int),  intent(in   ), optional :: do_explicit_diffusion
      
      ! Conservative div and EB stuff
      real(ar)  ::    &
           &  divc(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2), &
           & optmp(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2), &
           &  delm(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2), &
           &  mask(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2)

      ! FIXME -- pulled this out of bc_mod.f90 for now
      ! Flags for periodic boundary conditions
      logical :: cyclic_x = .true.
      logical :: cyclic_y = .true.

      real(ar), allocatable          :: divdiff_w(:,:)
      integer(c_int)                 :: i, j, n, nbr(-1:1,-1:1)
      integer(c_int)                 :: nwalls
      real(ar)                       :: idx, idy
      logical                        :: is_viscous

      idx = one / dx(1)
      idy = one / dx(2)

      ! Check number of ghost cells
      if (ng < 5) call amrex_abort( "compute_divop(): ng must be >= 5")

      ! Check if we are computing divergence for viscous term
      if ( present(eta) ) then
         is_viscous = .true.
      else
         is_viscous = .false.
      end if

      if ( abs(dx(1) - dx(2)) > epsilon(0.0_ar) ) then
         call amrex_abort("compute_divop(): grid spacing must be uniform")
      end if

      !
      ! Allocate arrays to host viscous wall fluxes
      !
      nwalls = 0
      if (is_viscous) then
         do j = lo(2)-2, hi(2)+2
            do i = lo(1)-2, hi(1)+2
               if (is_single_valued_cell(flags(i,j))) nwalls = nwalls + 1
            end do
         end do
         ! should this be SDIM instead of 3?
         allocate( divdiff_w(3,nwalls) )
         divdiff_w = zero
      end if

      !
      ! Array "mask" is used to sever the link to ghost cells when the BCs are not periodic
      ! It is set to 1 when the a cell can be used in computations, 0 otherwise
      !
      do j = lo(2)-2, hi(2)+2
         do i = lo(1)-2, hi(1)+2
            if ( ( .not. cyclic_x .and. (i < domlo(1) .or. i > domhi(1)) ) .or. &
                 ( .not. cyclic_y .and. (j < domlo(2) .or. j > domhi(2)) ) ) then
               mask(i,j) = zero
            else
               mask(i,j) = one
            end if
         end do
      end do

      !
      ! We use the EB algorithmm to compute the divergence at cell centers
      !
      ncomp_loop: do n = 1, nc

         !
         ! Step 1: compute conservative divergence on stencil (lo-2,hi+2)
         !
         compute_divc: block
           real(ar) :: fxp, fxm, fyp, fym, fzp, fzm
           integer  :: iwall

           iwall = 0

           divc = zero

           do j = lo(2)-2, hi(2)+2
              do i = lo(1)-2, hi(1)+2

                 if (is_covered_cell(flags(i,j))) then

                    divc(i,j) = 0.0d0 !my_huge

                 else if (is_single_valued_cell(flags(i,j))) then

                    call get_neighbor_cells( flags(i,j), nbr )

                    ! interp_to_face_centroid returns the proper flux multiplied
                    ! by the face area
                    fxp = interp_to_face_centroid( i+1, j, 1, fx, fxlo, n,  &
                         & afrac_x, axlo, cent_x, cxlo, nbr )

                    fxm = interp_to_face_centroid( i  , j, 1, fx, fxlo, n,  &
                         & afrac_x, axlo, cent_x, cxlo, nbr )

                    fyp = interp_to_face_centroid( i, j+1, 2, fy, fylo, n,  &
                         & afrac_y, aylo, cent_y, cylo, nbr )

                    fym = interp_to_face_centroid( i, j, 2, fy, fylo, n,  &
                         & afrac_y, aylo, cent_y, cylo, nbr )

                    divc(i,j) = ( ( fxp - fxm ) * idx + &
                         &        ( fyp - fym ) * idy ) / vfrac(i,j)
  !write(*,*) 'DEBUG in DIVOP ',i,j,  divc(i,j),fxp,fxm,fyp , fym
      

                    ! Add viscous wall fluxes (compute three components only
                    ! during the first pass, i.e. for n=1
                    iwall = iwall + 1
                    if (is_viscous) then
                       call amrex_abort("compute_divop(): not compatable with viscosity yet")
                       ! if (n==1) then
                       !    call compute_diff_wallflux(divdiff_w(:,iwall), &
                       !                               dx, i, j, k, &
                       !                               vel, vllo, vlhi, &
                       !                               eta, vflo, vfhi, &
                       !                               bcent, blo, bhi, &
                       !                               afrac_x, axlo, axhi, &
                       !                               afrac_y, aylo, ayhi, &
                       !                               afrac_z, azlo, azhi, &
                       !                               do_explicit_diffusion)
                       ! end if
                       ! divc(i,j,k) = divc(i,j,k) - divdiff_w(n,iwall) / (dx(n) * vfrac(i,j,k))
                    end if

                 else

                    divc(i,j) = (   fx(i+1,j  ,n) - fx(i,j,n) ) * idx  &
                         &      + ( fy(i  ,j+1,n) - fy(i,j,n) ) * idy
                         
   !write(*,*) 'DEBUG in DIVOP ',i,j,  divc(i,j),fx(i+1,j  ,n) , fx(i,j,n),fy(i  ,j+1,n) , fy(i,j,n)

                 end if

              end do
           end do
         end block compute_divc

         !
         ! Step 2: compute delta M ( mass gain or loss ) on (lo-1,hi+1)
         !
         optmp = zero
         block
           integer   :: ii, jj
           real(ar)  :: divnc, vtot
           real(ar)  :: vfrac_mask

           do j = lo(2)-1, hi(2)+1
              do i = lo(1)-1, hi(1)+1
                 if (is_single_valued_cell(flags(i,j))) then
                    divnc = zero
                    vtot  = zero
                    call get_neighbor_cells(flags(i,j),nbr)
                    do jj = -1, 1
                       do ii = -1, 1
                          ! Check if we have to include also cell (i,j) itself
                          if ( ( ii /= 0 .or. jj /= 0) &
                               .and. (nbr(ii,jj)==1) ) then
                             vfrac_mask  = vfrac(i+ii,j+jj) * mask(i+ii,j+jj)
                             vtot        = vtot  + vfrac_mask
                             divnc       = divnc + vfrac_mask * divc(i+ii,j+jj)
                          end if
                       end do
                    end do
                    divnc = divnc / vtot
                    optmp(i,j) = (one - vfrac(i,j)) * ( divnc - divc(i,j) )
                    !write(*,*) 'DEBUG in DIVOP ',i,j,vfrac(i,j) , divnc, divc(i,j)
                    delm(i,j) = - vfrac(i,j) * optmp(i,j)
                 else
                    delm(i,j) = zero
                 end if
                 
              end do
           end do

         end block

         !
         ! Step 3: redistribute excess/loss of mass
         !
         block
           real(ar) :: wtot
           integer  :: ii, jj

           do j = lo(2)-1, hi(2)+1
              do i = lo(1)-1, hi(1)+1
                 if (is_single_valued_cell(flags(i,j))) then
                    wtot = zero
                    call get_neighbor_cells(flags(i,j),nbr)
                    do jj = -1, 1
                       do ii = -1, 1
                          ! Check if we have to include also cell (i,j,k) itself
                          if ( ( ii /= 0 .or. jj /= 0 ) &
                               .and. (nbr(ii,jj)==1) ) then
                             wtot = wtot + vfrac(i+ii,j+jj) * mask(i+ii,j+jj)
                          end if
                       end do
                    end do

                    wtot = one/wtot

                    do jj = -1, 1
                       do ii = -1, 1
                          ! Check if we have to include also cell (i,j,k) itself
                          if ( ( ii /= 0 .or. jj /= 0) &
                               .and. (nbr(ii,jj)==1) ) then
                             optmp(i+ii,j+jj) = optmp(i+ii,j+jj) &
                                  + delm(i,j) * wtot &
                                  * mask(i+ii,j+jj)
                          end if
                       end do
                    end do
                 end if
              end do
           end do
         end block

         ! ****************************************************
         ! Resume the correct sign, AKA return the negative
         ! ****************************************************
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               div(i,j,n) = divc(i,j) + optmp(i,j)
               !write(*,*) 'DEBUG in DIVOP ',i,j,divc(i,j),optmp(i,j)
            end do
         end do


      end do ncomp_loop

      !
      ! Delete working arrays
      ! 
      if (is_viscous) deallocate(divdiff_w)

    end subroutine compute_divop

       !
   ! Returns the flux at the face centroid ALREADY multiplied by the face area
   !
   function interp_to_face_centroid ( i, j, dir, var, vlo,  n,  &
                                     afrac, alo, cent, clo, nbr )  result(ivar)

     use amrex_eb_util_module, only: amrex_eb_interpolate_to_face_centroid_per_cell

           ! Face indeces: these must be consistent with a staggered indexing
      ! and therefore consistent with the value of dir
      integer(c_int),  intent(in   ) :: i, j

      ! Direction of staggering (1=x, 2=y, 3=z): this specify how (i,j,k) must
      ! be interpreted, i.e. which staggered numbering the indexing refer to
      integer(c_int),  intent(in   ) :: dir

      ! The component to interpolate
      integer(c_int),  intent(in   ) :: n

      ! Array Bounds ( only start index )
      integer(c_int),  intent(in   ) :: vlo(SDIM), alo(SDIM), clo(SDIM)

      ! Arrays
      real(ar),        intent(in   ) ::           &
           &   var(vlo(1):, vlo(2):,1:), &
           & afrac(alo(1):, alo(2):),    &
           &  cent(clo(1):, clo(2):,1:)

      ! Neighbors information
      integer(c_int),  intent(in   ) :: nbr(-1:1,-1:1)

      ! Output: the interpolated value
      real(ar)                       :: ivar

      ! Local variables
      real(ar)                       :: fracx, fracy


      ivar = amrex_eb_interpolate_to_face_centroid_per_cell( &
           i, j, dir, var, vlo,  n, afrac, alo, cent, clo, nbr )
      
      !
      ! Return the flux multiplied by the face area
      !
      ivar = ivar * afrac(i,j)

    end function interp_to_face_centroid

  end module divop_mod
