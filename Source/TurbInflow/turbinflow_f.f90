module turbinflow_module

  use amrex_fort_module, only: amrex_real

  implicit none

  logical, save :: turbinflow_initialized = .false.

  integer, save :: lenfname
  integer, save :: iturbfile(128)

  integer, save :: npboxcells(3)
  real(kind=amrex_real), save :: pboxlo(3), dx(3), dxinv(3)
  real(kind=amrex_real), pointer, save :: plane_time(:)
  integer, save :: turb_ncomp = 3

  integer, parameter :: nplane = 32
  real(kind=amrex_real), pointer, save :: sdata(:,:,:,:)
  real(kind=amrex_real), save :: szlo=0.d0, szhi=0.d0
  integer, save :: izlo, izhi

  integer, save :: isswirltype  ! 1=swirl, 0=periodic
  real(kind=amrex_real), save :: units_conversion = 100.d0  ! m --> cm & m/s --> cm/s

  private

  public :: turbinflow_initialized, init_turbinflow, get_turbstate, turb_ncomp

  interface
     subroutine getplane (filename, len, data, plane, ncomp, isswirltype) bind(c)
       use amrex_fort_module, only: amrex_real
       integer, intent(in) :: len, ncomp, isswirltype, plane
       integer, intent(in) :: filename(len)
       real(kind=amrex_real), intent(inout) :: data(*)
     end subroutine getplane
  end interface

!$omp threadprivate(turbinflow_initialized,plane_time,sdata,szlo,szhi,izlo,izhi)
  
contains

  subroutine init_turbinflow(turbfile, is_cgs)
    use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
    character (len=*), intent(in) :: turbfile
    logical, intent(in), optional :: is_cgs

    integer, parameter :: iunit = 20
    integer :: ierr, n, npts(3)
    real(kind=amrex_real) :: probsize(3), pboxsize(3), z
    integer :: is_periodic(3)

    !$omp critical
    if (present(is_cgs)) then
       if (is_cgs) then
          units_conversion = 100.d0
       else
          units_conversion = 1.d0
       end if
    end if

    lenfname = len_trim(turbfile)
    do n=1,lenfname
       iturbfile(n) = ichar(turbfile(n:n))
    end do

    ! read header
    open(iunit, file=trim(turbfile)//'/HDR', form='formatted', action='read', &
         status='old', iostat=ierr)
    
    if (ierr .ne. 0) then
       call bl_abort('Problem opening file: ' // trim(turbfile) // '/HDR')
    end if

    read(iunit,*) npts
    read(iunit,*) probsize
    read(iunit,*) is_periodic
    !read(iunit,*) turb_ncomp

    if (is_periodic(3) .eq. 0) then
       isswirltype = 1
       call amrex_allocate(plane_time,1,npts(3))
       do n=1,npts(3)
          read(iunit,*) plane_time(n)
       enddo
    else
       isswirltype = 0
    end if

    close(iunit)

    probsize = probsize * units_conversion

    dx = probsize / dble(npts-1)
    dxinv = 1.d0/dx

    pboxsize(1:2) = probsize(1:2) - 2.d0*dx(1:2)  ! because there is one ghost point on each side
    pboxsize(3) = probsize(3) ! no ghost point in z-direction

    npboxcells(1:2) = npts(1:2) - 3
    if (isswirltype .eq. 0) then
       npboxcells(3)   = npts(3) - 1
    else
       npboxcells(3)   = npts(3)
    endif
    ! period box covers -0.5*pboxsize(1) <= x <= 0.5*pboxsize(1)
    !                   -0.5*pboxsize(2) <= y <= 0.5*pboxsize(2)
    !                                  0 <= z <= pboxsize(3)
    pboxlo(1:2) = -0.5d0*pboxsize(1:2)
    pboxlo(3) = 0.d0
    !$omp end critical
    
    call amrex_allocate(sdata,1,npts(1),1,npts(2),1,nplane,1,turb_ncomp)

    if (isswirltype .eq. 0) then
       z = 0.d0
    else
       z = plane_time(1)
    endif
    call store_planes(z)

    turbinflow_initialized = .true.

  end subroutine init_turbinflow


  subroutine get_turbstate(lo1,lo2,hi1,hi2,x,y,z,s)
    integer, intent(in) :: lo1,lo2,hi1,hi2
    real(kind=amrex_real), intent(in) :: x(lo1:hi1), y(lo2:hi2)
    real(kind=amrex_real), intent(in) :: z
    real(kind=amrex_real), intent(out) :: s(lo1:hi1,lo2:hi2,turb_ncomp)

    integer :: i, j, k, n, i0, j0, k0, ii, jj
    real(kind=amrex_real) :: xx, yy, zz, zdata(0:2,0:2), ydata(0:2)
    real(kind=amrex_real) :: cx(0:2), cy(0:2), cz(0:2), zmin, zmax, z1, z2, z3

    if (.not. turbinflow_initialized) call bl_error("turbinflow module uninitialized")

    if (isswirltype.eq.1) then
       zmin = plane_time(izlo) + 0.5d0 * (plane_time(izlo+1) - plane_time(izlo))
       zmax = plane_time(izhi) - 0.5d0 * (plane_time(izhi) - plane_time(izhi-1))
       ! Special fix to allow interp up to edge of data
       if (izlo .eq. 1)             zmin = plane_time(izlo)
       if (izlo .eq. npboxcells(3)) zmax = plane_time(npboxcells(3))
       if (z.lt.zmin .or. z.gt.zmax) then
          call store_planes(z)
       end if
       do n = izlo,izhi-2
          zmin = plane_time(n  ) + 0.5d0 * (plane_time(n+1) - plane_time(n  ))
          zmax = plane_time(n+2) - 0.5d0 * (plane_time(n+2) - plane_time(n+1))

          ! Special fix to allow interp up to edge of data
          if (n .eq. 1)               zmin = plane_time(izlo)
          if (n+2 .eq. npboxcells(3)) zmax = plane_time(npboxcells(3))
          
          if (z.ge.zmin .and. z.le.zmax) then
             k0 = n - izlo + 1
             z1 = plane_time(n)
             z2 = plane_time(n+1)
             z3 = plane_time(n+2)
             
             cz(0) = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
             cz(1) = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
             cz(2) = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
          end if
       enddo
    else
       if (z.lt.szlo+0.5d0*dx(3) .or. z.gt.szhi-0.5d0*dx(3)) then
          call store_planes(z)
       end if
       zz = (z-szlo)*dxinv(3)
       k0 = nint(zz) - 1
       zz = zz - dble(k0)
       cz(0) = 0.5d0*(zz-1.d0)*(zz-2.d0)
       cz(1) = zz*(2.d0-zz)
       cz(2) = 0.5d0*zz*(zz-1.d0)
       k0 = k0 + 1 ! because it's Fortran
       k0 = min(max(k0,1),nplane-2)
    endif

    do n=1,turb_ncomp
       do j=lo2,hi2
          yy = (y(j)-pboxlo(2))*dxinv(2)
          j0 = nint(yy) - 1
          yy = yy - dble(j0)
          cy(0) = 0.5d0*(yy-1.d0)*(yy-2.d0)
          cy(1) = yy*(2.d0-yy)
          cy(2) = 0.5d0*yy*(yy-1.d0)
          j0 = modulo(j0, npboxcells(2)) + 2 ! +2 because Fortran index starts with 1 and there is a ghost point
          do i=lo1,hi1
             xx = (x(i)-pboxlo(1))*dxinv(1)
             i0 = nint(xx) - 1
             xx = xx - dble(i0)
             cx(0) = 0.5d0*(xx-1.d0)*(xx-2.d0)
             cx(1) = xx*(2.d0-xx)
             cx(2) = 0.5d0*xx*(xx-1.d0)
             i0 = modulo(i0, npboxcells(1)) + 2 ! +2 as j0
             
             do jj=0,2
                do ii=0,2
                   zdata(ii,jj) = cz(0)*sdata(i0+ii,j0+jj,k0  ,n) &
                        +         cz(1)*sdata(i0+ii,j0+jj,k0+1,n) &
                        +         cz(2)*sdata(i0+ii,j0+jj,k0+2,n)
                end do
             end do
             
             do ii=0,2
                ydata(ii) = cy(0)*zdata(ii,0) + cy(1)*zdata(ii,1) + cy(2)*zdata(ii,2)
             end do

             s(i,j,n) = cx(0)*ydata(0) + cx(1)*ydata(1) + cx(2)*ydata(2)

          end do
       end do
    end do

  end subroutine get_turbstate

  subroutine store_planes(z)
    real(kind=amrex_real), intent(in) :: z
    real(kind=amrex_real) :: zmin, zmax
    integer :: iplane, k, n, zimin, zimax

    if (isswirltype.eq.0) then
       izlo = nint(z*dxinv(3)) - 1
       izhi = izlo + nplane - 1
       szlo = izlo*dx(3)
       szhi = szlo + dble(nplane-1)*dx(3)
    else
       zimin = 1
       zimax = nplane
       szlo = -1.d0
       !do n=1,npboxcells(3)-nplane+1
       do n=npboxcells(3)-nplane+1,1,-1
          if (szlo .lt. 0) then
             zimin = n
             zimax = min(n+nplane-1,npboxcells(3))
             zmin = 0.5d0 * (plane_time(zimin) + plane_time(zimin+1))
             zmax = 0.5d0 * (plane_time(zimax) + plane_time(zimax-1))
             if (zimin .eq. 1)             zmin = plane_time(1)
             if (zimax .eq. npboxcells(3)) zmin = plane_time(npboxcells(3))
             if (z.ge.zmin .and. z.le.zmax) then
                izlo = n
                izhi = n + nplane - 1
                szlo = plane_time(zimin)
                szhi = plane_time(zimax)
             end if
          end if
       enddo       
       if (szlo .lt. 0.d0) call bl_pd_abort('Bad logic in plane reader')
    endif
    do n=1,turb_ncomp
       do iplane=1,nplane
          if (isswirltype .eq. 1) then
             k = iplane + izlo - 1
          else
             k = modulo(izlo+iplane-1, npboxcells(3)) + 1
          end if
          
          call getplane(iturbfile, lenfname, sdata(:,:,iplane,n), k, n, isswirltype)
          sdata(:,:,iplane,n) = sdata(:,:,iplane,n)*units_conversion
       end do
    end do
  end subroutine store_planes

end module turbinflow_module
