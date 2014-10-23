

module bpm
  use iso_c_binding
  implicit none

  ! if you change this, you MUST update BPM_F.H as well
  type, bind(c) :: bpm_parameters
     integer(c_int) :: type
     real(c_double) :: diameter
  end type bpm_parameters

contains

  subroutine bpm_compute_forces_2d(x1, v1, f1, x2, v2, f2, sconstant, slength) bind(c)
    real(c_double), intent(in   )        :: x1(2), v1(2), x2(2), v2(2)
    real(c_double), intent(inout)        :: f1(2), f2(2)
    real(c_double), intent(in   ), value :: sconstant, slength

    real(c_double) :: force(2), delta(2), distance

    ! compute spring force
    delta    = x2 - x1
    distance = sqrt(delta(1)**2 + delta(2)**2)
    if (abs(distance - slength) > 1.d-12) then
       force = -sconstant * (distance - slength) * delta / distance
    else
       force = 0
    end if

    print *, "bpm:: x1:", x1
    print *, "bpm:: x2:", x2
    print *, "bpm:: f: ", force

    f1 = f1 + force
  end subroutine bpm_compute_forces_2d

  subroutine bpm_spread_forces_2d(x1, force, params, lo, hi, ng, xlo, xhi, dx, f, divf) bind(c)
    real(c_double),       intent(in   )        :: x1(2), force(2), xlo(2), xhi(2)
    type(bpm_parameters), intent(in   )        :: params
    real(c_double),       intent(in   ), value :: dx
    integer(c_int),       intent(in   )        :: lo(2), hi(2)
    integer(c_int),       intent(in   ), value :: ng
    real(c_double),       intent(inout)        :: f(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,2)
    real(c_double),       intent(inout)        :: divf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)

    real(c_double), allocatable, dimension(:,:) :: distx, disty, dist2, blob, dbdr

    real(c_double) :: a, b, x, y
    integer        :: p, i, j, c, cell(2)

    p = ceiling(params%diameter/2 / dx) + 1

    allocate(distx(-p:p,-p:p), disty(-p:p,-p:p), dist2(-p:p,-p:p), blob(-p:p,-p:p), dbdr(-p:p,-p:p))

    ! compute cell that contains particle
    cell = floor((x1 - xlo) / dx + lo)

    ! compute distances to cell centers
    do i = -p, p
       x = xlo(1) + dx * (cell(1) + i + 0.5d0)
       do j = -p, p
          y = xlo(2) + dx * (cell(2) + j + 0.5d0)
          distx(i,j) = x - x1(1)
          disty(i,j) = y - x1(2)
       end do
    end do

    dist2 = (distx**2 + disty**2) / params%diameter**2

    select case (params%type)
    case (1)
       ! B(r) = (6/pi) (1-x^2)^5
       a =   6.0d0 / (3.1415926535897931d0 * params%diameter**2)
       b = -60.0d0 / (3.1415926535897931d0 * params%diameter**4)
       where (dist2 <= 1.0d0)
          blob = a * (1.0d0 - dist2)**5
          dbdr = b * (1.0d0 - dist2)**4
       elsewhere
          blob = 0
          dbdr = 0
       end where
    case default
       stop "INVALID BLOB TYPE"
    end select

    ! print *, "bpm:: blob volume", sum(blob)*dx**2
    ! do i = -p, p
    !    print *, blob(i,:)
    ! end do

    ! print *, "bpm:: dbdr"
    ! do i = -p, p
    !    print *, dbdr(i,:)
    ! end do

    ! spread force
    do c = 1, 2
       do i = -p, p
          do j = -p, p
             f(cell(1)+i,cell(2)+j,c) = f(cell(1)+i,cell(2)+j,c) + force(c) * blob(i,j)
          end do
       end do
    end do

    ! spread divergence
    do i = -p, p
       do j = -p, p
          divf(cell(1)+i,cell(2)+j) = divf(cell(1)+i,cell(2)+j) + &
               (force(1) * distx(i,j) + force(2) * disty(i,j)) * dbdr(i,j)
       end do
    end do

  end subroutine bpm_spread_forces_2d

end module bpm
