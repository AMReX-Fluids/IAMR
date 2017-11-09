module iamr_mlmg_3d_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: iamr_mac_coef, iamr_mac_rhs

contains

  subroutine iamr_mac_coef (xlo, xhi, ylo, yhi, zlo, zhi, bx, bxlo, bxhi, by, bylo, byhi, &
       bz, bzlo, bzhi, rho, rlo, rhi, scale) &
       bind(c,name='iamr_mac_coef')
    integer, dimension(3), intent(in) :: xlo, xhi, ylo, yhi, zlo, zhi, bxlo, bxhi, bylo, byhi, &
         bzlo, bzhi, rlo, rhi
    real(amrex_real), intent(in) :: scale
    real(amrex_real), intent(in) :: rho(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(inout) :: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
    real(amrex_real), intent(inout) :: by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
    real(amrex_real), intent(inout) :: bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))

    integer :: i, j, k

    do       k = xlo(3), xhi(3)
       do    j = xlo(2), xhi(2)
          do i = xlo(1), xhi(1)
             bx(i,j,k) = (2.d0*scale) / (rho(i-1,j,k)+rho(i,j,k))
          end do
       end do
    end do

    do       k = ylo(3), yhi(3)
       do    j = ylo(2), yhi(2)
          do i = ylo(1), yhi(1)
             by(i,j,k) = (2.d0*scale) / (rho(i,j-1,k)+rho(i,j,k))
          end do
       end do
    end do

    do       k = zlo(3), zhi(3)
       do    j = zlo(2), zhi(2)
          do i = zlo(1), zhi(1)
             bz(i,j,k) = (2.d0*scale) / (rho(i,j,k-1)+rho(i,j,k))
          end do
       end do
    end do
  end subroutine iamr_mac_coef


  subroutine iamr_mac_rhs (lo, hi, rhs, rlo, rhi, ux, uxlo, uxhi, uy, uylo, uyhi, &
       uz, uzlo, uzhi, dxinv) bind(c,name='iamr_mac_rhs')
    integer, dimension(3) :: lo, hi, rlo, rhi, uxlo, uxhi, uylo, uyhi, uzlo, uzhi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(in   ) :: ux (uxlo(1):uxhi(1),uxlo(2):uxhi(2),uxlo(3):uxhi(3))
    real(amrex_real), intent(in   ) :: uy (uylo(1):uyhi(1),uylo(2):uyhi(2),uylo(3):uyhi(3))
    real(amrex_real), intent(in   ) :: uz (uzlo(1):uzhi(1),uzlo(2):uzhi(2),uzlo(3):uzhi(3))

    integer :: i, j, k
    real(amrex_real) :: divu

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             divu = dxinv(1) * (ux(i+1,j,k) - ux(i,j,k))  &
                  + dxinv(2) * (uy(i,j+1,k) - uy(i,j,k))  &
                  + dxinv(3) * (uz(i,j,k+1) - uz(i,j,k))
             rhs(i,j,k) = rhs(i,j,k) - divu
          end do
       end do
    end do
  end subroutine iamr_mac_rhs

end module iamr_mlmg_3d_module
