module iamr_mlmg_2d_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: iamr_mac_coef, iamr_mac_rhs

contains

  subroutine iamr_mac_coef (xlo, xhi, ylo, yhi, bx, bxlo, bxhi, by, bylo, byhi, rho, rlo, rhi, scale) &
       bind(c,name='iamr_mac_coef')
    integer, dimension(2), intent(in) :: xlo, xhi, ylo, yhi, bxlo, bxhi, bylo, byhi, rlo, rhi
    real(amrex_real), intent(in) :: scale
    real(amrex_real), intent(in) :: rho(rlo(1):rhi(1),rlo(2):rhi(2))
    real(amrex_real), intent(inout) :: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2))
    real(amrex_real), intent(inout) :: by(bylo(1):byhi(1),bylo(2):byhi(2))

    integer :: i, j

    do    j = xlo(2), xhi(2)
       do i = xlo(1), xhi(1)
          bx(i,j) = (2.d0*scale) / (rho(i-1,j)+rho(i,j))
       end do
    end do

    do    j = ylo(2), yhi(2)
       do i = ylo(1), yhi(1)
          by(i,j) = (2.d0*scale) / (rho(i,j-1)+rho(i,j))
       end do
    end do
  end subroutine iamr_mac_coef


  subroutine iamr_mac_rhs (lo, hi, rhs, rlo, rhi, ux, uxlo, uxhi, uy, uylo, uyhi, &
       ax, axlo, axhi, ay, aylo, ayhi, vol, vlo, vhi) bind(c,name='iamr_mac_rhs')
    integer, dimension(2) :: lo, hi, rlo, rhi, uxlo, uxhi, uylo, uyhi, axlo, axhi, aylo, ayhi, vlo, vhi
    real(amrex_real), intent(inout) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2))
    real(amrex_real), intent(in   ) :: ux (uxlo(1):uxhi(1),uxlo(2):uxhi(2))
    real(amrex_real), intent(in   ) :: uy (uylo(1):uyhi(1),uylo(2):uyhi(2))
    real(amrex_real), intent(in   ) :: ax (axlo(1):axhi(1),axlo(2):axhi(2))
    real(amrex_real), intent(in   ) :: ay (aylo(1):ayhi(1),aylo(2):ayhi(2))
    real(amrex_real), intent(in   ) :: vol( vlo(1): vhi(1), vlo(2): vhi(2))

    integer :: i, j
    real(amrex_real) :: divu

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          divu = (ax(i+1,j)*ux(i+1,j) - ax(i,j)*ux(i,j)   &
               +  ay(i,j+1)*uy(i,j+1) - ay(i,j)*uy(i,j)) / vol(i,j)
          rhs(i,j) = rhs(i,j) - divu
       end do
    end do
  end subroutine iamr_mac_rhs

end module iamr_mlmg_2d_module
