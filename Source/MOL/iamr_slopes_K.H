#ifndef IAMR_SLOPES_K_H_
#define IAMR_SLOPES_K_H_

#include <AMReX_FArrayBox.H>
// #include <AMReX_BCRec.H>


namespace {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_xslope (int i, int j, int k, int n,
                           amrex::Array4<amrex::Real const> const& vcc) noexcept
{
    amrex::Real dl = 2.0*(vcc(i  ,j,k,n) - vcc(i-1,j,k,n));
    amrex::Real dr = 2.0*(vcc(i+1,j,k,n) - vcc(i  ,j,k,n));
    amrex::Real dc = 0.5*(vcc(i+1,j,k,n) - vcc(i-1,j,k,n));
    amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_xslope_extdir (int i, int j, int k, int n,
                                  amrex::Array4<amrex::Real const> const& vcc,
                                  bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real dl = 2.0*(vcc(i  ,j,k,n) - vcc(i-1,j,k,n));
    amrex::Real dr = 2.0*(vcc(i+1,j,k,n) - vcc(i  ,j,k,n));
    amrex::Real dc = 0.5*(vcc(i+1,j,k,n) - vcc(i-1,j,k,n));
    if (edlo and i == domlo) {
        dc = (vcc(i+1,j,k,n)+3.0*vcc(i,j,k,n)-4.0*vcc(i-1,j,k,n))/3.0;
    } else if (edhi and i == domhi) {
        dc = (4.0*vcc(i+1,j,k,n)-3.0*vcc(i,j,k,n)-vcc(i-1,j,k,n))/3.0;
    }
    amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_yslope (int i, int j, int k, int n,
                           amrex::Array4<amrex::Real const> const& vcc) noexcept
{
    amrex::Real dl = 2.0*(vcc(i,j  ,k,n) - vcc(i,j-1,k,n));
    amrex::Real dr = 2.0*(vcc(i,j+1,k,n) - vcc(i,j  ,k,n));
    amrex::Real dc = 0.5*(vcc(i,j+1,k,n) - vcc(i,j-1,k,n));
    amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_yslope_extdir (int i, int j, int k, int n,
                                  amrex::Array4<amrex::Real const> const& vcc,
                                  bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real dl = 2.0*(vcc(i,j  ,k,n) - vcc(i,j-1,k,n));
    amrex::Real dr = 2.0*(vcc(i,j+1,k,n) - vcc(i,j  ,k,n));
    amrex::Real dc = 0.5*(vcc(i,j+1,k,n) - vcc(i,j-1,k,n));
    if (edlo and j == domlo) {
        dc = (vcc(i,j+1,k,n)+3.0*vcc(i,j,k,n)-4.0*vcc(i,j-1,k,n))/3.0;
    } else if (edhi and j == domhi) {
        dc = (4.0*vcc(i,j+1,k,n)-3.0*vcc(i,j,k,n)-vcc(i,j-1,k,n))/3.0;
    }
    amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

#if (AMREX_SPACEDIM==3)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_zslope (int i, int j, int k, int n,
                           amrex::Array4<amrex::Real const> const& vcc) noexcept
{
    amrex::Real dl = 2.0*(vcc(i,j,k  ,n) - vcc(i,j,k-1,n));
    amrex::Real dr = 2.0*(vcc(i,j,k+1,n) - vcc(i,j,k  ,n));
    amrex::Real dc = 0.5*(vcc(i,j,k+1,n) - vcc(i,j,k-1,n));
    amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_zslope_extdir (int i, int j, int k, int n,
                                  amrex::Array4<amrex::Real const> const& vcc,
                                  bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real dl = 2.0*(vcc(i,j,k  ,n) - vcc(i,j,k-1,n));
    amrex::Real dr = 2.0*(vcc(i,j,k+1,n) - vcc(i,j,k  ,n));
    amrex::Real dc = 0.5*(vcc(i,j,k+1,n) - vcc(i,j,k-1,n));
    if (edlo and k == domlo) {
        dc = (vcc(i,j,k+1,n)+3.0*vcc(i,j,k,n)-4.0*vcc(i,j,k-1,n))/3.0;
    } else if (edhi and k == domhi) {
        dc = (4.0*vcc(i,j,k+1,n)-3.0*vcc(i,j,k,n)-vcc(i,j,k-1,n))/3.0;
    }
    amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}
#endif

}


#endif
