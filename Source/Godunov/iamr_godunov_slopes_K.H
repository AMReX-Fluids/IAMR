#ifndef IAMR_CONVECTION_K_H_
#define IAMR_CONVECTION_K_H_

#include <AMReX_Slopes_K.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EB_slopes_K.H>
#else
#include <AMReX_FArrayBox.H>
#endif

namespace {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_ho_xslope (int i, int j, int k, int n,
                              amrex::Array4<amrex::Real const> const& q) noexcept
{
    int order = 4; 
    return amrex_calc_xslope(i,j,k,n,order,q);
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_ho_xslope_extdir (int i, int j, int k, int n,
                                     amrex::Array4<amrex::Real const> const& q,
                                     bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    int order = 4; 
    return amrex_calc_xslope_extdir(i,j,k,n,order,q,edlo,edhi,domlo,domhi);
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_ho_yslope (int i, int j, int k, int n,
                              amrex::Array4<amrex::Real const> const& q) noexcept
{
    int order = 4; 
    return amrex_calc_yslope(i,j,k,n,order,q);
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_ho_yslope_extdir (int i, int j, int k, int n,
                                     amrex::Array4<amrex::Real const> const& q,
                                     bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    int order = 4; 
    return amrex_calc_yslope_extdir(i,j,k,n,order,q,edlo,edhi,domlo,domhi);
}

#if (AMREX_SPACEDIM==3)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_ho_zslope (int i, int j, int k, int n,
                              amrex::Array4<amrex::Real const> const& q) noexcept
{
    int order = 4; 
    return amrex_calc_zslope(i,j,k,n,order,q);
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_ho_zslope_extdir (int i, int j, int k, int n,
                                     amrex::Array4<amrex::Real const> const& q,
                                     bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    int order = 4; 
    return amrex_calc_zslope_extdir(i,j,k,n,order,q,edlo,edhi,domlo,domhi);
}
#endif
}
#endif
