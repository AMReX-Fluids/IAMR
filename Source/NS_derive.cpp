#include <NavierStokesBase.H>
#include "NS_derive.H"


using namespace amrex;



void der_vel_avg (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real time, const int* /*bcrec*/, int level)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(datfab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(datfab.nComp() >= BL_SPACEDIM*2);
    AMREX_ASSERT(ncomp == BL_SPACEDIM*2);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);


    amrex::Real inv_time = NavierStokesBase::time_avg[level];

    amrex::ParallelFor(bx, BL_SPACEDIM, [inv_time,der,in_dat]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        der(i,j,k,n) = in_dat(i,j,k,n) * inv_time;
        der(i,j,k,n+BL_SPACEDIM) = sqrt(in_dat(i,j,k,n+BL_SPACEDIM) * inv_time);
    });

}

