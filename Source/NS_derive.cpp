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
    amrex::Real inv_time;
    amrex::Real inv_time_fluct;

    if (NavierStokesBase::time_avg[level] == 0){
      inv_time = 1.0;
    }else{
      inv_time = 1.0 / NavierStokesBase::time_avg[level];
    }

    if (NavierStokesBase::time_avg_fluct[level] == 0){
      inv_time_fluct = 1.0;
    }else{
      inv_time_fluct = 1.0 / NavierStokesBase::time_avg_fluct[level];
    }

    amrex::ParallelFor(bx, BL_SPACEDIM, [inv_time,inv_time_fluct,der,in_dat]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        der(i,j,k,n) = in_dat(i,j,k,n) * inv_time;
        der(i,j,k,n+BL_SPACEDIM) = sqrt(in_dat(i,j,k,n+BL_SPACEDIM) * inv_time_fluct);
    });

}

