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
    AMREX_ASSERT(datfab.nComp() >= BL_SPACEDIM);
    AMREX_ASSERT(ncomp == BL_SPACEDIM);
    auto const in_dat = datfab.array();
    auto          der = derfab.array(dcomp);


    amrex::ParallelFor(bx, BL_SPACEDIM,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        der(i,j,k,n) = in_dat(i,j,k,n);
    });
    amrex::Print() << " " << std::endl;
    amrex::Print() << "DEBUG in der_vel_avg level = " << level << std::endl;
    amrex::Print() << "DEBUG in der_vel_avg time_avg = " << NavierStokesBase::time_avg[level] << std::endl;

    amrex::Print() << " " << std::endl;

}

