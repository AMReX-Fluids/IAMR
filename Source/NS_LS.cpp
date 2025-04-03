
#include <NS_LS.H>

#include <AMReX_ParmParse.H>
#include <AMReX_TagBox.H>
#include <AMReX_Utility.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_FillPatchUtil.H>
#include <iamr_constants.H>

using namespace amrex;

// ls related

Real calculate_eps (const Geometry& geom, int epsilon)
{
    const Real* dx    = geom.CellSize();
    Real dxmin        = dx[0];
    for (int d=1; d<AMREX_SPACEDIM; ++d) {
        dxmin = std::min(dxmin,dx[d]);
    }
    Real eps = epsilon * dxmin;

    return eps;
}

void
phi_to_heavi(const Geometry& geom, int epsilon, MultiFab& phi, MultiFab& heaviside)
{

    Print() << "In the phi_to_heavi " << std::endl;
    
    const Real pi     = 3.141592653589793238462643383279502884197;
    Real eps = calculate_eps(geom, epsilon);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();
        auto const& phifab   = phi.array(mfi);
        auto const& heavifab = heaviside.array(mfi);
        amrex::ParallelFor(bx, [phifab, heavifab, pi, eps]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {

            if (phifab(i,j,k) > eps) {
                heavifab(i,j,k) = 1.0;
            } else if (phifab(i,j,k) > -eps) {
                heavifab(i,j,k) = 0.5 * (1.0 + phifab(i,j,k) / eps + 1.0 / pi * std::sin(phifab(i,j,k) * pi / eps));
            } else {
                heavifab(i,j,k) = 0.0;
            }

        });
    }
}

void
heavi_to_rhoormu(MultiFab& heaviside, Real var1, Real var2, MultiFab& outmf)
{

    amrex::Print() << "In the heavi_to_rhomu " << std::endl;
    
    BL_ASSERT(heaviside.nComp() == outmf.nComp());

    int ncomp = outmf.nComp();
    int ngrow = outmf.nGrow();

    // build heaviside_temp because we might need to smooth the heaviside function
    MultiFab heaviside_temp(heaviside.boxArray(), heaviside.DistributionMap(), ncomp, ngrow);
    MultiFab::Copy(heaviside_temp, heaviside, 0, 0, ncomp, ngrow);

    // if( smoothforrhomu==1 && (parent->levelSteps(0)%lev0step_of_smoothforrhomu == 0) ){
    // smooth_sf(heaviside_temp);
    // }

    heaviside_temp.mult(var1-var2, 0, ncomp, ngrow);

    MultiFab rtmp(heaviside.boxArray(), heaviside.DistributionMap(), ncomp, ngrow);
    rtmp.setVal(var2, 0, ncomp, ngrow);

    MultiFab::Add(heaviside_temp, rtmp, 0, 0, ncomp, ngrow);
    MultiFab::Copy(outmf, heaviside_temp, 0, 0, ncomp, ngrow);

}

Real calculate_eps_one (const Geometry& geom, int reinit_levelset)
{

    Print() << "In the calculate_eps_one " << std::endl;

    Real eps_one = 0.0;

    const Real* dx    = geom.CellSize();
    Real dxmax        = dx[0];
    for (int d=1; d<AMREX_SPACEDIM; ++d) {
        dxmax = std::max(dxmax,dx[d]);
    }

    if (reinit_levelset == 2)
    {
        // improved conservative level set, JCP, 2007
        eps_one = 0.5 * std::pow(dxmax, 0.9);
    }
    else if (reinit_levelset == 3)
    {
        // multiUQ, JCP, 2019/CPC, 2021
        eps_one = 0.875 * dxmax;
    }

    return eps_one;
}

Real calculate_eps_two (const Geometry& geom, int reinit_levelset)
{

    Print() << "In the calculate_eps_two " << std::endl;

    Real eps_two = 0.0;

    const Real* dx    = geom.CellSize();
    Real dxmax        = dx[0];
    for (int d=1; d<AMREX_SPACEDIM; ++d) {
        dxmax = std::max(dxmax,dx[d]);
    }

    if (reinit_levelset == 2)
    {
        // improved conservative level set, JCP, 2007
        eps_two = 0.0;
    }
    else if (reinit_levelset == 3)
    {
        // multiUQ, JCP, 2019/CPC, 2021
        eps_two = 0.125 * dxmax;
    }

    return eps_two;
}

void levelset_diffcomp (Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM>& phi_cc_grad,
                        MultiFab& phi_ctime,
                        MultiFab& phi1,
                        MultiFab& phi2, 
                        int epsG,
                        int epsG2)
{

    Print() << "In the levelset_diffcomp " << std::endl;

    // Step 1: calculate the comp term phi1

    // Step 2: calculate the diff terms phi2


}
