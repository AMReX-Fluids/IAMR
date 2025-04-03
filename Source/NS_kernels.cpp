
#include <NS_kernels.H>

#include <AMReX_ParmParse.H>
#include <AMReX_TagBox.H>
#include <AMReX_Utility.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_FillPatchUtil.H>
#include <iamr_constants.H>

using namespace amrex;

// ls related
// Assume all gts are filled
// phi: cc variables, at least 2 gts
// phi_cc_grad: cc variables, 1 gt
void cc_to_cc_grad(Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM>& phi_cc_grad, const MultiFab& phi, const Geometry& geom,
                   int normalize)
{

    Print() << "In the cc_to_cc_grad " << std::endl;

    const BoxArray& ba = phi.boxArray();
    const DistributionMapping& dm = phi.DistributionMap();
    const int ncomp = phi.nComp();
    const int ngrow = phi.nGrow();
    BL_ASSERT(ngrow>=2);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        phi_cc_grad[idim] = std::make_unique<MultiFab>(ba, dm, ncomp, 1);
    }

    // 1. Calculate the gradient
    AMREX_D_TERM(const Real dxi = static_cast<Real>(geom.InvCellSize(0));,
                 const Real dyi = static_cast<Real>(geom.InvCellSize(1));,
                 const Real dzi = static_cast<Real>(geom.InvCellSize(2)););
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi, TilingIfNotGPU());  mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(1); // 1 gt
        const auto& s = phi.array(mfi);
        AMREX_D_TERM(const auto& gx = phi_cc_grad[0]->array(mfi);,
                     const auto& gy = phi_cc_grad[1]->array(mfi);,
                     const auto& gz = phi_cc_grad[2]->array(mfi););

        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            gx(i,j,k,n) = 0.5 * dxi*(s(i+1,j,k,n) - s(i-1,j,k,n));
        });
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            gy(i,j,k,n) = 0.5 * dyi*(s(i,j+1,k,n) - s(i,j-1,k,n));
        });
#if (AMREX_SPACEDIM == 3)
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            gz(i,j,k,n) = 0.5 * dzi*(s(i,j,k+1,n) - s(i,j,k-1,n));
        });
#endif
    }

    // 2. Normalize the gradient if needed
    if (normalize == 1) {
        Real eps = 1.e-10;
        FArrayBox  magfb;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(phi, TilingIfNotGPU());  mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(1); // 1 gt
            magfb.resize(bx,ncomp);
            Elixir eli = magfb.elixir();
            const auto& mag = magfb.array();
            AMREX_D_TERM(const auto& gx = phi_cc_grad[0]->array(mfi);,
                         const auto& gy = phi_cc_grad[1]->array(mfi);,
                         const auto& gz = phi_cc_grad[2]->array(mfi););
#if (AMREX_SPACEDIM == 2)            
            ParallelFor(bx, ncomp, [gx, gy, mag]
            AMREX_GPU_DEVICE (int i, int j, int k, int n ) noexcept
            {
                mag(i,j,k,n) = std::sqrt(gx(i,j,k,n) *  gx(i,j,k,n) + gy(i,j,k,n) *  gy(i,j,k,n));
            });
            ParallelFor(bx, ncomp, [gx, gy, mag, eps]
            AMREX_GPU_DEVICE (int i, int j, int k, int n ) noexcept
            {
                gx(i,j,k,n) = gx(i,j,k,n) / (mag(i,j,k,n) + eps);
                gy(i,j,k,n) = gy(i,j,k,n) / (mag(i,j,k,n) + eps);
            });
#endif
#if (AMREX_SPACEDIM == 3)
            ParallelFor(bx, ncomp, [gx, gy, gz, mag]
            AMREX_GPU_DEVICE (int i, int j, int k, int n ) noexcept
            {
                mag(i,j,k,n) = std::sqrt(gx(i,j,k,n) *  gx(i,j,k,n) + gy(i,j,k,n) *  gy(i,j,k,n) + gz(i,j,k,n) *  gz(i,j,k,n));
            });
            ParallelFor(bx, ncomp, [gx, gy, gz, mag, eps]
            AMREX_GPU_DEVICE (int i, int j, int k, int n ) noexcept
            {
                gx(i,j,k,n) = gx(i,j,k,n) / (mag(i,j,k,n) + eps);
                gy(i,j,k,n) = gy(i,j,k,n) / (mag(i,j,k,n) + eps);
                gz(i,j,k,n) = gz(i,j,k,n) / (mag(i,j,k,n) + eps);
            });
#endif
        }
    }

    // 2. Normalize the gradient based on phi_max
    if (normalize == 2) {

    }

//     const Real* dx    = geom.CellSize();

}

// Assume all gts are filled
// phi_cc_grad: cc variables, at least 1 gt
// phi_cc_div: cc variables, 0 gt, it has been defined outside this function
void cc_grad_to_cc_div(MultiFab& phi_cc_div,
                  Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM>& phi_cc_grad,
                  const Geometry& geom)
{

    Print() << "In the cc_grad_to_cc_div " << std::endl;

    const int ncomp = phi_cc_grad[0]->nComp();
    BL_ASSERT(ncomp==phi_cc_div.nComp());
    const int ngrow = phi_cc_grad[0]->nGrow();
    BL_ASSERT(ngrow>=1);

    // 1. Calculate the divergence
    AMREX_D_TERM(const Real dxi = static_cast<Real>(geom.InvCellSize(0));,
                 const Real dyi = static_cast<Real>(geom.InvCellSize(1));,
                 const Real dzi = static_cast<Real>(geom.InvCellSize(2)););
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi_cc_div, TilingIfNotGPU());  mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox(); // 0 gt
        const auto& s = phi_cc_div.array(mfi);
        AMREX_D_TERM(const auto& gx = phi_cc_grad[0]->array(mfi);,
                     const auto& gy = phi_cc_grad[1]->array(mfi);,
                     const auto& gz = phi_cc_grad[2]->array(mfi););

#if (AMREX_SPACEDIM == 2)            
        ParallelFor(bx, ncomp, [s, gx, gy, dxi, dyi]
        AMREX_GPU_DEVICE (int i, int j, int k, int n ) noexcept
        {
            s(i,j,k,n) = 0.5 * dxi*(gx(i+1,j,k,n) - gx(i-1,j,k,n)) +
                         0.5 * dyi*(gy(i,j+1,k,n) - gy(i,j-1,k,n));
        });
#endif
#if (AMREX_SPACEDIM == 3)
        ParallelFor(bx, ncomp, [s, gx, gy, gz, dxi, dyi, dzi]
        AMREX_GPU_DEVICE (int i, int j, int k, int n ) noexcept
        {
            s(i,j,k,n) = 0.5 * dxi*(gx(i+1,j,k,n) - gx(i-1,j,k,n)) +
                         0.5 * dyi*(gy(i,j+1,k,n) - gy(i,j-1,k,n)) +
                         0.5 * dzi*(gz(i,j,k+1,n) - gz(i,j,k-1,n));
        });
#endif
    }

}

// Assume all gts are filled
// phi: cc variables, at least 2 gts
// phi_cc_lap: cc variables, 1 gt, it has been defined outside this function
void cc_to_cc_lap(MultiFab& phi_cc_lap, MultiFab& phi, const Geometry& geom)
{

    Print() << "In the cc_to_cc_lap " << std::endl;

    const int ncomp = phi.nComp();
    BL_ASSERT(ncomp==phi_cc_lap.nComp());
    const int ngrow = phi.nGrow();
    BL_ASSERT(ngrow>=2);

    // 1. Calculate the divergence
    AMREX_D_TERM(const Real dxi = static_cast<Real>(geom.InvCellSize(0));,
                 const Real dyi = static_cast<Real>(geom.InvCellSize(1));,
                 const Real dzi = static_cast<Real>(geom.InvCellSize(2)););
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi_cc_lap, TilingIfNotGPU());  mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(1); // 1 gt
        const auto& s = phi_cc_lap.array(mfi);
        const auto& s_in = phi.array(mfi);

#if (AMREX_SPACEDIM == 2)            
        ParallelFor(bx, ncomp, [s, s_in, dxi, dyi]
        AMREX_GPU_DEVICE (int i, int j, int k, int n ) noexcept
        {
            s(i,j,k,n) = dxi*dxi*(s_in(i+1,j,k,n) - 2.0 * s_in(i,j,k,n) + s_in(i-1,j,k,n)) +
                         dyi*dyi*(s_in(i,j+1,k,n) - 2.0 * s_in(i,j,k,n) + s_in(i,j-1,k,n));
        });
#endif
#if (AMREX_SPACEDIM == 3)
        ParallelFor(bx, ncomp, [s, s_in, dxi, dyi, dzi]
        AMREX_GPU_DEVICE (int i, int j, int k, int n ) noexcept
        {
            s(i,j,k,n) = dxi*dxi*(s_in(i+1,j,k,n) - 2.0 * s_in(i,j,k,n) + s_in(i-1,j,k,n)) +
                         dyi*dyi*(s_in(i,j+1,k,n) - 2.0 * s_in(i,j,k,n) + s_in(i,j-1,k,n)) +
                         dzi*dzi*(s_in(i,j,k+1,n) - 2.0 * s_in(i,j,k,n) + s_in(i,j,k-1,n));
        });
#endif
    }

}

void cc_to_fc(Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM>& phi_fc, const MultiFab& phi, const Geometry& geom)
{

    Print() << "In the cc_to_fc " << std::endl;

    const BoxArray& ba = phi.boxArray();
    const DistributionMapping& dm = phi.DistributionMap();
    const int ncomp = phi.nComp();
    const int ngrow = phi.nGrow();

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        BoxArray nba = convert(ba,IntVect::TheDimensionVector(idim));
        phi_fc[idim] = std::make_unique<MultiFab>(nba, dm, ncomp, 0);
    }
    average_cellcenter_to_face(GetArrOfPtrs(phi_fc), phi, geom, ncomp);

}