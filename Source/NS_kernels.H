
#ifndef NS_KERNELS_H
#define NS_KERNELS_H

#include <AMReX_MultiFab.H>

using namespace amrex;

AMREX_INLINE AMREX_GPU_DEVICE
Real minmod(Real alpha_1, Real beta_1)
{
    if (alpha_1 * beta_1 <= 0.0) {
        return 0.0;
    }
    else if (fabs(alpha_1) <= fabs(beta_1)) {
        return alpha_1;
    }
    else {
        return beta_1;
    }
}

void cc_to_cc_grad(Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM>& phi_cc_grad, const MultiFab& phi, const Geometry& geom, int normalize);
void cc_grad_to_cc_div(MultiFab& phi_cc_div, Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM>& phi_cc_grad, const Geometry& geom);
void cc_to_cc_lap(MultiFab& phi_cc_lap, MultiFab& phi, const Geometry& geom);
void cc_to_fc(Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM>& phi_fc, const MultiFab& phi, const Geometry& geom);
#endif
