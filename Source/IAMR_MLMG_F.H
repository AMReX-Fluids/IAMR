#ifndef IAMR_MLMG_F_H_
#define IAMR_MLMG_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void iamr_mac_coef (const int* xlo, const int* xhi,
                        const int* ylo, const int* yhi,
#if (AMREX_SPACEDIM == 3)
                        const int* zlo, const int* zhi,
#endif
                        amrex_real* bx, const int* bxlo, const int* bxhi,
                        amrex_real* by, const int* bylo, const int* byhi,
#if (AMREX_SPACEDIM == 3)
                        amrex_real* bz, const int* bzlo, const int* bzhi,
#endif
                        const amrex_real* rho, const int* rlo, const int* rhi,
                        const amrex_real* scale);

    void iamr_mac_rhs (const int* lo, const int* hi,
                       amrex_real* rhs, const int* rlo, const int* rhi,
                       const amrex_real* ux, const int* uxlo, const int* uxhi,
                       const amrex_real* uy, const int* uylo, const int* uyhi,
#if (AMREX_SPACEDIM == 3)
                       const amrex_real* uz, const int* uzlo, const int* uzhi,
#endif
#if (AMREX_SPACEDIM == 2)
                       const amrex_real* ax, const int* axlo, const int* axhi,
                       const amrex_real* ay, const int* aylo, const int* ayhi,
                       const amrex_real* vol, const int* vlo, const int* vhi);
#else
                       const amrex_real* dxinv);
#endif
                       

#ifdef __cplusplus
}
#endif

#endif
