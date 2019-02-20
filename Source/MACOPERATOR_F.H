
#ifndef _MACOPERATOR_F_H_
#define _MACOPERATOR_F_H_

#ifndef BL_LANG_FORT

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C"
{
#endif
#if (BL_SPACEDIM == 2)
    void maccoef (amrex::Real* cx, ARLIM_P(cxlo), ARLIM_P(cxhi),
                       amrex::Real* cy, ARLIM_P(cylo), ARLIM_P(cyhi),
                       const amrex::Real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                       const amrex::Real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
                       const amrex::Real* rho, ARLIM_P(rlo), ARLIM_P(rhi),
                       const int* lo, const int* hi, const int* vbxhi,
                       const amrex::Real* dx);

    void macrhs  (const amrex::Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                       const amrex::Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                       const amrex::Real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                       const amrex::Real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
                       const amrex::Real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                       amrex::Real* rhs, ARLIM_P(rlo), ARLIM_P(rhi),
                       const int* lo, const int* hi,
                       const amrex::Real* scale);

    void macupdate (const int *init,
                         amrex::Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                         amrex::Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                         const amrex::Real* phi, ARLIM_P(p_lo), ARLIM_P(p_hi),
                         const amrex::Real* rho, ARLIM_P(rlo), ARLIM_P(rhi),
                         const int* lo, const int* hi, const int* vbxhi,
                         const amrex::Real* dx, const amrex::Real* scale);
#endif

#if (BL_SPACEDIM == 3)
    void maccoef (amrex::Real* cx, ARLIM_P(cxlo), ARLIM_P(cxhi),
                       amrex::Real* cy, ARLIM_P(cylo), ARLIM_P(cyhi),
                       amrex::Real* cz, ARLIM_P(czlo), ARLIM_P(czhi),
                       const amrex::Real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                       const amrex::Real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
                       const amrex::Real* az, ARLIM_P(azlo), ARLIM_P(azhi),
                       const amrex::Real* rho, ARLIM_P(rlo), ARLIM_P(rhi),
                       const int* lo, const int* hi, const int* vbxhi,
                       const amrex::Real* dx);

    void macrhs  (const amrex::Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                       const amrex::Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                       const amrex::Real* uz, ARLIM_P(uzlo), ARLIM_P(uzhi),
                       const amrex::Real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                       const amrex::Real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
                       const amrex::Real* az, ARLIM_P(azlo), ARLIM_P(azhi),
                       const amrex::Real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                       amrex::Real* rhs, ARLIM_P(rlo), ARLIM_P(rhi),
                       const int* lo, const int* hi,
                       const amrex::Real* scale);

    void macupdate (const int *init,
                         amrex::Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                         amrex::Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                         amrex::Real* uz, ARLIM_P(uzlo), ARLIM_P(uzhi),
                         const amrex::Real* phi, ARLIM_P(p_lo), ARLIM_P(p_hi),
                         const amrex::Real* rho, ARLIM_P(rlo), ARLIM_P(rhi),
                         const int* lo, const int* hi, const int* vbxhi,
                         const amrex::Real* dx, const amrex::Real* scale); 

#endif

    void macsyncrhs (amrex::Real* rhs, ARLIM_P(rlo), ARLIM_P(rhi),
                          const int* lo, const int* hi,
                          const amrex::Real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                          const amrex::Real* rhs_scale);
#ifdef __cplusplus
}
#endif
#endif

#endif /*_MACOPERATOR_F_H_*/
