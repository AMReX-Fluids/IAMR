
#ifndef _NAVIERSTOKESBASE_F_H_
#define _NAVIERSTOKESBASE_F_H_

#ifndef BL_LANG_FORT

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef AMREX_USE_EB
void fort_set_body_state(const int* lo, const int* hi,
                    amrex::Real* S, const int* Slo, const int* Shi,
                    const int* mask, const int* mlo, const int* mhi,
                    const amrex::Real* b, const int* nc, const int* bval);
#endif

#ifdef __cplusplus
}
#endif

#endif
#endif /*_NAVIERSTOKESBASE_F_H_*/