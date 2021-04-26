
#ifndef _PROJECTION_F_H_
#define _PROJECTION_F_H_ 

#ifndef BL_LANG_FORT

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C" {
#endif
    
#if (BL_SPACEDIM < 3)
    void radmpyscal (const int* lo, const int* hi,
		     amrex::Real*a, ARLIM_P(alo), ARLIM_P(ahi),
		     const int* domlo, const int* domhi,
		     const amrex::Real* radius, const int* rlo,
		     const int* rhi);

    void radmpyvel (const int* lo, const int* hi,
		    amrex::Real*a, ARLIM_P(alo), ARLIM_P(ahi),
		    const int* domlo, const int* domhi,
		    const amrex::Real* radius, const int* rlo,
		    const int* rhi, const int* ndim);
    
    void fort_raddiv (const int* lo, const int* hi,
		      amrex::Real*a, ARLIM_P(alo), ARLIM_P(ahi),
		      const int* domlo, const int* domhi,
		      const amrex::Real* radius, const int* rlo,
		      const int* rhi, const amrex::Real* bogus_value);
#endif

#ifdef __cplusplus
}
#endif

#endif
#endif /*_PROJECTION_F_H_*/
