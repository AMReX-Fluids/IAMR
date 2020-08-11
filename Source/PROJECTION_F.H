
#ifndef _PROJECTION_F_H_
#define _PROJECTION_F_H_ 

#ifndef BL_LANG_FORT

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C" {
#endif

    void accel_to_vel (const int* lo, const int* hi,
                            const amrex::Real* uold, ARLIM_P(u_lo), ARLIM_P(u_hi),
                            const amrex::Real* dt,
                            const amrex::Real* unew, ARLIM_P(u_lo), ARLIM_P(u_hi) );
    
    void proj_update (const int* lo,   const int* hi,
                           const int *nvar, 
                           const amrex::Real* unew, ARLIM_P(u_lo), ARLIM_P(u_hi),
                           const  amrex::Real* dt,
                           const amrex::Real* uold, ARLIM_P(u_lo), ARLIM_P(u_hi) );
    
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

    void anelcoeffmpy (const int* lo, const int* hi,
		       amrex::Real*a, ARLIM_P(grid_lo), ARLIM_P(grid_hi),
		       const int* domlo, const int* domhi,
		       const amrex::Real* anelcoeff, const int* anel_lo,
		       const int* anel_hi,
		       const amrex::Real* bogus_value, const int* mult);
    
    void hgc2n(const int* nghost, 
                    ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                    const amrex::Real* dat, const amrex::Real* rcen,
                    ARLIM_P(rhs_lo), ARLIM_P(rhs_hi),
                    amrex::Real* rhsdat,
                    const int* domlo, const int* domhi, 
                    const amrex::Real* hr,
                    const int* is_rz, const int* imax);
    
    void hgn2c(const int* isrz, const int* rweighted, 
                    ARLIM_P(nodelo),ARLIM_P(nodehi), 
                    const amrex::Real* nodedat,
                    ARLIM_P(cclo), ARLIM_P(cchi), const int* lo, const int* hi,
                    amrex::Real* ccdat);


    void filtrhs (const amrex::Real* phi , ARLIM_P(p_lo), ARLIM_P(p_hi),
                       const amrex::Real* sig , ARLIM_P(s_lo), ARLIM_P(s_hi),
                       amrex::Real* divgp, ARLIM_P(dglo), ARLIM_P(dghi),
                       const int* lo, const int* hi,
                       const amrex::Real* dx, const amrex::Real* mult, int* is_rz);
#ifdef __cplusplus
}
#endif

#endif
#endif /*_PROJECTION_F_H_*/
