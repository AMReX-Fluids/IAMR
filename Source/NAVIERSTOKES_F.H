
#ifndef _NAVIERSTOKES_F_H_
#define _NAVIERSTOKES_F_H_

#ifndef BL_LANG_FORT

#include <AMReX_ArrayLim.H>
#include <AMReX_filcc_f.H>

#ifdef __cplusplus
extern "C" {
#endif
    void filcc (const amrex::Real * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
                     const int * domlo, const int * domhi,
                     const amrex::Real * dx_crse, const amrex::Real * xlo, 
                     const int * bc);
    void filcc_tile (ARLIM_P(lo), ARLIM_P(hi),
                          const amrex::Real * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
                          const int * domlo, const int * domhi,
                          const amrex::Real * dx_crse, const amrex::Real * xlo, 
                          const int * bc);
#if (BL_SPACEDIM == 2)
    void gradp (const amrex::Real* p, ARLIM_P(p_lo), ARLIM_P(p_hi),
                     const amrex::Real* gp, ARLIM_P(gp_lo), ARLIM_P(gp_hi),
                     const int* lo, const int* hi, const amrex::Real* dx ,
                     const int* is_full);
#else
    void gradp (const amrex::Real* p, ARLIM_P(p_lo), ARLIM_P(p_hi),
                     const amrex::Real* gp, ARLIM_P(gp_lo), ARLIM_P(gp_hi),
                     const int* lo, const int* hi, const amrex::Real* dx);
#endif
    
    void fort_putdown (amrex::Real* crse, ARLIM_P(crse_lo), ARLIM_P(crse_hi),
                       const amrex::Real* fine, ARLIM_P(fine_lo), ARLIM_P(fine_hi),
                       const int* lo, const int* hi,
                       const int* ratio);
    
    void incrmult (amrex::Real* a, ARLIM_P(a_lo), ARLIM_P(a_hi),
                        const amrex::Real* b, ARLIM_P(b_lo), ARLIM_P(b_hi),
                        const int* lo, const int* hi,
                        const amrex::Real* alpha);

#if (BL_SPACEDIM == 3)
    void sumturb (const amrex::Real* dat, const amrex::Real* pres,
		       ARLIM_P(dat_lo), ARLIM_P(dat_hi),
		       ARLIM_P(pres_lo), ARLIM_P(pres_hi),
                       ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       const amrex::Real* delta, amrex::Real* turb,
		       int *ksize, int *turbVars);

#if SUMJET
    void sumjet (const amrex::Real* dat, const amrex::Real* pres,
		      ARLIM_P(dat_lo), ARLIM_P(dat_hi),
		      ARLIM_P(pres_lo), ARLIM_P(pres_hi),
                      ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                      const amrex::Real* delta, amrex::Real* turb,
		      int *levRsize, int *levKsize, int *rsize, int *ksize,
		      int *jetVars, int* numSplit,
		      const amrex::Real* xlo, const amrex::Real* xhi);
#endif
#endif

#if (BL_SPACEDIM == 2)
    void summass (const amrex::Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
		  ARLIM_P(grid_lo), ARLIM_P(grid_hi),
		  const amrex::Real* delta, amrex::Real* mass,
		  const amrex::Real* radius,
		  const int* irlo, const int* irhi,
		  const int* rz_flag);
  
    void summass_eb (const amrex::Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
		     ARLIM_P(grid_lo), ARLIM_P(grid_hi),
		     const amrex::Real* vf, ARLIM_P(vf_lo), ARLIM_P(vf_hi),
		     const amrex::Real* delta,
		     amrex::Real* mass,
		     const amrex::Real* radius,
		     const int* irlo, const int* irhi,
		     const int* rz_flag);

  
    void summass_cyl (const amrex::Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
		      ARLIM_P(grid_lo), ARLIM_P(grid_hi),
		      const amrex::Real* delta, amrex::Real* mass,
		      const amrex::Real* radius,
		      const int* irlo, const int* irhi,
		      const int* rz_flag, const amrex::Real* plo,
		      const amrex::Real* vws_dz, const amrex::Real* vws_Rcyl);
#else
    void summass (const amrex::Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
		  ARLIM_P(grid_lo), ARLIM_P(grid_hi),
		  const amrex::Real* delta, amrex::Real* mass);

    void summass_eb (const amrex::Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
		     ARLIM_P(grid_lo), ARLIM_P(grid_hi),
		     const amrex::Real* vf, ARLIM_P(vf_lo), ARLIM_P(vf_hi),
		     const amrex::Real* delta,
		     amrex::Real* mass);
		  
    void summass_cyl (const amrex::Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
		      ARLIM_P(grid_lo), ARLIM_P(grid_hi),
		      const amrex::Real* delta, const amrex::Real* plo, const amrex::Real* vws_dz,
		      const amrex::Real* vws_Rcyl, amrex::Real* mass);
#endif

    void fort_maxval (const amrex::Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                      ARLIM_P(grid_lo), ARLIM_P(grid_hi), amrex::Real* maxval);


    void cen2edg(const int* lo, const int* hi,
                      ARLIM_P(clo), ARLIM_P(chi), const amrex::Real* cfabdat,
                      ARLIM_P(elo), ARLIM_P(ehi), amrex::Real* efabdat,
                      const int* nc, const int* dir, const int* isharm);

    void edge_interp(const int* flo, const int* fhi,
                          const int* nc, const int* refRatio, const int* dir,
                          amrex::Real*       fine, ARLIM_P(flo), ARLIM_P(fhi));
    
    void pc_edge_interp(const int* lo, const int* hi,
                                const int* nc, const int* refRatio, const int* dir,
                                const amrex::Real* crse, ARLIM_P(clo), ARLIM_P(chi),
                                amrex::Real*       fine, ARLIM_P(flo), ARLIM_P(fhi));



#ifdef __cplusplus
}
#endif
    
#endif
#endif /*_NAVIERSTOKES_F_H_*/
