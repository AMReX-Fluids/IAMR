
#ifndef _NAVIERSTOKES_F_H_
#define _NAVIERSTOKES_F_H_

#ifdef BL_LANG_FORT
#    define FORT_FILCC       filcc
#    define FORT_GRADP       gradp
#    define FORT_PUTDOWN     putdown
#    define FORT_INCRMULT    incrmult
#    define FORT_MAXVAL      maxval
#    define FORT_SUMMASS     summass
#    define FORT_SUMMASS_CYL summass_cyl
#    define FORT_CEN2EDG     cen2edg
#    define FORT_EDGE_INTERP edge_interp
#    define FORT_PC_CF_EDGE_INTERP pc_edge_interp
#    define FORT_BINGHAM 	 bingham
#    define FORT_SET_NS_PARAMS  setnsparams
#if (BL_SPACEDIM == 3)
#    define FORT_SUMTURB     sumturb
#ifdef SUMJET
#    define FORT_SUMJET      sumjet
#endif
#endif
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_FILCC       FILCC
#    define FORT_GRADP       GRADP
#    define FORT_PUTDOWN     PUTDOWN
#    define FORT_INCRMULT    INCRMULT
#    define FORT_MAXVAL      MAXVAL
#    define FORT_SUMMASS     SUMMASS
#    define FORT_SUMMASS_CYL SUMMASS_CYL
#    define FORT_CEN2EDG     CEN2EDG
#    define FORT_EDGE_INTERP EDGE_INTERP
#    define FORT_PC_CF_EDGE_INTERP PC_EDGE_INTERP
#    define FORT_HOEXTRAPTOCC HOEXTRAPTOCC
#    define FORT_BINGHAM 	 BINGHAM
#    define FORT_SET_NS_PARAMS  SETNSPARAMS
#if (BL_SPACEDIM == 3)
#    define FORT_SUMTURB     SUMTURB
#ifdef SUMJET
#    define FORT_SUMJET      SUMJET
#endif
#endif
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_FILCC       filcc
#    define FORT_GRADP       gradp
#    define FORT_PUTDOWN     putdown
#    define FORT_INCRMULT    incrmult
#    define FORT_MAXVAL      maxval
#    define FORT_SUMMASS     summass
#    define FORT_SUMMASS_CYL summass_cyl
#    define FORT_CEN2EDG     cen2edg
#    define FORT_EDGE_INTERP edge_interp
#    define FORT_PC_CF_EDGE_INTERP pc_edge_interp
#    define FORT_HOEXTRAPTOCC hoextraptocc
#    define FORT_BINGHAM 	 bingham
#    define FORT_SET_NS_PARAMS  setnsparams
#if (BL_SPACEDIM == 3)
#    define FORT_SUMTURB     sumturb
#ifdef SUMJET
#    define FORT_SUMJET      sumjet
#endif
#endif
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_FILCC       filcc_
#    define FORT_GRADP       gradp_
#    define FORT_PUTDOWN     putdown_
#    define FORT_INCRMULT    incrmult_
#    define FORT_MAXVAL      maxval_
#    define FORT_SUMMASS     summass_
#    define FORT_SUMMASS_CYL summass_cyl_
#    define FORT_CEN2EDG     cen2edg_
#    define FORT_EDGE_INTERP edge_interp_
#    define FORT_PC_CF_EDGE_INTERP pc_edge_interp_
#    define FORT_HOEXTRAPTOCC hoextraptocc_
#    define FORT_BINGHAM 	 bingham_
#    define FORT_SET_NS_PARAMS  setnsparams_
#if (BL_SPACEDIM == 3)
#    define FORT_SUMTURB     sumturb_
#ifdef SUMJET
#    define FORT_SUMJET      sumjet_
#endif
#endif
#  endif

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C" {
#endif
    void FORT_FILCC (const amrex::Real * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
                     const int * domlo, const int * domhi,
                     const amrex::Real * dx_crse, const amrex::Real * xlo, 
                     const int * bc);
#if (BL_SPACEDIM == 2)
    void FORT_GRADP (const amrex::Real* p, ARLIM_P(p_lo), ARLIM_P(p_hi),
                     const amrex::Real* gp, ARLIM_P(gp_lo), ARLIM_P(gp_hi),
                     const int* lo, const int* hi, const amrex::Real* dx ,
                     const int* is_full);
#else
    void FORT_GRADP (const amrex::Real* p, ARLIM_P(p_lo), ARLIM_P(p_hi),
                     const amrex::Real* gp, ARLIM_P(gp_lo), ARLIM_P(gp_hi),
                     const int* lo, const int* hi, const amrex::Real* dx);
#endif
    
    void FORT_PUTDOWN (amrex::Real* crse, ARLIM_P(crse_lo), ARLIM_P(crse_hi),
                       const amrex::Real* fine, ARLIM_P(fine_lo), ARLIM_P(fine_hi),
                       const int* lo, const int* hi,
                       const int* ratio);
    
    void FORT_INCRMULT (amrex::Real* a, ARLIM_P(a_lo), ARLIM_P(a_hi),
                        const amrex::Real* b, ARLIM_P(b_lo), ARLIM_P(b_hi),
                        const int* lo, const int* hi,
                        const amrex::Real* alpha);

#if (BL_SPACEDIM == 3)
    void FORT_SUMTURB (const amrex::Real* dat, const amrex::Real* pres,
		       ARLIM_P(dat_lo), ARLIM_P(dat_hi),
		       ARLIM_P(pres_lo), ARLIM_P(pres_hi),
                       ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       const amrex::Real* delta, amrex::Real* turb,
		       int *ksize, int *turbVars);

#if SUMJET
    void FORT_SUMJET (const amrex::Real* dat, const amrex::Real* pres,
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
    void FORT_SUMMASS (const amrex::Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                       ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       const amrex::Real* delta, amrex::Real* mass,
                       const amrex::Real* radius,
                       const int* irlo, const int* irhi,
                       const int* rz_flag);

    void FORT_SUMMASS_CYL (const amrex::Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                           const amrex::Real* delta, amrex::Real* mass,
                           const amrex::Real* radius,
                           const int* irlo, const int* irhi,
                           const int* rz_flag, const amrex::Real* plo,
                           const amrex::Real* vws_dz, const amrex::Real* vws_Rcyl);
#else
    void FORT_SUMMASS (const amrex::Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                       ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       const amrex::Real* delta, amrex::Real* mass);

    void FORT_SUMMASS_CYL (const amrex::Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                           const amrex::Real* delta, const amrex::Real* plo, const amrex::Real* vws_dz,
                           const amrex::Real* vws_Rcyl, amrex::Real* mass);
#endif

    void FORT_MAXVAL (const amrex::Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                      ARLIM_P(grid_lo), ARLIM_P(grid_hi), amrex::Real* maxval);


    void FORT_CEN2EDG(const int* lo, const int* hi,
                      ARLIM_P(clo), ARLIM_P(chi), const amrex::Real* cfabdat,
                      ARLIM_P(elo), ARLIM_P(ehi), amrex::Real* efabdat,
                      const int* nc, const int* dir, const int* isharm);

    void FORT_EDGE_INTERP(const int* flo, const int* fhi,
                          const int* nc, const int* refRatio, const int* dir,
                          amrex::Real*       fine, ARLIM_P(flo), ARLIM_P(fhi));
    
    void FORT_PC_CF_EDGE_INTERP(const int* lo, const int* hi,
                                const int* nc, const int* refRatio, const int* dir,
                                const amrex::Real* crse, ARLIM_P(clo), ARLIM_P(chi),
                                amrex::Real*       fine, ARLIM_P(flo), ARLIM_P(fhi));

    void FORT_HOEXTRAPTOCC (const amrex::Real* q, ARLIM_P(q_lo), ARLIM_P(q_hi),
                     const int* dom_lo, const int* dom_hi,
                     const amrex::Real* dx ,const amrex::Real* xlo);

    void FORT_BINGHAM(amrex::Real* visc, ARLIM_P(visc_lo), ARLIM_P(visc_hi),
		              const amrex::Real* vel, ARLIM_P(vel_lo), ARLIM_P(vel_hi),
                      const int* lo, const int* hi, const int* domain_lo, const int* domain_hi,
		              const amrex::Real* delta, const int* vel_bc);

   void FORT_SET_NS_PARAMS(const amrex::Real& vel_visc_coef,
                        const amrex::Real& yield_stress, 
                        const amrex::Real& reg_param);
#ifdef __cplusplus
}
#endif
    
#endif
#endif /*_NAVIERSTOKES_F_H_*/
