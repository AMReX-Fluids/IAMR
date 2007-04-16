#ifndef _NAVIERSTOKES_F_H_
#define _NAVIERSTOKES_F_H_

/*
** $Id: NAVIERSTOKES_F.H,v 1.15 2007-04-16 18:14:53 aaspden Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_FILCC       filcc
#    define FORT_GRADP       gradp
#    define FORT_AVGDOWN     avgdown
#    define FORT_PUTDOWN     putdown
#    define FORT_TESTINJECT  testinject
#    define FORT_INCRMULT    incrmult
#    define FORT_SUMTURB     sumturb
#    define FORT_SUMMASS     summass
#    define FORT_SUMMASSSQ   summasssq
#    define FORT_CEN2EDG     cen2edg
#    define FORT_EDGE_INTERP edge_interp
#    define FORT_PC_CF_EDGE_INTERP pc_edge_interp
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_FILCC       FILCC
#    define FORT_GRADP       GRADP
#    define FORT_AVGDOWN     AVGDOWN
#    define FORT_PUTDOWN     PUTDOWN
#    define FORT_TESTINJECT  TESTINJECT
#    define FORT_INCRMULT    INCRMULT
#    define FORT_SUMTURB     SUMTURB
#    define FORT_SUMMASS     SUMMASS
#    define FORT_SUMMASSSQ   SUMMASSSQ
#    define FORT_CEN2EDG     CEN2EDG
#    define FORT_EDGE_INTERP EDGE_INTERP
#    define FORT_PC_CF_EDGE_INTERP PC_EDGE_INTERP
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_FILCC       filcc
#    define FORT_GRADP       gradp
#    define FORT_AVGDOWN     avgdown
#    define FORT_PUTDOWN     putdown
#    define FORT_TESTINJECT  testinject
#    define FORT_INCRMULT    incrmult
#    define FORT_SUMTURB     sumturb
#    define FORT_SUMMASS     summass
#    define FORT_SUMMASSSQ   summasssq
#    define FORT_CEN2EDG     cen2edg
#    define FORT_EDGE_INTERP edge_interp
#    define FORT_PC_CF_EDGE_INTERP pc_edge_interp
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_FILCC       filcc_
#    define FORT_GRADP       gradp_
#    define FORT_AVGDOWN     avgdown_
#    define FORT_PUTDOWN     putdown_
#    define FORT_TESTINJECT  testinject_
#    define FORT_INCRMULT    incrmult_
#    define FORT_SUMTURB     sumturb_
#    define FORT_SUMMASS     summass_
#    define FORT_SUMMASSSQ   summasssq_
#    define FORT_CEN2EDG     cen2edg_
#    define FORT_EDGE_INTERP edge_interp_
#    define FORT_PC_CF_EDGE_INTERP pc_edge_interp_
#  endif

#include <ArrayLim.H>

extern "C" {
    void FORT_FILCC (const Real * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
                     const int * domlo, const int * domhi,
                     const Real * dx_crse, const Real * xlo, 
                     const int * bc);
#if (BL_SPACEDIM == 2)
    void FORT_GRADP (const Real* p, ARLIM_P(p_lo), ARLIM_P(p_hi),
                     const Real* gp, ARLIM_P(gp_lo), ARLIM_P(gp_hi),
                     const int* lo, const int* hi, const Real* dx ,
                     const int* is_full);
#else
    void FORT_GRADP (const Real* p, ARLIM_P(p_lo), ARLIM_P(p_hi),
                     const Real* gp, ARLIM_P(gp_lo), ARLIM_P(gp_hi),
                     const int* lo, const int* hi, const Real* dx);
#endif
    
    void FORT_AVGDOWN (const Real* crse, ARLIM_P(clo),     ARLIM_P(chi),
                       const int* nvar,
                       const Real* fine, ARLIM_P(fine_lo), ARLIM_P(fine_hi),
                       const Real* cvol, ARLIM_P(cv_lo),   ARLIM_P(cv_hi), 
                       const Real* fvol, ARLIM_P(fv_lo),   ARLIM_P(fv_hi),
                       const int* lo, const int* hi, const int* ratio);
    
    void FORT_PUTDOWN (Real* crse, ARLIM_P(crse_lo), ARLIM_P(crse_hi),
                       const Real* fine, ARLIM_P(fine_lo), ARLIM_P(fine_hi),
                       const int* lo, const int* hi,
                       const int* ratio);
    
    void FORT_TESTINJECT (Real* crse, ARLIM_P(crse_lo), ARLIM_P(crse_hi),
                          const Real* fine, ARLIM_P(fine_lo), ARLIM_P(fine_hi),
                          const int* lo, const int* hi,
                          const int* ratio);
    
    void FORT_INCRMULT (Real* a, ARLIM_P(a_lo), ARLIM_P(a_hi),
                        const Real* b, ARLIM_P(b_lo), ARLIM_P(b_hi),
                        const int* lo, const int* hi,
                        const Real* alpha);

    void FORT_SUMTURB (const Real* dat, const Real* pres,
		       ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                       ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       const Real* delta, Real* turb,
		       int *ksize, int *turbVars);

#if (BL_SPACEDIM == 2)
    void FORT_SUMMASS (const Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                       ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       const Real* delta, Real* mass,
                       const Real* radius,
                       const int* irlo, const int* irhi,
                       const int* rz_flag, Real* tmp);

    void FORT_SUMMASS (const Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                       ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       const Real* delta, Real* mass,
                       const Real* radius,
                       const int* irlo, const int* irhi,
                       const int* rz_flag, Real* tmp);
#else
    void FORT_SUMMASS (const Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                       ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       const Real* delta, Real* mass, Real* tmp);

    void FORT_SUMMASSSQ (const Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                         ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                         const Real* delta, Real* mass, Real* tmp);
#endif

    void FORT_CEN2EDG(const int* lo, const int* hi,
                      ARLIM_P(clo), ARLIM_P(chi), const Real* cfabdat,
                      ARLIM_P(elo), ARLIM_P(ehi), Real* efabdat,
                      const int* nc, const int* dir, const int* isharm);

    void FORT_EDGE_INTERP(const int* flo, const int* fhi,
                          const int* nc, const int* refRatio, const int* dir,
                          Real*       fine, ARLIM_P(flo), ARLIM_P(fhi));
    
    void FORT_PC_CF_EDGE_INTERP(const int* lo, const int* hi,
                                const int* nc, const int* refRatio, const int* dir,
                                const Real* crse, ARLIM_P(clo), ARLIM_P(chi),
                                Real*       fine, ARLIM_P(flo), ARLIM_P(fhi));
}
    
#endif
#endif /*_NAVIERSTOKES_F_H_*/
