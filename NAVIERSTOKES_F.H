#ifndef _NAVIERSTOKES_F_H_
#define _NAVIERSTOKES_F_H_

/*
** $Id: NAVIERSTOKES_F.H,v 1.5 1997-12-14 23:35:59 lijewski Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_FILCC       filcc
#    define FORT_GRADP       gradp
#    define FORT_AVGDOWN     avgdown
#    define FORT_PUTDOWN     putdown
#    define FORT_TESTINJECT  testinject
#    define FORT_INCRMULT    incrmult
#    define FORT_TEST_UMAC_PERIODIC    testumacperiodic
#    define FORT_SUMMASS     summass
#    define FORT_SUMMASSSQ   summasssq
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_FILCC       FILCC
#    define FORT_GRADP       GRADP
#    define FORT_AVGDOWN     AVGDOWN
#    define FORT_PUTDOWN     PUTDOWN
#    define FORT_TESTINJECT  TESTINJECT
#    define FORT_INCRMULT    INCRMULT
#    define FORT_TEST_UMAC_PERIODIC    TESTUMACPERIODIC
#    define FORT_SUMMASS     SUMMASS
#    define FORT_SUMMASSSQ   SUMMASSSQ
#  else
#    define FORT_FILCC       filcc_
#    define FORT_GRADP       gradp_
#    define FORT_AVGDOWN     avgdown_
#    define FORT_PUTDOWN     putdown_
#    define FORT_TESTINJECT  testinject_
#    define FORT_INCRMULT    incrmult_
#    define FORT_TEST_UMAC_PERIODIC    testumacperiodic_
#    define FORT_SUMMASS     summass_
#    define FORT_SUMMASSSQ   summasssq_
#  endif

#include <ArrayLim.H>

extern "C" {
    void FORT_FILCC (const Real * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
                     const int * domlo, const int * domhi,
                     const Real * dx_crse, const Real * xlo, 
                     const int * bc);

    void FORT_GRADP (const Real* p, ARLIM_P(p_lo), ARLIM_P(p_hi),
                     const Real* gp, ARLIM_P(gp_lo), ARLIM_P(gp_hi),
                     const int* lo, const int* hi, const Real* dx );
    
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

    void FORT_TEST_UMAC_PERIODIC( const int* lo, const int* hi,
                                  const Real* u, ARLIM_P(u_lo), ARLIM_P(u_hi),
                                  int *xper, Real *udiff,
                                  const Real* v, ARLIM_P(v_lo), ARLIM_P(v_hi),
                                  int *yper, Real *vdiff,
#if (BL_SPACEDIM == 3)                                                   
                                  const Real* w, ARLIM_P(w_lo), ARLIM_P(w_hi),
                                  int *zper, Real *wdiff,
#endif
                                  int *level );
    
    
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
}
    
#endif
#endif /*_NAVIERSTOKES_F_H_*/
