#ifndef _NAVIERSTOKES_F_H_
#define _NAVIERSTOKES_F_H_  $Id: NAVIERSTOKES_F.H,v 1.2 1997-07-17 22:01:04 vince Exp $

#ifdef BL_LANG_FORT
#    define FORT_FILCC       filcc
#    define FORT_GRADP       gradp
#    define FORT_AVGDOWN     avgdown
#    define FORT_PUTDOWN     putdown
#    define FORT_TESTINJECT  testinject
#    define FORT_INCRMULT    incrmult
#    define FORT_TEST_UMAC_PERIODIC    test_umac_periodic
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
#    define FORT_TEST_UMAC_PERIODIC    TEST_UMAC_PERIODIC
#    define FORT_SUMMASS     SUMMASS
#    define FORT_SUMMASSSQ   SUMMASSSQ
#  else
#    define FORT_FILCC       filcc_
#    define FORT_GRADP       gradp_
#    define FORT_AVGDOWN     avgdown_
#    define FORT_PUTDOWN     putdown_
#    define FORT_TESTINJECT  testinject_
#    define FORT_INCRMULT    incrmult_
#    define FORT_TEST_UMAC_PERIODIC    test_umac_periodic_
#    define FORT_SUMMASS     summass_
#    define FORT_SUMMASSSQ   summasssq_
#  endif

#include <ArrayLim.H>

extern "C" {
    void FORT_FILCC (const REAL * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
                     const int * domlo, const int * domhi,
                     const REAL * dx_crse, const REAL * xlo, 
                     const int * bc);

    void FORT_GRADP (const REAL* p, ARLIM_P(p_lo), ARLIM_P(p_hi),
                     const REAL* gp, ARLIM_P(gp_lo), ARLIM_P(gp_hi),
                     const int* lo, const int* hi, const REAL* dx );
    
    void FORT_AVGDOWN (const REAL* crse, ARLIM_P(clo),     ARLIM_P(chi),
                       const int* nvar,
                       const REAL* fine, ARLIM_P(fine_lo), ARLIM_P(fine_hi),
                       const REAL* cvol, ARLIM_P(cv_lo),   ARLIM_P(cv_hi), 
                       const REAL* fvol, ARLIM_P(fv_lo),   ARLIM_P(fv_hi),
                       const int* lo, const int* hi, const int* ratio);
    
    void FORT_PUTDOWN (REAL* crse, ARLIM_P(crse_lo), ARLIM_P(crse_hi),
                       const REAL* fine, ARLIM_P(fine_lo), ARLIM_P(fine_hi),
                       const int* lo, const int* hi,
                       const int* ratio);
    
    void FORT_TESTINJECT (REAL* crse, ARLIM_P(crse_lo), ARLIM_P(crse_hi),
                          const REAL* fine, ARLIM_P(fine_lo), ARLIM_P(fine_hi),
                          const int* lo, const int* hi,
                          const int* ratio);
    
    void FORT_INCRMULT (REAL* a, ARLIM_P(a_lo), ARLIM_P(a_hi),
                        const REAL* b, ARLIM_P(b_lo), ARLIM_P(b_hi),
                        const int* lo, const int* hi,
                        const REAL* alpha);

    void FORT_TEST_UMAC_PERIODIC( const int* lo, const int* hi,
                                  const REAL* u, ARLIM_P(u_lo), ARLIM_P(u_hi),
                                  int *xper, REAL *udiff,
                                  const REAL* v, ARLIM_P(v_lo), ARLIM_P(v_hi),
                                  int *yper, REAL *vdiff,
#if (BL_SPACEDIM == 3)                                                   
                                  const REAL* w, ARLIM_P(w_lo), ARLIM_P(w_hi),
                                  int *zper, REAL *wdiff,
#endif
                                  int *level );
    
    
#if (BL_SPACEDIM == 2)
    void FORT_SUMMASS (const REAL* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                       ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       const REAL* delta, REAL* mass,
                       const REAL* radius,
                       const int* irlo, const int* irhi,
                       const int* rz_flag, REAL* tmp);

    void FORT_SUMMASS (const REAL* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                       ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       const REAL* delta, REAL* mass,
                       const REAL* radius,
                       const int* irlo, const int* irhi,
                       const int* rz_flag, REAL* tmp);
#else
    void FORT_SUMMASS (const REAL* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                       ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       const REAL* delta, REAL* mass, REAL* tmp);

    void FORT_SUMMASSSQ (const REAL* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                         ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                         const REAL* delta, REAL* mass, REAL* tmp);
#endif
}
    
#endif
#endif

