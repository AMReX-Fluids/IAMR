#ifndef _PROJECTION_F_H_
#define _PROJECTION_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_FILTERP       filterp
#    define FORT_VEL_TO_ACCEL  vel_to_accel
#    define FORT_ACCEL_TO_VEL  accel_to_vel
#    define FORT_PROJ_UPDATE   proj_update
#    define FORT_RADMPY        radmpy
#    define FORT_RADDIV        raddiv
#    define FORT_HGC2N         hgc2n
#    define FORT_HGN2C         hgn2c
#    define FORT_HGPHIBC       hgphibc
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_FILTERP       FILTERP
#    define FORT_VEL_TO_ACCEL  VEL_TO_ACCEL
#    define FORT_ACCEL_TO_VEL  ACCEL_TO_VEL
#    define FORT_PROJ_UPDATE   PROJ_UPDATE
#    define FORT_RADMPY        RADMPY
#    define FORT_RADDIV        RADDIV
#    define FORT_HGC2N         HGC2N
#    define FORT_HGN2C         HGN2C
#    define FORT_HGPHIBC       HGPHIBC
#  else
#    define FORT_FILTERP       filterp_
#    define FORT_VEL_TO_ACCEL  vel_to_accel_
#    define FORT_ACCEL_TO_VEL  accel_to_vel_
#    define FORT_PROJ_UPDATE   proj_update_
#    define FORT_RADMPY        radmpy_
#    define FORT_RADDIV        raddiv_
#    define FORT_HGC2N         hgc2n_
#    define FORT_HGN2C         hgn2c_
#    define FORT_HGPHIBC       hgphibc_
#  endif

#include <ArrayLim.H>

extern "C" {
    
    void FORT_FILTERP(REAL * pnew, REAL * scratch, ARLIM_P(p_lo), ARLIM_P(p_hi),
                      const int* lo, const int* hi, const REAL * dx,
                      const REAL * fac);

    void FORT_VEL_TO_ACCEL (const int* lo, const int* hi,
                            const REAL* unew, ARLIM_P(u_lo), ARLIM_P(u_hi),
                            const REAL* uold, ARLIM_P(u_lo), ARLIM_P(u_hi),
                            const REAL* dt );
    
    void FORT_ACCEL_TO_VEL (const int* lo, const int* hi,
                            const REAL* uold, ARLIM_P(u_lo), ARLIM_P(u_hi),
                            const REAL* dt,
                            const REAL* unew, ARLIM_P(u_lo), ARLIM_P(u_hi) );
    
    void FORT_PROJ_UPDATE (const int* lo,   const int* hi,
                           const int *nvar, const int *ngrow, 
                           const REAL* unew, ARLIM_P(u_lo), ARLIM_P(u_hi),
                           const  REAL* dt,
                           const REAL* uold, ARLIM_P(u_lo), ARLIM_P(u_hi) );
    
    void FORT_RADMPY (REAL*a, ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                      const int* nghost,
                      const REAL* row, const int* nrow,
                      const int* n);
    
    void FORT_RADDIV (REAL*a, ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       const int* nghost,
                      const REAL* row, const int* nrow,
                      const int* n);
    
    void FORT_HGC2N(const int* nghost, 
                    ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                    const REAL* dat, const REAL* rcen,
                    ARLIM_P(rhs_lo), ARLIM_P(rhs_hi),
                    REAL* rhsdat,
                    const int* domlo,
                    const int* domhi, const int& lowfix,
                    const int& hifix, const REAL* hr,
                    const int* extrap_edges, 
                    const int* extrap_corners, const int* is_rz);
    
    void FORT_HGN2C(const int* isrz, const int* rweighted, ARLIM_P(nodelo),ARLIM_P(nodehi), 
                    const REAL* nodedat,
                    ARLIM_P(cclo), ARLIM_P(cchi), const int* lo, const int* hi,
                    REAL* ccdat);

    void FORT_HGPHIBC(ARLIM_P(u_lo), ARLIM_P(u_hi), const REAL* u,
                      ARLIM_P(d_lo), ARLIM_P(d_hi), const REAL* divcor,
                      ARLIM_P(rho_lo), ARLIM_P(rho_hi), const REAL* rho,
                      ARLIM_P(r_lo), ARLIM_P(r_hi), const REAL* rcen,
                      const REAL* hx,
                      ARLIM_P(phi_lo), ARLIM_P(phi_hi), const REAL* phi);
    
};

#endif
#endif

