#ifndef _PROJECTION_F_H_
#define _PROJECTION_F_H_ 

/*
** $Id: PROJECTION_F.H,v 1.18 2004-09-08 20:58:54 almgren Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_FILTERP       filterp
#    define FORT_VEL_TO_ACCEL  veltoaccel
#    define FORT_ACCEL_TO_VEL  acceltovel
#    define FORT_PROJ_UPDATE   projupdate
#    define FORT_RADMPY        radmpy
#    define FORT_RADDIV        raddiv
#    define FORT_ANELCOEFFMPY  anelcoeffmpy
#    define FORT_HGC2N         hgc2n
#    define FORT_HGN2C         hgn2c
#    define FORT_FILTRHS       filtrhs
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_FILTERP       FILTERP
#    define FORT_VEL_TO_ACCEL  VELTOACCEL
#    define FORT_ACCEL_TO_VEL  ACCELTOVEL
#    define FORT_PROJ_UPDATE   PROJUPDATE
#    define FORT_RADMPY        RADMPY
#    define FORT_RADDIV        RADDIV
#    define FORT_ANELCOEFFMPY  ANELCOEFFMPY
#    define FORT_HGC2N         HGC2N
#    define FORT_HGN2C         HGN2C
#    define FORT_FILTRHS       FILTRHS
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_FILTERP       filterp
#    define FORT_VEL_TO_ACCEL  veltoaccel
#    define FORT_ACCEL_TO_VEL  acceltovel
#    define FORT_PROJ_UPDATE   projupdate
#    define FORT_RADMPY        radmpy
#    define FORT_RADDIV        raddiv
#    define FORT_ANELCOEFFMPY  anelcoeffmpy
#    define FORT_HGC2N         hgc2n
#    define FORT_HGN2C         hgn2c
#    define FORT_FILTRHS       filtrhs
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_FILTERP       filterp_
#    define FORT_VEL_TO_ACCEL  veltoaccel_
#    define FORT_ACCEL_TO_VEL  acceltovel_
#    define FORT_PROJ_UPDATE   projupdate_
#    define FORT_RADMPY        radmpy_
#    define FORT_RADDIV        raddiv_
#    define FORT_ANELCOEFFMPY  anelcoeffmpy_
#    define FORT_HGC2N         hgc2n_
#    define FORT_HGN2C         hgn2c_
#    define FORT_FILTRHS       filtrhs_
#  endif

#include <ArrayLim.H>

extern "C" {
    
    void FORT_FILTERP(Real * pnew, Real * scratch, ARLIM_P(p_lo), ARLIM_P(p_hi),
                      const int* lo, const int* hi, const Real * dx,
                      const Real * fac,
#if (BL_SPACEDIM == 2)
                      const int * wrap_around_x, const int * wrap_around_y);
#elif (BL_SPACEDIM == 3)
                      const int * wrap_around_x, const int * wrap_around_y,
                      const int * wrap_around_z);
#endif

    void FORT_VEL_TO_ACCEL (const int* lo, const int* hi,
                            const Real* unew, ARLIM_P(u_lo), ARLIM_P(u_hi),
                            const Real* uold, ARLIM_P(u_lo), ARLIM_P(u_hi),
                            const Real* dt );
    
    void FORT_ACCEL_TO_VEL (const int* lo, const int* hi,
                            const Real* uold, ARLIM_P(u_lo), ARLIM_P(u_hi),
                            const Real* dt,
                            const Real* unew, ARLIM_P(u_lo), ARLIM_P(u_hi) );
    
    void FORT_PROJ_UPDATE (const int* lo,   const int* hi,
                           const int *nvar, const int *ngrow, 
                           const Real* unew, ARLIM_P(u_lo), ARLIM_P(u_hi),
                           const  Real* dt,
                           const Real* uold, ARLIM_P(u_lo), ARLIM_P(u_hi) );
    
    void FORT_RADMPY (Real*a, ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                      const int* domlo, const int* domhi,
                      const int* nghost,
                      const Real* radius, const int* nrad,
                      const Real* bogus_value);
    
    void FORT_RADDIV (Real*a, ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                      const int* domlo, const int* domhi,
                      const int* nghost,
                      const Real* radius, const int* nrad,
                      const Real* bogus_value);

    void FORT_ANELCOEFFMPY (Real*a, ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                            const int* domlo, const int* domhi,
                            const int* nghost,
                            const Real* anelcoeff, const int* nr,
                            const Real* bogus_value, const int* mult);
    
    void FORT_HGC2N(const int* nghost, 
                    ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                    const Real* dat, const Real* rcen,
                    ARLIM_P(rhs_lo), ARLIM_P(rhs_hi),
                    Real* rhsdat,
                    const int* domlo, const int* domhi, 
                    const Real* hr,
                    const int* is_rz, const int* imax);
    
    void FORT_HGN2C(const int* isrz, const int* rweighted, 
                    ARLIM_P(nodelo),ARLIM_P(nodehi), 
                    const Real* nodedat,
                    ARLIM_P(cclo), ARLIM_P(cchi), const int* lo, const int* hi,
                    Real* ccdat);


    void FORT_FILTRHS (const Real* phi , ARLIM_P(p_lo), ARLIM_P(p_hi),
                       const Real* sig , ARLIM_P(s_lo), ARLIM_P(s_hi),
                       Real* divgp, ARLIM_P(dglo), ARLIM_P(dghi),
                       const int* lo, const int* hi,
                       const Real* dx, const Real* mult, int* is_rz);
}

#endif
#endif /*_PROJECTION_F_H_*/
