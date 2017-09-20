
#ifndef _PROJECTION_F_H_
#define _PROJECTION_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_VEL_TO_ACCEL  veltoaccel
#    define FORT_ACCEL_TO_VEL  acceltovel
#    define FORT_PROJ_UPDATE   projupdate
#    define FORT_RADMPYSCAL    radmpyscal
#    define FORT_RADMPYVEL     radmpyvel
#    define FORT_RADDIV        raddiv
#    define FORT_ANELCOEFFMPY  anelcoeffmpy
#    define FORT_HGC2N         hgc2n
#    define FORT_HGN2C         hgn2c
#    define FORT_FILTRHS       filtrhs
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_VEL_TO_ACCEL  VELTOACCEL
#    define FORT_ACCEL_TO_VEL  ACCELTOVEL
#    define FORT_PROJ_UPDATE   PROJUPDATE
#    define FORT_RADMPYSCAL    RADMPYSCAL
#    define FORT_RADMPYVEL     RADMPVEL
#    define FORT_RADDIV        RADDIV
#    define FORT_ANELCOEFFMPY  ANELCOEFFMPY
#    define FORT_HGC2N         HGC2N
#    define FORT_HGN2C         HGN2C
#    define FORT_FILTRHS       FILTRHS
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_VEL_TO_ACCEL  veltoaccel
#    define FORT_ACCEL_TO_VEL  acceltovel
#    define FORT_PROJ_UPDATE   projupdate
#    define FORT_RADMPYSCAL    radmpyscal
#    define FORT_RADMPYVEL     radmpyvel
#    define FORT_RADDIV        raddiv
#    define FORT_ANELCOEFFMPY  anelcoeffmpy
#    define FORT_HGC2N         hgc2n
#    define FORT_HGN2C         hgn2c
#    define FORT_FILTRHS       filtrhs
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_VEL_TO_ACCEL  veltoaccel_
#    define FORT_ACCEL_TO_VEL  acceltovel_
#    define FORT_PROJ_UPDATE   projupdate_
#    define FORT_RADMPYSCAL    radmpyscal_
#    define FORT_RADMPYVEL     radmpyvel_
#    define FORT_RADDIV        raddiv_
#    define FORT_ANELCOEFFMPY  anelcoeffmpy_
#    define FORT_HGC2N         hgc2n_
#    define FORT_HGN2C         hgn2c_
#    define FORT_FILTRHS       filtrhs_
#  endif

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C" {
#endif

    void FORT_VEL_TO_ACCEL (const int* lo, const int* hi,
                            const amrex::Real* unew, ARLIM_P(u_lo), ARLIM_P(u_hi),
                            const amrex::Real* uold, ARLIM_P(u_lo), ARLIM_P(u_hi),
                            const amrex::Real* dt );
    
    void FORT_ACCEL_TO_VEL (const int* lo, const int* hi,
                            const amrex::Real* uold, ARLIM_P(u_lo), ARLIM_P(u_hi),
                            const amrex::Real* dt,
                            const amrex::Real* unew, ARLIM_P(u_lo), ARLIM_P(u_hi) );
    
    void FORT_PROJ_UPDATE (const int* lo,   const int* hi,
                           const int *nvar, const int *ngrow, 
                           const amrex::Real* unew, ARLIM_P(u_lo), ARLIM_P(u_hi),
                           const  amrex::Real* dt,
                           const amrex::Real* uold, ARLIM_P(u_lo), ARLIM_P(u_hi) );
    
#if (BL_SPACEDIM < 3)
    void FORT_RADMPYSCAL (amrex::Real*a, ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                          const int* domlo, const int* domhi,
                          const int* nghost,
                          const amrex::Real* radius, const int* nrad,
                          const amrex::Real* bogus_value);

    void FORT_RADMPYVEL (amrex::Real*a, ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                         const int* domlo, const int* domhi,
                         const int* nghost,
                         const amrex::Real* radius, const int* nrad,
                         const int* ndim);
    
    void FORT_RADDIV (amrex::Real*a, ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                      const int* domlo, const int* domhi,
                      const int* nghost,
                      const amrex::Real* radius, const int* nrad,
                      const amrex::Real* bogus_value);
#endif

    void FORT_ANELCOEFFMPY (amrex::Real*a, ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                            const int* domlo, const int* domhi,
                            const int* nghost,
                            const amrex::Real* anelcoeff, const int* nr,
                            const amrex::Real* bogus_value, const int* mult);
    
    void FORT_HGC2N(const int* nghost, 
                    ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                    const amrex::Real* dat, const amrex::Real* rcen,
                    ARLIM_P(rhs_lo), ARLIM_P(rhs_hi),
                    amrex::Real* rhsdat,
                    const int* domlo, const int* domhi, 
                    const amrex::Real* hr,
                    const int* is_rz, const int* imax);
    
    void FORT_HGN2C(const int* isrz, const int* rweighted, 
                    ARLIM_P(nodelo),ARLIM_P(nodehi), 
                    const amrex::Real* nodedat,
                    ARLIM_P(cclo), ARLIM_P(cchi), const int* lo, const int* hi,
                    amrex::Real* ccdat);


    void FORT_FILTRHS (const amrex::Real* phi , ARLIM_P(p_lo), ARLIM_P(p_hi),
                       const amrex::Real* sig , ARLIM_P(s_lo), ARLIM_P(s_hi),
                       amrex::Real* divgp, ARLIM_P(dglo), ARLIM_P(dghi),
                       const int* lo, const int* hi,
                       const amrex::Real* dx, const amrex::Real* mult, int* is_rz);
#ifdef __cplusplus
}
#endif

#endif
#endif /*_PROJECTION_F_H_*/
