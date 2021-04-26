#ifndef _NAVIERSTOKES_F_H_
#define _NAVIERSTOKES_F_H_

#ifndef BL_LANG_FORT

#include <AMReX_ArrayLim.H>
#include <AMReX_filcc_f.H>

#ifdef __cplusplus
extern "C" {
#endif


    void cen2edg(const int* lo, const int* hi, 
                      ARLIM_P(clo), ARLIM_P(chi), const amrex::Real* cfabdat,
                      ARLIM_P(elo), ARLIM_P(ehi), amrex::Real* efabdat,
                      const int* nc, const int* dir, const int* isharm);

    void fort_maxval (const amrex::Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                      ARLIM_P(grid_lo), ARLIM_P(grid_hi), amrex::Real* maxval);

#ifdef __cplusplus
}
#endif
    
#endif
#endif /*_NAVIERSTOKES_F_H_*/
