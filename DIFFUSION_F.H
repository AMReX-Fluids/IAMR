#ifndef _DIFFUSION_F_H
#define _DIFFUSION_F_H  $Id: DIFFUSION_F.H,v 1.3 1997-09-29 17:06:24 lijewski Exp $

#ifdef BL_LANG_FORT
#    define FORT_VISCFLUX    viscflux
#    define FORT_VISCSYNCFLUX viscsyncflux
#    define FORT_HOOPSRC     hoopsrc
#    define FORT_HOOPRHS     hooprhs
#    define FORT_TENSOR_HOOPSRC     tensorhoopsrc
#    define FORT_TENSOR_HOOPRHS     tensorhooprhs
#    define FORT_SETALPHA    setalpha
#    define FORT_SET_TENSOR_ALPHA    settensoralpha
#    define FORT_DIV_VARMU_SI    divvarmusi
#    define FORT_DIV_MU_SI    divmuusii
#    define FORT_CCTOEDGE        cctoedge
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_VISCFLUX    VISCFLUX
#    define FORT_VISCSYNCFLUX VISCSYNCFLUX
#    define FORT_HOOPSRC     HOOPSRC
#    define FORT_HOOPRHS     HOOPRHS
#    define FORT_TENSOR_HOOPSRC     TENSORHOOPSRC
#    define FORT_TENSOR_HOOPRHS     TENSORHOOPRHS
#    define FORT_SETALPHA    SETALPHA
#    define FORT_SET_TENSOR_ALPHA    SETTENSORALPHA
#    define FORT_DIV_VARMU_SI    DIVVARMUSI
#    define FORT_DIV_MU_SI    DIVMUUSII
#    define FORT_CCTOEDGE        CCTOEDGE
#  else
#    define FORT_VISCFLUX    viscflux_
#    define FORT_VISCSYNCFLUX viscsyncflux_
#    define FORT_HOOPSRC     hoopsrc_
#    define FORT_HOOPRHS     hooprhs_
#    define FORT_TENSOR_HOOPSRC     tensorhoopsrc_
#    define FORT_TENSOR_HOOPRHS     tensorhooprhs_
#    define FORT_SETALPHA    setalpha_
#    define FORT_SET_TENSOR_ALPHA    settensoralpha_
#    define FORT_DIV_VARMU_SI    divvarmusi_
#    define FORT_DIV_MU_SI    divmuusii_
#    define FORT_CCTOEDGE        cctoedge_
#  endif

#include <ArrayLim.H>

extern "C" {
#if (BL_SPACEDIM == 2)
    void FORT_VISCFLUX (const REAL* s_old, const REAL* s_new,
                        ARLIM_P(slo), ARLIM_P(shi),
                        const int*  lo, const int*  hi,
                        const REAL* xflux, ARLIM_P(xf_lo), ARLIM_P(xf_hi),
                        const REAL* yflux, ARLIM_P(yf_lo), ARLIM_P(yf_hi),
                        const REAL* xarea, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
                        const REAL* yarea, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
                        const REAL* dx, const REAL * mult, const REAL* theta);
    
    void FORT_VISCSYNCFLUX (const REAL* ssync, ARLIM_P(slo), ARLIM_P(shi),
                            const int*  lo, const int*  hi,
                              const REAL* xflux, ARLIM_P(xflo), ARLIM_P(xfhi),
                            const REAL* yflux, ARLIM_P(yflo), ARLIM_P(yfhi),
                            const REAL* xarea, ARLIM_P(axlo), ARLIM_P(axhi),
		             const REAL* yarea, ARLIM_P(aylo), ARLIM_P(ayhi),
                            const REAL* dx, const REAL * mult);
    
#endif
#if (BL_SPACEDIM == 3)
    void FORT_VISCFLUX (const REAL* s_old, const REAL* s_new,
                        ARLIM_P(s_lo), ARLIM_P(s_hi),
                        const int*  lo, const int*  hi,
                        const REAL* xflux, ARLIM_P(xf_lo), ARLIM_P(xf_hi),
                        const REAL* yflux, ARLIM_P(yf_lo), ARLIM_P(yf_hi),
                        const REAL* zflux, ARLIM_P(zf_lo), ARLIM_P(zf_hi),
                        const REAL* xarea, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
                        const REAL* yarea, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
                        const REAL* zarea, ARLIM_P(az_lo), ARLIM_P(az_hi),
                        const REAL* dx, const REAL * mult, const REAL* theta);
    
    void FORT_VISCSYNCFLUX (const REAL* ssync, ARLIM_P(slo), ARLIM_P(shi),
                            const int*  lo, const int*  hi,
                            const REAL* xflux, ARLIM_P(xflo), ARLIM_P(xfhi),
                            const REAL* yflux, ARLIM_P(yflo), ARLIM_P(yfhi),
                            const REAL* zflux, ARLIM_P(zflo), ARLIM_P(zfhi),
                            const REAL* xarea, ARLIM_P(axlo), ARLIM_P(axhi),
                            const REAL* yarea, ARLIM_P(aylo), ARLIM_P(ayhi),
                            const REAL* zarea, ARLIM_P(azlo), ARLIM_P(azhi),
                            const REAL* dx, const REAL * mult);
    
#endif
    void FORT_HOOPSRC (ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       REAL *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                       const REAL *u, ARLIM_P(u_lo), ARLIM_P(u_hi),
                       const REAL *r, const REAL* mu);
    
    void FORT_HOOPRHS (REAL *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                       REAL *u,   ARLIM_P(u_lo), ARLIM_P(u_hi),
                       const REAL* r, const REAL* b, 
                       const REAL *vol, ARLIM_P(vol_lo), ARLIM_P(vol_hi));

    void FORT_SETALPHA (REAL *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                        const int* lo, const int* hi,
                        const REAL* r, ARLIM_P(r_lo), ARLIM_P(r_hi),
                        const REAL* b,
                        const REAL *vfab, ARLIM_P(v_lo), ARLIM_P(v_hi),
                        REAL *denfab, ARLIM_P(den_lo), ARLIM_P(den_hi),
                        int* usehoop, int* useden);
    
    void FORT_SET_TENSOR_ALPHA (REAL *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                        const int* lo, const int* hi,
                        const REAL* r, ARLIM_P(r_lo), ARLIM_P(r_hi),
                        const REAL* b,
                        const REAL *vfab, ARLIM_P(v_lo), ARLIM_P(v_hi),
                        REAL *denfab, ARLIM_P(den_lo), ARLIM_P(den_hi),
                        const REAL *betax, ARLIM_P(betax_lo), ARLIM_P(betax_hi),
                        const REAL *betay, ARLIM_P(betay_lo), ARLIM_P(betay_hi),
                        const int* isrz);
    
    void FORT_TENSOR_HOOPRHS (const int* xvel_comp, 
                       REAL *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                       REAL *u,   ARLIM_P(u_lo), ARLIM_P(u_hi),
                       const REAL* r, const REAL* b, 
                       const REAL *vol, ARLIM_P(vol_lo), ARLIM_P(vol_hi),
                       const REAL *betax, ARLIM_P(betax_lo), ARLIM_P(betax_hi),
                       const REAL *betay, ARLIM_P(betay_lo), ARLIM_P(betay_hi));

    void FORT_TENSOR_HOOPSRC (const int* xvel_comp,
                       ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       REAL *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                       const REAL *u, ARLIM_P(u_lo), ARLIM_P(u_hi),
                       const REAL *r,
                       const REAL *betax, ARLIM_P(betax_lo), ARLIM_P(betax_hi),
                       const REAL *betay, ARLIM_P(betay_lo), ARLIM_P(betay_hi));
    
    void FORT_CCTOEDGE(const int* lo, const int* hi,
        ARLIM_P(blo), ARLIM_P(bhi), const REAL* betadat,
        ARLIM_P(betaxlo), ARLIM_P(betaxhi), REAL* betax,
        ARLIM_P(betaylo), ARLIM_P(betayhi), REAL* betay
#if (BL_SPACEDIM==3)
        ,
        ARLIM_P(betazlo), ARLIM_P(betazhi), REAL* betaz
#endif
      );


    void FORT_DIV_MU_SI(const int* lo, const int* hi,
          const REAL* dx, const REAL* mu,
          ARLIM_P(divunlo), ARLIM_P(divunhi), const REAL* divun_dat,
          ARLIM_P(divmusilo), ARLIM_P(divmusihi), const REAL* divmusi);
    
    void FORT_DIV_VARMU_SI(const int* lo, const int* hi,
          const REAL* dx, 
          ARLIM_P(divunlo), ARLIM_P(divunhi), const REAL* divun_dat,
          ARLIM_P(betaxlo), ARLIM_P(betaxhi), const REAL* betax,
          ARLIM_P(betaylo), ARLIM_P(betayhi), const REAL* betay,
#if (BL_SPACEDIM==3)
          ARLIM_P(betazlo), ARLIM_P(betazhi), const REAL* betaz,
#endif
          ARLIM_P(divmusilo), ARLIM_P(divmusihi), const REAL* divmusi);
    
}

#endif
#endif

