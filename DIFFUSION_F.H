#ifndef _DIFFUSION_F_H
#define _DIFFUSION_F_H  $Id: DIFFUSION_F.H,v 1.4 1997-10-08 20:15:26 car Exp $

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
    void FORT_VISCFLUX (const Real* s_old, const Real* s_new,
                        ARLIM_P(slo), ARLIM_P(shi),
                        const int*  lo, const int*  hi,
                        const Real* xflux, ARLIM_P(xf_lo), ARLIM_P(xf_hi),
                        const Real* yflux, ARLIM_P(yf_lo), ARLIM_P(yf_hi),
                        const Real* xarea, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
                        const Real* yarea, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
                        const Real* dx, const Real * mult, const Real* theta);
    
    void FORT_VISCSYNCFLUX (const Real* ssync, ARLIM_P(slo), ARLIM_P(shi),
                            const int*  lo, const int*  hi,
                              const Real* xflux, ARLIM_P(xflo), ARLIM_P(xfhi),
                            const Real* yflux, ARLIM_P(yflo), ARLIM_P(yfhi),
                            const Real* xarea, ARLIM_P(axlo), ARLIM_P(axhi),
		             const Real* yarea, ARLIM_P(aylo), ARLIM_P(ayhi),
                            const Real* dx, const Real * mult);
    
#endif
#if (BL_SPACEDIM == 3)
    void FORT_VISCFLUX (const Real* s_old, const Real* s_new,
                        ARLIM_P(s_lo), ARLIM_P(s_hi),
                        const int*  lo, const int*  hi,
                        const Real* xflux, ARLIM_P(xf_lo), ARLIM_P(xf_hi),
                        const Real* yflux, ARLIM_P(yf_lo), ARLIM_P(yf_hi),
                        const Real* zflux, ARLIM_P(zf_lo), ARLIM_P(zf_hi),
                        const Real* xarea, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
                        const Real* yarea, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
                        const Real* zarea, ARLIM_P(az_lo), ARLIM_P(az_hi),
                        const Real* dx, const Real * mult, const Real* theta);
    
    void FORT_VISCSYNCFLUX (const Real* ssync, ARLIM_P(slo), ARLIM_P(shi),
                            const int*  lo, const int*  hi,
                            const Real* xflux, ARLIM_P(xflo), ARLIM_P(xfhi),
                            const Real* yflux, ARLIM_P(yflo), ARLIM_P(yfhi),
                            const Real* zflux, ARLIM_P(zflo), ARLIM_P(zfhi),
                            const Real* xarea, ARLIM_P(axlo), ARLIM_P(axhi),
                            const Real* yarea, ARLIM_P(aylo), ARLIM_P(ayhi),
                            const Real* zarea, ARLIM_P(azlo), ARLIM_P(azhi),
                            const Real* dx, const Real * mult);
    
#endif
    void FORT_HOOPSRC (ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       Real *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                       const Real *u, ARLIM_P(u_lo), ARLIM_P(u_hi),
                       const Real *r, const Real* mu);
    
    void FORT_HOOPRHS (Real *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                       Real *u,   ARLIM_P(u_lo), ARLIM_P(u_hi),
                       const Real* r, const Real* b, 
                       const Real *vol, ARLIM_P(vol_lo), ARLIM_P(vol_hi));

    void FORT_SETALPHA (Real *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                        const int* lo, const int* hi,
                        const Real* r, ARLIM_P(r_lo), ARLIM_P(r_hi),
                        const Real* b,
                        const Real *vfab, ARLIM_P(v_lo), ARLIM_P(v_hi),
                        Real *denfab, ARLIM_P(den_lo), ARLIM_P(den_hi),
                        int* usehoop, int* useden);
    
    void FORT_SET_TENSOR_ALPHA (Real *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                        const int* lo, const int* hi,
                        const Real* r, ARLIM_P(r_lo), ARLIM_P(r_hi),
                        const Real* b,
                        const Real *vfab, ARLIM_P(v_lo), ARLIM_P(v_hi),
                        Real *denfab, ARLIM_P(den_lo), ARLIM_P(den_hi),
                        const Real *betax, ARLIM_P(betax_lo), ARLIM_P(betax_hi),
                        const Real *betay, ARLIM_P(betay_lo), ARLIM_P(betay_hi),
                        const int* isrz);
    
    void FORT_TENSOR_HOOPRHS (const int* xvel_comp, 
                       Real *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                       Real *u,   ARLIM_P(u_lo), ARLIM_P(u_hi),
                       const Real* r, const Real* b, 
                       const Real *vol, ARLIM_P(vol_lo), ARLIM_P(vol_hi),
                       const Real *betax, ARLIM_P(betax_lo), ARLIM_P(betax_hi),
                       const Real *betay, ARLIM_P(betay_lo), ARLIM_P(betay_hi));

    void FORT_TENSOR_HOOPSRC (const int* xvel_comp,
                       ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       Real *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                       const Real *u, ARLIM_P(u_lo), ARLIM_P(u_hi),
                       const Real *r,
                       const Real *betax, ARLIM_P(betax_lo), ARLIM_P(betax_hi),
                       const Real *betay, ARLIM_P(betay_lo), ARLIM_P(betay_hi));
    
    void FORT_CCTOEDGE(const int* lo, const int* hi,
        ARLIM_P(blo), ARLIM_P(bhi), const Real* betadat,
        ARLIM_P(betaxlo), ARLIM_P(betaxhi), Real* betax,
        ARLIM_P(betaylo), ARLIM_P(betayhi), Real* betay
#if (BL_SPACEDIM==3)
        ,
        ARLIM_P(betazlo), ARLIM_P(betazhi), Real* betaz
#endif
      );


    void FORT_DIV_MU_SI(const int* lo, const int* hi,
          const Real* dx, const Real* mu,
          ARLIM_P(divunlo), ARLIM_P(divunhi), const Real* divun_dat,
          ARLIM_P(divmusilo), ARLIM_P(divmusihi), const Real* divmusi);
    
    void FORT_DIV_VARMU_SI(const int* lo, const int* hi,
          const Real* dx, 
          ARLIM_P(divunlo), ARLIM_P(divunhi), const Real* divun_dat,
          ARLIM_P(betaxlo), ARLIM_P(betaxhi), const Real* betax,
          ARLIM_P(betaylo), ARLIM_P(betayhi), const Real* betay,
#if (BL_SPACEDIM==3)
          ARLIM_P(betazlo), ARLIM_P(betazhi), const Real* betaz,
#endif
          ARLIM_P(divmusilo), ARLIM_P(divmusihi), const Real* divmusi);
    
}

#endif
#endif

