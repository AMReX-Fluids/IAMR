
#ifndef _DIFFUSION_F_H
#define _DIFFUSION_F_H

#ifdef BL_LANG_FORT
#    define FORT_VISCSYNCFLUX     viscsyncflux
#    define FORT_HOOPSRC          hoopsrc
#    define FORT_HOOPRHS          hooprhs
#    define FORT_TENSOR_HOOPSRC   tensorhoopsrc
#    define FORT_TENSOR_HOOPRHS   tensorhooprhs
#    define FORT_SETALPHA         setalpha
#    define FORT_SET_TENSOR_ALPHA settensoralpha
#    define FORT_DIV_VARMU_SI     divvarmusi
#    define FORT_DIV_MU_SI        divmuusii
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_VISCSYNCFLUX     VISCSYNCFLUX
#    define FORT_HOOPSRC          HOOPSRC
#    define FORT_HOOPRHS          HOOPRHS
#    define FORT_TENSOR_HOOPSRC   TENSORHOOPSRC
#    define FORT_TENSOR_HOOPRHS   TENSORHOOPRHS
#    define FORT_SETALPHA         SETALPHA
#    define FORT_SET_TENSOR_ALPHA SETTENSORALPHA
#    define FORT_DIV_VARMU_SI     DIVVARMUSI
#    define FORT_DIV_MU_SI        DIVMUUSII
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_VISCSYNCFLUX     viscsyncflux
#    define FORT_HOOPSRC          hoopsrc
#    define FORT_HOOPRHS          hooprhs
#    define FORT_TENSOR_HOOPSRC   tensorhoopsrc
#    define FORT_TENSOR_HOOPRHS   tensorhooprhs
#    define FORT_SETALPHA         setalpha
#    define FORT_SET_TENSOR_ALPHA settensoralpha
#    define FORT_DIV_VARMU_SI     divvarmusi
#    define FORT_DIV_MU_SI        divmuusii
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_VISCSYNCFLUX     viscsyncflux_
#    define FORT_HOOPSRC          hoopsrc_
#    define FORT_HOOPRHS          hooprhs_
#    define FORT_TENSOR_HOOPSRC   tensorhoopsrc_
#    define FORT_TENSOR_HOOPRHS   tensorhooprhs_
#    define FORT_SETALPHA         setalpha_
#    define FORT_SET_TENSOR_ALPHA settensoralpha_
#    define FORT_DIV_VARMU_SI     divvarmusi_
#    define FORT_DIV_MU_SI        divmuusii_
#  endif

#include <AMReX_ArrayLim.H>

extern "C"
{
#if (BL_SPACEDIM == 2)
    void FORT_VISCSYNCFLUX (const Real* ssync, ARLIM_P(slo), ARLIM_P(shi),
                            const int*  xlo, const int*  xhi,
                            const int*  ylo, const int*  yhi,
                            const Real* xflux, ARLIM_P(xflo), ARLIM_P(xfhi),
                            const Real* yflux, ARLIM_P(yflo), ARLIM_P(yfhi),
                            const Real* xarea, ARLIM_P(axlo), ARLIM_P(axhi),
                            const Real* yarea, ARLIM_P(aylo), ARLIM_P(ayhi),
                            const Real* dx, const Real * mult);
    
#endif
#if (BL_SPACEDIM == 3)
    void FORT_VISCSYNCFLUX (const Real* ssync, ARLIM_P(slo), ARLIM_P(shi),
                            const int*  xlo, const int*  xhi,
                            const int*  ylo, const int*  yhi,
                            const int*  zlo, const int*  zhi,
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
                       const Real *u,   ARLIM_P(u_lo), ARLIM_P(u_hi),
                       const Real* r, const Real* b, 
                       const Real *vol, ARLIM_P(vol_lo), ARLIM_P(vol_hi));

    void FORT_SETALPHA (Real *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                        const int* lo, const int* hi,
                        const Real* r, ARLIM_P(r_lo), ARLIM_P(r_hi),
                        const Real* b,
                        const Real *vfab, ARLIM_P(v_lo), ARLIM_P(v_hi),
                        const Real *denfab, ARLIM_P(den_lo), ARLIM_P(den_hi),
                        int* usehoop, int* useden);
    
    void FORT_SET_TENSOR_ALPHA (Real *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                        const int* lo, const int* hi,
                        const Real* r, ARLIM_P(r_lo), ARLIM_P(r_hi),
                        const Real* b,
                        const Real *vfab, ARLIM_P(v_lo), ARLIM_P(v_hi),
                        const Real *denfab, ARLIM_P(den_lo), ARLIM_P(den_hi),
                        const Real *betax, ARLIM_P(betax_lo), ARLIM_P(betax_hi),
                        const Real *betay, ARLIM_P(betay_lo), ARLIM_P(betay_hi),
#if (BL_SPACEDIM ==3)
                        const Real *betaz, ARLIM_P(betaz_lo), ARLIM_P(betaz_hi),
#endif
                        const int* isrz);
    
    void FORT_TENSOR_HOOPRHS (const int* xvel_comp, 
                       Real *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                       const Real *u,   ARLIM_P(u_lo), ARLIM_P(u_hi),
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
#endif /*_DIFFUSION_F_H*/
