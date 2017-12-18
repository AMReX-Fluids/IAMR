
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

#ifdef __cplusplus
extern "C"
{
#endif
#if (BL_SPACEDIM == 2)
    void FORT_VISCSYNCFLUX (const amrex::Real* ssync, ARLIM_P(slo), ARLIM_P(shi),
                            const int*  xlo, const int*  xhi,
                            const int*  ylo, const int*  yhi,
                            const amrex::Real* xflux, ARLIM_P(xflo), ARLIM_P(xfhi),
                            const amrex::Real* yflux, ARLIM_P(yflo), ARLIM_P(yfhi),
                            const amrex::Real* xarea, ARLIM_P(axlo), ARLIM_P(axhi),
                            const amrex::Real* yarea, ARLIM_P(aylo), ARLIM_P(ayhi),
                            const amrex::Real* dx, const amrex::Real * mult);
    
#endif
#if (BL_SPACEDIM == 3)
    void FORT_VISCSYNCFLUX (const amrex::Real* ssync, ARLIM_P(slo), ARLIM_P(shi),
                            const int*  xlo, const int*  xhi,
                            const int*  ylo, const int*  yhi,
                            const int*  zlo, const int*  zhi,
                            const amrex::Real* xflux, ARLIM_P(xflo), ARLIM_P(xfhi),
                            const amrex::Real* yflux, ARLIM_P(yflo), ARLIM_P(yfhi),
                            const amrex::Real* zflux, ARLIM_P(zflo), ARLIM_P(zfhi),
                            const amrex::Real* xarea, ARLIM_P(axlo), ARLIM_P(axhi),
                            const amrex::Real* yarea, ARLIM_P(aylo), ARLIM_P(ayhi),
                            const amrex::Real* zarea, ARLIM_P(azlo), ARLIM_P(azhi),
                            const amrex::Real* dx, const amrex::Real * mult);
    
#endif
    void FORT_HOOPSRC (ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       amrex::Real *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                       const amrex::Real *u, ARLIM_P(u_lo), ARLIM_P(u_hi),
                       const amrex::Real *r, const amrex::Real* mu);
    
    void FORT_HOOPRHS (ARLIM_P(bx_lo), ARLIM_P(bx_hi),
		       amrex::Real *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                       const amrex::Real *u,   ARLIM_P(u_lo), ARLIM_P(u_hi),
                       const amrex::Real* r, const amrex::Real* b, 
                       const amrex::Real *vol, ARLIM_P(vol_lo), ARLIM_P(vol_hi));

    void FORT_SETALPHA (amrex::Real *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                        const int* lo, const int* hi,
                        const amrex::Real* r, ARLIM_P(r_lo), ARLIM_P(r_hi),
                        const amrex::Real* b,
                        const amrex::Real *vfab, ARLIM_P(v_lo), ARLIM_P(v_hi),
                        const amrex::Real *denfab, ARLIM_P(den_lo), ARLIM_P(den_hi),
                        int* usehoop, int* useden);
    
    void FORT_SET_TENSOR_ALPHA (amrex::Real *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                        const int* lo, const int* hi,
                        const amrex::Real* r, ARLIM_P(r_lo), ARLIM_P(r_hi),
                        const amrex::Real* b,
                        const amrex::Real *vfab, ARLIM_P(v_lo), ARLIM_P(v_hi),
                        const amrex::Real *denfab, ARLIM_P(den_lo), ARLIM_P(den_hi),
                        const amrex::Real *betax, ARLIM_P(betax_lo), ARLIM_P(betax_hi),
                        const amrex::Real *betay, ARLIM_P(betay_lo), ARLIM_P(betay_hi),
#if (BL_SPACEDIM ==3)
                        const amrex::Real *betaz, ARLIM_P(betaz_lo), ARLIM_P(betaz_hi),
#endif
                        const int* isrz);
    
    void FORT_TENSOR_HOOPRHS (const int* xvel_comp, 
		       ARLIM_P(bx_lo), ARLIM_P(bx_hi),
	       	       amrex::Real *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                       const amrex::Real *u,   ARLIM_P(u_lo), ARLIM_P(u_hi),
                       const amrex::Real* r, const amrex::Real* b, 
                       const amrex::Real *vol, ARLIM_P(vol_lo), ARLIM_P(vol_hi),
                       const amrex::Real *betax, ARLIM_P(betax_lo), ARLIM_P(betax_hi),
                       const amrex::Real *betay, ARLIM_P(betay_lo), ARLIM_P(betay_hi));

    void FORT_TENSOR_HOOPSRC (const int* xvel_comp,
                       ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       amrex::Real *fab, ARLIM_P(fab_lo), ARLIM_P(fab_hi),
                       const amrex::Real *u, ARLIM_P(u_lo), ARLIM_P(u_hi),
		       const amrex::Real *r,
                       const amrex::Real *betax, ARLIM_P(betax_lo), ARLIM_P(betax_hi),
                       const amrex::Real *betay, ARLIM_P(betay_lo), ARLIM_P(betay_hi));

    void FORT_DIV_MU_SI(const int* lo, const int* hi,
          const amrex::Real* dx, const amrex::Real* mu,
          ARLIM_P(divunlo), ARLIM_P(divunhi), const amrex::Real* divun_dat,
          ARLIM_P(divmusilo), ARLIM_P(divmusihi), const amrex::Real* divmusi);
    
    void FORT_DIV_VARMU_SI(const int* lo, const int* hi,
          const amrex::Real* dx, 
          ARLIM_P(divunlo), ARLIM_P(divunhi), const amrex::Real* divun_dat,
          ARLIM_P(betaxlo), ARLIM_P(betaxhi), const amrex::Real* betax,
          ARLIM_P(betaylo), ARLIM_P(betayhi), const amrex::Real* betay,
#if (BL_SPACEDIM==3)
          ARLIM_P(betazlo), ARLIM_P(betazhi), const amrex::Real* betaz,
#endif
          ARLIM_P(divmusilo), ARLIM_P(divmusihi), const amrex::Real* divmusi);
    
#ifdef __cplusplus
}
#endif

#endif
#endif /*_DIFFUSION_F_H*/
