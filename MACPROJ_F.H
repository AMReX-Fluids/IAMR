#ifndef _MACPROJ_F_H
#define _MACPROJ_F_H

/*
** $Id: MACPROJ_F.H,v 1.10 2004-09-08 21:02:25 almgren Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_MACDIV      macdiv
#    define FORT_SCALEAREA   scalearea
#    define FORT_UCORUPDT    ucorupdt
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_MACDIV      MACDIV
#    define FORT_SCALEAREA   SCALEAREA
#    define FORT_UCORUPDT    UCORUPDT
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_MACDIV      macdiv
#    define FORT_SCALEAREA   scalearea
#    define FORT_UCORUPDT    ucorupdt
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_MACDIV      macdiv_
#    define FORT_SCALEAREA   scalearea_
#    define FORT_UCORUPDT    ucorupdt_
#  endif

#include <ArrayLim.H>

extern "C" {
#if (BL_SPACEDIM == 2)
    void FORT_MACDIV (Real* dmac, ARLIM_P(dlo), ARLIM_P(dhi),
                      const int* lo, const int* hi,
                      const Real* ux, ARLIM_P(ux_lo), ARLIM_P(ux_hi),
                      const Real* uy, ARLIM_P(uy_lo), ARLIM_P(uy_hi),
                      const Real* ax, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
                      const Real* ay, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
                      const Real* vol,ARLIM_P(vol_lo), ARLIM_P(vol_hi));

    void FORT_SCALEAREA (const Real* ax, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
                         const Real* ay, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
                         const Real* anel_coeff, 
                         const int* anel_coeff_lo, const int* anel_coeff_hi,
                         const int* lo, const int* hi, const int* mult);
#endif
#if (BL_SPACEDIM == 3)
    void FORT_MACDIV (Real* dmac, ARLIM_P(dlo), ARLIM_P(dhi),
                      const int* lo, const int* hi,
                      const Real* ux, ARLIM_P(ux_lo), ARLIM_P(ux_hi),
                      const Real* uy, ARLIM_P(uy_lo), ARLIM_P(uy_hi),
                      const Real* uz, ARLIM_P(uz_lo), ARLIM_P(uz_hi),
                      const Real* ax, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
                      const Real* ay, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
                      const Real* az, ARLIM_P(az_lo), ARLIM_P(az_hi),
                      const Real* vol,ARLIM_P(vol_lo), ARLIM_P(vol_hi));

    void FORT_SCALEAREA (const Real* ax, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
                         const Real* ay, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
                         const Real* az, ARLIM_P(az_lo), ARLIM_P(az_hi),
                         const Real* anel_coeff, 
                         const int* anel_coeff_lo, const int* anel_coeff_hi,
                         const int* lo, const int* hi, const int* mult);
#endif

    void FORT_UCORUPDT (Real* snew, const Real* sold,
                        ARLIM_P(s_lo), ARLIM_P(s_hi),
                        const int* nvar,
                        const int* lo, const int* hi,
                        const Real* dmac, ARLIM_P(dmac_lo), ARLIM_P(dmac_hi),
                        const Real* mult);

}

#endif
#endif /*_MACPROJ_F_H*/
