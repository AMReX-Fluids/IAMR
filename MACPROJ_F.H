#ifndef _MACPROJ_F_H
#define _MACPROJ_F_H

/*
** $Id: MACPROJ_F.H,v 1.4 1997-12-14 23:35:59 lijewski Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_MACDIV      macdiv
#    define FORT_UCORUPDT    ucorupdt
#    define FORT_MACPHIBC    macphibc
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_MACDIV      MACDIV
#    define FORT_UCORUPDT    UCORUPDT
#    define FORT_MACPHIBC    MACPHIBC
#  else
#    define FORT_MACDIV      macdiv_
#    define FORT_UCORUPDT    ucorupdt_
#    define FORT_MACPHIBC    macphibc_
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
#endif

    void FORT_UCORUPDT (Real* snew, const Real* sold,
                        ARLIM_P(s_lo), ARLIM_P(s_hi),
                        const int* nvar,
                        const int* lo, const int* hi,
                        const Real* dmac, ARLIM_P(dmac_lo), ARLIM_P(dmac_hi),
                        const Real* mult);

    void FORT_MACPHIBC(ARLIM_P(u_lo),ARLIM_P(u_hi),const Real* u,
                       ARLIM_P(div_lo),ARLIM_P(div_hi),const Real* divcor,
                       ARLIM_P(sig_lo),ARLIM_P(sig_hi),const Real* sigmax,
                       ARLIM_P(rcen_lo),ARLIM_P(rcen_hi),const Real* rcen, 
                       ARLIM_P(rhalf_lo),ARLIM_P(rhalf_hi),const Real* rhalf, 
                       const Real* hx,
                       ARLIM_P(phi_lo),ARLIM_P(phi_hi),const Real* phi);
}

#endif
#endif /*_MACPROJ_F_H*/
