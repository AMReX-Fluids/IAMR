#ifndef _MACPROJ_F_H
#define _MACPROJ_F_H  $Id: MACPROJ_F.H,v 1.1 1997-07-08 23:08:12 vince Exp $

#ifdef BL_LANG_FORT
#    define FORT_MACDIV      macdiv
#    define FORT_UCORUPDT    ucorupdt
#    define FORT_MACPHIBC     macphibc
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
    void FORT_MACDIV (REAL* dmac, ARLIM_P(dlo), ARLIM_P(dhi),
                      const int* lo, const int* hi,
                      const REAL* ux, ARLIM_P(ux_lo), ARLIM_P(ux_hi),
                      const REAL* uy, ARLIM_P(uy_lo), ARLIM_P(uy_hi),
                      const REAL* ax, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
                      const REAL* ay, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
                      const REAL* vol,ARLIM_P(vol_lo), ARLIM_P(vol_hi));
#endif
#if (BL_SPACEDIM == 3)
    void FORT_MACDIV (REAL* dmac, ARLIM_P(dlo), ARLIM_P(dhi),
                      const int* lo, const int* hi,
                      const REAL* ux, ARLIM_P(ux_lo), ARLIM_P(ux_hi),
                      const REAL* uy, ARLIM_P(uy_lo), ARLIM_P(uy_hi),
                      const REAL* uz, ARLIM_P(uz_lo), ARLIM_P(uz_hi),
                      const REAL* ax, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
                      const REAL* ay, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
                      const REAL* az, ARLIM_P(az_lo), ARLIM_P(az_hi),
                      const REAL* vol,ARLIM_P(vol_lo), ARLIM_P(vol_hi));
#endif

    void FORT_UCORUPDT (REAL* snew, const REAL* sold,
                        ARLIM_P(s_lo), ARLIM_P(s_hi),
                        const int* nvar,
                        const int* lo, const int* hi,
                        const REAL* dmac, ARLIM_P(dmac_lo), ARLIM_P(dmac_hi),
                        const REAL* mult);

    void FORT_MACPHIBC(ARLIM_P(u_lo),ARLIM_P(u_hi),const REAL* u,
                       ARLIM_P(div_lo),ARLIM_P(div_hi),const REAL* divcor,
                       ARLIM_P(sig_lo),ARLIM_P(sig_hi),const REAL* sigmax,
                       ARLIM_P(rcen_lo),ARLIM_P(rcen_hi),const REAL* rcen, 
                       ARLIM_P(rhalf_lo),ARLIM_P(rhalf_hi),const REAL* rhalf, 
                       const REAL* hx,
                       ARLIM_P(phi_lo),ARLIM_P(phi_hi),const REAL* phi);
};

#endif
#endif


