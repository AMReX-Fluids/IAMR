
#ifndef _MACPROJ_F_H
#define _MACPROJ_F_H

#ifdef BL_LANG_FORT
#    define FORT_MACDIV      macdiv
#    define FORT_SCALEAREA   scalearea
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_MACDIV      MACDIV
#    define FORT_SCALEAREA   SCALEAREA
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_MACDIV      macdiv
#    define FORT_SCALEAREA   scalearea
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_MACDIV      macdiv_
#    define FORT_SCALEAREA   scalearea_
#  endif

#include <AMReX_ArrayLim.H>

extern "C" {
#if (BL_SPACEDIM == 2)
    void FORT_MACDIV (amrex::Real* dmac, ARLIM_P(dlo), ARLIM_P(dhi),
                      const int* lo, const int* hi,
                      const amrex::Real* ux, ARLIM_P(ux_lo), ARLIM_P(ux_hi),
                      const amrex::Real* uy, ARLIM_P(uy_lo), ARLIM_P(uy_hi),
                      const amrex::Real* ax, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
                      const amrex::Real* ay, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
                      const amrex::Real* vol,ARLIM_P(vol_lo), ARLIM_P(vol_hi));

    void FORT_SCALEAREA (const amrex::Real* ax, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
                         const amrex::Real* ay, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
                         const amrex::Real* anel_coeff, 
                         const int* anel_coeff_lo, const int* anel_coeff_hi,
                         const int* lo, const int* hi, const int* mult);
#endif
#if (BL_SPACEDIM == 3)
    void FORT_MACDIV (amrex::Real* dmac, ARLIM_P(dlo), ARLIM_P(dhi),
                      const int* lo, const int* hi,
                      const amrex::Real* ux, ARLIM_P(ux_lo), ARLIM_P(ux_hi),
                      const amrex::Real* uy, ARLIM_P(uy_lo), ARLIM_P(uy_hi),
                      const amrex::Real* uz, ARLIM_P(uz_lo), ARLIM_P(uz_hi),
                      const amrex::Real* ax, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
                      const amrex::Real* ay, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
                      const amrex::Real* az, ARLIM_P(az_lo), ARLIM_P(az_hi),
                      const amrex::Real* vol,ARLIM_P(vol_lo), ARLIM_P(vol_hi));

    void FORT_SCALEAREA (const amrex::Real* ax, ARLIM_P(ax_lo), ARLIM_P(ax_hi),
                         const amrex::Real* ay, ARLIM_P(ay_lo), ARLIM_P(ay_hi),
                         const amrex::Real* az, ARLIM_P(az_lo), ARLIM_P(az_hi),
                         const amrex::Real* anel_coeff, 
                         const int* anel_coeff_lo, const int* anel_coeff_hi,
                         const int* lo, const int* hi, const int* mult);
#endif

}

#endif
#endif /*_MACPROJ_F_H*/
