#ifndef _MACOPERATOR_F_H_
#define _MACOPERATOR_F_H_  "$Id: MACOPERATOR_F.H,v 1.2 1997-07-17 22:01:03 vince Exp $"

#ifdef BL_LANG_FORT
#    define FORT_MACCOEF       maccoef
#    define FORT_MACRHS        macrhs
#    define FORT_MACUPDATE     macupdate
#    define FORT_MACSYNCRHS    macsyncrhs
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_MACCOEF       MACCOEF
#    define FORT_MACRHS        MACRHS
#    define FORT_MACUPDATE     MACUPDATE
#    define FORT_MACSYNCRHS    MACSYNCRHS
#  else
#    define FORT_MACCOEF       maccoef_
#    define FORT_MACRHS        macrhs_
#    define FORT_MACUPDATE     macupdate_
#    define FORT_MACSYNCRHS    macsyncrhs_
#endif

#include <ArrayLim.H>

extern "C"
{
#if (BL_SPACEDIM == 2)
    void FORT_MACCOEF (REAL* cx, ARLIM_P(cxlo), ARLIM_P(cxhi),
		       REAL* cy, ARLIM_P(cylo), ARLIM_P(cyhi),
		       const REAL* ax, ARLIM_P(axlo), ARLIM_P(axhi),
		       const REAL* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
		       const REAL* rho, ARLIM_P(rlo), ARLIM_P(rhi),
		       const int* lo, const int* hi,
		       const REAL* dx);

    void FORT_MACRHS  (const REAL* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
		       const REAL* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
		       const REAL* ax, ARLIM_P(axlo), ARLIM_P(axhi),
		       const REAL* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
		       const REAL* vol, ARLIM_P(vlo), ARLIM_P(vhi),
		       REAL* rhs, ARLIM_P(rlo), ARLIM_P(rhi),
		       const int* lo, const int* hi,
		       const REAL* scale);

    void FORT_MACUPDATE (const int *init,
                         REAL* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
			 REAL* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
			 const REAL* phi, ARLIM_P(p_lo), ARLIM_P(p_hi),
			 const REAL* rho, ARLIM_P(rlo), ARLIM_P(rhi),
			 const REAL* dx, const REAL* scale);
#endif

#if (BL_SPACEDIM == 3)
    void FORT_MACCOEF (REAL* cx, ARLIM_P(cxlo), ARLIM_P(cxhi),
                       REAL* cy, ARLIM_P(cylo), ARLIM_P(cyhi),
                       REAL* cz, ARLIM_P(czlo), ARLIM_P(czhi),
                       const REAL* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                       const REAL* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
                       const REAL* az, ARLIM_P(azlo), ARLIM_P(azhi),
                       const REAL* rho, ARLIM_P(rlo), ARLIM_P(rhi),
                       const int* lo, const int* hi,
                       const REAL* dx);

    void FORT_MACRHS  (const REAL* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                       const REAL* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                       const REAL* uz, ARLIM_P(uzlo), ARLIM_P(uzhi),
                       const REAL* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                       const REAL* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
                       const REAL* az, ARLIM_P(azlo), ARLIM_P(azhi),
                       const REAL* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                       REAL* rhs, ARLIM_P(rlo), ARLIM_P(rhi),
                       const int* lo, const int* hi,
                       const REAL* scale);

    void FORT_MACUPDATE ( const int *init,
                          REAL* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                          REAL* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                          REAL* uz, ARLIM_P(uzlo), ARLIM_P(uzhi),
                          const REAL* phi, ARLIM_P(p_lo), ARLIM_P(p_hi),
                          const REAL* rho, ARLIM_P(rlo), ARLIM_P(rhi),
                          const REAL* dx, const REAL* scale); 

#endif

    void FORT_MACSYNCRHS (REAL* rhs, ARLIM_P(rlo), ARLIM_P(rhi),
			  const int* lo, const int* hi,
			  const REAL* vol, ARLIM_P(vlo), ARLIM_P(vhi),
			  const REAL* dx, const REAL* rhs_scale);
}
#endif

#endif
