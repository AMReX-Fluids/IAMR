#ifndef _VISCOPERATOR_F_H_
#define _VISCOPERATOR_F_H_  "$Id: VISCOPERATOR_F.H,v 1.4 1997-12-11 23:30:25 lijewski Exp $"

#ifdef BL_LANG_FORT
#    define FORT_VISCEXTRAP   viscextrap
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_VISCEXTRAP    VISCEXTRAP
#  else
#    define FORT_VISCEXTRAP    viscextrap_
#endif

#include <ArrayLim.H>

extern "C"
{
    void FORT_VISCEXTRAP (Real* visc, ARLIM_P(visclo), ARLIM_P(isclo),
                          const int* lo, const int* hi,
                          const int* comp);
}
#endif

#endif
