#ifndef _VISCOPERATOR_F_H_
#define _VISCOPERATOR_F_H_

/*
** $Id: VISCOPERATOR_F.H,v 1.6 1999-08-25 18:04:05 propp Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_VISCEXTRAP   viscextrap
#endif

#ifdef BL_LANG_CC
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_VISCEXTRAP    VISCEXTRAP
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_VISCEXTRAP    viscextrap
#  elif defined(BL_FORT_USE_UNDERSCORE)
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

#endif /*_VISCOPERATOR_F_H_*/
