#ifndef _SYNCREG_F_H_
#define _SYNCREG_F_H_

/*
** $Id: SYNCREG_F.H,v 1.11 2002-04-26 22:10:02 lijewski Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_SRCRSEREG   srcrsereg
#    define FORT_MAKEMASK    makemask
#    define FORT_CONVERTMASK convmask
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_SRCRSEREG   SRCRSEREG
#    define FORT_MAKEMASK    MAKEMASK
#    define FORT_CONVERTMASK CONVMASK
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_SRCRSEREG   srcrsereg
#    define FORT_MAKEMASK    makemask
#    define FORT_CONVERTMASK convmask
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_SRCRSEREG   srcrsereg_
#    define FORT_MAKEMASK    makemask_
#    define FORT_CONVERTMASK convmask_
#  endif

#include <ArrayLim.H>

extern "C" 
{
    void FORT_SRCRSEREG (const Real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                         Real* crse, ARLIM_P(clo), ARLIM_P(chi),
                         const int* lo, const int* hi,
                         const int* dir, const int* ratio);

    void FORT_MAKEMASK (Real* mask , ARLIM_P(mlo), ARLIM_P(mhi),
                        Real* cells, ARLIM_P(clo), ARLIM_P(chi));

    void FORT_CONVERTMASK (Real* mask , ARLIM_P(mlo), ARLIM_P(mhi));

}
#endif

#endif /*_SYNCREG_F_H_*/
