#ifndef _SYNCREG_F_H_
#define _SYNCREG_F_H_

/*
** $Id: SYNCREG_F.H,v 1.6 1999-06-30 22:28:41 almgren Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_SRDIVU      srdivu
#    define FORT_SRDGPHI     srdgphi
#    define FORT_SRCRSEREG   srcrsereg
#    define FORT_MAKEMASK    makemask
#    define FORT_CONVERTMASK convmask
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_SRDIVU      SRDIVU
#    define FORT_SRDGPHI     SRDGPHI
#    define FORT_SRCRSEREG   SRCRSEREG
#    define FORT_MAKEMASK    MAKEMASK
#    define FORT_CONVERTMASK CONVMASK
#  else
#    define FORT_SRDIVU      srdivu_
#    define FORT_SRDGPHI     srdgphi_
#    define FORT_SRCRSEREG   srcrsereg_
#    define FORT_MAKEMASK    makemask_
#    define FORT_CONVERTMASK convmask_
#  endif

#include <ArrayLim.H>

extern "C" 
{
    void FORT_SRDIVU (const Real* u, ARLIM_P(ulo), ARLIM_P(uhi),
                      Real* divu, ARLIM_P(dlo), ARLIM_P(dhi),
                      const int* lo, const int* hi,
                      const Real* dx, const Real* mult, int * is_rz);

    void FORT_SRDGPHI (const Real* phi, ARLIM_P(p_lo), ARLIM_P(p_hi),
                       const Real* sigma, ARLIM_P(slo), ARLIM_P(shi),
                       Real* divgp, ARLIM_P(dlo), ARLIM_P(dhi),
                       const int* lo, const int* hi,
                       const Real* dx, const Real* mult, int * is_rz);

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
