#ifndef _SYNCREG_F_H_
#define _SYNCREG_F_H_  $Id: SYNCREG_F.H,v 1.1 1997-07-08 23:08:15 vince Exp $

#ifdef BL_LANG_FORT
#    define FORT_SRDIVU      srdivu
#    define FORT_SRDGPPHI    srdgpphi
#    define FORT_SRDGPHI     srdgphi
#    define FORT_SRCRSEREG   srcrsereg
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_SRDIVU      SRDIVU
#    define FORT_SRDGPPHI    SRDGPPHI
#    define FORT_SRDGPHI     SRDGPHI
#    define FORT_SRCRSEREG   SRCRSEREG
#  else
#    define FORT_SRDIVU      srdivu_
#    define FORT_SRDGPPHI    srdgpphi_
#    define FORT_SRDGPHI     srdgphi_
#    define FORT_SRCRSEREG   srcrsereg_
#  endif

#include <ArrayLim.H>

extern "C" 
{
    void FORT_SRDIVU (const REAL* u, ARLIM_P(ulo), ARLIM_P(uhi),
                      REAL* divu, ARLIM_P(dlo), ARLIM_P(dhi),
                      const int* lo, const int* hi,
                      const REAL* dx, const REAL* mult, int * is_rz);

    void FORT_SRDGPPHI (const REAL* phi, ARLIM_P(p_lo), ARLIM_P(p_hi),
                        const REAL* sigma, ARLIM_P(slo), ARLIM_P(shi),
                        REAL* divgp, ARLIM_P(dglo), ARLIM_P(dghi),
                        const int* lo, const int* hi,
                        const REAL* dx, const REAL* mult, int * is_rz);

    void FORT_SRDGPHI (const REAL* phi, ARLIM_P(p_lo), ARLIM_P(p_hi),
                       const REAL* sigma, ARLIM_P(slo), ARLIM_P(shi),
                       REAL* divgp, ARLIM_P(dlo), ARLIM_P(dhi),
                       const int* lo, const int* hi,
                       const REAL* dx, const REAL* mult, int * is_rz);

    void FORT_SRCRSEREG (const REAL* fine, ARLIM_P(flo), ARLIM_P(fhi),
                         REAL* crse, ARLIM_P(clo), ARLIM_P(chi),
                         const int* lo, const int* hi,
                         const int* dir, const int* ratio);

};
#endif

#endif


