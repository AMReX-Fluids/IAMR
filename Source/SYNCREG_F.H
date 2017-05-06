
#ifndef _SYNCREG_F_H_
#define _SYNCREG_F_H_

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

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C" 
{
#endif
    void FORT_SRCRSEREG (const amrex::Real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                         amrex::Real* crse, ARLIM_P(clo), ARLIM_P(chi),
                         const int* lo, const int* hi,
                         const int* dir, const int* ratio);

    void FORT_MAKEMASK (amrex::Real* mask , ARLIM_P(mlo), ARLIM_P(mhi),
                        amrex::Real* cells, ARLIM_P(clo), ARLIM_P(chi));

    void FORT_CONVERTMASK (amrex::Real* mask , ARLIM_P(mlo), ARLIM_P(mhi));

#ifdef __cplusplus
}
#endif
#endif

#endif /*_SYNCREG_F_H_*/
