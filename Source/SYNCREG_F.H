
#ifndef _SYNCREG_F_H_
#define _SYNCREG_F_H_

#ifndef BL_LANG_FORT

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C" 
{
#endif
    void srcrsereg (const amrex::Real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                         amrex::Real* crse, ARLIM_P(clo), ARLIM_P(chi),
                         const int* lo, const int* hi,
                         const int* dir, const int* ratio);

    void makemask (amrex::Real* mask , ARLIM_P(mlo), ARLIM_P(mhi),
                        amrex::Real* cells, ARLIM_P(clo), ARLIM_P(chi));

    void convertmask (amrex::Real* mask , ARLIM_P(mlo), ARLIM_P(mhi));

#ifdef __cplusplus
}
#endif
#endif

#endif /*_SYNCREG_F_H_*/
