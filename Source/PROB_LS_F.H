#ifndef _PROB_LS_F_H_
#define _PROB_LS_F_H_

#ifdef BL_LANG_FORT
#    define FORT_PHIFILL     phifill
#    define FORT_INITPHI     initphi
#    define FORT_LSERROR     lserror
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_PHIFILL     PHIFILL
#    define FORT_INITPHI     INITPHI
#    define FORT_LSERROR     LSERROR
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_PHIFILL     phifill
#    define FORT_INITPHI     initphi
#    define FORT_LSERROR     lserror
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_PHIFILL     phifill_
#    define FORT_INITPHI     initphi_
#    define FORT_LSERROR     lserror_
#  endif

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
  extern "C" {
#endif


    void FORT_INITPHI(const int* lo,const int* hi,amrex::Real * phi,
                      ARLIM_P(lo), ARLIM_P(hi),
                      const amrex::Real * dx,const amrex::Real * x_lo,
                      int* type, ARLIM_P(type_lo), ARLIM_P(type_hi));

     void FORT_PHIFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
	                          const int* dom_lo, const int* dom_hi,
	                          const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

     void FORT_LSERROR   (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                          const int* tagval, const int* clearval,
                          amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                          const int* lo, const int* hi, const int* nvar,
                          const int* domain_lo, const int* domain_hi,
                          const amrex::Real* dx, const amrex::Real* xlo,
                          const amrex::Real* prob_lo, const amrex::Real* time,
                          const int* level);
#ifdef __cplusplus
  }
#endif
#endif
#endif 
