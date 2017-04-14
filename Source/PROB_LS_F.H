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

  extern "C" {


    void FORT_INITPHI(const int* lo,const int* hi,Real * phi,
                      ARLIM_P(lo), ARLIM_P(hi),
                      const Real * dx,const Real * x_lo,
                      int* type, ARLIM_P(type_lo), ARLIM_P(type_hi));

     void FORT_PHIFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
	                          const int* dom_lo, const int* dom_hi,
	                          const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

     void FORT_LSERROR   (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                          const int* tagval, const int* clearval,
                          Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                          const int* lo, const int* hi, const int* nvar,
                          const int* domain_lo, const int* domain_hi,
                          const Real* dx, const Real* xlo,
                          const Real* prob_lo, const Real* time,
                          const int* level);
  }
#endif
#endif 
