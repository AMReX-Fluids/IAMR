#ifndef _PROB_F_H_
#define _PROB_F_H_

/*
** $Id: PROB_F.H,v 1.7 1999-02-25 17:35:38 car Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_INITDATA    initdata
#    define FORT_DENFILL     denfill
#    define FORT_ADVFILL     advfill
#    define FORT_XVELFILL    xvelfill
#    define FORT_YVELFILL    yvelfill
#    define FORT_ZVELFILL    zvelfill
#    define FORT_PRESFILL    presfill
#    define FORT_TEMPFILL    tempfill
#    define FORT_DIVUFILL    divufill
#    define FORT_DSDTFILL    dsdtfill
#    define FORT_DENERROR    denerror
#    define FORT_ADVERROR    adverror
#    define FORT_TEMPERROR   temperror
#    define FORT_MVERROR     mverror

#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_INITDATA    INITDATA
#    define FORT_DENFILL     DENFILL
#    define FORT_ADVFILL     ADVFILL
#    define FORT_XVELFILL    XVELFILL
#    define FORT_YVELFILL    YVELFILL
#    define FORT_ZVELFILL    ZVELFILL
#    define FORT_PRESFILL    PRESFILL
#    define FORT_TEMPFILL    TEMPFILL
#    define FORT_DIVUFILL    DIVUFILL
#    define FORT_DSDTFILL    DSDTFILL
#    define FORT_DENERROR    DENERROR
#    define FORT_ADVERROR    ADVERROR
#    define FORT_TEMPERROR   TEMPERROR
#    define FORT_MVERROR     MVERROR
#  else
#    define FORT_INITDATA    initdata_
#    define FORT_DENFILL     denfill_
#    define FORT_ADVFILL     advfill_
#    define FORT_XVELFILL    xvelfill_
#    define FORT_YVELFILL    yvelfill_
#    define FORT_ZVELFILL    zvelfill_
#    define FORT_PRESFILL    presfill_
#    define FORT_TEMPFILL    tempfill_
#    define FORT_DIVUFILL    divufill_
#    define FORT_DSDTFILL    dsdtfill_
#    define FORT_DENERROR    denerror_
#    define FORT_ADVERROR    adverror_
#    define FORT_TEMPERROR   temperror_
#    define FORT_MVERROR     mverror_
#  endif

#include "ArrayLim.H"

  extern "C" {

     void FORT_INITDATA (const int* level, const Real* time,
                         const int* lo, const int* hi,
                         const int* nscal,
                         Real* vel, Real* scal, 
                         ARLIM_P(s_lo), ARLIM_P(s_hi),
                         Real* press,
                         ARLIM_P(p_lo), ARLIM_P(p_hi),
                         const Real* dx,
                         const Real* xlo, const Real* xhi);

     void FORT_DENFILL  (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

     void FORT_ADVFILL  (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi), 
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

     void FORT_XVELFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

     void FORT_YVELFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

#if (BL_SPACEDIM == 3)
     void FORT_ZVELFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);
#endif
     void FORT_PRESFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

     void FORT_TEMPFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi), 
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

     void FORT_DIVUFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

     void FORT_DSDTFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

     void FORT_DENERROR  (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                          const int* tagval, const int* clearval,
                          Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                          const int* lo, const int* hi, const int* nvar,
                          const int* domain_lo, const int* domain_hi,
                          const Real* dx, const Real* xlo,
                          const Real* prob_lo, const Real* time,
                          const int* level);
 
     void FORT_ADVERROR  (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                          const int* tagval, const int* clearval,
                          Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                          const int* lo, const int* hi, const int* nvar,
                          const int* domain_lo, const int* domain_hi,
                          const Real* dx, const Real* xlo,
                          const Real* prob_lo, const Real* time,
                          const int* level);
 
     void FORT_TEMPERROR  (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                           const int* tagval, const int* clearval,
                           Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                           const int* lo, const int* hi, const int* nvar,
                           const int* domain_lo, const int* domain_hi,
                           const Real* dx, const Real* xlo,
                           const Real* prob_lo, const Real* time,
                           const int* level);

     void FORT_MVERROR   (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                          const int* tagval, const int* clearval,
                          Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                          const int* lo, const int* hi, const int* nvar,
                          const int* domain_lo, const int* domain_hi,
                          const Real* dx, const Real* xlo,
                          const Real* prob_lo, const Real* time,
                          const int* level);

  }
#endif
#endif /*_PROB_F_H_*/
