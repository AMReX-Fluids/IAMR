#ifndef _PROB_NS_F_H_
#define _PROB_NS_F_H_

#ifndef BL_LANG_FORT

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
  extern "C" {
#endif

     void FORT_INITDATA (const int* level, const amrex::Real* time,
                         const int* lo, const int* hi,
                         const int* nscal,
                         amrex::Real* vel, amrex::Real* scal, 
                         ARLIM_P(s_lo), ARLIM_P(s_hi),
                         amrex::Real* press,
                         ARLIM_P(p_lo), ARLIM_P(p_hi),
                         const amrex::Real* dx,
                         const amrex::Real* xlo, const amrex::Real* xhi);

     void FORT_XVELFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);
     
     void FORT_YVELFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);
 
     void FORT_ZVELFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

     void FORT_DENFILL  (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

     void FORT_ADVFILL  (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi), 
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

     void FORT_ADV2FILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi), 
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

     void FORT_PRESFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

     void FORT_TEMPFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi), 
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

     void FORT_DIVUFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

     void FORT_DSDTFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

     void FORT_AVERAGE_EDGE_STATES(      BL_FORT_FAB_ARG_ANYD(Vel),
                                   const BL_FORT_FAB_ARG_ANYD(umacx),
                                   const BL_FORT_FAB_ARG_ANYD(umacy),
#if (AMREX_SPACEDIM==3)
                                   const BL_FORT_FAB_ARG_ANYD(umacz),
#endif
                                   const int* getForceVerbose);

#ifdef __cplusplus
  }
#endif
#endif

#endif /*_PROB_NS_F_H_*/
