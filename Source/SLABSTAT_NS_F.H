
#ifndef _SLABSTAT_NS_F_H_
#define _SLABSTAT_NS_F_H_

#ifdef BL_LANG_FORT
#    define FORT_NS_BASICSTATS_NCTRAC    ns_basicstats_nctrac
#    define FORT_NS_BASICSTATS_CTRAC     ns_basicstats_ctrac
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_NS_BASICSTATS_NCTRAC    NS_BASICSTATS_NCTRAC
#    define FORT_NS_BASICSTATS_CTRAC     NS_BASICSTATS_CTRAC
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_NS_BASICSTATS_NCTRAC    ns_basicstats_nctrac
#    define FORT_NS_BASICSTATS_CTRAC     ns_basicstats_ctrac
#  else
#    define FORT_NS_BASICSTATS_NCTRAC    ns_basicstats_nctrac_
#    define FORT_NS_BASICSTATS_CTRAC     ns_basicstats_ctrac_
#  endif
#  endif

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
     void FORT_NS_BASICSTATS_NCTRAC (const amrex::Real* src, 
                                       ARLIM_P(src_lo), ARLIM_P(src_hi),
                                     const int* nsrc,
                                     amrex::Real* dst, 
                                       ARLIM_P(dst_lo), ARLIM_P(dst_hi),
                                     const int* ndst,
                                     const amrex::Real* dt,
                                     const amrex::Real* dx);

     void FORT_NS_BASICSTATS_CTRAC (const amrex::Real* src, 
                                      ARLIM_P(src_lo), ARLIM_P(src_hi),
                                    const int* nsrc,
                                    amrex::Real* dst, 
                                      ARLIM_P(dst_lo), ARLIM_P(dst_hi),
                                    const int* ndst,
                                    const amrex::Real* dt,
                                    const amrex::Real* dx);
#ifdef __cplusplus
}
#endif

#endif

#endif /*_SLABSTAT_NS_F_H_*/
