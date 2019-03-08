
#ifndef _SLABSTAT_NS_F_H_
#define _SLABSTAT_NS_F_H_

#ifndef BL_LANG_FORT

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
     void ns_basicstats_nctrac (const amrex::Real* src, 
                                       ARLIM_P(src_lo), ARLIM_P(src_hi),
                                     const int* nsrc,
                                     amrex::Real* dst, 
                                       ARLIM_P(dst_lo), ARLIM_P(dst_hi),
                                     const int* ndst,
                                     const amrex::Real* dt,
                                     const amrex::Real* dx);

     void ns_basicstats_ctrac (const amrex::Real* src, 
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
