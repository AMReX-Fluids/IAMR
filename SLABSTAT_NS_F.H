#ifndef _SLABSTAT_NS_F_H_
#define _SLABSTAT_NS_F_H_

/*
** $Id: SLABSTAT_NS_F.H,v 1.1 2000-05-22 19:57:43 sstanley Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_NS_BASICSTATS_NCTRAC    ns_basicstats_nctrac
#    define FORT_NS_BASICSTATS_CTRAC     ns_basicstats_ctrac
#endif

#ifdef BL_LANG_CC
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

extern "C"
{
     void FORT_NS_BASICSTATS_NCTRAC (const Real* src, 
                                       ARLIM_P(src_lo), ARLIM_P(src_hi),
                                     const int* nsrc,
                                     Real* dst, 
                                       ARLIM_P(dst_lo), ARLIM_P(dst_hi),
                                     const int* ndst,
                                     const Real* dt,
                                     const Real* dx);

     void FORT_NS_BASICSTATS_CTRAC (const Real* src, 
                                      ARLIM_P(src_lo), ARLIM_P(src_hi),
                                    const int* nsrc,
                                    Real* dst, 
                                      ARLIM_P(dst_lo), ARLIM_P(dst_hi),
                                    const int* ndst,
                                    const Real* dt,
                                    const Real* dx);
}

#endif

#endif /*_SLABSTAT_NS_F_H_*/
