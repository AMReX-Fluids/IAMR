#ifdef BL_LANG_FORT
c $Id: DERIVE_F.H,v 1.1 1997-07-08 23:08:11 vince Exp $
#endif

#ifdef BL_LANG_FORT
#    define FORT_DERKENG     derkeng
#    define FORT_DERMGVORT   dermgvort
#    define FORT_DERVORT     dervort
#    define FORT_DERMVEL     dermvel
#    define FORT_DERLOGS     derlogs
#    define FORT_DERDVRHO    derdvrho
#    define FORT_DERMPRHO    dermprho
#    define FORT_DERLGRHODUST   derlgrc
#    define FORT_DERMGDIVU   dermgdivu
#    define FORT_GRADP_DIR   gradpdir
#    define FORT_DERGRDPX    dergrdpx
#    define FORT_DERGRDPY    dergrdpy
#if (BL_SPACEDIM == 3)
#    define FORT_DERGRDPZ    dergrdpz
#endif
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_DERKENG     DERKENG
#    define FORT_DERMGVORT   DERMGVORT
#    define FORT_DERVORT     DERVORT
#    define FORT_DERMVEL     DERMVEL
#    define FORT_DERLOGS     DERLOGS
#    define FORT_DERDVRHO    DERDVRHO
#    define FORT_DERMPRHO    DERMPRHO
#    define FORT_DERLGRHODUST   DERLGRC
#    define FORT_DERMGDIVU   DERMGDIVU
#    define FORT_GRADP_DIR   GRADPDIR
#    define FORT_DERGRDPX    DERGRDPX
#    define FORT_DERGRDPY    DERGRDPY
#if (BL_SPACEDIM == 3)
#    define FORT_DERGRDPZ    DERGRDPZ
#endif
#  else
#    define FORT_GRADP_DIR   gradpdir_
#    define FORT_DERKENG     derkeng_
#    define FORT_DERMGVORT   dermgvort_
#    define FORT_DERVORT     dervort_
#    define FORT_DERMVEL     dermvel_
#    define FORT_DERLOGS     derlogs_
#    define FORT_DERDVRHO    derdvrho_
#    define FORT_DERMPRHO    dermprho_
#    define FORT_DERLGRHODUST   derlgrc_
#    define FORT_DERMGDIVU   dermgdivu_
#    define FORT_GRADP_DIR   gradpdir_
#    define FORT_DERGRDPX    dergrdpx_
#    define FORT_DERGRDPY    dergrdpy_
#if (BL_SPACEDIM == 3)
#    define FORT_DERGRDPZ    dergrdpz_
#endif
#  endif

#include <ArrayLim.H>

  extern "C" {
     void FORT_DERKENG    (REAL* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const REAL* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const REAL* delta, const REAL* xlo,
			   const REAL* time, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERMGVORT  (REAL* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const REAL* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const REAL* delta, const REAL* xlo,
			   const REAL* time, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERVORT    (REAL* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const REAL* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const REAL* delta, const REAL* xlo,
			   const REAL* time, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERMVEL    (REAL* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const REAL* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const REAL* delta, const REAL* xlo,
			   const REAL* time, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERLOGS   (REAL* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const REAL* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const REAL* delta, const REAL* xlo,
			   const REAL* time, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERMPRHO   (REAL* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const REAL* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const REAL* delta, const REAL* xlo,
			   const REAL* time, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERDVRHO   (REAL* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const REAL* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const REAL* delta, const REAL* xlo,
			   const REAL* time, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERLGRHODUST   (REAL* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const REAL* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const REAL* delta, const REAL* xlo,
			   const REAL* time, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERMGDIVU  (REAL* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const REAL* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const REAL* delta, const REAL* xlo,
			   const REAL* time, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_GRADP_DIR (const REAL* p, ARLIM_P(p_lo), ARLIM_P(p_hi),
                          const REAL* gp, ARLIM_P(gp_lo), ARLIM_P(gp_hi),
                          const int* lo, const int* hi, const int* dir,
                          const REAL* dx );

     void FORT_DERGRDPX   (REAL* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const REAL* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const REAL* delta, const REAL* xlo,
			   const REAL* time, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERGRDPY   (REAL* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const REAL* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const REAL* delta, const REAL* xlo,
			   const REAL* time, const int* bcrec, 
                           const int* level, const int* grid_no) ;

#if (BL_SPACEDIM == 3)
     void FORT_DERGRDPZ   (REAL* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const REAL* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const REAL* delta, const REAL* xlo,
			   const REAL* time, const int* bcrec, 
                           const int* level, const int* grid_no) ;
#endif
  };
#endif
