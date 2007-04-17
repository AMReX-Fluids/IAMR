#ifndef _DERIVE_F_H
#define _DERIVE_F_H

/*
** $Id: DERIVE_F.H,v 1.18 2007-04-17 16:33:15 aaspden Exp $
*/

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
#    define FORT_DERAVGPRES  deravgpres
#    define FORT_DERGRDP     dergrdp
#    define FORT_DERGRDPX    dergrdpx
#    define FORT_DERGRDPY    dergrdpy
#if (BL_SPACEDIM == 3)
#    define FORT_DERGRDPZ    dergrdpz
#ifdef DO_IAMR_FORCE
#    define FORT_DERFORCING  derforcing
#    define FORT_DERFORCEX   derforcex
#    define FORT_DERFORCEY   derforcey
#    define FORT_DERFORCEZ   derforcez
#    define FORT_DERTURBVARS derturbvars
#    define FORT_DERPRESVARS derpresvars
#endif
#endif
#else
#  if defined(BL_FORT_USE_UPPERCASE)
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
#    define FORT_DERAVGPRES  DERAVGPRES
#    define FORT_DERGRDP    DERGRDP
#    define FORT_DERGRDPX    DERGRDPX
#    define FORT_DERGRDPY    DERGRDPY
#if (BL_SPACEDIM == 3)
#    define FORT_DERGRDPZ    DERGRDPZ
#ifdef DO_IAMR_FORCE
#    define FORT_DERFORCING  DERFORCING
#    define FORT_DERFORCEX   DERFORCEX
#    define FORT_DERFORCEY   DERFORCEY
#    define FORT_DERFORCEZ   DERFORCEZ
#    define FORT_DERTURBVARS DERTURBVARS
#    define FORT_DERPRESVARS DERPRESVARS
#endif
#endif
#  elif defined (BL_FORT_USE_LOWERCASE)
#    define FORT_GRADP_DIR   gradpdir
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
#    define FORT_DERAVGPRES  deravgpres
#    define FORT_DERGRDP     dergrdp
#    define FORT_DERGRDPX    dergrdpx
#    define FORT_DERGRDPY    dergrdpy
#if (BL_SPACEDIM == 3)
#    define FORT_DERGRDPZ    dergrdpz
#ifdef DO_IAMR_FORCE
#    define FORT_DERFORCING  derforcing
#    define FORT_DERFORCEX   derforcex
#    define FORT_DERFORCEY   derforcey
#    define FORT_DERFORCEZ   derforcez
#    define FORT_DERTURBVARS derturbvars
#    define FORT_DERPRESVARS derpresvars
#endif
#endif
#  elif defined (BL_FORT_USE_UNDERSCORE)
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
#    define FORT_DERAVGPRES  deravgpres_
#    define FORT_DERGRDP     dergrdp_
#    define FORT_DERGRDPX    dergrdpx_
#    define FORT_DERGRDPY    dergrdpy_
#if (BL_SPACEDIM == 3)
#    define FORT_DERGRDPZ    dergrdpz_
#ifdef DO_IAMR_FORCE
#    define FORT_DERFORCING  derforcing_
#    define FORT_DERFORCEX   derforcex_
#    define FORT_DERFORCEY   derforcey_
#    define FORT_DERFORCEZ   derforcez_
#    define FORT_DERTURBVARS derturbvars_
#    define FORT_DERPRESVARS derpresvars_
#endif
#endif
#  endif

#include <ArrayLim.H>

extern "C" {
     void FORT_DERKENG    (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERMGVORT  (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERVORT    (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERMVEL    (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERLOGS   (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERMPRHO   (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERDVRHO   (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERLGRHODUST   (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERMGDIVU  (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_GRADP_DIR (const Real* p, ARLIM_P(p_lo), ARLIM_P(p_hi),
                          const Real* gp, ARLIM_P(gp_lo), ARLIM_P(gp_hi),
                          const int* lo, const int* hi, const int* dir,
                          const Real* dx );

     void FORT_DERAVGPRES (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERGRDP    (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERGRDPX   (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERGRDPY   (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

#if (BL_SPACEDIM == 3)
     void FORT_DERGRDPZ   (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;
#ifdef DO_IAMR_FORCE
     void FORT_DERFORCING (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERFORCEX  (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERFORCEY  (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERFORCEZ  (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;
#endif
     void FORT_DERTURBVARS(Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERPRESVARS(Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

#endif
}
#endif
#endif
