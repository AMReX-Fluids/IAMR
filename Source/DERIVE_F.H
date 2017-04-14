
#ifndef _DERIVE_F_H
#define _DERIVE_F_H

#ifdef BL_LANG_FORT
#    define FORT_DERMODGRADRHO     dermodgradrho
#    define FORT_DERUDOTLAPU     derudotlapu
#    define FORT_DERKENG     derkeng
#    define FORT_DERMGVORT   dermgvort
#if (BL_SPACEDIM == 3)
#    define FORT_DERVORTX    dervortx
#    define FORT_DERVORTY    dervorty
#    define FORT_DERVORTZ    dervortz
#    define FORT_DERDMAG     derdmag
#endif
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
#    define FORT_DERRADVEL   derradvel
#    define FORT_DERAZIVEL   derazivel
#    define FORT_DERXVELROT  derxvelrot
#    define FORT_DERYVELROT  deryvelrot
#    define FORT_DERMAGVELROT dermagvelrot
#    define FORT_DERMAGVORTROT dermagvortrot
#ifdef DO_IAMR_FORCE
#    define FORT_DERFORCING  derforcing
#    define FORT_DERFORCEX   derforcex
#    define FORT_DERFORCEY   derforcey
#    define FORT_DERFORCEZ   derforcez
#endif
#    define FORT_DERTURBVARS derturbvars
#    define FORT_DERPRESVARS derpresvars
#ifdef SUMJET
#    define FORT_DERJETVARS  derjetvars
#    define FORT_DERJETPRESVARS derjetpresvars
#endif
#endif
#    define FORT_DERNULL     dernull
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_DERMODGRADRHO     DERMODGRADRHO
#    define FORT_DERUDOTLAPU     DERUDOTLAPU
#    define FORT_DERKENG     DERKENG
#    define FORT_DERMGVORT   DERMGVORT
#if (BL_SPACEDIM == 3)
#    define FORT_DERVORTX    DERVORTX
#    define FORT_DERVORTY    DERVORTY
#    define FORT_DERVORTZ    DERVORTZ
#    define FORT_DERDMAG     DERDMAG
#endif
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
#    define FORT_DERRADVEL   DERRADVEL
#    define FORT_DERAZIVEL   DERAZIVEL
#    define FORT_DERXVELROT  DERXVELROT
#    define FORT_DERYVELROT  DERYVELROT
#    define FORT_DERMAGVELROT DERMAGVELROT
#    define FORT_DERMAGVORTROT DERMAGVORTROT
#ifdef DO_IAMR_FORCE
#    define FORT_DERFORCING  DERFORCING
#    define FORT_DERFORCEX   DERFORCEX
#    define FORT_DERFORCEY   DERFORCEY
#    define FORT_DERFORCEZ   DERFORCEZ
#endif
#    define FORT_DERTURBVARS DERTURBVARS
#    define FORT_DERPRESVARS DERPRESVARS
#ifdef SUMJET
#    define FORT_DERJETVARS  DERJETVARS
#    define FORT_DERJETPRESVARS DERJETPRESVARS
#endif
#endif
#    define FORT_DERNULL     DERNULL
#  elif defined (BL_FORT_USE_LOWERCASE)
#    define FORT_GRADP_DIR   gradpdir
#    define FORT_DERMODGRADRHO     dermodgradrho
#    define FORT_DERUDOTLAPU     derudotlapu
#    define FORT_DERKENG     derkeng
#    define FORT_DERMGVORT   dermgvort
#if (BL_SPACEDIM == 3)
#    define FORT_DERVORTX    dervortx
#    define FORT_DERVORTY    dervorty
#    define FORT_DERVORTZ    dervortz
#    define FORT_DERDMAG     derdmag
#endif
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
#    define FORT_DERRADVEL   derradvel
#    define FORT_DERAZIVEL   derazivel
#    define FORT_DERXVELROT  derxvelrot
#    define FORT_DERYVELROT  deryvelrot
#    define FORT_DERMAGVELROT dermagvelrot
#    define FORT_DERMAGVORTROT dermagvortrot
#ifdef DO_IAMR_FORCE
#    define FORT_DERFORCING  derforcing
#    define FORT_DERFORCEX   derforcex
#    define FORT_DERFORCEY   derforcey
#    define FORT_DERFORCEZ   derforcez
#endif
#    define FORT_DERTURBVARS derturbvars
#    define FORT_DERPRESVARS derpresvars
#ifdef SUMJET
#    define FORT_DERJETVARS  derjetvars
#    define FORT_DERJETPRESVARS derjetpresvars
#endif
#endif
#    define FORT_DERNULL     dernull
#  elif defined (BL_FORT_USE_UNDERSCORE)
#    define FORT_GRADP_DIR   gradpdir_
#    define FORT_DERMODGRADRHO     dermodgradrho_
#    define FORT_DERUDOTLAPU     derudotlapu_
#    define FORT_DERKENG     derkeng_
#    define FORT_DERMGVORT   dermgvort_
#if (BL_SPACEDIM == 3)
#    define FORT_DERVORTX    dervortx_
#    define FORT_DERVORTY    dervorty_
#    define FORT_DERVORTZ    dervortz_
#    define FORT_DERDMAG     derdmag_
#endif
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
#    define FORT_DERRADVEL   derradvel_
#    define FORT_DERAZIVEL   derazivel_
#    define FORT_DERXVELROT  derxvelrot_
#    define FORT_DERYVELROT  deryvelrot_
#    define FORT_DERMAGVELROT dermagvelrot_
#    define FORT_DERMAGVORTROT dermagvortrot_
#ifdef DO_IAMR_FORCE
#    define FORT_DERFORCING  derforcing_
#    define FORT_DERFORCEX   derforcex_
#    define FORT_DERFORCEY   derforcey_
#    define FORT_DERFORCEZ   derforcez_
#endif
#    define FORT_DERTURBVARS derturbvars_
#    define FORT_DERPRESVARS derpresvars_
#ifdef SUMJET
#    define FORT_DERJETVARS  derjetvars_
#    define FORT_DERJETPRESVARS derjetpresvars_
#endif
#endif
#    define FORT_DERNULL     dernull_
#  endif

#include <AMReX_ArrayLim.H>

extern "C" {
     void FORT_DERMODGRADRHO    (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERUDOTLAPU    (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

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
#if (BL_SPACEDIM == 3)

     void FORT_DERVORTX  (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERVORTY  (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERVORTZ  (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERDMAG     (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

#endif

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

     void FORT_DERRADVEL (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERAZIVEL (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERXVELROT (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERYVELROT (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERMAGVELROT (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERMAGVORTROT (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
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

#ifdef SUMJET
     void FORT_DERJETVARS(Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                          const int* nvar, const Real* compdat,
			  ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			  const int* lo, const int* hi,
			  const int* domain_lo, const int* domain_hi,
			  const Real* delta, const Real* xlo,
			  const Real* time, const Real* dt, const int* bcrec, 
                          const int* level, const int* grid_no) ;

     void FORT_DERJETPRESVARS(Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                              const int* nvar, const Real* compdat,
			      ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			      const int* lo, const int* hi,
			      const int* domain_lo, const int* domain_hi,
			      const Real* delta, const Real* xlo,
			      const Real* time, const Real* dt, const int* bcrec, 
			      const int* level, const int* grid_no) ;
#endif
#endif

     void FORT_DERNULL(Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                       const int* nvar, const Real* compdat,
		       ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
		       const int* lo, const int* hi,
       		       const int* domain_lo, const int* domain_hi,
		       const Real* delta, const Real* xlo,
		       const Real* time, const Real* dt, const int* bcrec, 
                       const int* level, const int* grid_no) ;
}
#endif
#endif
