
#ifndef _DERIVE_F_H
#define _DERIVE_F_H

#ifndef BL_LANG_FORT

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C" {
#endif
     void dermodgradrho    (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void derudotlapu    (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void derkeng    (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void dermgvort  (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;
#if (BL_SPACEDIM == 3)

     void dervortx  (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void dervorty  (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void dervortz  (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void derdmag     (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

#endif

     void dervort    (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void dermvel    (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void derlogs   (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void dermprho   (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void derdvrho   (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void derlgrhodust   (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void dermgdivu  (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void gradp_dir (const amrex::Real* p, ARLIM_P(p_lo), ARLIM_P(p_hi),
                          const amrex::Real* gp, ARLIM_P(gp_lo), ARLIM_P(gp_hi),
                          const int* lo, const int* hi, const int* dir,
                          const amrex::Real* dx );

     void deravgpres (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void dergrdp    (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void dergrdpx   (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void dergrdpy   (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

#if (BL_SPACEDIM == 3)
     void dergrdpz   (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void derradvel (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void derazivel (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void derxvelrot (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void deryvelrot (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void dermagvelrot (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void dermagvortrot (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

#ifdef DO_IAMR_FORCE
     void derforcing (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void derforcex  (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void derforcey  (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void derforcez  (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;
#endif

     void derturbvars(amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void derpresvars(amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

#ifdef SUMJET
     void derjetvars(amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                          const int* nvar, const amrex::Real* compdat,
			  ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			  const int* lo, const int* hi,
			  const int* domain_lo, const int* domain_hi,
			  const amrex::Real* delta, const amrex::Real* xlo,
			  const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                          const int* level, const int* grid_no) ;

     void derjetpresvars(amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                              const int* nvar, const amrex::Real* compdat,
			      ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			      const int* lo, const int* hi,
			      const int* domain_lo, const int* domain_hi,
			      const amrex::Real* delta, const amrex::Real* xlo,
			      const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
			      const int* level, const int* grid_no) ;
#endif
#endif

     void dernull(amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                       const int* nvar, const amrex::Real* compdat,
		       ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
		       const int* lo, const int* hi,
       		       const int* domain_lo, const int* domain_hi,
		       const amrex::Real* delta, const amrex::Real* xlo,
		       const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                       const int* level, const int* grid_no) ;
#ifdef __cplusplus
}
#endif
#endif
#endif
