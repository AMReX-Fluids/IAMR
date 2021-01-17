
#ifndef _DERIVE_F_H
#define _DERIVE_F_H

#ifndef BL_LANG_FORT

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C" {
#endif
     void dermodgradrho (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                         const int* nvar, const amrex::Real* compdat,
                         ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
                         const int* lo, const int* hi,
                         const int* domain_lo, const int* domain_hi,
                         const amrex::Real* delta, const amrex::Real* xlo,
                         const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                         const int* level, const int* grid_no) ;

     void derudotlapu (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                       const int* nvar, const amrex::Real* compdat,
                       ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
                       const int* lo, const int* hi,
                       const int* domain_lo, const int* domain_hi,
                       const amrex::Real* delta, const amrex::Real* xlo,
                       const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                       const int* level, const int* grid_no) ;

     void derkeng (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                   const int* nvar, const amrex::Real* compdat,
                   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
                   const int* lo, const int* hi,
                   const int* domain_lo, const int* domain_hi,
                   const amrex::Real* delta, const amrex::Real* xlo,
                   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                   const int* level, const int* grid_no) ;

     void dermgvort (BL_FORT_FAB_ARG_ANYD(der),const int* nvar,
                     const BL_FORT_FAB_ARG_ANYD(data),const int* ncomp,
                     AMREX_ARLIM_P(vlo), AMREX_ARLIM_P(vhi),
                     AMREX_ARLIM_P(dlo), AMREX_ARLIM_P(dhi),
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

     void dermprho (BL_FORT_FAB_ARG_3D(der),const int* nvar,
                    const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
                    const int* lo, const int* hi,
                    const int* domain_lo, const int* domain_hi,
                    const amrex::Real* delta, const amrex::Real* xlo,
                    const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
                    const int* level, const int* grid_no) ;

     void derdvrho (BL_FORT_FAB_ARG_3D(der),const int* nvar,
                    const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
                    const int* lo, const int* hi,
                    const int* domain_lo, const int* domain_hi,
                    const amrex::Real* delta, const amrex::Real* xlo,
                    const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
                    const int* level, const int* grid_no) ;

     void derlgrhodust (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                        const int* nvar, const amrex::Real* compdat,
                        ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
                        const int* lo, const int* hi,
                        const int* domain_lo, const int* domain_hi,
                        const amrex::Real* delta, const amrex::Real* xlo,
                        const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                        const int* level, const int* grid_no) ;

     void dermgdivu (BL_FORT_FAB_ARG_ANYD(der),const int* nvar,
                     const BL_FORT_FAB_ARG_ANYD(data),const int* ncomp,
                     AMREX_ARLIM_P(vlo), AMREX_ARLIM_P(vhi),
                     AMREX_ARLIM_P(dlo), AMREX_ARLIM_P(dhi),
                     const amrex::Real* delta, const amrex::Real* xlo,
                     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
                     const int* level, const int* grid_no) ;

     void deravgpres (BL_FORT_FAB_ARG_ANYD(der),const int* nvar,
                      const BL_FORT_FAB_ARG_ANYD(data),const int* ncomp,
                      AMREX_ARLIM_P(vlo), AMREX_ARLIM_P(vhi),
                      AMREX_ARLIM_P(dlo), AMREX_ARLIM_P(dhi),
                      const amrex::Real* delta, const amrex::Real* xlo,
                      const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
                      const int* level, const int* grid_no) ;

     void dergrdp (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                   const int* nvar, const amrex::Real* compdat,
                   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
                   const int* lo, const int* hi,
                   const int* domain_lo, const int* domain_hi,
                   const amrex::Real* delta, const amrex::Real* xlo,
                   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                   const int* level, const int* grid_no) ;


#if (BL_SPACEDIM == 3)
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

     void derpresvars(amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

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
