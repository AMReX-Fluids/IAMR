#ifndef _DERIVE_F_H
#define _DERIVE_F_H

/*
** $Id: DERIVE_F.H,v 1.15 2006-09-08 21:21:02 almgren Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_DERGRDPX    dergrdpx
#    define FORT_GRADP_DIR   gradpdir
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_DERGRDPX    DERGRDPX
#    define FORT_GRADP_DIR   GRADPDIR
#  elif defined (BL_FORT_USE_LOWERCASE)
#    define FORT_DERGRDPX    dergrdpx
#    define FORT_GRADP_DIR   gradpdir
#  elif defined (BL_FORT_USE_UNDERSCORE)
#    define FORT_DERGRDPX    dergrdpx_
#    define FORT_GRADP_DIR   gradpdir_
#  endif

#include <ArrayLim.H>

extern "C" {

     void FORT_DERGRDPX   (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
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

}
#endif
#endif
