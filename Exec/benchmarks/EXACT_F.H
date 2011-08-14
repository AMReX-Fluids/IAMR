#ifndef _FORT_INTERFACE_
#define _FORT_INTERFACE_

#ifdef BL_LANG_FORT
c
c  -------------------------
c  ::: Fortran Interface :::
c  -------------------------
c
#  define   FORT_VISCBENCH        viscbench
#else
//
// ---------------------------------------
// ::: C Interface to Fortran Routines :::
// ---------------------------------------
//

#ifdef BL_FORT_USE_UNDERSCORE
#  define   FORT_VISCBENCH        viscbench_
#else
#  define   FORT_VISCBENCH        VISCBENCH
#endif

//
// Use ARLIM in the following.
//

extern "C"
{
#if (BL_SPACEDIM==3)
  void FORT_VISCBENCH(const Real *time, const Real *nu, const int *unifdir,
                      const int *lo, const int *hi, 
                      const int *ncomp,
                      const Real *state, ARLIM_P(lo), ARLIM_P(hi),
                      const Real *delta,
                      const Real *xlo, const Real *xhi);
#else
  void FORT_VISCBENCH(const Real *time, const Real *nu, 
                      const int *lo, const int *hi, 
                      const int *ncomp,
                      const Real *state, ARLIM_P(lo), ARLIM_P(hi),
                      const Real *delta,
                      const Real *xlo, const Real *xhi);
#endif

}
#endif
#endif
