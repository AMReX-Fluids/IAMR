#ifndef _FORT_INTERFACE_
#define _FORT_INTERFACE_

#ifdef BL_LANG_FORT
c
c  -------------------------
c  ::: Fortran Interface :::
c  -------------------------
c
#  define   FORT_VISCBENCH        viscbench
#endif



#ifdef BL_LANG_CC
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
  void FORT_VISCBENCH(const Real *time, const Real *nu,
                      const int *lo, const int *hi, 
                      const int *ncomp,
                      const Real *state, ARLIM_P(lo), ARLIM_P(hi),
                      const Real *delta,
                      const Real *xlo, const Real *xhi);

}
#endif
#endif
