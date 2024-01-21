#ifndef IAMR_FORT_INTERFACE_
#define IAMR_FORT_INTERFACE_

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
// Use AMREX_ARLIM in the following.
//

extern "C"
{
#if (BL_SPACEDIM==3)
  void FORT_VISCBENCH(const amrex::Real *time, const amrex::Real *nu, const int *unifdir,
                      const int *lo, const int *hi,
                      const int *ncomp,
                      const amrex::Real *state, AMREX_ARLIM_P(lo), AMREX_ARLIM_P(hi),
                      const amrex::Real *delta,
                      const amrex::Real *xlo, const amrex::Real *xhi);
#else
  void FORT_VISCBENCH(const amrex::Real *time, const amrex::Real *nu,
                      const int *lo, const int *hi,
                      const int *ncomp,
                      const amrex::Real *state, AMREX_ARLIM_P(lo), AMREX_ARLIM_P(hi),
                      const amrex::Real *delta,
                      const amrex::Real *xlo, const amrex::Real *xhi);
#endif

}
#endif
#endif
