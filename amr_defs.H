#ifndef _AMR_DEFS_H_
#define _AMR_DEFS_H_

//
// $Id: amr_defs.H,v 1.5 1997-10-08 20:15:48 car Exp $
//

#if 0
#if (defined BL_IRIX || defined BL_OSF1 || \
     defined BL_Linux || defined BL_ARCH_CRAY)
#  ifndef SYSV
#    define SYSV
#  endif
#endif
#endif

#include <MultiFab.H>

//#undef SPACEDIM

//#define BOX_NODE IndexType::NODE
//#define BOX_CELL IndexType::CELL

//#define zerovect   IntVect::TheZeroVector()
//#define unitvect   IntVect::TheUnitVector()
//#define nodevect   IntVect::TheNodeVector()
//#define cellvect   IntVect::TheCellVector()

#if (BL_SPACEDIM == 1)
#  define intRS const int&
#  define intS  const int&, const int&
#  define dimlist(a) a.smallEnd(0), a.bigEnd(0)
#  define RealRS const Real&
#  define RealPS Real*
#elif (BL_SPACEDIM == 2)
#  define intRS const int&, const int&
#  define intS  const int&, const int&, const int&, const int&
#  define dimlist(a) a.smallEnd(0), a.bigEnd(0), a.smallEnd(1), a.bigEnd(1)
#  define RealRS const Real&, const Real&
#  define RealPS Real*, Real*
#elif (BL_SPACEDIM == 3)
#  define intRS const int&, const int&, const int&
#  define intS  const int&, const int&, const int&, const int&, \
                const int&, const int&
#  define dimlist(a) a.smallEnd(0), a.bigEnd(0), a.smallEnd(1), a.bigEnd(1), \
                     a.smallEnd(2), a.bigEnd(2)
#  define RealRS const Real&, const Real&, const Real&
#  define RealPS Real*, Real*, Real*
#endif

// index type helper functions:

inline IntVect
type(const FArrayBox& f) 
{
  return f.box().type();
}

inline IndexType::CellIndex 
type(const FArrayBox& f, int idir) 
{
  return f.box().type(idir);
}

inline IntVect 
type(const FabArray<Real,FArrayBox>& r) 
{
  assert(r.length() > 0 && r.boxArray().ok());
  return r.box(0).type();
}

inline IndexType::CellIndex 
type(const FabArray<Real,FArrayBox>& r, int idir) 
{
  assert(r.length() > 0 && r.boxArray().ok());
  return r.box(0).type(idir);
}

// norm helper functions:

inline Real
mfnorm(const MultiFab& mf)
{
    Real r = 0;
    // PARALLEL -- REDUCTION
    // for (i = 0; i < mf.length(); i++) 
    for ( ConstMultiFabIterator mfi(mf); mfi.isValid(); ++mfi ) 
    {
	Real s = mfi().norm(0);
	r = (r > s) ? r : s;
    }
    ParallelDescriptor::ReduceRealMax(r);
    return r;
}

inline Real 
pmfnorm(const PArray<MultiFab>& pmf) 
{
  Real r = 0;
  for (int i = 0; i < pmf.length(); i++) 
  {
    Real s = mfnorm(pmf[i]);
    r = (r > s) ? r : s;
  }
  return r;
}

inline Real
pmfmin(const PArray<MultiFab>& pmf) 
{
  Real r = pmf[0].min(0);
  for (int i = 1; i < pmf.length(); i++) 
  {
    Real s = pmf[i].min(0);
    r = (r < s) ? r : s;
  }
  return r;
}

inline Real 
pmfmax(const PArray<MultiFab>& pmf) 
{
  Real r = pmf[0].max(0);
  for (int i = 1; i < pmf.length(); i++) {
    Real s = pmf[i].max(0);
    r = (r > s) ? r : s;
  }
  return r;
}

#endif
