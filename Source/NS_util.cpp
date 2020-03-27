#include <NS_util.H>

#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

  Vector<Real>
  VectorMax(const Vector<const MultiFab *>& mfs,
            const IntVect&                  tilesize,
            int                             sComp,
            int                             nComp,
            int                             nGrow)
  {
    Vector<Real> gMax(nComp);

    for (int n=0; n<nComp; ++n) {
        gMax[n] = mfs[0]->max(sComp+n,nGrow);
    }

    for (int i=1; i<mfs.size(); ++i) {

      for (int n=0; n<nComp; ++n) {
        gMax[n] = std::max(gMax[n], mfs[i]->max(sComp+n,nGrow));
      }
    }

    return gMax;
  }

  Vector<Real>
  VectorMaxAbs(const Vector<const MultiFab *>& mfs,
               const IntVect&                  tilesize,
               int                             sComp,
               int                             nComp,
               int                             nGrow)
  {
    Vector<Real> gMaxAbs(nComp);

    for (int n=0; n<nComp; ++n) {
        gMaxAbs[n] = mfs[0]->norm0(sComp+n,nGrow,false, true);
    }

    for (int i=1; i<mfs.size(); ++i) {

      for (int n=0; n<nComp; ++n) {
          gMaxAbs[n] = std::max(gMaxAbs[n], mfs[i]->norm0(sComp+n,nGrow,false,true));
      }
    }

    return gMaxAbs;
  }

  Vector<Real>
  VectorMin(const Vector<const MultiFab *>& mfs,
            const IntVect&                  tilesize,
            int                             sComp,
            int                             nComp,
            int                             nGrow)
  {
    Vector<Real> gMin(nComp);

    for (int n=0; n<nComp; ++n) {
      gMin[n] = mfs[0]->min(sComp+n,nGrow);
    }

    for (int i=1; i<mfs.size(); ++i) {

      for (int n=0; n<nComp; ++n) {
        gMin[n] = std::min(gMin[n], mfs[i]->min(sComp+n,nGrow));
      }
    }

    return gMin;
  }

}
