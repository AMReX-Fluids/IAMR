#include <NS_util.H>

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
  int nThreads = 1;
#ifdef _OPENMP
#pragma omp parallel
  nThreads = omp_get_num_threads();
#endif

  Vector<Vector<Real>> gMax(nComp,Vector<Real>(nThreads,-1));

  // Write each thread value into array over [comp][thread]
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    Vector<int> first(nComp,1);
    int threadID = omp_get_thread_num();
    
    for (int i=0; i<mfs.size(); ++i) {

      const MultiFab& mf = *mfs[i];
    
      AMREX_ALWAYS_ASSERT(mf.nComp() >= sComp + nComp);

      for (MFIter mfi(mf,tilesize); mfi.isValid(); ++mfi) {

        for (int n=0; n<nComp; ++n) {

          auto my_max = mf[mfi].max(mfi.growntilebox(nGrow),sComp+n);

          auto& the_max = gMax[n][threadID];

          if (first[n]==1) {
            the_max = my_max;
            first[n] = 0;
          } else {
            the_max = std::max(my_max, the_max);
          }

        }
      }
    }
  }

  // Find max of each comp over thread array
  Vector<Real> globalMax(nComp);
  for (int n=0; n<nComp; ++n) {

    globalMax[n] = *max_element(std::begin(gMax[n]), std::end(gMax[n]));

  }

  // MPI reduce max
  ParallelDescriptor::ReduceRealMax(&(globalMax[0]),globalMax.size());

  return globalMax;
}

Vector<Real>
VectorMin(const Vector<const MultiFab *>& mfs,
          const IntVect&                  tilesize,
          int                             sComp,
          int                             nComp,
          int                             nGrow)
{
  int nThreads = 1;
#ifdef _OPENMP
#pragma omp parallel
  nThreads = omp_get_num_threads();
#endif

  Vector<Vector<Real>> gMin(nComp,Vector<Real>(nThreads,-1));

  // Write each thread value into array over [comp][thread]
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    Vector<int> first(nComp,1);
    int threadID = omp_get_thread_num();
    
    for (int i=0; i<mfs.size(); ++i) {

      const MultiFab& mf = *mfs[i];
    
      AMREX_ALWAYS_ASSERT(mf.nComp() >= sComp + nComp);

      for (MFIter mfi(mf,tilesize); mfi.isValid(); ++mfi) {

        for (int n=0; n<nComp; ++n) {

          auto my_min = mf[mfi].min(mfi.growntilebox(nGrow),sComp+n);

          auto& the_min = gMin[n][threadID];

          if (first[n]==1) {
            the_min = my_min;
            first[n] = 0;
          } else {
            the_min = std::min(my_min, the_min);
          }
        }
      }
    }
  }

  // Find min of each comp over thread array
  Vector<Real> globalMin(nComp);
  for (int n=0; n<nComp; ++n) {

    globalMin[n] = *min_element(std::begin(gMin[n]), std::end(gMin[n]));

  }

  // MPI reduce min
  ParallelDescriptor::ReduceRealMin(&(globalMin[0]),globalMin.size());

  return globalMin;
}

Vector<Real>
VectorMaxAbs(const Vector<const MultiFab *>& mfs,
             const IntVect&                  tilesize,
             int                             sComp,
             int                             nComp,
             int                             nGrow)
{
  int nThreads = 1;
#ifdef _OPENMP
#pragma omp parallel
  nThreads = omp_get_num_threads();
#endif

  Vector<Vector<Real>> gMax(nComp,Vector<Real>(nThreads,-1));

  // Write each thread value into array over [comp][thread]
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    Vector<int> first(nComp,1);
    int threadID = omp_get_thread_num();
    
    for (int i=0; i<mfs.size(); ++i) {

      const MultiFab& mf = *mfs[i];
    
      AMREX_ALWAYS_ASSERT(mf.nComp() >= sComp + nComp);

      for (MFIter mfi(mf,tilesize); mfi.isValid(); ++mfi) {

        for (int n=0; n<nComp; ++n) {

          auto my_max = mf[mfi].maxabs(mfi.growntilebox(nGrow),sComp+n);

          auto& the_max = gMax[n][threadID];

          if (first[n]==1) {
            the_max = my_max;
            first[n] = 0;
          } else {
            the_max = std::max(my_max, the_max);
          }

        }
      }
    }
  }

  // Find max of each comp over thread array
  Vector<Real> globalMax(nComp,-1);
  for (int n=0; n<nComp; ++n) {

    const auto& v = gMax[n];

    for (int i=0; i<v.size(); ++i) {

      auto my_max = std::abs(v[i]);

      auto& the_max = globalMax[n];
      
      the_max = the_max < 0 ? my_max : max(my_max, the_max);

    }
  }

  // MPI reduce max
  ParallelDescriptor::ReduceRealMax(&(globalMax[0]),globalMax.size());

  return globalMax;
}

}
