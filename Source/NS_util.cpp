#include <NS_util.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

/*
  Useful functions to get min, max and maxabs vectors over a range of components in a vector of MultiFabs,
  and to do so with a single reduction call. Note: these are extremely messy.  There is probably a nicer
  way to do this...maybe someone will volunteer to make these better?
 */
  
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
  Vector<Vector<int>> first(nComp,Vector<int>(nThreads,1)); // A flag per thread/comp to avoid initializing with an arbitrary number

  // Write each thread value into array over [comp][thread]
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    int threadID = omp_get_thread_num();
    
    for (int i=0; i<mfs.size(); ++i) {

      const MultiFab& mf = *mfs[i];
    
      AMREX_ALWAYS_ASSERT(mf.nComp() >= sComp + nComp);
      AMREX_ALWAYS_ASSERT(mf.nGrow() >= nGrow);

      for (MFIter mfi(mf,tilesize); mfi.isValid(); ++mfi) {

        const Box& box = mfi.growntilebox(nGrow);
        
        for (int n=0; n<nComp; ++n) {

          auto my_max = mf[mfi].max(box,sComp+n);

          auto& the_max = gMax[n][threadID];

          if (first[n][threadID]==1) {
            the_max = my_max;
            first[n][threadID] = 0;
          } else {
            the_max = std::max(my_max, the_max);
          }

        }
      }
    }
  }

  // Find max of each comp over thread array
  Vector<int> gFirst(nComp,1); // A flag per comp to avoid initializing final answers with an arbitrary number
  Vector<Real> globalMax(nComp);
  for (int n=0; n<nComp; ++n) {
    for (int i=0; i<nThreads; ++i) {
      if (first[n][i] != 1) {
        if (gFirst[n] == 1) {
          globalMax[n] = gMax[n][i];
          gFirst[n] = 0;
        }
        else
        {
          globalMax[n] = std::max(globalMax[n],gMax[n][i]);
        }
      }
    }
    AMREX_ALWAYS_ASSERT(gFirst[n] != 1); // This means that we never picked up a useful number
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
  Vector<Vector<int>> first(nComp,Vector<int>(nThreads,1));

  // Write each thread value into array over [comp][thread]
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    int threadID = omp_get_thread_num();
    
    for (int i=0; i<mfs.size(); ++i) {

      const MultiFab& mf = *mfs[i];
    
      AMREX_ALWAYS_ASSERT(mf.nComp() >= sComp + nComp);
      AMREX_ALWAYS_ASSERT(mf.nGrow() >= nGrow);

      for (MFIter mfi(mf,tilesize); mfi.isValid(); ++mfi) {

        const Box& box = mfi.growntilebox(nGrow);
        
        for (int n=0; n<nComp; ++n) {

          auto my_min = mf[mfi].min(box,sComp+n);

          auto& the_min = gMin[n][threadID];

          if (first[n][threadID]==1) {
            the_min = my_min;
            first[n][threadID] = 0;
          } else {
            the_min = std::min(my_min, the_min);
          }
        }
      }
    }
  }

  // Find min of each comp over thread array
  Vector<int> gFirst(nComp,1);
  Vector<Real> globalMin(nComp);
  for (int n=0; n<nComp; ++n) {
    for (int i=0; i<nThreads; ++i) {
      if (first[n][i] != 1) {
        if (gFirst[n] == 1) {
          globalMin[n] = gMin[n][i];
          gFirst[n] = 0;
        }
        else
        {
          globalMin[n] = std::min(globalMin[n],gMin[n][i]);
        }
      }
    }
    AMREX_ALWAYS_ASSERT(gFirst[n] != 1); // This means that we never picked up a useful number
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
  Vector<Vector<int>> first(nComp,Vector<int>(nThreads,1));

  // Write each thread value into array over [comp][thread]
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    int threadID = omp_get_thread_num();
    
    for (int i=0; i<mfs.size(); ++i) {

      const MultiFab& mf = *mfs[i];
    
      AMREX_ALWAYS_ASSERT(mf.nComp() >= sComp + nComp);
      AMREX_ALWAYS_ASSERT(mf.nGrow() >= nGrow);

      for (MFIter mfi(mf,tilesize); mfi.isValid(); ++mfi) {

        const Box& box = mfi.growntilebox(nGrow);
        
        for (int n=0; n<nComp; ++n) {

          auto my_max = mf[mfi].maxabs(box,sComp+n);

          auto& the_max = gMax[n][threadID];

          if (first[n][threadID]==1) {
            the_max = my_max;
            first[n][threadID] = 0;
          } else {
            the_max = std::max(my_max, the_max);
          }

        }
      }
    }
  }

  // Find max of each comp over thread array
  Vector<int> gFirst(nComp,1);
  Vector<Real> globalMax(nComp);
  for (int n=0; n<nComp; ++n) {
    for (int i=0; i<nThreads; ++i) {
      if (first[n][i] != 1) {
        if (gFirst[n] == 1) {
          globalMax[n] = gMax[n][i];
          gFirst[n] = 0;
        }
        else
        {
          globalMax[n] = std::max(globalMax[n],gMax[n][i]);
        }
      }
    }
    AMREX_ALWAYS_ASSERT(gFirst[n] != 1); // This means that we never picked up a useful number
  }

  // MPI reduce max
  ParallelDescriptor::ReduceRealMax(&(globalMax[0]),globalMax.size());

  return globalMax;
}

}
