#include <NS_util.H>

#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

/*
  Useful functions to get min, max and maxabs vectors over a range of
  components in a vector of MultiFabs, and to do so with a single reduction call.
  Note: these are extremely messy.  There is probably a nicer way to do
  this...maybe someone will volunteer to make these better?
*/

// User-defined reductions in OpenMP:
//
// MyType:
// The non-native type which will be reduced using an OpenMP user-defined
// reduction (UDR).
struct MyType {
  amrex::Real val;
};

// The 'combiner' portion of the UDR. The first arg contains the result from an
// arbitrary number of reductions of the array which have already occured; this
// corresponds to the special variable 'omp_out' which is defined by OpenMP.
// The second arg is a new element to combine in the reduction; this
// corresponds to the special variable 'omp_in' which is defined by OpenMP.
void MyType_max_func(MyType *reduced_t, MyType *t_new) {
  reduced_t->val = reduced_t->val > t_new->val ? reduced_t->val : t_new->val;
  return;
}

// The 'initializer' portion of the UDR. This is executed prior to the
// 'combiner'. Since we are computing a max, and since we have chosen that all
// of the values are positive, we start with a negative number to guarantee
// that this will be less than any value in the array being reduced.
void MyType_max_init(MyType *t) {
  t->val = std::numeric_limits<amrex::Real>::lowest();
  return;
}

// This is the UDR 'declaration'. Note the usage of the three special OpenMP
// variables 'omp_out', 'omp_in', and 'omp_priv'.
#pragma omp declare reduction(my_max_func: MyType:                      \
                              MyType_max_func(&omp_out,&omp_in))        \
                              initializer(MyType_max_init(&omp_priv))

// Define another combiner/initializer for a UDR for computing min.
void MyType_min_func(MyType *reduced_t, MyType *t_new) {
  reduced_t->val = reduced_t->val < t_new->val ? reduced_t->val : t_new->val;
  return;
}

void MyType_min_init(MyType *t) {
  t->val = std::numeric_limits<amrex::Real>::max();
  return;
}

#pragma omp declare reduction(my_min_func: MyType:                      \
                              MyType_min_func(&omp_out,&omp_in))        \
                              initializer(MyType_min_init(&omp_priv))

namespace amrex {




  Vector<Real>
  VectorMax(const Vector<const MultiFab *>& mfs,
            const IntVect&                  tilesize,
            int                             sComp,
            int                             nComp,
            int                             nGrow)
  {
    Vector<Real> gMax(nComp);

    for (int i=0; i<mfs.size(); ++i) {

      const MultiFab& mf = *mfs[i];
    
      AMREX_ALWAYS_ASSERT(mf.nComp() >= sComp + nComp);

      for (int n=0; n<nComp; ++n) {

        MyType the_max;
        the_max.val = std::numeric_limits<Real>::lowest();
#ifdef _OPENMP
#pragma omp parallel reduction(my_max_func: the_max)
#endif
        for (MFIter mfi(mf,tilesize); mfi.isValid(); ++mfi) {

          const Box& bx = mfi.growntilebox(nGrow);
        
          auto my_max = mf[mfi].max(bx,sComp+n);

          the_max.val = std::max(my_max, the_max.val);
        }

        gMax[n] = the_max.val;
      }
    }

    // MPI reduce max
    ParallelDescriptor::ReduceRealMax(&(gMax[0]),nComp);

    return gMax;
  }

  Vector<Real>
  VectorMaxAbs(const Vector<const MultiFab *>& mfs,
               const IntVect&                  tilesize,
               int                             sComp,
               int                             nComp,
               int                             nGrow)
  {
    Vector<Real> gMaxabs(nComp);

    for (int i=0; i<mfs.size(); ++i) {

      const MultiFab& mf = *mfs[i];
    
      AMREX_ALWAYS_ASSERT(mf.nComp() >= sComp + nComp);

      for (int n=0; n<nComp; ++n) {

        MyType the_maxabs;
        the_maxabs.val = std::numeric_limits<Real>::lowest();
#ifdef _OPENMP
#pragma omp parallel reduction(my_max_func: the_maxabs)
#endif
        for (MFIter mfi(mf,tilesize); mfi.isValid(); ++mfi) {

          const Box& bx = mfi.growntilebox(nGrow);
        
          auto my_maxabs = mf[mfi].maxabs(bx,sComp+n);

          the_maxabs.val = std::max(my_maxabs, the_maxabs.val);
        }

        gMaxabs[n] = the_maxabs.val;
      }
    }

    // MPI reduce maxabs
    ParallelDescriptor::ReduceRealMax(&(gMaxabs[0]),nComp);

    return gMaxabs;
  }



  Vector<Real>
  VectorMin(const Vector<const MultiFab *>& mfs,
            const IntVect&                  tilesize,
            int                             sComp,
            int                             nComp,
            int                             nGrow)
  {
    Vector<Real> gMin(nComp);

    for (int i=0; i<mfs.size(); ++i) {

      const MultiFab& mf = *mfs[i];
    
      AMREX_ALWAYS_ASSERT(mf.nComp() >= sComp + nComp);

      for (int n=0; n<nComp; ++n) {

        MyType the_min;
        the_min.val = std::numeric_limits<Real>::max();
#ifdef _OPENMP
#pragma omp parallel reduction(my_min_func: the_min)
#endif
        for (MFIter mfi(mf,tilesize); mfi.isValid(); ++mfi) {

          const Box& bx = mfi.growntilebox(nGrow);
        
          auto my_min = mf[mfi].min(bx,sComp+n);

          the_min.val = std::min(my_min, the_min.val);
        }

        gMin[n] = the_min.val;
      }
    }

    // MPI reduce min
    ParallelDescriptor::ReduceRealMin(&(gMin[0]),nComp);

    return gMin;
  }
}
