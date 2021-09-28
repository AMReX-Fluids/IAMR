#ifndef NSB_K_H_
#define NSB_K_H_

#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <cmath>


AMREX_GPU_HOST_DEVICE
inline
void
set_body_state_k(int i, int j, int k,
                 int ncomp,
                 amrex::Real                   const b [],
                 int                           bval,
                 amrex::Array4<int>            const& mask,
                 amrex::Array4<amrex::Real>    const& state) noexcept
{
   using namespace amrex::literals;

   if ( mask(i,j,k) == bval ) {
      for (int n = 0; n < ncomp; n++) {
         state(i,j,k,n) = b[n];
      }
   }

}

#endif
