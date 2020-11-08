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
//
// Do linear in dir, pc transverse to dir, leave alone the fine values
// lining up with coarse edges--assume these have been set to hold the
// values you want to interpolate to the rest.
//
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_interp_x (int i, int j, int k, int n, amrex::Array4<amrex::Real> const& fine,
	       amrex::GpuArray<int,AMREX_SPACEDIM>   const& ratio) noexcept
{
    int rmod = i%ratio[0];

    if ( rmod != 0 )
    {
      int ii = amrex::coarsen(i,ratio[0]);
      amrex::Real const w = static_cast<amrex::Real>(i-ii*ratio[0]) * (amrex::Real(1.)/ratio[0]);
      rmod = std::abs(rmod);
      fine(i,j,k,n) = (amrex::Real(1.)-w) * fine(i-rmod,j,k,n) + w * fine(i+(ratio[0]-rmod),j,k,n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_interp_y (int i, int j, int k, int n, amrex::Array4<amrex::Real> const& fine,
	       amrex::GpuArray<int,AMREX_SPACEDIM>   const& ratio) noexcept
{
    int rmod = j%ratio[1];

    if ( rmod != 0 )
    {
      int jj = amrex::coarsen(j,ratio[1]);
      amrex::Real const w = static_cast<amrex::Real>(j-jj*ratio[1]) * (amrex::Real(1.)/ratio[1]);
      rmod = std::abs(rmod);
      fine(i,j,k,n) = (amrex::Real(1.)-w) * fine(i,j-rmod,k,n) + w * fine(i,j+(ratio[1]-rmod),k,n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_interp_z (int i, int j, int k, int n, amrex::Array4<amrex::Real> const& fine,
	       amrex::GpuArray<int,AMREX_SPACEDIM>   const& ratio) noexcept
{
    int rmod = k%ratio[2];

    if ( rmod != 0 )
    {
      int kk = amrex::coarsen(k,ratio[2]);
      amrex::Real const w = static_cast<amrex::Real>(k-kk*ratio[2]) * (amrex::Real(1.)/ratio[2]);
      rmod = std::abs(rmod);
      fine(i,j,k,n) = (amrex::Real(1.)-w) * fine(i,j,k-rmod,n) + w * fine(i,j,k+(ratio[2]-rmod),n);
    }
}
#endif
