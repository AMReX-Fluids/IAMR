#ifndef NSB_K_H_
#define NSB_K_H_

#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE
inline
void
edg2cen_average(int i, int j, int k,
                D_DECL(amrex::Array4<const amrex::Real> const& efab_x,
                       amrex::Array4<const amrex::Real> const& efab_y,
                       amrex::Array4<const amrex::Real> const& efab_z),
                amrex::Array4<      amrex::Real> const& cfab) noexcept
{      
   using namespace amrex::literals;
   
   D_TERM(cfab(i,j,k,0) = 0.5_rt * ( efab_x(i,j,k) + efab_x(i+1,j,k) );,
          cfab(i,j,k,1) = 0.5_rt * ( efab_y(i,j,k) + efab_y(i,j+1,k) );,
          cfab(i,j,k,2) = 0.5_rt * ( efab_z(i,j,k) + efab_z(i,j,k+1) ););
   
}
#endif
