#ifndef NSB_K_H_
#define NSB_K_H_

#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE
inline
void
edg2cen_average(int i, int j, int k,
                AMREX_D_DECL(amrex::Array4<const amrex::Real> const& efab_x,
                             amrex::Array4<const amrex::Real> const& efab_y,
                             amrex::Array4<const amrex::Real> const& efab_z),
                amrex::Array4<      amrex::Real> const& cfab) noexcept
{
   using namespace amrex::literals;

   AMREX_D_TERM(cfab(i,j,k,0) = 0.5_rt * ( efab_x(i,j,k) + efab_x(i+1,j,k) );,
                cfab(i,j,k,1) = 0.5_rt * ( efab_y(i,j,k) + efab_y(i,j+1,k) );,
                cfab(i,j,k,2) = 0.5_rt * ( efab_z(i,j,k) + efab_z(i,j,k+1) ););

}

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

AMREX_GPU_HOST_DEVICE
inline
void
calc_gradp(int i, int j, int k,
           amrex::Real                                  const space_scale,
           amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const dxinv,
           amrex::Array4<const amrex::Real>             const& p,
           amrex::Array4<amrex::Real>                   const& gp) noexcept
{
   gp(i,j,k,0) = space_scale * dxinv[0] * ( AMREX_D_TERM(  p(i+1,j  ,k  )-p(i  ,j  ,k  ),
                                                         + p(i+1,j+1,k  )-p(i  ,j+1,k  ),
                                                         + p(i+1,j  ,k+1)-p(i  ,j  ,k+1)
                                                         + p(i+1,j+1,k+1)-p(i  ,j+1,k+1)));
   gp(i,j,k,1) = space_scale * dxinv[1] * ( AMREX_D_TERM(  p(i  ,j+1,k  )-p(i  ,j  ,k  ),
                                                         + p(i+1,j+1,k  )-p(i+1,j  ,k  ),
                                                         + p(i  ,j+1,k+1)-p(i  ,j  ,k+1)
                                                         + p(i+1,j+1,k+1)-p(i+1,j  ,k+1)));
#if ( AMREX_SPACEDIM == 3 )
   gp(i,j,k,2) = space_scale * dxinv[2] *             (  p(i  ,j  ,k+1)-p(i  ,j  ,k  )
                                                       + p(i+1,j  ,k+1)-p(i+1,j  ,k  )
                                                       + p(i  ,j+1,k+1)-p(i  ,j+1,k  )
                                                       + p(i+1,j+1,k+1)-p(i+1,j+1,k  ));
#endif
}

AMREX_GPU_HOST_DEVICE
inline
void
srcrsereg_k(amrex::Box const& bx,
            amrex::Box const& fbx,
            int                                   dir,
            amrex::GpuArray<int,AMREX_SPACEDIM>   const &ratio,
            amrex::Array4<amrex::Real>            const& fine,
            amrex::Array4<amrex::Real>            const& crse) noexcept
{
   using namespace amrex::literals;

   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);
   const auto flo = amrex::lbound(fbx);
   const auto fhi = amrex::ubound(fbx);
   int idlo[3] = {D_DECL(lo.x,lo.y,lo.z)};
   int idhi[3] = {D_DECL(hi.x,hi.y,hi.z)};

   // Get the "other" dimension(s)
#if ( AMREX_SPACEDIM == 2 )
   int dim1 = ( dir == 0 ) ? 1 : 0;
   int d1lo = idlo[dim1];
   int d1hi = idhi[dim1];
#elif ( AMREX_SPACEDIM == 3 )
   const amrex::Real fourthirds = 4.0_rt/3.0_rt;
   const amrex::Real threefourths = 3.0_rt/4.0_rt;
   int idflo[3] = {D_DECL(flo.x,flo.y,flo.z)};
   int idfhi[3] = {D_DECL(fhi.x,fhi.y,fhi.z)};
   int dim1 = ( dir != 0 ) ? 0 : 1;
   int dim2 = ( dir != 0 ) ? ( ( dir == 2 ) ? 1 : 2 ) : 2;
   int d1lo = idlo[dim1];
   int d1hi = idhi[dim1];
   int d2lo = idlo[dim2];
   int d2hi = idhi[dim2];
   int d1lo_f = idflo[dim1];
   int d1hi_f = idfhi[dim1];
   int d2lo_f = idflo[dim2];
   int d2hi_f = idfhi[dim2];
#endif

   // Index array
   int idx[3] = {0};
   
   // 2D : Set crse lo box edge in dir to 0.0
   // 3D : Set crse lo box face in dir to 0.0
   int dc = idlo[dir];
   idx[dir] = dc;
   for (int d1c = d1lo; d1c <= d1hi; ++d1c) {
#if ( AMREX_SPACEDIM == 2 )
      idx[dim1] = d1c;
      crse(idx[0],idx[1],idx[2]) = 0.0_rt;
#elif ( AMREX_SPACEDIM == 3 )
      for (int d2c = d2lo; d2c <= d2hi; ++d2c) {
         idx[dim1] = d1c;
         idx[dim2] = d2c;
         crse(idx[0],idx[1],idx[2]) = 0.0_rt;
      }
#endif
   }

   int df  = dc * ratio[dir];
   idx[dir] = df;

#if ( AMREX_SPACEDIM == 2 )
   // Scale fine by 1/2 on corners
   idx[dim1] = ratio[dim1] * d1lo;
   fine(idx[0],idx[1],idx[2]) *= 0.5_rt;
   idx[dim1] = ratio[dim1] * d1hi;
   fine(idx[0],idx[1],idx[2]) *= 0.5_rt;

#elif ( AMREX_SPACEDIM == 3 )
   // Scale fine by 1/2 on edges
   idx[dim1] = ratio[dim1] * d1lo;
   for (int d2f = d2lo_f; d2f <= d2hi_f; ++d2f) {
      idx[dim2] = d2f;
      fine(idx[0],idx[1],idx[2]) *= 0.5_rt;
   }
   idx[dim1] = ratio[dim1] * d1hi;
   for (int d2f = d2lo_f; d2f <= d2hi_f; ++d2f) {
      idx[dim2] = d2f;
      fine(idx[0],idx[1],idx[2]) *= 0.5_rt;
   }
   idx[dim2] = ratio[dim2] * d2lo;
   for (int d1f = d1lo_f; d1f <= d1hi_f; ++d1f) {
      idx[dim1] = d1f;
      fine(idx[0],idx[1],idx[2]) *= 0.5_rt;
   }
   idx[dim2] = ratio[dim2] * d2hi;
   for (int d1f = d1lo_f; d1f <= d1hi_f; ++d1f) {
      idx[dim1] = d1f;
      fine(idx[0],idx[1],idx[2]) *= 0.5_rt;
   }

   // Scale fine by 1/3 on corners
   idx[dim1] = ratio[dim1] * d1lo;
   idx[dim2] = ratio[dim2] * d2lo;
   fine(idx[0],idx[1],idx[2]) *= fourthirds;
   idx[dim1] = ratio[dim1] * d1hi;
   idx[dim2] = ratio[dim2] * d2lo;
   fine(idx[0],idx[1],idx[2]) *= fourthirds;
   idx[dim1] = ratio[dim1] * d1lo;
   idx[dim2] = ratio[dim2] * d2hi;
   fine(idx[0],idx[1],idx[2]) *= fourthirds;
   idx[dim1] = ratio[dim1] * d1hi;
   idx[dim2] = ratio[dim2] * d2hi;
   fine(idx[0],idx[1],idx[2]) *= fourthirds;
#endif

   // Fine index array
   int idx_f[3] = {0};
   int idx_f2[3] = {0};

   idx[dir] = dc;
   idx_f[dir] = dc * ratio[dir];
   idx_f2[dir] = idx_f[dir];

#if ( AMREX_SPACEDIM == 2 )
   amrex::Real denom = 1.0_rt / ( ratio[dir] * ratio[dim1] * ratio[dim1] );
   for (int nd1 = 0; nd1 < ratio[dim1]; nd1++) {
      amrex::Real coeff = (ratio[dim1] - nd1) * denom;
      if (nd1 == 0) coeff *= 0.5_rt;
      for (int d1c = d1lo; d1c <= d1hi; ++d1c) {
         int d1f = ratio[dim1] * d1c;
         idx[dim1] = d1c;
         idx_f[dim1] = d1f + nd1;
         idx_f2[dim1] = d1f - nd1;
         crse(idx[0],idx[1],idx[2]) += coeff * (fine(idx_f[0],idx_f[1],idx_f[2]) + fine(idx_f2[0],idx_f2[1],idx_f2[2]) );
      }
   }
#elif ( AMREX_SPACEDIM == 3 )
   int idx_f3[3] = {0};
   idx_f3[dir] = idx_f[dir];
   int idx_f4[3] = {0};
   idx_f4[dir] = idx_f[dir];
   amrex::Real denom = 1.0_rt / ( ratio[dir] * ratio[dim1] * ratio[dim1]
                                             * ratio[dim2] * ratio[dim2] );
   for (int nd2 = 0; nd2 < ratio[dim2]; nd2++) {
      for (int nd1 = 0; nd1 < ratio[dim1]; nd1++) {
         amrex::Real coeff = (ratio[dim1] - nd1) * (ratio[dim2] - nd2) * denom;
         if (nd2 == 0) coeff *= 0.5_rt;
         if (nd1 == 0) coeff *= 0.5_rt;
         for (int d2c = d2lo; d2c <= d2hi; ++d2c) {
            int d2f = ratio[dim2] * d2c;
            idx[dim2] = d2c;
            for (int d1c = d1lo; d1c <= d1hi; ++d1c) {
               int d1f = ratio[dim1] * d1c;
               idx[dim1] = d1c;
               idx_f[dim1] = d1f + nd1;
               idx_f[dim2] = d2f + nd2;
               idx_f2[dim1] = d1f + nd1;
               idx_f2[dim2] = d2f - nd2;
               idx_f3[dim1] = d1f - nd1;
               idx_f3[dim2] = d2f + nd2;
               idx_f4[dim1] = d1f - nd1;
               idx_f4[dim2] = d2f - nd2;
               crse(idx[0],idx[1],idx[2]) += coeff * (  fine(idx_f[0],idx_f[1],idx_f[2]) + fine(idx_f3[0],idx_f3[1],idx_f3[2])
                                                      + fine(idx_f2[0],idx_f2[1],idx_f2[2]) + fine(idx_f4[0],idx_f4[1],idx_f4[2]) );
            }
         }
      }
   }
#endif

   idx[dir] = dc * ratio[dir];

#if ( AMREX_SPACEDIM == 2 )
   // Scale back fine by 2 on corners
   idx[dim1] = ratio[dim1] * d1lo;
   fine(idx[0],idx[1],idx[2]) *= 2.0_rt;
   idx[dim1] = ratio[dim1] * d1hi;
   fine(idx[0],idx[1],idx[2]) *= 2.0_rt;

#elif ( AMREX_SPACEDIM == 3 )
   // Scale back fine by 2 on edges
   idx[dim1] = ratio[dim1] * d1lo;
   for (int d2f = d2lo_f; d2f <= d2hi_f; ++d2f) {
      idx[dim2] = d2f;
      fine(idx[0],idx[1],idx[2]) *= 2.0_rt;
   }
   idx[dim1] = ratio[dim1] * d1hi;
   for (int d2f = d2lo_f; d2f <= d2hi_f; ++d2f) {
      idx[dim2] = d2f;
      fine(idx[0],idx[1],idx[2]) *= 2.0_rt;
   }
   idx[dim2] = ratio[dim2] * d2lo;
   for (int d1f = d1lo_f; d1f <= d1hi_f; ++d1f) {
      idx[dim1] = d1f;
      fine(idx[0],idx[1],idx[2]) *= 2.0_rt;
   }
   idx[dim2] = ratio[dim2] * d2hi;
   for (int d1f = d1lo_f; d1f <= d1hi_f; ++d1f) {
      idx[dim1] = d1f;
      fine(idx[0],idx[1],idx[2]) *= 2.0_rt;
   }

   // Scale back fine by 1/3 on corners
   idx[dim1] = ratio[dim1] * d1lo;
   idx[dim2] = ratio[dim2] * d2lo;
   fine(idx[0],idx[1],idx[2]) *= threefourths;
   idx[dim1] = ratio[dim1] * d1hi;
   idx[dim2] = ratio[dim2] * d2lo;
   fine(idx[0],idx[1],idx[2]) *= threefourths;
   idx[dim1] = ratio[dim1] * d1lo;
   idx[dim2] = ratio[dim2] * d2hi;
   fine(idx[0],idx[1],idx[2]) *= threefourths;
   idx[dim1] = ratio[dim1] * d1hi;
   idx[dim2] = ratio[dim2] * d2hi;
   fine(idx[0],idx[1],idx[2]) *= threefourths;
#endif
}
#endif
