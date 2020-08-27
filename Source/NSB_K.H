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
edge_interp_k(amrex::Box const& bx,
              int                                   ncomp,
              int                                   dir,
              amrex::GpuArray<int,AMREX_SPACEDIM>   const& ratio,
              amrex::Array4<amrex::Real>            const& fine) noexcept
{
   using namespace amrex::literals;

   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);
   int idlo[3] = {D_DECL(lo.x,lo.y,lo.z)};
   int idhi[3] = {D_DECL(hi.x,hi.y,hi.z)};

   // Get the "other" dimension(s)
#if ( AMREX_SPACEDIM == 2 )
   int dim1 = ( dir == 0 ) ? 1 : 0;
   int d1lo = idlo[dim1];
   int d1hi = idhi[dim1];
#elif ( AMREX_SPACEDIM == 3 )
   int dim1 = ( dir != 0 ) ? 0 : 1;
   int dim2 = ( dir != 0 ) ? ( ( dir == 2 ) ? 1 : 2 ) : 2;
   int d1lo = idlo[dim1];
   int d1hi = idhi[dim1];
   int d2lo = idlo[dim2];
   int d2hi = idhi[dim2];
#endif

   int idx[3] = {0};
   int idx_1[3] = {0};
   int idx_2[3] = {0};

#if ( AMREX_SPACEDIM == 2 )
   for (int comp = 0; comp < ncomp; comp++) {
      for (int d1 = d1lo; d1 <= d1hi; d1 += ratio[dim1]) {
         idx_1[dim1] = d1;
         idx_2[dim1] = d1;
         for (int d = idlo[dir]; d <= idhi[dir] - ratio[dir]; d += ratio[dir]) {
            idx_1[dir] = d + ratio[dir];
            idx_2[dir] = d;
            amrex::Real df = fine(idx_1[0],idx_1[1],idx_1[2],comp) - fine(idx_2[0],idx_2[1],idx_2[2],comp);
            for (int nd = 1; nd < ratio[dir]; nd++) {
               amrex::Real val = fine(idx_2[0],idx_2[1],idx_2[2],comp) + df * nd / ratio[dir];
               for (int nd1 = d1; nd1 < d1+ratio[dim1]; nd1++) {
                  idx[dir] = d + nd;
                  idx[dim1] = nd1;
                  fine(idx[0],idx[1],idx[2],comp) = val;
               }
            }
         }
      }
   }
#elif ( AMREX_SPACEDIM == 3 )
   for (int comp = 0; comp < ncomp; comp++) {
      for (int d2 = d2lo; d2 <= d2hi; d2 += ratio[dim2]) {
         idx_1[dim2] = d2;
         idx_2[dim2] = d2;
         for (int d1 = d1lo; d1 <= d1hi; d1 += ratio[dim1]) {
            idx_1[dim1] = d1;
            idx_2[dim1] = d1;
            for (int d = idlo[dir]; d <= idhi[dir] - ratio[dir]; d += ratio[dir]) {
               idx_1[dir] = d + ratio[dir];
               idx_2[dir] = d;
               amrex::Real df = fine(idx_1[0],idx_1[1],idx_1[2],comp) - fine(idx_2[0],idx_2[1],idx_2[2],comp);
               for (int nd = 1; nd < ratio[dir]; nd++) {
                  amrex::Real val = fine(idx_2[0],idx_2[1],idx_2[2],comp) + df * nd / ratio[dir];
                  for (int nd1 = d1; nd1 < d1+ratio[dim1]; nd1++) {
                     for (int nd2 = d2; nd2 < d2+ratio[dim2]; nd2++) {
                        idx[dir] = d + nd;
                        idx[dim1] = nd1;
                        idx[dim2] = nd2;
                        fine(idx[0],idx[1],idx[2],comp) = val;
                     }
                  }
               }
            }
         }
      }
   }
#endif
}
#endif
