#ifndef IAMR_EB_SLOPES_K_H_
#define IAMR_EB_SLOPES_K_H_

#include <AMReX_Array4.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_BCRec.H>


namespace {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
iamr_slopes_eb (D_DECL(int i, int j, int k), int n,
                amrex::Array4<amrex::Real const> const& state,
                amrex::Array4<amrex::Real const> const& ccent,
                amrex::Array4<amrex::EBCellFlag const> const& flag) noexcept
{
#if (AMREX_SPACEDIM==2)
    const int dim_stencil = 9;
    const int k  = 0;
    const int kk = 0;
#else
    const int dim_stencil = 27;
#endif

    amrex::Real AtA[AMREX_SPACEDIM][AMREX_SPACEDIM];
    amrex::Real Atb[AMREX_SPACEDIM];

    {
        // A and du are local to this block, so they get deleted at the
        // end of it once Atb and AtA are filled
        amrex::Real A[dim_stencil][AMREX_SPACEDIM];
        amrex::Real du[dim_stencil];

        {
            int lc = 0;

#if (AMREX_SPACEDIM==3)
            for (int kk = -1; kk<=1; kk++)
#endif
            {
                for (int jj = -1; jj<=1; jj++)
                {
                    for(int ii = -1; ii<=1; ii++)
                    {
                        if( flag(i,j,k).isConnected(ii,jj,kk) and
                            not (ii==0 and jj==0 and kk==0))
                        {
                            // Not multiplying by dx to be consistent with how the
                            // slope is stored. Also not including the global shift
                            // wrt plo or i,j,k. We only need relative distance.
                            D_TERM( A[lc][0] = ii + ccent(i+ii,j+jj,k+kk,0) - ccent(i,j,k,0);,
                                    A[lc][1] = jj + ccent(i+ii,j+jj,k+kk,1) - ccent(i,j,k,1);,
                                    A[lc][2] = kk + ccent(i+ii,j+jj,k+kk,2) - ccent(i,j,k,2););

                            du[lc] = state(i+ii,j+jj,k+kk,n) - state(i,j,k,n);
                        }
                        else
                        {
                            D_TERM( A[lc][0] = 0.0;,
                                    A[lc][1] = 0.0;,
                                    A[lc][2] = 0.0; );

                            du[lc] = 0.0;
                        }

                        lc++;
                    }
                }
            }
        }


        //
        // Fill Atb and AtA
        //
        for (int jj = 0 ; jj<AMREX_SPACEDIM; ++jj)
        {
            for(int ii = 0; ii<AMREX_SPACEDIM; ++ii)
            {
                AtA[ii][jj] = 0.0;
            }
            Atb[jj] = 0.0;
        }


        for (int lc = 0; lc < dim_stencil; ++lc)
        {
            // Matrix is symmetric: just fill (and the use) only the upper
            //half of it
            for(int jj = 0; jj<AMREX_SPACEDIM; ++jj)
            {
                for(int ii = 0; ii <= jj; ++ii)
                {
                    AtA[ii][jj] += A[lc][ii] * A[lc][jj];
                }

                Atb[jj] += A[lc][jj] * du[lc];
            }
        }

    }

    amrex::Real detAtA =
#if (AMREX_SPACEDIM==3)
        AtA[0][0]*(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[1][2]) -
        AtA[0][1]*(AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
        AtA[0][2]*(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);
#else
        AtA[0][0]*AtA[1][1] - AtA[1][0]*AtA[1][0];
#endif

    D_TERM( amrex::Real xs = 0.0;,
            amrex::Real ys = 0.0;,
            amrex::Real zs = 0.0;);

    // X direction
    if( flag(i  ,j,k).isSingleValued() or
       (flag(i-1,j,k).isSingleValued() or not flag(i,j,k).isConnected(-1,0,0)) or
       (flag(i+1,j,k).isSingleValued() or not flag(i,j,k).isConnected( 1,0,0)))
    {

        amrex::Real detAtA_x =
#if (AMREX_SPACEDIM==3)
            Atb[0]   *(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[1][2]) -
            AtA[0][1]*(Atb[1] *  AtA[2][2] - AtA[1][2]*Atb[2]   ) +
            AtA[0][2]*(Atb[1] *  AtA[2][1] - AtA[1][1]*Atb[2]   );
#else
            AtA[1][1]*Atb[0] - AtA[0][1]*Atb[1];
#endif
        // Slope at centroid of (i,j,k)
        xs = detAtA_x / detAtA;
    }
    else
    {
        amrex::Real du_xl = 2.0*(state(i  ,j,k,n) - state(i-1,j,k,n));
        amrex::Real du_xr = 2.0*(state(i+1,j,k,n) - state(i  ,j,k,n));
        amrex::Real du_xc = 0.5*(state(i+1,j,k,n) - state(i-1,j,k,n));

        amrex::Real xslope = amrex::min(std::abs(du_xl), std::abs(du_xc), std::abs(du_xr));
        xslope = (du_xr*du_xl > 0.0) ? xslope : 0.0;
        xs = (du_xc > 0.0) ? xslope : -xslope;
    }

    // Y direction
    if (flag(i,j  ,k).isSingleValued() or
       (flag(i,j-1,k).isSingleValued() or not flag(i,j,k).isConnected(0,-1,0)) or
       (flag(i,j+1,k).isSingleValued() or not flag(i,j,k).isConnected(0, 1,0)))
    {

        amrex::Real detAtA_y =
#if (AMREX_SPACEDIM==3)
            AtA[0][0]*(Atb[1]  * AtA[2][2] - AtA[1][2]*Atb[2]   ) -
            Atb[0] *  (AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
            AtA[0][2]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]);
#else
            AtA[0][0]*Atb[1] - AtA[0][1]*Atb[0];
#endif

            // Slope at centroid of (i,j,k)
            ys = detAtA_y / detAtA;

    }
    else
    {
      amrex::Real du_yl = 2.0*(state(i,j  ,k,n) - state(i,j-1,k,n));
      amrex::Real du_yr = 2.0*(state(i,j+1,k,n) - state(i,j  ,k,n));
      amrex::Real du_yc = 0.5*(state(i,j+1,k,n) - state(i,j-1,k,n));

      amrex::Real yslope = amrex::min(std::abs(du_yl), std::abs(du_yc), std::abs(du_yr));
      yslope = (du_yr*du_yl > 0.0) ? yslope : 0.0;
      ys = (du_yc > 0.0) ? yslope : -yslope;
    }

#if (AMREX_SPACEDIM==3)
    // Z direction
    if (flag(i,j,k  ).isSingleValued() or
       (flag(i,j,k-1).isSingleValued() or not flag(i,j,k).isConnected(0,0,-1)) or
       (flag(i,j,k+1).isSingleValued() or not flag(i,j,k).isConnected(0,0, 1)))
    {
        amrex::Real detAtA_z =
            AtA[0][0]*(AtA[1][1]*Atb[2]    - Atb[1]   *AtA[1][2]) -
            AtA[0][1]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]) +
            Atb[0]   *(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

        zs = detAtA_z / detAtA;

    }
    else
    {
        amrex::Real du_zl = 2.0*(state(i,j,k  ,n) - state(i,j,k-1,n));
        amrex::Real du_zr = 2.0*(state(i,j,k+1,n) - state(i,j,k  ,n));
        amrex::Real du_zc = 0.5*(state(i,j,k+1,n) - state(i,j,k-1,n));

        amrex::Real zslope = amrex::min(std::abs(du_zl), std::abs(du_zc), std::abs(du_zr));
        zslope = (du_zr*du_zl > 0.0) ? zslope : 0.0;
        zs = (du_zc > 0.0) ? zslope : -zslope;
    }
#endif
    return {D_DECL(xs,ys,zs)};
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
iamr_slopes_extdir_eb (D_DECL(int i, int j, int k), int n,
                       amrex::Array4<amrex::Real const> const& state,
                       amrex::Array4<amrex::Real const> const& ccent,
                       amrex::Array4<amrex::EBCellFlag const> const& flag,
                       amrex::BCRec const*  d_bcrec,
                       amrex::Box const&  domain   ) noexcept
{
#if (AMREX_SPACEDIM==2)
    const int k = 0;
#endif

    // First get EB-aware slope that doesn't know about extdir
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> eb_slope = iamr_slopes_eb(D_DECL(i,j,k),n,state,ccent,flag);
    D_TERM( amrex::Real xslope = eb_slope[0];,
            amrex::Real yslope = eb_slope[1];,
            amrex::Real zslope = eb_slope[2];);

    // Now correct only those cells which are not affected by EB but are by extdir
    bool extlo = d_bcrec[n].lo(0) == amrex::BCType::ext_dir || d_bcrec[n].lo(0) == amrex::BCType::hoextrap;
    bool exthi = d_bcrec[n].hi(0) == amrex::BCType::ext_dir || d_bcrec[n].hi(0) == amrex::BCType::hoextrap;

    if ( extlo  and i == domain.smallEnd(0) )
    {
        if( !( flag(i  ,j,k).isSingleValued() or
              (flag(i-1,j,k).isSingleValued() or not flag(i,j,k).isConnected(-1,0,0)) or
              (flag(i+1,j,k).isSingleValued() or not flag(i,j,k).isConnected( 1,0,0))) )
        {
           amrex::Real dl = 2.0*(state(i  ,j,k,n) - state(i-1,j,k,n));
           amrex::Real dr = 2.0*(state(i+1,j,k,n) - state(i  ,j,k,n));
           amrex::Real dc =
              (state(i+1,j,k,n)+3.0*state(i,j,k,n)-4.0*state(i-1,j,k,n))/3.0;

           amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           xslope = (dc > 0.0) ? slope : -slope;
        }

    }
    else if ( exthi and i == domain.bigEnd(0) )
    {
        if( !( flag(i  ,j,k).isSingleValued() or
              (flag(i-1,j,k).isSingleValued() or not flag(i,j,k).isConnected(-1,0,0)) or
              (flag(i+1,j,k).isSingleValued() or not flag(i,j,k).isConnected( 1,0,0))) )
        {
           amrex::Real dl = 2.0*(state(i  ,j,k,n) - state(i-1,j,k,n));
           amrex::Real dr = 2.0*(state(i+1,j,k,n) - state(i  ,j,k,n));
           amrex::Real dc =
              (4.0*state(i+1,j,k,n)-3.0*state(i,j,k,n)-state(i-1,j,k,n))/3.0;

           amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           xslope = (dc > 0.0) ? slope : -slope;
        }
    }

    extlo = d_bcrec[n].lo(1) == amrex::BCType::ext_dir || d_bcrec[n].lo(1) == amrex::BCType::hoextrap;
    exthi = d_bcrec[n].hi(1) == amrex::BCType::ext_dir || d_bcrec[n].hi(1) == amrex::BCType::hoextrap;

    if ( extlo and j == domain.smallEnd(1))
    {
        if( !( flag(i,j  ,k).isSingleValued() or
              (flag(i,j-1,k).isSingleValued() or not flag(i,j,k).isConnected(0,-1,0)) or
              (flag(i,j+1,k).isSingleValued() or not flag(i,j,k).isConnected(0, 1,0))) )
        {
           amrex::Real dl = 2.0*(state(i,j  ,k,n) - state(i,j-1,k,n));
           amrex::Real dr = 2.0*(state(i,j+1,k,n) - state(i,j  ,k,n));
           amrex::Real dc =
              (state(i,j+1,k,n)+3.0*state(i,j,k,n)-4.0*state(i,j-1,k,n))/3.0;

           amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           yslope = (dc > 0.0) ? slope : -slope;
        }
    }
    else if ( exthi and j == domain.bigEnd(1))
    {
        if( !(flag(i,j  ,k).isSingleValued() or
           (flag(i,j-1,k).isSingleValued() or not flag(i,j,k).isConnected(0,-1,0)) or
           (flag(i,j+1,k).isSingleValued() or not flag(i,j,k).isConnected(0, 1,0))) )
        {
           amrex::Real dl = 2.0*(state(i,j  ,k,n) - state(i,j-1,k,n));
           amrex::Real dr = 2.0*(state(i,j+1,k,n) - state(i,j  ,k,n));
           amrex::Real dc =
              (4.0*state(i,j+1,k,n)-3.0*state(i,j,k,n)-state(i,j-1,k,n))/3.0;

           amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           yslope = (dc > 0.0) ? slope : -slope;
        }
    }

#if (AMREX_SPACEDIM==3)
    extlo = d_bcrec[n].lo(2) == amrex::BCType::ext_dir || d_bcrec[n].lo(2) == amrex::BCType::hoextrap;
    exthi = d_bcrec[n].hi(2) == amrex::BCType::ext_dir || d_bcrec[n].hi(2) == amrex::BCType::hoextrap;

    if ( extlo and k == domain.smallEnd(2) )
    {
        if( !( flag(i,j,k  ).isSingleValued() or
              (flag(i,j,k-1).isSingleValued() or not flag(i,j,k).isConnected(0,0,-1)) or
              (flag(i,j,k+1).isSingleValued() or not flag(i,j,k).isConnected(0,0, 1))) )
        {
           amrex::Real dl = 2.0*(state(i,j,k  ,n) - state(i,j,k-1,n));
           amrex::Real dr = 2.0*(state(i,j,k+1,n) - state(i,j,k  ,n));
           amrex::Real dc =
              (state(i,j,k+1,n)+3.0*state(i,j,k,n)-4.0*state(i,j,k-1,n))/3.0;

           amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           zslope = (dc > 0.0) ? slope : -slope;
        }
    }
    else if ( exthi and k  == domain.bigEnd(2))
    {
        if( !( flag(i,j,k  ).isSingleValued() or
             (flag(i,j,k-1).isSingleValued() or not flag(i,j,k).isConnected(0,0,-1)) or
             (flag(i,j,k+1).isSingleValued() or not flag(i,j,k).isConnected(0,0, 1))) )
        {
           amrex::Real dl = 2.0*(state(i,j,k  ,n) - state(i,j,k-1,n));
           amrex::Real dr = 2.0*(state(i,j,k+1,n) - state(i,j,k  ,n));
           amrex::Real dc =
              (4.0*state(i,j,k+1,n)-3.0*state(i,j,k,n)-state(i,j,k-1,n))/3.0;

           amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           zslope = (dc > 0.0) ? slope : -slope;
        }
    }
#endif
    return {D_DECL(xslope,yslope,zslope)};
}


}


#endif
