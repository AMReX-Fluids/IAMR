#ifndef IAMR_EDGE_STATE_MOL_K_H_
#define IAMR_EDGE_STATE_MOL_K_H_

#include <iamr_slopes_K.H>
#include <iamr_constants.H>

namespace {


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_xedge_state_mol_extdir ( int i, int j, int k, int n,
                                          amrex::Array4<amrex::Real const> const& q,
                                          amrex::Array4<amrex::Real const> const& umac,
                                          bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real  qs;

    if (edlo and i <= domlo)
    {
        qs = q(domlo-1,j,k,n);
    }
    else if ( edhi and i >= domhi+1)
    {
        qs = q(domhi+1,j,k,n);
    }
    else
    {
        amrex::Real qpls = q(i  ,j,k,n) - 0.5 * iamr_xslope_extdir( i  , j, k, n, q, edlo, edhi, domlo, domhi );
        amrex::Real qmns = q(i-1,j,k,n) + 0.5 * iamr_xslope_extdir( i-1, j, k, n, q, edlo, edhi, domlo, domhi );

        if ( umac(i,j,k) > small_vel)
        {
            qs = qmns;
        }
        else if ( umac(i,j,k) < -small_vel)
        {
            qs = qpls;
        }
        else
        {
            qs = 0.5*(qmns+qpls);
        }
    }

    return qs;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_xedge_state_mol ( int i, int j, int k, int n,
                                   amrex::Array4<amrex::Real const> const& q,
                                   amrex::Array4<amrex::Real const> const& umac ) noexcept
{
    amrex::Real qs;
    amrex::Real qpls = q(i  ,j,k,n) - 0.5 * iamr_xslope( i  , j, k, n, q );
    amrex::Real qmns = q(i-1,j,k,n) + 0.5 * iamr_xslope( i-1, j, k, n, q );

    if ( umac(i,j,k) > small_vel)
    {
        qs = qmns;
    }
    else if ( umac(i,j,k) < -small_vel)
    {
        qs = qpls;
    }
    else
    {
        qs = 0.5*(qmns+qpls);
    }

    return qs;
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_yedge_state_mol_extdir ( int i, int j, int k, int n,
                                          amrex::Array4<amrex::Real const> const& q,
                                          amrex::Array4<amrex::Real const> const& vmac,
                                          bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real qs;

    if (edlo and j <= domlo)
    {
        qs = q(i,domlo-1,k,n);
    }
    else if ( edhi and j >= domhi+1)
    {
        qs = q(i,domhi+1,k,n);
    }
    else
    {
        amrex::Real qpls = q(i,j  ,k,n) - 0.5 * iamr_yslope_extdir( i, j  , k, n, q, edlo, edhi, domlo, domhi );
        amrex::Real qmns = q(i,j-1,k,n) + 0.5 * iamr_yslope_extdir( i, j-1, k, n, q, edlo, edhi, domlo, domhi );

        if ( vmac(i,j,k) > small_vel)
        {
            qs = qmns;
        }
        else if ( vmac(i,j,k) < -small_vel)
        {
            qs = qpls;
        }
        else
        {
            qs = 0.5*(qmns+qpls);
        }
    }

    return qs;
}



AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_yedge_state_mol ( int i, int j, int k, int n,
                                   amrex::Array4<amrex::Real const> const& q,
                                   amrex::Array4<amrex::Real const> const& vmac ) noexcept
{
    amrex::Real qs;
    amrex::Real qpls = q(i,j  ,k,n) - 0.5 * iamr_yslope( i, j  , k, n, q );
    amrex::Real qmns = q(i,j-1,k,n) + 0.5 * iamr_yslope( i, j-1, k, n, q );

    if ( vmac(i,j,k) > small_vel)
    {
        qs = qmns;
    }
    else if ( vmac(i,j,k) < -small_vel)
    {
        qs = qpls;
    }
    else
    {
        qs = 0.5*(qmns+qpls);
    }

    return qs;

}

#if (AMREX_SPACEDIM==3)

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_zedge_state_mol_extdir ( int i, int j, int k, int n,
                                          amrex::Array4<amrex::Real const> const& q,
                                          amrex::Array4<amrex::Real const> const& wmac,
                                          bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real qs;

    if (edlo and k <= domlo)
    {
        qs = q(i,j,domlo-1,n);
    }
    else if ( edhi and k >= domhi+1)
    {
        qs = q(i,j,domhi+1,n);
    }
    else
    {
        amrex::Real qpls = q(i,j,k  ,n) - 0.5 * iamr_zslope_extdir( i, j, k  , n, q, edlo, edhi, domlo, domhi );
        amrex::Real qmns = q(i,j,k-1,n) + 0.5 * iamr_zslope_extdir( i, j, k-1, n, q, edlo, edhi, domlo, domhi );

        if ( wmac(i,j,k) > small_vel)
        {
            qs = qmns;
        }
        else if ( wmac(i,j,k) < -small_vel)
        {
            qs = qpls;
        }
        else
        {
            qs = 0.5*(qmns+qpls);
        }
    }

    return qs;
}



AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_zedge_state_mol ( int i, int j, int k, int n,
                                   amrex::Array4<amrex::Real const> const& q,
                                   amrex::Array4<amrex::Real const> const& wmac ) noexcept
{
    amrex::Real qs;
    amrex::Real qpls = q(i,j,k  ,n) - 0.5 * iamr_zslope( i, j, k  , n, q );
    amrex::Real qmns = q(i,j,k-1,n) + 0.5 * iamr_zslope( i, j, k-1, n, q );

    if ( wmac(i,j,k) > small_vel)
    {
        qs = qmns;
    }
    else if ( wmac(i,j,k) < -small_vel)
    {
        qs = qpls;
    }
    else
    {
        qs = 0.5*(qmns+qpls);
    }

    return qs;

}

#endif


}

#endif
