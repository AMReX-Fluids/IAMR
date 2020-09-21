#ifndef IAMR_EDGE_STATE_MOL_K_H_
#define IAMR_EDGE_STATE_MOL_K_H_

#include <AMReX_Slopes_K.H>
#include <iamr_constants.H>

namespace {


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_xedge_state_mol_extdir ( int i, int j, int k, int n,
                                          amrex::Array4<amrex::Real const> const& q,
                                          amrex::Array4<amrex::Real const> const& umac,
					  amrex::BCRec const* const d_bcrec,
                                          int domlo, int domhi) noexcept
{
    amrex::Real  qs;

    //slope order
    int order = 2;
    
    bool edlo = (d_bcrec[n].lo(0) == amrex::BCType::ext_dir);
    bool edhi = (d_bcrec[n].hi(0) == amrex::BCType::ext_dir);
 
    bool extdir_or_ho_lo = (d_bcrec[n].lo(0) == amrex::BCType::ext_dir) || d_bcrec[n].lo(0) == amrex::BCType::hoextrap;
    bool extdir_or_ho_hi = (d_bcrec[n].hi(0) == amrex::BCType::ext_dir) || d_bcrec[n].lo(0) == amrex::BCType::hoextrap;

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
        amrex::Real qpls = q(i  ,j,k,n) - 0.5 *
	  amrex_calc_xslope_extdir( i  , j, k, n, order, q, extdir_or_ho_lo, extdir_or_ho_hi, domlo, domhi );
        amrex::Real qmns = q(i-1,j,k,n) + 0.5 *
	  amrex_calc_xslope_extdir( i-1, j, k, n, order, q, extdir_or_ho_lo, extdir_or_ho_hi, domlo, domhi );

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
    //slope order
    int order = 2;

    amrex::Real qs;
    amrex::Real qpls = q(i  ,j,k,n) - 0.5 * amrex_calc_xslope( i  , j, k, n, order, q );
    amrex::Real qmns = q(i-1,j,k,n) + 0.5 * amrex_calc_xslope( i-1, j, k, n, order, q );

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
					  amrex::BCRec const* const d_bcrec,
					  int domlo, int domhi) noexcept
{
    //slope order
    int order = 2;

    bool edlo = (d_bcrec[n].lo(1) == amrex::BCType::ext_dir);
    bool edhi = (d_bcrec[n].hi(1) == amrex::BCType::ext_dir);
 
    bool extdir_or_ho_lo = (d_bcrec[n].lo(1) == amrex::BCType::ext_dir) || d_bcrec[n].lo(1) == amrex::BCType::hoextrap;
    bool extdir_or_ho_hi = (d_bcrec[n].hi(1) == amrex::BCType::ext_dir) || d_bcrec[n].lo(1) == amrex::BCType::hoextrap;

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
        amrex::Real qpls = q(i,j  ,k,n) - 0.5 *
	  amrex_calc_yslope_extdir( i, j  , k, n, order, q, extdir_or_ho_lo, extdir_or_ho_hi, domlo, domhi );
        amrex::Real qmns = q(i,j-1,k,n) + 0.5 *
	  amrex_calc_yslope_extdir( i, j-1, k, n, order, q, extdir_or_ho_lo, extdir_or_ho_hi, domlo, domhi );

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
    //slope order
    int order = 2;

    amrex::Real qs;
    amrex::Real qpls = q(i,j  ,k,n) - 0.5 * amrex_calc_yslope( i, j  , k, n, order, q );
    amrex::Real qmns = q(i,j-1,k,n) + 0.5 * amrex_calc_yslope( i, j-1, k, n, order, q );

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
					  amrex::BCRec const* const d_bcrec,
					  int domlo, int domhi) noexcept
{
    //slope order
    int order = 2;

    amrex::Real qs;

    bool edlo = (d_bcrec[n].lo(2) == amrex::BCType::ext_dir);
    bool edhi = (d_bcrec[n].hi(2) == amrex::BCType::ext_dir);
 
    bool extdir_or_ho_lo = (d_bcrec[n].lo(2) == amrex::BCType::ext_dir) || d_bcrec[n].lo(2) == amrex::BCType::hoextrap;
    bool extdir_or_ho_hi = (d_bcrec[n].hi(2) == amrex::BCType::ext_dir) || d_bcrec[n].lo(2) == amrex::BCType::hoextrap;

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
        amrex::Real qpls = q(i,j,k  ,n) - 0.5 *
	  amrex_calc_zslope_extdir( i, j, k  , n, order, q, extdir_or_ho_lo, extdir_or_ho_hi, domlo, domhi );
        amrex::Real qmns = q(i,j,k-1,n) + 0.5 *
	  amrex_calc_zslope_extdir( i, j, k-1, n, order, q, extdir_or_ho_lo, extdir_or_ho_hi, domlo, domhi );

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
    //slope order
    int order = 2;

    amrex::Real qs;
    amrex::Real qpls = q(i,j,k  ,n) - 0.5 * amrex_calc_zslope( i, j, k  , n, order, q );
    amrex::Real qmns = q(i,j,k-1,n) + 0.5 * amrex_calc_zslope( i, j, k-1, n, order, q );

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
