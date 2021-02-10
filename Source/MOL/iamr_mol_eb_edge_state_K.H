#ifndef IAMR_EB_EDGE_STATE_MOL_K_H_
#define IAMR_EB_EDGE_STATE_MOL_K_H_

#include <AMReX_EB_slopes_K.H>
#include <AMReX_BCRec.H>
#include <iamr_constants.H>


namespace {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_eb_xedge_state_mol_extdir ( D_DECL(int i, int j, int k), int n,
                                             amrex::Array4<amrex::Real const> const& q,
                                             amrex::Array4<amrex::Real const> const& umac,
					     AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
							  amrex::Array4<amrex::Real const> const& fcy,
							  amrex::Array4<amrex::Real const> const& fcz),
                                             amrex::Array4<amrex::Real const> const& ccc,
                                             amrex::Array4<amrex::EBCellFlag const> const& flag,
                                             amrex::BCRec const* const d_bcrec,
                                             amrex::Box const&  domain ) noexcept
{

#if (AMREX_SPACEDIM==2)
    const int k = 0;
#endif

    amrex::Real qs;
    int  domlo = domain.smallEnd(0);
    int  domhi = domain.bigEnd(0);
    bool extlo = d_bcrec[n].lo(0) == amrex::BCType::ext_dir; 
    bool exthi = d_bcrec[n].hi(0) == amrex::BCType::ext_dir; 

    if ( extlo and i <= domlo)
    {
        qs = q(domlo-1,j,k,n);
    }
    else if ( exthi and i >= domhi+1)
    {
        qs = q(domhi+1,j,k,n);
    }
    else
    {
        const int domain_ilo = domain.smallEnd(0);
	const int domain_ihi = domain.bigEnd(0);
	const int domain_jlo = domain.smallEnd(1);
	const int domain_jhi = domain.bigEnd(1);
#if (AMREX_SPACEDIM == 3)
	const int domain_klo = domain.smallEnd(2);
	const int domain_khi = domain.bigEnd(2);
#endif

        amrex::Real yf = fcx(i,j,k,0); // local (y,z) of centroid of z-face we are extrapolating to
#if (AMREX_SPACEDIM==3)
        amrex::Real zf = fcx(i,j,k,1);
#endif

        D_TERM(amrex::Real xc = ccc(i,j,k,0);, // centroid of cell (i,j,k)
               amrex::Real yc = ccc(i,j,k,1);,
               amrex::Real zc = ccc(i,j,k,2););

        D_TERM(amrex::Real delta_x = 0.5 + xc;,
               amrex::Real delta_y = yf  - yc;,
               amrex::Real delta_z = zf  - zc;);

        amrex::Real cc_qmax = amrex::max(q(i,j,k,n),q(i-1,j,k,n));
        amrex::Real cc_qmin = amrex::min(q(i,j,k,n),q(i-1,j,k,n));

	AMREX_D_TERM(bool extdir_or_ho_ilo = (d_bcrec[n].lo(0) == amrex::BCType::ext_dir) or
		                             (d_bcrec[n].lo(0) == amrex::BCType::hoextrap);,
		     bool extdir_or_ho_jlo = (d_bcrec[n].lo(1) == amrex::BCType::ext_dir) or
                                             (d_bcrec[n].lo(1) == amrex::BCType::hoextrap);,
		     bool extdir_or_ho_klo = (d_bcrec[n].lo(2) == amrex::BCType::ext_dir) or
                                             (d_bcrec[n].lo(2) == amrex::BCType::hoextrap););

	AMREX_D_TERM(bool extdir_or_ho_ihi = (d_bcrec[n].hi(0) == amrex::BCType::ext_dir) or
                                             (d_bcrec[n].hi(0) == amrex::BCType::hoextrap);,
		     bool extdir_or_ho_jhi = (d_bcrec[n].hi(1) == amrex::BCType::ext_dir) or
                                             (d_bcrec[n].hi(1) == amrex::BCType::hoextrap);,
		     bool extdir_or_ho_khi = (d_bcrec[n].hi(2) == amrex::BCType::ext_dir) or
                                             (d_bcrec[n].hi(2) == amrex::BCType::hoextrap););

	const auto& slopes_eb_hi =
	  amrex_calc_slopes_extdir_eb(i, j, k, n, q, ccc,
				      AMREX_D_DECL(fcx,fcy,fcz), flag,
				      AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
				      AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
				      AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
				      AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

#if (AMREX_SPACEDIM==3)
        amrex::Real qpls = q(i  ,j,k,n) - delta_x * slopes_eb_hi[0]
                                        + delta_y * slopes_eb_hi[1]
                                        + delta_z * slopes_eb_hi[2];
#else
        amrex::Real qpls = q(i  ,j,k,n) - delta_x * slopes_eb_hi[0]
                                        + delta_y * slopes_eb_hi[1];
#endif

        qpls = amrex::max(amrex::min(qpls, cc_qmax), cc_qmin);

        D_TERM(xc = ccc(i-1,j,k,0);, // centroid of cell (i-1,j,k)
               yc = ccc(i-1,j,k,1);,
               zc = ccc(i-1,j,k,2););

        D_TERM( delta_x = 0.5 - xc;,
                delta_y = yf  - yc;,
                delta_z = zf  - zc;);

        // Compute slopes of component "n" of q
        const auto& slopes_eb_lo =
	  amrex_calc_slopes_extdir_eb(i-1, j, k, n, q, ccc,
				      AMREX_D_DECL(fcx,fcy,fcz), flag,
				      AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo), 
				      AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi), 
				      AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo), 
				      AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

#if (AMREX_SPACEDIM==3)
        amrex::Real qmns = q(i-1,j,k,n) + delta_x * slopes_eb_lo[0]
                                        + delta_y * slopes_eb_lo[1]
                                        + delta_z * slopes_eb_lo[2];
#else
        amrex::Real qmns = q(i-1,j,k,n) + delta_x * slopes_eb_lo[0]
                                        + delta_y * slopes_eb_lo[1];
#endif

        qmns = amrex::max(amrex::min(qmns, cc_qmax), cc_qmin);

        if (umac(i,j,k) > small_vel)
        {
            qs = qmns;
        }
        else if (umac(i,j,k) < - small_vel)
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
amrex::Real iamr_eb_xedge_state_mol ( D_DECL(int i, int j, int k), int n,
                                      amrex::Array4<amrex::Real const> const& q,
                                      amrex::Array4<amrex::Real const> const& umac,
                                      amrex::Array4<amrex::Real const> const& fcx,
                                      amrex::Array4<amrex::Real const> const& ccc,
                                      amrex::Array4<amrex::EBCellFlag const> const& flag ) noexcept
{
#if (AMREX_SPACEDIM==2)
    const int k = 0;
#endif

    // local (y,z) of centroid of x-face we are extrapolating to
    amrex::Real yf = fcx(i,j,k,0);
#if (AMREX_SPACEDIM==3)
    amrex::Real zf = fcx(i,j,k,1);
#endif

    D_TERM(amrex::Real xc = ccc(i,j,k,0);, // centroid of cell (i,j,k)
           amrex::Real yc = ccc(i,j,k,1);,
           amrex::Real zc = ccc(i,j,k,2););

    D_TERM(amrex::Real delta_x = 0.5 + xc;,
           amrex::Real delta_y = yf  - yc;,
           amrex::Real delta_z = zf  - zc;);

    amrex::Real cc_qmax = amrex::max(q(i,j,k,n),q(i-1,j,k,n));
    amrex::Real cc_qmin = amrex::min(q(i,j,k,n),q(i-1,j,k,n));

    // Compute slopes of component "n" of q
    const auto& slopes_eb_hi = amrex_calc_slopes_eb(i, j, k, n, q, ccc, flag);

#if (AMREX_SPACEDIM==3)
    amrex::Real qpls = q(i  ,j,k,n) - delta_x * slopes_eb_hi[0]
                                    + delta_y * slopes_eb_hi[1]
                                    + delta_z * slopes_eb_hi[2];
#else
    amrex::Real qpls = q(i  ,j,k,n) - delta_x * slopes_eb_hi[0]
                                    + delta_y * slopes_eb_hi[1];
#endif

    qpls = amrex::max(amrex::min(qpls, cc_qmax), cc_qmin);

    D_TERM(xc = ccc(i-1,j,k,0);, // centroid of cell (i-1,j,k)
           yc = ccc(i-1,j,k,1);,
           zc = ccc(i-1,j,k,2););

    D_TERM(delta_x = 0.5 - xc;,
           delta_y = yf  - yc;,
           delta_z = zf  - zc;);

    // Compute slopes of component "n" of q
    const auto& slopes_eb_lo = amrex_calc_slopes_eb(i-1, j, k, n, q, ccc, flag);

#if (AMREX_SPACEDIM==3)
    amrex::Real qmns = q(i-1,j,k,n) + delta_x * slopes_eb_lo[0]
                                    + delta_y * slopes_eb_lo[1]
                                    + delta_z * slopes_eb_lo[2];
#else
    amrex::Real qmns = q(i-1,j,k,n) + delta_x * slopes_eb_lo[0]
                                    + delta_y * slopes_eb_lo[1];
#endif

    qmns = amrex::max(amrex::min(qmns, cc_qmax), cc_qmin);


    amrex::Real qs;

    if (umac(i,j,k) > small_vel)
    {
        qs = qmns;
    }
    else if (umac(i,j,k) < - small_vel)
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
amrex::Real iamr_eb_yedge_state_mol_extdir ( D_DECL(int i, int j, int k), int n,
                                             amrex::Array4<amrex::Real const> const& q,
                                             amrex::Array4<amrex::Real const> const& vmac,
					     AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                             amrex::Array4<amrex::Real const> const& fcy,
							  amrex::Array4<amrex::Real const> const& fcz),
                                             amrex::Array4<amrex::Real const> const& ccc,
                                             amrex::Array4<amrex::EBCellFlag const> const& flag,
                                             amrex::BCRec const* const d_bcrec,
                                             amrex::Box const&  domain ) noexcept
{
#if (AMREX_SPACEDIM==2)
    const int k = 0;
#endif

    amrex::Real qs;
    int domlo = domain.smallEnd(1);
    int domhi = domain.bigEnd(1);
    bool extlo = d_bcrec[n].lo(1) == amrex::BCType::ext_dir;
    bool exthi = d_bcrec[n].hi(1) == amrex::BCType::ext_dir;

    if ( extlo and j <= domlo)
    {
        qs = q(i,domlo-1,k,n);
    }
    else if (exthi and j >= domhi+1)
    {
        qs = q(i,domhi+1,k,n);
    }
    else
    {
        const int domain_ilo = domain.smallEnd(0);
	const int domain_ihi = domain.bigEnd(0);
	const int domain_jlo = domain.smallEnd(1);
	const int domain_jhi = domain.bigEnd(1);
#if (AMREX_SPACEDIM == 3)
	const int domain_klo = domain.smallEnd(2);
	const int domain_khi = domain.bigEnd(2);
#endif
	
	AMREX_D_TERM(bool extdir_or_ho_ilo = (d_bcrec[n].lo(0) == amrex::BCType::ext_dir) or
		                             (d_bcrec[n].lo(0) == amrex::BCType::hoextrap);,
		     bool extdir_or_ho_jlo = (d_bcrec[n].lo(1) == amrex::BCType::ext_dir) or
                                             (d_bcrec[n].lo(1) == amrex::BCType::hoextrap);,
		     bool extdir_or_ho_klo = (d_bcrec[n].lo(2) == amrex::BCType::ext_dir) or
                                             (d_bcrec[n].lo(2) == amrex::BCType::hoextrap););

	AMREX_D_TERM(bool extdir_or_ho_ihi = (d_bcrec[n].hi(0) == amrex::BCType::ext_dir) or
                                             (d_bcrec[n].hi(0) == amrex::BCType::hoextrap);,
		     bool extdir_or_ho_jhi = (d_bcrec[n].hi(1) == amrex::BCType::ext_dir) or
                                             (d_bcrec[n].hi(1) == amrex::BCType::hoextrap);,
		     bool extdir_or_ho_khi = (d_bcrec[n].hi(2) == amrex::BCType::ext_dir) or
                                             (d_bcrec[n].hi(2) == amrex::BCType::hoextrap););

        amrex::Real xf = fcy(i,j,k,0); // local (x,z) of centroid of z-face we are extrapolating to
#if (AMREX_SPACEDIM==3)
        amrex::Real zf = fcy(i,j,k,1);
#endif

        D_TERM(amrex::Real xc = ccc(i,j,k,0);, // centroid of cell (i,j,k)
               amrex::Real yc = ccc(i,j,k,1);,
               amrex::Real zc = ccc(i,j,k,2););

        D_TERM(amrex::Real delta_x = xf  - xc;,
               amrex::Real delta_y = 0.5 + yc;,
               amrex::Real delta_z = zf  - zc;);

        amrex::Real cc_qmax = amrex::max(q(i,j,k,n),q(i,j-1,k,n));
        amrex::Real cc_qmin = amrex::min(q(i,j,k,n),q(i,j-1,k,n));

        // Compute slopes of component "n" of q
        const auto& slopes_eb_hi =
	  amrex_calc_slopes_extdir_eb(i, j, k, n, q, ccc,
				      AMREX_D_DECL(fcx,fcy,fcz), flag,
				      AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo), 
				      AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi), 
				      AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo), 
				      AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

#if (AMREX_SPACEDIM==3)
        amrex::Real qpls = q(i  ,j,k,n) + delta_x * slopes_eb_hi[0]
                                        - delta_y * slopes_eb_hi[1]
                                        + delta_z * slopes_eb_hi[2];
#else
        amrex::Real qpls = q(i  ,j,k,n) + delta_x * slopes_eb_hi[0]
                                        - delta_y * slopes_eb_hi[1];
#endif

        qpls = amrex::max(amrex::min(qpls, cc_qmax), cc_qmin);

        D_TERM(xc = ccc(i,j-1,k,0);, // centroid of cell (i,j-1,k)
               yc = ccc(i,j-1,k,1);,
               zc = ccc(i,j-1,k,2););

        D_TERM(delta_x = xf  - xc;,
               delta_y = 0.5 - yc;,
               delta_z = zf  - zc;);

        // Compute slopes of component "n" of q
        const auto& slopes_eb_lo =
	  amrex_calc_slopes_extdir_eb(i, j-1, k, n, q, ccc,
				      AMREX_D_DECL(fcx,fcy,fcz), flag,
				      AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo), 
				      AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi), 
				      AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo), 
				      AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

#if (AMREX_SPACEDIM==3)
        amrex::Real qmns = q(i,j-1,k,n) + delta_x * slopes_eb_lo[0]
                                        + delta_y * slopes_eb_lo[1]
                                        + delta_z * slopes_eb_lo[2];
#else
        amrex::Real qmns = q(i,j-1,k,n) + delta_x * slopes_eb_lo[0]
                                        + delta_y * slopes_eb_lo[1];
#endif

        qmns = amrex::max(amrex::min(qmns, cc_qmax), cc_qmin);

        if (vmac(i,j,k) > small_vel)
        {
            qs = qmns;
        }
        else if (vmac(i,j,k) < - small_vel)
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
amrex::Real iamr_eb_yedge_state_mol ( D_DECL(int i, int j, int k), int n,
                                      amrex::Array4<amrex::Real const> const& q,
                                      amrex::Array4<amrex::Real const> const& vmac,
                                      amrex::Array4<amrex::Real const> const& fcy,
                                      amrex::Array4<amrex::Real const> const& ccc,
                                      amrex::Array4<amrex::EBCellFlag const> const& flag ) noexcept
{
#if (AMREX_SPACEDIM==2)
    const int k = 0;
#endif

    // local (x,z) of centroid of z-face we are extrapolating to
    amrex::Real xf = fcy(i,j,k,0);
#if (AMREX_SPACEDIM==3)
    amrex::Real zf = fcy(i,j,k,1);
#endif

    D_TERM(amrex::Real xc = ccc(i,j,k,0);, // centroid of cell (i,j,k)
           amrex::Real yc = ccc(i,j,k,1);,
           amrex::Real zc = ccc(i,j,k,2););

    D_TERM(amrex::Real delta_x = xf  - xc;,
           amrex::Real delta_y = 0.5 + yc;,
           amrex::Real delta_z = zf  - zc;);

    amrex::Real cc_qmax = amrex::max(q(i,j,k,n),q(i,j-1,k,n));
    amrex::Real cc_qmin = amrex::min(q(i,j,k,n),q(i,j-1,k,n));

    // Compute slopes of component "n" of q
    const auto& slopes_eb_hi = amrex_calc_slopes_eb(i, j, k, n, q, ccc, flag);

#if (AMREX_SPACEDIM==3)
    amrex::Real qpls = q(i  ,j,k,n) + delta_x * slopes_eb_hi[0]
                                    - delta_y * slopes_eb_hi[1]
                                    + delta_z * slopes_eb_hi[2];
#else
    amrex::Real qpls = q(i  ,j,k,n) + delta_x * slopes_eb_hi[0]
                                    - delta_y * slopes_eb_hi[1];
#endif

    qpls = amrex::max(amrex::min(qpls, cc_qmax), cc_qmin);

    D_TERM(xc = ccc(i,j-1,k,0);, // centroid of cell (i,j-1,k)
           yc = ccc(i,j-1,k,1);,
           zc = ccc(i,j-1,k,2););

    D_TERM(delta_x = xf  - xc;,
           delta_y = 0.5 - yc;,
           delta_z = zf  - zc;);

    // Compute slopes of component "n" of q
    const auto& slopes_eb_lo = amrex_calc_slopes_eb(i, j-1, k, n, q, ccc, flag);

#if (AMREX_SPACEDIM==3)
    amrex::Real qmns = q(i,j-1,k,n) + delta_x * slopes_eb_lo[0]
                                    + delta_y * slopes_eb_lo[1]
                                    + delta_z * slopes_eb_lo[2];
#else
    amrex::Real qmns = q(i,j-1,k,n) + delta_x * slopes_eb_lo[0]
                                    + delta_y * slopes_eb_lo[1];
#endif

    qmns = amrex::max(amrex::min(qmns, cc_qmax), cc_qmin);

    amrex::Real qs;

    if (vmac(i,j,k) > small_vel)
    {
        qs = qmns;
    }
    else if (vmac(i,j,k) < - small_vel)
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
amrex::Real iamr_eb_zedge_state_mol_extdir ( int i, int j, int k, int n,
                                             amrex::Array4<amrex::Real const> const& q,
                                             amrex::Array4<amrex::Real const> const& wmac,
					     AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
							  amrex::Array4<amrex::Real const> const& fcy,
							  amrex::Array4<amrex::Real const> const& fcz),
                                             amrex::Array4<amrex::Real const> const& ccc,
                                             amrex::Array4<amrex::EBCellFlag const> const& flag,
                                             amrex::BCRec const* const d_bcrec,
                                             amrex::Box const&  domain ) noexcept
{
    amrex::Real qs;
    int domlo = domain.smallEnd(2);
    int domhi = domain.bigEnd(2);
    bool extlo = d_bcrec[n].lo(2) == amrex::BCType::ext_dir;
    bool exthi = d_bcrec[n].hi(2) == amrex::BCType::ext_dir;

    if ( extlo and k <= domlo)
    {
        qs = q(i,j,domlo-1,n);
    }
    else if ( exthi and k >= domhi+1)
    {
        qs = q(i,j,domhi+1,n);
    }
    else
    {
        const int domain_ilo = domain.smallEnd(0);
	const int domain_ihi = domain.bigEnd(0);
	const int domain_jlo = domain.smallEnd(1);
	const int domain_jhi = domain.bigEnd(1);
#if (AMREX_SPACEDIM == 3)
	const int domain_klo = domain.smallEnd(2);
	const int domain_khi = domain.bigEnd(2);
#endif

	AMREX_D_TERM(bool extdir_or_ho_ilo = (d_bcrec[n].lo(0) == amrex::BCType::ext_dir) or
		                             (d_bcrec[n].lo(0) == amrex::BCType::hoextrap);,
		     bool extdir_or_ho_jlo = (d_bcrec[n].lo(1) == amrex::BCType::ext_dir) or
                                             (d_bcrec[n].lo(1) == amrex::BCType::hoextrap);,
		     bool extdir_or_ho_klo = (d_bcrec[n].lo(2) == amrex::BCType::ext_dir) or
                                             (d_bcrec[n].lo(2) == amrex::BCType::hoextrap););

	AMREX_D_TERM(bool extdir_or_ho_ihi = (d_bcrec[n].hi(0) == amrex::BCType::ext_dir) or
                                             (d_bcrec[n].hi(0) == amrex::BCType::hoextrap);,
		     bool extdir_or_ho_jhi = (d_bcrec[n].hi(1) == amrex::BCType::ext_dir) or
                                             (d_bcrec[n].hi(1) == amrex::BCType::hoextrap);,
		     bool extdir_or_ho_khi = (d_bcrec[n].hi(2) == amrex::BCType::ext_dir) or
                                             (d_bcrec[n].hi(2) == amrex::BCType::hoextrap););

        amrex::Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
        amrex::Real yf = fcz(i,j,k,1);

        amrex::Real xc = ccc(i,j,k,0); // centroid of cell (i,j,k)
        amrex::Real yc = ccc(i,j,k,1);
        amrex::Real zc = ccc(i,j,k,2);

        amrex::Real delta_x = xf  - xc;
        amrex::Real delta_y = yf  - yc;
        amrex::Real delta_z = 0.5 + zc;

        amrex::Real cc_qmax = amrex::max(q(i,j,k,n),q(i,j,k-1,n));
        amrex::Real cc_qmin = amrex::min(q(i,j,k,n),q(i,j,k-1,n));

        // Compute slopes of component "n" of q
        const auto& slopes_eb_hi =
	  amrex_calc_slopes_extdir_eb(i,j,k,n,q,ccc,
				      AMREX_D_DECL(fcx,fcy,fcz), flag,
				      AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo), 
				      AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi), 
				      AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo), 
				      AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

        amrex::Real qpls = q(i,j,k  ,n) + delta_x * slopes_eb_hi[0]
                                        + delta_y * slopes_eb_hi[1]
                                        - delta_z * slopes_eb_hi[2];

        qpls = amrex::max(amrex::min(qpls, cc_qmax), cc_qmin);

        xc = ccc(i,j,k-1,0); // centroid of cell (i,j,k-1)
        yc = ccc(i,j,k-1,1);
        zc = ccc(i,j,k-1,2);

        delta_x = xf  - xc;
        delta_y = yf  - yc;
        delta_z = 0.5 - zc;

        // Compute slopes of component "n" of q
        const auto& slopes_eb_lo =
	  amrex_calc_slopes_extdir_eb(i,j,k-1,n,q,ccc,
				      AMREX_D_DECL(fcx,fcy,fcz),flag,
				      AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo), 
				      AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi), 
				      AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo), 
				      AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

        amrex::Real qmns = q(i,j,k-1,n) + delta_x * slopes_eb_lo[0]
                                        + delta_y * slopes_eb_lo[1]
                                        + delta_z * slopes_eb_lo[2];

        qmns = amrex::max(amrex::min(qmns, cc_qmax), cc_qmin);

        if (wmac(i,j,k) > small_vel)
        {
            qs = qmns;
        }
        else if (wmac(i,j,k) < -small_vel)
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
amrex::Real iamr_eb_zedge_state_mol ( int i, int j, int k, int n,
                                      amrex::Array4<amrex::Real const> const& q,
                                      amrex::Array4<amrex::Real const> const& wmac,
                                      amrex::Array4<amrex::Real const> const& fcz,
                                      amrex::Array4<amrex::Real const> const& ccc,
                                      amrex::Array4<amrex::EBCellFlag const> const& flag ) noexcept
{

    amrex::Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
    amrex::Real yf = fcz(i,j,k,1);

    amrex::Real xc = ccc(i,j,k,0); // centroid of cell (i,j,k)
    amrex::Real yc = ccc(i,j,k,1);
    amrex::Real zc = ccc(i,j,k,2);

    amrex::Real delta_x = xf  - xc;
    amrex::Real delta_y = yf  - yc;
    amrex::Real delta_z = 0.5 + zc;

    amrex::Real cc_qmax = amrex::max(q(i,j,k,n),q(i,j,k-1,n));
    amrex::Real cc_qmin = amrex::min(q(i,j,k,n),q(i,j,k-1,n));

    // Compute slopes of component "n" of q
    const auto& slopes_eb_hi = amrex_calc_slopes_eb(i,j,k,n,q,ccc,flag);

    amrex::Real qpls = q(i,j,k  ,n) + delta_x * slopes_eb_hi[0]
                                    + delta_y * slopes_eb_hi[1]
                                    - delta_z * slopes_eb_hi[2];

    qpls = amrex::max(amrex::min(qpls, cc_qmax), cc_qmin);

    xc = ccc(i,j,k-1,0); // centroid of cell (i,j,k-1)
    yc = ccc(i,j,k-1,1);
    zc = ccc(i,j,k-1,2);

    delta_x = xf  - xc;
    delta_y = yf  - yc;
    delta_z = 0.5 - zc;

    // Compute slopes of component "n" of q
    const auto& slopes_eb_lo = amrex_calc_slopes_eb(i,j,k-1,n,q,ccc,flag);

    amrex::Real qmns = q(i,j,k-1,n) + delta_x * slopes_eb_lo[0]
                                    + delta_y * slopes_eb_lo[1]
                                    + delta_z * slopes_eb_lo[2];

    qmns = amrex::max(amrex::min(qmns, cc_qmax), cc_qmin);

    amrex::Real qs;

    if (wmac(i,j,k) > small_vel)
    {
        qs = qmns;
    }
    else if (wmac(i,j,k) < -small_vel)
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
