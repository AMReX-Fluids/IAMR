#include <Godunov.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_BCRec.H>
#include <iamr_convection_K.H>
#include <iamr_convection_K.H>
#include <NS_util.H>

using namespace amrex;

namespace {
    std::pair<bool,bool> has_extdir (BCRec const* bcrec, int ncomp, int dir)
    {
        std::pair<bool,bool> r{false,false};
        for (int n = 0; n < ncomp; ++n) {
            r.first = r.first or bcrec[n].lo(dir) == BCType::ext_dir;
            r.second = r.second or bcrec[n].hi(dir) == BCType::ext_dir;
        }
        return r;
    }
}


//
// Auxiliary namespace fluxes
//
namespace fluxes
{
//
// Compute upwind non-normal velocity
//
AMREX_GPU_HOST_DEVICE
Real
upwind(const Real velocity_minus,
       const Real velocity_plus,
       const Real u_edge)
{
    // Small value to protect against tiny velocities used in upwinding
    const Real small_velocity(1.e-10);

    if(std::abs(u_edge) < small_velocity)
        return .5*(velocity_minus+velocity_plus);

    return u_edge > 0 ? velocity_minus : velocity_plus;
}


//
// Compute fluxes on given REGULAR box
//
void
ComputeFluxesOnBox (const Box& a_bx,
                    D_DECL( Array4<Real> const& a_fx,
                            Array4<Real> const& a_fy,
                            Array4<Real> const& a_fz),
                    D_DECL( Array4<Real> const& a_edgeq_x,
                            Array4<Real> const& a_edgeq_y,
                            Array4<Real> const& a_edgeq_z),
                    Array4<Real const> const& a_q,
                    const int a_comp,
                    const int a_ncomp,
                    D_DECL( Array4<Real const> const& a_umac,
                            Array4<Real const> const& a_vmac,
                            Array4<Real const> const& a_wmac),
                    const Box&       a_domain,
                    const Vector<BCRec>& a_bcs,
                    int known_edgestate)
{

    constexpr Real small = 1.e-10;

    const int domain_ilo = a_domain.smallEnd(0);
    const int domain_ihi = a_domain.bigEnd(0);
    const int domain_jlo = a_domain.smallEnd(1);
    const int domain_jhi = a_domain.bigEnd(1);
#if (AMREX_SPACEDIM==3)
    const int domain_klo = a_domain.smallEnd(2);
    const int domain_khi = a_domain.bigEnd(2);
#endif

    D_TERM( const Box& ubx = amrex::surroundingNodes(a_bx,0);,
            const Box& vbx = amrex::surroundingNodes(a_bx,1);,
            const Box& wbx = amrex::surroundingNodes(a_bx,2););

    const auto bc = a_bcs.dataPtr();
    auto d_bcrec  = convertToDeviceVector(a_bcs);

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    auto extdir_lohi = has_extdir(a_bcs.dataPtr(), a_ncomp, 0);
    bool has_extdir_lo = extdir_lohi.first;
    bool has_extdir_hi = extdir_lohi.second;

    if ((has_extdir_lo and domain_ilo >= ubx.smallEnd(0)-1) or
        (has_extdir_hi and domain_ihi <= ubx.bigEnd(0)))
    {
        amrex::ParallelFor(ubx, a_ncomp, [d_bcrec,a_q,domain_ilo,domain_ihi,a_umac,small,a_fx, known_edgestate, a_edgeq_x]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            bool extdir_ilo = d_bcrec[n].lo(0) == BCType::ext_dir;
            bool extdir_ihi = d_bcrec[n].hi(0) == BCType::ext_dir;
            Real qs;
            if ( known_edgestate == 0) {
                if (extdir_ilo and i <= domain_ilo) {
                    qs = a_q(domain_ilo-1,j,k,n);
                } else if (extdir_ihi and i >= domain_ihi+1) {
                    qs = a_q(domain_ihi+1,j,k,n);
                } else {
                    Real qpls = a_q(i,j,k,n) - 0.5 * iamr_xslope_extdir
                        (i,j,k,n,a_q, extdir_ilo, extdir_ihi, domain_ilo, domain_ihi);
                    Real qmns = a_q(i-1,j,k,n) + 0.5 * iamr_xslope_extdir
                        (i-1,j,k,n,a_q, extdir_ilo, extdir_ihi, domain_ilo, domain_ihi);
                    if (a_umac(i,j,k) > small) {
                        qs = qmns;
                    } else if (a_umac(i,j,k) < -small) {
                        qs = qpls;
                    } else {
                        qs = 0.5*(qmns+qpls);
                    }
                }
                a_edgeq_x(i,j,k,n) = qs;
            }
            a_fx(i,j,k,n) =  a_edgeq_x(i,j,k,n) * a_umac(i,j,k);
        });
    }
    else
    {
        amrex::ParallelFor(ubx, a_ncomp, [a_q,a_umac,small,a_fx,known_edgestate,a_edgeq_x]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qpls = a_q(i  ,j,k,n) - 0.5 * iamr_xslope(i  ,j,k,n,a_q);
            Real qmns = a_q(i-1,j,k,n) + 0.5 * iamr_xslope(i-1,j,k,n,a_q);
            Real qs;
            if ( known_edgestate == 0)  {
                if (a_umac(i,j,k) > small) {
                    qs = qmns;
                } else if (a_umac(i,j,k) < -small) {
                    qs = qpls;
                } else {
                    qs = 0.5*(qmns+qpls);
                }

                a_edgeq_x(i,j,k,n) = qs;
            }
            a_fx(i,j,k,n) = a_edgeq_x(i,j,k,n) * a_umac(i,j,k);
        });
    }


    extdir_lohi = has_extdir(a_bcs.dataPtr(), a_ncomp, 1);
    has_extdir_lo = extdir_lohi.first;
    has_extdir_hi = extdir_lohi.second;
    if ((has_extdir_lo and domain_jlo >= vbx.smallEnd(1)-1) or
        (has_extdir_hi and domain_jhi <= vbx.bigEnd(1)))
    {
        amrex::ParallelFor(vbx, a_ncomp, [d_bcrec,a_q,domain_jlo,domain_jhi,a_vmac,small,a_fy,known_edgestate,a_edgeq_y]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qs;

            if ( known_edgestate == 0) {
                bool extdir_jlo = d_bcrec[n].lo(1) == BCType::ext_dir;
                bool extdir_jhi = d_bcrec[n].hi(1) == BCType::ext_dir;
                if (extdir_jlo and j <= domain_jlo) {
                    qs = a_q(i,domain_jlo-1,k,n);
                } else if (extdir_jhi and j >= domain_jhi+1) {
                    qs = a_q(i,domain_jhi+1,k,n);
                } else {
                    Real qpls = a_q(i,j,k,n) - 0.5 * iamr_yslope_extdir
                        (i,j,k,n,a_q, extdir_jlo, extdir_jhi, domain_jlo, domain_jhi);
                    Real qmns = a_q(i,j-1,k,n) + 0.5 * iamr_yslope_extdir
                        (i,j-1,k,n,a_q, extdir_jlo, extdir_jhi, domain_jlo, domain_jhi);
                    if (a_vmac(i,j,k) > small) {
                        qs = qmns;
                    } else if (a_vmac(i,j,k) < -small) {
                        qs = qpls;
                    } else {
                        qs = 0.5*(qmns+qpls);
                    }
                }

                a_edgeq_y(i,j,k,n) = qs;
            }

            a_fy(i,j,k,n) = a_edgeq_y(i,j,k,n) * a_vmac(i,j,k);
        });
    }
    else
    {
        amrex::ParallelFor(vbx, a_ncomp, [a_q,a_vmac,small,a_fy,known_edgestate,a_edgeq_y]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qs;
            if ( known_edgestate == 0 ) {
                Real qpls = a_q(i,j  ,k,n) - 0.5 * iamr_yslope(i,j  ,k,n,a_q);
                Real qmns = a_q(i,j-1,k,n) + 0.5 * iamr_yslope(i,j-1,k,n,a_q);

                if (a_vmac(i,j,k) > small) {
                    qs = qmns;
                } else if (a_vmac(i,j,k) < -small) {
                    qs = qpls;
                } else {
                    qs = 0.5*(qmns+qpls);
                }

                a_edgeq_y(i,j,k,n) = qs;
            }

            a_fy(i,j,k,n) = a_edgeq_y(i,j,k,n) * a_vmac(i,j,k);
        });
    }


#if ( AMREX_SPACEDIM ==3 )

    extdir_lohi = has_extdir(a_bcs.dataPtr(), a_ncomp, 2);
    has_extdir_lo = extdir_lohi.first;
    has_extdir_hi = extdir_lohi.second;
    if ((has_extdir_lo and domain_klo >= wbx.smallEnd(2)-1) or
        (has_extdir_hi and domain_khi <= wbx.bigEnd(2)))
    {
        amrex::ParallelFor(wbx, a_ncomp, [d_bcrec,a_q,domain_klo,domain_khi,a_wmac,small,a_fz, a_edgeq_z,known_edgestate]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qs;

            if ( known_edgestate == 0 ) {
                bool extdir_klo = d_bcrec[n].lo(2) == BCType::ext_dir;
                bool extdir_khi = d_bcrec[n].hi(2) == BCType::ext_dir;

                if (extdir_klo and k <= domain_klo) {
                    qs = a_q(i,j,domain_klo-1,n);
                } else if (extdir_khi and k >= domain_khi+1) {
                    qs = a_q(i,j,domain_khi+1,n);
                } else {
                    Real qpls = a_q(i,j,k,n) - 0.5 * iamr_zslope_extdir
                        (i,j,k,n,a_q, extdir_klo, extdir_khi, domain_klo, domain_khi);
                    Real qmns = a_q(i,j,k-1,n) + 0.5 * iamr_zslope_extdir(
                        i,j,k-1,n,a_q, extdir_klo, extdir_khi, domain_klo, domain_khi);
                    if (a_wmac(i,j,k) > small) {
                        qs = qmns;
                    } else if (a_wmac(i,j,k) < -small) {
                        qs = qpls;
                    } else {
                        qs = 0.5*(qmns+qpls);
                    }
                }
                a_edgeq_z(i,j,k,n) = qs;
            }
            a_fz(i,j,k,n) = a_edgeq_z(i,j,k,n) * a_wmac(i,j,k);
        });
    }
    else
    {
        amrex::ParallelFor(wbx, a_ncomp, [a_q,a_wmac,small,a_fz, a_edgeq_z, known_edgestate]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qs;

            if ( known_edgestate == 0 ) {
                Real qpls = a_q(i,j,k  ,n) - 0.5 * iamr_zslope(i,j,k  ,n,a_q);
                Real qmns = a_q(i,j,k-1,n) + 0.5 * iamr_zslope(i,j,k-1,n,a_q);

                if (a_wmac(i,j,k) > small) {
                    qs = qmns;
                } else if (a_wmac(i,j,k) < -small) {
                    qs = qpls;
                } else {
                    qs = 0.5*(qmns+qpls);
                }

                a_edgeq_z(i,j,k,n) = qs;
            }
            a_fz(i,j,k,n) = a_edgeq_z(i,j,k,n) * a_wmac(i,j,k);
        });
    }

#endif
}

//
// Compute fluxes on given EB box
//
void
ComputeFluxesOnEBBox ( Box const& a_bx,
                       D_DECL( Array4<Real> const& a_fx,
                               Array4<Real> const& a_fy,
                               Array4<Real> const& a_fz),
                       D_DECL( Array4<Real> const& a_edgeq_x,
                               Array4<Real> const& a_edgeq_y,
                               Array4<Real> const& a_edgeq_z),
                       Array4<Real> const& a_q,
                       const int a_comp,
                       const int a_ncomp,
                       D_DECL( Array4<Real const> const& a_umac,
                               Array4<Real const> const& a_vmac,
                               Array4<Real const> const& a_wmac),
                       const Box&       a_domain,
                       const Vector<BCRec>& a_bcs,
                       D_DECL( Array4<Real const> const& a_fcx,
                               Array4<Real const> const& a_fcy,
                               Array4<Real const> const& a_fcz),
                       Array4<Real const> const& a_ccc,
                       Array4<EBCellFlag const> const& a_flag,
                       int known_edgestate)
{
    constexpr Real small = 1.e-10;

    const int domain_ilo = a_domain.smallEnd(0);
    const int domain_ihi = a_domain.bigEnd(0);
    const int domain_jlo = a_domain.smallEnd(1);
    const int domain_jhi = a_domain.bigEnd(1);
#if (AMREX_SPACEDIM==3)
    const int domain_klo = a_domain.smallEnd(2);
    const int domain_khi = a_domain.bigEnd(2);
#endif

    D_TERM( const Box& ubx = amrex::surroundingNodes(a_bx,0);,
            const Box& vbx = amrex::surroundingNodes(a_bx,1);,
            const Box& wbx = amrex::surroundingNodes(a_bx,2););

    const auto bc = a_bcs.dataPtr();
    auto d_bcrec  = convertToDeviceVector(a_bcs);


    // ****************************************************************************
    // Predict to x-faces
    // ****************************************************************************

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    auto extdir_lohi = has_extdir(a_bcs.dataPtr(), a_ncomp, 0);
    bool has_extdir_lo = extdir_lohi.first;
    bool has_extdir_hi = extdir_lohi.second;

    if ((has_extdir_lo and domain_ilo >= ubx.smallEnd(0)-1) or
        (has_extdir_hi and domain_ihi <= ubx.bigEnd(0)))
    {
        amrex::ParallelFor(ubx, a_ncomp,
        [d_bcrec,domain_ilo,domain_jlo,domain_klo,domain_ihi,domain_jhi,domain_khi,
         a_q,a_ccc,a_fcx,a_flag,a_umac,small,a_fx, known_edgestate, a_edgeq_x]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           bool extdir_ilo = d_bcrec[n].lo(0) == BCType::ext_dir;
           bool extdir_ihi = d_bcrec[n].hi(0) == BCType::ext_dir;
           bool extdir_jlo = d_bcrec[n].lo(1) == BCType::ext_dir;
           bool extdir_jhi = d_bcrec[n].hi(1) == BCType::ext_dir;
           bool extdir_klo = d_bcrec[n].lo(2) == BCType::ext_dir;
           bool extdir_khi = d_bcrec[n].hi(2) == BCType::ext_dir;
           Real qs;

           if (a_flag(i,j,k).isConnected(-1,0,0))
           {
               if (extdir_ilo and i <= domain_ilo) {
                   qs = a_q(domain_ilo-1,j,k,n);
               } else if (extdir_ihi and i >= domain_ihi+1) {
                   qs = a_q(domain_ihi+1,j,k,n);
               } else {

                   Real yf = a_fcx(i,j,k,0); // local (y,z) of centroid of z-face we are extrapolating to
                   Real zf = a_fcx(i,j,k,1);

                   Real xc = a_ccc(i,j,k,0); // centroid of cell (i,j,k)
                   Real yc = a_ccc(i,j,k,1);
                   Real zc = a_ccc(i,j,k,2);

                   Real delta_x = 0.5 + xc;
                   Real delta_y = yf  - yc;
                   Real delta_z = zf  - zc;

                   Real cc_qmax = std::max(a_q(i,j,k,n),a_q(i-1,j,k,n));
                   Real cc_qmin = std::min(a_q(i,j,k,n),a_q(i-1,j,k,n));

                   // Compute slopes of component "n" of q
                   const auto& slopes_eb_hi = iamr_slopes_extdir_eb(i,j,k,n,a_q,a_ccc,a_flag,
                                              extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                              extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                              extdir_klo, extdir_khi, domain_klo, domain_khi);

                   Real qpls = a_q(i  ,j,k,n) - delta_x * slopes_eb_hi[0]
                                              + delta_y * slopes_eb_hi[1]
                                              + delta_z * slopes_eb_hi[2];

                   qpls = std::max(std::min(qpls, cc_qmax), cc_qmin);

                   xc = a_ccc(i-1,j,k,0); // centroid of cell (i-1,j,k)
                   yc = a_ccc(i-1,j,k,1);
                   zc = a_ccc(i-1,j,k,2);

                   delta_x = 0.5 - xc;
                   delta_y = yf  - yc;
                   delta_z = zf  - zc;

                   // Compute slopes of component "n" of q
                   const auto& slopes_eb_lo = iamr_slopes_extdir_eb(i-1,j,k,n,a_q,a_ccc,a_flag,
                                              extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                              extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                              extdir_klo, extdir_khi, domain_klo, domain_khi);

                   Real qmns = a_q(i-1,j,k,n) + delta_x * slopes_eb_lo[0]
                                              + delta_y * slopes_eb_lo[1]
                                              + delta_z * slopes_eb_lo[2];

                   qmns = std::max(std::min(qmns, cc_qmax), cc_qmin);

                    if (a_umac(i,j,k) > small) {
                        qs = qmns;
                    } else if (a_umac(i,j,k) < -small) {
                        qs = qpls;
                    } else {
                        qs = 0.5*(qmns+qpls);
                    }
               }

               a_fx(i,j,k,n) = a_umac(i,j,k) * qs;

           } else {
               a_fx(i,j,k,n) = 0.0;
           }
        });
    }
    else
    {
        amrex::ParallelFor(ubx, a_ncomp,
        [a_q,a_ccc,a_fcx,a_flag,a_umac,small,a_fx]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           Real qs;
           if (a_flag(i,j,k).isConnected(-1,0,0))
           {
               Real yf = a_fcx(i,j,k,0); // local (y,z) of centroid of z-face we are extrapolating to
               Real zf = a_fcx(i,j,k,1);

               Real xc = a_ccc(i,j,k,0); // centroid of cell (i,j,k)
               Real yc = a_ccc(i,j,k,1);
               Real zc = a_ccc(i,j,k,2);

               Real delta_x = 0.5 + xc;
               Real delta_y = yf  - yc;
               Real delta_z = zf  - zc;

               Real cc_qmax = std::max(a_q(i,j,k,n),a_q(i-1,j,k,n));
               Real cc_qmin = std::min(a_q(i,j,k,n),a_q(i-1,j,k,n));

               // Compute slopes of component "n" of q
               const auto& slopes_eb_hi = iamr_slopes_eb(i,j,k,n,a_q,a_ccc,a_flag);

               Real qpls = a_q(i  ,j,k,n) - delta_x * slopes_eb_hi[0]
                                          + delta_y * slopes_eb_hi[1]
                                          + delta_z * slopes_eb_hi[2];

               qpls = std::max(std::min(qpls, cc_qmax), cc_qmin);

               xc = a_ccc(i-1,j,k,0); // centroid of cell (i-1,j,k)
               yc = a_ccc(i-1,j,k,1);
               zc = a_ccc(i-1,j,k,2);

               delta_x = 0.5 - xc;
               delta_y = yf  - yc;
               delta_z = zf  - zc;

               // Compute slopes of component "n" of q
               const auto& slopes_eb_lo = iamr_slopes_eb(i-1,j,k,n,a_q,a_ccc,a_flag);

               Real qmns = a_q(i-1,j,k,n) + delta_x * slopes_eb_lo[0]
                                          + delta_y * slopes_eb_lo[1]
                                          + delta_z * slopes_eb_lo[2];

               qmns = std::max(std::min(qmns, cc_qmax), cc_qmin);

               if (a_umac(i,j,k) > small) {
                    qs = qmns;
                } else if (a_umac(i,j,k) < -small) {
                    qs = qpls;
                } else {
                    qs = 0.5*(qmns+qpls);
                }

                a_fx(i,j,k,n) = a_umac(i,j,k) * qs;

           } else {
               a_fx(i,j,k,n) = 0.0;
           }
        });
    }


    // ****************************************************************************
    // Predict to y-faces
    // ****************************************************************************
    extdir_lohi = has_extdir(a_bcs.dataPtr(), a_ncomp, 1);
    has_extdir_lo = extdir_lohi.first;
    has_extdir_hi = extdir_lohi.second;
    if ((has_extdir_lo and domain_jlo >= vbx.smallEnd(1)-1) or
        (has_extdir_hi and domain_jhi <= vbx.bigEnd(1)))
    {
        amrex::ParallelFor(vbx, a_ncomp,
        [d_bcrec,domain_ilo,domain_jlo,domain_klo,domain_ihi,domain_jhi,domain_khi,
         a_q,a_ccc,a_fcy,a_flag,a_vmac,small,a_fy]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qs;
            if (a_flag(i,j,k).isConnected(0,-1,0))
            {
                bool extdir_ilo = d_bcrec[n].lo(0) == BCType::ext_dir;
                bool extdir_ihi = d_bcrec[n].hi(0) == BCType::ext_dir;
                bool extdir_jlo = d_bcrec[n].lo(1) == BCType::ext_dir;
                bool extdir_jhi = d_bcrec[n].hi(1) == BCType::ext_dir;
                bool extdir_klo = d_bcrec[n].lo(2) == BCType::ext_dir;
                bool extdir_khi = d_bcrec[n].hi(2) == BCType::ext_dir;

                if (extdir_jlo and j <= domain_jlo) {
                    qs = a_q(i,domain_jlo-1,k,n);
                } else if (extdir_jhi and j >= domain_jhi+1) {
                    qs = a_q(i,domain_jhi+1,k,n);
                } else {

                   Real xf = a_fcy(i,j,k,0); // local (x,z) of centroid of z-face we are extrapolating to
                   Real zf = a_fcy(i,j,k,1);

                   Real xc = a_ccc(i,j,k,0); // centroid of cell (i,j,k)
                   Real yc = a_ccc(i,j,k,1);
                   Real zc = a_ccc(i,j,k,2);

                   Real delta_x = xf  - xc;
                   Real delta_y = 0.5 + yc;
                   Real delta_z = zf  - zc;

                   Real cc_qmax = std::max(a_q(i,j,k,n),a_q(i,j-1,k,n));
                   Real cc_qmin = std::min(a_q(i,j,k,n),a_q(i,j-1,k,n));

                   // Compute slopes of component "n" of q
                   const auto& slopes_eb_hi = iamr_slopes_extdir_eb(i,j,k,n,a_q,a_ccc,a_flag,
                                              extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                              extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                              extdir_klo, extdir_khi, domain_klo, domain_khi);

                   Real qpls = a_q(i,j  ,k,n) + delta_x * slopes_eb_hi[0]
                                              - delta_y * slopes_eb_hi[1]
                                              + delta_z * slopes_eb_hi[2];

                   qpls = std::max(std::min(qpls, cc_qmax), cc_qmin);

                   xc = a_ccc(i,j-1,k,0); // centroid of cell (i-1,j,k)
                   yc = a_ccc(i,j-1,k,1);
                   zc = a_ccc(i,j-1,k,2);

                   delta_x = xf  - xc;
                   delta_y = 0.5 - yc;
                   delta_z = zf  - zc;

                   // Compute slopes of component "n" of q
                   const auto& slopes_eb_lo = iamr_slopes_extdir_eb(i,j-1,k,n,a_q,a_ccc,a_flag,
                                              extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                              extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                              extdir_klo, extdir_khi, domain_klo, domain_khi);

                   Real qmns = a_q(i,j-1,k,n) + delta_x * slopes_eb_lo[0]
                                              + delta_y * slopes_eb_lo[1]
                                              + delta_z * slopes_eb_lo[2];

                   qmns = std::max(std::min(qmns, cc_qmax), cc_qmin);

                    if (a_vmac(i,j,k) > small) {
                        qs = qmns;
                    } else if (a_vmac(i,j,k) < -small) {
                        qs = qpls;
                    } else {
                        qs = 0.5*(qmns+qpls);
                    }
                }

                a_fy(i,j,k,n) = a_vmac(i,j,k) * qs;

           } else {
                a_fy(i,j,k,n) = 0.0;
           }
        });
    }
    else
    {
        amrex::ParallelFor(vbx, a_ncomp,
        [a_q,a_ccc,a_fcy,a_flag,a_vmac,small,a_fy]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qs;
            if (a_flag(i,j,k).isConnected(0,-1,0)) {

               Real xf = a_fcy(i,j,k,0); // local (x,z) of centroid of z-face we are extrapolating to
               Real zf = a_fcy(i,j,k,1);

               Real xc = a_ccc(i,j,k,0); // centroid of cell (i,j,k)
               Real yc = a_ccc(i,j,k,1);
               Real zc = a_ccc(i,j,k,2);

               Real delta_x = xf  - xc;
               Real delta_y = 0.5 + yc;
               Real delta_z = zf  - zc;

               Real cc_qmax = std::max(a_q(i,j,k,n),a_q(i,j-1,k,n));
               Real cc_qmin = std::min(a_q(i,j,k,n),a_q(i,j-1,k,n));

               // Compute slopes of component "n" of q
               const auto& slopes_eb_hi = iamr_slopes_eb(i,j,k,n,a_q,a_ccc,a_flag);

               Real qpls = a_q(i,j  ,k,n) + delta_x * slopes_eb_hi[0]
                                          - delta_y * slopes_eb_hi[1]
                                          + delta_z * slopes_eb_hi[2];

               qpls = std::max(std::min(qpls, cc_qmax), cc_qmin);

               xc = a_ccc(i,j-1,k,0); // centroid of cell (i-1,j,k)
               yc = a_ccc(i,j-1,k,1);
               zc = a_ccc(i,j-1,k,2);

               delta_x = xf  - xc;
               delta_y = 0.5 - yc;
               delta_z = zf  - zc;

               // Compute slopes of component "n" of q
               const auto& slopes_eb_lo = iamr_slopes_eb(i,j-1,k,n,a_q,a_ccc,a_flag);

               Real qmns = a_q(i,j-1,k,n) + delta_x * slopes_eb_lo[0]
                                          + delta_y * slopes_eb_lo[1]
                                          + delta_z * slopes_eb_lo[2];

               qmns = std::max(std::min(qmns, cc_qmax), cc_qmin);

               if (a_vmac(i,j,k) > small) {
                   qs = qmns;
               } else if (a_vmac(i,j,k) < -small) {
                   qs = qpls;
               } else {
                   qs = 0.5*(qmns+qpls);
               }

               a_fy(i,j,k,n) = a_vmac(i,j,k) * qs;

           } else {
               a_fy(i,j,k,n) = 0.0;
           }
        });
    }



    //
    // ===================== Z =====================
    //
#if ( AMREX_SPACEDIM == 3 )


    // ****************************************************************************
    // Predict to z-faces
    // ****************************************************************************

    extdir_lohi =  has_extdir(a_bcs.dataPtr(), a_ncomp, 2);
    has_extdir_lo = extdir_lohi.first;
    has_extdir_hi = extdir_lohi.second;
    if ((has_extdir_lo and domain_klo >= wbx.smallEnd(2)-1) or
        (has_extdir_hi and domain_khi <= wbx.bigEnd(2)))
    {
        amrex::ParallelFor(wbx, a_ncomp,
        [d_bcrec,domain_ilo,domain_jlo,domain_klo,domain_ihi,domain_jhi,domain_khi,
         a_q,a_ccc,a_fcz,a_flag,a_wmac,small,a_fz]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (a_flag(i,j,k).isConnected(0,0,-1)) {

                bool extdir_ilo = d_bcrec[n].lo(0) == BCType::ext_dir;
                bool extdir_ihi = d_bcrec[n].hi(0) == BCType::ext_dir;
                bool extdir_jlo = d_bcrec[n].lo(1) == BCType::ext_dir;
                bool extdir_jhi = d_bcrec[n].hi(1) == BCType::ext_dir;
                bool extdir_klo = d_bcrec[n].lo(2) == BCType::ext_dir;
                bool extdir_khi = d_bcrec[n].hi(2) == BCType::ext_dir;

                Real qs;
                if (extdir_klo and k <= domain_klo) {
                    qs = a_q(i,j,domain_klo-1,n);
                } else if (extdir_khi and k >= domain_khi+1) {
                    qs = a_q(i,j,domain_khi+1,n);
                } else {

                    Real xf = a_fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
                    Real yf = a_fcz(i,j,k,1);

                    Real xc = a_ccc(i,j,k,0); // centroid of cell (i,j,k)
                    Real yc = a_ccc(i,j,k,1);
                    Real zc = a_ccc(i,j,k,2);

                    Real delta_x = xf  - xc;
                    Real delta_y = yf  - yc;
                    Real delta_z = 0.5 + zc;

                    Real cc_qmax = std::max(a_q(i,j,k,n),a_q(i,j,k-1,n));
                    Real cc_qmin = std::min(a_q(i,j,k,n),a_q(i,j,k-1,n));

                    // Compute slopes of component "n" of q
                    const auto& slopes_eb_hi = iamr_slopes_extdir_eb(i,j,k,n,a_q,a_ccc,a_flag,
                                               extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                               extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                               extdir_klo, extdir_khi, domain_klo, domain_khi);

                    Real qpls = a_q(i,j,k  ,n) + delta_x * slopes_eb_hi[0]
                                               + delta_y * slopes_eb_hi[1]
                                               - delta_z * slopes_eb_hi[2];

                    qpls = std::max(std::min(qpls, cc_qmax), cc_qmin);

                    xc = a_ccc(i,j,k-1,0); // centroid of cell (i,j,k-1)
                    yc = a_ccc(i,j,k-1,1);
                    zc = a_ccc(i,j,k-1,2);

                    delta_x = xf  - xc;
                    delta_y = yf  - yc;
                    delta_z = 0.5 - zc;

                    // Compute slopes of component "n" of q
                    const auto& slopes_eb_lo = iamr_slopes_extdir_eb(i,j,k-1,n,a_q,a_ccc,a_flag,
                                               extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                               extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                               extdir_klo, extdir_khi, domain_klo, domain_khi);

                    Real qmns = a_q(i,j,k-1,n) + delta_x * slopes_eb_lo[0]
                                               + delta_y * slopes_eb_lo[1]
                                               + delta_z * slopes_eb_lo[2];

                    qmns = std::max(std::min(qmns, cc_qmax), cc_qmin);

                    if (a_wmac(i,j,k) > small) {
                        qs = qmns;
                    } else if (a_wmac(i,j,k) < -small) {
                        qs = qpls;
                    } else {
                        qs = 0.5*(qmns+qpls);
                    }
                }

                a_fz(i,j,k,n) = a_wmac(i,j,k) * qs;

           } else {
                a_fz(i,j,k,n) = 0.0;
           }
        });
    }
    else
    {
        amrex::ParallelFor(wbx, a_ncomp,
        [a_q,a_ccc,a_fcz,a_flag,a_wmac,small,a_fz]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (a_flag(i,j,k).isConnected(0,0,-1)) {
                Real qs;

                Real xf = a_fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
                Real yf = a_fcz(i,j,k,1);

                Real xc = a_ccc(i,j,k,0); // centroid of cell (i,j,k)
                Real yc = a_ccc(i,j,k,1);
                Real zc = a_ccc(i,j,k,2);

                Real delta_x = xf  - xc;
                Real delta_y = yf  - yc;
                Real delta_z = 0.5 + zc;

                Real cc_qmax = std::max(a_q(i,j,k,n),a_q(i,j,k-1,n));
                Real cc_qmin = std::min(a_q(i,j,k,n),a_q(i,j,k-1,n));

                // Compute slopes of component "n" of q
                const auto& slopes_eb_hi = iamr_slopes_eb(i,j,k,n,a_q,a_ccc,a_flag);

                Real qpls = a_q(i,j,k  ,n) + delta_x * slopes_eb_hi[0]
                                           + delta_y * slopes_eb_hi[1]
                                           - delta_z * slopes_eb_hi[2];

                qpls = std::max(std::min(qpls, cc_qmax), cc_qmin);

                xc = a_ccc(i,j,k-1,0); // centroid of cell (i,j,k-1)
                yc = a_ccc(i,j,k-1,1);
                zc = a_ccc(i,j,k-1,2);

                delta_x = xf  - xc;
                delta_y = yf  - yc;
                delta_z = 0.5 - zc;

                // Compute slopes of component "n" of q
                const auto& slopes_eb_lo = iamr_slopes_eb(i,j,k-1,n,a_q,a_ccc,a_flag);

                Real qmns = a_q(i,j,k-1,n) + delta_x * slopes_eb_lo[0]
                                           + delta_y * slopes_eb_lo[1]
                                           + delta_z * slopes_eb_lo[2];

                qmns = std::max(std::min(qmns, cc_qmax), cc_qmin);

                if (a_wmac(i,j,k) > small) {
                    qs = qmns;
                } else if (a_wmac(i,j,k) < -small) {
                    qs = qpls;
                } else {
                    qs = 0.5*(qmns+qpls);
                }

                a_fz(i,j,k,n) = a_wmac(i,j,k) * qs;

           } else {
                a_fz(i,j,k,n) = 0.0;
           }
        });
    }

#endif
}


} // End of namespace "fluxes"

//
// Compute the three components of the convection term
//
void
Godunov::ComputeFluxes(  D_DECL(MultiFab& a_fx,
                                MultiFab& a_fy,
                                MultiFab& a_fz),
                         D_DECL(MultiFab& edgestate_x,
                                MultiFab& edgestate_y,
                                MultiFab& edgestate_z),
                         MultiFab& a_state,
                         const int a_comp,
                         const int a_ncomp,
                         D_DECL( const MultiFab& a_xsl,
                                 const MultiFab& a_ysl,
                                 const MultiFab& a_zsl),
                         const int a_sl_comp,
                         D_DECL( const MultiFab& a_umac,
                                 const MultiFab& a_vmac,
                                 const MultiFab& a_wmac),
                         const Geometry& a_geom,
                         const Vector<BCRec>& a_bcs,
                         int known_edgestate)
{

    AMREX_ALWAYS_ASSERT(a_state.hasEBFabFactory());
    AMREX_ALWAYS_ASSERT(a_state.ixType().cellCentered());
    AMREX_ALWAYS_ASSERT(a_bcs.size() == a_ncomp );

    // For now use 4 ghost nodes
    const int nghost(4);

    Box domain(a_geom.Domain());

    auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(a_state.Factory());

    // Get EB geometric info
    Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
    Array< const MultiCutFab*,AMREX_SPACEDIM> facecent;

    areafrac  =   ebfactory.getAreaFrac();
    facecent  =   ebfactory.getFaceCent();

    // Initialize fluxes
    D_TERM(a_fx.setVal(COVERED_VAL);,
           a_fy.setVal(COVERED_VAL);,
           a_fz.setVal(COVERED_VAL););

    for (MFIter mfi(a_state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox ();

        const EBFArrayBox&  state_fab = static_cast<EBFArrayBox const&>(a_state[mfi]);
        const EBCellFlagFab&    flags = state_fab.getEBCellFlagFab();
        Array4<EBCellFlag const> const& flag = flags.const_array();

        if (flags.getType(amrex::grow(bx,0)) != FabType::covered )
        {
            // No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
            if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular )
            {
                fluxes::ComputeFluxesOnBox( bx, D_DECL(a_fx.array(mfi), a_fy.array(mfi), a_fz.array(mfi)),
                                            D_DECL(edgestate_x.array(mfi), edgestate_y.array(mfi), edgestate_z.array(mfi)),
                                            a_state.array(mfi), a_comp, a_ncomp,
                                            D_DECL(a_umac.array(mfi), a_vmac.array(mfi), a_wmac.array(mfi)), domain, a_bcs, known_edgestate);

            }
            else
            {
                Array4<Real const> ccc = ebfactory.getCentroid().const_array(mfi);
                fluxes::ComputeFluxesOnEBBox(bx, D_DECL(a_fx.array(mfi), a_fy.array(mfi), a_fz.array(mfi)),
                                             D_DECL(edgestate_x.array(mfi), edgestate_y.array(mfi), edgestate_z.array(mfi)),
                                             a_state.array(mfi), a_comp, a_ncomp,
                                             D_DECL(a_umac.array(mfi), a_vmac.array(mfi), a_wmac.array(mfi)), domain, a_bcs,
                                             D_DECL(facecent[0]->array(mfi), facecent[1]->array(mfi), facecent[2]->array(mfi)),
                                             ccc, flag, known_edgestate);
            }
        }

    }

    // MR: incflo does not have this: should it be added?
    a_fx.FillBoundary(a_geom.periodicity());
    a_fy.FillBoundary(a_geom.periodicity());
#if ( AMREX_SPACEDIM == 3 )
    a_fz.FillBoundary(a_geom.periodicity());
#endif
}


//
// Compute the three components of the convection term -- MAC sync version
//
void
Godunov::ComputeSyncFluxes(  D_DECL(MultiFab& a_fx,
                                    MultiFab& a_fy,
                                    MultiFab& a_fz),
                             D_DECL(MultiFab& edgestate_x,
                                    MultiFab& edgestate_y,
                                    MultiFab& edgestate_z),
                             MultiFab& a_state,
                             const int a_comp,
                             const int a_ncomp,
                             D_DECL( const MultiFab& a_xsl,
                                     const MultiFab& a_ysl,
                                     const MultiFab& a_zsl),
                             const int a_sl_comp,
                             D_DECL( const MultiFab& a_umac,
                                     const MultiFab& a_vmac,
                                     const MultiFab& a_wmac),
                             D_DECL( const MultiFab& a_ucorr,
                                     const MultiFab& a_vcorr,
                                     const MultiFab& a_wcorr),
                             const Geometry& a_geom,
                             const Vector<BCRec>& a_bcs,
                             int known_edgestate)
{

    AMREX_ALWAYS_ASSERT(a_state.hasEBFabFactory());
    AMREX_ALWAYS_ASSERT(a_state.ixType().cellCentered());
    AMREX_ALWAYS_ASSERT(a_bcs.size() == a_ncomp );

    // For now use 4 ghost nodes
    const int nghost(4);

    Box domain(a_geom.Domain());

    auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(a_state.Factory());

    // Get EB geometric info
    Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
    Array< const MultiCutFab*,AMREX_SPACEDIM> facecent;

    areafrac  =   ebfactory.getAreaFrac();
    facecent  =   ebfactory.getFaceCent();

    //
    // For the time being we compute the fluxes in two steps
    //
    // 1) Compute fluxes using umac: this will give us the correct edge states
    // 2) Recompute the fluxes with ucorr and the edge state computed at 1)
    //
    // We could do 1) and 2) in one pass by modifying the flux functions but for now
    // we keep it simple and cross that bridge if we bump into efficiency issues

    //
    // Step 1: compute edge values using umac for upwinding -- only if we do not know edge state yet
    //
    if (!known_edgestate)
    {
        for (MFIter mfi(a_state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.tilebox ();

            const EBFArrayBox&  state_fab = static_cast<EBFArrayBox const&>(a_state[mfi]);
            const EBCellFlagFab&    flags = state_fab.getEBCellFlagFab();
            Array4<EBCellFlag const> const& flag = flags.const_array();

            if (flags.getType(amrex::grow(bx,0)) != FabType::covered )
            {
                // No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
                if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular )
                {
                    fluxes::ComputeFluxesOnBox( bx, D_DECL(a_fx.array(mfi), a_fy.array(mfi), a_fz.array(mfi)),
                                                D_DECL(edgestate_x.array(mfi), edgestate_y.array(mfi), edgestate_z.array(mfi)),
                                                a_state.array(mfi), a_comp, a_ncomp,
                                                D_DECL(a_umac.array(mfi), a_vmac.array(mfi), a_wmac.array(mfi)), domain, a_bcs, 0);

                }
                else
                {
                    Array4<Real const> ccc = ebfactory.getCentroid().const_array(mfi);
                    fluxes::ComputeFluxesOnEBBox(bx, D_DECL(a_fx.array(mfi), a_fy.array(mfi), a_fz.array(mfi)),
                                                 D_DECL(edgestate_x.array(mfi), edgestate_y.array(mfi), edgestate_z.array(mfi)),
                                                 a_state.array(mfi), a_comp, a_ncomp,
                                                 D_DECL(a_umac.array(mfi), a_vmac.array(mfi), a_wmac.array(mfi)), domain, a_bcs,
                                                 D_DECL(facecent[0]->array(mfi), facecent[1]->array(mfi), facecent[2]->array(mfi)),
                                                 ccc, flag, 0);
                }
            }

        }
    }


    // Initialize (or reset if know_edgestate=0) fluxes
    D_TERM(a_fx.setVal(COVERED_VAL);,
           a_fy.setVal(COVERED_VAL);,
           a_fz.setVal(COVERED_VAL););

    //
    // Step 2: re-advect with Ucorr
    //
    for (MFIter mfi(a_state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox ();

        const EBFArrayBox&  state_fab = static_cast<EBFArrayBox const&>(a_state[mfi]);
        const EBCellFlagFab&    flags = state_fab.getEBCellFlagFab();
        Array4<EBCellFlag const> const& flag = flags.const_array();

        if (flags.getType(amrex::grow(bx,0)) != FabType::covered )
        {
            // No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
            if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular )
            {
                fluxes::ComputeFluxesOnBox( bx, D_DECL(a_fx.array(mfi), a_fy.array(mfi), a_fz.array(mfi)),
                                            D_DECL(edgestate_x.array(mfi), edgestate_y.array(mfi), edgestate_z.array(mfi)),
                                            a_state.array(mfi), a_comp, a_ncomp,
                                            D_DECL(a_ucorr.array(mfi), a_vcorr.array(mfi), a_wcorr.array(mfi)), domain, a_bcs, 1);

            }
            else
            {
                Array4<Real const> ccc = ebfactory.getCentroid().const_array(mfi);
                fluxes::ComputeFluxesOnEBBox(bx, D_DECL(a_fx.array(mfi), a_fy.array(mfi), a_fz.array(mfi)),
                                             D_DECL(edgestate_x.array(mfi), edgestate_y.array(mfi), edgestate_z.array(mfi)),
                                             a_state.array(mfi), a_comp, a_ncomp,
                                             D_DECL(a_ucorr.array(mfi), a_vcorr.array(mfi), a_wcorr.array(mfi)), domain, a_bcs,
                                             D_DECL(facecent[0]->array(mfi), facecent[1]->array(mfi), facecent[2]->array(mfi)),
                                             ccc, flag, 1);
            }
        }

    }

    // MR: incflo does not have this: should it be added?
    a_fx.FillBoundary(a_geom.periodicity());
    a_fy.FillBoundary(a_geom.periodicity());
#if ( AMREX_SPACEDIM == 3 )
    a_fz.FillBoundary(a_geom.periodicity());
#endif
}
