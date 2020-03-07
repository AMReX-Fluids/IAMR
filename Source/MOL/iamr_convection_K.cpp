#include <Godunov.H>
#include <iamr_slopes_K.H>
#include <iamr_constants.H>

using namespace amrex;

void Godunov::predict_vels_on_faces ( D_DECL( Box const& a_ubx,
                                              Box const& a_vbx,
                                              Box const& a_wbx) ,
                                      D_DECL( Array4<Real> const& a_u,
                                              Array4<Real> const& a_v,
                                              Array4<Real> const& a_w),
                                      Array4<Real const> const& a_vcc,
                                      const Geometry&  a_geom,
                                      const Vector<BCRec>& a_bcs )
{

    const Box& domain_box = a_geom.Domain();
    const int  domain_ilo = domain_box.smallEnd(0);
    const int  domain_ihi = domain_box.bigEnd(0);
    const int  domain_jlo = domain_box.smallEnd(1);
    const int  domain_jhi = domain_box.bigEnd(1);
#if (AMREX_SPACEDIM==3)
    const int  domain_klo = domain_box.smallEnd(2);
    const int  domain_khi = domain_box.bigEnd(2);
#endif

    const auto bc = a_bcs.dataPtr();
    bool extdir_ilo = (bc[0].lo(0) == BCType::ext_dir);
    bool extdir_ihi = (bc[0].hi(0) == BCType::ext_dir);
    bool extdir_jlo = (bc[1].lo(1) == BCType::ext_dir);
    bool extdir_jhi = (bc[1].lo(1) == BCType::ext_dir);
#if (AMREX_SPACEDIM==3)
    bool extdir_klo = (bc[2].lo(2) == BCType::ext_dir);
    bool extdir_khi = (bc[2].lo(2) == BCType::ext_dir);
#endif

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    if ((extdir_ilo and domain_ilo >= a_ubx.smallEnd(0)-1) or
        (extdir_ihi and domain_ihi <= a_ubx.bigEnd(0)))
    {
        amrex::ParallelFor(a_ubx, [a_vcc,extdir_ilo,extdir_ihi,domain_ilo,domain_ihi,a_u]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real upls = a_vcc(i,j,k,0) - 0.5 * iamr_xslope_extdir
                (i,j,k,0,a_vcc, extdir_ilo, extdir_ihi, domain_ilo, domain_ihi);
            Real umns = a_vcc(i-1,j,k,0) + 0.5 * iamr_xslope_extdir
                (i-1,j,k,0,a_vcc, extdir_ilo, extdir_ihi, domain_ilo, domain_ihi);
            if (umns < 0.0 and upls > 0.0) {
                a_u(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (upls + umns);
                if (std::abs(avg) < small_vel) {
                    a_u(i,j,k) = 0.0;
                } else if (avg > 0.0) {
                    a_u(i,j,k) = umns;
                } else {
                    a_u(i,j,k) = upls;
                }
            }

            if (extdir_ilo and i == domain_ilo) {
                a_u(i,j,k) = a_vcc(i-1,j,k,0);
            } else if (extdir_ihi and i == domain_ihi+1) {
                a_u(i,j,k) = a_vcc(i,j,k,0);
            }
        });
    }
    else
    {
        amrex::ParallelFor(a_ubx, [a_vcc,a_u]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real upls = a_vcc(i  ,j,k,0) - 0.5 * iamr_xslope(i  ,j,k,0,a_vcc);
            Real umns = a_vcc(i-1,j,k,0) + 0.5 * iamr_xslope(i-1,j,k,0,a_vcc);
            if (umns < 0.0 and upls > 0.0) {
                a_u(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (upls + umns);
                if (std::abs(avg) < small_vel) {
                    a_u(i,j,k) = 0.0;
                } else if (avg > 0.0) {
                    a_u(i,j,k) = umns;
                } else {
                    a_u(i,j,k) = upls;
                }
            }
        });
    }

    if ((extdir_jlo and domain_jlo >= a_vbx.smallEnd(1)-1) or
        (extdir_jhi and domain_jhi <= a_vbx.bigEnd(1)))
    {
        amrex::ParallelFor(a_vbx, [a_vcc,extdir_jlo,extdir_jhi,domain_jlo,domain_jhi,a_v]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real vpls = a_vcc(i,j,k,1) - 0.5 * iamr_yslope_extdir
                (i,j,k,1,a_vcc, extdir_jlo, extdir_jhi, domain_jlo, domain_jhi);
            Real vmns = a_vcc(i,j-1,k,1) + 0.5 * iamr_yslope_extdir
                (i,j-1,k,1,a_vcc, extdir_jlo, extdir_jhi, domain_jlo, domain_jhi);
            if (vmns < 0.0 and vpls > 0.0) {
                a_v(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (vpls + vmns);
                if (std::abs(avg) < small_vel) {
                    a_v(i,j,k) = 0.0;
                } else if (avg > 0.0) {
                    a_v(i,j,k) = vmns;
                } else {
                    a_v(i,j,k) = vpls;
                }
            }

            if (extdir_jlo and j == domain_jlo) {
                a_v(i,j,k) = a_vcc(i,j-1,k,1);
            } else if (extdir_jhi and j == domain_jhi+1) {
                a_v(i,j,k) = a_vcc(i,j,k,1);
            }
        });
    }
    else
    {
        amrex::ParallelFor(a_vbx, [a_vcc,a_v]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real vpls = a_vcc(i,j  ,k,1) - 0.5 * iamr_yslope(i,j  ,k,1,a_vcc);
            Real vmns = a_vcc(i,j-1,k,1) + 0.5 * iamr_yslope(i,j-1,k,1,a_vcc);
            if (vmns < 0.0 and vpls > 0.0) {
                a_v(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (vpls + vmns);
                if (std::abs(avg) < small_vel) {
                    a_v(i,j,k) = 0.0;
                } else if (avg > 0.0) {
                    a_v(i,j,k) = vmns;
                } else {
                    a_v(i,j,k) = vpls;
                }
            }
        });
    }

#if (AMREX_SPACEDIM==3)
    if ((extdir_klo and domain_klo >= a_wbx.smallEnd(2)-1) or
        (extdir_khi and domain_khi <= a_wbx.bigEnd(2)))
    {
        amrex::ParallelFor(a_wbx, [a_vcc,extdir_klo,extdir_khi,domain_klo,domain_khi,a_w]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real wpls = a_vcc(i,j,k,2) - 0.5 * iamr_zslope_extdir
                (i,j,k,2,a_vcc, extdir_klo, extdir_khi, domain_klo, domain_khi);
            Real wmns = a_vcc(i,j,k-1,2) + 0.5 * iamr_zslope_extdir(
                i,j,k-1,2,a_vcc, extdir_klo, extdir_khi, domain_klo, domain_khi);
            if (wmns < 0.0 and wpls > 0.0) {
                a_w(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (wpls + wmns);
                if (std::abs(avg) < small_vel) {
                    a_w(i,j,k) = 0.0;
                } else if (avg > 0.0) {
                    a_w(i,j,k) = wmns;
                } else {
                    a_w(i,j,k) = wpls;
                }
            }

            if (extdir_klo and k == domain_klo) {
                a_w(i,j,k) = a_vcc(i,j,k-1,2);
            } else if (extdir_khi and k == domain_khi+1) {
                a_w(i,j,k) = a_vcc(i,j,k,2);
            }
        });
    }
    else
    {
        amrex::ParallelFor(a_wbx, [a_vcc,a_w]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real wpls = a_vcc(i,j,k  ,2) - 0.5 * iamr_zslope(i,j,k  ,2,a_vcc);
            Real wmns = a_vcc(i,j,k-1,2) + 0.5 * iamr_zslope(i,j,k-1,2,a_vcc);
            if (wmns < 0.0 and wpls > 0.0) {
                a_w(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (wpls + wmns);
                if (std::abs(avg) < small_vel) {
                    a_w(i,j,k) = 0.0;
            } else if (avg > 0.0) {
                    a_w(i,j,k) = wmns;
                } else {
                    a_w(i,j,k) = wpls;
                }
            }
        });
    }
#endif
}

// void Godunov::compute_convective_rate (int lev, Box const& bx, int ncomp,
//                                        Array4<Real> const& dUdt,
//                                        Array4<Real const> const& fx,
//                                        Array4<Real const> const& fy,
//                                        Array4<Real const> const& fz)
// {
//     const auto dxinv = Geom(lev).InvCellSizeArray();
//     amrex::ParallelFor(bx, ncomp,
//     [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
//     {
//         dUdt(i,j,k,n) = dxinv[0] * (fx(i,j,k,n) - fx(i+1,j,k,n))
//             +           dxinv[1] * (fy(i,j,k,n) - fy(i,j+1,k,n))
//             +           dxinv[2] * (fz(i,j,k,n) - fz(i,j,k+1,n));
//     });
// }
