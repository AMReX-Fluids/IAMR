#include <iamr_convection_K.H>
#include <Godunov.H>

using namespace amrex;

#ifdef AMREX_USE_EB
void Godunov::predict_vels_on_faces_eb ( Box const& a_ccbx,
                                         D_DECL( Box const& a_ubx,
                                                 Box const& a_vbx,
                                                 Box const& a_wbx ),
                                         D_DECL( Array4<Real> const& a_u,
                                                 Array4<Real> const& a_v,
                                                 Array4<Real> const& a_w ),
                                         Array4<Real const> const& a_vcc,
                                         Array4<EBCellFlag const> const& a_flag,
                                         D_DECL( Array4<Real const> const& a_fcx,
                                                 Array4<Real const> const& a_fcy,
                                                 Array4<Real const> const& a_fcz ),
                                         Array4<Real const> const& a_ccc,
                                         const Geometry&  a_geom,
                                         const Vector<BCRec>& a_bcs )
{
    constexpr Real small = 1.e-10;

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

    Real small_vel = 1.e-10;

    // ****************************************************************************
    // Predict to x-faces
    // ****************************************************************************
    if ((extdir_ilo and domain_ilo >= a_ubx.smallEnd(0)-1) or
        (extdir_ihi and domain_ihi <= a_ubx.bigEnd(0)))
    {
        amrex::ParallelFor(Box(a_ubx),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (a_flag(i,j,k).isConnected(-1,0,0))
            {
               Real yf = a_fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
               Real zf = a_fcx(i,j,k,1);

               Real xc = a_ccc(i,j,k,0); // centroid of cell (i,j,k)
               Real yc = a_ccc(i,j,k,1);
               Real zc = a_ccc(i,j,k,2);

               Real delta_x = 0.5 + xc;
               Real delta_y = yf  - yc;
               Real delta_z = zf  - zc;

               Real cc_umax = std::max(a_vcc(i,j,k,0), a_vcc(i-1,j,k,0));
               Real cc_umin = std::min(a_vcc(i,j,k,0), a_vcc(i-1,j,k,0));

               // Compute slopes of component "0" of a_vcc
               const auto& slopes_eb_hi = iamr_slopes_extdir_eb(i,j,k,0,a_vcc,a_ccc,a_flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real upls = a_vcc(i  ,j,k,0) - delta_x * slopes_eb_hi[0]
                                            + delta_y * slopes_eb_hi[1]
                                            + delta_z * slopes_eb_hi[2];

               upls = std::max(std::min(upls, cc_umax), cc_umin);

               xc = a_ccc(i-1,j,k,0); // centroid of cell (i-1,j,k)
               yc = a_ccc(i-1,j,k,1);
               zc = a_ccc(i-1,j,k,2);

               delta_x = 0.5 - xc;
               delta_y = yf  - yc;
               delta_z = zf  - zc;

               // Compute slopes of component "0" of a_vcc
               const auto& slopes_eb_lo = iamr_slopes_extdir_eb(i-1,j,k,0,a_vcc,a_ccc,a_flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real umns = a_vcc(i-1,j,k,0) + delta_x * slopes_eb_lo[0]
                                            + delta_y * slopes_eb_lo[1]
                                            + delta_z * slopes_eb_lo[2];

               umns = std::max(std::min(umns, cc_umax), cc_umin);

               if ( umns < 0.0 && upls > 0.0 ) {
                  a_u(i,j,k) = 0.0;
               } else {
                  Real avg = 0.5 * ( upls + umns );
                  if ( std::abs(avg) <  small_vel) { a_u(i,j,k) = 0.0;
                  } else if (avg >= 0)             { a_u(i,j,k) = umns;
                  } else                           { a_u(i,j,k) = upls;
                  }
               }

               if (extdir_ilo and i == domain_ilo) {
                   a_u(i,j,k) = a_vcc(i-1,j,k,0);
               } else if (extdir_ihi and i == domain_ihi+1) {
                   a_u(i,j,k) = a_vcc(i,j,k,0);
               }

            } else {
               a_u(i,j,k) = 0.0;
            }
        });
    }
    else
    {
        amrex::ParallelFor(Box(a_ubx),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (a_flag(i,j,k).isConnected(-1,0,0))
            {
               Real yf = a_fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
               Real zf = a_fcx(i,j,k,1);

               Real xc = a_ccc(i,j,k,0); // centroid of cell (i,j,k)
               Real yc = a_ccc(i,j,k,1);
               Real zc = a_ccc(i,j,k,2);

               Real delta_x = 0.5 + xc;
               Real delta_y = yf  - yc;
               Real delta_z = zf  - zc;

               Real cc_umax = std::max(a_vcc(i,j,k,0), a_vcc(i-1,j,k,0));
               Real cc_umin = std::min(a_vcc(i,j,k,0), a_vcc(i-1,j,k,0));

               // Compute slopes of component "0" of a_vcc
               const auto slopes_eb_hi = iamr_slopes_eb(i,j,k,0,a_vcc,a_ccc,a_flag);

               Real upls = a_vcc(i  ,j,k,0) - delta_x * slopes_eb_hi[0]
                                            + delta_y * slopes_eb_hi[1]
                                            + delta_z * slopes_eb_hi[2];

               upls = std::max(std::min(upls, cc_umax), cc_umin);

               xc = a_ccc(i-1,j,k,0); // centroid of cell (i-1,j,k)
               yc = a_ccc(i-1,j,k,1);
               zc = a_ccc(i-1,j,k,2);

               delta_x = 0.5 - xc;
               delta_y = yf  - yc;
               delta_z = zf  - zc;

               // Compute slopes of component "0" of a_vcc
               const auto& slopes_eb_lo = iamr_slopes_eb(i-1,j,k,0,a_vcc,a_ccc,a_flag);

               Real umns = a_vcc(i-1,j,k,0) + delta_x * slopes_eb_lo[0]
                                            + delta_y * slopes_eb_lo[1]
                                            + delta_z * slopes_eb_lo[2];

               umns = std::max(std::min(umns, cc_umax), cc_umin);

               if ( umns < 0.0 && upls > 0.0 ) {
                  a_u(i,j,k) = 0.0;
               } else {
                  Real avg = 0.5 * ( upls + umns );
                  if ( std::abs(avg) <  small_vel) { a_u(i,j,k) = 0.0;
                  } else if (avg >= 0)             { a_u(i,j,k) = umns;
                  } else                           { a_u(i,j,k) = upls;
                  }
               }

            } else {
               a_u(i,j,k) = 0.0;
            }
        });
    }

    // ****************************************************************************
    // Predict to y-faces
    // ****************************************************************************
    if ((extdir_jlo and domain_jlo >= a_vbx.smallEnd(1)-1) or
        (extdir_jhi and domain_jhi <= a_vbx.bigEnd(1)))
    {
        amrex::ParallelFor(Box(a_vbx),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (a_flag(i,j,k).isConnected(0,-1,0))
            {
               Real xf = a_fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
               Real zf = a_fcy(i,j,k,1);

               Real xc = a_ccc(i,j,k,0); // centroid of cell (i,j,k)
               Real yc = a_ccc(i,j,k,1);
               Real zc = a_ccc(i,j,k,2);

               Real delta_x = xf  - xc;
               Real delta_y = 0.5 + yc;
               Real delta_z = zf  - zc;

               Real cc_vmax = std::max(a_vcc(i,j,k,1), a_vcc(i,j-1,k,1));
               Real cc_vmin = std::min(a_vcc(i,j,k,1), a_vcc(i,j-1,k,1));

               // Compute slopes of component "1" of a_vcc
               const auto& slopes_eb_hi = iamr_slopes_extdir_eb(i,j,k,1,a_vcc,a_ccc,a_flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real vpls = a_vcc(i,j  ,k,1) + delta_x * slopes_eb_hi[0]
                                            - delta_y * slopes_eb_hi[1]
                                            + delta_z * slopes_eb_hi[2];

               vpls = std::max(std::min(vpls, cc_vmax), cc_vmin);

               xc = a_ccc(i,j-1,k,0); // centroid of cell (i,j-1,k)
               yc = a_ccc(i,j-1,k,1);
               zc = a_ccc(i,j-1,k,2);

               delta_x = xf  - xc;
               delta_y = 0.5 - yc;
               delta_z = zf  - zc;

               // Compute slopes of component "1" of a_vcc
               const auto& slopes_eb_lo = iamr_slopes_extdir_eb(i,j-1,k,1,a_vcc,a_ccc,a_flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real vmns = a_vcc(i,j-1,k,1) + delta_x * slopes_eb_lo[0]
                                          + delta_y * slopes_eb_lo[1]
                                          + delta_z * slopes_eb_lo[2];

               vmns = std::max(std::min(vmns, cc_vmax), cc_vmin);

               if ( vmns < 0.0 && vpls > 0.0 ) {
                  a_v(i,j,k) = 0.0;
               } else {
                  Real avg = 0.5 * ( vpls + vmns );
                  if ( std::abs(avg) <  small_vel) { a_v(i,j,k) = 0.0;
                  } else if (avg >= 0)             { a_v(i,j,k) = vmns;
                  } else                           { a_v(i,j,k) = vpls;
                  }
               }

               if (extdir_jlo and j == domain_jlo) {
                   a_v(i,j,k) = a_vcc(i,j-1,k,1);
               } else if (extdir_jhi and j == domain_jhi+1) {
                   a_v(i,j,k) = a_vcc(i,j,k,1);
               }

            } else {
               a_v(i,j,k) = 0.0;
            }
        });
    }
    else
    {
        amrex::ParallelFor(Box(a_vbx),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (a_flag(i,j,k).isConnected(0,-1,0))
            {
               Real xf = a_fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
               Real zf = a_fcy(i,j,k,1);

               Real xc = a_ccc(i,j,k,0); // centroid of cell (i,j,k)
               Real yc = a_ccc(i,j,k,1);
               Real zc = a_ccc(i,j,k,2);

               Real delta_x = xf  - xc;
               Real delta_y = 0.5 + yc;
               Real delta_z = zf  - zc;

               Real cc_vmax = std::max(a_vcc(i,j,k,1), a_vcc(i,j-1,k,1));
               Real cc_vmin = std::min(a_vcc(i,j,k,1), a_vcc(i,j-1,k,1));

               // Compute slopes of component "1" of a_vcc
               const auto slopes_eb_hi = iamr_slopes_eb(i,j,k,1,a_vcc,a_ccc,a_flag);

               Real vpls = a_vcc(i,j  ,k,1) + delta_x * slopes_eb_hi[0]
                                            - delta_y * slopes_eb_hi[1]
                                            + delta_z * slopes_eb_hi[2];

               vpls = std::max(std::min(vpls, cc_vmax), cc_vmin);

               xc = a_ccc(i,j-1,k,0); // centroid of cell (i,j-1,k)
               yc = a_ccc(i,j-1,k,1);
               zc = a_ccc(i,j-1,k,2);

               delta_x = xf  - xc;
               delta_y = 0.5 - yc;
               delta_z = zf  - zc;

               // Compute slopes of component "1" of a_vcc
               const auto& slopes_eb_lo = iamr_slopes_eb(i,j-1,k,1,a_vcc,a_ccc,a_flag);

               Real vmns = a_vcc(i,j-1,k,1) + delta_x * slopes_eb_lo[0]
                                            + delta_y * slopes_eb_lo[1]
                                            + delta_z * slopes_eb_lo[2];

               vmns = std::max(std::min(vmns, cc_vmax), cc_vmin);

               if ( vmns < 0.0 && vpls > 0.0 ) {
                  a_v(i,j,k) = 0.0;
               } else {
                  Real avg = 0.5 * ( vpls + vmns );
                  if ( std::abs(avg) <  small_vel) { a_v(i,j,k) = 0.0;
                  } else if (avg >= 0)             { a_v(i,j,k) = vmns;
                  } else                           { a_v(i,j,k) = vpls;
                  }
               }

            } else {
               a_v(i,j,k) = 0.0;
            }
        });
    }
#if (AMREX_SPACEDIM==3)
    // ****************************************************************************
    // Predict to z-faces
    // ****************************************************************************
    if ((extdir_klo and domain_klo >= a_wbx.smallEnd(2)-1) or
        (extdir_khi and domain_khi <= a_wbx.bigEnd(2)))
    {
        amrex::ParallelFor(Box(a_wbx),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (a_flag(i,j,k).isConnected(0,0,-1))
            {
               Real xf = a_fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
               Real yf = a_fcz(i,j,k,1);

               Real xc = a_ccc(i,j,k,0); // centroid of cell (i,j,k)
               Real yc = a_ccc(i,j,k,1);
               Real zc = a_ccc(i,j,k,2);

               Real delta_x = xf  - xc;
               Real delta_y = yf  - yc;
               Real delta_z = 0.5 + zc;

               Real cc_wmax = std::max(a_vcc(i,j,k,2), a_vcc(i,j,k-1,2));
               Real cc_wmin = std::min(a_vcc(i,j,k,2), a_vcc(i,j,k-1,2));

               // Compute slopes of component "2" of a_vcc
               const auto& slopes_eb_hi = iamr_slopes_extdir_eb(i,j,k,2,a_vcc,a_ccc,a_flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real wpls = a_vcc(i,j,k  ,2) + delta_x * slopes_eb_hi[0]
                                            + delta_y * slopes_eb_hi[1]
                                            - delta_z * slopes_eb_hi[2];

               wpls = std::max(std::min(wpls, cc_wmax), cc_wmin);

               xc = a_ccc(i,j,k-1,0); // centroid of cell (i,j,k-1)
               yc = a_ccc(i,j,k-1,1);
               zc = a_ccc(i,j,k-1,2);

               delta_x = xf  - xc;
               delta_y = yf  - yc;
               delta_z = 0.5 - zc;

               // Compute slopes of component "2" of a_vcc
               const auto& slopes_eb_lo = iamr_slopes_extdir_eb(i,j,k-1,2,a_vcc,a_ccc,a_flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real wmns = a_vcc(i,j,k-1,2) + delta_x * slopes_eb_lo[0]
                                            + delta_y * slopes_eb_lo[1]
                                            + delta_z * slopes_eb_lo[2];

               wmns = std::max(std::min(wmns, cc_wmax), cc_wmin);

               if ( wmns < 0.0 && wpls > 0.0 ) {
                  a_w(i,j,k) = 0.0;
               } else {
                  Real avg = 0.5 * ( wpls + wmns );
                  if ( std::abs(avg) <  small_vel) { a_w(i,j,k) = 0.0;
                  } else if (avg >= 0)             { a_w(i,j,k) = wmns;
                  } else                           { a_w(i,j,k) = wpls;
                  }
               }

                if (extdir_klo and k == domain_klo) {
                    a_w(i,j,k) = a_vcc(i,j,k-1,2);
                } else if (extdir_khi and k == domain_khi+1) {
                    a_w(i,j,k) = a_vcc(i,j,k,2);
                }

            } else {
               a_w(i,j,k) = 0.0;
            }
        });
    }
    else
    {
        amrex::ParallelFor(Box(a_wbx),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (a_flag(i,j,k).isConnected(0,0,-1))
            {
               Real xf = a_fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
               Real yf = a_fcz(i,j,k,1);

               Real xc = a_ccc(i,j,k,0); // centroid of cell (i,j,k)
               Real yc = a_ccc(i,j,k,1);
               Real zc = a_ccc(i,j,k,2);

               Real delta_x = xf  - xc;
               Real delta_y = yf  - yc;
               Real delta_z = 0.5 + zc;

               Real cc_wmax = std::max(a_vcc(i,j,k,2), a_vcc(i,j,k-1,2));
               Real cc_wmin = std::min(a_vcc(i,j,k,2), a_vcc(i,j,k-1,2));

               // Compute slopes of component "2" of a_vcc
               const auto slopes_eb_hi = iamr_slopes_eb(i,j,k,2,a_vcc,a_ccc,a_flag);

               Real wpls = a_vcc(i,j,k  ,2) + delta_x * slopes_eb_hi[0]
                                            + delta_y * slopes_eb_hi[1]
                                            - delta_z * slopes_eb_hi[2];

               wpls = std::max(std::min(wpls, cc_wmax), cc_wmin);

               xc = a_ccc(i,j,k-1,0); // centroid of cell (i,j,k-1)
               yc = a_ccc(i,j,k-1,1);
               zc = a_ccc(i,j,k-1,2);

               delta_x = xf  - xc;
               delta_y = yf  - yc;
               delta_z = 0.5 - zc;

               // Compute slopes of component "2" of a_vcc
               const auto& slopes_eb_lo = iamr_slopes_eb(i,j,k-1,2,a_vcc,a_ccc,a_flag);

               Real wmns = a_vcc(i,j,k-1,2) + delta_x * slopes_eb_lo[0]
                                            + delta_y * slopes_eb_lo[1]
                                            + delta_z * slopes_eb_lo[2];

               wmns = std::max(std::min(wmns, cc_wmax), cc_wmin);

               if ( wmns < 0.0 && wpls > 0.0 ) {
                  a_w(i,j,k) = 0.0;
               } else {
                  Real avg = 0.5 * ( wpls + wmns );
                  if ( std::abs(avg) <  small_vel) { a_w(i,j,k) = 0.0;
                  } else if (avg >= 0)             { a_w(i,j,k) = wmns;
                  } else                           { a_w(i,j,k) = wpls;
                  }
               }

            } else {
               a_w(i,j,k) = 0.0;
            }
        });
    }
#endif
}

// void Godunov::compute_convective_rate_eb (int lev, Box const& bx, int ncomp,
//                                          Array4<Real> const& dUdt,
//                                          Array4<Real const> const& fx,
//                                          Array4<Real const> const& fy,
//                                          Array4<Real const> const& fz,
//                                          Array4<EBCellFlag const> const& flag,
//                                          Array4<Real const> const& vfrac,
//                                          Array4<Real const> const& apx,
//                                          Array4<Real const> const& apy,
//                                          Array4<Real const> const& apz)
// {
//     const auto dxinv = Geom(lev).InvCellSizeArray();
//     const Box dbox = Geom(lev).growPeriodicDomain(2);
//     amrex::ParallelFor(bx, ncomp,
//     [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
//     {
//         if (!dbox.contains(IntVect(i,j,k)) or flag(i,j,k).isCovered()) {
//             dUdt(i,j,k,n) = 0.0;
//         } else if (flag(i,j,k).isRegular()) {
//             dUdt(i,j,k,n) = dxinv[0] * (fx(i,j,k,n) - fx(i+1,j,k,n))
//                 +           dxinv[1] * (fy(i,j,k,n) - fy(i,j+1,k,n))
//                 +           dxinv[2] * (fz(i,j,k,n) - fz(i,j,k+1,n));
//         } else {
//             dUdt(i,j,k,n) = (1.0/vfrac(i,j,k)) *
//                 ( dxinv[0] * (apx(i,j,k)*fx(i,j,k,n) - apx(i+1,j,k)*fx(i+1,j,k,n))
//                 + dxinv[1] * (apy(i,j,k)*fy(i,j,k,n) - apy(i,j+1,k)*fy(i,j+1,k,n))
//                 + dxinv[2] * (apz(i,j,k)*fz(i,j,k,n) - apz(i,j,k+1)*fz(i,j,k+1,n)) );
//         }
//     });
// }

// void Godunov::redistribute_eb (int lev, Box const& bx, int ncomp,
//                               Array4<Real> const& dUdt,
//                               Array4<Real const> const& dUdt_in,
//                               Array4<Real> const& scratch,
//                               Array4<EBCellFlag const> const& flag,
//                               Array4<Real const> const& vfrac)
// {
//     const Box dbox = Geom(lev).growPeriodicDomain(2);

//     Array4<Real> tmp(scratch, 0);
//     Array4<Real> delm(scratch, ncomp);
//     Array4<Real> wgt(scratch, 2*ncomp);

//     Box const& bxg1 = amrex::grow(bx,1);
//     Box const& bxg2 = amrex::grow(bx,2);

//     // xxxxx TODO: more weight options
//     amrex::ParallelFor(bxg2,
//     [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//     {
//         wgt(i,j,k) = (dbox.contains(IntVect(i,j,k))) ? 1.0 : 0.0;
//     });

//     amrex::ParallelFor(bxg1, ncomp,
//     [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
//     {
//         if (flag(i,j,k).isSingleValued()) {
//             Real vtot = 0.0;
//             Real divnc = 0.0;
//             for (int kk = -1; kk <= 1; ++kk) {
//             for (int jj = -1; jj <= 1; ++jj) {
//             for (int ii = -1; ii <= 1; ++ii) {
//                 if ((ii != 0 or jj != 0 or kk != 0) and
//                     flag(i,j,k).isConnected(ii,jj,kk) and
//                     dbox.contains(IntVect(i+ii,j+jj,k+kk)))
//                 {
//                     Real vf = vfrac(i+ii,j+jj,k+kk);
//                     vtot += vf;
//                     divnc += vf * dUdt_in(i+ii,j+jj,k+kk,n);
//                 }
//             }}}
//             divnc /= (vtot + 1.e-80);
//             Real optmp = (1.0-vfrac(i,j,k))*(divnc-dUdt_in(i,j,k,n));
//             tmp(i,j,k,n) = optmp;
//             delm(i,j,k,n) = -vfrac(i,j,k)*optmp;
//         } else {
//             tmp(i,j,k,n) = 0.0;
//         }
//     });

//     amrex::ParallelFor(bxg1 & dbox, ncomp,
//     [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
//     {
//         if (flag(i,j,k).isSingleValued()) {
//             Real wtot = 0.0;
//             for (int kk = -1; kk <= 1; ++kk) {
//             for (int jj = -1; jj <= 1; ++jj) {
//             for (int ii = -1; ii <= 1; ++ii) {
//                 if ((ii != 0 or jj != 0 or kk != 0) and
//                     flag(i,j,k).isConnected(ii,jj,kk))
//                 {
//                     wtot += vfrac(i+ii,j+jj,k+kk) * wgt(i+ii,j+jj,k+kk);
//                 }
//             }}}
//             wtot = 1.0/(wtot+1.e-80);

//             Real dtmp = delm(i,j,k,n) * wtot;
//             for (int kk = -1; kk <= 1; ++kk) {
//             for (int jj = -1; jj <= 1; ++jj) {
//             for (int ii = -1; ii <= 1; ++ii) {
//                 if ((ii != 0 or jj != 0 or kk != 0) and
//                     bx.contains(IntVect(i+ii,j+jj,k+kk)) and
//                     flag(i,j,k).isConnected(ii,jj,kk))
//                 {
//                     Gpu::Atomic::Add(&tmp(i+ii,j+jj,k+kk,n), dtmp*wgt(i+ii,j+jj,k+kk));
//                 }
//             }}}
//         }
//     });

//     amrex::ParallelFor(bx, ncomp,
//     [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
//     {
//         dUdt(i,j,k,n) = dUdt_in(i,j,k,n) + tmp(i,j,k,n);
//     });
// }
#endif
