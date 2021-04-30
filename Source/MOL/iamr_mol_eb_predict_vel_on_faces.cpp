#include <iamr_constants.H>
#include <iamr_mol.H>
#include <AMReX_EB_slopes_K.H>

using namespace amrex;

namespace {
    std::pair<bool,bool> has_extdir_or_ho (BCRec const* bcrec, int n, int dir)
    {
        std::pair<bool,bool> r{false,false};
        r.first = r.first
             or (bcrec[n].lo(dir) == BCType::ext_dir)
             or (bcrec[n].lo(dir) == BCType::hoextrap);
        r.second = r.second
             or (bcrec[n].hi(dir) == BCType::ext_dir)
             or (bcrec[n].hi(dir) == BCType::hoextrap);
        return r;
    }
}

#ifdef AMREX_USE_EB
void
MOL::EB_PredictVelOnFaces ( AMREX_D_DECL( Box const& ubx,
                                          Box const& vbx,
                                          Box const& wbx ),
                            AMREX_D_DECL( Array4<Real> const& u,
                                          Array4<Real> const& v,
                                          Array4<Real> const& w ),
                            Array4<Real const> const& vcc,
                            Array4<EBCellFlag const> const& flag,
                            D_DECL( Array4<Real const> const& fcx,
                                    Array4<Real const> const& fcy,
                                    Array4<Real const> const& fcz ),
                            Array4<Real const> const& ccc,
                            Array4<Real const> const& vfrac,
                            const Geometry&  geom,
                            Vector<BCRec> const& h_bcrec,
                            const BCRec* d_bcrec )

{

    const Box& domain_box = geom.Domain();
    AMREX_D_TERM(
        const int domain_ilo = domain_box.smallEnd(0);
        const int domain_ihi = domain_box.bigEnd(0);,
        const int domain_jlo = domain_box.smallEnd(1);
        const int domain_jhi = domain_box.bigEnd(1);,
        const int domain_klo = domain_box.smallEnd(2);
        const int domain_khi = domain_box.bigEnd(2););

    int order = 2;

    // ****************************************************************************
    // Decide whether the stencil at each cell might need to see values that
    //     live on face centroids rather than cell centroids, i.e.
    //     are at a domain boundary with ext_dir or hoextrap boundary conditions
    // ****************************************************************************


    int n_for_xbc = 0;

    auto extdir_lohi_x_for_u = has_extdir_or_ho(h_bcrec.data(), n_for_xbc, static_cast<int>(Direction::x));
    bool has_extdir_or_ho_lo_x_for_u = extdir_lohi_x_for_u.first;
    bool has_extdir_or_ho_hi_x_for_u = extdir_lohi_x_for_u.second;

    auto extdir_lohi_y_for_u = has_extdir_or_ho(h_bcrec.data(), n_for_xbc, static_cast<int>(Direction::y));
    bool has_extdir_or_ho_lo_y_for_u = extdir_lohi_y_for_u.first;
    bool has_extdir_or_ho_hi_y_for_u = extdir_lohi_y_for_u.second;

#if (AMREX_SPACEDIM == 3)
    auto extdir_lohi_z_for_u = has_extdir_or_ho(h_bcrec.data(), n_for_xbc, static_cast<int>(Direction::z));
    bool has_extdir_or_ho_lo_z_for_u = extdir_lohi_z_for_u.first;
    bool has_extdir_or_ho_hi_z_for_u = extdir_lohi_z_for_u.second;
#endif



    // ****************************************************************************
    // Predict to x-faces
    // ****************************************************************************
    if ((has_extdir_or_ho_lo_x_for_u and domain_ilo >= ubx.smallEnd(0)-1) or
        (has_extdir_or_ho_hi_x_for_u and domain_ihi <= ubx.bigEnd(0)    ) or
        (has_extdir_or_ho_lo_y_for_u and domain_jlo >= vbx.smallEnd(1)-1) or
        (has_extdir_or_ho_hi_y_for_u and domain_jhi <= vbx.bigEnd(1)    )
#if (AMREX_SPACEDIM == 2)
        )
#elif (AMREX_SPACEDIM == 3)
        or
        (has_extdir_or_ho_lo_z_for_u and domain_klo >= wbx.smallEnd(2)-1) or
        (has_extdir_or_ho_hi_z_for_u and domain_khi <= wbx.bigEnd(2)    ) )
#endif
    {
        amrex::ParallelFor(Box(ubx),
        [u,vcc,flag,ccc,vfrac,d_bcrec,
         AMREX_D_DECL(fcx,fcy,fcz),
         AMREX_D_DECL(domain_ilo,domain_jlo,domain_klo),
         AMREX_D_DECL(domain_ihi,domain_jhi,domain_khi),
         order]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real u_val(0);

            AMREX_D_TERM(bool extdir_or_ho_ilo = (d_bcrec[0].lo(0) == BCType::ext_dir) or
                                                 (d_bcrec[0].lo(0) == BCType::hoextrap);,
                         bool extdir_or_ho_jlo = (d_bcrec[0].lo(1) == BCType::ext_dir) or
                                                 (d_bcrec[0].lo(1) == BCType::hoextrap);,
                         bool extdir_or_ho_klo = (d_bcrec[0].lo(2) == BCType::ext_dir) or
                                                 (d_bcrec[0].lo(2) == BCType::hoextrap););

            AMREX_D_TERM(bool extdir_or_ho_ihi = (d_bcrec[0].hi(0) == BCType::ext_dir) or
                                                 (d_bcrec[0].hi(0) == BCType::hoextrap);,
                         bool extdir_or_ho_jhi = (d_bcrec[0].hi(1) == BCType::ext_dir) or
                                                 (d_bcrec[0].hi(1) == BCType::hoextrap);,
                         bool extdir_or_ho_khi = (d_bcrec[0].hi(2) == BCType::ext_dir) or
                                                 (d_bcrec[0].hi(2) == BCType::hoextrap););

            if (flag(i,j,k).isConnected(-1,0,0))
            {
               Real yf = fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
               Real zf = fcx(i,j,k,1);
#endif
               AMREX_D_TERM(Real delta_x = 0.5 + ccc(i,j,k,0);,
                            Real delta_y = yf  - ccc(i,j,k,1);,
                            Real delta_z = zf  - ccc(i,j,k,2););

               Real vcc_mns = vcc(i-1,j,k,0);
               Real vcc_pls = vcc(i,j,k,0);

               Real cc_umax = amrex::max(vcc_pls, vcc_mns);
               Real cc_umin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "0" of vcc
               const auto& slopes_eb_hi = amrex_lim_slopes_extdir_eb(i,j,k,0,vcc,ccc,vfrac,
                                          AMREX_D_DECL(fcx,fcy,fcz), flag,
                                          AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                          AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                          AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                          AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                          order);

#if (AMREX_SPACEDIM == 3)
               Real upls = vcc_pls - delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1]
                                   + delta_z * slopes_eb_hi[2];
#else
               Real upls = vcc_pls - delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1];
#endif
               upls = amrex::max(amrex::min(upls, cc_umax), cc_umin);

               AMREX_D_TERM(delta_x = 0.5 - ccc(i-1,j,k,0);,
                            delta_y = yf  - ccc(i-1,j,k,1);,
                            delta_z = zf  - ccc(i-1,j,k,2););

               // Compute slopes of component "0" of vcc
               const auto& slopes_eb_lo = amrex_lim_slopes_extdir_eb(i-1,j,k,0,vcc,ccc,vfrac,
                                          AMREX_D_DECL(fcx,fcy,fcz), flag,
                                          AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                          AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                          AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                          AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                          order);

#if (AMREX_SPACEDIM == 3)
               Real umns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];
#else
               Real umns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1];
#endif
               umns = amrex::max(amrex::min(umns, cc_umax), cc_umin);

               if ( umns >= 0.0 or upls <= 0.0 ) {
                  Real avg = 0.5 * ( upls + umns );

                  if (avg >= small_vel) {
                    u_val = umns;
                  }
                  else if (avg <= -small_vel) {
                    u_val = upls;
                  }
               }

               if (i == domain_ilo && (d_bcrec[0].lo(0) == BCType::ext_dir)) {
                   u_val = vcc_mns;
               } else if (i == domain_ihi+1 && (d_bcrec[0].hi(0) == BCType::ext_dir)) {
                   u_val = vcc_pls;
               }
            }

            u(i,j,k) = u_val;
        });
    }
    else
    {
        amrex::ParallelFor(Box(ubx),
        [u,vcc,flag,AMREX_D_DECL(fcx,fcy,fcz),ccc,vfrac,order]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real u_val(0);

            if (flag(i,j,k).isConnected(-1,0,0))
            {
               Real yf = fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
               Real zf = fcx(i,j,k,1);
#endif

               AMREX_D_TERM(Real delta_x = 0.5 + ccc(i,j,k,0);,
                            Real delta_y = yf  - ccc(i,j,k,1);,
                            Real delta_z = zf  - ccc(i,j,k,2););

               const Real vcc_mns = vcc(i-1,j,k,0);
               const Real vcc_pls = vcc(i,j,k,0);

               Real cc_umax = amrex::max(vcc_pls, vcc_mns);
               Real cc_umin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "0" of vcc
               const auto slopes_eb_hi = amrex_lim_slopes_eb(i,j,k,0,vcc,ccc,vfrac,AMREX_D_DECL(fcx,fcy,fcz),flag,order);

#if (AMREX_SPACEDIM == 3)
               Real upls = vcc_pls - delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1]
                                   + delta_z * slopes_eb_hi[2];
#else
               Real upls = vcc_pls - delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1];
#endif
               upls = amrex::max(amrex::min(upls, cc_umax), cc_umin);

               AMREX_D_TERM(delta_x = 0.5 - ccc(i-1,j,k,0);,
                            delta_y = yf  - ccc(i-1,j,k,1);,
                            delta_z = zf  - ccc(i-1,j,k,2););

               // Compute slopes of component "0" of vcc
               const auto& slopes_eb_lo = amrex_lim_slopes_eb(i-1,j,k,0,vcc,ccc,vfrac,AMREX_D_DECL(fcx,fcy,fcz),flag,order);

#if (AMREX_SPACEDIM == 3)
               Real umns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];
#else
               Real umns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1];
#endif
               umns = amrex::max(amrex::min(umns, cc_umax), cc_umin);

               if ( umns >= 0.0 or upls <= 0.0 ) {
                  Real avg = 0.5 * ( upls + umns );

                  if (avg >= small_vel) {
                    u_val = umns;
                  }
                  else if (avg <= -small_vel) {
                    u_val = upls;
                  }
               }
            }

            u(i,j,k) = u_val;
        });
    }

    // ****************************************************************************
    // Predict to y-faces
    // ****************************************************************************

    int n_for_ybc = 1;

    auto extdir_lohi_x_for_v = has_extdir_or_ho(h_bcrec.data(), n_for_ybc, static_cast<int>(Direction::x));
    bool has_extdir_or_ho_lo_x_for_v = extdir_lohi_x_for_v.first;
    bool has_extdir_or_ho_hi_x_for_v = extdir_lohi_x_for_v.second;

    auto extdir_lohi_y_for_v = has_extdir_or_ho(h_bcrec.data(), n_for_ybc, static_cast<int>(Direction::y));
    bool has_extdir_or_ho_lo_y_for_v = extdir_lohi_y_for_v.first;
    bool has_extdir_or_ho_hi_y_for_v = extdir_lohi_y_for_v.second;

#if (AMREX_SPACEDIM == 3)
    auto extdir_lohi_z_for_v = has_extdir_or_ho(h_bcrec.data(), n_for_ybc, static_cast<int>(Direction::z));
    bool has_extdir_or_ho_lo_z_for_v = extdir_lohi_z_for_v.first;
    bool has_extdir_or_ho_hi_z_for_v = extdir_lohi_z_for_v.second;
#endif

    if ((has_extdir_or_ho_lo_x_for_v and domain_ilo >= ubx.smallEnd(0)-1) or
        (has_extdir_or_ho_hi_x_for_v and domain_ihi <= ubx.bigEnd(0)    ) or
        (has_extdir_or_ho_lo_y_for_v and domain_jlo >= vbx.smallEnd(1)-1) or
        (has_extdir_or_ho_hi_y_for_v and domain_jhi <= vbx.bigEnd(1)    )
#if (AMREX_SPACEDIM == 2)
        )
#elif (AMREX_SPACEDIM == 3)
        or
        (has_extdir_or_ho_lo_z_for_v and domain_klo >= wbx.smallEnd(2)-1) or
        (has_extdir_or_ho_hi_z_for_v and domain_khi <= wbx.bigEnd(2)    ) )
#endif
    {
        amrex::ParallelFor(Box(vbx),
        [v,vcc,flag,ccc,vfrac,d_bcrec,
         AMREX_D_DECL(fcx,fcy,fcz),
         AMREX_D_DECL(domain_ilo,domain_jlo,domain_klo),
         AMREX_D_DECL(domain_ihi,domain_jhi,domain_khi),
         order]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real v_val(0);

            AMREX_D_TERM(bool extdir_or_ho_ilo = (d_bcrec[1].lo(0) == BCType::ext_dir) or
                                                 (d_bcrec[1].lo(0) == BCType::hoextrap);,
                         bool extdir_or_ho_jlo = (d_bcrec[1].lo(1) == BCType::ext_dir) or
                                                 (d_bcrec[1].lo(1) == BCType::hoextrap);,
                         bool extdir_or_ho_klo = (d_bcrec[1].lo(2) == BCType::ext_dir) or
                                                 (d_bcrec[1].lo(2) == BCType::hoextrap););

            AMREX_D_TERM(bool extdir_or_ho_ihi = (d_bcrec[1].hi(0) == BCType::ext_dir) or
                                                 (d_bcrec[1].hi(0) == BCType::hoextrap);,
                         bool extdir_or_ho_jhi = (d_bcrec[1].hi(1) == BCType::ext_dir) or
                                                 (d_bcrec[1].hi(1) == BCType::hoextrap);,
                         bool extdir_or_ho_khi = (d_bcrec[1].hi(2) == BCType::ext_dir) or
                                                 (d_bcrec[1].hi(2) == BCType::hoextrap););

            if (flag(i,j,k).isConnected(0,-1,0))
            {
               Real xf = fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
               Real zf = fcy(i,j,k,1);
#endif

               AMREX_D_TERM(Real delta_x = xf  - ccc(i,j,k,0);,
                            Real delta_y = 0.5 + ccc(i,j,k,1);,
                            Real delta_z = zf  - ccc(i,j,k,2););

               const Real vcc_mns = vcc(i,j-1,k,1);
               const Real vcc_pls = vcc(i,j,k,1);

               Real cc_vmax = amrex::max(vcc_pls, vcc_mns);
               Real cc_vmin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "1" of vcc
               const auto& slopes_eb_hi = amrex_lim_slopes_extdir_eb(i,j,k,1,vcc,ccc,vfrac,
                                          AMREX_D_DECL(fcx,fcy,fcz), flag,
                                          AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                          AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                          AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                          AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                          order);

#if (AMREX_SPACEDIM == 3)
               Real vpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   - delta_y * slopes_eb_hi[1]
                                   + delta_z * slopes_eb_hi[2];
#else
               Real vpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   - delta_y * slopes_eb_hi[1];
#endif

               vpls = amrex::max(amrex::min(vpls, cc_vmax), cc_vmin);

               AMREX_D_TERM(delta_x = xf  - ccc(i,j-1,k,0);,
                            delta_y = 0.5 - ccc(i,j-1,k,1);,
                            delta_z = zf  - ccc(i,j-1,k,2););

               // Compute slopes of component "1" of vcc
               const auto& slopes_eb_lo = amrex_lim_slopes_extdir_eb(i,j-1,k,1,vcc,ccc,vfrac,
                                          AMREX_D_DECL(fcx,fcy,fcz), flag,
                                          AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                          AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                          AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                          AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                          order);

#if (AMREX_SPACEDIM == 3)
               Real vmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];
#else
               Real vmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1];
#endif

               vmns = amrex::max(amrex::min(vmns, cc_vmax), cc_vmin);

               if ( vmns >= 0.0 or vpls <= 0.0 ) {
                  Real avg = 0.5 * ( vpls + vmns );

                  if (avg >= small_vel) {
                    v_val = vmns;
                  }
                  else if (avg <= -small_vel) {
                    v_val = vpls;
                  }
               }

               if (j == domain_jlo && (d_bcrec[1].lo(1) == BCType::ext_dir)) {
                   v_val = vcc_mns;
               }
               else if (j == domain_jhi+1 && (d_bcrec[1].hi(1) == BCType::ext_dir)) {
                   v_val = vcc_pls;
               }
            }

            v(i,j,k) = v_val;
        });
    }
    else
    {
        amrex::ParallelFor(Box(vbx),
        [v,vcc,flag,AMREX_D_DECL(fcx,fcy,fcz),ccc,vfrac,order]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real v_val(0);

            if (flag(i,j,k).isConnected(0,-1,0))
            {
               Real xf = fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
               Real zf = fcy(i,j,k,1);
#endif

               AMREX_D_TERM(Real delta_x = xf  - ccc(i,j,k,0);,
                            Real delta_y = 0.5 + ccc(i,j,k,1);,
                            Real delta_z = zf  - ccc(i,j,k,2););

               const Real vcc_mns = vcc(i,j-1,k,1);
               const Real vcc_pls = vcc(i,j,k,1);

               Real cc_vmax = amrex::max(vcc_pls, vcc_mns);
               Real cc_vmin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "1" of vcc
               const auto slopes_eb_hi = amrex_lim_slopes_eb(i,j,k,1,vcc,ccc,vfrac,AMREX_D_DECL(fcx,fcy,fcz),flag,order);

#if (AMREX_SPACEDIM == 3)
               Real vpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   - delta_y * slopes_eb_hi[1]
                                   + delta_z * slopes_eb_hi[2];
#else
               Real vpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   - delta_y * slopes_eb_hi[1];
#endif

               vpls = amrex::max(amrex::min(vpls, cc_vmax), cc_vmin);

               AMREX_D_TERM(delta_x = xf  - ccc(i,j-1,k,0);,
                            delta_y = 0.5 - ccc(i,j-1,k,1);,
                            delta_z = zf  - ccc(i,j-1,k,2););

               // Compute slopes of component "1" of vcc
               const auto& slopes_eb_lo = amrex_lim_slopes_eb(i,j-1,k,1,vcc,ccc,vfrac,AMREX_D_DECL(fcx,fcy,fcz),flag,order);

#if (AMREX_SPACEDIM == 3)
               Real vmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];
#else
               Real vmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1];
#endif

               vmns = amrex::max(amrex::min(vmns, cc_vmax), cc_vmin);

               if ( vmns >= 0.0 or vpls <= 0.0 ) {
                  Real avg = 0.5 * ( vpls + vmns );

                  if (avg >= small_vel) {
                    v_val = vmns;
                  }
                  else if (avg <= -small_vel) {
                    v_val = vpls;
                  }
               }
            }

            v(i,j,k) = v_val;
        });
    }

#if (AMREX_SPACEDIM == 3)
    // ****************************************************************************
    // Predict to z-faces
    // ****************************************************************************

    int n_for_zbc = 2;

    auto extdir_lohi_x_for_w = has_extdir_or_ho(h_bcrec.data(), n_for_zbc, static_cast<int>(Direction::x));
    bool has_extdir_or_ho_lo_x_for_w = extdir_lohi_x_for_w.first;
    bool has_extdir_or_ho_hi_x_for_w = extdir_lohi_x_for_w.second;

    auto extdir_lohi_y_for_w = has_extdir_or_ho(h_bcrec.data(), n_for_zbc, static_cast<int>(Direction::y));
    bool has_extdir_or_ho_lo_y_for_w = extdir_lohi_y_for_w.first;
    bool has_extdir_or_ho_hi_y_for_w = extdir_lohi_y_for_w.second;

#if (AMREX_SPACEDIM == 3)
    auto extdir_lohi_z_for_w = has_extdir_or_ho(h_bcrec.data(), n_for_zbc, static_cast<int>(Direction::z));
    bool has_extdir_or_ho_lo_z_for_w = extdir_lohi_z_for_w.first;
    bool has_extdir_or_ho_hi_z_for_w = extdir_lohi_z_for_w.second;
#endif

    if ((has_extdir_or_ho_lo_x_for_w and domain_ilo >= ubx.smallEnd(0)-1) or
        (has_extdir_or_ho_hi_x_for_w and domain_ihi <= ubx.bigEnd(0)    ) or
        (has_extdir_or_ho_lo_y_for_w and domain_jlo >= vbx.smallEnd(1)-1) or
        (has_extdir_or_ho_hi_y_for_w and domain_jhi <= vbx.bigEnd(1)    )
#if (AMREX_SPACEDIM == 2)
        )
#elif (AMREX_SPACEDIM == 3)
        or
        (has_extdir_or_ho_lo_z_for_w and domain_klo >= wbx.smallEnd(2)-1) or
        (has_extdir_or_ho_hi_z_for_w and domain_khi <= wbx.bigEnd(2)    ) )
#endif
    {
        amrex::ParallelFor(Box(wbx),
        [w,vcc,flag,ccc,vfrac,d_bcrec,
         AMREX_D_DECL(fcx,fcy,fcz),
         domain_ilo,domain_ihi,domain_jlo,domain_jhi,domain_klo,domain_khi,
         order]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real w_val(0);

            bool extdir_or_ho_ilo = (d_bcrec[2].lo(0) == BCType::ext_dir) or
                                    (d_bcrec[2].lo(0) == BCType::hoextrap);
            bool extdir_or_ho_ihi = (d_bcrec[2].hi(0) == BCType::ext_dir) or
                                    (d_bcrec[2].hi(0) == BCType::hoextrap);

            bool extdir_or_ho_jlo = (d_bcrec[2].lo(1) == BCType::ext_dir) or
                                    (d_bcrec[2].lo(1) == BCType::hoextrap);
            bool extdir_or_ho_jhi = (d_bcrec[2].hi(1) == BCType::ext_dir) or
                                    (d_bcrec[2].hi(1) == BCType::hoextrap);

            bool extdir_or_ho_klo = (d_bcrec[2].lo(2) == BCType::ext_dir) or
                                    (d_bcrec[2].lo(2) == BCType::hoextrap);
            bool extdir_or_ho_khi = (d_bcrec[2].hi(2) == BCType::ext_dir) or
                                    (d_bcrec[2].hi(2) == BCType::hoextrap);

            if (flag(i,j,k).isConnected(0,0,-1))
            {
               Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
               Real yf = fcz(i,j,k,1);

               Real delta_x = xf  - ccc(i,j,k,0);
               Real delta_y = yf  - ccc(i,j,k,1);
               Real delta_z = 0.5 + ccc(i,j,k,2);

               const Real vcc_mns = vcc(i,j,k-1,2);
               const Real vcc_pls = vcc(i,j,k,2);

               Real cc_wmax = amrex::max(vcc_pls, vcc_mns);
               Real cc_wmin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "2" of vcc
               const auto& slopes_eb_hi = amrex_lim_slopes_extdir_eb(i,j,k,2,vcc,ccc,vfrac,
                                          AMREX_D_DECL(fcx,fcy,fcz), flag,
                                          AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                          AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                          AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                          AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                          order);

               Real wpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1]
                                   - delta_z * slopes_eb_hi[2];

               wpls = amrex::max(amrex::min(wpls, cc_wmax), cc_wmin);

               delta_x = xf  - ccc(i,j,k-1,0);
               delta_y = yf  - ccc(i,j,k-1,1);
               delta_z = 0.5 - ccc(i,j,k-1,2);

               // Compute slopes of component "2" of vcc
               const auto& slopes_eb_lo = amrex_lim_slopes_extdir_eb(i,j,k-1,2,vcc,ccc,vfrac,
                                          AMREX_D_DECL(fcx,fcy,fcz), flag,
                                          AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                          AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                          AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                          AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                          order);

               Real wmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];

               wmns = amrex::max(amrex::min(wmns, cc_wmax), cc_wmin);

               if ( wmns >= 0.0 or wpls <= 0.0 ) {
                  Real avg = 0.5 * ( wpls + wmns );

                  if (avg >= small_vel) {
                    w_val = wmns;
                  }
                  else if (avg <= -small_vel) {
                    w_val = wpls;
                  }
               }

                if (k == domain_klo && (d_bcrec[2].lo(2) == BCType::ext_dir)) {
                    w_val = vcc_mns;
                }
                else if (k == domain_khi+1 && (d_bcrec[2].hi(2) == BCType::ext_dir)) {
                    w_val = vcc_pls;
                }
            }

            w(i,j,k) = w_val;
        });
    }
    else
    {
        amrex::ParallelFor(Box(wbx),
        [w,vcc,flag,AMREX_D_DECL(fcx,fcy,fcz),ccc,vfrac,order]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real w_val(0);

            if (flag(i,j,k).isConnected(0,0,-1))
            {
               Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
               Real yf = fcz(i,j,k,1);

               Real delta_x = xf  - ccc(i,j,k,0);
               Real delta_y = yf  - ccc(i,j,k,1);
               Real delta_z = 0.5 + ccc(i,j,k,2);

               const Real vcc_mns = vcc(i,j,k-1,2);
               const Real vcc_pls = vcc(i,j,k,2);

               Real cc_wmax = amrex::max(vcc_pls, vcc_mns);
               Real cc_wmin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "2" of vcc
               const auto slopes_eb_hi = amrex_lim_slopes_eb(i,j,k,2,vcc,ccc,vfrac,AMREX_D_DECL(fcx,fcy,fcz),flag,order);

               Real wpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1]
                                   - delta_z * slopes_eb_hi[2];

               wpls = amrex::max(amrex::min(wpls, cc_wmax), cc_wmin);

               delta_x = xf  - ccc(i,j,k-1,0);
               delta_y = yf  - ccc(i,j,k-1,1);
               delta_z = 0.5 - ccc(i,j,k-1,2);

               // Compute slopes of component "2" of vcc
               const auto& slopes_eb_lo = amrex_lim_slopes_eb(i,j,k-1,2,vcc,ccc,vfrac,AMREX_D_DECL(fcx,fcy,fcz),flag,order);

               Real wmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];

               wmns = amrex::max(amrex::min(wmns, cc_wmax), cc_wmin);

               if ( wmns >= 0.0 or wpls <= 0.0 ) {
                  Real avg = 0.5 * ( wpls + wmns );

                  if (avg >= small_vel) {
                    w_val = wmns;
                  }
                  else if (avg <= -small_vel) {
                    w_val = wpls;
                  }
               }
            }

            w(i,j,k) = w_val;
        });
    }
#endif
}
#endif


// #ifdef AMREX_USE_EB
// void
// MOL::EB_PredictVelOnFaces ( Box const& a_ccbx,
//                             D_DECL( Box const& a_ubx,
//                                     Box const& a_vbx,
//                                     Box const& a_wbx ),
//                             D_DECL( Array4<Real> const& a_u,
//                                     Array4<Real> const& a_v,
//                                     Array4<Real> const& a_w ),
//                             Array4<Real const> const& a_vcc,
//                             Array4<EBCellFlag const> const& a_flag,
//                             D_DECL( Array4<Real const> const& a_fcx,
//                                     Array4<Real const> const& a_fcy,
//                                     Array4<Real const> const& a_fcz ),
//                             Array4<Real const> const& a_ccc,
//                             const Geometry&  a_geom,
//                             const BCRec* bc )
// {
//     const Box& domain_box = a_geom.Domain();
//     const int  domain_ilo = domain_box.smallEnd(0);
//     const int  domain_ihi = domain_box.bigEnd(0);
//     const int  domain_jlo = domain_box.smallEnd(1);
//     const int  domain_jhi = domain_box.bigEnd(1);
// #if (AMREX_SPACEDIM==3)
//     const int  domain_klo = domain_box.smallEnd(2);
//     const int  domain_khi = domain_box.bigEnd(2);
// #endif

//     bool extdir_or_ho_ilo = (bc[0].lo(0) == BCType::ext_dir) || (bc[0].lo(0) == BCType::hoextrap);
//     bool extdir_or_ho_ihi = (bc[0].hi(0) == BCType::ext_dir) || (bc[0].hi(0) == BCType::hoextrap);
//     bool extdir_or_ho_jlo = (bc[1].lo(1) == BCType::ext_dir) || (bc[1].lo(1) == BCType::hoextrap);
//     bool extdir_or_ho_jhi = (bc[1].lo(1) == BCType::ext_dir) || (bc[1].lo(1) == BCType::hoextrap);
// #if (AMREX_SPACEDIM==3)
//     bool extdir_or_ho_klo = (bc[2].lo(2) == BCType::ext_dir) || (bc[2].lo(2) == BCType::hoextrap);
//     bool extdir_or_ho_khi = (bc[2].lo(2) == BCType::ext_dir) || (bc[2].lo(2) == BCType::hoextrap);
// #endif

//     bool extdir_ilo = (bc[0].lo(0) == BCType::ext_dir);
//     bool extdir_ihi = (bc[0].hi(0) == BCType::ext_dir);
//     bool extdir_jlo = (bc[1].lo(1) == BCType::ext_dir);
//     bool extdir_jhi = (bc[1].lo(1) == BCType::ext_dir);
// #if (AMREX_SPACEDIM==3)
//     bool extdir_klo = (bc[2].lo(2) == BCType::ext_dir);
//     bool extdir_khi = (bc[2].lo(2) == BCType::ext_dir);
// #endif

//     // At an ext_dir boundary, the boundary value is on the face, not cell center.

//     // ****************************************************************************
//     // Predict to x-faces
//     // ****************************************************************************
//     if ((extdir_or_ho_ilo and domain_ilo >= a_ubx.smallEnd(0)-1) or
//         (extdir_or_ho_ihi and domain_ihi <= a_ubx.bigEnd(0)))
//     {
//         amrex::ParallelFor(Box(a_ubx),
//         [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//         {
//             if (a_flag(i,j,k).isConnected(-1,0,0))
//             {
// 	      if (extdir_ilo and i == domain_ilo) {
// 		  a_u(i,j,k) = a_vcc(i-1,j,k,0);
//                } else if (extdir_ihi and i == domain_ihi+1) {
//                    a_u(i,j,k) = a_vcc(i,j,k,0);
//                } else {

//                Real yf = a_fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
// #if (AMREX_SPACEDIM==3)
//                Real zf = a_fcx(i,j,k,1);
// #endif

//                D_TERM(Real xc = a_ccc(i,j,k,0);, // centroid of cell (i,j,k)
//                       Real yc = a_ccc(i,j,k,1);,
//                       Real zc = a_ccc(i,j,k,2););

//                D_TERM(Real delta_x = 0.5 + xc;,
//                       Real delta_y = yf  - yc;,
//                       Real delta_z = zf  - zc;);

//                Real cc_umax = amrex::max(a_vcc(i,j,k,0), a_vcc(i-1,j,k,0));
//                Real cc_umin = amrex::min(a_vcc(i,j,k,0), a_vcc(i-1,j,k,0));

//                // Compute slopes of component "0" of a_vcc
// 		const auto& slopes_eb_hi =
// 		  amrex_calc_slopes_extdir_eb(i,j,k, 0, a_vcc, a_ccc,
// 					      AMREX_D_DECL(a_fcx,a_fcy,a_fcz), a_flag,
// 					      AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
// 					      AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
// 					      AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
// 					      AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

// #if (AMREX_SPACEDIM==3)
//                Real upls = a_vcc(i  ,j,k,0) - delta_x * slopes_eb_hi[0]
//                                             + delta_y * slopes_eb_hi[1]
//                                             + delta_z * slopes_eb_hi[2];
// #else
//                Real upls = a_vcc(i  ,j,k,0) - delta_x * slopes_eb_hi[0]
//                                             + delta_y * slopes_eb_hi[1];
// #endif

//                upls = amrex::max(amrex::min(upls, cc_umax), cc_umin);

//                D_TERM(xc = a_ccc(i-1,j,k,0);, // centroid of cell (i-1,j,k)
//                       yc = a_ccc(i-1,j,k,1);,
//                       zc = a_ccc(i-1,j,k,2););

//                D_TERM(delta_x = 0.5 - xc;,
//                       delta_y = yf  - yc;,
//                       delta_z = zf  - zc;);

//                // Compute slopes of component "0" of a_vcc
// 		const auto& slopes_eb_lo =
// 		  amrex_calc_slopes_extdir_eb(i-1,j,k,0,a_vcc,a_ccc,
// 					      AMREX_D_DECL(a_fcx,a_fcy,a_fcz), a_flag,
// 					      AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
// 					      AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
// 					      AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
// 					      AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

// #if (AMREX_SPACEDIM==3)
//                Real umns = a_vcc(i-1,j,k,0) + delta_x * slopes_eb_lo[0]
//                                             + delta_y * slopes_eb_lo[1]
//                                             + delta_z * slopes_eb_lo[2];
// #else
//                Real umns = a_vcc(i-1,j,k,0) + delta_x * slopes_eb_lo[0]
//                                             + delta_y * slopes_eb_lo[1];
// #endif

//                umns = amrex::max(amrex::min(umns, cc_umax), cc_umin);

//                if ( umns < 0.0 && upls > 0.0 ) {
//                   a_u(i,j,k) = 0.0;
//                } else {
//                   Real avg = 0.5 * ( upls + umns );
//                   if ( std::abs(avg) <  small_vel) { a_u(i,j,k) = 0.0;
//                   } else if (avg >= 0)             { a_u(i,j,k) = umns;
//                   } else                           { a_u(i,j,k) = upls;
//                   }
//                }
//                }
//             } else {
//                a_u(i,j,k) = 0.0;
//             }
//         });
//     }
//     else
//     {
//         amrex::ParallelFor(Box(a_ubx),
//         [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//         {
//             if (a_flag(i,j,k).isConnected(-1,0,0))
//             {
//                Real yf = a_fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
// #if (AMREX_SPACEDIM==3)
//                Real zf = a_fcx(i,j,k,1);
// #endif

//                D_TERM(Real xc = a_ccc(i,j,k,0);, // centroid of cell (i,j,k)
//                       Real yc = a_ccc(i,j,k,1);,
//                       Real zc = a_ccc(i,j,k,2););

//                D_TERM(Real delta_x = 0.5 + xc;,
//                       Real delta_y = yf  - yc;,
//                       Real delta_z = zf  - zc;);

//                Real cc_umax = amrex::max(a_vcc(i,j,k,0), a_vcc(i-1,j,k,0));
//                Real cc_umin = amrex::min(a_vcc(i,j,k,0), a_vcc(i-1,j,k,0));

//                // Compute slopes of component "0" of a_vcc
//                const auto slopes_eb_hi = amrex_calc_slopes_eb(i,j,k,0,a_vcc,a_ccc,a_flag);

// #if (AMREX_SPACEDIM==3)
//                Real upls = a_vcc(i  ,j,k,0) - delta_x * slopes_eb_hi[0]
//                                             + delta_y * slopes_eb_hi[1]
//                                             + delta_z * slopes_eb_hi[2];
// #else
//                Real upls = a_vcc(i  ,j,k,0) - delta_x * slopes_eb_hi[0]
//                                             + delta_y * slopes_eb_hi[1];
// #endif

//                upls = amrex::max(amrex::min(upls, cc_umax), cc_umin);

//                D_TERM(xc = a_ccc(i-1,j,k,0);, // centroid of cell (i-1,j,k)
//                       yc = a_ccc(i-1,j,k,1);,
//                       zc = a_ccc(i-1,j,k,2););

//                D_TERM(delta_x = 0.5 - xc;,
//                       delta_y = yf  - yc;,
//                       delta_z = zf  - zc;);

//                // Compute slopes of component "0" of a_vcc
//                const auto& slopes_eb_lo =
// 		 amrex_calc_slopes_eb(i-1,j,k,0,a_vcc,a_ccc,a_flag);

// #if (AMREX_SPACEDIM==3)
//                Real umns = a_vcc(i-1,j,k,0) + delta_x * slopes_eb_lo[0]
//                                             + delta_y * slopes_eb_lo[1]
//                                             + delta_z * slopes_eb_lo[2];
// #else
//                Real umns = a_vcc(i-1,j,k,0) + delta_x * slopes_eb_lo[0]
//                                             + delta_y * slopes_eb_lo[1];
// #endif

//                umns = amrex::max(amrex::min(umns, cc_umax), cc_umin);

//                if ( umns < 0.0 && upls > 0.0 ) {
//                   a_u(i,j,k) = 0.0;
//                } else {
//                   Real avg = 0.5 * ( upls + umns );
//                   if ( std::abs(avg) <  small_vel) { a_u(i,j,k) = 0.0;
//                   } else if (avg >= 0)             { a_u(i,j,k) = umns;
//                   } else                           { a_u(i,j,k) = upls;
//                   }
//                }

//             } else {
//                a_u(i,j,k) = 0.0;
//             }
//         });
//     }

//     // ****************************************************************************
//     // Predict to y-faces
//     // ****************************************************************************
//     if ((extdir_or_ho_jlo and domain_jlo >= a_vbx.smallEnd(1)-1) or
//         (extdir_or_ho_jhi and domain_jhi <= a_vbx.bigEnd(1)))
//     {
//         amrex::ParallelFor(Box(a_vbx),
//         [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//         {
//             if (a_flag(i,j,k).isConnected(0,-1,0))
//             {
//                if (extdir_jlo and j == domain_jlo) {
//                    a_v(i,j,k) = a_vcc(i,j-1,k,1);
//                } else if (extdir_jhi and j == domain_jhi+1) {
//                    a_v(i,j,k) = a_vcc(i,j,k,1);
//                } else {

//                Real xf = a_fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
// #if (AMREX_SPACEDIM==3)
//                Real zf = a_fcy(i,j,k,1);
// #endif

//                D_TERM(Real xc = a_ccc(i,j,k,0);, // centroid of cell (i,j,k)
//                       Real yc = a_ccc(i,j,k,1);,
//                       Real zc = a_ccc(i,j,k,2););

//                D_TERM(Real delta_x = xf  - xc;,
//                       Real delta_y = 0.5 + yc;,
//                       Real delta_z = zf  - zc;);

//                Real cc_vmax = amrex::max(a_vcc(i,j,k,1), a_vcc(i,j-1,k,1));
//                Real cc_vmin = amrex::min(a_vcc(i,j,k,1), a_vcc(i,j-1,k,1));

//                // Compute slopes of component "1" of a_vcc
// 		 const auto& slopes_eb_hi =
// 		   amrex_calc_slopes_extdir_eb(i,j,k,1,a_vcc,a_ccc,
// 					       AMREX_D_DECL(a_fcx,a_fcy,a_fcz),a_flag,
// 					       AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
// 					       AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
// 					       AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
// 					       AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));


// #if (AMREX_SPACEDIM==3)
//                Real vpls = a_vcc(i,j  ,k,1) + delta_x * slopes_eb_hi[0]
//                                             - delta_y * slopes_eb_hi[1]
//                                             + delta_z * slopes_eb_hi[2];
// #else
//                Real vpls = a_vcc(i,j  ,k,1) + delta_x * slopes_eb_hi[0]
//                                             - delta_y * slopes_eb_hi[1];
// #endif

//                vpls = amrex::max(amrex::min(vpls, cc_vmax), cc_vmin);

//                D_TERM(xc = a_ccc(i,j-1,k,0);, // centroid of cell (i,j-1,k)
//                       yc = a_ccc(i,j-1,k,1);,
//                       zc = a_ccc(i,j-1,k,2););

//                D_TERM(delta_x = xf  - xc;,
//                       delta_y = 0.5 - yc;,
//                       delta_z = zf  - zc;);

//                // Compute slopes of component "1" of a_vcc
// 		 const auto& slopes_eb_lo =
// 		   amrex_calc_slopes_extdir_eb(i,j-1,k,1,a_vcc,a_ccc,
// 					       AMREX_D_DECL(a_fcx,a_fcy,a_fcz),a_flag,
// 					       AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
// 					       AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
// 					       AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
// 					       AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

// #if (AMREX_SPACEDIM==3)
//                Real vmns = a_vcc(i,j-1,k,1) + delta_x * slopes_eb_lo[0]
//                                             + delta_y * slopes_eb_lo[1]
//                                             + delta_z * slopes_eb_lo[2];
// #else
//                Real vmns = a_vcc(i,j-1,k,1) + delta_x * slopes_eb_lo[0]
//                                             + delta_y * slopes_eb_lo[1];
// #endif

//                vmns = amrex::max(amrex::min(vmns, cc_vmax), cc_vmin);

//                if ( vmns < 0.0 && vpls > 0.0 ) {
//                   a_v(i,j,k) = 0.0;
//                } else {
//                   Real avg = 0.5 * ( vpls + vmns );
//                   if ( std::abs(avg) <  small_vel) { a_v(i,j,k) = 0.0;
//                   } else if (avg >= 0)             { a_v(i,j,k) = vmns;
//                   } else                           { a_v(i,j,k) = vpls;
//                   }
//                }
//                }
//             } else {
//                a_v(i,j,k) = 0.0;
//             }
//         });
//     }
//     else
//     {
//         amrex::ParallelFor(Box(a_vbx),
//         [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//         {
//             if (a_flag(i,j,k).isConnected(0,-1,0))
//             {
//                Real xf = a_fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
// #if (AMREX_SPACEDIM==3)
//                Real zf = a_fcy(i,j,k,1);
// #endif

//                D_TERM(Real xc = a_ccc(i,j,k,0);, // centroid of cell (i,j,k)
//                       Real yc = a_ccc(i,j,k,1);,
//                       Real zc = a_ccc(i,j,k,2););

//                D_TERM(Real delta_x = xf  - xc;,
//                       Real delta_y = 0.5 + yc;,
//                       Real delta_z = zf  - zc;);

//                Real cc_vmax = amrex::max(a_vcc(i,j,k,1), a_vcc(i,j-1,k,1));
//                Real cc_vmin = amrex::min(a_vcc(i,j,k,1), a_vcc(i,j-1,k,1));

//                // Compute slopes of component "1" of a_vcc
//                const auto slopes_eb_hi =
// 		 amrex_calc_slopes_eb(i,j,k,1,a_vcc,a_ccc,a_flag);

// #if (AMREX_SPACEDIM==3)
//                Real vpls = a_vcc(i,j  ,k,1) + delta_x * slopes_eb_hi[0]
//                                             - delta_y * slopes_eb_hi[1]
//                                             + delta_z * slopes_eb_hi[2];
// #else
//                Real vpls = a_vcc(i,j  ,k,1) + delta_x * slopes_eb_hi[0]
//                                             - delta_y * slopes_eb_hi[1];
// #endif

//                vpls = amrex::max(amrex::min(vpls, cc_vmax), cc_vmin);

//                D_TERM(xc = a_ccc(i,j-1,k,0);,// centroid of cell (i,j-1,k)
//                       yc = a_ccc(i,j-1,k,1);,
//                       zc = a_ccc(i,j-1,k,2););

//                D_TERM(delta_x = xf  - xc;,
//                       delta_y = 0.5 - yc;,
//                       delta_z = zf  - zc;);

//                // Compute slopes of component "1" of a_vcc
//                const auto& slopes_eb_lo =
// 		 amrex_calc_slopes_eb(i,j-1,k,1,a_vcc,a_ccc,a_flag);

// #if (AMREX_SPACEDIM==3)
//                Real vmns = a_vcc(i,j-1,k,1) + delta_x * slopes_eb_lo[0]
//                                             + delta_y * slopes_eb_lo[1]
//                                             + delta_z * slopes_eb_lo[2];
// #else
//                Real vmns = a_vcc(i,j-1,k,1) + delta_x * slopes_eb_lo[0]
//                                             + delta_y * slopes_eb_lo[1];
// #endif

//                vmns = amrex::max(amrex::min(vmns, cc_vmax), cc_vmin);

//                if ( vmns < 0.0 && vpls > 0.0 ) {
//                   a_v(i,j,k) = 0.0;
//                } else {
//                   Real avg = 0.5 * ( vpls + vmns );
//                   if ( std::abs(avg) <  small_vel) { a_v(i,j,k) = 0.0;
//                   } else if (avg >= 0)             { a_v(i,j,k) = vmns;
//                   } else                           { a_v(i,j,k) = vpls;
//                   }
//                }

//             } else {
//                a_v(i,j,k) = 0.0;
//             }
//         });
//     }

// #if (AMREX_SPACEDIM==3)
//     // ****************************************************************************
//     // Predict to z-faces
//     // ****************************************************************************
//     if ((extdir_or_ho_klo and domain_klo >= a_wbx.smallEnd(2)-1) or
//         (extdir_or_ho_khi and domain_khi <= a_wbx.bigEnd(2)))
//     {
//         amrex::ParallelFor(Box(a_wbx),
//         [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//         {
//             if (a_flag(i,j,k).isConnected(0,0,-1))
//             {
// 	        if (extdir_or_ho_klo and k == domain_klo) {
//                     a_w(i,j,k) = a_vcc(i,j,k-1,2);
//                 } else if (extdir_or_ho_khi and k == domain_khi+1) {
//                     a_w(i,j,k) = a_vcc(i,j,k,2);
//                 } else {

//                Real xf = a_fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
//                Real yf = a_fcz(i,j,k,1);

//                Real xc = a_ccc(i,j,k,0); // centroid of cell (i,j,k)
//                Real yc = a_ccc(i,j,k,1);
//                Real zc = a_ccc(i,j,k,2);

//                Real delta_x = xf  - xc;
//                Real delta_y = yf  - yc;
//                Real delta_z = 0.5 + zc;

//                Real cc_wmax = amrex::max(a_vcc(i,j,k,2), a_vcc(i,j,k-1,2));
//                Real cc_wmin = amrex::min(a_vcc(i,j,k,2), a_vcc(i,j,k-1,2));

//                // Compute slopes of component "2" of a_vcc
// 		  const auto& slopes_eb_hi =
// 		    amrex_calc_slopes_extdir_eb(i,j,k,2,a_vcc,a_ccc,
// 						AMREX_D_DECL(a_fcx,a_fcy,a_fcz),a_flag,
// 						AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
// 						AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
// 						AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
// 						AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));


//                Real wpls = a_vcc(i,j,k  ,2) + delta_x * slopes_eb_hi[0]
//                                             + delta_y * slopes_eb_hi[1]
//                                             - delta_z * slopes_eb_hi[2];

//                wpls = amrex::max(amrex::min(wpls, cc_wmax), cc_wmin);

//                xc = a_ccc(i,j,k-1,0); // centroid of cell (i,j,k-1)
//                yc = a_ccc(i,j,k-1,1);
//                zc = a_ccc(i,j,k-1,2);

//                delta_x = xf  - xc;
//                delta_y = yf  - yc;
//                delta_z = 0.5 - zc;

//                // Compute slopes of component "2" of a_vcc
// 		  const auto& slopes_eb_lo =
// 		    amrex_calc_slopes_extdir_eb(i,j,k-1,2,a_vcc,a_ccc,
// 						AMREX_D_DECL(a_fcx,a_fcy,a_fcz),a_flag,
// 						AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
// 						AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
// 						AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
// 						AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi));

//                Real wmns = a_vcc(i,j,k-1,2) + delta_x * slopes_eb_lo[0]
//                                             + delta_y * slopes_eb_lo[1]
//                                             + delta_z * slopes_eb_lo[2];

//                wmns = amrex::max(amrex::min(wmns, cc_wmax), cc_wmin);

//                if ( wmns < 0.0 && wpls > 0.0 ) {
//                   a_w(i,j,k) = 0.0;
//                } else {
//                   Real avg = 0.5 * ( wpls + wmns );
//                   if ( std::abs(avg) <  small_vel) { a_w(i,j,k) = 0.0;
//                   } else if (avg >= 0)             { a_w(i,j,k) = wmns;
//                   } else                           { a_w(i,j,k) = wpls;
//                   }
//                }
//                 }
//             } else {
//                a_w(i,j,k) = 0.0;
//             }
//         });
//     }
//     else
//     {
//         amrex::ParallelFor(Box(a_wbx),
//         [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//         {
//             if (a_flag(i,j,k).isConnected(0,0,-1))
//             {
//                Real xf = a_fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
//                Real yf = a_fcz(i,j,k,1);

//                Real xc = a_ccc(i,j,k,0); // centroid of cell (i,j,k)
//                Real yc = a_ccc(i,j,k,1);
//                Real zc = a_ccc(i,j,k,2);

//                Real delta_x = xf  - xc;
//                Real delta_y = yf  - yc;
//                Real delta_z = 0.5 + zc;

//                Real cc_wmax = amrex::max(a_vcc(i,j,k,2), a_vcc(i,j,k-1,2));
//                Real cc_wmin = amrex::min(a_vcc(i,j,k,2), a_vcc(i,j,k-1,2));

//                // Compute slopes of component "2" of a_vcc
//                const auto slopes_eb_hi = amrex_calc_slopes_eb(i,j,k,2,a_vcc,a_ccc,a_flag);

//                Real wpls = a_vcc(i,j,k  ,2) + delta_x * slopes_eb_hi[0]
//                                             + delta_y * slopes_eb_hi[1]
//                                             - delta_z * slopes_eb_hi[2];

//                wpls = amrex::max(amrex::min(wpls, cc_wmax), cc_wmin);

//                xc = a_ccc(i,j,k-1,0); // centroid of cell (i,j,k-1)
//                yc = a_ccc(i,j,k-1,1);
//                zc = a_ccc(i,j,k-1,2);

//                delta_x = xf  - xc;
//                delta_y = yf  - yc;
//                delta_z = 0.5 - zc;

//                // Compute slopes of component "2" of a_vcc
//                const auto& slopes_eb_lo = amrex_calc_slopes_eb(i,j,k-1,2,a_vcc,a_ccc,a_flag);

//                Real wmns = a_vcc(i,j,k-1,2) + delta_x * slopes_eb_lo[0]
//                                             + delta_y * slopes_eb_lo[1]
//                                             + delta_z * slopes_eb_lo[2];

//                wmns = amrex::max(amrex::min(wmns, cc_wmax), cc_wmin);

//                if ( wmns < 0.0 && wpls > 0.0 ) {
//                   a_w(i,j,k) = 0.0;
//                } else {
//                   Real avg = 0.5 * ( wpls + wmns );
//                   if ( std::abs(avg) <  small_vel) { a_w(i,j,k) = 0.0;
//                   } else if (avg >= 0)             { a_w(i,j,k) = wmns;
//                   } else                           { a_w(i,j,k) = wpls;
//                   }
//                }

//             } else {
//                a_w(i,j,k) = 0.0;
//             }
//         });
//     }
// #endif
// }


// #endif
