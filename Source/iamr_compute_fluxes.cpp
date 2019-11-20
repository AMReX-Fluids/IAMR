#include <Godunov.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_BCRec.H>

using namespace amrex;


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
                    D_DECL( FArrayBox& a_fx,
                            FArrayBox& a_fy,
                            FArrayBox& a_fz),
                    const FArrayBox& a_state,
                    const int a_comp,
                    const int a_ncomp,
                    D_DECL( const FArrayBox& a_xsl,
                            const FArrayBox& a_ysl,
                            const FArrayBox& a_zsl),
                    const int a_sl_comp,
                    D_DECL( const FArrayBox& a_umac,
                            const FArrayBox& a_vmac,
                            const FArrayBox& a_wmac),
                    const Box&       a_domain,
                    const Vector<BCRec>& a_bcs )
{

    const Dim3 domlo = amrex::lbound(a_domain);
    const Dim3 domhi = amrex::ubound(a_domain);

    const auto& state = a_state.array();

    D_TERM( const auto& fx  = a_fx.array();,
            const auto& fy  = a_fy.array();,
            const auto& fz  = a_fz.array(););

    D_TERM( const auto& u   = a_umac.array();,
            const auto& v   = a_vmac.array();,
            const auto& w   = a_wmac.array(););

    D_TERM( const auto& xsl = a_xsl.array();,
            const auto& ysl = a_ysl.array();,
            const auto& zsl = a_zsl.array(););

    D_TERM( const Box ubx = amrex::surroundingNodes(a_bx,0);,
            const Box vbx = amrex::surroundingNodes(a_bx,1);,
            const Box wbx = amrex::surroundingNodes(a_bx,2););

    const auto bc = a_bcs.dataPtr();

    AMREX_FOR_4D(ubx, a_ncomp, i, j, k, n,
    {
        Real state_w(0)  ;
        Real state_mns(0);
        Real state_pls(0);

        //
        // West face
        //
        // In the case of inflow we are using the prescribed Dirichlet value
        // This is from incflo: In the case of PINF, POUT we are using the upwind value
        if ( (i == domlo.x) and (bc[n].lo(0) == BCType::ext_dir) )
        {
            state_w = state(i-1,j,k,a_comp+n);
        }
        else if ( (i == domhi.x+1) and (bc[n].hi(0) == BCType::ext_dir) )
        {
            state_w = state(i,j,k,a_comp+n);
        }
        else
        {
            state_pls = state(i  ,j,k,a_comp+n) - .5*xsl(i  ,j,k,a_sl_comp+n);
            state_mns = state(i-1,j,k,a_comp+n) + .5*xsl(i-1,j,k,a_sl_comp+n);
            state_w   = upwind( state_mns, state_pls, u(i,j,k) );
        }

        fx(i,j,k,n) = u(i,j,k) * state_w;
    });

    AMREX_FOR_4D(vbx, a_ncomp, i, j, k, n,
    {
        Real state_s(0);
        Real state_mns(0);
        Real state_pls(0);

        //
        // South face
        //
        // In the case of inflow we are using the prescribed Dirichlet value
        // This is from incflo: In the case of PINF, POUT we are using the upwind value
        if ( (j == domlo.y) and (bc[n].lo(1) == BCType::ext_dir) )
        {
            state_s = state(i,j-1,k,a_comp+n);
        }
        else if ( (j == domhi.y+1) and (bc[n].hi(1) == BCType::ext_dir) )
        {
            state_s = state(i,j,k,a_comp+n);
        }
        else
        {
            state_pls = state(i,j  ,k,a_comp+n) - .5*ysl(i,j  ,k,a_sl_comp+n);
            state_mns = state(i,j-1,k,a_comp+n) + .5*ysl(i,j-1,k,a_sl_comp+n);
            state_s   = upwind( state_mns, state_pls, v(i,j,k) );
        }

        fy(i,j,k,n) = v(i,j,k) * state_s;
    });

#if ( AMREX_SPACEDIM ==3 )
    AMREX_FOR_4D(wbx, a_ncomp, i, j, k, n,
    {
        Real state_b(0);
        Real state_mns(0);
        Real state_pls(0);

        //
        // Bottom face
        //
        // In the case of inflow we are using the prescribed Dirichlet value
        // This is from incflo: In the case of PINF, POUT we are using the upwind value
        if ( (k == domlo.z) and (bc[n].lo(2) == BCType::ext_dir) )
        {
            state_b = state(i,j,k-1,a_comp+n);
        }
        else if ( (k == domhi.z+1) and (bc[n].lo(2) == BCType::ext_dir) )
        {
            state_b = state(i,j,k,a_comp+n);
        }
        else
        {
            state_pls = state(i,j,k  ,a_comp+n) - .5*zsl(i,j,k  ,a_sl_comp+n);
            state_mns = state(i,j,k-1,a_comp+n) + .5*zsl(i,j,k-1,a_sl_comp+n);
            state_b   = upwind( state_mns, state_pls, w(i,j,k) );
        }

        fz(i,j,k,n) = w(i,j,k) * state_b;
    });
#endif
}

//
// Compute fluxes on given EB box
//
void
ComputeFluxesOnEBBox (const Box& a_bx,
                      D_DECL( FArrayBox& a_fx,
                              FArrayBox& a_fy,
                              FArrayBox& a_fz),
                      const FArrayBox& a_state,
                      const int a_comp,
                      const int a_ncomp,
                      D_DECL( const FArrayBox& a_xsl,
                              const FArrayBox& a_ysl,
                              const FArrayBox& a_zsl),
                      const int a_sl_comp,
                      D_DECL( const FArrayBox& a_umac,
                              const FArrayBox& a_vmac,
                              const FArrayBox& a_wmac),
                      const Box&       a_domain,
                      const Vector<BCRec>& a_bcs,
                      D_DECL( const FArrayBox& a_afracx,
                              const FArrayBox& a_afracy,
                              const FArrayBox& a_afracz),
                      D_DECL( const FArrayBox& a_face_centx,
                              const FArrayBox& a_face_centy,
                              const FArrayBox& a_face_centz),
                      const IArrayBox& a_cc_mask,
                      const EBCellFlagFab& a_flags)
{
    const Dim3 domlo = amrex::lbound(a_domain);
    const Dim3 domhi = amrex::ubound(a_domain);

    const auto& state = a_state.array();

    D_TERM( const auto& fx  = a_fx.array();,
            const auto& fy  = a_fy.array();,
            const auto& fz  = a_fz.array(););

    D_TERM( const auto& u   = a_umac.array();,
            const auto& v   = a_vmac.array();,
            const auto& w   = a_wmac.array(););

    D_TERM( const auto& xsl = a_xsl.array();,
            const auto& ysl = a_ysl.array();,
            const auto& zsl = a_zsl.array(););

    D_TERM( const Box ubx = amrex::surroundingNodes(a_bx,0);,
            const Box vbx = amrex::surroundingNodes(a_bx,1);,
            const Box wbx = amrex::surroundingNodes(a_bx,2););

    D_TERM( const Box ubx_grown = amrex::surroundingNodes(amrex::grow(a_bx,1),0);,
            const Box vbx_grown = amrex::surroundingNodes(amrex::grow(a_bx,1),1);,
            const Box wbx_grown = amrex::surroundingNodes(amrex::grow(a_bx,1),2););

    D_TERM( const auto& areafrac_x = a_afracx.array();,
            const auto& areafrac_y = a_afracy.array();,
            const auto& areafrac_z = a_afracz.array(););

    D_TERM( FArrayBox s_on_x_face(ubx_grown, a_ncomp);,
            FArrayBox s_on_y_face(vbx_grown, a_ncomp);,
            FArrayBox s_on_z_face(wbx_grown, a_ncomp););

    D_TERM( const auto& fcx_fab = a_face_centx.array();,
            const auto& fcy_fab = a_face_centy.array();,
            const auto& fcz_fab = a_face_centz.array(););

    // These lines ensure that the temporary Fabs above aren't destroyed
    //   before we're done with them when running with GPUs
    D_TERM( Elixir eli_x = s_on_x_face.elixir();,
            Elixir eli_y = s_on_y_face.elixir();,
            Elixir eli_z = s_on_z_face.elixir(););

    D_TERM( const auto& sx = s_on_x_face.array();,
            const auto& sy = s_on_y_face.array();,
            const auto& sz = s_on_z_face.array(););

    const auto& ccm_fab = a_cc_mask.const_array();

    const auto bc = a_bcs.dataPtr();

    //
    // ===================== X =====================
    //
    AMREX_FOR_4D(ubx_grown, a_ncomp, i, j, k, n,
    {
        Real upls(0.0);
        Real umns(0.0);

        if( areafrac_x(i,j,k) > 0 )
        {
            if ( (i == domlo.x) and (bc[n].lo(0) == BCType::ext_dir) )
            {
                sx(i,j,k,n) = state(domlo.x-1,j,k,a_comp+n);
            }
            else if ( (i == domhi.x+1) and (bc[n].hi(0) == BCType::ext_dir) )
            {
                sx(i,j,k,n) = state(domhi.x+1,j,k,a_comp+n);
            }
            else
            {
                upls = state(i  ,j,k,a_comp+n) - .5*xsl(i  ,j,k,a_sl_comp+n);
                umns = state(i-1,j,k,a_comp+n) + .5*xsl(i-1,j,k,a_sl_comp+n);

                sx(i,j,k,n) = upwind( umns, upls, u(i,j,k) );
            }
        }
        else
        {
            sx(i,j,k,n) = COVERED_VAL;
        }
    });

    // Interpolate to face centroid
    AMREX_FOR_4D(ubx, a_ncomp, i, j, k, n,
    {
        if( areafrac_x(i,j,k) > 0 )
        {
            int jj = j + static_cast<int>(std::copysign(1.0, fcx_fab(i,j,k,0)));
#if (AMREX_SPACEDIM==3)
            int kk = k + static_cast<int>(std::copysign(1.0, fcx_fab(i,j,k,1)));
#else
            int kk = k;
#endif

            Real fracy = (ccm_fab(i-1,jj,k) || ccm_fab(i,jj,k)) ? std::abs(fcx_fab(i,j,k,0)) : 0.0;
#if (AMREX_SPACEDIM==3)
            Real fracz = (ccm_fab(i-1,j,kk) || ccm_fab(i,j,kk)) ? std::abs(fcx_fab(i,j,k,1)) : 0.0;
#else
            Real fracz(0.0);
#endif

            Real s_on_x_centroid = (1.0-fracy)*(1.0-fracz)*sx(i, j,k ,n)+
                fracy *(1.0-fracz)*sx(i,jj,k ,n)+
                fracz *(1.0-fracy)*sx(i, j,kk,n)+
                fracy *     fracz *sx(i,jj,kk,n);

            fx(i,j,k,n) = u(i,j,k) * s_on_x_centroid;
        }
        else
        {
            fx(i,j,k,n) = COVERED_VAL;
        }
    });


    //
    // ===================== Y =====================
    //
    AMREX_FOR_4D(vbx_grown, a_ncomp, i, j, k, n,
    {
        Real vpls(0.0);
        Real vmns(0.0);

        if( areafrac_y(i,j,k) > 0 )
        {
            if ( (j == domlo.y) and (bc[n].lo(1) == BCType::ext_dir) )
            {
                sy(i,j,k,n) = state(i,domlo.y-1,k,a_comp+n);
            }
            else if ( (j == domhi.y+1) and (bc[n].hi(1) == BCType::ext_dir) )
            {
                sy(i,j,k,n) = state(i,domhi.y+1,k,a_comp+n);
            }
            else
            {
                vpls = state(i,j  ,k,a_comp+n) - .5*ysl(i,j  ,k,a_sl_comp+n);
                vmns = state(i,j-1,k,a_comp+n) + .5*ysl(i,j-1,k,a_sl_comp+n);

                sy(i,j,k,n) = upwind( vmns, vpls, v(i,j,k) );
            }
        }
        else
        {
            sy(i,j,k,n) = COVERED_VAL;
        }
    });

    AMREX_FOR_4D(vbx, a_ncomp, i, j, k, n,
    {
        if ( areafrac_y(i,j,k) > 0 )
        {
            int ii = i + static_cast<int>(std::copysign(1.0,fcy_fab(i,j,k,0)));
#if (AMREX_SPACEDIM==3)
            int kk = k + static_cast<int>(std::copysign(1.0,fcy_fab(i,j,k,1)));
#else
            int kk = k;
#endif

            Real fracx = (ccm_fab(ii,j-1,k) || ccm_fab(ii,j,k)) ? std::abs(fcy_fab(i,j,k,0)) : 0.0;
#if (AMREX_SPACEDIM==3)
            Real fracz = (ccm_fab(i,j-1,kk) || ccm_fab(i,j,kk)) ? std::abs(fcy_fab(i,j,k,1)) : 0.0;
#else
            Real fracz(0.0);
#endif

            Real s_on_y_centroid = (1.0-fracx)*(1.0-fracz)*sy(i ,j,k ,n)+
                fracx *(1.0-fracz)*sy(ii,j,k ,n)+
                fracz *(1.0-fracx)*sy(i ,j,kk,n)+
                fracx *     fracz *sy(ii,j,kk,n);
            fy(i,j,k,n) = v(i,j,k) * s_on_y_centroid;
        }
        else
        {
            fy(i,j,k,n) = COVERED_VAL;
        }
    });


    //
    // ===================== Z =====================
    //
#if ( AMREX_SPACEDIM == 3 )
    AMREX_FOR_4D(wbx_grown, a_ncomp, i, j, k, n,
    {
        Real wpls(0.0);
        Real wmns(0.0);

        if( areafrac_z(i,j,k) > 0 )
        {
            if ( (k == domlo.z) and (bc[n].lo(2) == BCType::ext_dir) )
            {
                sz(i,j,k,n) = state(i,j,domlo.z-1,a_comp+n);
            }
            else if ( (k == domhi.z+1) and (bc[n].hi(2) == BCType::ext_dir) )
            {
                sz(i,j,k,n) = state(i,j,domhi.z+1,a_comp+n);
            }
            else
            {
                wpls = state(i,j,k  ,a_comp+n) - .5*zsl(i,j,k  ,a_sl_comp+n);
                wmns = state(i,j,k-1,a_comp+n) + .5*zsl(i,j,k-1,a_sl_comp+n);

                sz(i,j,k,n) = upwind( wmns, wpls, w(i,j,k) );
            }
        }
        else
        {
            sz(i,j,k,n) = COVERED_VAL;
        }
    });

    AMREX_FOR_4D(wbx, a_ncomp, i, j, k, n,
    {
        if( areafrac_z(i,j,k) > 0 )
        {
            int ii = i + static_cast<int>(std::copysign(1.0,fcz_fab(i,j,k,0)));
            int jj = j + static_cast<int>(std::copysign(1.0,fcz_fab(i,j,k,1)));

            Real fracx = (ccm_fab(ii,j,k-1) || ccm_fab(ii,j,k)) ? std::abs(fcz_fab(i,j,k,0)) : 0.0;
            Real fracy = (ccm_fab(i,jj,k-1) || ccm_fab(i,jj,k)) ? std::abs(fcz_fab(i,j,k,1)) : 0.0;

            Real s_on_z_centroid = (1.0-fracx)*(1.0-fracy)*sz(i ,j ,k,n)+
                fracx *(1.0-fracy)*sz(ii,j ,k,n)+
                fracy *(1.0-fracx)*sz(i ,jj,k,n)+
                fracx *     fracy *sz(ii,jj,k,n);

            fz(i,j,k,n) = w(i,j,k) * s_on_z_centroid;
        }
        else
        {
            fz(i,j,k,n) = COVERED_VAL;
        }
    });
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
                         const Vector<BCRec>& a_bcs )
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

    // Create cc_mask
    iMultiFab cc_mask(a_state.boxArray(), a_state.DistributionMap(), 1, 1);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        std::vector< std::pair<int,Box> > isects;
        const std::vector<IntVect>& pshifts = a_geom.periodicity().shiftIntVect();
        const BoxArray& ba = cc_mask.boxArray();
        for (MFIter mfi(cc_mask); mfi.isValid(); ++mfi)
        {
            Array4<int> const& fab = cc_mask.array(mfi);

            const Box& bx = mfi.fabbox();
            for (const auto& iv : pshifts)
            {
                ba.intersections(bx+iv, isects);
                for (const auto& is : isects)
                {
                    const Box& b = is.second-iv;
                    AMREX_FOR_3D ( b, i, j, k,
                    {
                        fab(i,j,k) = 1;
                    });
                }
            }
            // NOTE: here we do not need host-device synchronization since it
            // is already included in the MFIter destructor
        }
    }

    // Initialize fluxes
    D_TERM(a_fx.setVal(COVERED_VAL);,
           a_fy.setVal(COVERED_VAL);,
           a_fz.setVal(COVERED_VAL););

    for (MFIter mfi(a_state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox ();

        const EBFArrayBox&  state_fab = static_cast<EBFArrayBox const&>(a_state[mfi]);
        const EBCellFlagFab&    flags = state_fab.getEBCellFlagFab();

        if (flags.getType(amrex::grow(bx,0)) != FabType::covered )
        {
            // No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
            if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular )
            {
                fluxes::ComputeFluxesOnBox( bx, D_DECL(a_fx[mfi], a_fy[mfi], a_fz[mfi]), a_state[mfi], a_comp, a_ncomp,
                                            D_DECL(a_xsl[mfi], a_ysl[mfi], a_zsl[mfi]), a_sl_comp,
                                            D_DECL(a_umac[mfi], a_vmac[mfi], a_wmac[mfi]), domain, a_bcs);

            }
            else
            {
                fluxes::ComputeFluxesOnEBBox(bx, D_DECL(a_fx[mfi], a_fy[mfi], a_fz[mfi]), a_state[mfi], a_comp, a_ncomp,
                                             D_DECL(a_xsl[mfi], a_ysl[mfi], a_zsl[mfi]), a_sl_comp,
                                             D_DECL(a_umac[mfi], a_vmac[mfi], a_wmac[mfi]), domain, a_bcs,
                                             D_DECL((*areafrac[0])[mfi], (*areafrac[1])[mfi], (*areafrac[2])[mfi]),
                                             D_DECL((*facecent[0])[mfi], (*facecent[1])[mfi], (*facecent[2])[mfi]),
                                             cc_mask[mfi], flags);
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
