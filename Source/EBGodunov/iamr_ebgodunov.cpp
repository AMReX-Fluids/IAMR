#include <iamr_ebgodunov.H>
#include <iamr_godunov.H>
#include <iamr_redistribution.H>
// #include <NS_util.H>

using namespace amrex;


void
EBGodunov::ComputeAofs ( MultiFab& aofs, const int aofs_comp, const int ncomp,
                         MultiFab const& state, const int state_comp,
                         AMREX_D_DECL( MultiFab const& umac,
                                       MultiFab const& vmac,
                                       MultiFab const& wmac),
                         AMREX_D_DECL( MultiFab& xedge,
                                       MultiFab& yedge,
                                       MultiFab& zedge),
                         const int  edge_comp,
                         const bool known_edgestate,
                         AMREX_D_DECL( MultiFab& xfluxes,
                                       MultiFab& yfluxes,
                                       MultiFab& zfluxes),
                         int fluxes_comp,
                         MultiFab const& fq,
                         const int fq_comp,
                         MultiFab const& divu,
                         Vector<BCRec> const& h_bc,
                         BCRec const* d_bc,
                         Geometry const& geom,
                         Gpu::DeviceVector<int>& iconserv,
                         const Real dt,
                         const bool is_velocity,
                         std::string redistribution_type)
{
    BL_PROFILE("EBGodunov::ComputeAofs()");

    AMREX_ALWAYS_ASSERT(state.hasEBFabFactory());

    for (int n = 0; n < ncomp; n++)
       if (!iconserv[n]) amrex::Abort("EBGodunov does not support non-conservative form");

    auto const& ebfact= dynamic_cast<EBFArrayBoxFactory const&>(state.Factory());
    auto const& flags = ebfact.getMultiEBCellFlagFab();
    auto const& fcent = ebfact.getFaceCent();
    auto const& ccent = ebfact.getCentroid();
    auto const& vfrac = ebfact.getVolFrac();
    auto const& areafrac = ebfact.getAreaFrac();

    //FIXME - check on adding tiling here
    for (MFIter mfi(aofs); mfi.isValid(); ++mfi)
    {

        const Box& bx   = mfi.tilebox();

        auto const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
        bool regular = (flagfab.getType(amrex::grow(bx,2)) == FabType::regular);

        // Get handlers to Array4
        //
        AMREX_D_TERM( const auto& fx = xfluxes.array(mfi,fluxes_comp);,
                      const auto& fy = yfluxes.array(mfi,fluxes_comp);,
                      const auto& fz = zfluxes.array(mfi,fluxes_comp););

        AMREX_D_TERM( const auto& xed = xedge.array(mfi,edge_comp);,
                      const auto& yed = yedge.array(mfi,edge_comp);,
                      const auto& zed = zedge.array(mfi,edge_comp););

        AMREX_D_TERM( const auto& u = umac.const_array(mfi);,
                      const auto& v = vmac.const_array(mfi);,
                      const auto& w = wmac.const_array(mfi););

        if (regular)   // Plain Godunov
        {
            if (not known_edgestate)
            {
                Godunov::ComputeEdgeState( bx, ncomp,
                                           state.array(mfi,state_comp),
                                           AMREX_D_DECL( xed, yed, zed ),
                                           AMREX_D_DECL( u, v, w ),
                                           divu.array(mfi),
                                           fq.array(mfi,fq_comp),
                                           geom, dt, d_bc,
                                           iconserv.data(),
                                           false,
                                           false,
                                           is_velocity );
            }

            Godunov::ComputeFluxes( bx, AMREX_D_DECL( fx, fy, fz ),
                                    AMREX_D_DECL( u, v, w ),
                                    AMREX_D_DECL( xed, yed, zed ),
                                    geom, ncomp );

            Godunov::ComputeDivergence( bx,
                                        aofs.array(mfi,aofs_comp),
                                        AMREX_D_DECL( fx, fy, fz ),
                                        AMREX_D_DECL( xed, yed, zed ),
                                        AMREX_D_DECL( u, v, w ),
                                        ncomp, geom, iconserv.data() );
        }
        else     // EB Godunov
        {


            Box gbx = bx;
            // We need 3 if we are doing state redistribution
            if (redistribution_type == "StateRedist" ||
                redistribution_type == "MergeRedist")
                gbx.grow(3);
            else if (redistribution_type == "FluxRedist")
                gbx.grow(2);
            else if (redistribution_type == "NoRedist")
                gbx.grow(1);
            else
                amrex::Abort("Dont know this redistribution type");

            AMREX_D_TERM(Array4<Real const> const& fcx = fcent[0]->const_array(mfi);,
                         Array4<Real const> const& fcy = fcent[1]->const_array(mfi);,
                         Array4<Real const> const& fcz = fcent[2]->const_array(mfi););

            AMREX_D_TERM(Array4<Real const> const& apx = areafrac[0]->const_array(mfi);,
                         Array4<Real const> const& apy = areafrac[1]->const_array(mfi);,
                         Array4<Real const> const& apz = areafrac[2]->const_array(mfi););

            Array4<Real const> const& ccent_arr = ccent.const_array(mfi);
            Array4<Real const> const& vfrac_arr = vfrac.const_array(mfi);
            auto const& flags_arr  = flags.const_array(mfi);

            int ngrow = 4;

            if (redistribution_type=="StateRedist")
                ++ngrow;

            FArrayBox tmpfab(amrex::grow(bx,ngrow),  (4*AMREX_SPACEDIM + 2)*ncomp);
            Elixir    eli = tmpfab.elixir();


            if (not known_edgestate)
            {
                EBGodunov::ComputeEdgeState( gbx, ncomp,
                                             state.array(mfi,state_comp),
                                             AMREX_D_DECL( xed, yed, zed ),
                                             AMREX_D_DECL( u, v, w ),
                                             divu.array(mfi),
                                             fq.array(mfi,fq_comp),
                                             geom, dt, h_bc, d_bc,
                                             iconserv.data(),
                                             tmpfab.dataPtr(),
                                             flags_arr,
                                             AMREX_D_DECL( apx, apy, apz ),
                                             vfrac_arr,
                                             AMREX_D_DECL( fcx, fcy, fcz ),
                                             ccent_arr,
                                             is_velocity );
            }

            EBGodunov::ComputeFluxes( gbx, AMREX_D_DECL( fx, fy, fz ),
                                      AMREX_D_DECL( u, v, w ),
                                      AMREX_D_DECL( xed, yed, zed ),
                                      AMREX_D_DECL( apx, apy, apz ),
                                      geom, ncomp, flags_arr );

            // div at ncomp*3 to make space for the 3 redistribute temporaries
            Array4<Real> divtmp_arr = tmpfab.array(ncomp*3);

            EBGodunov::ComputeDivergence( gbx,
                                          divtmp_arr,
                                          AMREX_D_DECL( fx, fy, fz ),
                                          vfrac_arr, ncomp, geom );

            Array4<Real> scratch = tmpfab.array(0);
            Redistribution::Apply( bx, ncomp, aofs.array(mfi, aofs_comp), divtmp_arr,
                                   state.const_array(mfi, state_comp), scratch, flags_arr,
                                   AMREX_D_DECL(apx,apy,apz), vfrac_arr,
                                   AMREX_D_DECL(fcx,fcy,fcz), ccent_arr, geom, dt,
                                   redistribution_type );

            // Change sign because for EB we computed -div
	    auto const& aofs_arr = aofs.array(mfi, aofs_comp); 
            amrex::ParallelFor(bx, ncomp, [aofs_arr]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            { aofs_arr( i, j, k, n ) *=  - 1.0; });
	 }

        // Note this sync is needed since ComputeEdgeState() contains temporaries
	// Not sure it's really needed when known_edgestate==true
        Gpu::streamSynchronize();  // otherwise we might be using too much memory
    }


}



void
EBGodunov::ComputeSyncAofs ( MultiFab& aofs, const int aofs_comp, const int ncomp,
                             MultiFab const& state, const int state_comp,
                             AMREX_D_DECL( MultiFab const& umac,
                                           MultiFab const& vmac,
                                           MultiFab const& wmac),
                             AMREX_D_DECL( MultiFab const& ucorr,
                                           MultiFab const& vcorr,
                                           MultiFab const& wcorr),
                             AMREX_D_DECL( MultiFab& xedge,
                                           MultiFab& yedge,
                                           MultiFab& zedge),
                             const int  edge_comp,
                             const bool known_edgestate,
                             AMREX_D_DECL( MultiFab& xfluxes,
                                           MultiFab& yfluxes,
                                           MultiFab& zfluxes),
                             int fluxes_comp,
                             MultiFab const& fq,
                             const int fq_comp,
                             MultiFab const& divu,
                             Vector<BCRec> const& h_bc,
                             BCRec const* d_bc,
                             Geometry const& geom,
                             Gpu::DeviceVector<int>& iconserv,
                             const Real dt,
                             const bool is_velocity,
                             std::string redistribution_type )
{
    BL_PROFILE("EBGodunov::ComputeAofs()");

    AMREX_ALWAYS_ASSERT(state.hasEBFabFactory());

    for (int n = 0; n < ncomp; n++)
       if (!iconserv[n]) amrex::Abort("EBGodunov does not support non-conservative form");

    auto const& ebfact= dynamic_cast<EBFArrayBoxFactory const&>(state.Factory());
    auto const& flags = ebfact.getMultiEBCellFlagFab();
    auto const& fcent = ebfact.getFaceCent();
    auto const& ccent = ebfact.getCentroid();
    auto const& vfrac = ebfact.getVolFrac();
    auto const& areafrac = ebfact.getAreaFrac();


    //FIXME - check on adding tiling here
    for (MFIter mfi(aofs); mfi.isValid(); ++mfi)
    {

        const Box& bx   = mfi.tilebox();

        auto const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
        bool regular = (flagfab.getType(amrex::grow(bx,2)) == FabType::regular);

        //
        // Get handlers to Array4
        //
        AMREX_D_TERM( const auto& fx = xfluxes.array(mfi,fluxes_comp);,
                      const auto& fy = yfluxes.array(mfi,fluxes_comp);,
                      const auto& fz = zfluxes.array(mfi,fluxes_comp););

        AMREX_D_TERM( const auto& xed = xedge.array(mfi,edge_comp);,
                      const auto& yed = yedge.array(mfi,edge_comp);,
                      const auto& zed = zedge.array(mfi,edge_comp););

        AMREX_D_TERM( const auto& uc = ucorr.const_array(mfi);,
                      const auto& vc = vcorr.const_array(mfi);,
                      const auto& wc = wcorr.const_array(mfi););

        if (regular) // Plain Godunov
        {
            if (not known_edgestate)
            {
                AMREX_D_TERM( const auto& u = umac.const_array(mfi);,
                              const auto& v = vmac.const_array(mfi);,
                              const auto& w = wmac.const_array(mfi););

                Godunov::ComputeEdgeState( bx, ncomp,
                                           state.array(mfi,state_comp),
                                           AMREX_D_DECL( xed, yed, zed ),
                                           AMREX_D_DECL( u, v, w ),
                                           divu.array(mfi),
                                           fq.array(mfi,fq_comp),
                                           geom, dt, d_bc,
                                           iconserv.data(),
                                           false,
                                           false,
                                           is_velocity );
            }

            Godunov::ComputeFluxes( bx, AMREX_D_DECL( fx, fy, fz ),
                                    AMREX_D_DECL( uc, vc, wc ),
                                    AMREX_D_DECL( xed, yed, zed ),
                                    geom, ncomp );


            Godunov::ComputeSyncDivergence( bx,
                                            aofs.array(mfi,aofs_comp),
                                            AMREX_D_DECL( fx, fy, fz ),
                                            ncomp, geom );
        }
        else  // EB Godunov
        {
            Box gbx = bx;
            gbx.grow(2);

            AMREX_D_TERM(Array4<Real const> const& fcx = fcent[0]->const_array(mfi);,
                         Array4<Real const> const& fcy = fcent[1]->const_array(mfi);,
                         Array4<Real const> const& fcz = fcent[2]->const_array(mfi););

            AMREX_D_TERM(Array4<Real const> const& apx = areafrac[0]->const_array(mfi);,
                         Array4<Real const> const& apy = areafrac[1]->const_array(mfi);,
                         Array4<Real const> const& apz = areafrac[2]->const_array(mfi););

            Array4<Real const> const& ccent_arr = ccent.const_array(mfi);
            Array4<Real const> const& vfrac_arr = vfrac.const_array(mfi);
            auto const& flags_arr  = flags.const_array(mfi);

            int ngrow = 4;
            FArrayBox tmpfab(amrex::grow(bx,ngrow),  (4*AMREX_SPACEDIM + 2)*ncomp);
            Elixir    eli = tmpfab.elixir();


            if (not known_edgestate)
            {
                AMREX_D_TERM( const auto& u = umac.const_array(mfi);,
                              const auto& v = vmac.const_array(mfi);,
                              const auto& w = wmac.const_array(mfi););

                EBGodunov::ComputeEdgeState( gbx, ncomp,
                                             state.array(mfi,state_comp),
                                             AMREX_D_DECL( xed, yed, zed ),
                                             AMREX_D_DECL( u, v, w ),
                                             divu.array(mfi),
                                             fq.array(mfi,fq_comp),
                                             geom, dt, h_bc, d_bc,
                                             iconserv.data(),
                                             tmpfab.dataPtr(),
                                             flags_arr,
                                             AMREX_D_DECL( apx, apy, apz ),
                                             vfrac_arr,
                                             AMREX_D_DECL( fcx, fcy, fcz ),
                                             ccent_arr,
                                             is_velocity );
            }

            EBGodunov::ComputeFluxes( gbx, AMREX_D_DECL( fx, fy, fz ),
                                      AMREX_D_DECL( uc, vc, wc ),
                                      AMREX_D_DECL( xed, yed, zed ),
                                      AMREX_D_DECL( apx, apy, apz ),
                                      geom, ncomp, flags_arr );

            // div at ncomp*3 to make space for the 3 redistribute temporaries
            Array4<Real> divtmp_arr = tmpfab.array(ncomp*3);
            Array4<Real> divtmp_redist_arr = tmpfab.array(ncomp*4);

            EBGodunov::ComputeDivergence( gbx,
                                          divtmp_arr,
                                          AMREX_D_DECL( fx, fy, fz ),
                                          vfrac_arr, ncomp, geom );

            Array4<Real> scratch = tmpfab.array(0);

            Redistribution::Apply( bx, ncomp, divtmp_redist_arr, divtmp_arr,
                                   state.const_array(mfi, state_comp), scratch, flags_arr,
                                   AMREX_D_DECL(apx,apy,apz), vfrac_arr,
                                   AMREX_D_DECL(fcx,fcy,fcz), ccent_arr, geom, dt,
                                   redistribution_type );

            // Subtract contribution to sync aofs -- sign of divergence is aofs is opposite
            // of sign to div as computed by EBGOdunov::ComputeDivergence, thus it must be subtracted.
            auto const& aofs_arr = aofs.array(mfi, aofs_comp);
            amrex::ParallelFor(bx, ncomp, [aofs_arr, divtmp_redist_arr]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                aofs_arr( i, j, k, n ) += -divtmp_redist_arr( i, j, k, n );
            });

        }

	// Note this sync is needed since ComputeEdgeState() contains temporaries
	// Not sure it's really needed when known_edgestate==true
        Gpu::streamSynchronize();  // otherwise we might be using too much memory
    }

}


void
EBGodunov::ComputeFluxes ( Box const& bx,
                           AMREX_D_DECL( Array4<Real> const& fx,
                                         Array4<Real> const& fy,
                                         Array4<Real> const& fz),
                           AMREX_D_DECL( Array4<Real const> const& umac,
                                         Array4<Real const> const& vmac,
                                         Array4<Real const> const& wmac),
                           AMREX_D_DECL( Array4<Real const> const& xed,
                                         Array4<Real const> const& yed,
                                         Array4<Real const> const& zed),
                           AMREX_D_DECL( Array4<Real const> const& apx,
                                         Array4<Real const> const& apy,
                                         Array4<Real const> const& apz),
                           Geometry const& geom, const int ncomp,
                           Array4<EBCellFlag const> const& flag)
{

    const auto dx = geom.CellSizeArray();

    GpuArray<Real,AMREX_SPACEDIM> area;

#if ( AMREX_SPACEDIM == 3 )
    area[0] = dx[1]*dx[2];
    area[1] = dx[0]*dx[2];
    area[2] = dx[0]*dx[1];
#else
    area[0] = dx[1];
    area[1] = dx[0];
#endif

    //
    //  X flux
    //
    const Box& xbx = amrex::surroundingNodes(bx,0);

    amrex::ParallelFor(xbx, ncomp, [fx, umac, xed, area, apx, flag]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isConnected(-1,0,0))
            fx(i,j,k,n) = xed(i,j,k,n) * umac(i,j,k) * apx(i,j,k) * area[0];
        else
            fx(i,j,k,n) = 0.;
    });

    //
    //  y flux
    //
    const Box& ybx = amrex::surroundingNodes(bx,1);

    amrex::ParallelFor(ybx, ncomp, [fy, vmac, yed, area, apy, flag]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isConnected(0,-1,0))
            fy(i,j,k,n) = yed(i,j,k,n) * vmac(i,j,k) * apy(i,j,k) * area[1];
        else
            fy(i,j,k,n) = 0.;
    });

#if (AMREX_SPACEDIM==3)
    //
    //  z flux
    //
    const Box& zbx = amrex::surroundingNodes(bx,2);

    amrex::ParallelFor(zbx, ncomp, [fz, wmac, zed, area, apz, flag]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isConnected(0,0,-1))
            fz(i,j,k,n) = zed(i,j,k,n) * wmac(i,j,k) * apz(i,j,k) * area[2];
        else
            fz(i,j,k,n) = 0.;
    });
#endif

}



void
EBGodunov::ComputeDivergence ( Box const& bx,
                               Array4<Real> const& div,
                               AMREX_D_DECL( Array4<Real const> const& fx,
                                             Array4<Real const> const& fy,
                                             Array4<Real const> const& fz),
                               Array4<Real const> const& vfrac,
                               const int ncomp, Geometry const& geom )
{

    const auto dxinv = geom.InvCellSizeArray();

#if (AMREX_SPACEDIM==3)
    Real qvol = dxinv[0] * dxinv[1] * dxinv[2];
#else
    Real qvol = dxinv[0] * dxinv[1];
#endif

    // Return -div because reinitialization algo operates on it
    // instead of operatin on div
    amrex::ParallelFor(bx, ncomp,[=]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if ( vfrac(i,j,k) > 0.)
        {
            div(i,j,k,n) =  - qvol/vfrac(i,j,k) *
                (
                         fx(i+1,j,k,n) -  fx(i,j,k,n)
                       + fy(i,j+1,k,n) -  fy(i,j,k,n)
#if (AMREX_SPACEDIM==3)
                       + fz(i,j,k+1,n) -  fz(i,j,k,n)
#endif
                );
        }
        else
        {
            div(i,j,k,n) = 0.0;
        }

    });
}
