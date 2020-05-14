#include <iamr_mol.H>
#include <iamr_constants.H>
#include <NS_util.H>
#ifdef AMREX_USE_EB
#include <AMReX_MultiCutFab.H>
#endif

using namespace amrex;


void
MOL::ComputeAofs ( MultiFab& aofs, int aofs_comp, int ncomp,
                   MultiFab const& state, int state_comp,
                   D_DECL( MultiFab const& umac,
                           MultiFab const& vmac,
                           MultiFab const& wmac),
                   D_DECL( MultiFab& xedge,
                           MultiFab& yedge,
                           MultiFab& zedge),
                   int  edge_comp,
                   bool known_edgestate,
                   D_DECL( MultiFab& xfluxes,
                           MultiFab& yfluxes,
                           MultiFab& zfluxes),
                   int fluxes_comp,
                   Vector<BCRec> const& bcs,
                   Geometry const&  geom )
{
    BL_PROFILE("MOL::ComputeAofs()");

    AMREX_ALWAYS_ASSERT(aofs.nComp()  >= aofs_comp  + ncomp);
    AMREX_ALWAYS_ASSERT(state.nComp() >= state_comp + ncomp);
    D_TERM( AMREX_ALWAYS_ASSERT(xedge.nComp() >= edge_comp  + ncomp);,
            AMREX_ALWAYS_ASSERT(yedge.nComp() >= edge_comp  + ncomp);,
            AMREX_ALWAYS_ASSERT(zedge.nComp() >= edge_comp  + ncomp););
    D_TERM( AMREX_ALWAYS_ASSERT(xfluxes.nComp() >= fluxes_comp  + ncomp);,
            AMREX_ALWAYS_ASSERT(yfluxes.nComp() >= fluxes_comp  + ncomp);,
            AMREX_ALWAYS_ASSERT(zfluxes.nComp() >= fluxes_comp  + ncomp););
    D_TERM( AMREX_ALWAYS_ASSERT(xfluxes.nGrow() == xedge.nGrow());,
            AMREX_ALWAYS_ASSERT(yfluxes.nGrow() == yedge.nGrow());,
            AMREX_ALWAYS_ASSERT(zfluxes.nGrow() == zedge.nGrow()););

#ifdef AMREX_USE_EB
    // We need at least two ghost nodes for redistribution
    AMREX_ALWAYS_ASSERT(state.hasEBFabFactory());
    D_TERM( AMREX_ALWAYS_ASSERT(xedge.nGrow() >= 2 );,
            AMREX_ALWAYS_ASSERT(yedge.nGrow() >= 2 );,
            AMREX_ALWAYS_ASSERT(zedge.nGrow() >= 2 ););
    auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(state.Factory());
#endif

    Box  const& domain = geom.Domain();

    Gpu::DeviceVector<BCRec> bcs_device = convertToDeviceVector(bcs);
    const BCRec* bcs_ptr = bcs_device.dataPtr();

    MFItInfo mfi_info;

    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling(IntVect(D_DECL(1024,1024,1024))).SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aofs,mfi_info); mfi.isValid(); ++mfi)
    {
        auto const& bx = mfi.tilebox();

        D_TERM( Array4<Real> fx = xfluxes.array(mfi,fluxes_comp);,
                Array4<Real> fy = yfluxes.array(mfi,fluxes_comp);,
                Array4<Real> fz = zfluxes.array(mfi,fluxes_comp););

#ifdef AMREX_USE_EB
        // Initialize covered cells
        auto const& flagfab = ebfactory.getMultiEBCellFlagFab()[mfi];
        auto const& flag    = flagfab.const_array();

        if (flagfab.getType(bx) == FabType::covered)
        {
            auto const& aofs_arr = aofs.array(mfi, aofs_comp);
            amrex::ParallelFor(bx, ncomp, [aofs_arr] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                aofs_arr( i, j, k, n ) = covered_val;
            });

            const Box&  xbx = amrex::surroundingNodes(bx,0);
            amrex::ParallelFor(xbx, ncomp, [fx] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                fx( i, j, k, n ) = 0.0;
            });

            const Box&  ybx = amrex::surroundingNodes(bx,1);
            amrex::ParallelFor(ybx, ncomp, [fy] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                fy( i, j, k, n ) = 0.0;
            });

#if (AMREX_SPACEDIM==3)
            const Box&  zbx = amrex::surroundingNodes(bx,2);
            amrex::ParallelFor(zbx, ncomp, [fz] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                fz( i, j, k, n ) = 0.0;
            });
#endif

        }
        else
#endif
        {
            D_TERM( Array4<Real> xed = xedge.array(mfi,edge_comp);,
                    Array4<Real> yed = yedge.array(mfi,edge_comp);,
                    Array4<Real> zed = zedge.array(mfi,edge_comp););

            D_TERM( Array4<Real const> u = umac.const_array(mfi);,
                    Array4<Real const> v = vmac.const_array(mfi);,
                    Array4<Real const> w = wmac.const_array(mfi););

#ifdef AMREX_USE_EB
            bool regular = flagfab.getType(amrex::grow(bx,1)) == FabType::regular;

            if (!regular)
            {
                D_TERM( Array4<Real const> fcx = ebfactory.getFaceCent()[0]->const_array(mfi);,
                        Array4<Real const> fcy = ebfactory.getFaceCent()[1]->const_array(mfi);,
                        Array4<Real const> fcz = ebfactory.getFaceCent()[2]->const_array(mfi););

                Array4<Real const> ccc = ebfactory.getCentroid().const_array(mfi);

                // Grown box on which to compute the fluxes and divergence.
                // We need at least two ghost nodes for redistribution
                Box gbx = amrex::grow(bx,xedge.nGrow());

                // Compute edge state if needed
                if (!known_edgestate)
                {
		    Array4<Real const> const q = state.const_array(mfi,state_comp);

		    EB_ComputeEdgeState( gbx, D_DECL(xed,yed,zed), q, ncomp,
                                         D_DECL(u,v,w), domain, bcs,
                                         D_DECL(fcx,fcy,fcz), ccc, flag );
                }

                // Compute fluxes
                EB_ComputeFluxes(gbx, D_DECL(fx,fy,fz), D_DECL(u,v,w), D_DECL(xed,yed,zed), ncomp, flag );

                //
                // Compute divergence and redistribute
                //
                FArrayBox    divtmp(gbx,ncomp);
                Elixir       diveli = divtmp.elixir();
                Array4<Real> divtmp_arr = divtmp.array();

                auto vfrac = ebfactory.getVolFrac().const_array(mfi);

                D_TERM( auto apx = ebfactory.getAreaFrac()[0]->const_array(mfi);,
                        auto apy = ebfactory.getAreaFrac()[1]->const_array(mfi);,
                        auto apz = ebfactory.getAreaFrac()[2]->const_array(mfi); );

                // Compute conservative divergence
                EB_ComputeDivergence(gbx, divtmp_arr, D_DECL(fx,fy,fz), ncomp, geom, flag, vfrac, D_DECL(apx,apy,apz));

                // Redistribute
                Redistribute(bx, ncomp, aofs.array(mfi, aofs_comp), divtmp_arr, flag, vfrac, geom);

                // Weight fluxes by area and copy them  to output
                EB_AreaWeightFluxes(bx, D_DECL(fx,fy,fz), D_DECL(apx,apy,apz), ncomp, geom.CellSize());

            }
            else
#endif
            {
                Box gbx = amrex::grow(bx,xedge.nGrow());

                // Compute edge state if needed
                if (!known_edgestate)
                {
                    Array4<Real const> const q = state.const_array(mfi,state_comp);
                    ComputeEdgeState( gbx, D_DECL( xed, yed, zed ), q, ncomp,
                                      D_DECL( u, v, w ), domain, bcs );

                }

                // Compute fluxes
                ComputeFluxes(gbx, D_DECL(fx,fy,fz), D_DECL(u,v,w), D_DECL(xed,yed,zed), ncomp );

                // Compute divergence
                ComputeDivergence(bx, aofs.array(mfi, aofs_comp), D_DECL(fx,fy,fz), ncomp, geom);

                // Weight fluxes by area and copy them  to output
                AreaWeightFluxes(bx, D_DECL(fx,fy,fz), ncomp, geom.CellSize());

            }
        }
    }
}




void
MOL::ComputeSyncAofs ( MultiFab& aofs, int aofs_comp, int ncomp,
                       MultiFab const& state, int state_comp,
                       D_DECL( MultiFab const& umac,
                               MultiFab const& vmac,
                               MultiFab const& wmac),
                       D_DECL( MultiFab const& ucorr,
                               MultiFab const& vcorr,
                               MultiFab const& wcorr),
                       D_DECL( MultiFab& xedge,
                               MultiFab& yedge,
                               MultiFab& zedge),
                       int  edge_comp,
                       bool known_edgestate,
                       D_DECL( MultiFab& xfluxes,
                               MultiFab& yfluxes,
                               MultiFab& zfluxes),
                       int fluxes_comp,
                       Vector<BCRec> const& bcs,
                       Geometry const&  geom )
{
    BL_PROFILE("MOL::ComputeSyncAofs()");

    AMREX_ALWAYS_ASSERT(state.nComp() >= state_comp + ncomp);

#ifdef AMREX_USE_EB
    AMREX_ALWAYS_ASSERT(state.hasEBFabFactory());
    auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(state.Factory());
#endif

    Box  const& domain = geom.Domain();

    Gpu::DeviceVector<BCRec> bcs_device = convertToDeviceVector(bcs);
    const BCRec* bcs_ptr = bcs_device.dataPtr();

    MFItInfo mfi_info;

    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling(IntVect(D_DECL(1024,1024,1024))).SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aofs,mfi_info); mfi.isValid(); ++mfi)
    {
        auto const& bx = mfi.tilebox();

        D_TERM( Array4<Real> fx = xfluxes.array(mfi,fluxes_comp);,
                Array4<Real> fy = yfluxes.array(mfi,fluxes_comp);,
                Array4<Real> fz = zfluxes.array(mfi,fluxes_comp););

#ifdef AMREX_USE_EB
        // Initialize covered cells
        auto const& flagfab = ebfactory.getMultiEBCellFlagFab()[mfi];
        auto const& flag    = flagfab.const_array();

        if (flagfab.getType(bx) == FabType::covered)
        {
            auto const& aofs_arr = aofs.array(mfi, aofs_comp);
            amrex::ParallelFor(bx, ncomp, [aofs_arr] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                aofs_arr( i, j, k, n ) = covered_val;
            });

            const Box&  xbx = amrex::surroundingNodes(bx,0);
            amrex::ParallelFor(xbx, ncomp, [fx] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                fx( i, j, k, n ) = 0.0;
            });

            const Box&  ybx = amrex::surroundingNodes(bx,1);
            amrex::ParallelFor(ybx, ncomp, [fy] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                fy( i, j, k, n ) = 0.0;
            });

#if (AMREX_SPACEDIM==3)
            const Box&  zbx = amrex::surroundingNodes(bx,2);
            amrex::ParallelFor(zbx, ncomp, [fz] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                fz( i, j, k, n ) = 0.0;
            });
#endif

        }
        else
#endif
        {
            D_TERM( Array4<Real> xed = xedge.array(mfi,edge_comp);,
                    Array4<Real> yed = yedge.array(mfi,edge_comp);,
                    Array4<Real> zed = zedge.array(mfi,edge_comp););

            D_TERM( Array4<Real const> uc = ucorr.const_array(mfi);,
                    Array4<Real const> vc = vcorr.const_array(mfi);,
                    Array4<Real const> wc = wcorr.const_array(mfi););

#ifdef AMREX_USE_EB
            bool regular = flagfab.getType(amrex::grow(bx,1)) == FabType::regular;

            if (!regular)
            {
                D_TERM( Array4<Real const> fcx = ebfactory.getFaceCent()[0]->const_array(mfi);,
                        Array4<Real const> fcy = ebfactory.getFaceCent()[1]->const_array(mfi);,
                        Array4<Real const> fcz = ebfactory.getFaceCent()[2]->const_array(mfi););

                Array4<Real const> ccc = ebfactory.getCentroid().const_array(mfi);

                // Grown box on which to compute the fluxes and divergence.
                // We need at least two ghost nodes for redistribution
                Box gbx = amrex::grow(bx,2);

                // Compute edge state if needed
                if (!known_edgestate)
                {
		    Array4<Real const> const q = state.const_array(mfi,state_comp);

		    D_TERM( Array4<Real const> u = umac.const_array(mfi);,
			    Array4<Real const> v = vmac.const_array(mfi);,
			    Array4<Real const> w = wmac.const_array(mfi););

                    EB_ComputeEdgeState( gbx, D_DECL(xed,yed,zed), q, ncomp,
                                         D_DECL(u,v,w), domain, bcs,
                                         D_DECL(fcx,fcy,fcz), ccc, flag );
                }

                // Compute fluxes
                EB_ComputeFluxes(gbx, D_DECL(fx,fy,fz), D_DECL(uc,vc,wc), D_DECL(xed,yed,zed), ncomp, flag );

                //
                // Compute divergence and redistribute
                //
                FArrayBox    divtmp(gbx,2*ncomp);
                Elixir       diveli = divtmp.elixir();
                Array4<Real> divtmp_arr        = divtmp.array(0);
                Array4<Real> divtmp_redist_arr = divtmp.array(ncomp);

                auto vfrac = ebfactory.getVolFrac().const_array(mfi);

                D_TERM( auto apx = ebfactory.getAreaFrac()[0]->const_array(mfi);,
                        auto apy = ebfactory.getAreaFrac()[1]->const_array(mfi);,
                        auto apz = ebfactory.getAreaFrac()[2]->const_array(mfi); );

                // Compute conservative divergence
                EB_ComputeDivergence(gbx, divtmp_arr, D_DECL(fx,fy,fz), ncomp, geom, flag, vfrac, D_DECL(apx,apy,apz));

                // Redistribute
                Redistribute(bx, ncomp, divtmp_redist_arr, divtmp_arr, flag, vfrac, geom);

                // Sum contribution to sync aofs
                auto const& aofs_arr = aofs.array(mfi, aofs_comp);

                amrex::ParallelFor(bx, ncomp, [aofs_arr, divtmp_redist_arr]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    aofs_arr( i, j, k, n ) += divtmp_redist_arr( i, j, k, n );
                });

                // Weight fluxes by area and copy them  to output
                EB_AreaWeightFluxes(bx, D_DECL(fx,fy,fz), D_DECL(apx,apy,apz), ncomp, geom.CellSize());

            }
            else
#endif
            {
                // Compute edge state if needed
                if (!known_edgestate)
                {
		    Array4<Real const> const q = state.const_array(mfi,state_comp);

		    D_TERM( Array4<Real const> u = umac.const_array(mfi);,
			    Array4<Real const> v = vmac.const_array(mfi);,
			    Array4<Real const> w = wmac.const_array(mfi););

                    ComputeEdgeState( bx, D_DECL( xed, yed, zed ), q, ncomp,
                                      D_DECL( u, v, w ), domain, bcs );

                }

                // Compute fluxes
                ComputeFluxes(bx, D_DECL(fx,fy,fz), D_DECL(uc,vc,wc), D_DECL(xed,yed,zed), ncomp );

                // Compute divergence
                FArrayBox    divtmp(bx,ncomp);
                Array4<Real> divtmp_arr = divtmp.array(0);
                ComputeDivergence(bx, divtmp_arr, D_DECL(fx,fy,fz), ncomp, geom);

                // Sum contribution to sync aofs
                auto const& aofs_arr = aofs.array(mfi, aofs_comp);

                amrex::ParallelFor(bx, ncomp, [aofs_arr, divtmp_arr]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    aofs_arr( i, j, k, n ) += divtmp_arr( i, j, k, n );
                });

                // Weight fluxes by area and copy them  to output
                AreaWeightFluxes(bx, D_DECL(fx,fy,fz), ncomp, geom.CellSize());

            }
        }
    }
}


void
MOL::ComputeFluxes ( Box const& bx,
                     D_DECL( Array4<Real> const& fx,
                             Array4<Real> const& fy,
                             Array4<Real> const& fz),
                     D_DECL( Array4<Real const> const& umac,
                             Array4<Real const> const& vmac,
                             Array4<Real const> const& wmac),
                     D_DECL( Array4<Real const> const& xedge,
                             Array4<Real const> const& yedge,
                             Array4<Real const> const& zedge),
                     int ncomp )
{
    //
    //  X flux
    //
    const Box& xbx = amrex::surroundingNodes(bx,0);

    amrex::ParallelFor(xbx, ncomp, [fx, umac, xedge]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fx(i,j,k,n) = xedge(i,j,k,n) * umac(i,j,k);
    });

    //
    //  y flux
    //
    const Box& ybx = amrex::surroundingNodes(bx,1);

    amrex::ParallelFor(ybx, ncomp, [fy, vmac, yedge]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fy(i,j,k,n) = yedge(i,j,k,n) * vmac(i,j,k);
    });

#if (AMREX_SPACEDIM==3)
    //
    //  z flux
    //
    const Box& zbx = amrex::surroundingNodes(bx,2);

    amrex::ParallelFor(zbx, ncomp, [fz, wmac, zedge]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fz(i,j,k,n) = zedge(i,j,k,n) * wmac(i,j,k);
    });

#endif

}


void
MOL::ComputeDivergence ( Box const& bx,
                         Array4<Real> const& div,
                         D_DECL( Array4<Real const> const& fx,
                                 Array4<Real const> const& fy,
                                 Array4<Real const> const& fz),
                         int ncomp, Geometry const& geom )
{
    const auto dxinv = geom.InvCellSizeArray();
    amrex::ParallelFor(bx, ncomp, [ div, D_DECL(fx,fy,fz), dxinv]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        div(i,j,k,n) = dxinv[0] * (fx(i+1,j,k,n) - fx(i,j,k,n))
            +          dxinv[1] * (fy(i,j+1,k,n) - fy(i,j,k,n))
#if (AMREX_SPACEDIM==3)
            +          dxinv[2] * (fz(i,j,k+1,n) - fz(i,j,k,n))
#endif
            ;
    });

}


void
MOL::AreaWeightFluxes( Box const& bx,
                       D_DECL( Array4<Real> const& fx,
                               Array4<Real> const& fy,
                               Array4<Real> const& fz),
                       int ncomp,
                       const Real* dx )
{

    Real  area[AMREX_SPACEDIM];
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

    amrex::ParallelFor(xbx, ncomp, [fx,area]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fx(i,j,k,n) *= area[0];
    });

    //
    //  y flux
    //
    const Box& ybx = amrex::surroundingNodes(bx,1);

    amrex::ParallelFor(ybx, ncomp, [fy,area]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fy(i,j,k,n) *= area[1];
    });

#if (AMREX_SPACEDIM==3)
    //
    //  z flux
    //
    const Box& zbx = amrex::surroundingNodes(bx,2);

    amrex::ParallelFor(zbx, ncomp, [fz,area]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fz(i,j,k,n) *= area[2];
    });

#endif

}

//
// ============================  EB-only routines ============================
//

#ifdef AMREX_USE_EB
void
MOL::EB_ComputeFluxes ( Box const& bx,
                        D_DECL( Array4<Real> const& fx,
                                Array4<Real> const& fy,
                                Array4<Real> const& fz),
                        D_DECL( Array4<Real const> const& umac,
                                Array4<Real const> const& vmac,
                                Array4<Real const> const& wmac),
                        D_DECL( Array4<Real const> const& xedge,
                                Array4<Real const> const& yedge,
                                Array4<Real const> const& zedge),
                        int ncomp,
                        Array4<EBCellFlag const> const& flag)
{
    //
    //  X flux
    //
    const Box& xbx = amrex::surroundingNodes(bx,0);

    amrex::ParallelFor(xbx, ncomp, [fx, umac, xedge, flag]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isConnected(-1,0,0))
        {
            fx(i,j,k,n) = xedge(i,j,k,n) * umac(i,j,k);
        }
        else
        {
            fx(i,j,k,n) = 0.0;
        }
    });

    //
    //  y flux
    //
    const Box& ybx = amrex::surroundingNodes(bx,1);

    amrex::ParallelFor(ybx, ncomp, [fy, vmac, yedge, flag]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isConnected(0,-1,0))
        {
            fy(i,j,k,n) = yedge(i,j,k,n) * vmac(i,j,k);
        }
        else
        {
            fy(i,j,k,n) = 0.0;
        }
    });

#if (AMREX_SPACEDIM == 3)

    //
    //  z flux
    //
    const Box& zbx = amrex::surroundingNodes(bx,2);

    amrex::ParallelFor(zbx, ncomp, [fz, wmac, zedge, flag]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isConnected(0,0,-1))
        {
            fz(i,j,k,n) = zedge(i,j,k,n) * wmac(i,j,k);
        }
        else
        {
            fz(i,j,k,n) = 0.0;
        }
    });

#endif

}


void
MOL::EB_ComputeDivergence ( Box const& bx,
                            Array4<Real> const& div,
                            D_DECL( Array4<Real const> const& fx,
                                    Array4<Real const> const& fy,
                                    Array4<Real const> const& fz),
                            int ncomp, Geometry const& geom,
                            Array4<EBCellFlag const> const& flag,
                            Array4<Real const> const& vfrac,
                            D_DECL( Array4<Real const> const& apx,
                                    Array4<Real const> const& apy,
                                    Array4<Real const> const& apz ) )
{
    const auto dxinv = geom.InvCellSizeArray();
    const auto dbox  = geom.growPeriodicDomain(2);
    amrex::ParallelFor(bx, ncomp, [ div, D_DECL(fx,fy,fz), dxinv, dbox, flag, vfrac, D_DECL(apx,apy,apz) ]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (!dbox.contains(IntVect(D_DECL(i,j,k))) or flag(i,j,k).isCovered())
        {
            div(i,j,k,n) = 0.0;
        }
        else if (flag(i,j,k).isRegular())
        {
            div(i,j,k,n) = dxinv[0] * (fx(i+1,j,k,n) - fx(i,j,k,n))
                +          dxinv[1] * (fy(i,j+1,k,n) - fy(i,j,k,n))
#if (AMREX_SPACEDIM==3)
                +          dxinv[2] * (fz(i,j,k+1,n) - fz(i,j,k,n))
#endif
                ;
        }
        else
        {
            div(i,j,k,n) = (1.0/vfrac(i,j,k)) *
                (       dxinv[0] * (apx(i+1,j,k)*fx(i+1,j,k,n) - apx(i,j,k)*fx(i,j,k,n))
                      + dxinv[1] * (apy(i,j+1,k)*fy(i,j+1,k,n) - apy(i,j,k)*fy(i,j,k,n))
#if (AMREX_SPACEDIM==3)
                      + dxinv[2] * (apz(i,j,k+1)*fz(i,j,k+1,n) - apz(i,j,k)*fz(i,j,k,n))
#endif
                    );

        }
    });

}

void
MOL::Redistribute (  Box const& bx, int ncomp,
                     Array4<Real> const& div,
                     Array4<Real const> const& div_in,
                     Array4<EBCellFlag const> const& flag,
                     Array4<Real const> const& vfrac,
                     Geometry const& geom )
{
    const Box dbox = geom.growPeriodicDomain(2);

    Box const& bxg1 = amrex::grow(bx,1);
    Box const& bxg2 = amrex::grow(bx,2);

    // Temporaries
    FArrayBox scratch(bxg2,2*ncomp+1);
    Elixir eli = scratch.elixir();
    Array4<Real>  tmp  = scratch.array(0);
    Array4<Real>  delm = scratch.array(ncomp);
    Array4<Real>  wgt  = scratch.array(2*ncomp);


    // Weight by EB volume fraction
    amrex::ParallelFor(bxg2, [wgt,dbox,vfrac]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        wgt(i,j,k) = (dbox.contains(IntVect(D_DECL(i,j,k)))) ? vfrac(i,j,k) : 0.0;
    });

    amrex::ParallelFor(bxg1, ncomp, [flag, dbox, vfrac, div_in, tmp, delm, wgt]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isSingleValued()) {
            Real vtot = 0.0;
            Real divnc = 0.0;
#if (AMREX_SPACEDIM==3)
            for (int kk = -1; kk <= 1; ++kk)
#else
                const int kk = 0;
#endif
            {
                for (int jj = -1; jj <= 1; ++jj)
                {
                    for (int ii = -1; ii <= 1; ++ii)
                    {
                        if ( (ii != 0 or jj != 0 or kk != 0) and
                             flag(i,j,k).isConnected(ii,jj,kk) and
                             dbox.contains(IntVect(D_DECL(i+ii,j+jj,k+kk))))
                        {
                            Real wted_vf = vfrac(i+ii,j+jj,k+kk) * wgt(i+ii,j+jj,k+kk);
                            vtot += wted_vf;
                            divnc += wted_vf * div_in(i+ii,j+jj,k+kk,n);
                        }
                    }
                }
            }
            divnc /= (vtot + 1.e-80);
            Real optmp = (1.0-vfrac(i,j,k))*(divnc-div_in(i,j,k,n));
            tmp(i,j,k,n) = optmp;
            delm(i,j,k,n) = -vfrac(i,j,k)*optmp;
        }
        else
        {
            tmp(i,j,k,n) = 0.0;
        }
    });

    amrex::ParallelFor(bxg1 & dbox, ncomp, [flag, vfrac, wgt, bx, tmp, delm]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isSingleValued())
        {
            Real wtot = 0.0;
#if (AMREX_SPACEDIM==3)
            for (int kk = -1; kk <= 1; ++kk)
#else
                const int kk = 0;
#endif
            {
                for (int jj = -1; jj <= 1; ++jj)
                {
                    for (int ii = -1; ii <= 1; ++ii)
                    {
                        if ((ii != 0 or jj != 0 or kk != 0) and
                            flag(i,j,k).isConnected(ii,jj,kk))
                        {
                            wtot += vfrac(i+ii,j+jj,k+kk) * wgt(i+ii,j+jj,k+kk);
                        }
                    }
                }
            }

            wtot = 1.0/(wtot+1.e-80);

            Real dtmp = delm(i,j,k,n) * wtot;
#if (AMREX_SPACEDIM==3)
            for (int kk = -1; kk <= 1; ++kk)
#endif
            {
                for (int jj = -1; jj <= 1; ++jj)
                {
                    for (int ii = -1; ii <= 1; ++ii)
                    {
                        if ((ii != 0 or jj != 0 or kk != 0) and
                            bx.contains(IntVect(D_DECL(i+ii,j+jj,k+kk))) and
                            flag(i,j,k).isConnected(ii,jj,kk))
                        {
                            Gpu::Atomic::Add(&tmp(i+ii,j+jj,k+kk,n), dtmp*wgt(i+ii,j+jj,k+kk));
                        }
                    }
                }
            }
        }
    });

    amrex::ParallelFor(bx, ncomp, [div, div_in, tmp]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        div(i,j,k,n) = div_in(i,j,k,n) + tmp(i,j,k,n);
    });
}


void
MOL::EB_AreaWeightFluxes ( Box const& bx,
                           D_DECL( Array4<Real> const& fx,
                                   Array4<Real> const& fy,
                                   Array4<Real> const& fz),
                           D_DECL( Array4<Real const> const& apx,
                                   Array4<Real const> const& apy,
                                   Array4<Real const> const& apz ),
                           int ncomp,
                           const Real* dx )
{

    Real  area[AMREX_SPACEDIM];
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

    amrex::ParallelFor(xbx, ncomp, [fx,area,apx]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fx(i,j,k,n) *= area[0]*apx(i,j,k);
    });

    //
    //  y flux
    //
    const Box& ybx = amrex::surroundingNodes(bx,1);

    amrex::ParallelFor(ybx, ncomp, [fy,area,apy]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fy(i,j,k,n) *= area[1]*apy(i,j,k);
    });

#if (AMREX_SPACEDIM==3)
    //
    //  z flux
    //
    const Box& zbx = amrex::surroundingNodes(bx,2);

    amrex::ParallelFor(zbx, ncomp, [fz,area,apz]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fz(i,j,k,n) *= area[2]*apz(i,j,k);
    });

#endif

}




#endif
