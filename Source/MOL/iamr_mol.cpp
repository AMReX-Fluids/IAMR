#include <iamr_mol.H>
#include <iamr_constants.H>
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
    AMREX_ALWAYS_ASSERT(aofs.nGrow() == 0);
    D_TERM( AMREX_ALWAYS_ASSERT(xfluxes.nGrow() == xedge.nGrow());,
            AMREX_ALWAYS_ASSERT(yfluxes.nGrow() == yedge.nGrow());,
            AMREX_ALWAYS_ASSERT(zfluxes.nGrow() == zedge.nGrow()););
    if ( !known_edgestate )
      // To compute edge states, need at least 2 more ghost cells in state than in
      //  xedge 
      AMREX_ALWAYS_ASSERT(state.nGrow() >= xedge.nGrow()+2);

#ifdef AMREX_USE_EB
    AMREX_ALWAYS_ASSERT(state.hasEBFabFactory());
    // We need at least two ghost nodes for redistribution
    D_TERM( AMREX_ALWAYS_ASSERT(xedge.nGrow() >= 2 );,
            AMREX_ALWAYS_ASSERT(yedge.nGrow() >= 2 );,
            AMREX_ALWAYS_ASSERT(zedge.nGrow() >= 2 ););

    auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(state.Factory());
#endif

    Box  const& domain = geom.Domain();

    MFItInfo mfi_info;

    if (Gpu::notInLaunchRegion())  mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aofs,mfi_info); mfi.isValid(); ++mfi)
    {
        auto const& bx = mfi.tilebox();
	int ng_f = xfluxes.nGrow();
	auto const& gtbx = mfi.growntilebox(ng_f);
	D_TERM( const Box& xbx = mfi.grownnodaltilebox(0,ng_f);,
		const Box& ybx = mfi.grownnodaltilebox(1,ng_f);,
		const Box& zbx = mfi.grownnodaltilebox(2,ng_f); );

        D_TERM( Array4<Real> fx = xfluxes.array(mfi,fluxes_comp);,
                Array4<Real> fy = yfluxes.array(mfi,fluxes_comp);,
                Array4<Real> fz = zfluxes.array(mfi,fluxes_comp););

	D_TERM( Array4<Real> xed = xedge.array(mfi,edge_comp);,
		Array4<Real> yed = yedge.array(mfi,edge_comp);,
		Array4<Real> zed = zedge.array(mfi,edge_comp););

#ifdef AMREX_USE_EB
        // Initialize covered cells
        auto const& flagfab = ebfactory.getMultiEBCellFlagFab()[mfi];
        auto const& flag    = flagfab.const_array();

        if (flagfab.getType(gtbx) == FabType::covered)
        {
            auto const& aofs_arr = aofs.array(mfi, aofs_comp);
            amrex::ParallelFor(bx, ncomp, [aofs_arr] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                aofs_arr( i, j, k, n ) = covered_val;
            });

            amrex::ParallelFor(xbx, ncomp, [fx,xed] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                fx( i, j, k, n ) = 0.0;
		xed( i, j, k, n ) = 0.0;
            });

            amrex::ParallelFor(ybx, ncomp, [fy,yed] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                fy( i, j, k, n ) = 0.0;
		yed( i, j, k, n ) = 0.0;
            });

#if (AMREX_SPACEDIM==3)
            amrex::ParallelFor(zbx, ncomp, [fz,zed] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                fz( i, j, k, n ) = 0.0;
		zed( i, j, k, n ) = 0.0;
            });
#endif

        }
        else
#endif
        {
            D_TERM( Array4<Real const> u = umac.const_array(mfi);,
                    Array4<Real const> v = vmac.const_array(mfi);,
                    Array4<Real const> w = wmac.const_array(mfi););

	    // Grown box on which to compute the edge states and fluxes for regular boxes
	    Box gbx = mfi.growntilebox(ng_f);

	    Box tmpbox = amrex::surroundingNodes(gbx);
	    // Space for fluxes
	    int tmpcomp = ncomp*AMREX_SPACEDIM;

#ifdef AMREX_USE_EB
	    // PeleLM needs valid flux on all ng_f cells. If !known_edgestate, need 2 additional
	    //  cells in state to compute the slopes needed to compute the edge state.
	    int halo = known_edgestate ? 0 : ng_f+2;
            bool regular = flagfab.getType(amrex::grow(bx,halo)) == FabType::regular;
	    if (!regular) {
	      // Grown box on which to compute the edge states and fluxes for EB containing boxes
	      // need at least 2 filled ghost cells all around for redistribution
	      gbx = amrex::grow(bx,ng_f);
	      tmpbox = amrex::surroundingNodes(gbx);
	      // Not sure if we really need 3(incflo) here or 2
	      int ng_diff = 3-ng_f;
	      if ( ng_diff>0 )  
		tmpbox.grow(ng_diff);

	      // Add space for the temporaries needed by Redistribute
#if (AMREX_SPACEDIM == 3)	      
	      tmpcomp += ncomp;
#else
	      tmpcomp += 2*ncomp;
#endif
	    }
#endif
	    FArrayBox tmpfab(tmpbox, tmpcomp);
	    Elixir eli = tmpfab.elixir();

	    D_TERM( Array4<Real> fxtmp = tmpfab.array(0);,
		    Array4<Real> fytmp = tmpfab.array(ncomp);,
		    Array4<Real> fztmp = tmpfab.array(ncomp*2););
	
#ifdef AMREX_USE_EB
            if (!regular)
            {
                D_TERM( Array4<Real const> fcx = ebfactory.getFaceCent()[0]->const_array(mfi);,
                        Array4<Real const> fcy = ebfactory.getFaceCent()[1]->const_array(mfi);,
                        Array4<Real const> fcz = ebfactory.getFaceCent()[2]->const_array(mfi););

                Array4<Real const> ccc = ebfactory.getCentroid().const_array(mfi);

                // Compute edge state if needed
                if (!known_edgestate)
                {
		    Array4<Real const> const q = state.const_array(mfi,state_comp);

		    EB_ComputeEdgeState( gbx, D_DECL(xed,yed,zed), q, ncomp,
                                         D_DECL(u,v,w), domain, bcs,
                                         D_DECL(fcx,fcy,fcz), ccc, flag );
                }

                // Compute fluxes
                EB_ComputeFluxes(gbx, D_DECL(fxtmp,fytmp,fztmp), D_DECL(u,v,w),
				 D_DECL(xed,yed,zed), ncomp, flag );

                //
                // Compute divergence and redistribute
                //
		// div at ncomp*3 to make space for the 3 redistribute temporaries
                Array4<Real> divtmp_arr = tmpfab.array(ncomp*3);

                auto vfrac = ebfactory.getVolFrac().const_array(mfi);

                D_TERM( auto apx = ebfactory.getAreaFrac()[0]->const_array(mfi);,
                        auto apy = ebfactory.getAreaFrac()[1]->const_array(mfi);,
                        auto apz = ebfactory.getAreaFrac()[2]->const_array(mfi); );

                // Compute conservative divergence
                // Redistribute needs 2 ghost cells in div
	        Box g2bx = amrex::grow(bx,2);
                EB_ComputeDivergence(g2bx, divtmp_arr, D_DECL(fxtmp,fytmp,fztmp), ncomp,
				     geom, flag, vfrac, D_DECL(apx,apy,apz));

                // Weight fluxes by area
		// Do this before Redistribute so that we can reuse the tmp space
                EB_AreaWeightFluxes(D_DECL(xbx, ybx, zbx), D_DECL(fx,fy,fz),
				    D_DECL(fxtmp,fytmp,fztmp), D_DECL(apx,apy,apz),
				    ncomp, geom.CellSize());

                // Redistribute
		Array4<Real> scratch = tmpfab.array(0);
                Redistribute(bx, ncomp, aofs.array(mfi, aofs_comp), divtmp_arr, scratch,
			     flag, vfrac, geom);

            }
            else
#endif
            {
                // Compute edge state if needed
                if (!known_edgestate)
                {
                    Array4<Real const> const q = state.const_array(mfi,state_comp);
                    ComputeEdgeState( gbx, D_DECL( xed, yed, zed ), q, ncomp,
                                      D_DECL( u, v, w ), domain, bcs );

                }

                // Compute fluxes
                ComputeFluxes(gbx, D_DECL(fxtmp,fytmp,fztmp), D_DECL(u,v,w),
			      D_DECL(xed,yed,zed), ncomp );

                // Compute divergence
                ComputeDivergence(bx, aofs.array(mfi, aofs_comp), D_DECL(fxtmp,fytmp,fztmp),
				  ncomp, geom);

                // Weight fluxes by area and copy them  to output
                AreaWeightFluxes(D_DECL(xbx, ybx, zbx), D_DECL(fx,fy,fz),
				 D_DECL(fxtmp,fytmp,fztmp), ncomp, geom.CellSize());

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
    AMREX_ALWAYS_ASSERT(aofs.nComp()  >= aofs_comp  + ncomp);
    D_TERM( AMREX_ALWAYS_ASSERT(xedge.nComp() >= edge_comp  + ncomp);,
            AMREX_ALWAYS_ASSERT(yedge.nComp() >= edge_comp  + ncomp);,
            AMREX_ALWAYS_ASSERT(zedge.nComp() >= edge_comp  + ncomp););
    D_TERM( AMREX_ALWAYS_ASSERT(xfluxes.nComp() >= fluxes_comp  + ncomp);,
            AMREX_ALWAYS_ASSERT(yfluxes.nComp() >= fluxes_comp  + ncomp);,
            AMREX_ALWAYS_ASSERT(zfluxes.nComp() >= fluxes_comp  + ncomp););
    D_TERM( AMREX_ALWAYS_ASSERT(xfluxes.nGrow() == xedge.nGrow());,
            AMREX_ALWAYS_ASSERT(yfluxes.nGrow() == yedge.nGrow());,
            AMREX_ALWAYS_ASSERT(zfluxes.nGrow() == zedge.nGrow()););
    if ( !known_edgestate )
      // To compute edge states, need at least 2 more ghost cells in state than in
      //  xedge 
      AMREX_ALWAYS_ASSERT(state.nGrow() >= xedge.nGrow()+2);

#ifdef AMREX_USE_EB
    AMREX_ALWAYS_ASSERT(state.hasEBFabFactory());
    // We need at least two ghost nodes for redistribution
    D_TERM( AMREX_ALWAYS_ASSERT(xedge.nGrow() >= 2 );,
            AMREX_ALWAYS_ASSERT(yedge.nGrow() >= 2 );,
            AMREX_ALWAYS_ASSERT(zedge.nGrow() >= 2 ););

    auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(state.Factory());
#endif

    Box  const& domain = geom.Domain();

    MFItInfo mfi_info;

    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(aofs,mfi_info); mfi.isValid(); ++mfi)
    {
        auto const& bx = mfi.tilebox();
	D_TERM( const Box& xbx = mfi.nodaltilebox(0);,
		const Box& ybx = mfi.nodaltilebox(1);,
		const Box& zbx = mfi.nodaltilebox(2); );

        D_TERM( Array4<Real> fx = xfluxes.array(mfi,fluxes_comp);,
                Array4<Real> fy = yfluxes.array(mfi,fluxes_comp);,
                Array4<Real> fz = zfluxes.array(mfi,fluxes_comp););
	
	D_TERM( Array4<Real> xed = xedge.array(mfi,edge_comp);,
		Array4<Real> yed = yedge.array(mfi,edge_comp);,
		Array4<Real> zed = zedge.array(mfi,edge_comp););

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

            amrex::ParallelFor(xbx, ncomp, [fx,xed] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                fx( i, j, k, n ) = 0.0;
		xed( i, j, k, n ) = 0.0;
            });

            amrex::ParallelFor(ybx, ncomp, [fy,yed] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                fy( i, j, k, n ) = 0.0;
		yed( i, j, k, n ) = 0.0;
            });

#if (AMREX_SPACEDIM==3)
            amrex::ParallelFor(zbx, ncomp, [fz,zed] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                fz( i, j, k, n ) = 0.0;
		zed( i, j, k, n ) = 0.0;
            });
#endif
        }
        else
#endif
        {
            D_TERM( Array4<Real const> uc = ucorr.const_array(mfi);,
                    Array4<Real const> vc = vcorr.const_array(mfi);,
                    Array4<Real const> wc = wcorr.const_array(mfi););

	    Box tmpbox = amrex::surroundingNodes(bx);
	    int tmpcomp = ncomp*(AMREX_SPACEDIM+1);
#ifdef AMREX_USE_EB
	    Box gbx = bx;
	    // Need 2 grow cells in state to compute the slopes needed to compute the edge state.
	    int halo = known_edgestate ? 0 : 2;
            bool regular = flagfab.getType(amrex::grow(bx,halo)) == FabType::regular;
	    if (!regular) {
	      // Grown box on which to compute the fluxes and divergence.
	      gbx.grow(2);
	      tmpbox.grow(3);

	      // Add space for the temporaries needed by Redistribute
#if (AMREX_SPACEDIM == 3)	      
	      tmpcomp += ncomp;
#else
	      tmpcomp += 2*ncomp;
#endif
	    }
#endif
	    FArrayBox tmpfab(tmpbox, tmpcomp);
	    Elixir eli = tmpfab.elixir();

	    D_TERM( Array4<Real> fxtmp = tmpfab.array(0);,
		    Array4<Real> fytmp = tmpfab.array(ncomp);,
		    Array4<Real> fztmp = tmpfab.array(ncomp*2););

#ifdef AMREX_USE_EB
            if (!regular)
            {
                D_TERM( Array4<Real const> fcx = ebfactory.getFaceCent()[0]->const_array(mfi);,
                        Array4<Real const> fcy = ebfactory.getFaceCent()[1]->const_array(mfi);,
                        Array4<Real const> fcz = ebfactory.getFaceCent()[2]->const_array(mfi););

                Array4<Real const> ccc = ebfactory.getCentroid().const_array(mfi);

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
                EB_ComputeFluxes(gbx, D_DECL(fxtmp,fytmp,fztmp), D_DECL(uc,vc,wc),
				 D_DECL(xed,yed,zed), ncomp, flag );

                //
                // Compute divergence and redistribute
                //
		// div at ncomp*3 to make space for the 3 redistribute temporaries
                Array4<Real> divtmp_arr = tmpfab.array(ncomp*3);
                Array4<Real> divtmp_redist_arr = tmpfab.array(ncomp*4);

                auto vfrac = ebfactory.getVolFrac().const_array(mfi);

                D_TERM( auto apx = ebfactory.getAreaFrac()[0]->const_array(mfi);,
                        auto apy = ebfactory.getAreaFrac()[1]->const_array(mfi);,
                        auto apz = ebfactory.getAreaFrac()[2]->const_array(mfi); );

                // Compute conservative divergence
                EB_ComputeDivergence(gbx, divtmp_arr, D_DECL(fxtmp,fytmp,fztmp), ncomp,
				     geom, flag, vfrac, D_DECL(apx,apy,apz));

		// Weight fluxes by area
		// Do this before Redistribute so that we can reuse the tmp space
                EB_AreaWeightFluxes(D_DECL(xbx, ybx, zbx), D_DECL(fx,fy,fz),
				    D_DECL(fxtmp,fytmp,fztmp), D_DECL(apx,apy,apz),
				    ncomp, geom.CellSize());

                // Redistribute
		Array4<Real> scratch = tmpfab.array(0);
                Redistribute(bx, ncomp, divtmp_redist_arr, divtmp_arr, scratch,
			     flag, vfrac, geom);

                // Sum contribution to sync aofs
                auto const& aofs_arr = aofs.array(mfi, aofs_comp);

                amrex::ParallelFor(bx, ncomp, [aofs_arr, divtmp_redist_arr]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    aofs_arr( i, j, k, n ) += divtmp_redist_arr( i, j, k, n );
                });

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
                ComputeFluxes(bx, D_DECL(fxtmp,fytmp,fztmp), D_DECL(uc,vc,wc),
			      D_DECL(xed,yed,zed), ncomp );

                // Compute divergence
                Array4<Real> divtmp_arr = tmpfab.array(ncomp*AMREX_SPACEDIM);
                ComputeDivergence(bx, divtmp_arr, D_DECL(fxtmp,fytmp,fztmp), ncomp, geom);

                // Sum contribution to sync aofs
                auto const& aofs_arr = aofs.array(mfi, aofs_comp);

                amrex::ParallelFor(bx, ncomp, [aofs_arr, divtmp_arr]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    aofs_arr( i, j, k, n ) += divtmp_arr( i, j, k, n );
                });

                // Weight fluxes by area and copy them  to output
                AreaWeightFluxes(D_DECL(xbx, ybx, zbx), D_DECL(fx,fy,fz),
				 D_DECL(fxtmp,fytmp,fztmp), ncomp, geom.CellSize());

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
MOL::AreaWeightFluxes( D_DECL( Box const& xbx,
			       Box const& ybx,
			       Box const& zbx),
                       D_DECL( Array4<Real> const& fx,
                               Array4<Real> const& fy,
                               Array4<Real> const& fz),
		       D_DECL( Array4<Real> const& fxtmp,
                               Array4<Real> const& fytmp,
                               Array4<Real> const& fztmp),
                       int ncomp,
                       const Real* dx )
{
    Array<Real,AMREX_SPACEDIM> area;
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
    amrex::ParallelFor(xbx, ncomp, [fx,fxtmp,area]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fx(i,j,k,n) = fxtmp(i,j,k,n)*area[0];
    });

    //
    //  y flux
    //
    amrex::ParallelFor(ybx, ncomp, [fy,fytmp,area]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fy(i,j,k,n) = fytmp(i,j,k,n)*area[1];
    });

#if (AMREX_SPACEDIM==3)
    //
    //  z flux
    //
    amrex::ParallelFor(zbx, ncomp, [fz,fztmp,area]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fz(i,j,k,n) = fztmp(i,j,k,n)*area[2];
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
		     Array4<Real> const& scratch,
                     Array4<EBCellFlag const> const& flag,
                     Array4<Real const> const& vfrac,
                     Geometry const& geom )
{
    const Box dbox = geom.growPeriodicDomain(2);

    Box const& bxg1 = amrex::grow(bx,1);
    Box const& bxg2 = amrex::grow(bx,2);

    // Temporaries
    Array4<Real>  tmp(scratch, 0);
    Array4<Real>  delm(scratch, ncomp);
    Array4<Real>  wgt(scratch, 2*ncomp);

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
MOL::EB_AreaWeightFluxes ( D_DECL( Box const& xbx,
				   Box const& ybx,
				   Box const& zbx),
                           D_DECL( Array4<Real> const& fx,
                                   Array4<Real> const& fy,
                                   Array4<Real> const& fz),
                           D_DECL( Array4<Real> const& fxtmp,
                                   Array4<Real> const& fytmp,
                                   Array4<Real> const& fztmp),
                           D_DECL( Array4<Real const> const& apx,
                                   Array4<Real const> const& apy,
                                   Array4<Real const> const& apz ),
                           int ncomp,
                           const Real* dx )
{
    Array<Real,AMREX_SPACEDIM> area;
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
    amrex::ParallelFor(xbx, ncomp, [fx,fxtmp,area,apx]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fx(i,j,k,n) = fxtmp(i,j,k,n)*area[0]*apx(i,j,k);
    });

    //
    //  y flux
    //
    amrex::ParallelFor(ybx, ncomp, [fy,fytmp,area,apy]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fy(i,j,k,n) = fytmp(i,j,k,n)*area[1]*apy(i,j,k);
    });

#if (AMREX_SPACEDIM==3)
    //
    //  z flux
    //
    amrex::ParallelFor(zbx, ncomp, [fz,fztmp,area,apz]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fz(i,j,k,n) = fztmp(i,j,k,n)*area[2]*apz(i,j,k);
    });

#endif

}




#endif
