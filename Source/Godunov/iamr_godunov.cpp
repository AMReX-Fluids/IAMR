#include <iamr_godunov.H>
#include <NS_util.H>

using namespace amrex;


void
Godunov::ComputeAofs ( MultiFab& aofs, const int aofs_comp, const int ncomp,
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
                       BCRec const* d_bc,
                       Geometry const& geom,
                       Gpu::DeviceVector<int>& iconserv,
                       const Real dt,
                       const bool use_ppm,
                       const bool use_forces_in_trans,
                       const bool is_velocity  )
{
    BL_PROFILE("Godunov::ComputeAofs()");

#if (AMREX_SPACEDIM==2)
    MultiFab* volume;
    MultiFab* area[AMREX_SPACEDIM];

    if ( geom.IsRZ() )
    {
      const DistributionMapping& dmap = aofs.DistributionMap();
      const BoxArray& grids = aofs.boxArray();
      const int ngrow_vol = aofs.nGrow();

      volume = new MultiFab(grids,dmap,1,ngrow_vol);
      geom.GetVolume(*volume);

      const int ngrow_area = xfluxes.nGrow();

      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
      {
	BoxArray edge_ba(grids);
	area[dir] = new MultiFab(edge_ba.surroundingNodes(dir),dmap,1,ngrow_area);
	geom.GetFaceArea(*area[dir],dir);
      }
    }
#endif

    //FIXME - check on adding tiling here
    for (MFIter mfi(aofs); mfi.isValid(); ++mfi)
    {

        const Box& bx   = mfi.tilebox();

        //
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

        if (not known_edgestate)
        {
            ComputeEdgeState( bx, ncomp,
                              state.array(mfi,state_comp),
                              AMREX_D_DECL( xed, yed, zed ),
                              AMREX_D_DECL( u, v, w ),
                              divu.array(mfi),
                              fq.array(mfi,fq_comp),
                              geom, dt, d_bc,
                              iconserv.data(),
                              use_ppm,
                              use_forces_in_trans,
                              is_velocity );
        }

#if (AMREX_SPACEDIM == 2)
	if ( geom.IsRZ() )
	{
	  const auto& areax = area[0]->array(mfi);
	  const auto& areay = area[1]->array(mfi);
	  const auto& vol   = volume->array(mfi);

	  ComputeFluxes_rz( bx, AMREX_D_DECL( fx, fy, fz ),
			    AMREX_D_DECL( u, v, w ),
			    AMREX_D_DECL( xed, yed, zed ),
			    areax, areay,
			    ncomp );


	  ComputeDivergence_rz( bx,
				aofs.array(mfi,aofs_comp),
				AMREX_D_DECL( fx, fy, fz ),
				AMREX_D_DECL( xed, yed, zed ),
				AMREX_D_DECL( u, v, w ),
				areax, areay, vol,
				ncomp, iconserv.data() );

	}
	else
#endif
	{
	  ComputeFluxes( bx, AMREX_D_DECL( fx, fy, fz ),
			 AMREX_D_DECL( u, v, w ),
			 AMREX_D_DECL( xed, yed, zed ),
			 geom, ncomp );

	  ComputeDivergence( bx,
			     aofs.array(mfi,aofs_comp),
			     AMREX_D_DECL( fx, fy, fz ),
			     AMREX_D_DECL( xed, yed, zed ),
			     AMREX_D_DECL( u, v, w ),
			     ncomp, geom, iconserv.data() );
	}

	//
	// NOTE this sync cannot protect temporaries in ComputeEdgeState, ComputeFluxes
	// or ComputeDivergence, since functions have their own scope. As soon as the
	// CPU hits the end of the function, it will call the destructor for all
	// temporaries created in that function.
	//
        Gpu::streamSynchronize();  // otherwise we might be using too much memory
    }

}



void
Godunov::ComputeSyncAofs ( MultiFab& aofs, const int aofs_comp, const int ncomp,
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
                           BCRec const* d_bc,
                           Geometry const& geom,
                           Gpu::DeviceVector<int>& iconserv,
                           const Real dt,
                           const bool use_ppm,
                           const bool use_forces_in_trans,
                           const bool is_velocity  )
{
    BL_PROFILE("Godunov::ComputeAofs()");

#if (AMREX_SPACEDIM==2)
    MultiFab* volume;
    MultiFab* area[AMREX_SPACEDIM];

    if ( geom.IsRZ() )
    {
      const DistributionMapping& dmap = aofs.DistributionMap();
      const BoxArray& grids = aofs.boxArray();
      const int ngrow_vol = aofs.nGrow();

      volume = new MultiFab(grids,dmap,1,ngrow_vol);
      geom.GetVolume(*volume);

      const int ngrow_area = xfluxes.nGrow();

      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
      {
	BoxArray edge_ba(grids);
	area[dir] = new MultiFab(edge_ba.surroundingNodes(dir),dmap,1,ngrow_area);
	geom.GetFaceArea(*area[dir],dir);
      }
    }
#endif

    //FIXME - check on adding tiling here
    for (MFIter mfi(aofs); mfi.isValid(); ++mfi)
    {

        const Box& bx   = mfi.tilebox();

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

        if (not known_edgestate)
        {

            AMREX_D_TERM( const auto& u = umac.const_array(mfi);,
                          const auto& v = vmac.const_array(mfi);,
                          const auto& w = wmac.const_array(mfi););

            ComputeEdgeState( bx, ncomp,
                              state.array(mfi,state_comp),
                              AMREX_D_DECL( xed, yed, zed ),
                              AMREX_D_DECL( u, v, w ),
                              divu.array(mfi),
                              fq.array(mfi,fq_comp),
                              geom, dt, d_bc,
                              iconserv.data(),
                              use_ppm,
                              use_forces_in_trans,
                              is_velocity );
        }

#if (AMREX_SPACEDIM == 2)
	if ( geom.IsRZ() )
	{
	  const auto& areax = area[0]->array(mfi);
	  const auto& areay = area[1]->array(mfi);
	  const auto& vol   = volume->array(mfi);

	  ComputeFluxes_rz( bx, AMREX_D_DECL( fx, fy, fz ),
			    AMREX_D_DECL( uc, vc, wc ),
			    AMREX_D_DECL( xed, yed, zed ),
			    areax, areay,
			    ncomp );


	  ComputeSyncDivergence_rz( bx,
				    aofs.array(mfi,aofs_comp),
				    AMREX_D_DECL( fx, fy, fz ),
				    vol,
				    ncomp );

	}
	else
#endif
	{
	  ComputeFluxes( bx, AMREX_D_DECL( fx, fy, fz ),
			 AMREX_D_DECL( uc, vc, wc ),
			 AMREX_D_DECL( xed, yed, zed ),
			 geom, ncomp );


	  ComputeSyncDivergence( bx,
				 aofs.array(mfi,aofs_comp),
				 AMREX_D_DECL( fx, fy, fz ),
				 ncomp, geom );
	}

        Gpu::streamSynchronize();  // otherwise we might be using too much memory
    }

}


void
Godunov::ComputeFluxes ( Box const& bx,
                         AMREX_D_DECL( Array4<Real> const& fx,
                                       Array4<Real> const& fy,
                                       Array4<Real> const& fz),
                         AMREX_D_DECL( Array4<Real const> const& umac,
                                       Array4<Real const> const& vmac,
                                       Array4<Real const> const& wmac),
                         AMREX_D_DECL( Array4<Real const> const& xed,
                                       Array4<Real const> const& yed,
                                       Array4<Real const> const& zed),
                         Geometry const& geom, const int ncomp )
{

    const auto dx = geom.CellSizeArray();

    GpuArray<Real,AMREX_SPACEDIM> area;
#if ( AMREX_SPACEDIM == 3 )
    area[0] = dx[1]*dx[2];
    area[1] = dx[0]*dx[2];
    area[2] = dx[0]*dx[1];
#else
    area[0] = dx[0];
    area[1] = dx[1];
#endif

    //
    //  X flux
    //
    const Box& xbx = amrex::surroundingNodes(bx,0);

    amrex::ParallelFor(xbx, ncomp, [fx, umac, xed, area]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fx(i,j,k,n) = xed(i,j,k,n) * umac(i,j,k) * area[0];
    });

    //
    //  y flux
    //
    const Box& ybx = amrex::surroundingNodes(bx,1);

    amrex::ParallelFor(ybx, ncomp, [fy, vmac, yed, area]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fy(i,j,k,n) = yed(i,j,k,n) * vmac(i,j,k) * area[1];
    });

#if (AMREX_SPACEDIM==3)
    //
    //  z flux
    //
    const Box& zbx = amrex::surroundingNodes(bx,2);

    amrex::ParallelFor(zbx, ncomp, [fz, wmac, zed, area]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fz(i,j,k,n) = zed(i,j,k,n) * wmac(i,j,k) * area[2];
    });

#endif

}



void
Godunov::ComputeDivergence ( Box const& bx,
                             Array4<Real> const& div,
                             AMREX_D_DECL( Array4<Real const> const& fx,
                                           Array4<Real const> const& fy,
                                           Array4<Real const> const& fz),
                             AMREX_D_DECL( Array4<Real const> const& xed,
                                           Array4<Real const> const& yed,
                                           Array4<Real const> const& zed),
                             AMREX_D_DECL( Array4<Real const> const& umac,
                                           Array4<Real const> const& vmac,
                                           Array4<Real const> const& wmac),
                             const int ncomp, Geometry const& geom,
                             int const* iconserv )
{
    const auto dxinv = geom.InvCellSizeArray();

#if (AMREX_SPACEDIM==3)
    Real qvol = dxinv[0] * dxinv[1] * dxinv[2];
#else
    Real qvol = dxinv[0] * dxinv[1];
#endif

    amrex::ParallelFor(bx, ncomp,[=]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (iconserv[n])
        {
            div(i,j,k,n) =  qvol *
                (
                         fx(i+1,j,k,n) -  fx(i,j,k,n)
                       + fy(i,j+1,k,n) -  fy(i,j,k,n)
#if (AMREX_SPACEDIM==3)
                       + fz(i,j,k+1,n) -  fz(i,j,k,n)
#endif
                );
	    if ( i ==1 && j==25){
	      std::cout<<"fx "<<fx(i+1,j,k,n)
		       <<" "<<fx(i,j,k,n)
		       <<"\nfy "<<fy(i,j+1,k,n)
		       <<" "<<fy(i,j,k,n)<<std::endl;
	    }
        }
        else
        {
	    div(i,j,k,n) = 0.5*dxinv[0]*( umac(i+1,j,k  ) +  umac(i,j,k  ))
                *                       (  xed(i+1,j,k,n) -   xed(i,j,k,n))
                +          0.5*dxinv[1]*( vmac(i,j+1,k  ) +  vmac(i,j,k  ))
                *                       (  yed(i,j+1,k,n) -   yed(i,j,k,n))
#if (AMREX_SPACEDIM==3)
                +          0.5*dxinv[2]*( wmac(i,j,k+1  ) +  wmac(i,j,k  ))
                *                       (  zed(i,j,k+1,n) -   zed(i,j,k,n))
#endif
                ;
       }

    });
}


void
Godunov::ComputeSyncDivergence ( Box const& bx,
                                 Array4<Real> const& div,
                                 AMREX_D_DECL( Array4<Real const> const& fx,
                                               Array4<Real const> const& fy,
                                               Array4<Real const> const& fz),
                                 const int ncomp, Geometry const& geom )
{

    const auto dxinv = geom.InvCellSizeArray();

#if (AMREX_SPACEDIM==3)
    Real qvol = dxinv[0] * dxinv[1] * dxinv[2];
#else
    Real qvol = dxinv[0] * dxinv[1];
#endif

    amrex::ParallelFor(bx, ncomp,[=]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        div(i,j,k,n) +=  qvol * (
              fx(i+1,j,k,n) -  fx(i,j,k,n)
            + fy(i,j+1,k,n) -  fy(i,j,k,n)
#if (AMREX_SPACEDIM==3)
            + fz(i,j,k+1,n) -  fz(i,j,k,n)
#endif
            );
    });
}


void
Godunov::ComputeFluxes_rz ( Box const& bx,
			    Array4<Real> const& fx,
			    Array4<Real> const& fy,
			    Array4<Real const> const& umac,
			    Array4<Real const> const& vmac,
			    Array4<Real const> const& xed,
			    Array4<Real const> const& yed,
			    Array4<Real const> const& areax,
			    Array4<Real const> const& areay,
			    const int ncomp )
{
    //
    //  X flux
    //
    const Box& xbx = amrex::surroundingNodes(bx,0);

    amrex::ParallelFor(xbx, ncomp, [fx, umac, xed, areax]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      fx(i,j,k,n) = xed(i,j,k,n) * umac(i,j,k) * areax(i,j,k);
    });

    //
    //  y flux
    //
    const Box& ybx = amrex::surroundingNodes(bx,1);

    amrex::ParallelFor(ybx, ncomp, [fy, vmac, yed, areay]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fy(i,j,k,n) = yed(i,j,k,n) * vmac(i,j,k) * areay(i,j,k);
    });
}



void
Godunov::ComputeDivergence_rz ( Box const& bx,
				Array4<Real> const& div,
				Array4<Real const> const& fx,
				Array4<Real const> const& fy,
				Array4<Real const> const& xed,
				Array4<Real const> const& yed,
				Array4<Real const> const& umac,
				Array4<Real const> const& vmac,
				Array4<Real const> const& areax,
				Array4<Real const> const& areay,
				Array4<Real const> const& vol,
				const int ncomp,
				int const* iconserv )
{

    amrex::ParallelFor(bx, ncomp,[=]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (iconserv[n])
        {
	    div(i,j,k,n) = ( fx(i+1,j,k,n) -  fx(i,j,k,n) +
			     fy(i,j+1,k,n) -  fy(i,j,k,n)
			   ) / vol(i,j,k) ;
        }
        else
        {
	    //
	    // Compute  (U dot grad)c as ( -c(div U) + div(cU) )
	    //
	    const Real divux = ( areax(i+1,j,k)*umac(i+1,j,k) -
				 areax(i,  j,k)*umac(i,  j,k)
			       ) / vol(i,j,k);
	    const Real divuy = ( areay(i,j+1,k)*vmac(i,j+1,k) -
				 areay(i,j  ,k)*vmac(i,j  ,k)
			       ) / vol(i,j,k);
	    div(i,j,k,n) =  ( fx(i+1,j,k,n) - fx(i,j,k,n) +
			      fy(i,j+1,k,n) - fy(i,j,k,n) ) / vol(i,j,k)
	                  - ( divux*0.5*(xed(i+1,j,k,n) + xed(i,j,k,n)) +
			      divuy*0.5*(yed(i,j+1,k,n) + yed(i,j,k,n)) );
	}
    });
}

void
Godunov::ComputeSyncDivergence_rz ( Box const& bx,
				    Array4<Real> const& div,
				    Array4<Real const> const& fx,
				    Array4<Real const> const& fy,
				    Array4<Real const> const& vol,
				    const int ncomp )
{
    amrex::ParallelFor(bx, ncomp,[=]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      div(i,j,k,n) +=  ( fx(i+1,j,k,n) -  fx(i,j,k,n) +
			 fy(i,j+1,k,n) -  fy(i,j,k,n)
		       ) / vol(i,j,k);
    });
}
