#include <iamr_godunov.H>


using namespace amrex;


void
godunov::ComputeAofs ( MultiFab& aofs, const int aofs_comp, const int ncomp,
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
                       Vector<BCRec> const& bcs,
                       Geometry const& geom,
                       Gpu::DeviceVector<int>& iconserv,
                       const Real dt,
                       const bool use_ppm,
                       const bool use_forces_in_trans,
                       const bool is_velocity  )
{
    BL_PROFILE("Godunov::ComputeAofs()");

    for (MFIter mfi(aofs,true); mfi.isValid(); ++mfi)
    {

        const Box& bx = mfi.tilebox();

        FArrayBox tmpfab(amrex::grow(bx,1), 14*ncomp);

        compute_godunov_advection( bx, ncomp,
                                   aofs.array(mfi,aofs_comp),
                                   state.array(mfi,state_comp),
                                   AMREX_D_DECL( xfluxes.array(mfi,fluxes_comp),
                                                 yfluxes.array(mfi,fluxes_comp),
                                                 zfluxes.array(mfi,fluxes_comp)),
                                   AMREX_D_DECL( xedge.array(mfi,edge_comp),
                                                 yedge.array(mfi,edge_comp),
                                                 zedge.array(mfi,edge_comp)),
                                   known_edgestate,
                                   AMREX_D_DECL(umac.array(mfi),
                                                vmac.array(mfi),
                                                wmac.array(mfi)),
                                   divu.array(mfi),
                                   fq.array(mfi,fq_comp),
                                   geom, dt, &bcs[0],
                                   iconserv.data(),
                                   tmpfab.dataPtr(),
                                   use_ppm,
                                   is_velocity,
                                   use_forces_in_trans);
    }

}
