#include <Godunov.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <AMReX_BCRec.H>
#include <iamr_edge_state_mol_K.H>
#if AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <iamr_eb_edge_state_mol_K.H>
#endif
#include <iamr_mol.H>
#include <NS_util.H>

using namespace amrex;

namespace {
    std::pair<bool,bool> has_extdir (BCRec const* bcrec, int ncomp, int dir)
    {
        std::pair<bool,bool> r{false,false};
        for (int n = 0; n < ncomp; ++n)
        {
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
                    const int a_ncomp,
                    D_DECL( Array4<Real const> const& a_umac,
                            Array4<Real const> const& a_vmac,
                            Array4<Real const> const& a_wmac),
                    const Box&       a_domain,
                    const Vector<BCRec>& a_bcs,
                    int known_edgestate)
{
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

    // MOL::ComputeEdgeState(a_bx, a_edgeq_x, a_edgeq_y, a_edgeq_z, a_q, a_ncomp,
    //                       a_umac, a_vmac, a_wmac, a_domain, a_bcs );

    const auto bc = a_bcs.dataPtr();
    auto d_bcrec  = convertToDeviceVector(a_bcs);

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    auto extdir_lohi = has_extdir(a_bcs.dataPtr(), a_ncomp, 0);
    bool has_extdir_lo = extdir_lohi.first;
    bool has_extdir_hi = extdir_lohi.second;

    amrex::ParallelFor(ubx, a_ncomp, [a_q,a_umac,a_fx,known_edgestate,a_edgeq_x]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        a_fx(i,j,k,n) = a_edgeq_x(i,j,k,n) * a_umac(i,j,k);
    });

    amrex::ParallelFor(vbx, a_ncomp, [a_q,a_vmac,a_fy,known_edgestate,a_edgeq_y]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        a_fy(i,j,k,n) = a_edgeq_y(i,j,k,n) * a_vmac(i,j,k);
    });


#if ( AMREX_SPACEDIM ==3 )
        amrex::ParallelFor(wbx, a_ncomp, [a_q,a_wmac,a_fz, a_edgeq_z, known_edgestate]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            a_fz(i,j,k,n) = a_edgeq_z(i,j,k,n) * a_wmac(i,j,k);
        });
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
    amrex::ParallelFor(ubx, a_ncomp,
    [a_q,a_ccc,a_fcx,a_flag,a_umac,a_fx,known_edgestate, a_edgeq_x]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (a_flag(i,j,k).isConnected(-1,0,0))
        {
            a_fx(i,j,k,n) = a_umac(i,j,k) * a_edgeq_x(i,j,k,n);
        }
        else
        {
            a_fx(i,j,k,n) = 0.0;
        }
    });

    // ****************************************************************************
    // Predict to y-faces
    // ****************************************************************************
    extdir_lohi = has_extdir(a_bcs.dataPtr(), a_ncomp, 1);
    amrex::ParallelFor(vbx, a_ncomp,
    [a_q,a_ccc,a_fcy,a_flag,a_vmac,a_fy,known_edgestate,a_edgeq_y]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Real qs;
        if (a_flag(i,j,k).isConnected(0,-1,0))
        {
            a_fy(i,j,k,n) = a_vmac(i,j,k) * a_edgeq_y(i,j,k,n);
        }
        else
        {
            a_fy(i,j,k,n) = 0.0;
        }
    });


    //
    // ===================== Z =====================
    //
#if ( AMREX_SPACEDIM == 3 )


    // ****************************************************************************
    // Predict to z-faces
    // ****************************************************************************

    extdir_lohi =  has_extdir(a_bcs.dataPtr(), a_ncomp, 2);
    amrex::ParallelFor(wbx, a_ncomp,
    [a_q,a_ccc,a_fcz,a_flag,a_wmac,a_fz,known_edgestate,a_edgeq_z]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (a_flag(i,j,k).isConnected(0,0,-1))
        {
            a_fz(i,j,k,n) = a_wmac(i,j,k) * a_edgeq_z(i,j,k,n);
        }
        else
        {
            a_fz(i,j,k,n) = 0.0;
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
                MOL::ComputeEdgeState(bx, D_DECL(edgestate_x.array(mfi), edgestate_y.array(mfi), edgestate_z.array(mfi)),
                                      a_state.array(mfi,a_comp), a_ncomp,
                                      D_DECL(a_umac.array(mfi), a_vmac.array(mfi), a_wmac.array(mfi)),
                                      domain, a_bcs );

                MOL::ComputeFluxes( bx, D_DECL(a_fx.array(mfi), a_fy.array(mfi), a_fz.array(mfi)),
                                    D_DECL(a_umac.array(mfi), a_vmac.array(mfi), a_wmac.array(mfi)),
                                    D_DECL(edgestate_x.array(mfi), edgestate_y.array(mfi), edgestate_z.array(mfi)), a_ncomp );
            }
            else
            {
                Box tmpbox  = amrex::surroundingNodes(bx);
                Array4<Real const> ccc = ebfactory.getCentroid().const_array(mfi);

                MOL::EB_ComputeEdgeState(tmpbox,
                                         D_DECL(edgestate_x.array(mfi), edgestate_y.array(mfi), edgestate_z.array(mfi)),
                                         a_state.array(mfi,a_comp), a_ncomp,
                                         D_DECL(a_umac.array(mfi), a_vmac.array(mfi), a_wmac.array(mfi)),
                                         domain, a_bcs,
                                         D_DECL(facecent[0]->array(mfi), facecent[1]->array(mfi), facecent[2]->array(mfi)),
                                         ccc, flag );

                MOL::EB_ComputeFluxes( bx, D_DECL(a_fx.array(mfi), a_fy.array(mfi), a_fz.array(mfi)),
                                       D_DECL(a_umac.array(mfi), a_vmac.array(mfi), a_wmac.array(mfi)),
                                       D_DECL(edgestate_x.array(mfi), edgestate_y.array(mfi), edgestate_z.array(mfi)),
                                       a_ncomp, flag );
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
                                                a_state.array(mfi,a_comp), a_ncomp,
                                                D_DECL(a_umac.array(mfi), a_vmac.array(mfi), a_wmac.array(mfi)), domain, a_bcs, 0);

                }
                else
                {
                    Array4<Real const> ccc = ebfactory.getCentroid().const_array(mfi);
                    fluxes::ComputeFluxesOnEBBox(bx, D_DECL(a_fx.array(mfi), a_fy.array(mfi), a_fz.array(mfi)),
                                                 D_DECL(edgestate_x.array(mfi), edgestate_y.array(mfi), edgestate_z.array(mfi)),
                                                 a_state.array(mfi,a_comp), a_ncomp,
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
                                            a_state.array(mfi,a_comp), a_ncomp,
                                            D_DECL(a_ucorr.array(mfi), a_vcorr.array(mfi), a_wcorr.array(mfi)), domain, a_bcs, 1);

            }
            else
            {
                Array4<Real const> ccc = ebfactory.getCentroid().const_array(mfi);
                fluxes::ComputeFluxesOnEBBox(bx, D_DECL(a_fx.array(mfi), a_fy.array(mfi), a_fz.array(mfi)),
                                             D_DECL(edgestate_x.array(mfi), edgestate_y.array(mfi), edgestate_z.array(mfi)),
                                             a_state.array(mfi,a_comp), a_ncomp,
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
