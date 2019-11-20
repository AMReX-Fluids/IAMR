#include <Godunov.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBFArrayBox.H>

using namespace amrex;


//
// Auxiliar namespace for Riemann solver
//
namespace riemann
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real
    solver(Real umns, Real upls)
    {

        if ( umns < 0.0 && upls > 0.0 )
        {
            return 0.0;
        }
        else
        {
            Real small_vel(1.e-10);
            Real avg = 0.5 * ( upls + umns );

            if ( std::abs(avg) < small_vel )
            {
                return 0.0;
            }
            else if ( avg >= 0 )
            {
                return umns;
            }
            else
            {
                return upls;
            }

        }
    }
}


//
// Compute upwinded FC velocities by extrapolating CC values in SPACE ONLY
// This is NOT a Godunov type extrapolation: there is NO dependence on time!
// The resulting FC velocities are computed at the CENTROID of the face.
//
void
Godunov::ExtrapVelToFaces ( const MultiFab&  a_vel,
                            D_DECL( MultiFab& a_umac,
                                    MultiFab& a_vmac,
                                    MultiFab& a_wmac ),
                            D_DECL( MultiFab& a_xslopes,
                                    MultiFab& a_yslopes,
                                    MultiFab& a_zslopes),
                            const Geometry&  a_geom,
                            const amrex::Vector<amrex::BCRec>& a_bcs )
{
    BL_PROFILE("Godunov::ExtrapVelToFaces");
    AMREX_ASSERT(a_vel.hasEBFabFactory());

    Box domain(a_geom.Domain());

    const Dim3 domlo = amrex::lbound(domain);
    const Dim3 domhi = amrex::ubound(domain);

    auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(a_vel.Factory());

    //
    // Create mask
    //
    iMultiFab cc_mask(a_vel.boxArray(), a_vel.DistributionMap(), 1, 1);

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
                    AMREX_FOR_3D( b, i, j, k,
                    {
                        fab(i,j,k) = 1;
                    });
                }
            }
            // NOTE: here we do not need host-device synchronization since it
            // is already included in the MFIter destructor
        }
    }

    // Get EB geometric info
    Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
    Array< const MultiCutFab*,AMREX_SPACEDIM> facecent;

    areafrac = ebfactory.getAreaFrac();
    facecent = ebfactory.getFaceCent();

    Real  huge_vel = 1.e100;

    // ****************************************************************************
    // We will store the left and right states in arrays for interpolation
    // ****************************************************************************

    MultiFab upls(a_umac.boxArray(), a_umac.DistributionMap(), 1, 1, MFInfo(), ebfactory);
    MultiFab umns(a_umac.boxArray(), a_umac.DistributionMap(), 1, 1, MFInfo(), ebfactory);

    MultiFab vpls(a_vmac.boxArray(), a_vmac.DistributionMap(), 1, 1, MFInfo(), ebfactory);
    MultiFab vmns(a_vmac.boxArray(), a_vmac.DistributionMap(), 1, 1, MFInfo(), ebfactory);

#if (AMREX_SPACEDIM ==3)
    MultiFab wpls(a_wmac.boxArray(), a_wmac.DistributionMap(), 1, 1, MFInfo(), ebfactory);
    MultiFab wmns(a_wmac.boxArray(), a_wmac.DistributionMap(), 1, 1, MFInfo(), ebfactory);
#endif
    // We need this just to avoid FPE (eg for outflow faces)
    upls.setVal(COVERED_VAL);
    umns.setVal(COVERED_VAL);
    vpls.setVal(COVERED_VAL);
    vmns.setVal(COVERED_VAL);
#if (AMREX_SPACEDIM ==3)
    wpls.setVal(COVERED_VAL);
    wmns.setVal(COVERED_VAL);
#endif

    // Shifts to define staggered boxes
#if ( AMREX_SPACEDIM == 2 )
    IntVect e_x(1,0);
    IntVect e_y(0,1);
#elif ( AMREX_SPACEDIM == 3 )
    IntVect e_x(1,0,0);
    IntVect e_y(0,1,0);
    IntVect e_z(0,0,1);
#endif

    // Get ptr to BCRec
    const auto bc = a_bcs.dataPtr();

    // ****************************************************************************
    // Predict to face centers
    // ****************************************************************************
    for (MFIter mfi(a_vel,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Tilebox
        Box  bx = mfi.tilebox();

        Box ubx = mfi.tilebox(e_x);
        Box vbx = mfi.tilebox(e_y);
#if ( AMREX_SPACEDIM == 3 )
        Box wbx = mfi.tilebox(e_z);
#endif

        Box ubx_grown = mfi.growntilebox(e_x);
        Box vbx_grown = mfi.growntilebox(e_y);
#if ( AMREX_SPACEDIM == 3 )
        Box wbx_grown = mfi.growntilebox(e_z);
#endif

        const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>(a_vel[mfi]);
        const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

        // ===================  COVERED BRANCH =================== //
        if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
        {
             const auto& umac_fab = a_umac[mfi].array();
             const auto& vmac_fab = a_vmac[mfi].array();

             AMREX_FOR_3D(ubx, i, j, k, { umac_fab(i,j,k) = COVERED_VAL; });
             AMREX_FOR_3D(vbx, i, j, k, { vmac_fab(i,j,k) = COVERED_VAL; });

#if (AMREX_SPACEDIM == 3)
             const auto& wmac_fab = a_wmac[mfi].array();
             AMREX_FOR_3D(wbx, i, j, k, { wmac_fab(i,j,k) = COVERED_VAL; });
#endif
        }
        // ===================  REGULAR BRANCH =================== //
        else if (flags.getType(amrex::grow(bx,1)) == FabType::regular )
        {
            // Cell-centered velocity
            const auto& ccvel_fab = a_vel.array(mfi);

            // Cell-centered slopes
            D_TERM( const auto& xslopes_fab = a_xslopes.array(mfi);,
                    const auto& yslopes_fab = a_yslopes.array(mfi);,
                    const auto& zslopes_fab = a_zslopes.array(mfi););

            // Face-centered velocity components
            D_TERM( const auto& umac_fab = a_umac.array(mfi);,
                    const auto& vmac_fab = a_vmac.array(mfi);,
                    const auto& wmac_fab = a_wmac.array(mfi););

            // Face-centered left and right states
            D_TERM( const auto& upls_fab = upls.array(mfi);,
                    const auto& vpls_fab = vpls.array(mfi);,
                    const auto& wpls_fab = wpls.array(mfi););

            D_TERM( const auto& umns_fab = umns.array(mfi);,
                    const auto& vmns_fab = vmns.array(mfi);,
                    const auto& wmns_fab = wmns.array(mfi););

            // No cut cells in tile + 1-cell witdh halo -> use non-eb routine
            AMREX_FOR_3D(ubx, i, j, k,
            {
                if ( (i == domlo.x) and (bc[0].lo(0) == BCType::ext_dir) )
                {
                    umac_fab(i,j,k) = ccvel_fab(domlo.x-1,j,k,0);
                }
                else if ( (i == domhi.x+1) and (bc[0].hi(0) == BCType::ext_dir) )
                {
                    umac_fab(i,j,k) = ccvel_fab(domhi.x+1,j,k,0);
                }
                else
                {
                    upls_fab(i,j,k) = ccvel_fab(i  ,j,k,0) - 0.5 * xslopes_fab(i  ,j,k,0);
                    umns_fab(i,j,k) = ccvel_fab(i-1,j,k,0) + 0.5 * xslopes_fab(i-1,j,k,0);
                    umac_fab(i,j,k) = riemann::solver( umns_fab(i,j,k), upls_fab(i,j,k) );
                }
            });

            AMREX_FOR_3D(vbx, i, j, k,
            {
                if ( (j == domlo.y) and (bc[1].lo(1) == BCType::ext_dir) )
                {
                    vmac_fab(i,j,k) = ccvel_fab(i,domlo.y-1,k,1);
                }
                else if ( (j == domhi.y+1) and (bc[1].hi(1) == BCType::ext_dir) )
                {
                    vmac_fab(i,j,k) = ccvel_fab(i,domhi.y+1,k,1);
                }
                else
                {
                    vpls_fab(i,j,k) = ccvel_fab(i,j  ,k,1) - 0.5 * yslopes_fab(i,j  ,k,1);
                    vmns_fab(i,j,k) = ccvel_fab(i,j-1,k,1) + 0.5 * yslopes_fab(i,j-1,k,1);
                    vmac_fab(i,j,k) = riemann::solver( vmns_fab(i,j,k), vpls_fab(i,j,k) );
                }
            });
#if (AMREX_SPACEDIM == 3 )
            AMREX_FOR_3D(wbx, i, j, k,
            {
                if ( (k == domlo.z) and (bc[2].lo(2) == BCType::ext_dir) )
                {
                    wmac_fab(i,j,k) = ccvel_fab(i,j,domlo.z-1,2);
                }
                else if ( (k == domhi.z+1) and (bc[2].lo(2) == BCType::ext_dir) )
                {
                    wmac_fab(i,j,k) = ccvel_fab(i,j,domhi.z+1,2);
                }
                else
                {
                    wpls_fab(i,j,k) = ccvel_fab(i,j,k  ,2) - 0.5 * zslopes_fab(i,j,k  ,2);
                    wmns_fab(i,j,k) = ccvel_fab(i,j,k-1,2) + 0.5 * zslopes_fab(i,j,k-1,2);
                    wmac_fab(i,j,k) = riemann::solver( wmns_fab(i,j,k), wpls_fab(i,j,k) );
                }
            });
#endif
            Gpu::synchronize();
        }
        // ===================  CUT-CELLS BRANCH =================== //
        else
        {
            //
            // For cut-cells we only computed cell-centered left and right values
            // and later on we average those to the centroid and apply upwinding to them
            //

            // Cell-centered velocity
            const auto& ccvel_fab = a_vel.array(mfi);

            // Cell-centered slopes
            D_TERM( const auto& xslopes_fab = a_xslopes.array(mfi);,
                    const auto& yslopes_fab = a_yslopes.array(mfi);,
                    const auto& zslopes_fab = a_zslopes.array(mfi););

            // Face-centered left and right states
            D_TERM( const auto& upls_fab = upls.array(mfi);,
                    const auto& vpls_fab = vpls.array(mfi);,
                    const auto& wpls_fab = wpls.array(mfi););

            D_TERM( const auto& umns_fab = umns.array(mfi);,
                    const auto& vmns_fab = vmns.array(mfi);,
                    const auto& wmns_fab = wmns.array(mfi););

            // Face-centered areas
            D_TERM( const auto& apx_fab = areafrac[0]->array(mfi);,
                    const auto& apy_fab = areafrac[1]->array(mfi);,
                    const auto& apz_fab = areafrac[2]->array(mfi););

            // If near a Dirichlet boundary, define both upls and umns to be
            // the Dirichlet value so that we do not need to modify the the interpolation
            // to face centroid that follows.
            AMREX_FOR_3D(ubx_grown, i, j, k,
            {
                if (apx_fab(i,j,k) > 0.0)
                {
                    if ( (i <= domlo.x) and (bc[0].lo(0) == BCType::ext_dir) )
                    {
                        upls_fab(i,j,k) = ccvel_fab(domlo.x-1,j,k,0);
                        umns_fab(i,j,k) = ccvel_fab(domlo.x-1,j,k,0);
                    }
                    else if ( (i >= domhi.x+1) and (bc[0].hi(0) == BCType::ext_dir) )
                    {
                        upls_fab(i,j,k) = ccvel_fab(domhi.x+1,j,k,0);
                        umns_fab(i,j,k) = ccvel_fab(domhi.x+1,j,k,0);
                    }
                    else
                    {
                        upls_fab(i,j,k) = ccvel_fab(i  ,j,k,0) - 0.5 * xslopes_fab(i  ,j,k,0);
                        umns_fab(i,j,k) = ccvel_fab(i-1,j,k,0) + 0.5 * xslopes_fab(i-1,j,k,0);
                    }
                }
            });

            AMREX_FOR_3D(vbx_grown, i, j, k,
            {
                if (apy_fab(i,j,k) > 0.0)
                {
                    if ( (j <= domlo.y) and (bc[1].lo(1) == BCType::ext_dir) )
                    {
                        vpls_fab(i,j,k) = ccvel_fab(i,domlo.y-1,k,1);
                        vmns_fab(i,j,k) = ccvel_fab(i,domlo.y-1,k,1);
                    }
                    else if ( (j >= domhi.y+1) and (bc[1].hi(1) == BCType::ext_dir) )
                    {
                        vpls_fab(i,j,k) = ccvel_fab(i,domhi.y+1,k,1);
                        vmns_fab(i,j,k) = ccvel_fab(i,domhi.y+1,k,1);
                    }
                    else
                    {
                        vpls_fab(i,j,k) = ccvel_fab(i,j  ,k,1) - 0.5 * yslopes_fab(i,j  ,k,1);
                        vmns_fab(i,j,k) = ccvel_fab(i,j-1,k,1) + 0.5 * yslopes_fab(i,j-1,k,1);
                    }
                }
            });

#if (AMREX_SPACEDIM == 3 )
            AMREX_FOR_3D(wbx_grown, i, j, k,
            {
                if (apz_fab(i,j,k) > 0.0)
                {
                    if ( (k <= domlo.z) and (bc[2].lo(2) == BCType::ext_dir) )
                    {
                        wpls_fab(i,j,k) = ccvel_fab(i,j,domlo.z-1,2);
                        wmns_fab(i,j,k) = ccvel_fab(i,j,domlo.z-1,2);
                    }
                    else if ( (k >= domhi.z+1) and (bc[2].lo(2) == BCType::ext_dir) )
                    {
                        wpls_fab(i,j,k) = ccvel_fab(i,j,domhi.z+1,2);
                        wmns_fab(i,j,k) = ccvel_fab(i,j,domhi.z+1,2);
                    }
                    else
                    {
                        wpls_fab(i,j,k) = ccvel_fab(i,j,k  ,2) - 0.5 * zslopes_fab(i,j,k  ,2);
                        wmns_fab(i,j,k) = ccvel_fab(i,j,k-1,2) + 0.5 * zslopes_fab(i,j,k-1,2);
                    }
                }
            });
#endif
            Gpu::synchronize();

        } // Cut cells

    } // MFIter

    // ****************************************************************************
    // Make sure to fill face-centered values outside our grid before interpolating
    // ****************************************************************************
    upls.FillBoundary(a_geom.periodicity());
    umns.FillBoundary(a_geom.periodicity());
    vpls.FillBoundary(a_geom.periodicity());
    vmns.FillBoundary(a_geom.periodicity());
#if (AMREX_SPACEDIM == 3 )
    wpls.FillBoundary(a_geom.periodicity());
    wmns.FillBoundary(a_geom.periodicity());
#endif

    // ****************************************************************************
    // Do interpolation to centroids -- only for cut cells
    // ****************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(a_vel,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Tilebox
        Box  bx = mfi.tilebox();
        Box ubx = mfi.tilebox(e_x);
        Box vbx = mfi.tilebox(e_y);
#if ( AMREX_SPACEDIM == 3 )
        Box wbx = mfi.tilebox(e_z);
#endif

        // Check efficiently if this tile contains any eb stuff
        const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>(a_vel[mfi]);
        const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

        if ( !(flags.getType(amrex::grow(bx,0)) == FabType::covered ||
               flags.getType(amrex::grow(bx,1)) == FabType::regular ) )
        {
            // Face-centered velocity components
            D_TERM( const auto& umac_fab = a_umac.array(mfi);,
                    const auto& vmac_fab = a_vmac.array(mfi);,
                    const auto& wmac_fab = a_wmac.array(mfi););

            // Face-centered left and right states
            D_TERM( const auto& upls_fab = upls.array(mfi);,
                    const auto& vpls_fab = vpls.array(mfi);,
                    const auto& wpls_fab = wpls.array(mfi););

            D_TERM( const auto& umns_fab = umns.array(mfi);,
                    const auto& vmns_fab = vmns.array(mfi);,
                    const auto& wmns_fab = wmns.array(mfi););

            // Face-centered areas
            D_TERM( const auto& apx_fab = areafrac[0]->array(mfi);,
                    const auto& apy_fab = areafrac[1]->array(mfi);,
                    const auto& apz_fab = areafrac[2]->array(mfi););

            // Face centroids
            D_TERM( const auto& fcx_fab = facecent[0]->array(mfi);,
                    const auto& fcy_fab = facecent[1]->array(mfi);,
                    const auto& fcz_fab = facecent[2]->array(mfi););

            const auto& ccm_fab = cc_mask.const_array(mfi);

            AMREX_FOR_3D(ubx, i, j, k,
            {
                if (apx_fab(i,j,k) == 0.0)
                {
                    umac_fab(i,j,k) = huge_vel;
                }
                else if (apx_fab(i,j,k) < 1.0)
                {
                    int jj = j + static_cast<int>(std::copysign(1.0, fcx_fab(i,j,k,0)));
#if (AMREX_SPACEDIM == 3)
                    int kk = k + static_cast<int>(std::copysign(1.0, fcx_fab(i,j,k,1)));
#else
                    int kk = k;
#endif
                    Real fracy = (ccm_fab(i-1,jj,k) || ccm_fab(i,jj,k)) ? std::abs(fcx_fab(i,j,k,0)) : 0.0;
#if (AMREX_SPACEDIM == 3)
                    Real fracz = (ccm_fab(i-1,j,kk) || ccm_fab(i,j,kk)) ? std::abs(fcx_fab(i,j,k,1)) : 0.0;
#else
                    Real fracz = 0.0;
#endif
                    Real upls_on_centroid = (1.0-fracy)*(1.0-fracz)*upls_fab(i, j,k )+
                        fracy *(1.0-fracz)*upls_fab(i,jj,k )+
                        fracz *(1.0-fracy)*upls_fab(i, j,kk)+
                        fracy *     fracz *upls_fab(i,jj,kk);
                    Real umns_on_centroid = (1.0-fracy)*(1.0-fracz)*umns_fab(i, j,k )+
                        fracy *(1.0-fracz)*umns_fab(i,jj,k )+
                        fracz *(1.0-fracy)*umns_fab(i, j,kk)+
                        fracy *     fracz *umns_fab(i,jj,kk);

                    umac_fab(i,j,k) = riemann::solver(umns_on_centroid, upls_on_centroid);
                }
                else
                {
                    umac_fab(i,j,k) = riemann::solver(umns_fab(i,j,k), upls_fab(i,j,k));
                }
            });

            AMREX_FOR_3D(vbx, i, j, k,
            {
                if (apy_fab(i,j,k) == 0.0)
                {
                    vmac_fab(i,j,k) = huge_vel;
                }
                else if (apy_fab(i,j,k) < 1.0)
                {
                    int ii = i + static_cast<int>(std::copysign(1.0,fcy_fab(i,j,k,0)));
#if (AMREX_SPACEDIM == 3)
                    int kk = k + static_cast<int>(std::copysign(1.0,fcy_fab(i,j,k,1)));
#else
                    int kk = k;
#endif
                    Real fracx = (ccm_fab(ii,j-1,k ) || ccm_fab(ii,j,k)) ? std::abs(fcy_fab(i,j,k,0)) : 0.0;
#if (AMREX_SPACEDIM == 3)
                    Real fracz = (ccm_fab( i,j-1,kk) || ccm_fab(i,j,kk)) ? std::abs(fcy_fab(i,j,k,1)) : 0.0;
#else
                    Real fracz = 0.0;
#endif

                    Real vpls_on_centroid = (1.0-fracx)*(1.0-fracz)*vpls_fab(i ,j,k )+
                        fracx *(1.0-fracz)*vpls_fab(ii,j,k )+
                        fracz *(1.0-fracx)*vpls_fab(i ,j,kk)+
                        fracx *     fracz *vpls_fab(ii,j,kk);
                    Real vmns_on_centroid = (1.0-fracx)*(1.0-fracz)*vmns_fab(i ,j,k )+
                        fracx *(1.0-fracz)*vmns_fab(ii,j,k )+
                        fracz *(1.0-fracx)*vmns_fab(i ,j,kk)+
                        fracx *     fracz *vmns_fab(ii,j,kk);

                    vmac_fab(i,j,k) = riemann::solver(vmns_on_centroid, vpls_on_centroid);
                }
                else
                {
                    vmac_fab(i,j,k) = riemann::solver(vmns_fab(i,j,k), vpls_fab(i,j,k));
                }
            });

#if (AMREX_SPACEDIM == 3)
             AMREX_FOR_3D(wbx, i, j, k,
             {
                if (apz_fab(i,j,k) == 0.0)

                    wmac_fab(i,j,k) = huge_vel;

                else if (apz_fab(i,j,k) < 1.0) {

                    int ii = i + static_cast<int>(std::copysign(1.0,fcz_fab(i,j,k,0)));
                    int jj = j + static_cast<int>(std::copysign(1.0,fcz_fab(i,j,k,1)));

                    Real fracx = (ccm_fab(ii,j,k-1) || ccm_fab(ii,j,k)) ? std::abs(fcz_fab(i,j,k,0)) : 0.0;
                    Real fracy = (ccm_fab(i,jj,k-1) || ccm_fab(i,jj,k)) ? std::abs(fcz_fab(i,j,k,1)) : 0.0;

                    Real wpls_on_centroid = (1.0-fracx)*(1.0-fracy)*wpls_fab(i ,j ,k)+
                                                 fracx *(1.0-fracy)*wpls_fab(ii,j ,k)+
                                                 fracy *(1.0-fracx)*wpls_fab(i ,jj,k)+
                                                 fracx *     fracy *wpls_fab(ii,jj,k);
                    Real wmns_on_centroid = (1.0-fracx)*(1.0-fracy)*wmns_fab(i ,j ,k)+
                                                 fracx *(1.0-fracy)*wmns_fab(ii,j ,k)+
                                                 fracy *(1.0-fracx)*wmns_fab(i ,jj,k)+
                                                 fracx *     fracy *wmns_fab(ii,jj,k);

                    wmac_fab(i,j,k) = riemann::solver(wmns_on_centroid, wpls_on_centroid);
                }
                else
                {
                    wmac_fab(i,j,k) = riemann::solver(wmns_fab(i,j,k), wpls_fab(i,j,k));
                }

             });
#endif
             Gpu::synchronize();

          } // Cut cells
       } // MFIter
}
