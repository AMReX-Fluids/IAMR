#include <iamr_mol.H>
#include <NS_util.H>
#include <iamr_mol_edge_state_K.H>

#ifdef AMREX_USE_EB
#include <iamr_mol_eb_edge_state_K.H>
#endif


using namespace amrex;

namespace {
    std::pair<bool,bool> has_extdir_or_ho (BCRec const* bcrec, int ncomp, int dir)
    {
        std::pair<bool,bool> r{false,false};
        for (int n = 0; n < ncomp; ++n)
        {
            r.first = r.first or bcrec[n].lo(dir) == BCType::ext_dir
                              or bcrec[n].lo(dir) == BCType::hoextrap;
            r.second = r.second or bcrec[n].hi(dir) == BCType::ext_dir
                                or bcrec[n].hi(dir) == BCType::hoextrap;
        }
        return r;
    }
}



//
// Compute edge state on REGULAR box
//
void
MOL::ComputeEdgeState (const Box& bx,
                       D_DECL( Array4<Real> const& xedge,
                               Array4<Real> const& yedge,
                               Array4<Real> const& zedge),
                       Array4<Real const> const& q,
                       const int ncomp,
                       D_DECL( Array4<Real const> const& umac,
                               Array4<Real const> const& vmac,
                               Array4<Real const> const& wmac),
                       const Box&       domain,
                       const Vector<BCRec>& bcs,
                       const        BCRec * d_bcrec_ptr)
{
    const int domain_ilo = domain.smallEnd(0);
    const int domain_ihi = domain.bigEnd(0);
    const int domain_jlo = domain.smallEnd(1);
    const int domain_jhi = domain.bigEnd(1);
#if (AMREX_SPACEDIM==3)
    const int domain_klo = domain.smallEnd(2);
    const int domain_khi = domain.bigEnd(2);
#endif

    D_TERM( const Box& ubx = amrex::surroundingNodes(bx,0);,
            const Box& vbx = amrex::surroundingNodes(bx,1);,
            const Box& wbx = amrex::surroundingNodes(bx,2););

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    auto extdir_lohi = has_extdir_or_ho(bcs.dataPtr(), ncomp, 0);
    bool has_extdir_or_ho_lo = extdir_lohi.first;
    bool has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo and domain_ilo >= ubx.smallEnd(0)-1) or
        (has_extdir_or_ho_hi and domain_ihi <= ubx.bigEnd(0)))
    {
        amrex::ParallelFor(ubx, ncomp, [d_bcrec_ptr,q,domain_ilo,domain_ihi,umac,xedge]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            xedge(i,j,k,n) = iamr_xedge_state_mol_extdir( i, j, k, n, q, umac,
                                                          d_bcrec_ptr,
                                                          domain_ilo, domain_ihi );
        });
    }
    else
    {
        amrex::ParallelFor(ubx, ncomp, [q,umac,xedge]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            xedge(i,j,k,n) = iamr_xedge_state_mol( i, j, k, n, q, umac);
        });
    }


    extdir_lohi = has_extdir_or_ho(bcs.dataPtr(), ncomp, 1);
    has_extdir_or_ho_lo = extdir_lohi.first;
    has_extdir_or_ho_hi = extdir_lohi.second;
    if ((has_extdir_or_ho_lo and domain_jlo >= vbx.smallEnd(1)-1) or
        (has_extdir_or_ho_hi and domain_jhi <= vbx.bigEnd(1)))
    {
        amrex::ParallelFor(vbx, ncomp, [d_bcrec_ptr,q,domain_jlo,domain_jhi,vmac,yedge]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            yedge(i,j,k,n) = iamr_yedge_state_mol_extdir( i, j, k, n, q, vmac,
                                                          d_bcrec_ptr,
                                                          domain_jlo, domain_jhi );
        });
    }
    else
    {
        amrex::ParallelFor(vbx, ncomp, [q,vmac,yedge]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            yedge(i,j,k,n) = iamr_yedge_state_mol( i, j, k, n, q, vmac );
        });
    }


#if ( AMREX_SPACEDIM ==3 )

    extdir_lohi = has_extdir_or_ho(bcs.dataPtr(), ncomp, 2);
    has_extdir_or_ho_lo = extdir_lohi.first;
    has_extdir_or_ho_hi = extdir_lohi.second;
    if ((has_extdir_or_ho_lo and domain_klo >= wbx.smallEnd(2)-1) or
        (has_extdir_or_ho_hi and domain_khi <= wbx.bigEnd(2)))
    {
        amrex::ParallelFor(wbx, ncomp, [d_bcrec_ptr,q,domain_klo,domain_khi,wmac,zedge]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            zedge(i,j,k,n) = iamr_zedge_state_mol_extdir( i, j, k, n, q, wmac,
                                                          d_bcrec_ptr,
                                                          domain_klo, domain_khi );
        });
    }
    else
    {
        amrex::ParallelFor(wbx, ncomp, [q,wmac,zedge]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            zedge(i,j,k,n) = iamr_zedge_state_mol( i, j, k, n, q, wmac );
        });
    }

#endif
}


#ifdef AMREX_USE_EB
//
// Compute edge state on EB box
//
void
MOL::EB_ComputeEdgeState ( Box const& bx,
                           AMREX_D_DECL( Array4<Real> const& xedge,
                                         Array4<Real> const& yedge,
                                         Array4<Real> const& zedge),
                           Array4<Real const> const& q,
                           const int ncomp,
                           AMREX_D_DECL( Array4<Real const> const& umac,
                                         Array4<Real const> const& vmac,
                                         Array4<Real const> const& wmac),
                           const Box&       domain,
                           const Vector<BCRec>& bcs,
			   const        BCRec * d_bcrec_ptr,
                           AMREX_D_DECL( Array4<Real const> const& fcx,
                                         Array4<Real const> const& fcy,
                                         Array4<Real const> const& fcz),
                           Array4<Real const> const& ccc,
                           Array4<Real const> const& vfrac,
                           Array4<EBCellFlag const> const& flag)
{

    int order = 2;

    AMREX_D_TERM( const Box& ubx = amrex::surroundingNodes(bx,0);,
                  const Box& vbx = amrex::surroundingNodes(bx,1);,
                  const Box& wbx = amrex::surroundingNodes(bx,2););

    AMREX_D_TERM(
        const int domain_ilo = domain.smallEnd(0);
        const int domain_ihi = domain.bigEnd(0);,
        const int domain_jlo = domain.smallEnd(1);
        const int domain_jhi = domain.bigEnd(1);,
        const int domain_klo = domain.smallEnd(2);
        const int domain_khi = domain.bigEnd(2););

    // ****************************************************************************
    // Decide whether the stencil at each cell might need to see values that
    //     live on face centroids rather than cell centroids, i.e.
    //     are at a domain boundary with ext_dir or hoextrap boundary conditions
    // ****************************************************************************

    auto extdir_lohi_x = has_extdir_or_ho(bcs.dataPtr(), ncomp, static_cast<int>(Direction::x));
    bool has_extdir_or_ho_lo_x = extdir_lohi_x.first;
    bool has_extdir_or_ho_hi_x = extdir_lohi_x.second;

    auto extdir_lohi_y = has_extdir_or_ho(bcs.dataPtr(), ncomp, static_cast<int>(Direction::y));
    bool has_extdir_or_ho_lo_y = extdir_lohi_y.first;
    bool has_extdir_or_ho_hi_y = extdir_lohi_y.second;

#if (AMREX_SPACEDIM == 3)
    auto extdir_lohi_z = has_extdir_or_ho(bcs.dataPtr(), ncomp, static_cast<int>(Direction::z));
    bool has_extdir_or_ho_lo_z = extdir_lohi_z.first;
    bool has_extdir_or_ho_hi_z = extdir_lohi_z.second;
#endif


    if ((has_extdir_or_ho_lo_x and domain_ilo >= ubx.smallEnd(0)-1) or
        (has_extdir_or_ho_hi_x and domain_ihi <= ubx.bigEnd(0)    ) or
        (has_extdir_or_ho_lo_y and domain_jlo >= vbx.smallEnd(1)-1) or
        (has_extdir_or_ho_hi_y and domain_jhi <= vbx.bigEnd(1)    )
#if (AMREX_SPACEDIM == 2)
        )
#elif (AMREX_SPACEDIM == 3)
        or
        (has_extdir_or_ho_lo_z and domain_klo >= wbx.smallEnd(2)-1) or
        (has_extdir_or_ho_hi_z and domain_khi <= wbx.bigEnd(2)    ) )
#endif
    {

        // ****************************************************************************
        // Predict to x-faces
        // ****************************************************************************
        amrex::ParallelFor(ubx, ncomp, [d_bcrec_ptr,q,ccc,AMREX_D_DECL(fcx,fcy,fcz),flag,umac, xedge, domain, vfrac, order]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           if (flag(i,j,k).isConnected(-1,0,0))
           {
               xedge(i,j,k,n) = iamr_eb_xedge_state_mol_extdir( AMREX_D_DECL(i, j, k), n, q, umac,
								AMREX_D_DECL(fcx,fcy,fcz), ccc, vfrac,
                                                                flag, d_bcrec_ptr, domain, order );
           }
           else
           {
               xedge(i,j,k,n) = 0.0;
           }
        });

        // ****************************************************************************
        // Predict to y-faces
        // ****************************************************************************
        amrex::ParallelFor(vbx, ncomp, [d_bcrec_ptr,q,ccc,AMREX_D_DECL(fcx,fcy,fcz),flag,vmac,yedge,domain,vfrac,order]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (flag(i,j,k).isConnected(0,-1,0))
            {
                yedge(i,j,k,n) = iamr_eb_yedge_state_mol_extdir( AMREX_D_DECL(i, j, k), n, q, vmac,
								 AMREX_D_DECL(fcx,fcy,fcz), ccc, vfrac,
                                                                 flag, d_bcrec_ptr, domain, order );
            }
            else
            {
                yedge(i,j,k,n) = 0.0;
            }
        });

#if ( AMREX_SPACEDIM == 3 )
        // ****************************************************************************
        // Predict to z-faces
        // ****************************************************************************
        amrex::ParallelFor(wbx, ncomp, [d_bcrec_ptr,q,ccc,AMREX_D_DECL(fcx,fcy,fcz),flag,wmac,zedge,domain,vfrac,order]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (flag(i,j,k).isConnected(0,0,-1))
            {
                zedge(i,j,k,n) = iamr_eb_zedge_state_mol_extdir( i, j, k, n, q, wmac,
								 AMREX_D_DECL(fcx,fcy,fcz), ccc, vfrac,
                                                                 flag, d_bcrec_ptr, domain, order );
            }
            else
            {
                zedge(i,j,k,n) = 0.0;
            }
        });
#endif
    }
    else // We assume below that the stencil does not need to use hoextrap or extdir boundaries
    {
        // ****************************************************************************
        // Predict to x-faces
        // ****************************************************************************
        amrex::ParallelFor(ubx, ncomp, [q,ccc,AMREX_D_DECL(fcx,fcy,fcz),flag, umac, xedge, vfrac, order]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           if (flag(i,j,k).isConnected(-1,0,0))
           {
               xedge(i,j,k,n) = iamr_eb_xedge_state_mol( AMREX_D_DECL(i, j, k), n, q, umac,
                                                         AMREX_D_DECL(fcx,fcy,fcz), ccc, vfrac,
                                                         flag, order);
           }
           else
           {
               xedge(i,j,k,n) = 0.0;
           }
        });

        // ****************************************************************************
        // Predict to y-faces
        // ****************************************************************************
        amrex::ParallelFor(vbx, ncomp, [q,ccc,AMREX_D_DECL(fcx,fcy,fcz),flag, vmac, yedge, vfrac, order]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (flag(i,j,k).isConnected(0,-1,0))
            {
                yedge(i,j,k,n) = iamr_eb_yedge_state_mol( AMREX_D_DECL(i, j, k), n, q, vmac,
                                                          AMREX_D_DECL(fcx,fcy,fcz), ccc, vfrac,
                                                          flag, order );
            }
            else
            {
                yedge(i,j,k,n) = 0.0;
            }
        });

#if ( AMREX_SPACEDIM == 3 )
        // ****************************************************************************
        // Predict to z-faces
        // ****************************************************************************
        amrex::ParallelFor(wbx, ncomp, [q,ccc,AMREX_D_DECL(fcx,fcy,fcz),flag, wmac, zedge, vfrac, order]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (flag(i,j,k).isConnected(0,0,-1))
            {
                zedge(i,j,k,n) = iamr_eb_zedge_state_mol( AMREX_D_DECL(i, j, k), n, q, wmac,
                                                          AMREX_D_DECL(fcx,fcy,fcz), ccc, vfrac,
                                                          flag, order );
            }
            else
            {
                zedge(i,j,k,n) = 0.0;
            }
        });
#endif
    }

}

#endif //End of AMREX_USE_EB
