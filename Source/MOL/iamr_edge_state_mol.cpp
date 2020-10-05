#include <iamr_mol.H>
#include <NS_util.H>
#include <iamr_edge_state_mol_K.H>

#if AMREX_USE_EB
#include <iamr_eb_edge_state_mol_K.H>
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
			   const        BCRec * d_bcrec_ptr,
                           D_DECL( Array4<Real const> const& fcx,
                                   Array4<Real const> const& fcy,
                                   Array4<Real const> const& fcz),
                           Array4<Real const> const& ccc,
                      Array4<EBCellFlag const> const& flag)
{

    D_TERM( const Box& ubx = amrex::surroundingNodes(bx,0);,
            const Box& vbx = amrex::surroundingNodes(bx,1);,
            const Box& wbx = amrex::surroundingNodes(bx,2););

    // ****************************************************************************
    // Predict to x-faces
    // ****************************************************************************

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    auto extdir_lohi = has_extdir_or_ho(bcs.dataPtr(), ncomp, 0);
    if ((extdir_lohi.first  and domain.smallEnd(0) >= ubx.smallEnd(0)-1) or
        (extdir_lohi.second and domain.bigEnd(0)  <= ubx.bigEnd(0)))
    {
      amrex::ParallelFor(ubx, ncomp, [d_bcrec_ptr,q,ccc,AMREX_D_DECL(fcx,fcy,fcz),flag,umac, xedge, domain]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           if (flag(i,j,k).isConnected(-1,0,0))
           {
               xedge(i,j,k,n) = iamr_eb_xedge_state_mol_extdir( D_DECL(i, j, k), n, q, umac,
								D_DECL(fcx,fcy,fcz),
								ccc, flag, d_bcrec_ptr, domain );
           }
           else
           {
               xedge(i,j,k,n) = covered_val;
           }
        });
    }
    else
    {
        amrex::ParallelFor(ubx, ncomp, [q,ccc,fcx,flag,umac,xedge]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           if (flag(i,j,k).isConnected(-1,0,0))
           {
               xedge(i,j,k,n) = iamr_eb_xedge_state_mol( D_DECL(i, j, k), n, q, umac, fcx, ccc, flag );
           }
           else
           {
               xedge(i,j,k,n) = covered_val;
           }
        });
    }


    // ****************************************************************************
    // Predict to y-faces
    // ****************************************************************************
    extdir_lohi = has_extdir_or_ho(bcs.dataPtr(), ncomp, 1);
    if ((extdir_lohi.first  and domain.smallEnd(1) >= vbx.smallEnd(1)-1) or
        (extdir_lohi.second and domain.bigEnd(1)   <= vbx.bigEnd(1)))
    {
        amrex::ParallelFor(vbx, ncomp, [d_bcrec_ptr,q,ccc,AMREX_D_DECL(fcx,fcy,fcz),flag,vmac,yedge,domain]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (flag(i,j,k).isConnected(0,-1,0))
            {
                yedge(i,j,k,n) = iamr_eb_yedge_state_mol_extdir( D_DECL(i, j, k), n, q, vmac,
								 AMREX_D_DECL(fcx,fcy,fcz), ccc,
                                                                 flag, d_bcrec_ptr, domain );
            }
            else
            {
                yedge(i,j,k,n) = covered_val;
            }
        });
    }
    else
    {
        amrex::ParallelFor(vbx, ncomp, [q,ccc,fcy,flag,vmac,yedge]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (flag(i,j,k).isConnected(0,-1,0))
            {
                yedge(i,j,k,n) = iamr_eb_yedge_state_mol( D_DECL(i, j, k), n, q, vmac, fcy, ccc, flag );
            }
            else
            {
                yedge(i,j,k,n) = covered_val;
            }
        });
    }



#if ( AMREX_SPACEDIM == 3 )

    // ****************************************************************************
    // Predict to z-faces
    // ****************************************************************************
    extdir_lohi =  has_extdir_or_ho(bcs.dataPtr(), ncomp, 2);
    if ((extdir_lohi.first  and domain.smallEnd(2) >= wbx.smallEnd(2)-1) or
        (extdir_lohi.second and domain.bigEnd(2)   <= wbx.bigEnd(2)))
    {
        amrex::ParallelFor(wbx, ncomp, [d_bcrec_ptr,q,ccc,AMREX_D_DECL(fcx,fcy,fcz),flag,wmac,zedge,domain]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (flag(i,j,k).isConnected(0,0,-1))
            {
                zedge(i,j,k,n) = iamr_eb_zedge_state_mol_extdir( i, j, k, n, q, wmac,
								 AMREX_D_DECL(fcx,fcy,fcz), ccc,
                                                                 flag, d_bcrec_ptr, domain );
            }
            else
            {
                zedge(i,j,k,n) = covered_val;
            }
        });
    }
    else
    {
        amrex::ParallelFor(wbx, ncomp, [q,ccc,fcz,flag,wmac,zedge]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (flag(i,j,k).isConnected(0,0,-1))
            {
                zedge(i,j,k,n) = iamr_eb_zedge_state_mol( i, j, k, n, q, wmac, fcz, ccc, flag );
            }
            else
            {
                zedge(i,j,k,n) = covered_val;
            }
        });
    }

#endif
}

#endif //End of AMREX_USE_EB
