#include <NS_util.H>
#include <iamr_mol.H>
#ifdef AMREX_USE_EB
#include <AMReX_MultiCutFab.H>
#endif

using namespace amrex;

//
// Compute upwinded FC velocities by extrapolating CC values in SPACE ONLY
// This is NOT a Godunov type extrapolation: there is NO dependence on time!
// The resulting FC velocities are computed at the CENTROID of the face.
//
void
MOL::ExtrapVelToFaces ( const MultiFab&  a_vel,
                        D_DECL( MultiFab& a_umac,
                                MultiFab& a_vmac,
                                MultiFab& a_wmac ),
                        const Geometry&  a_geom,
                        const amrex::Vector<amrex::BCRec>& a_bcs )
{
    BL_PROFILE("Godunov::ExtrapVelToFaces");

#ifdef AMREX_USE_EB
    AMREX_ASSERT(a_vel.hasEBFabFactory());
    auto const& fact  = dynamic_cast<EBFArrayBoxFactory const&>(a_vel.Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();
    auto const& fcent = fact.getFaceCent();
    auto const& ccent = fact.getCentroid();
#endif

    Box const& domain = a_geom.Domain();
    Gpu::DeviceVector<BCRec> bcs_device = convertToDeviceVector(a_bcs);
    const BCRec* bcs_ptr = bcs_device.dataPtr();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(a_vel, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            D_TERM( Box const& ubx = mfi.nodaltilebox(0);,
                    Box const& vbx = mfi.nodaltilebox(1);,
                    Box const& wbx = mfi.nodaltilebox(2););

            D_TERM( Array4<Real> const& u = a_umac.array(mfi);,
                    Array4<Real> const& v = a_vmac.array(mfi);,
                    Array4<Real> const& w = a_wmac.array(mfi););

            Array4<Real const> const& vcc = a_vel.const_array(mfi);

#ifdef AMREX_USE_EB
            Box const& bx = mfi.tilebox();
            EBCellFlagFab const& flagfab = flags[mfi];
            Array4<EBCellFlag const> const& flagarr = flagfab.const_array();
            auto const typ = flagfab.getType(amrex::grow(bx,2));
            if (typ == FabType::covered)
            {
                amrex::ParallelFor(ubx, [u] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { u(i,j,k) = 0.0; });
                amrex::ParallelFor(vbx, [v] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { v(i,j,k) = 0.0; });
#if (AMREX_SPACEDIM==3)
                amrex::ParallelFor(wbx, [w] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { w(i,j,k) = 0.0; });
#endif
            }
            else if (typ == FabType::singlevalued)
            {
                D_TERM( Array4<Real const> const& fcx = fcent[0]->const_array(mfi);,
                        Array4<Real const> const& fcy = fcent[1]->const_array(mfi);,
                        Array4<Real const> const& fcz = fcent[2]->const_array(mfi););

                Array4<Real const> const& ccc = ccent.const_array(mfi);

                MOL::EB_PredictVelOnFaces(bx,D_DECL(ubx,vbx,wbx),D_DECL(u,v,w),vcc,flagarr,
                                          D_DECL(fcx,fcy,fcz),ccc,a_geom, bcs_ptr);
            }
            else
#endif
            {
                MOL::PredictVelOnFaces(D_DECL(ubx,vbx,wbx),D_DECL(u,v,w),vcc,a_geom,bcs_ptr);
            }

            MOL::SetMacBCs(domain,D_DECL(ubx,vbx,wbx),D_DECL(u,v,w),vcc,a_bcs);
        }
    }

}

void
MOL::SetMacBCs ( Box const& a_domain,
                 D_DECL( Box const& a_ubx,
                         Box const& a_vbx,
                         Box const& a_wbx ) ,
                 D_DECL( Array4<Real> const& a_u,
                         Array4<Real> const& a_v,
                         Array4<Real> const& a_w ),
                 Array4<Real const> const& a_vel,
                 Vector<BCRec> const& a_bcs )
{
    int idim = 0;
    if (a_bcs[idim].lo(idim) == BCType::ext_dir and
        a_domain.smallEnd(idim) == a_ubx.smallEnd(idim))
    {
        amrex::ParallelFor(amrex::bdryLo(a_ubx,idim),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a_u(i,j,k) = a_vel(i-1,j,k,0);
        });
    }

    if (a_bcs[idim].hi(idim) == BCType::ext_dir and
        a_domain.bigEnd(idim)+1 == a_ubx.bigEnd(idim))
    {
        amrex::ParallelFor(amrex::bdryHi(a_ubx,idim),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a_u(i,j,k) = a_vel(i,j,k,0);
        });
    }

    idim = 1;
    if (a_bcs[idim].lo(idim) == BCType::ext_dir and
        a_domain.smallEnd(idim) == a_vbx.smallEnd(idim))
    {
        amrex::ParallelFor(amrex::bdryLo(a_vbx,idim),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a_v(i,j,k) = a_vel(i,j-1,k,1);
        });
    }

    if (a_bcs[idim].hi(idim) == BCType::ext_dir and
        a_domain.bigEnd(idim)+1 == a_vbx.bigEnd(idim))
    {
        amrex::ParallelFor(amrex::bdryHi(a_vbx,idim),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a_v(i,j,k) = a_vel(i,j,k,1);
        });
    }

#if (AMREX_SPACEDIM==3)
    idim = 2;
    if (a_bcs[idim].lo(idim) == BCType::ext_dir and
        a_domain.smallEnd(idim) == a_wbx.smallEnd(idim))
    {
        amrex::ParallelFor(amrex::bdryLo(a_wbx,idim),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a_w(i,j,k) = a_vel(i,j,k-1,2);
        });
    }

    if (a_bcs[idim].hi(idim) == BCType::ext_dir and
        a_domain.bigEnd(idim)+1 == a_wbx.bigEnd(idim))
    {
        amrex::ParallelFor(amrex::bdryHi(a_wbx,idim),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a_w(i,j,k) = a_vel(i,j,k,2);
        });
    }
#endif
}
