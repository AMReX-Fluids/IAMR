#include <Godunov.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBFArrayBox.H>
#include <NS_util.H>

using namespace amrex;

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

#ifdef AMREX_USE_EB
    AMREX_ASSERT(a_vel.hasEBFabFactory());
    auto const& fact  = dynamic_cast<EBFArrayBoxFactory const&>(a_vel.Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();
    auto const& fcent = fact.getFaceCent();
    auto const& ccent = fact.getCentroid();
#endif

    Box const& domain = a_geom.Domain();
    Gpu::DeviceVector<BCRec> bcs_device = convertToDeviceVector(a_bcs);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(a_vel, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& ubx = mfi.nodaltilebox(0);
            Box const& vbx = mfi.nodaltilebox(1);
            Box const& wbx = mfi.nodaltilebox(2);
            Array4<Real> const& u = a_umac.array(mfi);
            Array4<Real> const& v = a_vmac.array(mfi);
            Array4<Real> const& w = a_wmac.array(mfi);
            Array4<Real const> const& vcc = a_vel.const_array(mfi);
#ifdef AMREX_USE_EB
            Box const& bx = mfi.tilebox();
            EBCellFlagFab const& flagfab = flags[mfi];
            Array4<EBCellFlag const> const& flagarr = flagfab.const_array();
            auto const typ = flagfab.getType(amrex::grow(bx,1));
            if (typ == FabType::covered)
            {
                amrex::ParallelFor(ubx, vbx, wbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { u(i,j,k) = 0.0; },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { v(i,j,k) = 0.0; },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { w(i,j,k) = 0.0; });
            }
            else if (typ == FabType::singlevalued)
            {
                Array4<Real const> const& fcx = fcent[0]->const_array(mfi);
                Array4<Real const> const& fcy = fcent[1]->const_array(mfi);
                Array4<Real const> const& fcz = fcent[2]->const_array(mfi);
                Array4<Real const> const& ccc = ccent.const_array(mfi);
                predict_vels_on_faces_eb(bx,ubx,vbx,wbx,u,v,w,vcc,flagarr,fcx,fcy,fcz,ccc,a_geom,a_bcs);
            }
            else
#endif
            {
                predict_vels_on_faces(ubx,vbx,wbx,u,v,w,vcc,a_geom,a_bcs);
            }

            iamr_set_mac_bcs(domain,ubx,vbx,wbx,u,v,w,vcc,a_bcs);
        }
    }

}
