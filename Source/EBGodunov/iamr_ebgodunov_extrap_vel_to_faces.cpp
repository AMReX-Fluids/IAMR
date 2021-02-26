#include <iamr_ebgodunov_plm.H>
#include <iamr_godunov_plm.H>
#include <iamr_ebgodunov.H>
#include <iamr_godunov.H>
#include <iamr_godunov_K.H>
#include <iamr_ebgodunov_transverse_2D_K.H>
//#include <iomanip>

using namespace amrex;
void
EBGodunov::ExtrapVelToFaces ( MultiFab const& vel,
                              MultiFab const& vel_forces,
                              AMREX_D_DECL(MultiFab& u_mac,
                                           MultiFab& v_mac,
                                           MultiFab& w_mac),
                              Vector<BCRec> const& h_bcrec,
                              BCRec  const* d_bcrec,
                              Geometry& geom,
                              Real l_dt)
{
    BL_PROFILE("EBGodunov::ExtrapVelToFaces()");
    AMREX_ALWAYS_ASSERT(vel.hasEBFabFactory());

    Box const& domain = geom.Domain();
    const Real* dx    = geom.CellSize();

    auto const& ebfact= dynamic_cast<EBFArrayBoxFactory const&>(vel.Factory());
    auto const& flags = ebfact.getMultiEBCellFlagFab();
    auto const& fcent = ebfact.getFaceCent();
    auto const& ccent = ebfact.getCentroid();
    auto const& vfrac = ebfact.getVolFrac();
    auto const& areafrac = ebfact.getAreaFrac();

    // Since we don't fill all the ghost cells in the mac vel arrays
    // we need to initialize to something which won't make the code crash
    AMREX_D_TERM( u_mac.setVal(1.e40);,
                  v_mac.setVal(1.e40);,
                  w_mac.setVal(1.e40););

    const int ncomp = AMREX_SPACEDIM;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox scratch;
        for (MFIter mfi(vel,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Box const& bxg2 = amrex::grow(bx,2);

            EBCellFlagFab const& flagfab = flags[mfi];
            Array4<EBCellFlag const> const& flagarr = flagfab.const_array();

            AMREX_D_TERM( Array4<Real> const& a_umac = u_mac.array(mfi);,
                          Array4<Real> const& a_vmac = v_mac.array(mfi);,
                          Array4<Real> const& a_wmac = w_mac.array(mfi););

            Array4<Real const> const& a_vel = vel.const_array(mfi);
            Array4<Real const> const& a_f = vel_forces.const_array(mfi);

            // Not all the arrays have exactly the size of bxg3 but none are bigger
            // In 2-d:
            //  8*ncomp are:  Imx, Ipx, Imy, Ipy, xlo/xhi, ylo/yhi
            //  2       are:  u_ad, v_ad
            // In 3-d:
            // 12*ncomp are:  Imx, Ipx, Imy, Ipy, Imz, Ipz, xlo/xhi, ylo/yhi, zlo/zhi
            //  3       are:  u_ad, v_ad, w_ad
            Box const& bxg3 = amrex::grow(bx,3);
            scratch.resize(bxg3, 4*AMREX_SPACEDIM*ncomp + AMREX_SPACEDIM);
            Real* p  = scratch.dataPtr();

            AMREX_D_TERM(Box const& xbx = mfi.nodaltilebox(0);,
                         Box const& ybx = mfi.nodaltilebox(1);,
                         Box const& zbx = mfi.nodaltilebox(2));;

            AMREX_D_TERM(Box xebx(Box(bx).grow(1).surroundingNodes(0));,
                         Box yebx(Box(bx).grow(1).surroundingNodes(1));,
                         Box zebx(Box(bx).grow(1).surroundingNodes(2)));

#if (AMREX_SPACEDIM == 2)
            Box xebx_g2(Box(bx).grow(1).grow(1,1).surroundingNodes(0));
            Box yebx_g2(Box(bx).grow(1).grow(0,1).surroundingNodes(1));
#else
            Box xebx_g2(Box(bx).grow(1).grow(1,1).grow(2,1).surroundingNodes(0));
            Box yebx_g2(Box(bx).grow(1).grow(0,1).grow(2,1).surroundingNodes(1));
            Box zebx_g2(Box(bx).grow(1).grow(0,1).grow(1,1).surroundingNodes(2));
#endif

            Array4<Real> Imx = makeArray4(p,bxg2,ncomp);
            p +=         Imx.size();
            Array4<Real> Ipx = makeArray4(p,bxg2,ncomp);
            p +=         Ipx.size();
            Array4<Real> Imy = makeArray4(p,bxg2,ncomp);
            p +=         Imy.size();
            Array4<Real> Ipy = makeArray4(p,bxg2,ncomp);
            p +=         Ipy.size();

            Array4<Real> u_ad = makeArray4(p,xebx_g2,1);
            p +=         u_ad.size();
            Array4<Real> v_ad = makeArray4(p,yebx_g2,1);
            p +=         v_ad.size();

#if (AMREX_SPACEDIM == 3)
            Array4<Real> Imz = makeArray4(p,bxg2,ncomp);
            p +=         Imz.size();
            Array4<Real> Ipz = makeArray4(p,bxg2,ncomp);
            p +=         Ipz.size();

            Array4<Real> w_ad = makeArray4(p,zebx_g2,1);
            p +=         w_ad.size();
#endif

            // This tests on covered cells just in the box itself
            if (flagfab.getType(bx) == FabType::covered)
            {
                // We shouldn't need to zero these

            // This tests on only regular cells including two rows of ghost cells
            }
            else if (flagfab.getType(amrex::grow(bx,2)) == FabType::regular)
            {

                PLM::PredictVelOnXFace( xebx_g2, AMREX_SPACEDIM, Imx, Ipx, a_vel, a_vel,
                                         geom, l_dt, h_bcrec, d_bcrec);

                PLM::PredictVelOnYFace( yebx_g2, AMREX_SPACEDIM, Imy, Ipy, a_vel, a_vel,
                                        geom, l_dt, h_bcrec, d_bcrec);

#if ( AMREX_SPACEDIM == 3 )
                PLM::PredictVelOnZFace( zebx_g2, AMREX_SPACEDIM, Imz, Ipz, a_vel, a_vel,
                                        geom, l_dt, h_bcrec, d_bcrec);
#endif


                bool local_use_forces_in_trans = false;
                Godunov::ComputeAdvectiveVel( AMREX_D_DECL(Box(u_ad), Box(v_ad), Box(w_ad)),
                                              AMREX_D_DECL(u_ad, v_ad, w_ad),
                                              AMREX_D_DECL(Imx, Imy, Imz),
                                              AMREX_D_DECL(Ipx, Ipy, Ipz),
                                              a_vel, a_f, domain, l_dt, d_bcrec,
                                              local_use_forces_in_trans);

                Godunov::ExtrapVelToFacesOnBox( bx, ncomp,
                                                AMREX_D_DECL(xbx, ybx, zbx),
                                                AMREX_D_DECL(a_umac, a_vmac, a_wmac),
                                                a_vel,
                                                AMREX_D_DECL(u_ad, v_ad, w_ad),
                                                AMREX_D_DECL(Imx, Imy, Imz),
                                                AMREX_D_DECL(Ipx, Ipy, Ipz),
                                                a_f, domain, dx, l_dt, d_bcrec,
                                                local_use_forces_in_trans, p);
            }
            else
            {

                AMREX_D_TERM(Array4<Real const> const& fcx = fcent[0]->const_array(mfi);,
                             Array4<Real const> const& fcy = fcent[1]->const_array(mfi);,
                             Array4<Real const> const& fcz = fcent[2]->const_array(mfi););

                Array4<Real const> const& ccent_arr = ccent.const_array(mfi);
                Array4<Real const> const& vfrac_arr = vfrac.const_array(mfi);

                EBPLM::PredictVelOnXFace( xebx_g2, Imx, Ipx, a_vel, a_vel,
                                          flagarr, vfrac_arr,
                                          AMREX_D_DECL(fcx,fcy,fcz),ccent_arr,
                                          geom, l_dt, h_bcrec, d_bcrec );

                EBPLM::PredictVelOnYFace( yebx_g2, Imy, Ipy, a_vel, a_vel,
                                          flagarr, vfrac_arr,
                                          AMREX_D_DECL(fcx,fcy,fcz),ccent_arr,
                                          geom, l_dt, h_bcrec, d_bcrec );

#if (AMREX_SPACEDIM == 3)
                EBPLM::PredictVelOnZFace( zebx_g2, Imz, Ipz, a_vel, a_vel,
                                          flagarr, vfrac_arr,
                                          AMREX_D_DECL(fcx,fcy,fcz),ccent_arr,
                                          geom, l_dt, h_bcrec, d_bcrec );
#endif

                EBGodunov::ComputeAdvectiveVel( AMREX_D_DECL(xebx_g2, yebx_g2, zebx_g2),
                                                AMREX_D_DECL(u_ad, v_ad, w_ad),
                                                AMREX_D_DECL(Imx, Imy, Imz),
                                                AMREX_D_DECL(Ipx, Ipy, Ipz),
                                                a_vel, flagarr, domain, d_bcrec);

                AMREX_D_TERM(Array4<Real const> const& apx = areafrac[0]->const_array(mfi);,
                             Array4<Real const> const& apy = areafrac[1]->const_array(mfi);,
                             Array4<Real const> const& apz = areafrac[2]->const_array(mfi););

                EBGodunov::ExtrapVelToFacesOnBox( bx, ncomp,
                                                  AMREX_D_DECL(xbx,ybx,zbx),
                                                  AMREX_D_DECL(xebx_g2,yebx_g2,zebx_g2),
                                                  AMREX_D_DECL(a_umac, a_vmac, a_wmac),
                                                  a_vel,
                                                  AMREX_D_DECL(u_ad, v_ad, w_ad),
                                                  AMREX_D_DECL(Imx, Imy, Imz),
                                                  AMREX_D_DECL(Ipx, Ipy, Ipz),
                                                  a_f,
                                                  domain, dx, l_dt, d_bcrec,
                                                  flagarr,
                                                  AMREX_D_DECL(apx, apy, apz),
#if (AMREX_SPACEDIM == 3)
                                                  vfrac_arr,
#endif
                                                  AMREX_D_DECL(fcx, fcy, fcz),
                                                  p);
            }

            Gpu::streamSynchronize();  // otherwise we might be using too much memory
        }
    }
}

void
EBGodunov::ComputeAdvectiveVel ( AMREX_D_DECL(Box const& xbx,
                                              Box const& ybx,
                                              Box const& zbx),
                                 AMREX_D_DECL(Array4<Real> const& u_ad,
                                              Array4<Real> const& v_ad,
                                              Array4<Real> const& w_ad),
                                 AMREX_D_DECL(Array4<Real const> const& Imx,
                                              Array4<Real const> const& Imy,
                                              Array4<Real const> const& Imz),
                                 AMREX_D_DECL(Array4<Real const> const& Ipx,
                                              Array4<Real const> const& Ipy,
                                              Array4<Real const> const& Ipz),
                                 Array4<Real const> const& vel,
                                 Array4<EBCellFlag const> const& flag,
                                 const Box& domain,
                                 BCRec  const* pbc)
  {
    const Dim3 dlo = amrex::lbound(domain);
    const Dim3 dhi = amrex::ubound(domain);

    amrex::ParallelFor(AMREX_D_DECL(xbx, ybx, zbx),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // We only care about x-velocity on x-faces here
        if (flag(i,j,k).isConnected(-1,0,0))
        {
            constexpr int n = 0;

            Real lo = Ipx(i-1,j,k,n);
            Real hi = Imx(i  ,j,k,n);

            auto bc = pbc[n];
            SetTransTermXBCs(i, j, k, n, vel, lo, hi, bc.lo(0), bc.hi(0), dlo.x, dhi.x, true);

            Real st = ( (lo+hi) >= 0.) ? lo : hi;
            bool ltm = ( (lo <= 0. && hi >= 0.) || (amrex::Math::abs(lo+hi) < small_vel) );
            u_ad(i,j,k) = ltm ? 0. : st;
        } else {
            u_ad(i,j,k) = 0.;
        }
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // We only care about y-velocity on y-faces here
        if (flag(i,j,k).isConnected(0,-1,0))
        {
            constexpr int n = 1;

            Real lo = Ipy(i,j-1,k,n);
            Real hi = Imy(i,j  ,k,n);

            auto bc = pbc[n];
            SetTransTermYBCs(i, j, k, n, vel, lo, hi, bc.lo(1), bc.hi(1), dlo.y, dhi.y, true);

            Real st = ( (lo+hi) >= 0.) ? lo : hi;
            bool ltm = ( (lo <= 0. && hi >= 0.) || (amrex::Math::abs(lo+hi) < small_vel) );
            v_ad(i,j,k) = ltm ? 0. : st;
        } else {
            v_ad(i,j,k) = 0.;
        }
#if (AMREX_SPACEDIM == 3)
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (flag(i,j,k).isConnected(0,0,-1))
        {
            // We only care about z-velocity on z-faces here
            constexpr int n = 2;

            Real lo = Ipz(i,j,k-1,n);
            Real hi = Imz(i,j,k  ,n);

            auto bc = pbc[n];
            SetTransTermZBCs(i, j, k, n, vel, lo, hi, bc.lo(2), bc.hi(2), dlo.z, dhi.z, true);

            Real st = ( (lo+hi) >= 0.) ? lo : hi;
            bool ltm = ( (lo <= 0. && hi >= 0.) || (amrex::Math::abs(lo+hi) < small_vel) );
            w_ad(i,j,k) = ltm ? 0. : st;
        } else {
            w_ad(i,j,k) = 0.;
        }
#endif
    });
}
