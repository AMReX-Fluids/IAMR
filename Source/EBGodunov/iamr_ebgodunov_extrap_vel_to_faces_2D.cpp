//#include <NS_util.H>
#include <iamr_ebgodunov_plm.H>
#include <iamr_godunov_plm.H>
// #include <iamr_godunov.H>
#include <iamr_ebgodunov.H>
#include <iamr_godunov.H>
#include <iamr_godunov_K.H>
#include <iamr_ebgodunov_transverse_2D_K.H>
//#include <iomanip>

// using namespace amrex;
// void
// EBGodunov::ExtrapVelToFaces ( MultiFab const& vel,
//                               MultiFab const& vel_forces,
//                               MultiFab& u_mac, MultiFab& v_mac,
//                               Vector<BCRec> const& h_bcrec,
//                               BCRec  const* d_bcrec,
//                               Geometry& geom,
//                               Real l_dt)
// {
//     Box const& domain = geom.Domain();
//     const Real* dx    = geom.CellSize();

//     auto const& ebfact= dynamic_cast<EBFArrayBoxFactory const&>(vel.Factory());
//     auto const& flags = ebfact.getMultiEBCellFlagFab();
//     auto const& fcent = ebfact.getFaceCent();
//     auto const& ccent = ebfact.getCentroid();
//     auto const& vfrac = ebfact.getVolFrac();
//     auto const& areafrac = ebfact.getAreaFrac();

//     const int ncomp = AMREX_SPACEDIM;
// #ifdef _OPENMP
// #pragma omp parallel if (Gpu::notInLaunchRegion())
// #endif
//     {
//         FArrayBox scratch;
//         for (MFIter mfi(vel,TilingIfNotGPU()); mfi.isValid(); ++mfi)
//         {
//             Box const& bx = mfi.tilebox();
//             Box const& bxg1 = amrex::grow(bx,1);
//             Box const& xbx = mfi.nodaltilebox(0);
//             Box const& ybx = mfi.nodaltilebox(1);

//             EBCellFlagFab const& flagfab = flags[mfi];
//             Array4<EBCellFlag const> const& flagarr = flagfab.const_array();

//             Array4<Real> const& a_umac = u_mac.array(mfi);
//             Array4<Real> const& a_vmac = v_mac.array(mfi);

//             Array4<Real const> const& a_vel = vel.const_array(mfi);
//             Array4<Real const> const& a_f = vel_forces.const_array(mfi);

//             scratch.resize(bxg1, ncomp*(4*AMREX_SPACEDIM)+AMREX_SPACEDIM);
//             Real* p = scratch.dataPtr();

//             Array4<Real> Imx = makeArray4(p,bxg1,ncomp);
//             p +=         Imx.size();
//             Array4<Real> Ipx = makeArray4(p,bxg1,ncomp);
//             p +=         Ipx.size();
//             Array4<Real> Imy = makeArray4(p,bxg1,ncomp);
//             p +=         Imy.size();
//             Array4<Real> Ipy = makeArray4(p,bxg1,ncomp);
//             p +=         Ipy.size();
//             Array4<Real> u_ad = makeArray4(p,Box(bx).grow(1,1).surroundingNodes(0),1);
//             p +=         u_ad.size();
//             Array4<Real> v_ad = makeArray4(p,Box(bx).grow(0,1).surroundingNodes(1),1);
//             p +=         v_ad.size();

//             // This tests on covered cells just in the box itself
//             if (flagfab.getType(bx) == FabType::covered)
//             {
//                 // We shouldn't need to zero these

//             // This tests on only regular cells including two rows of ghost cells
//             }
//             else if (flagfab.getType(amrex::grow(bx,2)) == FabType::regular)
//             {

//                 PLM::PredictVelOnXFace( bx, AMREX_SPACEDIM, Imx, Ipx, a_vel, a_vel,
//                                          geom, l_dt, h_bcrec, d_bcrec);

//                 PLM::PredictVelOnYFace( bx, AMREX_SPACEDIM, Imy, Ipy, a_vel, a_vel,
//                                         geom, l_dt, h_bcrec, d_bcrec);

//                 bool local_use_forces_in_trans = false;
//                 Godunov::ComputeAdvectiveVel( Box(u_ad), Box(v_ad),
//                                               u_ad, v_ad,
//                                               Imx, Imy, Ipx, Ipy,
//                                               a_vel, a_f, domain, l_dt, d_bcrec,
//                                               local_use_forces_in_trans);

//                 Godunov::ExtrapVelToFacesOnBox( bx, ncomp, xbx, ybx,
//                                                 a_umac, a_vmac, a_vel,
//                                                 u_ad, v_ad,
//                                                 Imx, Imy, Ipx, Ipy,
//                                                 a_f, domain, dx, l_dt, d_bcrec,
//                                                 local_use_forces_in_trans, p);

//             }
//             else
//             {

//                 AMREX_D_TERM(Array4<Real const> const& fcx = fcent[0]->const_array(mfi);,
//                              Array4<Real const> const& fcy = fcent[1]->const_array(mfi);,
//                              Array4<Real const> const& fcz = fcent[2]->const_array(mfi););

//                 Array4<Real const> const& ccent_arr = ccent.const_array(mfi);
//                 Array4<Real const> const& vfrac_arr = vfrac.const_array(mfi);

//                 EBPLM::PredictVelOnXFace( bx, Imx, Ipx, a_vel, a_vel,
//                                           flagarr, vfrac_arr,
//                                           AMREX_D_DECL(fcx,fcy,fcz),ccent_arr,
//                                           geom, l_dt, h_bcrec, d_bcrec );

//                 EBPLM::PredictVelOnYFace( bx, Imy, Ipy, a_vel, a_vel,
//                                           flagarr, vfrac_arr,
//                                           AMREX_D_DECL(fcx,fcy,fcz),ccent_arr,
//                                           geom, l_dt, h_bcrec, d_bcrec );

//                 EBGodunov::ComputeAdvectiveVel(Box(u_ad), Box(v_ad),
//                                                  u_ad, v_ad,
//                                                  Imx, Imy, Ipx, Ipy, a_vel,
//                                                  flagarr, domain, d_bcrec);

//                 AMREX_D_TERM(Array4<Real const> const& apx = areafrac[0]->const_array(mfi);,
//                              Array4<Real const> const& apy = areafrac[1]->const_array(mfi);,
//                              Array4<Real const> const& apz = areafrac[2]->const_array(mfi););

//                 EBGodunov::ExtrapVelToFacesOnBox( bx, ncomp, xbx, ybx,
//                                                   a_umac, a_vmac,
//                                                   a_vel, u_ad, v_ad,
//                                                   Imx, Imy, Ipx, Ipy, a_f,
//                                                   domain, dx, l_dt, d_bcrec,
//                                                   flagarr, apx,apy,vfrac_arr,
//                                                   fcx,fcy,p);
//             }

//             Gpu::streamSynchronize();  // otherwise we might be using too much memory
//         }
//     }
// }

// void
// EBGodunov::ComputeAdvectiveVel ( Box const& xbx, Box const& ybx,
//                                  Array4<Real> const& u_ad,
//                                  Array4<Real> const& v_ad,
//                                  Array4<Real const> const& Imx,
//                                  Array4<Real const> const& Imy,
//                                  Array4<Real const> const& Ipx,
//                                  Array4<Real const> const& Ipy,
//                                  Array4<Real const> const& vel,
//                                  Array4<EBCellFlag const> const& flag,
//                                  const Box& domain,
//                                  BCRec  const* pbc)
//   {
//     const Dim3 dlo = amrex::lbound(domain);
//     const Dim3 dhi = amrex::ubound(domain);

//     amrex::ParallelFor(xbx, ybx, //zbx,
//     [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//     {
//         // We only care about x-velocity on x-faces here
//         if (flag(i,j,k).isConnected(-1,0,0))
//         {
//             constexpr int n = 0;

//             Real lo = Ipx(i-1,j,k,n);
//             Real hi = Imx(i  ,j,k,n);

//             auto bc = pbc[n];
//             SetTransTermXBCs(i, j, k, n, vel, lo, hi, bc.lo(0), bc.hi(0), dlo.x, dhi.x, true);

//             Real st = ( (lo+hi) >= 0.) ? lo : hi;
//             bool ltm = ( (lo <= 0. && hi >= 0.) || (amrex::Math::abs(lo+hi) < small_vel) );
//             u_ad(i,j,k) = ltm ? 0. : st;
//         } else {
//             u_ad(i,j,k) = 0.;
//         }
//     },
//     [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//     {
//         // We only care about y-velocity on y-faces here
//         if (flag(i,j,k).isConnected(0,-1,0))
//         {
//             constexpr int n = 1;

//             Real lo = Ipy(i,j-1,k,n);
//             Real hi = Imy(i,j  ,k,n);

//             auto bc = pbc[n];
//             SetTransTermYBCs(i, j, k, n, vel, lo, hi, bc.lo(1), bc.hi(1), dlo.y, dhi.y, true);

//             Real st = ( (lo+hi) >= 0.) ? lo : hi;
//             bool ltm = ( (lo <= 0. && hi >= 0.) || (amrex::Math::abs(lo+hi) < small_vel) );
//             v_ad(i,j,k) = ltm ? 0. : st;
//         } else {
//             v_ad(i,j,k) = 0.;
//         }
//     });
// }

void
EBGodunov::ExtrapVelToFacesOnBox (Box const& bx, int ncomp,
                                  Box const& xbx, Box const& ybx,
                                  Box const& xebx, Box const& yebx,
                                  Array4<Real> const& qx,
                                  Array4<Real> const& qy,
                                  Array4<Real const> const& q,
                                  Array4<Real const> const& u_ad,
                                  Array4<Real const> const& v_ad,
                                  Array4<Real> const& Imx,
                                  Array4<Real> const& Imy,
                                  Array4<Real> const& Ipx,
                                  Array4<Real> const& Ipy,
                                  Array4<Real const> const& f,
                                  const Box& domain, const Real* dx_arr,
                                  Real l_dt, BCRec  const* pbc,
                                  Array4<EBCellFlag const> const& flag,
                                  Array4<Real const> const& apx,
                                  Array4<Real const> const& apy,
                                  Array4<Real const> const& fcx,
                                  Array4<Real const> const& fcy,
                                  Real* p)
{
    const Dim3 dlo = amrex::lbound(domain);
    const Dim3 dhi = amrex::ubound(domain);

    Real dx = dx_arr[0];
    Real dy = dx_arr[1];

    Array4<Real> xlo = makeArray4(p, xebx, ncomp);
    p += xlo.size();
    Array4<Real> xhi = makeArray4(p, xebx, ncomp);
    p += xhi.size();
    Array4<Real> ylo = makeArray4(p, yebx, ncomp);
    p += ylo.size();
    Array4<Real> yhi = makeArray4(p, yebx, ncomp);
    p += yhi.size();

    //
    // Upwind the states that result from plm_x and plm_y
    //
    amrex::ParallelFor(
        xebx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo = Ipx(i-1,j,k,n);
            Real hi = Imx(i  ,j,k,n);

            Real uad = u_ad(i,j,k);
            auto bc = pbc[n];

            SetTransTermXBCs(i, j, k, n, q, lo, hi, bc.lo(0), bc.hi(0), dlo.x, dhi.x, true);

            xlo(i,j,k,n) = lo;
            xhi(i,j,k,n) = hi;

            Real st = (uad >= 0.) ? lo : hi;
            Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
            Imx(i, j, k, n) = fu*st + (1.0 - fu) *0.5 * (hi + lo); // store xedge
        },
        yebx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo = Ipy(i,j-1,k,n);
            Real hi = Imy(i,j  ,k,n);

            Real vad = v_ad(i,j,k);
            auto bc = pbc[n];

            SetTransTermYBCs(i, j, k, n, q, lo, hi, bc.lo(1), bc.hi(1), dlo.y, dhi.y, true);

            ylo(i,j,k,n) = lo;
            yhi(i,j,k,n) = hi;

            Real st = (vad >= 0.) ? lo : hi;
            Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
            Imy(i, j, k, n) = fu*st + (1.0 - fu)*0.5*(hi + lo); // store yedge
        });


    // We can reuse the space in Ipy

    //
    // X-Flux
    //
    Box const xbxtmp = Box(xbx).enclosedCells().grow(0,1);
    Array4<Real> yhat = makeArray4(Ipy.dataPtr(), Box(v_ad), 1);
    // Add du/dy term to u on x-faces
    // Start with {ylo,yhi} --> upwind using vad to {yhat}

    amrex::ParallelFor(Box(yhat),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (flag(i,j,k).isConnected(0,-1,0))
        {
            constexpr int n = 0;
            const auto bc = pbc[n];
            Real l_yzlo, l_yzhi;

            l_yzlo = ylo(i,j,k,n);
            l_yzhi = yhi(i,j,k,n);
            Real vad = v_ad(i,j,k);
            SetTransTermYBCs(i, j, k, n, q, l_yzlo, l_yzhi, bc.lo(1), bc.hi(1), dlo.y, dhi.y, true);

            Real st = (vad >= 0.) ? l_yzlo : l_yzhi;
            Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
            yhat(i,j,k) = fu*st + (1.0 - fu) * 0.5 * (l_yzhi + l_yzlo);
        } else {
            yhat(i,j,k) = 0.0;
        }
    });


    //
    // Define du/dy and add (v du/dy) to u on x-faces
    //
    amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (flag(i,j,k).isConnected(-1,0,0))
        {
        constexpr int n = 0;
        auto bc = pbc[n];

        // stl is on the left  side of the lo-x side of cell (i,j)
        // sth is on the right side of the lo-x side of cell (i,j)
        Real stl = xlo(i,j,k,n);
        Real sth = xhi(i,j,k,n);

        Real trans_y;

        // Left side of interface
        {
            int ic = i-1;
            if (flag(ic,j,k).isRegular())
            {
                // For full cells this is the transverse term
                stl += - (0.25*l_dt/dy)*(v_ad(ic,j+1,k)+v_ad(ic,j,k))*
                                        (yhat(ic,j+1,k)-yhat(ic,j,k));
                stl += 0.5 * l_dt * f(ic,j,k,n);

            } else {

                // If either y-face is covered then don't include any dt-based terms
                if (apy(ic,j,k) > 0.0 && apy(ic,j+1,k) > 0.0)
                {
                    create_transverse_terms_for_xface(ic,j,k,v_ad,yhat,apy,fcy,trans_y,dy);

                    stl += -0.5 * l_dt * trans_y;
                    stl +=  0.5 * l_dt * f(ic,j,k,n);
                }
            }
        }

        // Right side of interface
        {
            int ic = i;
            if (flag(ic,j,k).isRegular())
            {

                // For full cells this is the transverse term
                sth += - (0.25*l_dt/dy)*(v_ad(ic,j+1,k)+v_ad(ic,j,k))*
                                        (yhat(ic,j+1,k)-yhat(ic,j,k));
                sth +=  0.5 * l_dt * f(ic,j,k,n);

            } else {

                // If either y-face is covered then don't include any dt-based terms
                if (apy(ic,j,k) > 0.0 && apy(ic,j+1,k) > 0.0)
                {
                    create_transverse_terms_for_xface(ic,j,k,v_ad,yhat,apy,fcy,trans_y,dy);

                    sth += -0.5 * l_dt * trans_y;
                    sth +=  0.5 * l_dt * f(ic,j,k,n);
                }
            }
        }

        SetXEdgeBCs(i, j, k, n, q, stl, sth, bc.lo(0), dlo.x, bc.hi(0), dhi.x, true);

        // Prevent backflow
        if ( (i==dlo.x) and (bc.lo(0) == BCType::foextrap || bc.lo(0) == BCType::hoextrap) )
        {
#ifndef ALLOWXINFLOW
            sth = amrex::min(sth,0.);
#endif
            stl = sth;
        }
        if ( (i==dhi.x+1) and (bc.hi(0) == BCType::foextrap || bc.hi(0) == BCType::hoextrap) )
        {
#ifndef ALLOWXINFLOW
             stl = amrex::max(stl,0.);
#endif
             sth = stl;
        }
        Real st = ( (stl+sth) >= 0.) ? stl : sth;
        bool ltm = ( (stl <= 0. && sth >= 0.) || (amrex::Math::abs(stl+sth) < small_vel) );
        qx(i,j,k) = ltm ? 0. : st;

        } else {
            qx(i,j,k) = 0.;
        }
    });

    //
    // Y-Flux
    //
    Box const ybxtmp = Box(ybx).enclosedCells().grow(1,1);
    Array4<Real> xhat = makeArray4(Ipy.dataPtr(), Box(u_ad), 1);

    // Start with {xlo,xhi} --> upwind using uad to {xhat}
    amrex::ParallelFor(Box(xhat),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (flag(i,j,k).isConnected(-1,0,0))
        {
            constexpr int n = 1;
            const auto bc = pbc[n];
            Real l_xzlo, l_xzhi;

            l_xzlo = xlo(i,j,k,n);
            l_xzhi = xhi(i,j,k,n);

            Real uad = u_ad(i,j,k);
            SetTransTermXBCs(i, j, k, n, q, l_xzlo, l_xzhi, bc.lo(0), bc.hi(0), dlo.x, dhi.x, true);

            Real st = (uad >= 0.) ? l_xzlo : l_xzhi;
            Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
            xhat(i,j,k) = fu*st + (1.0 - fu) * 0.5 * (l_xzhi + l_xzlo);
        } else {
            xhat(i,j,k) = 0.0;
        }
    });

    //
    // Define dv/dx and add (u dv/dx) to v on y-faces
    //
    amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (flag(i,j,k).isConnected(0,-1,0))
        {
        constexpr int n = 1;
        auto bc = pbc[n];

        // stl is on the low  side of the lo-y side of cell (i,j)
        // sth is on the high side of the lo-y side of cell (i,j)
        Real stl = ylo(i,j,k,n);
        Real sth = yhi(i,j,k,n);

        Real trans_x;

        // d/dx computed in (i,j-1)
        {
            int jc = j-1;
            if (flag(i,jc,k).isRegular())
            {
                // For full cells this is the transverse term
                stl += - (0.25*l_dt/dx)*(u_ad(i+1,jc,k)+u_ad(i,jc,k))*
                                        (xhat(i+1,jc,k)-xhat(i,jc,k));
                stl += 0.5 * l_dt * f(i,jc,k,n);

            } else {

                // If either x-face is covered then don't include any dt-based terms
                if (apx(i,jc,k) > 0.0 && apx(i+1,jc,k) > 0.0)
                {
                    create_transverse_terms_for_yface(i,jc,k,u_ad,xhat,apx,fcx,trans_x,dx);

                    stl += -0.5 * l_dt * trans_x;
                    stl +=  0.5 * l_dt * f(i,jc,k,n);
                }
            }
        }

        // d/dx computed in (i,j)
        {
            int jc = j;
            if (flag(i,jc,k).isRegular())
            {
                // For full cells this is the transverse term
                sth += - (0.25*l_dt/dx)*(u_ad(i+1,jc,k)+u_ad(i,jc,k))*
                                        (xhat(i+1,jc,k)-xhat(i,jc,k));
                sth += 0.5 * l_dt * f(i,jc,k,n);

            } else {

                // If either x-face is covered then don't include any dt-based terms
                if (apx(i,jc,k) > 0.0 && apx(i+1,jc,k) > 0.0)
                {
                    create_transverse_terms_for_yface(i,jc,k,u_ad,xhat,apx,fcx,trans_x,dx);

                    sth += -0.5 * l_dt * trans_x;
                    sth +=  0.5 * l_dt * f(i,jc,k,n);
                }
            }
        }

        SetYEdgeBCs(i, j, k, n, q, stl, sth, bc.lo(1), dlo.y, bc.hi(1), dhi.y, true);

        // Prevent backflow
        if ( (j==dlo.y) and (bc.lo(1) == BCType::foextrap || bc.lo(1) == BCType::hoextrap) )
        {
#ifndef ALLOWYINFLOW
            sth = amrex::min(sth,0.);
#endif
            stl = sth;
        }
        if ( (j==dhi.y+1) and (bc.hi(1) == BCType::foextrap || bc.hi(1) == BCType::hoextrap) )
        {
#ifndef ALLOWYINFLOW
            stl = amrex::max(stl,0.);
#endif
            sth = stl;
        }

        Real st = ( (stl+sth) >= 0.) ? stl : sth;
        bool ltm = ( (stl <= 0. && sth >= 0.) || (amrex::Math::abs(stl+sth) < small_vel) );
        qy(i,j,k) = ltm ? 0. : st;

        } else {
            qy(i,j,k) = 0.;
        }
    });

}
