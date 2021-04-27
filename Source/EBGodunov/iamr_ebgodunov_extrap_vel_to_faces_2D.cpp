#include <iamr_ebgodunov_plm.H>
#include <iamr_godunov_plm.H>
#include <iamr_ebgodunov.H>
#include <iamr_godunov.H>
#include <iamr_godunov_K.H>
#include <iamr_ebgodunov_transverse_2D_K.H>

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
