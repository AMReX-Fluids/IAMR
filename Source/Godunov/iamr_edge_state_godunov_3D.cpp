#include <iamr_plm_godunov.H>
#include <iamr_ppm_godunov.H>
#include <iamr_godunov.H>
#include <iamr_godunov_K.H>


using namespace amrex;

void
Godunov::ComputeEdgeState (Box const& bx, int ncomp,
                           Array4<Real const> const& q,
                           Array4<Real> const& xedge,
                           Array4<Real> const& yedge,
                           Array4<Real> const& zedge,
                           Array4<Real const> const& umac,
                           Array4<Real const> const& vmac,
                           Array4<Real const> const& wmac,
                           Array4<Real const> const& divu,
                           Array4<Real const> const& fq,
                           Geometry geom,
                           Real l_dt,
                           BCRec const* pbc, int const* iconserv,
                           bool use_ppm,
                           bool use_forces_in_trans,
                           bool is_velocity)
{
    Box const& xbx = amrex::surroundingNodes(bx,0);
    Box const& ybx = amrex::surroundingNodes(bx,1);
    Box const& zbx = amrex::surroundingNodes(bx,2);

    Box const& bxg1 = amrex::grow(bx,1);

    FArrayBox tmpfab(amrex::grow(bx,1),  (4*AMREX_SPACEDIM + 2)*ncomp);
    Elixir tmpeli = tmpfab.elixir();
    Real* p   = tmpfab.dataPtr();

    Box xebox = Box(xbx).grow(1,1).grow(2,1);
    Box yebox = Box(ybx).grow(0,1).grow(2,1);
    Box zebox = Box(zbx).grow(0,1).grow(1,1);

    const Real dx = geom.CellSize(0);
    const Real dy = geom.CellSize(1);
    const Real dz = geom.CellSize(2);

    Real dtdx = l_dt/dx;
    Real dtdy = l_dt/dy;
    Real dtdz = l_dt/dz;

    Box const& domain = geom.Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);
    const auto dxinv = geom.InvCellSizeArray();

    Array4<Real> Imx = makeArray4(p, bxg1, ncomp);
    p +=         Imx.size();
    Array4<Real> Ipx = makeArray4(p, bxg1, ncomp);
    p +=         Ipx.size();
    Array4<Real> Imy = makeArray4(p, bxg1, ncomp);
    p +=         Imy.size();
    Array4<Real> Ipy = makeArray4(p, bxg1, ncomp);
    p +=         Ipy.size();
    Array4<Real> Imz = makeArray4(p, bxg1, ncomp);
    p +=         Imz.size();
    Array4<Real> Ipz = makeArray4(p, bxg1, ncomp);
    p +=         Ipz.size();
    Array4<Real> xlo = makeArray4(p, xebox, ncomp);
    p +=         xlo.size();
    Array4<Real> xhi = makeArray4(p, xebox, ncomp);
    p +=         xhi.size();
    Array4<Real> ylo = makeArray4(p, yebox, ncomp);
    p +=         ylo.size();
    Array4<Real> yhi = makeArray4(p, yebox, ncomp);
    p +=         yhi.size();
    Array4<Real> zlo = makeArray4(p, zebox, ncomp);
    p +=         zlo.size();
    Array4<Real> zhi = makeArray4(p, zebox, ncomp);
    p +=         zhi.size();
    Array4<Real> xyzlo = makeArray4(p, bxg1, ncomp);
    p +=         xyzlo.size();
    Array4<Real> xyzhi = makeArray4(p, bxg1, ncomp);
    p +=         xyzhi.size();

    // Use PPM to generate Im and Ip */
    if (use_ppm)
    {
        amrex::ParallelFor(bxg1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            PPM::PredictStateOnXFace(i, j, k, n, l_dt, dx, Imx(i,j,k,n), Ipx(i,j,k,n),
                                   q, umac, pbc[n], dlo.x, dhi.x);
            PPM::PredictStateOnYFace(i, j, k, n, l_dt, dy, Imy(i,j,k,n), Ipy(i,j,k,n),
                                   q, vmac, pbc[n], dlo.y, dhi.y);
            PPM::PredictStateOnZFace(i, j, k, n, l_dt, dz, Imz(i,j,k,n), Ipz(i,j,k,n),
                                   q, wmac, pbc[n], dlo.z, dhi.z);
        });
    // Use PLM to generate Im and Ip */
    }
    else
    {

        amrex::ParallelFor(xebox, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            PLM::PredictStateOnXFace(i, j, k, n, l_dt, dx, Imx(i,j,k,n), Ipx(i-1,j,k,n),
                                     q, umac(i,j,k), pbc[n], dlo.x, dhi.x, is_velocity);
        });

        amrex::ParallelFor(yebox, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            PLM::PredictStateOnYFace(i, j, k, n, l_dt, dy, Imy(i,j,k,n), Ipy(i,j-1,k,n),
                                     q, vmac(i,j,k), pbc[n], dlo.y, dhi.y, is_velocity);
        });
        amrex::ParallelFor(zebox, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            PLM::PredictStateOnZFace(i, j, k, n, l_dt, dz, Imz(i,j,k,n), Ipz(i,j,k-1,n),
                                     q, wmac(i,j,k), pbc[n], dlo.z, dhi.z, is_velocity);
        });
    }


    amrex::ParallelFor(
    xebox, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Real uad = umac(i,j,k);
        Real fux = (amrex::Math::abs(uad) < small_vel)? 0. : 1.;
        bool uval = uad >= 0.;
        Real lo = Ipx(i-1,j,k,n);
        Real hi = Imx(i  ,j,k,n);

        if (use_forces_in_trans)
        {
            lo += (iconserv[n]) ? -0.5*l_dt*q(i-1,j,k,n)*divu(i-1,j,k) : 0.;
            hi += (iconserv[n]) ? -0.5*l_dt*q(i  ,j,k,n)*divu(i  ,j,k) : 0.;
            lo += 0.5*l_dt*fq(i-1,j,k,n);
            hi += 0.5*l_dt*fq(i  ,j,k,n);
        }

        auto bc = pbc[n];

        SetTransTermXBCs(i, j, k, n, q, lo, hi, uad, bc.lo(0), bc.hi(0), dlo.x, dhi.x, is_velocity);
        xlo(i,j,k,n) = lo;
        xhi(i,j,k,n) = hi;
        Real st = (uval) ? lo : hi;
        Imx(i,j,k,n) = fux*st + (1. - fux)*0.5*(hi + lo);

    },
    yebox, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Real vad = vmac(i,j,k);
        Real fuy = (amrex::Math::abs(vad) < small_vel)? 0. : 1.;
        bool vval = vad >= 0.;
        Real lo = Ipy(i,j-1,k,n);
        Real hi = Imy(i,j  ,k,n);

        if (use_forces_in_trans)
        {
            lo += (iconserv[n]) ? -0.5*l_dt*q(i,j-1,k,n)*divu(i,j-1,k) : 0.;
            hi += (iconserv[n]) ? -0.5*l_dt*q(i,j  ,k,n)*divu(i,j  ,k) : 0.;
            lo += 0.5*l_dt*fq(i,j-1,k,n);
            hi += 0.5*l_dt*fq(i,j  ,k,n);
        }

        auto bc = pbc[n];

        SetTransTermYBCs(i, j, k, n, q, lo, hi, vad, bc.lo(1), bc.hi(1), dlo.y, dhi.y, is_velocity);

        ylo(i,j,k,n) = lo;
        yhi(i,j,k,n) = hi;
        Real st = (vval) ? lo : hi;
        Imy(i,j,k,n) = fuy*st + (1. - fuy)*0.5*(hi + lo);
    },
    zebox, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Real wad = wmac(i,j,k);
        Real fuz = (amrex::Math::abs(wad) < small_vel) ? 0. : 1.;
        bool wval = wad >= 0.;
        Real lo = Ipz(i,j,k-1,n);
        Real hi = Imz(i,j,k  ,n);

        if (use_forces_in_trans)
        {
            lo += (iconserv[n]) ? -0.5*l_dt*q(i,j,k-1,n)*divu(i,j,k-1) : 0.;
            hi += (iconserv[n]) ? -0.5*l_dt*q(i,j,k  ,n)*divu(i,j,k  ) : 0.;
            lo += 0.5*l_dt*fq(i,j,k-1,n);
            hi += 0.5*l_dt*fq(i,j,k  ,n);
        }

        auto bc = pbc[n];

        SetTransTermZBCs(i, j, k, n, q, lo, hi, wad, bc.lo(2), bc.hi(2), dlo.z, dhi.z, is_velocity);

        zlo(i,j,k,n) = lo;
        zhi(i,j,k,n) = hi;
        Real st = (wval) ? lo : hi;
        Imz(i,j,k,n) = fuz*st + (1. - fuz)*0.5*(hi + lo);
    }
    );



    //
    // x-direction
    //
    Box const& xbxtmp = amrex::grow(bx,0,1);
    Array4<Real> yzlo = makeArray4(xyzlo.dataPtr(), amrex::surroundingNodes(xbxtmp,1), ncomp);
    Array4<Real> zylo = makeArray4(xyzhi.dataPtr(), amrex::surroundingNodes(xbxtmp,2), ncomp);
    amrex::ParallelFor(
    Box(zylo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_zylo, l_zyhi;
        AddCornerCoupleTermZY(l_zylo, l_zyhi,
                              i, j, k, n, l_dt, dy, iconserv[n],
                              zlo(i,j,k,n), zhi(i,j,k,n),
                              q, divu, vmac, Imy);

        Real wad = wmac(i,j,k);
        SetTransTermZBCs(i, j, k, n, q, l_zylo, l_zyhi, wad, bc.lo(2), bc.hi(2), dlo.z, dhi.z, is_velocity);

        Real st = (wad >= 0.) ? l_zylo : l_zyhi;
        Real fu = (amrex::Math::abs(wad) < small_vel) ? 0.0 : 1.0;
        zylo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_zyhi + l_zylo);
    },
    Box(yzlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_yzlo, l_yzhi;
        AddCornerCoupleTermYZ(l_yzlo, l_yzhi,
                              i, j, k, n, l_dt, dz, iconserv[n],
                              ylo(i,j,k,n), yhi(i,j,k,n),
                              q, divu, wmac, Imz);

        Real vad = vmac(i,j,k);
        SetTransTermYBCs(i, j, k, n, q, l_yzlo, l_yzhi, vad, bc.lo(1), bc.hi(1), dlo.y, dhi.y, is_velocity);

        Real st = (vad >= 0.) ? l_yzlo : l_yzhi;
        Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
        yzlo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_yzhi + l_yzlo);
    });


    //
    amrex::ParallelFor(xbx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Real stl, sth;
        if (iconserv[n])
        {
            stl = xlo(i,j,k,n) - (0.5*dtdy)*(yzlo(i-1,j+1,k  ,n)*vmac(i-1,j+1,k  )
                                           - yzlo(i-1,j  ,k  ,n)*vmac(i-1,j  ,k  ))
                               - (0.5*dtdz)*(zylo(i-1,j  ,k+1,n)*wmac(i-1,j  ,k+1)
                                           - zylo(i-1,j  ,k  ,n)*wmac(i-1,j  ,k  ))
                + (0.5*dtdy)*q(i-1,j,k,n)*(vmac(i-1,j+1,k  ) - vmac(i-1,j,k))
                + (0.5*dtdz)*q(i-1,j,k,n)*(wmac(i-1,j  ,k+1) - wmac(i-1,j,k));

            sth = xhi(i,j,k,n) - (0.5*dtdy)*(yzlo(i,j+1,k  ,n)*vmac(i,j+1,k  )
                                           - yzlo(i,j  ,k  ,n)*vmac(i,j  ,k  ))
                               - (0.5*dtdz)*(zylo(i,j  ,k+1,n)*wmac(i,j  ,k+1)
                                           - zylo(i,j  ,k  ,n)*wmac(i,j  ,k  ))
                + (0.5*dtdy)*q(i,j,k,n)*(vmac(i,j+1,k  ) - vmac(i,j,k))
                + (0.5*dtdz)*q(i,j,k,n)*(wmac(i,j  ,k+1) - wmac(i,j,k));
        }
        else
        {
            stl = xlo(i,j,k,n) - (0.25*dtdy)*(vmac(i-1,j+1,k  ) + vmac(i-1,j,k)) *
                                             (yzlo(i-1,j+1,k,n) - yzlo(i-1,j,k,n))
                               - (0.25*dtdz)*(wmac(i-1,j,k+1  ) + wmac(i-1,j,k))*
                                             (zylo(i-1,j,k+1,n) - zylo(i-1,j,k,n));

            sth = xhi(i,j,k,n) - (0.25*dtdy)*(vmac(i,j+1,k  ) + vmac(i,j,k))*
                                             (yzlo(i,j+1,k,n) - yzlo(i,j,k,n))
                               - (0.25*dtdz)*(wmac(i,j,k+1  ) + wmac(i,j,k))*
                                             (zylo(i,j,k+1,n) - zylo(i,j,k,n));
        }

        if (not use_forces_in_trans)
        {
            stl += (iconserv[n]) ? -0.5*l_dt*q(i-1,j,k,n)*divu(i-1,j,k) : 0.;
            sth += (iconserv[n]) ? -0.5*l_dt*q(i  ,j,k,n)*divu(i  ,j,k) : 0.;
            stl += 0.5*l_dt*fq(i-1,j,k,n);
            sth += 0.5*l_dt*fq(i  ,j,k,n);
        }

        auto bc = pbc[n];
        SetXEdgeBCs(i, j, k, n, q, stl, sth, umac, bc.lo(0), dlo.x, bc.hi(0), dhi.x, is_velocity);

        if ( (i==dlo.x) and (bc.lo(0) == BCType::foextrap || bc.lo(0) == BCType::hoextrap) )
        {
#ifndef ALLOWXINFLOW
            if ( umac(i,j,k) >= 0. and n==XVEL and is_velocity )  sth = amrex::min(sth,0.);
#endif
            stl = sth;
        }
        if ( (i==dhi.x+1) and (bc.hi(0) == BCType::foextrap || bc.hi(0) == BCType::hoextrap) )
        {
#ifndef ALLOWXINFLOW
            if ( umac(i,j,k) <= 0. and n==XVEL and is_velocity ) stl = amrex::max(stl,0.);
#endif
             sth = stl;
        }

        Real temp = (umac(i,j,k) >= 0.) ? stl : sth;
        temp = (amrex::Math::abs(umac(i,j,k)) < small_vel) ? 0.5*(stl + sth) : temp;
        xedge(i,j,k,n) = temp;
    });

    //
    // y-direction
    //
    Box const& ybxtmp = amrex::grow(bx,1,1);
    Array4<Real> xzlo = makeArray4(xyzlo.dataPtr(), amrex::surroundingNodes(ybxtmp,0), ncomp);
    Array4<Real> zxlo = makeArray4(xyzhi.dataPtr(), amrex::surroundingNodes(ybxtmp,2), ncomp);
    amrex::ParallelFor(
    Box(xzlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_xzlo, l_xzhi;
        AddCornerCoupleTermXZ(l_xzlo, l_xzhi,
                              i, j, k, n, l_dt, dz, iconserv[n],
                              xlo(i,j,k,n),  xhi(i,j,k,n),
                              q, divu, wmac, Imz);

        Real uad = umac(i,j,k);
        SetTransTermXBCs(i, j, k, n, q, l_xzlo, l_xzhi, uad, bc.lo(0), bc.hi(0), dlo.x, dhi.x, is_velocity);

        Real st = (uad >= 0.) ? l_xzlo : l_xzhi;
        Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
        xzlo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_xzhi + l_xzlo);
    },
    Box(zxlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_zxlo, l_zxhi;
        AddCornerCoupleTermZX(l_zxlo, l_zxhi,
                              i, j, k, n, l_dt, dx, iconserv[n],
                              zlo(i,j,k,n), zhi(i,j,k,n),
                              q, divu, umac, Imx);

        Real wad = wmac(i,j,k);
        SetTransTermZBCs(i, j, k, n, q, l_zxlo, l_zxhi, wad, bc.lo(2), bc.hi(2), dlo.z, dhi.z, is_velocity);

        Real st = (wad >= 0.) ? l_zxlo : l_zxhi;
        Real fu = (amrex::Math::abs(wad) < small_vel) ? 0.0 : 1.0;
        zxlo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_zxhi + l_zxlo);
    });

    //
    amrex::ParallelFor(ybx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Real stl, sth;
        if (iconserv[n])
        {
            stl = ylo(i,j,k,n) - (0.5*dtdx)*(xzlo(i+1,j-1,k  ,n)*umac(i+1,j-1,k  )
                                           - xzlo(i  ,j-1,k  ,n)*umac(i  ,j-1,k  ))
                               - (0.5*dtdz)*(zxlo(i  ,j-1,k+1,n)*wmac(i  ,j-1,k+1)
                                           - zxlo(i  ,j-1,k  ,n)*wmac(i  ,j-1,k  ))
                + (0.5*dtdx)*q(i,j-1,k,n)*(umac(i+1,j-1,k  ) - umac(i,j-1,k))
                + (0.5*dtdz)*q(i,j-1,k,n)*(wmac(i  ,j-1,k+1) - wmac(i,j-1,k));

            sth = yhi(i,j,k,n) - (0.5*dtdx)*(xzlo(i+1,j,k  ,n)*umac(i+1,j,k  )
                                           - xzlo(i  ,j,k  ,n)*umac(i  ,j,k  ))
                               - (0.5*dtdz)*(zxlo(i  ,j,k+1,n)*wmac(i  ,j,k+1)
                                           - zxlo(i  ,j,k  ,n)*wmac(i  ,j,k  ))
                + (0.5*dtdx)*q(i,j,k,n)*(umac(i+1,j,k  ) - umac(i,j,k))
                + (0.5*dtdz)*q(i,j,k,n)*(wmac(i  ,j,k+1) - wmac(i,j,k));
        }
        else
        {
            stl = ylo(i,j,k,n) - (0.25*dtdx)*(umac(i+1,j-1,k    ) + umac(i,j-1,k))*
                                             (xzlo(i+1,j-1,k  ,n) - xzlo(i,j-1,k,n))
                               - (0.25*dtdz)*(wmac(i  ,j-1,k+1  ) + wmac(i,j-1,k))*
                                             (zxlo(i  ,j-1,k+1,n) - zxlo(i,j-1,k,n));

            sth = yhi(i,j,k,n) - (0.25*dtdx)*(umac(i+1,j,k  ) + umac(i,j,k))*
                                             (xzlo(i+1,j,k,n) - xzlo(i,j,k,n))
                               - (0.25*dtdz)*(wmac(i,j,k+1  ) + wmac(i,j,k))*
                                             (zxlo(i,j,k+1,n) - zxlo(i,j,k,n));
        }


        if (not use_forces_in_trans)
        {
            stl += (iconserv[n]) ? -0.5*l_dt*q(i,j-1,k,n)*divu(i,j-1,k) : 0.;
            sth += (iconserv[n]) ? -0.5*l_dt*q(i,j  ,k,n)*divu(i,j  ,k) : 0.;
            stl += 0.5*l_dt*fq(i,j-1,k,n);
            sth += 0.5*l_dt*fq(i,j  ,k,n);
        }

        auto bc = pbc[n];
        SetYEdgeBCs(i, j, k, n, q, stl, sth, vmac, bc.lo(1), dlo.y, bc.hi(1), dhi.y, is_velocity);

        if ( (j==dlo.y) and (bc.lo(1) == BCType::foextrap || bc.lo(1) == BCType::hoextrap) )
        {
#ifndef ALLOWYINFLOW
            if ( vmac(i,j,k) >= 0. and n==YVEL and is_velocity ) sth = amrex::min(sth,0.);
#endif
            stl = sth;
        }
        if ( (j==dhi.y+1) and (bc.hi(1) == BCType::foextrap || bc.hi(1) == BCType::hoextrap) )
        {
#ifndef ALLOWYINFLOW
            if ( vmac(i,j,k) <= 0. and n==YVEL and is_velocity ) stl = amrex::max(stl,0.);
#endif
            sth = stl;
        }

        Real temp = (vmac(i,j,k) >= 0.) ? stl : sth;
        temp = (amrex::Math::abs(vmac(i,j,k)) < small_vel) ? 0.5*(stl + sth) : temp;
        yedge(i,j,k,n) = temp;
    });

    //
    // z-direcion
    //
    Box const& zbxtmp = amrex::grow(bx,2,1);
    Array4<Real> xylo = makeArray4(xyzlo.dataPtr(), amrex::surroundingNodes(zbxtmp,0), ncomp);
    Array4<Real> yxlo = makeArray4(xyzhi.dataPtr(), amrex::surroundingNodes(zbxtmp,1), ncomp);
    amrex::ParallelFor(
    Box(xylo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_xylo, l_xyhi;
        AddCornerCoupleTermXY(l_xylo, l_xyhi,
                              i, j, k, n, l_dt, dy, iconserv[n],
                              xlo(i,j,k,n), xhi(i,j,k,n),
                              q, divu, vmac, Imy);

        Real uad = umac(i,j,k);
        SetTransTermXBCs(i, j, k, n, q, l_xylo, l_xyhi, uad, bc.lo(0), bc.hi(0), dlo.x, dhi.x, is_velocity);

        Real st = (uad >= 0.) ? l_xylo : l_xyhi;
        Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
        xylo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_xyhi + l_xylo);
    },
    Box(yxlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_yxlo, l_yxhi;
        AddCornerCoupleTermYX(l_yxlo, l_yxhi,
                              i, j, k, n, l_dt, dx, iconserv[n],
                              ylo(i,j,k,n), yhi(i,j,k,n),
                              q, divu, umac, Imx);

        Real vad = vmac(i,j,k);
        SetTransTermYBCs(i, j, k, n, q, l_yxlo, l_yxhi, vad, bc.lo(1), bc.hi(1), dlo.y, dhi.y, is_velocity);

        Real st = (vad >= 0.) ? l_yxlo : l_yxhi;
        Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
        yxlo(i,j,k,n) = fu*st + (1.0 - fu) * 0.5 * (l_yxhi + l_yxlo);
    });
    //

    amrex::ParallelFor(zbx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Real stl, sth;

        if (iconserv[n])
        {
            stl = zlo(i,j,k,n) - (0.5*dtdx)*(xylo(i+1,j  ,k-1,n)*umac(i+1,j  ,k-1)
                                           - xylo(i  ,j  ,k-1,n)*umac(i  ,j  ,k-1))
                               - (0.5*dtdy)*(yxlo(i  ,j+1,k-1,n)*vmac(i  ,j+1,k-1)
                                           - yxlo(i  ,j  ,k-1,n)*vmac(i  ,j  ,k-1))
                + (0.5*dtdx)*q(i,j,k-1,n)*(umac(i+1,j,k-1) -umac(i,j,k-1))
                + (0.5*dtdy)*q(i,j,k-1,n)*(vmac(i,j+1,k-1) -vmac(i,j,k-1));

            sth = zhi(i,j,k,n) - (0.5*dtdx)*(xylo(i+1,j  ,k,n)*umac(i+1,j  ,k)
                                           - xylo(i  ,j  ,k,n)*umac(i  ,j  ,k))
                               - (0.5*dtdy)*(yxlo(i  ,j+1,k,n)*vmac(i  ,j+1,k)
                                           - yxlo(i  ,j  ,k,n)*vmac(i  ,j  ,k))
                + (0.5*dtdx)*q(i,j,k,n)*(umac(i+1,j,k) -umac(i,j,k))
                + (0.5*dtdy)*q(i,j,k,n)*(vmac(i,j+1,k) -vmac(i,j,k));
        }
        else
        {
            stl = zlo(i,j,k,n) - (0.25*dtdx)*(umac(i+1,j  ,k-1  ) + umac(i,j,k-1))*
                                             (xylo(i+1,j  ,k-1,n) - xylo(i,j,k-1,n))
                               - (0.25*dtdy)*(vmac(i  ,j+1,k-1  ) + vmac(i,j,k-1))*
                                             (yxlo(i  ,j+1,k-1,n) - yxlo(i,j,k-1,n));

            sth = zhi(i,j,k,n) - (0.25*dtdx)*(umac(i+1,j  ,k  ) + umac(i,j,k))*
                                             (xylo(i+1,j  ,k,n) - xylo(i,j,k,n))
                               - (0.25*dtdy)*(vmac(i  ,j+1,k  ) + vmac(i,j,k))*
                                             (yxlo(i  ,j+1,k,n) - yxlo(i,j,k,n));
        }

        if (not use_forces_in_trans)
        {
            stl += (iconserv[n]) ? -0.5*l_dt*q(i,j,k-1,n)*divu(i,j,k-1) : 0.;
            sth += (iconserv[n]) ? -0.5*l_dt*q(i,j,k  ,n)*divu(i,j,k  ) : 0.;
            stl += 0.5*l_dt*fq(i,j,k-1,n);
            sth += 0.5*l_dt*fq(i,j,k  ,n);
        }

        auto bc = pbc[n];
        SetZEdgeBCs(i, j, k, n, q, stl, sth, wmac, bc.lo(2), dlo.z, bc.hi(2), dhi.z, is_velocity);

        if ( (k==dlo.z) and (bc.lo(2) == BCType::foextrap || bc.lo(2) == BCType::hoextrap) )
        {
#ifndef ALLOWZINFLOW
            if ( wmac(i,j,k) >= 0. and n==ZVEL and is_velocity ) sth = amrex::min(sth,0.);
#endif
            stl = sth;
        }
        if ( (k==dhi.z+1) and (bc.hi(2) == BCType::foextrap || bc.hi(2) == BCType::hoextrap) )
        {
#ifndef ALLOWZINFLOW
            if ( wmac(i,j,k) <= 0. and n==ZVEL and is_velocity ) stl = amrex::max(stl,0.);
#endif
            sth = stl;
        }

        Real temp = (wmac(i,j,k) >= 0.) ? stl : sth;
        temp = (amrex::Math::abs(wmac(i,j,k)) < small_vel) ? 0.5*(stl + sth) : temp;
        zedge(i,j,k,n) = temp;
    });

}
