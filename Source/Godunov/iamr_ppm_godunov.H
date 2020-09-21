#ifndef IAMR_PPM_GODUNOV_H
#define IAMR_PPM_GODUNOV_H

#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <iamr_constants_godunov.H>


namespace {
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real vanLeer(const amrex::Real a, const amrex::Real b, const amrex::Real c)
{
    constexpr amrex::Real small_qty_sq = (1.e-10)*(1.e-10);

    amrex::Real dsc = 0.5e0*(b - c);
    amrex::Real dsl = 2.e0*(a - c);
    amrex::Real dsr = 2.e0*(b - a);
    return (dsl*dsr > small_qty_sq) ?
        std::copysign(1., dsc)*amrex::min(amrex::Math::abs(dsc),amrex::min(amrex::Math::abs(dsl), amrex::Math::abs(dsr))) : 0.;
}
}


namespace PPM {


void PredictVelOnFaces (amrex::Box const& bx, int ncomp,
                        AMREX_D_DECL(amrex::Array4<amrex::Real> const& Imx,
                                     amrex::Array4<amrex::Real> const& Imy,
                                     amrex::Array4<amrex::Real> const& Imz),
                        AMREX_D_DECL(amrex::Array4<amrex::Real> const& Ipx,
                                     amrex::Array4<amrex::Real> const& Ipy,
                                     amrex::Array4<amrex::Real> const& Ipz),
                        amrex::Array4<amrex::Real const> const& q,
                        amrex::Array4<amrex::Real const> const& vel,
                        amrex::Geometry geom,
                        amrex::Real dt,
                        amrex::BCRec const* d_bcrec);


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetXBCs ( const int i, const int j, const int k, const int n,
               amrex::Real &sm, amrex::Real &sp,
               amrex::Real &sedge1, amrex::Real &sedge2,
               const amrex::Array4<const amrex::Real> &s,
               const int bclo, const int bchi,
               const int domlo, const int domhi)
{
    using namespace amrex;

    if (bclo == BCType::ext_dir || bclo == BCType::hoextrap)
    {
        if (i == domlo)
        {
            sedge2 = -0.2*s(domlo-1,j,k,n) + 0.75*s(domlo,j, k, n)
                   +  0.5*s(domlo+1,j,k,n) - 0.05*s(domlo+2,j,k,n);
            sedge2 = amrex::max(sedge2, amrex::min(s(domlo+1,j,k,n), s(domlo,j,k,n)));
            sedge2 = amrex::min(sedge2, amrex::max(s(domlo+1,j,k,n), s(domlo,j,k,n)));

            sm = s(domlo-1,j,k,n);
            sp = sedge2;

        } else if (i == domlo+1) {

            sedge1 = -0.2*s(domlo-1,j,k,n) + 0.75*s(domlo  ,j,k,n)
                   +  0.5*s(domlo+1,j,k,n) - 0.05*s(domlo+2,j,k,n);
            sedge1 = amrex::max(sedge1, amrex::min(s(domlo+1,j,k,n), s(domlo,j,k,n)));
            sedge1 = amrex::min(sedge1, amrex::max(s(domlo+1,j,k,n), s(domlo,j,k,n)));

            sp = sedge2;
            sm = sedge1;

            if ( (sp - s(domlo+1,j,k,n))*(s(domlo+1,j,k,n) - sm) <= 0.)
            {
                sp = s(domlo+1,j,k,n);
                sm = s(domlo+1,j,k,n);
            }
            else if(amrex::Math::abs(sp - s(domlo+1,j,k,n)) >= 2.*amrex::Math::abs(sm - s(domlo+1,j,k,n)))
                sp = 3.*s(domlo+1,j,k,n) - 2.*sm;
            else if(amrex::Math::abs(sm - s(domlo+1,j,k,n)) >= 2.*amrex::Math::abs(sp - s(domlo+1,j,k,n)))
                sm = 3.*s(domlo+1,j,k,n) - 2.*sp;
        }
    }

    if (bchi == BCType::ext_dir || bchi == BCType::hoextrap)
    {
        if (i == domhi)
        {
            sedge1 =  -0.2*s(domhi+1,j,k,n) + 0.75*s(domhi,j,k, n)
                      +0.5*s(domhi-1,j,k,n) - 0.05*s(domhi-2,j,k,n);
            sedge1 = amrex::max(sedge1, amrex::min(s(domhi-1,j,k,n), s(domhi,j,k,n)));
            sedge1 = amrex::min(sedge1, amrex::max(s(domhi-1,j,k,n), s(domhi,j,k,n)));

            sp = s(domhi+1,j,k,n);
            sm = sedge1;

        } else if (i == domhi-1) {

            sedge2 = -0.2*s(domhi+1,j,k,n) + 0.75*s(domhi  ,j,k,n)
                     +0.5*s(domhi-1,j,k,n) - 0.05*s(domhi-2,j,k,n);
            sedge2 = amrex::max(sedge2, amrex::min(s(domhi-1,j,k,n), s(domhi,j,k,n)));
            sedge2 = amrex::min(sedge2, amrex::max(s(domhi-1,j,k,n), s(domhi,j,k,n)));

            sp = sedge2;
            sm = sedge1;

            if( (sp - s(domhi-1,j,k,n))*(s(domhi-1,j,k,n) - sm) <= 0.)
            {
                sp = s(domhi-1,j,k,n);
                sm = s(domhi-1,j,k,n);
            }
            else if(amrex::Math::abs(sp - s(domhi-1,j,k,n)) >= 2.*amrex::Math::abs(sm - s(domhi-1,j,k,n)))
                sp = 3.*s(domhi-1,j,k,n) - 2.*sm;
            else if(amrex::Math::abs(sm - s(domhi-1,j,k,n)) >= 2.*amrex::Math::abs(sp - s(domhi-1,j,k,n)))
                sm = 3.*s(domhi-1,j,k,n) - 2.*sp;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetYBCs ( const int i, const int j, const int k, const int n,
               amrex::Real &sm, amrex::Real &sp,
               amrex::Real &sedge1, amrex::Real &sedge2,
               const amrex::Array4<const amrex::Real> &s,
               const int bclo, const int bchi,
               const int domlo, const int domhi)
{
    using namespace amrex;

    if (bclo == BCType::ext_dir || bclo == BCType::hoextrap)
    {
        if (j == domlo)
        {
            sedge2 = -0.2*s(i,domlo-1,k,n) + 0.75*s(i,domlo  ,k,n)
                     +0.5*s(i,domlo+1,k,n) - 0.05*s(i,domlo+2,k,n);
            sedge2 = amrex::max(sedge2, amrex::min(s(i,domlo+1,k,n), s(i,domlo,k,n)));
            sedge2 = amrex::min(sedge2, amrex::max(s(i,domlo+1,k,n), s(i,domlo,k,n)));

            sm = s(i,domlo-1,k,n);
            sp = sedge2;

        } else if (j == domlo+1) {

            sedge1 = -0.2*s(i,domlo-1,k,n) + 0.75*s(i,domlo  ,k,n)
                     +0.5*s(i,domlo+1,k,n) - 0.05*s(i,domlo+2,k,n);
            sedge1 = amrex::max(sedge1, amrex::min(s(i,domlo+1,k,n), s(i,domlo,k,n)));
            sedge1 = amrex::min(sedge1, amrex::max(s(i,domlo+1,k,n), s(i,domlo,k,n)));

            sp = sedge2;
            sm = sedge1;

            if ( (sp - s(i,domlo+1,k,n))*(s(i,domlo+1,k,n) - sm) <= 0.)
            {
                sp = s(i,domlo+1,k,n);
                sm = s(i,domlo+1,k,n);
            }
            else if(amrex::Math::abs(sp - s(i,domlo+1,k,n)) >= 2.*amrex::Math::abs(sm - s(i,domlo+1,k,n)))
                sp = 3.*s(i,domlo+1,k,n) - 2.*sm;
            else if(amrex::Math::abs(sm - s(i,domlo+1,k,n)) >= 2.*amrex::Math::abs(sp - s(i,domlo+1,k,n)))
                sm = 3.*s(i,domlo+1,k,n) - 2.*sp;
        }
    }

    if (bchi == BCType::ext_dir || bchi == BCType::hoextrap)
    {
        if (j == domhi)
        {
            sedge1 =  -0.2*s(i,domhi+1,k,n) + 0.75*s(i,domhi  ,k,n)
                      +0.5*s(i,domhi-1,k,n) - 0.05*s(i,domhi-2,k,n);
            sedge1 = amrex::max(sedge1, amrex::min(s(i,domhi-1,k,n), s(i,domhi,k,n)));
            sedge1 = amrex::min(sedge1, amrex::max(s(i,domhi-1,k,n), s(i,domhi,k,n)));

            sp = s(i,domhi+1,k, n);
            sm = sedge1;

        } else if (j == domhi-1) {

            sedge2 =  -0.2*s(i,domhi+1,k,n) + 0.75*s(i,domhi  ,k,n)
                      +0.5*s(i,domhi-1,k,n) - 0.05*s(i,domhi-2,k,n);
            sedge2 = amrex::max(sedge2, amrex::min(s(i,domhi-1,k,n), s(i,domhi,k,n)));
            sedge2 = amrex::min(sedge2, amrex::max(s(i,domhi-1,k,n), s(i,domhi,k,n)));

            sp = sedge2;
            sm = sedge1;

            if( (sp - s(i,domhi-1,k,n))*(s(i,domhi-1,k,n) - sm) <= 0.){
                sp = s(i,domhi-1,k,n);
                sm = s(i,domhi-1,k,n);
            }
            else if(amrex::Math::abs(sp - s(i,domhi-1,k,n)) >= 2.*amrex::Math::abs(sm - s(i,domhi-1,k,n)))
                sp = 3.*s(i,domhi-1,k,n) - 2.*sm;

            else if(amrex::Math::abs(sm - s(i,domhi-1,k,n)) >= 2.*amrex::Math::abs(sp - s(i,domhi-1,k,n)))
                sm = 3.*s(i,domhi-1,k,n) - 2.*sp;
        }
    }
}

#if (AMREX_SPACEDIM==3)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetZBCs ( const int i, const int j, const int k, const int n,
               amrex::Real &sm, amrex::Real &sp,
               amrex::Real &sedge1, amrex::Real &sedge2,
               const amrex::Array4<const amrex::Real> &s,
               const int bclo, const int bchi,
               const int domlo, const int domhi)
{
    using namespace amrex;

    if (bclo == BCType::ext_dir || bclo == BCType::hoextrap)
    {

        if (k == domlo)
        {
            sedge2 = -0.2*s(i,j,domlo-1,n) + 0.75*s(i,j,domlo  ,n)
                     +0.5*s(i,j,domlo+1,n) - 0.05*s(i,j,domlo+2,n);
            sedge2 = amrex::max(sedge2, amrex::min(s(i,j,domlo+1,n), s(i,j,domlo,n)));
            sedge2 = amrex::min(sedge2, amrex::max(s(i,j,domlo+1,n), s(i,j,domlo,n)));

            sm = s(i,j,domlo-1,n);
            sp = sedge2;

        } else if (k == domlo+1) {

            sedge1 = -0.2*s(i,j,domlo-1,n) + 0.75*s(i,j,domlo  ,n)
                     +0.5*s(i,j,domlo+1,n) - 0.05*s(i,j,domlo+2,n);
            sedge1 = amrex::max(sedge1, amrex::min(s(i,j,domlo+1,n), s(i,j,domlo,n)));
            sedge1 = amrex::min(sedge1, amrex::max(s(i,j,domlo+1,n), s(i,j,domlo,n)));

            sp = sedge2;
            sm = sedge1;

            if ( (sp - s(i,j,domlo+1,n))*(s(i,j,domlo+1,n) - sm) <= 0. )
            {
                sp = s(i,j,domlo+1,n);
                sm = s(i,j,domlo+1,n);
            }
            else if(amrex::Math::abs(sp - s(i,j,domlo+1,n)) >= 2.*amrex::Math::abs(sm - s(i,j,domlo+1,n)))
                sp = 3.*s(i,j,domlo+1,n) - 2.*sm;

            else if(amrex::Math::abs(sm - s(i,j,domlo+1,n)) >= 2.*amrex::Math::abs(sp - s(i,j,domlo+1,n)))
                sm = 3.*s(i,j,domlo+1,n) - 2.*sp;
        }
    }

    if (bchi == BCType::ext_dir || bchi == BCType::hoextrap)
    {
        if (k == domhi)
        {
            sedge1 =  -0.2*s(i,j,domhi+1,n) + 0.75*s(i,j,domhi  ,n)
                      +0.5*s(i,j,domhi-1,n) - 0.05*s(i,j,domhi-2,n);
            sedge1 = amrex::max(sedge1, amrex::min(s(i,j,domhi-1,n), s(i,j,domhi,n)));
            sedge1 = amrex::min(sedge1, amrex::max(s(i,j,domhi-1,n), s(i,j,domhi,n)));

            sp = s(i,j,domhi+1,n);
            sm = sedge1;

        } else if (k == domhi-1) {

            sedge2 = -0.2*s(i,j,domhi+1,n) + 0.75*s(i,j,domhi  ,n)
                     +0.5*s(i,j,domhi-1,n) - 0.05*s(i,j,domhi-2,n);
            sedge2 = amrex::max(sedge2, amrex::min(s(i,j,domhi-1,n), s(i,j,domhi,n)));
            sedge2 = amrex::min(sedge2, amrex::max(s(i,j,domhi-1,n), s(i,j,domhi,n)));

            sp = sedge2;
            sm = sedge1;

            if ( (sp - s(i,j,domhi-1,n))*(s(i,j,domhi-1,n) - sm) <= 0. )
            {
                sp = s(i,j,domhi-1,n);
                sm = s(i,j,domhi-1,n);
            }
            else if(amrex::Math::abs(sp - s(i,j,domhi-1,n)) >= 2.*amrex::Math::abs(sm - s(i,j,domhi-1,n)))
                sp = 3.*s(i,j,domhi-1,n) - 2.*sm;

            else if(amrex::Math::abs(sm - s(i,j,domhi-1,n)) >= 2.*amrex::Math::abs(sp - s(i,j,domhi-1,n)))
                sm = 3.*s(i,j,domhi-1,n) - 2.*sp;

        }
    }
}
#endif


// Right now only ppm type 1 is supported on GPU
// This version is called before the MAC projection, when we use the cell-centered velocity
//      for upwinding
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void PredictVelOnXFace ( const int i, const int j, const int k, const int n,
                         const amrex::Real dtdx, const amrex::Real v_ad,
                         const amrex::Array4<const amrex::Real> &S,
                         const amrex::Array4<amrex::Real> &Im,
                         const amrex::Array4<amrex::Real> &Ip,
                         const amrex::BCRec bc, const int domlo, const int domhi)
{
    using namespace amrex;

    constexpr amrex::Real sixth = 1.0/6.0;

    amrex::Real sedge1, sedge2, sp, sm;

    amrex::Real sm2 = S(i-2,j,k,n);
    amrex::Real sm1 = S(i-1,j,k,n);
    amrex::Real s0  = S(i  ,j,k,n);
    amrex::Real sp1 = S(i+1,j,k,n);
    amrex::Real sp2 = S(i+2,j,k,n);

    amrex::Real d1 = vanLeer(s0,sp1,sm1);
    amrex::Real d2 = vanLeer(sm1,s0,sm2);

    sedge1 = 0.5e0*(s0 + sm1) - sixth*(d1 - d2);
    sedge1 = amrex::min(amrex::max(sedge1, amrex::min(s0, sm1)),amrex::max(s0,sm1));

    d1 = vanLeer(sp1,sp2,s0);
    d2 = vanLeer(s0,sp1,sm1);

    sedge2 = 0.5e0*(sp1 + s0) - sixth*(d1 - d2);
    sedge2 = amrex::min(amrex::max(sedge2, amrex::min(s0, sp1)),amrex::max(s0,sp1));

    sm = sedge1;
    sp = sedge2;

    if ((sedge2-s0)*(s0-sedge1) < 0.e0)
    {
        sp = s0;
        sm = s0;

    } else if (amrex::Math::abs(sedge2-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge1-s0))
        sp = 3.0*s0 - 2.0*sedge1;

      else if (amrex::Math::abs(sedge1-S(i,j,k,n)) >=  2.0*amrex::Math::abs(sedge2-s0))
        sm = 3.0*s0 - 2.0*sedge2;

    SetXBCs(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(0), bc.hi(0), domlo, domhi);

    amrex::Real s6 = 6.0*s0 - 3.0*(sm + sp);

    amrex::Real sigma = amrex::Math::abs(v_ad)*dtdx;

    if (v_ad > small_vel)
    {
        Ip(i,j,k,n) = sp - (0.5*sigma)*((sp-sm) - (1.0 - 2.0/3.0*sigma)*s6);
        Im(i,j,k,n) = S(i,j,k,n);
    }
    else if (v_ad < -small_vel)
    {
        Ip(i,j,k,n) = S(i,j,k,n);
        Im(i,j,k,n) = sm + (0.5*sigma)*((sp-sm) + (1.0 - 2.0/3.0*sigma)*s6);
    } else
    {
        Ip(i,j,k,n) = S(i,j,k,n);
        Im(i,j,k,n) = S(i,j,k,n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void PredictVelOnYFace ( const int i, const int j, const int k, const int n,
                         const amrex::Real dtdy, const amrex::Real v_ad,
                         const amrex::Array4<const amrex::Real> &S,
                         const amrex::Array4<amrex::Real> &Im,
                         const amrex::Array4<amrex::Real> &Ip,
                         const amrex::BCRec bc, const int domlo, const int domhi)
{
    using namespace amrex;

    constexpr amrex::Real sixth = 1.0/6.0;

    amrex::Real sedge1, sedge2, sp, sm;

    amrex::Real sm2 = S(i,j-2,k,n);
    amrex::Real sm1 = S(i,j-1,k,n);
    amrex::Real s0  = S(i,j  ,k,n);
    amrex::Real sp1 = S(i,j+1,k,n);
    amrex::Real sp2 = S(i,j+2,k,n);

    amrex::Real d1 = vanLeer(s0,sp1,sm1);
    amrex::Real d2 = vanLeer(sm1,s0,sm2);

    sedge1 = 0.5e0*(s0 + sm1) - sixth*(d1 - d2);
    sedge1 = amrex::min(amrex::max(sedge1, amrex::min(s0, sm1)),amrex::max(s0,sm1));

    d1 = vanLeer(sp1,sp2,s0);
    d2 = vanLeer(s0,sp1,sm1);

    sedge2 = 0.5e0*(sp1 + s0) - sixth*(d1 - d2);
    sedge2 = amrex::min(amrex::max(sedge2, amrex::min(s0, sp1)),amrex::max(s0,sp1));

    sm = sedge1;
    sp = sedge2;

    if ((sedge2-s0)*(s0-sedge1) < 0.e0)
    {
        sp = s0;
        sm = s0;

    } else if (amrex::Math::abs(sedge2-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge1-s0))
        sp = 3.0*s0 - 2.0*sedge1;

      else if (amrex::Math::abs(sedge1-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge2-s0))
        sm = 3.0*s0 - 2.0*sedge2;

    SetYBCs(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(1), bc.hi(1), domlo, domhi);

    amrex::Real s6 = 6.0*s0- 3.0*(sm + sp);

    amrex::Real sigma = amrex::Math::abs(v_ad)*dtdy;

    if (v_ad > small_vel)
    {
        Ip(i,j,k,n) = sp - (0.5*sigma)*((sp-sm) - (1.0 - 2.0/3.0*sigma)*s6);
        Im(i,j,k,n) = S(i,j,k,n);
    }
    else if (v_ad < -small_vel)
    {
        Ip(i,j,k,n) = S(i,j,k,n);
        Im(i,j,k,n) = sm + (0.5*sigma)*((sp-sm) + (1.0 - 2.0/3.0*sigma)*s6);
    } else
    {
        Ip(i,j,k,n) = S(i,j,k,n);
        Im(i,j,k,n) = S(i,j,k,n);
    }
}

#if (AMREX_SPACEDIM==3)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void PredictVelOnZFace ( const int i, const int j, const int k, const int n,
                         const amrex::Real dtdz, const amrex::Real v_ad,
                         const amrex::Array4<const amrex::Real> &S,
                         const amrex::Array4<amrex::Real> &Im,
                         const amrex::Array4<amrex::Real> &Ip,
                         const amrex::BCRec bc, const int domlo, const int domhi)
{
    using namespace amrex;

    constexpr amrex::Real sixth = 1.0/6.0;

    amrex::Real sedge1, sedge2, sp, sm;

    amrex::Real sm2 = S(i,j,k-2,n);
    amrex::Real sm1 = S(i,j,k-1,n);
    amrex::Real s0  = S(i,j,k  ,n);
    amrex::Real sp1 = S(i,j,k+1,n);
    amrex::Real sp2 = S(i,j,k+2,n);

    amrex::Real d1 = vanLeer(s0,sp1,sm1);
    amrex::Real d2 = vanLeer(sm1,s0,sm2);

    sedge1 = 0.5e0*(s0 + sm1) - sixth*(d1 - d2);
    sedge1 = amrex::min(amrex::max(sedge1, amrex::min(s0, sm1)),amrex::max(s0,sm1));

    d1 = vanLeer(sp1,sp2,s0);
    d2 = vanLeer(s0,sp1,sm1);

    sedge2 = 0.5e0*(sp1 + s0) - sixth*(d1 - d2);
    sedge2 = amrex::min(amrex::max(sedge2, amrex::min(s0, sp1)),amrex::max(s0,sp1));

    sm = sedge1;
    sp = sedge2;

    if ((sedge2-s0)*(s0-sedge1) < 0.e0)
    {
        sp = s0;
        sm = s0;

    } else if (amrex::Math::abs(sedge2-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge1-s0))
        sp = 3.0*s0 - 2.0*sedge1;

      else if (amrex::Math::abs(sedge1-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge2-s0))
        sm = 3.0*s0 - 2.0*sedge2;

    SetZBCs(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(2), bc.hi(2), domlo, domhi);

    amrex::Real s6 = 6.0*s0- 3.0*(sm + sp);

    amrex::Real sigma = amrex::Math::abs(v_ad)*dtdz;

    if (v_ad > small_vel)
    {
        Ip(i,j,k,n) = sp - (0.5*sigma)*((sp-sm) - (1.0 - 2.0/3.0*sigma)*s6);
        Im(i,j,k,n) = S(i,j,k,n);
    }
    else if (v_ad < -small_vel)
    {
        Ip(i,j,k,n) = S(i,j,k,n);
        Im(i,j,k,n) = sm + (0.5*sigma)*((sp-sm) + (1.0 - 2.0/3.0*sigma)*s6);
    } else
    {
        Ip(i,j,k,n) = S(i,j,k,n);
        Im(i,j,k,n) = S(i,j,k,n);
    }
}

#endif

// Right now only ppm type 1 is supported on GPU
// This version is called after the MAC projection, when we use the MAC-projected velocity
//      for upwinding
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void PredictStateOnXFace ( const int i, const int j, const int k, const int n,
                           const amrex::Real dt, const amrex::Real dx,
                           amrex::Real& Im, amrex::Real& Ip,
                           const amrex::Array4<const amrex::Real> &S,
                           const amrex::Array4<const amrex::Real> &vel_edge,
                           const amrex::BCRec bc,
                           const int domlo, const int domhi)
{
    {
    using namespace amrex;

    constexpr Real sixth = 1.0/6.0;

    amrex::Real sedge1, sedge2, sp, sm;

    amrex::Real sm2 = S(i-2,j,k,n);
    amrex::Real sm1 = S(i-1,j,k,n);
    amrex::Real s0  = S(i  ,j,k,n);
    amrex::Real sp1 = S(i+1,j,k,n);
    amrex::Real sp2 = S(i+2,j,k,n);

    amrex::Real d1 = vanLeer(s0,sp1,sm1);
    amrex::Real d2 = vanLeer(sm1,s0,sm2);

    sedge1 = 0.5e0*(s0 + sm1) - sixth*(d1 - d2);
    sedge1 = amrex::min(amrex::max(sedge1, amrex::min(s0, sm1)),amrex::max(s0,sm1));

    d1 = vanLeer(sp1,sp2,s0);
    d2 = vanLeer(s0,sp1,sm1);

    sedge2 = 0.5e0*(sp1 + s0) - sixth*(d1 - d2);
    sedge2 = amrex::min(amrex::max(sedge2, amrex::min(s0, sp1)),amrex::max(s0,sp1));

    sm = sedge1;
    sp = sedge2;

    if ((sedge2-s0)*(s0-sedge1) < 0.e0)
    {
        sp = s0;
        sm = s0;
    }
    else if (amrex::Math::abs(sedge2-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge1-s0))
      sp = 3.0*s0 - 2.0*sedge1;

    else if (amrex::Math::abs(sedge1-S(i,j,k,n)) >=  2.0*amrex::Math::abs(sedge2-s0))
      sm = 3.0*s0 - 2.0*sedge2;

    SetXBCs(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(0), bc.hi(0), domlo, domhi);

    Real s6 = 6.0*s0 - 3.0*(sm + sp);

    Real sigmap = amrex::Math::abs(vel_edge(i+1,j,k))*dt/dx;
    Real sigmam = amrex::Math::abs(vel_edge(i  ,j,k))*dt/dx;

    if (vel_edge(i+1,j,k) > small_vel)
        Ip = sp - (0.5*sigmap)*((sp - sm) - (1.e0 -2.e0/3.e0*sigmap)*s6);
    else
        Ip = S(i,j,k,n);

    if(vel_edge(i,j,k) < -small_vel)
        Im = sm + (0.5*sigmam)*((sp-sm) + (1.e0 - 2.e0/3.e0*sigmam)*s6);
    else
        Im = S(i,j,k,n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void PredictStateOnYFace ( const int i, const int j, const int k, const int n,
                           const amrex::Real dt, const amrex::Real dx,
                           amrex::Real& Im, amrex::Real& Ip,
                           const amrex::Array4<const amrex::Real> &S,
                           const amrex::Array4<const amrex::Real> &vel_edge,
                           const amrex::BCRec bc,
                           const int domlo, const int domhi)
{
    {
    using namespace amrex;

    constexpr amrex::Real sixth = 1.0/6.0;

    Real sedge1, sedge2, sp, sm;

    Real sm2 = S(i,j-2,k,n);
    Real sm1 = S(i,j-1,k,n);
    Real s0  = S(i,j  ,k,n);
    Real sp1 = S(i,j+1,k,n);
    Real sp2 = S(i,j+2,k,n);

    Real d1 = vanLeer(s0,sp1,sm1);
    Real d2 = vanLeer(sm1,s0,sm2);

    sedge1 = 0.5e0*(s0 + sm1) - sixth*(d1 - d2);
    sedge1 = amrex::min(amrex::max(sedge1, amrex::min(s0, sm1)),amrex::max(s0,sm1));

    d1 = vanLeer(sp1,sp2,s0);
    d2 = vanLeer(s0,sp1,sm1);

    sedge2 = 0.5e0*(sp1 + s0) - sixth*(d1 - d2);
    sedge2 = amrex::min(amrex::max(sedge2, amrex::min(s0, sp1)),amrex::max(s0,sp1));

    sm = sedge1;
    sp = sedge2;

    if ((sedge2-s0)*(s0-sedge1) < 0.e0)
    {
        sp = s0;
        sm = s0;
    }
    else if (amrex::Math::abs(sedge2-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge1-s0))
        sp = 3.0*s0 - 2.0*sedge1;

    else if (amrex::Math::abs(sedge1-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge2-s0))
        sm = 3.0*s0 - 2.0*sedge2;

    SetYBCs(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(1), bc.hi(1), domlo, domhi);

    amrex::Real s6 = 6.0*s0- 3.0*(sm + sp);

    Real sigmap = amrex::Math::abs(vel_edge(i,j+1,k))*dt/dx;
    Real sigmam = amrex::Math::abs(vel_edge(i,j  ,k))*dt/dx;

    if (vel_edge(i,j+1,k) > small_vel)
        Ip = sp - (0.5*sigmap)*((sp - sm) - (1.e0 -2.e0/3.e0*sigmap)*s6);
    else
        Ip = S(i,j,k,n);

    if (vel_edge(i,j,k) < -small_vel)
        Im = sm + (0.5*sigmam)*((sp-sm) + (1.e0 - 2.e0/3.e0*sigmam)*s6);
    else
        Im = S(i,j,k,n);
    }
}


#if (AMREX_SPACEDIM==3)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void PredictStateOnZFace ( const int i, const int j, const int k, const int n,
                           const amrex::Real dt, const amrex::Real dx,
                           amrex::Real& Im, amrex::Real& Ip,
                           const amrex::Array4<const amrex::Real> &S,
                           const amrex::Array4<const amrex::Real> &vel_edge,
                           const amrex::BCRec bc,
                           const int domlo, const int domhi)
{
    {
    using namespace amrex;

    constexpr amrex::Real sixth = 1.0/6.0;

    Real sedge1, sedge2, sp, sm;

    Real sm2 = S(i,j,k-2,n);
    Real sm1 = S(i,j,k-1,n);
    Real s0  = S(i,j,k  ,n);
    Real sp1 = S(i,j,k+1,n);
    Real sp2 = S(i,j,k+2,n);

    Real d1 = vanLeer(s0,sp1,sm1);
    Real d2 = vanLeer(sm1,s0,sm2);

    sedge1 = 0.5e0*(s0 + sm1) - sixth*(d1 - d2);
    sedge1 = amrex::min(amrex::max(sedge1, amrex::min(s0, sm1)),amrex::max(s0,sm1));

    d1 = vanLeer(sp1,sp2,s0);
    d2 = vanLeer(s0,sp1,sm1);

    sedge2 = 0.5e0*(sp1 + s0) - sixth*(d1 - d2);
    sedge2 = amrex::min(amrex::max(sedge2, amrex::min(s0, sp1)),amrex::max(s0,sp1));

    sm = sedge1;
    sp = sedge2;

    if ((sedge2-s0)*(s0-sedge1) < 0.e0)
    {
        sp = s0;
        sm = s0;
    }
    else if (amrex::Math::abs(sedge2-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge1-s0))
        sp = 3.0*s0 - 2.0*sedge1;

    else if (amrex::Math::abs(sedge1-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge2-s0))
        sm = 3.0*s0 - 2.0*sedge2;

    SetZBCs(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(2), bc.hi(2), domlo, domhi);

    Real s6 = 6.0*s0- 3.0*(sm + sp);
    Real sigmap = amrex::Math::abs(vel_edge(i,j,k+1))*dt/dx;
    Real sigmam = amrex::Math::abs(vel_edge(i,j,k  ))*dt/dx;

    if(vel_edge(i,j,k+1) > small_vel)
        Ip = sp - (0.5*sigmap)*((sp-sm) - (1.e0 -2.e0/3.e0*sigmap)*s6);
    else
        Ip = S(i,j,k,n);

    if(vel_edge(i,j,k) < -small_vel)
        Im = sm + (0.5*sigmam)*((sp-sm) + (1.e0 - 2.e0/3.e0*sigmam)*s6);
    else
        Im = S(i,j,k,n);
    }
}
#endif


}

#endif
