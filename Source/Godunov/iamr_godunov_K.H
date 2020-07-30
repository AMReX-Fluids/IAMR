#ifndef IAMR_GODUNOV_K_H
#define IAMR_GODUNOV_K_H

#include <AMReX_Gpu.H>
// #include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Array.H>
#include <iomanip>
#include <iamr_constants_godunov.H>

/* This header file contains the inlined __host__ __device__ functions required for
   the scalar advection routines for 3D Godunov. It also contains function declarations
   for controlling host functions.  */

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetTransTermXBCs (const int i, const int j, const int k, const int n,
                       const amrex::Array4<const amrex::Real> &s,
                       amrex::Real &lo,
                       amrex::Real &hi,
                       amrex::Real &uedge,
                       const int bclo, const int bchi,
                       const int domlo, const int domhi,
                       const bool is_velocity )
{
    using namespace amrex;


    // Low X
    if (i <= domlo)
    {
        if (bclo==BCType::ext_dir)
        {
            lo = s(domlo-1,j,k,n);
            if( n == XVEL and is_velocity ) hi=lo;
        }
        else if(bclo == BCType::foextrap || bclo == BCType::hoextrap || bclo == BCType::reflect_even)
        {
            lo = hi;
        }
        else if(bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
    // High X
    else if (i > domhi)
    {
        if (bchi==BCType::ext_dir)
        {
            hi = s(domhi+1,j,k,n) ;
            if( n ==XVEL and is_velocity ) lo=hi;
        }
        else if (bchi == BCType::foextrap || bchi == BCType::hoextrap || bchi == BCType::reflect_even)
        {
            hi = lo;
        }
        else if(bchi == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
    else
        return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetTransTermYBCs ( const int i, const int j, const int k, const int n,
                        const amrex::Array4<const amrex::Real> &s,
                        amrex::Real &lo,
                        amrex::Real &hi,
                        amrex::Real vedge,
                        const int bclo, const int bchi,
                        const int domlo, const int domhi,
                        const bool is_velocity )
{
    using namespace amrex;


    // Low Y
    if (j <= domlo)
    {
        if (bclo==BCType::ext_dir)
        {
            lo = s(i,domlo-1,k,n);
            if ( n == YVEL and is_velocity ) hi = lo;
        }
        else if (bclo == BCType::foextrap || bclo == BCType::hoextrap || bclo == BCType::reflect_even)
        {
            lo = hi;
        }
        else if(bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
    // High Y
    else if (j > domhi)
    {
        if (bchi==BCType::ext_dir)
        {
            hi = s(i,domhi+1,k,n);
            if( n == YVEL and is_velocity ) lo = hi ;
        }
        else if (bchi == BCType::foextrap || bchi == BCType::hoextrap || bchi == BCType::reflect_even)
        {
            hi = lo;
        }
        else if (bchi == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
    else
        return;
}

#if (AMREX_SPACEDIM==3)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetTransTermZBCs(const int i, const int j, const int k, const int n,
                      const amrex::Array4<const amrex::Real> &s,
                      amrex::Real &lo,
                      amrex::Real &hi,
                      amrex::Real wedge,
                      const int bclo, const int bchi,
                      const int domlo, const int domhi,
                      const bool is_velocity)
{
    using namespace amrex;


    // Low Z
    if (k <= domlo)
    {
        if (bclo==BCType::ext_dir)
        {
            lo =s(i,j,domlo-1,n);
            if ( n == ZVEL and is_velocity ) hi = lo;
        }
        else if (bclo == BCType::foextrap || bclo == BCType::hoextrap || bclo == BCType::reflect_even)
        {
            lo = hi;
        }
        else if(bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
    // High Z
    else if (k > domhi)
    {
        if (bchi==BCType::ext_dir)
        {
            hi = s(i,j,domhi+1,n);
            if ( n == ZVEL and is_velocity ) lo = hi ;
        }
        else if (bchi == BCType::foextrap || bchi == BCType::hoextrap || bchi == BCType::reflect_even)
        {
            hi = lo;
        }
        else if (bchi == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
    else
        return;
}
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetXEdgeBCs( const int i, const int j, const int k, const int n,
                  const amrex::Array4<const amrex::Real> &s,
                  amrex::Real &lo, amrex::Real &hi,
                  const amrex::Array4< const amrex::Real> &uedge,
                  const int bclo, const int domlo,
                  const int bchi, const int domhi,
                  const bool is_velocity )
{
    using namespace amrex;

    if (i == domlo)
    {
        if (bclo == BCType::ext_dir)
        {
            lo = s(domlo-1, j, k, n);
            if ( n==XVEL and is_velocity ) hi=lo;
        }
        else if ( bclo == BCType::foextrap || bclo == BCType::hoextrap ||
                  bclo == BCType::reflect_even )
        {
            lo = hi;
        }
        else if (bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = hi;
        }
    }
    else if (i == domhi+1)
    {
        if (bchi == BCType::ext_dir)
        {
            hi = s(domhi+1, j, k, n);
            if ( n==XVEL and is_velocity ) lo = hi;
        }
        else if ( bchi == BCType::foextrap || bchi == BCType::hoextrap ||
                  bchi == BCType::reflect_even )
        {
            hi = lo;
        }
        else if (bchi == BCType::reflect_odd)
        {
            lo = 0.;
            hi = lo;
        }
    }
    else
        return;
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetYEdgeBCs ( const int i, const int j, const int k, const int n,
                   const amrex::Array4<const amrex::Real> &s,
                   amrex::Real &lo, amrex::Real &hi,
                   const amrex::Array4< const amrex::Real> &vedge,
                   const int bclo, const int domlo,
                   const int bchi, const int domhi,
                   const bool is_velocity )
{
    using namespace amrex;

    if (j == domlo)
    {
        if (bclo == BCType::ext_dir)
        {
            lo = s(i, domlo-1, k, n);
            if ( n==YVEL and is_velocity ) hi=lo;
        }
        else if ( bclo == BCType::foextrap || bclo == BCType::hoextrap ||
                  bclo == BCType::reflect_even )
        {
            lo = hi;
        }
        else if (bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = hi;
        }
    }
    else  if (j == domhi+1)
    {
        if (bchi == BCType::ext_dir)
        {
            hi = s(i, domhi+1, k, n);
            if ( n==YVEL and is_velocity ) lo=hi;
        }
        else if ( bchi == BCType::foextrap || bchi == BCType::hoextrap ||
                  bchi == BCType::reflect_even )
        {
            hi = lo;
        }
        else if(bchi == BCType::reflect_odd)
        {
            lo = 0.;
            hi = lo;
        }
    }
    else
        return;
}


//                                                                             //
//   =============================== 3D ONLY ===============================   //
//                                                                             //
#if (AMREX_SPACEDIM==3)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetZEdgeBCs ( const int i, const int j, const int k, const int n,
                   const amrex::Array4<const amrex::Real> &s,
                   amrex::Real &lo, amrex::Real &hi,
                   const amrex::Array4< const amrex::Real> &wedge,
                   const int bclo, const int domlo,
                   const int bchi, const int domhi,
                   const bool is_velocity )
{
    using namespace amrex;

    if (k == domlo)
    {
        if (bclo == BCType::ext_dir)
        {
            lo = s(i, j, domlo-1, n);
            if ( n==ZVEL and is_velocity ) hi=lo;
        }
        else if ( bclo == BCType::foextrap || bclo == BCType::hoextrap ||
                  bclo == BCType::reflect_even )
        {
            lo = hi;
        }
        else if(bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = hi;
        }
    }
    else if (k == domhi+1)
    {
        if (bchi == BCType::ext_dir)
        {
            hi = s(i,j,domhi+1, n);
            if ( n==ZVEL and is_velocity ) lo=hi;
        }
        else if ( bchi == BCType::foextrap || bchi == BCType::hoextrap ||
                  bchi == BCType::reflect_even )
        {
            hi = lo;
        }
        else if(bchi == BCType::reflect_odd)
        {
            lo = 0.;
            hi = lo;
        }
    }
    else
        return;
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void AddCornerCoupleTermYX ( amrex::Real& lo1, amrex::Real& hi1,
                             int i, int j, int k, int n, amrex::Real dt, amrex::Real dx,
                             bool iconserv, amrex::Real lo, amrex::Real hi,
                             amrex::Array4<amrex::Real const> const& s,
                             amrex::Array4<amrex::Real const> const& divu_cc,
                             amrex::Array4<amrex::Real const> const& mac,
                             amrex::Array4<amrex::Real const> const& state )
{
    if (iconserv)
    {
        lo1 = lo - dt/(3.e0*dx)*(state(i+1,j-1,k,n)*mac(i+1,j-1,k)
                               - state(i  ,j-1,k,n)*mac(i  ,j-1,k))
            - dt/(3.e0)   *s(i,j-1,k,n)*divu_cc(i,j-1,k)
            + dt/(3.e0*dx)*s(i,j-1,k,n)*(mac(i+1,j-1,k)-mac(i,j-1,k));

        hi1  = hi - dt/(3.e0*dx)*(state(i+1,j,k,n)*mac(i+1,j,k)
                                - state(i  ,j,k,n)*mac(i  ,j,k))
            - dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
            + dt/(3.e0*dx)*s(i,j,k,n)*(mac(i+1,j,k) - mac(i,j,k));
    }
    else
    {
        lo1 = lo - dt/(6.e0*dx)*(mac  (i+1,j-1,k)   +   mac(i,j-1,k))
                              * (state(i+1,j-1,k,n) - state(i,j-1,k,n));
        hi1 = hi - dt/(6.e0*dx)*(mac  (i+1,j  ,k)   +   mac(i,j  ,k))
                              * (state(i+1,j  ,k,n) - state(i,j  ,k,n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void AddCornerCoupleTermZX ( amrex::Real& lo1, amrex::Real& hi1,
                             int i, int j, int k, int n, amrex::Real dt, amrex::Real dx,
                             bool iconserv, amrex::Real lo, amrex::Real hi,
                             amrex::Array4<amrex::Real const> const& s,
                             amrex::Array4<amrex::Real const> const& divu_cc,
                             amrex::Array4<amrex::Real const> const& mac,
                             amrex::Array4<amrex::Real const> const& state )
{
    if (iconserv)
    {
        lo1 = lo - dt/(3.e0*dx)*(state(i+1,j,k-1,n)*mac(i+1,j,k-1)
                               - state(i  ,j,k-1,n)*mac(i  ,j,k-1))
            - dt/(3.e0)   *s(i,j,k-1,n)*divu_cc(i,j,k-1)
            + dt/(3.e0*dx)*s(i,j,k-1,n)*(mac(i+1,j,k-1)-mac(i,j,k-1));

        hi1 = hi - dt/(3.e0*dx)*(state(i+1,j,k,n)*mac(i+1,j,k)
                               - state(i  ,j,k,n)*mac(i  ,j,k))
            - dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
            + dt/(3.e0*dx)*s(i,j,k,n)*(mac(i+1,j,k) - mac(i,j,k));
    }
    else
    {
        lo1 = lo - dt/(6.e0*dx)*(  mac(i+1,j,k-1)   +   mac(i,j,k-1))
                              * (state(i+1,j,k-1,n) - state(i,j,k-1,n));
        hi1 = hi - dt/(6.e0*dx)*(  mac(i+1,j,k  )   +   mac(i,j,k))
                              * (state(i+1,j,k  ,n) - state(i,j,k  ,n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void AddCornerCoupleTermXY ( amrex::Real& lo1, amrex::Real& hi1,
                             int i, int j, int k, int n, amrex::Real dt, amrex::Real dy,
                             bool iconserv, amrex::Real lo, amrex::Real hi,
                             amrex::Array4<amrex::Real const> const& s,
                             amrex::Array4<amrex::Real const> const& divu_cc,
                             amrex::Array4<amrex::Real const> const& mac,
                             amrex::Array4<amrex::Real const> const& state )
{
    if (iconserv)
    {
        lo1 = lo - dt/(3.e0*dy)*(state(i-1,j+1,k,n)*mac(i-1,j+1,k)
                               - state(i-1,j  ,k,n)*mac(i-1,j  ,k))
            - dt/3.e0     *s(i-1,j,k,n)*divu_cc(i-1,j,k)
            + dt/(3.e0*dy)*s(i-1,j,k,n)*(mac(i-1,j+1,k)-mac(i-1,j,k));

        hi1 = hi - dt/(3.e0*dy)*(state(i,j+1,k,n)*mac(i,j+1,k)
                               - state(i,j  ,k,n)*mac(i,j  ,k))
            - dt/3.e0     *s(i,j,k,n)*divu_cc(i,j,k)
            + dt/(3.e0*dy)*s(i,j,k,n)*(mac(i,j+1,k)-mac(i,j,k));
    }
    else
    {
        lo1 = lo - dt/(6.e0*dy)*(  mac(i-1,j+1,k)   +   mac(i-1,j,k))
                              * (state(i-1,j+1,k,n) - state(i-1,j,k,n));
        hi1 = hi - dt/(6.e0*dy)*(mac  (i  ,j+1,k)   +   mac(i  ,j,k))
                              * (state(i  ,j+1,k,n) - state(i  ,j,k,n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void AddCornerCoupleTermZY ( amrex::Real& lo1, amrex::Real& hi1,
                             int i, int j, int k, int n, amrex::Real dt, amrex::Real dy,
                             bool iconserv, amrex::Real lo, amrex::Real hi,
                             amrex::Array4<amrex::Real const> const& s,
                             amrex::Array4<amrex::Real const> const& divu_cc,
                             amrex::Array4<amrex::Real const> const& mac,
                             amrex::Array4<amrex::Real const> const& state )
{
    if (iconserv)
    {
        lo1 = lo - dt/(3.e0*dy)*(state(i,j+1,k-1,n)*mac(i,j+1,k-1)
                               - state(i,j  ,k-1,n)*mac(i,j  ,k-1))
            - dt/(3.e0)   *s(i,j,k-1,n)*divu_cc(i,j,k-1)
            + dt/(3.e0*dy)*s(i,j,k-1,n)*(mac(i,j+1,k-1)-mac(i,j,k-1));

        hi1 = hi - dt/(3.e0*dy)*(state(i,j+1,k,n)*mac(i,j+1,k)
                               - state(i,j  ,k,n)*mac(i,j  ,k))
            - dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
            + dt/(3.e0*dy)*s(i,j,k,n)*(mac(i,j+1,k) - mac(i,j,k));
    }
    else
    {
        lo1 = lo - dt/(6.e0*dy)*(  mac(i,j+1,k-1)   +   mac(i,j,k-1))
                              * (state(i,j+1,k-1,n) - state(i,j,k-1,n));
        hi1 = hi - dt/(6.e0*dy)*(  mac(i,j+1,k)     +   mac(i,j,k))
                              * (state(i,j+1,k  ,n) - state(i,j,k  ,n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void AddCornerCoupleTermXZ ( amrex::Real& lo1, amrex::Real& hi1,
                             int i, int j, int k, int n, amrex::Real dt, amrex::Real dz,
                             bool iconserv, amrex::Real lo, amrex::Real hi,
                             amrex::Array4<amrex::Real const> const& s,
                             amrex::Array4<amrex::Real const> const& divu_cc,
                             amrex::Array4<amrex::Real const> const& mac,
                             amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv)
    {
        lo1 = lo - dt/(3.e0*dz)*(state(i-1,j,k+1,n)*mac(i-1,j,k+1)
                               - state(i-1,j,k  ,n)*mac(i-1,j,k))
            - dt/(3.e0)   *s(i-1,j,k,n)*divu_cc(i-1,j,k)
            + dt/(3.e0*dz)*s(i-1,j,k,n)*(mac(i-1,j,k+1) - mac(i-1,j,k));

        hi1 = hi - dt/(3.e0*dz)*(state(i,j,k+1,n)*mac(i,j,k+1)
                               - state(i,j,k  ,n)*mac(i,j,k))
            - dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
            + dt/(3.e0*dz)*s(i,j,k,n)*(mac(i,j,k+1) - mac(i,j,k));
    }
    else
    {
        lo1 = lo - dt/(6.e0*dz)*(  mac(i-1,j,k+1)   +   mac(i-1,j,k))
                              * (state(i-1,j,k+1,n) - state(i-1,j,k,n));
        hi1 = hi - dt/(6.e0*dz)*(  mac(i  ,j,k+1)   +   mac(i  ,j,k))
                              * (state(i  ,j,k+1,n) - state(i  ,j,k,n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void AddCornerCoupleTermYZ ( amrex::Real& lo1, amrex::Real& hi1,
                             int i, int j, int k, int n, amrex::Real dt, amrex::Real dz,
                             bool iconserv, amrex::Real lo, amrex::Real hi,
                             amrex::Array4<amrex::Real const> const& s,
                             amrex::Array4<amrex::Real const> const& divu_cc,
                             amrex::Array4<amrex::Real const> const& mac,
                             amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv)
    {
        lo1 = lo - dt/(3.e0*dz)*(state(i,j-1,k+1,n)*mac(i,j-1,k+1)
                               - state(i,j-1,k  ,n)*mac(i,j-1,k))
            - dt/(3.e0)   *s(i,j-1,k,n)*divu_cc(i,j-1,k)
            + dt/(3.e0*dz)*s(i,j-1,k,n)*(mac(i,j-1,k+1)-mac(i,j-1,k));

        hi1 = hi - dt/(3.e0*dz)*(state(i,j,k+1,n)*mac(i,j,k+1)
                               - state(i,j,k  ,n)*mac(i,j,k))
            - dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
            + dt/(3.e0*dz)*s(i,j,k,n)*(mac(i,j,k+1) - mac(i,j,k));
    }
    else
    {
        lo1 = lo - dt/(6.e0*dz)*(  mac(i,j-1,k+1)   +   mac(i,j-1,k))
                              * (state(i,j-1,k+1,n) - state(i,j-1,k,n));
        hi1 = hi - dt/(6.e0*dz)*(  mac(i,j  ,k+1)   +   mac(i,j  ,k))
                              * (state(i,j  ,k+1,n) - state(i,j  ,k,n));
    }
}
#endif
}
#endif
