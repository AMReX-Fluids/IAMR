#ifndef EBGODUNOV_CORNER_COUPLE_H
#define EBGODUNOV_CORNER_COUPLE_H

#include <AMReX_Gpu.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Array.H>

/* This header file contains the inlined __host__ __device__ functions required for
   the scalar advection routines for 3D Godunov. It also contains function declarations
   for controlling host functions.  */

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void EBGodunov_corner_couple_yx (amrex::Real& lo1, amrex::Real& hi1,
                                 int i, int j, int k, int n, amrex::Real dt, amrex::Real dx,
                                 bool iconserv, amrex::Real lo, amrex::Real hi,
                                 amrex::Array4<amrex::Real const> const& s,
                                 amrex::Array4<amrex::Real const> const& divu_cc,
                                 amrex::Array4<amrex::Real const> const& apx,
                                 amrex::Array4<amrex::Real const> const& apy,
                                 amrex::Array4<amrex::Real const> const& apz,
                                 amrex::Array4<amrex::Real const> const& vfrac,
                                 amrex::Array4<amrex::Real const> const& umac,
                                 amrex::Array4<amrex::Real const> const& state)
{
    lo1 = lo;
    hi1 = hi;

    // Modify state on y-faces with x-derivatives to be used for computing state on z-faces
    if (apy(i,j,k) > 0.)
    {
        if (apx(i+1,j-1,k) > 0. && apx(i,j-1,k) > 0.)
        {
            if (iconserv)
            {
                lo1 += ( - dt/(3.0*dx)*(apx(i+1,j-1,k)*state(i+1,j-1,k,n)*umac(i+1,j-1,k)
                                      - apx(i  ,j-1,k)*state(i  ,j-1,k,n)*umac(i  ,j-1,k))
                         + dt/(3.0*dx)*s(i,j-1,k,n)*(apx(i+1,j-1,k)*umac(i+1,j-1,k)
                                                    -apx(i  ,j-1,k)*umac(i  ,j-1,k)) ) / vfrac(i,j-1,k)
                         - dt/(3.0)   *s(i,j-1,k,n)*divu_cc(i,j-1,k);
            } else {

                lo1 += - dt/(6.0*dx)*( umac(i+1,j-1,k)   +  umac(i,j-1,k))
                                    *(state(i+1,j-1,k,n) - state(i,j-1,k,n));
            }
        }

        if (apx(i+1,j,k) > 0. && apx(i,j,k) > 0.)
        {
            if (iconserv)
            {
                hi1 += ( - dt/(3.0*dx)*(apx(i+1,j,k)*state(i+1,j,k,n)*umac(i+1,j,k)
                                      - apx(i  ,j,k)*state(i  ,j,k,n)*umac(i  ,j,k))
                         + dt/(3.0*dx)*s(i,j,k,n)*(apx(i+1,j,k)*umac(i+1,j,k) - 
                                                   apx(i  ,j,k)*umac(i  ,j,k)) ) / vfrac(i,j,k)
                         - dt/(3.0)   *s(i,j,k,n)*divu_cc(i,j,k);
            } else {

                hi1 += - dt/(6.0*dx)*( umac(i+1,j  ,k)   +  umac(i,j  ,k))
                                    *(state(i+1,j  ,k,n) - state(i,j  ,k,n));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void EBGodunov_corner_couple_zx (amrex::Real& lo1, amrex::Real& hi1,
                                 int i, int j, int k, int n, amrex::Real dt, amrex::Real dx,
                                 bool iconserv, amrex::Real lo, amrex::Real hi,
                                 amrex::Array4<amrex::Real const> const& s,
                                 amrex::Array4<amrex::Real const> const& divu_cc,
                                 amrex::Array4<amrex::Real const> const& apx,
                                 amrex::Array4<amrex::Real const> const& apy,
                                 amrex::Array4<amrex::Real const> const& apz,
                                 amrex::Array4<amrex::Real const> const& vfrac,
                                 amrex::Array4<amrex::Real const> const& umac,
                                 amrex::Array4<amrex::Real const> const& state)
{
    lo1 = lo;
    hi1 = hi;

    // Modify state on z-faces with x-derivatives to be used for computing state on y-faces
    if (apz(i,j,k) > 0.)
    {
        if (apx(i+1,j,k-1) > 0. && apx(i,j,k-1) > 0.)
        {
            if (iconserv)
            {
                lo1 += ( - dt/(3.0*dx)*(apx(i+1,j,k-1)*state(i+1,j,k-1,n)*umac(i+1,j,k-1)
                                      - apx(i  ,j,k-1)*state(i  ,j,k-1,n)*umac(i  ,j,k-1))
                         + dt/(3.0*dx)*s(i,j,k-1,n)*(apx(i+1,j,k-1)*umac(i+1,j,k-1)
                                                    -apx(i  ,j,k-1)*umac(i  ,j,k-1)) ) / vfrac(i,j,k-1) 
                         - dt/(3.0)   *s(i,j,k-1,n)*divu_cc(i,j,k-1);
            } else {

                lo1 += - dt/(6.0*dx)*( umac(i+1,j,k-1)   +  umac(i,j,k-1))
                                    *(state(i+1,j,k-1,n) - state(i,j,k-1,n));

            }
        }
        if (apx(i+1,j,k) > 0. && apx(i,j,k) > 0.)
        {
            if (iconserv)
            {
                hi1 += ( - dt/(3.0*dx)*(apx(i+1,j,k)*state(i+1,j,k,n)*umac(i+1,j,k)
                                      - apx(i  ,j,k)*state(i  ,j,k,n)*umac(i  ,j,k))
                         + dt/(3.0*dx)*s(i,j,k,n)*(apx(i+1,j,k)*umac(i+1,j,k) 
                                                  -apx(i  ,j,k)*umac(i  ,j,k)) ) / vfrac(i,j,k)
                         - dt/(3.0)   *s(i,j,k,n)*divu_cc(i,j,k);
            } else {

                hi1 += - dt/(6.0*dx)*( umac(i+1,j,k  )   +  umac(i,j,k))
                                    *(state(i+1,j,k  ,n) - state(i,j,k  ,n));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void EBGodunov_corner_couple_xy (amrex::Real& lo1, amrex::Real& hi1,
                                 int i, int j, int k, int n, amrex::Real dt, amrex::Real dy,
                                 bool iconserv, amrex::Real lo, amrex::Real hi,
                                 amrex::Array4<amrex::Real const> const& s,
                                 amrex::Array4<amrex::Real const> const& divu_cc,
                                 amrex::Array4<amrex::Real const> const& apx,
                                 amrex::Array4<amrex::Real const> const& apy,
                                 amrex::Array4<amrex::Real const> const& apz,
                                 amrex::Array4<amrex::Real const> const& vfrac,
                                 amrex::Array4<amrex::Real const> const& vmac,
                                 amrex::Array4<amrex::Real const> const& state)
{
    lo1 = lo;
    hi1 = hi;

    // Modify state on x-faces with y-derivatives to be used for computing state on z-faces
    if (apx(i,j,k) > 0.)
    {
        if (apy(i-1,j+1,k) > 0. && apy(i-1,j,k) > 0.)
        {
            if (iconserv)
            {
                lo1 += ( - dt/(3.0*dy)*(apy(i-1,j+1,k)*state(i-1,j+1,k,n)*vmac(i-1,j+1,k)
                                      - apy(i-1,j  ,k)*state(i-1,j  ,k,n)*vmac(i-1,j  ,k))
                         + dt/(3.0*dy)*s(i-1,j,k,n)*(apy(i-1,j+1,k)*vmac(i-1,j+1,k)
                                                    -apy(i-1,j  ,k)*vmac(i-1,j  ,k)) ) /vfrac(i-1,j,k)
                         - dt/3.0     *s(i-1,j,k,n)*divu_cc(i-1,j,k);
            } else {

                lo1 += - dt/(6.0*dy)*( vmac(i-1,j+1,k)   +  vmac(i-1,j,k))
                                    *(state(i-1,j+1,k,n) - state(i-1,j,k,n));
            }
        }
    
        if (apy(i,j+1,k) > 0. && apy(i,j,k) > 0.)
        {
            if (iconserv)
            {
                hi1 += ( - dt/(3.0*dy)*(apy(i,j+1,k)*state(i,j+1,k,n)*vmac(i,j+1,k)
                                       -apy(i,j  ,k)*state(i,j  ,k,n)*vmac(i,j  ,k))
                         + dt/(3.0*dy)*s(i,j,k,n)*(apy(i,j+1,k)*vmac(i,j+1,k)
                                                  -apy(i,j  ,k)*vmac(i,j  ,k)) ) /vfrac(i,j,k)
                         - dt/3.0     *s(i,j,k,n)*divu_cc(i,j,k);
            } else {

                hi1 += - dt/(6.0*dy)*( vmac(i  ,j+1,k)   +  vmac(i  ,j,k))
                                    *(state(i  ,j+1,k,n) - state(i  ,j,k,n));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void EBGodunov_corner_couple_zy (amrex::Real& lo1, amrex::Real& hi1,
                                 int i, int j, int k, int n, amrex::Real dt, amrex::Real dy,
                                 bool iconserv, amrex::Real lo, amrex::Real hi,
                                 amrex::Array4<amrex::Real const> const& s,
                                 amrex::Array4<amrex::Real const> const& divu_cc,
                                 amrex::Array4<amrex::Real const> const& apx,
                                 amrex::Array4<amrex::Real const> const& apy,
                                 amrex::Array4<amrex::Real const> const& apz,
                                 amrex::Array4<amrex::Real const> const& vfrac,
                                 amrex::Array4<amrex::Real const> const& vmac,
                                 amrex::Array4<amrex::Real const> const& state)
{
    lo1 = lo;
    hi1 = hi;

    // Modify state on z-faces with y-derivatives to be used for computing state on x-faces
    if (apz(i,j,k) > 0.)
    {
        if (apy(i,j+1,k-1) > 0. && apy(i,j,k-1) > 0.)
        {
            if (iconserv)
            {
                lo1 += ( - dt/(3.0*dy)*(apy(i,j+1,k-1)*state(i,j+1,k-1,n)*vmac(i,j+1,k-1)
                                      - apy(i,j  ,k-1)*state(i,j  ,k-1,n)*vmac(i,j  ,k-1))
                         + dt/(3.0*dy)*s(i,j,k-1,n)*(apy(i,j+1,k-1)*vmac(i,j+1,k-1)
                                                    -apy(i,j  ,k-1)*vmac(i,j  ,k-1)) ) /vfrac(i,j,k-1)
                         - dt/(3.0)   *s(i,j,k-1,n)*divu_cc(i,j,k-1);
            } else {

                lo1 += - dt/(6.0*dy)*( vmac(i,j+1,k-1)   +  vmac(i,j,k-1))
                                    *(state(i,j+1,k-1,n) - state(i,j,k-1,n));
            }
        }

        if (apy(i,j+1,k) > 0. && apy(i,j,k) > 0.)
        {
            if (iconserv)
            {
                hi1 += ( - dt/(3.0*dy)*(apy(i,j+1,k)*state(i,j+1,k,n)*vmac(i,j+1,k)
                                      - apy(i,j  ,k)*state(i,j  ,k,n)*vmac(i,j  ,k))
                         + dt/(3.0*dy)*s(i,j,k,n)*(apy(i,j+1,k)*vmac(i,j+1,k) 
                                                  -apy(i,j  ,k)*vmac(i,j  ,k)) ) /vfrac(i,j,k)
                         - dt/(3.0)   *s(i,j,k,n)*divu_cc(i,j,k);
            } else {

                hi1 += - dt/(6.0*dy)*( vmac(i,j+1,k)     +  vmac(i,j,k))
                                    *(state(i,j+1,k  ,n) - state(i,j,k  ,n));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void EBGodunov_corner_couple_xz (amrex::Real& lo1, amrex::Real& hi1,
                                 int i, int j, int k, int n, amrex::Real dt, amrex::Real dz,
                                 bool iconserv, amrex::Real lo, amrex::Real hi,
                                 amrex::Array4<amrex::Real const> const& s,
                                 amrex::Array4<amrex::Real const> const& divu_cc,
                                 amrex::Array4<amrex::Real const> const& apx,
                                 amrex::Array4<amrex::Real const> const& apy,
                                 amrex::Array4<amrex::Real const> const& apz,
                                 amrex::Array4<amrex::Real const> const& vfrac,
                                 amrex::Array4<amrex::Real const> const& wmac,
                                 amrex::Array4<amrex::Real const> const& state)
{
    lo1 = lo;
    hi1 = hi;

    // Modify state on x-faces with z-derivatives to be used for computing state on y-faces
    if (apx(i,j,k) > 0.)
    {
        if (apz(i-1,j,k+1) > 0. && apz(i-1,j,k) > 0.)
        {
            if (iconserv)
            {
                lo1 += ( - dt/(3.0*dz)*(apz(i-1,j,k+1)*state(i-1,j,k+1,n)*wmac(i-1,j,k+1)
                                      - apz(i-1,j,k  )*state(i-1,j,k  ,n)*wmac(i-1,j,k  ))
                         + dt/(3.0*dz)*s(i-1,j,k,n)*(apz(i-1,j,k+1)*wmac(i-1,j,k+1) 
                                                    -apz(i-1,j,k  )*wmac(i-1,j,k  )) ) / vfrac(i-1,j,k)
                         - dt/(3.0)   *s(i-1,j,k,n)*divu_cc(i-1,j,k);
            } else {

                lo1 += - dt/(6.0*dz)*( wmac(i-1,j,k+1)   +  wmac(i-1,j,k))
                                    *(state(i-1,j,k+1,n) - state(i-1,j,k,n));
            }
        }

        if (apz(i,j,k+1) > 0. && apz(i,j,k) > 0.)
        {
            if (iconserv)
            {
                hi1 += ( - dt/(3.0*dz)*(apz(i,j,k+1)*state(i,j,k+1,n)*wmac(i,j,k+1)
                                      - apz(i,j,k  )*state(i,j,k  ,n)*wmac(i,j,k  ))
                         + dt/(3.0*dz)*s(i,j,k,n)*(apz(i,j,k+1)*wmac(i,j,k+1) 
                                                  -apz(i,j,k  )*wmac(i,j,k  )) ) / vfrac(i,j,k)
                         - dt/(3.0)   *s(i,j,k,n)*divu_cc(i,j,k);
            } else {

                hi1 += - dt/(6.0*dz)*( wmac(i  ,j,k+1)   +  wmac(i  ,j,k))
                                    *(state(i  ,j,k+1,n) - state(i  ,j,k,n));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void EBGodunov_corner_couple_yz (amrex::Real& lo1, amrex::Real& hi1,
                                 int i, int j, int k, int n, amrex::Real dt, amrex::Real dz,
                                 bool iconserv, amrex::Real lo, amrex::Real hi,
                                 amrex::Array4<amrex::Real const> const& s,
                                 amrex::Array4<amrex::Real const> const& divu_cc,
                                 amrex::Array4<amrex::Real const> const& apx,
                                 amrex::Array4<amrex::Real const> const& apy,
                                 amrex::Array4<amrex::Real const> const& apz,
                                 amrex::Array4<amrex::Real const> const& vfrac,
                                 amrex::Array4<amrex::Real const> const& wmac,
                                 amrex::Array4<amrex::Real const> const& state)
{
    lo1 = lo;
    hi1 = hi;

    // Modify state on y-faces with z-derivatives to be used for computing state on x-faces
    if (apy(i,j,k) > 0.)
    {
        if (apz(i,j-1,k+1) > 0. && apz(i,j-1,k) > 0.)
        {
            if (iconserv)
            {
                lo1 += ( - dt/(3.0*dz)*(apz(i,j-1,k+1)*state(i,j-1,k+1,n)*wmac(i,j-1,k+1)
                                      - apz(i,j-1,k  )*state(i,j-1,k  ,n)*wmac(i,j-1,k))
                         + dt/(3.0*dz)*s(i,j-1,k,n)*(apz(i,j-1,k+1)*wmac(i,j-1,k+1)- 
                                                     apz(i,j-1,k  )*wmac(i,j-1,k  )) ) /vfrac(i,j-1,k)
                         - dt/(3.0)   *s(i,j-1,k,n)*divu_cc(i,j-1,k);
            } else {

                lo1 += - dt/(6.0*dz)*( wmac(i,j-1,k+1)   +  wmac(i,j-1,k))
                                    *(state(i,j-1,k+1,n) - state(i,j-1,k,n));
            }
        }

        if (apz(i,j,k+1) > 0. && apz(i,j,k) > 0.)
        {
            if (iconserv)
            {
                hi1 += ( - dt/(3.0*dz)*(apz(i,j,k+1)*state(i,j,k+1,n)*wmac(i,j,k+1)
                                      - apz(i,j,k  )*state(i,j,k  ,n)*wmac(i,j,k))
                         + dt/(3.0*dz)*s(i,j,k,n)*(apz(i,j,k+1)*wmac(i,j,k+1) 
                                                  -apz(i,j,k)*wmac(i,j,k)) ) / vfrac(i,j,k) 
                         - dt/(3.0)   *s(i,j,k,n)*divu_cc(i,j,k);

            } else {

                hi1 += - dt/(6.0*dz)*( wmac(i,j  ,k+1)   +  wmac(i,j  ,k))
                                    *(state(i,j  ,k+1,n) - state(i,j  ,k,n));
            }
        }
    }

}
}
#endif
