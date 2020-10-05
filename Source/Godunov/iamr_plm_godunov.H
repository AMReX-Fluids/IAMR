#ifndef IAMR_PLM_GODUNOV_H
#define IAMR_PLM_GODUNOV_H

#include <iamr_slopes_godunov_K.H>
#include <AMReX_Gpu.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Array.H>
#include <iomanip>
#include <iamr_constants_godunov.H>

namespace PLM {

void PredictVelOnXFace ( amrex::Box const& bx, int ncomp,
                         amrex::Array4<amrex::Real> const& Imx, amrex::Array4<amrex::Real> const& Ipx,
                         amrex::Array4<amrex::Real const> const& q,
                         amrex::Array4<amrex::Real  const> const& vcc,
                         const amrex::Geometry& geom,
                         amrex::Real dt,
                         amrex::Vector<amrex::BCRec> const& h_bcrec,
                         amrex::BCRec const* d_bcrec );

void PredictVelOnYFace ( amrex::Box const& bx, int ncomp,
                         amrex::Array4<amrex::Real> const& Imy, amrex::Array4<amrex::Real> const& Ipy,
                         amrex::Array4<amrex::Real const> const& q,
                         amrex::Array4<amrex::Real  const> const& vcc,
                         const amrex::Geometry& geom,
                         amrex::Real dt,
                         amrex::Vector<amrex::BCRec> const& h_bcrec,
                         amrex::BCRec const* d_bcrec );

#if (AMREX_SPACEDIM==3)
void PredictVelOnZFace ( amrex::Box const& bx, int ncomp,
                         amrex::Array4<amrex::Real> const& Imz, amrex::Array4<amrex::Real> const& Ipz,
                         amrex::Array4<amrex::Real const> const& q,
                         amrex::Array4<amrex::Real  const> const& vcc,
                         const amrex::Geometry& geom,
                         amrex::Real dt,
                         amrex::Vector<amrex::BCRec> const& h_bcrec,
                         amrex::BCRec const* d_bcrec );
#endif


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void PredictStateOnXFace ( const int i, const int j, const int k, const int n,
                           const amrex::Real dt, const amrex::Real dx,
                           amrex::Real& Im, amrex::Real& Ip,
                           const amrex::Array4<const amrex::Real> &S,
                           const amrex::Real& umac,
                           const amrex::BCRec bc,
                           const int domain_ilo, const int domain_ihi,
                           const bool is_velocity )
{
    using namespace amrex;
    {
        bool extdir_or_ho_ilo = (bc.lo(0) == BCType::ext_dir) or
                                (bc.lo(0) == BCType::hoextrap);
        bool extdir_or_ho_ihi = (bc.hi(0) == BCType::ext_dir) or
                                (bc.hi(0) == BCType::hoextrap);

        Real upls, umns;

        if (i == domain_ilo && (bc.lo(0) == BCType::ext_dir))
        {
            umns = S(i-1,j,k,n);

            if ( n==XVEL and is_velocity )
            {
              upls = S(i-1,j,k,n);
            }
            else
            {
                upls = S(i  ,j,k,n) + 0.5 * (-1.0 - umac * dt/dx) *
                    iamr_ho_xslope_extdir(i  ,j,k,n,S, extdir_or_ho_ilo, extdir_or_ho_ihi, domain_ilo, domain_ihi);
            }

        }
        else if (i == domain_ihi+1 && (bc.hi(0) == BCType::ext_dir))
        {
            upls = S(i  ,j,k,n);

            if ( n==XVEL and is_velocity )
            {
                umns = S(i,j,k,n);
            }
            else
            {
                umns = S(i-1,j,k,n) + 0.5 * ( 1.0 - umac * dt/dx) *
                    iamr_ho_xslope_extdir(i-1,j,k,n,S, extdir_or_ho_ilo, extdir_or_ho_ihi, domain_ilo, domain_ihi);
            }
        }
        else
        {
	    // Note that we still call the "extdir version" here because interior cells one
	    // away from the boundary will still see the boundary condition in the 4th order
	    // slope
            upls = S(i  ,j,k,n) + 0.5 * (-1.0 - umac * dt/dx) *
                iamr_ho_xslope_extdir(i  ,j,k,n,S, extdir_or_ho_ilo, extdir_or_ho_ihi, domain_ilo, domain_ihi);
            umns = S(i-1,j,k,n) + 0.5 * ( 1.0 - umac * dt/dx) *
                iamr_ho_xslope_extdir(i-1,j,k,n,S, extdir_or_ho_ilo, extdir_or_ho_ihi, domain_ilo, domain_ihi);
        }

        Ip = umns;
        Im = upls;
    }
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void PredictStateOnYFace ( const int i, const int j, const int k, const int n,
                           const amrex::Real dt, const amrex::Real dy,
                           amrex::Real& Im, amrex::Real& Ip,
                           const amrex::Array4<const amrex::Real> &S,
                           const amrex::Real& vmac,
                           const amrex::BCRec bc,
                           const int domain_jlo, const int domain_jhi,
                           const bool is_velocity )
{
    using namespace amrex;
    {
        bool extdir_or_ho_jlo = (bc.lo(1) == BCType::ext_dir) or
                                (bc.lo(1) == BCType::hoextrap);
        bool extdir_or_ho_jhi = (bc.hi(1) == BCType::ext_dir) or
                                (bc.hi(1) == BCType::hoextrap);

        Real vpls, vmns;

        if (j == domain_jlo && (bc.lo(1) == BCType::ext_dir))
        {
            vmns = S(i,j-1,k,n);
            if ( n==YVEL and is_velocity )
            {
                vpls = S(i,j-1,k,n);
            }
            else
            {
                vpls = S(i,j  ,k,n) + 0.5 * (-1.0 - vmac * dt/dy) *
                    iamr_ho_yslope_extdir(i,j  ,k,n,S, extdir_or_ho_jlo, extdir_or_ho_jhi, domain_jlo, domain_jhi);
            }
        }
        else if (j == domain_jhi+1 && (bc.hi(1) == BCType::ext_dir))
        {
            vpls = S(i,j  ,k,n);
            if ( n==YVEL and is_velocity )
            {
                vmns = S(i,j  ,k,n);
            }
            else
            {
                vmns = S(i,j-1,k,n) + 0.5 * ( 1.0 - vmac * dt/dy) *
                    iamr_ho_yslope_extdir(i,j-1,k,n,S, extdir_or_ho_jlo, extdir_or_ho_jhi, domain_jlo, domain_jhi);
            }
        }
        else
        {
	    // Note that we still call the "extdir version" here because interior cells one
	    // away from the boundary will still see the boundary condition in the 4th order
	    // slope
            vpls = S(i,j  ,k,n) + 0.5 * (-1.0 - vmac * dt/dy) *
                iamr_ho_yslope_extdir(i,j  ,k,n,S, extdir_or_ho_jlo, extdir_or_ho_jhi, domain_jlo, domain_jhi);
            vmns = S(i,j-1,k,n) + 0.5 * ( 1.0 - vmac * dt/dy) *
                iamr_ho_yslope_extdir(i,j-1,k,n,S, extdir_or_ho_jlo, extdir_or_ho_jhi, domain_jlo, domain_jhi);
        }

        Ip = vmns;
        Im = vpls;
    }
}

#if (AMREX_SPACEDIM==3)

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void PredictStateOnZFace ( const int i, const int j, const int k, const int n,
                           const amrex::Real dt, const amrex::Real dz,
                           amrex::Real& Im, amrex::Real& Ip,
                           const amrex::Array4<const amrex::Real> &S,
                           const amrex::Real& wmac,
                           const amrex::BCRec bc,
                           const int domain_klo, const int domain_khi,
                           const bool is_velocity )
{
    using namespace amrex;
    {
        bool extdir_or_ho_klo = (bc.lo(2) == BCType::ext_dir) or
                                (bc.lo(2) == BCType::hoextrap);
        bool extdir_or_ho_khi = (bc.hi(2) == BCType::ext_dir) or
                                (bc.hi(2) == BCType::hoextrap);

        Real wpls, wmns;

        if (k == domain_klo && (bc.lo(2) == BCType::ext_dir))
        {
            wmns = S(i,j,k-1,n);
            if ( n == ZVEL and is_velocity )
            {
                wpls = S(i,j,k-1,n);
            }
            else
            {
                wpls = S(i,j,k  ,n) + 0.5 * (-1.0 - wmac * dt/dz) *
                    iamr_ho_zslope_extdir(i,j,k  ,n,S, extdir_or_ho_klo, extdir_or_ho_khi, domain_klo, domain_khi);
            }
        }
        else if (k == domain_khi+1 && (bc.hi(2) == BCType::ext_dir))
        {
            wpls = S(i,j,k  ,n);
            if ( n == ZVEL and is_velocity )
            {
                wmns = S(i,j,k  ,n);
            }
            else
            {
                wmns = S(i,j,k-1,n) + 0.5 * ( 1.0 - wmac * dt/dz) *
                    iamr_ho_zslope_extdir(i,j,k-1,n,S, extdir_or_ho_klo, extdir_or_ho_khi, domain_klo, domain_khi);
            }
        }
        else
        {
	    // Note that we still call the "extdir version" here because interior cells one
	    // away from the boundary will still see the boundary condition in the 4th order
	    // slope
            wpls = S(i,j,k  ,n) + 0.5 * (-1.0 - wmac * dt/dz) *
                iamr_ho_zslope_extdir(i,j,k  ,n,S, extdir_or_ho_klo, extdir_or_ho_khi, domain_klo, domain_khi);
            wmns = S(i,j,k-1,n) + 0.5 * ( 1.0 - wmac * dt/dz) *
                iamr_ho_zslope_extdir(i,j,k-1,n,S, extdir_or_ho_klo, extdir_or_ho_khi, domain_klo, domain_khi);
        }

        Ip = wmns;
        Im = wpls;
    }
}
#endif

}
#endif
