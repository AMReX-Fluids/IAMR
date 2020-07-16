#ifndef IAMR_CONVECTION_K_H_
#define IAMR_CONVECTION_K_H_

#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>
#else
#include <AMReX_FArrayBox.H>
#endif

namespace {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_ho_xslope (int i, int j, int k, int n,
                              amrex::Array4<amrex::Real const> const& q) noexcept
{
    amrex::Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp;
    amrex::Real qm, qp, qi;
    qi = q(i,j,k,n);
    qm = q(i-1,j,k,n);
    qp = q(i+1,j,k,n);

    dlft = qm - q(i-2,j,k,n);
    drgt = qi - qm;
    dcen = 0.5*(dlft+drgt);
    dsgn = std::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
    dfm = dsgn*amrex::min(dlim, std::abs(dcen));

    dlft = qp - qi;
    drgt = q(i+2,j,k,n) - qp;
    dcen = 0.5*(dlft+drgt);
    dsgn = std::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
    dfp = dsgn*amrex::min(dlim, std::abs(dcen));

    dlft = qi - qm;
    drgt = qp - qi;
    dcen = 0.5*(dlft+drgt);
    dsgn = std::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;

    dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);

    return dsgn*amrex::min(dlim, std::abs(dtemp));

}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_ho_xslope_extdir (int i, int j, int k, int n,
                                     amrex::Array4<amrex::Real const> const& q,
                                     bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp, dlimsh, dsgnsh;
    amrex::Real qm, qp, qi;
    qi = q(i,j,k,n);
    qm = q(i-1,j,k,n);
    qp = q(i+1,j,k,n);

    dlft = qm - q(i-2,j,k,n);
    drgt = qi - qm;
    dcen = 0.5*(dlft+drgt);
    dsgn = std::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
    dfm = dsgn*amrex::min(dlim, std::abs(dcen));

    dlft = qp - qi;
    drgt = q(i+2,j,k,n) - qp;
    dcen = 0.5*(dlft+drgt);
    dsgn = std::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
    dfp = dsgn*amrex::min(dlim, std::abs(dcen));

    dlft = qi - qm;
    drgt = qp - qi;
    dcen = 0.5*(dlft+drgt);
    dsgn = std::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;

    dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);

    if (edlo and i == domlo) {
       dtemp  = -16./15.*q(i-1,j,k,n) + .5*q(i,j,k,n) + 2./3.*q(i+1,j,k,n) -  0.1*q(i+2,j,k,n);
       dlft = 2.*(q(i  ,j,k,n)-q(i-1,j,k,n));
       drgt = 2.*(q(i+1,j,k,n)-q(i  ,j,k,n));
       dlim = (dlft*drgt >= 0.0) ? amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
       dsgn = std::copysign(1.e0, dtemp);
    } else if (edlo and i == domlo+1) {
       dfm  = -16./15.*q(domlo-1,j,k,n) + .5*q(domlo,j,k,n) + 2./3.*q(domlo+1,j,k,n) -  0.1*q(domlo+2,j,k,n);
       dlft = 2.*(q(domlo  ,j,k,n)-q(domlo-1,j,k,n));
       drgt = 2.*(q(domlo+1,j,k,n)-q(domlo  ,j,k,n));
       dlimsh = (dlft*drgt >= 0.0) ? amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
       dsgnsh = std::copysign(1.e0, dfm);
       dfm = dsgnsh*amrex::min(dlimsh, std::abs(dfm));
       dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);
    }

    if (edhi and i == domhi) {
       dtemp  = 16./15.*q(i+1,j,k,n) - .5*q(i,j,k,n) - 2./3.*q(i-1,j,k,n) +  0.1*q(i-2,j,k,n);
       dlft = 2.*(q(i  ,j,k,n)-q(i-1,j,k,n));
       drgt = 2.*(q(i+1,j,k,n)-q(i  ,j,k,n));
       dlim = (dlft*drgt >= 0.0) ? amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
       dsgn = std::copysign(1.e0, dtemp);
    } else if (edhi and i == domhi-1) {
       dfp  = 16./15.*q(domhi+1,j,k,n) - .5*q(domhi,j,k,n) - 2./3.*q(domhi-1,j,k,n) +  0.1*q(domhi-2,j,k,n);
       dlft = 2.*(q(domhi  ,j,k,n)-q(domhi-1,j,k,n));
       drgt = 2.*(q(domhi+1,j,k,n)-q(domhi  ,j,k,n));
       dlimsh = (dlft*drgt >= 0.0) ? amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
       dsgnsh = std::copysign(1.e0, dfp);
       dfp = dsgnsh*amrex::min(dlimsh, std::abs(dfp));
       dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);
    }

    return dsgn*amrex::min(dlim, std::abs(dtemp));
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_ho_yslope (int i, int j, int k, int n,
                              amrex::Array4<amrex::Real const> const& q) noexcept
{
    amrex::Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp;
    amrex::Real qm, qp, qj;
    qj = q(i,j,k,n);
    qm = q(i,j-1,k,n);
    qp = q(i,j+1,k,n);

    dlft = qm - q(i,j-2,k,n);
    drgt = qj - qm;
    dcen = 0.5*(dlft+drgt);
    dsgn = std::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
    dfm = dsgn*amrex::min(dlim, std::abs(dcen));

    dlft = qp - qj;
    drgt = q(i,j+2,k,n) - qp;
    dcen = 0.5*(dlft+drgt);
    dsgn = std::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
    dfp = dsgn*amrex::min(dlim, std::abs(dcen));

    dlft = qj - qm;
    drgt = qp - qj;
    dcen = 0.5*(dlft+drgt);
    dsgn = std::copysign(1.e0, dcen);

    dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);
    return dsgn*amrex::min(dlim, std::abs(dtemp));
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_ho_yslope_extdir (int i, int j, int k, int n,
                                     amrex::Array4<amrex::Real const> const& q,
                                     bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp, dlimsh,dsgnsh;
    amrex::Real qm, qp, qj;
    qj = q(i,j,k,n);
    qm = q(i,j-1,k,n);
    qp = q(i,j+1,k,n);

    dlft = qm - q(i,j-2,k,n);
    drgt = qj - qm;
    dcen = 0.5*(dlft+drgt);
    dsgn = std::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
    dfm = dsgn*amrex::min(dlim, std::abs(dcen));

    dlft = qp - qj;
    drgt = q(i,j+2,k,n) - qp;
    dcen = 0.5*(dlft+drgt);
    dsgn = std::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
    dfp = dsgn*amrex::min(dlim, std::abs(dcen));

    dlft = qj - qm;
    drgt = qp - qj;
    dcen = 0.5*(dlft+drgt);
    dsgn = std::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;

    dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);

    if (edlo and j == domlo) {
       dtemp  = -16./15.*q(i,j-1,k,n) + .5*q(i,j,k,n) + 2./3.*q(i,j+1,k,n) -  0.1*q(i,j+2,k,n);
       dlft = 2.*(q(i  ,j,k,n)-q(i,j-1,k,n));
       drgt = 2.*(q(i,j+1,k,n)-q(i  ,j,k,n));
       dlim = (dlft*drgt >= 0.0) ? amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
       dsgn = std::copysign(1.e0, dtemp);
    } else if (edlo and j == domlo+1) {
       dfm  = -16./15.*q(i,domlo-1,k,n) + .5*q(i,domlo,k,n) + 2./3.*q(i,domlo+1,k,n) -  0.1*q(i,domlo+2,k,n);
       dlft = 2.*(q(i  ,domlo,k,n)-q(i,domlo-1,k,n));
       drgt = 2.*(q(i,domlo+1,k,n)-q(i  ,domlo,k,n));
       dlimsh = (dlft*drgt >= 0.0) ? amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
       dsgnsh = std::copysign(1.e0, dfm);
       dfm = dsgnsh*amrex::min(dlimsh, std::abs(dfm));
       dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);
    }

    if (edhi and j == domhi) {
       dtemp  = 16./15.*q(i,j+1,k,n) - .5*q(i,j,k,n) - 2./3.*q(i,j-1,k,n) +  0.1*q(i,j-2,k,n);
       dlft = 2.*(q(i  ,j,k,n)-q(i,j-1,k,n));
       drgt = 2.*(q(i,j+1,k,n)-q(i  ,j,k,n));
       dlim = (dlft*drgt >= 0.0) ? amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
       dsgn = std::copysign(1.e0, dtemp);
    } else if (edhi and j == domhi-1) {
       dfp  = 16./15.*q(i,domhi+1,k,n) - .5*q(i,domhi,k,n) - 2./3.*q(i,domhi-1,k,n) +  0.1*q(i,domhi-2,k,n);
       dlft = 2.*(q(i  ,domhi,k,n)-q(i,domhi-1,k,n));
       drgt = 2.*(q(i,domhi+1,k,n)-q(i  ,domhi,k,n));
       dlimsh = (dlft*drgt >= 0.0) ? amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
       dsgnsh = std::copysign(1.e0, dfp);
       dfp = dsgnsh*amrex::min(dlimsh, std::abs(dfp));
       dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);
    }

    return dsgn*amrex::min(dlim, std::abs(dtemp));

}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_ho_zslope (int i, int j, int k, int n,
                              amrex::Array4<amrex::Real const> const& q) noexcept
{
    amrex::Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp;
    amrex::Real qm, qp, qk;
    qk = q(i,j,k,n);
    qm = q(i,j,k-1,n);
    qp = q(i,j,k+1,n);

    dlft = qm - q(i,j,k-2,n);
    drgt = qk - qm;
    dcen = 0.5*(dlft+drgt);
    dsgn = std::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
    dfm = dsgn*amrex::min(dlim, std::abs(dcen));

    dlft = qp - qk;
    drgt = q(i,j,k+2,n) - qp;
    dcen = 0.5*(dlft+drgt);
    dsgn = std::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
    dfp = dsgn*amrex::min(dlim, std::abs(dcen));

    dlft = qk - qm;
    drgt = qp - qk;
    dcen = 0.5*(dlft+drgt);
    dsgn = std::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;

    dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);
    return dsgn*amrex::min(dlim, std::abs(dtemp));

}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real iamr_ho_zslope_extdir (int i, int j, int k, int n,
                                     amrex::Array4<amrex::Real const> const& q,
                                     bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp, dlimsh, dsgnsh;
    amrex::Real qm, qp, qk;
    qk = q(i,j,k,n);
    qm = q(i,j,k-1,n);
    qp = q(i,j,k+1,n);

    dlft = qm - q(i,j,k-2,n);
    drgt = qk - qm;
    dcen = 0.5*(dlft+drgt);
    dsgn = std::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
    dfm = dsgn*amrex::min(dlim, std::abs(dcen));

    dlft = qp - qk;
    drgt = q(i,j,k+2,n) - qp;
    dcen = 0.5*(dlft+drgt);
    dsgn = std::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
    dfp = dsgn*amrex::min(dlim, std::abs(dcen));

    dlft = qk - qm;
    drgt = qp - qk;
    dcen = 0.5*(dlft+drgt);
    dsgn = std::copysign(1.e0, dcen);
    dlim = (dlft*drgt >= 0.0) ? 2.0*amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;

    dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);

    if (edlo and k == domlo) {
       dtemp  = -16./15.*q(i,j,k-1,n) + .5*q(i,j,k,n) + 2./3.*q(i,j,k+1,n) -  0.1*q(i,j,k+2,n);
       dlft = 2.*(q(i  ,j,k,n)-q(i,j,k-1,n));
       drgt = 2.*(q(i,j,k+1,n)-q(i  ,j,k,n));
       dlim = (dlft*drgt >= 0.0) ? amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
       dsgn = std::copysign(1.e0, dtemp);
    } else if (edlo and k == domlo+1) {
       dfm  = -16./15.*q(i,j,domlo-1,n) + .5*q(i,j,domlo,n) + 2./3.*q(i,j,domlo+1,n) -  0.1*q(i,j,domlo+2,n);
       dlft = 2.*(q(i  ,j,domlo,n)-q(i,j,domlo-1,n));
       drgt = 2.*(q(i,j,domlo+1,n)-q(i  ,j,domlo,n));
       dlimsh = (dlft*drgt >= 0.0) ? amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
       dsgnsh = std::copysign(1.e0, dfm);
       dfm = dsgnsh*amrex::min(dlimsh, std::abs(dfm));
       dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);
    }

    if (edhi and k == domhi) {
       dtemp  = 16./15.*q(i,j,k+1,n) - .5*q(i,j,k,n) - 2./3.*q(i,j,k-1,n) +  0.1*q(i,j,k-2,n);
       dlft = 2.*(q(i  ,j,k,n)-q(i,j,k-1,n));
       drgt = 2.*(q(i,j,k+1,n)-q(i  ,j,k,n));
       dlim = (dlft*drgt >= 0.0) ? amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
       dsgn = std::copysign(1.e0, dtemp);
    } else if (edhi and k == domhi-1) {
       dfp  = 16./15.*q(i,j,domhi+1,n) - .5*q(i,j,domhi,n) - 2./3.*q(i,j,domhi-1,n) +  0.1*q(i,j,domhi-2,n);
       dlft = 2.*(q(i  ,j,domhi,n)-q(i,j,domhi-1,n));
       drgt = 2.*(q(i,j,domhi+1,n)-q(i  ,j,domhi,n));
       dlimsh = (dlft*drgt >= 0.0) ? amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0;
       dsgnsh = std::copysign(1.e0, dfp);
       dfp = dsgnsh*amrex::min(dlimsh, std::abs(dfp));
       dtemp  = 4.0/3.0*dcen - 1.0/6.0*(dfp + dfm);
    }

    return dsgn*amrex::min(dlim, std::abs(dtemp));
}

}
#endif
