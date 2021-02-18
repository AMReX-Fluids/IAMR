#ifndef IAMR_PLM_EBGODUNOV_H
#define IAMR_PLM_EBGODUNOV_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiCutFab.H>

// #include <iamr_slopes_godunov_K.H>
// #include <AMReX_Gpu.H>
// #include <AMReX_FArrayBox.H>
// #include <AMReX_BCRec.H>
// #include <AMReX_BC_TYPES.H>
// #include <AMReX_Array.H>
// #include <iomanip>
// #include <iamr_constants_godunov.H>

namespace EBPLM {

void PredictVelOnXFace ( amrex::Box const& bx,
                         amrex::Array4<amrex::Real> const& Imx,
                         amrex::Array4<amrex::Real> const& Ipx,
                         amrex::Array4<amrex::Real const> const& q,
                         amrex::Array4<amrex::Real  const> const& vcc,
                         amrex::Array4<amrex::EBCellFlag const> const& flag,
                         amrex::Array4<amrex::Real const> const& vfrac,
                         AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                      amrex::Array4<amrex::Real const> const& fcy,
                                      amrex::Array4<amrex::Real const> const& fcz),
                         amrex::Array4<amrex::Real const> const& ccc,
                         amrex::Geometry& geom,
                         amrex::Real dt,
                         amrex::Vector<amrex::BCRec> const& h_bcrec,
                         amrex::BCRec const* d_bcrec);

void PredictVelOnYFace ( amrex::Box const& bx,
                         amrex::Array4<amrex::Real> const& Imy,
                         amrex::Array4<amrex::Real> const& Ipy,
                         amrex::Array4<amrex::Real const> const& q,
                         amrex::Array4<amrex::Real  const> const& vcc,
                         amrex::Array4<amrex::EBCellFlag const> const& flag,
                         amrex::Array4<amrex::Real const> const& vfrac,
                         AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                      amrex::Array4<amrex::Real const> const& fcy,
                                      amrex::Array4<amrex::Real const> const& fcz),
                         amrex::Array4<amrex::Real const> const& ccc,
                         amrex::Geometry& geom,
                         amrex::Real dt,
                         amrex::Vector<amrex::BCRec> const& h_bcrec,
                         amrex::BCRec const* d_bcrec);

#if (AMREX_SPACEDIM==3)
void PredictVelOnZFace ( amrex::Box const& bx,
                         amrex::Array4<amrex::Real> const& Imz,
                         amrex::Array4<amrex::Real> const& Ipz,
                         amrex::Array4<amrex::Real const> const& q,
                         amrex::Array4<amrex::Real  const> const& vcc,
                         amrex::Array4<amrex::EBCellFlag const> const& flag,
                         amrex::Array4<amrex::Real const> const& vfrac,
                         AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                      amrex::Array4<amrex::Real const> const& fcy,
                                      amrex::Array4<amrex::Real const> const& fcz),
                         amrex::Array4<amrex::Real const> const& ccc,
                         amrex::Geometry& geom,
                         amrex::Real dt,
                         amrex::Vector<amrex::BCRec> const& h_bcrec,
                         amrex::BCRec const* d_bcrec);
#endif



void PredictStateOnXFace ( amrex::Box const& bx, int ncomp,
                           amrex::Array4<amrex::Real> const& Imx, amrex::Array4<amrex::Real> const& Ipx,
                           amrex::Array4<amrex::Real const> const& q,
                           amrex::Array4<amrex::Real  const> const& umac,
                           amrex::Array4<amrex::EBCellFlag const> const& flag,
                           amrex::Array4<amrex::Real const> const& vfrac,
                           AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                        amrex::Array4<amrex::Real const> const& fcy,
                                        amrex::Array4<amrex::Real const> const& fcz),
                           amrex::Array4<amrex::Real const> const& ccc,
                           amrex::Geometry const& geom,
                           amrex::Real dt,
                           amrex::Vector<amrex::BCRec> const& h_bcrec,
                           amrex::BCRec const* d_bcrec, bool is_velocity);

void PredictStateOnYFace ( amrex::Box const& bx, int ncomp,
                           amrex::Array4<amrex::Real> const& Imy, amrex::Array4<amrex::Real> const& Ipy,
                           amrex::Array4<amrex::Real const> const& q,
                           amrex::Array4<amrex::Real  const> const& vmac,
                           amrex::Array4<amrex::EBCellFlag const> const& flag,
                           amrex::Array4<amrex::Real const> const& vfrac,
                           AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                        amrex::Array4<amrex::Real const> const& fcy,
                                        amrex::Array4<amrex::Real const> const& fcz),
                           amrex::Array4<amrex::Real const> const& ccc,
                           amrex::Geometry const& geom,
                           amrex::Real dt,
                           amrex::Vector<amrex::BCRec> const& h_bcrec,
                           amrex::BCRec const* d_bcrec, bool is_velocity);

#if (AMREX_SPACEDIM == 3)
void PredictStateOnZFace ( amrex::Box const& bx, int ncomp,
                           amrex::Array4<amrex::Real> const& Imz, amrex::Array4<amrex::Real> const& Ipz,
                           amrex::Array4<amrex::Real const> const& q,
                           amrex::Array4<amrex::Real  const> const& wmac,
                           amrex::Array4<amrex::EBCellFlag const> const& flag,
                           amrex::Array4<amrex::Real const> const& vfrac,
                           AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                        amrex::Array4<amrex::Real const> const& fcy,
                                        amrex::Array4<amrex::Real const> const& fcz),
                           amrex::Array4<amrex::Real const> const& ccc,
                           amrex::Geometry const& geom,
                           amrex::Real dt,
                           amrex::Vector<amrex::BCRec> const& h_bcrec,
                           amrex::BCRec const* d_bcrec, bool is_velocity);
#endif

}
#endif
