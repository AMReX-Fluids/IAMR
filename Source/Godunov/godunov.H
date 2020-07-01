//
//  Godunov.H
//

#ifndef godunov_H
#define godunov_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_BCRec.H>

namespace godunov {

    constexpr amrex::Real small_vel = 1e-8;

    void predict_godunov ( amrex::MultiFab const& vel,
                           amrex::MultiFab const& vel_forces,
                           D_DECL( amrex::MultiFab& u_mac,
                                   amrex::MultiFab& v_mac,
                                   amrex::MultiFab& w_mac ),
                           const amrex::Vector<amrex::BCRec> & h_bcrec,
                           const amrex::Geometry& geom, amrex::Real l_dt,
                           bool use_ppm, bool use_forces_in_trans);


    void predict_plm_x (amrex::Box const& bx, int ncomp,
                        amrex::Array4<amrex::Real> const& Imx, amrex::Array4<amrex::Real> const& Ipx,
                        amrex::Array4<amrex::Real const> const& q,
                        amrex::Array4<amrex::Real  const> const& vcc,
                        const amrex::Geometry& geom,
                        amrex::Real dt,
                        amrex::Vector<amrex::BCRec> const& h_bcrec,
                        amrex::BCRec const* d_bcrec);

    void predict_plm_y (amrex::Box const& bx, int ncomp,
                        amrex::Array4<amrex::Real> const& Imy, amrex::Array4<amrex::Real> const& Ipy,
                        amrex::Array4<amrex::Real const> const& q,
                        amrex::Array4<amrex::Real  const> const& vcc,
                        const amrex::Geometry& geom,
                        amrex::Real dt,
                        amrex::Vector<amrex::BCRec> const& h_bcrec,
                        amrex::BCRec const* d_bcrec);

    void predict_plm_z (amrex::Box const& bx, int ncomp,
                        amrex::Array4<amrex::Real> const& Imz, amrex::Array4<amrex::Real> const& Ipz,
                        amrex::Array4<amrex::Real const> const& q,
                        amrex::Array4<amrex::Real  const> const& vcc,
                        const amrex::Geometry& geom,
                        amrex::Real dt,
                        amrex::Vector<amrex::BCRec> const& h_bcrec,
                        amrex::BCRec const* d_bcrec);

    void predict_ppm (amrex::Box const& bx, int ncomp,
                      amrex::Array4<amrex::Real> const& Imx,
                      amrex::Array4<amrex::Real> const& Ipx,
                      amrex::Array4<amrex::Real> const& Imy,
                      amrex::Array4<amrex::Real> const& Ipy,
                      amrex::Array4<amrex::Real> const& Imz,
                      amrex::Array4<amrex::Real> const& Ipz,
                      amrex::Array4<amrex::Real const> const& q,
                      amrex::Array4<amrex::Real const> const& vel,
                      amrex::Geometry geom,
                      amrex::Real dt,
                      amrex::BCRec const* d_bcrec);

    void make_trans_velocities (amrex::Box const& xbx,
                                amrex::Box const& ybx, amrex::Box const& zbx,
                                amrex::Array4<amrex::Real> const& u_ad,
                                amrex::Array4<amrex::Real> const& v_ad,
                                amrex::Array4<amrex::Real> const& w_ad,
                                amrex::Array4<amrex::Real const> const& Imx,
                                amrex::Array4<amrex::Real const> const& Ipx,
                                amrex::Array4<amrex::Real const> const& Imy,
                                amrex::Array4<amrex::Real const> const& Ipy,
                                amrex::Array4<amrex::Real const> const& Imz,
                                amrex::Array4<amrex::Real const> const& Ipz,
                                amrex::Array4<amrex::Real const> const& vel,
                                amrex::Array4<amrex::Real const> const& f,
                                const amrex::Box& domain,
                                amrex::Real dt,
                                amrex::BCRec const* d_bcrec,
                                bool use_forces_in_trans);

    void predict_godunov_on_box (amrex::Box const& bx, int ncomp,
                                 amrex::Box const& xbx, amrex::Box const& ybx, amrex::Box const& zbx,
                                 amrex::Array4<amrex::Real> const& qx,
                                 amrex::Array4<amrex::Real> const& qy,
                                 amrex::Array4<amrex::Real> const& qz,
                                 amrex::Array4<amrex::Real const> const& q,
                                 amrex::Array4<amrex::Real const> const& u_ad,
                                 amrex::Array4<amrex::Real const> const& v_ad,
                                 amrex::Array4<amrex::Real const> const& w_ad,
                                 amrex::Array4<amrex::Real> const& Imx,
                                 amrex::Array4<amrex::Real> const& Ipx,
                                 amrex::Array4<amrex::Real> const& Imy,
                                 amrex::Array4<amrex::Real> const& Ipy,
                                 amrex::Array4<amrex::Real> const& Imz,
                                 amrex::Array4<amrex::Real> const& Ipz,
                                 amrex::Array4<amrex::Real const> const& f,
                                 const amrex::Box& domain,
                                 const amrex::Real* dx,
                                 amrex::Real dt,
                                 amrex::BCRec const* d_bcrec,
                                 bool use_forces_in_trans,
                                 amrex::Real* p);

    void compute_godunov_advection (int lev, amrex::Box const& bx, int ncomp,
                                    amrex::Array4<amrex::Real> const& dqdt,
                                    amrex::Array4<amrex::Real const> const& q,
                                    amrex::Array4<amrex::Real const> const& umac,
                                    amrex::Array4<amrex::Real const> const& vmac,
                                    amrex::Array4<amrex::Real const> const& wmac,
                                    amrex::Array4<amrex::Real const> const& fq,
                                    amrex::Vector<amrex::Geometry> geom,
                                    amrex::Real dt,
                                    amrex::BCRec const* d_bcrec,
                                    int const* iconserv,
                                    amrex::Real* p, bool use_ppm);

} // namespace godunov

#endif /* Godunov_H */
