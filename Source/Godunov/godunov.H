//
//  Godunov.H
//

#ifndef godunov_H
#define godunov_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_BCRec.H>

namespace godunov {

    void ExtrapVelToFaces ( amrex::MultiFab const& vel,
                            amrex::MultiFab const& vel_forces,
                            AMREX_D_DECL( amrex::MultiFab& u_mac,
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
                      D_DECL(amrex::Array4<amrex::Real> const& Imx,
                             amrex::Array4<amrex::Real> const& Imy,
                             amrex::Array4<amrex::Real> const& Imz),
                      D_DECL(amrex::Array4<amrex::Real> const& Ipx,
                             amrex::Array4<amrex::Real> const& Ipy,
                             amrex::Array4<amrex::Real> const& Ipz),
                      amrex::Array4<amrex::Real const> const& q,
                      amrex::Array4<amrex::Real const> const& vel,
                      amrex::Geometry geom,
                      amrex::Real dt,
                      amrex::BCRec const* d_bcrec);

    void make_trans_velocities (D_DECL(amrex::Box const& xbx,
                                       amrex::Box const& ybx,
                                       amrex::Box const& zbx),
                                D_DECL(amrex::Array4<amrex::Real> const& u_ad,
                                       amrex::Array4<amrex::Real> const& v_ad,
                                       amrex::Array4<amrex::Real> const& w_ad),
                                D_DECL(amrex::Array4<amrex::Real const> const& Imx,
                                       amrex::Array4<amrex::Real const> const& Imy,
                                       amrex::Array4<amrex::Real const> const& Imz),
                                D_DECL(amrex::Array4<amrex::Real const> const& Ipx,
                                       amrex::Array4<amrex::Real const> const& Ipy,
                                       amrex::Array4<amrex::Real const> const& Ipz),
                                amrex::Array4<amrex::Real const> const& vel,
                                amrex::Array4<amrex::Real const> const& f,
                                const amrex::Box& domain,
                                amrex::Real dt,
                                amrex::BCRec const* d_bcrec,
                                bool use_forces_in_trans);

    void predict_godunov_on_box (amrex::Box const& bx, int ncomp,
                                 D_DECL(amrex::Box const& xbx,
                                        amrex::Box const& ybx,
                                        amrex::Box const& zbx),
                                 D_DECL(amrex::Array4<amrex::Real> const& qx,
                                        amrex::Array4<amrex::Real> const& qy,
                                        amrex::Array4<amrex::Real> const& qz),
                                 amrex::Array4<amrex::Real const> const& q,
                                 D_DECL(amrex::Array4<amrex::Real const> const& u_ad,
                                        amrex::Array4<amrex::Real const> const& v_ad,
                                        amrex::Array4<amrex::Real const> const& w_ad),
                                 D_DECL(amrex::Array4<amrex::Real> const& Imx,
                                        amrex::Array4<amrex::Real> const& Imy,
                                        amrex::Array4<amrex::Real> const& Imz),
                                 D_DECL(amrex::Array4<amrex::Real> const& Ipx,
                                        amrex::Array4<amrex::Real> const& Ipy,
                                        amrex::Array4<amrex::Real> const& Ipz),
                                 amrex::Array4<amrex::Real const> const& f,
                                 const amrex::Box& domain,
                                 const amrex::Real* dx,
                                 amrex::Real dt,
                                 amrex::BCRec const* d_bcrec,
                                 bool use_forces_in_trans,
                                 amrex::Real* p);

    void compute_godunov_advection (amrex::Box const& bx, int ncomp,
                                    amrex::Array4<amrex::Real> const& dqdt,
                                    amrex::Array4<amrex::Real const> const& q,
                                    AMREX_D_DECL(amrex::Array4<amrex::Real> const& q_x,
                                                 amrex::Array4<amrex::Real> const& q_y,
                                                 amrex::Array4<amrex::Real> const& q_z),
                                    AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                                 amrex::Array4<amrex::Real const> const& vmac,
                                                 amrex::Array4<amrex::Real const> const& wmac),
                                    amrex::Array4<amrex::Real const> const& divu,
                                    amrex::Array4<amrex::Real const> const& fq,
                                    amrex::Geometry geom,
                                    amrex::Real dt,
                                    amrex::BCRec const* d_bcrec,
                                    int const* iconserv,
                                    amrex::Real* p, bool use_ppm, bool is_velocity,
                                    bool use_forces_in_trans);

} // namespace godunov

#endif /* Godunov_H */
