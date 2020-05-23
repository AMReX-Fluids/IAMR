#ifndef IAMR_MOL_H_
#define IAMR_MOL_H_

#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

namespace MOL {

    void ComputeAofs ( amrex::MultiFab& aofs, int aofs_comp, int ncomp,
                       amrex::MultiFab const& state, int state_comp,
                       D_DECL( amrex::MultiFab const& umac,
                               amrex::MultiFab const& vmac,
                               amrex::MultiFab const& wmac),
                       D_DECL( amrex::MultiFab& xedge,
                               amrex::MultiFab& yedge,
                               amrex::MultiFab& zedge),
                       int  edge_comp,
                       bool known_edgestate,
                       D_DECL( amrex::MultiFab& xfluxes,
                               amrex::MultiFab& yfluxes,
                               amrex::MultiFab& zfluxes),
                       int fluxes_comp,
                       amrex::Vector<amrex::BCRec> const& bcs,
                       amrex::Geometry const& geom);

    void ComputeSyncAofs ( amrex::MultiFab& aofs, int aofs_comp, int ncomp,
                           amrex::MultiFab const& state, int state_comp,
                           D_DECL( amrex::MultiFab const& umac,
                                   amrex::MultiFab const& vmac,
                                   amrex::MultiFab const& wmac),
                           D_DECL( amrex::MultiFab const& ucorr,
                                   amrex::MultiFab const& vcorr,
                                   amrex::MultiFab const& wcorr),
                           D_DECL( amrex::MultiFab& xedge,
                                   amrex::MultiFab& yedge,
                                   amrex::MultiFab& zedge),
                           int  edge_comp,
                           bool known_edgestate,
                           D_DECL( amrex::MultiFab& xfluxes,
                                   amrex::MultiFab& yfluxes,
                                   amrex::MultiFab& zfluxes),
                           int fluxes_comp,
                           amrex::Vector<amrex::BCRec> const& bcs,
                           amrex::Geometry const& geom);


    void ComputeEdgeState (amrex::Box const& a_bx,
                           D_DECL( amrex::Array4<amrex::Real> const& xedge,
                                   amrex::Array4<amrex::Real> const& yedge,
                                   amrex::Array4<amrex::Real> const& zedge),
                           amrex::Array4<amrex::Real const> const& q,
                           const int ncomp,
                           D_DECL( amrex::Array4<amrex::Real const> const& umac,
                                   amrex::Array4<amrex::Real const> const& vmac,
                                   amrex::Array4<amrex::Real const> const& wmac),
                           amrex::Box const&       domain,
                           amrex::Vector<amrex::BCRec> const& bcs);

    void ComputeFluxes ( amrex::Box const& bx,
                         D_DECL( amrex::Array4<amrex::Real> const& fx,
                                 amrex::Array4<amrex::Real> const& fy,
                                 amrex::Array4<amrex::Real> const& fz),
                         D_DECL( amrex::Array4<amrex::Real const> const& umac,
                                 amrex::Array4<amrex::Real const> const& vmac,
                                 amrex::Array4<amrex::Real const> const& wmac),
                         D_DECL( amrex::Array4<amrex::Real const> const& xedge,
                                 amrex::Array4<amrex::Real const> const& yedge,
                                 amrex::Array4<amrex::Real const> const& zedge),
                         int ncomp );

    void ComputeDivergence ( amrex::Box const& bx,
                             amrex::Array4<amrex::Real> const& div,
                             D_DECL( amrex::Array4<amrex::Real const> const& fx,
                                     amrex::Array4<amrex::Real const> const& fy,
                                     amrex::Array4<amrex::Real const> const& fz),
                             int ncomp, amrex::Geometry const& geom );

    void ExtrapVelToFaces ( const amrex::MultiFab&  vel,
                            D_DECL(amrex::MultiFab& umac,
                                   amrex::MultiFab& vmac,
                                   amrex::MultiFab& wmac ),
                            const amrex::Geometry&  a_geom,
                            const amrex::Vector<amrex::BCRec>& a_bcs );

    void PredictVelOnFaces ( D_DECL( amrex::Box const& ubx,
                                     amrex::Box const& vbx,
                                     amrex::Box const& wbx ),
                             D_DECL( amrex::Array4<amrex::Real> const& u,
                                     amrex::Array4<amrex::Real> const& v,
                                     amrex::Array4<amrex::Real> const& w ),
                             amrex::Array4<amrex::Real const > const& vcc,
                             const amrex::Geometry&  geom,
                             const amrex::BCRec* bcs );

    void SetMacBCs ( amrex::Box const& domain,
                     D_DECL( amrex::Box const& ubx,
                             amrex::Box const& vbx,
                             amrex::Box const& wbx ),
                     D_DECL( amrex::Array4<amrex::Real> const& u,
                             amrex::Array4<amrex::Real> const& v,
                             amrex::Array4<amrex::Real> const& w ),
                     amrex::Array4<amrex::Real const > const& vel,
                     const amrex::Vector<amrex::BCRec>& bcs );


    void AreaWeightFluxes ( D_DECL( const amrex::Box& xbx,
				    const amrex::Box& ybx,
				    const amrex::Box& zbx),
                            D_DECL( amrex::Array4<amrex::Real> const& fx,
                                    amrex::Array4<amrex::Real> const& fy,
                                    amrex::Array4<amrex::Real> const& fz),
                            D_DECL( amrex::Array4<amrex::Real> const& fxtmp,
                                    amrex::Array4<amrex::Real> const& fytmp,
                                    amrex::Array4<amrex::Real> const& fztmp),
                            int ncomp,
                            const amrex::Real* dx );

#ifdef AMREX_USE_EB
    void EB_ComputeEdgeState ( amrex::Box const& bx,
                               D_DECL( amrex::Array4<amrex::Real> const& xedge,
                                       amrex::Array4<amrex::Real> const& yedge,
                                       amrex::Array4<amrex::Real> const& zedge),
                               amrex::Array4<amrex::Real const> const& q,
                               const int ncomp,
                               D_DECL( amrex::Array4<amrex::Real const> const& umac,
                                       amrex::Array4<amrex::Real const> const& vmac,
                                       amrex::Array4<amrex::Real const> const& wmac),
                               amrex::Box const&       domain,
                               amrex::Vector<amrex::BCRec> const& bcs,
                               D_DECL( amrex::Array4<amrex::Real const> const& fcx,
                                       amrex::Array4<amrex::Real const> const& fcy,
                                       amrex::Array4<amrex::Real const> const& fcz),
                               amrex::Array4<amrex::Real const> const& ccc,
                               amrex::Array4<amrex::EBCellFlag const> const& flag);

    void EB_ComputeFluxes ( amrex::Box const& bx,
                            D_DECL( amrex::Array4<amrex::Real> const& fx,
                                    amrex::Array4<amrex::Real> const& fy,
                                    amrex::Array4<amrex::Real> const& fz),
                            D_DECL( amrex::Array4<amrex::Real const> const& umac,
                                    amrex::Array4<amrex::Real const> const& vmac,
                                    amrex::Array4<amrex::Real const> const& wmac),
                            D_DECL( amrex::Array4<amrex::Real const> const& xedge,
                                    amrex::Array4<amrex::Real const> const& yedge,
                                    amrex::Array4<amrex::Real const> const& zedge),
                            int ncomp,
                            amrex::Array4<amrex::EBCellFlag const> const& flag);

    void EB_ComputeDivergence ( amrex::Box const& bx,
                                amrex::Array4<amrex::Real> const& div,
                                D_DECL( amrex::Array4<amrex::Real const> const& fx,
                                        amrex::Array4<amrex::Real const> const& fy,
                                        amrex::Array4<amrex::Real const> const& fz),
                                int ncomp, amrex::Geometry const& geom,
                                amrex::Array4<amrex::EBCellFlag const> const& flag,
                                amrex::Array4<amrex::Real const> const& vfrac,
                                D_DECL( amrex::Array4<amrex::Real const> const& apx,
                                        amrex::Array4<amrex::Real const> const& apy,
                                        amrex::Array4<amrex::Real const> const& apz ) );

    void EB_PredictVelOnFaces ( amrex::Box const& ccbx,
                                D_DECL( amrex::Box const& ubx,
                                        amrex::Box const& vbx,
                                        amrex::Box const& wbx ),
                                D_DECL( amrex::Array4<amrex::Real> const& u,
                                        amrex::Array4<amrex::Real> const& v,
                                        amrex::Array4<amrex::Real> const& w ),
                                amrex::Array4<amrex::Real const > const& vcc,
                                amrex::Array4<amrex::EBCellFlag const> const& flag,
                                D_DECL( amrex::Array4<amrex::Real const> const& fcx,
                                        amrex::Array4<amrex::Real const> const& fcy,
                                        amrex::Array4<amrex::Real const> const& fcz),
                                amrex::Array4<amrex::Real const> const& ccc,
                                const amrex::Geometry&  geom,
                                const amrex::BCRec* bcs );

    void Redistribute (  amrex::Box const& bx, int ncomp,
                         amrex::Array4<amrex::Real> const& div,
                         amrex::Array4<amrex::Real const> const& div_in,
			 amrex::Array4<amrex::Real> const& scratch,
                         amrex::Array4<amrex::EBCellFlag const> const& flag,
                         amrex::Array4<amrex::Real const> const& vfrac,
                         amrex::Geometry const& geom );

    void EB_AreaWeightFluxes ( D_DECL( const amrex::Box& xbx,
				       const amrex::Box& ybx,
				       const amrex::Box& zbx),
                               D_DECL( amrex::Array4<amrex::Real> const& fx,
                                       amrex::Array4<amrex::Real> const& fy,
                                       amrex::Array4<amrex::Real> const& fz),
                               D_DECL( amrex::Array4<amrex::Real> const& fxtmp,
                                       amrex::Array4<amrex::Real> const& fytmp,
                                       amrex::Array4<amrex::Real> const& fztmp),
                               D_DECL( amrex::Array4<amrex::Real const> const& apx,
                                       amrex::Array4<amrex::Real const> const& apy,
                                       amrex::Array4<amrex::Real const> const& apz ),
                               int ncomp,
                               const amrex::Real* dx );

#endif

}


#endif
