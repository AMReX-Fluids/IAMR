
#ifndef IAMR_MacProj_H_
#define IAMR_MacProj_H_

#include <AMReX_Amr.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_MLLinOp.H>
#include <AMReX_YAFluxRegister.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBFluxRegister.H>
#endif

//
// A useful enumeration of the forms for the advective terms
//
enum AdvectionForm {Conservative = 0, NonConservative};

class MacProj
{
public:

    MacProj (amrex::Amr*   parent,
             int    _finest_level,
             amrex::BCRec* _phys_bc,
             int    /*not used*/ );

    void install_level (int       level,
                        amrex::AmrLevel* level_data);

    void setup (int level);
    void cleanup (int level);
    //
    // The level advance mac projection.
    //
    void mac_project (int             level,
                      amrex::MultiFab*       u_mac,
                      amrex::MultiFab&       S,
                      amrex::Real            dt,
                      amrex::Real            prev_time,
                      const amrex::MultiFab& divu,
                      int             have_divu,
                      const amrex::BCRec&   density_math_bc,
                      bool            increment_vel_register = true );

    //
    // The sync solve.
    //
    void mac_sync_solve (int       level,
                         amrex::Real      dt,
                         amrex::MultiFab& rho_half,
                         const amrex::BCRec&  rho_math_bc,
                         amrex::IntVect&  fine_ratio,
                         amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>& Ucorr,
                         amrex::MultiFab* Rhs_increment = nullptr);

    //
    // Computing the sync tendency.
    //
    void mac_sync_compute (int                          level,
                           amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>& Ucorr,
                           amrex::MultiFab&             Vsync,
                           amrex::MultiFab&             Ssync,
                           amrex::Vector<AdvectionForm>& advectionType,
                           amrex::Real                  prev_time,
                           amrex::Real                  dt,
                           int                          num_state_comps,
                           amrex::Real                  be_cn_theta,
                           int                          do_mom_diff,
                           bool                         update_fluxreg = true);

    void mac_sync_compute (int                          level,
                           amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>& Ucorr,
                           amrex::MultiFab&             Sync,
                           int                          comp,
                           int                          Sync_indx,
                           amrex::MultiFab* const*      edgestate,
                           int                          edge_comp,
                           amrex::Real                  dt,
                           bool                         update_fluxreg = true);

    static void mlmg_mac_solve (amrex::Amr* parent, const amrex::MultiFab* cphi,
                         const amrex::BCRec& phys_bc,
                         const amrex::BCRec& density_math_bc,
                         int level, amrex::Real mac_tol, amrex::Real mac_abs_tol,
                         amrex::Real rhs_scale,
                         const amrex::MultiFab &rho, const amrex::MultiFab &Rhs,
                         amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>& u_mac,
                         amrex::MultiFab *mac_phi,
                         amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>& fluxes);

    static void set_mac_solve_bc (amrex::Array<amrex::MLLinOp::BCType,AMREX_SPACEDIM>& mlmg_lobc,
                           amrex::Array<amrex::MLLinOp::BCType,AMREX_SPACEDIM>& mlmg_hibc,
                           const amrex::BCRec& phys_bc, const amrex::Geometry& geom);

    static void Initialize ();
    static void Finalize ();

    void test_umac_periodic (int level, amrex::MultiFab* u_mac) const;

    //
    // Test the divergence constraint.
    //
    void check_div_cond (int      level,
                         amrex::MultiFab U_edge[]) const;
    //
    // Boundary conditions.
    //
    void set_outflow_bcs (int level,
                          amrex::MultiFab* mac_phi,
                          const amrex::MultiFab* u_mac,
                          const amrex::MultiFab& S,
                          const  amrex::MultiFab& divu);
    //
    // Pointers to amr,amrlevel.
    //
    amrex::Amr*             parent;
    amrex::Vector<amrex::AmrLevel*> LevelData;
    //
    // Boundary condition objects.
    //
    amrex::BCRec*                phys_bc;
    //
    // MAC sync correction and solution.
    //
    amrex::Vector<std::unique_ptr<amrex::MultiFab    > > mac_phi_crse;
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > mac_reg;
    //
    // Parameters.
    //
    int        finest_level;
    int        finest_level_allocated;

    static int  verbose;
    static int  do_outflow_bcs;
    static amrex::Real mac_tol;
    static amrex::Real mac_abs_tol;
    static amrex::Real mac_sync_tol;
    static int  check_umac_periodicity;
    static int fix_mac_sync_rhs;

    //
    // Options for MacProjector
    // NOTE: IAMR uses a different max_order default than amrex::MacProjector,
    // which uses a default of 3
    //
    static int max_order;
    static int agglomeration;
    static int consolidation;
    static int max_fmg_iter;

};

#endif
