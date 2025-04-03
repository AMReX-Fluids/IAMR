
#ifndef IAMR_NavierStokesBase_H_
#define IAMR_NavierStokesBase_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrLevel.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_BLFort.H>
#include <Diffusion.H>
#include <AMReX_ErrorList.H>
#include <MacProj.H>
#include <Projection.H>
#include <SyncRegister.H>
#include <AMReX_Utility.H>

#ifdef AMREX_PARTICLES
#include <AMReX_AmrParticles.H>
#endif

#ifdef AMREX_USE_EB
#include <AMReX_EBCellFlag.H>
#endif

//
// Parameters used when we call BuildMask
//
// interior  : interior cells (i.e., valid cells)
// covered   : ghost cells covered by valid cells of this FabArray
//             (including periodically shifted valid cells)
// notcovered: ghost cells not covered by valid cells
//             (including ghost cells outside periodic boundaries where the
//             periodically shifted cells don't exist at this level)
// physbnd   : boundary cells outside the domain (excluding periodic boundaries)
static constexpr int level_mask_interior   = 0; // valid cells
static constexpr int level_mask_covered    = 1; // ghost cells covered by valid cells of this level
static constexpr int level_mask_notcovered = 2; // ghost cells not covered
static constexpr int level_mask_physbnd    = 3; // outside domain

//
// "Divu_Type" means S, where divergence U = S
// "Dsdt_Type" means pd S/pd t, where S is as above
//
// Determine what you want in the state -- Divu, Dsdt -- in
// NavierStokes::variableSetUp in NS_setup.cpp.
//
enum StateType {State_Type=0, Press_Type, Gradp_Type};

#if (AMREX_SPACEDIM == 2)
enum StateNames  { Xvel=0, Yvel, Density};
#else
enum StateNames  { Xvel=0, Yvel, Zvel, Density};
#endif

enum PressureNames { Pressure=0 };
#if (AMREX_SPACEDIM == 2)
enum GradpNames  { Gradpx=0, Gradpy};
#else
enum GradpNames  { Gradpx=0, Gradpy, Gradpz};
#endif
enum DivuNames { Divu=0};
enum DsDtNames { Dsdt=0};

class NavierStokesBase
    :
    public amrex::AmrLevel
{
public:

    friend class Projection;
    friend class MacProj;
    friend class Diffusion;

    ~NavierStokesBase () override;

    NavierStokesBase (NavierStokesBase const&) = delete;
    NavierStokesBase (NavierStokesBase &&) = delete;
    NavierStokesBase& operator= (NavierStokesBase const&) = delete;
    NavierStokesBase& operator= (NavierStokesBase &&) = delete;

    void define_workspace ();

    ////////////////////////////////////////////////////////////////////////////
    //    AmrLevel virtual functions                                          //
    ////////////////////////////////////////////////////////////////////////////

    void checkPoint (const std::string& dir,
                     std::ostream&      os,
                     amrex::VisMF::How         how = amrex::VisMF::OneFilePerCPU,
                     bool               dump_old = true) override;

    void computeInitialDt (int                   finest_level,
                           int                   sub_cycle,
                           amrex::Vector<int>&           n_cycle,
                           const amrex::Vector<amrex::IntVect>& ref_ratio,
                           amrex::Vector<amrex::Real>&          dt_level,
                           amrex::Real                  stop_time) override;

    void computeNewDt (int                   finest_level,
                       int                   sub_cycle,
                       amrex::Vector<int>&           n_cycle,
                       const amrex::Vector<amrex::IntVect>& ref_ratio,
                       amrex::Vector<amrex::Real>&          dt_min,
                       amrex::Vector<amrex::Real>&          dt_level,
                       amrex::Real                  stop_time,
                       int                   post_regrid_flag) override;
    //
    // Error estimation for regridding.
    //
    void errorEst (amrex::TagBoxArray& tb,
                   int          clearval,
                   int          tagval,
                   amrex::Real         time,
                   int          n_error_buf = 0,
                   int          ngrow = 0) override;
    //
    // This function fills a new level n with the best
    // level n and coarser data available (usually during regrid).
    //
    void init (amrex::AmrLevel& old) override;
    //
    // Init data on this level after regridding if old level
    // did not exist previously.
    //
    void init () override;
    //
    // Called in grid_places after other tagging routines to modify
    //   the list of tagged points
    //
    void manual_tags_placement (amrex::TagBoxArray&    tags,
                                const amrex::Vector<amrex::IntVect>& bf_lev) override;
    //
    // Estimate the end of the simulation for amrLevel.
    //
    int okToContinue () override;
    //
    // Check whether simulation has reached steady state, i.e. whether the
    // change from last iteration is below a prescribed threshold steady_tol.
    //
    virtual int steadyState ();
    //
    // Build any additional data structures after regrid.
    //
    void post_regrid (int lbase, int new_finest) override;
    //
    // Build any additional data structures after restart.
    //
    void post_restart () override;
    //
    // Sync state and pressure at the end of a composite timestep.
    //
    void post_timestep (int iteration) override;
    //
    void restart (amrex::Amr&          papa,
                  std::istream& is,
                  bool          bReadSpecial = false) override;
    //
    // Old checkpoint file may have different number of state types than the
    // new source code.
    //
    void set_state_in_checkpoint (amrex::Vector<int>& state_in_checkpoint) override;
    //
    // Set time levels of state data.
    //
    void setTimeLevel (amrex::Real time,
                       amrex::Real dt_old,
                       amrex::Real dt_new) override;
    //
    // A string written as the first item in writePlotFile() at level zero.
    // It is so we can distinguish between different types of plot files.
    //
    // For NavierStokes it has the form: NavierStokes-Vnnn
    //
    std::string thePlotFileType () const override;

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokesBase public functions                                   //
    ////////////////////////////////////////////////////////////////////////////

    const amrex::MultiFab* Area () const { return area; }
    const amrex::MultiFab& Area (int dir) const { return area[dir]; }
    const amrex::MultiFab& Volume () const { return volume; }
    //
    // Get rho at time.
    // Must be one of AmrOldTime, AmrHalfTime or AmrNewTime.
    // Aborts if not one of the three above times.
    //
    const amrex::MultiFab& get_rho (amrex::Real time);
    //
    // Get divU from the state data.
    //
    amrex::MultiFab* getDivCond (int ngrow, amrex::Real time);
    //
    // Get pressure gradient data, fill-patch ghost cells.
    //
    void computeGradP(amrex::Real time);

#ifdef AMREX_USE_EB
    const amrex::MultiFab* VolFrac () const { return volfrac; }
#endif
    //
    // Boundary condition access functions.
    //
    amrex::Vector<int> fetchBCArray (int State_Type, const amrex::Box& bx,
                                     int scomp, int ncomp);

    amrex::Vector<amrex::BCRec> fetchBCArray (int State_Type, int scomp, int ncomp);

    //
    // ls related
    //
    void fill_allgts (amrex::MultiFab& mf, int type, int scomp, int ncomp, amrex::Real time);
    void reinit ();
    void reinitialization_sussman (amrex::Real dt, int loop_iter);
    void reinitialization_consls (amrex::Real dt, int loop_iter, amrex::Real epsG, amrex::Real epsG2);
    void phi_to_sgn0 (amrex::MultiFab& phi);
    amrex::MultiFab& get_phi_half_time ();
    void rk_first_reinit (amrex::MultiFab& phi_ctime,
                          amrex::MultiFab& phi2,
                          amrex::MultiFab& phi3,
                          amrex::MultiFab& sgn0,
                          amrex::MultiFab& G0,
                          amrex::Real dt, 
                          amrex::MultiFab& phi_ori);
    void rk_second_reinit (amrex::MultiFab& phi_ctime,
                          amrex::MultiFab& phi2,
                          amrex::MultiFab& phi3,
                          amrex::MultiFab& sgn0,
                          amrex::MultiFab& G0,
                          amrex::Real dt,
                          amrex::MultiFab& phi_ori);
    void mass_fix (amrex::MultiFab& phi_ctime,
                          amrex::MultiFab& phi_original,
                          amrex::MultiFab& phi2,
                          amrex::MultiFab& phi3,
                          amrex::MultiFab& ld,
                          amrex::MultiFab& lambdad,
                          amrex::MultiFab& deltafunc_alias,
                          amrex::Real delta_t,
                          int  loop_iter);

    amrex::Vector<amrex::BCRec> const& get_bcrec_velocity () const noexcept { return m_bcrec_velocity; }
    amrex::Vector<amrex::BCRec> const& get_bcrec_scalars () const noexcept { return m_bcrec_scalars; }

    amrex::BCRec const* get_bcrec_velocity_d_ptr () const noexcept {
        return m_bcrec_velocity_d.dataPtr(); }
    amrex::BCRec const* get_bcrec_scalars_d_ptr () const noexcept {
        return m_bcrec_scalars_d.dataPtr(); }

    //
    // Select appropriate AMReX average_down() based on EB/non-EB and dimensionality
    // for use in derived classes and Projection sync fns and initial vel proj
    //
    void average_down(const amrex::MultiFab& S_fine, amrex::MultiFab& S_crse,
              int scomp, int ncomp);

    void calc_mut_LES(amrex::MultiFab* mu_LES[AMREX_SPACEDIM], amrex::Real time);

    void time_average(amrex::Real& time_avg, amrex::Real& time_avg_fluct, amrex::Real& dt_avg, const amrex::Real& dt_level);

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokesBase public static functions                            //
    ////////////////////////////////////////////////////////////////////////////

    //
    // Returns the value of "gravity" for use in the projection outflow bcs.
    //
    static amrex::Real getGravity () { return gravity; }

    //
    // Returns the value of "gravity" for use in the projection outflow bcs.
    //
#ifdef AMREX_USE_EB
    static std::string const& getRedistType () { return redistribution_type; }
#endif

    static int DoTrac2() { return NavierStokesBase::do_trac2; }

    static bool GodunovUseForcesInTrans () {return godunov_use_forces_in_trans;}

    //////////////////////////////////////////////////////////////////

    NavierStokesBase ();

    NavierStokesBase (amrex::Amr&            papa,
                      int                    lev,
                      const amrex::Geometry& level_geom,
                      const amrex::BoxArray& bl,
                      const amrex::DistributionMapping& dm,
                      amrex::Real            time);

    //////////////////////////////////////////////////////////////////
    //    NavierStokesBase functions                                //
    //////////////////////////////////////////////////////////////////

    void advance_cleanup (int iteration,int ncycle);

    static void diffuse_scalar_setup (int sigma, int& rho_flag);
    //
    // Get the forcing term.
    //
    virtual void getForce (amrex::FArrayBox&       force,
                           const amrex::Box&       bx,
                           int                     scomp,
                           int                     ncomp,
                           amrex::Real             time,
                           const amrex::FArrayBox& State,
                           const amrex::FArrayBox& Aux,
                           int                     auxScomp,
                           const amrex::MFIter&    mfi);

    auto& getAdvFluxReg () {
        AMREX_ASSERT(advflux_reg);
        return *advflux_reg;
    }
    auto& getAdvFluxReg (int lev) {
        return getLevel(lev).getAdvFluxReg();
    }

    amrex::FluxRegister& getViscFluxReg () {
        AMREX_ASSERT(viscflux_reg);
        return *viscflux_reg;
    }
    amrex::FluxRegister& getViscFluxReg (int lev) {
        return getLevel(lev).getViscFluxReg();
    }
    //
    // Get rho at time n+1/2
    //
    amrex::MultiFab& get_rho_half_time ();
    //
    // Get dSdt from the state data.
    //
    [[nodiscard]] amrex::MultiFab* getDsdt (int ngrow, amrex::Real time);
    //
    void incrRhoAvg (amrex::Real alpha);
    void incrRhoAvg (const amrex::MultiFab& rho_incr,
                     int             sComp,
                     amrex::Real            alpha);

    void initial_velocity_diffusion_update (amrex::Real dt);
    //
    // Compute level projection.
    //
    void level_projector (amrex::Real dt,
                          amrex::Real time,
                          int  iteration);
    //
    // Compute the level sync correction.
    //
    void level_sync (int crse_iteration);
    //
    // Impose divergence constraint on MAC velocities.
    //
    void mac_project (amrex::Real      time,
                      amrex::Real      dt,
                      amrex::MultiFab& S_old,
                      amrex::MultiFab* divu,
                      int       ngrow,
                      bool      increment_vel_register);
    //
    // Make rho at time n.
    //
    void make_rho_prev_time ();
    //
    // Make rho at time n+1.
    //
    void make_rho_curr_time ();
    //
    // This function estimates the initial timesteping used by the model.
    //
    void post_init_estDT (amrex::Real&        dt_init,
                          amrex::Vector<int>&  nc_save,
                          amrex::Vector<amrex::Real>& dt_save,
                          amrex::Real         stop_time);
    //
    // Ensure state is consistent, i.e. velocity field is nondivergent,
    // coarse level data are averages of fine level data, pressure is zero.
    //
    void post_init_state ();
    //
    // Interpolate cell-centered cync correction from coarse to fine.
    //
    enum SyncInterpType
    {
        PC_T,
        CellCons_T,
        CellConsLin_T,
        CellConsProt_T
    };
    //
    void SyncInterp (amrex::MultiFab&      CrseSync,
                     int            c_lev,
                     amrex::MultiFab&      FineSync,
                     int            f_lev,
                     amrex::IntVect&       ratio,
                     int            src_comp,
                     int            dest_comp,
                     int            num_comp,
                     int            increment,
                     amrex::Real           dt_clev,
                     int**          bc_orig_qty,
                     SyncInterpType which_interp = CellCons_T,
                     int            state_comp   = -1);
    //
    // Bilinear interpolate nodal pressures from coarse to fine.
    //
    void SyncProjInterp (amrex::MultiFab& phi,
                         int       c_lev,
                         amrex::MultiFab& P_new,
                         amrex::MultiFab& P_old,
                         int       f_lev,
                         amrex::IntVect&  ratio);

    void sync_setup (amrex::MultiFab*& DeltaSsync);
    static void sync_cleanup (amrex::MultiFab*& DeltaSsync);
    //
    // Grow by nGrow and fillpatch the MAC-projected velocities.
    //
    void create_umac_grown (int nGrow,
                            const amrex::MultiFab *a_divu);
    //
    // Predict face velocities (before projection), estimate current cfl
    // from t^n cell-centered velocities (plus cfl and change_max)
    //
    amrex::Real predict_velocity (amrex::Real dt);
    //
    // Advect velocities.
    //
    void velocity_advection (amrex::Real dt);
    //
    // Update velocities (or momentum) with advection and diffusion terms.
    //
    void velocity_update (amrex::Real dt);
    //
    // Update velocities (or momentum) with advection terms.
    //
    void velocity_advection_update (amrex::Real dt);
    //
    // Average down State_Type and Press_Type.
    //
    void avgDown_StatePress ();

    void LES_setDomainBC (std::array<amrex::LinOpBCType,AMREX_SPACEDIM>& mlmg_lobc,
                          std::array<amrex::LinOpBCType,AMREX_SPACEDIM>& mlmg_hibc,
                          int src_comp);


    //////////////////////////////////////////////////////////////////
    //    NavierStokesBase virtual functions                        //
    //////////////////////////////////////////////////////////////////

    virtual void advance_setup (amrex::Real time,
                                amrex::Real dt,
                                int  iteration,
                                int  ncycle);

    virtual void avgDown () = 0; // Average down for all the state types.

    virtual void calc_divu (amrex::Real      time,
                            amrex::Real      dt,
                            amrex::MultiFab& fab) = 0;
    //
    // Calculate nonuniform viscosity and diffusivity
    //
    virtual void calcViscosity (amrex::Real time,
                                amrex::Real dt,
                                int  iteration,
                                int  ncycle) = 0;
    virtual void calcDiffusivity (amrex::Real time) = 0;
    //
    // Create the RHS for the MAC projection.
    //
    virtual void create_mac_rhs (amrex::MultiFab& rhs, int nGrow, amrex::Real time, amrex::Real dt);

    //
    // Timestep estimation functions follow.
    //
    virtual amrex::Real estTimeStep ();

    virtual void getViscosity (amrex::MultiFab*  viscosity[AMREX_SPACEDIM],
                               amrex::Real time) = 0;
    //
    // Compute viscous terms.
    //
    virtual void getViscTerms (amrex::MultiFab& visc_terms,
                               int              src_comp,
                               int              num_comp,
                               amrex::Real      time) = 0;
    //
    virtual void mac_sync () = 0;
    //
    virtual void reflux () = 0;
    //
    // Reset time levels for the initial iterations.
    //
    virtual void resetState (amrex::Real time,
                             amrex::Real dt_old,
                             amrex::Real dt_new);

    virtual void scalar_advection_update (amrex::Real dt,
                                          int  first_scalar,
                                          int  last_scalar);
    //
    virtual void sum_integrated_quantities () = 0;

    virtual void velocity_diffusion_update (amrex::Real dt) = 0;

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokesBase static functions                                   //
    ////////////////////////////////////////////////////////////////////////////

    //
    // Cleanup data descriptors at end of run.
    //
    static void variableCleanUp ();
    //
    // Read input file
    //
    static void Initialize ();
    static void Finalize ();
    //
    // Check for RZ geometry
    //
    static void read_geometry ();

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokesBase particle functions                                 //
    ////////////////////////////////////////////////////////////////////////////

#ifdef AMREX_PARTICLES
    static amrex::AmrTracerParticleContainer* theNSPC ();
    static void read_particle_params ();

    void initParticleData ();
    void post_restart_particle ();

    void post_timestep_particle (int iteration);

    virtual int timestamp_num_extras () { return 0; }
    virtual void timestamp_add_extras (int lev, amrex::Real time, amrex::MultiFab& mf) { }

    std::unique_ptr<amrex::MultiFab> ParticleDerive (const std::string& name, amrex::Real time, int ngrow);
    void ParticleDerive (const std::string& name, amrex::Real time, amrex::MultiFab& mf, int dcomp);
#endif

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokesBase diagnostics functions                              //
    ////////////////////////////////////////////////////////////////////////////
    void printMaxVel (bool new_data = true);

    void printMaxGp (bool new_data = true);

    void printMaxValues (bool new_data = true);

    ////////////////////////////////////////////////////////////////////////////

    void buildMetrics (); // 1-D metric arrays for RZ
    //
    void calc_dsdt (amrex::Real      time,
            amrex::Real      dt,
            amrex::MultiFab& dsdt);
    //
    NavierStokesBase& getLevel (int lev) {
        return dynamic_cast<NavierStokesBase&> ( parent->getLevel(lev) );
    }
    //
    SyncRegister& getSyncReg () {
        AMREX_ASSERT(sync_reg);
        return *sync_reg;
    }
    //
    static void getOutFlowFaces (amrex::Vector<amrex::Orientation>& outFaces);
    //
    // Get state data.
    //
    [[nodiscard]] amrex::MultiFab* getState (int  ngrow,
                                             int  state_indx,
                                             int  strt_comp,
                                             int  num_comp,
                                             amrex::Real time);
    //
    // Compile p_avg in advance.
    //
    void incrPAvg ();
    void initOldFromNew (int type, int lev = -1); // Initialize old with new
    //
    // Compile rho_avg in advance.
    //
    void initRhoAvg (amrex::Real alpha);
    //
    static void init_additional_state_types ();
    //
    amrex::Real initialTimeStep ();

    //
    // Correct a conservatively-advected scalar for under-over shoots.
    //
    void ConservativeScalMinMax ( amrex::MultiFab& Snew, int snew_comp, int new_density_comp,
                                  amrex::MultiFab const& Sold, int sold_comp, int old_density_comp );

    //
    // Correct a convectively-advected scalar for under-over shoots.
    //
    void ConvectiveScalMinMax ( amrex::MultiFab& Snew, int snew_comp,
                                amrex::MultiFab const& Sold, int sold_comp );
    //
    // Floor small values of states to be extrapolated
    //
    static void floor(amrex::MultiFab& mf);

#ifdef AMREX_USE_EB
    void init_eb (const amrex::Geometry& level_geom, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    void initialize_eb2_structs();
    void define_body_state();

    void InitialRedistribution ();

    static bool ebInitialized();
    static bool eb_initialized;
    static bool body_state_set;
    static bool      no_eb_in_domain;
    static std::vector<amrex::Real> body_state;
#endif

    void ComputeAofs (int comp, int ncomp,
                      amrex::MultiFab const& state,
                      int S_comp,
                      amrex::MultiFab const& forcing_term,
                      amrex::MultiFab const& divu,
                      bool is_velocity, amrex::Real dt);

    void ComputeAofs (amrex::MultiFab& advc, int a_comp,
                      int state_indx, int ncomp,
                      amrex::MultiFab const& S, int S_comp,
                      amrex::MultiFab const* forcing, int f_comp,
                      amrex::MultiFab const* divu,
                      amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& cfluxes, int flux_comp,
                      amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& edgestate, int edge_comp,
                      bool known_edge_state,
                      bool is_velocity, amrex::Real dt,
                      bool is_sync, amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> const& U_corr,
                      bool do_crse_add, bool do_fine_add);
//
// Data
//

    // volume and area are intentionally without EB knowledge
    amrex::MultiFab volume;
    amrex::MultiFab area[AMREX_SPACEDIM];

    std::unique_ptr<amrex::iMultiFab> coarse_fine_mask;

#ifdef AMREX_USE_EB
    void set_body_state(amrex::MultiFab& S);

    const amrex::MultiFab* volfrac;
    amrex::FabArray<amrex::BaseFab<int>> ebmask;

    std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> areafrac;
    std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> facecent;

#endif


    //
    // MAC edge velocities.
    //
    amrex::MultiFab* u_mac = nullptr;
    int  umac_n_grow;
    //
    // Advective update terms.
    //
    amrex::MultiFab* aofs = nullptr;
    //
    // BC info for advection
    //
    amrex::Vector<amrex::BCRec> m_bcrec_velocity;
    amrex::Vector<amrex::BCRec> m_bcrec_scalars;
    // NOTE that DeviceVector reverts to amrex::Vector if !AMREX_USE_GPU
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_velocity_d;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_scalars_d;

    std::unique_ptr<Diffusion> diffusion;
    //
    // arrays for variable viscosity and diffusivity
    //
    amrex::MultiFab *diffn_cc, *diffnp1_cc;
    amrex::MultiFab *viscn_cc, *viscnp1_cc;
    //
    // Sync update quantities.
    //
    // Average rho over a composite timestep.
    // Used only in the multlilevel sync projection
    //
    amrex::MultiFab rho_avg;
    //
    // Average p over a composite timestep.
    // Used only to average fine pressure data to coarse.
    //
    amrex::MultiFab p_avg;

    amrex::MultiFab Vsync;    // Velocity sync update storage
    amrex::MultiFab Ssync;    // Scalar sync update storage
    //
    // Density at time n+1/2 (used in advance).
    //
    amrex::MultiFab rho_half;
    amrex::MultiFab* rho_qtime = nullptr;
    amrex::MultiFab* rho_tqtime = nullptr;
    //
    // Density at prev_time used in advance().
    //
    amrex::MultiFab rho_ptime;
    //
    // Density at cur_time used in advance().
    //
    amrex::MultiFab rho_ctime;
    //
    // Data structure used to compute RHS for sync project.
    //
    SyncRegister* sync_reg = nullptr;
    //
    // Data structures to store advective and viscous refluxing quantities
    // on coarse-fine grid interfaces.
    //
#ifdef AMREX_USE_EB
    std::unique_ptr<amrex::EBFluxRegister> advflux_reg;
#else
    std::unique_ptr<amrex::YAFluxRegister> advflux_reg;
#endif
    amrex::FluxRegister* viscflux_reg = nullptr;
    //
    // Radii for r-z calculations.
    //
    amrex::Vector< amrex::Vector<amrex::Real> > radius;
    //
    // Static objects.
    //
    static amrex::BCRec phys_bc;
    static Projection*  projector;
    static MacProj*     mac_projector;
    //
    // Internal parameters for timestepping.
    //
    static amrex::Real init_shrink;   // reduction factor of first estimated timestep
    static int  init_iter;            // # of timestep iterations for initial pressure
    static int  init_vel_iter;        // # of iterations for initial velocity projection
    static amrex::Real cfl;           // desired maximum cfl
    static amrex::Real change_max;    // maximum change in dt over a timestep
    static amrex::Real fixed_dt;      // set > 0 to specify dt
    static amrex::Real init_dt;       // set > 0 to specify initial dt at coarsest level
    static bool stop_when_steady;     // set to true if simulation should stop at steady-state
    static amrex::Real steady_tol;    // tolerance for assuming steady-state has been reached
    static int  initial_iter;         // flag for initial pressure iterations
    static int  initial_step;         // flag for initial iterations
    static amrex::Real dt_cutoff;     // minimum dt allowed
    static int  sum_interval;         // number of timesteps for conservation stats
    //
    // Internal parameters for options.
    //
    static int          radius_grow;
    static int          verbose;
    static amrex::Real  gravity;
    static int          NUM_SCALARS;      // number of non-velocity components
    static int          NUM_STATE;        // total number of state variables
#ifdef AMREX_USE_EB
    static const int    NUM_GROW=4;       // number of ghost cells
    static int          refine_cutcells;
#else
    //
    // ls related
    // change to two gts
    // 
    static const int    NUM_GROW=2;       // number of ghost cells
#endif
    //
    // Controls over how the convective and diffusive terms are calculated
    // for each state variable.
    //
    static amrex::Vector<AdvectionForm> advectionType;
    static amrex::Vector<DiffusionForm> diffusionType;
    //
    // Viscosity parameters.
    //
    static amrex::Vector<int>  is_diffusive;      // does variable diffuse?
    static amrex::Vector<amrex::Real> visc_coef;  // const coef viscosity terms
    static amrex::Real        visc_tol;
    static amrex::Real        visc_abs_tol;
    static amrex::Real        be_cn_theta;
    //
    // Internal switches.
    //
    static int  Tracer;
    static int  Tracer2;
    static int  Temp;
    static int  do_trac2;
    static int  do_temp;
    static int  do_cons_trac;
    static int  do_cons_trac2;
    static int  do_sync_proj;
    static int  do_reflux;
    static int  do_mac_proj;
    static int  do_refine_outflow;          // The following three items control
    static int  do_derefine_outflow;        //   how grids are treated at
    static int  Nbuf_outflow;               //   outflow boundaries.
    static int  do_denminmax;               // The code for these was in NavierStokes.cpp,
    static int  do_scalminmax;              //   but the flags were not declared or read in.
    static int  getForceVerbose;            // Does exactly what it says on the tin
    //
    // LES parameters
    //
    static int  do_LES;             // To activate LES modelling
    static int  getLESVerbose;      // Does exactly what it says on the tin
    static std::string LES_model;   // To choose a LES model
    static amrex::Real smago_Cs_cst;       // A parameter for the Smagorinsky model, usually set to 0.18
    static amrex::Real sigma_Cs_cst;       // A parameter for the Sigma model, usually set to 1.5
    //
    // Parameters for averaging
    //
    static amrex::Vector<amrex::Real> time_avg;
    static amrex::Vector<amrex::Real> time_avg_fluct;
    static amrex::Vector<amrex::Real> dt_avg;
    static int avg_interval;
    static int compute_fluctuations;
    //
    // Members for non-zero divu.
    //
    static int  additional_state_types_initialized;
    static int  Divu_Type;
    static int  Dsdt_Type;
    static int  Average_Type;
    static int  num_state_type;
    static int  have_divu;
    static int  have_dsdt;
    static int  do_init_vort_proj;
    static int  do_init_proj;
    //
    // Control velocity vs momentum update
    //
    static int  do_mom_diff;

    //
    // skip_level_projector
    //
    static int  skip_level_projector;

    static int isolver;

    //
    // ls related
    //
    static int  do_phi;
    static int  phicomp;
    static int  epsilon;
    amrex::MultiFab phi_half;
    amrex::MultiFab phi_ptime;
    amrex::MultiFab phi_ctime;
    amrex::MultiFab heaviside;

    static amrex::Real mu_a;
    static amrex::Real rho_a;
    static amrex::Real mu_w;
    static amrex::Real rho_w;

    static int do_reinit;
    static int lev0step_of_reinit;
    static int number_of_reinit;
    static int reinit_levelset;

    // for reinit_sussman  
    amrex::MultiFab phi_original;
    amrex::MultiFab sgn0;

    static int  do_cons_phi;
    static int  prescribed_vel;

    // for cons_levelset
    static int do_cons_levelset;

    //
    // diffused ib
    //
    static int  do_diffused_ib;
    amrex::MultiFab phi_nodal;
    amrex::MultiFab pvf;
    static int  advect_and_update_scalar;
    static amrex::Real fluid_rho;

// #ifdef AMREX_PARTICLES
//     mParticle *mParticle_obj;
// #endif
    int nghost_state () const;
    static constexpr int nghost_force () { return 1; }

    //
    // static integer constants
    //
    static const int GEOM_GROW  = 1;
#ifdef AMREX_USE_EB
    // This must be consistent with what nghost_force returns for
    // godunov
    static const int gradp_grow = 4;
#else
    static const int gradp_grow = 1;
#endif

    //
    // Controls advection scheme
    //
    static std::string advection_scheme;
    static bool godunov_use_forces_in_trans;

#ifdef AMREX_USE_EB
    // Controls redistribution
    static std::string redistribution_type;
#endif

  //
  // For restart, is GradP in checkpoint file
  //
  static int gradp_in_checkpoint;

  static int average_in_checkpoint;

#ifdef AMREX_PARTICLES
  //
  // Check in case different particle methods are being used
  //
  static bool do_nspc;
  //
  // Put particle info in plotfile (using ParticleContainer::Checkpoint)?
  //
  static bool particles_in_plotfile;
#endif
};

#endif
