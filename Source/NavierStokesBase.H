
#ifndef _NavierStokesBase_H_
#define _NavierStokesBase_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrLevel.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_BLFort.H>
#include <Diffusion.H>
#include <AMReX_ErrorList.H>
#include <AMReX_FluxRegister.H>
#include <MacProj.H>
#include <Projection.H>
#include <SyncRegister.H>
#include <AMReX_Utility.H>

#ifdef AMREX_PARTICLES
#include <AMReX_AmrParticles.H>
#endif

#ifdef AMREX_USE_EB
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBFluxRegister.H>
#endif

//
// "Divu_Type" means S, where divergence U = S
// "Dsdt_Type" means pd S/pd t, where S is as above
//
// Determine what you want in the state -- Divu, Dsdt -- in
// NavierStokes::variableSetUp in NS_setup.cpp.
//
enum StateType {State_Type=0, Press_Type, Gradp_Type};

//
// Note: enumerated value NUM_STATE_TYPE no longer defined
// use num_state_type (static int member of NavierStokes)
//
#if (BL_SPACEDIM == 2)
enum StateNames  { Xvel=0, Yvel, Density};
#else
enum StateNames  { Xvel=0, Yvel, Zvel, Density};
#endif

enum PressureNames { Pressure=0 };
#if (BL_SPACEDIM == 2)
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

    virtual ~NavierStokesBase ();

    ////////////////////////////////////////////////////////////////////////////
    //    AmrLevel virtual functions                                          //
    ////////////////////////////////////////////////////////////////////////////

    virtual void checkPoint (const std::string& dir,
                             std::ostream&      os,
                             amrex::VisMF::How         how = amrex::VisMF::OneFilePerCPU,
                             bool               dump_old = true) override;

    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   amrex::Vector<int>&           n_cycle,
                                   const amrex::Vector<amrex::IntVect>& ref_ratio,
                                   amrex::Vector<amrex::Real>&          dt_level,
                                   amrex::Real                  stop_time) override;

    virtual void computeNewDt (int                   finest_level,
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
    virtual void errorEst (amrex::TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           amrex::Real         time,
                           int          n_error_buf = 0,
                           int          ngrow = 0) override;
    //
    // This function fills a new level n with the best
    // level n and coarser data available (usually during regrid).
    //
    virtual void init (amrex::AmrLevel& old) override;
    //
    // Init data on this level after regridding if old level
    // did not exist previously.
    //
    virtual void init () override;
    //
    // Called in grid_places after other tagging routines to modify
    //   the list of tagged points
    //
    virtual void manual_tags_placement (amrex::TagBoxArray&    tags,
					const amrex::Vector<amrex::IntVect>& bf_lev) override;
    //
    // Estimate the end of the simulation for amrLevel.
    //
    virtual int okToContinue () override;
    //
    // Check whether simulation has reached steady state, i.e. whether the
    // change from last iteration is below a prescribed threshold steady_tol.
    //
    virtual int steadyState ();
    //
    // Build any additional data structures after regrid.
    //
    virtual void post_regrid (int lbase, int new_finest) override;
    //
    // Build any additional data structures after restart.
    //
    virtual void post_restart () override;
    //
    // Sync state and pressure at the end of a composite timestep.
    //
    virtual void post_timestep (int iteration) override;
    //
    virtual void restart (amrex::Amr&          papa,
                          std::istream& is,
                          bool          bReadSpecial = false) override;
    //
    // Old checkpoint file may have different number of state types than the
    // new source code.
    //
    virtual void set_state_in_checkpoint (amrex::Vector<int>& state_in_checkpoint) override;
    //
    // Set time levels of state data.
    //
    virtual void setTimeLevel (amrex::Real time,
                               amrex::Real dt_old,
                               amrex::Real dt_new) override;
    //
    // A string written as the first item in writePlotFile() at level zero.
    // It is so we can distinguish between different types of plot files.
    //
    // For NavierStokes it has the form: NavierStokes-Vnnn
    //
    virtual std::string thePlotFileType () const override;

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
    void computeGradP(int clev, int flev);

#ifdef AMREX_USE_EB
    const amrex::MultiFab* VolFrac () const { return volfrac; }
#endif
    //
    // Boundary condition access functions.
    //
    amrex::Vector<int> fetchBCArray (int State_Type, const amrex::Box& bx,
				     int scomp, int ncomp);

    amrex::Vector<amrex::BCRec> fetchBCArray (int State_Type, int scomp, int ncomp);

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

    void calc_mut_LES(amrex::MultiFab*  diffusivity[BL_SPACEDIM], amrex::Real time);

    void time_average(amrex::Real& time_avg, amrex::Real& time_avg_fluct, amrex::Real& dt_avg, const amrex::Real& dt_level);

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokesBase public static functions                            //
    ////////////////////////////////////////////////////////////////////////////

    //
    // Returns the value of "gravity" for use in the projection outflow bcs.
    //
    static amrex::Real getGravity () { return gravity; }

    static int DoTrac2();

    static bool GodunovUsePPM () {return godunov_use_ppm;};
    static bool GodunovUseForcesInTrans () {return godunov_use_forces_in_trans;};
    static int  GodunovHypgrow () {return godunov_hyp_grow;};

    //////////////////////////////////////////////////////////////////
  
    NavierStokesBase ();

    NavierStokesBase (amrex::Amr&            papa,
		      int             lev,
		      const amrex::Geometry& level_geom,
		      const amrex::BoxArray& bl,
                      const amrex::DistributionMapping& dm,
		      amrex::Real            time);

    //////////////////////////////////////////////////////////////////
    //    NavierStokesBase functions                                //
    //////////////////////////////////////////////////////////////////

    void advance_cleanup (int iteration,int ncycle);

    void diffuse_scalar_setup (int sigma, int& rho_flag);
    //
    // Get the forcing term.
    //
    virtual void getForce (amrex::FArrayBox&       force,
                           const amrex::Box&       bx,
                           int                     strt_comp,
                           int                     num_comp,
                           const amrex::Real       time,
                           const amrex::FArrayBox& Vel,
                           const amrex::FArrayBox& Scal,
                           int                     scalScomp);
    //
    amrex::FluxRegister& getAdvFluxReg () {
        BL_ASSERT(advflux_reg);
        return *advflux_reg;
    }
    //
    amrex::FluxRegister& getAdvFluxReg (int lev) {
        return getLevel(lev).getAdvFluxReg();
    }
    //
    amrex::FluxRegister& getViscFluxReg () {
        BL_ASSERT(viscflux_reg);
        return *viscflux_reg;
    }
    //
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
    amrex::MultiFab* getDsdt (int ngrow, amrex::Real time);
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
    void sync_cleanup (amrex::MultiFab*& DeltaSsync);
    //
    // Grow by nGrow and fillpatch the MAC-projected velocities.
    //
    void create_umac_grown (int nGrow);
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


    amrex::Real volWgtSum (const std::string& name,
                    amrex::Real               time);

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
    virtual void calcViscosity (const amrex::Real time,
				const amrex::Real dt,
				const int  iteration,
				const int  ncycle) = 0;
    virtual void calcDiffusivity (const amrex::Real time) = 0;
    //
    // Create the RHS for the MAC projection.
    //
    virtual void create_mac_rhs (amrex::MultiFab& rhs, int nGrow, amrex::Real time, amrex::Real dt);

    //
    // Timestep estimation functions follow.
    //
    virtual amrex::Real estTimeStep ();

    virtual void getViscosity (amrex::MultiFab*  viscosity[BL_SPACEDIM],
			       const amrex::Real time) = 0;
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

    static void variableCleanUp (); // Cleanup data descriptors at end of run.
    static void Initialize ();   // Read input file common to IAMR and PeleLM
    static void Initialize_specific ();   // Read input file specific to IAMR
    static void Finalize ();
    static void read_geometry (); // Set metrics

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
        BL_ASSERT(sync_reg);
        return *sync_reg;
    }
    //
    void getOutFlowFaces (amrex::Vector<amrex::Orientation>& outFaces);
    //
    // Get state data.
    //
    amrex::MultiFab* getState (int  ngrow,
			       int  state_indx,
			       int  strt_comp,
			       int  num_comp,
			       amrex::Real time);
    //
    // Compile p_avg in advance.
    //
    void incrPAvg ();
    void initOldFromNew (int type); // Initialize old with new
    //
    // Compile rho_avg in advance.
    //
    void initRhoAvg (amrex::Real alpha);
    //
    void init_additional_state_types ();
    //
    amrex::Real initialTimeStep ();

    //
    // Correct a conservatively-advected scalar for under-over shoots.
    //
    void ConservativeScalMinMax ( amrex::MultiFab& Snew, const int snew_comp, const int new_density_comp,
                                  amrex::MultiFab const& Sold, const int sold_comp, const int old_density_comp );

    //
    // Correct a convectively-advected scalar for under-over shoots.
    //
    void ConvectiveScalMinMax ( amrex::MultiFab& Snew, const int snew_comp,
                                amrex::MultiFab const& Sold, const int sold_comp );
    //
    // Floor small values of states to be extrapolated
    //
    void floor(amrex::MultiFab& mf);

#ifdef AMREX_USE_EB

    void init_eb (const amrex::Geometry& level_geom, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    void initialize_eb2_structs();
    void define_body_state();


    static bool ebInitialized();
    static bool eb_initialized;
    static bool body_state_set;
    static bool      no_eb_in_domain;
    static std::vector<amrex::Real> body_state;

#endif

//
// Data
//

    // volume and area are intentionally without EB knowledge
    amrex::MultiFab volume;
    amrex::MultiFab area[BL_SPACEDIM];

#ifdef AMREX_USE_EB

    void set_body_state(amrex::MultiFab& S);

  //fixme - don't really need to carry these around anymore....
    const amrex::MultiFab* volfrac;
    const amrex::MultiCutFab* bndrycent;
    amrex::FabArray<amrex::BaseFab<int>> ebmask;

    std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> areafrac;
    std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> facecent;

#endif

    //
    // MAC edge velocities.
    //
    amrex::MultiFab* u_mac;
    int  umac_n_grow;
    //
    // Advective update terms.
    //
    amrex::MultiFab* aofs;
    //
    // BC info for advection
    //
    amrex::Vector<amrex::BCRec> m_bcrec_velocity;
    amrex::Vector<amrex::BCRec> m_bcrec_scalars;
    // NOTE that DeviceVector reverts to amrex::Vector if !AMREX_USE_GPU
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_velocity_d;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_scalars_d;
  
    Diffusion* diffusion;
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
    amrex::MultiFab* rho_qtime;
    amrex::MultiFab* rho_tqtime;
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
    SyncRegister* sync_reg;
    //
    // Data structures to store advective and viscous refluxing quantities
    // on coarse-fine grid interfaces.
    //
    amrex::FluxRegister* advflux_reg;
    amrex::FluxRegister* viscflux_reg;
    //
    // Radii for r-z calculations.
    //
    amrex::Vector< amrex::Vector<amrex::Real> > radius;
    //
    // Intersection of current amrex::BoxArray with that from before regrid.
    //
    amrex::BoxArray old_intersect_new;
    //
    // Static objects.
    //
    static amrex::ErrorList   err_list;
    static amrex::BCRec       phys_bc;
    static Projection* projector;
    static MacProj*    mac_projector;
    //
    // Internal parameters for timestepping.
    //
    static amrex::Real init_shrink;   // reduction factor of first esimated timestep
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
  // fixme --- not sure that NUM_GROW or refine_cutcells are still used...
    static const int    NUM_GROW=4;       // number of ghost cells
    static int          refine_cutcells;
#else
    static const int    NUM_GROW=1;       // number of ghost cells
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
    static int  modify_reflux_normal_vel;
    static int  do_mac_proj;
    static int  do_refine_outflow;          // The following three items control
    static int  do_derefine_outflow;        //   how grids are treated at
    static int  Nbuf_outflow;               //   outflow boundaries.
    static int  do_denminmax;               // The code for these was in NavierStokes.cpp,
    static int  do_scalminmax;              //   but the flags were not declared or read in.
    static int  do_density_ref;             // Flags added for refinement control
    static int  do_tracer_ref;
    static int  do_tracer2_ref;
    static int  do_vorticity_ref;
    static int  do_temp_ref;
    static int  do_scalar_update_in_order;  // Flags to allow evaluation of source terms
    static amrex::Vector<int> scalarUpdateOrder;
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
    static amrex::Real divu_relax_factor;
    static int  do_init_vort_proj;
    static int  do_init_proj;
    //
    // Running statistics controls
    //
    static amrex::Real volWgtSum_sub_origin_x;
    static amrex::Real volWgtSum_sub_origin_y;
    static amrex::Real volWgtSum_sub_origin_z;
    static amrex::Real volWgtSum_sub_Rcyl;
    static amrex::Real volWgtSum_sub_dx;
    static amrex::Real volWgtSum_sub_dy;
    static amrex::Real volWgtSum_sub_dz;
    //
    // Control velocity vs momentum update
    //
    static int  do_mom_diff;
    //
    // static integer constants
    //
  // FIXME?  do these need to be adjusted for EB??? see NUM_GROW above
    static const int GEOM_GROW  = 1;
    static const int PRESS_GROW = 1;
    static const int DIVU_GROW  = 1;
    static const int DSDT_GROW  = 1;
#ifdef AMREX_USE_EB
    // MOL doesn't require gradp ghost cells
    static const int gradp_grow = 0;
#else
    static const int gradp_grow = 1;
#endif
    //
    // Controls Godunov solver
    //
    static bool godunov_use_ppm;
    static bool godunov_use_forces_in_trans;
#ifdef AMREX_USE_EB
    static const int  godunov_hyp_grow = 4;
#else
    static const int  godunov_hyp_grow = 3;
#endif
  //
  // For restart, is GradP in checkpoint file 
  //
  static int gradp_in_checkpoint;

  static int average_in_checkpoint;

#ifdef AMREX_PARTICLES
  //
  // Put particle info in plotfile (using ParticleContainer::Checkpoint)?
  //
  static bool particles_in_plotfile;
#endif
};

#endif
