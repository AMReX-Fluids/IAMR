
#ifndef _NavierStokesBase_H_
#define _NavierStokesBase_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AmrLevel.H>
#include <BL_CXX11.H>
#include <BC_TYPES.H>
#include <BLFort.H>
#include <Diffusion.H>
#include <ErrorList.H>
#include <FluxRegister.H>
#include <Godunov.H>
#include <MacProj.H> 
#include <Projection.H>
#include <SyncRegister.H>
#include <Utility.H>
#include <ViscBndry.H>

#ifdef PARTICLES
#include <AmrParticles.H>
//
// The "BL_SPACEDIM" bytes in m_data[] are to hold temporary positions or velocities
//      in AdvectWithUmac().
//
typedef AmrParticleContainer<BL_SPACEDIM> NSParticleContainer;
#endif

#ifdef USE_BPM
#include <BPM.H>
#endif

//
// "Divu_Type" means S, where divergence U = S
// "Dsdt_Type" means pd S/pd t, where S is as above
//
// Determine what you want in the state -- Divu, Dsdt -- in
// NavierStokes::variableSetUp in NS_setup.cpp.
//
enum StateType {State_Type=0, Press_Type};

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
enum DivuNames { Divu=0};
enum DsDtNames { Dsdt=0};
enum DpDtNames { Dpdt=0};

class NavierStokesBase
    :
    public AmrLevel
{
public:

    friend class Projection;
    friend class MacProj;
    friend class Diffusion;
    
    virtual ~NavierStokesBase ();

    ////////////////////////////////////////////////////////////////////////////
    //    AmrLevel virtual functions                                          //
    ////////////////////////////////////////////////////////////////////////////

    virtual void allocOldData () BL_OVERRIDE;

    virtual void checkPoint (const std::string& dir,
                             std::ostream&      os,
                             VisMF::How         how = VisMF::OneFilePerCPU,
                             bool               dump_old = true) BL_OVERRIDE;

    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   Array<int>&           n_cycle,
                                   const Array<IntVect>& ref_ratio,
                                   Array<Real>&          dt_level,
                                   Real                  stop_time) BL_OVERRIDE;

    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               Array<int>&           n_cycle,
                               const Array<IntVect>& ref_ratio,
                               Array<Real>&          dt_min,
                               Array<Real>&          dt_level,
                               Real                  stop_time,
                               int                   post_regrid_flag) BL_OVERRIDE;
    //
    // Error estimation for regridding.
    //
    virtual void errorEst (TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           Real         time,
                           int          n_error_buf = 0,
                           int          ngrow = 0) BL_OVERRIDE;
    //
    // This function fills a new level n with the best
    // level n and coarser data available (usually during regrid).
    //
    virtual void init (AmrLevel& old) BL_OVERRIDE;
    //
    // Init data on this level after regridding if old level
    // did not exist previously.
    //
    virtual void init () BL_OVERRIDE;
    //
    // Called in grid_places after other tagging routines to modify
    //   the list of tagged points
    //
    virtual void manual_tags_placement (TagBoxArray&    tags, 
					Array<IntVect>& bf_lev) BL_OVERRIDE;
    //
    // Estimate the end of the simulation for amrLevel.
    //
    virtual int okToContinue () BL_OVERRIDE;
    //
    // Build any additional data structures after regrid.
    //
    virtual void post_regrid (int lbase, int new_finest) BL_OVERRIDE;
    //
    // Build any additional data structures after restart.
    //
    virtual void post_restart () BL_OVERRIDE;
    //
    // Sync state and pressure at the end of a composite timestep.
    //
    virtual void post_timestep (int iteration) BL_OVERRIDE;
    //
    virtual void restart (Amr&          papa,
                          std::istream& is,
                          bool          bReadSpecial = false) BL_OVERRIDE;
    //
    // Set time levels of state data.
    //
    virtual void setTimeLevel (Real time,
                               Real dt_old,
                               Real dt_new) BL_OVERRIDE;
    //
    // A string written as the first item in writePlotFile() at level zero.
    // It is so we can distinguish between different types of plot files.
    //
    // For NavierStokes it has the form: NavierStokes-Vnnn
    //
    virtual std::string thePlotFileType () const BL_OVERRIDE;

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokesBase public functions                                   //
    ////////////////////////////////////////////////////////////////////////////

    const MultiFab* Area () const { return area; }
    const MultiFab& Area (int dir) const { return area[dir]; }
    const MultiFab& Volume () const { return volume; }
    //
    // Get rho at time.
    // Must be one of AmrOldTime, AmrHalfTime or AmrNewTime.
    // Aborts if not one of the three above times.
    //
    const MultiFab& get_rho (Real time);
    //
    // Get divU from the state data.
    //
    MultiFab* getDivCond (int ngrow, Real time);
    //
    // Get pressure gradient data via fill-patching.
    //
    void getGradP (MultiFab& gp, Real time);

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokesBase public static functions                            //
    ////////////////////////////////////////////////////////////////////////////

    //
    // Returns the value of "gravity" for use in the projection outflow bcs.
    //
    static Real getGravity () { return gravity; }

    static int DoTrac2();

protected:

    NavierStokesBase ();

    NavierStokesBase (Amr&            papa,
		      int             lev,
		      const Geometry& level_geom,
		      const BoxArray& bl,
		      Real            time);

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokesBase protected functions                                //
    ////////////////////////////////////////////////////////////////////////////

    void advance_cleanup (int iteration,int ncycle);

    void diffuse_scalar_setup (int sigma, int& rho_flag);
    //
    // Get the forcing term.
    //
#ifdef BOUSSINESQ
    void getForce (FArrayBox&       force,
                   int              gridno,
                   int              ngrow,
                   int              strt_comp,
                   int              num_comp,
		   const Real       time,
                   const FArrayBox& Scal);
#else
#ifdef MOREGENGETFORCE
    void getForce (FArrayBox&       force,
                   int              gridno,
                   int              ngrow,
                   int              strt_comp,
                   int              num_comp,
		   const Real       time,
                   const FArrayBox& Vel,
                   const FArrayBox& Scal,
		   int              scalScomp);
#else
    void getForce (FArrayBox&       force,
                   int              gridno,
                   int              ngrow,
                   int              strt_comp,
                   int              num_comp,
#ifdef GENGETFORCE
		   const Real       time,
#endif		 
                   const FArrayBox& Rho,
                   int              RComp = 0);
#endif		 
#endif		 
    //
    FluxRegister& getAdvFluxReg () {
        BL_ASSERT(advflux_reg);
        return *advflux_reg;
    }
    //
    FluxRegister& getAdvFluxReg (int lev) {
        return getLevel(lev).getAdvFluxReg();
    }
    //
    FluxRegister& getViscFluxReg () {
        BL_ASSERT(viscflux_reg);
        return *viscflux_reg;
    }
    //
    FluxRegister& getViscFluxReg (int lev) {
        return getLevel(lev).getViscFluxReg();
    }
    //
    // Get rho at time n+1/2
    //
    MultiFab& get_rho_half_time ();
    //
    // Get dSdt from the state data.
    //
    MultiFab* getDsdt (int ngrow, Real time);
    //
    void incrRhoAvg (Real alpha);
    void incrRhoAvg (const MultiFab& rho_incr,
                     int             sComp,
                     Real            alpha);

    void initial_velocity_diffusion_update (Real dt);

    //
    // Inject fine pressure nodes down onto coarse nodes.
    //
    void injectDown (const Box&       ovlp,
                     FArrayBox&       Pcrse,
                     const FArrayBox& Pfine,
                     IntVect&         fine_ratio);
    //
    // Compute level projection.
    //
    void level_projector (Real dt,
                          Real time,
                          int  iteration);
    //
    // Compute the level sync correction.
    //
    void level_sync (int crse_iteration);
    //
    // Impose divergence constraint on MAC velocities.
    //
    void mac_project (Real      time,
                      Real      dt,
                      MultiFab& S_old, 
                      MultiFab* divu,
                      int       have_divu,
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
    void post_init_estDT (Real&        dt_init,
			  Array<int>&  nc_save,
			  Array<Real>& dt_save,
			  Real         stop_time);
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
    void SyncInterp (MultiFab&      CrseSync,
                     int            c_lev,
                     MultiFab&      FineSync,
                     int            f_lev,
                     IntVect&       ratio,
                     int            src_comp,
                     int            dest_comp,
                     int            num_comp,
                     int            increment,
                     Real           dt_clev,
                     int**          bc_orig_qty,
                     SyncInterpType which_interp = CellCons_T,
                     int            state_comp   = -1);
    //
    // Bilinear interpolate nodal pressures from coarse to fine.
    //
    void SyncProjInterp (MultiFab& phi,
                         int       c_lev,
                         MultiFab& P_new,
                         MultiFab& P_old,
                         int       f_lev,
                         IntVect&  ratio,
                         bool      first_crse_step_after_initial_iters,
                         Real      cur_crse_pres_time,
                         Real      prev_crse_pres_time);

    void sync_setup (MultiFab*& DeltaSsync);
    void sync_cleanup (MultiFab*& DeltaSsync);

    //
    // Advect velocities.
    //
    void velocity_advection (Real dt);
    //
    // Update velocities (or momentum) after advection.
    //
    void velocity_update (Real dt);

    void velocity_advection_update (Real dt);

    Real volWgtSum (const std::string& name,
                    Real               time);

#if (BL_SPACEDIM == 3)
    void sum_turbulent_quantities ();
    void TurbSum (Real time, Real *tmp, int  ksize, int turbVars);
#ifdef SUMJET
    void sum_jet_quantities ();
    void JetSum   (Real time, 
		   Real *jetData,
		   int  levRsize,
		   int  levKsize,
		   int  rsize,
		   int  ksize,
		   int  jetVars);
#endif
#endif

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokesBase protected virtual functions                        //
    ////////////////////////////////////////////////////////////////////////////

    virtual void advance_setup (Real time,
				Real dt,
				int  iteration,
				int  ncycle);

    virtual void avgDown () = 0; // Average down for all the state types.
    //
    virtual void calc_divu (Real      time,
			    Real      dt,
			    MultiFab& fab) = 0;
    virtual void calc_dsdt (Real      time,
			    Real      dt,
			    MultiFab& dsdt) = 0;
    //
    // Calculate nonuniform viscosity and diffusivity
    //
    virtual void calcViscosity (const Real time,
				const Real dt,
				const int  iteration,
				const int  ncycle) = 0;
    virtual void calcDiffusivity (const Real time) = 0;
    //
    // Create the RHS for the MAC projection.
    //
    virtual void create_mac_rhs (MultiFab& rhs, int nGrow, Real time, Real dt);

    //
    // Timestep estimation functions follow.
    //
    virtual Real estTimeStep ();

    virtual void getViscosity (MultiFab* viscosity[BL_SPACEDIM],
			       const Real time) = 0;
    //
    // Compute viscous terms.
    //
    virtual void getViscTerms (MultiFab& visc_terms,
			       int       src_comp, 
			       int       num_comp,
			       Real      time) = 0;
    //
    virtual void mac_sync () = 0;
    //
    virtual void reflux () = 0;
    //
    // Reset time levels for the initial iterations.
    //
    virtual void resetState (Real time,
	                     Real dt_old,
		             Real dt_new);

    virtual void scalar_advection_update (Real dt,
					  int  first_scalar,
					  int  last_scalar);
    //
    virtual void sum_integrated_quantities () = 0;

    virtual void velocity_diffusion_update (Real dt) = 0;

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokesBase protected static functions                         //
    ////////////////////////////////////////////////////////////////////////////

    static void variableCleanUp (); // Cleanup data descriptors at end of run.
    static void Initialize ();   // Read input file
    static void Finalize ();
    static void read_geometry (); // Set metrics

#ifdef PARTICLES
    static NSParticleContainer* theNSPC ();
    static void read_particle_params ();

    void initParticleData ();
    void post_restart_particle ();

    void post_timestep_particle (int iteration);

    virtual int timestamp_num_extras () { return 0; }
    virtual void timestamp_add_extras (int lev, Real time, MultiFab& mf) { }

    MultiFab* ParticleDerive(const std::string& name, Real time, int ngrow);
    void ParticleDerive (const std::string& name, Real time, MultiFab& mf, int dcomp);
#endif

private:

    void buildMetrics (); // 1-D metric arrays for RZ
    //
    // Compute dpdt
    //
    void calcDpdt();
    //
    // Grow by nGrow and fillpatch the MAC-projected velocities.
    //
    void create_umac_grown (int nGrow);
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
    void getOutFlowFaces (Array<Orientation>& outFaces);
    //
    // Get state data.
    //
    MultiFab* getState (int  ngrow,
                        int  state_indx,
                        int  strt_comp,
                        int  num_comp,
                        Real time);
    //
    // Compile p_avg in advance.
    //
    void incrPAvg ();
    void initOldPress (); // Initialize old pressure with new
    void zeroNewPress (); // Set new pressure to zero
    void zeroOldPress (); // Set old pressure to zero
    //
    // Compile rho_avg in advance.
    //
    void initRhoAvg (Real alpha);
    //
    void init_additional_state_types ();
    //
    Real initialTimeStep ();
    //
    void SetGodunov () { if (godunov == 0) godunov = new Godunov(); };

//
// Data
//
protected:

    PArray<FluxBoxes> raii_fbs;
    PArray<MultiFab>  raii_mfs;

    MultiFab volume;
    MultiFab area[BL_SPACEDIM];

    //
    // MAC edge velocities.
    //
    MultiFab* u_mac;
    //
    // Advective update terms.
    //
    MultiFab* aofs;

    Diffusion* diffusion;
    //
    // MultiFab arrays for variable viscosity and diffusivity
    //
    MultiFab *viscn_cc, *viscnp1_cc;
    MultiFab *diffn_cc, *diffnp1_cc;
    //
    // Sync update quantities.
    //
    // Average rho over a composite timestep.
    // Used only in the multlilevel sync projection
    //
    MultiFab rho_avg;
    //
    // Average p over a composite timestep.
    // Used only to average fine pressure data to coarse.
    //
    MultiFab p_avg;
    
    MultiFab Vsync;    // Velocity sync update storage
    MultiFab Ssync;    // Scalar sync update storage
    //
    // Density at time n+1/2 (used in advance).
    //
    MultiFab rho_half;
    MultiFab* rho_qtime;
    MultiFab* rho_tqtime;
    //
    // Density at prev_time used in advance().
    //
    MultiFab rho_ptime;
    //
    // Density at cur_time used in advance().
    //
    MultiFab rho_ctime;
    //
    // Data structure used to compute RHS for sync project.
    //
    SyncRegister* sync_reg;
    //
    // Data structures to store advective and viscous refluxing quantities 
    // on coarse-fine grid interfaces.
    //
    FluxRegister* advflux_reg;
    FluxRegister* viscflux_reg;
    //
    // Radii for r-z calculations.
    //
    Array< Array<Real> > radius;
    //
    // Flag for doing first step after regrid
    //
    bool is_first_step_after_regrid;
    //
    // Intersection of current BoxArray with that from before regrid.
    //
    BoxArray old_intersect_new;
    //
    // Controls for particle subcycling
    //
    int  umac_n_grow;
    //
    // Static objects.
    //
    static Godunov*    godunov;
    static ErrorList   err_list;
    static BCRec       phys_bc;
    static Projection* projector;
    static MacProj*    mac_projector;
    //
    // Internal parameters for timestepping.
    //
    static Real init_shrink;   // reduction factor of first esimated timestep
    static int  init_iter;     // # of timestep iterations for initial pressure
    static Real cfl;           // desired maximum cfl
    static Real change_max;    // maximum change in dt over a timestep
    static Real fixed_dt;      // set > 0 to specify dt
    static int  initial_iter;  // flag for initial pressure iterations
    static int  initial_step;  // flag for initial iterations
    static Real dt_cutoff;     // minimum dt allowed
    static int  sum_interval;  // number of timesteps for conservation stats
    static int  turb_interval; // number of timesteps for turbulence stats
    static int  jet_interval;  // number of timesteps for turbulence stats for the round jet
    static int  jet_interval_split;  // cell divisor for each dimension in numerical integration 
    //
    // Internal parameters for options.
    //
    static int        radius_grow;
    static int        verbose;
    static Real       gravity;
    static int        NUM_SCALARS;      // number of non-velocity components
    static int        NUM_STATE;        // total number of state variables
    //
    // Controls over how the convective and diffusive terms are calculated
    // for each state variable.
    //
    static Array<AdvectionForm> advectionType;
    static Array<DiffusionForm> diffusionType;
    //
    // Viscosity parameters.
    //
    static Array<int>  is_diffusive;    // does variable diffuse?
    static Array<Real> visc_coef;       // const coef viscosity terms
    static Real        visc_tol;
    static Real        visc_abs_tol;
    static Real        be_cn_theta;
    static int         variable_vel_visc;  // variable viscosity flag
    static int         variable_scal_diff; // variable scalar diffusion flag
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
    static int  do_MLsync_proj;
    static int  do_reflux;
    static int  modify_reflux_normal_vel;
    static int  do_mac_proj;
    static int  do_divu_sync;
    static int  do_refine_outflow;          // The following three items control
    static int  do_derefine_outflow;        //   how grids are treated at 
    static int  Nbuf_outflow;               //   outflow boundaries.
    static int  do_denminmax;               // The code for these was in NavierStokes.cpp, 
    static int  do_scalminmax;              //   but the flags were not declared or read in.
    static int  do_density_ref;             // Flags added for refinement control
    static int  do_tracer_ref;
    static int  do_tracer2_ref;
    static int  do_vorticity_ref;
    static int  do_scalar_update_in_order;  // Flags to allow evaluation of source terms
    static Array<int> scalarUpdateOrder;
    static int  getForceVerbose;            // Does exactly what it says on the tin
    //
    // Member when pressure defined at points in time rather than interval
    //
    static int  Dpdt_Type;
    //
    // Members for non-zero divu.
    //
    static int  additional_state_types_initialized;
    static int  Divu_Type;
    static int  Dsdt_Type;
    static int  num_state_type;
    static int  have_divu;
    static int  have_dsdt;
    static Real divu_relax_factor;
    static int  S_in_vel_diffusion;
    static int  do_init_vort_proj;
    static int  do_init_proj;
    //
    // Running statistics controls
    //
    static int  do_running_statistics;
    static Real volWgtSum_sub_origin_x;
    static Real volWgtSum_sub_origin_y;
    static Real volWgtSum_sub_origin_z;
    static Real volWgtSum_sub_Rcyl;
    static Real volWgtSum_sub_dx;
    static Real volWgtSum_sub_dy;
    static Real volWgtSum_sub_dz;
    //
    // Control velocity vs momentum update
    //
    static int  do_mom_diff;
    static int  predict_mom_together;
    static bool def_harm_avg_cen2edge;
    //
    // static integer constants
    //
    static const int GEOM_GROW  = 1;
    static const int PRESS_GROW = 1;
    static const int DIVU_GROW  = 1;
    static const int DSDT_GROW  = 1;
};

#endif
