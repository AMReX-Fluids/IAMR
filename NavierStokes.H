
#ifndef _NavierStokes_H_
#define _NavierStokes_H_

//
// $Id: NavierStokes.H,v 1.96 2007-11-07 20:58:27 aaspden Exp $
//

#include <LevelBld.H>
#include <BC_TYPES.H>
#include <AmrLevel.H>
#include <ErrorList.H>
#include <SyncRegister.H>
#include <FluxRegister.H>
#include <ViscBndry.H>
#include <Projection.H>
#include <Godunov.H>
#include <Diffusion.H>
#include <MacProj.H> 

//
// "Divu_Type" means S, where divergence U = S
// "Dsdt_Type" means pd S/pd t, where S is as above
//
// Determine what you want in the state -- Divu, Dsdt -- in
// NavierStokes::variableSetUp in NS_setup.C. Look at the example in 
// the _tube_flow problem dependent directory
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

class NSBld
    :
    public LevelBld
{
    virtual void variableSetUp ();
    virtual void variableCleanUp ();
    virtual AmrLevel *operator() ();
    virtual AmrLevel *operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& level_geom,
                                  const BoxArray& ba,
                                  Real            time);
};

class NavierStokes
    :
    public AmrLevel 
{
public:

    friend class Projection;
    friend class MacProj;
    friend class Diffusion;

    NavierStokes ();

    NavierStokes (Amr&            papa,
                  int             lev,
                  const Geometry& level_geom,
                  const BoxArray& bl,
                  Real            time);

    virtual ~NavierStokes ();

    virtual void restart (Amr&          papa,
                          std::istream& is,
                          bool          bReadSpecial = false);
    //
    // Define data descriptors.
    //
    static void variableSetUp ();
    //
    // Cleanup data descriptors at end of run.
    //
    static void variableCleanUp ();
    //
    // Init grid data at problem start-up.
    //
    virtual void initData ();
    virtual void initDataOtherTypes () {}
    //
    // Reset time levels for the initial iterations.
    //
    virtual void resetState (Real time,
                             Real dt_old,
                             Real dt_new);
    //
    // Set time levels of state data.
    //
    virtual void setTimeLevel (Real time,
                               Real dt_old,
                               Real dt_new);
    //
    // This function fills a new level n with the best
    // level n and coarser data available (usually during regrid).
    //
    virtual void init (AmrLevel& old);
    //
    // Init data on this level after regridding if old level
    // did not exist previously.
    //
    virtual void init ();

    virtual void allocOldData ();

    virtual void removeOldData ();

    void init_additional_state_types ();
    //
    // Error estimation for regridding.
    //
    virtual void errorEst (TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           Real         time,
                           int          n_error_buf = 0,
                           int          ngrow = 0);
    //
    // Integrate derived quantities over domain.
    //
    Real sumDerive (const std::string& name,
                    Real           time);

    Real MaxVal (const std::string& name,
                 Real           time);

    Real volWgtSum (const std::string& name,
                    Real           time);

#if (BL_SPACEDIM == 3)
    void TurbSum   (Real time, 
		    Real *tmp,
		    int  ksize,
		    int  turbVars);
#ifdef SUMJET
    void JetSum   (Real time, 
		   Real *jetData,
		   int  levRsize,
		   int  levKsize,
		   int  rsize,
		   int  ksize,
		   int  jetVars);
#endif
#endif
    //
    // A string written as the first item in writePlotFile() at level zero.
    // It is so we can distinguish between different types of plot files.
    //
    // For NavierStokes it has the form: NavierStokes-Vnnn
    //
    virtual std::string thePlotFileType () const;

    //
    // Write plot file stuff to specified directory.
    //
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&  os,
                                VisMF::How     how);
    //
    // Modify list of variables to be plotted
    //
    virtual void setPlotVariables();

    //
    // Timestep estimation functions follow ...
    //
    virtual Real estTimeStep ();

    Real initialTimeStep ();

    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   Array<int>&           n_cycle,
                                   const Array<IntVect>& ref_ratio,
                                   Array<Real>&          dt_level,
                                   Real                  stop_time);

    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               Array<int>&           n_cycle,
                               const Array<IntVect>& ref_ratio,
                               Array<Real>&          dt_min,
                               Array<Real>&          dt_level,
                               Real                  stop_time);
    //
    // This function estimates the initial timesteping used by the model.
    //
   void post_init_estDT (Real&        dt_init,
                         Array<int>&  nc_save,
                         Array<Real>& dt_save,
                         Real         stop_time);
    //
    // Estimate the end of the simulation for amrLevel.
    //
    virtual int okToContinue ();
    //
    // Sync state and pressure at the end of a composite timestep.
    //
    virtual void post_timestep (int iteration);
    //
    // Build any additional data structures after regrid.
    //
    virtual void post_regrid (int lbase, int new_finest);
    //
    // Build any additional data structures after restart.
    //
    virtual void post_restart ();
    //
    // Insure state, and pressure are consistent.
    //
    virtual void post_init (Real stop_time);
    //
    // Advance grids at this level in time.
    //
    virtual Real advance (Real time,
                          Real dt,
                          int  iteration,
                          int  ncycle);
    //
    // Get divU from the state data.
    //
    MultiFab* getDivCond (int  ngrow,
                          Real time);
    //
    // Public so MacProj can access it.
    //
    virtual void SetGodunov ();

    //
    // Returns the value of "gravity" for use in the projection outflow bcs.
    //
    Real getGravity () 
    {
      return gravity;
    }

protected:

    void buildMetrics ();      // 1-D metric arrays for RZ

    void initOldPress (); // Initialize old pressure with new
    void zeroNewPress (); // Set new pressure to zero
    void zeroOldPress (); // Set old pressure to zero

    static void read_params ();   // Read input file
    static void read_geometry (); // Set metrics

    virtual void sum_integrated_quantities ();

#if (BL_SPACEDIM == 3)
    virtual void sum_turbulent_quantities ();
#ifdef SUMJET
    virtual void sum_jet_quantities ();
#endif
#endif
    //
    // Setup, clean up for a level timestep.
    //
    virtual void advance_setup (Real time,
                                Real dt,
                                int  iteration,
                                int  ncycle);

    void advance_cleanup (Real dt,
                          int  iteration,
                          int  ncycle);
    //
    // Compute divergent mac velocities, estimate best timestep
    // from t^n data, and compile cfl number.
    //
    Real predict_velocity (Real  dt,
                           Real& comp_cfl);

    //
    // Create the RHS for the MAC projection.
    //
    virtual MultiFab* create_mac_rhs (Real time, Real dt);

    //
    // Create a grown-by-1 RHS for the MAC projection, which
    //   is identical to mac_rhs in all the valid regions.
    //
    MultiFab* create_mac_rhs_grown (int nGrow, Real time, Real dt);

    //
    // Impose divergence constraint on MAC velocities.
    //
    void mac_project (Real time, Real dt, MultiFab& S_old, 
                      MultiFab* divu, int have_divu);
    //
    // Grow by 1 and fillpatch the MAC-projected velocities.
    //
    void create_umac_grown ();
    //
    // Advect velocities.
    //
    void velocity_advection (Real dt);
    //
    // Advect scalars.
    //
    void scalar_advection (Real dt,
                           int  first_scalar,
                           int  last_scalar);
    //
    // Update scalars, (viscous solve in scalar_update).
    //
    void scalar_update (Real dt,
                        int  first_scalar,
                        int  last_scalar);

    void scalar_advection_update (Real dt,
                                  int  first_scalar,
                                  int  last_scalar);

    void scalar_diffusion_update (Real dt,
                                  int  first_scalar,
                                  int  last_scalar);

    virtual void diffuse_scalar_setup (Real        dt,
                                       int         sigma,
                                       int*        rho_flag,
                                       MultiFab*&  delta_rhs,
                                       MultiFab*&  alpha,
                                       MultiFab**& diffn,
                                       MultiFab**& diffnp1); 
    //
    // Make rho at time n.
    //
    virtual void make_rho_prev_time ();
    //
    // Make rho at time n+1.
    //
    virtual void make_rho_curr_time ();
    //
    // Get rho at time n+1/2 -- don't delete the pointer !!!
    //
    MultiFab* get_rho_half_time ();
    //
    // Get rho at time.
    // Must be one of AmrOldTime, AmrHalfTime or AmrNewTime.
    // Aborts if not one of the three above times.
    //
    const MultiFab& get_rho (Real time);
    //
    // Update velocities (or momentum) after advection.
    //
    void velocity_update (Real dt);

    void velocity_advection_update (Real dt);

    virtual void velocity_diffusion_update (Real dt);

    void initial_velocity_diffusion_update (Real dt);

    virtual void diffuse_velocity_setup (Real       dt,
                                         MultiFab*& delta_rhs,
                                         MultiFab**& viscn,
                                         MultiFab**& viscnp1); 
    //
    // Compute level projection.
    //
    void level_projector (Real dt,
                          Real time,
                          int  iteration);
    //
    // Compute dpdt
    //
    void calcDpdt();

    //
    // Ensure state is consistent, i.e. velocity field is nondivergent,
    // coarse level data are averages of fine level data, pressure is zero.
    //
    virtual void post_init_state ();
    //
    // Initialize the pressure by iterating the initial timestep.
    //
    virtual void post_init_press (Real&        dt_init,
                                  Array<int>&  nc_save,
                                  Array<Real>& dt_save);
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
    // Two virtual hooks for anelastic SyncInterp.
    //
    virtual void ScaleCrseSyncInterp (FArrayBox& cdata,
                                      int        c_lev,
                                      int        num_comp) {}

    virtual void reScaleFineSyncInterp (FArrayBox& fdata,
                                        int        f_lev,
                                        int        num_comp) {}
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
    //
    // Average a fine multifab down onto a coarse one.
    //
    void avgDown (const BoxArray& grids,
                  const BoxArray& fgrids,
                  MultiFab&       S_crse,
                  MultiFab&       S_fine,
                  MultiFab&       volume,
                  MultiFab&       fvolume,
                  int             c_level,
                  int             f_level,
                  int             strt_comp,
                  int             num_comp,
                  const IntVect&  fratio);
    //
    // Average fine down to coarse in the ovlp intersection.
    //
    void avgDown (const FArrayBox& fine_fab,
                  const FArrayBox& crse_fab, 
                  const FArrayBox& fine_vol,
                  const FArrayBox& crse_vol,
                  int              f_level,
                  int              c_level,
                  const Box&       ovlp,
                  int              start_comp,
                  int              num_comp,
                  const IntVect&   fratio);
    //
    // Average fine down to coarse in the ovlp intersection.
    //
    static void avgDown_doit (const FArrayBox& fine_fab,
                              const FArrayBox& crse_fab, 
                              const FArrayBox& fine_vol,
                              const FArrayBox& crse_vol,
                              int              f_level,
                              int              c_level,
                              const Box&       ovlp,
                              int              start_comp,
                              int              num_comp,
                              const IntVect&   fratio);
    //
    // Inject fine pressure nodes down onto coarse nodes.
    //
    void injectDown (const Box&       ovlp,
                     FArrayBox&       Pcrse,
                     const FArrayBox& Pfine,
                     IntVect&         fine_ratio);
    //
    // Test for consistency between fine dirichlet and coarse nodes.
    //
    void testInject (const Box&       ovlp,
                     FArrayBox&       Pcrse,
                     const FArrayBox& Pfine,
                     IntVect&         fine_ratio);
    //
    // Compute the level sync correction.
    //
    virtual void level_sync (int crse_iteration);
    //
    // Compute the mac sync correction.
    //
    virtual void mac_sync ();
    virtual void sync_setup (MultiFab*& DeltaSsync);
    virtual void sync_cleanup (MultiFab*& DeltaSsync);

    //
    // Reflux function.
    //
    virtual void reflux ();

    virtual void avgDown (); // Average down for all the state types.
    void avgDown (int comp); // Average down for a single StateType scalar
    //
    // Compile rho_avg in advance.
    //
    void initRhoAvg (Real alpha);
    void incrRhoAvg (Real alpha);
    void incrRhoAvg (const MultiFab& rho_incr,
                     int             sComp,
                     Real            alpha);
    //
    // Compile p_avg in advance.
    //
    void incrPAvg ();
    //
    // Virtual function get fluxes out of the advection routines.
    //
    virtual void pullFluxes (int        gridno,
                             int        start_ind,
                             int        ncomp,
                             FArrayBox& xflux,
                             FArrayBox& yflux,
                             FArrayBox& zflux,
                             Real       dt);
    //
    // Compute viscous terms.
    //
    virtual void getViscTerms (MultiFab& visc_terms,
                               int       src_comp, 
                               int       num_comp,
                               Real      time);
    //
    // Calculate divU, which sets them to zero by default.
    //
    virtual void calc_divu (Real      time,
                            Real      dt,
                            MultiFab& fab);
    //
    // Calculate dSdT, which sets them to zero by default.
    //
    virtual void calc_dsdt (Real      time,
                            Real      dt,
                            MultiFab& fab);
    //
    // Get the forcing term.
    //
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
                   const FArrayBox& Rho);
#endif		 
    //
    // Get state data.
    //
    MultiFab* getState (int  ngrow,
                        int  state_indx,
                        int  strt_comp,
                        int  num_comp,
                        Real time);
    //
    // Get pressure gradient data via fill-patching.
    //
    void getGradP (MultiFab& gp,
                   Real      time);
    //
    // Get dSdt from the state data.
    //
    MultiFab* getDsdt (int  ngrow,
                       Real time);
    //
    // Fill ghost cells of state.
    //
    void FillStateBndry (Real time,
                         int  state_indx,
                         int  src_comp, 
                         int  num_comp); 
    //
    // Calculate nonuniform viscosity and diffusivity
    //
    virtual void calcViscosity (const Real time,
                                const Real dt,
                                const int  iteration,
                                const int  ncycle);

    virtual void calcDiffusivity (const Real time,
                                  const Real dt,
                                  const int  iteration,
                                  const int  ncycle,
                                  const int  src_comp = 0,
                                  const int  num_comp = 1);

    virtual void getViscosity (MultiFab*  viscosity[BL_SPACEDIM],
                               const Real time);

    virtual void getDiffusivity (MultiFab*  diffusivity[BL_SPACEDIM],
                                 const Real time,
                                 const int  state_comp,
                                 const int  dst_comp,
                                 const int  num_comp);

    void center_to_edge_plain (const FArrayBox& ccfab,
                               FArrayBox&       ecfab,
                               int              sComp,
                               int              dComp,
                               int              nComp);
    //
    // Called in grid_places after other tagging routines to modify
    //   the list of tagged points
    //
    virtual void manual_tags_placement (TagBoxArray&    tags, 
                                        Array<IntVect>& bf_lev);
    NavierStokes& getLevel (int lev)
    {
        return *(NavierStokes*) &parent->getLevel(lev);
    }

    Diffusion& getDiffusion ()
    {
        BL_ASSERT(diffusion);
        return *diffusion;
    }

    SyncRegister& getSyncReg ()
    {
        BL_ASSERT(sync_reg);
        return *sync_reg;
    }

    FluxRegister& getAdvFluxReg ()
    {
        BL_ASSERT(advflux_reg);
        return *advflux_reg;
    }

    FluxRegister& getAdvFluxReg (int lev)
    {
        return getLevel(lev).getAdvFluxReg();
    }

    FluxRegister& getViscFluxReg ()
    {
        BL_ASSERT(viscflux_reg);
        return *viscflux_reg;
    }

    FluxRegister& getViscFluxReg (int lev)
    {
        return getLevel(lev).getViscFluxReg();
    }

protected:
    //
    // MAC edge velocities.
    //
    MultiFab* u_mac;
    MultiFab* u_macG;
    MultiFab* u_corr;
    //
    // Advective update terms.
    //
    MultiFab* aofs;

    //
    // RHS for the MAC projection.
    //
    MultiFab* mac_rhs;

    Diffusion* diffusion;
    //
    // MultiFab arrays for variable viscosity and diffusivity
    //
    MultiFab *viscn_cc, *viscnp1_cc;
    MultiFab *diffn_cc, *diffnp1_cc;
    //
    // Volume and area fractions.
    //
    MultiFab volume;
    MultiFab area[BL_SPACEDIM];
    //
    // Sync update quantities.
    //
    // Average rho over a composite timestep.
    // Used only in the multlilevel sync projection
    //
    MultiFab* rho_avg;
    //
    // Average p over a composite timestep.
    // Used only to average fine pressure data to coarse.
    //
    MultiFab* p_avg;
    
    MultiFab* Vsync;    // Velocity sync update storage
    MultiFab* Ssync;    // Scalar sync update storage
    //
    // Density at time n+1/2 (used in advance).
    //
    MultiFab* rho_half;
    MultiFab* rho_qtime;
    MultiFab* rho_tqtime;
    //
    // Density at prev_time used in advance().
    //
    MultiFab* rho_ptime;
    //
    // Density at cur_time used in advance().
    //
    MultiFab* rho_ctime;
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
    //
    // Control velocity vs momentum update
    //
    static int  do_mom_diff;
    static int  predict_mom_together;
    static bool def_harm_avg_cen2edge;

};

#endif /*_NavierStokes_H_*/
