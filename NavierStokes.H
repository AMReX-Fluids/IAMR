//BL_COPYRIGHT_NOTICE

#ifndef _NavierStokes_H_
#define _NavierStokes_H_

//
// $Id: NavierStokes.H,v 1.31 1999-02-25 17:35:37 car Exp $
//

#include <LevelBld.H>
#include <BC_TYPES.H>
#include <AmrLevel.H>
#include <ErrorList.H>
#include <SyncRegister.H>
#include <FluxRegister.H>
#include <ViscBndry.H>
#include <Projection.H>
#include <MacProj.H> 
#include <Godunov.H>
#include <Diffusion.H>

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

    NavierStokes ();

    NavierStokes (Amr&            papa,
                  int             lev,
                  const Geometry& level_geom,
                  const BoxArray& bl,
                  Real            time);

    virtual ~NavierStokes ();

    virtual void restart (Amr&     papa,
                          istream& is,
                          bool     bReadSpecial = false);
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
    Real sumDerive (const aString& name,
                    Real           time);

    Real volWgtSum (const aString& name,
                    Real           time);
    //
    // A string written as the first item in writePlotFile() at level zero.
    // It is so we can distinguish between different types of plot files.
    //
    // For NavierStokes it has the form: NavierStokes-Vnnn
    //
    virtual aString thePlotFileType () const;
    //
    // Write plot file stuff to specified directory.
    //
    virtual void writePlotFile (const aString& dir,
                                ostream&       os,
                                VisMF::How     how);
    //
    // Returns count of additional components of StateData MultiFabs
    // derived classes want to write out to plotfiles.
    //
    virtual int additionalStateDataItemsCount () const;
    //
    // Returns names of the additional components of StateData MultiFabs
    // derived classes want to write out to plotfiles.
    //
    virtual Array<aString> additionalStateDataItemsNames () const;
    //
    // Pointers to additional StateData MultiFabs derived
    // classes want to write out to plotfiles.
    //
    virtual Array<const MultiFab*> additionalStateDataItems () const;
    //
    // For each MultiFab in additionalStateDataItems(), this is a map
    // of those components of the respective MultiFab that should be written.
    //
    virtual Array< vector<int> > additionalStateDataItemsMap () const;
    //
    // Timestep estimation functions follow ...
    //
    Real estTimeStep ();

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
    virtual void post_timestep ();
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
    // Boundary condition functions.
    //
    void SaveOldBoundary (Real time); // Save ghost cells back into.
    void SaveNewBoundary (Real time); // State (old/new) multifabs.
    //
    // Get divU from the state data.
    //
    void getDivCond (FArrayBox& fab,
                     int        ngrow,
                     Real       time);

    void getDivCond (FArrayBox& fab,
                     int        gridno,
                     int        ngrow,
                     Real       time);

    MultiFab* getDivCond (int  ngrow,
                          Real time);
    //
    // Public so MacProj can access it.
    //
    virtual void SetGodunov ();

protected:

    void buildMetrics ();      // 1-D metric arrays for RZ
    virtual void buildRho0 (); // Build atmospheric rho0

    void initOldPress (); // Initialize old pressure with new
    void zeroNewPress (); // Set new pressure to zero
    void zeroOldPress (); // Set old pressure to zero

    static void read_params ();   // Read input file
    static void read_geometry (); // Set metrics

    virtual void sum_integrated_quantities ();
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

    void test_umac_periodic ();  // Ensure level 0 periodic umac are periodic.
    //
    // Make mac velocities divergence free.
    //
    void mac_project (Real dt);
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

    void diffuse_scalar_setup (int  sigma,
                               int* rho_flag); 
    //
    // Make rho at time n+1/2.
    //
    void makerhonph (Real dt);
    //
    // Update velocities after advection.
    //
    void velocity_update (Real dt);

    void velocity_advection_update (Real dt);

    virtual void velocity_diffusion_update (Real dt);

    void initial_velocity_diffusion_update (Real dt);

    void diffuse_velocity_setup (Real       dt,
                                 MultiFab*& delta_rhs); 
    //
    // Compute level projection.
    //
    void level_projector (Real dt,
                          Real time,
                          int  iteration);
    //
    // Insure state is consistent, ie velocity field is nondivergent,
    // coarse levels are fine level averages, pressure is zero.
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
    void SyncInterp (MultiFab& CrseSync,
                     int       c_lev,
                     MultiFab& FineSync,
                     int       f_lev,
                     IntVect&  ratio,
                     int       src_comp,
                     int       dest_comp,
                     int       num_comp,
                     int       increment,
                     Real      dt_clev,
                     int**     bc_orig_qty,
                     int       which_interp = 0);
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
                         int       f_lev,
                         IntVect&  ratio);
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
                  IntVect&        fratio);
    //
    // Average fine down to coarse in the ovlp intersection.
    //
    virtual void avgDown (const FArrayBox& fine_fab,
                          const FArrayBox& crse_fab, 
                          const FArrayBox& fine_vol,
                          const FArrayBox& crse_vol,
                          int              f_level,
                          int              c_level,
                          const Box&       ovlp,
                          int              start_comp,
                          int              num_comp,
                          IntVect&         fratio);
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
    virtual void level_sync ();
    //
    // Compute the mac sync correction.
    //
    virtual void mac_sync ();
    //
    // Reflux function.
    //
    virtual void reflux ();

    void avgDown ();         // Average down for all the state types.
    void avgDown (int comp); // Average down forfor a single scalar.
    //
    // Compile rho_avg in advance.
    //
    void initRhoAvg (Real alpha);
    void incrRhoAvg (Real alpha);
    //
    // Compile p_avg in advance.
    //
    void incrPAvg (int  iteration,
                   Real alpha); 
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
    // Compute the gradient of divU - dSdt.
    //
    void compute_grad_divu_minus_s (Real      time,
                                    MultiFab* grad_divu_minus_s,
                                    int       scaleRhoDivDt);
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
    // Virtual function to allow arbitrary forcing.
    //
    virtual void getForce (FArrayBox& force,
                           int        gridno,
                           int        ngrow,
                           int        strt_comp,
                           int        num_comp,
                           Real       time);

    void getForce (FArrayBox&       force,
                   int              gridno,
                   int              ngrow,
                   int              strt_comp,
                   int              num_comp,
                   const FArrayBox& Rho);
    //
    // Get state data.
    //
    void getState (FArrayBox& fab,
                   int        gridno,
                   int        ngrow,
                   int        strt_comp,
                   int        num_comp,
                   Real       time);

    void getState (FArrayBox& fab,
                   int        gridno,
                   int        ngrow,
                   int        state_indx,
                   int        strt_comp,
                   int        num_comp,
                   Real       time);

    MultiFab* getState (int  ngrow,
                        int  state_indx,
                        int  strt_comp,
                        int  num_comp,
                        Real time);
    
    void getState (FArrayBox&  fab,
                   int         gridno,
                   int         ngrow, 
                   Real        time,
                   int         have_state,
                   int         state_indx,
                   Array<Box>& unfilled,
                   int         start_comp = 0,
                   int         num_comp = 1);
    //
    // Get pressure gradient data via fill-patching.
    //
    void getGradP (FArrayBox& gp,
                   int        gridno,
                   int        ngrow,
                   Real       time);
    //
    // Given a nodal pressure P compute the pressure gradient at the
    // contained cell centers.
    //
    void getGradP (FArrayBox& p_fab,
                   FArrayBox& gp,
                   const Box& grd,
                   int        ngrow);
    //
    // Get dSdt from the state data.
    //
    void getDsdt (FArrayBox& fab,
                  int        gridno,
                  int        ngrow,
                  Real       time);

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
    void calcViscosity(Real time,
                       const MultiFab& state,
                       MultiFab* visc[BL_SPACEDIM]);

    void calcDiffusivity(Real time,
                         const MultiFab& state,
                         MultiFab* diff[BL_SPACEDIM],
                         int sComp = 0,
                         int num_comp = 1);

    void getViscosity(MultiFab* viscosity[BL_SPACEDIM],
                      Real time);

    void getDiffusivity(MultiFab* diffusivity[BL_SPACEDIM],
                        Real time,
                        int sComp);


    NavierStokes& getLevel (int lev)
    {
        return *(NavierStokes*) &parent->getLevel(lev);
    }

    Diffusion& getDiffusion ()
    {
        assert(diffusion);
        return *diffusion;
    }

    SyncRegister& getSyncReg ()
    {
        assert(sync_reg);
        return *sync_reg;
    }

    FluxRegister& getAdvFluxReg ()
    {
        assert(advflux_reg);
        return *advflux_reg;
    }

    FluxRegister& getAdvFluxReg (int lev)
    {
        return getLevel(lev).getAdvFluxReg();
    }

    FluxRegister& getViscFluxReg ()
    {
        assert(viscflux_reg);
        return *viscflux_reg;
    }

    FluxRegister& getViscFluxReg (int lev)
    {
        return getLevel(lev).getViscFluxReg();
    }

protected:
    //
    // Mac edge velocities.
    //
    MultiFab* u_mac;
    //
    // Advective tendencies.
    //
    MultiFab* aofs;

    Diffusion* diffusion;
    //
    // MultiFab arrays for variable viscosity and diffusivity
    //
    MultiFab **viscn, **viscnp1;
    MultiFab **diffn, **diffnp1;
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
    PArray<Real> radius;
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
    //
    // Internal parameters for options.
    //
    static int        radius_grow;
    static int        verbose;
    static Real       gravity;
    static int        NUM_SCALARS;      // number of non-velocity components
    static int        NUM_STATE;        // total number of state variables
    static Array<int> is_conservative;  // conservative or convective?
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
    static int Temp;
    static int do_temp;
    static int do_sync_proj;
    static int do_MLsync_proj;
    static int do_reflux;
    static int do_mac_proj;
    static int check_umac_periodicity;
    //
    // Members for non-zero divu.
    //
    static int  additional_state_types_initialized;
    static int  Divu_Type;
    static int  Dsdt_Type;
    static int  num_state_type;
    static int  have_divu;
    static int  have_dsdt;
    static Real divu_minus_s_factor;
    static Real divu_relax_factor;
    static int  S_in_vel_diffusion;
    static int  do_init_vort_proj;
};
#endif /*_NavierStokes_H_*/
