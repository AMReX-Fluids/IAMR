#ifndef _NavierStokes_H_
#define _NavierStokes_H_ $Id: NavierStokes.H,v 1.3 1997-09-19 18:07:58 vince Exp $

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

// "Divu_Type" means S, where divergence U = S
// "Dsdt_Type" means pd S/pd t, where S is as above

// Determine what you want in the state -- Divu, Dsdt -- in
// NavierStokes::variableSetUp in NS_setup.C. Look at the example in 
// the _tube_flow problem dependent directory

enum StateType {State_Type=0, Press_Type}; 
// note: enumerated value NUM_STATE_TYPE no longer defined
// use num_state_type (static int member of NavierStokes)

#if (BL_SPACEDIM == 2)
enum StateNames  { Xvel=0, Yvel, Density};
#else
enum StateNames  { Xvel=0, Yvel, Zvel, Density};
#endif
enum PressureNames { Pressure=0 };
enum DivuNames { Divu=0};
enum DsDtNames { Dsdt=0};

class NSBld : public LevelBld
{
    virtual void variableSetUp();
    virtual void variableCleanUp();
    virtual AmrLevel *operator()();
    virtual AmrLevel *operator()(Amr &papa, int lev,
                                 const Geometry &level_geom,
                                 const BoxArray &ba, REAL time);
};

class NavierStokes : public AmrLevel {

public:   // ===================== public member functions ===

    friend class Projection;
    friend class MacProj;

    //=================================================================
    // Initialization functions follow
    //=================================================================
    NavierStokes();
    NavierStokes(Amr& papa, int lev, const Geometry &level_geom,
		 const BoxArray& bl, REAL time);
    virtual ~NavierStokes();

    virtual void restart(Amr& papa, istream &is);

    // define data descriptors
    static void variableSetUp();

    // cleanup data descriptors at end of run
    static void variableCleanUp();

    // init grid data at problem start-up
    virtual void initData();
    virtual void initDataOtherTypes() {return;} ;

    // reset time levels for the initial iterations
    virtual void resetState( REAL time, REAL dt );
    
    // set time levels of state data
    virtual void setTimeLevel(REAL time, REAL dt);
    
    // This function fills a new level n with the best
    // level n and coarser data available (usually during regrid)
    virtual void init(AmrLevel &old);

    // init data on this level after regridding if old level
    // did not exist previously
    virtual void init();

    virtual void allocOldData();
    virtual void removeOldData();

    void init_additional_state_types();

    //=================================================================
    // Diagnostics and IO functions
    //=================================================================

    // error estimation for regridding
    virtual void errorEst(TagBoxArray &tb, int clearval, int tagval,
			  REAL time);

    // integrate derived quantities over domain
    REAL sumDerive(const aString& name, REAL time);
    REAL volWgtSum(const aString& name, REAL time);

    virtual void checkPoint(ostream& os);

    // write plot file
    virtual void writePlotFile(ostream &os);
    virtual void incrPlotFile_ndata(int* n_data) { return;} ;
    virtual void writeOtherPlotFileNames(ostream &os) {return;} ;
    virtual void writeOtherPlotFileGridData(ostream &os, int igrid) 
      {return;} ;

    //=================================================================
    // Timestep estimation functions follow
    //=================================================================

    REAL estTimeStep();
    REAL initialTimeStep();
    virtual void computeInitialDt(int finest_level, int sub_cycle,
				  Array<int>& n_cycle,
				  const Array<IntVect>& ref_ratio,
				  Array<REAL>& dt_level);
    virtual void computeNewDt(int finest_level, int sub_cycle,
			      Array<int>& n_cycle,
			      const Array<IntVect>& ref_ratio,
			      Array<REAL>& dt_min,
			      Array<REAL>& dt_level,
	                      REAL stop_time);

    // this function estimates the initial timesteping used by the model
    void post_init_estDT( REAL &dt_init,
                          Array<int> &nc_save,
                          Array<REAL> &dt_save );

    // estimate the end of the simulation for amrLevel
    virtual int okToContinue();

    //=================================================================
    // the main hooks into amr and amrLevel
    //=================================================================
    
    // sync state and pressure at the end of a composite timestep
    virtual void post_timestep();
    
    // build any additional data structures after regrid
    virtual void post_regrid(int lbase, int new_finest);
    
    // build any additional data structures after restart
    virtual void post_restart();

    // insure state, and pressure are consistent
    virtual void post_init();

    // advance grids at this level in time
    virtual REAL advance(REAL time, REAL dt, int iteration, int ncycle);

    //============ public for debugging purposes
    // mac edge velocities
    MultiFab *u_mac;

    // advective tendencies
    MultiFab *aofs;

    // boundary condition functions
    void SaveOldBoundary( REAL time ); // save ghost cells back into 
    void SaveNewBoundary( REAL time ); // state (old/new) multifabs
    
    // get divU from the state data
    void getDivCond(FARRAYBOX& fab, int ngrow, REAL time);
    void getDivCond(FARRAYBOX& fab, int gridno, int ngrow, REAL time);

    // public so MacProj can access it
    virtual void SetGodunov();
    static Godunov *godunov;

protected: // ===================== protected member functions ===

    void buildMetrics();  // 1-D metric arrays for RZ
    virtual void buildRho0(); // build atmospheric rho0

    void initOldPress(); // initialize old pressure with new
    void zeroNewPress(); // set new pressure to zero
    void zeroOldPress(); // set old pressure to zero

    static void read_params();   // read input file
    static void read_geometry(); // set metrics

    virtual void sum_integrated_quantities();

    //=================================================================
    // ADVANCE FUNCTIONS
    //=================================================================

    // setup, clean up for a level timestep
    virtual void advance_setup( REAL time, REAL dt, int iteration, int ncycle);
    void advance_cleanup(REAL dt, int iteration, int ncycle);
    
    // compute divergent mac velocities, estimate best timestep
    // from t^n data, and compile cfl number
    REAL predict_velocity(REAL dt, REAL &comp_cfl );
    void test_umac_periodic();  // make sure level 0 periodic umac are
                                // periodic
    
    // make mac velocities divergence free
    void mac_project(REAL dt);

    // advect velocities
    void velocity_advection(REAL dt);

    // advect scalars
    void scalar_advection(REAL dt, int first_scalar, int last_scalar);

    // update scalars, (viscous solve in scalar_update)
    void scalar_update(REAL dt, int first_scalar, int last_scalar);
    void scalar_advection_update(REAL dt, int first_scalar, int last_scalar);
    void scalar_diffusion_update(REAL dt, int first_scalar, int last_scalar);
    void diffuse_scalar_setup(int sigma, int* rho_flag); 

    // make rho at time n+1/2
    void makerhonph(REAL dt);
    
    // update velocities after advection
    void velocity_update(REAL dt);
    void velocity_advection_update(REAL dt);
    virtual void velocity_diffusion_update(REAL dt);
    void initial_velocity_diffusion_update(REAL dt);
    void diffuse_velocity_setup(REAL dt, MultiFab*& delta_rhs); 

    // compute level projection
    void level_projector(REAL dt, REAL time, int iteration);

    //=================================================================
    // MULTILEVEL SYNC FUNCTIONS
    //=================================================================

    // insure state is consistent, ie velocity field is nondivergent,
    // coarse levels are fine level averages, pressure is zero
    virtual void post_init_state();
    
    // initialize the pressure by iterating the initial timestep
    virtual void post_init_press( REAL &dt_init,
                                  Array<int> &nc_save,
                                  Array<REAL> &dt_save );

    // interpolate cell-centered cync correction from coarse to fine
    void SyncInterp( MultiFab &CrseSync, int c_lev,
                     MultiFab &FineSync, int f_lev, IntVect& ratio,
                     int src_comp, int dest_comp, int num_comp,
                     int increment, REAL dt_clev, int** bc_orig_qty,
                     int which_interp=0);

    // two virtual hooks for anelastic SyncInterp
    virtual void ScaleCrseSyncInterp( FARRAYBOX &cdata,
                                      int c_lev, int num_comp ) {}
    virtual void reScaleFineSyncInterp( FARRAYBOX &fdata,
                                        int f_lev, int num_comp ) {}

    // bilinear interpolate nodal pressures from coarse to fine
    void SyncProjInterp( MultiFab &phi,   int c_lev,
                         MultiFab &P_new, int f_lev,
                         IntVect& ratio);
    
    // average a fine multifab down onto a coarse one
    void avgDown( const BoxArray &grids, const BoxArray &fgrids,
                  MultiFab &S_crse, MultiFab &S_fine,
                  MultiFab &volume, MultiFab &fvolume,
                  int c_level,            int f_level,
                  int strt_comp, int num_comp, IntVect& fratio );
 
    // average fine down to coarse in the ovlp intersection
    virtual void avgDown( const FARRAYBOX& fine_fab, const FARRAYBOX& crse_fab, 
                          const FARRAYBOX& fine_vol, const FARRAYBOX& crse_vol,
                          int f_level,               int c_level,
                          const BOX &ovlp,
                          int start_comp, int num_comp, IntVect& fratio );

    // inject fine pressure nodes down onto coarse nodes
    void injectDown( const BOX &ovlp, FARRAYBOX &Pcrse,
                     const FARRAYBOX &Pfine, IntVect &fine_ratio );

    // test for consistency between fine dirichlet and coarse nodes
    void testInject( const BOX &ovlp, FARRAYBOX &Pcrse,
                     const FARRAYBOX &Pfine, IntVect &fine_ratio );

    // compute the level sync correction
    virtual void level_sync();
    
    // compute the mac sync correction
    virtual void mac_sync();
    
    // reflux function
    virtual void reflux();

    void avgDown();         // average down for all the state types
    void avgDown(int comp); // for a single scalar

    // compile rho_avg in advance
    void initRhoAvg(REAL alpha);
    void incrRhoAvg(REAL alpha);

    // compile p_avg in advance
    void incrPAvg( int iteration, REAL alpha); 

    //=================================================================
    // ACCESS FUNCTIONS FOLLOW
    //=================================================================
    
    // virtual function get fluxes out of the advection routines
    virtual void pullFluxes( int level, int gridno, int start_ind, int ncomp,
                             FARRAYBOX &xflux,
                             FARRAYBOX &yflux,
                             FARRAYBOX &zflux, REAL dt );

    // compute the gradient of divU - dSdt
    void compute_grad_divu_minus_s(REAL time, MultiFab* grad_divu_minus_s,
                                   int scaleRhoDivDt);
    
    // compute viscous terms
    virtual void getViscTerms(MultiFab& visc_terms, int src_comp, 
                              int num_comp, REAL time);

    // calculate divU, which sets them to zero by default
    virtual void calc_divu(REAL time, REAL dt, MultiFab& fab);

    // calculate dSdT, which sets them to zero by default
    virtual void calc_dsdt(REAL time, REAL dt, MultiFab& fab);

    // virtual function to allow arbitrary forcing
    virtual void getForce(FARRAYBOX& force, int gridno, int ngrow,
                          int strt_comp, int num_comp, REAL time);
    
    // get state data
    void getState(FARRAYBOX& fab,int gridno, int ngrow,
                  int strt_comp, int num_comp, REAL time);

    void getState(FARRAYBOX& fab,int gridno, int ngrow,
                  int state_indx, int strt_comp, int num_comp, REAL time);
    
    void getState(FARRAYBOX& fab, int gridno, int ngrow, 
                  REAL time, int have_state, int state_indx,
                  //BoxAssoc& assoc, Array<BOX>& unfilled,
                  Array<BOX>& unfilled,
                  int num_comp = 1, int start_comp = 0);
    
    // get pressure gradient data via fill-patching
    void getGradP(FARRAYBOX& gp, int gridno, int ngrow, REAL time);

    // given a nodal pressure P compute the pressure gradient at the
    // contained cell centers
    void getGradP(FARRAYBOX &p_fab, FARRAYBOX& gp,
                  const BOX &grd, int ngrow );
    
    // get dSdt from the state data
    void getDsdt(FARRAYBOX& fab, int gridno, int ngrow, REAL time);

    // fill ghost cells of state
    void FillStateBndry(REAL time, int state_indx, int src_comp, 
                            int num_comp); 


    NavierStokes& getLevel(int lev)
    {return *(NavierStokes*) &parent->getLevel(lev); }

    Diffusion& getDiffusion() {
	assert( diffusion );
	return *diffusion;
    }

    SyncRegister& getSyncReg() {
	assert( sync_reg );
	return *sync_reg;
    }

    FluxRegister& getAdvFluxReg() {
	assert( advflux_reg );
	return *advflux_reg;
    }

    FluxRegister& getAdvFluxReg(int lev) {
	return getLevel(lev).getAdvFluxReg();
    }

    FluxRegister& getViscFluxReg() {
	assert( viscflux_reg );
	return *viscflux_reg;
    }

    FluxRegister& getViscFluxReg(int lev) {
	return getLevel(lev).getViscFluxReg();
    }

protected: // ===================== protected data members ===

    // ----------------------- static objects 
    // In the end, there can be only one. (Louis Howellism)
    static ErrorList   err_list;
    static BCRec       phys_bc;
    static Projection *projector;
    static MacProj    *mac_projector;
    Diffusion         *diffusion; // not static
    
    // volume and area fractions
    MultiFab volume;
    MultiFab area[BL_SPACEDIM];

    // ======sync update quantities
    // average rho over a composite timestep
    // used only in the multlilevel sync projection
    MultiFab *rho_avg;

    // average p over a composite timestep
    // used only to average fine pressure data to coarse
    MultiFab *p_avg;
    
    MultiFab *Vsync;    // velocity sync update storage
    MultiFab *Ssync;    // scalar sync update storage
    
    // density at time n+1/2 (used in advance)
    MultiFab  *rho_half;

    // data structure used to compute RHS for sync project
    SyncRegister *sync_reg;
    
    // data structures to store advective and viscous refluxing quantities 
    // on coarse-fine grid interfaces
    FluxRegister *advflux_reg;
    FluxRegister *viscflux_reg;

    // Radii for r-z calculations
    PArray<REAL> radius;
    
    // box associations (nearest neighbor lists) used by MultiFab
    // for fast filpatch operations
    //BoxAssoc hyp_assoc;       // 3 zone wide cell centered association
    //BoxAssoc divu_assoc;     // 1 zone wide cell centered association
    //BoxAssoc dsdt_assoc;     // 1 zone wide cell centered association
    
    // for fast filpatch, we precompute the subbox of a grid that
    // must be interpolated from coarser data.
    Array<BOX>  hyp_unfilled;      // unfilled region for Godunov step
    Array<BOX>  cc1_unfilled;      // unfilled region for Godunov step
    Array<BOX>  divu_unfilled;     // unfilled region for divu fill
    Array<BOX>  dsdt_unfilled;     // unfilled region for dsdt fill

    // build all the (*_unfilled) arrays
    void buildUnfilledRegions();

    // ----------------------- static internal parameters for timestepping
    static REAL init_shrink;   // reduction factor of first esimated timestep
    static int  init_iter;     // # of timestep iterations for initial pressure
    static REAL cfl;           // desired maximum cfl
    static REAL change_max;    // maximum change in dt over a timestep
    static REAL fixed_dt;      // set > 0 to specify dt
    static int  initial_iter;  // flag for initial pressure iterations
    static int  initial_step;  // flag for initial iterations
    static REAL dt_cutoff;     // minimum dt allowed
    static int sum_interval;   // number of timesteps for conservation stats
    

    // ----------------------- static internal parameters for options
    
    static int radius_grow;
    static int  verbose;
    static REAL gravity;
    static int  NUM_SCALARS;       // number of non-velocity components
    static int  NUM_STATE;         // total number of state variables
    static Array<int> is_conservative;  // conservative or convective?

    // ----------------------- viscosity parameters
    static Array<int> is_diffusive;     // does variable diffuse?
    static Array<REAL> visc_coef;       // const coef viscosity terms
    static REAL visc_tol;
    static REAL visc_abs_tol;
    static REAL be_cn_theta;

    // ----------------------- internal switches
    static int  do_sync_proj;
    static int  do_MLsync_proj;
    static int  do_reflux;
    static int  do_mac_proj;
    static int  do_diffusion;

    // ----------------------- static new members for non-zero divu
    static int additional_state_types_initialized;
    static int Divu_Type;
    static int Dsdt_Type;
    static int num_state_type;
    static int have_divu;
    static int have_dsdt;
    static REAL divu_minus_s_factor;
    static REAL divu_relax_factor;
    static int S_in_vel_diffusion;
};
#endif
