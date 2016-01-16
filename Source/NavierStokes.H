
#ifndef _NavierStokes_H_
#define _NavierStokes_H_

#include <NavierStokesBase.H>
#include <LevelBld.H>

class NSBld
    :
    public LevelBld
{
    virtual void variableSetUp () BL_OVERRIDE;
    virtual void variableCleanUp () BL_OVERRIDE;
    virtual AmrLevel *operator() () BL_OVERRIDE;
    virtual AmrLevel *operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& level_geom,
                                  const BoxArray& ba,
                                  Real            time) BL_OVERRIDE;
};

class NavierStokes
    :
    public NavierStokesBase
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

    //
    // Define data descriptors.
    //
    static void variableSetUp ();
    //
    // Cleanup data descriptors at end of run.
    //
    static void variableCleanUp ();

    //
    // Define error estimators
    //
    static void error_setup ();

    //
    // Init grid data at problem start-up.
    //
    virtual void initData () BL_OVERRIDE;
    void initDataOtherTypes () { }

    virtual void allocOldData () BL_OVERRIDE;

    virtual void removeOldData () BL_OVERRIDE;
    //
    // Integrate derived quantities over domain.
    //
    Real sumDerive (const std::string& name,
                    Real           time);

    Real MaxVal (const std::string& name,
                 Real           time);

#if (BL_SPACEDIM == 3)
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
    // Write plot file stuff to specified directory.
    //
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&  os,
                                VisMF::How     how) BL_OVERRIDE;
    //
    // Modify list of variables to be plotted
    //
    virtual void setPlotVariables() BL_OVERRIDE;
    //
    // Returns a MultiFab containing the derived data for this level.
    // The user is responsible for deleting this pointer when done
    // with it.  If ngrow>0 the MultiFab is built on the appropriately
    // grown BoxArray.
    //
    virtual MultiFab* derive (const std::string& name,
                              Real               time,
                              int                ngrow) BL_OVERRIDE;
    //
    // This version of derive() fills the dcomp'th component of mf with the derived quantity.
    //
    virtual void derive (const std::string& name,
                         Real               time,
                         MultiFab&          mf,
                         int                dcomp) BL_OVERRIDE;
    //
    // Insure state, and pressure are consistent.
    //
    virtual void post_init (Real stop_time) BL_OVERRIDE;
    //
    // Advance grids at this level in time.
    //
    virtual Real advance (Real time,
                          Real dt,
                          int  iteration,
                          int  ncycle) BL_OVERRIDE;

protected:

    static void Initialize ();   // Read input file
    static void Finalize ();

    virtual void sum_integrated_quantities () BL_OVERRIDE;

#if (BL_SPACEDIM == 3)
#ifdef SUMJET
    void sum_jet_quantities ();
#endif
#endif

    //
    // Compute divergent mac velocities, estimate best timestep
    // from t^n data, and compile cfl number.
    //
    Real predict_velocity (Real dt, Real& comp_cfl);
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

    void scalar_diffusion_update (Real dt,
                                  int  first_scalar,
                                  int  last_scalar);

    virtual void velocity_diffusion_update (Real dt) BL_OVERRIDE;

    void diffuse_velocity_setup (Real       dt,
				 MultiFab*& delta_rhs,
				 MultiFab**& viscn,
				 MultiFab**& viscnp1); 
    //
    // Initialize the pressure by iterating the initial timestep.
    //
    void post_init_press (Real&        dt_init,
			  Array<int>&  nc_save,
			  Array<Real>& dt_save);
    //
    // Compute the mac sync correction.
    //
    virtual void mac_sync () BL_OVERRIDE;

    //
    // Reflux function.
    //
    virtual void reflux () BL_OVERRIDE;

    virtual void avgDown () BL_OVERRIDE; // Average down for all the state types.
    void avgDown (int comp); // Average down for a single StateType scalar
    //
    // get fluxes out of the advection routines.
    //
    void pullFluxes (int        gridno,
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
			       Real      time) BL_OVERRIDE;
    //
    // Calculate divU, which sets them to zero by default.
    //
    virtual void calc_divu (Real      time,
			    Real      dt,
			    MultiFab& fab) BL_OVERRIDE;
    //
    // Calculate dSdT, which sets them to zero by default.
    //
    virtual void calc_dsdt (Real      time,
			    Real      dt,
			    MultiFab& fab) BL_OVERRIDE;
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
				const int  ncycle) BL_OVERRIDE;

    virtual void calcDiffusivity (const Real time) BL_OVERRIDE;

    virtual void getViscosity (MultiFab* viscosity[BL_SPACEDIM],
			       const Real time) BL_OVERRIDE;

    void getDiffusivity (MultiFab*  diffusivity[BL_SPACEDIM],
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
    void manual_tags_placement (TagBoxArray&    tags, 
				Array<IntVect>& bf_lev);

    NavierStokes& getLevel (int lev) { 
	return dynamic_cast<NavierStokes&> ( parent->getLevel(lev) );
    }
};

#endif /*_NavierStokes_H_*/
