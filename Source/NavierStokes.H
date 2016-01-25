
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

    ////////////////////////////////////////////////////////////////////////////
    //    AmrLevel virtual functions                                          //
    ////////////////////////////////////////////////////////////////////////////
    //
    // Init grid data at problem start-up.
    //
    virtual void initData () BL_OVERRIDE;
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


    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokes public static functions                                //
    ////////////////////////////////////////////////////////////////////////////

    //
    // Define data descriptors.
    //
    static void variableSetUp ();
    //
    // Cleanup data descriptors at end of run.
    //
    static void variableCleanUp ();

protected:

    ////////////////////////////////////////////////////////////////////////////
    //    Overriding Virtual Functions in NavierStokesBase                    //
    ////////////////////////////////////////////////////////////////////////////

    virtual void avgDown () BL_OVERRIDE; // Average down for all the state types.
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
    // Calculate nonuniform viscosity and diffusivity
    //
    virtual void calcViscosity (const Real time,
				const Real dt,
				const int  iteration,
				const int  ncycle) BL_OVERRIDE;
    //
    virtual void calcDiffusivity (const Real time) BL_OVERRIDE;
    //
    virtual void getViscosity (MultiFab*  viscosity[BL_SPACEDIM],
			       const Real time) BL_OVERRIDE;
    //
    // Compute viscous terms.
    //
    virtual void getViscTerms (MultiFab& visc_terms,
			       int       src_comp, 
			       int       num_comp,
			       Real      time) BL_OVERRIDE;
    //
    // Compute the mac sync correction.
    //
    virtual void mac_sync () BL_OVERRIDE;
    //
    // Reflux function.
    //
    virtual void reflux () BL_OVERRIDE;
    //
    virtual void sum_integrated_quantities () BL_OVERRIDE;
    //
    virtual void velocity_diffusion_update (Real dt) BL_OVERRIDE;

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokes protected static functions                             //
    ////////////////////////////////////////////////////////////////////////////

    static void Initialize ();   // Read input file
    static void Finalize ();

private:

    void avgDown (int comp); // Average down for a single StateType scalar

    void center_to_edge_plain (const FArrayBox& ccfab,
                               FArrayBox&       ecfab,
                               int              sComp,
                               int              dComp,
                               int              nComp);
    //
    void diffuse_velocity_setup (Real       dt,
				 MultiFab*& delta_rhs,
				 MultiFab**& viscn,
				 MultiFab**& viscnp1); 

    void getDiffusivity (MultiFab*  diffusivity[BL_SPACEDIM],
			 const Real time,
			 const int  state_comp,
			 const int  dst_comp,
			 const int  num_comp);

    NavierStokes& getLevel (int lev) { 
	return dynamic_cast<NavierStokes&> ( parent->getLevel(lev) );
    }

    Real MaxVal (const std::string& name,
                 Real               time);
    //
    // Initialize the pressure by iterating the initial timestep.
    //
    void post_init_press (Real&        dt_init,
			  Array<int>&  nc_save,
			  Array<Real>& dt_save);
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
    //
    void scalar_diffusion_update (Real dt,
                                  int  first_scalar,
                                  int  last_scalar);
    //
    // Define error estimators
    //
    static void error_setup ();

};

#endif /*_NavierStokes_H_*/
