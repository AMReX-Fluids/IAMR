
#ifndef _NavierStokes_H_
#define _NavierStokes_H_

#include <NavierStokesBase.H>
#include <AMReX_LevelBld.H>

class NSBld
    :
    public LevelBld
{
    virtual void variableSetUp () override;
    virtual void variableCleanUp () override;
    virtual AmrLevel *operator() () override;
    virtual AmrLevel *operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& level_geom,
                                  const BoxArray& ba,
                                  Real            time) override;
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
    virtual void initData () override;
    //
    // Write plot file stuff to specified directory.
    //
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&  os,
                                VisMF::How     how) override;
    //
    // Modify list of variables to be plotted
    //
    virtual void setPlotVariables() override;
    //
    // Returns a MultiFab containing the derived data for this level.
    // The user is responsible for deleting this pointer when done
    // with it.  If ngrow>0 the MultiFab is built on the appropriately
    // grown BoxArray.
    //
    virtual MultiFab* derive (const std::string& name,
                              Real               time,
                              int                ngrow) override;
    //
    // This version of derive() fills the dcomp'th component of mf with the derived quantity.
    //
    virtual void derive (const std::string& name,
                         Real               time,
                         MultiFab&          mf,
                         int                dcomp) override;
    //
    // Insure state, and pressure are consistent.
    //
    virtual void post_init (Real stop_time) override;
    //
    // Advance grids at this level in time.
    //
    virtual Real advance (Real time,
                          Real dt,
                          int  iteration,
                          int  ncycle) override;


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

    virtual void avgDown () override; // Average down for all the state types.
    //
    // Calculate divU, which sets them to zero by default.
    //
    virtual void calc_divu (Real      time,
			    Real      dt,
			    MultiFab& fab) override;
    //
    // Calculate dSdT, which sets them to zero by default.
    //
    virtual void calc_dsdt (Real      time,
			    Real      dt,
			    MultiFab& fab) override;
    //
    // Calculate nonuniform viscosity and diffusivity
    //
    virtual void calcViscosity (const Real time,
				const Real dt,
				const int  iteration,
				const int  ncycle) override;
    //
    virtual void calcDiffusivity (const Real time) override;
    //
    virtual void getViscosity (MultiFab*  viscosity[BL_SPACEDIM],
			       const Real time) override;
    //
    // Compute viscous terms.
    //
    virtual void getViscTerms (MultiFab& visc_terms,
			       int       src_comp, 
			       int       num_comp,
			       Real      time) override;
    //
    // Compute the mac sync correction.
    //
    virtual void mac_sync () override;
    //
    // Reflux function.
    //
    virtual void reflux () override;
    //
    virtual void sum_integrated_quantities () override;
    //
    virtual void velocity_diffusion_update (Real dt) override;

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
