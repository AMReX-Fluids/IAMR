
#ifndef _NavierStokes_H_
#define _NavierStokes_H_

#include <NavierStokesBase.H>
#include <AMReX_LevelBld.H>
#include "Utilities.H"

class NSBld
    :
    public amrex::LevelBld
{
    virtual void variableSetUp () override;
    virtual void variableCleanUp () override;
    virtual amrex::AmrLevel *operator() () override;
    virtual amrex::AmrLevel *operator() (amrex::Amr&            papa,
                                         int             lev,
                                         const amrex::Geometry& level_geom,
                                         const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm,
                                         amrex::Real            time) override;
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

    NavierStokes (amrex::Amr&            papa,
                  int             lev,
                  const amrex::Geometry& level_geom,
                  const amrex::BoxArray& bl,
                  const amrex::DistributionMapping& dm,
                  amrex::Real            time);

    virtual ~NavierStokes ();

    ////////////////////////////////////////////////////////////////////////////
    //    AmrLevel virtual functions                                          //
    ////////////////////////////////////////////////////////////////////////////
    //
    // Init grid data at problem start-up.
    //
    virtual void initData () override;
    //
    // Problem initialization functions
    //
#include "prob_init.H"
    //
    // Write plot file stuff to specified directory.
    //
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&  os,
                                amrex::VisMF::How     how) override;
    //
    // Modify list of variables to be plotted
    //
    virtual void setPlotVariables() override;
    //
    // Returns a amrex::MultiFab containing the derived data for this level.
    // The user is responsible for deleting this pointer when done
    // with it.  If ngrow>0 the amrex::MultiFab is built on the appropriately
    // grown BoxArray.
    //
    virtual std::unique_ptr<amrex::MultiFab>
    derive (const std::string& name, amrex::Real time, int ngrow) override;
    //
    // This version of derive() fills the dcomp'th component of mf with the derived quantity.
    //
    virtual void derive (const std::string& name,
                         amrex::Real        time,
                         amrex::MultiFab&   mf,
                         int                dcomp) override;
    //
    // Insure state, and pressure are consistent.
    //
    virtual void post_init (amrex::Real stop_time) override;
    //
    // Advance grids at this level in time.
    //
    virtual amrex::Real advance (amrex::Real time,
                          amrex::Real dt,
                          int  iteration,
                          int  ncycle) override;

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokesBase public virtual functions                           //
    ////////////////////////////////////////////////////////////////////////////

    //
    // Compute the mac sync correction.
    //
    virtual void mac_sync () override;
    //
    // Reflux function.
    //
    virtual void reflux () override;
    //
    // Calculate divU, which sets them to zero by default.
    //
    virtual void calc_divu (amrex::Real      time,
			    amrex::Real      dt,
			    amrex::MultiFab& fab) override;

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokes public static functions                                //
    ////////////////////////////////////////////////////////////////////////////

    //
    // Define data descriptors.
    //
    static void variableSetUp ();

    // FIXME? - size m_bc_values better.
    // Max possible number of scalars is 4: Density, 2 Tracers, Temperature
    static constexpr int NUM_STATE_MAX = AMREX_SPACEDIM+4;
    static amrex::GpuArray<amrex::GpuArray<amrex::Real, NUM_STATE_MAX>, AMREX_SPACEDIM*2>
    //static auto
    get_bc_values () { return m_bc_values; }

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokes public functions                                       //
    ////////////////////////////////////////////////////////////////////////////

    //
    // Advect scalars.
    //
    void scalar_advection (amrex::Real dt,
                           int  first_scalar,
                           int  last_scalar);
    //
    void scalar_diffusion_update (amrex::Real dt,
                                  int  first_scalar,
                                  int  last_scalar);

protected:

    ////////////////////////////////////////////////////////////////////////////
    //    Overriding Virtual Functions in NavierStokesBase                    //
    ////////////////////////////////////////////////////////////////////////////

    virtual void avgDown () override; // Average down for all the state types.
    //
    // Calculate nonuniform viscosity and diffusivity
    //
    virtual void calcViscosity (const amrex::Real time,
				const amrex::Real dt,
				const int  iteration,
				const int  ncycle) override;
    //
    virtual void calcDiffusivity (const amrex::Real time) override;
    //
    virtual void getViscosity (amrex::MultiFab*  viscosity[BL_SPACEDIM],
			       const amrex::Real time) override;
    //
    // Compute viscous terms.
    //
    virtual void getViscTerms (amrex::MultiFab& visc_terms,
			       int       src_comp, 
			       int       num_comp,
			       amrex::Real      time) override;
    //
    virtual void sum_integrated_quantities () override;
    //
    virtual void velocity_diffusion_update (amrex::Real dt) override;

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokes protected static functions                             //
    ////////////////////////////////////////////////////////////////////////////

    virtual void errorEst (amrex::TagBoxArray& tags,
                           int                 clearval,
                           int                 tagval,
                           amrex::Real         time,
                           int                 n_error_buf,
                           int                 ngrow) override;

    static void Initialize ();   // Read input file
    static void Initialize_specific ();   // Read input file specific to IAMR
    static void Finalize ();

private:

    void getDiffusivity (amrex::MultiFab*  diffusivity[BL_SPACEDIM],
			 const amrex::Real time,
			 const int  state_comp,
			 const int  dst_comp,
			 const int  num_comp);

    NavierStokes& getLevel (int lev) { 
	return dynamic_cast<NavierStokes&> ( parent->getLevel(lev) );
    }

    amrex::Real MaxVal (const std::string& name,
                 amrex::Real               time);
    //
    // Initialize the pressure by iterating the initial timestep.
    //
    void post_init_press (amrex::Real&        dt_init,
			  amrex::Vector<int>&  nc_save,
			  amrex::Vector<amrex::Real>& dt_save);
    //
    // Update scalars, (viscous solve in scalar_update).
    //
    void scalar_update (amrex::Real dt,
                        int  first_scalar,
                        int  last_scalar);
    //
    // Define error estimators
    //
    static void error_setup ();

    ///////////////////
    // Private data  //
    ///////////////////

    //
    // Holds info for dynamically generated error estimators
    //
    static amrex::Vector<amrex::AMRErrorTag> errtags;

    //
    // Hold external Dirichlet BC values
    //
    // FIXME? - dymanic size. And fix in NavierStokes.cpp NS_bcfill as well!
    // Max possible number of scalars is 4: Density, 2 Tracers, Temperature
    //static amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM+4>,
    static amrex::GpuArray<amrex::GpuArray<amrex::Real, NUM_STATE_MAX>,
    			   AMREX_SPACEDIM*2> m_bc_values;
};

#endif /*_NavierStokes_H_*/
