
#ifndef IAMR_NavierStokes_H_
#define IAMR_NavierStokes_H_

#include <NavierStokesBase.H>
#include <AMReX_LevelBld.H>
#include "Utilities.H"

class NSBld
    :
    public amrex::LevelBld
{
    void variableSetUp () override;
    void variableCleanUp () override;
    amrex::AmrLevel *operator() () override;
    amrex::AmrLevel *operator() (amrex::Amr&            papa,
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

    ~NavierStokes () override = default;

    NavierStokes (NavierStokes const&) = delete;
    NavierStokes (NavierStokes &&) = delete;
    NavierStokes& operator= (NavierStokes const&) = delete;
    NavierStokes& operator= (NavierStokes &&) = delete;

    ////////////////////////////////////////////////////////////////////////////
    //    AmrLevel virtual functions                                          //
    ////////////////////////////////////////////////////////////////////////////
    //
    // Init grid data at problem start-up.
    //
    void initData () override;
    //
    // Problem initialization functions
    //
#include "prob_init.H"
    //
    // Build/fill any additional data structures after restart.
    //
    void post_restart () override;
    //
    // Do pre-plotfile work
    //
    void writePlotFilePre (const std::string& dir,
                           std::ostream&      os) override;
    //
    // Do post-plotfile work
    //
    void writePlotFilePost (const std::string& dir,
                            std::ostream&      os) override;
    //
    // Modify list of variables to be plotted
    //
    void setPlotVariables() override;
    //
    // Returns a amrex::MultiFab containing the derived data for this level.
    // The user is responsible for deleting this pointer when done
    // with it.  If ngrow>0 the amrex::MultiFab is built on the appropriately
    // grown BoxArray.
    //
    std::unique_ptr<amrex::MultiFab>
    derive (const std::string& name, amrex::Real time, int ngrow) override;
    //
    // This version of derive() fills the dcomp'th component of mf with the derived quantity.
    //
    void derive (const std::string& name,
                 amrex::Real        time,
                 amrex::MultiFab&   mf,
                 int                dcomp) override;
    //
    // Insure state, and pressure are consistent.
    //
    void post_init (amrex::Real stop_time) override;
    //
    // Advance grids at this level in time.
    //
    amrex::Real advance (amrex::Real time,
                         amrex::Real dt,
                         int  iteration,
                         int  ncycle) override;

    // ls related
    amrex::Real advance_semistaggered_twophase_ls (amrex::Real time,
                         amrex::Real dt,
                         int  iteration,
                         int  ncycle);
#ifdef AMREX_PARTICLES
    amrex::Real advance_semistaggered_fsi_diffusedib (amrex::Real time,
                         amrex::Real dt,
                         int  iteration,
                         int  ncycle);
#endif

    amrex::Real advance_semistaggered_twophase_phasefield (amrex::Real time,
                         amrex::Real dt,
                         int  iteration,
                         int  ncycle);

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokesBase public virtual functions                           //
    ////////////////////////////////////////////////////////////////////////////

    //
    // Compute the mac sync correction.
    //
    void mac_sync () override;
    //
    // Reflux function.
    //
    void reflux () override;
    //
    // Calculate divU, which sets them to zero by default.
    //
    void calc_divu (amrex::Real      time,
                    amrex::Real      dt,
                    amrex::MultiFab& divu) override;

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokes public static functions                                //
    ////////////////////////////////////////////////////////////////////////////

    //
    // Define data descriptors.
    //
    static void variableSetUp ();

    //
    // ls related
    //
    // Max possible number of scalars is 5: Density, 2 Tracers, Temperature, ls
    static constexpr int NUM_STATE_MAX = AMREX_SPACEDIM+5;
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
                           int  fscalar,
                           int  lscalar);
    //
    void scalar_diffusion_update (amrex::Real dt,
                                  int  first_scalar,
                                  int  last_scalar);

protected:

    ////////////////////////////////////////////////////////////////////////////
    //    Overriding Virtual Functions in NavierStokesBase                    //
    ////////////////////////////////////////////////////////////////////////////

    void avgDown () override; // Average down for all the state types.
    //
    // Calculate nonuniform viscosity and diffusivity
    //
    void calcViscosity (amrex::Real time,
                        amrex::Real dt,
                        int  iteration,
                        int  ncycle) override;
    //
    void calcDiffusivity (amrex::Real time) override;
    //
    void getViscosity (amrex::MultiFab*  viscosity[AMREX_SPACEDIM],
                       amrex::Real time) override;
    //
    // Compute viscous terms.
    //
    void getViscTerms (amrex::MultiFab& visc_terms,
                       int       src_comp,
                       int       ncomp,
                       amrex::Real      time) override;
    //
    void sum_integrated_quantities () override;
    //
    void velocity_diffusion_update (amrex::Real dt) override;

    ////////////////////////////////////////////////////////////////////////////
    //    NavierStokes protected static functions                             //
    ////////////////////////////////////////////////////////////////////////////

    void errorEst (amrex::TagBoxArray& tags,
                   int                 clearval,
                   int                 tagval,
                   amrex::Real         time,
                   int                 n_error_buf,
                   int                 ngrow) override;

    static void Initialize ();   // Read input file
    static void Initialize_bcs ();
    static void Initialize_diffusivities ();
    static void Finalize ();

private:

    static void getDiffusivity (amrex::MultiFab*  diffusivity[AMREX_SPACEDIM],
                                amrex::Real time,
                                int  state_comp,
                                int  dst_comp,
                                int  ncomp);

    NavierStokes& getLevel (int lev) {
    return dynamic_cast<NavierStokes&> ( parent->getLevel(lev) );
    }

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
    // Runtime parameters
    //
    // For EB, set all covered cells to zero (1, default), or show whatever is
    // in covered cells in plotfile (0).
    static int set_plot_coveredCell_val;

    //
    // Holds info for dynamically generated error estimators
    //
    static amrex::Vector<amrex::AMRErrorTag> errtags;

    //
    // Hold external Dirichlet BC values
    //
    //
    // ls related
    //
    // Max possible number of scalars is 5: Density, 2 Tracers, Temperature, ls
    static amrex::GpuArray<amrex::GpuArray<amrex::Real, NUM_STATE_MAX>,
                   AMREX_SPACEDIM*2> m_bc_values;
};

#endif /*_NavierStokes_H_*/
