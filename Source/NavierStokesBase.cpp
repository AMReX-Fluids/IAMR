
#include <AMReX_ParmParse.H>
#include <AMReX_TagBox.H>
#include <AMReX_Utility.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBInterpolater.H>
#include <AMReX_EBFArrayBox.H>
#include <iamr_mol.H>
#endif

#include <NavierStokesBase.H>
#include <NAVIERSTOKES_F.H>
#include <NAVIERSTOKESBASE_F.H>

#include <PROB_NS_F.H>

//fixme, for writesingle level plotfile
#include<AMReX_PlotFileUtil.H>

using namespace amrex;

Godunov*    NavierStokesBase::godunov       = 0;
ErrorList   NavierStokesBase::err_list;
BCRec       NavierStokesBase::phys_bc;
Projection* NavierStokesBase::projector     = 0;
MacProj*    NavierStokesBase::mac_projector = 0;

Real NavierStokesBase::init_shrink        = 1.0;
int  NavierStokesBase::init_iter          = 2;
int  NavierStokesBase::init_vel_iter      = 1;
Real NavierStokesBase::cfl                = 0.8;
Real NavierStokesBase::change_max         = 1.1;
Real NavierStokesBase::init_dt            = -1.0;
Real NavierStokesBase::fixed_dt           = -1.0;
bool NavierStokesBase::stop_when_steady   = false;
Real NavierStokesBase::steady_tol         = 1.0e-10;
int  NavierStokesBase::initial_iter       = false;
int  NavierStokesBase::initial_step       = false;
Real NavierStokesBase::dt_cutoff          = 0.0;
int  NavierStokesBase::sum_interval       = -1;
int  NavierStokesBase::turb_interval      = -1;
int  NavierStokesBase::jet_interval       = -1;
int  NavierStokesBase::jet_interval_split = 2;

int  NavierStokesBase::radius_grow = 1;
int  NavierStokesBase::verbose     = 0;
Real NavierStokesBase::gravity     = 0.0;
int  NavierStokesBase::NUM_SCALARS = 0;
int  NavierStokesBase::NUM_STATE   = 0;

Vector<AdvectionForm> NavierStokesBase::advectionType;
Vector<DiffusionForm> NavierStokesBase::diffusionType;

Vector<int>  NavierStokesBase::is_diffusive;
Vector<Real> NavierStokesBase::visc_coef;
Real        NavierStokesBase::visc_tol           = 1.0e-10;
Real        NavierStokesBase::visc_abs_tol       = 1.0e-10;
Real        NavierStokesBase::be_cn_theta        = 0.5;

int         NavierStokesBase::Tracer                    = -1;
int         NavierStokesBase::Tracer2                   = -1;
int         NavierStokesBase::Temp                      = -1;
int         NavierStokesBase::do_trac2                  = 0;
int         NavierStokesBase::do_temp                   = 0;
int         NavierStokesBase::do_cons_trac              = 0;
int         NavierStokesBase::do_cons_trac2             = 0;
int         NavierStokesBase::do_sync_proj              = 1;
int         NavierStokesBase::do_reflux                 = 1;
int         NavierStokesBase::modify_reflux_normal_vel  = 0;
int         NavierStokesBase::do_mac_proj               = 1;
int         NavierStokesBase::do_refine_outflow         = 0;
int         NavierStokesBase::do_derefine_outflow       = 1;
int         NavierStokesBase::Nbuf_outflow              = 1;
int         NavierStokesBase::do_denminmax              = 0;
int         NavierStokesBase::do_scalminmax             = 0;
int         NavierStokesBase::do_density_ref            = 0;
int         NavierStokesBase::do_tracer_ref             = 0;
int         NavierStokesBase::do_tracer2_ref            = 0;
int         NavierStokesBase::do_vorticity_ref          = 0;
int         NavierStokesBase::do_temp_ref               = 0;
int         NavierStokesBase::do_scalar_update_in_order = 0;
Vector<int>  NavierStokesBase::scalarUpdateOrder;
int         NavierStokesBase::getForceVerbose           = 0;
int         NavierStokesBase::do_LES                    = 0;
int         NavierStokesBase::getLESVerbose             = 0;
std::string NavierStokesBase::LES_model                 = "Smagorinsky";
Real        NavierStokesBase::smago_Cs_cst              = 0.18;
Real        NavierStokesBase::sigma_Cs_cst              = 1.5;



int  NavierStokesBase::Dpdt_Type = -1;

int  NavierStokesBase::additional_state_types_initialized = 0;
int  NavierStokesBase::Divu_Type                          = -1;
int  NavierStokesBase::Dsdt_Type                          = -1;
int  NavierStokesBase::num_state_type                     = 2;
int  NavierStokesBase::have_divu                          = 0;
int  NavierStokesBase::have_dsdt                          = 0;
Real NavierStokesBase::divu_relax_factor                  = 0.0;
int  NavierStokesBase::S_in_vel_diffusion                 = 0;
int  NavierStokesBase::do_init_vort_proj                  = 0;
int  NavierStokesBase::do_init_proj                       = 1;

int  NavierStokesBase::do_running_statistics  = 0;
Real NavierStokesBase::volWgtSum_sub_origin_x = 0;
Real NavierStokesBase::volWgtSum_sub_origin_y = 0;
Real NavierStokesBase::volWgtSum_sub_origin_z = 0;
Real NavierStokesBase::volWgtSum_sub_Rcyl     = -1;
Real NavierStokesBase::volWgtSum_sub_dx       = -1;
Real NavierStokesBase::volWgtSum_sub_dy       = -1;
Real NavierStokesBase::volWgtSum_sub_dz       = -1;

int  NavierStokesBase::do_mom_diff            = 0;
int  NavierStokesBase::predict_mom_together   = 0;
bool NavierStokesBase::def_harm_avg_cen2edge  = false;

#ifdef AMREX_USE_EB
int          NavierStokesBase::refine_cutcells     = 1;
bool         NavierStokesBase::eb_initialized      = false;
bool         NavierStokesBase::no_eb_in_domain     = true;
bool         NavierStokesBase::body_state_set      = false;
std::vector<Real> NavierStokesBase::body_state;
#endif

namespace
{
    bool initialized = false;
    int  dump_plane  = -1;
    std::string dump_plane_name("SLABS/vel-");
    bool benchmarking = false;
}

#ifdef AMREX_PARTICLES
namespace
{
    //
    // Name of subdirectory in chk???? holding checkpointed particles.
    //
    const std::string the_ns_particle_file_name("Particles");
    //
    // There's really only one of these.
    //
    AmrTracerParticleContainer* NSPC = 0;

    std::string      timestamp_dir                   ("Timestamps");
    std::vector<int> timestamp_indices;
    std::string      particle_init_file;
    std::string      particle_restart_file;
    std::string      particle_output_file;
    bool             restart_from_nonparticle_chkfile = false;
    int              pverbose                         = 2;
}

AmrTracerParticleContainer* NavierStokesBase::theNSPC () { return NSPC; }
#endif

int NavierStokesBase::DoTrac2() {return NavierStokesBase::do_trac2;}
//
BL_FORT_PROC_DECL(BL_NS_DOTRAC2,bl_ns_dotrac2)(int* dotrac2)
{
    *dotrac2 = NavierStokesBase::DoTrac2();
}

NavierStokesBase::NavierStokesBase ()
{
    rho_qtime    = 0;
    rho_tqtime   = 0;
    sync_reg     = 0;
    advflux_reg  = 0;
    viscflux_reg = 0;
    u_mac        = 0;
    aofs         = 0;
    diffusion    = 0;


    if (!additional_state_types_initialized)
        init_additional_state_types();
}

NavierStokesBase::NavierStokesBase (Amr&            papa,
				    int             lev,
				    const Geometry& level_geom,
				    const BoxArray& bl,
                                    const DistributionMapping& dm,
				    Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dm,time)
{
    if(!additional_state_types_initialized) {
        init_additional_state_types();
    }

    const BoxArray& P_grids = state[Press_Type].boxArray();
    //
    // Alloc old_time pressure.
    //
    state[Press_Type].allocOldData();
    //
    // Alloc space for density and temporary pressure variables.
    //
    if (level > 0)
    {
        rho_avg.define(grids,dmap,1,1,MFInfo(),Factory());
        p_avg.define(P_grids,dmap,1,0,MFInfo(),Factory());
    }

    //
    // rho_half is passed into level_project to be used as sigma in the MLMG
    // solve, but MLMG doesn't copy any ghost cells, it fills what it needs itself.
    // does rho_half still need any ghost cells?
    rho_half.define (grids,dmap,1,1,MFInfo(),Factory());
    rho_ptime.define(grids,dmap,1,1,MFInfo(),Factory());
    rho_ctime.define(grids,dmap,1,1,MFInfo(),Factory());
    rho_qtime  = 0;
    rho_tqtime = 0;
    //
    // Build metric coefficients for RZ calculations.
    // Build volume and areas.
    //
    buildMetrics();


#ifdef AMREX_USE_EB

    init_eb(level_geom, bl, dm);

    //fixme? not 100% sure this is the right place
    gradp.reset(new MultiFab(grids,dmap,BL_SPACEDIM,1, MFInfo(), Factory()));
    gradp->setVal(0.);

    //FIXME --- this fn is really similar to restart()... work on that later
#endif

    //
    // Set up reflux registers.
    //
    sync_reg = 0;
    if (level > 0 && do_sync_proj)
    {
        sync_reg = new SyncRegister(grids,dmap,crse_ratio);
    }
    advflux_reg  = 0;
    viscflux_reg = 0;
    if (level > 0 && do_reflux)
    {
        advflux_reg  = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
        viscflux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
    }
    //
    // Initialize work multifabs.
    //
    u_mac   = 0;
    aofs    = 0;

    //
    // Set up the level projector.
    //
    if (projector == 0)
    {
        projector = new Projection(parent,&phys_bc,do_sync_proj,
                                   parent->finestLevel(),radius_grow);
    }
    projector->install_level(level,this,&radius);
    //
    // Set up the godunov box.
    //
    SetGodunov();
    //
    // Set up diffusion.
    //
    diffusion = new Diffusion(parent,this,
                              (level > 0) ? getLevel(level-1).diffusion : 0,
                              NUM_STATE,viscflux_reg,is_diffusive,visc_coef);
    //
    // Allocate the storage for variable viscosity and diffusivity
    //
    viscn   = fb_viscn.define(this,1,0);
    viscnp1 = fb_viscnp1.define(this,1,0);
    diffn   = fb_diffn.define(this,NUM_STATE-Density-1,0);
    diffnp1 = fb_diffnp1.define(this,NUM_STATE-Density-1,0);

    //
    // Set up the mac projector.
    //
    if (mac_projector == 0)
    {
        mac_projector = new MacProj(parent,parent->finestLevel(),
                                    &phys_bc,radius_grow);
    }
    mac_projector->install_level(level,this);
}

NavierStokesBase::~NavierStokesBase ()
{
    delete rho_qtime;
    delete rho_tqtime;
    delete sync_reg;
    delete advflux_reg;
    delete viscflux_reg;
    delete [] u_mac;

    if (mac_projector != 0)
        mac_projector->cleanup(level);
    //
    // Remove the arrays for variable viscosity and diffusivity
    // and delete the Diffusion object
    //
    fb_viscn.clear();
    fb_viscnp1.clear();
    fb_diffn.clear();
    fb_diffnp1.clear();

    delete diffusion;
}

void
NavierStokesBase::allocOldData ()
{
    bool init_pres = !(state[Press_Type].hasOldData());
    AmrLevel::allocOldData();
    if (init_pres)
        initOldPress();
}

void
NavierStokesBase::variableCleanUp ()
{
    desc_lst.clear();
    derive_lst.clear();

    delete godunov;
    godunov = 0;

    err_list.clear();

    delete projector;
    projector = 0;

    delete mac_projector;
    mac_projector = 0;

#ifdef AMREX_PARTICLES
    delete NSPC;
    NSPC = 0;
#endif
}

void
NavierStokesBase::Initialize ()
{
    if (initialized) return;

    ParmParse pp("ns");

    pp.query("dump_plane",dump_plane);

    pp.query("benchmarking",benchmarking);

    pp.query("v",verbose);


    //
    // Get timestepping parameters.
    //
    pp.get("cfl",cfl);
    pp.query("init_iter",init_iter);
    pp.query("init_vel_iter",init_vel_iter);
    pp.query("init_shrink",init_shrink);
    pp.query("dt_cutoff",dt_cutoff);
    pp.query("change_max",change_max);
    pp.query("fixed_dt",fixed_dt);
    pp.query("init_dt", init_dt);
    pp.query("stop_when_steady",stop_when_steady);
    pp.query("steady_tol",steady_tol);
    pp.query("sum_interval",sum_interval);
    pp.query("turb_interval",turb_interval);
    pp.query("jet_interval",jet_interval);
    pp.query("jet_interval_split",jet_interval_split);
    pp.query("gravity",gravity);
    //
    // Get run options.
    //
    pp.query("do_temp",                  do_temp          );
    pp.query("do_trac2",                 do_trac2         );
    pp.query("do_cons_trac",             do_cons_trac     );
    pp.query("do_cons_trac2",            do_cons_trac2    );
    int initial_do_sync_proj =           do_sync_proj;
    pp.query("do_sync_proj",             do_sync_proj     );
    pp.query("do_reflux",                do_reflux        );
    pp.query("modify_reflux_normal_vel", modify_reflux_normal_vel);
    pp.query("do_init_vort_proj",        do_init_vort_proj);
    pp.query("do_init_proj",             do_init_proj     );
    pp.query("do_mac_proj",              do_mac_proj      );
    pp.query("do_denminmax",             do_denminmax     );
    pp.query("do_scalminmax",            do_scalminmax    );
    pp.query("do_density_ref",           do_density_ref   );
    pp.query("do_tracer_ref",            do_tracer_ref    );
    pp.query("do_tracer2_ref",           do_tracer2_ref   );
    pp.query("do_vorticity_ref",         do_vorticity_ref );
    pp.query("do_temp_ref",              do_temp_ref      );

    pp.query("visc_tol",visc_tol);
    pp.query("visc_abs_tol",visc_abs_tol);

    if (modify_reflux_normal_vel)
        amrex::Abort("modify_reflux_normal_vel is no longer supported");

    pp.query("getForceVerbose",          getForceVerbose  );
    pp.query("do_LES",                   do_LES  );
    pp.query("getLESVerbose",            getLESVerbose  );
    pp.query("LES_model",                LES_model  );
    pp.query("smago_Cs_cst",             smago_Cs_cst  );
    pp.query("sigma_Cs_cst",             sigma_Cs_cst  );

#ifdef AMREX_USE_EB
    pp.query("refine_cutcells", refine_cutcells);
#endif

    pp.query("do_scalar_update_in_order",do_scalar_update_in_order );
    if (do_scalar_update_in_order) {
	    const int n_scalar_update_order_vals = pp.countval("scalar_update_order");
	    scalarUpdateOrder.resize(n_scalar_update_order_vals);
	    int got_scalar_update_order = pp.queryarr("scalar_update_order",scalarUpdateOrder,0,n_scalar_update_order_vals);
    }

    // Don't let init_shrink be greater than 1
    if (init_shrink > 1.0)
        amrex::Abort("NavierStokesBase::Initialize(): init_shrink cannot be greater than 1");

    pp.query("divu_relax_factor",divu_relax_factor);
    pp.query("S_in_vel_diffusion",S_in_vel_diffusion);
    if ( S_in_vel_diffusion ){
#ifdef AMREX_USE_EB
      // Currently, we should use the TensorOp to compute the divU terms in divtau.
      // The code is still present to use the source term S instead of a numerically
      // computed divu, however, divmusi terms isn't EB-aware.
      // Perhaps one day a comparision would be interesting.
      amrex::Abort("S_in_vel_diffusion not currently supported.\n");
#else
      amrex::Warning("WARNING: S_in_vel_diffusion is probably not what you want anymore. \nSuggested option is now to set S_in_vel_diffusion=0 to allow the tensor diffusion solver to compute divU.");
#endif
    }
    pp.query("be_cn_theta",be_cn_theta);
    if (be_cn_theta > 1.0 || be_cn_theta < .5)
        amrex::Abort("NavierStokesBase::Initialize(): Must have be_cn_theta <= 1.0 && >= .5");
    //
    // Set parameters dealing with how grids are treated at outflow boundaries.
    //
    pp.query("do_refine_outflow",do_refine_outflow);
    pp.query("do_derefine_outflow",do_derefine_outflow);
    if (do_derefine_outflow == 1 && do_refine_outflow == 1)
      amrex::Abort("NavierStokesBase::Initialize(): Cannot have both do_refine_outflow==1 and do_derefine_outflow==1");

    pp.query("Nbuf_outflow",Nbuf_outflow);
    BL_ASSERT(Nbuf_outflow >= 0);
    BL_ASSERT(!(Nbuf_outflow <= 0 && do_derefine_outflow == 1));

    //
    // Check whether we are doing running statistics.
    //
    pp.query("do_running_statistics",do_running_statistics);

    // If dx,dy,dz,Rcyl<0 (default) the volWgtSum is computed over the entire domain
    pp.query("volWgtSum_sub_origin_x",volWgtSum_sub_origin_x);
    pp.query("volWgtSum_sub_origin_y",volWgtSum_sub_origin_y);
    pp.query("volWgtSum_sub_origin_z",volWgtSum_sub_origin_z);
    pp.query("volWgtSum_sub_Rcyl",volWgtSum_sub_Rcyl);
    pp.query("volWgtSum_sub_dx",volWgtSum_sub_dx);
    pp.query("volWgtSum_sub_dy",volWgtSum_sub_dy);
    pp.query("volWgtSum_sub_dz",volWgtSum_sub_dz);


    // Are we going to do velocity or momentum update?

    pp.query("do_mom_diff",do_mom_diff);
    pp.query("predict_mom_together",predict_mom_together);

    if (do_mom_diff == 0 && predict_mom_together == 1)
    {
      amrex::Print() << "MAKES NO SENSE TO HAVE DO_MOM_DIFF=0 AND PREDICT_MOM_TOGETHER=1\n";
      exit(0);
    }

    pp.query("harm_avg_cen2edge", def_harm_avg_cen2edge);

#ifdef AMREX_PARTICLES
    read_particle_params ();
#endif

    amrex::ExecOnFinalize(NavierStokesBase::Finalize);

    initialized = true;
}

// The following Initialize_specific is dedicated to read and set data
// only specific for IAMR, because it conflicts with PeleLM.
// PeleLM calls NavierStokesBase::Initialize() and its own PelelM::Initialize_specific ()
void
NavierStokesBase::Initialize_specific ()
{
    ParmParse pp("ns");

    Vector<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
    pp.getarr("lo_bc",lo_bc,0,BL_SPACEDIM);
    pp.getarr("hi_bc",hi_bc,0,BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        phys_bc.setLo(i,lo_bc[i]);
        phys_bc.setHi(i,hi_bc[i]);
    }

    read_geometry();
    //
    // Check phys_bc against possible periodic geometry
    // if periodic, must have internal BC marked.
    //
    if (DefaultGeometry().isAnyPeriodic())
    {
        //
        // Do idiot check.  Periodic means interior in those directions.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (DefaultGeometry().isPeriodic(dir))
            {
                if (lo_bc[dir] != Interior)
                {
                    std::cerr << "NavierStokesBase::variableSetUp:periodic in direction "
                              << dir
                              << " but low BC is not Interior\n";
                    amrex::Abort("NavierStokesBase::Initialize()");
                }
                if (hi_bc[dir] != Interior)
                {
                    std::cerr << "NavierStokesBase::variableSetUp:periodic in direction "
                              << dir
                              << " but high BC is not Interior\n";
                    amrex::Abort("NavierStokesBase::Initialize()");
                }
            }
        }
    }

    {
        //
        // Do idiot check.  If not periodic, should be no interior.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (!DefaultGeometry().isPeriodic(dir))
            {
              if (lo_bc[dir] == Interior)
              {
                  std::cerr << "NavierStokesBase::variableSetUp:Interior bc in direction "
                            << dir
                            << " but not defined as periodic\n";
                  amrex::Abort("NavierStokesBase::Initialize()");
              }
              if (hi_bc[dir] == Interior)
              {
                  std::cerr << "NavierStokesBase::variableSetUp:Interior bc in direction "
                            << dir
                            << " but not defined as periodic\n";
                  amrex::Abort("NavierStokesBase::Initialize()");
              }
            }
        }
    }

    //
    // Read viscous/diffusive parameters and array of viscous/diffusive coeffs.
    // NOTE: at this point, we dont know number of state variables
    //       so just read all values listed.
    //

    const int n_vel_visc_coef   = pp.countval("vel_visc_coef");
    const int n_temp_cond_coef  = pp.countval("temp_cond_coef");
    const int n_scal_diff_coefs = pp.countval("scal_diff_coefs");

    if (n_vel_visc_coef != 1)
        amrex::Abort("NavierStokesBase::Initialize(): Only one vel_visc_coef allowed");

    if (do_temp && n_temp_cond_coef != 1)
        amrex::Abort("NavierStokesBase::Initialize(): Only one temp_cond_coef allowed");

    int n_visc = BL_SPACEDIM + 1 + n_scal_diff_coefs;
    if (do_temp)
        n_visc++;
    visc_coef.resize(n_visc);
    is_diffusive.resize(n_visc);

    pp.get("vel_visc_coef",visc_coef[0]);
    for (int i = 1; i < BL_SPACEDIM; i++)
      visc_coef[i] = visc_coef[0];
    //
    // Here we set the coefficient for density, which does not diffuse.
    //
    visc_coef[Density] = -1;
    //
    // Set the coefficients for the scalars, but temperature.
    //
    Vector<Real> scal_diff_coefs(n_scal_diff_coefs);
    pp.getarr("scal_diff_coefs",scal_diff_coefs,0,n_scal_diff_coefs);

    int scalId = Density;

    // Will need to add more lines when more variables are added
    Tracer = Density+1;
    if (do_trac2)
	    Tracer2 = Density+2;

    for (int i = 0; i < n_scal_diff_coefs; i++)
    {
        visc_coef[++scalId] = scal_diff_coefs[i];
    }
    //
    // Set the coefficient for temperature.
    //
    if (do_temp)
    {
	    Temp = ++scalId;
	    pp.get("temp_cond_coef",visc_coef[Temp]);
    }
}


void
NavierStokesBase::Finalize ()
{
    initialized = false;
}

void
NavierStokesBase::read_geometry ()
{
#if (BL_SPACEDIM == 2)
    //
    // Must load coord here because Geometry hasn't read it in yet.
    //
    ParmParse pp("geometry");

    int coord;
    pp.get("coord_sys",coord);

    if ((Geometry::CoordType) coord == Geometry::RZ && phys_bc.lo(0) != Symmetry)
    {
        phys_bc.setLo(0,Symmetry);
	amrex::Print() << "\nWarning: Setting phys_bc at xlo to Symmetry\n\n";
    }
#endif
}

void
NavierStokesBase::advance_setup (Real time,
                                 Real dt,
	                         int  iteration,
                                 int  ncycle)
{
    BL_PROFILE("NavierStokesBase::advance_setup()");

    const int finest_level = parent->finestLevel();

#ifdef AMREX_USE_EB
    // incflo now uses: use_godunov ? 4 : 3;
    umac_n_grow = 4;
#else
    umac_n_grow = 1;
#endif

#ifdef AMREX_PARTICLES
    if (ncycle > umac_n_grow)
      umac_n_grow = ncycle;
#endif

    mac_projector->setup(level);
    //
    // Why are they defined here versus the constructor?
    //
    if (level < finest_level)
    {
        if (Vsync.empty())
            Vsync.define(grids,dmap,BL_SPACEDIM,1,MFInfo(),Factory());
        if (Ssync.empty())
	  Ssync.define(grids,dmap,NUM_STATE-BL_SPACEDIM,1,MFInfo(),Factory());
        Vsync.setVal(0);
        Ssync.setVal(0);
    }
    //
    // Set reflux registers to zero.
    //
    if (do_reflux && level < finest_level)
    {
        getAdvFluxReg(level+1).setVal(0);
        getViscFluxReg(level+1).setVal(0);
    }
    //
    // Alloc space for edge velocities (normal comp only).
    //
    if (u_mac == 0 || u_mac[0].nGrow() < umac_n_grow)
    {
	if (u_mac != 0) delete [] u_mac;

        u_mac = new MultiFab[BL_SPACEDIM];

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
	    const BoxArray& edgeba = getEdgeBoxArray(dir);
            u_mac[dir].define(edgeba,dmap,1,umac_n_grow,MFInfo(),Factory());
            u_mac[dir].setVal(1.e40);
        }
    }
    //
    // Alloc MultiFab to hold advective update terms.
    //
    BL_ASSERT(aofs == 0);
    // NOTE: nghost=0 for aofs appears to work. mfix also uses no ghost cells.
    aofs = new MultiFab(grids,dmap,NUM_STATE,0,MFInfo(),Factory());

#ifdef AMREX_USE_EB
    //Slopes in x-direction
    m_xslopes.define(grids, dmap, AMREX_SPACEDIM, Godunov::hypgrow(), MFInfo(), Factory());
    m_xslopes.setVal(0.);
    // Slopes in y-direction
    m_yslopes.define(grids, dmap, AMREX_SPACEDIM, Godunov::hypgrow(), MFInfo(), Factory());
    m_yslopes.setVal(0.);
#if (AMREX_SPACEDIM > 2)
    // Slopes in z-direction
    m_zslopes.define(grids, dmap, AMREX_SPACEDIM, Godunov::hypgrow(), MFInfo(), Factory());
    m_zslopes.setVal(0.);
#endif
#endif

    //
    // Set rho_avg.
    //
    if (!initial_step && level > 0 && iteration == 1)
        initRhoAvg(0.5/Real(ncycle));
    //
    // Set up state multifabs for the advance.
    //
    for (int k = 0; k < num_state_type; k++)
    {
	bool has_old_data = state[k].hasOldData();
        state[k].allocOldData();
	if (! has_old_data) state[k].oldData().setVal(0.0);
        state[k].swapTimeLevels(dt);
    }

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
    {
        const Real new_press_time = .5 * (state[State_Type].prevTime() +
                                          state[State_Type].curTime());
        state[Press_Type].setNewTimeLevel(new_press_time);
    }

    make_rho_prev_time();

    // refRatio==4 is not currently supported
    //
    // If refRatio==4 to the next level coarser, and we're going to diffuse
    // scalars as SoverRho, we're going to need rho at 1/4 and 3/4 time there.
    // Make these things if need be.
    //
    // if (level > 0)
    // {
    //     bool needs_rho4 = false;

    //     if (parent->nCycle(level) == 4)
    //         for (int i = 0; i < NUM_STATE && !needs_rho4; ++i)
    //             needs_rho4 = (diffusionType[i] == Laplacian_SoverRho);

    //     if (needs_rho4)
    //     {
    //         NavierStokesBase& clevel = getLevel(level-1);
    //         const BoxArray& cgrids = clevel.boxArray();
    //         const DistributionMapping& cdmap = clevel.DistributionMap();
    //         const Real      ptime  = clevel.state[State_Type].prevTime();
    //         const Real      ctime  = clevel.state[State_Type].curTime();

    //         if (clevel.rho_qtime == 0)
    //         {
    //             const Real qtime = ptime + 0.25*(ctime-ptime);
    //             clevel.rho_qtime = new MultiFab(cgrids,cdmap,1,1);
    //             FillPatch(clevel,*(clevel.rho_qtime),1,qtime,State_Type,Density,1,0);
    //         }
    //         if (clevel.rho_tqtime == 0)
    //         {
    //             const Real tqtime = ptime + 0.75*(ctime-ptime);
    //             clevel.rho_tqtime = new MultiFab(cgrids,cdmap,1,1);
    // 		FillPatch(clevel, *(clevel.rho_tqtime), 1, tqtime, State_Type, Density, 1, 0);
    //         }
    //    }
    // }
}

//
// Clean up after the advance function.
//
void
NavierStokesBase::advance_cleanup (int iteration, int ncycle)
{
    delete aofs;
    aofs = 0;
}

void
NavierStokesBase::buildMetrics ()
{
    //
    // We "should" only need radius when we're RZ, but some 2-D code is written to
    // access it first and then "use" if if RZ.  It's easier to just always build
    // it for 2D than try to fix the underlying Fortran calls that take radius.
    //
#if (BL_SPACEDIM == 2)
    radius.resize(grids.size());

    const Real dxr = geom.CellSize()[0];

    for (int i = 0; i < grids.size(); i++)
    {
        const int ilo = grids[i].smallEnd(0)-radius_grow;
        const int ihi = grids[i].bigEnd(0)+radius_grow;
        const int len = ihi - ilo + 1;

        radius[i].resize(len);

        RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());

        const Real xlo = gridloc.lo(0) + (0.5 - radius_grow)*dxr;
        for (int j = 0; j < len; j++)
            radius[i][j] = xlo + j*dxr;
    }
#endif

    // fixme? for now, volume and area are intentionally without EB knowledge
    volume.clear();
    volume.define(grids,dmap,1,GEOM_GROW);
    geom.GetVolume(volume);

    for (int dir = 0; dir < BL_SPACEDIM; ++dir)
    {
        area[dir].clear();
	area[dir].define(getEdgeBoxArray(dir),dmap,1,GEOM_GROW);
        geom.GetFaceArea(area[dir],dir);
    }

#ifdef AMREX_USE_EB
    // make sure dx == dy == dz
    const Real* dx = geom.CellSize();
    Print()<<"dx = "<<dx[0]<<" "<<dx[1]<<" "<<dx[2]<<" \n";
    for (int i = 1; i < BL_SPACEDIM; i++){
      if (std::abs(dx[i]-dx[i-1]) > 1.e-12*dx[0])
        amrex::Abort("EB requires dx == dy (== dz)\n");
    }

    const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());
    volfrac = &(ebfactory.getVolFrac());
    areafrac = ebfactory.getAreaFrac();


    //fixme? assume will need this part cribbed from CNS
    // level_mask.clear();
    // level_mask.define(grids,dmap,1,1);
    // level_mask.BuildMask(geom.Domain(), geom.periodicity(),
    //                      level_mask_covered,
    //                      level_mask_notcovered,
    //                      level_mask_physbnd,
    //                      level_mask_interior);

#endif
}

void
NavierStokesBase::calcDpdt ()
{
    BL_ASSERT(state[Press_Type].descriptor()->timeType() == StateDescriptor::Point);

    MultiFab&  new_press   = get_new_data(Press_Type);
    MultiFab&  old_press   = get_old_data(Press_Type);
    MultiFab&  dpdt        = get_new_data(Dpdt_Type);
    const Real dt_for_dpdt = state[Press_Type].curTime()-state[Press_Type].prevTime();

    if (dt_for_dpdt != 0.0)
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(dpdt,true); mfi.isValid(); ++mfi)
        {
            const Box& vbx     = mfi.tilebox();
            FArrayBox& dpdtfab = dpdt[mfi];
            // FIXME MSD: Combine these
            dpdtfab.copy<RunOn::Host>(new_press[mfi],vbx,0,vbx,0,1);
            dpdtfab.minus<RunOn::Host>(old_press[mfi],vbx,0,0,1);
            dpdtfab.divide<RunOn::Host>(dt_for_dpdt,vbx,0,1);
        }
    }
    else
    {
        dpdt.setVal(0);
    }
}

void
NavierStokesBase::checkPoint (const std::string& dir,
			      std::ostream&      os,
			      VisMF::How         how,
			      bool               dump_old)
{
    AmrLevel::checkPoint(dir, os, how, dump_old);

#ifdef AMREX_PARTICLES
    if (level == 0)
    {
        if (NSPC != 0)
            NSPC->Checkpoint(dir,the_ns_particle_file_name);
    }
#endif

# ifdef AMREX_USE_EB
// Need to add gradp in the checkpoint
   std::string LevelDir, FullPath;
   LevelDirectoryNames(dir, LevelDir, FullPath);
   std::string gradp_mf_fullpath = FullPath + "/gradp";
   VisMF::Write(*gradp,gradp_mf_fullpath,how);
#endif
}

void
NavierStokesBase::computeInitialDt (int                   finest_level,
				    int                   sub_cycle,
				    Vector<int>&           n_cycle,
				    const Vector<IntVect>& ref_ratio,
				    Vector<Real>&          dt_level,
				    Real                  stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    int i;

    Real dt_0    = 1.0e+100;
    int n_factor = 1;
    ///TODO/DEBUG: This will need to change for optimal subcycling.
    for (i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor   *= n_cycle[i];
        dt_0        = std::min(dt_0,n_factor*dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    if (stop_time >= 0.0)
    {
        const Real eps      = 0.0001*dt_0;
        const Real cur_time = state[State_Type].curTime();
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        n_factor   *= n_cycle[i];
        dt_level[i] = dt_0/( (Real)n_factor );
    }
}

void
NavierStokesBase::computeNewDt (int                   finest_level,
				int                   sub_cycle,
				Vector<int>&           n_cycle,
				const Vector<IntVect>& ref_ratio,
				Vector<Real>&          dt_min,
				Vector<Real>&          dt_level,
				Real                  stop_time,
				int                   post_regrid_flag)
{
    //
    // We are at the start of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;

    int i;

    Real dt_0     = 1.0e+100;
    int  n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        NavierStokesBase& adv_level = getLevel(i);
        dt_min[i] = std::min(dt_min[i],adv_level.estTimeStep());
    }

    if (fixed_dt <= 0.0)
    {
       if (post_regrid_flag == 1)
       {
          //
          // Limit dt's by pre-regrid dt
          //
          for (i = 0; i <= finest_level; i++)
          {
              dt_min[i] = std::min(dt_min[i],dt_level[i]);
          }
       }
       else
       {
          //
          // Limit dt's by change_max * old dt
          //
          for (i = 0; i <= finest_level; i++)
          {
	    if (verbose)
                 if (dt_min[i] > change_max*dt_level[i])
                 {
                     amrex::Print() << "NavierStokesBase::compute_new_dt : limiting dt at level "
				    << i << '\n';
                     amrex::Print() << " ... new dt computed: " << dt_min[i]
				    << '\n';
                     amrex::Print() << " ... but limiting to: "
				    << change_max * dt_level[i] << " = " << change_max
				    << " * " << dt_level[i] << '\n';
                 }
             dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
          }
       }
    }

    //
    // Find the minimum over all levels
    //
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0      = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps      = 0.0001*dt_0;
    const Real cur_time = state[State_Type].curTime();
    if (stop_time >= 0.0)
    {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    //
    // Set dt at each level of refinement
    //
    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        n_factor   *= n_cycle[i];
        dt_level[i] = dt_0/( (Real)n_factor );
    }
}

void
NavierStokesBase::create_mac_rhs (MultiFab& rhs, int nGrow, Real time, Real dt)
{
    BL_PROFILE("NavierStokesBase::create_mac_rhs()");

    BL_ASSERT(rhs.nGrow()>=nGrow);
    BL_ASSERT(rhs.boxArray()==grids);

    const int sCompDivU = 0;
    const int nCompDivU = 1;
    const int sCompDsdt = 0;
    const int nCompDsdt = 1;

    if (have_divu)
    {
	FillPatch(*this,rhs,nGrow,time,Divu_Type,sCompDivU,nCompDivU,sCompDivU);
    }
    else
    {
        rhs.setVal(0.);
    }

    if (have_dsdt)
    {
	FillPatchIterator fpi(*this,rhs,nGrow,time,Dsdt_Type,sCompDsdt,nCompDsdt);
	const MultiFab& mf = fpi.get_mf();
	MultiFab::Saxpy(rhs, 0.5*dt, mf, 0, sCompDsdt, nCompDsdt, nGrow);
    }
}

void
NavierStokesBase::create_umac_grown (int nGrow)
{
    BL_PROFILE("NavierStokesBase::create_umac_grown()");

    if (level > 0)
    {
        BoxList bl = amrex::GetBndryCells(grids,nGrow);

        BoxArray f_bnd_ba(std::move(bl));

        BoxArray c_bnd_ba = f_bnd_ba; c_bnd_ba.coarsen(crse_ratio);

        c_bnd_ba.maxSize(32);

        f_bnd_ba = c_bnd_ba; f_bnd_ba.refine(crse_ratio);

        for (int n = 0; n < BL_SPACEDIM; ++n)
        {
            //
            // crse_src & fine_src must have same parallel distribution.
            // We'll use the KnapSack distribution for the fine_src_ba.
            // Since fine_src_ba should contain more points, this'll lead
            // to a better distribution.
            //
            BoxArray crse_src_ba(c_bnd_ba), fine_src_ba(f_bnd_ba);

            crse_src_ba.surroundingNodes(n);
            fine_src_ba.surroundingNodes(n);

            const int N = fine_src_ba.size();

            std::vector<long> wgts(N);

#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < N; i++)
                wgts[i] = fine_src_ba[i].numPts();

            DistributionMapping dm;
            // This DM won't be put into the cache.
            dm.KnapSackProcessorMap(wgts,ParallelDescriptor::NProcs());

	    // FIXME
	    // Declaring in this way doesn't work. I think it's because the box arrays
	    // have been changed and each src box is not completely contained within a
	    // single box in the Factory's BA
	    // For now, coarse-fine boundary doesn't intersect EB, so should be okay...
            // MultiFab crse_src(crse_src_ba, dm, 1, 0, MFInfo(), getLevel(level-1).Factory());
            // MultiFab fine_src(fine_src_ba, dm, 1, 0, MFInfo(), Factory());
            MultiFab crse_src(crse_src_ba, dm, 1, 0);
            MultiFab fine_src(fine_src_ba, dm, 1, 0);

            crse_src.setVal(1.e200);
            fine_src.setVal(1.e200);
            //
            // We want to fill crse_src from lower level u_mac including u_mac's grow cells.
            //
	    const MultiFab& u_macLL = getLevel(level-1).u_mac[n];
	    crse_src.copy(u_macLL,0,0,1,1,0);

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(crse_src); mfi.isValid(); ++mfi)
            {
                const int  nComp = 1;
                const Box& box   = crse_src[mfi].box();
                const int* rat   = crse_ratio.getVect();
                pc_edge_interp(box.loVect(), box.hiVect(), &nComp, rat, &n,
                                       crse_src[mfi].dataPtr(),
                                       ARLIM(crse_src[mfi].loVect()),
                                       ARLIM(crse_src[mfi].hiVect()),
                                       fine_src[mfi].dataPtr(),
                                       ARLIM(fine_src[mfi].loVect()),
                                       ARLIM(fine_src[mfi].hiVect()));
            }
            crse_src.clear();
            //
            // Replace pc-interpd fine data with preferred u_mac data at
            // this level u_mac valid only on surrounding faces of valid
            // region - this op will not fill grow region.
            //
            fine_src.copy(u_mac[n]);
            //
            // Interpolate unfilled grow cells using best data from
            // surrounding faces of valid region, and pc-interpd data
            // on fine edges overlaying coarse edges.
            //
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(fine_src); mfi.isValid(); ++mfi)
            {
                const int  nComp = 1;
                const Box& fbox  = fine_src[mfi].box();
                const int* rat   = crse_ratio.getVect();
                edge_interp(fbox.loVect(), fbox.hiVect(), &nComp, rat, &n,
                                 fine_src[mfi].dataPtr(),
                                 ARLIM(fine_src[mfi].loVect()),
                                 ARLIM(fine_src[mfi].hiVect()));
            }

	    MultiFab u_mac_save(u_mac[n].boxArray(),u_mac[n].DistributionMap(),1,0,MFInfo(),Factory());
	    u_mac_save.copy(u_mac[n]);
	    u_mac[n].copy(fine_src,0,0,1,0,nGrow);
	    u_mac[n].copy(u_mac_save);
        }
    }
    //
    // Now we set the boundary data
    // FillBoundary fills grow cells that overlap valid regions.
    // HOEXTRAPTOCC fills outside of domain cells.
    //
    const Real* xlo = geom.ProbLo(); //these aren't actually used by the FORT method
    const Real* dx  = geom.CellSize();

    Box domain_box = geom.Domain();
    for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
        if (geom.isPeriodic(idim)) {
            domain_box.grow(idim,nGrow);
        }
    }

    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
        Box dm = domain_box;
        dm.surroundingNodes(n);
        const int*  lo  = dm.loVect();
        const int*  hi  = dm.hiVect();

        //
        // HOEXTRAPTOCC isn't threaded.  OMP over calls to it.
        //

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(u_mac[n]); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = u_mac[n][mfi];
            amrex_hoextraptocc(BL_TO_FORTRAN_ANYD(fab),lo,hi,dx,xlo);
        }
        // call FillBoundary to make sure that fine/fine grow cells are valid
	u_mac[n].FillBoundary(geom.periodicity());
    }
}

void
NavierStokesBase::diffuse_scalar_setup (int sigma, int& rho_flag)
{

    rho_flag = Diffusion::set_rho_flag(diffusionType[sigma]);
}

void
NavierStokesBase::errorEst (TagBoxArray& tags,
			    int          clearval,
			    int          tagval,
			    Real         time,
			    int          n_error_buf,
			    int          ngrow)
{
    const int*  domain_lo = geom.Domain().loVect();
    const int*  domain_hi = geom.Domain().hiVect();
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    for (int j = 0; j < err_list.size(); j++)
    {
        auto mf = derive(err_list[j].name(), time, err_list[j].nGrow());
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
        {
	          const Box&  vbx     = mfi.tilebox();
            RealBox     gridloc = RealBox(vbx,geom.CellSize(),geom.ProbLo());
            Vector<int>  itags   = tags[mfi].tags();
            int*        tptr    = itags.dataPtr();
            const int*  tlo     = tags[mfi].box().loVect();
            const int*  thi     = tags[mfi].box().hiVect();
            const int*  lo      = vbx.loVect();
            const int*  hi      = vbx.hiVect();
            const Real* xlo     = gridloc.lo();
            FArrayBox&  fab     = (*mf)[mfi];
            Real*       dat     = fab.dataPtr();
            const int*  dlo     = fab.box().loVect();
            const int*  dhi     = fab.box().hiVect();
            const int   ncomp   = fab.nComp();

            err_list[j].errFunc()(tptr, ARLIM(tlo), ARLIM(thi), &tagval,
                                  &clearval, dat, ARLIM(dlo), ARLIM(dhi),
                                  lo,hi, &ncomp, domain_lo, domain_hi,
                                  dx, xlo, prob_lo, &time, &level);
            //
            // Don't forget to set the tags in the TagBox.
            //
            tags[mfi].tags(itags);
        }
    }

#ifdef AMREX_USE_EB
    // Enforce that the EB not cross the coarse-fine boundary
    const auto& ebfactory = dynamic_cast<amrex::EBFArrayBoxFactory const&>(Factory());
    if ( !ebfactory.isAllRegular() )
    {
      //
      // FIXME - For now, always refine cut cells
      //   Later, figure out a slick way to check if EB and CFB cross
      //   and allow !refine_cutcells
      //
      if (!refine_cutcells) amrex::Abort("For now, cutcells must always exist at finest level.");

      // Refine on cut cells
      if (refine_cutcells) // or if EB and CBF cross
      {
        const MultiFab& S_new = get_new_data(State_Type);
        amrex::TagCutCells(tags, S_new);
      }
    }
#endif
}


//
// Estimate the maximum allowable timestep at a cell center.
//
Real
NavierStokesBase::estTimeStep ()
{
    BL_PROFILE("NavierStokesBase::estTimeStep()");

    if (fixed_dt > 0.0)
    {
        Real factor = 1.0;

        if (!(level == 0))
        {
            int ratio = 1;
            for (int lev = 1; lev <= level; lev++)
            {
                ratio *= parent->nCycle(lev);
            }
            factor = 1.0/double(ratio);
        }

        return factor*fixed_dt;
    }

    const Real  small         = 1.0e-8;
    const int   n_grow        = 0;
    Real        estdt         = 1.0e+20;

    const Real  cur_pres_time = state[Press_Type].curTime();
    MultiFab&   U_new         = get_new_data(State_Type);

    Vector<Real> u_max(AMREX_SPACEDIM);
    Vector<Real> f_max(AMREX_SPACEDIM);

#ifdef AMREX_USE_EB
    MultiFab& Gp = getGradP();
    Gp.FillBoundary(geom.periodicity());
#else
    MultiFab Gp(grids,dmap,BL_SPACEDIM,1);
    getGradP(Gp, cur_pres_time);
#endif

    //
    // Find local max of velocity
    //
    u_max = U_new.norm0({AMREX_D_DECL(0,1,2)},0,true,true);

    //
    // Compute forcing terms: in this case this means external forces and grad(p)
    // Viscous terms not included since Crack-Nicholson is unconditionally stable
    // so no need to account for explicit part of viscous term
    //
    MultiFab tforces(grids,dmap,AMREX_SPACEDIM,n_grow,MFInfo(),Factory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox tforces_fab;
        for (MFIter mfi(rho_ctime,true); mfi.isValid(); ++mfi)
        {
            const auto& bx = mfi.tilebox();
            const auto  cur_time = state[State_Type].curTime();

            if (getForceVerbose)
                amrex::Print() << "---" << '\n'
                               << "H - est Time Step:" << '\n'
                               << "Calling getForce..." << '\n';

            getForce(tforces_fab,bx,n_grow,Xvel,AMREX_SPACEDIM,cur_time,U_new[mfi],U_new[mfi],Density);
            tforces[mfi].copy<RunOn::Host>(tforces_fab,bx,0,bx,0,AMREX_SPACEDIM);
        }
    }

    MultiFab::Subtract(tforces, Gp, 0, 0, AMREX_SPACEDIM, 0);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim )
        MultiFab::Divide(tforces, rho_ctime, 0, idim, 1, 0);

    //
    // Find local max of tforces
    //
    f_max = tforces.norm0({AMREX_D_DECL(0,1,2)},0,true,true);

    //
    // Compute local estdt
    //
    const Real* dx = geom.CellSize();

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (u_max[idim] > small)
        {
            estdt = std::min(estdt, dx[idim]/u_max[idim]);
        }

        if (f_max[idim] > small)
        {
            estdt = std::min(estdt, std::sqrt(2.0*dx[idim]/f_max[idim]));
        }
    }

    //
    // Reduce estimated dt by CFL factor and find global min
    //
    ParallelDescriptor::ReduceRealMin(estdt);

    if ( estdt < 1.0e+20) {
      //
      // timestep estimation successful
      //
      estdt = estdt * cfl;
    }
    else if (init_dt > 0 ) {
      //
      // use init_dt, scale for amr level
      //
      Real factor = 1.0;

      if (!(level == 0))
      {
	int ratio = 1;
	for (int lev = 1; lev <= level; lev++)
	{
	  ratio *= parent->nCycle(lev);
	}
	factor = 1.0/double(ratio);
      }

      estdt = factor*init_dt;
    }
    else {
      Print()<<"\nNavierStokesBase::estTimeStep() failed to provide a good timestep "
	     <<"(probably because initial velocity field is zero with no external forcing).\n"
	     <<"Use ns.init_dt to provide a reasonable timestep on coarsest level.\n"
	     <<"Note that ns.init_shrink will be applied to init_dt."<<std::endl;
      amrex::Abort("\n");
    }

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        ParallelDescriptor::ReduceRealMax(u_max.dataPtr(), AMREX_SPACEDIM, IOProc);

	amrex::Print() << "estTimeStep :: \n" << "LEV = " << level << " UMAX = ";
	for (int k = 0; k < AMREX_SPACEDIM; k++)
        {
            amrex::Print() << u_max[k] << "  ";
        }
	amrex::Print() << '\n';

	if (getForceVerbose){
	  ParallelDescriptor::ReduceRealMax(f_max.dataPtr(), AMREX_SPACEDIM, IOProc);
	  amrex::Print() << "        FMAX = ";
	  for (int k = 0; k < AMREX_SPACEDIM; k++)
	  {
            amrex::Print() << f_max[k] << "  ";
	  }
	  amrex::Print() << '\n';
	}
	Print()<<"estimated timestep: dt = "<<estdt<<std::endl;
    }

    return estdt;
}

const MultiFab&
NavierStokesBase::get_rho (Real time)
{
    const TimeLevel whichTime = which_time(State_Type,time);

    if (whichTime == AmrOldTime)
    {
        return rho_ptime;
    }
    else if (whichTime == AmrNewTime)
    {
        return rho_ctime;
    }
    else if (whichTime == Amr1QtrTime)
    {
        BL_ASSERT(rho_qtime);
        return *rho_qtime;
    }
    else if (whichTime == Amr3QtrTime)
    {
        BL_ASSERT(rho_tqtime);
        return *rho_tqtime;
    }
    else if (whichTime == AmrHalfTime)
    {
        return get_rho_half_time();
    }
    else
    {
        amrex::Error("NavierStokesBase::get_rho(): bad time");

        return rho_ptime; // Got to return something to shut up compiler.
    }
}

MultiFab&
NavierStokesBase::get_rho_half_time ()
{
    //
    // Fill it in when needed ...
    //
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(rho_half,true); mfi.isValid(); ++mfi)
    {
        FArrayBox& rhofab = rho_half[mfi];
	const Box& bx = mfi.growntilebox();

        // FIXME MSD: Combine these
        rhofab.copy<RunOn::Host>(rho_ptime[mfi],bx,0,bx,0,1);
        rhofab.plus<RunOn::Host>(rho_ctime[mfi],bx,0,0);
        rhofab.mult<RunOn::Host>(.5,bx);

    }

    return rho_half;
}

//
// Fill patch divU.
//
MultiFab*
NavierStokesBase::getDivCond (int ngrow, Real time)
{
    MultiFab* divu = 0;

    if (!have_divu)
    {
        divu = new MultiFab(grids,dmap,1,ngrow,MFInfo(),Factory());

        divu->setVal(0);
    }
    else
    {
        divu = getState(ngrow,Divu_Type,0,1,time);
    }

    return divu;
}

//
// Fill patch dSdt.
//
MultiFab*
NavierStokesBase::getDsdt (int ngrow, Real time)
{
    MultiFab* dsdt = 0;

    if (!(have_dsdt && have_divu))
    {
        dsdt = new MultiFab(grids,dmap,1,ngrow,MFInfo(),Factory());

        dsdt->setVal(0);
    }
    else
    {
        dsdt = getState(ngrow,Dsdt_Type,0,1,time);
    }

    return dsdt;
}


void
NavierStokesBase::getGradP (MultiFab& gp, Real      time)
{
    BL_PROFILE("NavierStokesBase::getGradP()");

    const int   NGrow = gp.nGrow();
    MultiFab&   P_old = get_old_data(Press_Type);
    const Real* dx    = geom.CellSize();

    if (level > 0 && state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
    {
        //
        // We want to be sure the intersection of old and new grids is
        // entirely contained within gp.boxArray()
        //
        BL_ASSERT(gp.boxArray() == grids);

        {
            const BoxArray& pBA = state[Press_Type].boxArray();
            MultiFab pMF(pBA,dmap,1,NGrow);

            if (time == getLevel(level-1).state[Press_Type].prevTime() ||
                time == getLevel(level-1).state[Press_Type].curTime())
            {
                FillCoarsePatch(pMF,0,time,Press_Type,0,1,NGrow);
            }
            else
            {
                Real crse_time;

                if (time > getLevel(level-1).state[State_Type].prevTime())
                {
                    crse_time = getLevel(level-1).state[Press_Type].curTime();
                }
                else
                {
                    crse_time = getLevel(level-1).state[Press_Type].prevTime();
                }

                FillCoarsePatch(pMF,0,crse_time,Press_Type,0,1,NGrow);

                MultiFab dpdtMF(pBA,dmap,1,NGrow);

                FillCoarsePatch(dpdtMF,0,time,Dpdt_Type,0,1,NGrow);

                Real dt_temp = time - crse_time;

                dpdtMF.mult(dt_temp,0,1,NGrow);

                pMF.plus(dpdtMF,0,1,NGrow);
            }
#ifdef _OPENMP
#pragma omp parallel
#endif
	    for (MFIter mfi(gp, true); mfi.isValid(); ++mfi)
            {
              const Box& bx=mfi.growntilebox();
              Projection::getGradP(pMF[mfi],gp[mfi],bx,dx);
            }
        }
        //
        // We've now got good coarse data everywhere in gp.
        //
        MultiFab gpTmp(gp.boxArray(),gp.DistributionMap(),1,NGrow);

	{

	  FillPatchIterator P_fpi(*this,P_old,NGrow,time,Press_Type,0,1);
	  MultiFab& pMF = P_fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel
#endif
	  for (MFIter mfi(gpTmp, true); mfi.isValid(); ++mfi)
	  {
      	    const Box& bx=mfi.growntilebox();
	    Projection::getGradP(pMF[mfi],gpTmp[mfi],bx,dx);
	  }
	}
        //
        // Now must decide which parts of gpTmp to copy to gp.
        //
        const int M = old_intersect_new.size();

        BoxArray fineBA(M);

        for (int j = 0; j < M; j++)
        {
            Box bx = old_intersect_new[j];

            for (int i = 0; i < BL_SPACEDIM; i++)
            {
                if (!geom.isPeriodic(i))
                {
                    if (bx.smallEnd(i) == geom.Domain().smallEnd(i))
                        bx.growLo(i,NGrow);
                    if (bx.bigEnd(i) == geom.Domain().bigEnd(i))
                        bx.growHi(i,NGrow);
                }
            }

            fineBA.set(j,bx);
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
	  std::vector< std::pair<int,Box> > isects;

	  for (MFIter mfi(gpTmp,true); mfi.isValid(); ++mfi)
	  {
            fineBA.intersections(mfi.growntilebox(),isects);

            FArrayBox&       gpfab    =    gp[mfi];
            const FArrayBox& gptmpfab = gpTmp[mfi];

            for (int ii = 0, N = isects.size(); ii < N; ii++)
            {
                gpfab.copy<RunOn::Host>(gptmpfab,isects[ii].second);
            }
	  }
	}

	gp.EnforcePeriodicity(geom.periodicity());
    }
    else
    {

        FillPatchIterator P_fpi(*this,P_old,NGrow,time,Press_Type,0,1);
	MultiFab& pMF = P_fpi.get_mf();
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(gp, true); mfi.isValid(); ++mfi)
        {
	  BL_ASSERT(amrex::grow(grids[mfi.index()],NGrow) == gp[mfi].box());

	  Projection::getGradP(pMF[mfi],gp[mfi],mfi.growntilebox(),dx);
	}
    }
}

//
// Fill patch a state component.
//
MultiFab*
NavierStokesBase::getState (int  ngrow,
			    int  state_idx,
			    int  scomp,
			    int  ncomp,
			    Real time)
{
    BL_PROFILE("NavierStokesBase::getState()");

    MultiFab* mf = new MultiFab(state[state_idx].boxArray(),
                                state[state_idx].DistributionMap(),
                                ncomp,ngrow,MFInfo(),Factory());

    FillPatch(*this,*mf,ngrow,time,state_idx,scomp,ncomp,0);

    return mf;
}

void
NavierStokesBase::getOutFlowFaces (Vector<Orientation>& outFaces)
{
    outFaces.resize(0);
    for (int idir = 0; idir < BL_SPACEDIM; idir++)
    {
        if (phys_bc.lo(idir) == Outflow)
        {
            const int len = outFaces.size();
            outFaces.resize(len+1);
            outFaces[len] = Orientation(idir,Orientation::low);
        }

        if (phys_bc.hi(idir) == Outflow)
        {
            const int len = outFaces.size();
            outFaces.resize(len+1);
            outFaces[len] = Orientation(idir,Orientation::high);
        }
    }
}

void
NavierStokesBase::incrPAvg ()
{
    //
    // Increment p_avg with 1/ncycle times current pressure
    //
    MultiFab& P_new = get_new_data(Press_Type);

    Real alpha = 1.0/Real(parent->nCycle(level));

    MultiFab::Saxpy(p_avg,alpha,P_new,0,0,1,0);
}

void
NavierStokesBase::initRhoAvg (Real alpha)
{
    const MultiFab& S_new = get_new_data(State_Type);

    // Set to a ridiculous number just for debugging -- shouldn't need this otherwise
    rho_avg.setVal(1.e200);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter rho_avgmfi(rho_avg,true); rho_avgmfi.isValid(); ++rho_avgmfi)
    {
//  	const Box& bx = rho_avgmfi.growntilebox();
    	const Box& bx = rho_avgmfi.tilebox();
        FArrayBox& rhoavgfab = rho_avg[rho_avgmfi];

    	rhoavgfab.copy<RunOn::Host>(S_new[rho_avgmfi],bx,Density,bx,0,1);
        rhoavgfab.mult<RunOn::Host>(alpha,bx);
    }
}

void
NavierStokesBase::incrRhoAvg(const MultiFab& rho_incr,
                         int             sComp,
                         Real            alpha)
{
    MultiFab::Saxpy(rho_avg,alpha,rho_incr,sComp,0,1,0);
}

void
NavierStokesBase::incrRhoAvg (Real alpha)
{
    const MultiFab& S_new = get_new_data(State_Type);
    incrRhoAvg(S_new,Density,alpha);
}

//
// Fills a new level n with best level n and coarser data available.
//
void
NavierStokesBase::init (AmrLevel &old)
{
    NavierStokesBase* oldns = (NavierStokesBase*) &old;
    const Real    dt_new    = parent->dtLevel(level);
    const Real    cur_time  = oldns->state[State_Type].curTime();
    const Real    prev_time = oldns->state[State_Type].prevTime();
    const Real    dt_old    = cur_time - prev_time;
    MultiFab&     S_new     = get_new_data(State_Type);
    MultiFab&     P_new     = get_new_data(Press_Type);
    MultiFab&     P_old     = get_old_data(Press_Type);

    setTimeLevel(cur_time,dt_old,dt_new);

    const Real cur_pres_time = state[Press_Type].curTime();
    //
    // Get best state and pressure data.
    //
    FillPatch(old,S_new,0,cur_time,State_Type,0,NUM_STATE);
    //
    // Note: we don't need to worry here about using FillPatch because
    //       it will automatically use the "old dpdt" to interpolate,
    //       since we haven't yet defined a new pressure at the lower level.
    //
    {
	FillPatchIterator fpi(old,P_new,0,cur_pres_time,Press_Type,0,1);
	const MultiFab& mf_fpi = fpi.get_mf();
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(mf_fpi,true); mfi.isValid(); ++mfi)
	{
	  const Box& vbx  = mfi.tilebox();
	  const FArrayBox& pfab = mf_fpi[mfi];

	  P_old[mfi].copy<RunOn::Host>(pfab,vbx,0,vbx,0,1);
	  P_new[mfi].copy<RunOn::Host>(pfab,vbx,0,vbx,0,1);
	}
    }

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
    {
        MultiFab& Dpdt_new = get_new_data(Dpdt_Type);
	FillPatch(old,Dpdt_new,0,cur_pres_time,Dpdt_Type,0,1);
    }
    //
    // Get best divu and dSdt data.
    //
    if (have_divu)
    {
        MultiFab& Divu_new = get_new_data(Divu_Type);
	FillPatch(old,Divu_new,0,cur_time,Divu_Type,0,1);

        if (have_dsdt)
        {
            MultiFab& Dsdt_new = get_new_data(Dsdt_Type);
	    FillPatch(old,Dsdt_new,0,cur_time,Dsdt_Type,0,1);
        }
    }

    old_intersect_new          = amrex::intersect(grids,oldns->boxArray());
    is_first_step_after_regrid = true;
}

void
NavierStokesBase::init ()
{
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& P_new = get_new_data(Press_Type);

    BL_ASSERT(level > 0);

    const Vector<Real>& dt_amr = parent->dtLevel();
    Vector<Real>        dt_new(level+1);

    for (int lev = 0; lev < level; lev++)
        dt_new[lev] = dt_amr[lev];
    //
    // Guess new dt from new data (interpolated from coarser level).
    //
    const Real dt = dt_new[level-1]/Real(parent->MaxRefRatio(level-1));
    dt_new[level] = dt;

    parent->setDtLevel(dt_new);
    //
    // Compute dt based on old data.
    //
    NavierStokesBase& old   = getLevel(level-1);
    const Real    cur_time  = old.state[State_Type].curTime();
    const Real    prev_time = old.state[State_Type].prevTime();
    const Real    dt_old    = (cur_time-prev_time)/Real(parent->MaxRefRatio(level-1));

    setTimeLevel(cur_time,dt_old,dt);

    Real cur_pres_time = state[Press_Type].curTime();
    //
    // Get best coarse state and pressure data.
    //
    FillCoarsePatch(S_new,0,cur_time,State_Type,0,NUM_STATE);
    FillCoarsePatch(P_new,0,cur_pres_time,Press_Type,0,1);

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
        FillCoarsePatch(get_new_data(Dpdt_Type),0,cur_time,Dpdt_Type,0,1);

    initOldPress();

    //
    // Get best coarse divU and dSdt data.
    //
    if (have_divu)
    {
        FillCoarsePatch(get_new_data(Divu_Type),0,cur_time,Divu_Type,0,1);
        if (have_dsdt)
            FillCoarsePatch(get_new_data(Dsdt_Type),0,cur_time,Dsdt_Type,0,1);
    }
    old_intersect_new = grids;
}

void
NavierStokesBase::init_additional_state_types ()
{
    additional_state_types_initialized = 1;
    //
    // Set "Temp" from user's variable setup.
    //
    int dummy_State_Type;
    int have_temp = isStateVariable("temp", dummy_State_Type, Temp);
    have_temp &= (dummy_State_Type == State_Type);
    BL_ASSERT((do_temp && have_temp)  ||  (!do_temp && !have_temp));

    int _Divu = -1;
    int dummy_Divu_Type;
    //
    // FIXME-- have divu was already set to zero...
    // also check have_dsdt
    have_divu = 0;
    have_divu = isStateVariable("divu", dummy_Divu_Type, _Divu);
    have_divu = have_divu && dummy_Divu_Type == Divu_Type;
    if (verbose)
    {
        amrex::Print() << "NavierStokesBase::init_additional_state_types()::have_divu = "
                  << have_divu << '\n';
    }
    if (have_divu && _Divu!=Divu)
    {
        amrex::Print() << "divu must be 0-th Divu_Type component in the state\n";

        amrex::Abort("NavierStokesBase::init_additional_state_types()");
    }

    int _Dsdt = -1;
    int dummy_Dsdt_Type;
    have_dsdt = 0;
    have_dsdt = isStateVariable("dsdt", dummy_Dsdt_Type, _Dsdt);
    have_dsdt = have_dsdt && dummy_Dsdt_Type==Dsdt_Type;
    if (verbose)
    {
        amrex::Print() << "NavierStokesBase::init_additional_state_types()::have_dsdt = "
		       << have_dsdt << '\n';
    }
    if (have_dsdt && _Dsdt!=Dsdt)
    {
        amrex::Print() << "dsdt must be 0-th Dsdt_Type component in the state\n";

        amrex::Abort("NavierStokesBase::init_additional_state_types()");
    }
    if (have_dsdt && !have_divu)
    {
        amrex::Print() << "Must have divu in order to have dsdt\n";

        amrex::Abort("NavierStokesBase::init_additional_state_types()");
    }

    num_state_type = desc_lst.size();
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "NavierStokesBase::init_additional_state_types: num_state_type = "
		       << num_state_type << '\n';
    }
}

Real
NavierStokesBase::initialTimeStep ()
{
    Real returnDt = init_shrink*estTimeStep();

    amrex::Print() << "Multiplying dt by init_shrink: dt = "
		   << returnDt << '\n';
    return returnDt;
}

//
// Since the pressure solver always stores its estimate of the
// pressure solver in Pnew, we need to copy it to Pold at the start.
//
void
NavierStokesBase::initOldPress ()
{
    MultiFab& P_new = get_new_data(Press_Type);
    MultiFab& P_old = get_old_data(Press_Type);

    MultiFab::Copy(P_old, P_new, 0, 0, P_old.nComp(), P_old.nGrow());
}

void
NavierStokesBase::zeroNewPress ()
{
    get_new_data(Press_Type).setVal(0);
}

void
NavierStokesBase::zeroOldPress ()
{
    get_old_data(Press_Type).setVal(0);
}

void
NavierStokesBase::level_projector (Real dt,
				   Real time,
				   int  iteration)
{
    BL_PROFILE_REGION_START("R::NavierStokesBase::level_projector()");
    BL_PROFILE("NavierStokesBase::level_projector()");

    BL_ASSERT(iteration > 0);

    MultiFab& U_old = get_old_data(State_Type);
    MultiFab& U_new = get_new_data(State_Type);
    MultiFab& P_old = get_old_data(Press_Type);
    MultiFab& P_new = get_new_data(Press_Type);

    SyncRegister* crse_ptr = 0;

    if (level < parent->finestLevel() && do_sync_proj)
    {
        crse_ptr = &(getLevel(level+1).getSyncReg());
    }

    int        crse_dt_ratio  = (level > 0) ? parent->nCycle(level) : -1;
    const Real cur_pres_time  = state[Press_Type].curTime();
    const Real prev_pres_time = state[Press_Type].prevTime();

    projector->level_project(level,time,dt,cur_pres_time,prev_pres_time,
                             geom,U_old,U_new,P_old,P_new,
                             get_rho_half_time(),crse_ptr,sync_reg,
                             crse_dt_ratio,iteration,have_divu);

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
        calcDpdt();

    BL_PROFILE_REGION_STOP("R::NavierStokesBase::level_projector()");
}

void
NavierStokesBase::level_sync (int crse_iteration)
{
    BL_PROFILE_REGION_START("R::NavierStokesBase::level_sync()");
    BL_PROFILE("NavierStokesBase::level_sync()");

    const Real*     dx            = geom.CellSize();
    IntVect         ratio         = parent->refRatio(level);
    const int       finest_level  = parent->finestLevel();
    int             crse_dt_ratio = parent->nCycle(level);
    Real            dt            = parent->dtLevel(level);
    MultiFab&       pres          = get_new_data(Press_Type);
    MultiFab&       vel           = get_new_data(State_Type);
    SyncRegister&   rhs_sync_reg  = getLevel(level+1).getSyncReg();
    SyncRegister*   crsr_sync_ptr = 0;
    NavierStokesBase&   fine_level    = getLevel(level+1);
    MultiFab&       pres_fine     = fine_level.get_new_data(Press_Type);
    MultiFab&       vel_fine      = fine_level.get_new_data(State_Type);
    const BoxArray& finegrids     = vel_fine.boxArray();
    const DistributionMapping& finedmap = vel_fine.DistributionMap();

    if (level > 0)
        crsr_sync_ptr = &(getLevel(level).getSyncReg());
    //
    // Get boundary conditions.
    //
    const int N = grids.size();

    Vector<int*>         sync_bc(N);
    Vector< Vector<int> > sync_bc_array(N);

    for (int i = 0; i < N; i++)
    {
        sync_bc_array[i] = getBCArray(State_Type,i,Xvel,BL_SPACEDIM);
        sync_bc[i] = sync_bc_array[i].dataPtr();
    }

    //
    // Multilevel or single-level sync projection.
    //
    MultiFab& Rh = get_rho_half_time();
    MultiFab cc_rhs_crse, cc_rhs_fine;

    cc_rhs_crse.define(    grids,    dmap,1,1,MFInfo(),           Factory());
    cc_rhs_fine.define(finegrids,finedmap,1,1,MFInfo(),fine_level.Factory());
    cc_rhs_crse.setVal(0);
    cc_rhs_fine.setVal(0);

    MultiFab&         v_fine    = fine_level.get_new_data(State_Type);
    MultiFab&       rho_fine    = fine_level.rho_avg;
    const Geometry& crse_geom   = parent->Geom(level);
    const BoxArray& P_finegrids = pres_fine.boxArray();
    const DistributionMapping& P_finedmap = pres_fine.DistributionMap();

    MultiFab phi(P_finegrids,P_finedmap,1,1,MFInfo(),fine_level.Factory());
    MultiFab V_corr(finegrids,finedmap,BL_SPACEDIM,1,MFInfo(),fine_level.Factory());

    V_corr.setVal(0);
    //
    // If periodic, enforce periodicity on Vsync.
    //
    if (crse_geom.isAnyPeriodic()) {
      Vsync.FillBoundary(0, BL_SPACEDIM, crse_geom.periodicity());
    }
    //
    // Interpolate Vsync to fine grid correction in Vcorr.
    //
    SyncInterp(Vsync, level, V_corr, level+1, ratio,
	       0, 0, BL_SPACEDIM, 0 , dt, sync_bc.dataPtr());
    //
    // The multilevel projection.  This computes the projection and
    // adds in its contribution to levels (level) and (level+1).
    //
    Real  cur_crse_pres_time = state[Press_Type].curTime();
    Real prev_crse_pres_time = state[Press_Type].prevTime();

    NavierStokesBase& fine_lev   = getLevel(level+1);
    Real  cur_fine_pres_time = fine_lev.state[Press_Type].curTime();
    Real prev_fine_pres_time = fine_lev.state[Press_Type].prevTime();

    bool first_crse_step_after_initial_iters =
      (prev_crse_pres_time > state[State_Type].prevTime());

    bool pressure_time_is_interval =
      (state[Press_Type].descriptor()->timeType() == StateDescriptor::Interval);
    projector->MLsyncProject(level,pres,vel,cc_rhs_crse,
			     pres_fine,v_fine,cc_rhs_fine,
			     Rh,rho_fine,Vsync,V_corr,
			     phi,&rhs_sync_reg,crsr_sync_ptr,
			     dt,ratio,crse_iteration,crse_dt_ratio,
			     geom,pressure_time_is_interval,
			     first_crse_step_after_initial_iters,
			     cur_crse_pres_time,prev_crse_pres_time,
			     cur_fine_pres_time,prev_fine_pres_time);
    cc_rhs_crse.clear();
    cc_rhs_fine.clear();
    //
    // Correct pressure and velocities after the projection.
    //
    const int Nf = finegrids.size();

    ratio = IntVect::TheUnitVector();

    Vector<int*>         fine_sync_bc(Nf);
    Vector< Vector<int> > fine_sync_bc_array(Nf);

    for (int i = 0; i < Nf; i++)
    {
      fine_sync_bc_array[i] = getLevel(level+1).getBCArray(State_Type,
							   i,
							   Xvel,
							   BL_SPACEDIM);
      fine_sync_bc[i] = fine_sync_bc_array[i].dataPtr();
    }

    for (int lev = level+2; lev <= finest_level; lev++)
    {
      ratio                 *= parent->refRatio(lev-1);
      NavierStokesBase& flev = getLevel(lev);
      MultiFab&     P_new    = flev.get_new_data(Press_Type);
      MultiFab&     P_old    = flev.get_old_data(Press_Type);
      MultiFab&     U_new    = flev.get_new_data(State_Type);

      SyncInterp(V_corr, level+1, U_new, lev, ratio,
		 0, 0, BL_SPACEDIM, 1 , dt, fine_sync_bc.dataPtr());
      SyncProjInterp(phi, level+1, P_new, P_old, lev, ratio,
		     first_crse_step_after_initial_iters,
		     cur_crse_pres_time, prev_crse_pres_time);
    }

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
      calcDpdt();

    BL_PROFILE_REGION_STOP("R::NavierStokesBase::level_sync()");
}

void
NavierStokesBase::make_rho_prev_time ()
{
    const Real prev_time = state[State_Type].prevTime();

    FillPatch(*this,rho_ptime,1,prev_time,State_Type,Density,1,0);

#ifdef AMREX_USE_EB
    EB_set_covered(rho_ptime,COVERED_VAL);
#endif
}

void
NavierStokesBase::make_rho_curr_time ()
{
    const Real curr_time = state[State_Type].curTime();

    FillPatch(*this,rho_ctime,1,curr_time,State_Type,Density,1,0);

#ifdef AMREX_USE_EB
    EB_set_covered(rho_ctime,COVERED_VAL);
#endif
}

void
NavierStokesBase::mac_project (Real      time,
			       Real      dt,
			       MultiFab& Sold,
			       MultiFab* divu,
			       int       ngrow,
			       bool      increment_vel_register)
{
    BL_PROFILE_REGION_START("R::NavierStokesBase::mac_project()");
    BL_PROFILE("NavierStokesBase::mac_project()");

    if (verbose) amrex::Print() << "... mac_projection\n";

    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real strt_time = ParallelDescriptor::second();

    Vector<BCRec> density_math_bc = fetchBCArray(State_Type,Density,1);

    mac_projector->mac_project(level,u_mac,Sold,dt,time,*divu,have_divu,
                               density_math_bc[0], increment_vel_register);

    create_umac_grown(ngrow);

    if (verbose)
    {
        Real run_time    = ParallelDescriptor::second() - strt_time;
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	amrex::Print() << "NavierStokesBase:mac_project(): lev: "
		       << level
		       << ", time: " << run_time << '\n';
    }
    BL_PROFILE_REGION_STOP("R::NavierStokesBase::mac_project()");
}

void
NavierStokesBase::manual_tags_placement (TagBoxArray&    tags,
					 const Vector<IntVect>& bf_lev)
{
    Vector<Orientation> outFaces;
    getOutFlowFaces(outFaces);
    if (outFaces.size()>0)
    {
        for (int i=0; i<outFaces.size(); ++i)
        {
            const Orientation& outFace = outFaces[i];
            const int oDir = outFace.coordDir();
            const Box& crse_domain = amrex::coarsen(geom.Domain(),bf_lev[level]);
            const int mult = (outFace.isLow() ? +1 : -1);
            if (do_refine_outflow)
            {
                //
                // Refine entire outflow boundary if new boxes within grid_tol
                // from outflow
                //
                const int grid_tol = 1;

                Box outflowBox = amrex::adjCell(crse_domain,outFace,grid_tol);

                outflowBox.shift(oDir,mult*grid_tol);

                //
                // Only refine if there are already tagged cells in the outflow
                // region
                //
                bool hasTags = false;
                for (MFIter tbi(tags); !hasTags && tbi.isValid(); ++tbi)
                    if (tags[tbi].numTags(outflowBox) > 0)
                        hasTags = true;
		ParallelAllReduce::Or(hasTags, ParallelContext::CommunicatorSub());
                if (hasTags)
                    tags.setVal(BoxArray(&outflowBox,1),TagBox::SET);
	    }
            else if (do_derefine_outflow)
            {
                const int np = parent->nProper();
                //
                // Calculate the number of level 0 cells to be left uncovered
                // at the outflow.  The convoluted logic allows for the fact that
                // the number of uncovered cells must be a multiple of the level
                // blocking factor.  So, when calculating the number of coarse
                // cells below, we always round the division up.
                //
                int N_coarse_cells = Nbuf_outflow / bf_lev[0][oDir];
                if (Nbuf_outflow % bf_lev[0][oDir] != 0)
                    N_coarse_cells++;

                int N_level_cells = N_coarse_cells * bf_lev[0][oDir];

                //
                // Adjust this to get the number of cells to be left uncovered at
                // levels higher than 0
                //
                for (int j = 1; j <= level; ++j)
                {
                    /*** Calculate the minimum cells at this level ***/

                    const int rat = (parent->refRatio(j-1))[oDir];
                    N_level_cells = N_level_cells * rat + np;

                    /*** Calculate the required number of coarse cells ***/

                    N_coarse_cells = N_level_cells / bf_lev[j][oDir];
                    if (N_level_cells % bf_lev[j][oDir] != 0)
                        N_coarse_cells++;

                    /*** Calculate the corresponding number of level cells ***/

                    N_level_cells = N_coarse_cells * bf_lev[j][oDir];
                }
                //
                // Untag the cells near the outflow
                //
                if (N_coarse_cells > 0)
                {
                    //
                    // Generate box at the outflow and grow it in all directions
                    // other than the outflow.  This forces outflow cells in the
                    // ghostcells in directions other that oDir to be cleared.
                    //
                    Box outflowBox = amrex::adjCell(crse_domain, outFace, 1);
                    for (int dir = 0; dir < BL_SPACEDIM; dir++)
                        if (dir != oDir) outflowBox.grow(dir, 1);
                    //
                    // Now, grow the box into the domain (opposite direction as
                    // outFace) the number of cells we need to clear.
                    //
                    if (outFace.isLow())
                        outflowBox.growHi(oDir, N_coarse_cells);
                    else
                        outflowBox.growLo(oDir, N_coarse_cells);

                    tags.setVal(BoxArray(&outflowBox,1),TagBox::CLEAR);
                }
            }
        }
    }
}

int
NavierStokesBase::okToContinue ()
{
	//
	// Check that dt is OK across AMR levels
	//
  	int okLevel = (level > 0) ? true : (parent->dtLevel(0) > dt_cutoff);

	if (stop_when_steady)
		//
		// If stop_when_steady is enabled, also check that we haven't reached
		// steady-state.
		//
		return (okLevel && !steadyState());
	else
	  	return okLevel;
}

int
NavierStokesBase::steadyState()
{
    if (!get_state_data(State_Type).hasOldData()) {
        return false; // If nothing to compare to, must not yet be steady :)
    }

    Real        max_change    = 0.0;
    MultiFab&   U_old         = get_old_data(State_Type);
    MultiFab&   U_new         = get_new_data(State_Type);

	//
	// Estimate the maximum change in velocity magnitude since previous
	// iteration
	//
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(max:max_change)
#endif
    for (MFIter Rho_mfi(rho_ctime,true); Rho_mfi.isValid(); ++Rho_mfi)
    {
      const Box& bx=Rho_mfi.tilebox();
      Real change = godunov->maxchng_velmag(U_new[Rho_mfi],U_old[Rho_mfi],bx);

      max_change = std::max(change, max_change);
    }

    ParallelDescriptor::ReduceRealMax(max_change);

	//
	// System is classified as steady if the maximum change is smaller than
	// prescribed tolerance
	//
    bool steady = max_change < steady_tol;

    if (verbose)
    {
        amrex::Print() << "steadyState :: \n" << "LEV = " << level
                       << " MAX_CHANGE = " << max_change << std::endl;

        if (steady)
        {
            amrex::Print()
                << "System reached steady-state, stopping simulation."
                << std::endl;
        }
    }

    return steady;
}

//
// This function estimates the initial timesteping used by the model.
//
void
NavierStokesBase::post_init_estDT (Real&        dt_init,
				   Vector<int>&  nc_save,
				   Vector<Real>& dt_save,
				   Real         stop_time)
{
    const Real strt_time    = state[State_Type].curTime();
    const int  finest_level = parent->finestLevel();

    dt_init = 1.0e+100;

    int  n_factor;
    for (int k = 0; k <= finest_level; k++)
    {
        nc_save[k] = parent->nCycle(k);
        dt_save[k] = getLevel(k).initialTimeStep();

        n_factor   = 1;
        for (int m = finest_level; m > k; m--)
             n_factor *= parent->nCycle(m);
        dt_init    = std::min( dt_init, dt_save[k]/((Real) n_factor) );
    }

    Vector<Real> dt_level(finest_level+1,dt_init);
    Vector<int>  n_cycle(finest_level+1,1);

    Real dt0 = dt_save[0];
    n_factor = 1;
    for (int k = 0; k <= finest_level; k++)
    {
        n_factor *= nc_save[k];
        dt0       = std::min(dt0,n_factor*dt_save[k]);
    }

    if (stop_time >= 0.0)
    {
        const Real eps = 0.0001*dt0;
        if ((strt_time + dt0) > (stop_time - eps))
            dt0 = stop_time - strt_time;
    }

    n_factor = 1;
    for (int k = 0; k <= finest_level; k++)
    {
        n_factor  *= nc_save[k];
        dt_save[k] = dt0/( (Real) n_factor);
    }
    //
    // Hack.
    //
    parent->setDtLevel(dt_level);
    parent->setNCycle(n_cycle);
    for (int k = 0; k <= finest_level; k++)
    {
        getLevel(k).setTimeLevel(strt_time,dt_init,dt_init);
    }
}

//
// This function ensures that the state is initially consistent
// with respect to the divergence condition and fields are initially consistent
//
void
NavierStokesBase::post_init_state ()
{
    const int finest_level = parent->finestLevel();
    const Real divu_time   = have_divu ? state[Divu_Type].curTime()
                                       : state[Press_Type].curTime();

    //
    // Make sure we're not trying to use ref_ratio=4
    // MLMG does not support rr=4 yet
    //
    // Derived class PeleLM seems to also use this function , so it's a
    // convienient place for the test, even though it's not the most
    // logical place to put the check
    int maxlev = parent->maxLevel();
    for (int i = 0; i<maxlev; i++)
    {
      const int rr = parent->MaxRefRatio(i);
      if (rr == 4)
      {
	amrex::Abort("Refinement ratio of 4 not currently supported.\n");
      }
    }

    if (do_init_vort_proj)
    {
        //
	// NOTE: this assumes have_divu == 0.
	// Only used if vorticity is used to initialize the velocity field.
        //
        BL_ASSERT(!(projector == 0));

	if (verbose) amrex::Print() << "calling initialVorticityProject" << std::endl;

	projector->initialVorticityProject(0);

	if (verbose) amrex::Print() << "done calling initialVorticityProject" << std::endl;
    }

    if (do_init_proj && projector)
    {
      //
      // Do sync project to define divergence free velocity field.
      //
      if (verbose) amrex::Print() << "calling initialVelocityProject" << std::endl;

      projector->initialVelocityProject(0,divu_time,have_divu,init_vel_iter);

      if (verbose) amrex::Print() << "done calling initialVelocityProject" << std::endl;
    }

    NavierStokesBase::initial_step = true;
    //
    // Average velocity and scalar data down from finer levels
    // so that conserved data is consistant between levels.
    //
    for (int k = finest_level-1; k>= 0; k--)
    {
      getLevel(k).avgDown();
    }
    make_rho_curr_time();

    if (do_init_proj && projector && (std::fabs(gravity)) > 0.){
      //
      // Do projection to establish initially hydrostatic pressure field.
      //
      if (verbose) amrex::Print() << "calling initialPressureProject" << std::endl;

      projector->initialPressureProject(0);

      if (verbose) amrex::Print() << "done calling initialPressureProject" << std::endl;
    }
    // make sure there's not NANs in old pressure field
    // end up with P_old = P_new as is the case when exiting initialPressureProject
    if(!do_init_proj){
      MultiFab& p_old=get_old_data(Press_Type);
      MultiFab& p_new=get_new_data(Press_Type);
      MultiFab::Copy(p_old, p_new, 0, 0, 1, p_new.nGrow());
    }

}

//
// Build any additional data structures after regrid.
//
void
NavierStokesBase::post_regrid (int lbase,
			       int new_finest)
{
#ifdef AMREX_PARTICLES
    if (NSPC && level == lbase)
    {
        NSPC->Redistribute(lbase);
    }
#endif
}

//
// Build any additional data structures after restart.
//
void
NavierStokesBase::post_restart ()
{
    make_rho_prev_time();
    make_rho_curr_time();

#ifdef AMREX_PARTICLES
    post_restart_particle ();
#endif
}

//
// Integration cycle on fine level grids is complete .
// post_timestep() is responsible for syncing levels together.
//
// The registers used for level syncing are initialized in the
// coarse level advance and incremented in the fine level advance.
// These quantities are described in comments above advance_setup.
//
void
NavierStokesBase::post_timestep (int crse_iteration)
{

  BL_PROFILE("NavierStokesBase::post_timestep()");

    const int finest_level = parent->finestLevel();

#ifdef AMREX_PARTICLES
    post_timestep_particle (crse_iteration);
#endif

    if (level == parent->finestLevel())
    {
        delete [] u_mac;
        u_mac = 0;
    }

    if (do_reflux && level < finest_level)
        reflux();

    if (level < finest_level)
        avgDown();

    if (do_mac_proj && level < finest_level)
        mac_sync();

    if (do_sync_proj && (level < finest_level))
        level_sync(crse_iteration);
    //
    // Test for conservation.
    //
    if (level==0 && sum_interval>0 && (parent->levelSteps(0)%sum_interval == 0))
    {
        sum_integrated_quantities();
    }
#if (BL_SPACEDIM==3)
    //
    // Derive turbulent statistics
    //
    if (level==0 && turb_interval>0 && (parent->levelSteps(0)%turb_interval == 0))
    {
        sum_turbulent_quantities();
    }
#ifdef SUMJET
    //
    // Derive turbulent statistics for the round jet
    //
    if (level==0 && jet_interval>0 && (parent->levelSteps(0)%jet_interval == 0))
    {
        sum_jet_quantities();
    }
#endif
#endif

    if (level > 0) incrPAvg();

    old_intersect_new          = grids;
    is_first_step_after_regrid = false;

    if (level == 0 && dump_plane >= 0)
    {
        Box bx = geom.Domain();

        BL_ASSERT(bx.bigEnd(BL_SPACEDIM-1) >= dump_plane);

        bx.setSmall(BL_SPACEDIM-1, dump_plane);
        bx.setBig  (BL_SPACEDIM-1, dump_plane);

        BoxArray ba(bx);
        DistributionMapping dm{ba};

        MultiFab mf(ba, dm, BL_SPACEDIM, 0, MFInfo(), Factory());

        mf.copy(get_new_data(State_Type), Xvel, 0, BL_SPACEDIM);

        if (ParallelDescriptor::MyProc() == mf.DistributionMap()[0])
        {
            char buf[64];
            sprintf(buf, "%14.12e", state[State_Type].curTime());

            std::string name(dump_plane_name);
            name += buf;
            name += ".fab";

            std::ofstream ofs;
            ofs.open(name.c_str(),std::ios::out|std::ios::trunc|std::ios::binary);
            if (!ofs.good())
                amrex::FileOpenFailed(name);

            mf[0].writeOn(ofs);
        }
    }
}

//
// Reset the time levels to time (time) and timestep dt.
// This is done at the start of the timestep in the pressure iteration section.
//
void
NavierStokesBase::resetState (Real time,
			      Real dt_old,
			      Real dt_new)
{
    //
    // Reset state types.
    //
    state[State_Type].reset();
    state[State_Type].setTimeLevel(time,dt_old,dt_new);

    initOldPress();
    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Interval)
    {
        state[Press_Type].setTimeLevel(time-dt_old,dt_old,dt_new);
    }
    else if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
    {
        state[Press_Type].setTimeLevel(time-.5*dt_old,dt_old,dt_old);
        state[Dpdt_Type].setTimeLevel(time-dt_old,dt_old,dt_old);
    }
    //
    // Reset state types for divu not equal to zero.
    //
    if (have_divu)
    {
        state[Divu_Type].reset();
        state[Divu_Type].setTimeLevel(time,dt_old,dt_new);
        if (have_dsdt)
        {
            //
            // Dont do this, we want to improve dsdt with press iters
            // but we do need to make sure time is set correctly..
            // state[Dsdt_Type].reset();
            state[Dsdt_Type].setTimeLevel(time,dt_old,dt_new);
        }
    }
}

void
NavierStokesBase::restart (Amr&          papa,
                       std::istream& is,
                       bool          bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

#ifdef AMREX_USE_EB
    amrex::Warning("Restart not tested with EB yet.");
#endif
    //
    // Build metric coefficients for RZ calculations.
    // Build volume and areas.
    //
    buildMetrics();

    if (projector == 0)
    {
        projector = new Projection(parent,&phys_bc,do_sync_proj,
                                   parent->finestLevel(),radius_grow);
    }
    projector->install_level(level, this, &radius);
    //
    // Set the godunov box.
    //
    SetGodunov();

    if (mac_projector == 0)
    {
        mac_projector = new MacProj(parent,parent->finestLevel(),
                                    &phys_bc,radius_grow);
    }
    mac_projector->install_level(level,this);

    const BoxArray& P_grids = state[Press_Type].boxArray();
#ifdef AMREX_USE_EB
    init_eb(parent->Geom(level), grids, dmap);

    //fixme? not 100% sure this is the right place
    // note --- this fn is really similar to constructor
    //  need to make sure gradp is getting properly filled for this restart case?
    //  incflo style advection does not use Gp in tracing states to edges,
    //  don't think Gp is needed until the vel update after the projection.
    // But ultimately we need gradp in the checkpoint file.
    gradp.reset(new MultiFab(grids,dmap,BL_SPACEDIM,1, MFInfo(), Factory()));

    std::string file=papa.theRestartFile();
    std::string LevelDir, FullPath;
    LevelDirectoryNames(file, LevelDir, FullPath);
    std::string gradp_mf_fullpath = FullPath;
    gradp_mf_fullpath += "/gradp";
    const char *faHeader = 0;
    VisMF::Read(*gradp, gradp_mf_fullpath, faHeader);
#endif

    //
    // Alloc space for density and temporary pressure variables.
    //
    if (level > 0)
    {
        rho_avg.define(grids,dmap,1,1,MFInfo(),Factory());
        p_avg.define(P_grids,dmap,1,0,MFInfo(),Factory());
    }
    //FIXME see similar stuff in constructor
    rho_half.define (grids,dmap,1,1,MFInfo(),Factory());
    rho_ptime.define(grids,dmap,1,1,MFInfo(),Factory());
    rho_ctime.define(grids,dmap,1,1,MFInfo(),Factory());
    rho_qtime  = 0;
    rho_tqtime = 0;

    BL_ASSERT(sync_reg == 0);
    if (level > 0 && do_sync_proj)
    {
        sync_reg = new SyncRegister(grids,dmap,crse_ratio);
    }
    BL_ASSERT(advflux_reg == 0);
    if (level > 0 && do_reflux)
    {
        advflux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
    }
    BL_ASSERT(viscflux_reg == 0);
    if (level > 0 && do_reflux)
    {
        viscflux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
    }

    if (level < parent->finestLevel())
    {
        Vsync.define(grids,dmap,BL_SPACEDIM,1,MFInfo(),Factory());
        Ssync.define(grids,dmap,NUM_STATE-BL_SPACEDIM,1,MFInfo(),Factory());
    }

    diffusion = new Diffusion(parent, this,
                              (level > 0) ? getLevel(level-1).diffusion : 0,
                              NUM_STATE, viscflux_reg,is_diffusive, visc_coef);
    //
    // Allocate the storage for variable viscosity and diffusivity
    //
    viscn   = fb_viscn.define(this,1,0);
    viscnp1 = fb_viscnp1.define(this,1,0);
    diffn   = fb_diffn.define(this,NUM_STATE-Density-1,0);
    diffnp1 = fb_diffnp1.define(this,NUM_STATE-Density-1,0);

    is_first_step_after_regrid = false;
    old_intersect_new          = grids;
}

void
NavierStokesBase::scalar_advection_update (Real dt,
					   int  first_scalar,
					   int  last_scalar)
{
    BL_PROFILE("NavierStokesBase::scalar_advection_update()");

    MultiFab&  S_old     = get_old_data(State_Type);
    MultiFab&  S_new     = get_new_data(State_Type);
    MultiFab&  Aofs      = *aofs;

    const Real prev_time = state[State_Type].prevTime();


    //
    // Compute inviscid estimate of scalars.
    // (do rho separate, as we do not have rho at new time yet)
    //
    int sComp = first_scalar;

    if (sComp == Density)
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
{
      FArrayBox  tforces;
      for (MFIter S_oldmfi(S_old,true); S_oldmfi.isValid(); ++S_oldmfi)
      {

	    const Box& bx = S_oldmfi.tilebox();
            tforces.resize(bx,1);
            tforces.setVal<RunOn::Host>(0);
            godunov->Add_aofs_tf(S_old[S_oldmfi],S_new[S_oldmfi],Density,1,
                                 Aofs[S_oldmfi],Density,tforces,0,bx,dt);
      }
}
        //
        // Call ScalMinMax to avoid overshoots in density.
        //
      if (do_denminmax)
      {
	    //
            // Must do FillPatch here instead of MF iterator because we need the
            // boundary values in the old data (especially at inflow)
            //
            const int index_new_s   = Density;
            const int index_new_rho = Density;
            const int index_old_s   = index_new_s   - Density;
            const int index_old_rho = index_new_rho - Density;

            FillPatchIterator S_fpi(*this,S_old,1,prev_time,State_Type,Density,1);
            MultiFab& Smf=S_fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel
#endif
{
            Vector<int> state_bc;
	    for (MFIter mfi(Smf,true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                state_bc = fetchBCArray(State_Type,bx,Density,1);
                godunov->ConservativeScalMinMax(Smf[mfi],S_new[mfi],
                                                index_old_s, index_old_rho,
                                                index_new_s, index_new_rho,
                                                state_bc.dataPtr(),bx);
            }
}
      }
      ++sComp;
    }

    if (sComp <= last_scalar)
    {
        const MultiFab& rho_halftime = get_rho_half_time();
#ifdef _OPENMP
#pragma omp parallel
#endif
{
        FArrayBox  tforces;

        for (MFIter Rho_mfi(rho_halftime,true); Rho_mfi.isValid(); ++Rho_mfi)
        {
            const Box& bx = Rho_mfi.tilebox();

            for (int sigma = sComp; sigma <= last_scalar; sigma++)
            {
               // Need to do some funky half-time stuff
               if (getForceVerbose)
                  amrex::Print() << "---" << '\n' << "E - scalar advection update (half time):" << '\n';

               // Average the mac face velocities to get cell centred velocities
               const Real halftime = 0.5*(state[State_Type].curTime()+state[State_Type].prevTime());
               FArrayBox Vel(amrex::grow(bx,0),AMREX_SPACEDIM);
               FORT_AVERAGE_EDGE_STATES(BL_TO_FORTRAN_ANYD(Vel),
                                        BL_TO_FORTRAN_ANYD(u_mac[0][Rho_mfi]),
                                        BL_TO_FORTRAN_ANYD(u_mac[1][Rho_mfi]),
#if (AMREX_SPACEDIM==3)
                                        BL_TO_FORTRAN_ANYD(u_mac[2][Rho_mfi]),
#endif
                                        &getForceVerbose);

               //
               // Average the new and old time to get Crank-Nicholson half time approximation.
               //
               FArrayBox Scal(amrex::grow(bx,0),NUM_SCALARS);
               Scal.copy<RunOn::Host>(S_old[Rho_mfi],bx,Density,bx,0,NUM_SCALARS);
               Scal.plus<RunOn::Host>(S_new[Rho_mfi],bx,Density,0,NUM_SCALARS);
               Scal.mult<RunOn::Host>(0.5,bx);

               if (getForceVerbose) amrex::Print() << "Calling getForce..." << '\n';
                  getForce(tforces,bx,0,sigma,1,halftime,Vel,Scal,0);

               godunov->Add_aofs_tf(S_old[Rho_mfi],S_new[Rho_mfi],sigma,1,
                                    Aofs[Rho_mfi],sigma,tforces,0,bx,dt);
            }
        }
}
    }
    //
    // Call ScalMinMax to avoid overshoots in the scalars.
    //

    if ( do_scalminmax && (sComp <= last_scalar) )
    {
        const int num_scalars = last_scalar - Density + 1;
        //
        // Must do FillPatch here instead of MF iterator because we need the
        // boundary values in the old data (especially at inflow).
        //

        FillPatchIterator S_fpi(*this,S_old,1,prev_time,State_Type,Density,num_scalars);
        MultiFab& Smf=S_fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel
#endif
{
        Vector<int> state_bc;
	for (MFIter mfi(Smf,true); mfi.isValid();++mfi)
        {
            const Box& bx = mfi.tilebox();
            for (int sigma = sComp; sigma <= last_scalar; sigma++)
            {
                const int index_new_s   = sigma;
                const int index_new_rho = Density;
                const int index_old_s   = index_new_s   - Density;
                const int index_old_rho = index_new_rho - Density;

                state_bc = fetchBCArray(State_Type,bx,sigma,1);
                if (advectionType[sigma] == Conservative)
                {
		    godunov->ConservativeScalMinMax(Smf[mfi],S_new[mfi],
                                                    index_old_s, index_old_rho,
                                                    index_new_s, index_new_rho,
                                                    state_bc.dataPtr(),bx);
                }
                else if (advectionType[sigma] == NonConservative)
                {
                    godunov->ConvectiveScalMinMax(Smf[mfi],S_new[mfi],index_old_s,sigma,
                                                  state_bc.dataPtr(),bx);
                }
            }
        }
}
    }
}

//
// Set the time levels to time (time) and timestep dt.
//
void
NavierStokesBase::setTimeLevel (Real time,
				Real dt_old,
				Real dt_new)
{
    state[State_Type].setTimeLevel(time,dt_old,dt_new);

    if (have_divu)
    {
        state[Divu_Type].setTimeLevel(time,dt_old,dt_new);
        if (have_dsdt)
        {
            state[Dsdt_Type].setTimeLevel(time,dt_old,dt_new);
        }
    }

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Interval)
    {
        state[Press_Type].setTimeLevel(time-dt_old,dt_old,dt_old);
    }
    else if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
    {
        state[Press_Type].setTimeLevel(time-.5*dt_old,dt_old,dt_old);
        state[Dpdt_Type].setTimeLevel(time-dt_old,dt_old,dt_old);
    }
}

void
NavierStokesBase::sync_setup (MultiFab*& DeltaSsync)
{
    BL_ASSERT(DeltaSsync == 0);

    int nconserved = Godunov::how_many(advectionType, Conservative,
                                       BL_SPACEDIM, NUM_STATE-BL_SPACEDIM);

    if (nconserved > 0 && level < parent->finestLevel())
    {
        DeltaSsync = new MultiFab(grids, dmap, nconserved, 1, MFInfo(), Factory());
        DeltaSsync->setVal(0,1);
    }
}

void
NavierStokesBase::sync_cleanup (MultiFab*& DeltaSsync)
{
    delete DeltaSsync;

    DeltaSsync = 0;
}

//
// Helper function for NavierStokesBase::SyncInterp().
//
static
void
set_bc_new (int*            bc_new,
            int             n,
            int             src_comp,
            const int*      clo,
            const int*      chi,
            const int*      cdomlo,
            const int*      cdomhi,
            const BoxArray& cgrids,
            int**           bc_orig_qty)

{
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        int bc_index = (n+src_comp)*(2*BL_SPACEDIM) + dir;
        bc_new[bc_index]             = INT_DIR;
        bc_new[bc_index+BL_SPACEDIM] = INT_DIR;

        if (clo[dir] < cdomlo[dir] || chi[dir] > cdomhi[dir])
        {
            for (int crse = 0, N = cgrids.size(); crse < N; crse++)
            {
		const Box& bx = cgrids[crse];
                const int* c_lo = bx.loVect();
                const int* c_hi = bx.hiVect();

                if (clo[dir] < cdomlo[dir] && c_lo[dir] == cdomlo[dir])
                    bc_new[bc_index] = bc_orig_qty[crse][bc_index];
                if (chi[dir] > cdomhi[dir] && c_hi[dir] == cdomhi[dir])
                    bc_new[bc_index+BL_SPACEDIM] = bc_orig_qty[crse][bc_index+BL_SPACEDIM];
            }
        }
    }
}

//
// FIXME filcc eventually needs attention, see comments below
//
// Interpolate A cell centered Sync correction from a
// coarse level (c_lev) to a fine level (f_lev).
//
// This routine interpolates the num_comp components of CrseSync
// (starting at src_comp) and either increments or puts the result into
// the num_comp components of FineSync (starting at dest_comp)
// The components of bc_orig_qty corespond to the quantities of CrseSync.
//
void
NavierStokesBase::SyncInterp (MultiFab&      CrseSync,
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
			      SyncInterpType which_interp,
			      int            state_comp)
{
    BL_PROFILE("NavierStokesBase::SyncInterp()");

    BL_ASSERT(which_interp >= 0 && which_interp <= 5);

    Interpolater* interpolater = 0;

#ifdef AMREX_USE_EB
    switch (which_interp)
    {
      // As with the non-EB case, both of these point to the same interpolater
    case CellCons_T:     interpolater = &eb_cell_cons_interp;    break;
    case CellConsLin_T:  interpolater = &eb_lincc_interp;        break;
    default:
      amrex::Abort("NavierStokesBase::SyncInterp(): EB currently requires Cell Conservative interpolater. \n");
    }
#else
    switch (which_interp)
    {
    case PC_T:           interpolater = &pc_interp;           break;
    case CellCons_T:     interpolater = &cell_cons_interp;    break;
    case CellConsLin_T:  interpolater = &lincc_interp;        break;
    case CellConsProt_T: interpolater = &protected_interp;    break;
    default:
        amrex::Abort("NavierStokesBase::SyncInterp(): how did this happen \n");
    }
#endif

    NavierStokesBase& fine_level = getLevel(f_lev);
    const BoxArray& fgrids     = fine_level.boxArray();
    const DistributionMapping& fdmap = fine_level.DistributionMap();
    const Geometry& fgeom      = parent->Geom(f_lev);
    const BoxArray& cgrids     = getLevel(c_lev).boxArray();
    const Geometry& cgeom      = parent->Geom(c_lev);
    const Real*     dx_crse    = cgeom.CellSize();
    Box             cdomain    = amrex::coarsen(fgeom.Domain(),ratio);
    const int*      cdomlo     = cdomain.loVect();
    const int*      cdomhi     = cdomain.hiVect();
    const int       N          = fgrids.size();

    BoxArray cdataBA(N);

    for (int i = 0; i < N; i++)
        cdataBA.set(i,interpolater->CoarseBox(fgrids[i],ratio));
    //
    // Note: The boxes in cdataBA may NOT be disjoint !!!
    //
#ifdef AMREX_USE_EB
    // I am unsure of EBSupport and ng (set to zero here)
    auto factory = makeEBFabFactory(cgeom,cdataBA,fdmap,{0,0,0},EBSupport::basic);
    MultiFab cdataMF(cdataBA,fdmap,num_comp,0,MFInfo(),*factory);
#else
    //    ,MFInfo(),getLevel(c_lev).Factory());
    MultiFab cdataMF(cdataBA,fdmap,num_comp,0);
#endif



    // Coarse box could expand beyond the extent of fine box depending on the interpolation type, so initialize here
    cdataMF.setVal(0);

    cdataMF.copy(CrseSync, src_comp, 0, num_comp);
    //
    // Set physical boundary conditions in cdataMF.
    //
    //////////
    // Should be fine for EB for now, since EB doesn't intersect Phys BC
    // Although calling filcc_tile is not what we really want to be doing,
    // there must be something higher level
    // Maybe setPhysBoundaryValues? see example below
    //
    // MultiFab test(cdataBA,fdmap,num_comp,0,MFInfo(),Factory());
    // // replace divu with test and compare with cdataMF?
    // divu_old.FillBoundary();
    // Real prev_time = amr_level.get_state_data(Divu_Type).prevTime();
    // for (MFIter mfi(divu_new); mfi.isValid(); ++mfi)
    // {
    //   amr_level.setPhysBoundaryValues(divu_old[mfi],Divu_Type,prev_time,0,0,1);
    // }
    ///////

    // tiling may not be needed here, but what the hey

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      int* bc_new = new int[2*BL_SPACEDIM*(src_comp+num_comp)];

      for (MFIter mfi(cdataMF,true); mfi.isValid(); ++mfi)
      {

        FArrayBox&  cdata   = cdataMF[mfi];
        const int*  clo     = cdata.loVect();
        const int*  chi     = cdata.hiVect();
	const Box&  bx      = mfi.tilebox();
        RealBox     gridloc = RealBox(bx,fine_level.geom.CellSize(),fine_level.geom.ProbLo());
        const int*  lo      = bx.loVect();
        const int*  hi      = bx.hiVect();
        const Real* xlo     = gridloc.lo();

        for (int n = 0; n < num_comp; n++)
        {
          set_bc_new(bc_new,n,src_comp,lo,hi,cdomlo,cdomhi,cgrids,bc_orig_qty);

	  filcc_tile(ARLIM(lo),ARLIM(hi),
		     cdata.dataPtr(n), ARLIM(clo), ARLIM(chi),
		     cdomlo, cdomhi, dx_crse, xlo,
		     &(bc_new[2*BL_SPACEDIM*(n+src_comp)]));
        }
      }
      delete [] bc_new;
    }

    cdataMF.EnforcePeriodicity(cgeom.periodicity());
    //
    // Interpolate from cdataMF to fdata and update FineSync.
    // Note that FineSync and cdataMF will have the same distribution
    // since the length of their BoxArrays are equal.
    //

    MultiFab* fine_stateMF = 0;
    if (interpolater == &protected_interp)
    {
        fine_stateMF = &(getLevel(f_lev).get_new_data(State_Type));
    }


#ifdef AMREX_USE_EB
    //FIXME?
    // there's currently no way to reassign the EBCellFlagFab for a EBFArrayBox,
    // so the non-EB strategy of creating one FAB and resizing it for MFiters doens't work
    const FabArray<EBCellFlagFab>& flags = dynamic_cast<EBFArrayBoxFactory const&>(getLevel(f_lev).Factory()).getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
#ifndef AMREX_USE_EB
      FArrayBox    fdata;
#endif
      Vector<BCRec> bc_interp(num_comp);
      int* bc_new = new int[2*BL_SPACEDIM*(src_comp+num_comp)];

      for (MFIter mfi(FineSync,true); mfi.isValid(); ++mfi)
      {
	  FArrayBox& cdata = cdataMF[mfi];
	  // cdataMF has no ghost cells
	  const Box&  fbx     = mfi.tilebox();
	  const Box cbx = interpolater->CoarseBox(fbx,ratio);
	  const int* clo   = cbx.loVect();
	  const int* chi   = cbx.hiVect();

#ifdef AMREX_USE_EB
	  EBFArrayBox fdata(flags[mfi],fbx,num_comp,FineSync[mfi].arena());
#else
	  fdata.resize(fbx, num_comp);
#endif
	  //
	  // Set the boundary condition array for interpolation.
	  //
	  for (int n = 0; n < num_comp; n++)
	  {
	      set_bc_new(bc_new,n,src_comp,clo,chi,cdomlo,cdomhi,cgrids,bc_orig_qty);
	  }

	  for (int n = 0; n < num_comp; n++)
	  {
	      for (int dir = 0; dir < BL_SPACEDIM; dir++)
	      {
		  int bc_index = (n+src_comp)*(2*BL_SPACEDIM) + dir;
		  bc_interp[n].setLo(dir,bc_new[bc_index]);
		  bc_interp[n].setHi(dir,bc_new[bc_index+BL_SPACEDIM]);
	      }
	  }

	  //        ScaleCrseSyncInterp(cdata, c_lev, num_comp);

	  interpolater->interp(cdata,0,fdata,0,num_comp,fbx,ratio,
			       cgeom,fgeom,bc_interp,src_comp,State_Type,RunOn::Cpu);

	  //        reScaleFineSyncInterp(fdata, f_lev, num_comp);

	  if (increment)
	  {
	      fdata.mult<RunOn::Host>(dt_clev);

	      if (interpolater == &protected_interp)
	      {
		  cdata.mult<RunOn::Host>(dt_clev,cbx);
		  FArrayBox& fine_state = (*fine_stateMF)[mfi];
		  interpolater->protect(cdata,0,fdata,0,fine_state,state_comp,
					num_comp,fbx,ratio,
					cgeom,fgeom,bc_interp,RunOn::Cpu);

		  Real dt_clev_inv = 1./dt_clev;
		  cdata.mult<RunOn::Host>(dt_clev_inv,cbx);
	      }

	      FineSync[mfi].plus<RunOn::Host>(fdata,fbx,0,dest_comp,num_comp);
	  }
	  else
	  {
	      FineSync[mfi].copy<RunOn::Host>(fdata,fbx,0,fbx,dest_comp,num_comp);
	  }
      }
    delete [] bc_new;
    }
}

//
// Interpolate sync pressure correction to a finer level.
//
void
NavierStokesBase::SyncProjInterp (MultiFab& phi,
				  int       c_lev,
				  MultiFab& P_new,
				  MultiFab& P_old,
				  int       f_lev,
				  IntVect&  ratio,
				  bool      first_crse_step_after_initial_iters,
				  Real      cur_crse_pres_time,
				  Real      prev_crse_pres_time)
{
    BL_PROFILE("NavierStokesBase:::SyncProjInterp()");

    const BoxArray& P_grids = P_new.boxArray();
    const int       N       = P_grids.size();

    BoxArray crse_ba(N);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        crse_ba.set(i,node_bilinear_interp.CoarseBox(P_grids[i],ratio));

    // None  of these 3 are actually used by node_bilinear_interp()
    Vector<BCRec> bc(BL_SPACEDIM);
    const Geometry& fgeom   = parent->Geom(f_lev);
    const Geometry& cgeom   = parent->Geom(c_lev);

#ifdef AMREX_USE_EB
    // I am unsure of EBSupport and ng (set to 1 here)
    // need 1 ghost cell to use EB_set_covered on nodal MF
    // Factory is always CC, regardless of status of crse_ba
    auto factory = makeEBFabFactory(cgeom,crse_ba,P_new.DistributionMap(),{1,1,1},EBSupport::basic);
    MultiFab     crse_phi(crse_ba,P_new.DistributionMap(),1,0,MFInfo(),*factory);

#else
    MultiFab     crse_phi(crse_ba,P_new.DistributionMap(),1,0);
#endif

    crse_phi.setVal(1.e200);
    crse_phi.copy(phi,0,0,1);

#ifdef AMREX_USE_EB
    // FIXME -
    // For now, just zero covered fine cells. Better interpolation to come...
    ///
    EB_set_covered(crse_phi,0.);
#endif

    NavierStokesBase& fine_lev        = getLevel(f_lev);
    const Real    cur_fine_pres_time  = fine_lev.state[Press_Type].curTime();
    const Real    prev_fine_pres_time = fine_lev.state[Press_Type].prevTime();

    if (state[Press_Type].descriptor()->timeType() ==
	StateDescriptor::Point && first_crse_step_after_initial_iters)
    {
        const Real time_since_zero  = cur_crse_pres_time - prev_crse_pres_time;
        const Real dt_to_prev_time  = prev_fine_pres_time - prev_crse_pres_time;
        const Real dt_to_cur_time   = cur_fine_pres_time - prev_crse_pres_time;
        const Real cur_mult_factor  = dt_to_cur_time / time_since_zero;
        const Real prev_mult_factor = dt_to_prev_time / dt_to_cur_time;

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
	  FArrayBox     fine_phi;

	  for (MFIter mfi(P_new,true); mfi.isValid(); ++mfi)
          {
	    const Box&  fbx     = mfi.tilebox();

            fine_phi.resize(fbx,1);
            fine_phi.setVal<RunOn::Host>(1.e200);
            node_bilinear_interp.interp(crse_phi[mfi],0,fine_phi,0,1,
                                        fine_phi.box(),ratio,cgeom,fgeom,bc,
                                        0,Press_Type,RunOn::Cpu);

            fine_phi.mult<RunOn::Host>(cur_mult_factor);
            P_new[mfi].plus<RunOn::Host>(fine_phi,fbx,0,0);
            fine_phi.mult<RunOn::Host>(prev_mult_factor);
            P_old[mfi].plus<RunOn::Host>(fine_phi,fbx,0,0);
	  }
	}
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
      {
        FArrayBox     fine_phi;

	for (MFIter mfi(P_new,true); mfi.isValid(); ++mfi)
        {
	    const Box&  fbx     = mfi.tilebox();

            fine_phi.resize(fbx,1);
            fine_phi.setVal<RunOn::Host>(1.e200);
            node_bilinear_interp.interp(crse_phi[mfi],0,fine_phi,0,1,
                                        fine_phi.box(),ratio,cgeom,fgeom,bc,
                                        0,Press_Type,RunOn::Cpu);

            P_new[mfi].plus<RunOn::Host>(fine_phi,fbx,0,0);
            P_old[mfi].plus<RunOn::Host>(fine_phi,fbx,0,0);
        }
      }
    }
#ifdef AMREX_USE_EB
    // FIXME? - this can probably go after new interpolation is implemented
    EB_set_covered(P_new,0.);
    EB_set_covered(P_old,0.);
#endif

}

std::string
NavierStokesBase::thePlotFileType () const
{
    //
    // Increment this whenever the writePlotFile() format changes.
    //
    static const std::string the_plot_file_type("NavierStokes-V1.1");

    return the_plot_file_type;
}

//
// This routine advects the velocities
//
void
NavierStokesBase::velocity_advection (Real dt)
{
    BL_PROFILE("NavierStokesBase::velocity_advection()");

    if (verbose)
    {
        if (do_mom_diff == 0)
        {
            amrex::Print() << "... advect velocities\n";
        }
        else
        {
            if (predict_mom_together == 0)
            {
                amrex::Print() << "Must set predict_mom_together == 1 in NavierStokesBase." << '\n';
                exit(0);
            }
            amrex::Print() << "... advect momenta\n";
        }
    }

    const int   finest_level   = parent->finestLevel();
    const Real  prev_time      = state[State_Type].prevTime();

    //
    // Compute viscosity components.
    //
#ifdef AMREX_USE_EB
    MultiFab& Gp = getGradP();
    Gp.FillBoundary(geom.periodicity());
#else
    MultiFab Gp(grids,dmap,BL_SPACEDIM,1);
    getGradP(Gp, state[Press_Type].prevTime());
#endif

    MultiFab visc_terms(grids,dmap,AMREX_SPACEDIM,1,MFInfo(),Factory());

    // No need to compute this is we are using EB because we will
    // not use Godunov.
#ifndef AMREX_USE_EB
    if (be_cn_theta != 1.0)
    {
        getViscTerms(visc_terms,Xvel,AMREX_SPACEDIM,prev_time);
    }
    else
    {
        visc_terms.setVal(0.,1);
    }
#endif

    MultiFab divu_fp(grids,dmap,1,1,MFInfo(),Factory());
    create_mac_rhs(divu_fp,1,prev_time,dt);

    MultiFab fluxes[BL_SPACEDIM];

    if (do_reflux)
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        {
            const BoxArray& ba = getEdgeBoxArray(i);
            fluxes[i].define(ba, dmap, AMREX_SPACEDIM, 0, MFInfo(),Factory());
        }
    }

    //
    // Compute the advective forcing.
    //
    {
        FillPatchIterator U_fpi(*this,visc_terms,Godunov::hypgrow(),prev_time,State_Type,Xvel,AMREX_SPACEDIM);
        MultiFab& Umf=U_fpi.get_mf();

#ifndef AMREX_USE_EB
        //
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>  NON-EB ALGORITHM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        //
        FillPatchIterator Rho_fpi(*this,visc_terms,Godunov::hypgrow(),prev_time,State_Type,Density,1);
        FillPatchIterator S_fpi(*this,visc_terms,1,prev_time,State_Type,Density,NUM_SCALARS);

        MultiFab& Rmf=Rho_fpi.get_mf();
        MultiFab& Smf=S_fpi.get_mf();

        const Real* dx = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            Vector<int> bndry[AMREX_SPACEDIM];
            FArrayBox tforces;
            FArrayBox S;
            FArrayBox cfluxes[AMREX_SPACEDIM];
            FArrayBox edgstate[AMREX_SPACEDIM];

            for (MFIter U_mfi(Umf,true); U_mfi.isValid(); ++U_mfi)
            {
                const Box& bx=U_mfi.tilebox();

                if (getForceVerbose)
                {
                    amrex::Print() << "---" << '\n'
                                   << "B - velocity advection:" << '\n'
                                   << "Calling getForce..." << '\n';
                }

                getForce(tforces,bx,1,Xvel,BL_SPACEDIM,prev_time,Umf[U_mfi],Smf[U_mfi],0);

                godunov->Sum_tf_gp_visc(tforces,visc_terms[U_mfi],Gp[U_mfi],rho_ptime[U_mfi]);

                D_TERM(bndry[0] = fetchBCArray(State_Type,bx,0,1);,
                       bndry[1] = fetchBCArray(State_Type,bx,1,1);,
                       bndry[2] = fetchBCArray(State_Type,bx,2,1););

                for (int d=0; d<AMREX_SPACEDIM; ++d)
                {
                    const Box& ebx = amrex::surroundingNodes(bx,d);
                    cfluxes[d].resize(ebx,BL_SPACEDIM+1);
                    edgstate[d].resize(ebx,BL_SPACEDIM+1);
                }

                //
                // Loop over the velocity components.
                //
                S.resize(grow(bx,Godunov::hypgrow()),BL_SPACEDIM);
                S.copy<RunOn::Host>(Umf[U_mfi],0,0,BL_SPACEDIM);

                FArrayBox& divufab = divu_fp[U_mfi];
                FArrayBox& aofsfab = (*aofs)[U_mfi];

                D_TERM(FArrayBox& u_mac_fab0 = u_mac[0][U_mfi];,
                       FArrayBox& u_mac_fab1 = u_mac[1][U_mfi];,
                       FArrayBox& u_mac_fab2 = u_mac[2][U_mfi];);

                for (int comp = 0 ; comp < BL_SPACEDIM ; comp++ )
                {
                    int use_conserv_diff = (advectionType[comp] == Conservative) ? true : false;

                    if (do_mom_diff == 1)
                    {
                        S.mult<RunOn::Host>(Rmf[U_mfi],S.box(),S.box(),0,comp,1);
                        tforces.mult<RunOn::Host>(rho_ptime[U_mfi],tforces.box(),tforces.box(),0,comp,1);
                    }

                    // WARNING: FPU argument is not used because FPU is by default in AdvectState
                    godunov->AdvectState(bx, dx, dt,
                                         area[0][U_mfi], u_mac_fab0, cfluxes[0],
                                         area[1][U_mfi], u_mac_fab1, cfluxes[1],
#if (AMREX_SPACEDIM == 3)
                                         area[2][U_mfi], u_mac_fab2, cfluxes[2],
#endif
                                         Umf[U_mfi], S, tforces, divufab, comp,
                                         aofsfab,comp,use_conserv_diff,
                                         comp,bndry[comp].dataPtr(),FPU,volume[U_mfi]);

                    if (do_reflux)
                    {
                        for (int d = 0; d < AMREX_SPACEDIM; d++)
                        {
                            const Box& ebx = U_mfi.nodaltilebox(d);
                            fluxes[d][U_mfi].copy<RunOn::Host>(cfluxes[d],ebx,0,ebx,comp,1);
                        }
                    }
                }
            } // end of MFIter
        } // end OMP region
#else
        //
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>  EB ALGORITHM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        //

         MultiFab cfluxes[AMREX_SPACEDIM];
         MultiFab edgstate[AMREX_SPACEDIM];
         int nghost(2);

         for (int i(0); i < AMREX_SPACEDIM; i++)
         {
             const BoxArray& ba = getEdgeBoxArray(i);
             cfluxes[i].define(ba, dmap, AMREX_SPACEDIM, nghost, MFInfo(), Umf.Factory());
             edgstate[i].define(ba, dmap, AMREX_SPACEDIM, nghost, MFInfo(), Umf.Factory());
         }

         Vector<BCRec> math_bcs(AMREX_SPACEDIM);
         math_bcs = fetchBCArray(State_Type, Xvel, AMREX_SPACEDIM);

         MOL::ComputeAofs(*aofs, Xvel, AMREX_SPACEDIM, Umf, 0,
                          D_DECL(u_mac[0],u_mac[1],u_mac[2]),
                          D_DECL(edgstate[0],edgstate[1],edgstate[2]), 0, false,
                          D_DECL(cfluxes[0],cfluxes[1],cfluxes[2]), 0,
                          math_bcs, geom  );

	 // don't think this is needed here any more. Godunov sets covered vals now...
         EB_set_covered(*aofs, 0.);

         if (do_reflux)
         {
             for (int d(0); d < AMREX_SPACEDIM; d++)
                 MultiFab::Copy(fluxes[d], cfluxes[d], 0, 0, AMREX_SPACEDIM, 0 );

         }

#endif
    } //end scope of FillPatchIter

    if (do_reflux)
    {
        if (level > 0 )
	{
            for (int d = 0; d < BL_SPACEDIM; d++)
                advflux_reg->FineAdd(fluxes[d],d,0,0,BL_SPACEDIM,dt);
	}
        if(level < finest_level)
	{
            for (int i = 0; i < BL_SPACEDIM; i++)
                getAdvFluxReg(level+1).CrseInit(fluxes[i],i,0,0,BL_SPACEDIM,-dt);
	}
    }
}

//
// This subroutine updates the velocity field before the level projection.
//
// At this point in time, all we know is u^n, rho^n+1/2, and the
// general forcing terms at t^n, and after solving in this routine
// viscous forcing at t^n+1/2.  Except for a simple buoyancy term,
// b = -rho^n+1/2 g, it is usually not possible to estimate more
// general forcing terms at t^n+1/2.  Since the default getForce, handles
// this case automatically, F_new and F_old have been replaced by a single
// tforces FArrayBox.
//
// We assume that if one component of velocity is viscous that all must be.
//

void
NavierStokesBase::velocity_update (Real dt)
{
    BL_PROFILE("NavierStokesBase::velocity_update()");

    if (verbose)
    {
      if (do_mom_diff == 0)
      {
	amrex::Print() << "... update velocities \n";
      }
      else
      {
	amrex::Print() << "... update momenta \n";
      }
    }

    velocity_advection_update(dt);

    if (!initial_iter)
        velocity_diffusion_update(dt);
    else
        initial_velocity_diffusion_update(dt);

    MultiFab&  S_new     = get_new_data(State_Type);

    for (int sigma = 0; sigma < BL_SPACEDIM; sigma++)
    {
       if (S_new.contains_nan(sigma,1,0))
       {
	 amrex::Print() << "New velocity " << sigma << " contains Nans" << '\n';
	 exit(0);
       }
    }
}

void
NavierStokesBase::velocity_advection_update (Real dt)
{
    BL_PROFILE("NavierStokesBase::velocity_advection_update()");

    MultiFab&  U_old          = get_old_data(State_Type);
    MultiFab&  U_new          = get_new_data(State_Type);
    MultiFab&  Aofs           = *aofs;
    const Real prev_pres_time = state[Press_Type].prevTime();

#ifdef AMREX_USE_EB
    MultiFab& Gp=*gradp;
    Gp.FillBoundary(geom.periodicity());
    if (do_mom_diff == 1)
      amrex::Abort("NavierStokesBase::velocity_advection_update(): do_mom_diff==1 not currently working with EB.");
    // Changing Gp to a pointer causes problems for non-EB ...
    // Think issue was need to use MultiFab::Copy; plain Copy was going someplace funny
    // MultiFab* Gp;
      // if (do_mom_diff == 1){
      // 	//need to make a copy of gradp
      // 	Gp->define(grids,dmap,BL_SPACEDIM,1);
      // 	Copy(*Gp,*gradp,0,0,3,0);
      // }
      // else{
      // 	Gp = gradp.get();
      // }
#else
      MultiFab Gp(grids,dmap,BL_SPACEDIM,1);
      getGradP(Gp, prev_pres_time);
#endif

    MultiFab& halftime = get_rho_half_time();
#ifdef _OPENMP
#pragma omp parallel
#endif
{

    FArrayBox  tforces, S;
    for (MFIter Rhohalf_mfi(halftime,true); Rhohalf_mfi.isValid(); ++Rhohalf_mfi)
    {
        const int i = Rhohalf_mfi.index();
        const Box& bx = Rhohalf_mfi.tilebox();

        //
        // Need to do some funky half-time stuff.
        //
        if (getForceVerbose)
           amrex::Print() << "---" << '\n' << "F - velocity advection update (half time):" << '\n';
        //
        // Average the mac face velocities to get cell centred velocities.
        //
        FArrayBox Vel(amrex::grow(bx,0),BL_SPACEDIM);
        FORT_AVERAGE_EDGE_STATES(BL_TO_FORTRAN_ANYD(Vel),
                                 BL_TO_FORTRAN_ANYD(u_mac[0][Rhohalf_mfi]),
                                 BL_TO_FORTRAN_ANYD(u_mac[1][Rhohalf_mfi]),
#if (AMREX_SPACEDIM==3)
                                 BL_TO_FORTRAN_ANYD(u_mac[2][Rhohalf_mfi]),
#endif
                                 &getForceVerbose);
        //
        // Average the new and old time to get Crank-Nicholson half time approximation.
        //
        FArrayBox Scal(amrex::grow(bx,0),NUM_SCALARS);
        Scal.copy<RunOn::Host>(U_old[Rhohalf_mfi],bx,Density,bx,0,NUM_SCALARS);
        Scal.plus<RunOn::Host>(U_new[Rhohalf_mfi],bx,Density,0,NUM_SCALARS);
        Scal.mult<RunOn::Host>(0.5,bx,0,NUM_SCALARS);

        if (getForceVerbose) amrex::Print() << "Calling getForce..." << '\n';
        const Real half_time = 0.5*(state[State_Type].prevTime()+state[State_Type].curTime());
        getForce(tforces,bx,0,Xvel,BL_SPACEDIM,half_time,Vel,Scal,0);

        //
        // Do following only at initial iteration--per JBB.
        //
        if (initial_iter && is_diffusive[Xvel])
            tforces.setVal<RunOn::Host>(0);

	//why is this on a grown box?
	   const Box& sbx = Rhohalf_mfi.growntilebox();
        S.resize(sbx,BL_SPACEDIM);
        S.copy<RunOn::Host>(U_old[Rhohalf_mfi],sbx,0,sbx,0,BL_SPACEDIM);

        if (do_mom_diff == 1)
        {
            for (int d = 0; d < BL_SPACEDIM; d++)
            {
                Gp[Rhohalf_mfi].mult<RunOn::Host>(halftime[i],bx,0,d,1);
                tforces.mult<RunOn::Host>(halftime[i],bx,0,d,1);
                S.mult<RunOn::Host>(rho_ptime[Rhohalf_mfi],bx,0,d,1);
            }
        }

        godunov->Add_aofs_tf_gp(S,U_new[Rhohalf_mfi],Aofs[Rhohalf_mfi],tforces,
                                Gp[Rhohalf_mfi],halftime[i],bx,dt);
        if (do_mom_diff == 1)
        {
            for (int d = 0; d < BL_SPACEDIM; d++)
                U_new[Rhohalf_mfi].divide<RunOn::Host>(rho_ctime[Rhohalf_mfi],bx,0,d,1);
        }
    }
}

    for (int sigma = 0; sigma < BL_SPACEDIM; sigma++)
    {
       if (U_old.contains_nan(sigma,1,0))
       {
	 amrex::Print() << "Old velocity " << sigma << " contains Nans" << '\n';
       }
       if (U_new.contains_nan(sigma,1,0))
       {
         for(MFIter mfi(U_new); mfi.isValid(); ++mfi){
             if(U_new[mfi].contains_nan<RunOn::Host>(mfi.validbox(),sigma,1)){
		            FArrayBox tmp(mfi.validbox(),1);
                tmp.copy<RunOn::Host>(U_new[mfi],mfi.validbox(),sigma,mfi.validbox(),0,1);

		}
         }
	 amrex::Print() << "New velocity " << sigma << " contains Nans" << '\n';
       }
    }
}

void
NavierStokesBase::initial_velocity_diffusion_update (Real dt)
{
    //
    // Do following only at initial iteration.
    //
    if (is_diffusive[Xvel])
    {
        MultiFab&  U_old          = get_old_data(State_Type);
        MultiFab&  U_new          = get_new_data(State_Type);
        MultiFab&  Rh             = get_rho_half_time();
        const Real prev_time      = state[State_Type].prevTime();
        const int  xvel           = Xvel;

        int              ng(1);
	MultiFab visc_terms(grids,dmap,AMREX_SPACEDIM,ng,MFInfo(),Factory());
        MultiFab    tforces(grids,dmap,AMREX_SPACEDIM,ng,MFInfo(),Factory());

        //
        // Get grad(p)
        //
#ifdef AMREX_USE_EB
        MultiFab& Gp = *gradp;
#else
        const Real prev_pres_time = state[Press_Type].prevTime();
        MultiFab         Gp(grids,dmap,AMREX_SPACEDIM,ng,MFInfo(),Factory());
        getGradP(Gp, prev_pres_time);
#endif

        //
        // Compute additional forcing terms
        //
        tforces.setVal(0.0,1);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            FArrayBox tforces_fab;
            for (MFIter mfi(tforces,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const auto& bx = mfi.tilebox();

                if (getForceVerbose)
                {
                    amrex::Print() << "---" << '\n'
                                   << "G - initial velocity diffusion update:" << '\n'
                                   << "Calling getForce..." << '\n';
                }

                getForce(tforces_fab,bx,0,Xvel,AMREX_SPACEDIM,prev_time,U_old[mfi],U_old[mfi],Density);
                tforces[mfi].copy<RunOn::Host>(tforces_fab,bx,0,bx,0,AMREX_SPACEDIM);
            }
        }

        //
        // Compute viscous terms
        //
	if (be_cn_theta != 1.0)
        {
	    getViscTerms(visc_terms,Xvel,AMREX_SPACEDIM,prev_time);
        }
        else
	{
	    visc_terms.setVal(0);
	}

        //
        // Assemble RHS
        //

        // Set tforces = tforces + visc - Gp
        MultiFab::Add(tforces, visc_terms, 0, 0, AMREX_SPACEDIM, 0);
        MultiFab::Subtract(tforces, Gp, 0, 0, AMREX_SPACEDIM, 0);

        // If NOT solving for momentum, divide tforces for Rho
        if (!(do_mom_diff==1))
        {
            for (int idim(0); idim < AMREX_SPACEDIM; ++idim )
                MultiFab::Divide(tforces, Rh, 0, idim, 1, 0);
        }

        // Subtract convective term from tforces
        // WARNING: I think it is assumed that the convective term
        // may or may not account for rho depending on do_mom_diff
        // so no need to account for rho here
        MultiFab::Subtract(tforces, *aofs, 0, xvel, AMREX_SPACEDIM, 0);

        if (do_mom_diff==1)
        {
            // If solving for momentum, set U_new=dt*tforces
            // Then set u_new += rho_old*u_old = rho_old*u_old + dt*tforce
            // Finally divide u_new by rho_new to recover velocities
            MultiFab::Copy(U_new, tforces, 0, Xvel, AMREX_SPACEDIM, 0);
            U_new.mult(dt);
            for (int idim(0); idim < AMREX_SPACEDIM; ++idim )
            {
                MultiFab::AddProduct(U_new, rho_ptime, 0, U_old, xvel+idim, xvel+idim, 1, 0);
                MultiFab::Divide(U_new, rho_ctime, 0, xvel+idim, 1, 0);
            }
        }
        else
        {
            // Then set u_new = u_old + dt*tforces
            MultiFab::LinComb(U_new, 1.0, U_old, xvel, dt, tforces, 0, xvel, AMREX_SPACEDIM, 0);
        }
    }

}

Real
NavierStokesBase::volWgtSum (const std::string& name,
			     Real               time)
{
    Real        sum = 0.0;
    const Real* dx  = geom.CellSize();
    auto        mf  = derive(name,time,0);
    BoxArray    baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    std::vector< std::pair<int,Box> > isects;

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            baf.intersections(grids[mfi.index()],isects);

            for (int ii = 0, N = isects.size(); ii < N; ii++)
                fab.setVal<RunOn::Host>(0,isects[ii].second,0,fab.nComp());
        }
        Real        s;
        const Real* dat = fab.dataPtr();
        const int*  dlo = fab.loVect();
        const int*  dhi = fab.hiVect();
        const Box&  grdbx = grids[mfi.index()];
        const int*  lo  = grdbx.loVect();
        const int*  hi  = grdbx.hiVect();
#ifdef AMREX_USE_EB
        if ( (volWgtSum_sub_dz > 0 && volWgtSum_sub_Rcyl > 0) || Geom().IsRZ() )
	  amrex::Abort("EB volWgtSum currently only works over entire cartesian domain.");

	const Real* vf   = (*volfrac)[mfi].dataPtr();
        const int*  vflo = (*volfrac)[mfi].loVect();
        const int*  vfhi = (*volfrac)[mfi].hiVect();
#endif

#if (BL_SPACEDIM == 2)
        int   rz_flag = Geom().IsRZ() ? 1 : 0;
        Real* rad     = &radius[mfi.index()][0];
        int   irlo    = lo[0]-radius_grow;
        int   irhi    = hi[0]+radius_grow;
        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
        if (volWgtSum_sub_dz > 0 && volWgtSum_sub_Rcyl > 0)
        {
            const Real* plo = geom.ProbLo();
            summass_cyl(dat,ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),
                             dx,&s,rad,&irlo,&irhi,&rz_flag,plo,
                             &volWgtSum_sub_dz,&volWgtSum_sub_Rcyl);
        }
        else
        {
#ifdef AMREX_USE_EB
	    summass_eb(dat,ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),
		       vf,ARLIM(vflo),ARLIM(vfhi),dx,&s,rad,&irlo,&irhi,&rz_flag);
#else
	    summass(dat,ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),
		    dx,&s,rad,&irlo,&irhi,&rz_flag);
#endif
	}
#endif

#if (BL_SPACEDIM == 3)
        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
        if (volWgtSum_sub_dz > 0 && volWgtSum_sub_Rcyl > 0)
        {
            const Real* plo = geom.ProbLo();
            summass_cyl(dat,ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),
                             dx,plo,&volWgtSum_sub_dz,&volWgtSum_sub_Rcyl,&s);
        }
        else
        {
#ifdef AMREX_USE_EB
	    summass_eb(dat,ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),
		       vf,ARLIM(vflo),ARLIM(vfhi),dx,&s);
#else
	    summass(dat,ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),dx,&s);
#endif
        }
#endif
        sum += s;
    }

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

#if (BL_SPACEDIM == 3)
void
NavierStokesBase::sum_turbulent_quantities ()
{
    Real time = state[State_Type].curTime();
    const int finestLevel = parent->finestLevel();
    const Real *dx = parent->Geom(finestLevel).CellSize();
    const int ksize(parent->Geom(finestLevel).Domain().length(2));
    const int turbVars(33);
    int refRatio(1);

    Real* turb = new Real[turbVars*ksize];

    for (int i=0; i<turbVars*ksize; i++) turb[i]=0;

    for (int lev = finestLevel; lev >= 0; lev--)
    {
	const int levKsize(parent->Geom(lev).Domain().length(2));

	Real* levTurb = new Real[turbVars*levKsize];

	for (int i=0; i<turbVars*levKsize; i++) levTurb[i]=0;

        NavierStokesBase& ns_level = getLevel(lev);
	ns_level.TurbSum(time,levTurb,levKsize,turbVars);

	if (lev<finestLevel)  refRatio *= parent->refRatio(lev)[2];
	else                  refRatio  = 1;

	for (int l=0, k=0; l<levKsize; l++)
	    for (int r=0; r<refRatio; r++, k++)
		for (int v=0; v<turbVars; v++)
		    turb[k*turbVars+v] += levTurb[l*turbVars+v];

	delete [] levTurb;
    }

    ParallelDescriptor::ReduceRealSum(&turb[0], ksize*turbVars, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
        std::string DirPath = "TurbData";
        if (!amrex::UtilCreateDirectory(DirPath, 0755))
            amrex::CreateDirectoryFailed(DirPath);

        const int steps = parent->levelSteps(0);
        FILE *file;

        std::string filename = amrex::Concatenate("TurbData/TurbData_", steps, 4);
        filename += ".dat";

        file = fopen(filename.c_str(),"w");
        for (int k=0; k<ksize; k++)
        {
            fprintf(file,"%e ",dx[2]*(0.5+(double)k));
            for (int v=0; v<turbVars; v++)
                fprintf(file,"%e ",turb[k*turbVars+v]);
            fprintf(file,"\n");
        }
        fclose(file);
    }

    delete [] turb;
}

void
NavierStokesBase::TurbSum (Real time, Real *turb, int ksize, int turbVars)
{
    const Real* dx = geom.CellSize();

    const int turbGrow(0);
    const int presGrow(0);
    auto turbMF = derive("TurbVars",time,turbGrow);
    auto presMF = derive("PresVars",time,presGrow);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    std::vector< std::pair<int,Box> > isects;

    for (MFIter turbMfi(*turbMF), presMfi(*presMF);
	 turbMfi.isValid() && presMfi.isValid();
	 ++turbMfi, ++presMfi)
    {
	FArrayBox& turbFab = (*turbMF)[turbMfi];
	FArrayBox& presFab = (*presMF)[presMfi];

        if (level < parent->finestLevel())
        {
            baf.intersections(grids[turbMfi.index()],isects);

            for (int ii = 0, N = isects.size(); ii < N; ii++)
            {
                presFab.setVal<RunOn::Host>(0,isects[ii].second,0,presMF->nComp());
                turbFab.setVal<RunOn::Host>(0,isects[ii].second,0,turbMF->nComp());
            }
        }
    }

    turbMF->FillBoundary(0,turbMF->nComp(), geom.periodicity());
    presMF->FillBoundary(0,presMF->nComp(), geom.periodicity());

    for (MFIter turbMfi(*turbMF), presMfi(*presMF);
	 turbMfi.isValid() && presMfi.isValid();
	 ++turbMfi, ++presMfi)
    {
	FArrayBox& turbFab = (*turbMF)[turbMfi];
	FArrayBox& presFab = (*presMF)[presMfi];

        const Real* turbData = turbFab.dataPtr();
        const Real* presData = presFab.dataPtr();
        const int*  dlo = turbFab.loVect();
        const int*  dhi = turbFab.hiVect();
        const int*  plo = presFab.loVect();
        const int*  phi = presFab.hiVect();
	const Box& grdbx = grids[turbMfi.index()];
        const int*  lo  = grdbx.loVect();
        const int*  hi  = grdbx.hiVect();

        sumturb(turbData,presData,ARLIM(dlo),ARLIM(dhi),ARLIM(plo),ARLIM(phi),ARLIM(lo),ARLIM(hi),
		     dx,turb,&ksize,&turbVars);
   }
}

#ifdef SUMJET
void
NavierStokesBase::JetSum (Real time, Real *jetData, int levRsize,  int levKsize,  int rsize,  int ksize, int jetVars)
{
    const Real* dx = geom.CellSize();

    const int turbGrow(0);
    const int presGrow(0);

    auto turbMF = derive("JetVars",time,turbGrow);
    auto presMF = derive("JetPresVars",time,presGrow);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    std::vector< std::pair<int,Box> > isects;

    for (MFIter turbMfi(*turbMF), presMfi(*presMF);
	 turbMfi.isValid() && presMfi.isValid();
	 ++turbMfi, ++presMfi)
    {
	FArrayBox& turbFab = (*turbMF)[turbMfi];
	FArrayBox& presFab = (*presMF)[presMfi];

        if (level < parent->finestLevel())
        {
            baf.intersections(grids[turbMfi.index()],isects);

            for (int ii = 0, N = isects.size(); ii < N; ii++)
            {
                presFab.setVal(0,isects[ii].second,0,presMF->nComp());
                turbFab.setVal(0,isects[ii].second,0,turbMF->nComp());
            }
        }
    }

    turbMF->FillBoundary(0,turbMF->nComp(), geom.periodicity());
    presMF->FillBoundary(0,presMF->nComp(), geom.periodicity());

    for (MFIter turbMfi(*turbMF), presMfi(*presMF);
	 turbMfi.isValid() && presMfi.isValid();
	 ++turbMfi, ++presMfi)
    {
	FArrayBox& turbFab = (*turbMF)[turbMfi];
	FArrayBox& presFab = (*presMF)[presMfi];

        RealBox     gridloc  = RealBox(grids[turbMfi.index()],geom.CellSize(),geom.ProbLo());
        const Real* turbData = turbFab.dataPtr();
        const Real* presData = presFab.dataPtr();
        const int*  dlo = turbFab.loVect();
        const int*  dhi = turbFab.hiVect();
        const int*  plo = presFab.loVect();
        const int*  phi = presFab.hiVect();
        const int*  lo  = grids[turbMfi.index()].loVect();
        const int*  hi  = grids[turbMfi.index()].hiVect();

        sumjet(turbData,presData,ARLIM(dlo),ARLIM(dhi),ARLIM(plo),ARLIM(phi),ARLIM(lo),ARLIM(hi),
		    dx,jetData,&levRsize,&levKsize,&rsize,&ksize,&jetVars,&jet_interval_split,
		    gridloc.lo(),gridloc.hi());
    }
}

void
NavierStokesBase::sum_jet_quantities ()
{
    Real time = state[State_Type].curTime();
    const int finestLevel = parent->finestLevel();
    const Real *dx = parent->Geom(finestLevel).CellSize();
    const int isize(parent->Geom(finestLevel).Domain().length(0));
    const int ksize(parent->Geom(finestLevel).Domain().length(2));
    const int rsize=isize>>1;
    const int jetVars(104);

    amrex::Print() << "NavierStokesBase::sum_jet_quantities():" << '\n'
		   << "   jetVars: " << jetVars << '\n'
		   << "   rsize  : " << rsize << '\n'
		   << "   ksize  : " << ksize << '\n';

    Real* jetData = new Real[jetVars*ksize*rsize];

    for (int i=0; i<jetVars*ksize*rsize; i++) jetData[i]=0;

    for (int lev = finestLevel; lev >= 0; lev--)
    {
	const int levIsize(parent->Geom(lev).Domain().length(0));
	const int levKsize(parent->Geom(lev).Domain().length(2));
	const int levRsize(levIsize>>1);

        NavierStokesBase& ns_level = getLevel(lev);
	ns_level.JetSum(time,jetData,levRsize,levKsize,rsize,ksize,jetVars);
    }

    ParallelDescriptor::ReduceRealSum(&jetData[0], ksize*rsize*jetVars, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "      Creating JetData..." << '\n';
        std::string DirPath = "JetData";
        if (!amrex::UtilCreateDirectory(DirPath, 0755))
            amrex::CreateDirectoryFailed(DirPath);

        const int steps = parent->levelSteps(0);
        FILE *file;
        std::string filename;

	Vector<Real> r(rsize);
	for (int i=0; i<rsize; i++)
	    r[i] = dx[0]*(0.5+(double)i);
	Vector<Real> z(ksize);
	for (int k=0; k<ksize; k++)
	    z[k] = dx[2]*(0.5+(double)k);

#if 0
        filename  = amrex::Concatenate("JetData/JetData_", steps, 4);
        filename += "_r.dat";

	file = fopen(filename.c_str(),"w");
	for (int i=0; i<rsize; i++)
	    fprintf(file,"%e ",r[i]);
	fclose(file);

        filename  = amrex::Concatenate("JetData/JetData_", steps, 4);
        filename += "_z.dat";

	file = fopen(filename.c_str(),"w");
	for (int k=0; k<ksize; k++)
	    fprintf(file,"%e ",dx[2]*(0.5+(double)k));
	fclose(file);

	for (int v=0; v<jetVars; v++) {

            filename  = amrex::Concatenate("JetData/JetData_", steps, 4);
            filename += amrex::Concatenate(filename + "_v", v, 4);
            filename += ".dat";

	    file = fopen(filename.c_str(),"w");
	    for (int k=0; k<ksize; k++) {
		for (int i=0; i<rsize; i++) {
		    fprintf(file,"%e ",jetData[(k*rsize+i)*jetVars+v]);
		}
		fprintf(file,"\n");
	    }
	    fclose(file);
	    amrex::Print() << "   ...done." << '\n';
	}
#else
	std::string FullPath = amrex::Concatenate("JetData/JD", steps, 4);

	if (!amrex::UtilCreateDirectory(FullPath, 0755))
	    amrex::CreateDirectoryFailed(FullPath);

        filename = FullPath;
        filename += '/';
        filename += "data.bin";

	file=fopen(filename.c_str(),"w");
	fwrite(&time,sizeof(double),1,file);
	fwrite(&rsize,sizeof(int),1,file);
	fwrite(&ksize,sizeof(int),1,file);
	fwrite(&jetVars,sizeof(int),1,file);
	fwrite(r.dataPtr(),sizeof(Real),rsize,file);
	fwrite(z.dataPtr(),sizeof(Real),ksize,file);
	fwrite(jetData,sizeof(Real),jetVars*rsize*ksize,file);
	fclose(file);
#endif
    }

    delete [] jetData;
}
#endif // SUMJET

#endif  // (BL_SPACEDIM == 3)

#ifdef AMREX_PARTICLES

void
NavierStokesBase::read_particle_params ()
{
    ParmParse ppp("particles");
    //
    // The directory in which to store timestamp files.
    //
    ppp.query("timestamp_dir", timestamp_dir);
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!amrex::UtilCreateDirectory(timestamp_dir, 0755))
            amrex::CreateDirectoryFailed(timestamp_dir);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (int nc = ppp.countval("timestamp_indices"))
    {
        timestamp_indices.resize(nc);

        ppp.getarr("timestamp_indices", timestamp_indices, 0, nc);
    }

    ppp.query("pverbose",pverbose);
    //
    // Used in initData() on startup to read in a file of particles.
    //
    ppp.query("particle_init_file", particle_init_file);
    //
    // Used in post_restart() to read in a file of particles.
    //
    ppp.query("particle_restart_file", particle_restart_file);
    //
    // This must be true the first time you try to restart from a checkpoint
    // that was written with USE_PARTICLES=FALSE; i.e. one that doesn't have
    // the particle checkpoint stuff (even if there are no active particles).
    // Otherwise the code will fail when trying to read the checkpointed particles.
    //
    ppp.query("restart_from_nonparticle_chkfile", restart_from_nonparticle_chkfile);
    //
    // Used in post_restart() to write out the file of particles.
    //
    ppp.query("particle_output_file", particle_output_file);
}

void
NavierStokesBase::initParticleData ()
{
    if (level == 0)
    {
        if (NSPC == 0)
        {
            NSPC = new AmrTracerParticleContainer(parent);
        }

        NSPC->SetVerbose(pverbose);

        if (!particle_init_file.empty())
        {
            NSPC->InitFromAsciiFile(particle_init_file,0);
        }
    }
}

void
NavierStokesBase::post_restart_particle ()
{
    if (level == 0)
    {
        BL_ASSERT(NSPC == 0);

        NSPC = new AmrTracerParticleContainer(parent);

        NSPC->SetVerbose(pverbose);
        //
        // We want to be able to add new particles on a restart.
        // As well as the ability to write the particles out to an ascii file.
        //
        if (!restart_from_nonparticle_chkfile)
        {
            NSPC->Restart(parent->theRestartFile(), the_ns_particle_file_name);
        }

        if (!particle_restart_file.empty())
        {
            NSPC->InitFromAsciiFile(particle_restart_file,0);
        }

        if (!particle_output_file.empty())
        {
            NSPC->WriteAsciiFile(particle_output_file);
        }
    }
}

void
NavierStokesBase::post_timestep_particle (int crse_iteration)
{
    const int ncycle = parent->nCycle(level);
    const int finest_level = parent->finestLevel();
    //
    // Don't redistribute/timestamp on the final subiteration except on the coarsest grid.
    //
    if (NSPC != 0 && (crse_iteration < ncycle || level == 0))
    {
        const Real curr_time = state[State_Type].curTime();

	int ngrow = (level == 0) ? 0 : crse_iteration;

        NSPC->Redistribute(level, finest_level, ngrow);

        if (!timestamp_dir.empty())
        {
            std::string basename = timestamp_dir;

            if (basename[basename.length()-1] != '/') basename += '/';

            basename += "Timestamp";

	    static bool first = true;
	    static int n, nextras;
	    static std::vector<int> tindices;

	    if (first)
	    {
		first = false;

		n = timestamp_indices.size();
		nextras = timestamp_num_extras();

		int sz = n + nextras;
		tindices.reserve(sz);

		for (int i = 0; i < sz; ++i) {
		    tindices.push_back(i);
		}
	    }

            for (int lev = level; lev <= finest_level; lev++)
            {
                if (NSPC->NumberOfParticlesAtLevel(lev) <= 0) continue;

		int ng = (lev == level) ? ngrow+1 : 1;

		AmrLevel& amr_level = parent->getLevel(lev);
		MultiFab& S_new = amr_level.get_new_data(State_Type);

		MultiFab tmf;

		if (tindices.size() > 0)
		{
		    tmf.define(S_new.boxArray(), S_new.DistributionMap(), tindices.size(), ng, MFInfo(), Factory());

		    if (n > 0)
		    {
		      FillPatchIterator fpi(parent->getLevel(lev), S_new,
                                            ng, curr_time, State_Type, 0, NUM_STATE);
                      const MultiFab& S = fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel
#endif
                      for (MFIter mfi(tmf,true); mfi.isValid(); ++mfi)
                      {
                        FArrayBox& tfab = tmf[mfi];
                        const FArrayBox& sfab = S[mfi];
                        const Box& box = mfi.growntilebox();
                        for (int i = 0; i < n; ++i)
                        {
                          tfab.copy(sfab, box, timestamp_indices[i], box, i, 1);
                        }
                      }
		    }

		    if (nextras > 0)
		    {
			timestamp_add_extras(lev, curr_time, tmf);
		    }
		}

		NSPC->Timestamp(basename, tmf, lev, curr_time, tindices);
            }
        }
    }
}

std::unique_ptr<MultiFab>
NavierStokesBase::ParticleDerive (const std::string& name,
				  Real               time,
				  int                ngrow)
{
    if (name == "particle_count" || name == "total_particle_count") {
	int ncomp = 1;
	const DeriveRec* rec = derive_lst.get(name);
	if (rec)
	{
	    ncomp = rec->numDerive();
	}

        MultiFab* ret = new MultiFab(grids, dmap, ncomp, ngrow, MFInfo(), Factory());
	ParticleDerive(name,time,*ret,0);
	return std::unique_ptr<MultiFab>{ret};
    }
    else {
	return AmrLevel::derive(name, time, ngrow);
    }
}

void
NavierStokesBase::ParticleDerive (const std::string& name,
				  Real               time,
				  MultiFab&          mf,
				  int                dcomp)
{
    if (NSPC == 0 || !(name == "particle_count" || name == "total_particle_count"))
    {
        AmrLevel::derive(name,time,mf,dcomp);
    }
    else {
	if (name == "particle_count")
	{
	    MultiFab temp_dat(grids,dmap,1,0, MFInfo(), Factory());
	    temp_dat.setVal(0);
	    NSPC->Increment(temp_dat,level);
	    MultiFab::Copy(mf,temp_dat,0,dcomp,1,0);
	}
	else if (name == "total_particle_count")
	{
	    //
	    // We want the total particle count at this level or higher.
	    //
	    ParticleDerive("particle_count",time,mf,dcomp);

	    IntVect trr(D_DECL(1,1,1));

	    for (int lev = level+1; lev <= parent->finestLevel(); lev++)
	    {
		BoxArray ba = parent->boxArray(lev);

		MultiFab temp_dat(ba,parent->DistributionMap(lev),1,0,MFInfo(),Factory());

		trr *= parent->refRatio(lev-1);

		ba.coarsen(trr);

                // FIXME? Won't work because ba has been coarsened. But don't actually need
                //  facotry here anyway...
		//MultiFab ctemp_dat(ba,parent->DistributionMap(lev),1,0,MFInfo(),Factory());
		MultiFab ctemp_dat(ba,parent->DistributionMap(lev),1,0);

		temp_dat.setVal(0);
		ctemp_dat.setVal(0);

		NSPC->Increment(temp_dat,lev);

#ifdef _OPENMP
#pragma omp parallel
#endif
		for (MFIter mfi(temp_dat,true); mfi.isValid(); ++mfi)
		{
		    const FArrayBox& ffab =  temp_dat[mfi];
		    FArrayBox&       cfab = ctemp_dat[mfi];
		    const Box&       fbx  = mfi.tilebox();

		    BL_ASSERT(cfab.box() == amrex::coarsen(fbx,trr));

		    for (IntVect p = fbx.smallEnd(); p <= fbx.bigEnd(); fbx.next(p))
		    {
		        const Real val = ffab(p);
			if (val > 0)
			    cfab(amrex::coarsen(p,trr)) += val;
		    }
		}

		temp_dat.clear();

		MultiFab dat(grids,dmap,1,0,MFInfo(),Factory());
		dat.setVal(0);
		dat.copy(ctemp_dat);

		MultiFab::Add(mf,dat,0,dcomp,1,0);
	    }
	}
	else
	{
	    amrex::Abort("NavierStokesBase::ParticleDerive: how did this happen?");
	}
    }
}

#endif  // AMREX_PARTICLES

// Boundary condition access function.
Vector<int>
NavierStokesBase::fetchBCArray (int State_Type, const Box& bx, int scomp, int ncomp)
{
    Vector<int> bc(2*BL_SPACEDIM*ncomp);
    BCRec bcr;
    const StateDescriptor* stDesc;
    const Box& domain = geom.Domain();

    for (int n = 0; n < ncomp; n++)
    {
      stDesc=state[State_Type].descriptor();
      setBC(bx,domain,stDesc->getBC(scomp+n),bcr);

      const int* b_rec = bcr.vect();
      for (int m = 0; m < 2*BL_SPACEDIM; m++)
	bc[2*BL_SPACEDIM*n + m] = b_rec[m];
    }

    return bc;
}

Vector<BCRec>
NavierStokesBase::fetchBCArray (int State_Type, int scomp, int ncomp)
{
    Vector<BCRec> bc(ncomp);
    const StateDescriptor* stDesc;
    const Box& domain = geom.Domain();

    for (int n(0); n < ncomp; ++n)
    {
      stDesc=state[State_Type].descriptor();
      setBC(domain,domain,stDesc->getBC(scomp+n), bc[n] );
    }

    return bc;
}

void
NavierStokesBase::average_down(const MultiFab& S_fine, MultiFab& S_crse,
			       int scomp, int ncomp)
{
  //
  // Choose the appropriate AMReX average_down() based on
  // whether EB or non-EB, and dimensionality
  //

#ifdef AMREX_USE_EB

  // FIXME?
  // Assume we want EB to behave the same as non-EB in regards to dimensionality
  // Not sure why we'd want 2D to be different than 3D
  // Note that 3D volume weighting doesn't exist for non-EB
#if (AMREX_SPACEDIM == 3)
    // no volume weighting
    amrex::EB_average_down(S_fine, S_crse, scomp, ncomp, fine_ratio);
#else
    // volume weighting
    amrex::EB_average_down(S_fine, S_crse, this->getLevel(level+1).Volume(),
			   *(this->getLevel(level+1).VolFrac()),
			   scomp, ncomp, fine_ratio);
#endif

#else
    // non-EB aware, uses volume weighting for 1D,2D but no volume weighting for 3D
    amrex::average_down(S_fine, S_crse,
			this->getLevel(level+1).geom, this->getLevel(level).geom,
			scomp, ncomp, fine_ratio);
#endif
}


//
//  Diagnostics functions
//
void
NavierStokesBase::printMaxVel (bool new_data)
{

    MultiFab& S = new_data? get_new_data(State_Type) : get_old_data(State_Type);

#if (AMREX_SPACEDIM==3)
    amrex::Print() << "max(abs(u/v/w))  = "
#else
        amrex::Print() << "max(abs(u/v))  = "
#endif
                   << S.norm0( Xvel,   0, false, true )
                   << "  "
		   << S.norm0( Xvel+1, 0, false, true )
#if (AMREX_SPACEDIM==3)
                   << "  "
                   << S.norm0( Xvel+2, 0, false, true )
#endif
                   << std::endl;
}


void
NavierStokesBase::printMaxGp (bool new_data)
{
#ifdef AMREX_USE_EB
    MultiFab& Gp = getGradP();
#else
    MultiFab Gp(grids,dmap,BL_SPACEDIM,1);
    const Real time = new_data ? state[Press_Type].curTime() : state[Press_Type].prevTime();
    getGradP(Gp, time);
#endif
    MultiFab& P  = new_data? get_new_data(Press_Type) : get_old_data(Press_Type);

#if (AMREX_SPACEDIM==3)
    amrex::Print() << "max(abs(gpx/gpy/gpz/p)) = "
#else
    amrex::Print() << "max(abs(gpx/gpy/p)) = "
#endif
                   << Gp.norm0( 0, 0, false, true )
                   << "  "
		   << Gp.norm0( 1, 0, false, true )
#if (AMREX_SPACEDIM==3)
                   << "  "
                   << Gp.norm0( 2, 0, false, true )
#endif
                   << "  "
                   << P.norm0(0, 0, false, true )
                   << std::endl;
}

void
NavierStokesBase::printMaxValues (bool new_data)
{
    printMaxVel(new_data);
    printMaxGp(new_data);
}
