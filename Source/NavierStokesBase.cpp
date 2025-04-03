
#include <AMReX_ParmParse.H>
#include <AMReX_TagBox.H>
#include <AMReX_Utility.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_FillPatchUtil.H>
#include <NavierStokesBase.H>
#include <NSB_K.H>
#include <NS_util.H>
#include <iamr_constants.H>
#include <NS_LS.H>
#include <NS_kernels.H>

#include <hydro_godunov.H>
#include <hydro_bds.H>
#include <hydro_utils.H>

#ifdef AMREX_PARTICLES
#include <DiffusedIB.H>
#endif
#ifdef AMREX_USE_EB
#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBInterpolater.H>
#include <AMReX_EBFArrayBox.H>
#include <hydro_ebgodunov.H>
#include <AMReX_EB_Redistribution.H>
#include <AMReX_EBMultiFabUtil_C.H>
#endif
#ifdef AMREX_USE_TURBULENT_FORCING
#include <TurbulentForcing_params.H>
#endif


using namespace amrex;

//
// Set external dirichlet BC to zero
//
struct HomExtDirFill
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real /*time*/,
                     const BCRec* bcr, const int /*bcomp*/,
                     const int /*orig_comp*/) const
    {
       const int* domlo = geom.Domain().loVect();
       const int* domhi = geom.Domain().hiVect();
       for (int n = 0; n < numcomp; n++ ) {
          const int* bc = bcr[n].data();
          for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
             if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {
                dest(iv, dcomp+n) = 0.0;
             }
             if ((bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and (iv[idir] > domhi[idir])) {
                dest(iv, dcomp+n) = 0.0;
             }
          }
       }
    }
};

//
// A dummy function because FillPatch requires something to exist for filling dirichlet boundary conditions,
// even if we know we cannot have an ext_dir BC.
// u_mac BCs are only either periodic (BCType::int_dir) or first order extrapolation (FOEXTRAP).
//
struct umacFill
{
    AMREX_GPU_DEVICE
    void operator()(
       const amrex::IntVect& /*iv*/,
       amrex::Array4<amrex::Real> const& /*dummy*/,
       const int /*dcomp*/,
       const int numcomp,
       amrex::GeometryData const& /*geom*/,
       const amrex::Real /*time*/,
       const amrex::BCRec* bcr,
       const int bcomp,
       const int /*orig_comp*/) const
    {
        // Abort if this function is expected to fill an ext_dir BC.
        for (int n = bcomp; n < bcomp+numcomp; ++n) {
            const amrex::BCRec& bc = bcr[n];
            if ( AMREX_D_TERM(   bc.lo(0) == amrex::BCType::ext_dir || bc.hi(0) == amrex::BCType::ext_dir,
                              || bc.lo(1) == amrex::BCType::ext_dir || bc.hi(1) == amrex::BCType::ext_dir,
                              || bc.lo(2) == amrex::BCType::ext_dir || bc.hi(2) == amrex::BCType::ext_dir ) ) {
               amrex::Abort("NavierStokesBase::umacFill: umac should not have BCType::ext_dir");
            }
        }
    }
};


BCRec       NavierStokesBase::phys_bc;
Projection* NavierStokesBase::projector     = nullptr;
MacProj*    NavierStokesBase::mac_projector = nullptr;

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
int         NavierStokesBase::do_mac_proj               = 1;
int         NavierStokesBase::do_refine_outflow         = 0;
int         NavierStokesBase::do_derefine_outflow       = 1;
int         NavierStokesBase::Nbuf_outflow              = 1;
int         NavierStokesBase::do_denminmax              = 0;
int         NavierStokesBase::do_scalminmax             = 0;
int         NavierStokesBase::getForceVerbose           = 0;
int         NavierStokesBase::do_LES                    = 0;
int         NavierStokesBase::getLESVerbose             = 0;
std::string NavierStokesBase::LES_model                 = "Smagorinsky";
Real        NavierStokesBase::smago_Cs_cst              = 0.18;
Real        NavierStokesBase::sigma_Cs_cst              = 1.5;

amrex::Vector<amrex::Real> NavierStokesBase::time_avg;
amrex::Vector<amrex::Real> NavierStokesBase::time_avg_fluct;
amrex::Vector<amrex::Real> NavierStokesBase::dt_avg;
int  NavierStokesBase::avg_interval                    = 0;
int  NavierStokesBase::compute_fluctuations            = 0;
int  NavierStokesBase::additional_state_types_initialized = 0;
//
// "Divu_Type" means S, where divergence U = S
// "Dsdt_Type" means pd S/pd t, where S is as above
//
int  NavierStokesBase::Divu_Type                          = -1;
int  NavierStokesBase::Dsdt_Type                          = -1;
int  NavierStokesBase::Average_Type                       = -1;
int  NavierStokesBase::num_state_type                     = 2;
int  NavierStokesBase::have_divu                          = 0;
int  NavierStokesBase::have_dsdt                          = 0;
int  NavierStokesBase::do_init_vort_proj                  = 0;
int  NavierStokesBase::do_init_proj                       = 1;

int  NavierStokesBase::do_mom_diff            = 0;

std::string  NavierStokesBase::advection_scheme = "Godunov_PLM";

bool NavierStokesBase::godunov_use_forces_in_trans = false;

#ifdef AMREX_USE_EB
int          NavierStokesBase::refine_cutcells     = 1;
bool         NavierStokesBase::eb_initialized      = false;
bool         NavierStokesBase::no_eb_in_domain     = true;
bool         NavierStokesBase::body_state_set      = false;
std::vector<Real> NavierStokesBase::body_state;
std::string  NavierStokesBase::redistribution_type = "StateRedist";
#endif

//
// For restart, is GradP in checkpoint file
//
int NavierStokesBase::gradp_in_checkpoint = -1;

// is Average in checkpoint file
int NavierStokesBase::average_in_checkpoint = -1;

//
// skip_level_projector
//
int NavierStokesBase::skip_level_projector = 0;

int NavierStokesBase::isolver  = 0;

//
// ls related
//
int NavierStokesBase::do_phi   = 0;
int NavierStokesBase::phicomp  = 0;
int NavierStokesBase::epsilon  = 2;

Real NavierStokesBase::mu_a      = 0.00018;
Real NavierStokesBase::mu_w      = 0.01;
Real NavierStokesBase::rho_a     = 0.0012;
Real NavierStokesBase::rho_w     = 1.0;

int NavierStokesBase::do_reinit          = 0;
int NavierStokesBase::lev0step_of_reinit = 1;
int NavierStokesBase::number_of_reinit   = 4;
int NavierStokesBase::reinit_levelset    = 1;

int NavierStokesBase::do_cons_phi        = 0;
int NavierStokesBase::prescribed_vel     = 0;

int NavierStokesBase::do_cons_levelset   = 0;

//
// diffused ib
//
int NavierStokesBase::do_diffused_ib   = 0;
int NavierStokesBase::advect_and_update_scalar   = 1;
Real NavierStokesBase::fluid_rho   = 1.0;

namespace
{
    bool initialized = false;
    int  dump_plane  = -1;
    std::string dump_plane_name("SLABS/vel-");
    bool benchmarking = false;
}

#ifdef AMREX_PARTICLES
bool NavierStokesBase::do_nspc = true;
bool NavierStokesBase::particles_in_plotfile = false;

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
    int              pverbose                         = 0;
}

AmrTracerParticleContainer* NavierStokesBase::theNSPC () { return NSPC; }
#endif

NavierStokesBase::NavierStokesBase ()
{
    if (!additional_state_types_initialized) {
        init_additional_state_types();
    }
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

    //
    // 2/2022 - Only allow RZ if there's no visc.
    //   MLMG Tensor solver does not currently support RZ
    //
    if ( level_geom.IsRZ() )
    {
#ifdef AMREX_USE_EB
      amrex::Abort("Embedded boundaries with RZ geometry is not currently supported.");
#endif
        for ( int n = 0; n < AMREX_SPACEDIM; n++ ) {
            if ( visc_coef[n] > 0 ) {
                amrex::Abort("RZ geometry with viscosity is not currently supported. To use set ns.vel_visc_coef=0");
            }
        }
    }

    if(!additional_state_types_initialized) {
        init_additional_state_types();
    }

    //
    // Alloc old_time pressure.
    //
    state[Press_Type].allocOldData();
    state[Gradp_Type].allocOldData();

    define_workspace();
}

void NavierStokesBase::define_workspace()
{
    //
    // Alloc space for density and temporary pressure variables.
    //
    if (level > 0)
    {
        rho_avg.define(grids,dmap,1,1,MFInfo(),Factory());

        const BoxArray& P_grids = state[Press_Type].boxArray();
        p_avg.define(P_grids,dmap,1,0,MFInfo(),Factory());
    }

    rho_half.define (grids,dmap,1,1,MFInfo(),Factory());
    rho_ptime.define(grids,dmap,1,1,MFInfo(),Factory());
    rho_ctime.define(grids,dmap,1,1,MFInfo(),Factory());
    rho_qtime  = nullptr;
    rho_tqtime = nullptr;

    //
    // Build metric coefficients for RZ calculations.
    // Build volume and areas.
    //
    buildMetrics();

    //
    // ls related
    // 2 ghost cells
    //
    if (do_phi) {
        phi_half.define(grids,dmap,1,2,MFInfo(),Factory());
        phi_ptime.define(grids,dmap,1,2,MFInfo(),Factory());
        phi_ctime.define(grids,dmap,1,2,MFInfo(),Factory());
        heaviside.define(grids,dmap,1,2,MFInfo(),Factory());
        sgn0.define(grids,dmap,1,2,MFInfo(),Factory());
        phi_original.define(grids,dmap,1,2,MFInfo(),Factory());
    }

    //
    // diffused ib
    //
    if (do_diffused_ib) {
        const BoxArray& nba = amrex::convert(grids,IntVect::TheNodeVector());
        phi_nodal.define(nba,dmap,1,2,MFInfo(),Factory());
        pvf.define(grids,dmap,1,2,MFInfo(),Factory());
#ifdef AMREX_PARTICLES
        amrex::Print() << "check level " << level << " " << Particles::ParticleFinestLevel() << std::endl;
        if (level == Particles::ParticleFinestLevel()) {
            if(!Particles::isInitial){
                //particles
                Particles::init_particle(gravity, geom.CellSizeArray()[0]);
            }
            //largra
            Particles::create_particles(geom, dmap, grids); // Class constructor
        }
#endif
    }

    //
    // Set up reflux registers.
    //
        sync_reg = nullptr;
    viscflux_reg = nullptr;

    AMREX_ASSERT(sync_reg == nullptr);
    if (level > 0 && do_sync_proj)
    {
        sync_reg = new SyncRegister(grids,dmap,crse_ratio);
    }

    if (level > 0 && do_reflux)
    {
#ifdef AMREX_USE_EB
        advflux_reg  = std::make_unique<EBFluxRegister>();
#else
        advflux_reg  = std::make_unique<YAFluxRegister>();
#endif

        NavierStokesBase& clevel = getLevel(level-1);
        advflux_reg->define(grids, clevel.boxArray(), dmap, clevel.DistributionMap(),
                            parent->Geom(level),                 parent->Geom(level-1),
                            parent->refRatio(level-1),level,NUM_STATE);

#ifndef AMREX_USE_EB
        if (parent->Geom(level-1).IsRZ()) {
            advflux_reg->setCrseVolume(&(clevel.Volume()));
        }
#endif

        AMREX_ASSERT(viscflux_reg == nullptr);
        viscflux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
    }

    //
    // Set up the level projector.
    //
    if (projector == nullptr)
    {
        projector = new Projection(parent,&phys_bc,do_sync_proj,
                                   parent->finestLevel(),radius_grow);
    }
    projector->install_level(level,this,&radius);

    //
    // Set up the mac projector.
    //
    if (mac_projector == nullptr)
    {
        mac_projector = new MacProj(parent,parent->finestLevel(),
                                    &phys_bc,radius_grow);
    }
    mac_projector->install_level(level,this);

    //
    // Set up diffusion.
    //
    diffusion = std::make_unique<Diffusion>(parent,this,
                                            (level > 0) ? getLevel(level-1).diffusion.get()
                                                        : nullptr,
                                            NUM_STATE,viscflux_reg,is_diffusive);
    //
    // Allocate the storage for variable viscosity and diffusivity
    //
    diffn_cc = new MultiFab(grids, dmap, NUM_STATE-Density-1, 1, MFInfo(), Factory());
    diffnp1_cc = new MultiFab(grids, dmap, NUM_STATE-Density-1, 1, MFInfo(), Factory());
    viscn_cc = new MultiFab(grids, dmap, 1, 1, MFInfo(), Factory());
    viscnp1_cc = new MultiFab(grids, dmap, 1, 1, MFInfo(), Factory());

    //
    // Initialize BCRec for use with advection
    //
    m_bcrec_velocity.resize(AMREX_SPACEDIM);
    m_bcrec_velocity = fetchBCArray(State_Type,Xvel,AMREX_SPACEDIM);

    m_bcrec_velocity_d.resize(AMREX_SPACEDIM);
    m_bcrec_velocity_d = convertToDeviceVector(m_bcrec_velocity);

    m_bcrec_scalars.resize(NUM_SCALARS);
    m_bcrec_scalars = fetchBCArray(State_Type,Density,NUM_SCALARS);

    m_bcrec_scalars_d.resize(NUM_SCALARS);
    m_bcrec_scalars_d = convertToDeviceVector(m_bcrec_scalars);

}

NavierStokesBase::~NavierStokesBase ()
{
    delete rho_qtime;
    delete rho_tqtime;
    delete sync_reg;
    delete viscflux_reg;
    delete [] u_mac;

    if (mac_projector != nullptr) {
        mac_projector->cleanup(level);
    }
    //
    // Remove the arrays for variable viscosity and diffusivity
    // and delete the Diffusion object
    //
    delete viscn_cc;
    delete viscnp1_cc;
    delete diffn_cc;
    delete diffnp1_cc;

    diffusion.reset();
}

void
NavierStokesBase::variableCleanUp ()
{
    desc_lst.clear();
    derive_lst.clear();

    delete projector;
    projector = nullptr;

    delete mac_projector;
    mac_projector = nullptr;

#ifdef AMREX_PARTICLES
    delete NSPC;
    NSPC = nullptr;
#endif
}

void
NavierStokesBase::Initialize ()
{
    if (initialized) {
        return;
    }

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
    pp.query("gravity",gravity);
    //
    // Get run options.
    //
    pp.query("do_temp",                  do_temp          );
    pp.query("do_trac2",                 do_trac2         );
    pp.query("do_cons_trac",             do_cons_trac     );
    pp.query("do_cons_trac2",            do_cons_trac2    );
    pp.query("do_sync_proj",             do_sync_proj     );
    pp.query("do_reflux",                do_reflux        );
    pp.query("do_init_vort_proj",        do_init_vort_proj);
    pp.query("do_init_proj",             do_init_proj     );
    pp.query("do_mac_proj",              do_mac_proj      );
    pp.query("do_denminmax",             do_denminmax     );
    pp.query("do_scalminmax",            do_scalminmax    );

    if ( pp.contains("do_temp_ref") ||
         pp.contains("do_density_ref") ||
         pp.contains("do_tracer_ref") ||
         pp.contains("do_tracer2_ref") ||
         pp.contains("do_vorticity_ref") ) {
        amrex::Abort("ns.do_*_ref no longer supported. Refinement now implemented using refinement_indicators. For help, see UsersGuide or examples in /Exec");
    }

    pp.query("visc_tol",visc_tol);
    pp.query("visc_abs_tol",visc_abs_tol);

    pp.query("getForceVerbose",          getForceVerbose  );
    pp.query("do_LES",                   do_LES  );
    pp.query("getLESVerbose",            getLESVerbose  );
    pp.query("LES_model",                LES_model  );
    pp.query("smago_Cs_cst",             smago_Cs_cst  );
    pp.query("sigma_Cs_cst",             sigma_Cs_cst  );

    pp.query("avg_interval",             avg_interval  );
    pp.query("compute_fluctuations",     compute_fluctuations  );

#ifdef AMREX_USE_EB
    pp.query("refine_cutcells", refine_cutcells);
#endif

    int do_scalar_update_in_order = 0;
    pp.query("do_scalar_update_in_order",do_scalar_update_in_order );
    if (do_scalar_update_in_order) {
        amrex::Abort("NavierStokesBase::Initialize(): do_scalar_update_in_order no longer supported. If needed, please open issue on github.");
    }

    // Don't let init_shrink be greater than 1
    if (init_shrink > 1.0) {
        amrex::Abort("NavierStokesBase::Initialize(): init_shrink cannot be greater than 1");
    }

    pp.query("be_cn_theta",be_cn_theta);
    if (be_cn_theta > 1.0 || be_cn_theta < .5) {
        amrex::Abort("NavierStokesBase::Initialize(): Must have be_cn_theta <= 1.0 && >= .5");
    }
    //
    // Set parameters dealing with how grids are treated at outflow boundaries.
    //
    pp.query("do_refine_outflow",do_refine_outflow);
    pp.query("do_derefine_outflow",do_derefine_outflow);
    if (do_derefine_outflow == 1 && do_refine_outflow == 1) {
      amrex::Abort("NavierStokesBase::Initialize(): Cannot have both do_refine_outflow==1 and do_derefine_outflow==1");
    }

    pp.query("Nbuf_outflow",Nbuf_outflow);
    AMREX_ASSERT(Nbuf_outflow >= 0);
    AMREX_ASSERT(!(Nbuf_outflow <= 0 && do_derefine_outflow == 1));

    // Provide error message for depreciated volume weighted sum over a sub-domain.
    // NSB only ever supported cylinder sub-domains, so check for that one.
    if (pp.contains("volWgtSum_sub_dz")) {
        Abort("Computing volume weighted sum over sub-domains is no longer supported. If desired, submit an issue on github");
    };

    // Are we going to do velocity or momentum update?
    pp.query("do_mom_diff",do_mom_diff);

#ifdef AMREX_PARTICLES
    read_particle_params ();
#endif

    //
    // Get checkpoint info
    //
    pp.query("gradp_in_checkpoint", gradp_in_checkpoint);
    pp.query("avg_in_checkpoint",   average_in_checkpoint);

    //
    // Get advection scheme options
    //
    if ( pp.contains("use_godunov") ) {
        Abort("ns.use_godunov is depreciated. Please use ns.advection_scheme instead. Options are Godunov_PLM (default), Godunov_PPM, or BDS");
    }

    pp.query("advection_scheme", advection_scheme);
    if ( advection_scheme == "MOL" ) {
        Abort("MOL advection scheme is no longer supported. Current options are Godunov_PLM (default), Godunov_PPM, or BDS");
    }
    if (advection_scheme != "Godunov_PLM" && advection_scheme != "Godunov_PPM" && advection_scheme != "BDS") {
        Abort("Invalid advection_scheme. Options are Godunov_PLM, Godunov_PPM, BDS");
    }

    ParmParse pp2("godunov");
    pp2.query("use_forces_in_trans", godunov_use_forces_in_trans);

#ifdef AMREX_USE_EB
    //
    // Advection scheme restrictions
    //
    if ( advection_scheme == "Godunov_PPM" || advection_scheme == "BDS") {
        amrex::Abort("This advection_scheme is not implemented for EB. Please use Godunov_PLM (default)");
    }
    if ( godunov_use_forces_in_trans ) {
        amrex::Abort("use_forces_in_trans not implemented within EB Godunov. Set godunov.use_forces_in_trans=0.");
    }

    //
    // Redistribution
    //
    pp.query("redistribution_type", redistribution_type);
    if (redistribution_type != "NoRedist" &&
        redistribution_type != "FluxRedist" &&
        redistribution_type != "StateRedist" ) {
        amrex::Abort("redistribution type must be NoRedist, FluxRedist, or StateRedist");
    }
#endif

    //
    // composite time advancement
    //
    pp.query("skip_level_projector", skip_level_projector);

    //
    // ls related
    //
    pp.query("do_phi", do_phi);
    if (do_phi) {
        pp.query("epsilon", epsilon);
        pp.query("mu_a", mu_a);
        pp.query("mu_w", mu_w);
        pp.query("rho_a", rho_a);
        pp.query("rho_w", rho_w);

        pp.query("do_reinit", do_reinit);
        pp.query("lev0step_of_reinit", lev0step_of_reinit);
        pp.query("number_of_reinit", number_of_reinit);
        pp.query("reinit_levelset", reinit_levelset);

        pp.query("do_cons_phi", do_cons_phi);

    }

    pp.query("prescribed_vel", prescribed_vel);
    pp.query("isolver", isolver);
    pp.query("do_diffused_ib", do_diffused_ib);
    advect_and_update_scalar = !(do_diffused_ib == 1);
    pp.query("fluid_rho", fluid_rho);

    amrex::ExecOnFinalize(NavierStokesBase::Finalize);

#ifdef AMREX_PARTICLES
    Particles::Initialize();
#endif

    initialized = true;
}

void
NavierStokesBase::Finalize ()
{
    initialized = false;
}

void
NavierStokesBase::read_geometry ()
{
#if (AMREX_SPACEDIM == 2)
    //
    // Must load coord here because Geometry hasn't read it in yet.
    //
    ParmParse pp("geometry");

    int coord;
    pp.get("coord_sys",coord);

    if ((Geometry::CoordType) coord == Geometry::RZ && phys_bc.lo(0) != PhysBCType::symmetry)
    {
        phys_bc.setLo(0,PhysBCType::symmetry);
        amrex::Print() << "\nWarning: Setting phys_bc at xlo to PhysBCType::symmetry\n\n";
    }
#endif
}

void
NavierStokesBase::advance_setup (Real /*time*/,
                                 Real dt,
                                 int  iteration,
                                 int  ncycle)
{
    BL_PROFILE("NavierStokesBase::advance_setup()");

    const int finest_level = parent->finestLevel();

    // Same for EB vs not.
    umac_n_grow = 1;

#ifdef AMREX_PARTICLES
    if (ncycle > umac_n_grow && NSPC) {
        umac_n_grow = ncycle;
    }
#endif

    mac_projector->setup(level);
    //
    // Why are they defined here versus the constructor?
    //
    if (level < finest_level)
    {
#ifdef AMREX_USE_EB
        int ng_sync = (redistribution_type == "StateRedist") ? nghost_state() : 1;
#else
        int ng_sync = 1;
#endif

        if (Vsync.empty()) {
            Vsync.define(grids,dmap,AMREX_SPACEDIM,ng_sync,MFInfo(),Factory());
        }
        if (Ssync.empty()) {
            Ssync.define(grids,dmap,NUM_STATE-AMREX_SPACEDIM,ng_sync,MFInfo(),Factory());
        }
        Vsync.setVal(0);
        Ssync.setVal(0);
    }
    //
    // Set reflux registers to zero.
    //
    if (do_reflux && level < finest_level)
    {
        getAdvFluxReg(level+1).reset();
        getViscFluxReg(level+1).setVal(0.);
    }
    //
    // Alloc space for edge velocities (normal comp only).
    //
    if (u_mac == nullptr || u_mac[0].nGrow() < umac_n_grow)
    {
        delete [] u_mac;

        u_mac = new MultiFab[AMREX_SPACEDIM];

        for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            const BoxArray& edgeba = getEdgeBoxArray(dir);
            u_mac[dir].define(edgeba,dmap,1,umac_n_grow,MFInfo(),Factory());
            u_mac[dir].setVal(1.e40);
        }
    }
    //
    // Alloc MultiFab to hold advective update terms.
    //
    AMREX_ASSERT(aofs == nullptr);
    aofs = new MultiFab(grids,dmap,NUM_STATE,0,MFInfo(),Factory());

    //
    // Set rho_avg.
    //
    if (!initial_step && level > 0 && iteration == 1) {
        initRhoAvg(0.5/Real(ncycle));
    }
    //
    // Set up state multifabs for the advance.
    //
    for (int k = 0; k < num_state_type; k++)
    {
        bool has_old_data = state[k].hasOldData();
        // does nothing if old_data!=null
        state[k].allocOldData();
        if (! has_old_data) {
            state[k].oldData().setVal(0.0);
        }
        // swaps pointers-- reuses space, but doesn't leave new with good data.
        state[k].swapTimeLevels(dt);
    }

    //
    // ls related
    // fill the gts of old state data in the beginning
    // 
    if (do_phi) {
        // amrex::Print() << "1 " << std::endl;
        const Real prev_time = state[State_Type].prevTime();
        MultiFab&  S_old    = get_old_data(State_Type);
        int nScomp = S_old.nComp();
        fill_allgts(S_old,State_Type,0,nScomp,prev_time);
    }

    make_rho_prev_time();

    //
    // ls related
    // update the rho_ptime
    // 
    if (do_phi) {
        // amrex::Print() << "2 " << std::endl;
        MultiFab&  S_old    = get_old_data(State_Type);
        MultiFab::Copy(phi_ptime, S_old, phicomp, 0, 1, S_old.nGrow()); 
        phi_to_heavi(geom, epsilon, phi_ptime, heaviside);
        heavi_to_rhoormu(heaviside, rho_w, rho_a, rho_ptime);
        MultiFab::Copy(S_old, rho_ptime, 0, Density, 1, rho_ptime.nGrow());

        MultiFab outmf_mu_ptime(grids, dmap, 1, 1, MFInfo(), Factory());
        heavi_to_rhoormu(heaviside, mu_w, mu_a, outmf_mu_ptime);
        MultiFab::Copy(*viscn_cc, outmf_mu_ptime, 0, 0, 1, 1);
    }

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
    //             FillPatch(clevel, *(clevel.rho_tqtime), 1, tqtime, State_Type, Density, 1, 0);
    //         }
    //    }
    // }
}

//
// Clean up after the advance function.
//
void
NavierStokesBase::advance_cleanup (int /*iteration*/, int /*ncycle*/)
{
    delete aofs;
    aofs = nullptr;
}

void
NavierStokesBase::buildMetrics ()
{
    //
    // We "should" only need radius when we're RZ, but some 2-D code is written to
    // access it first and then "use" if if RZ.  It's easier to just always build
    // it for 2D than try to fix the underlying Fortran calls that take radius.
    //
#if (AMREX_SPACEDIM == 2)
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
        for (int j = 0; j < len; j++) {
            radius[i][j] = xlo + j*dxr;
        }
    }
#endif

    // volume and area are intentionally without EB knowledge
    volume.clear();
    volume.define(grids,dmap,1,GEOM_GROW);
    geom.GetVolume(volume);

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {
        area[dir].clear();
        area[dir].define(getEdgeBoxArray(dir),dmap,1,GEOM_GROW);
        geom.GetFaceArea(area[dir],dir);
    }

#ifdef AMREX_USE_EB
    // make sure dx == dy == dz
    const Real* dx = geom.CellSize();
    for (int i = 1; i < AMREX_SPACEDIM; i++){
        if (std::abs(dx[i]-dx[i-1]) > 1.e-12*dx[0]){
            Print()<<"dx = "
                   <<AMREX_D_TERM(dx[0], <<" "<<dx[1], <<" "<<dx[2])
                   <<std::endl;
            amrex::Abort("EB requires dx == dy (== dz)\n");
        }
    }

    const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());
    volfrac = &(ebfactory.getVolFrac());
    areafrac = ebfactory.getAreaFrac();

#endif
}

//
// Default dSdt is set to zero.
//
void
NavierStokesBase::calc_dsdt (Real      /*time*/,
                         Real      dt,
                         MultiFab& dsdt)
{
    if (have_divu && have_dsdt)
    {
        // Don't think we need this here. Instead, code will use FillPatch to
        // fill ghosts.
        //dsdt.setVal(0);

        if (do_temp) {
            MultiFab& Divu_new = get_new_data(Divu_Type);
            MultiFab& Divu_old = get_old_data(Divu_Type);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(dsdt,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box&  bx      = mfi.tilebox();
                auto const& div_new = Divu_new.array(mfi);
                auto const& div_old = Divu_old.array(mfi);
                auto const& dsdtarr = dsdt.array(mfi);

                amrex::ParallelFor(bx, [div_new, div_old, dsdtarr, dt]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    dsdtarr(i,j,k) = ( div_new(i,j,k) - div_old(i,j,k) )/ dt;
                });
            }
        }
        else
        {
            dsdt.setVal(0);
        }
    }
}

void
NavierStokesBase::checkPoint (const std::string& dir,
                              std::ostream&      os,
                              VisMF::How         how,
                              bool               dump_old)
{
    AmrLevel::checkPoint(dir, os, how, dump_old);

    if (avg_interval > 0)
    {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

        if (ParallelDescriptor::IOProcessor())
        {
            std::ofstream TImeAverageFile;
            TImeAverageFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
            std::string TAFileName(dir + "/TimeAverage");
            TImeAverageFile.open(TAFileName.c_str(), std::ofstream::out |
                                 std::ofstream::trunc |
                                 std::ofstream::binary);

            if( !TImeAverageFile.good()) {
                amrex::FileOpenFailed(TAFileName);
            }

            TImeAverageFile.precision(17);

            // write out title line
            TImeAverageFile << "Writing time_average to checkpoint\n";

            TImeAverageFile << NavierStokesBase::time_avg[level] << "\n";
            TImeAverageFile << NavierStokesBase::time_avg_fluct[level] << "\n";
        }
    }

#ifdef AMREX_PARTICLES
    if (level == 0)
    {
        if (NSPC != 0)
            NSPC->Checkpoint(dir,the_ns_particle_file_name);
    }
#endif
}

void
NavierStokesBase::computeInitialDt (int                   finest_level,
                                    int                   /*sub_cycle*/,
                                    Vector<int>&           n_cycle,
                                    const Vector<IntVect>& /*ref_ratio*/,
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
                                int                   /*sub_cycle*/,
                                Vector<int>&           n_cycle,
                                const Vector<IntVect>& /*ref_ratio*/,
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

    AMREX_ASSERT(rhs.nGrow()>=nGrow);
    AMREX_ASSERT(rhs.boxArray()==grids);

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
NavierStokesBase::create_umac_grown (int nGrow,
                                     const MultiFab* a_divu)
{
    if ( nGrow <= 0 ) { return; }
    if ( nGrow > 1 )  { Print()<<"\n\nWARNING!\n  NSB::create_umac_grown currently only enforces the divergnece constraint on 1 ghost cell, but nGrow > 1\n\n"; }


    Array<MultiFab*, AMREX_SPACEDIM> u_mac_fine;
    AMREX_D_TERM(u_mac_fine[0] = &u_mac[0];,
                 u_mac_fine[1] = &u_mac[1];,
                 u_mac_fine[2] = &u_mac[2];);

    Geometry *fine_geom = &geom;

    //Grab the velocity phys bc fill function from the StateData StateDescriptor
    AMREX_D_TERM(StateDataPhysBCFunct fine_bndry_func_x(get_state_data(State_Type),0,geom);,
                 StateDataPhysBCFunct fine_bndry_func_y(get_state_data(State_Type),1,geom);,
                 StateDataPhysBCFunct fine_bndry_func_z(get_state_data(State_Type),2,geom););

    Array<StateDataPhysBCFunct,AMREX_SPACEDIM> fbndyFuncArr = {AMREX_D_DECL(fine_bndry_func_x,
                                                                            fine_bndry_func_y,
                                                                            fine_bndry_func_z)};


    if ( level == 0 )
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            //
            // BDS needs physical BCs filled.
            // Godunov needs periodic and coarse-fine ghosts filled (and handles
            // physical BCs internally).
            //
            Real fake_time = 0.;
            amrex::FillPatchSingleLevel(u_mac[idim], IntVect(nGrow), fake_time,
                                        {u_mac_fine[idim]}, {fake_time},
                                        0, 0, 1, geom,
                                        fbndyFuncArr[idim], 0);
        }
    }
    else // level > 0
    {
        Array<MultiFab*, AMREX_SPACEDIM> u_mac_crse;
        AMREX_D_TERM(u_mac_crse[0] = &getLevel(level-1).u_mac[0];,
                     u_mac_crse[1] = &getLevel(level-1).u_mac[1];,
                     u_mac_crse[2] = &getLevel(level-1).u_mac[2];);

        Geometry *crse_geom = &getLevel(level-1).geom;

        //
        // First interpolate, ignoring divergence constraint. Then correct
        // the 1-cell wide halo of ghosts cells we need to enforce the
        // constraint.
        //

        // Divergence preserving interp -- This is for case of MAC solve on
        // composite grid; doesn't make sense to use it here.
        //Interpolater* mapper = &face_divfree_interp;
        // Use linear interpolation, which matches up with old create umac grown
        Interpolater* mapper = &face_linear_interp;

        // This never actually gets used because the FaceLinear Interpolator doesn't use it.
        // Inside FillPatchTwoLevels, it's only use is that it's passed to the interpolator.
        // Fill it with the correct thing anyway.
        Array<Vector<BCRec>,AMREX_SPACEDIM> bcrecArr = {AMREX_D_DECL(m_bcrec_velocity,
                                                                     m_bcrec_velocity,
                                                                     m_bcrec_velocity)};
        Array<int, AMREX_SPACEDIM> bc_idx = {AMREX_D_DECL(0,1,2)};

        // Grab the velocity phys bc fill function from the StateData StateDescriptor
        // This will use the BCRec for velocity stored in the StateDescriptor
        AMREX_D_TERM(
            StateDataPhysBCFunct crse_bndry_func_x(getLevel(level-1).get_state_data(State_Type),0,*crse_geom);,
            StateDataPhysBCFunct crse_bndry_func_y(getLevel(level-1).get_state_data(State_Type),1,*crse_geom);,
            StateDataPhysBCFunct crse_bndry_func_z(getLevel(level-1).get_state_data(State_Type),2,*crse_geom););
        Array<int, AMREX_SPACEDIM> bf_idx = {AMREX_D_DECL(0,0,0)};

        Array<StateDataPhysBCFunct,AMREX_SPACEDIM> cbndyFuncArr = {AMREX_D_DECL(crse_bndry_func_x,
                                                                                crse_bndry_func_y,
                                                                                crse_bndry_func_z)};

        // Use piecewise constant interpolation in time, so create fictitious variable for time
        Real fake_time = 0.;

        FillPatchTwoLevels(u_mac_fine, IntVect(nGrow), fake_time,
                           {u_mac_crse}, {fake_time},
                           {u_mac_fine}, {fake_time},
                           0, 0, 1,
                           *crse_geom, *fine_geom,
                           cbndyFuncArr, bf_idx, fbndyFuncArr, bf_idx,
                           crse_ratio, mapper, bcrecArr, bc_idx);

        //
        // Correct u_mac to enforce the divergence constraint in the ghost cells.
        // Do this by adjusting only the outer face (wrt the valid region) of the ghost
        // cell, i.e. for the hi-x face, adjust umac_x(i+1).
        // NOTE that this does not fill grid edges or corners.
        //

        // Need 2 ghost cells here so we can safely check the status of all faces of a
        // u_mac ghost cell
        if (coarse_fine_mask == nullptr) {
            coarse_fine_mask = std::make_unique<iMultiFab>(grids, dmap, 1, 2, MFInfo(), DefaultFabFactory<IArrayBox>());
            coarse_fine_mask->BuildMask(geom.Domain(), geom.periodicity(),
                                        level_mask_covered, level_mask_notcovered, level_mask_physbnd, level_mask_interior);
        }

        const GpuArray<Real,AMREX_SPACEDIM> dx = fine_geom->CellSizeArray();
        const GpuArray<Real,AMREX_SPACEDIM> dxinv = fine_geom->InvCellSizeArray();

        const bool is_rz = geom.IsRZ();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*coarse_fine_mask,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& tbx = mfi.tilebox();
            auto const& maskarr = coarse_fine_mask->const_array(mfi);
            auto const& divu = (a_divu) ? a_divu->const_array(mfi) : Array4<const Real> {};
            auto const& umac = u_mac_fine[0]->array(mfi);
            auto const& vmac = u_mac_fine[1]->array(mfi);
            auto const& wmac = (AMREX_SPACEDIM==3) ? u_mac_fine[2]->array(mfi) : Array4<Real> {};

            const auto& vol = (is_rz) ?  volume.const_array(mfi): Array4<Real> {};
            const auto&  ax = (is_rz) ? area[0].const_array(mfi): Array4<Real> {};
            const auto&  ay = (is_rz) ? area[1].const_array(mfi): Array4<Real> {};
#if (AMREX_SPACEDIM == 2)
            int ks =  0;
            int ke =  0;
#else
            int ks = -1;
            int ke =  1;
#endif

            AMREX_HOST_DEVICE_FOR_3D(mfi.growntilebox(1), i, j, k,
            {
                if ( maskarr(i,j,k) == level_mask_notcovered )
                {
                    //
                    // Leave cells on grid edges/corners unaltered.
                    // This correction scheme doesn't work for concave edges where a cell
                    // has faces that are all either valid or touching another ghost cell
                    // because then there's no "free" face to absorb the divergence constraint
                    // error.
                    // There are (>=) 1 case that are treatable, but we don't implement here:
                    // 1. Convex grid edges (cells that don't have any valid faces), e.g. by
                    //    dividing the divergence constraint error equally between each face
                    //    not touching another ghost cell
                    //

                    int count = 0;
                    for(int kk(ks); kk<=ke; kk++)
                    {
                        for(int jj(-1); jj<=1; jj++) {
                            for(int ii(-1); ii<=1; ii++) {
                                if ( Math::abs(ii)+Math::abs(jj)+Math::abs(kk) == 1 &&
                                     (maskarr(i+ii,j+jj,k+kk) == interior || maskarr(i+ii,j+jj,k+kk) == level_mask_covered) )
                                {
                                    count++;
                                }
                            }
                        }
                    }

                    if ( count == 1 )
                    {
                        Real div = (divu) ? divu(i,j,k) : 0.0;

                        if (is_rz)
                        {
                            Real dux = (ax(i+1,j,k)*umac(i+1,j,k) - ax(i,j,k)*umac(i,j,k));
                            Real duy = (ay(i,j+1,k)*vmac(i,j+1,k) - ay(i,j,k)*vmac(i,j,k));

                            // To avoid inconsistencies between boxes, we make sure to fix box
                            // corners (2D) or edges (3D) that are not grid corners/edges.
                            // The directional check ensures we only alter one face of these
                            // cells.
                            // It's unlikely there'd ever be a case of such ghost cells abutting the
                            // symmetry axis, but just in case, check here.
                            if ( i < tbx.smallEnd(0) && maskarr(i+1,j,k) != level_mask_notcovered && ax(i,j,k) != Real(0.0) )
                            {
                                umac(i,j,k) = (ax(i+1,j,k)*umac(i+1,j,k) + (duy - vol(i,j,k)*div))/ax(i,j,k);
                            }
                            else if ( i > tbx.bigEnd(0) && maskarr(i-1,j,k) != level_mask_notcovered )
                            {
                                umac(i+1,j,k) = (ax(i,j,k)*umac(i,j,k) - (duy - vol(i,j,k)*div))/ax(i+1,j,k);
                            }

                            if ( j < tbx.smallEnd(1) && maskarr(i,j+1,k) != level_mask_notcovered )
                            {
                                vmac(i,j,k) = (ay(i,j+1,k)*vmac(i,j+1,k) + (dux - vol(i,j,k)*div))/ay(i,j,k);
                            }
                            else if ( j > tbx.bigEnd(1) && maskarr(i,j-1,k) != level_mask_notcovered )
                            {
                                vmac(i,j+1,k) = (ay(i,j,k)*vmac(i,j,k) - (dux - vol(i,j,k)*div))/ay(i,j+1,k);
                            }
                        }
                        else
                        {
                            Real dux =          dxinv[0]*(umac(i+1,j,k) - umac(i,j,k));
                            Real duy =          dxinv[1]*(vmac(i,j+1,k) - vmac(i,j,k));
                            Real duz = (wmac) ? dxinv[2]*(wmac(i,j,k+1) - wmac(i,j,k)) : 0.0;

                            // To avoid inconsistencies between boxes, we make sure to fix box
                            // corners (2D) or edges (3D) that are not grid corners/edges.
                            // The directional check ensures we only alter one face of these
                            // cells.
                            if ( i < tbx.smallEnd(0) && maskarr(i+1,j,k) != level_mask_notcovered )
                            {
                                umac(i,j,k) = umac(i+1,j,k) + dx[0] * (duy + duz - div);
                            }
                            else if ( i > tbx.bigEnd(0) && maskarr(i-1,j,k) != level_mask_notcovered )
                            {
                                umac(i+1,j,k) = umac(i,j,k) - dx[0] * (duy + duz - div);
                            }

                            if ( j < tbx.smallEnd(1) && maskarr(i,j+1,k) != level_mask_notcovered )
                            {
                                vmac(i,j,k) = vmac(i,j+1,k) + dx[1] * (dux + duz - div);
                            }
                            else if ( j > tbx.bigEnd(1) && maskarr(i,j-1,k) != level_mask_notcovered )
                            {
                                vmac(i,j+1,k) = vmac(i,j,k) - dx[1] * (dux + duz - div);
                            }

                            if (wmac)
                            {
                                if ( k < tbx.smallEnd(2) && maskarr(i,j,k+1) != level_mask_notcovered )
                                {
                                    wmac(i,j,k) = wmac(i,j,k+1) + dx[2] * (dux + duy - div);
                                }
                                else if ( k > tbx.bigEnd(2) && maskarr(i,j,k-1) != level_mask_notcovered )
                                {
                                    wmac(i,j,k+1) = wmac(i,j,k) - dx[2] * (dux + duy - div);
                                }
                            }
                        }
                    }
                }
            });
        }
    }
}

void
NavierStokesBase::diffuse_scalar_setup (int sigma, int& rho_flag)
{
    rho_flag = Diffusion::set_rho_flag(diffusionType[sigma]);
}

void
NavierStokesBase::errorEst (TagBoxArray& tb,
                            int          /*clearval*/,
                            int          /*tagval*/,
                            Real         /*time*/,
                            int          /*n_error_buf*/,
                            int          /*ngrow*/)
{
#ifdef AMREX_USE_EB
    // Enforce that the EB not cross the coarse-fine boundary
    const auto& ebfactory = dynamic_cast<amrex::EBFArrayBoxFactory const&>(Factory());
    if ( !ebfactory.isAllRegular() )
    {
        if (!refine_cutcells) {
            amrex::Warning("Partially refined EB is still under development. This is not guaranteed to work! Please use ns.refine_cutcells=1 for now.");
        }

        // Refine on cut cells
        if (refine_cutcells)
        {
            const MultiFab& S_new = get_new_data(State_Type);
            amrex::TagCutCells(tb, S_new);
        }
    }
#else
    amrex::ignore_unused(tb);
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
            factor = 1.0/Real(ratio);
        }

        return factor*fixed_dt;
    }

    const Real  small         = 1.0e-8;
    Real        estdt         = 1.0e+20;

    MultiFab&   S_new         = get_new_data(State_Type);

    Vector<Real> u_max(AMREX_SPACEDIM);
    Vector<Real> f_max(AMREX_SPACEDIM);

    MultiFab& Gp = get_new_data(Gradp_Type);

    //
    // Find local max of velocity
    //
    u_max = S_new.norm0({AMREX_D_DECL(0,1,2)},0,true,true);

    //
    // Compute forcing terms: in this case this means external forces and grad(p)
    // Viscous terms not included since Crack-Nicholson is unconditionally stable
    // so no need to account for explicit part of viscous term
    //
    MultiFab tforces(grids,dmap,AMREX_SPACEDIM,0,MFInfo(),Factory());

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       const auto& bx          = mfi.tilebox();
       const auto  cur_time    = state[State_Type].curTime();
       auto& tforces_fab       = tforces[mfi];

       if (getForceVerbose) {
           amrex::Print() << "---" << '\n'
                          << "H - est Time Step:" << '\n'
                          << "Calling getForce..." << '\n';
       }
       getForce(tforces_fab,bx,0,AMREX_SPACEDIM,cur_time,S_new[mfi],S_new[mfi],Density,mfi);

       const auto& rho   = S_new.array(mfi,Density);
       const auto& gradp = Gp.array(mfi);
       const auto& force = tforces.array(mfi);
       amrex::ParallelFor(bx, [rho, gradp, force]
       AMREX_GPU_DEVICE(int i, int j, int k) noexcept
       {
          Real rho_inv = 1.0/rho(i,j,k);
          for (int n = 0; n < AMREX_SPACEDIM; n++) {
             force(i,j,k,n) -= gradp(i,j,k,n);
             force(i,j,k,n) *= rho_inv;
          }
       });
    }

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
         factor = 1.0/Real(ratio);
      }

      estdt = factor*init_dt;
    } else {
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

        if (getForceVerbose) {
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
        AMREX_ASSERT(rho_qtime);
        return *rho_qtime;
    }
    else if (whichTime == Amr3QtrTime)
    {
        AMREX_ASSERT(rho_tqtime);
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
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rho_half,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();
        auto const& rho_h = rho_half.array(mfi);
        auto const& rho_p = rho_ptime.array(mfi);
        auto const& rho_c = rho_ctime.array(mfi);
        amrex::ParallelFor(bx, [rho_h, rho_p, rho_c]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
           rho_h(i,j,k) = 0.5 * (rho_p(i,j,k) + rho_c(i,j,k));
        });
    }
    return rho_half;
}

//
// Fill patch divU.
//
MultiFab*
NavierStokesBase::getDivCond (int ngrow, Real time)
{
    MultiFab* divu = nullptr;

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
    MultiFab* dsdt = nullptr;

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

//
// Fill patch a state component.
//
MultiFab*
NavierStokesBase::getState (int  ngrow,
                            int  state_indx,
                            int  strt_comp,
                            int  num_comp,
                            Real time)
{
    BL_PROFILE("NavierStokesBase::getState()");

    auto* mf = new MultiFab(state[state_indx].boxArray(),
                            state[state_indx].DistributionMap(),
                            num_comp,ngrow,MFInfo(),Factory());

    FillPatch(*this,*mf,ngrow,time,state_indx,strt_comp,num_comp,0);

    return mf;
}

void
NavierStokesBase::getOutFlowFaces (Vector<Orientation>& outFaces)
{
    outFaces.resize(0);
    for (int idir = 0; idir < AMREX_SPACEDIM; idir++)
    {
        if (phys_bc.lo(idir) == PhysBCType::outflow)
        {
            auto len = outFaces.size();
            outFaces.resize(len+1);
            outFaces[len] = Orientation(idir,Orientation::low);
        }

        if (phys_bc.hi(idir) == PhysBCType::outflow)
        {
            auto len = outFaces.size();
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

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rho_avg,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.tilebox();
       auto const& rhoavg     = rho_avg.array(mfi);
       auto const& rho_new    = S_new.array(mfi,Density);
       amrex::ParallelFor(bx, [rhoavg,rho_new,alpha]
       AMREX_GPU_DEVICE(int i, int j, int k) noexcept
       {
          rhoavg(i,j,k) = rho_new(i,j,k) * alpha;
       });
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
    auto* oldns = dynamic_cast<NavierStokesBase*>(&old);
    const Real    dt_new    = parent->dtLevel(level);
    const Real    cur_time  = oldns->state[State_Type].curTime();
    const Real    prev_time = oldns->state[State_Type].prevTime();
    const Real    dt_old    = cur_time - prev_time;
    MultiFab&     S_new     = get_new_data(State_Type);
    MultiFab&     P_new     = get_new_data(Press_Type);
    MultiFab&     Gp_new    = get_new_data(Gradp_Type);

    setTimeLevel(cur_time,dt_old,dt_new);

    const Real cur_pres_time = state[Press_Type].curTime();
    //
    // Get best state and pressure data.
    //
    FillPatch(old,S_new,0,cur_time,State_Type,0,NUM_STATE);
    FillPatch(old,P_new,0,cur_pres_time,Press_Type,0,1);
    FillPatch(old,Gp_new,Gp_new.nGrow(),cur_pres_time,Gradp_Type,0,AMREX_SPACEDIM);

    if (avg_interval > 0){
      MultiFab& Save_new = get_new_data(Average_Type);
      FillPatch(old,Save_new,0,cur_time,Average_Type,0,AMREX_SPACEDIM*2);
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
}

//
// Fills a totally new level n with data interpolated from coarser level.
//
void
NavierStokesBase::init ()
{
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& P_new = get_new_data(Press_Type);
    MultiFab& Gp_new = get_new_data(Gradp_Type);

    AMREX_ASSERT(level > 0);

    const Vector<Real>& dt_amr = parent->dtLevel();
    Vector<Real>        dt_new(level+1);

    for (int lev = 0; lev < level; lev++) {
        dt_new[lev] = dt_amr[lev];
    }
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
    FillCoarsePatch(Gp_new,0,cur_pres_time,Gradp_Type,0,AMREX_SPACEDIM,Gp_new.nGrow());
    //
    // Get best coarse divU and dSdt data.
    //
    if (have_divu)
    {
        FillCoarsePatch(get_new_data(Divu_Type),0,cur_time,Divu_Type,0,1);
        if (have_dsdt)
            FillCoarsePatch(get_new_data(Dsdt_Type),0,cur_time,Dsdt_Type,0,1);
    }
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
    amrex::ignore_unused(have_temp);
    AMREX_ASSERT((do_temp && have_temp)  ||  (!do_temp && !have_temp));

    int i_Divu = -1;
    int dummy_Divu_Type;
    have_divu = 0;
    have_divu = isStateVariable("divu", dummy_Divu_Type, i_Divu);
    have_divu = have_divu && dummy_Divu_Type == Divu_Type;
    if (verbose)
    {
        amrex::Print() << "NavierStokesBase::init_additional_state_types()::have_divu = "
                  << have_divu << '\n';
    }
    if (have_divu && i_Divu!=Divu)
    {
        amrex::Print() << "divu must be 0-th Divu_Type component in the state\n";
        amrex::Abort("NavierStokesBase::init_additional_state_types()");
    }

    int i_Dsdt = -1;
    int dummy_Dsdt_Type;
    have_dsdt = 0;
    have_dsdt = isStateVariable("dsdt", dummy_Dsdt_Type, i_Dsdt);
    have_dsdt = have_dsdt && dummy_Dsdt_Type==Dsdt_Type;
    if (verbose)
    {
        amrex::Print() << "NavierStokesBase::init_additional_state_types()::have_dsdt = "
                       << have_dsdt << '\n';
    }
    if (have_dsdt && i_Dsdt!=Dsdt)
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

    if (verbose) amrex::Print() << "Multiplying dt by init_shrink: dt = "
                   << returnDt << '\n';
    return returnDt;
}

//
// Since the pressure solver always stores its estimate of the
// pressure solver in Pnew, we need to copy it to Pold at the start.
//
void
NavierStokesBase::initOldFromNew (int type, int lev)
{
    if ( lev < 0 ) {
      lev = level;
    }

    MultiFab& new_t = getLevel(lev).get_new_data(type);
    MultiFab& old_t = getLevel(lev).get_old_data(type);

    MultiFab::Copy(old_t, new_t, 0, 0, old_t.nComp(), old_t.nGrow());
}

void
NavierStokesBase::level_projector (Real dt,
                                   Real time,
                                   int  iteration)
{
    BL_PROFILE_REGION_START("R::NavierStokesBase::level_projector()");
    BL_PROFILE("NavierStokesBase::level_projector()");

    AMREX_ASSERT(iteration > 0);

    MultiFab& U_old = get_old_data(State_Type);
    MultiFab& U_new = get_new_data(State_Type);
    MultiFab& P_old = get_old_data(Press_Type);
    MultiFab& P_new = get_new_data(Press_Type);

    SyncRegister* crse_ptr = nullptr;

    if (level < parent->finestLevel() && do_sync_proj)
    {
        crse_ptr = &(getLevel(level+1).getSyncReg());
    }

    int        crse_dt_ratio  = (level > 0) ? parent->nCycle(level) : -1;
    const Real cur_pres_time  = state[Press_Type].curTime();

    projector->level_project(level,time,dt,cur_pres_time,
                             geom,U_old,U_new,P_old,P_new,
                             get_rho_half_time(),crse_ptr,sync_reg,
                             crse_dt_ratio,iteration,have_divu);

    BL_PROFILE_REGION_STOP("R::NavierStokesBase::level_projector()");
}

void
NavierStokesBase::level_sync (int crse_iteration)
{
    BL_PROFILE_REGION_START("R::NavierStokesBase::level_sync()");
    BL_PROFILE("NavierStokesBase::level_sync()");

    IntVect         ratio         = parent->refRatio(level);
    const int       finest_level  = parent->finestLevel();
    int             crse_dt_ratio = parent->nCycle(level);
    Real            dt            = parent->dtLevel(level);
    MultiFab&       pres          = get_new_data(Press_Type);
    MultiFab&       vel           = get_new_data(State_Type);
    SyncRegister&   rhs_sync_reg  = getLevel(level+1).getSyncReg();
    SyncRegister*   crsr_sync_ptr = nullptr;
    NavierStokesBase& fine_level  = getLevel(level+1);
    MultiFab&       pres_fine     = fine_level.get_new_data(Press_Type);
    MultiFab&       vel_fine      = fine_level.get_new_data(State_Type);
    const BoxArray& finegrids     = vel_fine.boxArray();
    const DistributionMapping& finedmap = vel_fine.DistributionMap();

    if (level > 0) {
        crsr_sync_ptr = &(getLevel(level).getSyncReg());
    }
    //
    // Get boundary conditions.
    //
    const auto N = int(grids.size());

    Vector<int*>         sync_bc(N);
    Vector< Vector<int> > sync_bc_array(N);

    for (int i = 0; i < N; i++)
    {
        sync_bc_array[i] = getBCArray(State_Type,i,Xvel,AMREX_SPACEDIM);
        sync_bc[i] = sync_bc_array[i].dataPtr();
    }

    //
    // Multilevel sync projection.
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
    MultiFab V_corr(finegrids,finedmap,AMREX_SPACEDIM,1,MFInfo(),fine_level.Factory());

    V_corr.setVal(0);
    //
    // If periodic, enforce periodicity on Vsync.
    //
    if (crse_geom.isAnyPeriodic()) {
      Vsync.FillBoundary(0, AMREX_SPACEDIM, crse_geom.periodicity());
    }
    //
    // Interpolate Vsync to fine grid correction in Vcorr.
    //
    SyncInterp(Vsync, level, V_corr, level+1, ratio,
               0, 0, AMREX_SPACEDIM, 0 , dt, sync_bc.dataPtr());
    //
    // The multilevel projection.  This computes the projection and
    // adds in its contribution to levels (level) and (level+1).
    //
    projector->MLsyncProject(level,pres,vel,cc_rhs_crse,
                             pres_fine,v_fine,cc_rhs_fine,
                             Rh,rho_fine,Vsync,V_corr,
                             phi,&rhs_sync_reg,crsr_sync_ptr,
                             dt,ratio,crse_iteration,crse_dt_ratio,
                             geom);
    cc_rhs_crse.clear();
    cc_rhs_fine.clear();
    //
    // Correct pressure and velocities after the projection.
    //
    const auto Nf = int(finegrids.size());

    ratio = IntVect::TheUnitVector();

    Vector<int*>         fine_sync_bc(Nf);
    Vector< Vector<int> > fine_sync_bc_array(Nf);

    for (int i = 0; i < Nf; i++)
    {
      fine_sync_bc_array[i] = getLevel(level+1).getBCArray(State_Type,
                                                           i,
                                                           Xvel,
                                                           AMREX_SPACEDIM);
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
                 0, 0, AMREX_SPACEDIM, 1 , dt, fine_sync_bc.dataPtr());
      SyncProjInterp(phi, level+1, P_new, P_old, lev, ratio);

      flev.computeGradP(flev.state[Gradp_Type].prevTime());
      flev.computeGradP(flev.state[Gradp_Type].curTime());

    }

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
                               MultiFab& S_old,
                               MultiFab* divu,
                               int       ngrow,
                               bool      increment_vel_register)
{
    BL_PROFILE_REGION_START("R::NavierStokesBase::mac_project()");
    BL_PROFILE("NavierStokesBase::mac_project()");

    if (verbose) {
        amrex::Print() << "... mac_projection\n";
    }

    if (verbose && benchmarking) {
        ParallelDescriptor::Barrier();
    }

    const Real strt_time = ParallelDescriptor::second();

    Vector<BCRec> density_math_bc = fetchBCArray(State_Type,Density,1);

    mac_projector->mac_project(level,u_mac,S_old,dt,time,*divu,have_divu,
                               density_math_bc[0], increment_vel_register);

    create_umac_grown(ngrow, divu);

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
    if (! outFaces.empty())
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
                bool hasTags = tags.hasTags(outflowBox);
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
                    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
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

    MultiFab&   u_old = get_old_data(State_Type);
    MultiFab&   u_new = get_new_data(State_Type);

        //
        // Estimate the maximum change in velocity magnitude since previous
        // iteration
        //
    Real max_change = 0.0;

    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    // Do not OpenMP-fy this loop for now
    // Unclear how to keep OpenMP and GPU implementation
    // from messing with each other
    for (MFIter mfi(u_old,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto& bx   = mfi.tilebox();
        const auto& uold = u_old[mfi].array();
        const auto& unew = u_new[mfi].array();

        reduce_op.eval(bx, reduce_data, [uold, unew]
        AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real uold_mag = 0.0;
            Real unew_mag = 0.0;
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
            {
                uold_mag += uold(i,j,k,d)*uold(i,j,k,d);
                unew_mag += unew(i,j,k,d)*unew(i,j,k,d);
            }

            uold_mag = std::sqrt(uold_mag);
            unew_mag = std::sqrt(unew_mag);

            return std::abs(unew_mag-uold_mag);
        });

        max_change = std::max(amrex::get<0>(reduce_data.value()),
                              max_change);
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

    if (do_init_vort_proj)
    {
        amrex::Abort("NavierStokesBase::post_init_state(): initialVorticityProject not tested with new Gradp!!! See comments Projection::initialVorticityProject.\n");
        //
        // NOTE: this assumes have_divu == 0.
        // Only used if vorticity is used to initialize the velocity field.
        //
        AMREX_ASSERT(projector != nullptr);

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
    // so that conserved data is consistent between levels.
    // This might not be the most efficient way of doing things
    // (since initialVelocityProject will average down vel, P and Gradp),
    // but it does ensure everything is averaged down for all cases
    // (e.g. initialVelocityProject doesn't get called or init_vel_iter<=0).
    //
    for (int k = finest_level-1; k>= 0; k--)
    {
      getLevel(k).avgDown();
    }

    if (do_init_proj && projector && (std::fabs(gravity)) > 0.){
      //
      // Do projection to establish initially hydrostatic pressure field.
      //
      if (verbose) amrex::Print() << "calling initialPressureProject" << std::endl;

      projector->initialPressureProject(0);

      if (verbose) amrex::Print() << "done calling initialPressureProject" << std::endl;
    }
    //
    // Make sure there's not NANs in old pressure field.
    // End up with P_old = P_new as is the case when exiting initialPressureProject
    //
    if(!do_init_proj)
    {
      for (int k = finest_level; k>= 0; k--)
      {
        initOldFromNew(Press_Type, k);
        initOldFromNew(Gradp_Type, k);
    }
    }
}

//
// Build any additional data structures after regrid.
//
void
NavierStokesBase::post_regrid (int lbase,
                               int /*new_finest*/)
{
#ifdef AMREX_PARTICLES
    if (NSPC && level == lbase)
    {
        NSPC->Redistribute(lbase);
    }
#else
    amrex::ignore_unused(lbase);
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

  if (avg_interval > 0){

    const int   finest_level = parent->finestLevel();
    NavierStokesBase::time_avg.resize(finest_level+1);
    NavierStokesBase::time_avg_fluct.resize(finest_level+1);
    NavierStokesBase::dt_avg.resize(finest_level+1);

    //
    // We assume that if Average_Type is not present, we have just activated
    // the start of averaging
    //
    if ( average_in_checkpoint==0 )
    {
      Print()<<"WARNING! Average not found in checkpoint file. Creating data"
             <<std::endl;

      Real cur_time = state[State_Type].curTime();
      Real prev_time = state[State_Type].prevTime();
      Real dt = cur_time - prev_time;
      state[Average_Type].define(geom.Domain(), grids, dmap, desc_lst[Average_Type],
                               cur_time, dt, Factory());

      MultiFab& Savg   = get_new_data(Average_Type);
      Savg.setVal(0.);
      state[Average_Type].allocOldData();
      MultiFab& Savg_old   = get_old_data(Average_Type);
      Savg_old.setVal(0.);

      NavierStokesBase::dt_avg[level]   = 0;
      NavierStokesBase::time_avg[level] = 0;
      NavierStokesBase::time_avg_fluct[level] = 0;


    }else{
      //
      // If Average_Type data were found, this means that we need to recover the
      // value of time_average
      //
      std::string line;
      std::string file=parent->theRestartFile();

      std::string File(file + "/TimeAverage");
      Vector<char> fileCharPtr;
      ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
      std::string fileCharPtrString(fileCharPtr.dataPtr());
      std::istringstream isp(fileCharPtrString, std::istringstream::in);

      // read in title line
      std::getline(isp, line);

      isp >> NavierStokesBase::time_avg[level];
      isp >> NavierStokesBase::time_avg_fluct[level];
      NavierStokesBase::dt_avg[level]   = 0;

    }
  }

#ifdef AMREX_USE_TURBULENT_FORCING
  //
  // Initialize data structures used for homogeneous isentropic forced turbulence.
  // Only need to do it once.
  if (level == 0)
      TurbulentForcing::init_turbulent_forcing(geom.ProbLoArray(),geom.ProbHiArray());
#endif

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
        u_mac = nullptr;
    }

    if (do_reflux && level < finest_level)
        reflux();

    //
    // Average everything down, including P and Gradp.
    // Even though the multilevel projections average down, only
    // single level projections have been done for current timestep.
    // The linearity of the average ensures that if we average down P
    // and Gp here, then we may simply add the incremental correction
    // (which get averaged down in amrex) during the sync projection.
    //
    // avgDown also updates rho_ctime since it's needed for rho_half,
    // which is used in the sync projection.
    //
    if (level < finest_level)
        avgDown();

    if (do_mac_proj && level < finest_level)
        mac_sync();

    // set Density to fluid_rho on all regions 
    if (do_diffused_ib) {
        MultiFab& S_new = get_new_data(State_Type);
        S_new.setVal(fluid_rho, Density, 1, S_new.nGrow());
        if (level < parent->finestLevel()) {
            auto&   fine_lev = getLevel(level+1);
            MultiFab& S_fine = fine_lev.get_new_data(State_Type);
            S_fine.setVal(fluid_rho, Density, 1, S_fine.nGrow());
        }
    }

    if (do_sync_proj && (level < finest_level))
        level_sync(crse_iteration);


    //
    // Test for conservation.
    //
    if (level==0 && sum_interval>0 && (parent->levelSteps(0)%sum_interval == 0))
    {
        sum_integrated_quantities();
    }

    if (level > 0) incrPAvg();

    // Copy pvf to Tracer before writing Tracer into plt files on the finest level,
    // or set Tracer to zero on the coarser/coarsest levels
    if (do_diffused_ib) {
        MultiFab& S_new = get_new_data(State_Type);
        if (level == parent->finestLevel()) {
            MultiFab::Copy(S_new, pvf, 0, Tracer, 1, pvf.nGrow()); // Note: the ghost cell region of pvf is zero. 
        }
        else {
            S_new.setVal(0.0, Tracer, 1, S_new.nGrow());
            // We need to average the Tracer here since we use it for refinement/de-refinement
            auto&   fine_lev = getLevel(level+1);
            MultiFab& S_fine = fine_lev.get_new_data(State_Type);
            average_down(S_fine, S_new, Tracer, 1);
        }
    }

    if (level == 0 && dump_plane >= 0)
    {
        Box bx = geom.Domain();

        AMREX_ASSERT(bx.bigEnd(AMREX_SPACEDIM-1) >= dump_plane);

        bx.setSmall(AMREX_SPACEDIM-1, dump_plane);
        bx.setBig  (AMREX_SPACEDIM-1, dump_plane);

        BoxArray ba(bx);
        DistributionMapping dm{ba};

        MultiFab mf(ba, dm, AMREX_SPACEDIM, 0, MFInfo(), Factory());

        mf.ParallelCopy(get_new_data(State_Type), Xvel, 0, AMREX_SPACEDIM);

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

    if (avg_interval > 0)
    {
      const amrex::Real dt_level = parent->dtLevel(level);
      time_average(time_avg[level], time_avg_fluct[level], dt_avg[level], dt_level);
    }

}

//
// Reset the time levels to time (time) and timestep dt.
// This is done at the end of the timestep in the pressure iteration section.
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

    //
    // Set P & gradP old = new. This way we retain new after
    // advance_setup() does swap(old,new).
    //
    initOldFromNew(Press_Type);
    state[Press_Type].setTimeLevel(time-dt_old,dt_old,dt_new);

    initOldFromNew(Gradp_Type);
    state[Gradp_Type].setTimeLevel(time-dt_old,dt_old,dt_new);
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

//
// Old checkpoint files may not have Gradp_Type and/or Average_Type.
//
void
NavierStokesBase::set_state_in_checkpoint (Vector<int>& state_in_checkpoint)
{
  //
  // Abort if any of the NSB::*_in_checkpoint variables haven't been set by user.
  //
  if ( gradp_in_checkpoint<0 || average_in_checkpoint<0 )
    Abort("\n\n   Checkpoint file is missing one or more state types. Set both\n ns.gradp_in_checkpoint and ns.avg_in_checkpoint to identify missing\n data. Set to 1 if present in checkpoint, 0 if not present. If unsure,\n try setting both to 0.\n\n If you just activated Time Averaging, you should add \n  ns.avg_in_checkpoint=0 ns.gradp_in_checkpoint=1 \n\n");

  //
  // Tell AmrLevel which types are in the checkpoint, so it knows what to copy.
  // state_in_checkpoint is initialized to all true.
  //
  if ( gradp_in_checkpoint==0 )
    state_in_checkpoint[Gradp_Type] = 0;

  if ( average_in_checkpoint==0 && avg_interval>0 )
    state_in_checkpoint[Average_Type] = 0;
}

void
NavierStokesBase::restart (Amr&          papa,
                           std::istream& is,
                           bool          bReadSpecial)
{
    Print()<<"\nWARNING! Note that you can't drop data from the checkpoint file.\n"
           <<" If your checkpoint file contains Average_Type, then your inputs\n"
           <<" must also specify ns.avg_interval>0.\n"<<std::endl;

    AmrLevel::restart(papa,is,bReadSpecial);

    if ( gradp_in_checkpoint==0 )
    {
      Print()<<"WARNING! GradP not found in checkpoint file. Recomputing from Pressure."
             <<std::endl;

      //
      // Compute GradP from the Pressure
      //
      computeGradP(state[Press_Type].curTime());
      computeGradP(state[Press_Type].prevTime());
    }

#ifdef AMREX_PARTICLES
    if(level == Particles::ParticleFinestLevel())
    {
        Particles::Restart(gravity, geom.CellSizeArray()[0],parent->levelSteps(0));
        ParallelDescriptor::Barrier();
    }
#endif

    define_workspace();
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
    // Do rho separately, as rho does not have forcing terms and is always conservative.
    //
    int sComp = first_scalar;

    if (sComp == Density)
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(S_old,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box&  bx = mfi.tilebox();
            const auto& Snew = S_new[mfi].array(Density);
            const auto& Sold = S_old[mfi].const_array(Density);
            const auto& advc = Aofs[mfi].const_array(Density);

            amrex::ParallelFor(bx, [ Snew, Sold, advc, dt]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Snew(i,j,k) = Sold(i,j,k) - dt * advc(i,j,k);
            });
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

            ConservativeScalMinMax(S_new, index_new_s, index_new_rho,
                                   Smf,   index_old_s, index_old_rho);

        }

        ++sComp;
    }

    //
    // Advective update of other scalars
    //
    if (sComp <= last_scalar)
    {
        //
        // Average mac face velocity to cell-centers for use in generating external
        // forcing term in getForce()
        //
        // NOTE that default getForce() does not use Vel or Scal, user must supply the
        // forcing function for that case.
        //
        MultiFab Vel(grids, dmap, AMREX_SPACEDIM, 0, MFInfo(), Factory());
#ifdef AMREX_USE_EB
        // This isn't quite right because it's face-centers to cell-centers
        // what's really wanted is face-centroid to cell-centroid
        EB_average_face_to_cellcenter(Vel, 0, Array<MultiFab const*,AMREX_SPACEDIM>{{AMREX_D_DECL(&u_mac[0],&u_mac[1],&u_mac[2])}});
#else
        average_face_to_cellcenter(Vel, 0, Array<MultiFab const*,AMREX_SPACEDIM>{{AMREX_D_DECL(&u_mac[0],&u_mac[1],&u_mac[2])}});
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            FArrayBox  tforces;

            for (MFIter mfi(S_old,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                if (getForceVerbose) {
                    amrex::Print() << "---" << '\n' << "E - scalar advection update (half time):" << '\n';
                }

                //
                // Create an estimate for time n+1 using advection term only; it's the best
                // we've got at this point. Then, average the new and old time to get
                // Crank-Nicholson half time approximation.
                //
                const Box& bx = mfi.tilebox();
                // getForce() expects all scalars to be present in Scal, and may need them
                // (particularly density) to make forcing term.
                FArrayBox Scal(bx,NUM_SCALARS);
                // Scal protected from early destruction by Gpu::synchronize at end of loop.
                const auto& Sn   = S_old[mfi].const_array(Density);
                const auto& Sarr = Scal.array();
                const auto& aofs_dens = Aofs[mfi].const_array(Density);
                // Create a local copy for lambda capture
                int numscal = NUM_SCALARS;

                amrex::ParallelFor(bx, [ Sn, Sarr, aofs_dens, dt, numscal]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    int n = 0;
                    // For density, we can create the Crank-Nicholson half-time approximation:
                    // Snew = Sold - dt*adv
                    // Shalftime = Sarr = (Snew + Sold)/2
                    Sarr(i,j,k,n) = Sn(i,j,k,n) - 0.5 * dt * aofs_dens(i,j,k,n);

                    // For other scalars, which may have diffusive or forcing terms, this is
                    // a safe choice.
                    for ( n = 1; n < numscal; n++ )
                    {
                        Sarr(i,j,k,n) = Sn(i,j,k,n);
                    }
                });

                const Real halftime = 0.5 * ( state[State_Type].curTime() +
                                              state[State_Type].prevTime() );
                FArrayBox& Vel_fab = Vel[mfi];

                int num_comp = last_scalar - sComp + 1;
                // Note that in general, num_comp != NUM_SCALAR
                tforces.resize(bx,num_comp);
                // tforces protected from early destruction by Gpu::synchronize at end of loop.
                getForce(tforces,bx,sComp,num_comp,halftime,Vel_fab,Scal,0,mfi);

                const auto& Snew = S_new[mfi].array(sComp);
                const auto& Sold = S_old[mfi].const_array(sComp);
                const auto& advc = Aofs[mfi].const_array(sComp);
                const auto& tf   = tforces.const_array();
                const auto& rho  = Scal.const_array();

                // Advection type
                amrex::Vector<int> iconserv_h;
                iconserv_h.resize(num_comp);
                for (int i = 0; i < num_comp; ++i) {
                    iconserv_h[i] = (advectionType[sComp+i] == Conservative) ? 1 : 0;
                }
                amrex::Gpu::DeviceVector<int> iconserv_d;
                iconserv_d.resize(num_comp);
                Gpu::copy(Gpu::hostToDevice, iconserv_h.begin(), iconserv_h.end(), iconserv_d.begin());
                const int* iconserv = iconserv_d.data();

                // Recall tforces is always density-weighted
                amrex::ParallelFor(bx, num_comp, [ Snew, Sold, advc, tf, dt, rho, iconserv]
                AMREX_GPU_DEVICE (int i, int j, int k, int n ) noexcept
                {
                    if ( iconserv[n] == 1 ) {
                        Snew(i,j,k,n) = Sold(i,j,k,n) + dt * ( - advc(i,j,k,n) + tf(i,j,k,n)  );
                    }
                    else {
                        Snew(i,j,k,n) = Sold(i,j,k,n) + dt * ( - advc(i,j,k,n) + tf(i,j,k,n)/rho(i,j,k) );
                    }
                });

                // Either need this synchronize here, or elixirs. Not sure if it matters which
                amrex::Gpu::synchronize();
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

        for (int sigma = sComp; sigma <= last_scalar; sigma++)
        {
            const int index_new_s   = sigma;
            const int index_new_rho = Density;
            const int index_old_s   = index_new_s   - Density;
            const int index_old_rho = index_new_rho - Density;

            if (advectionType[sigma] == Conservative)
            {
                ConservativeScalMinMax(S_new, index_new_s, index_new_rho,
                                       Smf,   index_old_s, index_old_rho);
            }
            else if (advectionType[sigma] == NonConservative)
            {
                ConvectiveScalMinMax(S_new, index_new_s, Smf, index_old_s);
            }
        }
    }

    //
    // Check the max of Snew and for NANs in solution
    //
    // static int count=0; count++;
    // VisMF::Write(S_new,"sn_"+std::to_string(count));
    // for (int sigma = first_scalar; sigma <= last_scalar; sigma++)
    // {
    //     std::cout << count <<" , comp = " << sigma << ", max(S_new)  = "
    //               << S_new.norm0( sigma, 0, false, true )
    //               << std::endl;
    // }

    // for (int sigma = first_scalar; sigma <= last_scalar; sigma++)
    // {
    //     if (S_old.contains_nan(sigma,1,0))
    //     {
    //         amrex::Print() << "SAU: Old scalar " << sigma << " contains Nans" << std::endl;

    //         IntVect mpt(AMREX_D_DECL(-100,100,-100));
    //         for (MFIter mfi(S_old); mfi.isValid(); ++mfi){
    //             if ( S_old[mfi].contains_nan<RunOn::Device>(mpt) )
    //                 amrex::Print() << " Nans at " << mpt << std::endl;
    //         }
    //     }
    //     if (S_new.contains_nan(sigma,1,0))
    //     {
    //         amrex::Print() << "SAU: New scalar " << sigma << " contains Nans" << std::endl;

    //         IntVect mpt(AMREX_D_DECL(-100,100,-100));
    //         for (MFIter mfi(S_new); mfi.isValid(); ++mfi){
    //             if ( S_new[mfi].contains_nan<RunOn::Device>(mpt) )
    //                 amrex::Print() << " Nans at " << mpt << std::endl;
    //         }
    //     }
    // }
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

    state[Press_Type].setTimeLevel(time-dt_old,dt_old,dt_old);

    state[Gradp_Type].setTimeLevel(time-dt_old,dt_old,dt_old);
}

void
NavierStokesBase::sync_setup (MultiFab*& DeltaSsync)
{
    AMREX_ASSERT(DeltaSsync == nullptr);

    int nconserved = 0;

    for (int comp = AMREX_SPACEDIM; comp < NUM_STATE; ++comp)
    {
        if (advectionType[comp] == Conservative)
            ++nconserved;
    }

    if (nconserved > 0 && level < parent->finestLevel())
    {
        DeltaSsync = new MultiFab(grids, dmap, nconserved, 0, MFInfo(), Factory());
        DeltaSsync->setVal(0);
    }
}

void
NavierStokesBase::sync_cleanup (MultiFab*& DeltaSsync)
{
    delete DeltaSsync;
    DeltaSsync = nullptr;
}

//
// Helper function for NavierStokesBase::SyncInterp().
//
static
void
set_bcrec_new (Vector<BCRec>  &bcrec,
               int             ncomp,
               int             src_comp,
               const Box&      box,
               const Box&      domain,
               const BoxArray& cgrids,
               int**           bc_orig_qty)

{
   for (int n = 0; n < ncomp; n++) {
      for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
      {
         int bc_index = (src_comp+n)*(2*AMREX_SPACEDIM) + dir;
         bcrec[n].setLo(dir,BCType::int_dir);
         bcrec[n].setHi(dir,BCType::int_dir);
         if ( ( box.smallEnd(dir) < domain.smallEnd(dir) ) ||
              ( box.bigEnd(dir)   > domain.bigEnd(dir) ) ) {
            for (int crse = 0; crse < cgrids.size(); crse++) {
               const Box& crsebx = cgrids[crse];
               if ( ( box.smallEnd(dir) < domain.smallEnd(dir) ) && ( crsebx.smallEnd(dir) == domain.smallEnd(dir) ) ) {
                  bcrec[n].setLo(dir,bc_orig_qty[crse][bc_index]);
               }
               if ( ( box.bigEnd(dir) > domain.bigEnd(dir) ) && ( crsebx.bigEnd(dir) == domain.bigEnd(dir) ) ) {
                  bcrec[n].setHi(dir,bc_orig_qty[crse][bc_index+AMREX_SPACEDIM]);
               }
            }
         }
      }
   }
}

//
// Interpolate A cell centered Sync correction from a
// coarse level (c_lev) to a fine level (f_lev).
//
// This routine interpolates the num_comp components of CrseSync
// (starting at src_comp) and either increments or puts the result into
// the num_comp components of FineSync (starting at dest_comp)
// The components of bc_orig_qty correspond to the quantities of CrseSync.
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

    Interpolater* interpolater = nullptr;

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

    NavierStokesBase& fine_level     = getLevel(f_lev);
    const BoxArray& fgrids           = fine_level.boxArray();
    const DistributionMapping& fdmap = fine_level.DistributionMap();
    const Geometry& fgeom            = parent->Geom(f_lev);
    const BoxArray& cgrids           = getLevel(c_lev).boxArray();
    const Geometry& cgeom            = parent->Geom(c_lev);
    Box             cdomain          = amrex::coarsen(fgeom.Domain(),ratio);
    const auto      N                = int(fgrids.size());

    BoxArray cdataBA(N);

    for (int i = 0; i < N; i++) {
        cdataBA.set(i,interpolater->CoarseBox(fgrids[i],ratio));
    }
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

    cdataMF.ParallelCopy(CrseSync, src_comp, 0, num_comp, cgeom.periodicity());

    //
    // Set physical boundary conditions in cdataMF.
    //
    //////////
    // Should be fine for EB for now, since EB doesn't intersect Phys BC
    // Not sure about what happens if EB intersects Phys BC
    ///////

    // tiling may not be needed here, but what the hey
    GpuBndryFuncFab<HomExtDirFill> gpu_bndry_func(HomExtDirFill{});
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(cdataMF,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       const Box& bx   = mfi.tilebox();
       FArrayBox& data = cdataMF[mfi];

       Vector<BCRec> bx_bcrec(num_comp);
       set_bcrec_new(bx_bcrec,num_comp,src_comp,bx,cdomain,cgrids,bc_orig_qty);
       gpu_bndry_func(bx,data,0,num_comp,cgeom,0.0,bx_bcrec,0,0);
    }

    //
    // Interpolate from cdataMF to fdata and update FineSync.
    // Note that FineSync and cdataMF will have the same distribution
    // since the length of their BoxArrays are equal.
    //
    MultiFab* fine_stateMF = nullptr;
    if (interpolater == &protected_interp)
    {
        fine_stateMF = &(getLevel(f_lev).get_new_data(State_Type));
    }


#ifdef AMREX_USE_EB
    const FabArray<EBCellFlagFab>& flags = dynamic_cast<EBFArrayBoxFactory const&>(getLevel(f_lev).Factory()).getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      for (MFIter mfi(FineSync,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         FArrayBox& cdata = cdataMF[mfi];
         const Box&  bx   = mfi.tilebox();
         const Box cbx    = interpolater->CoarseBox(bx,ratio);

#ifdef AMREX_USE_EB
         EBFArrayBox fdata(flags[mfi],bx,num_comp,FineSync[mfi].arena());
#else
         FArrayBox fdata(bx, num_comp);
#endif
         Elixir fdata_i = fdata.elixir();

         //
         // Set the boundary condition array for interpolation.
         //
         Vector<BCRec> bx_bcrec(num_comp);
         set_bcrec_new(bx_bcrec,num_comp,src_comp,cbx,cdomain,cgrids,bc_orig_qty);

         //ScaleCrseSyncInterp(cdata, c_lev, num_comp);

         interpolater->interp(cdata,0,fdata,0,num_comp,bx,ratio,
                              cgeom,fgeom,bx_bcrec,src_comp,State_Type,RunOn::Gpu);

         //reScaleFineSyncInterp(fdata, f_lev, num_comp);

         if (increment)
         {
            auto const& finedata    = fdata.array();
            auto const& coarsedata  = cdata.array();
            int scale_coarse = (interpolater == &protected_interp) ? 1 : 0;
            amrex::ParallelFor(bx, num_comp, [finedata,dt_clev]
            AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
            {
               finedata(i,j,k,n) *= dt_clev;
            });

            if ( scale_coarse ) {
                amrex::ParallelFor(cbx, num_comp, [coarsedata,dt_clev]
                AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
                {
                    coarsedata(i,j,k,n) *= dt_clev;
                });
            }


            if (interpolater == &protected_interp)
            {
               FArrayBox& fine_state = (*fine_stateMF)[mfi];
               interpolater->protect(cdata,0,fdata,0,fine_state,state_comp,
                                     num_comp,bx,ratio,
                                     cgeom,fgeom,bx_bcrec,RunOn::Gpu);
            }

            auto const& fsync       = FineSync.array(mfi,dest_comp);
            amrex::ParallelFor(bx, num_comp, [finedata,fsync]
            AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
            {
               fsync(i,j,k,n) += finedata(i,j,k,n);
            });

            if ( scale_coarse ) {
                amrex::ParallelFor(cbx, num_comp, [coarsedata,dt_clev]
                AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
                {
                    coarsedata(i,j,k,n) /= dt_clev;
                });
            }


         }
         else
         {
            auto const& finedata    = fdata.array();
            auto const& fsync       = FineSync.array(mfi,dest_comp);
            amrex::ParallelFor(bx, num_comp, [finedata,fsync]
            AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
            {
               fsync(i,j,k,n) = finedata(i,j,k,n);
            });
         }
       }
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
                                  IntVect&  ratio)
{
    BL_PROFILE("NavierStokesBase:::SyncProjInterp()");

    const BoxArray& P_grids = P_new.boxArray();
    const auto      N       = int(P_grids.size());

    BoxArray crse_ba(N);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        crse_ba.set(i,node_bilinear_interp.CoarseBox(P_grids[i],ratio));

    // None  of these 3 are actually used by node_bilinear_interp()
    Vector<BCRec> bc(AMREX_SPACEDIM);
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
    crse_phi.ParallelCopy(phi,0,0,1);

#ifdef AMREX_USE_EB
    EB_set_covered(crse_phi,0.);
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(P_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box&  bx     = mfi.tilebox();
        FArrayBox fine_phi(bx,1);
        Elixir fine_phi_i = fine_phi.elixir();
        node_bilinear_interp.interp(crse_phi[mfi],0,fine_phi,0,1,
                                    fine_phi.box(),ratio,cgeom,fgeom,bc,
                                    0,Press_Type,RunOn::Gpu);

        auto const& f_phi    = fine_phi.array();
        auto const& p_new    = P_new.array(mfi);
        auto const& p_old    = P_old.array(mfi);
        amrex::ParallelFor(bx, [f_phi, p_old, p_new]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
          p_new(i,j,k) += f_phi(i,j,k);
          p_old(i,j,k) += f_phi(i,j,k);
        });
    }

#ifdef AMREX_USE_EB
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
            amrex::Print() << "... advect momenta\n";
        }
    }

    const Real  prev_time      = state[State_Type].prevTime();

    std::unique_ptr<MultiFab> divu_fp(getDivCond(nghost_force(),prev_time));

    MultiFab forcing_term( grids, dmap, AMREX_SPACEDIM, nghost_force(), MFInfo(),Factory());
    forcing_term.setVal(0.0);

    FillPatchIterator U_fpi(*this,forcing_term, nghost_state(),prev_time,State_Type,Xvel,AMREX_SPACEDIM);
    MultiFab& Umf=U_fpi.get_mf();

    //
    // S_term is the state we are solving for: either velocity or momentum
    //
    MultiFab raii;
    MultiFab* S_term;
    if (do_mom_diff) {
        raii.define(grids, dmap, AMREX_SPACEDIM, nghost_state(), MFInfo(), Factory());
        S_term = &raii;
    } else {
        S_term = &Umf;
    }

    if (do_mom_diff)
    {
        FillPatchIterator Rho_fpi(*this,forcing_term,nghost_state(),prev_time,State_Type,Density,1);
        MultiFab& Rmf=Rho_fpi.get_mf();

        for (MFIter U_mfi(Umf,TilingIfNotGPU()); U_mfi.isValid(); ++U_mfi)
        {
            auto const state_bx = U_mfi.growntilebox(nghost_state());

            auto const& dens = Rmf.const_array(U_mfi); //Previous time, nghost_state() grow cells filled
            auto const& vel  = Umf.const_array(U_mfi);
            auto const& st   = S_term->array(U_mfi);

            amrex::ParallelFor(state_bx, AMREX_SPACEDIM, [ dens, vel, st ]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            { st(i,j,k,n) = vel(i,j,k,n) * dens(i,j,k); });
        }
    }

    MultiFab& Gp = get_old_data(Gradp_Type);

    FillPatchIterator S_fpi(*this,forcing_term,nghost_force(),prev_time,State_Type,Density,NUM_SCALARS);
    MultiFab& Smf=S_fpi.get_mf();

    // Get divu to time n+1/2
    {
        std::unique_ptr<MultiFab> dsdt(getDsdt(nghost_force(),prev_time));
        MultiFab::Saxpy(*divu_fp, 0.5*dt, *dsdt, 0, 0, 1, nghost_force());
    }

    MultiFab visc_terms(grids,dmap,AMREX_SPACEDIM,nghost_force(),MFInfo(),Factory());
    if (be_cn_theta != 1.0)
        getViscTerms(visc_terms,Xvel,AMREX_SPACEDIM,prev_time);
    else
        visc_terms.setVal(0.0);

    //
    // ls related
    // may add some surface tension subroutines later
    //

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter U_mfi(Umf,TilingIfNotGPU()); U_mfi.isValid(); ++U_mfi)
    {

        auto const force_bx = U_mfi.growntilebox(nghost_force()); // Box for forcing term

        if (getForceVerbose)
        {
            amrex::Print() << "---" << '\n'
                           << "B - velocity advection:" << '\n'
                           << "Calling getForce..." << '\n';
        }
        //
        // ls related
        // may consider the surface tension in a new getForce later
        //
        getForce(forcing_term[U_mfi],force_bx,Xvel,AMREX_SPACEDIM,
                 prev_time,Umf[U_mfi],Smf[U_mfi],0,U_mfi);

        //
        // Compute the total forcing.
        //
        auto const& tf   = forcing_term.array(U_mfi,Xvel);
        auto const& visc = visc_terms.const_array(U_mfi,Xvel);
        auto const& gp   = Gp.const_array(U_mfi);
        auto const& rho  = Smf.const_array(U_mfi); //Previous time, nghost_force() grow cells filled

        bool is_convective = do_mom_diff ? false : true;
        amrex::ParallelFor(force_bx, AMREX_SPACEDIM, [ tf, visc, gp, rho, is_convective]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            tf(i,j,k,n) = ( tf(i,j,k,n) + visc(i,j,k,n) - gp(i,j,k,n) );
            if (is_convective)
                tf(i,j,k,n) /= rho(i,j,k);
        });
    }

    ComputeAofs( Xvel, AMREX_SPACEDIM, *S_term, 0, forcing_term, *divu_fp, true, dt );
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

    for (int sigma = 0; sigma < AMREX_SPACEDIM; sigma++)
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
    MultiFab&  Gp    = get_old_data(Gradp_Type);
    MultiFab& Rh = get_rho_half_time();

    MultiFab Vel(grids, dmap, AMREX_SPACEDIM, 0, MFInfo(), Factory());
    //
    // Average mac face velocity to cell-centers for use in generating external
    // forcing term in getForce()
    // NOTE that default getForce() does not use Vel or Scal, user must supply the
    // forcing function for that case.
    //
#ifdef AMREX_USE_EB
    // This isn't quite right because it's face-centers to cell-centers
    // what's really wanted is face-centroid to cell-centroid
    EB_average_face_to_cellcenter(Vel, 0, Array<MultiFab const*,AMREX_SPACEDIM>{{AMREX_D_DECL(&u_mac[0],&u_mac[1],&u_mac[2])}});
#else
    average_face_to_cellcenter(Vel, 0, Array<MultiFab const*,AMREX_SPACEDIM>{{AMREX_D_DECL(&u_mac[0],&u_mac[1],&u_mac[2])}});
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
{
    FArrayBox  tforces, ScalFAB;

    for (MFIter mfi(Rh,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        FArrayBox& VelFAB = Vel[mfi];
        ScalFAB.resize(bx,NUM_SCALARS);
        Elixir scal_i = ScalFAB.elixir();

        //
        // Need to do some funky half-time stuff.
        //
        if (getForceVerbose)
           amrex::Print() << "---" << '\n' << "F - velocity advection update (half time):" << '\n';
        //
        // Average the new and old time to get Crank-Nicholson half time approximation.
        // Scalars always get updated before velocity (see NavierStokes::advance), so
        // this is guaranteed to be good.
        //
        auto const& scal = ScalFAB.array();
        auto const& scal_o = U_old.array(mfi,Density);
        auto const& scal_n = U_new.array(mfi,Density);
        amrex::ParallelFor(bx, NUM_SCALARS, [scal, scal_o, scal_n]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            scal(i,j,k,n) = 0.5 * ( scal_o(i,j,k,n) + scal_n(i,j,k,n) );
        });

        const Real half_time = 0.5*(state[State_Type].prevTime()+state[State_Type].curTime());
        tforces.resize(bx,AMREX_SPACEDIM);
        Elixir tf_i = tforces.elixir();
        getForce(tforces,bx,Xvel,AMREX_SPACEDIM,half_time,VelFAB,ScalFAB,0,mfi);

        //
        // Do following only at initial iteration--per JBB.
        //
        if (initial_iter && is_diffusive[Xvel]) {
           auto const& force  = tforces.array();
           amrex::ParallelFor(bx, AMREX_SPACEDIM, [force]
           AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
               force(i,j,k,n) = 0.0;
           });
        }

        // Update velocity
        auto const& vel_old  = U_old.array(mfi);
        auto const& vel_new  = U_new.array(mfi);
        auto const& gradp    = Gp.array(mfi);
        auto const& force    = tforces.array();
        auto const& advec    = Aofs.array(mfi);
        auto const& rho_old  = U_old.array(mfi, Density);
        auto const& rho_new  = U_new.array(mfi, Density);
        auto const& rho_Half = Rh.array(mfi);
        int mom_diff = do_mom_diff;
        amrex::ParallelFor(bx, AMREX_SPACEDIM, [vel_old,vel_new,gradp,force,advec,rho_old,rho_new,rho_Half,mom_diff,dt]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real velold = vel_old(i,j,k,n);

            if ( mom_diff ) {
               velold *= rho_old(i,j,k);
               vel_new(i,j,k,n) = velold - dt * advec(i,j,k,n)
                                         + dt * force(i,j,k,n)
                                         - dt * gradp(i,j,k,n);

               vel_new(i,j,k,n) /= rho_new(i,j,k);
            }
            else
            {
                vel_new(i,j,k,n) = velold - dt * advec(i,j,k,n)
                                          + dt * force(i,j,k,n) / rho_Half(i,j,k)
                                          - dt * gradp(i,j,k,n) / rho_Half(i,j,k);
            }
        });
    }
}

    for (int sigma = 0; sigma < AMREX_SPACEDIM; sigma++)
    {
       if (U_old.contains_nan(sigma,1,0))
       {
         amrex::Print() << "VAU: Old velocity " << sigma << " contains Nans" << std::endl;

         IntVect mpt(AMREX_D_DECL(-100,100,-100));
         for (MFIter mfi(U_old); mfi.isValid(); ++mfi){
           const Box& bx = mfi.tilebox();
           if ( U_old[mfi].contains_nan<RunOn::Device>(bx, sigma, 1, mpt) )
             amrex::Print() << " Nans at " << mpt << std::endl;
         }
       }
       if (U_new.contains_nan(sigma,1,0))
       {
         amrex::Print() << "VAU: New velocity " << sigma << " contains Nans" << std::endl;

         IntVect mpt(AMREX_D_DECL(-100,100,-100));
         for (MFIter mfi(U_new); mfi.isValid(); ++mfi){
           const Box& bx = mfi.tilebox();
           if ( U_new[mfi].contains_nan<RunOn::Device>(bx, sigma, 1, mpt) )
             amrex::Print() << " Nans at " << mpt << std::endl;
         }
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

        int   ngrow = 0;
        MultiFab visc_terms(grids,dmap,AMREX_SPACEDIM,ngrow,MFInfo(),Factory());
        MultiFab    tforces(grids,dmap,AMREX_SPACEDIM,ngrow,MFInfo(),Factory());

        //
        // Get grad(p)
        //
        MultiFab& Gp = get_old_data(Gradp_Type);

        //
        // Compute additional forcing terms
        //
        tforces.setVal(0.0);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(tforces,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const auto& bx = mfi.tilebox();
                  auto& tforces_fab = tforces[mfi];
            if (getForceVerbose)
            {
                amrex::Print() << "---" << '\n'
                               << "G - initial velocity diffusion update:" << '\n'
                               << "Calling getForce..." << '\n';
            }
            getForce(tforces_fab,bx,Xvel,AMREX_SPACEDIM,prev_time,U_old[mfi],U_old[mfi],Density,mfi);
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
          visc_terms.setVal(0.0);
        }

        //
        // Assemble RHS
        //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(tforces,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
           const Box& bx = mfi.tilebox();
           auto const& force   = tforces.array(mfi);
           auto const& viscT   = visc_terms.array(mfi);
           auto const& gradp   = Gp.array(mfi);
           auto const& rhohalf = Rh.array(mfi);
           auto const& rho_old = U_old.array(mfi,Density);
           auto const& rho_new = U_new.array(mfi,Density);
           auto const& vel_old = U_old.array(mfi,Xvel);
           auto const& vel_new = U_new.array(mfi,Xvel);
           auto const& advT    = aofs->array(mfi,Xvel);
           int mom_diff = do_mom_diff;
           amrex::ParallelFor(bx, AMREX_SPACEDIM, [force,viscT,gradp,rhohalf,advT,rho_old,rho_new,vel_old,vel_new,mom_diff,dt]
           AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
              // Set force += (visc - Gp) / rho_half - aofs
              force(i,j,k,n) += viscT(i,j,k,n) - gradp(i,j,k,n);
              if ( !mom_diff ) {
                 force(i,j,k,n) /= rhohalf(i,j,k);
              }
              force(i,j,k,n) -= advT(i,j,k,n);
              // if mom_diff : U_new = (force* dt + U_old * rho_old) / rho_new
              // else        : U_new = U_old + force* dt
              if ( mom_diff ) {
                 vel_new(i,j,k,n) = (force(i,j,k,n) * dt + vel_old(i,j,k,n) * rho_old(i,j,k)) / rho_new(i,j,k);
              } else {
                 vel_new(i,j,k,n) = vel_old(i,j,k,n) + force(i,j,k,n) * dt;
              }
           });
        }
    }
}

#ifdef AMREX_PARTICLES

void
NavierStokesBase::read_particle_params ()
{
    ParmParse ppp("particles");
    //
    // Ensure other particle methods aren't being used, like sprays
    //
    ppp.query("do_nspc_particles", do_nspc);
    if (!do_nspc) return;
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

    ppp.query("verbose",pverbose);
    if ( ppp.countname("pverbose") > 0) {
        amrex::Abort("particles.pverbose found in inputs. Please use particles.verbose");
    }
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
    //
    // Put particle info in plotfile (using ParticleContainer::Checkpoint)?
    //
    ppp.query("particles_in_plotfile", particles_in_plotfile);
}

void
NavierStokesBase::initParticleData ()
{

    if (!do_nspc) { return; }

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
    if (level == 0 && do_nspc)
    {
        AMREX_ASSERT(NSPC == 0);

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
                          tfab.copy<RunOn::Device>(sfab, box, timestamp_indices[i], box, i, 1);
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
    if ((name == "particle_count" || name == "total_particle_count") && do_nspc) {
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

            IntVect trr(AMREX_D_DECL(1,1,1));

            for (int lev = level+1; lev <= parent->finestLevel(); lev++)
            {
                BoxArray ba = parent->boxArray(lev);

                MultiFab temp_dat(ba,parent->DistributionMap(lev),1,0,MFInfo(),Factory());

                trr *= parent->refRatio(lev-1);

                ba.coarsen(trr);

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

                    AMREX_ASSERT(cfab.box() == amrex::coarsen(fbx,trr));

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
                dat.ParallelCopy(ctemp_dat);

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
    Vector<int> bc(2*AMREX_SPACEDIM*ncomp);
    BCRec bcr;
    const StateDescriptor* stDesc;
    const Box& domain = geom.Domain();

    for (int n = 0; n < ncomp; n++)
    {
        stDesc=state[State_Type].descriptor();
        setBC(bx,domain,stDesc->getBC(scomp+n),bcr);

        const int* b_rec = bcr.vect();
        for (int m = 0; m < 2*AMREX_SPACEDIM; m++) {
            bc[2*AMREX_SPACEDIM*n + m] = b_rec[m];
        }
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

//
// Compute gradient of P and fill ghost cells with FillPatch
//
void
NavierStokesBase::computeGradP(Real time)
{
    LPInfo info;
    info.setMaxCoarseningLevel(0);
    MLNodeLaplacian linop({geom}, {grids}, {dmap}, info, {&Factory()});
#ifdef AMREX_USE_EB
    linop.buildIntegral();
#endif

    // No call to set BCs because we're only calling compGrad(), which
    // doesn't use them. P already exists on surroundingNodes(Gp.validbox()),
    // and compGrad() does not fill ghost cells

    MultiFab& Press = get_data(Press_Type, time);
    MultiFab& Gp    = get_data(Gradp_Type, time);

    linop.compGrad(0, Gp, Press);

    // Now fill ghost cells
    FillPatch(*this,Gp,Gp.nGrow(),time,Gradp_Type,0,AMREX_SPACEDIM);
}

void
NavierStokesBase::avgDown_StatePress()
{
    auto&   fine_lev = getLevel(level+1);

    //
    // Average down the states at the new time.
    //
    MultiFab& S_crse = get_new_data(State_Type);
    MultiFab& S_fine = fine_lev.get_new_data(State_Type);

    average_down(S_fine, S_crse, 0, S_crse.nComp());

    //
    // Fill rho_ctime at the current and finer levels with the correct data.
    //
    for (int lev = level; lev <= parent->finestLevel(); lev++)
    {
        getLevel(lev).make_rho_curr_time();
    }

    //
    // Now average down pressure over time n-(n+1) interval.
    //
    MultiFab&       P_crse      = get_new_data(Press_Type);
    MultiFab&       P_fine_init = fine_lev.get_new_data(Press_Type);
    MultiFab&       P_fine_avg  = fine_lev.p_avg;
    MultiFab&       P_fine      = initial_step ? P_fine_init : P_fine_avg;

    // NOTE: this fills ghost cells, but amrex::average_down does not.
    amrex::average_down_nodal(P_fine,P_crse,fine_ratio);

    //
    // Average down Gradp
    //
    MultiFab& Gp_crse = get_new_data(Gradp_Type);
    MultiFab& Gp_fine = fine_lev.get_new_data(Gradp_Type);

    average_down(Gp_fine, Gp_crse, 0, Gp_crse.nComp());
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
#if (AMREX_SPACEDIM == 3)
    // Don't need volume weighting for EB -- EB is only used in Cartesian
    amrex::EB_average_down(S_fine, S_crse, scomp, ncomp, fine_ratio);
#else
    // Volume weighting
    amrex::EB_average_down(S_fine, S_crse, this->getLevel(level+1).Volume(),
                           *(this->getLevel(level+1).VolFrac()),
                           scomp, ncomp, fine_ratio);
#endif

#else
    //
    // non-EB aware, uses volume weighting for 1D,2D but no volume weighting for 3D
    //
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
    MultiFab& Gp = new_data? get_new_data(Gradp_Type) : get_old_data(Gradp_Type);
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


//
// Correct a conservatively-advected scalar for under-over shoots.
//
void
NavierStokesBase::ConservativeScalMinMax ( amrex::MultiFab&       Snew, const int snew_comp, const int new_density_comp,
                                           amrex::MultiFab const& Sold, const int sold_comp, const int old_density_comp )
{
    amrex::ignore_unused(this);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(Snew,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        const auto& bx = mfi.tilebox();

        const auto& sn   = Snew.array(mfi,snew_comp);
        const auto& so   = Sold.const_array(mfi,sold_comp);
        const auto& rhon = Snew.const_array(mfi,new_density_comp);
        const auto& rhoo = Sold.const_array(mfi,old_density_comp);
#ifdef AMREX_USE_EB
        const auto& ebfactory = dynamic_cast<amrex::EBFArrayBoxFactory const&>(Factory());
        const auto& vfrac = ebfactory.getVolFrac().const_array(mfi);
#endif

        amrex::ParallelFor(bx, [=]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real smn = std::numeric_limits<Real>::max();
            Real smx = std::numeric_limits<Real>::min();

#if (AMREX_SPACEDIM==3)
            int ks = -1;
            int ke =  1;
#else
            int ks = 0;
            int ke = 0;
#endif

            for (int kk = ks; kk <= ke; ++kk)
            {
                for (int jj = -1; jj <= 1; ++jj)
                {
                    for (int ii = -1; ii <= 1; ++ii)
                    {
#ifdef AMREX_USE_EB
                        if ( vfrac (i+ii,j+jj,k+kk) > 0. )
#endif
                        {
                            smn =  amrex::min(smn, so(i+ii,j+jj,k+kk)/rhoo(i+ii,j+jj,k+kk));
                            smx =  amrex::max(smx, so(i+ii,j+jj,k+kk)/rhoo(i+ii,j+jj,k+kk));
                        }
                    }
                }
            }
            sn(i,j,k) = amrex::min( amrex::max(sn(i,j,k)/rhon(i,j,k), smn), smx ) * rhon(i,j,k);
        });
    }
}

//
// Correct a convectively-advected  scalar for under-over shoots.
//
void
NavierStokesBase::ConvectiveScalMinMax ( amrex::MultiFab&       Snew, const int snew_comp,
                                         amrex::MultiFab const& Sold, const int sold_comp )
{
    amrex::ignore_unused(this);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(Snew,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        const auto& bx = mfi.tilebox();

        const auto& sn   = Snew.array(mfi,snew_comp);
        const auto& so   = Sold.const_array(mfi,sold_comp);
#ifdef AMREX_USE_EB
        const auto& ebfactory = dynamic_cast<amrex::EBFArrayBoxFactory const&>(Factory());
        const auto& vfrac = ebfactory.getVolFrac().const_array(mfi);
#endif

        amrex::ParallelFor(bx, [=]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real smn = std::numeric_limits<Real>::max();
            Real smx = std::numeric_limits<Real>::min();

#if (AMREX_SPACEDIM==3)
            int ks = -1;
            int ke =  1;
#else
            int ks = 0;
            int ke = 0;
#endif

            for (int kk = ks; kk <= ke; ++kk)
            {
                for (int jj = -1; jj <= 1; ++jj)
                {
                    for (int ii = -1; ii <= 1; ++ii)
                    {
#ifdef AMREX_USE_EB
                        if ( vfrac (i+ii,j+jj,k+kk) > 0. )
#endif
                        {
                            smn =  amrex::min(smn, so(i+ii,j+jj,k+kk));
                            smx =  amrex::max(smx, so(i+ii,j+jj,k+kk));
                        }
                    }
                }
            }
            sn(i,j,k) = amrex::min( amrex::max(sn(i,j,k), smn), smx );
        });
    }
}


//
// Predict the edge velocities which go into forming u_mac.  This
// function also returns an estimate of dt for use in variable timesteping.
//
Real
NavierStokesBase::predict_velocity (Real  dt)
{
   BL_PROFILE("NavierStokesBase::predict_velocity()");
   if (verbose) {
      amrex::Print() << "... predict edge velocities\n";
   }
   //
   // Get simulation parameters.
   //
   const int   nComp          = AMREX_SPACEDIM;
   const Real* dx             = geom.CellSize();
   const Real  prev_time      = state[State_Type].prevTime();
   const Real  prev_pres_time = state[Press_Type].prevTime();
   const Real  strt_time      = ParallelDescriptor::second();

   //
   // Compute viscous terms at level n.
   // Ensure reasonable values in 1 grow cell.  Here, do extrap for
   // c-f/phys boundary, since we have no interpolator fn, also,
   // preserve extrap for corners at periodic/non-periodic intersections.
   //
   MultiFab visc_terms(grids,dmap,nComp,nghost_force(),MFInfo(), Factory());

   FillPatchIterator U_fpi(*this,visc_terms,nghost_state(),prev_time,State_Type,Xvel,AMREX_SPACEDIM);
   MultiFab& Umf=U_fpi.get_mf();

   // Floor small values of states to be extrapolated
   floor(Umf);

   //
   // Compute "grid cfl number" based on cell-centered time-n velocities
   //
   auto umax = Umf.norm0({AMREX_D_DECL(0,1,2)},Umf.nGrow(), /*local = */false, /*ignore_covered = */true);
   Real cflmax = dt*umax[0]/dx[0];
   for (int d=1; d<AMREX_SPACEDIM; ++d) {
     cflmax = std::max(cflmax,dt*umax[d]/dx[d]);
   }
   Real tempdt = cflmax==0 ? change_max : std::min(change_max,cfl/cflmax);

   if ( advection_scheme == "Godunov_PLM" || advection_scheme == "Godunov_PPM" || advection_scheme == "BDS")
   {
       MultiFab& Gp = get_old_data(Gradp_Type);
       // FillPatch Gp here, as crse data has been updated
       // only really needed on the first step of the fine subcycle
       // because level_project at this level fills Gp ghost cells.
       // OR maybe better to not FP in level_proj for level>0 and do it
       // once the first time it's needed, which is presumably here...
       if ( level > 0 )
           FillPatch(*this,Gp,Gp.nGrow(),prev_pres_time,Gradp_Type,0,AMREX_SPACEDIM);

       if (be_cn_theta != 1.0)
       {
           getViscTerms(visc_terms,Xvel,nComp,prev_time);
       }
       else
       {
           visc_terms.setVal(0.0);
       }

       //
       // ls related
       // may add some surface tension subroutines later
       //

       FillPatchIterator S_fpi(*this,visc_terms,nghost_state(),prev_time,State_Type,Density,NUM_SCALARS);
       MultiFab& Smf=S_fpi.get_mf();

       MultiFab forcing_term( grids, dmap, AMREX_SPACEDIM, nghost_force() );

       //
       // Compute forcing
       //
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       {
           for (MFIter U_mfi(Umf,TilingIfNotGPU()); U_mfi.isValid(); ++U_mfi)
           {
               FArrayBox& Ufab = Umf[U_mfi];
               auto const  gbx = U_mfi.growntilebox(nghost_force());

               if (getForceVerbose) {
                   Print() << "---\nA - Predict velocity:\n Calling getForce...\n";
               }

               //
               // ls related
               // may consider the surface tension in a new getForce later
               //
               getForce(forcing_term[U_mfi],gbx,Xvel,AMREX_SPACEDIM,prev_time,Ufab,Smf[U_mfi],0,U_mfi);

               //
               // Compute the total forcing.
               //
               auto const& tf   = forcing_term.array(U_mfi,Xvel);
               auto const& visc = visc_terms.const_array(U_mfi,Xvel);
               auto const& gp   = Gp.const_array(U_mfi);
               auto const& rho  = Smf.const_array(U_mfi);

               amrex::ParallelFor(gbx, AMREX_SPACEDIM, [tf, visc, gp, rho]
               AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
               {
                   tf(i,j,k,n) = ( tf(i,j,k,n) + visc(i,j,k,n) - gp(i,j,k,n) ) / rho(i,j,k);
               });
           }
       }

#ifdef AMREX_USE_EB
       if (!EBFactory().isAllRegular())
       {
           EBGodunov::ExtrapVelToFaces( Umf, forcing_term,
                                        AMREX_D_DECL(u_mac[0], u_mac[1], u_mac[2]),
                                        m_bcrec_velocity, m_bcrec_velocity_d.dataPtr(),
                                        geom, dt );
       }
       else
#endif
       {
           bool godunov_use_ppm = ( advection_scheme == "Godunov_PPM" ? true : false );

           Godunov::ExtrapVelToFaces( Umf, forcing_term,
                                      AMREX_D_DECL(u_mac[0], u_mac[1], u_mac[2]),
                                      m_bcrec_velocity, m_bcrec_velocity_d.dataPtr(),
                                      geom, dt,
                                      godunov_use_ppm, godunov_use_forces_in_trans );
       }

   }
   else
   {
       Abort("NSB::predict_velocity: Unknown advection_scheme");
   }

   if (verbose > 1)
   {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

      ParallelDescriptor::ReduceRealMax(run_time,IOProc);

      Print() << "NavierStokesBase::predict_velocity(): lev: " << level
              << ", time: " << run_time << '\n';
   }

   return dt*tempdt;
}


//
// Floor small values of states to be extrapolated
//
void
NavierStokesBase::floor(MultiFab& mf){

  int ncomp = mf.nComp();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box gbx=mfi.growntilebox(mf.nGrow());
        auto const& fab_a = mf.array(mfi);
        AMREX_PARALLEL_FOR_4D ( gbx, ncomp, i, j, k, n,
        {
            auto& val = fab_a(i,j,k,n);
            val = amrex::Math::abs(val) > 1.e-20 ? val : 0;
        });
    }
}

int
NavierStokesBase::nghost_state () const
{
    amrex::ignore_unused(this);
#ifdef AMREX_USE_EB
  if (!EBFactory().isAllRegular())
  {
        return 4;
  }
  else
#endif
  {
    return 3;
  }
}

void
NavierStokesBase::ComputeAofs ( int comp, int ncomp,
                                MultiFab const& S,
                                int S_comp,
                                MultiFab const& forcing_term,
                                MultiFab const& divu,
                                bool is_velocity, Real dt)
{
    Array<MultiFab, AMREX_SPACEDIM> cfluxes;
    Array<MultiFab, AMREX_SPACEDIM> edgestate;

    //
    // Advection needs S to have 2-3 ghost cells.
    // Advection routines call slopes on cells i & i+1, and then
    // 2nd order slopes use i+/-1 => S needs 2 ghost cells (MOL)
    // 4th order slopes use i+/-2 => S needs 3 ghost cells (Godunov)
    //
    int nghost = 0;
    for (int i = 0; i < AMREX_SPACEDIM; ++i)
    {
        const BoxArray& ba = getEdgeBoxArray(i);
        cfluxes[i].define(ba, dmap, ncomp, nghost, MFInfo(), Factory());
        edgestate[i].define(ba, dmap, ncomp, nghost, MFInfo(), Factory());
    }

    bool do_crse_add = true;
    bool do_fine_add = true;

    ComputeAofs(*aofs, /*aofs_comp*/ comp, /*state_indx*/ comp, ncomp,
                S, S_comp,
                &forcing_term, /*forcing_term_comp*/ 0,
                &divu,
                cfluxes,  /*flux_comp*/ 0,
                edgestate, /*edge_comp*/ 0, /*known_edgestate*/ false,
                is_velocity, dt,
                /*is_sync*/ false, /*sync fluxing velocity Ucorr*/ {},
                do_crse_add, do_fine_add);
}

void
NavierStokesBase::ComputeAofs ( MultiFab& advc, int a_comp, // Advection term "Aofs" held here
                                int state_indx, // Index of first component in AmrLevel.state corresponding to quantity to advect
                                int ncomp,
                                MultiFab const& S, int S_comp, // State for computing edgestates, may have gotten massaged
                                                               // a bit compared to AmrLevel.state. Must have filled ghost cells.
                                MultiFab const* forcing, int f_comp,
                                MultiFab const* divu, // Constraint divu=Source, not div(Umac)
                                Array<MultiFab, AMREX_SPACEDIM>& cfluxes, int flux_comp,
                                Array<MultiFab, AMREX_SPACEDIM>& edgestate, int edge_comp,
                                bool known_edge_state,
                                bool is_velocity, Real dt,
                                bool is_sync, Array<MultiFab*, AMREX_SPACEDIM> const& U_corr,
                                bool do_crse_add, bool do_fine_add)
{
    BL_PROFILE("NSB::ComputeAofs_kernel");

    amrex::ignore_unused(do_fine_add);

    // Need U_corr to be defined for sync.
    AMREX_ASSERT( (is_sync && !U_corr.empty()) || !is_sync );

    // Advection type conservative or non?
    // Maybe there's something better than DeviceVector now...
    amrex::Gpu::DeviceVector<int> iconserv;
    Vector<int> iconserv_h;
    iconserv.resize(ncomp, 0);
    iconserv_h.resize(ncomp, 0);

    // Will we do any convective differencing?
    bool any_convective = false;
    for (int i = 0; i < ncomp; ++i) {
        iconserv_h[i] = (advectionType[state_indx+i] == Conservative) ? 1 : 0;
        if (!iconserv_h[i]) any_convective = true;
    }
    Gpu::copy(Gpu::hostToDevice,iconserv_h.begin(),iconserv_h.end(), iconserv.begin());
    int const* iconserv_ptr = iconserv.data();

    // As code is currently written, ComputeAofs is always called separately for
    // velocity vs scalars. May be called with an individual scalar.
    auto const& bcrec_h = fetchBCArray(State_Type, state_indx, ncomp);
    auto* const bcrec_d = is_velocity ? m_bcrec_velocity_d.dataPtr()
                                     : &m_bcrec_scalars_d.dataPtr()[state_indx-AMREX_SPACEDIM];

#ifdef AMREX_USE_EB
    auto const& ebfact= dynamic_cast<EBFArrayBoxFactory const&>(Factory());

    // Always need a temporary MF to hold advective update before redistribution.
    MultiFab update_MF(advc.boxArray(),advc.DistributionMap(),ncomp,3,MFInfo(),Factory());

    // Must initialize to zero because not all values may be set, e.g. outside the domain.
    update_MF.setVal(0.);
#endif

    //
    // Define some parameters for hydro routines
    //

    bool fluxes_are_area_weighted = true;

    // AMReX_Hydro only recognizes Godunov scheme and then uses ppm switch.
    std::string scheme   = (advection_scheme=="Godunov_PLM" || advection_scheme=="Godunov_PPM")
                           ? "Godunov" : advection_scheme;
    bool godunov_use_ppm = (advection_scheme == "Godunov_PPM") ? true : false ;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(advc,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx   = mfi.tilebox();

        const auto& S_arr = S.const_array(mfi, S_comp);
        AMREX_D_TERM( const auto& fx = cfluxes[0].array(mfi,flux_comp);,
                      const auto& fy = cfluxes[1].array(mfi,flux_comp);,
                      const auto& fz = cfluxes[2].array(mfi,flux_comp););
        AMREX_D_TERM( const auto& xed = edgestate[0].array(mfi,edge_comp);,
                      const auto& yed = edgestate[1].array(mfi,edge_comp);,
                      const auto& zed = edgestate[2].array(mfi,edge_comp););
        AMREX_D_TERM( const auto& umac = u_mac[0].const_array(mfi);,
                      const auto& vmac = u_mac[1].const_array(mfi);,
                      const auto& wmac = u_mac[2].const_array(mfi););
        AMREX_D_TERM( const auto& uflux = (is_sync) ? U_corr[0]->const_array(mfi) : u_mac[0].const_array(mfi);,
                      const auto& vflux = (is_sync) ? U_corr[1]->const_array(mfi) : u_mac[1].const_array(mfi);,
                      const auto& wflux = (is_sync) ? U_corr[2]->const_array(mfi) : u_mac[2].const_array(mfi););

#ifdef AMREX_USE_EB
        const auto& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
        const auto   fabtyp = flagfab.getType(bx);

        if (fabtyp == FabType::covered)
        {
            //
            // Set advection term and move on to next iteration.
            // Old ComputeAofs also set fluxes (=0) and edgestates (=covered_val) here.
            //
            auto const& aofs_arr = advc.array(mfi, a_comp);
            amrex::ParallelFor(bx, ncomp, [aofs_arr]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            { aofs_arr( i, j, k, n ) = COVERED_VAL;});

            continue;
        }
#endif

        //
        // 1. Compute fluxes
        //
        HydroUtils::ComputeFluxesOnBoxFromState(bx, ncomp, mfi, S_arr,
                                                AMREX_D_DECL(fx, fy, fz),
                                                AMREX_D_DECL(xed, yed, zed),
                                                known_edge_state,
                                                AMREX_D_DECL(umac, vmac, wmac), //used to create edge state
                                                AMREX_D_DECL(uflux, vflux, wflux), //used to create flux
                                                (divu) ? divu->const_array(mfi) : Array4<Real const>{},
                                                (forcing) ? forcing->const_array(mfi, f_comp) : Array4<Real const>{},
                                                geom, dt,
                                                bcrec_h, bcrec_d, iconserv_ptr,
#ifdef AMREX_USE_EB
                                                ebfact,
                                                /*values_on_eb_inflow*/ Array4<Real const> {},
#endif
                                                godunov_use_ppm, godunov_use_forces_in_trans,
                                                is_velocity, fluxes_are_area_weighted,
                                                scheme);


        //
        // 2. Get the right container to hold the advective update
        //
        FArrayBox* update_fab;
        int update_comp = 0;

#ifdef AMREX_USE_EB
        // Recall that for EB we always use a temporary MF for redistribution.
        update_fab = &update_MF[mfi];
#else
        FArrayBox tmp;
        if (is_sync)
        {
            // For the sync, we need a temporary FAB to hold the update because
            // we add the update to what's already in advc.
            tmp.resize(bx, ncomp, The_Async_Arena());
            update_fab = &tmp;
        }
        else
        {
            // Otherwise, we can overwrite advc.
            update_fab = &advc[mfi];
            update_comp = a_comp;
        }
#endif
        auto const& update_arr = update_fab->array(update_comp);


        //
        // 3. Compute flux divergence
        //
        // We compute -div here for consistency with the way we HAVE to do it for EB
        // (because redistribution operates on -div rather than div)
        Real mult = -1.0;
#ifdef AMREX_USE_EB
        const auto& vfrac  = ebfact.getVolFrac().const_array(mfi);

        if (fabtyp != FabType::regular)
        {
            HydroUtils::EB_ComputeDivergence( bx, update_arr,
                                              AMREX_D_DECL( fx, fy, fz ),
                                              vfrac,
                                              ncomp, geom,
                                              mult, fluxes_are_area_weighted);
        }
        else
#endif
        {
            HydroUtils::ComputeDivergence( bx, update_arr,
                                           AMREX_D_DECL( fx, fy, fz ),
                                           ncomp, geom,
                                           mult, fluxes_are_area_weighted);
        }

        //
        // 4. Formulate convective term, if needed
        //
        if (any_convective && !is_sync) // Recall sync is always a conservative update
        {
            // Compute div(u_mac) first
            FArrayBox div_umac(bx, 1, The_Async_Arena());
            auto const& divum_arr = div_umac.array();

#ifdef AMREX_USE_EB
            const GpuArray<Real,AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();

            if (fabtyp != FabType::regular)
            {
                // Need AMReX routine here. Hydro version takes a flux (which always has area-fraction
                // already included).
                bool already_on_centroids = true;
                AMREX_D_TERM(Array4<Real const> const& apx = ebfact.getAreaFrac()[0]->const_array(mfi);,
                             Array4<Real const> const& apy = ebfact.getAreaFrac()[1]->const_array(mfi);,
                             Array4<Real const> const& apz = ebfact.getAreaFrac()[2]->const_array(mfi));
                Array4<EBCellFlag const> const& flagarr = flagfab.const_array();
                AMREX_HOST_DEVICE_FOR_4D(bx,div_umac.nComp(),i,j,k,n,
                {
                    eb_compute_divergence(i,j,k,n,divum_arr,AMREX_D_DECL(umac,vmac,wmac),
                                          Array4<int const>{}, flagarr, vfrac,
                                          AMREX_D_DECL(apx,apy,apz),
                                          AMREX_D_DECL(Array4<Real const>{},
                                                       Array4<Real const>{},
                                                       Array4<Real const>{}),
                                          dxinv, already_on_centroids);
                });
            }
            else
#endif
            {
                HydroUtils::ComputeDivergence(bx, divum_arr, AMREX_D_DECL(umac,vmac,wmac),
                                              1, geom, Real(1.0), false );
            }

            HydroUtils::ComputeConvectiveTerm(bx, ncomp, mfi, S_arr,
                                              AMREX_D_DECL( xed, yed, zed ),
                                              div_umac.array(), update_arr,
                                              iconserv_ptr,
#ifdef AMREX_USE_EB
                                              ebfact,
#endif
                                              scheme);
        }


#ifndef AMREX_USE_EB
        //
        // non-EB step 5:
        //
        // Recall we computed -div above, because redistribution operates
        // on -div. Thus, we use -update here.
        //
        auto const& aofs_arr = advc.array(mfi,a_comp);
        if (is_sync)
        {
            amrex::ParallelFor(bx, ncomp, [aofs_arr, update_arr]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            { aofs_arr( i, j, k, n ) -= update_arr(i,j,k,n); });
        }
        else
        {
            amrex::ParallelFor(bx, ncomp, [aofs_arr, update_arr]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            { aofs_arr( i, j, k, n ) = -update_arr(i,j,k,n); });
        }
#endif
    }
    //
    // The non-EB computation is complete.
    //

    // EB step 5: Redistribute the advective update stashed in update_MF.
    //
#ifdef AMREX_USE_EB
    update_MF.FillBoundary(geom.periodicity());

    //
    // Define the "state" for StateRedistribution.
    //
    MultiFab  rstate_tmp;
    if (is_sync && redistribution_type == "StateRedist")
    {
        // For the sync, use the Sync data passed in via advc as the "state".
        // WARNING: This choice may lead to oversmoothing.
        //
        // Must create a temporary copy so we're not overwriting this "state" as
        // we go through the redistribution process.
        rstate_tmp.define(S.boxArray(),S.DistributionMap(),ncomp,S.nGrow(),
                          MFInfo(),ebfact);
        MultiFab::Copy(rstate_tmp,advc,a_comp,0,ncomp,S.nGrow());
    }
    MultiFab const* rstate = (is_sync && redistribution_type == "StateRedist")
                             ? &rstate_tmp : &S;
    int        rstate_comp = (is_sync && redistribution_type == "StateRedist")
                             ? 0           : S_comp;
#endif

    const auto& dx    = geom.CellSizeArray();

    // **********************************************************************************
    // We use this hack to allow CrseAdd/FineAdd to divide area-weighted fluxes by volume
    //    instead of needing to un-area-weight the fluxes then divide just by dx
    // **********************************************************************************
    AMREX_ALWAYS_ASSERT(fluxes_are_area_weighted);
    Real dx1 = dx[0];
    for (int dir = 1; dir < AMREX_SPACEDIM; ++dir) {
      dx1 *= dx[dir];
    }

    std::array<Real, AMREX_SPACEDIM> dxD = {{AMREX_D_DECL(dx1, dx1, dx1)}};
    const Real* dxDp = &(dxD[0]);

    // **********************************************************************************
    // Build mask to find the ghost cells we need to correct
    // **********************************************************************************
    if (coarse_fine_mask == nullptr) {
       coarse_fine_mask = std::make_unique<iMultiFab>(grids, dmap, 1, 2, MFInfo(), DefaultFabFactory<IArrayBox>());
       coarse_fine_mask->BuildMask(geom.Domain(), geom.periodicity(),
                   level_mask_covered, level_mask_notcovered, level_mask_physbnd, level_mask_interior);
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    // for (MFIter mfi(advc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    for (MFIter mfi(advc, false); mfi.isValid(); ++mfi)
    {
        AMREX_D_TERM( const auto& fx_fab = (cfluxes[0])[mfi];,
                      const auto& fy_fab = (cfluxes[1])[mfi];,
                      const auto& fz_fab = (cfluxes[2])[mfi];);

#ifdef AMREX_USE_EB
        auto const& bx = mfi.tilebox();

        auto const& flagfab   = ebfact.getMultiEBCellFlagFab()[mfi];
        auto const& flags_arr = flagfab.const_array();

        if (flagfab.getType(bx) != FabType::covered )
        {
            auto const& aofs_arr = advc.array(mfi, a_comp);
            auto const& update_arr = update_MF.array(mfi);

            FArrayBox         dm_as_fine(Box::TheUnitBox(),ncomp);
            FArrayBox   fab_drho_as_crse(Box::TheUnitBox(),ncomp);
            IArrayBox fab_rrflag_as_crse(Box::TheUnitBox());

            if (flagfab.getType(grow(bx,4)) != FabType::regular)
            {
                AMREX_D_TERM( auto apx = ebfact.getAreaFrac()[0]->const_array(mfi);,
                              auto apy = ebfact.getAreaFrac()[1]->const_array(mfi);,
                              auto apz = ebfact.getAreaFrac()[2]->const_array(mfi); );

                AMREX_D_TERM( Array4<Real const> fcx = ebfact.getFaceCent()[0]->const_array(mfi);,
                              Array4<Real const> fcy = ebfact.getFaceCent()[1]->const_array(mfi);,
                              Array4<Real const> fcz = ebfact.getFaceCent()[2]->const_array(mfi););

                Array4<Real const> ccent_arr = ebfact.getCentroid().const_array(mfi);
                Array4<Real const> const& vfrac_arr = ebfact.getVolFrac().const_array(mfi);

                // This is scratch space if calling StateRedistribute,
                //  but is used as the weights (here set to 1) if calling
                //  FluxRedistribute
                Box gbx = bx;

                if (redistribution_type == "StateRedist")
                    gbx.grow(3);
                else if (redistribution_type == "FluxRedist")
                    gbx.grow(2);

                int tmpfab_comp = (is_sync) ? ncomp*2 : ncomp;
                FArrayBox tmpfab(gbx, tmpfab_comp, The_Async_Arena());
                Array4<Real> scratch = tmpfab.array(0);
                if (redistribution_type == "FluxRedist")
                {
                    amrex::ParallelFor(Box(scratch),
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    { scratch(i,j,k) = 1.;});
                }

                Array4<Real> redist_arr = (is_sync) ? tmpfab.array(ncomp) : aofs_arr;

                EBFluxRegister* fr_as_crse = nullptr;
                if (do_reflux && level < parent->finestLevel()) {
                    NavierStokesBase& flevel = getLevel(level+1);
                    fr_as_crse = flevel.advflux_reg.get();
                }

                EBFluxRegister* fr_as_fine = nullptr;
                if (do_reflux && level > 0) {
                    fr_as_fine = advflux_reg.get();
                }

                int as_crse = (fr_as_crse != nullptr);
                int as_fine = (fr_as_fine != nullptr);

                FArrayBox* p_drho_as_crse = (fr_as_crse) ?
                        fr_as_crse->getCrseData(mfi) : &fab_drho_as_crse;
                const IArrayBox* p_rrflag_as_crse = (fr_as_crse) ?
                       fr_as_crse->getCrseFlag(mfi) : &fab_rrflag_as_crse;

                if (fr_as_fine) {
                    dm_as_fine.resize(amrex::grow(bx,1),ncomp,The_Async_Arena());
                    dm_as_fine.template setVal<RunOn::Device>(0.0);
                }

                if (redistribution_type == "FluxRedist") {
                    bool use_wts_in_divnc = true;
                    ApplyMLRedistribution( bx, ncomp, redist_arr, update_arr,
                                           rstate->const_array(mfi, rstate_comp), scratch, flags_arr,
                                           AMREX_D_DECL(apx,apy,apz), vfrac_arr,
                                           AMREX_D_DECL(fcx,fcy,fcz), ccent_arr, bcrec_d,
                                           geom, dt, redistribution_type,
                                           as_crse, p_drho_as_crse->array(), p_rrflag_as_crse->array(),
                                           as_fine, dm_as_fine.array(), coarse_fine_mask->const_array(mfi),
                                           level_mask_notcovered, use_wts_in_divnc);
                } else {
                    bool use_wts_in_divnc = true;
                    ApplyRedistribution( bx, ncomp, redist_arr, update_arr,
                                         rstate->const_array(mfi, rstate_comp), scratch, flags_arr,
                                         AMREX_D_DECL(apx,apy,apz), vfrac_arr,
                                         AMREX_D_DECL(fcx,fcy,fcz), ccent_arr, bcrec_d,
                                         geom, dt, redistribution_type, use_wts_in_divnc );
                }

                if (is_sync)
                {
                    amrex::ParallelFor(bx, ncomp, [aofs_arr, redist_arr]
                    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    { aofs_arr( i, j, k, n ) -= redist_arr( i, j, k, n ); });
                }
                else
                {
                    amrex::ParallelFor(bx, ncomp, [aofs_arr, redist_arr]
                    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    { aofs_arr( i, j, k, n ) =  -redist_arr( i, j, k, n ); });
                }
            }
            else // bx is EB regular
            {
                // Recall that we computed -div in previous MFIter, because
                // redistribution operates on -div. Thus, we use -update here.
                if (is_sync)
                {
                    amrex::ParallelFor(bx, ncomp, [aofs_arr, update_arr]
                    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    { aofs_arr( i, j, k, n ) -= update_arr(i,j,k,n); });
                }
                else
                {
                    amrex::ParallelFor(bx, ncomp, [aofs_arr, update_arr]
                    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    { aofs_arr( i, j, k, n ) = -update_arr(i,j,k,n); });
                }
            }

            AMREX_D_TERM(FArrayBox fx_fr_fab(fx_fab,amrex::make_alias,flux_comp,ncomp);,
                         FArrayBox fy_fr_fab(fy_fab,amrex::make_alias,flux_comp,ncomp);,
                         FArrayBox fz_fr_fab(fz_fab,amrex::make_alias,flux_comp,ncomp););

            // Now update the flux registers (inside test on AMREX_USE_EB)
            if ( do_reflux && do_crse_add && (level < parent->finestLevel()) ) {
              if (flagfab.getType(amrex::grow(bx,1)) == FabType::regular)
              {
                   getAdvFluxReg(level+1).CrseAdd(mfi,
                       {AMREX_D_DECL(&fx_fr_fab,&fy_fr_fab,&fz_fr_fab)},
                       dxDp, dt, 0, state_indx, ncomp, amrex::RunOn::Device);

              } else if (flagfab.getType(bx) != FabType::covered ) {
                   getAdvFluxReg(level + 1).CrseAdd(mfi,
                      {AMREX_D_DECL(&fx_fr_fab,&fy_fr_fab,&fz_fr_fab)},
                      dxDp, dt, (*volfrac)[mfi],
                      {AMREX_D_DECL(&(*areafrac[0])[mfi], &(*areafrac[1])[mfi], &(*areafrac[2])[mfi])},
                      0, state_indx, ncomp, amrex::RunOn::Device);
              }
            } // do_reflux && level < finest_level

            // This is a hack-y way of testing whether this ComputeAofs call
            // came from the mac_sync (do_crse_add = false)
            // or from the regular advance (do_crse_add = true).  When the call
            // comes from the mac_sync, the multiplier in FineAdd needs to have
            // the opposite sign
            Real sync_factor = do_crse_add ? 1.0 : -1.0;

            if ( do_reflux && do_fine_add && (level > 0)) {
              if (flagfab.getType(amrex::grow(bx,1)) == FabType::regular)
              {
                  advflux_reg->FineAdd(mfi,
                     {AMREX_D_DECL(&fx_fr_fab,&fy_fr_fab,&fz_fr_fab)},
                     dxDp, sync_factor*dt, 0, state_indx, ncomp, amrex::RunOn::Device);
              } else if (flagfab.getType(bx) != FabType::covered ) {
                  advflux_reg->FineAdd(mfi,
                     {AMREX_D_DECL(&fx_fr_fab,&fy_fr_fab,&fz_fr_fab)},
                     dxDp, sync_factor*dt, (*volfrac)[mfi],
                     {AMREX_D_DECL(&(*areafrac[0])[mfi], &(*areafrac[1])[mfi], &(*areafrac[2])[mfi])},
                     dm_as_fine, 0, state_indx, ncomp, amrex::RunOn::Device);
              }
            } // do_reflux && (level > 0)
        } // not covered
#else
        // This is a hack-y way of testing whether this ComputeAofs call
        // came from the mac_sync (do_crse_add = false)
        // or from the regular advance (do_crse_add = true).  When the call
        // comes from the mac_sync, the multiplier in FineAdd needs to have
        // the opposite sign
        Real sync_factor = do_crse_add ? 1.0 : -1.0;

        // Update the flux registers when no EB
        if ( do_reflux && (level < parent->finestLevel()) ) {
               getAdvFluxReg(level+1).CrseAdd(mfi,
                              {AMREX_D_DECL(&fx_fab,&fy_fab,&fz_fab)},
                              dxDp, sync_factor*dt, flux_comp, state_indx, ncomp, amrex::RunOn::Device);
        } // do_reflux && level < finest_level

        if ( do_reflux && (level > 0) ) {
              advflux_reg->FineAdd(mfi,
                              {AMREX_D_DECL(&fx_fab,&fy_fab,&fz_fab)},
                               dxDp, sync_factor*dt, flux_comp, state_indx, ncomp, amrex::RunOn::Device);
        } // do_reflux && (level > 0)
#endif
    } // mfi
}


#ifdef AMREX_USE_EB
void
NavierStokesBase::InitialRedistribution ()
{
    // Next we must redistribute the initial solution if we are going to use
    // MergeRedist or StateRedist redistribution schemes
    if ( redistribution_type != "StateRedist" ) {
        return;
    }

    if (verbose) {
      amrex::Print() << "Doing initial redistribution... " << std::endl;
    }

    // Initial data are set at new time step
    MultiFab& S_new = get_new_data(State_Type);
    // We must fill internal ghost values before calling redistribution
    // We also need any physical boundary conditions imposed if we are
    //    calling state redistribution (because that calls the slope routine)
    FillPatchIterator S_fpi(*this, S_new, nghost_state(), state[State_Type].curTime(),
                            State_Type, 0, NUM_STATE);
    MultiFab& Smf=S_fpi.get_mf();

    MultiFab tmp( grids, dmap, NUM_STATE, nghost_state(), MFInfo(), Factory() );
    MultiFab::Copy(tmp, Smf, 0, 0, NUM_STATE, nghost_state());

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto const& fact =  dynamic_cast<EBFArrayBoxFactory const&>(S_new.Factory());

        EBCellFlagFab const& flagfab = fact.getMultiEBCellFlagFab()[mfi];
        Array4<EBCellFlag const> const& flag = flagfab.const_array();

        if ( (flagfab.getType(bx) != FabType::covered) &&
             (flagfab.getType(amrex::grow(bx,4)) != FabType::regular) )
        {
            Array4<Real const> AMREX_D_DECL(fcx, fcy, fcz), ccc, vfrac, AMREX_D_DECL(apx, apy, apz);

            AMREX_D_TERM(fcx = fact.getFaceCent()[0]->const_array(mfi);,
                         fcy = fact.getFaceCent()[1]->const_array(mfi);,
                         fcz = fact.getFaceCent()[2]->const_array(mfi););

            ccc   = fact.getCentroid().const_array(mfi);

            AMREX_D_TERM(apx = fact.getAreaFrac()[0]->const_array(mfi);,
                         apy = fact.getAreaFrac()[1]->const_array(mfi);,
                         apz = fact.getAreaFrac()[2]->const_array(mfi););

            vfrac = fact.getVolFrac().const_array(mfi);

            ApplyInitialRedistribution( bx, AMREX_SPACEDIM,
                                        Smf.array(mfi), tmp.array(mfi),
                                        flag, AMREX_D_DECL(apx, apy, apz), vfrac,
                                        AMREX_D_DECL(fcx, fcy, fcz),
                                        ccc, m_bcrec_velocity_d.dataPtr(),
                                        geom, redistribution_type);
            ApplyInitialRedistribution( bx, NUM_SCALARS,
                                        Smf.array(mfi,Density), tmp.array(mfi,Density),
                                        flag, AMREX_D_DECL(apx, apy, apz), vfrac,
                                        AMREX_D_DECL(fcx, fcy, fcz),
                                        ccc,m_bcrec_scalars_d.dataPtr(),
                                        geom, redistribution_type);
        }
    }

    MultiFab::Copy(S_new, Smf, 0, 0, NUM_STATE, 0);
}
#endif

//
// ls related
// 
void
NavierStokesBase::fill_allgts(MultiFab& mf, int type, int scomp, int ncomp, Real time)
{
    // Fill phys bc ,ff bc, cf bc
    int ngrow = mf.nGrow();
    FillPatchIterator mf_fpi(*this, mf, ngrow, time, type, scomp, ncomp);
    MultiFab& mf_temp = mf_fpi.get_mf();
    // MultiFab::Copy(mfdst, mfsrc, sc, dc, nc, ng); // Copy from mfsrc to mfdst
    MultiFab::Copy(mf, mf_temp, 0, scomp, ncomp, ngrow);
}

void
NavierStokesBase::reinit()
{

    if(verbose) amrex::Print() << "In the NavierStokesBase::reinit() " << std::endl;

    const Real* dx    = geom.CellSize();
    Real dxmin        = dx[0];
    for (int d=1; d<AMREX_SPACEDIM; ++d) {
        dxmin = std::min(dxmin,dx[d]);
    }

    if(reinit_levelset==1) {
        BL_ASSERT(do_cons_levelset==0);

        // Step 1: copy phi_ctime to phi_original
        MultiFab::Copy(phi_original, phi_ctime, 0, 0, 1, phi_ctime.nGrow());

        // Step 2: get pseudo dt for the current level
        const Real coeff = 0.5;
        Real dtlevel = coeff * dxmin;
        amrex::Print() << "level " << level << " dtlevel " << dtlevel << std::endl;

        // Step 3: only reinitialize the ls function on the current level
        for (int k=1; k<=number_of_reinit; k++)
        {
            reinitialization_sussman(dtlevel, k);
        }

    }
    else {
        BL_ASSERT(do_cons_levelset==1);
        const Real coeff = 0.75;
        Real epsG  = calculate_eps_one(geom, reinit_levelset);
        Real epsG2 = calculate_eps_two(geom, reinit_levelset);
        Real dtlevel = coeff * std::min({dxmin, dxmin * dxmin / (4.0 * (epsG + epsG2))});

        for (int k=1; k<=number_of_reinit; k++)
        {
            reinitialization_consls(dtlevel, k, epsG, epsG2);
        }

    }
}

void
NavierStokesBase::reinitialization_consls (Real dt,
                       int  loop_iter, Real epsG, Real epsG2)
{
    
    if (verbose) amrex::Print() << "In the NavierStokesBase::reinitialization_consls() " << std::endl;
    if (verbose) amrex::Print() << "loop_iter " << loop_iter << std::endl;

    // Real Gconverge = 0.00000000000001;
    // Real resmax = 1.e10;
    // Real res = 0.0;

    // Step 1:
    MultiFab::Copy(phi_original, phi_ctime, 0, 0, 1, phi_ctime.nGrow());

    for (int i=0; i<4; i++) {

        // Step 2: 
        MultiFab::Add(phi_ctime, phi_original, 0, 0, 1, phi_ctime.nGrow());
        phi_ctime.mult(0.5, phi_ctime.nGrow());

        // Step 3: copy phi_ctime back to phi in S_new, fill phi's bc data in S_new, then
        // copy it to phi_ctime
        const Real cur_time = state[State_Type].curTime();
        MultiFab&  S_new = get_new_data(State_Type);
        MultiFab::Copy(S_new, phi_ctime, 0, phicomp, 1, phi_ctime.nGrow());
        fill_allgts(S_new,State_Type,phicomp,1,cur_time);
        MultiFab::Copy(phi_ctime, S_new, phicomp, 0, 1, phi_ctime.nGrow());

        // Step 4: calculate level set normal
        Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM> phi_normal;
        int normalize = 1;
        cc_to_cc_grad(phi_normal, phi_ctime, geom, normalize);

        // Step 5: calculate diffs (phi2) and comp (phi1), add them to diffs_comp (override phi_ctime)
        MultiFab phi1(grids,dmap,1,2);
        MultiFab phi2(grids,dmap,1,2);
        phi1.setVal(0.0); phi2.setVal(0.0);
        levelset_diffcomp(phi_normal, phi_ctime, phi1, phi2, epsG, epsG2);

        // Step 6: update the level set
        phi_ctime.mult(dt, 0);
        MultiFab::Add(phi_ctime, phi_original, 0, 0, 1, 0);

        // Step 7: same as the Step 3
        MultiFab::Copy(S_new, phi_ctime, 0, phicomp, 1, phi_ctime.nGrow());
        fill_allgts(S_new,State_Type,phicomp,1,cur_time);
        MultiFab::Copy(phi_ctime, S_new, phicomp, 0, 1, phi_ctime.nGrow());

    }

}

void
NavierStokesBase::reinitialization_sussman (Real dt,
                       int  loop_iter)
{
    
    if (verbose) amrex::Print() << "In the NavierStokesBase::reinitialization_sussman() " << std::endl;
    if (verbose) amrex::Print() << "loop_iter " << loop_iter << std::endl;

    // Step 1: get sgn, similiar to heaviside
    if (loop_iter==1) {
        phi_to_sgn0(phi_original);
    }

    // Step 2: RK2
    // phi2, phi3, and G0 have 2 ghost cells, initialize them as 0,
    // and these multifabs only influence the minmod function at the boundary;
    MultiFab phi2(grids,dmap,1,2);
    MultiFab phi3(grids,dmap,1,2);
    MultiFab G0(grids,dmap,1,2);
    phi2.setVal(0.0); phi3.setVal(0.0); G0.setVal(0.0);
    rk_first_reinit(phi_ctime, phi2, phi3, sgn0, G0, dt, phi_original);

    // Step 3: copy phi_ctime back to phi in S_new, fill phi's bc data in S_new, then
    // copy it to phi_ctime
    const Real cur_time = state[State_Type].curTime();
    MultiFab&  S_new = get_new_data(State_Type);
    MultiFab::Copy(S_new, phi_ctime, 0, phicomp, 1, phi_ctime.nGrow());
    fill_allgts(S_new,State_Type,phicomp,1,cur_time);
    MultiFab::Copy(phi_ctime, S_new, phicomp, 0, 1, phi_ctime.nGrow());

    // Step 4: RK2
    rk_second_reinit(phi_ctime, phi2, phi3, sgn0, G0, dt, phi_original);

    // Step 5: same as the Step 3
    MultiFab::Copy(S_new, phi_ctime, 0, phicomp, 1, phi_ctime.nGrow());
    fill_allgts(S_new,State_Type,phicomp,1,cur_time);
    MultiFab::Copy(phi_ctime, S_new, phicomp, 0, 1, phi_ctime.nGrow());

    // Step 6: Fix mass
    // Step 6-1: set inputs variables as 0.0
    MultiFab ld(grids,dmap,1,2);
    MultiFab lambdad(grids,dmap,1,2);
    MultiFab deltafunc(grids,dmap,1,2);
    phi2.setVal(0.0); phi3.setVal(0.0); ld.setVal(0.0); lambdad.setVal(0.0); deltafunc.setVal(0.0);

    // Step 6-2: mass_fix
    mass_fix(phi_ctime, phi_original, phi2, phi3, ld, lambdad, deltafunc, dt, loop_iter);

    // Step 6-3: same as Step 5
    MultiFab::Copy(S_new, phi_ctime, 0, phicomp, 1, phi_ctime.nGrow());
    fill_allgts(S_new,State_Type,phicomp,1,cur_time);
    MultiFab::Copy(phi_ctime, S_new, phicomp, 0, 1, phi_ctime.nGrow());
}

void
NavierStokesBase::phi_to_sgn0 (MultiFab& phi)
{

    if(verbose) amrex::Print() << "In the NavierStokesBase::phi_to_sgn0 " << std::endl;
    
    sgn0.setVal(0.0);

    const Real pi     = 3.141592653589793238462643383279502884197;
    Real eps = calculate_eps(geom, epsilon);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phifab  =  phi.array(mfi);
        auto const& sgn0fab = sgn0.array(mfi);
        amrex::ParallelFor(vbx, [phifab, sgn0fab, pi, eps]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {

            if (phifab(i,j,k) > eps) {
                sgn0fab(i,j,k) = 1.0;
            } else if (phifab(i,j,k) > -eps) {
                sgn0fab(i,j,k) = phifab(i,j,k) / eps + 1.0 / pi * std::sin(phifab(i,j,k) * pi / eps);
            } else {
                sgn0fab(i,j,k) = -1.0;
            }

        });
    }

}

void
NavierStokesBase::rk_first_reinit (MultiFab& phi_ctime,
                          MultiFab& phi2,
                          MultiFab& phi3,
                          MultiFab& sgn0,
                          MultiFab& G0,
                          Real delta_t,
                          MultiFab& phi_ori)
{

    if(verbose) amrex::Print() << "NavierStokesBase::rk_first_reinit " << std::endl;
    
    Real eps = calculate_eps(geom, epsilon);
    const GpuArray<Real,AMREX_SPACEDIM> dxGpu = geom.CellSizeArray();

    // MultiFab phi1_xface(amrex::convert(grids, IntVect(AMREX_D_DECL(1,0,0))), dmap, 1, 1);
    MultiFab phi1_face(amrex::convert(grids, IntVect(AMREX_D_DECL(1,1,1))), dmap, 1, 1); // A node-based mf actually
    phi1_face.setVal(0.0);
    
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi1_face,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.growntilebox();
        // std::cout << "bx " << bx <<std::endl;
        auto const& phifab   = phi_ctime.array(mfi);
        auto const& phi1fab  = phi1_face.array(mfi);
        amrex::ParallelFor(bx, [phi1fab, phifab, dxGpu]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            phi1fab(i,j,k) = ( phifab(i,j,k) - phifab(i-1,j,k) )/dxGpu[0];
        });
    }
    
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi2,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.growntilebox(1);
        auto const& phi1fab  = phi1_face.array(mfi);
        auto const& phi2fab  = phi2.array(mfi);
        amrex::ParallelFor(bx, [phi1fab, phi2fab, dxGpu]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            phi2fab(i,j,k) = ( phi1fab(i+1,j,k) - phi1fab(i,j,k) )/dxGpu[0];
        });
    }
    
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi3,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phi1fab = phi1_face.array(mfi);
        auto const& phi2fab = phi2.array(mfi);
        auto const& phi3fab = phi3.array(mfi);
        auto const& sgn0fab = sgn0.array(mfi);
        amrex::ParallelFor(vbx, [phi1fab, phi2fab, phi3fab, sgn0fab, dxGpu]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real dxm = phi1fab(i, j, k) + minmod(phi2fab(i - 1, j, k), phi2fab(i, j, k)) * dxGpu[0] / 2.0;
            Real dxp = phi1fab(i + 1, j, k) - minmod(phi2fab(i, j, k), phi2fab(i + 1, j, k)) * dxGpu[0] / 2.0;
            Real ddx = 0.0;
            if (dxp * sgn0fab(i, j, k) < 0.0 && dxm * sgn0fab(i, j, k) < -dxp * sgn0fab(i, j, k)) {
                ddx = dxp;
            } else if (dxm * sgn0fab(i, j, k) > 0.0 && dxp * sgn0fab(i, j, k) > -dxm * sgn0fab(i, j, k)) {
                ddx = dxm;
            } else {
                ddx = (dxp + dxm) / 2.0;
            }
            phi3fab(i, j, k) = pow(ddx, 2);
        });
    }
    
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi1_face,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.growntilebox();
        auto const& phifab   = phi_ctime.array(mfi);
        auto const& phi1fab  = phi1_face.array(mfi);
        amrex::ParallelFor(bx, [phifab, phi1fab, dxGpu]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            phi1fab(i,j,k) = ( phifab(i,j,k) - phifab(i,j-1,k) )/dxGpu[1];
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi2,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.growntilebox(1);
        auto const& phi1fab  = phi1_face.array(mfi);
        auto const& phi2fab  = phi2.array(mfi);
        amrex::ParallelFor(bx, [phi1fab, phi2fab, dxGpu]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            phi2fab(i,j,k) = ( phi1fab(i,j+1,k) - phi1fab(i,j,k) )/dxGpu[1];
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi3,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phi1fab = phi1_face.array(mfi);
        auto const& phi2fab = phi2.array(mfi);
        auto const& phi3fab = phi3.array(mfi);
        auto const& sgn0fab = sgn0.array(mfi);
        amrex::ParallelFor(vbx, [phi1fab, phi2fab, phi3fab, sgn0fab, dxGpu]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real dym = phi1fab(i, j, k) + minmod(phi2fab(i, j - 1, k), phi2fab(i, j, k)) * dxGpu[1] / 2.0;
            Real dyp = phi1fab(i, j + 1, k) - minmod(phi2fab(i, j, k), phi2fab(i, j + 1, k)) * dxGpu[1] / 2.0;
            Real ddy = 0.0;
            if (dyp * sgn0fab(i, j, k) < 0.0 && dym * sgn0fab(i, j, k) < -dyp * sgn0fab(i, j, k)) {
                ddy = dyp;
            } else if (dym * sgn0fab(i, j, k) > 0.0 && dyp * sgn0fab(i, j, k) > -dym * sgn0fab(i, j, k)) {
                ddy = dym;
            } else {
                ddy = (dyp + dym) / 2.0;
            }
            phi3fab(i, j, k) += pow(ddy, 2);
            if (AMREX_SPACEDIM==2) {
                phi3fab(i, j, k) = std::sqrt(phi3fab(i, j, k));
            }
        });
    }
    
#if (AMREX_SPACEDIM==3)

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi1_face,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.growntilebox();
        auto const& phifab   = phi_ctime.array(mfi);
        auto const& phi1fab  = phi1_face.array(mfi);
        amrex::ParallelFor(bx, [phifab, phi1fab, dxGpu]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            phi1fab(i,j,k) = ( phifab(i,j,k) - phifab(i,j,k-1) )/dxGpu[2];
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi2,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.growntilebox(1);
        auto const& phi1fab  = phi1_face.array(mfi);
        auto const& phi2fab  = phi2.array(mfi);
        amrex::ParallelFor(bx, [phi1fab, phi2fab, dxGpu]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            phi2fab(i,j,k) = ( phi1fab(i,j,k+1) - phi1fab(i,j,k) )/dxGpu[2];
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi3,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phi1fab = phi1_face.array(mfi);
        auto const& phi2fab = phi2.array(mfi);
        auto const& phi3fab = phi3.array(mfi);
        auto const& sgn0fab = sgn0.array(mfi);
        amrex::ParallelFor(vbx, [phi1fab, phi2fab, phi3fab, sgn0fab, dxGpu]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real dzm = phi1fab(i, j, k) + minmod(phi2fab(i, j, k - 1), phi2fab(i, j, k)) * dxGpu[2] / 2.0;
            Real dzp = phi1fab(i, j, k + 1) - minmod(phi2fab(i, j, k), phi2fab(i, j, k + 1)) * dxGpu[2] / 2.0;
            Real ddz = 0.0;
            if (dzp * sgn0fab(i, j, k) < 0.0 && dzm * sgn0fab(i, j, k) < -dzp * sgn0fab(i, j, k)) {
                ddz = dzp;
            } else if (dzm * sgn0fab(i, j, k) > 0.0 && dzp * sgn0fab(i, j, k) > -dzm * sgn0fab(i, j, k)) {
                ddz = dzm;
            } else {
                ddz = (dzp + dzm) / 2.0;
            }
            phi3fab(i, j, k) += pow(ddz, 2);
            phi3fab(i, j, k) = std::sqrt(phi3fab(i, j, k));
        });
    }

#endif

    MultiFab::Copy(G0, phi3, 0, 0, 1, 1); // 1 gt 
    G0.plus(-1.0, 1); // 1 gt

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi2,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phifab  = phi_ctime.array(mfi);
        auto const& phi2fab = phi2.array(mfi);
        auto const& phi3fab = phi3.array(mfi);
        auto const& G0fab   = G0.array(mfi);
        amrex::ParallelFor(vbx, [phifab, phi2fab, phi3fab, G0fab, delta_t]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real sgntmp = phifab(i, j, k) /
                std::sqrt( pow(phifab(i, j, k), 2) + pow(phi3fab(i, j, k)*2.0*delta_t,2) );
            phi2fab(i, j, k) = phifab(i,j,k) - delta_t*sgntmp*G0fab(i,j,k);
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi2,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phifab  = phi_ctime.array(mfi);
        auto const& phi2fab = phi2.array(mfi);
        amrex::ParallelFor(vbx, [phifab, phi2fab, eps]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if (phi2fab(i,j,k)*phifab(i,j,k) < 0.0) {
                if (std::abs(phifab(i,j,k)) <= eps) {
                    phi2fab(i, j, k) = phifab(i,j,k) * 0.1;
                }
                else {
                    amrex::Abort("sign change in rk1 bulk flow!");
                }
            }
        });
    }

    MultiFab::Copy(phi_ctime, phi2, 0, 0, 1, 0); // only copy interior cells
}

void
NavierStokesBase::rk_second_reinit (MultiFab& phi_ctime,
                          MultiFab& phi2,
                          MultiFab& phi3,
                          MultiFab& sgn0,
                          MultiFab& G0,
                          Real delta_t, 
                          MultiFab& phi_ori)
{

    if(verbose) amrex::Print() << "In the NavierStokesBase::rk_second_reinit " << std::endl;

    Real eps = calculate_eps(geom, epsilon);
    const GpuArray<Real,AMREX_SPACEDIM> dxGpu = geom.CellSizeArray();

    // MultiFab phi1_xface(amrex::convert(grids, IntVect(AMREX_D_DECL(1,0,0))), dmap, 1, 1);
    MultiFab phi1_face(amrex::convert(grids, IntVect(AMREX_D_DECL(1,1,1))), dmap, 1, 1); // A node-based mf actually
    phi1_face.setVal(0.0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi1_face,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.growntilebox();
        auto const& phifab   = phi_ctime.array(mfi);
        auto const& phi1fab  = phi1_face.array(mfi);
        amrex::ParallelFor(bx, [phifab, phi1fab, dxGpu]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            phi1fab(i,j,k) = ( phifab(i,j,k) - phifab(i-1,j,k) )/dxGpu[0];
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi2,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.growntilebox(1);
        auto const& phi1fab  = phi1_face.array(mfi);
        auto const& phi2fab  = phi2.array(mfi);
        amrex::ParallelFor(bx, [phi1fab, phi2fab, dxGpu]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            phi2fab(i,j,k) = ( phi1fab(i+1,j,k) - phi1fab(i,j,k) )/dxGpu[0];
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi3,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phi1fab = phi1_face.array(mfi);
        auto const& phi2fab = phi2.array(mfi);
        auto const& phi3fab = phi3.array(mfi);
        auto const& sgn0fab = sgn0.array(mfi);
        amrex::ParallelFor(vbx, [phi1fab, phi2fab, phi3fab, sgn0fab, dxGpu]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real dxm = phi1fab(i, j, k) + minmod(phi2fab(i - 1, j, k), phi2fab(i, j, k)) * dxGpu[0] / 2.0;
            Real dxp = phi1fab(i + 1, j, k) - minmod(phi2fab(i, j, k), phi2fab(i + 1, j, k)) * dxGpu[0] / 2.0;
            Real ddx = 0.0;
            if (dxp * sgn0fab(i, j, k) < 0.0 && dxm * sgn0fab(i, j, k) < -dxp * sgn0fab(i, j, k)) {
                ddx = dxp;
            } else if (dxm * sgn0fab(i, j, k) > 0.0 && dxp * sgn0fab(i, j, k) > -dxm * sgn0fab(i, j, k)) {
                ddx = dxm;
            } else {
                ddx = (dxp + dxm) / 2.0;
            }
            phi3fab(i, j, k) = pow(ddx, 2);
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi1_face,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.growntilebox();
        auto const& phifab   = phi_ctime.array(mfi);
        auto const& phi1fab  = phi1_face.array(mfi);
        amrex::ParallelFor(bx, [phifab, phi1fab, dxGpu]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            phi1fab(i,j,k) = ( phifab(i,j,k) - phifab(i,j-1,k) )/dxGpu[1];
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi2,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.growntilebox(1);
        auto const& phi1fab  = phi1_face.array(mfi);
        auto const& phi2fab  = phi2.array(mfi);
        amrex::ParallelFor(bx, [phi1fab, phi2fab, dxGpu]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            phi2fab(i,j,k) = ( phi1fab(i,j+1,k) - phi1fab(i,j,k) )/dxGpu[1];
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi3,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phi1fab = phi1_face.array(mfi);
        auto const& phi2fab = phi2.array(mfi);
        auto const& phi3fab = phi3.array(mfi);
        auto const& sgn0fab = sgn0.array(mfi);
        amrex::ParallelFor(vbx, [phi1fab, phi2fab, phi3fab, sgn0fab, dxGpu]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real dym = phi1fab(i, j, k) + minmod(phi2fab(i, j - 1, k), phi2fab(i, j, k)) * dxGpu[1] / 2.0;
            Real dyp = phi1fab(i, j + 1, k) - minmod(phi2fab(i, j, k), phi2fab(i, j + 1, k)) * dxGpu[1] / 2.0;
            Real ddy = 0.0;
            if (dyp * sgn0fab(i, j, k) < 0.0 && dym * sgn0fab(i, j, k) < -dyp * sgn0fab(i, j, k)) {
                ddy = dyp;
            } else if (dym * sgn0fab(i, j, k) > 0.0 && dyp * sgn0fab(i, j, k) > -dym * sgn0fab(i, j, k)) {
                ddy = dym;
            } else {
                ddy = (dyp + dym) / 2.0;
            }
            phi3fab(i, j, k) += pow(ddy, 2);
            if (AMREX_SPACEDIM==2) {
                phi3fab(i, j, k) = std::sqrt(phi3fab(i, j, k));
            }
        });
    }
    
#if (AMREX_SPACEDIM==3)

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi1_face,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.growntilebox();
        auto const& phifab   = phi_ctime.array(mfi);
        auto const& phi1fab  = phi1_face.array(mfi);
        amrex::ParallelFor(bx, [phifab, phi1fab, dxGpu]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            phi1fab(i,j,k) = ( phifab(i,j,k) - phifab(i,j,k-1) )/dxGpu[2];
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi2,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.growntilebox(1);
        auto const& phi1fab  = phi1_face.array(mfi);
        auto const& phi2fab  = phi2.array(mfi);
        amrex::ParallelFor(bx, [phi1fab, phi2fab, dxGpu]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            phi2fab(i,j,k) = ( phi1fab(i,j,k+1) - phi1fab(i,j,k) )/dxGpu[2];
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi3,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phi1fab = phi1_face.array(mfi);
        auto const& phi2fab = phi2.array(mfi);
        auto const& phi3fab = phi3.array(mfi);
        auto const& sgn0fab = sgn0.array(mfi);
        amrex::ParallelFor(vbx, [phi1fab, phi2fab, phi3fab, sgn0fab, dxGpu]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real dzm = phi1fab(i, j, k) + minmod(phi2fab(i, j, k - 1), phi2fab(i, j, k)) * dxGpu[2] / 2.0;
            Real dzp = phi1fab(i, j, k + 1) - minmod(phi2fab(i, j, k), phi2fab(i, j, k + 1)) * dxGpu[2] / 2.0;
            Real ddz = 0.0;
            if (dzp * sgn0fab(i, j, k) < 0.0 && dzm * sgn0fab(i, j, k) < -dzp * sgn0fab(i, j, k)) {
                ddz = dzp;
            } else if (dzm * sgn0fab(i, j, k) > 0.0 && dzp * sgn0fab(i, j, k) > -dzm * sgn0fab(i, j, k)) {
                ddz = dzm;
            } else {
                ddz = (dzp + dzm) / 2.0;
            }
            phi3fab(i, j, k) += pow(ddz, 2);
            phi3fab(i, j, k) = std::sqrt(phi3fab(i, j, k));
        });
    }

#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi2,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phifab  = phi_ctime.array(mfi);
        auto const& phi2fab = phi2.array(mfi);
        auto const& phi3fab = phi3.array(mfi);
        auto const& G0fab   = G0.array(mfi);
        amrex::ParallelFor(vbx, [phifab, phi2fab, phi3fab, G0fab, delta_t]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real sgntmp = phifab(i, j, k) /
                std::sqrt( pow(phifab(i, j, k), 2) + pow(phi3fab(i, j, k)*2.0*delta_t,2) );
            phi2fab(i, j, k) = phifab(i,j,k) - 0.5*delta_t*sgntmp*(
                               phi3fab(i,j,k) - 1.0 - G0fab(i,j,k));
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi2,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phifab  = phi_ctime.array(mfi);
        auto const& phi2fab = phi2.array(mfi);
        amrex::ParallelFor(vbx, [phifab, phi2fab, eps]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if (phi2fab(i,j,k)*phifab(i,j,k) < 0.0) {
                if (std::abs(phifab(i,j,k)) <= eps) {
                    phi2fab(i, j, k) = phifab(i,j,k) * 0.1;
                }
                else {
                    amrex::Abort("sign change in rk2 bulk flow!");
                }
            }
        });
    }

    MultiFab::Copy(phi_ctime, phi2, 0, 0, 1, 0); // only copy interior cells

    // {
    //   // amrex::Gpu::LaunchSafeGuard lsg(false); // no needed here!
    //   int idx = 0;
    //   amrex::Print() << "phi_ctime " << " " << phi_ctime.max(idx,0) << " " << phi_ctime.min(idx,0) << " " << phi_ctime.norm2(0) << "\n";
    // }
    // amrex::Abort("stop here");

}

void
NavierStokesBase::mass_fix (MultiFab& phi_ctime,
                          MultiFab& phi_original,
                          MultiFab& phi2,
                          MultiFab& phi3,
                          MultiFab& ld,
                          MultiFab& la,
                          MultiFab& deltafunc,
                          Real delta_t,
                          int loop_iter)
{
    if(verbose) amrex::Print() << "NavierStokesBase::mass_fix " << std::endl;
    
    const Real pi     = 3.141592653589793238462643383279502884197;
    Real eps = calculate_eps(geom, epsilon);
    Real tao = loop_iter * delta_t;

    // calculate ld and delta
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi_ctime,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();
        auto const& phifab    = phi_ctime.array(mfi);
        auto const& phiorifab = phi_original.array(mfi);
        auto const& deltafab  = deltafunc.array(mfi);
        auto const& ldfab     = ld.array(mfi);
        amrex::ParallelFor(bx, [phifab, phiorifab, deltafab, ldfab, pi, eps, tao]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            ldfab(i,j,k) = (phifab(i,j,k) - phiorifab(i,j,k)) / tao;
            if (phifab(i,j,k) > eps) {
                deltafab(i,j,k) = 0.0;
            } else if (phifab(i,j,k) > -eps) {
                deltafab(i,j,k) = 0.5 * (1.0 + std::cos(phiorifab(i,j,k) * pi / eps)) / eps;
            } else {
                deltafab(i,j,k) = 0.0;
            }
        });
    }

    // get grid value of numerator (phi_2) and denominator (phi_3)
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi2,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();
        auto const& phi2fab   = phi2.array(mfi);
        auto const& phi3fab   = phi3.array(mfi);
        auto const& deltafab  = deltafunc.array(mfi);
        auto const& ldfab     = ld.array(mfi);
        amrex::ParallelFor(bx, [phi2fab, phi3fab, deltafab, ldfab]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            phi2fab(i,j,k) = -1.0 * deltafab(i,j,k) * ldfab(i,j,k);
            phi3fab(i,j,k) = deltafab(i,j,k) * deltafab(i,j,k);
        });
    }

#if (AMREX_SPACEDIM==2)
    // numerical integration
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi2,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phi2fab   = phi2.array(mfi);
        auto const& lafab     = la.array(mfi);
        amrex::ParallelFor(vbx, [phi2fab, lafab]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real temp_n = 0.0;
            for(int jj=-1; jj<=1; jj++) {
                for(int ii=-1; ii<=1; ii++) {
                    temp_n += phi2fab(i+ii,j+jj,k);
                }
            }
            lafab(i,j,k) = 15.0 * phi2fab(i,j,k) + temp_n;
        });
    }
#else
    // numerical integration
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi2,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phi2fab   = phi2.array(mfi);
        auto const& lafab     = la.array(mfi);
        amrex::ParallelFor(vbx, [phi2fab, lafab]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real temp_n = 0.0;
            for(int kk=-1; kk<=1; kk++) {
                for(int jj=-1; jj<=1; jj++) {
                    for(int ii=-1; ii<=1; ii++) {
                        temp_n += phi2fab(i+ii,j+jj,k+kk);
                    }
                }
            }
            lafab(i,j,k) = 51.0 * phi2fab(i,j,k) + temp_n;
        });
    }
#endif

    // copy
    MultiFab::Copy(phi2, la, 0, 0, 1, 0); // only copy interior cells

#if (AMREX_SPACEDIM==2)
    // numerical integration
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi3,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phi3fab   = phi3.array(mfi);
        auto const& lafab     = la.array(mfi);
        amrex::ParallelFor(vbx, [phi3fab, lafab]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real temp_d = 0.0;
            for(int jj=-1; jj<=1; jj++) {
                for(int ii=-1; ii<=1; ii++) {
                    temp_d += phi3fab(i+ii,j+jj,k);
                }
            }
            lafab(i,j,k) = 15.0 * phi3fab(i,j,k) + temp_d;
        });
    }
#else
    // numerical integration
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi3,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phi3fab   = phi3.array(mfi);
        auto const& lafab     = la.array(mfi);
        amrex::ParallelFor(vbx, [phi3fab, lafab]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real temp_d = 0.0;
            for(int kk=-1; kk<=1; kk++) {
                for(int jj=-1; jj<=1; jj++) {
                    for(int ii=-1; ii<=1; ii++) {
                        temp_d += phi3fab(i+ii,j+jj,k+kk);
                    }
                }
            }
            lafab(i,j,k) = 51.0 * phi3fab(i,j,k) + temp_d;
        });
    }
#endif

    // copy
    MultiFab::Copy(phi3, la, 0, 0, 1, 0); // only copy interior cells

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi3,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phi2fab   = phi2.array(mfi);
        auto const& phi3fab   = phi3.array(mfi);
        auto const& lafab     = la.array(mfi);
        amrex::ParallelFor(vbx, [phi2fab, phi3fab, lafab]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            lafab(i,j,k) = phi2fab(i,j,k) / (phi3fab(i,j,k) + 1.e-12);
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi_ctime,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phifab    = phi_ctime.array(mfi);
        auto const& lafab     = la.array(mfi);
        auto const& deltafab  = deltafunc.array(mfi);
        amrex::ParallelFor(vbx, [phifab, lafab, deltafab, tao]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            phifab(i,j,k) += tao * lafab(i,j,k) * deltafab(i,j,k);
        });
    }

    // {
    //   // amrex::Gpu::LaunchSafeGuard lsg(false); // no needed here!
    //   int idx = 0;
    //   amrex::Print() << "phi_ctime " << " " << phi_ctime.max(idx,0) << " " << phi_ctime.min(idx,0) << " " << phi_ctime.norm2(0) << "\n";
    // }
    // amrex::Abort("stop here");

}

MultiFab&
NavierStokesBase::get_phi_half_time ()
{
    //
    // Fill it in when needed ...
    //
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi_half,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();
        auto const& phi_h = phi_half.array(mfi);
        auto const& phi_p = phi_ptime.array(mfi);
        auto const& phi_c = phi_ctime.array(mfi);
        amrex::ParallelFor(bx, [phi_h, phi_p, phi_c]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
           phi_h(i,j,k) = 0.5 * (phi_p(i,j,k) + phi_c(i,j,k));
        });
    }
    return phi_half;
}