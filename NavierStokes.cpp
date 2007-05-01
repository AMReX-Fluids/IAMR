//
// "Divu_Type" means S, where divergence U = S
// "Dsdt_Type" means pd S/pd t, where S is as above
//
#include <winstd.H>

#include <algorithm>
#include <vector>
#include <cstdio>
#include <cmath>

#include <CoordSys.H>
#include <Geometry.H>
#include <BoxDomain.H>
#include <ParmParse.H>
#include <ErrorList.H>
#include <FArrayBox.H>
#include <Godunov.H>
#include <Interpolater.H>
#include <NavierStokes.H>
#include <MultiGrid.H>
#include <ArrayLim.H>
#include <Utility.H>
#include <NAVIERSTOKES_F.H>
#include <Profiler.H>
#include <PROJECTION_F.H>
#include <PROB_NS_F.H>
#include <TagBox.H>
#include <VISCOPERATOR_F.H>

#define GEOM_GROW   1
#define HYP_GROW    3
#define PRESS_GROW  1
#define DIVU_GROW   1
#define DSDT_GROW   1

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

//
// Static objects.
//
ErrorList   NavierStokes::err_list;
BCRec       NavierStokes::phys_bc;
Projection *NavierStokes::projector     = 0;
MacProj    *NavierStokes::mac_projector = 0;
Godunov    *NavierStokes::godunov       = 0;

//
// Internal parameters.
//
int  NavierStokes::verbose      = 0;
Real NavierStokes::cfl          = 0.8;
Real NavierStokes::init_shrink  = 1.0;
Real NavierStokes::change_max   = 1.1;
Real NavierStokes::fixed_dt     = -1.0;
Real NavierStokes::dt_cutoff    = 0.0;
int  NavierStokes::init_iter    = 2;
Real NavierStokes::gravity      = 0.0;
int  NavierStokes::initial_step = false;
int  NavierStokes::initial_iter = false;
int  NavierStokes::radius_grow  = 1;
int  NavierStokes::sum_interval = -1;
int  NavierStokes::turb_interval= -1;
int  NavierStokes::NUM_SCALARS  = 0;
int  NavierStokes::NUM_STATE    = 0;

Array<AdvectionForm> NavierStokes::advectionType;
Array<DiffusionForm> NavierStokes::diffusionType;

bool NavierStokes::def_harm_avg_cen2edge = false;


//
// ----------------------- viscosity parameters.
//
Real NavierStokes::be_cn_theta  = 0.5;
Real NavierStokes::visc_tol     = 1.0e-10;  // tolerance for viscous solve
Real NavierStokes::visc_abs_tol = 1.0e-10;  // absolute tol. for visc solve
int  NavierStokes::variable_vel_visc  = 0;  // variable viscosity flag
int  NavierStokes::variable_scal_diff = 0;  // variable scalar diffusion flag

Array<int>  NavierStokes::is_diffusive;
Array<Real> NavierStokes::visc_coef;

//
// Internal switches.
//
int  NavierStokes::do_temp                    = 0;
int  NavierStokes::do_trac2                   = 0;
int  NavierStokes::Temp                       = -1;
int  NavierStokes::Tracer                     = -1; // AJA
int  NavierStokes::Tracer2                    = -1; // AJA
int  NavierStokes::do_sync_proj               = 1;
int  NavierStokes::do_MLsync_proj             = 1;
int  NavierStokes::do_reflux                  = 1;
int  NavierStokes::modify_reflux_normal_vel   = 0;
int  NavierStokes::do_mac_proj                = 1;
int  NavierStokes::do_init_vort_proj          = 0;
int  NavierStokes::do_init_proj               = 1;
int  NavierStokes::do_refine_outflow          = 0;
int  NavierStokes::do_derefine_outflow        = 1;
int  NavierStokes::Nbuf_outflow               = 1;
int  NavierStokes::do_running_statistics      = 0;
int  NavierStokes::do_denminmax               = 0;
int  NavierStokes::do_scalminmax              = 0;
int  NavierStokes::do_density_ref             = 0;
int  NavierStokes::do_tracer_ref              = 0;
int  NavierStokes::do_vorticity_ref           = 0;

int  NavierStokes::Dpdt_Type                  = -1;

int  NavierStokes::do_mom_diff                = 0;
int  NavierStokes::predict_mom_together       = 0;

//     
// New members for non-zero divu.
//
int  NavierStokes::additional_state_types_initialized = 0;
int  NavierStokes::Divu_Type                          = -1;
int  NavierStokes::Dsdt_Type                          = -1;
int  NavierStokes::have_divu                          = 0;
int  NavierStokes::have_dsdt                          = 0;
int  NavierStokes::S_in_vel_diffusion                 = 1;

Real NavierStokes::divu_relax_factor   = 0.0;
     
int  NavierStokes::num_state_type = 2;     // for backward compatibility

int  NavierStokes::do_divu_sync = 0;       // for debugging new correction to MLSP

static
BoxArray
GetBndryCells (const BoxArray& ba,
               int             ngrow,
               const Geometry& geom)
{
    //
    // First get list of all ghost cells.
    //
    BoxList gcells, bcells;

    for (int i = 0; i < ba.size(); ++i)
	gcells.join(BoxLib::boxDiff(BoxLib::grow(ba[i],ngrow),ba[i]));
    //
    // Now strip out intersections with original BoxArray.
    //
    for (BoxList::const_iterator it = gcells.begin(); it != gcells.end(); ++it)
    {
        std::vector< std::pair<int,Box> > isects = ba.intersections(*it);

        if (isects.empty())
            bcells.push_back(*it);
        else
        {
            //
            // Collect all the intersection pieces.
            //
            BoxList pieces;
            for (int i = 0; i < isects.size(); i++)
                pieces.push_back(isects[i].second);
            BoxList leftover = BoxLib::complementIn(*it,pieces);
            bcells.catenate(leftover);
        }
    }
    //
    // Now strip out overlaps.
    //
    gcells.clear();
    gcells = BoxLib::removeOverlap(bcells);
    bcells.clear();

    if (geom.isAnyPeriodic())
    {
        Array<IntVect> pshifts(27);

        const Box domain = geom.Domain();

        for (BoxList::const_iterator it = gcells.begin(); it != gcells.end(); ++it)
        {
            if (!domain.contains(*it))
            {
                //
                // Add in periodic ghost cells shifted to valid region.
                //
                geom.periodicShift(domain, *it, pshifts);

                for (int i = 0; i < pshifts.size(); i++)
                {
                    const Box shftbox = *it + pshifts[i];

                    bcells.push_back(domain & shftbox);
                }
            }
        }

        gcells.catenate(bcells);
    }

    return BoxArray(gcells);
}

void
NavierStokes::variableCleanUp ()
{
    desc_lst.clear();
    delete projector;
    projector = 0;
    delete mac_projector;
    mac_projector = 0;
    delete godunov;
    godunov = 0;
}

void
NavierStokes::read_geometry ()
{
#if (BL_SPACEDIM == 2)
    //
    // Must load coord here because CoordSys hasn't read it in yet.
    //
    ParmParse pp("geometry");

    int coord;
    pp.get("coord_sys",coord);

    if ((CoordSys::CoordType) coord == CoordSys::RZ && phys_bc.lo(0) != Symmetry)
    {
        phys_bc.setLo(0,Symmetry);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "\nWarning: Setting phys_bc at xlo to Symmetry\n\n";
    }
#endif
}

void
NavierStokes::read_params ()
{
    ParmParse pp("ns");

    pp.query("v",verbose);

    Array<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
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
    if (Geometry::isAnyPeriodic())
    {
        //
        // Do idiot check.  Periodic means interior in those directions.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (Geometry::isPeriodic(dir))
            {
                if (lo_bc[dir] != Interior)
                {
                    std::cerr << "NavierStokes::variableSetUp:periodic in direction "
                              << dir
                              << " but low BC is not Interior\n";
                    BoxLib::Abort("NavierStokes::read_params()");
                }
                if (hi_bc[dir] != Interior)
                {
                    std::cerr << "NavierStokes::variableSetUp:periodic in direction "
                              << dir
                              << " but high BC is not Interior\n";
                    BoxLib::Abort("NavierStokes::read_params()");
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
            if (!Geometry::isPeriodic(dir))
            {
              if (lo_bc[dir] == Interior)
              {
                  std::cerr << "NavierStokes::variableSetUp:Interior bc in direction "
                            << dir
                            << " but not defined as periodic\n";
                  BoxLib::Abort("NavierStokes::read_params()");
              }
              if (hi_bc[dir] == Interior)
              {
                  std::cerr << "NavierStokes::variableSetUp:Interior bc in direction "
                            << dir
                            << " but not defined as periodic\n";
                  BoxLib::Abort("NavierStokes::read_params()");
              }
            }
        }
    }
    //
    // Get timestepping parameters.
    //
    pp.get("cfl",cfl);
    pp.query("init_iter",init_iter);
    pp.query("init_shrink",init_shrink);
    pp.query("dt_cutoff",dt_cutoff);
    pp.query("change_max",change_max);
    pp.query("fixed_dt",fixed_dt);
    pp.query("sum_interval",sum_interval);
    pp.query("turb_interval",turb_interval);
    pp.query("gravity",gravity);

    //
    // Get run options.
    //
    pp.query("do_temp",                  do_temp          );
    pp.query("do_trac2",                 do_trac2         );
    int initial_do_sync_proj =           do_sync_proj;
    pp.query("do_sync_proj",             do_sync_proj     );
    pp.query("do_MLsync_proj",           do_MLsync_proj   );
    pp.query("do_reflux",                do_reflux        );
    pp.query("modify_reflux_normal_vel", modify_reflux_normal_vel);
    pp.query("do_init_vort_proj",        do_init_vort_proj);
    pp.query("do_init_proj",             do_init_proj     );
    pp.query("do_mac_proj",              do_mac_proj      );
    pp.query("do_divu_sync",             do_divu_sync     );
    pp.query("do_denminmax",             do_denminmax     );
    pp.query("do_scalminmax",            do_scalminmax    );
    pp.query("do_density_ref",           do_density_ref   );
    pp.query("do_tracer_ref",            do_tracer_ref    );
    pp.query("do_vorticity_ref",         do_vorticity_ref );

    //
    // Make sure we don't use divu_sync.
    //
    if (do_divu_sync)
        BoxLib::Error("do_divu_sync == 1 is the wrong setting");
    //
    // This test ensures that if the user toggles do_sync_proj,
    // the user has knowledge that do_MLsync_proj is meaningless.
    //
    if (do_MLsync_proj && !do_sync_proj && initial_do_sync_proj != do_sync_proj)
    {
        std::cout << "Mismatched options for NavierStokes\n"
                  << "do_MLsync_proj and do_sync_proj are inconsistent\n";

        BoxLib::Abort("NavierStokes::read_params()");
    }
    //
    // Read viscous/diffusive parameters and array of viscous/diffusive coeffs.
    // NOTE: at this point, we dont know number of state variables
    //       so just read all values listed.
    //
    pp.query("visc_tol",visc_tol);
    pp.query("visc_abs_tol",visc_abs_tol);
    pp.query("variable_vel_visc",variable_vel_visc);
    pp.query("variable_scal_diff",variable_scal_diff);

    const int n_vel_visc_coef   = pp.countval("vel_visc_coef");
    const int n_temp_cond_coef  = pp.countval("temp_cond_coef");
    const int n_scal_diff_coefs = pp.countval("scal_diff_coefs");

    if (n_vel_visc_coef != 1)
        BoxLib::Abort("NavierStokes::read_params(): Only one vel_visc_coef allowed");

    if (do_temp && n_temp_cond_coef != 1)
        BoxLib::Abort("NavierStokes::read_params(): Only one temp_cond_coef allowed");

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
    Array<Real> scal_diff_coefs(n_scal_diff_coefs);
    pp.getarr("scal_diff_coefs",scal_diff_coefs,0,n_scal_diff_coefs);

    int scalId = Density;

    // Will need to add more lines when more variables are added - AJA
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
    
    pp.query("divu_relax_factor",divu_relax_factor);
    pp.query("S_in_vel_diffusion",S_in_vel_diffusion);
    pp.query("be_cn_theta",be_cn_theta);
    if (be_cn_theta > 1.0 || be_cn_theta < .5)
        BoxLib::Abort("NavierStokes::read_params(): Must have be_cn_theta <= 1.0 && >= .5");
    //
    // Set parameters dealing with how grids are treated at outflow boundaries.
    //
    pp.query("do_refine_outflow",do_refine_outflow);
    pp.query("do_derefine_outflow",do_derefine_outflow);
    if (do_derefine_outflow) do_refine_outflow = 0;

    pp.query("Nbuf_outflow",Nbuf_outflow);
    BL_ASSERT(Nbuf_outflow >= 0);
    BL_ASSERT(!(Nbuf_outflow <= 0 && do_derefine_outflow == 1));

    //
    // Check whether we are doing running statistics.
    //
    pp.query("do_running_statistics",do_running_statistics);

    //
    // Are we going to do velocity or momentum update?
    //
    pp.query("do_mom_diff",do_mom_diff);
    pp.query("predict_mom_together",predict_mom_together);

    if (do_mom_diff == 0 && predict_mom_together == 1)
    {
      std::cout << "MAKES NO SENSE TO HAVE DO_MOM_DIFF=0 AND PREDICT_MOM_TOGETHER=1" << std::endl;
      exit(0);
    }

    pp.query("harm_avg_cen2edge", def_harm_avg_cen2edge);

}

NavierStokes::NavierStokes ()
{
    rho_avg      = 0;
    rho_half     = 0;
    rho_ptime    = 0;
    rho_ctime    = 0;
    rho_qtime    = 0;
    rho_tqtime   = 0;
    p_avg        = 0;
    Vsync        = 0;
    Ssync        = 0;
    sync_reg     = 0;
    advflux_reg  = 0;
    viscflux_reg = 0;
    u_mac        = 0;
    u_macG       = 0;
    u_corr       = 0;
    mac_rhs      = 0;
    aofs         = 0;
    diffusion    = 0;

    if (!additional_state_types_initialized)
        init_additional_state_types();
}

NavierStokes::NavierStokes (Amr&            papa,
                            int             lev,
                            const Geometry& level_geom,
                            const BoxArray& bl,
                            Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,time)
{
    if(!additional_state_types_initialized)
        init_additional_state_types();
    
    const BoxArray& P_grids = state[Press_Type].boxArray();
    //
    // Alloc old_time pressure.
    //
    state[Press_Type].allocOldData();
    //
    // Alloc space for density and temporary pressure variables.
    //
    p_avg   = 0;
    rho_avg = 0;
    if (level > 0)
    {
        rho_avg = new MultiFab(grids,1,1);
        p_avg   = new MultiFab(P_grids,1,0);
    }

    rho_half   = new MultiFab(grids,1,1);
    rho_ptime  = new MultiFab(grids,1,1);
    rho_ctime  = new MultiFab(grids,1,1);
    rho_qtime  = 0;
    rho_tqtime = 0;
    //
    // Build metric coefficients for RZ calculations.
    //
    buildMetrics();
    //
    // Set up reflux registers.
    //
    sync_reg = 0;
    if (level > 0 && do_sync_proj)
    {
        sync_reg = new SyncRegister(grids,crse_ratio,level);
    }
    advflux_reg  = 0;
    viscflux_reg = 0;
    if (level > 0 && do_reflux)
    {
        advflux_reg  = new FluxRegister(grids,crse_ratio,level,NUM_STATE);
        viscflux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);
    }
    //
    // Initialize work multifabs.
    //
    Vsync   = 0;
    Ssync   = 0;
    u_mac   = 0;
    u_macG  = 0;
    u_corr  = 0;
    mac_rhs = 0;
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
                              NUM_STATE,viscflux_reg,volume,area,
                              is_diffusive,visc_coef);
    //
    // Allocate the storage for variable viscosity and diffusivity
    //
    viscn_cc   = 0;
    viscnp1_cc = 0;
    if (variable_vel_visc) 
    {
        viscn_cc   = new MultiFab(grids, 1, 1);
        viscnp1_cc = new MultiFab(grids, 1, 1);
    }

    diffn_cc   = 0;
    diffnp1_cc = 0;
    if (variable_scal_diff) 
    {
        diffn_cc   = new MultiFab(grids, NUM_STATE-Density-1, 1);
        diffnp1_cc = new MultiFab(grids, NUM_STATE-Density-1, 1);
    }
    //
    // Set up the mac projector.
    //
    if (mac_projector == 0)
    {
        mac_projector = new MacProj(parent,parent->finestLevel(),
                                    &phys_bc,radius_grow);
    }
    mac_projector->install_level(level,this,volume,area,&radius);
}

NavierStokes::~NavierStokes ()
{
    delete rho_avg;
    delete p_avg;
    delete rho_half;
    delete rho_ptime;
    delete rho_ctime;
    delete rho_qtime;
    delete rho_tqtime;
    delete Vsync;
    delete Ssync;
    delete sync_reg;
    delete advflux_reg;
    delete viscflux_reg;
    delete [] u_mac;
    delete [] u_macG;
    delete [] u_corr;
    
    if (mac_projector != 0)
        mac_projector->cleanup(level);
    //
    // Remove the arrays for variable viscosity and diffusivity
    // and delete the Diffusion object
    //
    if (variable_vel_visc)
    {
        delete viscn_cc;
        delete viscnp1_cc;
    }

    if (variable_scal_diff)
    {
        delete diffn_cc;
        delete diffnp1_cc;
    }

    delete diffusion;
}

void
NavierStokes::init_additional_state_types ()
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
    have_divu = 0;
    have_divu = isStateVariable("divu", dummy_Divu_Type, _Divu);
    have_divu = have_divu && dummy_Divu_Type == Divu_Type;
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "NavierStokes::init_additional_state_types()::have_divu = "
                  << have_divu << '\n';
    }
    if (have_divu && _Divu!=Divu)
    {
        std::cout << "divu must be 0-th, Divu_Type component in the state\n";

        BoxLib::Abort("NavierStokes::init_additional_state_types()");
    }

    if (have_divu && do_sync_proj && !do_MLsync_proj) 
    {
        std::cout << "Must run the ML sync project if have_divu is true " << std::endl;
        std::cout << "  because the divu sync is only implemented in the " << std::endl;
        std::cout << "  multilevel sync (MLsyncProject), not in the single level " << std::endl;
        std::cout << "  (syncProject)." << std::endl;
        BoxLib::Abort("NavierStokes::init_additional_state_types()");
    }

    int _Dsdt = -1;
    int dummy_Dsdt_Type;
    have_dsdt = 0;
    have_dsdt = isStateVariable("dsdt", dummy_Dsdt_Type, _Dsdt);
    have_dsdt = have_dsdt && dummy_Dsdt_Type==Dsdt_Type;
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "NavierStokes::init_additional_state_types()::have_dsdt = "
                  << have_dsdt << '\n';
    }
    if (have_dsdt && _Dsdt!=Dsdt)
    {
        std::cout << "dsdt must be 0-th, Dsdt_Type component in the state\n";

        BoxLib::Abort("NavierStokes::init_additional_state_types()");
    }
    if (have_dsdt && !have_divu)
    {
        std::cout << "Must have divu in order to have dsdt\n";

        BoxLib::Abort("NavierStokes::init_additional_state_types()");
    }

    num_state_type = desc_lst.size();
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "NavierStokes::init_additional_state_types: num_state_type = "
                  << num_state_type << '\n';
    }
}

//
// Since the pressure solver always stores its estimate of the
// pressure solver in Pnew, we need to copy it to Pold at the start.
//
void
NavierStokes::initOldPress ()
{
    MultiFab& P_new = get_new_data(Press_Type);
    MultiFab& P_old = get_old_data(Press_Type);

    for (MFIter mfi(P_new); mfi.isValid(); ++mfi)
    {
        P_old[mfi].copy(P_new[mfi]);
    }
}

void
NavierStokes::zeroNewPress ()
{
    get_new_data(Press_Type).setVal(0);
}

void
NavierStokes::zeroOldPress ()
{
    get_old_data(Press_Type).setVal(0);
}

void
NavierStokes::allocOldData ()
{
    bool init_pres = !(state[Press_Type].hasOldData());

    for (int k = 0; k < num_state_type; k++)
    {
        state[k].allocOldData();
    }
    if (init_pres)
        initOldPress();
}

void
NavierStokes::removeOldData ()
{
    AmrLevel::removeOldData();
}

void
NavierStokes::SetGodunov()
{
    if (godunov == 0)
        godunov = new Godunov();
}

void
NavierStokes::restart (Amr&          papa,
                       std::istream& is,
                       bool          bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

    if (projector == 0)
    {
        projector = new Projection(parent,&phys_bc,do_sync_proj,
                                   parent->finestLevel(),radius_grow);
    }
    projector->install_level(level, this, &radius );
    //
    // Set the godunov box.
    //
    SetGodunov();
    
    if (mac_projector == 0)
    {
        mac_projector = new MacProj(parent,parent->finestLevel(),
                                    &phys_bc,radius_grow );
    }
    mac_projector->install_level(level,this,volume,area,&radius );

    rho_avg = 0;
    p_avg   = 0;
    const BoxArray& P_grids = state[Press_Type].boxArray();
    //
    // Alloc space for density and temporary pressure variables.
    //
    if (level > 0)
    {
        rho_avg = new MultiFab(grids,1,1);
        p_avg   = new MultiFab(P_grids,1,0);
    }
    rho_half   = new MultiFab(grids,1,1);
    rho_ptime  = new MultiFab(grids,1,1);
    rho_ctime  = new MultiFab(grids,1,1);
    rho_qtime  = 0;
    rho_tqtime = 0;
    //
    // Build metric coefficients for RZ calculations.
    //
    buildMetrics();

    BL_ASSERT(sync_reg == 0);
    if (level > 0 && do_sync_proj)
    {
        sync_reg = new SyncRegister(grids,crse_ratio,level);
    }
    BL_ASSERT(advflux_reg == 0);
    if (level > 0 && do_reflux)
    {
        advflux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);
    }
    BL_ASSERT(viscflux_reg == 0);
    if (level > 0 && do_reflux)
    {
        viscflux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);
    }

    BL_ASSERT(Vsync == 0);
    BL_ASSERT(Ssync == 0);
    if (level < parent->finestLevel())
    {
        Vsync = new MultiFab(grids,BL_SPACEDIM,1);
        Ssync = new MultiFab(grids,NUM_STATE-BL_SPACEDIM,1);
    }

    diffusion = new Diffusion(parent, this,
                              (level > 0) ? getLevel(level-1).diffusion : 0,
                              NUM_STATE, viscflux_reg, volume, area,
                              is_diffusive, visc_coef);
    //
    // Allocate the storage for variable viscosity and diffusivity
    //
    viscn_cc   = 0;
    viscnp1_cc = 0;
    if (variable_vel_visc)
    {
        viscn_cc   = new MultiFab(grids, 1, 1);
        viscnp1_cc = new MultiFab(grids, 1, 1);
    }

    diffn_cc   = 0;
    diffnp1_cc = 0;
    if (variable_scal_diff)
    {
        diffn_cc   = new MultiFab(grids, NUM_STATE-Density-1, 1);
        diffnp1_cc = new MultiFab(grids, NUM_STATE-Density-1, 1);
    }

    is_first_step_after_regrid = false;
    old_intersect_new          = grids;
}

void
NavierStokes::buildMetrics ()
{
    radius.resize(grids.size());

    const Real dxr = geom.CellSize()[0];

    for (int i = 0; i < grids.size(); i++)
    {
        const int ilo = grids[i].smallEnd(0)-radius_grow;
        const int ihi = grids[i].bigEnd(0)+radius_grow;
        const int len = ihi - ilo + 1;

        radius[i].resize(len);

        if (CoordSys::IsCartesian())
        {
            for (int j = 0; j < len; j++)
                radius[i][j] = 1.0;
        }
        else
        {
            const Real xlo = grid_loc[i].lo(0) + (0.5 - radius_grow)*dxr;
            for (int j = 0; j < len; j++)
                radius[i][j] = xlo + j*dxr;
        }
    }
    //
    // Build volume and face area arrays.
    //
    geom.GetVolume(volume,grids,GEOM_GROW);
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        geom.GetFaceArea(area[dir],grids,dir,GEOM_GROW);
    }
}

//
// Reset the time levels to time (time) and timestep dt.
// This is done at the start of the timestep in the pressure iteration section.
//

void
NavierStokes::resetState (Real time,
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

//
// Set the time levels to time (time) and timestep dt.
//
void
NavierStokes::setTimeLevel (Real time,
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

//
// This function initializes the State and Pressure with data.
//

void
NavierStokes::initData ()
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::initData()");
    //
    // Initialize the state and the pressure.
    //
    int         ns       = NUM_STATE - BL_SPACEDIM;
    const Real* dx       = geom.CellSize();
    MultiFab&   S_new    = get_new_data(State_Type);
    MultiFab&   P_new    = get_new_data(Press_Type);
    const Real  cur_time = state[State_Type].curTime();

    for (MFIter snewmfi(S_new); snewmfi.isValid(); ++snewmfi)
    {
        BL_ASSERT(grids[snewmfi.index()] == snewmfi.validbox());

        P_new[snewmfi].setVal(0);

        const int  i    = snewmfi.index();
        const int* lo   = snewmfi.validbox().loVect();
        const int* hi   = snewmfi.validbox().hiVect();
        const int* s_lo = S_new[snewmfi].loVect();
        const int* s_hi = S_new[snewmfi].hiVect();
        const int* p_lo = P_new[snewmfi].loVect();
        const int* p_hi = P_new[snewmfi].hiVect();

        FORT_INITDATA (&level,&cur_time,lo,hi,&ns,
                       S_new[snewmfi].dataPtr(Xvel),
                       S_new[snewmfi].dataPtr(BL_SPACEDIM),
                       ARLIM(s_lo), ARLIM(s_hi),
                       P_new[snewmfi].dataPtr(),
                       ARLIM(p_lo), ARLIM(p_hi),
                       dx,grid_loc[i].lo(),grid_loc[i].hi() );
    }

    make_rho_prev_time();
    make_rho_curr_time();
    //
    // Initialize other types.
    //
    initDataOtherTypes();
    //
    // Initialize divU and dSdt.
    //
    if (have_divu)
    {
        const Real dt       = 1.0;
        const Real dtin     = -1.0; // Dummy value denotes initialization.
        const Real cur_time = state[Divu_Type].curTime();
        MultiFab&  Divu_new = get_new_data(Divu_Type);

        state[State_Type].setTimeLevel(cur_time,dt,dt);

        calc_divu(cur_time,dtin,Divu_new);

        if (have_dsdt)
            get_new_data(Dsdt_Type).setVal(0);
    }

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point) 
    {
        get_new_data(Dpdt_Type).setVal(0);
    }

    is_first_step_after_regrid = false;
    old_intersect_new          = grids;
}

//
// Fills a new level n with best level n and coarser data available.
//

void
NavierStokes::init (AmrLevel &old)
{
    NavierStokes* oldns     = (NavierStokes*) &old;
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
    for (FillPatchIterator fpi(old,S_new,0,cur_time,State_Type,0,NUM_STATE);
        fpi.isValid();
         ++fpi)
    {
        S_new[fpi.index()].copy(fpi());
    }
    //
    // Note: we don't need to worry here about using FillPatch because
    //       it will automatically use the "old dpdt" to interpolate,
    //       since we haven't yet defined a new pressure at the lower level.
    //
    for (FillPatchIterator fpi(old,P_new,0,cur_pres_time,Press_Type,0,1);
         fpi.isValid();
         ++fpi)
    {
        P_old[fpi.index()].copy(fpi());
        P_new[fpi.index()].copy(fpi());
    }

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point) 
    {
        MultiFab& Dpdt_new = get_new_data(Dpdt_Type);

        for (FillPatchIterator fpi(old,Dpdt_new,0,cur_pres_time,Dpdt_Type,0,1);
             fpi.isValid();
             ++fpi)
        {
            Dpdt_new[fpi.index()].copy(fpi());
        }
    }
    //
    // Get best divu and dSdt data.
    //
    if (have_divu)
    {
        MultiFab& Divu_new = get_new_data(Divu_Type);
        
        for (FillPatchIterator fpi(old,Divu_new,0,cur_time,Divu_Type,0,1);
             fpi.isValid();
             ++fpi)
        {
            Divu_new[fpi.index()].copy(fpi());
        }

        if (have_dsdt)
        {
            MultiFab& Dsdt_new = get_new_data(Dsdt_Type);

            for (FillPatchIterator fpi(old,Dsdt_new,0,cur_time,Dsdt_Type,0,1);
                 fpi.isValid();
                 ++fpi)
            {
                Dsdt_new[fpi.index()].copy(fpi());
            }
        }
    }

    old_intersect_new          = BoxLib::intersect(grids,oldns->boxArray());
    is_first_step_after_regrid = true;
}

void
NavierStokes::init ()
{
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& P_new = get_new_data(Press_Type);
    MultiFab& P_old = get_old_data(Press_Type);

    BL_ASSERT(level > 0);

    const Array<Real>& dt_amr = parent->dtLevel();
    Array<Real>        dt_new(level+1);

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
    NavierStokes& old       = getLevel(level-1);
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

    for (MFIter mfi(P_new); mfi.isValid(); ++mfi)
    {
        P_old[mfi].copy(P_new[mfi]);
    }
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

//
// ADVANCE FUNCTIONS
//

//
// This function ensures that the multifab registers and boundary
// flux registers needed for syncing the composite grid
//
//     u_mac, umacG, Vsync, Ssync, rhoavg, fr_adv, fr_visc
//
// are initialized to zero.  In general these quantities
// along with the pressure sync registers (sync_reg) and
// advective velocity registers (mac_reg) are compiled by first
// setting them to the coarse value acquired during a coarse timestep
// and then incrementing in the fine values acquired during the
// subcycled fine timesteps.  This compilation procedure occurs in
// different parts for different quantities
//
// * u_mac is set in predict_velocity and mac_project.
// * fr_adv, fr_visc are set in velocity_advect and scalar_advect
// * Vsync, Ssync are set in subcycled calls to post_timestep
// * mac_reg is set in mac_project
// * sync_reg is set in level_project
// * rhoavg, pavg are set in advance_setup and advance
//
// After these quantities have been compiled during a coarse
// timestep and subcycled fine timesteps.  The post_timestep function
// uses them to sync the fine and coarse levels.  If the coarse level
// is not the base level, post_timestep modifies the next coarsest levels
// registers appropriately.
//
// Note :: There is a little ambiguity as to which level owns the
// boundary flux registers.  The Multifab registers are quantities
// sized by the coarse level BoxArray and belong to the coarse level.
// The fine levels own the boundary registers, since they are sized by
// the boundaries of the fine level BoxArray.
//

void
NavierStokes::advance_setup (Real time,
                             Real dt,
                             int  iteration,
                             int  ncycle)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::advance_setup()");

    const int finest_level = parent->finestLevel();

    mac_projector->setup(level);
    //
    // Why are they defined here versus the constructor?
    //
    if (level < finest_level)
    {
        if (Vsync == 0)
            Vsync = new MultiFab(grids,BL_SPACEDIM,1);
        if (Ssync == 0)
            Ssync = new MultiFab(grids,NUM_STATE-BL_SPACEDIM,1);
        Vsync->setVal(0);
        Ssync->setVal(0);
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
    if (u_mac == 0)
    {
        u_mac = new MultiFab[BL_SPACEDIM];

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            BoxArray edge_grids(grids);
            edge_grids.surroundingNodes(dir);
            u_mac[dir].define(edge_grids,1,0,Fab_allocate);
        }
    }
    if (u_macG == 0)
    {
        u_macG = new MultiFab[BL_SPACEDIM];

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            BoxArray edge_grids(grids);
            edge_grids.surroundingNodes(dir).grow(1);
            u_macG[dir].define(edge_grids,1,0,Fab_allocate);
            u_macG[dir].setVal(1.e40);
        }
    }
    if (u_corr == 0)
    {
        u_corr = new MultiFab[BL_SPACEDIM];

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            BoxArray edge_grids(grids);
            edge_grids.surroundingNodes(dir).grow(1);
            u_corr[dir].define(edge_grids,1,0,Fab_allocate);
        }
    }
    //
    // Alloc MultiFab to hold advective update terms.
    //
    BL_ASSERT(aofs == 0);
    aofs = new MultiFab(grids,NUM_STATE,0);
    //
    // Alloc MultiFab to hold RHS for MAC projection.
    //
    BL_ASSERT(mac_rhs == 0);
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
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point) 
    {
        const Real new_press_time = .5 * (state[State_Type].prevTime() +
                                          state[State_Type].curTime());
        state[Press_Type].setNewTimeLevel(new_press_time);
    }

    make_rho_prev_time();
    //
    // If refRatio==4 to the next level coarser, and we're going to diffuse
    // scalars as SoverRho, we're going to need rho at 1/4 and 3/4 time there.
    // Make these things if need be.
    //
    if (level > 0)
    {
        bool needs_rho4 = false;

        if (parent->nCycle(level) == 4)
            for (int i = 0; i < NUM_STATE && !needs_rho4; ++i)
                needs_rho4 = (diffusionType[i] == Laplacian_SoverRho);

        if (needs_rho4)
        {
            NavierStokes&   clevel = getLevel(level-1);
            const BoxArray& cgrids = clevel.boxArray();
            const Real      ptime  = clevel.state[State_Type].prevTime();
            const Real      ctime  = clevel.state[State_Type].curTime();

            if (clevel.rho_qtime == 0)
            {
                const Real qtime = ptime + 0.25*(ctime-ptime);
                clevel.rho_qtime = new MultiFab(cgrids,1,1);
                FillPatchIterator fpi(clevel,*(clevel.rho_qtime),
                                      1,qtime,State_Type,Density,1);
                for ( ; fpi.isValid(); ++fpi)
                    (*clevel.rho_qtime)[fpi.index()].copy(fpi());
            }
            if (clevel.rho_tqtime == 0)
            {
                const Real tqtime = ptime + 0.75*(ctime-ptime);
                clevel.rho_tqtime = new MultiFab(cgrids,1,1);
                FillPatchIterator fpi(clevel,*(clevel.rho_tqtime),
                                      1,tqtime,State_Type,Density,1);
                for ( ; fpi.isValid(); ++fpi)
                    (*clevel.rho_tqtime)[fpi.index()].copy(fpi());
            }
        }
    }
    //
    // Calculate the time N viscosity and diffusivity
    //   Note: The viscosity and diffusivity at time N+1 are 
    //         initialized here to the time N values just to
    //         have something reasonable.
    //
    const Real prev_time = state[State_Type].prevTime();

    if (variable_vel_visc)
    {
        calcViscosity(prev_time,dt,iteration,ncycle);

        for (MFIter np1Mfi(*viscnp1_cc); np1Mfi.isValid(); ++np1Mfi)
        {
            (*viscnp1_cc)[np1Mfi].copy((*viscn_cc)[np1Mfi],0,0,1);
        }
    }

    if (variable_scal_diff)
    {
        const int num_diff = NUM_STATE-Density-1;

        calcDiffusivity(prev_time,dt,iteration,ncycle,Density+1,num_diff);

        for (MFIter np1Mfi(*diffnp1_cc); np1Mfi.isValid(); ++np1Mfi)
        {
            (*diffnp1_cc)[np1Mfi].copy((*diffn_cc)[np1Mfi],0,0,num_diff);
        }
    }
}

//
// Clean up after the advance function.
//

void
NavierStokes::advance_cleanup (Real dt,
                               int  iteration,
                               int  ncycle)
{
    if (level == parent->finestLevel())
    {
        delete [] u_mac;
        u_mac = 0;

        delete [] u_macG;
        u_macG = 0;

        delete [] u_corr;
        u_corr = 0;
    }
    delete mac_rhs;
    mac_rhs = 0;
    delete aofs;
    aofs = 0;
}


//
// Compute a timestep at a level. Return largest safe timestep.
//

Real
NavierStokes::advance (Real time,
                       Real dt,
                       int  iteration,
                       int  ncycle)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::advance()");

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "Advancing grids at level " << level
                  << " : starting time = "       << time
                  << " with dt = "               << dt << '\n';
    }
    advance_setup(time,dt,iteration,ncycle);

    //
    // Compute traced states for normal comp of velocity at half time level.
    //
    Real dummy   = 0.0;
    Real dt_test = predict_velocity(dt,dummy);

    //
    // Do MAC projection and update edge velocities.
    //

    if (do_mac_proj) 
    {
      mac_rhs = create_mac_rhs(time,dt);
      MultiFab& S_old  = get_old_data(State_Type);
      mac_project(time,dt,S_old,mac_rhs,have_divu);
    }

    //
    // Advect velocities.
    //
    if (do_mom_diff == 0) 
        velocity_advection(dt);
    //
    // Advect scalars.
    //
    const int first_scalar = Density;
    const int last_scalar  = first_scalar + NUM_SCALARS - 1;
    scalar_advection(dt,first_scalar,last_scalar);
    //
    // Update Rho.
    //
    scalar_update(dt,first_scalar,first_scalar);

    make_rho_curr_time();
    //
    // Advect momenta after rho^(n+1) has been created.
    //
    if (do_mom_diff == 1) 
        velocity_advection(dt);
    //
    // Add the advective and other terms to get scalars at t^{n+1}.
    //
    scalar_update(dt,first_scalar+1,last_scalar);
    //
    // S appears in rhs of the velocity update, so we better do it now.
    //
    if (have_divu)
    {
        calc_divu(time+dt,dt,get_new_data(Divu_Type));
        if (have_dsdt)
        {
            calc_dsdt(time,dt,get_new_data(Dsdt_Type));
            if (initial_step)
                MultiFab::Copy(get_old_data(Dsdt_Type),
                               get_new_data(Dsdt_Type),0,0,1,0);
        }
    }
    //
    // Add the advective and other terms to get velocity at t^{n+1}.
    //
    velocity_update(dt);
    //
    // Clean up after the predicted value at t^n+1.
    // Estimate new timestep from umac cfl.
    //
    advance_cleanup(dt,iteration,ncycle);
    //
    // Increment rho average.
    //
    if (!initial_step && level > 0)
        incrRhoAvg((iteration==ncycle ? 0.5 : 1.0) / Real(ncycle));
    //
    // Do a level project to update the pressure and velocity fields.
    //
    if (!initial_step)
    {
        if (projector)
            level_projector(dt,time,iteration);
        if (level > 0 && iteration == 1)
           p_avg->setVal(0);
    }
    return dt_test;  // Return estimate of best new timestep.
}

MultiFab*
NavierStokes::create_mac_rhs (Real time, Real dt)
{
      MultiFab* mac_rhs = getDivCond(0,time);
      MultiFab* dsdt = getDsdt(0,time);

      for (MFIter mfi(*mac_rhs); mfi.isValid(); ++mfi)
      {
          (*dsdt)[mfi].mult(.5*dt);
          (*mac_rhs)[mfi].plus((*dsdt)[mfi]);
      }

      delete dsdt;
      return mac_rhs;
}

MultiFab*
NavierStokes::create_mac_rhs_grown (int  nGrow,
                                    Real time,
                                    Real dt)
{
    //
    // We do all this just to fill the ghost cells.
    //
    MultiFab* mac_rhs_grown = 0;

    if (nGrow > 0) 
    {
        mac_rhs_grown = getDivCond(nGrow,time);
        MultiFab* dsdt_grown = getDsdt(nGrow,time);
        for (MFIter mfi(*mac_rhs_grown); mfi.isValid(); ++mfi)
        {
            (*dsdt_grown)[mfi].mult(.5*dt);
            (*mac_rhs_grown)[mfi].plus((*dsdt_grown)[mfi]);
        }
        delete dsdt_grown;
    }
    //
    // Now we copy from the MultiFab mac_rhs which has no ghost cells.
    //
    for (MFIter mfi(*mac_rhs_grown); mfi.isValid(); ++mfi) 
    {
        (*mac_rhs_grown)[mfi].copy((*mac_rhs)[mfi],grids[mfi.index()]);
    }

    if (nGrow > 0) 
        mac_rhs_grown->FillBoundary();

    return mac_rhs_grown;
}

void
NavierStokes::mac_project (Real      time,
                           Real      dt,
                           MultiFab& Sold, 
                           MultiFab* divu,
                           int       have_divu)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::mac_project()");

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... mac_projection\n";

    const Real strt_time = ParallelDescriptor::second();

    mac_projector->mac_project(level,u_mac,Sold,dt,time,*divu,have_divu);

    create_umac_grown();

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_time = ParallelDescriptor::second() - strt_time;

    ParallelDescriptor::ReduceRealMax(run_time,IOProc);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "NavierStokes:mac_project(): lev: "
                  << level
                  << ", time: " << run_time << std::endl;
    }
}

void
NavierStokes::level_projector (Real dt,
                               Real time,
                               int  iteration)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::level_projector()");

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

    Array<int*>         sync_bc(grids.size());
    Array< Array<int> > sync_bc_array(grids.size());

    for (int i = 0; i < grids.size(); i++)
    {
        sync_bc_array[i] = getBCArray(State_Type,i,Xvel,BL_SPACEDIM);
        sync_bc[i]       = sync_bc_array[i].dataPtr();
    }

    int        crse_dt_ratio  = (level > 0) ? parent->nCycle(level) : -1;
    const Real cur_pres_time  = state[Press_Type].curTime();
    const Real prev_pres_time = state[Press_Type].prevTime();

    projector->level_project(level,time,dt,cur_pres_time,prev_pres_time,
                             geom,U_old,U_new,P_old,P_new,
                             get_rho_half_time(),crse_ptr,sync_reg,
                             crse_dt_ratio,sync_bc.dataPtr(),iteration,
                             have_divu,Divu_Type);

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
        calcDpdt();
}

void
NavierStokes::make_rho_prev_time ()
{
    const Real prev_time = state[State_Type].prevTime();

    FillPatchIterator fpi(*this,*rho_ptime,1,prev_time,State_Type,Density,1);

    for ( ; fpi.isValid(); ++fpi)
    {
        (*rho_ptime)[fpi.index()].copy(fpi());
    }
}

void
NavierStokes::make_rho_curr_time ()
{
    const Real curr_time = state[State_Type].curTime();

    FillPatchIterator fpi(*this,*rho_ctime,1,curr_time,State_Type,Density,1);

    for ( ; fpi.isValid(); ++fpi)
    {
        (*rho_ctime)[fpi.index()].copy(fpi());
    }
}

MultiFab*
NavierStokes::get_rho_half_time ()
{
    //
    // Fill it in when needed ...
    //
    for (MFIter mfi(*rho_half); mfi.isValid(); ++mfi)
    {
        (*rho_half)[mfi].copy((*rho_ptime)[mfi]);
        (*rho_half)[mfi] += (*rho_ctime)[mfi];
        (*rho_half)[mfi].mult(.5);
    }

    return rho_half;
}

const MultiFab&
NavierStokes::get_rho (Real time)
{
    const TimeLevel whichTime = which_time(State_Type,time);

    if (whichTime == AmrOldTime)
    {
        return *rho_ptime;
    }
    else if (whichTime == AmrNewTime)
    {
        return *rho_ctime;
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
        return *get_rho_half_time();
    }
    else
    {
        BoxLib::Error("NavierStokes::get_rho(): bad time");

        return *rho_ptime; // Got to return something to shut up compiler.
    }
}

//
// Predict the edge velocities which go into forming u_mac.  This
// function also returns an estimate of dt for use in variable timesteping.
//

Real
NavierStokes::predict_velocity (Real  dt,
                                Real& comp_cfl)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::predict_velocity()");

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... predict edge velocities\n";
    //
    // Get simulation parameters.
    //
    const int   nComp          = BL_SPACEDIM;
    const Real* dx             = geom.CellSize();
    const Real  prev_time      = state[State_Type].prevTime();
    const Real  prev_pres_time = state[Press_Type].prevTime();
    //
    // Compute viscous terms at level n.
    // Ensure reasonable values in 1 grow cell.  Here, do extrap for
    // c-f/phys boundary, since we have no interpolator fn, also,
    // preserve extrap for corners at periodic/non-periodic intersections.
    //
    MultiFab visc_terms(grids,nComp,1);

    if (be_cn_theta != 1.0)
    {
	getViscTerms(visc_terms,Xvel,nComp,prev_time);
    }
    else
    {
	visc_terms.setVal(0);
    }
    //
    // Set up the timestep estimation.
    //
    Real cflgrid,u_max[3];
    Real cflmax = 1.0e-10;
    comp_cfl    = (level == 0) ? cflmax : comp_cfl;

    FArrayBox tforces;
    Array<int> bndry[BL_SPACEDIM];

    MultiFab Gp(grids,BL_SPACEDIM,1);
    getGradP(Gp, prev_pres_time);
    
    FArrayBox* null_fab = 0;
    for (FillPatchIterator P_fpi(*this,get_old_data(Press_Type),1,prev_pres_time,Press_Type,0,1),
             U_fpi(*this,visc_terms,HYP_GROW,prev_time,State_Type,Xvel,BL_SPACEDIM);
         U_fpi.isValid() && P_fpi.isValid();
         ++U_fpi, ++P_fpi)
    {
        const int i = U_fpi.index();

        getForce(tforces,i,1,Xvel,BL_SPACEDIM,
#ifdef GENGETFORCE
		 prev_time,
#endif		 
		 (*rho_ptime)[i]);
        //
        // Test velocities, rho and cfl.
        //
        cflgrid  = godunov->test_u_rho(U_fpi(),(*rho_ptime)[i],grids[i],dx,dt,u_max);
        cflmax   = std::max(cflgrid,cflmax);
        comp_cfl = std::max(cflgrid,comp_cfl);
        //
        // Compute the total forcing.
        //
        godunov->Sum_tf_gp_visc(tforces,visc_terms[i],Gp[i],(*rho_ptime)[i]);

        D_TERM(bndry[0] = getBCArray(State_Type,i,0,1);,
               bndry[1] = getBCArray(State_Type,i,1,1);,
               bndry[2] = getBCArray(State_Type,i,2,1);)

        godunov->Setup(grids[i], dx, dt, 1,
                       *null_fab, bndry[0].dataPtr(),
                       *null_fab, bndry[1].dataPtr(),
#if (BL_SPACEDIM == 3)                         
                       *null_fab, bndry[2].dataPtr(),
#endif
                       U_fpi(), (*rho_ptime)[i], tforces);

        godunov->ComputeUmac(grids[i], dx, dt, 
                             u_mac[0][i], bndry[0].dataPtr(),
                             u_mac[1][i], bndry[1].dataPtr(),
#if (BL_SPACEDIM == 3)
                             u_mac[2][i], bndry[2].dataPtr(),
#endif
                             U_fpi(), tforces);
    }

    Real tempdt = std::min(change_max,cfl/cflmax);

    ParallelDescriptor::ReduceRealMin(tempdt);

    return dt*tempdt;
}

//
// This routine advects the velocities
//
void
NavierStokes::velocity_advection (Real dt)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::velocity_advection()");

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        if (do_mom_diff == 0) 
        {
            std::cout << "... advect velocities\n";
        }
        else
        {
            if (predict_mom_together == 0)
            {
                std::cout << "Must set predict_mom_together == 1 in NavierStokes." << std::endl;
                exit(0);
            }
            std::cout << "... advect momenta\n";
        }
    }

    const int   finest_level   = parent->finestLevel();
    const Real* dx             = geom.CellSize();
    const Real  prev_time      = state[State_Type].prevTime();
    const Real  prev_pres_time = state[Press_Type].prevTime();
    //
    // Compute viscosity components.
    //
    MultiFab visc_terms(grids,BL_SPACEDIM,1);

    if (be_cn_theta != 1.0)
        getViscTerms(visc_terms,Xvel,BL_SPACEDIM,prev_time);
    else
        visc_terms.setVal(0,1);

    Array<int> bndry[BL_SPACEDIM];
    FArrayBox xflux, yflux, zflux, divu, tforces;

    int nGrowF = 1;
    MultiFab* divu_fp = create_mac_rhs_grown(nGrowF,prev_time,dt);

#if 0
    MultiFab* divu_fp = getDivCond(nGrowF,prev_time);
    MultiFab* dsdt = getDsdt(nGrowF,prev_time);
    for (MFIter dsdtmfi(*dsdt); dsdtmfi.isValid(); ++dsdtmfi)
    {
       (*dsdt)[dsdtmfi].mult(.5*dt);
       (*divu_fp)[dsdtmfi].plus((*dsdt)[dsdtmfi]);
    }
    delete dsdt;
#endif

    MultiFab Gp(grids,BL_SPACEDIM,1);
    getGradP(Gp, prev_pres_time);
    //
    // Compute the advective forcing.
    //
    FArrayBox S;

    for (FillPatchIterator P_fpi(*this,get_old_data(Press_Type),1,prev_pres_time,Press_Type,0,1),
             U_fpi(*this,visc_terms,HYP_GROW,prev_time,State_Type,Xvel,BL_SPACEDIM),
             Rho_fpi(*this,visc_terms,HYP_GROW,prev_time,State_Type,Density,1);
         U_fpi.isValid() && P_fpi.isValid() && Rho_fpi.isValid(); 
         ++U_fpi, ++P_fpi, ++Rho_fpi)
    {
        //
        // Since all the MultiFabs are on same grid we'll just use indices.
        //
        const int i = U_fpi.index();

        getForce(tforces,i,1,Xvel,BL_SPACEDIM,
#ifdef GENGETFORCE
		 prev_time,
#endif		 
		 (*rho_ptime)[i]);

        godunov->Sum_tf_gp_visc(tforces,visc_terms[i],Gp[i],(*rho_ptime)[i]);

        D_TERM(bndry[0] = getBCArray(State_Type,i,0,1);,
               bndry[1] = getBCArray(State_Type,i,1,1);,
               bndry[2] = getBCArray(State_Type,i,2,1);)

        godunov->Setup(grids[i], dx, dt, 0,
                       xflux, bndry[0].dataPtr(), yflux, bndry[1].dataPtr(),
#if (BL_SPACEDIM == 3)                          
                       zflux, bndry[2].dataPtr(),
#endif
                       U_fpi(),(*rho_ptime)[i],tforces);
        //
        // Loop over the velocity components.
        //
        S.resize(U_fpi().box(),BL_SPACEDIM);
        S.copy(U_fpi(),0,0,BL_SPACEDIM);

        for (int comp = 0 ; comp < BL_SPACEDIM ; comp++ )
        {
            int use_conserv_diff = (advectionType[comp] == Conservative) 
                                                             ? true : false;
            if (do_mom_diff == 1)
            {
                S.mult(Rho_fpi(),S.box(),S.box(),0,comp,1);
                tforces.mult((*rho_ptime)[i],tforces.box(),tforces.box(),0,comp,1);
            }
            FArrayBox divu_dummy;
            godunov->AdvectState(grids[i], dx, dt, 
                                 area[0][i], u_mac[0][i], xflux,
                                 area[1][i], u_mac[1][i], yflux,
#if (BL_SPACEDIM == 3)                       
                                 area[2][i], u_mac[2][i], zflux,
#endif
                                 U_fpi(), S, tforces, (*divu_fp)[i], comp,
                                 (*aofs)[i],comp,use_conserv_diff,
                                 comp,bndry[comp].dataPtr(),PRE_MAC,volume[i]);
            //
            // Get fluxes for diagnostics and refluxing.
            //
            pullFluxes(i,comp,1,xflux,yflux,zflux,dt);
        }
    }
    delete divu_fp;
    //
    // pullFluxes() contains CrseInit() calls -- complete the process.
    //
    if (do_reflux && level < finest_level)
        getAdvFluxReg(level+1).CrseInitFinish();
}

//
// This routine advects the scalars
//

void
NavierStokes::scalar_advection (Real dt,
                                int  fscalar,
                                int  lscalar)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::scalar_advection()");

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... advect scalars\n";
    //
    // Get simulation parameters.
    //
    const int   num_scalars    = lscalar - fscalar + 1;
    const Real* dx             = geom.CellSize();
    const Real  prev_time      = state[State_Type].prevTime();
    const Real  prev_pres_time = state[Press_Type].prevTime();
    //
    // Get the viscous terms.
    //
    MultiFab visc_terms(grids,num_scalars,1);

    if (be_cn_theta != 1.0)
        getViscTerms(visc_terms,fscalar,num_scalars,prev_time);
    else
        visc_terms.setVal(0,1);
    //
    // Set up the grid loop.
    //
    FArrayBox xflux, yflux, zflux, tforces, tvelforces;

    MultiFab Gp(grids,BL_SPACEDIM,1);
    getGradP(Gp, prev_pres_time);

    MultiFab vel_visc_terms;

    const int use_forces_in_trans = godunov->useForcesInTrans();

    if (use_forces_in_trans)
    {
      vel_visc_terms.define(grids,BL_SPACEDIM,1,Fab_allocate);

      if (be_cn_theta != 1.0)
          getViscTerms(vel_visc_terms,Xvel,BL_SPACEDIM,prev_time);
      else
          vel_visc_terms.setVal(0,1);
    }
    Array<int> state_bc, bndry[BL_SPACEDIM];

    MultiFab* divu_fp = getDivCond(1,prev_time);

    int nGrowF = 1;
    MultiFab* dsdt = getDsdt(nGrowF,prev_time);
    for (MFIter dsdtmfi(*dsdt); dsdtmfi.isValid(); ++dsdtmfi)
    {
       (*dsdt)[dsdtmfi].mult(.5*dt);
       (*divu_fp)[dsdtmfi].plus((*dsdt)[dsdtmfi]);
    }
    delete dsdt;

    //
    // Compute the advective forcing.
    //
    for (FillPatchIterator P_fpi(*this,get_old_data(Press_Type),1,prev_pres_time,Press_Type,0,1),
             U_fpi(*this,visc_terms,HYP_GROW,prev_time,State_Type,Xvel,BL_SPACEDIM),
             S_fpi(*this,visc_terms,HYP_GROW,prev_time,State_Type,fscalar,num_scalars);
         U_fpi.isValid() && S_fpi.isValid() && P_fpi.isValid();
         ++U_fpi, ++S_fpi, ++P_fpi)
    {
        const int i = U_fpi.index();

        getForce(tforces,i,1,fscalar,num_scalars,
#ifdef GENGETFORCE
		 prev_time,
#endif		 
		 (*rho_ptime)[i]);
        
        if (use_forces_in_trans)
        {
            getForce(tvelforces,i,1,Xvel,BL_SPACEDIM,
#ifdef GENGETFORCE
		     prev_time,
#endif		 
		     (*rho_ptime)[i]);

            godunov->Sum_tf_gp_visc(tvelforces,vel_visc_terms[i],Gp[i],(*rho_ptime)[i]);
        }

        D_TERM(bndry[0] = getBCArray(State_Type,i,0,1);,
               bndry[1] = getBCArray(State_Type,i,1,1);,
               bndry[2] = getBCArray(State_Type,i,2,1);)

        godunov->Setup(grids[i], dx, dt, 0,
                       xflux, bndry[0].dataPtr(),
                       yflux, bndry[1].dataPtr(),
#if (BL_SPACEDIM == 3)                         
                       zflux, bndry[2].dataPtr(),
#endif
                       U_fpi(),(*rho_ptime)[i],tvelforces);
        //
        // Loop over the scalar components.
        //
        for (int comp = 0 ; comp < num_scalars ; comp++)
        {
            int state_ind = fscalar + comp;
            //
            // Compute total forcing.
            //
            int use_conserv_diff = (advectionType[state_ind] == Conservative)
                                                             ? true : false;
            AdvectionScheme adv_scheme = PRE_MAC;

            if (adv_scheme == PRE_MAC) {
              godunov->Sum_tf_divu_visc(S_fpi(),tforces,comp,1,visc_terms[i],
                                        comp,(*divu_fp)[i],(*rho_ptime)[i],
                                        use_conserv_diff);
            } else {
              FArrayBox junkDivu(tforces.box(),1);
              junkDivu.setVal(0.);
              godunov->Sum_tf_divu_visc(S_fpi(),tforces,comp,1,visc_terms[i],
                                        comp,junkDivu,(*rho_ptime)[i],
                                        use_conserv_diff);
            }
            //
            // Advect scalar.
            //
            state_bc = getBCArray(State_Type,i,state_ind,1);

            godunov->AdvectState(grids[i], dx, dt, 
                                 area[0][i], u_mac[0][i], xflux,
                                 area[1][i], u_mac[1][i], yflux,
#if (BL_SPACEDIM == 3)                        
                                 area[2][i], u_mac[2][i], zflux,
#endif
                                 U_fpi(),S_fpi(),tforces,(*divu_fp)[i],comp,
                                 (*aofs)[i],state_ind,use_conserv_diff,
                                 state_ind,state_bc.dataPtr(),adv_scheme,volume[i]);
            //
            // Get the fluxes for refluxing and diagnostic purposes.
            //
            pullFluxes(i, state_ind, 1, xflux, yflux, zflux, dt);
        }
    }
    delete divu_fp;
    //
    // pullFluxes() contains CrseInit() calls -- complete the process.
    //
    if (do_reflux && level < parent->finestLevel())
        getAdvFluxReg(level+1).CrseInitFinish();
}

//
// This subroutine updates the scalars, before the velocity update
// and the level projection
//
// AT this point in time, all we know is psi^n, rho^n+1/2, and the
// general forcing terms at t^n, and after solving in this routine
// viscous forcing at t^n+1/2.  Note, unless more complicated logic
// is invoked earlier, we do not have any estimate of general forcing
// terms at t^n+1/2.
//

void
NavierStokes::scalar_update (Real dt,
                             int  first_scalar,
                             int  last_scalar)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::scalar_update()");

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... update scalars\n";

    scalar_advection_update(dt, first_scalar, last_scalar);

    bool do_any_diffuse = false;
    for (int sigma = first_scalar; sigma <= last_scalar; sigma++)
        if (is_diffusive[sigma]) do_any_diffuse = true;

    if (do_any_diffuse)
      scalar_diffusion_update(dt, first_scalar, last_scalar);
}

void
NavierStokes::scalar_advection_update (Real dt,
                                       int  first_scalar,
                                       int  last_scalar)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::scalar_advection_update()");

    MultiFab&  S_old    = get_old_data(State_Type);
    MultiFab&  S_new    = get_new_data(State_Type);
    MultiFab&  Aofs     = *aofs;
    const Real halftime = 0.5*(state[State_Type].curTime()+state[State_Type].prevTime());
    const Real prev_time= state[State_Type].prevTime();
    Array<int> state_bc;
    FArrayBox  tforces;
    //
    // Compute inviscid estimate of scalars.
    // (do rho separate, as we do not have rho at new time yet)
    //
    int sComp = first_scalar;

    if (sComp == Density)
    {
        for (MFIter S_oldmfi(S_old); S_oldmfi.isValid(); ++S_oldmfi)
        {
            const int i = S_oldmfi.index();
            tforces.resize(grids[i],1);
            tforces.setVal(0);
            godunov->Add_aofs_tf(S_old[S_oldmfi],S_new[S_oldmfi],Density,1,
                                 Aofs[S_oldmfi],Density,tforces,0,grids[i],dt);
        }
#if 1
        // Call ScalMinMax to avoid overshoots in density
        if (do_denminmax) {
          // Must do FillPatch here instead of MF iterator because we need the
          //  boundary values in the old data (especially at inflow)

          int index_new_s   = Density;
          int index_new_rho = Density;
          int index_old_s   = index_new_s   - Density;
          int index_old_rho = index_new_rho - Density;

          for (FillPatchIterator
               S_fpi(*this,S_old,1,prev_time,State_Type,Density,1);
               S_fpi.isValid(); ++S_fpi)
           {
              const int i = S_fpi.index();
              state_bc = getBCArray(State_Type,i,Density,1);
              godunov->ConservativeScalMinMax(S_fpi(),S_new[S_fpi],
                                              index_old_s, index_old_rho,
                                              index_new_s, index_new_rho,
                                              state_bc.dataPtr(),grids[i]);
           }
        }
#endif
        ++sComp;
    }

    if (sComp <= last_scalar)
    {
        const MultiFab& rho_halftime = *get_rho_half_time();
        for (MFIter Rho_mfi(rho_halftime); Rho_mfi.isValid(); ++Rho_mfi)
        {
            const int i = Rho_mfi.index();
            for (int sigma = sComp; sigma <= last_scalar; sigma++)
            {
                getForce(tforces,i,0,sigma,1,
#ifdef GENGETFORCE
			 halftime,
#endif		 
			 rho_halftime[Rho_mfi]);
                godunov->Add_aofs_tf(S_old[Rho_mfi],S_new[Rho_mfi],sigma,1,
                                     Aofs[Rho_mfi],sigma,tforces,0,grids[i],dt);
            }
        }
    }

#if 1
    // Call ScalMinMax to avoid overshoots in the scalars 
    if ( do_scalminmax && (sComp <= last_scalar) )
    {
        int num_scalars = last_scalar - Density + 1;

        // Must do FillPatch here instead of MF iterator because we need the
        //  boundary values in the old data (especially at inflow)

        for (FillPatchIterator
             S_fpi(*this,S_old,1,prev_time,State_Type,Density,num_scalars);
             S_fpi.isValid(); ++S_fpi)
        {
            const int i = S_fpi.index();
            for (int sigma = sComp; sigma <= last_scalar; sigma++)
            {
              int index_new_s   = sigma;
              int index_new_rho = Density;
              int index_old_s   = index_new_s   - Density;
              int index_old_rho = index_new_rho - Density;
              state_bc = getBCArray(State_Type,i,sigma,1);
              if (advectionType[sigma] == Conservative)
              {
                  godunov->ConservativeScalMinMax(S_fpi(),S_new[S_fpi],
                                                  index_old_s, index_old_rho,
                                                  index_new_s, index_new_rho,
                                                  state_bc.dataPtr(),grids[i]);
              }
              else if (advectionType[sigma] == NonConservative)
              {
                  godunov->ConvectiveScalMinMax(S_fpi(),S_new[S_fpi],index_old_s,sigma,
                                                state_bc.dataPtr(),grids[i]);
              }
            }
        }
    }
#endif
}

void
NavierStokes::scalar_diffusion_update (Real dt,
                                       int  first_scalar,
                                       int  last_scalar)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::scalar_diffusion_update()");

    MultiFab** fluxSCn;
    MultiFab** fluxSCnp1;
    const int nGrow = 0;
    const int nComp = 1;
    diffusion->allocFluxBoxesLevel(fluxSCn,  nGrow,nComp);
    diffusion->allocFluxBoxesLevel(fluxSCnp1,nGrow,nComp);
    const MultiFab* Rh = get_rho_half_time();

    for (int sigma = first_scalar; sigma <= last_scalar; sigma++)
    {
        if (is_diffusive[sigma])
        {
            int        rho_flag    = 0;
            MultiFab*  delta_rhs   = 0;
            MultiFab*  alpha       = 0;
            MultiFab** cmp_diffn   = 0;
            MultiFab** cmp_diffnp1 = 0;

            if (variable_scal_diff)
            {
                Real diffTime = state[State_Type].prevTime();
                diffusion->allocFluxBoxesLevel(cmp_diffn, 0, 1);
                getDiffusivity(cmp_diffn, diffTime, sigma, 0, 1);

                diffTime = state[State_Type].curTime();
                diffusion->allocFluxBoxesLevel(cmp_diffnp1, 0, 1);
                getDiffusivity(cmp_diffnp1, diffTime, sigma, 0, 1);
            }

            diffuse_scalar_setup(dt, sigma, &rho_flag, 
                                 delta_rhs, alpha, cmp_diffn, cmp_diffnp1);

            diffusion->diffuse_scalar(dt,sigma,be_cn_theta,Rh,
                                      rho_flag,fluxSCn,fluxSCnp1,0,delta_rhs,
                                      alpha,cmp_diffn,cmp_diffnp1);
            if (variable_scal_diff)
            {
                diffusion->removeFluxBoxesLevel(cmp_diffn);
                diffusion->removeFluxBoxesLevel(cmp_diffnp1);
            }

            delete delta_rhs;
            delete alpha;
            //
            // Increment the viscous flux registers
            //
            if (do_reflux)
            {
                FArrayBox fluxtot;

                for (int d = 0; d < BL_SPACEDIM; d++)
                {
                    for (MFIter fmfi(*fluxSCn[d]); fmfi.isValid(); ++fmfi)
                    {
                        const Box& ebox = (*fluxSCn[d])[fmfi].box();
                        fluxtot.resize(ebox,nComp);
                        fluxtot.copy((*fluxSCn[d])[fmfi],ebox,0,ebox,0,nComp);
                        fluxtot.plus((*fluxSCnp1[d])[fmfi],ebox,0,0,nComp);
                        if (level < parent->finestLevel())
                            getLevel(level+1).getViscFluxReg().CrseInit(fluxtot,ebox,
                                                                        d,0,sigma,
                                                                        nComp,-dt);

                        if (level > 0)
                            getViscFluxReg().FineAdd((*fluxSCn[d])[fmfi],d,fmfi.index(),
                                                     0,sigma,nComp,dt);
                    }
                }
                if (level < parent->finestLevel())
                    getLevel(level+1).getViscFluxReg().CrseInitFinish();
            }
        }
    }
    diffusion->removeFluxBoxesLevel(fluxSCn);
    diffusion->removeFluxBoxesLevel(fluxSCnp1);
}

void
NavierStokes::diffuse_scalar_setup (Real        dt,
                                    int         sigma,
                                    int*        rho_flag,
                                    MultiFab*&  delta_rhs,
                                    MultiFab*&  alpha,
                                    MultiFab**& diffn,
                                    MultiFab**& diffnp1)
{
    (*rho_flag) = Diffusion::set_rho_flag(diffusionType[sigma]);
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
NavierStokes::velocity_update (Real dt)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::velocity_update()");

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        if (do_mom_diff == 0) 
        {
            std::cout << "... update velocities \n";
        }
        else
        {
            std::cout << "... update momenta \n";
        }
    }

    velocity_advection_update(dt);

    if (!initial_iter)
        velocity_diffusion_update(dt);
    else
        initial_velocity_diffusion_update(dt);
}

void
NavierStokes::velocity_advection_update (Real dt)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::velocity_advection_update()");

    FArrayBox  tforces;
    MultiFab&  U_old          = get_old_data(State_Type);
    MultiFab&  U_new          = get_new_data(State_Type);
    MultiFab&  P_old          = get_old_data(Press_Type);
    MultiFab&  Aofs           = *aofs;
    const Real prev_time      = state[State_Type].prevTime();
    const Real curr_time      = state[State_Type].curTime();
    const Real prev_pres_time = state[Press_Type].prevTime();
    const Real half_time      = 0.5*(curr_time+prev_time);

    MultiFab Gp(grids,BL_SPACEDIM,1);
    getGradP(Gp, prev_pres_time);

    MultiFab& halftime = *get_rho_half_time();
    MFIter    Rhohalf_mfi(halftime);
    FArrayBox S;

    for (FillPatchIterator P_fpi(*this,P_old,0,prev_pres_time,Press_Type,0,1);
         Rhohalf_mfi.isValid() && P_fpi.isValid();
         ++Rhohalf_mfi, ++P_fpi)
    {
        const int i = Rhohalf_mfi.index();

	getForce(tforces,i,0,Xvel,BL_SPACEDIM,
#ifdef GENGETFORCE
		 half_time,
#endif		 
		 halftime[i]);
        //
        // Do following only at initial iteration--per JBB.
        //
        if (initial_iter && is_diffusive[Xvel])
            tforces.setVal(0);

        S.resize(U_old[i].box(),BL_SPACEDIM);
        S.copy(U_old[i],0,0,BL_SPACEDIM);

        if (do_mom_diff == 1)
        {
            for (int d = 0; d < BL_SPACEDIM; d++)
            {
                Gp[i].mult(halftime[i],grids[i],grids[i],0,d,1);
                tforces.mult(halftime[i],grids[i],grids[i],0,d,1);
                S.mult((*rho_ptime)[i],grids[i],grids[i],0,d,1);
            }
        }

        godunov->Add_aofs_tf_gp(S,U_new[i],Aofs[i],tforces,
                                Gp[i],halftime[i],grids[i],dt);
        if (do_mom_diff == 1)
        {
            for (int d = 0; d < BL_SPACEDIM; d++)
                U_new[i].divide((*rho_ctime)[i],grids[i],grids[i],0,d,1);
        }
    }
}

void
NavierStokes::velocity_diffusion_update (Real dt)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::velocity_diffusion_update()");

    //
    // Compute the viscous forcing.
    // Do following except at initial iteration.
    //
    MultiFab& U_old = get_old_data(State_Type);
    MultiFab& U_new = get_new_data(State_Type);

    if (is_diffusive[Xvel])
    {
        int rho_flag = (do_mom_diff == 0) ? 1 : 3;

        MultiFab* delta_rhs = 0;
        if (S_in_vel_diffusion && have_divu)
        {
            delta_rhs = new MultiFab(grids,BL_SPACEDIM,0);
            delta_rhs->setVal(0);
        }

        MultiFab** loc_viscn   = 0;
        MultiFab** loc_viscnp1 = 0;

        if (variable_vel_visc)
        {
            Real viscTime = state[State_Type].prevTime();
            diffusion->allocFluxBoxesLevel(loc_viscn, 0, 1);
            getViscosity(loc_viscn, viscTime);

            viscTime = state[State_Type].curTime();
            diffusion->allocFluxBoxesLevel(loc_viscnp1, 0, 1);
            getViscosity(loc_viscnp1, viscTime);
        }

        diffuse_velocity_setup(dt, delta_rhs, loc_viscn, loc_viscnp1);

        diffusion->diffuse_velocity(dt,be_cn_theta,get_rho_half_time(),rho_flag,
                                    delta_rhs,loc_viscn,loc_viscnp1);

        if (variable_vel_visc)
        {
            diffusion->removeFluxBoxesLevel(loc_viscn);
            diffusion->removeFluxBoxesLevel(loc_viscnp1);
        }

        delete delta_rhs;
    }
}

void
NavierStokes::diffuse_velocity_setup (Real       dt,
                                      MultiFab*& delta_rhs,
                                      MultiFab**& viscn,
                                      MultiFab**& viscnp1)
{
    if (S_in_vel_diffusion && have_divu)
    {
        //
        // Include div mu S*I terms in rhs
        //  (i.e. make nonzero delta_rhs to add into RHS):
        //
        // The scalar and tensor solvers incorporate the relevant pieces of
        //  of Div(tau), provided the flow is divergence-free.  However, if
        //  Div(U) =/= 0, there is an additional piece not accounted for,
        //  which is of the form A.Div(U).  For constant viscosity, Div(tau)_i
        //  = Lapacian(U_i) + mu/3 d[Div(U)]/dx_i.  For mu not constant,
        //  Div(tau)_i = d[ mu(du_i/dx_j + du_j/dx_i) ]/dx_i - 2mu/3 d[Div(U)]/dx_i
        //
        // As a convenience, we treat this additional term as a "source" in
        // the diffusive solve, computing Div(U) in the "normal" way we
        // always do--via a call to calc_divu.  This routine computes delta_rhs
        // if necessary, and stores it as an auxilliary rhs to the viscous solves.
        // This is a little strange, but probably not bad.
        //
        const Real time = state[State_Type].prevTime();

        MultiFab divmusi(grids,BL_SPACEDIM,0);

        if (!variable_vel_visc)
        {
            diffusion->compute_divmusi(time,visc_coef[Xvel],divmusi);
            divmusi.mult((1./3.)*(1.0-be_cn_theta),0,BL_SPACEDIM,0);
            (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);

            diffusion->compute_divmusi(time+dt,visc_coef[Xvel],divmusi);
            divmusi.mult((1./3.)*be_cn_theta,0,BL_SPACEDIM,0);
            (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);
        }
        else
        {
            diffusion->compute_divmusi(time,viscn,divmusi);
            divmusi.mult((-2./3.)*(1.0-be_cn_theta),0,BL_SPACEDIM,0);
            (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);

            diffusion->compute_divmusi(time+dt,viscnp1,divmusi);
            divmusi.mult((-2./3.)*be_cn_theta,0,BL_SPACEDIM,0);
            (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);
        }
    }
}

void
NavierStokes::initial_velocity_diffusion_update (Real dt)
{
    //
    // Do following only at initial iteration.
    //
    if (is_diffusive[Xvel])
    {
        MultiFab&  U_old          = get_old_data(State_Type);
        MultiFab&  U_new          = get_new_data(State_Type);
        MultiFab&  Aofs           = *aofs;
        const int  nComp          = BL_SPACEDIM;
        const Real prev_time      = state[State_Type].prevTime();
        const Real prev_pres_time = state[Press_Type].prevTime();

        MultiFab Gp(grids,BL_SPACEDIM,1);
        getGradP(Gp, prev_pres_time);

	MultiFab visc_terms(grids,nComp,1);

	if (be_cn_theta != 1.0)
        {
	    getViscTerms(visc_terms,Xvel,nComp,prev_time);
        }
        else
	{
	    visc_terms.setVal(0);
	}

        FArrayBox tforces;
        FArrayBox S;
        //
        // Update U_new with viscosity.
        //
        MultiFab* Rh = get_rho_half_time();

        for (FillPatchIterator P_fpi(*this,get_old_data(Press_Type),0,prev_pres_time,Press_Type,0,1);
             P_fpi.isValid();
             ++P_fpi)
        {
            const int i = P_fpi.index();

            getForce(tforces,i,0,Xvel,BL_SPACEDIM,
#ifdef GENGETFORCE
		     prev_time,
#endif		 
		     (*rho_ptime)[i]);

            godunov->Sum_tf_gp_visc(tforces,visc_terms[i],Gp[i],(*Rh)[i]);

            S.resize(U_old[i].box(),BL_SPACEDIM);
            S.copy(U_old[i],0,0,BL_SPACEDIM);

            if (do_mom_diff == 1)
            {
                for (int d = 0; d < BL_SPACEDIM; d++)
                {
                    tforces.mult((*Rh)[i],grids[i],grids[i],0,d,1);
                    S.mult((*rho_ptime)[i],grids[i],grids[i],0,d,1);
                }
            }

            godunov->Add_aofs_tf(S,U_new[i],0,BL_SPACEDIM,Aofs[i],
                                 0,tforces,0,grids[i],dt);

            if (do_mom_diff == 1)
            {
                for (int d = 0; d < BL_SPACEDIM; d++)
                    U_new[i].divide((*rho_ctime)[i],grids[i],grids[i],0,d,1);
            }
        }
    }
}

void
NavierStokes::errorEst (TagBoxArray& tags,
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
    Array<int>  itags;

    for (int j = 0; j < err_list.size(); j++)
    {
        MultiFab* mf = derive(err_list[j].name(), time, err_list[j].nGrow());

        for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
        {
            itags             = tags[mfi.index()].tags();
            FArrayBox&  fab   = (*mf)[mfi];
            int*        tptr  = itags.dataPtr();
            const int*  tlo   = tags[mfi.index()].box().loVect();
            const int*  thi   = tags[mfi.index()].box().hiVect();
            const int*  lo    = mfi.validbox().loVect();
            const int*  hi    = mfi.validbox().hiVect();
            const Real* xlo   = grid_loc[mfi.index()].lo();
            Real*       dat   = (*mf)[mfi].dataPtr();
            const int*  dlo   = (*mf)[mfi].box().loVect();
            const int*  dhi   = (*mf)[mfi].box().hiVect();
            const int   ncomp = (*mf)[mfi].nComp();

            err_list[j].errFunc()(tptr, ARLIM(tlo), ARLIM(thi), &tagval,
                                  &clearval, dat, ARLIM(dlo), ARLIM(dhi),
                                  lo,hi, &ncomp, domain_lo, domain_hi,
                                  dx, xlo, prob_lo, &time, &level);
            //
            // Don't forget to set the tags in the TagBox.
            //
            tags[mfi.index()].tags(itags);
        }

        delete mf;
    }
}

Real
NavierStokes::sumDerive (const std::string& name, Real time)
{
    Real      sum = 0.0;
    MultiFab* mf  = derive(name,time,0);

    BL_ASSERT(!(mf == 0));

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = mf->get(mfi);

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }

        sum += fab.sum(0);
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
NavierStokes::volWgtSum (const std::string& name,
                         Real           time)
{
    Real        sum     = 0.0;
    int         rz_flag = CoordSys::IsRZ() ? 1 : 0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive(name,time,0);
    Array<Real> tmp;

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }
        Real        s;
        const Real* dat = fab.dataPtr();
        const int*  dlo = fab.loVect();
        const int*  dhi = fab.hiVect();
        const int*  lo  = grids[mfi.index()].loVect();
        const int*  hi  = grids[mfi.index()].hiVect();
        Real*       rad = &radius[mfi.index()][0];

        tmp.resize(hi[1]-lo[1]+1);

#if (BL_SPACEDIM == 2)
        int irlo  = lo[0]-radius_grow;
        int irhi  = hi[0]+radius_grow;
        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
        FORT_SUMMASS(dat,ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),
                     dx,&s,rad,&irlo,&irhi,&rz_flag,tmp.dataPtr());
#endif

#if (BL_SPACEDIM == 3)
        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
        FORT_SUMMASS(dat,ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),
                     dx,&s,tmp.dataPtr());
#endif
        sum += s;
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

void
NavierStokes::sum_integrated_quantities ()
{
    const int finest_level = parent->finestLevel();

    Real time = state[State_Type].curTime();
//    Real mass = 0.0;
//    Real trac = 0.0;
    Real energy = 0.0;
    Real forcing = 0.0;

    for (int lev = 0; lev <= finest_level; lev++)
    {
        NavierStokes& ns_level = getLevel(lev);
//        mass += ns_level.volWgtSum("density",time);
//        trac += ns_level.volWgtSum("tracer",time);
        energy += ns_level.volWgtSum("energy",time);
	if (BL_SPACEDIM==3)
	    forcing += ns_level.volWgtSum("forcing",time);
    }

    if (ParallelDescriptor::IOProcessor())
    {
        const int old_prec = std::cout.precision(12);
        std::cout << '\n';
//        std::cout << "TIME= " << time << " MASS= " << mass << '\n';
//        std::cout << "TIME= " << time << " TRAC= " << trac << '\n';
        std::cout << "TIME= " << time << " KENG= " << energy << '\n';
	if (BL_SPACEDIM==3)
	    std::cout << "TIME= " << time << " FORC= " << forcing << '\n';
        std::cout.precision(old_prec);
    }
}

void
NavierStokes::TurbSum (Real time, Real *turb, int ksize, int turbVars)
{
#if (BL_SPACEDIM==3)
    const Real* dx = geom.CellSize();

    const int turbGrow(1);
    const int presGrow(0);
    MultiFab* turbMF = derive("TurbVars",time,turbGrow);
    MultiFab* presMF = derive("PresVars",time,presGrow);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter turbMfi(*turbMF), presMfi(*presMF);
	 turbMfi.isValid() && presMfi.isValid();
	 ++turbMfi, ++presMfi)
    {
	FArrayBox& turbFab = (*turbMF)[turbMfi];
	FArrayBox& presFab = (*presMF)[presMfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[turbMfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
              presFab.setVal(0,isects[ii].second,0,presMF->nComp());
              turbFab.setVal(0,isects[ii].second,0,turbMF->nComp());
            }
        }
    }

    turbMF->FillBoundary(0,turbMF->nComp());
    geom.FillPeriodicBoundary(*turbMF,0,turbMF->nComp());

    presMF->FillBoundary(0,presMF->nComp());
    geom.FillPeriodicBoundary(*presMF,0,presMF->nComp());
    
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
        const int*  lo  = grids[turbMfi.index()].loVect();
        const int*  hi  = grids[turbMfi.index()].hiVect();

        FORT_SUMTURB(turbData,presData,ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),
		     dx,turb,&ksize,&turbVars);
    }

    delete turbMF;
    delete presMF;
#else
    BoxLib::Error("TurbSum not implemented in 2D");
#endif
}

void
NavierStokes::sum_turbulent_quantities ()
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
    
        NavierStokes& ns_level = getLevel(lev);
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
        if (!BoxLib::UtilCreateDirectory(DirPath, 0755))
            BoxLib::CreateDirectoryFailed(DirPath);

        const int steps = parent->levelSteps(0);
        FILE *file;
        char filename[256];
        sprintf(filename,"TurbData/TurbData_%04d.dat",steps);
        file = fopen(filename,"w");
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

    //Geometry   geom       = parent->Geom(finestLevel);
    //RealBox    probDomain = geom.ProbDomain();
    //Box        domain     = geom.Domain();
    //CoordSys   cs         = (CoordSys&) geom;
    //std::cout << "parent->Geom(finestLevel) = " << parent->Geom(finestLevel) << std::endl;
    //std::cout << "geom = " << geom << std::endl;
    //std::cout << "probDomain = " << probDomain << std::endl;
    //std::cout << "Domain = " << domain << std::endl;
    //std::cout << "cs = " << cs << std::endl;
    //const Real *dx = ((CoordSys&)parent->Geom(finestLevel)).CellSize();
    //
    //std::cout << "Outside level loop:" << std::endl;
    //std::cout << "ksize = " << ksize << std::endl;
    //std::cout << "dx = " << dx[0] << ", "  << dx[1] << ", "  << dx[2] << std::endl;
    //
    //const Real *levDx = ((CoordSys&)parent->Geom(lev)).CellSize();
    //    std::cout << "Inside level loop: level = " << lev << std::endl;
    //    std::cout << "levKsize = " << levKsize << std::endl;
    //    std::cout << "levDx = " << levDx[0] << ", "  << levDx[1] << ", "  << levDx[2] << std::endl;

void
NavierStokes::setPlotVariables()
{
    AmrLevel::setPlotVariables();
}

std::string
NavierStokes::thePlotFileType () const
{
    //
    // Increment this whenever the writePlotFile() format changes.
    //
    static const std::string the_plot_file_type("NavierStokes-V1.1");

    return the_plot_file_type;
}

void
NavierStokes::writePlotFile (const std::string& dir,
                             std::ostream&  os,
                             VisMF::How     how)
{
    if ( ! Amr::Plot_Files_Output() ) return;

    int i, n;
    //
    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    //
    std::vector<std::pair<int,int> > plot_var_map;
    for (int typ = 0; typ < desc_lst.size(); typ++)
        for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
            if (parent->isStatePlotVar(desc_lst[typ].name(comp)) &&
                desc_lst[typ].getType() == IndexType::TheCellType())
                plot_var_map.push_back(std::pair<int,int>(typ,comp));

    int num_derive = 0;
    std::list<std::string> derive_names;
    const std::list<DeriveRec>& dlist = derive_lst.dlist();

    for (std::list<DeriveRec>::const_iterator it = dlist.begin();
         it != dlist.end();
         ++it)
    {
        if (parent->isDerivePlotVar(it->name()))
	{
            derive_names.push_back(it->name());
            num_derive += it->numDerive();
	}
    }

    int n_data_items = plot_var_map.size() + num_derive;
    Real cur_time = state[State_Type].curTime();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0)
            BoxLib::Error("Must specify at least one valid data item to plot");

        os << n_data_items << '\n';

	//
	// Names of variables -- first state, then derived
	//
	for (i =0; i < plot_var_map.size(); i++)
        {
	    int typ  = plot_var_map[i].first;
	    int comp = plot_var_map[i].second;
	    os << desc_lst[typ].name(comp) << '\n';
        }

	for (std::list<std::string>::const_iterator it = derive_names.begin();
             it != derive_names.end();
             ++it)
        {
	    const DeriveRec* rec = derive_lst.get(*it);
	    for (i = 0; i < rec->numDerive(); i++)
                os << rec->variableName(i) << '\n';
        }
        os << BL_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbLo(i) << ' ';
        os << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbHi(i) << ' ';
        os << '\n';
        for (i = 0; i < f_lev; i++)
            os << parent->refRatio(i)[0] << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->Geom(i).Domain() << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->levelSteps(i) << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++)
                os << parent->Geom(i).CellSize()[k] << ' ';
            os << '\n';
        }
        os << (int) CoordSys::Coord() << '\n';
        os << "0\n"; // Write bndry data.
    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";
    char buf[64];
    sprintf(buf, "Level_%d", level);
    std::string Level = buf;
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
        FullPath += '/';
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
            BoxLib::CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (i = 0; i < grids.size(); ++i)
        {
            for (n = 0; n < BL_SPACEDIM; n++)
                os << grid_loc[i].lo(n) << ' ' << grid_loc[i].hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (n_data_items > 0)
        {
            std::string PathNameInHeader = Level;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }
    }
    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    // NOTE: we are assuming that each state variable has one component,
    // but a derived variable is allowed to have multiple components.
    int       cnt   = 0;
    int       ncomp = 1;
    const int nGrow = 0;
    MultiFab  plotMF(grids,n_data_items,nGrow);
    MultiFab* this_dat = 0;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (i = 0; i < plot_var_map.size(); i++)
    {
	int typ  = plot_var_map[i].first;
	int comp = plot_var_map[i].second;
	this_dat = &state[typ].newData();
	MultiFab::Copy(plotMF,*this_dat,comp,cnt,ncomp,nGrow);
	cnt+= ncomp;
    }
    //
    // Cull data from derived variables.
    // 
    Real plot_time;

    if (derive_names.size() > 0)
    {
	for (std::list<std::string>::const_iterator it = derive_names.begin();
             it != derive_names.end();
             ++it) 
	{
            if (*it == "avg_pressure" ||
                *it == "gradpx"       ||
                *it == "gradpy"       ||
                *it == "gradpz") 
            {
                if (state[Press_Type].descriptor()->timeType() == 
                    StateDescriptor::Interval) 
                {
                    plot_time = cur_time;
                }
                else
                {
                    int f_lev = parent->finestLevel();
                    plot_time = getLevel(f_lev).state[Press_Type].curTime();
                }
            }
            else
            {
                plot_time = cur_time;
            } 
	    const DeriveRec* rec = derive_lst.get(*it);
	    ncomp = rec->numDerive();
	    MultiFab* derive_dat = derive(*it,plot_time,nGrow);
	    MultiFab::Copy(plotMF,*derive_dat,0,cnt,ncomp,nGrow);
	    delete derive_dat;
	    cnt += ncomp;
	}
    }
    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(plotMF,TheFullPath,how,true);
}

Real
NavierStokes::estTimeStep ()
{
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

    const int   n_grow        = 0;
    Real        estdt         = 1.0e+20;
    const Real  cur_time      = state[State_Type].curTime();
    const Real  cur_pres_time = state[Press_Type].curTime();
    MultiFab&   P_new         = get_new_data(Press_Type);
    MultiFab&   U_new         = get_new_data(State_Type);
    const Real* dx            = geom.CellSize();

    Real gr_max[BL_SPACEDIM], u_max[BL_SPACEDIM] = {0};

    FArrayBox p_fab, tforces;
    MultiFab Gp(grids,BL_SPACEDIM,1);
    getGradP(Gp, cur_pres_time);

    for (MFIter Rho_mfi(*rho_ctime); Rho_mfi.isValid(); ++Rho_mfi)
    {
        const int i = Rho_mfi.index();
        //
        // Get the pressure.
        //
        p_fab.resize(BoxLib::surroundingNodes(grids[i]),1);
        p_fab.copy(P_new[i],p_fab.box());
        //
        // Get the velocity forcing.  For some reason no viscous forcing.
        //
        getForce(tforces,i,n_grow,Xvel,BL_SPACEDIM,
#ifdef GENGETFORCE
		 cur_time,
#endif		 
		 (*rho_ctime)[i]);
        tforces.minus(Gp[i],0,0,BL_SPACEDIM);
        //
        // Estimate the maximum allowable timestep from the Godunov box.
        //
        Real dt = godunov->estdt(U_new[i],tforces,(*rho_ctime)[i],grids[i],
                                 geom.CellSize(),cfl,gr_max);

        for (int k = 0; k < BL_SPACEDIM; k++)
        {
	    u_max[k] = std::max(u_max[k],gr_max[k]);
	}
	estdt = std::min(estdt,dt);
    }

    ParallelDescriptor::ReduceRealMin(estdt);

    if (verbose)
    {
        for (int k = 0; k < BL_SPACEDIM; k++)
        {
            ParallelDescriptor::ReduceRealMax(u_max[k]);
        }
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "estTimeStep :: \n" << "LEV = " << level << " UMAX = ";
            for (int k = 0; k < BL_SPACEDIM; k++)
                std::cout << u_max[k] << "  ";
            std::cout << '\n';
        }
    }

    return estdt;
}

Real
NavierStokes::initialTimeStep ()
{
    return init_shrink*estTimeStep();
}

void
NavierStokes::computeNewDt (int                   finest_level,
                            int                   sub_cycle,
                            Array<int>&           n_cycle,
                            const Array<IntVect>& ref_ratio,
                            Array<Real>&          dt_min,
                            Array<Real>&          dt_level,
                            Real                  stop_time)
{
    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0) return;
    //
    // For Navier Stokes compute the new dt based on current velocity field.
    //
    const int max_level = parent->maxLevel();

    n_cycle[0] = 1;
    for (int i = 1; i <= max_level; i++)
    {
        n_cycle[i] = sub_cycle ? parent->MaxRefRatio(i-1) : 1;
    }

    Real dt_0     = 1.0e+100;
    int  n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        dt_min[i] = std::min(dt_min[i],getLevel(i).estTimeStep());
        if (fixed_dt <= 0.0) 
          dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
        n_factor *= n_cycle[i];
        dt_0      = std::min(dt_0,n_factor*dt_min[i]);
    }

    const Real eps      = 0.0001*dt_0;
    const Real cur_time = state[State_Type].curTime();
    if (stop_time >= 0.0)
    {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }
    //
    // Adjust the time step to be able to output checkpoints at specific times.
    //
    const Real check_per = parent->checkPer();
    if (check_per > 0.0)
    {
        int a = int((cur_time + eps ) / check_per);
        int b = int((cur_time + dt_0) / check_per);
        if (a != b)
            dt_0 = b * check_per - cur_time;
    }
    //
    // Adjust the time step to be able to output plot files at specific times.
    //
    const Real plot_per = parent->plotPer();
    if (plot_per > 0.0)
    {
        int a = int((cur_time + eps ) / plot_per);
        int b = int((cur_time + dt_0) / plot_per);
        if (a != b)
            dt_0 = b * plot_per - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= max_level; i++)
    {
        n_factor   *= n_cycle[i];
        dt_level[i] = dt_0/( (Real)n_factor );
    }
}

void
NavierStokes::computeInitialDt (int                   finest_level,
                                int                   sub_cycle,
                                Array<int>&           n_cycle,
                                const Array<IntVect>& ref_ratio,
                                Array<Real>&          dt_level, 
                                Real                  stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;
    //
    // For Navier Stokes compute the new dt based on current velocity field.
    //
    const int max_level = parent->maxLevel();

    n_cycle[0] = 1;
    for (int i = 1; i <= max_level; i++)
    {
        n_cycle[i] = sub_cycle ? parent->MaxRefRatio(i-1) : 1;
    }

    Real dt_0    = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor   *= n_cycle[i];
        dt_0        = std::min(dt_0,n_factor*dt_level[i]);
    }

    if (stop_time >= 0.0)
    {
        const Real eps      = 0.0001*dt_0;
        const Real cur_time = state[State_Type].curTime();
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= max_level; i++)
    {
        n_factor   *= n_cycle[i];
        dt_level[i] = dt_0/( (Real)n_factor );
    }
}

//
// This function estimates the initial timesteping used by the model.
//

void
NavierStokes::post_init_estDT (Real&        dt_init,
                               Array<int>&  nc_save,
                               Array<Real>& dt_save,
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
 
    Array<Real> dt_level(finest_level+1,dt_init);
    Array<int>  n_cycle(finest_level+1,1);

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
// Fills in amrLevel okToContinue.
//

int
NavierStokes::okToContinue ()
{
    return (level > 0) ? true : (parent->dtLevel(0) > dt_cutoff);
}

//
// THE MAIN HOOKS INTO AMR AND AMRLEVEL
//

//
// Integration cycle on fine level grids is complete .
// post_timestep() is responsible for syncing levels together.
//
// The registers used for level syncing are initialized in the
// coarse level advance and incremented in the fine level advance.
// These quantities are described in comments above advance_setup.
//

void
NavierStokes::post_timestep (int crse_iteration)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::post_timestep()");

    const int finest_level = parent->finestLevel();

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
    //
    // Derive turbulent statistics
    //
    if (level==0 && turb_interval>0 && (parent->levelSteps(0)%turb_interval == 0) && BL_SPACEDIM==3)
    {
        sum_turbulent_quantities();
    }

    if (level > 0) incrPAvg();

    old_intersect_new          = grids;
    is_first_step_after_regrid = false;
}

//
// Build any additional data structures after restart.
//

void NavierStokes::post_restart()
{
    make_rho_prev_time();
    make_rho_curr_time();
}

//
// Build any additional data structures after regrid.
//

void
NavierStokes::post_regrid (int lbase,
                                int new_finest)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::post_regrid()");
    if (projector && level == lbase)
        projector->setFinestLevel(new_finest);
}

//
// Ensure state, and pressure are consistent.
//

void
NavierStokes::post_init (Real stop_time)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::post_init()");

    if (level > 0)
        //
        // Nothing to sync up at level > 0.
        //
        return;

    MultiFab&   P_new        = get_new_data(Press_Type);
    MultiFab&   P_old        = get_old_data(Press_Type);
    const int   finest_level = parent->finestLevel();
    Real        dt_init      = 0.0;
    Array<Real> dt_save(finest_level+1);
    Array<int>  nc_save(finest_level+1);
    //
    // Ensure state is consistent, i.e. velocity field is non-divergent,
    // Coarse levels are fine level averages, pressure is zero.
    //
    post_init_state();
    //
    // Estimate the initial timestepping.
    //
    post_init_estDT(dt_init, nc_save, dt_save, stop_time);
    //
    // Initialize the pressure by iterating the initial timestep.
    //
    post_init_press(dt_init, nc_save, dt_save);
    //
    // Compute the initial estimate of conservation.
    //
    if (sum_interval > 0)
        sum_integrated_quantities();
    //
    // Derive turbulent statistics
    //
    if (turb_interval > 0)
        sum_turbulent_quantities();
}

//
// MULTILEVEL SYNC FUNCTIONS
//

void
NavierStokes::initRhoAvg (Real alpha)
{
    MultiFab& S_new = get_new_data(State_Type);

    rho_avg->setVal(0);

    for (MFIter rho_avgmfi(*rho_avg); rho_avgmfi.isValid(); ++rho_avgmfi)
    {
        const int i = rho_avgmfi.index();
        (*rho_avg)[i].copy(S_new[i],S_new.box(i),Density,S_new.box(i),0,1);
        (*rho_avg)[i].mult(alpha);
    }
}

void
NavierStokes::incrRhoAvg(const MultiFab& rho_incr,
                         int             sComp,
                         Real            alpha)
{
    for (MFIter mfi(rho_incr); mfi.isValid(); ++mfi)
    {
        const int* lo      = mfi.validbox().loVect();
        const int* hi      = mfi.validbox().hiVect();
        const int* rlo     = rho_incr[mfi].loVect();
        const int* rhi     = rho_incr[mfi].hiVect();
        const Real* rhodat = rho_incr[mfi].dataPtr(sComp);
        const int* alo     = (*rho_avg)[mfi].loVect();
        const int* ahi     = (*rho_avg)[mfi].hiVect();
        Real* avgdat       = (*rho_avg)[mfi].dataPtr();

        FORT_INCRMULT(avgdat,ARLIM(alo),ARLIM(ahi),
                      rhodat,ARLIM(rlo),ARLIM(rhi),lo,hi,&alpha);
    }
}

void
NavierStokes::incrRhoAvg (Real alpha)
{
    const MultiFab& S_new = get_new_data(State_Type);
    incrRhoAvg(S_new,Density,alpha);
}

void
NavierStokes::incrPAvg ()
{
    //
    // Increment p_avg with 1/ncycle times current pressure
    //
    MultiFab& P_new = get_new_data(Press_Type);

    Real alpha = 1.0/Real(parent->nCycle(level));

    for (MFIter P_newmfi(P_new); P_newmfi.isValid(); ++P_newmfi)
    {
        const int*  lo     = P_newmfi.validbox().loVect();
        const int*  hi     = P_newmfi.validbox().hiVect();
        const int*  p_lo   = P_new[P_newmfi].loVect();
        const int*  p_hi   = P_new[P_newmfi].hiVect();
        const Real* pdat   = P_new[P_newmfi].dataPtr();
        const int*  alo    = (*p_avg)[P_newmfi].loVect();
        const int*  ahi    = (*p_avg)[P_newmfi].hiVect();
        Real*       avgdat = (*p_avg)[P_newmfi].dataPtr();

        FORT_INCRMULT(avgdat,ARLIM(alo),ARLIM(ahi),
                      pdat,ARLIM(p_lo),ARLIM(p_hi),lo,hi,&alpha);
    }
}

//
// This function ensures that the state is initially consistent
// with respect to the divergence condition and fields are initially consistent
//

void
NavierStokes::post_init_state ()
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::post_init_state()");

    const int finest_level = parent->finestLevel();
    const Real divu_time   = have_divu ? state[Divu_Type].curTime()
                                       : state[Press_Type].curTime();
    if (do_init_vort_proj)
    {
        //
	// NOTE: this assumes have_divu == 0.
	// Only used if vorticity is used to initialize the velocity field.
        //
        BL_ASSERT(!(projector == 0));
        
	projector->initialVorticityProject(0);
    }

    if (do_init_proj && projector)
        //
        // Do sync project to define divergence free velocity field.
        //
        projector->initialVelocityProject(0,divu_time,have_divu);

    NavierStokes::initial_step = true;
    //
    // Average velocity and scalar data down from finer levels
    // so that conserved data is consistant between levels.
    //
    for (int k = finest_level-1; k>= 0; k--)
    {
        getLevel(k).avgDown();
    }
    make_rho_curr_time();

    if (do_init_proj && projector && (std::abs(gravity)) > 0.)
        //
        // Do projection to establish initially hydrostatic pressure field.
        //
        projector->initialPressureProject(0);
}

//
// Initialize the pressure by iterating the initial timestep
//

void
NavierStokes::post_init_press (Real&        dt_init,
                               Array<int>&  nc_save,
                               Array<Real>& dt_save)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::post_init_press()");

    const Real strt_time       = state[State_Type].curTime();
    const int  finest_level    = parent->finestLevel();
    MultiFab&  P_new           = get_new_data(Press_Type);
    MultiFab&  P_old           = get_old_data(Press_Type);
    NavierStokes::initial_iter = true;
    //
    // Iterate over the advance function.
    //
    for (int iter = 0; iter < init_iter; iter++)
    {
        for (int k = 0; k <= finest_level; k++ )
        {
            getLevel(k).advance(strt_time,dt_init,1,1);
        }
        //
        // This constructs a guess at P, also sets p_old == p_new.
        //
        MultiFab** sig = new MultiFab*[finest_level+1];

        for (int k = 0; k <= finest_level; k++)
        {
            sig[k] = getLevel(k).get_rho_half_time();
        }
        if (projector)
            projector->initialSyncProject(0,sig,parent->dtLevel(0),
                                          strt_time,dt_init,have_divu);
        delete [] sig;

        for (int k = finest_level-1; k >= 0; k--)
        {
            getLevel(k).avgDown();
        }
        for (int k = 0; k <= finest_level; k++)
        {
            //
            // Reset state variables to initial time, but 
            // do not reset pressure variable, only pressure time.
            //
            getLevel(k).resetState(strt_time, dt_init, dt_init);
        }

        make_rho_curr_time();

        NavierStokes::initial_iter = false;
    }

    if (init_iter <= 0)
        NavierStokes::initial_iter = false; // Just being compulsive -- rbp.

    NavierStokes::initial_step = false;
    //
    // Re-instate timestep.
    //
    for (int k = 0; k <= finest_level; k++)
    {
        getLevel(k).setTimeLevel(strt_time,dt_save[k],dt_save[k]);
        if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
        {
            getLevel(k).state[Press_Type].setNewTimeLevel(.5*dt_init);
            getLevel(k).get_old_data(Dpdt_Type).setVal(0);
        }
    }

    parent->setDtLevel(dt_save);
    parent->setNCycle(nc_save);
}

//
// Helper function for NavierStokes::SyncInterp().
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
            for (int crse = 0; crse < cgrids.size(); crse++)
            {
                const int* c_lo = cgrids[crse].loVect();
                const int* c_hi = cgrids[crse].hiVect();

                if (clo[dir] < cdomlo[dir] && c_lo[dir] == cdomlo[dir])
                    bc_new[bc_index] = bc_orig_qty[crse][bc_index];
                if (chi[dir] > cdomhi[dir] && c_hi[dir] == cdomhi[dir])
                    bc_new[bc_index+BL_SPACEDIM] = bc_orig_qty[crse][bc_index+BL_SPACEDIM]; 
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
// The components of bc_orig_qty corespond to the quantities of CrseSync.
//

void
NavierStokes::SyncInterp (MultiFab&      CrseSync,
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
    BL_ASSERT(which_interp >= 0 && which_interp <= 5);

    Interpolater* interpolater = 0;

    switch (which_interp)
    {
    case PC_T:           interpolater = &pc_interp;           break;
    case CellCons_T:     interpolater = &cell_cons_interp;    break;
    case CellConsLin_T:  interpolater = &lincc_interp;        break;
    case CellConsProt_T: interpolater = &protected_interp;    break;
    default:
        BoxLib::Abort("NavierStokes::SyncInterp(): how did this happen");
    }

    NavierStokes&   fine_level = getLevel(f_lev);
    const BoxArray& fgrids     = fine_level.boxArray();
    const Geometry& fgeom      = parent->Geom(f_lev);
    const BoxArray& cgrids     = getLevel(c_lev).boxArray();
    const Geometry& cgeom      = parent->Geom(c_lev);
    const Real*     dx_crse    = cgeom.CellSize();
    Box             cdomain    = BoxLib::coarsen(fgeom.Domain(),ratio);
    const int*      cdomlo     = cdomain.loVect();
    const int*      cdomhi     = cdomain.hiVect();
    int*            bc_new     = new int[2*BL_SPACEDIM*(src_comp+num_comp)];

    BoxArray cdataBA(fgrids.size());

    for (int i = 0; i < fgrids.size(); i++)
        cdataBA.set(i,interpolater->CoarseBox(fgrids[i],ratio));
    //
    // Note: The boxes in cdataBA may NOT be disjoint !!!
    //
    MultiFab cdataMF(cdataBA,num_comp,0);

    cdataMF.setVal(0);

    cdataMF.copy(CrseSync, src_comp, 0, num_comp);
    //
    // Set physical boundary conditions in cdataMF.
    //
    for (MFIter mfi(cdataMF); mfi.isValid(); ++mfi)
    {
        int         i     = mfi.index();
        const Box&  grd   = fgrids[i];
        FArrayBox&  cdata = cdataMF[mfi];
        const Box&  cgrd  = cdata.box();
        const int*  clo   = cdata.loVect();
        const int*  chi   = cdata.hiVect();
        const Real* xlo   = fine_level.grid_loc[i].lo();

        for (int n = 0; n < num_comp; n++)
        {
            set_bc_new(bc_new,n,src_comp,clo,chi,cdomlo,cdomhi,cgrids,bc_orig_qty);

            FORT_FILCC(cdata.dataPtr(n), ARLIM(clo), ARLIM(chi),
                       cdomlo, cdomhi, dx_crse, xlo,
                       &(bc_new[2*BL_SPACEDIM*(n+src_comp)]));
        }
    }
    cgeom.FillPeriodicBoundary(cdataMF, 0, num_comp);
    //
    // Interpolate from cdataMF to fdata and update FineSync.
    // Note that FineSync and cdataMF will have the same distribution
    // since the length of their BoxArrays are equal.
    //
    FArrayBox    fdata;
    Array<BCRec> bc_interp(num_comp);

    MultiFab* fine_stateMF;
    if (interpolater == &protected_interp)
        fine_stateMF = &(getLevel(f_lev).get_new_data(State_Type));

    for (MFIter mfi(cdataMF); mfi.isValid(); ++mfi)
    {
        int        i     = mfi.index();
        FArrayBox& cdata = cdataMF[mfi];
        const int* clo   = cdata.loVect();
        const int* chi   = cdata.hiVect();

        fdata.resize(fgrids[i], num_comp);
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

        ScaleCrseSyncInterp(cdata, c_lev, num_comp);

        interpolater->interp(cdata,0,fdata,0,num_comp,fgrids[i],ratio,
                             cgeom,fgeom,bc_interp);

        reScaleFineSyncInterp(fdata, f_lev, num_comp);

        if (increment)
        {
            fdata.mult(dt_clev);

            if (interpolater == &protected_interp) {

              cdata.mult(dt_clev);
              FArrayBox& fine_state = (*fine_stateMF)[i];
              interpolater->protect(cdata,0,fdata,0,fine_state,state_comp,
                                    num_comp,fgrids[i],ratio,
                                    cgeom,fgeom,bc_interp);
              Real dt_clev_inv = 1./dt_clev;
              cdata.mult(dt_clev_inv);

            }
            
            FineSync[i].plus(fdata,0,dest_comp,num_comp);
        }
        else
        {
            FineSync[i].copy(fdata,0,dest_comp,num_comp);
        }
    }

    delete [] bc_new;
}

//
// Interpolate sync pressure correction to a finer level.
//

void
NavierStokes::SyncProjInterp (MultiFab& phi,
                              int       c_lev,
                              MultiFab& P_new,
                              MultiFab& P_old,
                              int       f_lev,
                              IntVect&  ratio,
                              bool      first_crse_step_after_initial_iters,
                              Real      cur_crse_pres_time,
                              Real      prev_crse_pres_time)
{
    const Geometry& fgeom   = parent->Geom(f_lev);
    const BoxArray& P_grids = P_new.boxArray();
    const Geometry& cgeom   = parent->Geom(c_lev);

    BoxArray crse_ba(P_grids.size());

    for (int i = 0; i < P_grids.size(); i++)
    {
        crse_ba.set(i,node_bilinear_interp.CoarseBox(P_grids[i],ratio));
    }

    Array<BCRec> bc(BL_SPACEDIM);
    MultiFab     crse_phi(crse_ba,1,0);

    crse_phi.setVal(1.e200);
    crse_phi.copy(phi,0,0,1);

    FArrayBox     fine_phi;
    NavierStokes& fine_lev            = getLevel(f_lev);
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

        for (MFIter mfi(crse_phi); mfi.isValid(); ++mfi)
        {
            fine_phi.resize(P_grids[mfi.index()],1);
            fine_phi.setVal(1.e200);
            node_bilinear_interp.interp(crse_phi[mfi],0,fine_phi,0,1,
                                        fine_phi.box(),ratio,cgeom,fgeom,bc);
            fine_phi.mult(cur_mult_factor);
            P_new[mfi.index()].plus(fine_phi);
            fine_phi.mult(prev_mult_factor);
            P_old[mfi.index()].plus(fine_phi);
        }
    }
    else 
    {
        for (MFIter mfi(crse_phi); mfi.isValid(); ++mfi)
        {
            fine_phi.resize(P_grids[mfi.index()],1);
            fine_phi.setVal(1.e200);
            node_bilinear_interp.interp(crse_phi[mfi],0,fine_phi,0,1,
                                        fine_phi.box(),ratio,cgeom,fgeom,bc);
            P_new[mfi.index()].plus(fine_phi);
            P_old[mfi.index()].plus(fine_phi);
        }
    }
}

//
// Averages a multifab of fine data down onto a multifab of coarse data.
//
// This should be an Amrlevel or Multifab function
//

void
NavierStokes::avgDown (const BoxArray& cgrids,
                       const BoxArray& fgrids,
                       MultiFab&       S_crse,
                       MultiFab&       S_fine,
                       MultiFab&       cvolume,
                       MultiFab&       fvolume,
                       int             c_level,
                       int             f_level,
                       int             scomp,
                       int             ncomp,
                       const IntVect&  fratio)
{
    BL_ASSERT(cgrids == S_crse.boxArray());
    BL_ASSERT(fgrids == S_fine.boxArray());
    BL_ASSERT(cvolume.boxArray() == cgrids);
    BL_ASSERT(fvolume.boxArray() == fgrids);
    BL_ASSERT(S_crse.nComp() == S_fine.nComp());
    BL_ASSERT(fvolume.nComp() == 1 && cvolume.nComp() == 1);

    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::avgDown()");
    //
    // Coarsen() the fine stuff on processors owning the fine data.
    //
    BoxArray crse_S_fine_BA(fgrids.size());

    for (int i = 0; i < fgrids.size(); ++i)
    {
        crse_S_fine_BA.set(i,BoxLib::coarsen(fgrids[i],fratio));
    }

    MultiFab crse_S_fine(crse_S_fine_BA,ncomp,0);
    MultiFab crse_fvolume(crse_S_fine_BA,1,0);

    crse_fvolume.copy(cvolume);

    for (MFIter mfi(S_fine); mfi.isValid(); ++mfi)
    {
        const int i = mfi.index();

        avgDown(S_fine[i],crse_S_fine[i],fvolume[i],crse_fvolume[i],
                f_level,c_level,crse_S_fine_BA[i],scomp,ncomp,fratio);
    }

    S_crse.copy(crse_S_fine,0,scomp,ncomp);
}

//
// Average fine down to coarse in the ovlp intersection.
//

void
NavierStokes::avgDown (const FArrayBox& fine_fab,
                       const FArrayBox& crse_fab, 
                       const FArrayBox& fine_vol,
                       const FArrayBox& crse_vol,
                       int              f_level,
                       int              c_level,
                       const Box&       ovlp,
                       int              scomp,
                       int              ncomp,
                       const IntVect&   fratio)
{
    avgDown_doit(fine_fab,crse_fab,fine_vol,crse_vol,
                 f_level,c_level,ovlp,scomp,ncomp,fratio);
}

//
// Actually average the data down (this is static)
//

void
NavierStokes::avgDown_doit (const FArrayBox& fine_fab,
                            const FArrayBox& crse_fab, 
                            const FArrayBox& fine_vol,
                            const FArrayBox& crse_vol,
                            int              f_level,
                            int              c_level,
                            const Box&       ovlp,
                            int              scomp,
                            int              ncomp,
                            const IntVect&   fratio)
{
//
//  NOTE: We copy from component scomp of the fine fab into component 0 of the crse fab
//        because the crse fab is a temporary which was made starting at comp 0, it is
//        not the actual state data.
//
    const int*  ovlo   = ovlp.loVect();
    const int*  ovhi   = ovlp.hiVect();
    const int*  flo    = fine_fab.loVect();
    const int*  fhi    = fine_fab.hiVect();
    const Real* f_dat  = fine_fab.dataPtr(scomp);
    const int*  fvlo   = fine_vol.loVect();
    const int*  fvhi   = fine_vol.hiVect();
    const Real* fv_dat = fine_vol.dataPtr();
    const int*  clo    = crse_fab.loVect();
    const int*  chi    = crse_fab.hiVect();
    const Real* c_dat  = crse_fab.dataPtr();
    const int*  cvlo   = crse_vol.loVect();
    const int*  cvhi   = crse_vol.hiVect();
    const Real* cv_dat = crse_vol.dataPtr();

    FORT_AVGDOWN(c_dat,ARLIM(clo),ARLIM(chi),&ncomp,
                 f_dat,ARLIM(flo),ARLIM(fhi),
                 cv_dat,ARLIM(cvlo),ARLIM(cvhi),
                 fv_dat,ARLIM(fvlo),ARLIM(fvhi),
                 ovlo,ovhi,fratio.getVect());
}

void
NavierStokes::level_sync (int crse_iteration)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::level_sync()");

    const Real*     dx            = geom.CellSize();
    IntVect         ratio         = parent->refRatio(level);
    const int       finest_level  = parent->finestLevel();
    int             crse_dt_ratio = parent->nCycle(level);
    Real            dt            = parent->dtLevel(level);
    const Real      half_time     = state[State_Type].prevTime() + 0.5*dt;
    MultiFab&       pres          = get_new_data(Press_Type);
    MultiFab&       vel           = get_new_data(State_Type);
    SyncRegister&   rhs_sync_reg  = getLevel(level+1).getSyncReg();
    SyncRegister*   crsr_sync_ptr = 0;
    NavierStokes&   fine_level    = getLevel(level+1);
    MultiFab&       pres_fine     = fine_level.get_new_data(Press_Type);
    MultiFab&       vel_fine      = fine_level.get_new_data(State_Type);
    const BoxArray& finegrids     = vel_fine.boxArray();
    
    if (level > 0)
        crsr_sync_ptr = &(getLevel(level).getSyncReg());
    //
    // Get boundary conditions.
    //
    Array<int*>         sync_bc(grids.size());
    Array< Array<int> > sync_bc_array(grids.size());

    for (int i = 0; i < grids.size(); i++)
    {
        sync_bc_array[i] = getBCArray(State_Type,i,Xvel,BL_SPACEDIM);
        sync_bc[i] = sync_bc_array[i].dataPtr();
    }

    MultiFab cc_rhs_crse(grids,1,1);
    MultiFab cc_rhs_fine(finegrids,1,1);
    cc_rhs_crse.setVal(0);
    cc_rhs_fine.setVal(0);

    MultiFab new_divu_crse(grids,1,0);
    MultiFab new_divu_fine(finegrids,1,0);
    //
    // At this point the Divu state data is what was used in the original
    // level project and has not been updated by avgDown or mac_sync.
    // We want to fill cc_rhs_crse and cc_rhs_fine with the difference
    // between the divu we now define using calc_divu and the divu which
    // is in the state data.
    // We are also copying the new computed value of divu into the Divu state.
    //
    if (do_sync_proj && have_divu && do_divu_sync == 1) 
    {
        const Real cur_time = state[Divu_Type].curTime();
        const Real dt_inv = 1.0 / dt;

        MultiFab& cur_divu_crse = get_new_data(Divu_Type);
        calc_divu(cur_time,dt,cc_rhs_crse);
        MultiFab::Copy(new_divu_crse,cc_rhs_crse,0,0,1,0);
        cc_rhs_crse.minus(cur_divu_crse,0,1,0);
        MultiFab::Copy(cur_divu_crse,new_divu_crse,0,0,1,0);
        cc_rhs_crse.mult(dt_inv,0,1,0);

        NavierStokes& fine_lev = getLevel(level+1);
        MultiFab& cur_divu_fine = fine_lev.get_new_data(Divu_Type);
        fine_lev.calc_divu(cur_time,dt,cc_rhs_fine);
        MultiFab::Copy(new_divu_fine,cc_rhs_fine,0,0,1,0);
        cc_rhs_fine.minus(cur_divu_fine,0,1,0);
        MultiFab::Copy(cur_divu_fine,new_divu_fine,0,0,1,0);
        cc_rhs_fine.mult(dt_inv,0,1,0);
        //
        // With new divu's, get new Dsdt, then average down.
        //
        calc_dsdt(cur_time, dt, get_new_data(Dsdt_Type));
        fine_lev.calc_dsdt(cur_time, dt/crse_dt_ratio,
                           fine_lev.get_new_data(Dsdt_Type));
        for (int k = level; k>= 0; k--)
        {
            NavierStokes&   flev     = getLevel(k+1);
            const BoxArray& fgrids   = flev.grids;
            MultiFab&       fvolume  = flev.volume;
          
            NavierStokes&   clev     = getLevel(k);
            const BoxArray& cgrids   = clev.grids;
            MultiFab&       cvolume  = clev.volume;
          
            IntVect&  fratio = clev.fine_ratio;
          
            NavierStokes::avgDown(cgrids, fgrids,
                                  clev.get_new_data(Divu_Type),
                                  flev.get_new_data(Divu_Type),
                                  cvolume, fvolume,
                                  k, k+1, 0, 1, fratio);

            NavierStokes::avgDown(cgrids, fgrids,
                                  clev.get_new_data(Dsdt_Type),
                                  flev.get_new_data(Dsdt_Type),
                                  cvolume, fvolume,
                                  k, k+1, 0, 1, fratio);
        }
    }
    //
    // Multilevel or single-level sync projection.
    //
    MultiFab* Rh = get_rho_half_time();

    if (do_MLsync_proj)
    {
        
        MultiFab&       vel_fine    = fine_level.get_new_data(State_Type);
        MultiFab&       rho_fine    = *fine_level.rho_avg;
        const Geometry& fine_geom   = parent->Geom(level+1);
        const Geometry& crse_geom   = parent->Geom(level);
        const BoxArray& P_finegrids = pres_fine.boxArray();

        MultiFab phi(P_finegrids,1,1);
        MultiFab V_corr(finegrids,BL_SPACEDIM,1);

        V_corr.setVal(0);
        //
        // If periodic, enforce periodicity on Vsync.
        //
        projector->EnforcePeriodicity(*Vsync, BL_SPACEDIM, grids, crse_geom);
        //
        // Interpolate Vsync to fine grid correction in Vcorr.
        //
        SyncInterp(*Vsync, level, V_corr, level+1, ratio,
                   0, 0, BL_SPACEDIM, 0 , dt, sync_bc.dataPtr());
        //
        // The multilevel projection.  This computes the projection and
        // adds in its contribution to levels (level) and (level+1).
        //
        Real  cur_crse_pres_time = state[Press_Type].curTime();
        Real prev_crse_pres_time = state[Press_Type].prevTime();

        NavierStokes& fine_lev   = getLevel(level+1);
        Real  cur_fine_pres_time = fine_lev.state[Press_Type].curTime();
        Real prev_fine_pres_time = fine_lev.state[Press_Type].prevTime();

        bool first_crse_step_after_initial_iters =
         (prev_crse_pres_time > state[State_Type].prevTime());

        bool pressure_time_is_interval = 
         (state[Press_Type].descriptor()->timeType() == StateDescriptor::Interval);
        projector->MLsyncProject(level,pres,vel,cc_rhs_crse,
                                 pres_fine,vel_fine,cc_rhs_fine,
                                 *Rh,rho_fine,Vsync,V_corr,
                                 phi,&rhs_sync_reg,crsr_sync_ptr,
                                 dx,dt,ratio,crse_iteration,crse_dt_ratio, 
                                 fine_geom,geom,pressure_time_is_interval,
                                 first_crse_step_after_initial_iters,
                                 cur_crse_pres_time,prev_crse_pres_time,
                                 cur_fine_pres_time,prev_fine_pres_time);
        //
        // Correct pressure and velocities after the projection.
        //
        ratio = IntVect::TheUnitVector();
        Array<int*>         fine_sync_bc(finegrids.size());
        Array< Array<int> > fine_sync_bc_array(finegrids.size());

        for (int i = 0; i < finegrids.size(); i++)
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
            NavierStokes& fine_lev = getLevel(lev);
            MultiFab&     P_new    = fine_lev.get_new_data(Press_Type);
            MultiFab&     P_old    = fine_lev.get_old_data(Press_Type);
            MultiFab&     U_new    = fine_lev.get_new_data(State_Type);

            SyncInterp(V_corr, level+1, U_new, lev, ratio,
                       0, 0, BL_SPACEDIM, 1 , dt, fine_sync_bc.dataPtr());
            SyncProjInterp(phi, level+1, P_new, P_old, lev, ratio,
                           first_crse_step_after_initial_iters,
                           cur_crse_pres_time, prev_crse_pres_time);
        }

        if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
            calcDpdt();
    }
    else if (do_sync_proj) 
    {
        MultiFab phi(pres.boxArray(),1,1);
        BoxArray sync_boxes = pres_fine.boxArray();
        sync_boxes.coarsen(ratio);
        //
        // The single level projection.  This computes the projection and
        // adds in its contribution to level (level).
        //
        projector->syncProject(level,pres,vel,Rh,Vsync,phi,
                               &rhs_sync_reg,crsr_sync_ptr,sync_boxes,
                               geom,dx,dt,crse_iteration,crse_dt_ratio);
        //
        // Correct pressure and velocities after the projection.
        //
        ratio = IntVect::TheUnitVector(); 

        const Real cur_crse_pres_time  = state[Press_Type].curTime();
        const Real prev_crse_pres_time = state[Press_Type].prevTime();

        bool first_crse_step_after_initial_iters =
         (prev_crse_pres_time > state[State_Type].prevTime());

        for (int lev = level+1; lev <= finest_level; lev++)
        {
            ratio                 *= parent->refRatio(lev-1);
            NavierStokes& fine_lev = getLevel(lev);
            MultiFab&     P_new    = fine_lev.get_new_data(Press_Type);
            MultiFab&     P_old    = fine_lev.get_old_data(Press_Type);
            MultiFab&     U_new    = fine_lev.get_new_data(State_Type);
            
            SyncInterp(*Vsync, level, U_new, lev, ratio,
                       0, 0, BL_SPACEDIM, 1 , dt, sync_bc.dataPtr());
            SyncProjInterp(phi, level, P_new, P_old, lev, ratio,
                           first_crse_step_after_initial_iters,
                           cur_crse_pres_time, prev_crse_pres_time);
        }

        if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
            calcDpdt();
    }
}

//
// The Mac Sync correction function
//

void
NavierStokes::mac_sync ()
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::mac_sync()");

    const int  numscal        = NUM_STATE - BL_SPACEDIM;
    const Real prev_time      = state[State_Type].prevTime();
    const Real prev_pres_time = state[Press_Type].prevTime();
    const Real dt             = parent->dtLevel(level);
    MultiFab*  DeltaSsync     = 0;// hold (Delta rho)*q for conserved quantities
    MultiFab*  Rh             = get_rho_half_time();

    sync_setup(DeltaSsync);
    //
    // Compute the u_mac for the correction.
    //
    mac_projector->mac_sync_solve(level,dt,Rh,fine_ratio);
    //
    // Update coarse grid state by adding correction from mac_sync solve
    // the correction is the advective tendency of the new velocities.
    //
    if (do_reflux)
    {
        MultiFab& S_new = get_new_data(State_Type);
        mac_projector->mac_sync_compute(level,u_mac,Vsync,Ssync,Rh,
                                        level > 0 ? &getAdvFluxReg(level) : 0,
                                        advectionType, prev_time,
                                        prev_pres_time,dt,
                                        NUM_STATE,be_cn_theta, 
                                        modify_reflux_normal_vel,
                                        do_mom_diff);
        //
        // The following used to be done in mac_sync_compute.  Ssync is
        //   the source for a rate of change to S over the time step, so
        //   Ssync*dt is the source to the actual sync amount.
        //
        Ssync->mult(dt,Ssync->nGrow());
        //
        // For all conservative variables Q (other than density)
        // express Q as rho*q and increment sync by -(sync_for_rho)*q
        // (See Pember, et. al., LBNL-41339, Jan. 1989)
        //
        FArrayBox delta_ssync;

        int iconserved = -1;
        for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
        {
            if (istate != Density && advectionType[istate] == Conservative)
            {
                iconserved++;

                for (MFIter Smfi(S_new); Smfi.isValid(); ++Smfi)
                {
                    const int  i   = Smfi.index();
                    const Box& grd = grids[i];

                    delta_ssync.resize(grd,1);
                    delta_ssync.copy(S_new[Smfi], grd, istate, grd, 0, 1);
                    delta_ssync.divide(S_new[Smfi], grd, Density, 0, 1);
                    delta_ssync.mult((*Ssync)[Smfi],grd,Density-BL_SPACEDIM,0,1);
                    (*DeltaSsync)[Smfi].copy(delta_ssync,grd,0,grd,iconserved,1);
                    (*Ssync)[Smfi].minus(delta_ssync,grd,0,istate-BL_SPACEDIM,1);
                }
            }
        }

        if (do_mom_diff == 1)
        {
            for (MFIter Vsyncmfi(*Vsync); Vsyncmfi.isValid(); ++Vsyncmfi)
            {
                const int i = Vsyncmfi.index();

                D_TERM((*Vsync)[i].divide((*rho_ctime)[i],rho_ctime->box(i),0,Xvel,1);,
                       (*Vsync)[i].divide((*rho_ctime)[i],rho_ctime->box(i),0,Yvel,1);,
                       (*Vsync)[i].divide((*rho_ctime)[i],rho_ctime->box(i),0,Zvel,1););
            }
        }
        //
        // Compute viscous sync.
        //
        if (is_diffusive[Xvel])
        {
            int rho_flag = (do_mom_diff == 0) ? 1 : 3;

            MultiFab** loc_viscn = 0;

            if (variable_vel_visc)
            {
                Real viscTime = state[State_Type].prevTime();
                diffusion->allocFluxBoxesLevel(loc_viscn, 0, 1);
                getViscosity(loc_viscn, viscTime);
            }

            diffusion->diffuse_Vsync(Vsync,dt,be_cn_theta,Rh,rho_flag,loc_viscn);

            if (variable_vel_visc)
            {
                diffusion->removeFluxBoxesLevel(loc_viscn);
            }
        }

        MultiFab** fluxSC        = 0;
        bool       any_diffusive = false;
        for (int sigma  = 0; sigma < numscal; sigma++)
            if (is_diffusive[BL_SPACEDIM+sigma])
                any_diffusive = true;

        if (any_diffusive)
            diffusion->allocFluxBoxesLevel(fluxSC,0,1);

        for (int sigma = 0; sigma<numscal; sigma++)
        {
            const int state_ind = BL_SPACEDIM + sigma;
            const int rho_flag  = Diffusion::set_rho_flag(diffusionType[state_ind]);

            if (is_diffusive[state_ind])
            {
                MultiFab** cmp_diffn=0;

                if (variable_scal_diff)
                {
                    Real diffTime = state[State_Type].prevTime();
                    diffusion->allocFluxBoxesLevel(cmp_diffn, 0, 1);
                    getDiffusivity(cmp_diffn, diffTime, BL_SPACEDIM+sigma,0,1);
                }

                diffusion->diffuse_Ssync(Ssync,sigma,dt,be_cn_theta,
                                         Rh,rho_flag,fluxSC,0,cmp_diffn);

                if (variable_scal_diff)
                {
                    diffusion->removeFluxBoxesLevel(cmp_diffn);
                }
                //
                // Increment the viscous flux registers
                //
                if (level > 0)
                {
                    for (int d = 0; d < BL_SPACEDIM; d++)
                    {
                        Real mult = dt;
                        MultiFab& fluxSCd = *fluxSC[d];
                        for (MFIter fmfi(fluxSCd); fmfi.isValid(); ++fmfi)
                            getViscFluxReg().FineAdd(fluxSCd[fmfi],d,
                                                     fmfi.index(),
                                                     0,state_ind,1,mult);
                    }
                }
            }
        }

        if (any_diffusive)
            diffusion->removeFluxBoxesLevel(fluxSC);
        //
        // For all conservative variables Q (other than density)
        // increment sync by (sync_for_rho)*q_presync.
        // (See Pember, et. al., LBNL-41339, Jan. 1989)
        //
        iconserved = -1;
        for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
        {
            if (istate != Density && advectionType[istate] == Conservative)
            {
                iconserved++;

                for (MFIter SsyncMfi(*Ssync); SsyncMfi.isValid(); ++SsyncMfi)
                {
                    const int i = SsyncMfi.index();
                    (*Ssync)[i].plus((*DeltaSsync)[i], grids[i],
                                     iconserved, istate-BL_SPACEDIM, 1);
                }
            }
        }
        //
        // Add the sync correction to the state.
        //
        for (int sigma  = 0; sigma < numscal; sigma++)
        {
            for (MFIter S_newmfi(S_new); S_newmfi.isValid(); ++S_newmfi)
            {
                S_new[S_newmfi].plus((*Ssync)[S_newmfi],S_newmfi.validbox(),
                                     sigma,BL_SPACEDIM+sigma,1);
            }
        }
        //
        // Update rho_ctime after rho is updated with Ssync.
        //
        make_rho_curr_time();
        if (level > 0) incrRhoAvg((*Ssync),Density-BL_SPACEDIM,1.0);
        //
        // Get boundary conditions.
        //
        Array<int*>         sync_bc(grids.size());
        Array< Array<int> > sync_bc_array(grids.size());

        for (int i = 0; i < grids.size(); i++)
        {
            sync_bc_array[i] = getBCArray(State_Type,i,Density,numscal);
            sync_bc[i]       = sync_bc_array[i].dataPtr();
        }
        //
        // Interpolate the sync correction to the finer levels,
        //  and update rho_ctime, rhoAvg at those levels.
        //
        IntVect    ratio = IntVect::TheUnitVector();
        const Real mult  = 1.0;
        for (int lev = level+1; lev <= parent->finestLevel(); lev++)
        {
            ratio                     *= parent->refRatio(lev-1);
            NavierStokes&     fine_lev = getLevel(lev);
            const BoxArray& fine_grids = fine_lev.boxArray();
            MultiFab sync_incr(fine_grids,numscal,0);
            sync_incr.setVal(0.0);

            SyncInterp(*Ssync,level,sync_incr,lev,ratio,0,0,
                       numscal,1,mult,sync_bc.dataPtr());

            MultiFab& S_new = fine_lev.get_new_data(State_Type);
            for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
                S_new[mfi].plus(sync_incr[mfi],fine_grids[mfi.index()],0,Density,numscal);

            fine_lev.make_rho_curr_time();
            fine_lev.incrRhoAvg(sync_incr,Density-BL_SPACEDIM,1.0);
        }
    }

    sync_cleanup(DeltaSsync);
}

void
NavierStokes::sync_setup (MultiFab*& DeltaSsync)
{
    BL_ASSERT(DeltaSsync == 0);

    int nconserved = Godunov::how_many(advectionType, Conservative,
                                       BL_SPACEDIM, NUM_STATE-BL_SPACEDIM);

    if (nconserved > 0 && level < parent->finestLevel())
    {
        DeltaSsync = new MultiFab(grids, nconserved, 1, Fab_allocate);
        DeltaSsync->setVal(0,1);
    }
}

void
NavierStokes::sync_cleanup (MultiFab*& DeltaSsync)
{
    delete DeltaSsync;

    DeltaSsync = 0;
}


//
// The reflux function
//

void
NavierStokes::reflux ()
{
    if (level == parent->finestLevel())
        return;

    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::reflux()");

    BL_ASSERT(do_reflux);

    MultiFab& S_crse = get_new_data(State_Type);
    //
    // First do refluxing step.
    //
    FluxRegister& fr_adv  = getAdvFluxReg(level+1);
    FluxRegister& fr_visc = getViscFluxReg(level+1);
    Real          dt_crse = parent->dtLevel(level);
    Real          scale   = 1.0/dt_crse;
    //
    // It is important, for do_mom_diff == 0, to do the viscous
    //   refluxing first, since this will be divided by rho_half
    //   before the advective refluxing is added.  In the case of
    //   do_mom_diff == 1, both components of the refluxing will
    //   be divided by rho^(n+1) in level_sync.
    //
    fr_visc.Reflux(*Vsync,volume,scale,0,0,BL_SPACEDIM,geom);
    fr_visc.Reflux(*Ssync,volume,scale,BL_SPACEDIM,0,NUM_STATE-BL_SPACEDIM,geom);

    const MultiFab* Rh = get_rho_half_time();

    if (do_mom_diff == 0)
    {
        for (MFIter Vsyncmfi(*Vsync); Vsyncmfi.isValid(); ++Vsyncmfi)
        {
            const int i = Vsyncmfi.index();

            D_TERM((*Vsync)[i].divide((*Rh)[i],Rh->box(i),0,Xvel,1);,
                   (*Vsync)[i].divide((*Rh)[i],Rh->box(i),0,Yvel,1);,
                   (*Vsync)[i].divide((*Rh)[i],Rh->box(i),0,Zvel,1););
        }
    }

    for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
    {
        if (advectionType[istate] == NonConservative)
        {
            const int sigma = istate -  BL_SPACEDIM;

            for (MFIter Ssyncmfi(*Ssync); Ssyncmfi.isValid(); ++Ssyncmfi)
            {
                const int i = Ssyncmfi.index();

                (*Ssync)[i].divide((*Rh)[i],Rh->box(i),0,sigma,1);
            }
        }
    }

    fr_adv.Reflux(*Vsync,volume,scale,0,0,BL_SPACEDIM,geom);
    fr_adv.Reflux(*Ssync,volume,scale,BL_SPACEDIM,0,NUM_STATE-BL_SPACEDIM,geom);

    const BoxArray& fine_boxes = getLevel(level+1).boxArray();
    const int       nfine      = fine_boxes.size();
    //
    // This is necessary in order to zero out the contribution to any
    // coarse grid cells which underlie fine grid cells.
    //
    for (int kf = 0; kf < nfine; kf++)
    {
        Box bf = BoxLib::coarsen(fine_boxes[kf],fine_ratio);

        for (MFIter Vsyncmfi(*Vsync); Vsyncmfi.isValid(); ++Vsyncmfi)
        {
            const int i = Vsyncmfi.index();

            BL_ASSERT(grids[i] == Vsyncmfi.validbox());

            Box bx = bf & Vsyncmfi.validbox();

            if (bx.ok())
            {
                (*Vsync)[i].setVal(0,bx,0,BL_SPACEDIM);
                (*Ssync)[i].setVal(0,bx,0,NUM_STATE-BL_SPACEDIM);
            }
        }
    }
}

//
// Average down a single state component.
//

void
NavierStokes::avgDown (int comp)
{
    if (level == parent->finestLevel())
        return;

    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::avgDown()");

    NavierStokes&   fine_lev = getLevel(level+1);
    const BoxArray& fgrids   = fine_lev.grids;
    MultiFab&       fvolume  = fine_lev.volume;
    MultiFab&       S_crse   = get_new_data(State_Type);
    MultiFab&       S_fine   = fine_lev.get_new_data(State_Type);

    avgDown(grids,fgrids,S_crse,S_fine,volume,fvolume,level,level+1,comp,1,fine_ratio);

    if (comp == Density) 
    {
        //
        // Fill rho_ctime at current and finer levels with the correct data.
        //
        for (int lev = level; lev <= parent->finestLevel(); lev++)
        {
            getLevel(lev).make_rho_curr_time();
        }
    }
}

//
// Inject fine pressure nodes down onto coarse nodes.
//

void
NavierStokes::injectDown (const Box&       ovlp,
                          FArrayBox&       Pcrse,
                          const FArrayBox& Pfine,
                          IntVect&         fine_ratio )
{
    const int*  ovlo  = ovlp.loVect();
    const int*  ovhi  = ovlp.hiVect();
    Real*       cpres = Pcrse.dataPtr();
    const int*  clo   = Pcrse.loVect();
    const int*  chi   = Pcrse.hiVect();
    const Real* fpres = Pfine.dataPtr();
    const int*  flo   = Pfine.loVect();
    const int*  fhi   = Pfine.hiVect();

    FORT_PUTDOWN(cpres,ARLIM(clo),ARLIM(chi),
                 fpres,ARLIM(flo),ARLIM(fhi),
                 ovlo,ovhi,fine_ratio.getVect());
}

//
// Test for consistency between fine and coarse nodes.
//

void
NavierStokes::testInject (const Box&       ovlp,
                          FArrayBox&       Pcrse,
                          const FArrayBox& Pfine,
                          IntVect&         fine_ratio)
{
    const int*  ovlo  = ovlp.loVect();
    const int*  ovhi  = ovlp.hiVect();
    Real*       cpres = Pcrse.dataPtr();
    const int*  clo   = Pcrse.loVect();
    const int*  chi   = Pcrse.hiVect();
    const Real* fpres = Pfine.dataPtr();
    const int*  flo   = Pfine.loVect();
    const int*  fhi   = Pfine.hiVect();

    FORT_TESTINJECT(cpres,ARLIM(clo),ARLIM(chi),
                    fpres,ARLIM(flo),ARLIM(fhi),
                    ovlo,ovhi,fine_ratio.getVect());
}

//
// Average fine information from the complete set of state types to coarse.
//

void
NavierStokes::avgDown ()
{
    if (level == parent->finestLevel())
        return;

    NavierStokes&   fine_lev = getLevel(level+1);
    const BoxArray& fgrids   = fine_lev.grids;
    MultiFab&       fvolume  = fine_lev.volume;
    //
    // Average down the states at the new time.
    //
    MultiFab& S_crse = get_new_data(State_Type);
    MultiFab& S_fine = fine_lev.get_new_data(State_Type);

    avgDown(grids,fgrids,S_crse,S_fine,volume,fvolume,level,level+1,0,S_crse.nComp(), fine_ratio);
    //
    // Now average down pressure over time n-(n+1) interval.
    //
    MultiFab&       P_crse      = get_new_data(Press_Type);
    MultiFab&       P_fine_init = fine_lev.get_new_data(Press_Type);
    MultiFab&       P_fine_avg  = *fine_lev.p_avg;
    MultiFab&       P_fine      = initial_step ? P_fine_init : P_fine_avg;
    const BoxArray& P_cgrids    = state[Press_Type].boxArray();
    const BoxArray& P_fgrids    = fine_lev.state[Press_Type].boxArray();

    BoxArray crse_P_fine_BA(P_fgrids.size());

    for (int i = 0; i < P_fgrids.size(); ++i)
    {
        crse_P_fine_BA.set(i,BoxLib::coarsen(P_fgrids[i],fine_ratio));
    }
    MultiFab crse_P_fine(crse_P_fine_BA,1,0);

    for (MFIter mfi(P_fine); mfi.isValid(); ++mfi)
    {
        const int i = mfi.index();

        injectDown(crse_P_fine_BA[i],crse_P_fine[i],P_fine[i],fine_ratio);
    }
    P_crse.copy(crse_P_fine);
    //
    // Next average down divu and dSdT at new time.
    //
    if (have_divu)
    {
        MultiFab& Divu_crse = get_new_data(Divu_Type);
        MultiFab& Divu_fine = fine_lev.get_new_data(Divu_Type);
        
        avgDown(grids,fgrids,Divu_crse,Divu_fine,volume,fvolume,level,level+1,0,1,fine_ratio);
    }
    if (have_dsdt)
    {
        MultiFab& Dsdt_crse = get_new_data(Dsdt_Type);
        MultiFab& Dsdt_fine = fine_lev.get_new_data(Dsdt_Type);
        
        avgDown(grids,fgrids,Dsdt_crse,Dsdt_fine,volume,fvolume,level,level+1,0,1,fine_ratio);
    }
    //
    // Fill rho_ctime at the current and finer levels with the correct data.
    //
    for (int lev = level; lev <= parent->finestLevel(); lev++)
    {
        getLevel(lev).make_rho_curr_time();
    }
}

//
// ACCESS FUNCTIONS FOLLOW
//

//
// Virtual access function for getting the advective flux out of the
// advection routines for diagnostics and refluxing.
//

void
NavierStokes::pullFluxes (int        i,
                          int        start_ind,
                          int        ncomp,
                          FArrayBox& xflux,
                          FArrayBox& yflux,
                          FArrayBox& zflux,
                          Real       dt)
{
    //
    // Add fluxes into the refluxing counters.
    //
    if (do_reflux)
    {
        if (level < parent->finestLevel())
        {
            FluxRegister& fr = getAdvFluxReg(level+1);
            if (!modify_reflux_normal_vel || start_ind != Xvel)
                fr.CrseInit(xflux,xflux.box(),0,0,start_ind,ncomp,-dt);
            if (!modify_reflux_normal_vel || start_ind != Yvel)
                fr.CrseInit(yflux,yflux.box(),1,0,start_ind,ncomp,-dt);
#if (BL_SPACEDIM == 3)                              
            if (!modify_reflux_normal_vel || start_ind != Zvel)
                fr.CrseInit(zflux,zflux.box(),2,0,start_ind,ncomp,-dt);
#endif
        }
        if (level > 0)
        {
            if (!modify_reflux_normal_vel || start_ind != Xvel)
                advflux_reg->FineAdd(xflux,0,i,0,start_ind,ncomp,dt);
            if (!modify_reflux_normal_vel || start_ind != Yvel)
                advflux_reg->FineAdd(yflux,1,i,0,start_ind,ncomp,dt);
#if (BL_SPACEDIM == 3)                                
            if (!modify_reflux_normal_vel || start_ind != Zvel)
                advflux_reg->FineAdd(zflux,2,i,0,start_ind,ncomp,dt);
#endif
        }
    }
}

//
// Virtual access function for getting the forcing terms for the
// velocities and scalars.  The base version computes a buoyancy.
//
// As NavierStokes is currently implemented.  Velocities are integrated
// according to the equation
//
//     ui_t + uj ui_j = S_ui        ===> tforces = rho S_ui
//
// and scalars psi where (psi = rho q) as
//
//     psi_t + (uj psi)_j = S_psi   ===> tforces = S_psi = rho S_q
//
// q is a concentration.  This function returns a rho weighted
// source term, which requires a division by rho in the predict_velocity
// and velocity_advection routines.
//
#ifdef GENGETFORCE
void
NavierStokes::getForce (FArrayBox&       force,
                        int              gridno,
                        int              ngrow,
                        int              scomp,
                        int              ncomp,
                        const Real       time,
                        const FArrayBox& Rho)
{
    BL_ASSERT(Rho.nComp() == 1);

    force.resize(BoxLib::grow(grids[gridno],ngrow),ncomp);

    BL_ASSERT(Rho.box().contains(force.box()));

    const Real* dx       = geom.CellSize();
    const Real  grav     = gravity;
    const int*  s_lo     = Rho.loVect();
    const int*  s_hi     = Rho.hiVect();
    const int*  f_lo     = force.loVect();
    const int*  f_hi     = force.hiVect();

    FORT_MAKEFORCE (&time,
		    force.dataPtr(),
		    Rho.dataPtr(),
		    ARLIM(f_lo), ARLIM(f_hi),
		    ARLIM(s_lo), ARLIM(s_hi),
		    dx,
		    grid_loc[gridno].lo(),
		    grid_loc[gridno].hi(),
		    &grav,&scomp,&ncomp);

}
#else
void
NavierStokes::getForce (FArrayBox&       force,
                        int              gridno,
                        int              ngrow,
                        int              scomp,
                        int              ncomp,
                        const FArrayBox& Rho)
{
    BL_ASSERT(Rho.nComp() == 1);

    force.resize(BoxLib::grow(grids[gridno],ngrow),ncomp);

    BL_ASSERT(Rho.box().contains(force.box()));

    const Real grav = gravity;

    for (int dc = 0; dc < ncomp; dc++)
    {
        const int sc = scomp + dc;
#if (BL_SPACEDIM == 2)
        if (sc == Yvel && std::abs(grav) > 0.001) 
#endif
#if (BL_SPACEDIM == 3)
        if (sc == Zvel && std::abs(grav) > 0.001) 
#endif
        {
            //
            // Set force to -rho*g.
            //
            force.copy(Rho,0,dc,1);
            force.mult(grav,dc,1);
        }
        else
        {
            force.setVal(0,dc);
        }
    }
}
// Generalised getForce
#endif

void
NavierStokes::getGradP (MultiFab& gp,
                        Real      time)
{
    const int NGrow = gp.nGrow();
    MultiFab& P_old = get_old_data(Press_Type);

    const Real* dx             = geom.CellSize();

    if (level > 0 && state[Press_Type].descriptor()->timeType() ==
                     StateDescriptor::Point)
    {
        //
        // We want to be sure the intersection of old and new grids is
        // entirely contained within gp.boxArray()
        //
        BL_ASSERT(gp.boxArray() == grids);

        {
            //
            // Build MultiFab whose valid region encompasses NGrow grow cells.
            // The valid region of the MultiFab will contain overlaps!
            //
            const BoxArray& pBA = state[Press_Type].boxArray();

            BoxArray ovlpBA(pBA.size());

            for (int j = 0; j < ovlpBA.size(); j++)
                ovlpBA.set(j,BoxLib::grow(pBA[j],NGrow));

            MultiFab pMF(ovlpBA,1,0);
            MultiFab dpdtMF(ovlpBA,1,0);

            if (time == getLevel(level-1).state[Press_Type].prevTime() || 
                time == getLevel(level-1).state[Press_Type].curTime())
            {
                FillCoarsePatch(pMF,0,time,Press_Type,0,1);
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
  
                FillCoarsePatch(pMF,0,crse_time,Press_Type,0,1);
  
                FillCoarsePatch(dpdtMF,0,time,Dpdt_Type,0,1);

                Real dt_temp = time - crse_time;

                dpdtMF.mult(dt_temp);

                pMF.plus(dpdtMF,0,1,0);
            }

            for (MFIter mfi(pMF); mfi.isValid(); ++mfi) 
            {
                Projection::getGradP(pMF[mfi],gp[mfi],gp[mfi].box(),dx);
            }
        }
        //
        // We've now got good coarse data everywhere in gp.
        // FillPatch temp version of gp having overlapping valid regions.
        //
        BoxArray ovlpBA(gp.boxArray().size());

        for (int j = 0; j < gp.boxArray().size(); j++)
            ovlpBA.set(j,BoxLib::grow(gp.boxArray()[j],NGrow));

        MultiFab gpTmp(ovlpBA,gp.nComp(),0);

        for (FillPatchIterator P_fpi(*this,P_old,NGrow,time,Press_Type,0,1);
             P_fpi.isValid();
             ++P_fpi) 
        {
            Projection::getGradP(P_fpi(),gpTmp[P_fpi],gpTmp[P_fpi].box(),dx);
        }
        //
        // Now must decide which parts of gpTmp to copy to gp.
        //
        BoxArray fineBA(old_intersect_new.size());

        for (int j = 0; j < old_intersect_new.size(); j++)
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

        for (MFIter mfi(gpTmp); mfi.isValid(); ++mfi) 
        {
            for (int j = 0; j < fineBA.size(); j++)
            {
                Box isect = fineBA[j] & gpTmp[mfi].box();

                if (isect.ok())
                {
                    gp[mfi].copy(gpTmp[mfi],isect);
                }
            }
        }

        geom.FillPeriodicBoundary(gp,true);
    }
    else
    {
        FillPatchIterator P_fpi(*this,P_old,NGrow,time,Press_Type,0,1);

        for ( ; P_fpi.isValid(); ++P_fpi) 
        {
            BL_ASSERT(BoxLib::grow(grids[P_fpi.index()],NGrow) == gp[P_fpi].box());

            Projection::getGradP(P_fpi(),gp[P_fpi],gp[P_fpi].box(),dx);
        }
    }
}

//
// Fill patch divU.
//

MultiFab*
NavierStokes::getDivCond (int  ngrow, 
                          Real time)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::getDivCond()");

    MultiFab* divu = 0;

    if (!have_divu)
    {
        divu = new MultiFab(grids,1,ngrow);

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
NavierStokes::getDsdt (int  ngrow, 
                       Real time)
{
    MultiFab* dsdt = 0;

    if (!(have_dsdt && have_divu))
    {
        dsdt = new MultiFab(grids,1,ngrow);

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
NavierStokes::getState (int  ngrow,
                        int  state_idx,
                        int  scomp,
                        int  ncomp, 
                        Real time)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::getState()");

    MultiFab* mf = new MultiFab(state[state_idx].boxArray(),ncomp,ngrow);

    FillPatchIterator fpi(*this,*mf,ngrow,time,state_idx,scomp,ncomp);

    for ( ; fpi.isValid(); ++fpi)
    {
        (*mf)[fpi.index()].copy(fpi());
    }

    return mf;
}

//
// Fills ghost cells of state:
//

void
NavierStokes::FillStateBndry (Real time,
                              int  state_idx,
                              int  src_comp, 
                              int  ncomp) 
{
    MultiFab& S = get_data(state_idx,time);

    if (S.nGrow() == 0)
        return;

    FillPatchIterator fpi(*this,S,S.nGrow(),time,state_idx,src_comp,ncomp);

    for ( ; fpi.isValid(); ++fpi)
    {
        //
        // Fill all ghost cells interior & exterior to valid region.
        //
        BoxList boxes = BoxLib::boxDiff(fpi().box(),grids[fpi.index()]);

        for (BoxList::iterator bli = boxes.begin(); bli != boxes.end(); ++bli)
        {
            S[fpi.index()].copy(fpi(),*bli,0,*bli,src_comp,ncomp);
        }
    }
}

//
// Default divU is set to zero.
//

void
NavierStokes::calc_divu (Real      time,
                         Real      dt,
                         MultiFab& divu)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::calc_divu()");

    if (have_divu)
    {
        divu.setVal(0);

        if (do_temp && visc_coef[Temp] > 0.0)
        {
            //
            // Compute Div(U) = Div(visc_cond_coef * Grad(T))/(c_p*rho*T)
            //
            getViscTerms(divu,Temp,1,time);

            const MultiFab&   rhotime = get_rho(time);
            FillPatchIterator temp_fpi(*this,divu,0,time,State_Type,Temp,1);
            MFIter            rho_mfi(rhotime);

            for ( ;
                  rho_mfi.isValid() && temp_fpi.isValid();
                  ++rho_mfi, ++temp_fpi)
            {
                const int i = rho_mfi.index();

                divu[i].divide(rhotime[i],divu.box(i),0,0,1);
                divu[i].divide(temp_fpi(),divu.box(i),0,0,1);
            }
            Real THERMO_cp_inv = 1.0 / 1004.6;
            divu.mult(THERMO_cp_inv);
        }
    }
}

//
// Default dSdt is set to zero.
//

void
NavierStokes::calc_dsdt (Real      time,
                         Real      dt,
                         MultiFab& dsdt)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::calc_dsdt()");

    if (have_divu && have_dsdt)
    {
        dsdt.setVal(0);

        if (do_temp)
        {
            MultiFab& Divu_new = get_new_data(Divu_Type);
            MultiFab& Divu_old = get_old_data(Divu_Type);

            for (MFIter mfi(dsdt); mfi.isValid(); ++mfi)
            {
                dsdt[mfi].copy(Divu_new[mfi],mfi.validbox(),0,mfi.validbox(),0,1);
                dsdt[mfi].minus(Divu_old[mfi],mfi.validbox(),0,0,1);
                dsdt[mfi].divide(dt,mfi.validbox(),0,1);
            }
        }
    }
}

void
NavierStokes::getViscTerms (MultiFab& visc_terms,
                            int       src_comp, 
                            int       ncomp,
                            Real      time)
{
    //
    // The logic below for selecting between scalar or tensor solves does 
    // not allow for calling NavierStokes::getViscTerms with src_comp=Yvel
    // or Zvel
    //
#ifndef NDEBUG
    if (src_comp<BL_SPACEDIM && (src_comp!=Xvel || ncomp<BL_SPACEDIM))
    {
        std::cout << "src_comp=" << src_comp << "   ncomp=" << ncomp << std::endl;
        BoxLib::Error("must call NavierStokes::getViscTerms with all three velocity components");
    }
#endif
    // 
    // Initialize all viscous terms to zero
    //
    const int nGrow = visc_terms.nGrow();
    visc_terms.setVal(0,0,ncomp,nGrow);
    //
    // 
    // Get Velocity Viscous Terms
    //
    if (src_comp == Xvel && is_diffusive[Xvel])
    {
        MultiFab** viscosity;

        if (variable_vel_visc)
        {
            diffusion->allocFluxBoxesLevel(viscosity, 0, 1);
            getViscosity(viscosity, time);

            diffusion->getTensorViscTerms(visc_terms,time,0,viscosity);
        }
        else
        {
            for (int icomp = Xvel; icomp < BL_SPACEDIM; icomp++)
            {
                int rho_flag = Diffusion::set_rho_flag(diffusionType[icomp]);

                diffusion->getViscTerms(visc_terms,src_comp,icomp,time,rho_flag);
            }
        }
        //
        // Add Div(u) term if desired, if this is velocity, and if Div(u) 
        // is nonzero.  If const-visc, term is mu.Div(u)/3, else 
        // it's -Div(mu.Div(u).I)*2/3
        //
        if (have_divu && S_in_vel_diffusion)
        {
            MultiFab divmusi(grids,BL_SPACEDIM,1);

            if (variable_vel_visc)
            {
                diffusion->compute_divmusi(time,viscosity,divmusi);
                divmusi.mult((-2./3.),0,BL_SPACEDIM,0);
            }
            else
            {
                diffusion->compute_divmusi(time,visc_coef[Xvel],divmusi);
                divmusi.mult((1./3.),0,BL_SPACEDIM,0);
            }

            visc_terms.plus(divmusi,Xvel,BL_SPACEDIM,0);
        }
        //
        // Clean up variable viscosity arrays
        //
        if (variable_vel_visc)
            diffusion->removeFluxBoxesLevel(viscosity);
    }
    //
    // Get Scalar Diffusive Terms
    //
    const int first_scal = (src_comp==Xvel) ? BL_SPACEDIM : src_comp;
    const int num_scal   = (src_comp==Xvel) ? ncomp-BL_SPACEDIM : ncomp;

    if (num_scal > 0)
    {
        for (int icomp = first_scal; icomp < first_scal+num_scal; icomp++)
        {
            if (is_diffusive[icomp])
            {
                int rho_flag = Diffusion::set_rho_flag(diffusionType[icomp]);

                MultiFab** cmp_diffn = 0;

                if (variable_scal_diff)
                {
                    diffusion->allocFluxBoxesLevel(cmp_diffn, 0, 1);
                    getDiffusivity(cmp_diffn, time, icomp, 0, 1);
                }

                diffusion->getViscTerms(visc_terms,src_comp,icomp,
                                        time,rho_flag,0,cmp_diffn);

                if (variable_scal_diff)
                {
                    diffusion->removeFluxBoxesLevel(cmp_diffn);
                }
            }
        }
    }
    //
    // Ensure consistent grow cells
    //    
    if (nGrow > 0)
    {
        for (MFIter mfi(visc_terms); mfi.isValid(); ++mfi)
        {
            FArrayBox& vt  = visc_terms[mfi];
            const Box& box = mfi.validbox();
            FORT_VISCEXTRAP(vt.dataPtr(),ARLIM(vt.loVect()),ARLIM(vt.hiVect()),
                            box.loVect(),box.hiVect(),&ncomp);
        }
        visc_terms.FillBoundary(0,ncomp);
        //
        // Note: this is a special periodic fill in that we want to
        // preserve the extrapolated grow values when periodic --
        // usually we preserve only valid data.  The scheme relies on
        // the fact that there is good data in the "non-periodic" grow cells.
        // ("good" data produced via VISCEXTRAP above)
        //
        geom.FillPeriodicBoundary(visc_terms,0,ncomp,true);
    }
}

//
// Functions for calculating the variable viscosity and diffusivity.
// These default to setting the variable viscosity and diffusivity arrays
// to the values in visc_coef and diff_coef.  These functions would
// need to be replaced in any class derived from NavierStokes that
// wants variable coefficients.
//
void 
NavierStokes::calcViscosity (const Real time, 
                             const Real dt,
                             const int  iteration,
                             const int  ncycle)
{
    const MultiFab& S = get_data(State_Type, time);
    //
    // Select time level to work with (N or N+1)
    //
    MultiFab* visc_cc;

    const TimeLevel whichTime = which_time(State_Type,time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    if (whichTime == AmrOldTime)               // time N
    {
        visc_cc = viscn_cc;
    }
    else if (whichTime == AmrNewTime)          // time N+1
    {
        visc_cc = viscnp1_cc;
    }
    //
    // Calculate viscosity
    //
    const int nGrow = visc_cc->nGrow();

    if (is_diffusive[Xvel]) 
    {
        if (visc_coef[Xvel] >= 0.0)
        {
            visc_cc->setVal(visc_coef[Xvel], 0, 1, nGrow);
        }
        else
        {
            BoxLib::Abort("NavierStokes::calcViscosity() : must have velocity visc_coef >= 0.0");
        }
    }
}

void 
NavierStokes::calcDiffusivity (const Real time, 
                               const Real dt,
                               const int  iteration,
                               const int  ncycle,
                               const int  src_comp, 
                               const int  ncomp)
{
    //
    // NOTE:  The component numbers passed into NavierStokes::calcDiffusivity
    //        correspond to the components in the state.  In the diffusivity 
    //        arrays, there is an offset since no diffusivity array
    //        is kept for the velocities or the density.  So, the scalar
    //        component Density+1 in the state corresponds to component
    //        0 in the arrays diffn and diffnp1.
    //
    BL_ASSERT(src_comp > Density);

    const MultiFab& S = get_data(State_Type, time);
    //
    // Select time level to work with (N or N+1)
    //
    MultiFab* diff_cc;

    const TimeLevel whichTime = which_time(State_Type,time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    if (whichTime == AmrOldTime)               // time N
    {
        diff_cc = diffn_cc;
    }
    else if (whichTime == AmrNewTime)          // time N+1
    {
        diff_cc = diffnp1_cc;
    }
    //
    // Calculate diffusivity
    //
    const int nGrow = diff_cc->nGrow();

    for (int comp=src_comp; comp<src_comp+ncomp; comp++)
    {
        int diff_comp = comp - Density - 1;

        if (is_diffusive[comp])
        {
            if (visc_coef[comp] >= 0.0)
            {
                diff_cc->setVal(visc_coef[comp], diff_comp, 1, nGrow);
            }
            else
            {
                BoxLib::Abort("NavierStokes::calcDiffusivity() : must have scalar diff_coefs >= 0.0");
            }
        }
    }
}


void 
NavierStokes::getViscosity (MultiFab* viscosity[BL_SPACEDIM],
                            const Real time)
{
    //
    // Select time level to work with (N or N+1)
    //
    MultiFab* visc_cc;

    const TimeLevel whichTime = which_time(State_Type,time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    if (whichTime == AmrOldTime)               // time N
    {
        visc_cc = viscn_cc;
    }
    else if (whichTime == AmrNewTime)          // time N+1
    {
        visc_cc = viscnp1_cc;
    }
    //
    // Fill edge-centered viscosity
    //
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        for (MFIter ecMfi(*viscosity[dir]); ecMfi.isValid(); ++ecMfi)
        {
            const int i = ecMfi.index();

            center_to_edge_plain((*visc_cc)[i],(*viscosity[dir])[i],0,0,1);
        }
    }
}

void 
NavierStokes::getDiffusivity (MultiFab* diffusivity[BL_SPACEDIM],
                              const Real time,
                              const int state_comp,
                              const int dst_comp,
                              const int ncomp)
{
    BL_ASSERT(state_comp > Density);
    //
    // Pick correct diffusivity component
    //
    int diff_comp = state_comp - Density - 1;
    //
    // Select time level to work with (N or N+1)
    //
    MultiFab* diff_cc;

    const TimeLevel whichTime = which_time(State_Type,time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    if (whichTime == AmrOldTime)               // time N
    {
        diff_cc = diffn_cc;
    }
    else if (whichTime == AmrNewTime)          // time N+1
    {
        diff_cc = diffnp1_cc;
    }
    //
    // Fill edge-centered diffusivities
    //
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        for (MFIter ecMfi(*diffusivity[dir]); ecMfi.isValid(); ++ecMfi)
        {
            center_to_edge_plain((*diff_cc)[ecMfi],(*diffusivity[dir])[ecMfi],
                                 diff_comp,dst_comp,ncomp);
        }
    }
}

void
NavierStokes::center_to_edge_plain (const FArrayBox& ccfab,
                                    FArrayBox&       ecfab,
                                    int              sComp,
                                    int              dComp,
                                    int              nComp)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::center_to_edge_plain()");

    //
    // This routine fills an edge-centered FAB from a cell-centered FAB.
    // It assumes that the data in all cells of the cell-centered FAB is
    // valid and totally ignores any concept of boundary conditions.  
    // It is assummed that the cell-centered FAB fully contains the 
    // edge-centered FAB.  If anything special needs to be done at boundaries, 
    // a varient of this routine needs to be written.  See 
    // HeatTransfer::center_to_edge_fancy().
    //
    const Box&      ccbox = ccfab.box();
    const Box&      ecbox = ecfab.box();
    const IndexType ixt   = ecbox.ixType();
    //
    // Get direction for interpolation to edges
    //
    int dir = -1;
    for (int d = 0; d < BL_SPACEDIM; d++)
        if (ixt.test(d))
            dir = d;
    //
    // Miscellanious checks
    //
    BL_ASSERT(!(ixt.cellCentered()) && !(ixt.nodeCentered()));
    BL_ASSERT(BoxLib::grow(ccbox,-BoxLib::BASISV(dir)).contains(BoxLib::enclosedCells(ecbox)));
    BL_ASSERT(sComp+nComp <= ccfab.nComp() && dComp+nComp <= ecfab.nComp());
    //
    // Shift cell-centered data to edges
    //
    Box fillBox = ccbox; 
    for (int d = 0; d < BL_SPACEDIM; d++)
        if (d != dir)
            fillBox.setRange(d, ecbox.smallEnd(d), ecbox.length(d));
    
    const int isharm = def_harm_avg_cen2edge;
    FORT_CEN2EDG(fillBox.loVect(), fillBox.hiVect(),
                 ARLIM(ccfab.loVect()), ARLIM(ccfab.hiVect()),
                 ccfab.dataPtr(sComp),
                 ARLIM(ecfab.loVect()), ARLIM(ecfab.hiVect()),
                 ecfab.dataPtr(dComp),
                 &nComp, &dir, &isharm);
}

//
// Note: this is a temporary function.  Eventually this will be moved to a
// boundary condition class.
//

static
void
getOutFlowFaces (Array<Orientation>& outFaces,
                 BCRec*              _phys_bc)
{
    outFaces.resize(0);
    for (int idir = 0; idir < BL_SPACEDIM; idir++)
    {
        if (_phys_bc->lo(idir) == Outflow)
        {
            const int len = outFaces.size();
            outFaces.resize(len+1);
            outFaces.set(len,Orientation(idir,Orientation::low));
        }

        if (_phys_bc->hi(idir) == Outflow)
        {
            const int len = outFaces.size();
            outFaces.resize(len+1);
            outFaces.set(len,Orientation(idir,Orientation::high));
        }
    }
}

void
NavierStokes::manual_tags_placement (TagBoxArray&    tags,
                                     Array<IntVect>& bf_lev)
{
    Array<Orientation> outFaces;
    getOutFlowFaces(outFaces,&phys_bc);
    if (outFaces.size()>0)
    {
        for (int i=0; i<outFaces.size(); ++i)
        {
            const Orientation& outFace = outFaces[i];
            const int oDir = outFace.coordDir();
            const Box& crse_domain = BoxLib::coarsen(geom.Domain(),bf_lev[level]);
            const int mult = (outFace.isLow() ? +1 : -1);
            if (do_refine_outflow)
            {
                //
                // Refine entire outflow boundary if new boxes within grid_tol
                // from outflow
                //
                const int grid_tol = 1;
                const Box outflowBox =
                    BoxLib::adjCell(crse_domain,outFace,grid_tol).shift(oDir,mult*grid_tol);
                //
                // Only refine if there are already tagged cells in the outflow
                // region
                //
                bool hasTags = false;
                for (MFIter tbi(tags); !hasTags && tbi.isValid(); ++tbi)
                    if (tags[tbi].numTags(outflowBox) > 0)
                        hasTags = true;
                
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
                for (int i = 1; i <= level; ++i)
                {
                    /*** Calculate the minimum cells at this level ***/
                    
                    const int rat = (parent->refRatio(i-1))[oDir];
                    N_level_cells = N_level_cells * rat + np;
                    
                    /*** Calculate the required number of coarse cells ***/
                    
                    N_coarse_cells = N_level_cells / bf_lev[i][oDir];
                    if (N_level_cells % bf_lev[i][oDir] != 0)
                        N_coarse_cells++;
                    
                    /*** Calculate the corresponding number of level cells ***/
                    
                    N_level_cells = N_coarse_cells * bf_lev[i][oDir];
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
                    Box outflowBox = BoxLib::adjCell(crse_domain, outFace, 1);
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

void
NavierStokes::calcDpdt ()
{
    BL_ASSERT(state[Press_Type].descriptor()->timeType() == StateDescriptor::Point);

    MultiFab&  new_press   = get_new_data(Press_Type);
    MultiFab&  old_press   = get_old_data(Press_Type);
    MultiFab&  dpdt        = get_new_data(Dpdt_Type);
    const Real dt_for_dpdt = state[Press_Type].curTime()-state[Press_Type].prevTime();

    if (dt_for_dpdt != 0.0) 
    {
        for (MFIter mfi(dpdt); mfi.isValid(); ++mfi)
        {
            dpdt[mfi].copy(new_press[mfi],mfi.validbox(),0,mfi.validbox(),0,1);
            dpdt[mfi].minus(old_press[mfi],mfi.validbox(),0,0,1);
            dpdt[mfi].divide(dt_for_dpdt,mfi.validbox(),0,1);
        }
    }
    else
    {
        dpdt.setVal(0);
    }
}

void
NavierStokes::create_umac_grown ()
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::create_umac_grown()");

    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
        u_macG[n].copy(u_mac[n]);
        u_macG[n].FillBoundary(0,1);
        geom.FillPeriodicBoundary(u_macG[n],0,1);
    }
        
    if (level > 0)
    {
        BoxArray f_bnd_ba = GetBndryCells(grids,1,geom);

        BoxArray c_bnd_ba = BoxArray(f_bnd_ba.size());

        for (int i = 0; i < f_bnd_ba.size(); ++i)
        {
            c_bnd_ba.set(i,Box(f_bnd_ba[i]).coarsen(crse_ratio));
            f_bnd_ba.set(i,Box(c_bnd_ba[i]).refine(crse_ratio));
        }
            
        const BoxArray& cgrids = getLevel(level-1).boxArray();
            
        for (int n = 0; n < BL_SPACEDIM; ++n)
        {
            MultiFab crseT(BoxArray(cgrids).surroundingNodes(n),1,0);
                
            crseT.setVal(1.e200);
            for (MFIter mfi(crseT); mfi.isValid(); ++mfi)
                crseT[mfi].copy(getLevel(level-1).u_mac[n][mfi]);
            crseT.FillBoundary(0,1);
            getLevel(level-1).geom.FillPeriodicBoundary(crseT,0,1);
            //
            // crse_src & fine_src must have same parallel distribution.
            // We'll use the KnapSack distribution for the fine_src_ba.
            // Since fine_src_ba should contain more points, this'll lead
            // to a better distribution.
            //
            BoxArray crse_src_ba(c_bnd_ba);
            BoxArray fine_src_ba(f_bnd_ba);

            crse_src_ba.surroundingNodes(n);
            fine_src_ba.surroundingNodes(n);

            std::vector<long> wgts(fine_src_ba.size());

            for (unsigned int i = 0; i < wgts.size(); i++)
            {
                wgts[i] = fine_src_ba[i].numPts();
            }
            DistributionMapping dm;
            //
            // This call doesn't invoke the MinimizeCommCosts() stuff.
            // There's very little to gain with these types of coverings
            // of trying to use SFC or anything else.
            // This also guarantees that these DMs won't be put into the
            // cache, as it's not representative of that used for more
            // usual MultiFabs.
            //
            dm.KnapSackProcessorMap(wgts,ParallelDescriptor::NProcs());

            MultiFab crse_src; crse_src.define(crse_src_ba, 1, 0, dm, Fab_allocate);
	    MultiFab fine_src; fine_src.define(fine_src_ba, 1, 0, dm, Fab_allocate);

            crse_src.setVal(1.e200);
            crse_src.copy(crseT);
            crse_src.FillBoundary(0,1);
            getLevel(level-1).geom.FillPeriodicBoundary(crse_src,0,1);

            for (MFIter mfi(crse_src); mfi.isValid(); ++mfi)
            {
                const int  nComp = 1;
                const Box  box   = crse_src[mfi].box();
                const int* rat   = crse_ratio.getVect();
                FORT_PC_CF_EDGE_INTERP(box.loVect(), box.hiVect(), &nComp, rat, &n,
                                       crse_src[mfi].dataPtr(),
                                       ARLIM(crse_src[mfi].loVect()),
                                       ARLIM(crse_src[mfi].hiVect()),
                                       fine_src[mfi].dataPtr(),
                                       ARLIM(fine_src[mfi].loVect()),
                                       ARLIM(fine_src[mfi].hiVect()));
            }
            //
            // Replace pc-interpd fine data with preferred u_mac data at
            // this level u_mac valid only on surrounding faces of valid
            // region - this op will not fill grow region.
            //
            fine_src.copy(u_mac[n]);

            for (MFIter mfi(fine_src); mfi.isValid(); ++mfi)
            {
                //
                // Interpolate unfilled grow cells using best data from
                // surrounding faces of valid region, and pc-interpd data
                // on fine edges overlaying coarse edges.
                //
                const int  nComp = 1;
                const Box& fbox  = fine_src[mfi.index()].box();
                const int* rat   = crse_ratio.getVect();
                FORT_EDGE_INTERP(fbox.loVect(), fbox.hiVect(), &nComp, rat, &n,
                                 fine_src[mfi].dataPtr(),
                                 ARLIM(fine_src[mfi].loVect()),
                                 ARLIM(fine_src[mfi].hiVect()));
            }

            u_macG[n].copy(fine_src);
            u_macG[n].copy(u_mac[n]);
            u_macG[n].FillBoundary(0,1);
            geom.FillPeriodicBoundary(u_macG[n],0,1);
        }
    }
}

