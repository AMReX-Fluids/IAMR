//
// $Id: NavierStokes.cpp,v 1.91 1998-09-29 19:56:01 lijewski Exp $
//
// "Divu_Type" means S, where divergence U = S
// "Dsdt_Type" means pd S/pd t, where S is as above
//

#include <Misc.H>
#include <CoordSys.H>
#include <Geometry.H>
#include <BoxDomain.H>
#include <RunStats.H>
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
#include <PROJECTION_F.H>
#include <PROB_F.H>

#ifdef BL_USE_NEW_HFILES
#include <vector>
#include <cstdio>
using std::vector;
using std::streampos;
#else
#include <vector.h>
#include <stdio.h>
#endif

#define GEOM_GROW 1
#define HYP_GROW 3
#define PRESS_GROW 1
#define DIVU_GROW 1
#define DSDT_GROW 1
#define bogus_value 1.e20

const char NL = '\n';

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

//
// Initialization functions follow
//

//
// Initialization of static members to default values.
//

//
// ----------------------- static objects.
//
ErrorList   NavierStokes::err_list;
BCRec       NavierStokes::phys_bc;
Projection *NavierStokes::projector     = 0;
MacProj    *NavierStokes::mac_projector = 0;
Godunov    *NavierStokes::godunov       = 0;

//
// ----------------------- internal parameters.
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
int  NavierStokes::NUM_SCALARS  = 0;
int  NavierStokes::NUM_STATE    = 0;
Array<int> NavierStokes::is_conservative;

//
// ----------------------- viscosity parameters.
//
Real NavierStokes::be_cn_theta  = 0.5;
Real NavierStokes::visc_tol     = 1.0e-10;  // tolerance for viscous solve
Real NavierStokes::visc_abs_tol = 1.0e-10;  // absolute tol. for visc solve
Array<int> NavierStokes::is_diffusive;
Array<Real> NavierStokes::visc_coef;

//
// ----------------------- internal switches.
//
int  NavierStokes::do_temp         = 0;
int  NavierStokes::Temp            = -1;
int  NavierStokes::do_sync_proj    = 1;
int  NavierStokes::do_MLsync_proj  = 1;
int  NavierStokes::do_reflux       = 1;
int  NavierStokes::do_mac_proj     = 1;

//     
// ------------------------ new members for non-zero divu.
//
int  NavierStokes::additional_state_types_initialized = 0;
int  NavierStokes::Divu_Type = -1;
int  NavierStokes::Dsdt_Type = -1;
int  NavierStokes::have_divu = 0;
int  NavierStokes::have_dsdt = 0;
int  NavierStokes::S_in_vel_diffusion = 1;


Real NavierStokes::divu_minus_s_factor = 0.0;
Real NavierStokes::divu_relax_factor   = 0.0;
     
int  NavierStokes::num_state_type = 2;     // for backward compatibility


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
    if ((CoordSys::CoordType) coord == CoordSys::RZ &&
            phys_bc.lo(0) != Symmetry)
    {
        phys_bc.setLo(0,Symmetry);
        if (ParallelDescriptor::IOProcessor())
        {
            cout << "\nWarning:  Setting phys_bc at xlo to Symmetry\n\n";
        }
    }
#endif
}

void
NavierStokes::read_params ()
{
    //
    // Read parameters from input file and command line.
    //
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
                    cerr << "NavierStokes::variableSetUp:periodic in direction";
                    cerr << dir << " but low BC is not Interior\n";
                    BoxLib::Abort("NavierStokes::read_params()");
                }
                if (hi_bc[dir] != Interior)
                {
                    cerr << "NavierStokes::variableSetUp:periodic in direction";
                    cerr << dir << " but high BC is not Interior\n";
                    BoxLib::Abort("NavierStokes::read_params()");
                }
            }
        }
    }
    else
    {
        //
        // Do idiot check.  If not periodic, should be no interior.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (lo_bc[dir] == Interior)
            {
                cerr << "NavierStokes::variableSetUp:interior bc in direction";
                cerr << dir << " but no periodic\n";
                BoxLib::Abort("NavierStokes::read_params()");
            }
            if (hi_bc[dir] == Interior)
            {
                cerr << "NavierStokes::variableSetUp:interior bc in direction";
                cerr << dir << " but no periodic\n";
                BoxLib::Abort("NavierStokes::read_params()");
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
    pp.query("gravity",gravity);
    //
    // Get run options.
    //
    pp.query("do_temp",        do_temp        );
    int initial_do_sync_proj = do_sync_proj;
    pp.query("do_sync_proj",   do_sync_proj   );
    pp.query("do_MLsync_proj", do_MLsync_proj );
    pp.query("do_reflux",      do_reflux      );
    pp.query("do_mac_proj",    do_mac_proj    );
    //
    // This test insures if the user toggles do_sync_proj,
    // the user has knowledge that do_MLsync_proj is meaningless.
    //
    if (do_MLsync_proj && !do_sync_proj && initial_do_sync_proj != do_sync_proj)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            cout << "mismatched options for NavierStokes\n";
            cout << "do_MLsync_proj and do_sync_proj are inconsistent\n";
        }
        BoxLib::Abort("NavierStokes::read_params()");
    }
    //
    // Read viscous/diffusive parameters and array of viscous/diffusive coeffs.
    // NOTE: at this point, we dont know number of state variables
    //       so just read all values listed.
    //
    pp.query("visc_tol",visc_tol);
    pp.query("visc_abs_tol",visc_abs_tol);

    int n_vel_visc_coef   = pp.countval("vel_visc_coef");
    int n_temp_cond_coef  = pp.countval("temp_cond_coef");
    int n_scal_diff_coefs = pp.countval("scal_diff_coefs");

    if (n_vel_visc_coef != 1)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            cout << "NavierStokes::read_params: you must input only one vel_visc_coef"
                 << endl;
        }
        exit(0);
    }

    if (do_temp && n_temp_cond_coef != 1)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            cout << "NavierStokes::read_params: you must input only one temp_cond_coef"
                 << endl;
        }
        exit(0);
    }

    int n_visc = BL_SPACEDIM + 1 + n_scal_diff_coefs;
    if (do_temp)
        n_visc++;
    visc_coef.resize(n_visc);
    is_diffusive.resize(n_visc);
 
    pp.get("vel_visc_coef",visc_coef[0]);
    for (int i = 1; i < BL_SPACEDIM; i++)
      visc_coef[i] = visc_coef[0];

    // Here we set the coefficient for density, which does not diffuse.
    visc_coef[Density] = -1;

    // Set the coefficients for the scalars, but temperature
    Array<REAL> scal_diff_coefs(n_scal_diff_coefs);
    pp.getarr("scal_diff_coefs",scal_diff_coefs,0,n_scal_diff_coefs);

    int scalId = Density;
    for (int i = 0; i < n_scal_diff_coefs; i++) {
      visc_coef[++scalId] = scal_diff_coefs[i];
    }

    // Set the coefficient for temperature
    if (do_temp)
    {
	Temp = ++scalId;
	pp.get("temp_cond_coef",visc_coef[Temp]);
    }
    
    pp.query("divu_minus_s_factor",divu_minus_s_factor);
    pp.query("divu_relax_factor",divu_relax_factor);
    pp.query("S_in_vel_diffusion",S_in_vel_diffusion);
    pp.query("be_cn_theta",be_cn_theta);
    if (be_cn_theta > 1.0 || be_cn_theta < .5)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            cout << "NavierStokes::read_params: must have be_cn_theta <= 1.0 && >= .5"
                 << NL;
        }
        BoxLib::Abort("NavierStokes::read_params()");
    }
}

NavierStokes::NavierStokes ()
    :
    radius(PArrayManage)
{
    rho_avg      = 0;
    rho_half     = 0;
    p_avg        = 0;
    Vsync        = 0;
    Ssync        = 0;
    sync_reg     = 0;
    advflux_reg  = 0;
    viscflux_reg = 0;
    u_mac        = 0;
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
    AmrLevel(papa,lev,level_geom,bl,time),
    radius(PArrayManage)
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
        rho_avg = new MultiFab(grids,1,1,Fab_allocate);
        p_avg   = new MultiFab(P_grids,1,0,Fab_allocate);
    }
    rho_half = new MultiFab(grids,1,1,Fab_allocate);
    //
    // Build metric coeficients for RZ calculations.
    //
    buildMetrics();
    //
    // Build base state rho for atmospheric calculations
    // this is probably an Atmosphere virtual function.
    //
    buildRho0();
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
    Vsync = 0;
    Ssync = 0;
    u_mac = 0;
    aofs  = 0;
    //
    // Set up the level projector.
    //
    if (do_MLsync_proj || do_sync_proj)
    {
        if (projector == 0)
        {
            projector = new Projection(parent,&phys_bc,
                                       do_sync_proj,parent->finestLevel(),
                                       radius_grow);
        }
        projector->install_level(level, this, &radius );
    }
    //
    // Set up the godunov box.
    //
    SetGodunov();
    //
    // Set up diffusion.
    //
    diffusion = new Diffusion(parent, this,
                              (level > 0) ? getLevel(level-1).diffusion : 0,
                              NUM_STATE, viscflux_reg, volume, area,
                              is_diffusive, visc_coef);
    //
    // Set up the mac projector.
    //
    if (mac_projector == 0)
    {
        mac_projector = new MacProj(parent, parent->finestLevel(),
                                    &phys_bc, radius_grow );
    }
    mac_projector->install_level(level, this, volume, area, &radius );
}

NavierStokes::~NavierStokes ()
{
    delete rho_avg;
    delete p_avg;
    delete rho_half;
    delete Vsync;
    delete Ssync;
    delete sync_reg;
    delete advflux_reg;
    delete viscflux_reg;
    delete [] u_mac;
    
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        area[dir].clearUnmanaged();
    }
    
    if (mac_projector != 0)
        mac_projector->cleanup(level);
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
    bool have_temp = isStateVariable("temp", dummy_State_Type, Temp);
    have_temp &= dummy_State_Type == State_Type;
    assert((do_temp && have_temp)  ||  (!do_temp && !have_temp));

    int _Divu = -1;
    int dummy_Divu_Type;
    have_divu = 0;
    have_divu = isStateVariable("divu", dummy_Divu_Type, _Divu);
    have_divu = have_divu && dummy_Divu_Type==Divu_Type;
    if (ParallelDescriptor::IOProcessor())
    {
      cout << "NavierStokes::init_additional_state_types()::have_divu = "
           << have_divu << NL;
    }
    if (have_divu && _Divu!=Divu)
    {
      if (ParallelDescriptor::IOProcessor())
      {
        cout << "NavierStokes::init_additional_state_types(): divu must be " <<
                "0-th, Divu_Type component in the state\n";
      }
      BoxLib::Abort("NavierStokes::init_additional_state_types()");
    }

    int _Dsdt = -1;
    int dummy_Dsdt_Type;
    have_dsdt = 0;
    have_dsdt = isStateVariable("dsdt", dummy_Dsdt_Type, _Dsdt);
    have_dsdt = have_dsdt && dummy_Dsdt_Type==Dsdt_Type;
    if (ParallelDescriptor::IOProcessor())
    {
      cout << "NavierStokes::init_additional_state_types()::have_dsdt = "
           << have_dsdt << NL;
    }
    if (have_dsdt && _Dsdt!=Dsdt)
    {
      if (ParallelDescriptor::IOProcessor())
      {
        cout << "NavierStokes::init_additional_state_types(): dsdt must be " <<
                "0-th, Dsdt_Type component in the state\n";
      }
      BoxLib::Abort("NavierStokes::init_additional_state_types()");
    }
    if (have_dsdt && !have_divu)
    {
      if (ParallelDescriptor::IOProcessor())
      {
        cout << "NavierStokes::init_additional_state_types(): "
             << "must have divu in order to have dsdt\n";
      }
      BoxLib::Abort("NavierStokes::init_additional_state_types()");
    }

    num_state_type = desc_lst.length();
    if (ParallelDescriptor::IOProcessor())
    {
      cout << "NavierStokes::init_additional_state_types: num_state_type = "
           << num_state_type << NL;
    }
}

void
NavierStokes::SaveOldBoundary (Real time)
{
    MultiFab& Sold = get_old_data(State_Type);

    FillPatchIterator Sold_fpi(*this,Sold,1,time,State_Type,0,NUM_STATE);

    for ( ; Sold_fpi.isValid(); ++Sold_fpi)
    {
        Sold[Sold_fpi.index()].copy(Sold_fpi(), 0, 0, NUM_STATE);
    }
}

void
NavierStokes::SaveNewBoundary (Real time)
{
    MultiFab& Snew = get_new_data(State_Type);

    FillPatchIterator Snew_fpi(*this,Snew,1,time,State_Type,0,NUM_STATE);

    for ( ; Snew_fpi.isValid(); ++Snew_fpi)
    {
        Snew[Snew_fpi.index()].copy(Snew_fpi(), 0, 0, NUM_STATE);
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

    for (MultiFabIterator mfi(P_new); mfi.isValid(); ++mfi)
    {
        DependentMultiFabIterator dmfi(mfi,P_old);
        dmfi().copy(mfi());
    }
}

void
NavierStokes::zeroNewPress ()
{
    get_new_data(Press_Type).setVal(0.0);
}

void
NavierStokes::zeroOldPress ()
{
    get_old_data(Press_Type).setVal(0.0);
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
    {
        initOldPress();
    }
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
    {
        godunov = new Godunov();
    }
}

void
NavierStokes::restart (Amr&     papa,
                       istream& is)
{
    AmrLevel::restart(papa,is);

    if (do_MLsync_proj || do_sync_proj)
    {
        if (projector == 0)
        {
            projector = new Projection(parent,&phys_bc,do_sync_proj,
                                       parent->finestLevel(),radius_grow);
        }
        projector->install_level(level, this, &radius );
    }
    //
    // Set the godunov box.
    //
    SetGodunov();
    
    if (mac_projector == 0)
    {
        mac_projector = new MacProj(parent, parent->finestLevel(),
                                    &phys_bc, radius_grow );
    }
    mac_projector->install_level(level, this, volume, area, &radius );

    rho_avg = 0;
    p_avg   = 0;
    const BoxArray& P_grids = state[Press_Type].boxArray();
    //
    // Alloc space for density and temporary pressure variables.
    //
    if (level > 0)
    {
        rho_avg = new MultiFab(grids,1,1,Fab_allocate);
        p_avg   = new MultiFab(P_grids,1,0,Fab_allocate);
    }
    rho_half = new MultiFab(grids,1,1,Fab_allocate);
    //
    // Build metric coeficients for RZ calculations.
    //
    buildMetrics();
    //
    // Build base state rho for atmospheric calculations.
    // This is probably an Atmosphere virtual function.
    //
    buildRho0();

    assert(sync_reg == 0);
    if (level > 0 && do_sync_proj)
    {
        sync_reg = new SyncRegister(grids,crse_ratio,level);
    }
    assert(advflux_reg == 0);
    if (level > 0 && do_reflux)
    {
        advflux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);
    }
    assert(viscflux_reg == 0);
    if (level > 0 && do_reflux)
    {
        viscflux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);
    }

    assert(Vsync == 0);
    assert(Ssync == 0);
    if (level < parent->finestLevel())
    {
        Vsync = new MultiFab(grids,BL_SPACEDIM,1,Fab_allocate);
        Ssync = new MultiFab(grids,NUM_STATE-BL_SPACEDIM,1,Fab_allocate);
    }

    diffusion = new Diffusion(parent, this,
                              (level > 0) ? getLevel(level-1).diffusion : 0,
                              NUM_STATE, viscflux_reg, volume, area,
                              is_diffusive, visc_coef);

}

//
// Build rho0 arrays as ones for the base class.  Empty function for now.
//

void NavierStokes::buildRho0 () {}

void
NavierStokes::buildMetrics ()
{
    radius.resize(grids.length());

    const Real dxr = geom.CellSize()[0];

    for (int i = 0; i < grids.length(); i++)
    {
        const int ilo = grids[i].smallEnd(0)-radius_grow;
        const int ihi = grids[i].bigEnd(0)+radius_grow;
        const int len = ihi - ilo + 1;
        Real* rad = new Real[len];
        radius.set(i,rad);
        if (CoordSys::IsCartesian())
        {
            for (int j = 0; j < len; j++)
            {
                rad[j] = 1.0;
            }
        }
        else
        {
            const Real xlo = grid_loc[i].lo(0) + (0.5 - radius_grow)*dxr;
            for (int j = 0; j < len; j++)
            {
                rad[j] = xlo + j*dxr;
            }
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
    // Reset state and pressure types.
    //
    state[State_Type].reset();
    state[State_Type].setTimeLevel(time,dt_old,dt_new);
    initOldPress();
    state[Press_Type].setTimeLevel(time-dt_old,dt_old,dt_old);
    //
    // Reset state types for divu not equal to zero.
    //
    if (have_divu)
    {
        state[Divu_Type].reset();
        state[Divu_Type].setTimeLevel(time,dt_old,dt_new);
        if (have_dsdt)
        {
            state[Dsdt_Type].reset();
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
    state[Press_Type].setTimeLevel(time-dt_old,dt_old,dt_old);
}

//
// This function initializes the State and Pressure with data.
//

void
NavierStokes::initData ()
{
    static const aString RunstatString("init_data");

    RunStats rs(RunstatString, Level());

    rs.start();
    //
    // Initialize the state and the pressure.
    //
    int ns          = NUM_STATE - BL_SPACEDIM;
    const Real* dx  = geom.CellSize();
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& P_new = get_new_data(Press_Type);
    Real cur_time   = state[State_Type].curTime();

    for (MultiFabIterator snewmfi(S_new); snewmfi.isValid(); ++snewmfi)
    {
        DependentMultiFabIterator pnewmfi(snewmfi, P_new);

        assert(grids[snewmfi.index()] == snewmfi.validbox());

        const int* lo   = snewmfi.validbox().loVect();
        const int* hi   = snewmfi.validbox().hiVect();
        const int* s_lo = snewmfi().loVect();
        const int* s_hi = snewmfi().hiVect();
        const int* p_lo = pnewmfi().loVect();
        const int* p_hi = pnewmfi().hiVect();
        pnewmfi().setVal(0.0);
        int i = snewmfi.index();
        FORT_INITDATA (&level,&cur_time,lo,hi,&ns,
                       snewmfi().dataPtr(Xvel),
                       snewmfi().dataPtr(BL_SPACEDIM),
                       ARLIM(s_lo), ARLIM(s_hi),
                       pnewmfi().dataPtr(),
                       ARLIM(p_lo), ARLIM(p_hi),
                       dx,grid_loc[i].lo(),grid_loc[i].hi() );
    }
    //
    // Initialize other types.
    //
    initDataOtherTypes();
    //
    // Initialize divU and dSdt.
    //
    if (have_divu)
    {
        const Real cur_time = state[Divu_Type].curTime();
        MultiFab& Divu_new  = get_new_data(Divu_Type);
        const Real dt       = 1.0;
        state[State_Type].setTimeLevel(cur_time,dt,dt);
        const Real dtin     = -1.0; // dummy value denotes initialization
        calc_divu(cur_time, dtin, Divu_new);
        if (have_dsdt)
        {
            get_new_data(Dsdt_Type).setVal(0);
        }
    }

    rs.end();
}

//
// Fills a new level n with best level n and coarser data available.
//

void
NavierStokes::init (AmrLevel &old)
{
    //
    // Get data pointers.
    //
    NavierStokes* oldns = (NavierStokes*) &old;
    MultiFab &S_new     = get_new_data(State_Type);
    MultiFab &P_new     = get_new_data(Press_Type);
    MultiFab &P_old     = get_old_data(Press_Type);
    //
    // Get time information.
    //
    const Real dt_new    = parent->dtLevel(level);
    const Real cur_time  = oldns->state[State_Type].curTime();
    const Real prev_time = oldns->state[State_Type].prevTime();
    const Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);
    const Real cur_pres_time = state[Press_Type].curTime();
    //
    // Get best state and pressure data.
    //
    for (FillPatchIterator snewfpi(old,S_new,0,cur_time,State_Type,0,NUM_STATE);
        snewfpi.isValid();
         ++snewfpi)
    {
        S_new[snewfpi.index()].copy(snewfpi());
    }

    for (FillPatchIterator pnewfpi(old,P_new,0,cur_pres_time,Press_Type,0,1);
         pnewfpi.isValid();
         ++pnewfpi)
    {
        P_new[pnewfpi.index()].copy(pnewfpi());
        P_old[pnewfpi.index()].copy(pnewfpi());
    }
    //
    // Get best divu and dSdt data.
    //
    if (have_divu)
    {
        MultiFab& Divu_new = get_new_data(Divu_Type);
        
        for (FillPatchIterator divunewfpi(old,Divu_new,0,cur_time,Divu_Type,0,1);
             divunewfpi.isValid();
             ++divunewfpi)
        {
            Divu_new[divunewfpi.index()].copy(divunewfpi());
        }

        if (have_dsdt)
        {
            MultiFab& Dsdt_new = get_new_data(Dsdt_Type);

            for (FillPatchIterator dsdtnewfpi(old,Dsdt_new,0,cur_time,
                                              Dsdt_Type,0,1);
                 dsdtnewfpi.isValid();
                 ++dsdtnewfpi)
            {
                Dsdt_new[dsdtnewfpi.index()].copy(dsdtnewfpi());
            }
        }
    }
}

//
// Inits the data on a new level that did not exist before regridding.
//

void
NavierStokes::init ()
{
    //
    // Get data pointers.
    //
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& P_new = get_new_data(Press_Type);
    MultiFab& P_old = get_old_data(Press_Type);
    //
    // Get time information.
    //
    assert(level > 0);

    const Array<Real>& dt_amr = parent->dtLevel();
    Array<Real> dt_new(level+1);
    for (int lev = 0; lev < level; lev++)
    {
        dt_new[lev] = dt_amr[lev];
    }
    const Real dt = dt_new[level-1]/(Real)parent->MaxRefRatio(level-1);
    dt_new[level] = dt;
    parent->setDtLevel(dt_new);

    NavierStokes& old    = getLevel(level-1);
    const Real cur_time  = old.state[State_Type].curTime();
    const Real prev_time = old.state[State_Type].prevTime();
    const Real dt_old    = (cur_time-prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);

    Real cur_pres_time = state[Press_Type].curTime();
    //
    // Get best coarse state and pressure data.
    //
    FillCoarsePatch(S_new,0,cur_time,State_Type,0,NUM_STATE);
    FillCoarsePatch(P_new,0,cur_pres_time,Press_Type,0,1);

    for (MultiFabIterator mfi(P_new); mfi.isValid(); ++mfi)
    {
        DependentMultiFabIterator dmfi(mfi,P_old);
        dmfi().copy(mfi());
    }
    //
    // Get best coarse divU and dSdt data.
    //
    if (have_divu)
    {
        FillCoarsePatch(get_new_data(Divu_Type),0,cur_time,Divu_Type,0,1);
        if (have_dsdt)
        {
            FillCoarsePatch(get_new_data(Dsdt_Type),0,cur_time,Dsdt_Type,0,1);
        }
    }
}

//
// ADVANCE FUNCTIONS
//

//
// This function ensures that the multifab registerss and boundary
// flux registers needed for syncing the composite grid
//
//     u_mac, Vsync, Ssync, rhoavg, fr_adv, fr_visc
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
    const int finest_level = parent->finestLevel();

    mac_projector->setup(level);
    //
    // Why are they defined here versus the constructor?
    //
    if (level < finest_level)
    {
        if (Vsync == 0)
        {
            Vsync = new MultiFab(grids,BL_SPACEDIM,1,Fab_allocate);
        }
        if (Ssync == 0)
        {
            Ssync = new MultiFab(grids,NUM_STATE-BL_SPACEDIM,1,Fab_allocate);
        }
        Vsync->setVal(0.0);
        Ssync->setVal(0.0);
    }
    //
    // Set reflux registers to zero.
    //
    if (do_reflux && level < finest_level)
    {
        getAdvFluxReg(level+1).setVal(0.0);
        getViscFluxReg(level+1).setVal(0.0);
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
    //
    // Alloc multifab to hold advective tendencies.
    //
    assert(aofs == 0);
    aofs = new MultiFab(grids,NUM_STATE,0,Fab_allocate);
    //
    // Set rho_avg.
    //
    if (!initial_step && (level > 0) && (iteration == 1))
    {
        const Real fratio = (Real) ncycle;
        const Real alpha = 0.5/fratio;
        initRhoAvg(alpha);
    }
    //
    // Set up state multifabs for the advance.
    //
    for (int k = 0; k < num_state_type; k++)
    {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    get_new_data(State_Type).setVal(bogus_value);

    if (level > 0 || geom.isAnyPeriodic())
    {
        //
        // This is neccessary so that diffusion works properly during the first
        // time step (Diffusion can call only AmrLevel::setPhysBndryValues() to
        // fill ghost cells.  That function does not know about ghost cells
        // over coarse values, or periodic bc's)
        //
        for (int k = 0; k < num_state_type; k++)
        {
            if (k != Press_Type)
            {
                MultiFab& new_state = get_new_data(k);
                MultiFab& old_state = get_old_data(k);
                FillStateBndry(time,k,0,old_state.nComp());
                FillStateBndry(time+dt,k,0,new_state.nComp());
            }
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
        //
        // Delete mac edge velocities.
        //
        delete [] u_mac;
        u_mac = 0;
    }
    //
    // Delete advective tendencies.
    //
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
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "Advancing grids at level " << level
             << " : starting time = " << time << " with dt = " << dt << NL;
    }
    advance_setup(time,dt,iteration,ncycle);

//     cout << "!!!!!!!!!!! we are here, reading data" << NL;
//     read_all( this );
//     FILE *dfile = fopen( "dt.dat", "rb" );
//     int nitems = fread( &dt, sizeof(dt), 1, dfile );
//     fclose(dfile);
//     cout << "!!!!!!!!!!! done reading data" << NL;
    
    SaveOldBoundary(time); 
    MultiFab& Snew  = get_new_data(State_Type);
    MultiFab& Sold  = get_old_data(State_Type);
    MultiFab& P_new = get_new_data(Press_Type);
    MultiFab& P_old = get_old_data(Press_Type);

    static const aString VelPredictStr("vel_predict");
    static const aString VelAdvectStr("vel_advect");
    static const aString ScalAdvectStr("scal_advect");
    static const aString VelUpdateStr("vel_update");
    static const aString ScalUpdateStr("scal_update");
    static const aString MacProjectStr("mac_project");

    RunStats vel_pred_stats(VelPredictStr, level);
    RunStats  vel_adv_stats(VelAdvectStr , level);
    RunStats scal_adv_stats(ScalAdvectStr, level);
    RunStats  vel_upd_stats(VelUpdateStr , level);
    RunStats scal_upd_stats(ScalUpdateStr, level);
    RunStats      mac_stats(MacProjectStr, level);
    //
    // ------------------ Advance starts here
    //
    // Compute traced states for normal comp of velocity at half time level.
    //
    vel_pred_stats.start();
    Real dt_test = 0.0, dummy = 0.0;
    dt_test = predict_velocity(dt,dummy);
    vel_pred_stats.end();
    //
    // Do MAC projection and update edge velocities.
    //
    mac_stats.start();
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "... mac_projection\n";
    }
    MultiFab* divu = getDivCond(0,time);
    MultiFab* dsdt = getDsdt(0,time);

    for (MultiFabIterator mfi(*divu); mfi.isValid(); ++mfi)
    {
        DependentMultiFabIterator dmfi(mfi,*dsdt);
        dmfi().mult(.5*dt);
        mfi().plus(dmfi());
    }
    delete dsdt;
    //
    // Compute mac velocities and maximum cfl number.
    //
    if (do_mac_proj)
        mac_projector->mac_project(level,u_mac,Sold,dt,time,*divu,have_divu);
    mac_stats.end();

    delete divu;
    //
    // Advect velocities.
    //
    vel_adv_stats.start();
    velocity_advection(dt);
    vel_adv_stats.end();
    //
    // Advect scalars.
    //
    const int first_scalar = Density;
    const int last_scalar  = first_scalar + NUM_SCALARS - 1;
    scal_adv_stats.start();
    scalar_advection(dt,first_scalar,last_scalar);
    scal_adv_stats.end();
    //
    // Add the advective and other terms to get scalars at t^{n+1}.
    //
    scal_upd_stats.start();
    scalar_update(dt,first_scalar,first_scalar);
    //
    // Compute rho at half time level including 1-zone bndry values.
    //
    makerhonph(dt);
    //
    // Add the advective and other terms to get scalars at t^{n+1}.
    //
    scalar_update(dt,first_scalar+1,last_scalar);
    scal_upd_stats.end();
    //
    // S appears in rhs of the velocity update, so we better do it now.
    //
    if (have_divu)
    {
        calc_divu(time+dt, dt, get_new_data(Divu_Type));
        if (have_dsdt)
        {
            calc_dsdt(time,dt,get_new_data(Dsdt_Type));
        }
    }
    //
    // Add the advective and other terms to get velocity at t^{n+1}.
    //
    vel_upd_stats.start();
    velocity_update(dt);
    vel_upd_stats.end();
    //
    // Clean up after the predicted value at t^n+1.
    // Estimate new timestep from umac cfl.
    //
    advance_cleanup(dt,iteration,ncycle);
    //
    // Increment rho average.
    //
    if (!initial_step && level > 0)
    {
        const Real fratio = (Real) ncycle;
        Real alpha = 1.0/fratio;
        if (iteration == ncycle)
            alpha = 0.5/fratio;
        incrRhoAvg(alpha);
    }
    //
    // Do a level project to update the pressure and velocity fields.
    //
    if (!initial_step)
    {
        if (projector)
            level_projector(dt,time,iteration);
        if (level > 0)
        {
            const Real alpha = 1.0/ (Real) ncycle;
            incrPAvg(iteration,alpha);
        }
    }
    //
    // Relax back to continuity constraint.
    //
    if (divu_relax_factor > 0.0 && !initial_step)
    {
        MultiFab* delta_U = new MultiFab (grids,BL_SPACEDIM,0,Fab_allocate);
        compute_grad_divu_minus_s(time+dt, delta_U, 0);
        Snew.plus(*delta_U, 0, BL_SPACEDIM, 0);
        delete delta_U;
    }

    return dt_test;  // Return estimate of best new timestep.
}

void
NavierStokes::level_projector (Real dt,
                               Real time,
                               int  iteration)
{
   if (iteration > 0)
   {
       static const aString RunstatString("level_project");

       RunStats lp_stats(RunstatString,level);

       lp_stats.start();

       MultiFab& U_old = get_old_data(State_Type);
       MultiFab& U_new = get_new_data(State_Type);
       MultiFab& P_old = get_old_data(Press_Type);
       MultiFab& P_new = get_new_data(Press_Type);

       SyncRegister* crse_ptr = 0;

       if (level < parent->finestLevel() && do_sync_proj)
       {
           crse_ptr = &(getLevel(level+1).getSyncReg());
       }

       Array<int*> sync_bc(grids.length());
       Array< Array<int> > sync_bc_array(grids.length());

       for (int i = 0; i < grids.length(); i++)
       {
           sync_bc_array[i] = getBCArray(State_Type,i,Xvel,BL_SPACEDIM);
           sync_bc[i]       = sync_bc_array[i].dataPtr();
       }

       MultiFab *dsdt = 0, *divuold = 0;

       if (have_divu)
       {
           FillStateBndry(time,Divu_Type,0,1);
           FillStateBndry(time+dt,Divu_Type,0,1);

           dsdt    = getDivCond(1,time+dt);
           divuold = getDivCond(1,time);

           dsdt->minus(*divuold,0,1,1);

           for (MultiFabIterator mfi(*dsdt); mfi.isValid(); ++mfi)
           {
               mfi().mult(1.0/dt,0,1);
           }
       }
       else
       {
           dsdt    = new MultiFab(grids,1,1,Fab_allocate);
           divuold = new MultiFab(grids,1,1,Fab_noallocate);
           dsdt->setVal(0.0);
       }

       int crse_dt_ratio  = (level > 0) ? parent->MaxRefRatio(level-1) : -1;
       const Real cur_pres_time = state[Press_Type].curTime();

       projector->level_project(level,dt,cur_pres_time,time,time+dt,geom,
                                U_old,U_new,P_old,P_new,rho_half,*dsdt,
                                crse_ptr,sync_reg,crse_dt_ratio,
                                sync_bc.dataPtr(),iteration,
                                divu_minus_s_factor,*divuold,have_divu);
       delete dsdt;
       delete divuold;

       lp_stats.end();
   }
   else
   {
       const Real cur_pres_time = state[Press_Type].curTime();
       MultiFab& P_old = get_old_data(Press_Type);
       BoxLib::Abort("NavierStokes::level_projector calling harmonic_project");
       projector->harmonic_project(level,dt,cur_pres_time,geom,P_old);
   }
}

void
NavierStokes::makerhonph (Real dt)
{
    const Real half_time = state[State_Type].prevTime() + 0.5*dt;

    FillPatchIterator rho_fpi(*this,*rho_half,1,half_time,State_Type,Density,1);

    for ( ; rho_fpi.isValid(); ++rho_fpi)
    {
        DependentMultiFabIterator dmfi(rho_fpi,*rho_half);
        dmfi().copy(rho_fpi());
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
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "... predict edge velocities\n";
    }
    //
    // Get simulation parameters.
    //
    const Real* dx            = geom.CellSize();
    const Real prev_time      = state[State_Type].prevTime();
    const Real prev_pres_time = state[Press_Type].prevTime();
    //
    // Compute viscous terms at level n.
    //
    MultiFab visc_terms(grids,BL_SPACEDIM,1,Fab_allocate);
    getViscTerms(visc_terms,Xvel,BL_SPACEDIM,prev_time);
    if (be_cn_theta == 1.0)
        visc_terms.setVal(0.0,1);
    //
    // Set up the timestep estimation.
    //
    Real cflgrid,u_max[3];
    Real cflmax = 1.0e-10;
    comp_cfl    = (level == 0) ? cflmax : comp_cfl;

    FArrayBox tforces, Gp;
    Array<int> bndry[BL_SPACEDIM];
    //
    // FillPatch'd state data.
    //
    FillPatchIterator P_fpi(*this,get_old_data(Press_Type),1,
                            prev_pres_time,Press_Type,0,1);

    FillPatchIterator U_fpi(*this,visc_terms,HYP_GROW,prev_time,
                            State_Type,Xvel,BL_SPACEDIM);

    FillPatchIterator Rho_fpi(*this,visc_terms,1,prev_time,
                              State_Type,Density,1);

    for ( ;
          U_fpi.isValid() && Rho_fpi.isValid() && P_fpi.isValid();
          ++U_fpi, ++Rho_fpi, ++P_fpi)
    {
        //
        // Since all the MultiFabs are based on same grid we'll use indices.
        //
        const int i = U_fpi.index();

        setForce(tforces,i,1,Xvel,BL_SPACEDIM,Rho_fpi());
        //
        // Test velocities, rho and cfl.
        //
        cflgrid  = godunov->test_u_rho(U_fpi(),Rho_fpi(),grids[i],dx,dt,u_max);
        cflmax   = Max(cflgrid,cflmax);
        comp_cfl = Max(cflgrid,comp_cfl);
        //
        // Compute the total forcing.
        //
        getGradP(P_fpi(),Gp,grids[i],1);

        godunov->Sum_tf_gp_visc(tforces,visc_terms[i],Gp,Rho_fpi());
        //
        // Set up the Godunov box.
        //
        D_TERM(bndry[0] = getBCArray(State_Type,i,0,1);,
               bndry[1] = getBCArray(State_Type,i,1,1);,
               bndry[2] = getBCArray(State_Type,i,2,1);)

        godunov->Setup(grids[i], dx, dt, 1,
                       u_mac[0][i], bndry[0].dataPtr(),
                       u_mac[1][i], bndry[1].dataPtr(),
#if (BL_SPACEDIM == 3)                         
                       u_mac[2][i], bndry[2].dataPtr(),
#endif
                       U_fpi(), Rho_fpi(), tforces);
        //
        // Predict the mac velocities.
        //
        godunov->ComputeUmac(grids[i], dx, dt, 
                             u_mac[0][i], bndry[0].dataPtr(),
                             u_mac[1][i], bndry[1].dataPtr(),
#if (BL_SPACEDIM == 3)
                             u_mac[2][i], bndry[2].dataPtr(),
#endif
                             U_fpi(), tforces);
    }

    if (level == 0 && geom.isAnyPeriodic())
    {
        test_umac_periodic();
    }
    //
    // Compute estimate of the timestep.
    //
    Real tempdt = Min(change_max,cfl/cflmax);

    ParallelDescriptor::ReduceRealMin(tempdt);

    return dt*tempdt;
}

//
// Test for periodic umac on a single level 0 grid.
//

void
NavierStokes::test_umac_periodic ()
{
    //
    // Error block.
    //
    if (grids.length() != 1)
        return;
    assert(level == 0);
    //
    // Get the bounds and grid size.
    //
    const Box &grd = grids[0];
    const int *lo  = grd.loVect();
    const int *hi  = grd.hiVect();
    //
    // Get the velocities.
    //
    int xperiod        = (geom.isPeriodic(0) ? 1 : 0);
    const int *ux_lo   = u_mac[0][0].loVect();
    const int *ux_hi   = u_mac[0][0].hiVect();
    const Real *ux_dat = u_mac[0][0].dataPtr();

    int yperiod        = (geom.isPeriodic(1) ? 1 : 0);
    const int *uy_lo   = u_mac[1][0].loVect();
    const int *uy_hi   = u_mac[1][0].hiVect();
    const Real *uy_dat = u_mac[1][0].dataPtr();
#if (BL_SPACEDIM == 3)
    int zperiod        = (geom.isPeriodic(2) ? 1 : 0);
    const int *uz_lo   = u_mac[2][0].loVect();
    const int *uz_hi   = u_mac[2][0].hiVect();
    const Real *uz_dat = u_mac[2][0].dataPtr();
#endif
    //
    // Call the fortran.
    //
    Real udiff = 0.0;
    Real vdiff = 0.0;
#if (BL_SPACEDIM == 3)                                             
    Real wdiff = 0.0;
#endif
    
    FORT_TEST_UMAC_PERIODIC(lo, hi, 
                            ux_dat, ARLIM(ux_lo), ARLIM(ux_hi),
                            &xperiod, &udiff,
                            uy_dat, ARLIM(uy_lo), ARLIM(uy_hi),
                            &yperiod, &vdiff,
#if (BL_SPACEDIM == 3)                                             
                            uz_dat, ARLIM(uz_lo), ARLIM(uz_hi),
                            &zperiod, &wdiff,
#endif
                            &level);
    
    if (level == 0 && udiff > 1.0e-10)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            cout << "!!!!!!!!!!!!!!!!!!!!!!!! udiff = " << udiff << NL;
        }
        BoxLib::Abort("Exiting.");
    }
    if (level == 0 && vdiff > 1.0e-10)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            cout << "!!!!!!!!!!!!!!!!!!!!!!!! vdiff = " << vdiff << NL;
        }
        BoxLib::Abort("Exiting.");
    }
#if (BL_SPACEDIM == 3)                                             
    if (level == 0 && vdiff > 1.0e-10)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            cout << "!!!!!!!!!!!!!!!!!!!!!!!! vdiff = " << vdiff << NL;
        }
        BoxLib::Abort("Exiting.");
    }
#endif
}

//
// This routine advects the velocities
//

void
NavierStokes::velocity_advection (Real dt)
{
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "... advect velocities\n";
    }
    //
    // Get simulation parameters.
    //
    const int finest_level    = parent->finestLevel();
    const Real *dx            = geom.CellSize();
    const Real prev_time      = state[State_Type].prevTime();
    const Real prev_pres_time = state[Press_Type].prevTime();
    //
    // Compute viscosity components.
    //
    MultiFab visc_terms(grids,BL_SPACEDIM,1,Fab_allocate);
    getViscTerms(visc_terms,Xvel,BL_SPACEDIM,prev_time);
    if (be_cn_theta == 1.0)
        visc_terms.setVal(0.0,1);

    Array<int> bndry[BL_SPACEDIM];
    FArrayBox xflux, yflux, zflux, divu, tforces, Gp;
    //
    // FillPatch'd state data.
    //
    FillPatchIterator P_fpi(*this,get_old_data(Press_Type),1,
                            prev_pres_time,Press_Type,0,1);

    FillPatchIterator U_fpi(*this,visc_terms,HYP_GROW,prev_time,
                            State_Type,Xvel,BL_SPACEDIM);

    FillPatchIterator Rho_fpi(*this,visc_terms,1,prev_time,
                              State_Type,Density,1);
    //
    // Compute the advective forcing.
    //
    for ( ;
          U_fpi.isValid() && Rho_fpi.isValid() && P_fpi.isValid();
          ++U_fpi, ++Rho_fpi, ++P_fpi)
    {
        //
        // Since all the MultiFabs are on same grid we'll just use indices.
        //
        const int i = U_fpi.index();

        setForce(tforces,i,1,Xvel,BL_SPACEDIM,Rho_fpi());
        //
        // Compute the total forcing.
        //
        getGradP(P_fpi(),Gp,grids[i],1);

        godunov->Sum_tf_gp_visc(tforces,visc_terms[i],Gp,Rho_fpi());
        //
        // Set up the workspace for the godunov Box.
        //
        D_TERM(bndry[0] = getBCArray(State_Type,i,0,1);,
               bndry[1] = getBCArray(State_Type,i,1,1);,
               bndry[2] = getBCArray(State_Type,i,2,1);)

        godunov->Setup(grids[i], dx, dt, 0,
                       xflux, bndry[0].dataPtr(),
                       yflux, bndry[1].dataPtr(),
#if (BL_SPACEDIM == 3)                          
                       zflux, bndry[2].dataPtr(),
#endif
                       U_fpi(), Rho_fpi(), tforces);
        //
        // Loop over the velocity components.
        //
        for (int comp = 0 ; comp < BL_SPACEDIM ; comp++ )
        {
            godunov->AdvectState(grids[i], dx, dt, 
                                 area[0][i], u_mac[0][i], xflux,
                                 area[1][i], u_mac[1][i], yflux,
#if (BL_SPACEDIM == 3)                       
                                 area[2][i], u_mac[2][i], zflux,
#endif
                                 U_fpi(), U_fpi(), tforces, comp,
                                 (*aofs)[i],    comp,
                                 is_conservative[comp],
                                 comp,
                                 bndry[comp].dataPtr(),
                                 volume[i]);
            //
            // Get fluxes for diagnostics and refluxing.
            //
            pullFluxes(i, comp, 1, xflux, yflux, zflux, dt);
        }
    }
    //
    // pullFluxes() contains CrseInit() calls. Got to complete the process.
    //
    if (do_reflux && level < finest_level)
    {
        FluxRegister& fr = getAdvFluxReg(level+1);
        fr.CrseInitFinish();
    }
}

//
// This routine advects the scalars
//

void
NavierStokes::scalar_advection (Real dt,
                                int  fscalar,
                                int  lscalar)
{
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "... advect scalars\n";
    }
    //
    // Get simulation parameters.
    //
    const int num_scalars     = lscalar - fscalar + 1;
    const Real* dx            = geom.CellSize();
    const Real prev_time      = state[State_Type].prevTime();
    const Real prev_pres_time = state[Press_Type].prevTime();
    //
    // Get the viscous terms.
    //
    MultiFab visc_terms(grids,num_scalars,1,Fab_allocate);
    getViscTerms(visc_terms,fscalar,num_scalars,prev_time);
    if (be_cn_theta == 1.0)
        visc_terms.setVal(0.0,1);
    //
    // Set up the grid loop.
    //
    FArrayBox xflux, yflux, zflux, tforces, tvelforces, Gp;

    MultiFab vel_visc_terms;

    int use_forces_in_trans = godunov->useForcesInTrans();

    if (use_forces_in_trans)
    {
      vel_visc_terms.define(grids,BL_SPACEDIM,1,Fab_allocate);
      getViscTerms(vel_visc_terms,Xvel,BL_SPACEDIM,prev_time);
      if (be_cn_theta == 1.0)
          vel_visc_terms.setVal(0.0,1);
    }
    Array<int> state_bc, bndry[BL_SPACEDIM];
    //
    // FillPatch'd state data.
    //
    MultiFab* divu_fp = getDivCond(1,prev_time);

    FillPatchIterator P_fpi(*this,get_old_data(Press_Type),1,
                            prev_pres_time,Press_Type,0,1);

    FillPatchIterator U_fpi(*this,visc_terms,HYP_GROW,prev_time,
                            State_Type,Xvel,BL_SPACEDIM);

    FillPatchIterator S_fpi(*this,visc_terms,HYP_GROW,prev_time,
                            State_Type,fscalar,num_scalars);

    FillPatchIterator Rho_fpi(*this,visc_terms,1,prev_time,
                              State_Type,Density,1);
    //
    // Compute the advective forcing.
    //
    for ( ;
          U_fpi.isValid() && S_fpi.isValid() &&
              Rho_fpi.isValid() && P_fpi.isValid();
          ++U_fpi, ++S_fpi, ++Rho_fpi, ++P_fpi)
    {
        //
        // Since all the MultiFabs are based on same grid just use indices.
        //
        const int i = U_fpi.index();

        setForce(tforces,i,1,fscalar,num_scalars,Rho_fpi());
        
        if (use_forces_in_trans)
        {
            setForce(tvelforces,i,1,Xvel,BL_SPACEDIM,Rho_fpi());

            getGradP(P_fpi(),Gp,grids[i],1);
            
            godunov->Sum_tf_gp_visc(tvelforces,vel_visc_terms[i],Gp,Rho_fpi());
        }
        //
        // Set up the workspace for the godunov Box.
        //
        D_TERM(bndry[0] = getBCArray(State_Type,i,0,1);,
               bndry[1] = getBCArray(State_Type,i,1,1);,
               bndry[2] = getBCArray(State_Type,i,2,1);)

        godunov->Setup(grids[i], dx, dt, 0,
                       xflux, bndry[0].dataPtr(),
                       yflux, bndry[1].dataPtr(),
#if (BL_SPACEDIM == 3)                         
                       zflux, bndry[2].dataPtr(),
#endif
                       U_fpi(), Rho_fpi(), tvelforces);
        //
        // Loop over the scalar components.
        //
        for (int comp = 0 ; comp < num_scalars ; comp++)
        {
            int state_ind = fscalar + comp;
            //
            // Compute total forcing.
            //
            godunov->Sum_tf_divu_visc(S_fpi(), tforces, comp, 1,
                                      visc_terms[i], comp,
                                      (*divu_fp)[i], Rho_fpi(),
                                      is_conservative[state_ind]);
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
                                 U_fpi(), S_fpi(), tforces, comp,
                                 (*aofs)[i],    state_ind,
                                 is_conservative[state_ind],
                                 state_ind,
                                 state_bc.dataPtr(),
                                 volume[i]);
            //
            // Get the fluxes for refluxing and diagnostic purposes.
            //
            pullFluxes(i, state_ind, 1, xflux, yflux, zflux, dt);
        }
    }

    delete divu_fp;
    //
    // pullFluxes() contains CrseInit() calls. Got to complete the process.
    //
    if (do_reflux && level < parent->finestLevel())
    {
        FluxRegister& fr = getAdvFluxReg(level+1);
        fr.CrseInitFinish();
    }
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
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "... update scalars\n";
    }
    scalar_advection_update(dt, first_scalar, last_scalar);
    scalar_diffusion_update(dt, first_scalar, last_scalar);
}

void
NavierStokes::scalar_advection_update (Real dt,
                                       int  first_scalar,
                                       int  last_scalar)
{
    //
    // Simulation parameters.
    //
    MultiFab& S_old      = get_old_data(State_Type);
    MultiFab& S_new      = get_new_data(State_Type);
    MultiFab& Aofs       = *aofs;
    const Real half_time = 0.5*(state[State_Type].curTime()+state[State_Type].prevTime());

    Array<int> state_bc;
    FArrayBox tforces;
    //
    // Loop over the desired scalars.
    //
    for (int sigma = first_scalar; sigma <= last_scalar; sigma++)
    {
        //
        // Compute inviscid estimate of scalars.
        //
        FillPatchIterator Rho_fpi(*this,S_old,0,half_time,State_Type,Density,1);

        for ( ; Rho_fpi.isValid(); ++Rho_fpi)
        {
            DependentMultiFabIterator S_oldmfi(Rho_fpi,S_old);
            DependentMultiFabIterator S_newmfi(Rho_fpi,S_new);
            DependentMultiFabIterator Aofsmfi(Rho_fpi,Aofs);

            int i = Rho_fpi.index();

            setForce(tforces,i,0,sigma,1,Rho_fpi());

            godunov->Add_aofs_tf(S_oldmfi(), S_newmfi(), sigma, 1, Aofsmfi(),
                                 sigma, tforces, 0, grids[i], dt);

	    if (is_conservative[sigma] ? false : true)
            {
                state_bc = getBCArray(State_Type,i,sigma,1);

                godunov->ScalMinMax(S_oldmfi(), S_newmfi(), sigma,
                                    state_bc.dataPtr(), grids[i]);
	    }
        }
    }
}

void
NavierStokes::scalar_diffusion_update (Real dt,
                                       int  first_scalar,
                                       int  last_scalar)
{
    //
    // Loop over the desired scalars.
    //
    for (int sigma = first_scalar; sigma <= last_scalar; sigma++)
    {
        //
        // Compute diffusion.
        //
        if (is_diffusive[sigma])
        {
            int rho_flag = 0;
            diffuse_scalar_setup(sigma, &rho_flag);
            diffusion->diffuse_scalar(dt,sigma,be_cn_theta,rho_half,rho_flag);
        }
    }    
}

void
NavierStokes::diffuse_scalar_setup (int  sigma,
                                    int* rho_flag) 
{
    (*rho_flag) = !is_conservative[sigma] ? 1 : 2;
}

//
// This subroutine updates the velocity field before the level projection.
//
// AT this point in time, all we know is u^n, rho^n+1/2, and the
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
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "... update velocities\n";
    }

    velocity_advection_update(dt);

    if (!initial_iter)
    {
      velocity_diffusion_update(dt);
    }
    else
    {
      initial_velocity_diffusion_update(dt);
    }
}

void
NavierStokes::velocity_advection_update (Real dt)
{
    FArrayBox Gp, tforces;
    //
    // Simulation parameters.
    //
    MultiFab& U_old           = get_old_data(State_Type);
    MultiFab& U_new           = get_new_data(State_Type);
    MultiFab& P_old           = get_old_data(Press_Type);
    MultiFab& Aofs            = *aofs;
    const Real cur_time       = state[State_Type].curTime();
    const Real prev_time      = state[State_Type].prevTime();
    const Real half_time      = 0.5*(prev_time+cur_time);
    const Real pres_prev_time = state[Press_Type].prevTime();

    FillPatchIterator P_fpi(*this,P_old,0,pres_prev_time,Press_Type,0,1);

    FillPatchIterator Rho_fpi(*this,U_old,0,half_time,State_Type,Density,1);

    for ( ; Rho_fpi.isValid() && P_fpi.isValid(); ++Rho_fpi, ++P_fpi)
    {
        //
        // Since all the MultiFabs are based on same grid we'll use indices.
        //
        const int i = Rho_fpi.index();
        //
        // Get the forcing terms.
        //
        getGradP(P_fpi(),Gp,grids[i],0);

	setForce(tforces,i,0,Xvel,BL_SPACEDIM,Rho_fpi());
        //
        // Do following only at initial iteration--per JBB.
        //
        if (initial_iter && is_diffusive[Xvel])
        {
            tforces.setVal(0);
        }

        godunov->Add_aofs_tf_gp(U_old[i],U_new[i],Aofs[i],tforces,
                                 Gp,Rho_fpi(),grids[i],dt);
    }
}

void
NavierStokes::velocity_diffusion_update (Real dt)
{
    //
    // Compute the viscous forcing.
    // Do following except at initial iteration--rbp, per jbb.
    //
    if (is_diffusive[Xvel])
    {
        MultiFab* delta_rhs = 0;
        diffuse_velocity_setup(dt, delta_rhs);
        diffusion->diffuse_velocity(dt, be_cn_theta, rho_half, 1, delta_rhs);
        delete delta_rhs;
    }
}

void
NavierStokes::diffuse_velocity_setup (Real       dt,
                                      MultiFab*& delta_rhs)
{
    if (S_in_vel_diffusion && have_divu)
    {
        Real time = state[State_Type].prevTime();
        delta_rhs = new MultiFab(grids,BL_SPACEDIM,0,Fab_allocate);
      
        MultiFab divmusi(grids,BL_SPACEDIM,0,Fab_allocate);

        diffusion->compute_divmusi(time,visc_coef[Xvel],divmusi);
        divmusi.mult((1./3.)*(1.0-be_cn_theta),0,BL_SPACEDIM,0);
        (*delta_rhs).copy(divmusi,0,0,BL_SPACEDIM);

        diffusion->compute_divmusi(time+dt,visc_coef[Xvel],divmusi);
        divmusi.mult((1./3.)*be_cn_theta,0,BL_SPACEDIM,0);
        (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);
    }
}

void
NavierStokes::initial_velocity_diffusion_update (Real dt)
{
    //
    // Simulation parameters.
    //
    MultiFab& U_old           = get_old_data(State_Type);
    MultiFab& U_new           = get_new_data(State_Type);
    MultiFab& Aofs            = *aofs;
    const Real prev_time      = state[State_Type].prevTime();
    const Real pres_prev_time = state[Press_Type].prevTime();
    //
    // Do following only at initial iteration -- rbp, per jbb.
    //
    if (is_diffusive[Xvel])
    {
        //
        // Get viscous forcing terms.
        //
        MultiFab visc_terms(grids,BL_SPACEDIM,1,Fab_allocate);
        getViscTerms(visc_terms,Xvel,BL_SPACEDIM,prev_time);
        if (be_cn_theta == 1.0)
            visc_terms.setVal(0.0,1);
        FArrayBox tforces, Gp;
        //
        // Update U_new with viscosity.
        //
        FillPatchIterator P_fpi(*this,get_old_data(Press_Type),0,
                                pres_prev_time,Press_Type,0,1);

        FillPatchIterator Rho_fpi(*this,visc_terms,0,prev_time,
                                  State_Type,Density,1);

        for ( ; Rho_fpi.isValid() && P_fpi.isValid(); ++Rho_fpi, ++P_fpi)
        {
            int i = Rho_fpi.index();

            setForce(tforces,i,0,Xvel,BL_SPACEDIM,Rho_fpi());

            getGradP(P_fpi(),Gp,grids[i],0);

            godunov->Sum_tf_gp_visc(tforces,visc_terms[i],Gp,(*rho_half)[i]);
            //
            // Compute the inviscid update.
            //
            godunov->Add_aofs_tf(U_old[i], U_new[i], 0, BL_SPACEDIM, Aofs[i],
                                 0, tforces, 0, grids[i], dt);
        }
    }
}

void
NavierStokes::errorEst (TagBoxArray& tags,
                        int          clearval,
                        int          tagval,
                        Real         time)
{
    const int*  domain_lo = geom.Domain().loVect();
    const int*  domain_hi = geom.Domain().hiVect();
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    for (int j = 0; j < err_list.length(); j++)
    {
        MultiFab* mf = derive(err_list[j]->name(), time, err_list[j]->nGrow());

        for (MultiFabIterator mfi(*mf); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab    = mfi();
            Array<int>  itags = tags[mfi.index()].tags();
            int*        tptr  = itags.dataPtr();
            const int*  tlo   = tags[mfi.index()].box().loVect();
            const int*  thi   = tags[mfi.index()].box().hiVect();
            const int*  lo    = mfi.validbox().loVect();
            const int*  hi    = mfi.validbox().hiVect();
            const Real* xlo   = grid_loc[mfi.index()].lo();
            Real*      dat    = mfi().dataPtr();
            const int* dlo    = mfi().box().loVect();
            const int* dhi    = mfi().box().hiVect();
            const int  ncomp  = mfi().nComp();

            err_list[j]->errFunc()(tptr, ARLIM(tlo), ARLIM(thi), &tagval,
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
NavierStokes::sumDerive (const aString& name, Real time)
{
    Real sum     = 0.0;
    MultiFab* mf = derive(name,time,0);

    assert(!(mf == 0));

    for (MultiFabIterator mfi(*mf); mfi.isValid(); ++mfi)
    {
        if (level < parent->finestLevel())
        {
            const BoxArray& f_box = parent->boxArray(level+1);

            for (int j = 0; j < f_box.length(); j++)
            {
                Box c_box = ::coarsen(f_box[j],fine_ratio);
                if (c_box.intersects(grids[mfi.index()]))
                    mfi().setVal(0.0,(c_box & grids[mfi.index()]),0);
            }
        }

        sum += mfi().sum(0);
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
NavierStokes::volWgtSum (const aString& name,
                         Real           time)
{
    Real sum       = 0.0;
    int rz_flag    = CoordSys::IsRZ() ? 1 : 0;
    const Real* dx = geom.CellSize();
    MultiFab* mf   = derive(name,time,0);

    for (MultiFabIterator mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = mfi();

        if (level < parent->finestLevel())
        {
            const BoxArray& f_box = parent->boxArray(level+1);

            for (int j = 0; j < f_box.length(); j++)
            {
                Box c_box = ::coarsen(f_box[j],fine_ratio);
                if (c_box.intersects(grids[mfi.index()]))
                    mfi().setVal(0.0,(c_box & grids[mfi.index()]),0);
            }
        }
        Real s;
        const Real* dat = mfi().dataPtr();
        const int* dlo  = mfi().loVect();
        const int* dhi  = mfi().hiVect();
        const int* lo   = grids[mfi.index()].loVect();
        const int* hi   = grids[mfi.index()].hiVect();
        Real* rad       = &radius[mfi.index()];

        Array<Real> tmp(hi[1]-lo[1]+1);

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

aString
NavierStokes::thePlotFileType () const
{
    //
    // Increment this whenever the writePlotFile() format changes.
    //
    static const aString the_plot_file_type("NavierStokes-V1.1");

    return the_plot_file_type;
}

void
NavierStokes::writePlotFile (const aString& dir,
                             ostream&       os,
                             VisMF::How     how)
{
    int i, n;

    bool PlotDivu = have_divu && parent->isPlotVar(desc_lst[Divu_Type].name(0));
    bool PlotDsdt = have_dsdt && parent->isPlotVar(desc_lst[Dsdt_Type].name(0));
    //
    // The list of indices of State to write to plotfile.
    //
    vector<int> idx_map;

    for (i = 0; i < NUM_STATE; i++)
        if (parent->isPlotVar(desc_lst[State_Type].name(i)))
            idx_map.push_back(i);

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        int n_data_items = idx_map.size();

        if (PlotDivu) n_data_items++;
        if (PlotDsdt) n_data_items++;

        if (n_data_items == 0)
            BoxLib::Error("Must specify at least one valid data item to plot");

        os << n_data_items << '\n';

        for (n = 0; n < idx_map.size(); n++)
            os << desc_lst[State_Type].name(idx_map[n]) << '\n';

        if (PlotDivu) os << desc_lst[Divu_Type].name(0) << '\n';
        if (PlotDsdt) os << desc_lst[Dsdt_Type].name(0) << '\n';

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
    //
    // Now write state data.
    //
    Real cur_time      = state[State_Type].curTime();
    MultiFab& cell_dat = state[State_Type].newData();
    MultiFab* divu_dat = 0;
    MultiFab* dsdt_dat = 0;

    if (PlotDivu) divu_dat = &state[Divu_Type].newData();
    if (PlotDsdt) dsdt_dat = &state[Dsdt_Type].newData();
    //
    // There may be up to three MultiFabs written out at each level.
    //
    static const aString BaseName[] =
    {
        aString("/Cell"), aString("/DivU"), aString("/DsDt")
    };
    //
    // Build the directory to hold the MultiFabs at this level.
    // The name is relative to the directory containing the Header file.
    //
    char buf[64];
    sprintf(buf, "Level_%d", level);
    aString Level = buf;
    //
    // Now for the full pathname of that directory.
    //
    aString FullPath = dir;
    if (!FullPath.isNull() && FullPath[FullPath.length()-1] != '/')
        FullPath += '/';
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!Utility::CreateDirectory(FullPath, 0755))
            Utility::CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.length() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (i = 0; i < cell_dat.boxArray().length(); ++i)
        {
            for (n = 0; n < BL_SPACEDIM; n++)
                os << grid_loc[i].lo(n) << ' ' << grid_loc[i].hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (idx_map.size() > 0)
        {
            aString PathNameInHeader = Level;
            PathNameInHeader += BaseName[0];
            os << PathNameInHeader << '\n';
        }
        if (PlotDivu)
        {
            aString PathNameInHeader = Level;
            PathNameInHeader += BaseName[1];
            os << PathNameInHeader << '\n';
        }
        if (PlotDsdt)
        {
            aString PathNameInHeader = Level;
            PathNameInHeader += BaseName[2];
            os << PathNameInHeader << '\n';
        }
    }
    //
    // Use the Full pathname when naming the MultiFab.
    //
    aString TheFullPath = FullPath;
    TheFullPath += BaseName[0];

    if (idx_map.size() == NUM_STATE)
    {
        RunStats::addBytes(VisMF::Write(cell_dat, TheFullPath, how, true));
    }
    else if (idx_map.size() > 0)
    {
        //
        // Make MultiFab containing copy of selected components.
        //
        // Note that we don't copy the ghost cells.
        //
        MultiFab mf(cell_dat.boxArray(), idx_map.size(), 0);

        for (MultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
        {
            DependentMultiFabIterator dmfi(mfi,cell_dat);

            for (i = 0; i < idx_map.size(); i++)
            {
                mfi().copy(dmfi(), idx_map[i], i, 1);
            }
        }
        RunStats::addBytes(VisMF::Write(mf, TheFullPath, how));
    }

    if (PlotDivu)
    {
        TheFullPath = FullPath;
        TheFullPath += BaseName[1];
        RunStats::addBytes(VisMF::Write(*divu_dat, TheFullPath, how, true));
    }
    if (PlotDsdt)
    {
        TheFullPath = FullPath;
        TheFullPath += BaseName[2];
        RunStats::addBytes(VisMF::Write(*dsdt_dat, TheFullPath, how, true));
    }
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
                ratio *= parent->MaxRefRatio(lev-1);
            }
            factor = 1.0/((double) ratio);
        }

        return factor*fixed_dt;
    }
    const int n_grow    = 0;
    Real estdt          = 1.0e+20;
    const Real cur_time = state[State_Type].curTime();
    MultiFab& P_new     = get_new_data(Press_Type);
    MultiFab& U_new     = get_new_data(State_Type);
    const Real* dx      = geom.CellSize();

    Real gr_max[BL_SPACEDIM], u_max[BL_SPACEDIM] = {0};

    FArrayBox p_fab, Gp, tforces;

    FillPatchIterator Rho_fpi(*this,U_new,n_grow,cur_time,State_Type,Density,1);

    for ( ; Rho_fpi.isValid(); ++Rho_fpi)
    {
        const int i = Rho_fpi.index();
        //
        // Get the pressure.
        //
        p_fab.resize(::surroundingNodes(grids[i]),1);
        p_fab.copy(P_new[i],p_fab.box());
        //
        // Get the velocity forcing.  For some reason no viscous forcing.
        //
        setForce(tforces,i,n_grow,Xvel,BL_SPACEDIM,Rho_fpi());
        getGradP(p_fab, Gp, grids[i], 0);
        tforces.minus(Gp,0,0,BL_SPACEDIM);
        //
        // Estimate the maximum allowable timestep from the Godunov box.
        //
        Real dt = godunov->estdt(U_new[i], tforces, Rho_fpi(),
                                 grids[i], geom.CellSize(), cfl, gr_max);

        for (int k = 0; k < BL_SPACEDIM; k++)
        {
	    u_max[k] = Max(u_max[k],gr_max[k]);
	}
	estdt = Min(estdt,dt);
    }
    //
    // Parallel reductions.
    //
    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        ParallelDescriptor::ReduceRealMax(u_max[k]);
    }
    ParallelDescriptor::ReduceRealMin(estdt);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "estTimeStep :: \n";
        cout << "LEV = " << level << " UMAX = ";
        for (int k = 0; k < BL_SPACEDIM; k++)
        {
            cout << u_max[k] << "  ";
        }
        cout << NL;
    }
    Real vel_max = u_max[0];
    for (int k = 1; k < BL_SPACEDIM; k++)
    {
        vel_max = Max(vel_max,u_max[k]);
    }
    if (vel_max < 1.0e-10)
    {
        const Real grav = Abs(gravity);
        vel_max         = (grav > 1.0e-5) ? grav : 1.0;
        Real dx_min     = dx[0];
        for (int k = 1; k < BL_SPACEDIM; k++)
        {
            dx_min = Min(dx_min,dx[k]);
        }
        estdt = cfl*dx_min/vel_max;
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
        dt_min[i] = getLevel(i).estTimeStep();
        dt_min[i] = Min(dt_min[i],change_max*dt_level[i]);
        n_factor *= n_cycle[i];
        dt_0      = Min(dt_0,n_factor*dt_min[i]);
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
        dt_level[i] = dt_0/( (float)n_factor );
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
        dt_0        = Min(dt_0,n_factor*dt_level[i]);
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
        dt_level[i] = dt_0/( (float)n_factor );
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
    const Real strt_time   = state[State_Type].curTime();
    const int finest_level = parent->finestLevel();

    dt_init = 1.0e+100;

    for (int k = 0; k <= finest_level; k++)
    {
        nc_save[k] = parent->nCycle(k);
        dt_save[k] = getLevel(k).initialTimeStep();
        dt_init    = Min(dt_init,dt_save[k]);
    }

    Array<Real> dt_level(finest_level+1,dt_init);
    Array<int> n_cycle(finest_level+1,1);

    Real dt0      = dt_save[0];
    int  n_factor = 1;
    for (int k = 0; k <= finest_level; k++)
    {
        n_factor *= nc_save[k];
        dt0       = Min(dt0,n_factor*dt_save[k]);
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
        dt_save[k] = dt0/( (float) n_factor);
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
NavierStokes::post_timestep ()
{
    const int finest_level = parent->finestLevel();

    MultiFab& pres = get_new_data(Press_Type);
    //
    // Reflux .
    //
    if (do_reflux && (level < finest_level))
    {
        reflux();
    }
    //
    // Average down.
    //
    if (level < finest_level)
    {
        avgDown();
        const Real dt        = parent->dtLevel(level);
        const Real prev_time = state[State_Type].prevTime();
        const Real half_time = prev_time + 0.5*dt;

        FillPatchIterator rho_halffpi(*this,*rho_half,1,half_time,
                                      State_Type,Density,1);

        for ( ; rho_halffpi.isValid(); ++rho_halffpi)
        {
            (*rho_half)[rho_halffpi.index()].copy(rho_halffpi());
        }
    }
    //
    // Mac Sync Correction.
    //
    if (do_mac_proj && (level < finest_level))
    {
        static const aString RunstatString("mac_sync");
        RunStats mac_sync_stats(RunstatString, level);
        mac_sync_stats.start();
        mac_sync();
        mac_sync_stats.end();
    }
    //
    // Level Sync Correction.
    //
    if (do_sync_proj && (level < finest_level))
    {
        level_sync();
    }
    //
    // Test for conservation.
    //
    if (level == 0       &&
        sum_interval > 0 &&
        (parent->levelSteps(0)%sum_interval == 0))
    {
        sum_integrated_quantities();
    }
}

//
// Build any additional data structures after restart.
//

void NavierStokes::post_restart() {}

//
// Build any additional data structures after regrid.
//

void
NavierStokes::post_regrid (int lbase,
                                int new_finest)
{
    if (projector && level == lbase)
        projector->setFinestLevel(new_finest);
}

//
// Ensure state, and pressure are consistent.
//

void
NavierStokes::post_init (Real stop_time)
{
    if (level > 0)
        //
        // Nothing to sync up at level > 0.
        //
        return;

    MultiFab& P_new = get_new_data(Press_Type);
    MultiFab& P_old = get_old_data(Press_Type);

    const int finest_level = parent->finestLevel();

    Real dt_init = 0.0;

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
    {
        sum_integrated_quantities();
    }
}

//
// MULTILEVEL SYNC FUNCTIONS
//

void
NavierStokes::initRhoAvg (Real alpha)
{
    MultiFab& S_new = get_new_data(State_Type);

    rho_avg->setVal(0);

    for (MultiFabIterator rho_avgmfi(*rho_avg); rho_avgmfi.isValid();
         ++rho_avgmfi)
    {
        DependentMultiFabIterator S_newmfi(rho_avgmfi, S_new);

        rho_avgmfi().copy(S_newmfi(),
                          S_newmfi.validbox(),
                          Density,
                          S_newmfi.validbox(),
                          0,
                          1);
        rho_avgmfi().mult(alpha);
    }
}

void
NavierStokes::incrRhoAvg (Real alpha)
{
    MultiFab& S_new = get_new_data(State_Type);

    for (MultiFabIterator S_newmfi(S_new); S_newmfi.isValid(); ++S_newmfi)
    {
        DependentMultiFabIterator rho_avgmfi(S_newmfi,*rho_avg);

        const int* lo      = S_newmfi.validbox().loVect();
        const int* hi      = S_newmfi.validbox().hiVect();
        const int* rlo     = S_newmfi().loVect();
        const int* rhi     = S_newmfi().hiVect();
        const Real* rhodat = S_newmfi().dataPtr(Density);
        const int* alo     = rho_avgmfi().loVect();
        const int* ahi     = rho_avgmfi().hiVect();
        Real* avgdat       = rho_avgmfi().dataPtr();

        FORT_INCRMULT(avgdat,ARLIM(alo),ARLIM(ahi),
                      rhodat,ARLIM(rlo),ARLIM(rhi),lo,hi,&alpha);
    }
}

void
NavierStokes::incrPAvg (int  iteration,
                        Real alpha)
{
    if (iteration == 1)
        //
        // Set p_avg to zero for the initial iteration.
        //
        p_avg->setVal(0.);
    //
    // Then increment p_avg with fine grid pressure.
    //
    MultiFab& P_new = get_new_data(Press_Type);

    for (MultiFabIterator P_newmfi(P_new); P_newmfi.isValid(); ++P_newmfi)
    {
        DependentMultiFabIterator p_avgmfi(P_newmfi,*p_avg);

        const int* lo    = P_newmfi.validbox().loVect();
        const int* hi    = P_newmfi.validbox().hiVect();
        const int* p_lo  = P_newmfi().loVect();
        const int* p_hi  = P_newmfi().hiVect();
        const Real* pdat = P_newmfi().dataPtr();
        const int* alo   = p_avgmfi().loVect();
        const int* ahi   = p_avgmfi().hiVect();
        Real* avgdat     = p_avgmfi().dataPtr();

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
    const int finest_level = parent->finestLevel();
    const Real divu_time   = have_divu ? state[Divu_Type].curTime()
                                       : state[Press_Type].curTime();
    if (projector)
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
    //
    // Zero pressure field.
    //
    for (int k = 0; k <= finest_level; k++)
    {
        NavierStokes& ns_level = getLevel(k);
        ns_level.zeroNewPress();
        ns_level.zeroOldPress();
    }
}

//
// Initialize the pressure by iterating the initial timestep
//

void
NavierStokes::post_init_press (Real&        dt_init,
                               Array<int>&  nc_save,
                               Array<Real>& dt_save )
{
    const Real strt_time   = state[State_Type].curTime();
    const int finest_level = parent->finestLevel();

    MultiFab& P_new = get_new_data(Press_Type);
    MultiFab& P_old = get_old_data(Press_Type);

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
            sig[k] = getLevel(k).rho_half;
        }
        if (projector)
        {
            projector->initialSyncProject(0,sig,parent->dtLevel(0),strt_time,
                                          dt_init,have_divu);
        }
        delete sig;

        for (int k = finest_level-1; k>= 0; k--)
        {
            getLevel(k).avgDown();
        }
        for (int k = 0; k <= finest_level; k++)
        {
            //
            // Reset state variables to initial time, do not reset
            // pressure variable.
            //
            getLevel(k).resetState(strt_time, dt_save[k], dt_save[k]);
        }
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
    }

    parent->setDtLevel(dt_save);
    parent->setNCycle(nc_save);
}

//
// Interpolate A cell centered Sync correction from a
// coarse level (c_lev) to a fine level (f_lev).
//

void
NavierStokes::SyncInterp (MultiFab& CrseSync,
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
                          int       which_interp)
{
    int i,n,dir;

    if (ParallelDescriptor::NProcs() > 1)
      BoxLib::Abort("NavierStokes::SyncInterp(): not implemented in parallel");
    //
    // This routine interpolates the num_comp components of CrseSync
    // (starting at src_comp) and either increments or puts the result into
    // the num_comp components of FineSync (starting at dest_comp)
    // The components of bc_orig_qty corespond to the quantities of CrseSync.
    //
    assert(which_interp >= 0 && which_interp <= 2);

    Interpolater* interpolater=0;
    interpolater = (which_interp == 0) ? &cell_cons_interp : interpolater;
    interpolater = (which_interp == 1) ? &pc_interp : interpolater;
    interpolater = (which_interp == 2) ? &unlimited_cc_interp : interpolater;

    Array<BCRec> bc_interp(num_comp);
    Array<IntVect> pshifts(27);

    FArrayBox cdata,fdata;
    //
    // Get fine parameters.
    //
    NavierStokes& fine_level = getLevel(f_lev);
    const BoxArray& fgrids   = fine_level.boxArray();
    const Geometry& fgeom    = parent->Geom(f_lev);
    int nfine                = fgrids.length();
    const Box& domain        = fgeom.Domain();
    //
    // Get coarse parameters.
    //
    NavierStokes& crse_level = getLevel(c_lev);
    const BoxArray& cgrids   = crse_level.boxArray();
    const Geometry& cgeom    = parent->Geom(c_lev);
    const Real* dx_crse      = cgeom.CellSize();

    const Box cdomain(coarsen(domain,ratio));
    const int* cdomlo = cdomain.loVect();
    const int* cdomhi = cdomain.hiVect();

    int* bc_new =  new int[2*BL_SPACEDIM*(src_comp+num_comp)];
    //
    // Loop over fine grids.
    //
    for (i = 0; i < nfine; i++)
    {
        //
        // Create storage for interpolation.
        //
        const Box& grd = fgrids[i];
        Box cgrd = interpolater->CoarseBox(grd,ratio);

        fdata.resize(grd,num_comp);
        cdata.resize(cgrd,num_comp);
        cdata.setVal(0);
        CrseSync.copy(cdata,src_comp,0,num_comp);

        const int* clo  = cdata.loVect();
        const int* chi  = cdata.hiVect();
        const Real* xlo = fine_level.grid_loc[i].lo();

        for (n = 0; n < num_comp; n++)
        {
            for (dir = 0; dir < BL_SPACEDIM; dir++)
            {
                int bc_index = (n+src_comp)*(2*BL_SPACEDIM) + dir;
                bc_new[bc_index]             = INT_DIR;
                bc_new[bc_index+BL_SPACEDIM] = INT_DIR;
 
                if (clo[dir] < cdomlo[dir] || chi[dir] > cdomhi[dir])
                {
                    for (int crse = 0; crse < cgrids.length(); crse++)
                    {
                        const int* crse_lo = cgrids[crse].loVect();
                        const int* crse_hi = cgrids[crse].hiVect();

                        if (clo[dir] < cdomlo[dir] &&
                            crse_lo[dir] == cdomlo[dir])
                        {
                            bc_new[bc_index] = bc_orig_qty[crse][bc_index];
                        }
                        if (chi[dir] > cdomhi[dir] &&
                            crse_hi[dir] == cdomhi[dir]) {
                            bc_new[bc_index+BL_SPACEDIM] = bc_orig_qty[crse][bc_index+BL_SPACEDIM]; 
                        }
                    }
                }
            }
            FORT_FILCC(cdata.dataPtr(n), ARLIM(clo), ARLIM(chi),
                       cdomlo, cdomhi, dx_crse, xlo,
                       &(bc_new[2*BL_SPACEDIM*(n+src_comp)]));
        }
        //
        // Fill in periodic images.
        //
        if (cgeom.isAnyPeriodic())
        {
            cgeom.periodicShift(cgeom.Domain(), cgrd, pshifts);

            for (int iiv = 0; iiv < pshifts.length(); iiv++)
            {
                cdata.shift(pshifts[iiv]);
                CrseSync.copy(cdata,src_comp,0,num_comp);
                cdata.shift(-pshifts[iiv]);
            }
        }
        //
        // Set the boundary condition array for interpolation.
        //
        for (n = 0; n < num_comp; n++)
        {
            for (int dir = 0; dir < BL_SPACEDIM; dir++)
            {
                int bc_index = (n+src_comp)*(2*BL_SPACEDIM) + dir;
                bc_interp[n].setLo(dir,bc_new[bc_index]);
                bc_interp[n].setHi(dir,bc_new[bc_index+BL_SPACEDIM]);
            }
        }
        //
        // Scale coarse interpolant for anelastic.
        //
        ScaleCrseSyncInterp(cdata, c_lev, num_comp);
        //
        // Compute the interpolated correction.
        //
        interpolater->interp(cdata,0,
                             fdata,0,num_comp,grd,ratio,
                             cgeom,fgeom,bc_interp);
        //
        // Rescale Fine interpolant for anelastic.
        //
        reScaleFineSyncInterp(fdata, f_lev, num_comp);
        //
        // Set Fine Sync equal to the correction or add it in.
        //
        if (increment)
        {
            fdata.mult(dt_clev);
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
                              int       f_lev,
                              IntVect&  ratio )
{
    if (ParallelDescriptor::NProcs() > 1)
        BoxLib::Abort("NavierStokes::SyncProjInterp(): not implemented in parallel");
    //
    // Because of phi.copy(crse_phi) within multifab loop below ...
    //
    // Get fine parameters.
    //
    NavierStokes& fine_level = getLevel(f_lev);
    const BoxArray& fgrids   = fine_level.boxArray();
    const Geometry& fgeom    = parent->Geom(f_lev);
    const BoxArray& P_grids  = P_new.boxArray();
    const int nfine          = fgrids.length();
    //
    // Get coarse parameters.
    //
    NavierStokes& crse_level = getLevel(c_lev);
    const Geometry& cgeom    = parent->Geom(c_lev);

    int i;
    FArrayBox fine_phi, crse_phi;
    Array<BCRec> bc(BL_SPACEDIM);
    //
    // Loop over the fine grids.
    //
    for (i = 0; i < nfine; i++)
    {
        //
        // Get space for bi/trilinear interpolation.
        //
        const Box& fbox = P_grids[i];
        const Box  cbox(node_bilinear_interp.CoarseBox(fbox,ratio));
        fine_phi.resize(fbox,1);
        crse_phi.resize(cbox,1);
        fine_phi.setVal(1.e30);
        crse_phi.setVal(1.e30);
        phi.copy(crse_phi);
        //
        // Compute the interpolated pressure sync correction.
        //
        node_bilinear_interp.interp(crse_phi,0,
                                    fine_phi,0,1,
                                    fbox,ratio,cgeom,fgeom,bc);
        //
        // Add the fine corrections to the state.
        //
        P_new[i].plus(fine_phi);
    }
}

//
// This function averages a multifab of fine data down onto
// a multifab of coarse data.
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
                       int             strt_comp,
                       int             num_comp,
                       IntVect&        fratio)
{
    assert(S_crse.nComp() == S_fine.nComp());

    MultiFabCopyDescriptor mfcd;

    MultiFabId mfidS_fine  = mfcd.RegisterFabArray(&S_fine);
    MultiFabId mfidFineVol = mfcd.RegisterFabArray(&fvolume);

    vector<FillBoxId> fillBoxIdList, fillBoxIdListVol;
    //
    // Loop over coarse grids and intersect coarse with fine.
    //
    for (ConstMultiFabIterator mfi(S_crse); mfi.isValid(); ++mfi)
    {
        assert(grids[mfi.index()] == mfi.validbox());

        for (int fine = 0; fine < fgrids.length(); fine++)
        {
            Box ovlp = ::coarsen(fgrids[fine],fratio) & mfi.validbox();

            if (ovlp.ok())
            {
              Box fine_ovlp = ::refine(ovlp,fratio);

              fillBoxIdList.push_back(mfcd.AddBox(mfidS_fine,
                                                  fine_ovlp,
                                                  0,
                                                  fine,
                                                  strt_comp,
                                                  0,
                                                  num_comp));

              assert(fillBoxIdList.back().box() == fine_ovlp);
              //
              // Also save the index of the coarse FAB needed filling.
              //
              fillBoxIdList.back().FabIndex(mfi.index());

              fillBoxIdListVol.push_back(mfcd.AddBox(mfidFineVol,
                                                     fine_ovlp,
                                                     0,
                                                     fine,
                                                     0,
                                                     0,
                                                     1));

              assert(fillBoxIdListVol.back().box() == fine_ovlp);
              //
              // Here we'll save the fine index so we can reconstruct `ovlp'.
              //
              fillBoxIdListVol.back().FabIndex(fine);
            }
        }
    }

    mfcd.CollectData();

    assert(fillBoxIdList.size() == fillBoxIdListVol.size());

    const int MyProc = ParallelDescriptor::MyProc();

    FArrayBox fine_fab, fine_vol;

    for (int i = 0; i < fillBoxIdList.size(); i++)
    {
        int c_idx = fillBoxIdList[i].FabIndex();
        int f_idx = fillBoxIdListVol[i].FabIndex();

        Box ovlp = ::coarsen(fgrids[f_idx],fratio) & cgrids[c_idx];

        assert(ovlp.ok());
        assert(S_crse.DistributionMap()[c_idx] == MyProc);
        assert(cvolume.DistributionMap()[c_idx] == MyProc);

        const FillBoxId& fbidFine    = fillBoxIdList[i];
        const FillBoxId& fbidFineVol = fillBoxIdListVol[i];

        fine_fab.resize(fbidFine.box(), num_comp);
        fine_vol.resize(fbidFineVol.box(), 1);

        mfcd.FillFab(mfidS_fine,  fbidFine,    fine_fab);
        mfcd.FillFab(mfidFineVol, fbidFineVol, fine_vol);

        const int* ovlo = ovlp.loVect();
        const int* ovhi = ovlp.hiVect();
        //
        // Get the fine data.
        //
        const int* flo     = fine_fab.loVect();
        const int* fhi     = fine_fab.hiVect();
        const Real* f_dat  = fine_fab.dataPtr();
        const int* fvlo    = fine_vol.loVect();
        const int* fvhi    = fine_vol.hiVect();
        const Real* fv_dat = fine_vol.dataPtr();
        //
        // Get the coarse data.
        //
        FArrayBox& cfab    = S_crse[c_idx];
        FArrayBox& cfabvol = cvolume[c_idx];
        const int* clo     = cfab.loVect();
        const int* chi     = cfab.hiVect();
        const Real* c_dat  = cfab.dataPtr(strt_comp);
        const int* cvlo    = cfabvol.loVect();
        const int* cvhi    = cfabvol.hiVect();
        const Real* cv_dat = cfabvol.dataPtr();

        FORT_AVGDOWN (c_dat,ARLIM(clo),ARLIM(chi),&num_comp,
                      f_dat,ARLIM(flo),ARLIM(fhi),
                      cv_dat,ARLIM(cvlo),ARLIM(cvhi),
                      fv_dat,ARLIM(fvlo),ARLIM(fvhi),
                      ovlo,ovhi,fratio.getVect());
    }
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
                       int              strt_comp,
                       int              num_comp,
                       IntVect&         fratio)
{
    const int* ovlo = ovlp.loVect();
    const int* ovhi = ovlp.hiVect();
    //
    // Get the fine data.
    //
    const int* flo            = fine_fab.loVect();
    const int* fhi            = fine_fab.hiVect();
    const Real* f_dat         = fine_fab.dataPtr(strt_comp);
    const int* fvlo           = fine_vol.loVect();
    const int* fvhi           = fine_vol.hiVect();
    const Real* fv_dat        = fine_vol.dataPtr();
    //
    // Get the coarse data.
    //
    const int* clo            = crse_fab.loVect();
    const int* chi            = crse_fab.hiVect();
    const Real* c_dat         = crse_fab.dataPtr(strt_comp);
    const int* cvlo           = crse_vol.loVect();
    const int* cvhi           = crse_vol.hiVect();
    const Real* cv_dat        = crse_vol.dataPtr();

    FORT_AVGDOWN(c_dat,ARLIM(clo),ARLIM(chi),&num_comp,
                 f_dat,ARLIM(flo),ARLIM(fhi),
                 cv_dat,ARLIM(cvlo),ARLIM(cvhi),
                 fv_dat,ARLIM(fvlo),ARLIM(fvhi),
                 ovlo,ovhi,fratio.getVect());
}

void
NavierStokes::level_sync ()
{
    //
    // Get parameters.
    //
    const Real* dx         = geom.CellSize();
    IntVect ratio          = parent->refRatio(level);
    int crse_dt_ratio      = parent->MaxRefRatio(level);
    const int finest_level = parent->finestLevel();
    Real dt                = parent->dtLevel(level);
    const Real prev_time   = state[State_Type].prevTime();
    const Real half_time   = prev_time + 0.5*dt;
    //
    // Get objects.
    //
    MultiFab& pres              = get_new_data(Press_Type);
    MultiFab& vel               = get_new_data(State_Type);
    SyncRegister& rhs_sync_reg  = getLevel(level+1).getSyncReg();
    SyncRegister* crsr_sync_ptr = 0;
    NavierStokes& fine_level    = getLevel(level+1);
    MultiFab& pres_fine         = fine_level.get_new_data(Press_Type);
    //
    // Get rho at t^{n+1/2}.
    //
    for (FillPatchIterator rho_halffpi(*this,*rho_half,1,half_time,
                                       State_Type,Density,1);
        rho_halffpi.isValid(); ++rho_halffpi)
    {
        (*rho_half)[rho_halffpi.index()].copy(rho_halffpi());
    }
    
    if (level > 0)
    {
        crsr_sync_ptr = &(getLevel(level).getSyncReg());
    }
    //
    // Get boundary conditions.
    //
    Array<int*> sync_bc(grids.length());
    Array< Array<int> > sync_bc_array(grids.length());

    for (int i = 0; i < grids.length(); i++)
    {
        sync_bc_array[i] = getBCArray(State_Type,i,Xvel,BL_SPACEDIM);
        sync_bc[i] = sync_bc_array[i].dataPtr();
    }
    //
    // Multilevel Sync projection or single level.
    //
    if (do_MLsync_proj)
    {
        
        MultiFab& vel_fine          = fine_level.get_new_data(State_Type);
        MultiFab& rho_fine          = *fine_level.rho_avg;
        const Geometry& fine_geom   = parent->Geom(level+1);
        const Geometry& crse_geom   = parent->Geom(level);
        const BoxArray& finegrids   = vel_fine.boxArray();
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
        projector->MLsyncProject(level,pres,vel, pres_fine,vel_fine,
                                 *rho_half, rho_fine, Vsync,V_corr,phi,
                                 &rhs_sync_reg,crsr_sync_ptr,sync_bc.dataPtr(),
                                 dx,dt,ratio,crse_dt_ratio, fine_geom,geom);
        //
        // Correct pressure and velocities after the projection.
        //
        ratio = IntVect::TheUnitVector();
        Array<int*> fine_sync_bc(finegrids.length());
        Array< Array<int> > fine_sync_bc_array(finegrids.length());

        for (int i = 0; i < finegrids.length(); i++)
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
            MultiFab& P_new        = fine_lev.get_new_data(Press_Type);
            MultiFab& U_new        = fine_lev.get_new_data(State_Type);

            SyncInterp(V_corr, level+1, U_new, lev, ratio,
                       0, 0, BL_SPACEDIM, 1 , dt, fine_sync_bc.dataPtr());
            SyncProjInterp(phi, level+1, P_new, lev, ratio);
        }
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
        projector->syncProject(level,pres,vel,rho_half,Vsync,phi,
                               &rhs_sync_reg,crsr_sync_ptr,sync_boxes,
                               sync_bc.dataPtr(),geom,dx,dt,crse_dt_ratio);
        //
        // Correct pressure and velocities after the projection.
        //
        ratio = IntVect::TheUnitVector(); 
        for (int lev = level+1; lev <= finest_level; lev++)
        {
            ratio                 *= parent->refRatio(lev-1);
            NavierStokes& fine_lev = getLevel(lev);
            MultiFab& P_new        = fine_lev.get_new_data(Press_Type);
            MultiFab& U_new        = fine_lev.get_new_data(State_Type);
            
            SyncInterp(*Vsync, level, U_new, lev, ratio,
                       0, 0, BL_SPACEDIM, 1 , dt, sync_bc.dataPtr());
            SyncProjInterp(phi, level, P_new, lev, ratio);
        }
    }
}

//
// The Mac Sync correction function
//

void
NavierStokes::mac_sync ()
{
    const int numscal = NUM_STATE - BL_SPACEDIM;
    
    const Real prev_time      = state[State_Type].prevTime();
    const Real prev_pres_time = state[Press_Type].prevTime();
    Real dt                   = parent->dtLevel(level);
    //
    // Compute the u_mac for the correction.
    //
    mac_projector->mac_sync_solve(level,u_mac,dt,rho_half,fine_ratio);
    //
    // Update coarse grid state by adding correction from mac_sync solve
    // the correction is the advective tendency of the new velocities.
    //
    if (do_reflux)
    {
        MultiFab& S_new = get_new_data(State_Type);
        mac_projector->mac_sync_compute(level,u_mac,Vsync,Ssync, rho_half,
                                        level > 0 ? &getAdvFluxReg(level) : 0,
                                        is_conservative, prev_time,
                                        prev_pres_time,dt,
                                        NUM_STATE,be_cn_theta);
        //
        // The following used to be done in mac_sync_compute.
        //
        Ssync->mult(dt,Ssync->nGrow());
        //
        // Compute viscous sync.
        //
        if (is_diffusive[Xvel])
        {
            diffusion->diffuse_Vsync(Vsync, dt, be_cn_theta, rho_half, 1);
        }
        for (int sigma  = 0; sigma < numscal; sigma++)
        {
            int do_viscsyncflux = 1;
            int rho_flag        = !is_conservative[BL_SPACEDIM+sigma] ? 1 : 2;
            if (is_diffusive[BL_SPACEDIM+sigma])
                diffusion->diffuse_Ssync(Ssync, sigma, dt, be_cn_theta,
                                         rho_half, rho_flag, do_viscsyncflux);
        }
        //
        // Add the sync correction to the state.
        //
        for (int sigma  = 0; sigma < numscal; sigma++)
        {
            for (MultiFabIterator S_newmfi(S_new); S_newmfi.isValid(); ++S_newmfi)
            {
                DependentMultiFabIterator Ssyncmfi(S_newmfi, *Ssync);

                S_newmfi().plus(Ssyncmfi(),
                                S_newmfi.validbox(),
                                sigma,
                                BL_SPACEDIM+sigma,
                                1);
            }
        }
        //
        // Get boundary conditions.
        //
        Array<int*> sync_bc(grids.length());
        Array< Array<int> > sync_bc_array(grids.length());

        for (int i = 0; i < grids.length(); i++)
        {
            sync_bc_array[i] = getBCArray(State_Type,i,Density,numscal);
            sync_bc[i]       = sync_bc_array[i].dataPtr();
        }
        //
        // Interpolate the sync correction to the finer levels.
        //
        IntVect ratio = IntVect::TheUnitVector();
        const Real mult = 1.0;
        for (int lev = level+1; lev <= parent->finestLevel(); lev++)
        {
            ratio                 *= parent->refRatio(lev-1);
            NavierStokes& fine_lev = getLevel(lev);
            MultiFab& S_new        = fine_lev.get_new_data(State_Type);
            SyncInterp(*Ssync, level, S_new, lev, ratio,
                       0, BL_SPACEDIM, numscal, 1 , mult, sync_bc.dataPtr());
        }
    }
}

//
// The reflux function
//

void
NavierStokes::reflux ()
{
    if (level == parent->finestLevel())
        return;

    assert(do_reflux);

    MultiFab& S_crse = get_new_data(State_Type);
    //
    // First do refluxing step.
    //
    FluxRegister& fr_adv  = getAdvFluxReg(level+1);
    FluxRegister& fr_visc = getViscFluxReg(level+1);

    Real dt_crse = parent->dtLevel(level);
    Real scale = 1.0/dt_crse;
    //
    // DONT ACTUALLY REFLUX HERE, JUST FILL VSYNC AND SSYNC
    // IMPORTANT TO DO VISCOUS FIRST BECAUSE OF DENSITY-WEIGHTING OF V_SYNC
    //
    fr_visc.Reflux(*Vsync,volume,scale,0,0,BL_SPACEDIM,geom);
    fr_visc.Reflux(*Ssync,volume,scale,BL_SPACEDIM,0,NUM_STATE-BL_SPACEDIM,geom);

    for (MultiFabIterator rho_halfmfi(*rho_half); rho_halfmfi.isValid();
        ++rho_halfmfi)
    {
        DependentMultiFabIterator Vsyncmfi(rho_halfmfi, *Vsync);

        Vsyncmfi().divide(rho_halfmfi(),rho_halfmfi.validbox(),0,Xvel,1);
        Vsyncmfi().divide(rho_halfmfi(),rho_halfmfi.validbox(),0,Yvel,1);
#if (BL_SPACEDIM == 3)
        Vsyncmfi().divide(rho_halfmfi(),rho_halfmfi.validbox(),0,Zvel,1);
#endif
    }

    for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
    {
        if (!is_conservative[istate])
        {
            const int sigma = istate -  BL_SPACEDIM;
            for (MultiFabIterator rho_halfmfi(*rho_half);
                 rho_halfmfi.isValid(); ++rho_halfmfi)
            {
                DependentMultiFabIterator Ssyncmfi(rho_halfmfi, *Ssync);

                Ssyncmfi().divide(rho_halfmfi(),rho_halfmfi.validbox(),0,sigma,1);
            }
        }
    }

    fr_adv.Reflux(*Vsync,volume,scale,0,0,BL_SPACEDIM,geom);
    fr_adv.Reflux(*Ssync,volume,scale,BL_SPACEDIM,0,NUM_STATE-BL_SPACEDIM,geom);

    const BoxArray& fine_boxes = getLevel(level+1).boxArray();
    const int nfine            = fine_boxes.length();
    //
    // This is necessary in order to zero out the contribution to any
    // coarse grid cells which underlie fine grid cells.
    //
    for (int kf = 0; kf < nfine; kf++)
    {
        Box bf = ::coarsen(fine_boxes[kf],fine_ratio);

        for (MultiFabIterator Vsyncmfi(*Vsync); Vsyncmfi.isValid(); ++Vsyncmfi)
        {
            DependentMultiFabIterator Ssyncmfi(Vsyncmfi, *Ssync);
            assert(grids[Vsyncmfi.index()] == Vsyncmfi.validbox());

            if (bf.intersects(Vsyncmfi.validbox()))
            {
                Box bx = bf & Vsyncmfi.validbox();
                Vsyncmfi().setVal(0.0,bx,0,BL_SPACEDIM);
                Ssyncmfi().setVal(0.0,bx,0,NUM_STATE-BL_SPACEDIM);
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

    NavierStokes& fine_lev = getLevel(level+1);
    const BoxArray& fgrids = fine_lev.grids;
    MultiFab& fvolume      = fine_lev.volume;

    MultiFab& S_crse = get_new_data(State_Type);
    MultiFab& S_fine = fine_lev.get_new_data(State_Type);

    avgDown(grids, fgrids, S_crse, S_fine, volume, fvolume,
            level,  level+1, comp, 1, fine_ratio);
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
    //
    // Get the coarse intersection bounds.
    //
    const int* ovlo = ovlp.loVect();
    const int* ovhi = ovlp.hiVect();
    //
    // Get the fine and coarse pressures.
    //
    Real* cpres       = Pcrse.dataPtr();
    const int* clo    = Pcrse.loVect();
    const int* chi    = Pcrse.hiVect();
    const Real* fpres = Pfine.dataPtr();
    const int* flo    = Pfine.loVect();
    const int* fhi    = Pfine.hiVect();

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
    //
    // Get the coarse intersection bounds.
    //
    const int* ovlo = ovlp.loVect();
    const int* ovhi = ovlp.hiVect();
    //
    // Get the fine and coarse pressures.
    //
    Real* cpres       = Pcrse.dataPtr();
    const int* clo    = Pcrse.loVect();
    const int* chi    = Pcrse.hiVect();
    const Real* fpres = Pfine.dataPtr();
    const int* flo    = Pfine.loVect();
    const int* fhi    = Pfine.hiVect();

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

    NavierStokes& fine_lev = getLevel(level+1);
    const BoxArray& fgrids = fine_lev.grids;
    MultiFab& fvolume      = fine_lev.volume;
    //
    // Average down the states at the new time.
    //
    MultiFab& S_crse = get_new_data(State_Type);
    MultiFab& S_fine = fine_lev.get_new_data(State_Type);

    avgDown(grids, fgrids, S_crse, S_fine, volume, fvolume,
            level,  level+1, 0, S_crse.nComp(), fine_ratio);
    //
    // Now average down pressure over time n-(n+1) interval.
    //
    MultiFab& P_crse      = get_new_data(Press_Type);
    MultiFab& P_fine_init = fine_lev.get_new_data(Press_Type);
    MultiFab& P_fine_avg  = *fine_lev.p_avg;
    MultiFab& P_fine      = initial_step ? P_fine_init : P_fine_avg;

    const BoxArray& P_cgrids = state[Press_Type].boxArray();
    const BoxArray& P_fgrids = fine_lev.state[Press_Type].boxArray();

    MultiFabCopyDescriptor mfcd;

    MultiFabId mfidP_fine = mfcd.RegisterFabArray(&P_fine);

    vector<FillBoxId> fillBoxIdList;
    //
    // Inject fine pressure nodes down onto coarse nodes.
    //
    for (MultiFabIterator mfi(P_crse); mfi.isValid(); ++mfi)
    {
        assert(P_cgrids[mfi.index()] == mfi.validbox());
        //
        // Loop over fine grids and periodic extensions.
        //
        for (int fine = 0; fine < fgrids.length(); fine++)
        {
            Box ovlp = ::coarsen(P_fgrids[fine],fine_ratio) & mfi.validbox();

            if (ovlp.ok())
            {
                Box fine_ovlp = ::refine(ovlp,fine_ratio);

                fillBoxIdList.push_back(mfcd.AddBox(mfidP_fine,
                                                    fine_ovlp,
                                                    0,
                                                    fine,
                                                    0,
                                                    0,
                                                    P_fine.nComp()));

                assert(fillBoxIdList.back().box() == fine_ovlp);
                //
                // I need to save both the fine and the coarse grid indices.
                // I'll try to stuff'm into the single integer place available.
                // This assumes that integers are >= 4 bytes in size.
                //
                assert(sizeof(int) >= 4);
                assert(fine < 0xFFFF);
                assert(mfi.index() < 0xFFFF);

                fillBoxIdList.back().FabIndex((fine << 16) | mfi.index());
            }
        }
    }

    mfcd.CollectData();

    FArrayBox fine_fab;

    const int MyProc = ParallelDescriptor::MyProc();

    for (int i = 0; i < fillBoxIdList.size(); i++)
    {
        const FillBoxId& fbidFine = fillBoxIdList[i];

        int c_idx = fbidFine.FabIndex() & 0xFFFF;
        int f_idx = (fbidFine.FabIndex() >> 16) & 0xFFFF;

        assert(c_idx >= 0 && c_idx < P_cgrids.length());
        assert(f_idx >= 0 && f_idx < P_fgrids.length());

        Box ovlp = ::coarsen(P_fgrids[f_idx],fine_ratio) & P_cgrids[c_idx];

        assert(ovlp.ok());
        assert(P_crse.DistributionMap()[c_idx] == MyProc);
        //
        // Inject fine down to coarse.
        //
        fine_fab.resize(fbidFine.box(), P_fine.nComp());
        mfcd.FillFab(mfidP_fine, fbidFine, fine_fab);
        injectDown(ovlp, P_crse[c_idx], fine_fab, fine_ratio);
    }
    //
    // Next average down divu and dSdT at new time.
    //
    if (have_divu)
    {
        MultiFab& Divu_crse = get_new_data(Divu_Type);
        MultiFab& Divu_fine = fine_lev.get_new_data(Divu_Type);
        
        avgDown(grids, fgrids, Divu_crse, Divu_fine,
                volume, fvolume, level, level+1, 0, 1, fine_ratio );
    }
    if (have_dsdt)
    {
        MultiFab& Dsdt_crse = get_new_data(Dsdt_Type);
        MultiFab& Dsdt_fine = fine_lev.get_new_data(Dsdt_Type);
        
        avgDown(grids, fgrids, Dsdt_crse, Dsdt_fine,
                volume, fvolume, level, level+1, 0, 1, fine_ratio );
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
            fr.CrseInit(xflux,xflux.box(),0,0,start_ind,ncomp,-dt);
            fr.CrseInit(yflux,yflux.box(),1,0,start_ind,ncomp,-dt);
#if (BL_SPACEDIM == 3)                              
            fr.CrseInit(zflux,zflux.box(),2,0,start_ind,ncomp,-dt);
#endif
        }
        if (level > 0)
        {
            advflux_reg->FineAdd(xflux,0,i,0,start_ind,ncomp,dt);
            advflux_reg->FineAdd(yflux,1,i,0,start_ind,ncomp,dt);
#if (BL_SPACEDIM == 3)                                
            advflux_reg->FineAdd(zflux,2,i,0,start_ind,ncomp,dt);
#endif
        }
    }
}

//
// Virtual access function for getting the forcing terms for the
// velocities and scalars.  The base version computes a buoyancy.
//
// As NavierStokes is currently implmented.  Velocities are integrated
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

void
NavierStokes::getForce (FArrayBox& force,
                        int        gridno,
                        int        ngrow,
                        int        strt_comp,
                        int        num_comp,
                        Real       time)
{
    BoxLib::Error("NavierStokes::getForce(): not implemented");
#if 0
    Box bx(grids[gridno]);
    bx.grow(ngrow);
    force.resize(bx,num_comp);

    const int* lo = bx.loVect();
    const int* hi = bx.hiVect();

    Real grav = Abs(gravity);
    FArrayBox rho;
    for (int dc = 0; dc < num_comp; dc++)
    {
        const int sc = strt_comp + dc;
#if (BL_SPACEDIM == 2)
        if (BL_SPACEDIM == 2 && sc == Yvel && grav > 0.001) 
#endif
#if (BL_SPACEDIM == 3)
        if (BL_SPACEDIM == 3 && sc == Zvel && grav > 0.001) 
#endif
        {
            //
            // Set force to -rho*g.
            //
            getState(rho,gridno,ngrow,Density,1,time);
            rho.mult(-grav);
            force.copy(rho,0,dc,1);
        }
        else
        {
            force.setVal(0.0,dc);
        }
    }
#endif
}

void
NavierStokes::setForce (FArrayBox&       force,
                        int              gridno,
                        int              ngrow,
                        int              strt_comp,
                        int              num_comp,
                        const FArrayBox& Rho)
{
    assert(Rho.nComp() == 1);

    force.resize(::grow(grids[gridno],ngrow),num_comp);

    assert(Rho.box() == force.box());

    const Real grav = Abs(gravity);

    for (int dc = 0; dc < num_comp; dc++)
    {
        const int sc = strt_comp + dc;
#if (BL_SPACEDIM == 2)
        if (sc == Yvel && grav > 0.001) 
#endif
#if (BL_SPACEDIM == 3)
        if (sc == Zvel && grav > 0.001) 
#endif
        {
            //
            // Set force to -rho*g.
            //
            force.copy(Rho,0,dc,1);
            force.mult(-grav,dc,1);
        }
        else
        {
            force.setVal(0.0,dc);
        }
    }
}

//
// Fill patch and then call the other GradP function.
//

void
NavierStokes::getGradP (FArrayBox& gp,
                        int        gridno,
                        int        ngrow,
                        Real       time)
{
    BoxLib::Error("NavierStokes::getGradP(): not implemented");
#if 0
    Box gpbx(grids[gridno]);
    gpbx.grow(ngrow);
    Box p_box(surroundingNodes(gpbx));
    FArrayBox p_fab(p_box,1);
    FillPatch(p_fab,0,time,Press_Type,0,1);
    getGradP( p_fab, gp, grids[gridno], ngrow );
#endif
}

//
// Given a nodal pressure P compute the pressure gradient at the
// contained cell centers.
//
// This function should live in Projection, but it made life easier
// to have it here.
//

void
NavierStokes::getGradP (FArrayBox& p_fab,
                        FArrayBox& gp,
                        const Box& grd,
                        int        ngrow)
{
    //
    // Size the pressure gradient storage.
    //
    gp.resize(::grow(grd,ngrow),BL_SPACEDIM);
    //
    // Test to see if p_fab contains gp.box().
    //
    assert(::enclosedCells(p_fab.box()).contains(gp.box()));

    const int*  plo    = p_fab.loVect();
    const int*  phi    = p_fab.hiVect();
    const int*  glo    = gp.box().loVect();
    const int*  ghi    = gp.box().hiVect();
    const Real* p_dat  = p_fab.dataPtr();
    const Real* gp_dat = gp.dataPtr();
    const Real* dx     = geom.CellSize();

    FORT_GRADP(p_dat,ARLIM(plo),ARLIM(phi),
               gp_dat,ARLIM(glo),ARLIM(ghi),glo,ghi,dx);
}

//
// Fill patch divU.
//

void
NavierStokes::getDivCond (FArrayBox& fab,
                          int        gridno,
                          int        ngrow, 
                          Real       time)
{
    BoxLib::Error("NavierStokes::getDivCond(): not implemented");
#if 0
    getState(fab, gridno, ngrow, time, have_divu, Divu_Type,
             //divu_assoc, divu_unfilled);
             divu_unfilled);
#endif
}

MultiFab*
NavierStokes::getDivCond (int  ngrow, 
                          Real time)
{
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

void
NavierStokes::getDsdt (FArrayBox& fab,
                       int        gridno,
                       int        ngrow, 
                       Real       time)
{
    BoxLib::Error("NavierStokes::getDsdt(): not implemented");
#if 0
    getState(fab, gridno, ngrow, time, (have_dsdt && have_divu), Dsdt_Type,
             //dsdt_assoc, dsdt_unfilled);
             dsdt_unfilled);
#endif
}

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

void
NavierStokes::getState (FArrayBox& fab,
                        int        gridno,
                        int        ngrow,
                        int        strt_comp,
                        int        num_comp,
                        Real       time)
{
    BoxLib::Error("NavierStokes::getState(1): not implemented");

    getState(fab,gridno,ngrow,State_Type,strt_comp,num_comp,time);
}

//
// Fill patch a state component.
//

void
NavierStokes::getState (FArrayBox& fab,
                        int        gridno,
                        int        ngrow,
                        int        state_indx,
                        int        strt_comp,
                        int        num_comp, 
                        Real       time)
{
    BoxLib::Error("NavierStokes::getState(2): not implemented");

#if 0
    Box bx(grids[gridno]);
    bx.grow(ngrow);
    fab.resize(bx,num_comp);

    if (ngrow == 1 && !Geometry::isAnyPeriodic())
    {
        //hyp_assoc.setCacheWidth(ngrow);
        FillPatch(fab,0,time,state_indx,strt_comp,num_comp,
                 //hyp_assoc,gridno,cc1_unfilled[gridno]);
                 cc1_unfilled[gridno]);

    }
    else if (ngrow == HYP_GROW && !Geometry::isAnyPeriodic())
    {
        //hyp_assoc.setCacheWidth(ngrow);
        FillPatch(fab,0,time,state_indx,strt_comp,num_comp,
                 //hyp_assoc,gridno,hyp_unfilled[gridno]);
                 hyp_unfilled[gridno]);
    }
    else
    {
        FillPatch(fab,0,time,state_indx,strt_comp,num_comp);
    }
#endif
}

//
// Fill patch a state component.
//

MultiFab*
NavierStokes::getState (int  ngrow,
                        int  state_idx,
                        int  strt_comp,
                        int  num_comp, 
                        Real time)
{
    MultiFab* mf = new MultiFab(state[state_idx].boxArray(),num_comp,ngrow);

    FillPatchIterator fpi(*this,*mf,ngrow,time,state_idx,strt_comp,num_comp);

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
                              int  num_comp) 
{
    MultiFab& S = get_data(state_idx,time);

    if (S.nGrow() == 0)
        return;

    FillPatchIterator fpi(*this,S,S.nGrow(),time,state_idx,src_comp,num_comp);

    for ( ; fpi.isValid(); ++fpi)
    {
        //
        // Fill all ghost cells interior & exterior to valid region.
        //
        BoxList boxes = ::boxDiff(fpi().box(),grids[fpi.index()]);

        for (BoxListIterator bli(boxes); bli; ++bli)
        {
            S[fpi.index()].copy(fpi(),
                                bli(),
                                0,
                                bli(),
                                src_comp,
                                num_comp);
        }
    }
}

//
// Fill patch a state component.  This is a driver for calling
// the amrLevel filpatch function.  amrLevel filpatch works
// by first interpolating coarse data (time, then space), overwriting
// coarse data with fine data, and letting physical boundary conditions
// take precedence over interpolated ghost cells
//
// NOTE:: It might be a bug to resize the fab, since it could
// be part of a multifab with unforeseen consequences
//

void
NavierStokes::getState (FArrayBox&  fab,
                        int         gridno,
                        int         ngrow, 
                        Real        time,
                        int         have_state,
                        int         state_idx,
                        Array<Box>& unfilled,
                        int         strt_comp,
                        int         num_comp)
{
    BoxLib::Error("NavierStokes::getState(3): not implemented");

#if 0
    //
    // Create the storage.
    //
    Box bx(grids[gridno]);
    bx.grow(ngrow);
    fab.resize(bx,num_comp);
    if (fab.box()!=bx)
    {
      cout << "NavierStokes::getState : fab.box()!=bx\n";
      cout << "bx        = " << bx << NL;
      cout << "fab.box() = " << fab.box() << NL;
      BoxLib::Abort("NavierStokes::getState(...)");
    }
    
    if(!have_state)
    {
        fab.setVal(0.0);  // Can't filpatch what we don't have.
    }
    else
    {
        fab.setVal(1.0e30); // For debugging only.
        
        if (ngrow == 1 && !Geometry::isAnyPeriodic())
        {
            //assoc.setCacheWidth(ngrow);
            FillPatch(fab,0,time,state_idx,strt_comp,num_comp,
                     //assoc,gridno,unfilled[gridno]);
                     unfilled[gridno]);
        }
        else
        {
            FillPatch(fab,0,time,state_idx,strt_comp,num_comp);
        }
    }
#endif
}

//
// Fill patch divU.
//

void
NavierStokes::getDivCond (FArrayBox& fab,
                          int        ngrow,
                          Real       time)
{
    BoxLib::Error("NavierStokes::getDivCond(): not implemented");
#if 0
    //
    // For NavierStokes, defaults divU = 0.
    //
    if (!have_divu)
    {
        fab.setVal(0.0);
    }
    else
    {
        fab.setVal(1.0e30); // For debugging only.
        FillPatch(fab,0,time,Divu_Type,0,1);
    }
#endif
}

//
// Default divU is set to zero.
//

void
NavierStokes::calc_divu (Real      time,
                         Real      dt,
                         MultiFab& divu)
{
    if (have_divu)
    {
        divu.setVal(0.0);

        if (do_temp && visc_coef[Temp] > 0.0)
        {
            //
            // Compute Div(U) = Div(cond.Grad(T))/(rho.T) assuming cond = k/cp.
            //
            getViscTerms(divu,Temp,1,time);

            FillPatchIterator rho_fpi(*this,divu,0,time,State_Type,Density,1);
            FillPatchIterator temp_fpi(*this,divu,0,time,State_Type,Temp,1);

            for ( ;
                  rho_fpi.isValid() && temp_fpi.isValid();
                  ++rho_fpi, ++temp_fpi)
            {
                int i = rho_fpi.index();

                FArrayBox& rho  = rho_fpi();
                FArrayBox& temp = temp_fpi();

                divu[i].divide(rho,divu.box(i),0,0,1);
                divu[i].divide(temp,divu.box(i),0,0,1);
            }
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
    if (have_divu && have_dsdt)
    {
        dsdt.setVal(0.0);

        if (do_temp)
        {
            MultiFab& Divu_new = get_new_data(Divu_Type);
            MultiFab& Divu_old = get_old_data(Divu_Type);

            for (MultiFabIterator mfi(dsdt); mfi.isValid(); ++mfi)
            {
                DependentMultiFabIterator dmfi_old(mfi, Divu_old);
                DependentMultiFabIterator dmfi_new(mfi, Divu_new);

                mfi().copy(dmfi_new(),mfi.validbox(),0,mfi.validbox(),0,1);
                mfi().minus(dmfi_old(),mfi.validbox(),0,0,1);
                mfi().divide(dt,mfi.validbox(),0,1);
            }
        }
    }
}

void
NavierStokes::getViscTerms (MultiFab& visc_terms,
                            int       src_comp, 
                            int       num_comp,
                            Real      time)
{
    for (int icomp = src_comp; icomp < src_comp+num_comp; icomp++)
    {
        int rho_flag = !is_conservative[icomp] ? 1 : 2;

        diffusion->getViscTerms(visc_terms,src_comp,icomp,time,rho_flag);
    }
}

void 
NavierStokes::compute_grad_divu_minus_s (Real      time,
                                         MultiFab* grad_divu_minus_s,
                                         int       scaleRhoDivDt)
{
    Real dt = parent->dtLevel(0);
    grad_divu_minus_s->setVal(0.0);
#if (BL_SPACEDIM == 2)
    MultiFab& S_new = get_new_data(State_Type);
    S_new.FillBoundary();
    setPhysBoundaryValues(State_Type,Xvel,BL_SPACEDIM,time);

    BoxArray node_grids = grids;
    node_grids.surroundingNodes();
    MultiFab S(node_grids,1,0,Fab_allocate);
    MultiFab Dv(node_grids,1,0,Fab_allocate);

    int is_rz = (CoordSys::IsRZ() ? 1 : 0);

    int nghost = 0;
    if (have_divu)
    {
        projector->put_divu_in_node_rhs(S, level, nghost, node_grids, time);
    }
    else
    {
        S.setVal(0.0);
    }

    const int nGrow = 1;

    MultiFab U(grids, BL_SPACEDIM, nGrow, Fab_allocate);

    for (FillPatchIterator Ufpi(*this,U,nGrow,time,State_Type,0,BL_SPACEDIM);
         Ufpi.isValid();
         ++Ufpi)
    {
        U[Ufpi.index()].copy(Ufpi());
    }

    if (is_rz)
    {
        for (int n = 0; n < BL_SPACEDIM; n++) 
            projector->radMult(level,U,n);
    }

    const Real* dx = geom.CellSize();
    projector->computeDV(Dv, U, Xvel, dx, is_rz);

    if (is_rz)
    {
        for (int n = 0; n < BL_SPACEDIM; n++) 
            projector->radDiv(level,U,n);
    }
    //
    // Fix DV (mult by 2) at walls.
    // Outflow, Inflow, and Symmetry should be ok already.
    //
    const Box& domain  = geom.Domain();
    Box node_domain    = ::surroundingNodes(domain);
    const int* phys_lo = phys_bc.lo();
    const int* phys_hi = phys_bc.hi();
    const int* dlo     = node_domain.loVect();
    const int* dhi     = node_domain.hiVect();

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        for (MultiFabIterator Dvmfi(Dv); Dvmfi.isValid(); ++Dvmfi)
        {
            if (phys_lo[dir] == SlipWall || phys_lo[dir] == NoSlipWall)
            {
                Box domlo(node_domain);
                domlo.setRange(dir,dlo[dir],1);

                if (domlo.intersects(node_grids[Dvmfi.index()]))
                {
                    Box bx = domlo & node_grids[Dvmfi.index()];
                    Dvmfi().mult(2.0,bx,0,1);
                }
            }
            if (phys_hi[dir] == SlipWall || phys_hi[dir] == NoSlipWall)
            {
                Box domhi(node_domain);
                domhi.setRange(dir,dhi[dir],1);

                if (domhi.intersects(node_grids[Dvmfi.index()]))
                {
                    Box bx = domhi & node_grids[Dvmfi.index()];
                    Dvmfi().mult(2.0,bx,0,1);
                }
            }
        }
    }
    MultiFab divu_minus_s(node_grids,1,0,Fab_allocate);

    for (MultiFabIterator divu_minus_smfi(divu_minus_s);
        divu_minus_smfi.isValid(); ++divu_minus_smfi)
    {
        DependentMultiFabIterator Dvmfi(divu_minus_smfi, Dv);
        DependentMultiFabIterator Smfi(divu_minus_smfi, S);
        //
        // Compute nodal divU-S.
        //
        divu_minus_smfi().copy(Dvmfi());
        divu_minus_smfi().minus(Smfi());
    }

    if (is_rz)
    {
        Box dom_lo(node_domain);
        dom_lo.setRange(0,dlo[0],1);
        for (MultiFabIterator divu_minus_smfi(divu_minus_s);
             divu_minus_smfi.isValid(); ++divu_minus_smfi)
        {
            assert(grids[divu_minus_smfi.index()] == divu_minus_smfi.validbox());
            Box grid = divu_minus_smfi.validbox();
            grid.surroundingNodes();
            if (grid.intersects(dom_lo))
            {
                //
                // TODO -- use left_strip ???
                //
                Box left_strip = (grid & dom_lo);
                divu_minus_smfi().setVal(0.0,grid,0,1);
            }
        }
    }

    for (MultiFabIterator divu_minus_smfi(divu_minus_s);
         divu_minus_smfi.isValid(); ++divu_minus_smfi)
    {
        DependentMultiFabIterator grad_divu_minus_smfi(divu_minus_smfi,
                                                       (*grad_divu_minus_s));
        assert(grids[divu_minus_smfi.index()] == divu_minus_smfi.validbox());
        //
        // Compute gradient of this field.
        //
        getGradP(divu_minus_smfi(), grad_divu_minus_smfi(), 
                 divu_minus_smfi.validbox(), 0);
    }

    if (is_rz)
    {
        //
        // What we just computed was grad (r * (dU - S)).  We now compute
        // grad (dU - S) by
        // grad (dU - S) = 1/r grad (r * (dU - S)) - 1/r grad (r) * (dU-S) (*)
        //
        for (int n = 0; n < BL_SPACEDIM; n++) 
            projector->radDiv(level,*grad_divu_minus_s,n);

        MultiFab divu_minus_s_cc(grids,1,0,Fab_allocate);

        for (MultiFabIterator divu_minus_smfi(divu_minus_s);
             divu_minus_smfi.isValid(); ++divu_minus_smfi)
        {
            DependentMultiFabIterator divu_minus_s_ccmfi(divu_minus_smfi,
                                                         divu_minus_s_cc);
            DEF_CLIMITS(divu_minus_smfi(),nodedat,nodelo,nodehi);
            DEF_LIMITS(divu_minus_s_ccmfi(),ccdat,cclo,cchi);
            int rweighted = 1;
            FORT_HGN2C(&is_rz,&rweighted, ARLIM(nodelo), ARLIM(nodehi),
                       nodedat, ARLIM(cclo), ARLIM(cchi), cclo, cchi, ccdat);
        }
        //
        // The first divide is because divu_minus_s_cc contains r*(dU-S)
        // the second divide account for the 1/r in expression (*) above.
        //
        projector->radDiv(level,divu_minus_s_cc,0);
        projector->radDiv(level,divu_minus_s_cc,0);
        //
        // Since grad r = (1,0), we subtract  1/r * (dU-S from first component.
        //
        for (MultiFabIterator divu_minus_s_ccmfi(divu_minus_s_cc);
             divu_minus_s_ccmfi.isValid(); ++divu_minus_s_ccmfi)
        {
            DependentMultiFabIterator grad_divu_minus_smfi(divu_minus_s_ccmfi,
                                                           (*grad_divu_minus_s));
            grad_divu_minus_smfi().minus(divu_minus_s_ccmfi(),0,0,1);
        }
    }

    for (MultiFabIterator grad_divu_minus_smfi(*grad_divu_minus_s);
         grad_divu_minus_smfi.isValid(); ++grad_divu_minus_smfi)
    {
        DependentMultiFabIterator rho_halfmfi(grad_divu_minus_smfi, (*rho_half));
        assert(grids[grad_divu_minus_smfi.index()] == grad_divu_minus_smfi.validbox());
        if (scaleRhoDivDt)
        {
            grad_divu_minus_smfi().mult(divu_relax_factor*dx[0]*dx[1]/
                                        parent->dtLevel(0));
            for (int n = 0; n < BL_SPACEDIM; n++)
                grad_divu_minus_smfi().mult(rho_halfmfi(),
                                            grad_divu_minus_smfi.validbox(),
                                            0,
                                            n,
                                            1);
        }
        else
        {
            grad_divu_minus_smfi().mult(divu_relax_factor*dx[0]*dx[1]*dt/
                                        parent->dtLevel(0));
        }
    }
#endif /*BL_SPACEDIM==2*/
}
