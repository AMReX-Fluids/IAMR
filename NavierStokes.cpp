//
// $Id: NavierStokes.cpp,v 1.23 1997-12-11 23:30:22 lijewski Exp $
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
#include <Godunov.H>
#include <Interpolater.H>
#include <NavierStokes.H>
#include <MultiGrid.H>
#include <ArrayLim.H>
#include <NAVIERSTOKES_F.H>
#include <PROJECTION_F.H>
#include <PROB_F.H>

#ifdef BL_USE_NEW_HFILES
using std::streampos;
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

//=================================================================
// Initialization functions follow
//=================================================================

// -------------------------------------------------------------
// initialization of static members to default values
// -------------------------------------------------------------

// ----------------------- static objects
ErrorList   NavierStokes::err_list;
BCRec       NavierStokes::phys_bc;
Projection *NavierStokes::projector     = 0;
MacProj    *NavierStokes::mac_projector = 0;
Godunov    *NavierStokes::godunov       = 0;


// ----------------------- internal parameters
int  NavierStokes::verbose      = false;
Real NavierStokes::cfl          = 0.8;
Real NavierStokes::init_shrink  = 0.5;
Real NavierStokes::change_max   = 1.1;
Real NavierStokes::fixed_dt     = -1.0;
Real NavierStokes::dt_cutoff    = 0.0;
int  NavierStokes::init_iter    = 2;
Real NavierStokes::gravity      = 0;
int  NavierStokes::initial_step = false;
int  NavierStokes::initial_iter = false;
int  NavierStokes::radius_grow  = 1;
int  NavierStokes::sum_interval = -1;
int  NavierStokes::NUM_SCALARS  = 0;
int  NavierStokes::NUM_STATE    = 0;
Array<int> NavierStokes::is_conservative;


// ----------------------- viscosity parameters
Real NavierStokes::be_cn_theta  = 0.5;
Real NavierStokes::visc_tol     = 1.0e-10;  // tolerance for viscous solve
Real NavierStokes::visc_abs_tol = 1.0e-10;  // absolute tol. for visc solve
Array<int> NavierStokes::is_diffusive;
Array<Real> NavierStokes::visc_coef;


// ----------------------- internal switches
int  NavierStokes::do_sync_proj    = 1;
int  NavierStokes::do_MLsync_proj  = 0;
int  NavierStokes::do_reflux       = 1;
int  NavierStokes::do_mac_proj     = 1;
int  NavierStokes::do_diffusion    = 1;
     
// ------------------------ new members for non-zero divu
int  NavierStokes::additional_state_types_initialized = 0;
int  NavierStokes::Divu_Type = -1;
int  NavierStokes::Dsdt_Type = -1;
int  NavierStokes::have_divu = 0;
int  NavierStokes::have_dsdt = 0;
int  NavierStokes::S_in_vel_diffusion = 1;


Real NavierStokes::divu_minus_s_factor = 0.0;
Real NavierStokes::divu_relax_factor   = 0.0;
     
int  NavierStokes::num_state_type = 2;     // for backward compatibility


// -------------------------------------------------------------
void NavierStokes::variableCleanUp()
{
    derive_lst.clear();
    desc_lst.clear();
    delete projector;
    projector = 0;
    delete mac_projector;
    mac_projector = 0;
    delete godunov;
    godunov = 0;
}

// -------------------------------------------------------------
void NavierStokes::read_geometry()
{
#if (BL_SPACEDIM == 2)
    // must load coord here because CoordSys hasn't read it in yet
    ParmParse pp("geometry");
    int coord;
    pp.get("coord_sys",coord);
    if( (CoordSys::CoordType) coord == CoordSys::RZ &&
            phys_bc.lo(0) != Symmetry) {
      phys_bc.setLo(0,Symmetry);
      cout << "\nWarning:  Setting phys_bc at xlo to Symmetry\n\n";
    }
#endif
}

// -------------------------------------------------------------

void
NavierStokes::read_params()
{
    // read parameters from input file and command line

    ParmParse pp("ns");

    if(ParallelDescriptor::IOProcessor()) {
      verbose = pp.contains("v");
    }

    Array<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
    pp.getarr("lo_bc",lo_bc,0,BL_SPACEDIM);
    pp.getarr("hi_bc",hi_bc,0,BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; i++) {
        phys_bc.setLo(i,lo_bc[i]);
        phys_bc.setHi(i,hi_bc[i]);
    }
  
    read_geometry();

    // check phys_bc against possible periodic geometry
    // if periodic, must have internal BC marked
    if( Geometry::isAnyPeriodic() ){
      // do idiot check.  Periodic means interior in those directions
      for(int dir=0; dir<BL_SPACEDIM; dir++ ){
        if( Geometry::isPeriodic(dir) ){
          if( lo_bc[dir] != Interior ){
            cerr << "NavierStokes::variableSetUp:periodic in direction";
            cerr << dir << " but low BC is not Interior\n";
            ParallelDescriptor::Abort("NavierStokes::read_params()");
          }
          if( hi_bc[dir] != Interior ){
            cerr << "NavierStokes::variableSetUp:periodic in direction";
            cerr << dir << " but high BC is not Interior\n";
            ParallelDescriptor::Abort("NavierStokes::read_params()");
          }
        }
      }
    } else {
      // do idiot check.  If not periodic, should be no interior
      for(int dir=0; dir<BL_SPACEDIM; dir++ ){
        if( lo_bc[dir] == Interior ){
          cerr << "NavierStokes::variableSetUp:interior bc in direction";
          cerr << dir << " but no periodic\n";
          ParallelDescriptor::Abort("NavierStokes::read_params()");
        }
        if( hi_bc[dir] == Interior ){
          cerr << "NavierStokes::variableSetUp:interior bc in direction";
          cerr << dir << " but no periodic\n";
          ParallelDescriptor::Abort("NavierStokes::read_params()");
        }
      }
    }
    
    // get timestepping parameters
    pp.get("init_shrink",init_shrink);
    pp.get("cfl",cfl);
    pp.get("init_iter",init_iter);
    pp.get("dt_cutoff",dt_cutoff);
    pp.query("change_max",change_max);
    pp.query("fixed_dt",fixed_dt);
    pp.query("sum_interval",sum_interval);
    pp.get("gravity",gravity);

    // get run options
    int initial_do_sync_proj = do_sync_proj;
    pp.query("do_sync_proj",   do_sync_proj   );
    pp.query("do_MLsync_proj", do_MLsync_proj );
    pp.query("do_reflux",      do_reflux      );
    pp.query("do_diffusion",   do_diffusion   );

    // this test insures if the user toggles do_sync_proj,
    // the user has knowledge that do_MLsync_proj is meaningless
    if ( do_MLsync_proj && !do_sync_proj &&
         initial_do_sync_proj != do_sync_proj ) {
        cout << "mismatched options for NavierStokes\n";
        cout << "do_MLsync_proj and do_sync_proj are inconsistent\n";
        ParallelDescriptor::Abort("NavierStokes::read_params()");
    }


    // read viscous parameters and array of viscous coeficients.
    // NOTE: at this point, we dont know number of state variables
    //       so just read all values listed
    pp.query("visc_tol",visc_tol);
    pp.query("visc_abs_tol",visc_abs_tol);

    int n_visc = pp.countval("visc_coef");
    visc_coef.resize(n_visc);
    is_diffusive.resize(n_visc);
    pp.getarr("visc_coef",visc_coef,0,n_visc);

    pp.query("divu_minus_s_factor",divu_minus_s_factor);
    pp.query("divu_relax_factor",divu_relax_factor);
    pp.query("S_in_vel_diffusion",S_in_vel_diffusion);
    pp.query("be_cn_theta",be_cn_theta);
    if(be_cn_theta > 1.0 || be_cn_theta < .5) {
      cout << "NavierStokes::read_params: must have be_cn_theta <= 1.0 && >= .5"
           << NL;
      ParallelDescriptor::Abort("NavierStokes::read_params()");
    }

}

// -------------------------------------------------------------
NavierStokes::NavierStokes()
    : AmrLevel(), radius(PArrayManage)
{
    rho_avg = 0;
    rho_half = 0;
    p_avg = 0;
    Vsync = 0;
    Ssync = 0;
    sync_reg = 0;
    advflux_reg = 0;
    viscflux_reg = 0;
    u_mac = 0;
    aofs = 0;

    diffusion = 0;

    if(!additional_state_types_initialized)
      init_additional_state_types();
}

// -------------------------------------------------------------
NavierStokes::NavierStokes(Amr& papa, int lev, const Geometry &level_geom,
                           const BoxArray& bl, Real time)
        : AmrLevel(papa,lev,level_geom,bl,time),
          radius(PArrayManage)
{
    if(!additional_state_types_initialized)
        init_additional_state_types();
    
    const BoxArray& P_grids = state[Press_Type].boxArray();
    
    // alloc old_time pressure
    state[Press_Type].allocOldData();
    
    // alloc space for density and temporary pressure variables
    p_avg   = NULL;
    rho_avg = NULL;
    if (level > 0) {
        rho_avg = new MultiFab(grids,1,1,Fab_allocate);
        p_avg   = new MultiFab(P_grids,1,0,Fab_allocate);
    }
    rho_half = new MultiFab(grids,1,1,Fab_allocate);
    
    // build association lists
    //hyp_assoc.define(grids,HYP_GROW);
    //hyp_assoc.setCacheWidth(1);
    
    //if(have_divu) {
        //const BoxArray& Divu_grids = state[Divu_Type].boxArray();
        //divu_assoc.define(Divu_grids,DIVU_GROW); // needed
        //divu_assoc.setCacheWidth(1);
        //if(have_dsdt) {
            //const BoxArray& Dsdt_grids = state[Dsdt_Type].boxArray();
            //dsdt_assoc.define(Dsdt_grids,DSDT_GROW); // needed
            //dsdt_assoc.setCacheWidth(1);
        //}
    //}
    
    // pre-determine unfilled regions for fast filpatch
    buildUnfilledRegions();
    
    // build metric coeficients for RZ calculations
    buildMetrics();
    
    // build base state rho for atmospheric calculations
    // this is probably an Atmosphere virtual function
    buildRho0();

    // set up reflux registers
    sync_reg = 0;
    if (level > 0 && do_sync_proj) {
        sync_reg = new SyncRegister(grids,crse_ratio,level);
    }

    advflux_reg = 0;
    viscflux_reg = 0;
    if (level > 0 && do_reflux) {
        advflux_reg  = new FluxRegister(grids,crse_ratio,level,NUM_STATE);
        viscflux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);
    }

    // initialize work multifabs
    Vsync = 0;
    Ssync = 0;
    u_mac = 0;
    aofs  = 0;

    // set up the level projector
    if ( do_MLsync_proj || do_sync_proj )
    {
        if (projector == 0)
        {
            projector = new Projection(parent,&phys_bc,
                                     do_sync_proj,parent->finestLevel(),
                                     radius_grow );
        }
        projector->install_level(level, this, &radius );
    }

    // set up the godunov box
    SetGodunov();
    
    // set up diffusion
    diffusion = new Diffusion(parent, this,
                              (level > 0) ? getLevel(level-1).diffusion
                              : (Diffusion*)NULL,
                              NUM_STATE, viscflux_reg, volume, area,
                              is_diffusive, visc_coef);

    // set up the mac projector
    if (mac_projector == 0) {
        int finest_level = parent->finestLevel();
        mac_projector = new MacProj( parent, finest_level,
                                     &phys_bc, radius_grow );
    }
    mac_projector->install_level(level, this,
                                 volume, area, &radius );
}

// -------------------------------------------------------------
NavierStokes::~NavierStokes()
{
    if (rho_avg != NULL) delete rho_avg;
    if (p_avg != NULL)   delete p_avg;
    delete rho_half;
    delete Vsync;
    delete Ssync;
    delete sync_reg;
    delete advflux_reg;
    delete viscflux_reg;
    delete [] u_mac;
    int ngrids = grids.length();
    
    // area is not PArrayManaged, must manually delete
    //for (int k = 0; k < ngrids; k++) {
        for (int dir = 0; dir < BL_SPACEDIM; dir++) {
            //FArrayBox *a = area[dir].remove(k);
            //delete a;
            area[dir].clearUnmanaged();
        }
    //}
    
    if (mac_projector != NULL)
        mac_projector->cleanup(level);
    delete diffusion;
}

// -------------------------------------------------------------
void NavierStokes::init_additional_state_types()
{

    additional_state_types_initialized = 1;

    int _Divu = -1;
    int dummy_Divu_Type;
    have_divu = 0;
    have_divu = isStateVariable("divu", dummy_Divu_Type, _Divu);
    have_divu = have_divu && dummy_Divu_Type==Divu_Type;
    if(ParallelDescriptor::IOProcessor()) {
      cout << "NavierStokes::init_additional_state_types()::have_divu = "
           << have_divu << NL;
    }
    if(have_divu && _Divu!=Divu) {
      if(ParallelDescriptor::IOProcessor()) {
        cout << "NavierStokes::init_additional_state_types(): divu must be " <<
                "0-th, Divu_Type component in the state\n";
      }
      ParallelDescriptor::Abort("NavierStokes::init_additional_state_types()");
    }

    int _Dsdt = -1;
    int dummy_Dsdt_Type;
    have_dsdt = 0;
    have_dsdt = isStateVariable("dsdt", dummy_Dsdt_Type, _Dsdt);
    have_dsdt = have_dsdt && dummy_Dsdt_Type==Dsdt_Type;
    if(ParallelDescriptor::IOProcessor()) {
      cout << "NavierStokes::init_additional_state_types()::have_dsdt = "
           << have_dsdt << NL;
    }
    if(have_dsdt && _Dsdt!=Dsdt) {
      if(ParallelDescriptor::IOProcessor()) {
        cout << "NavierStokes::init_additional_state_types(): dsdt must be " <<
                "0-th, Dsdt_Type component in the state\n";
      }
      ParallelDescriptor::Abort("NavierStokes::init_additional_state_types()");
    }
    if(have_dsdt && !have_divu) {
      if(ParallelDescriptor::IOProcessor()) {
        cout << "NavierStokes::init_additional_state_types(): "
             << "must have divu in order to have dsdt\n";
      }
      ParallelDescriptor::Abort("NavierStokes::init_additional_state_types()");
    }

    num_state_type = desc_lst.length();
    if(ParallelDescriptor::IOProcessor()) {
      cout << "NavierStokes::init_additional_state_types: num_state_type = "
           << num_state_type << NL;
    }


}

// -------------------------------------------------------------
void NavierStokes::SaveOldBoundary( Real time )
{
    MultiFab &Sold = get_old_data(State_Type);

#if (USEOLDFILLPATCH == 1)
    FArrayBox temp;
    for(int i = 0; i < grids.length(); i++) {
        getState(temp,i,1,0,NUM_STATE,time);
        Sold[i].copy(temp,0,0,NUM_STATE);
    }
#else
    int nGrow    = 1;
    int srcComp  = 0;
    int destComp = 0;
    int nComp    = NUM_STATE;
    for(FillPatchIterator Soldfpi(*this, Sold, nGrow,
                                  destComp, time, State_Type,
                                  srcComp, nComp);
        Soldfpi.isValid(); ++Soldfpi)
    {
      DependentMultiFabIterator Soldmfi(Soldfpi, Sold);
      Soldmfi().copy(Soldfpi(), srcComp, destComp, nComp);
    }
#endif
}

// -------------------------------------------------------------
void NavierStokes::SaveNewBoundary( Real time )
{
    MultiFab &Snew = get_new_data(State_Type);

#if (USEOLDFILLPATCH == 1)
    FArrayBox temp;
    for(int i = 0; i < grids.length(); i++) {
        getState(temp,i,1,0,NUM_STATE,time);
        Snew[i].copy(temp,0,0,NUM_STATE);
    }
#else
    int nGrow    = 1;
    int srcComp  = 0;
    int destComp = 0;
    int nComp    = NUM_STATE;
    for(FillPatchIterator Snewfpi(*this, Snew, nGrow,
                                  destComp, time, State_Type,
                                  srcComp, nComp);
        Snewfpi.isValid(); ++Snewfpi)
    {
      DependentMultiFabIterator Snewmfi(Snewfpi, Snew);
      Snewmfi().copy(Snewfpi(), srcComp, destComp, nComp);
    }
#endif
}


// -------------------------------------------------------------
// since the pressure solver always stores its estimate of the
// pressure solver in Pnew, we need to copy it to Pold at the start
// -------------------------------------------------------------

void NavierStokes::initOldPress()
{
    MultiFab &P_new = get_new_data(Press_Type);
    MultiFab &P_old = get_old_data(Press_Type);
    //int i;
    //for (i = 0; i < grids.length(); i++) {
        //P_old[i].copy(P_new[i]);
    //}
    P_old.copy(P_new);
}

// -------------------------------------------------------------
void NavierStokes::zeroNewPress()
{
    MultiFab &P_new = get_new_data(Press_Type);
    //int i;
    //for (i = 0; i < grids.length(); i++) {
        //P_new[i].setVal(0.0);
    //}
    P_new.setVal(0.0);
}

// -------------------------------------------------------------
void NavierStokes::zeroOldPress()
{
    MultiFab &P_old = get_old_data(Press_Type);
    //int i;
    //for (i = 0; i < grids.length(); i++) {
        //P_old[i].setVal(0.0);
    //}
    P_old.setVal(0.0);
}

// -------------------------------------------------------------
void NavierStokes::allocOldData()
{
    int init_pres = !(state[Press_Type].hasOldData());
    for (int k = 0; k < num_state_type; k++) state[k].allocOldData();
    if (init_pres) initOldPress();
}

// -------------------------------------------------------------
void NavierStokes::removeOldData()
{
    AmrLevel::removeOldData();
}


// create the godunov object
void NavierStokes::SetGodunov()
{
    if (godunov == 0) {
        godunov = new Godunov();
    }
}


// -------------------------------------------------------------
void
NavierStokes::restart(Amr& papa, istream& is)
{
    AmrLevel::restart(papa,is);

    if ( do_MLsync_proj || do_sync_proj )
    {
        if (projector == 0) {
          projector = new Projection(parent,&phys_bc,
                                     do_sync_proj,parent->finestLevel(),
                                     radius_grow );
        }
        projector->install_level(level, this, &radius );
    }

    // set the godunov box
    SetGodunov();
    
    if (mac_projector == 0) {
        int finest_level = parent->finestLevel();
        mac_projector = new MacProj( parent, finest_level,
                                     &phys_bc, radius_grow );
    }
    mac_projector->install_level(level, this,
                                 volume, area, &radius );

    rho_avg = NULL;
    p_avg = NULL;
    const BoxArray& P_grids = state[Press_Type].boxArray();

    // alloc space for density and temporary pressure variables
    if (level > 0) {
        rho_avg = new MultiFab(grids,1,1,Fab_allocate);
        p_avg   = new MultiFab(P_grids,1,0,Fab_allocate);
    }
    rho_half = new MultiFab(grids,1,1,Fab_allocate);

      // build association lists
    //hyp_assoc.define(grids,HYP_GROW);
    //hyp_assoc.setCacheWidth(1);

    //if(have_divu) {
      //const BoxArray& Divu_grids = state[Divu_Type].boxArray();
      //divu_assoc.define(Divu_grids,DIVU_GROW);
      //divu_assoc.setCacheWidth(1);
      //if(have_dsdt) {
        //const BoxArray& Dsdt_grids = state[Dsdt_Type].boxArray();
        //dsdt_assoc.define(Dsdt_grids,DSDT_GROW);
        //dsdt_assoc.setCacheWidth(1);
      //}
    //}

      // pre-determine unfilled regions for fast filpatch
    buildUnfilledRegions();

      // build metric coeficients for RZ calculations
    buildMetrics();

    // build base state rho for atmospheric calculations
    // this is probably an Atmosphere virtual function
    buildRho0();

    assert(sync_reg == 0);
    if (level > 0 && do_sync_proj) {
        sync_reg = new SyncRegister(grids,crse_ratio,level);
    }
    assert(advflux_reg == 0);
    if (level > 0 && do_reflux) {
        advflux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);
    }
    assert(viscflux_reg == 0);
    if (level > 0 && do_reflux) {
        viscflux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);
    }

    assert(Vsync == 0);
    assert(Ssync == 0);
    if (level < parent->finestLevel()) {
        Vsync = new MultiFab;
        Vsync->define(grids,BL_SPACEDIM,1,Fab_allocate);
        Ssync = new MultiFab;
        Ssync->define(grids,NUM_STATE-BL_SPACEDIM,1,Fab_allocate);
    }

    diffusion = new Diffusion(parent, this,
                              (level > 0) ? getLevel(level-1).diffusion
                                          : (Diffusion*)NULL,
                              NUM_STATE, viscflux_reg, volume, area,
                              is_diffusive, visc_coef);

}

// -------------------------------------------------------------
// build rho0 arrays as ones for the base class
// empty function for now
void NavierStokes::buildRho0()
{
}



// -------------------------------------------------------------
void
NavierStokes::buildMetrics()
{
    int ngrids = grids.length();
    radius.resize(ngrids);
    const Real* dx = geom.CellSize();
    Real dxr = dx[0];
    int i;
    for (i = 0; i < ngrids; i++) {
        const Box& b = grids[i];
        int ilo = b.smallEnd(0)-radius_grow;
        int ihi = b.bigEnd(0)+radius_grow;
        int len = ihi - ilo + 1;
        Real *rad = new Real[len];
        radius.set(i,rad);
        if (CoordSys::IsCartesian()) {
            for (int j = 0; j < len; j++) rad[j] = 1.0;
        } else {
            Real xlo = grid_loc[i].lo(0) + (0.5 - radius_grow)*dxr;
            for (int j = 0; j < len; j++) rad[j] = xlo + j*dxr;
        }
    }

    // build volume and face area arrays
    geom.GetVolume(volume,grids,GEOM_GROW);
    for (int dir = 0; dir < BL_SPACEDIM; dir++) {
        geom.GetFaceArea(area[dir],grids,dir,GEOM_GROW);
    }
}



// for fast filpatch, we precompute the subbox of a grid that
// must be interpolated from coarser data.  This routine precomputes
// these boxes and stores them in (*_unfilled) arrays
// 
void NavierStokes::buildUnfilledRegions()
{
    if (hyp_unfilled.ready()) {
        BoxLib::Error("buildUnfilledRegions: already built");
    }

    int numgrids = grids.length();
    hyp_unfilled.resize(numgrids);
    cc1_unfilled.resize(numgrids);

    int i, j;

    const Box& domain = geom.Domain();
    for (i = 0; i < numgrids; i++) {
        BoxDomain fd;
        Box sbox(grow(grids[i],HYP_GROW));
        sbox &= domain;
        fd.add(sbox);
        for (j = 0; j < grids.length(); j++) {
            const Box& gbox = grids[j];
            if (sbox.intersects(gbox)) fd.rmBox(gbox);
        }
        hyp_unfilled.set(i,fd.minimalBox());
        fd.clear();
    }

    for (i = 0; i < numgrids; i++) {
        BoxDomain fd;
        Box sbox(grow(grids[i],1));
        sbox &= domain;
        fd.add(sbox);
        for (j = 0; j < grids.length(); j++) {
            const Box& gbox = grids[j];
            if (sbox.intersects(gbox)) fd.rmBox(gbox);
        }
        cc1_unfilled.set(i,fd.minimalBox());
        fd.clear();
    }

    if(have_divu) {
      divu_unfilled.resize(numgrids);
      for (i = 0; i < numgrids; i++) {
        const Box& Divu_domain = state[Divu_Type].getDomain();
        IndexType divu_typ = Divu_domain.ixType();
        BoxDomain fd(divu_typ);
        const BoxArray& Divu_grids = state[Divu_Type].boxArray();
        Box sbox(grow(Divu_grids[i],DIVU_GROW));
        sbox &= Divu_domain;
        fd.add(sbox);
        for (j = 0; j < grids.length(); j++) {
            const Box& gbox = Divu_grids[j];
            if (sbox.intersects(gbox)) fd.rmBox(gbox);
        }
        divu_unfilled.set(i,fd.minimalBox());
        fd.clear();
      }
      if(have_dsdt) {
        dsdt_unfilled.resize(numgrids);
        for (i = 0; i < numgrids; i++) {
          const Box& Dsdt_domain = state[Dsdt_Type].getDomain();
          IndexType dsdt_typ = Dsdt_domain.ixType();
          BoxDomain fd(dsdt_typ);
          const BoxArray& Dsdt_grids = state[Dsdt_Type].boxArray();
          Box sbox(grow(Dsdt_grids[i],DSDT_GROW));
          sbox &= Dsdt_domain;
          fd.add(sbox);
          for (j = 0; j < grids.length(); j++) {
            const Box& gbox = Dsdt_grids[j];
            if (sbox.intersects(gbox)) fd.rmBox(gbox);
          }
          dsdt_unfilled.set(i,fd.minimalBox());
          fd.clear();
        }
      }
    }

    for (i = 0; i < numgrids; i++) {
        const Box& P_domain = state[Press_Type].getDomain();
        IndexType p_typ = P_domain.ixType();
        BoxDomain fd(p_typ);
        const BoxArray& P_grids = state[Press_Type].boxArray();
        Box sbox(grow(P_grids[i],PRESS_GROW));
        sbox &= P_domain;
        fd.add(sbox);
        for (j = 0; j < P_grids.length(); j++) {
            const Box& gbox = P_grids[j];
            if (sbox.intersects(gbox)) fd.rmBox(gbox);
        }
        fd.clear();
    }
}


// reset the time levels to time (time) and timestep dt.
// This is done at the start of the timestep in the pressure
// iteration section
void NavierStokes::resetState(Real time, Real dt_old, Real dt_new)
{
    // reset state and pressure types
    state[State_Type].reset();
    state[State_Type].setTimeLevel(time,dt_old,dt_new);
    initOldPress();
    state[Press_Type].setTimeLevel(time-dt_old,dt_old,dt_old);

    // reset state types for divu not equal to zero
    if (have_divu) {
        state[Divu_Type].reset();
        state[Divu_Type].setTimeLevel(time,dt_old,dt_new);
        if ( have_dsdt ) {
            state[Dsdt_Type].reset();
            state[Dsdt_Type].setTimeLevel(time,dt_old,dt_new);
        }
    }
}

// set the time levels to time (time) and timestep dt
void NavierStokes::setTimeLevel(Real time, Real dt_old, Real dt_new)
{
    state[State_Type].setTimeLevel(time,dt_old,dt_new);
    if(have_divu) {
      state[Divu_Type].setTimeLevel(time,dt_old,dt_new);
      if(have_dsdt) {
        state[Dsdt_Type].setTimeLevel(time,dt_old,dt_new);
      }
    }
    state[Press_Type].setTimeLevel(time-dt_old,dt_old,dt_old);
}



// This function initializes the State and Pressure with data
void NavierStokes::initData()
{
    // initialize the state and the pressure
    int ns = NUM_STATE - BL_SPACEDIM;
    const Real* dx = geom.CellSize();
    MultiFab &S_new = get_new_data(State_Type);
    MultiFab &P_new = get_new_data(Press_Type);
    Real cur_time = state[State_Type].curTime();
    //for (int i = 0; i < grids.length(); i++) {
    for(MultiFabIterator snewmfi(S_new); snewmfi.isValid(); ++snewmfi) {
      DependentMultiFabIterator pnewmfi(snewmfi, P_new);
      assert(grids[snewmfi.index()] == snewmfi.validbox());
        const int* lo = snewmfi.validbox().loVect();
        const int* hi = snewmfi.validbox().hiVect();
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

    // initialize other types
    initDataOtherTypes();

    // initialize divU and dSdt
    if(have_divu) {
      Real cur_time = state[Divu_Type].curTime();
      MultiFab &Divu_new = get_new_data(Divu_Type);
      Real dt = 1.0;
      state[State_Type].setTimeLevel(cur_time,dt,dt);
      Real dtin = -1.0; // dummy value denotes initialization
      calc_divu(cur_time, dtin, Divu_new);
      if(have_dsdt) {
        MultiFab &Dsdt_new = get_new_data(Dsdt_Type);
        Dsdt_new.setVal(0.0);
      }
    }
}



// This function fills a new level n with the best
// level n and coarser data available
void NavierStokes::init( AmrLevel &old)
{
    // get data pointers
    NavierStokes *oldns = (NavierStokes*) &old;
    MultiFab &S_new = get_new_data(State_Type);
    MultiFab &P_new = get_new_data(Press_Type);
    MultiFab &P_old = get_old_data(Press_Type);

    // get time information
    Real dt_new = parent->dtLevel(level);
    Real cur_time  = oldns->state[State_Type].curTime();
    Real prev_time = oldns->state[State_Type].prevTime();
    Real dt_old = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);
    Real cur_pres_time = state[Press_Type].curTime();

#if (USEOLDFILLPATCH == 1)
    int i;
    if(ParallelDescriptor::NProcs() > 1) {
      cerr << "Using old FillPatch with multiple processors.\n";
      ParallelDescriptor::Abort("NavierStokes::init( AmrLevel &old)");
    } else {
      cerr << "Using old FillPatch.\n";
    }
    // get best state and pressure data
    for (i = 0; i < grids.length(); i++) {
        const int* lo = grids[i].loVect();
        const int* hi = grids[i].hiVect();
        oldns->FillPatch(S_new[i],0,cur_time,State_Type,0,NUM_STATE);
        oldns->FillPatch(P_new[i],0,cur_pres_time,Press_Type,0,1);
        P_old[i].copy(P_new[i]);
    }

    // get best divu and dSdt data
    if (have_divu) {
      MultiFab &Divu_new = get_new_data(Divu_Type);
      for (i = 0; i < grids.length(); i++) {
        const int* lo = grids[i].loVect();
        const int* hi = grids[i].hiVect();
        oldns->FillPatch(Divu_new[i],0,cur_time,Divu_Type,0,1);
      }
      if (have_dsdt) {
        MultiFab &Dsdt_new = get_new_data(Dsdt_Type);
        for (i = 0; i < grids.length(); i++) {
          const int* lo = grids[i].loVect();
          const int* hi = grids[i].hiVect();
          oldns->FillPatch(Dsdt_new[i],0,cur_time,Dsdt_Type,0,1);
        }
      }
    }
#else
    int destComp = 0;
    int srcComp  = 0;
    int boxGrow  = 0;
    for(FillPatchIterator snewfpi(old, S_new, boxGrow, destComp, cur_time,
                              State_Type, srcComp, NUM_STATE);
        snewfpi.isValid();
        ++snewfpi)
    {
      DependentMultiFabIterator dsnewmfi(snewfpi, S_new);
      dsnewmfi().copy(snewfpi());
    }

    int numComp = 1;
    for(FillPatchIterator pnewfpi(old, P_new, boxGrow, destComp, cur_pres_time,
                              Press_Type, srcComp, numComp);
        pnewfpi.isValid();
        ++pnewfpi)
    {
      DependentMultiFabIterator dpnewmfi(pnewfpi, P_new);
      DependentMultiFabIterator dpoldmfi(pnewfpi, P_old);
      dpnewmfi().copy(pnewfpi());  // copy from the temp fill patched fab into new
      dpoldmfi().copy(pnewfpi());  // copy from the temp fill patched fab into old
    }

    // get best divu and dSdt data
    if (have_divu) {
      MultiFab &Divu_new = get_new_data(Divu_Type);
      int numComp = 1;
      for(FillPatchIterator divunewfpi(old, Divu_new, boxGrow, destComp, cur_time,
                                    Divu_Type, srcComp, numComp);
          divunewfpi.isValid();
          ++divunewfpi)
      {
        DependentMultiFabIterator ddivunewmfi(divunewfpi, Divu_new);
        ddivunewmfi().copy(divunewfpi());
      }

      if (have_dsdt) {
        MultiFab &Dsdt_new = get_new_data(Dsdt_Type);
        int numComp = 1;
        for(FillPatchIterator dsdtnewfpi(old, Dsdt_new, boxGrow, destComp, cur_time,
                                      Dsdt_Type, srcComp, numComp);
            dsdtnewfpi.isValid();
            ++dsdtnewfpi)
        {
          DependentMultiFabIterator ddsdtnewmfi(dsdtnewfpi, Dsdt_new);
          ddsdtnewmfi().copy(dsdtnewfpi());
        }
      }
    }
#endif
}


// this version inits the data on a new level that did not
// exist before regridding
void NavierStokes::init()
{
    // get data pointers
    MultiFab &S_new = get_new_data(State_Type);
    MultiFab &P_new = get_new_data(Press_Type);
    MultiFab &P_old = get_old_data(Press_Type);

    // get time information
    assert(level > 0);
    const Array<Real>& dt_amr = parent->dtLevel();
    Array<Real> dt_new(level+1);
    for (int lev = 0; lev < level; lev++) {
        dt_new[lev] = dt_amr[lev];
    }
    Real dt = dt_new[level-1]/(Real)parent->MaxRefRatio(level-1);
    dt_new[level] = dt;
    parent->setDtLevel(dt_new);

    NavierStokes& old = getLevel(level-1);
    Real cur_time  = old.state[State_Type].curTime();
    Real prev_time = old.state[State_Type].prevTime();
    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);

    Real cur_pres_time = state[Press_Type].curTime();
    
    // get best coarse state and pressure data
    int i;
    for (i = 0; i < grids.length(); i++) {
        FillCoarsePatch(S_new[i],0,cur_time,State_Type,0,NUM_STATE);
        FillCoarsePatch(P_new[i],0,cur_pres_time,Press_Type,0,1);
        P_old[i].copy(P_new[i]);
    }

    
    // get best coarse divU and dSdt data
    if (have_divu) {
      MultiFab &Divu_new = get_new_data(Divu_Type);
      for (i = 0; i < grids.length(); i++) {
        FillCoarsePatch(Divu_new[i],0,cur_time,Divu_Type,0,1);
      }
      if (have_dsdt) {
        MultiFab &Dsdt_new = get_new_data(Dsdt_Type);
        for (i = 0; i < grids.length(); i++) {
          FillCoarsePatch(Dsdt_new[i],0,cur_time,Dsdt_Type,0,1);
        }
      }
    }
}




//=================================================================
// ADVANCE FUNCTIONS
//=================================================================

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
// the boundaries of the fine level BoxArray
//----------------------------------------------------------------

// setup for the advance function
void NavierStokes::advance_setup( Real time, Real dt, int iteration, int ncycle)
{
    int finest_level = parent->finestLevel();
    mac_projector->setup(level);

    // why are they defined here versus the constructor?
    if (level < finest_level) {
        if (Vsync == 0) {
            //Vsync = new MultiFab;
            //Vsync->define(grids,BL_SPACEDIM,1,Fab_allocate);
            Vsync = new MultiFab(grids,BL_SPACEDIM,1,Fab_allocate);
        }
        if (Ssync == 0) {
            //Ssync = new MultiFab;
            //Ssync->define(grids,NUM_STATE-BL_SPACEDIM,1,Fab_allocate);
            Ssync = new MultiFab(grids,NUM_STATE-BL_SPACEDIM,1,Fab_allocate);
        }
        Vsync->setVal(0.0);
        Ssync->setVal(0.0);
    }

    // set reflux registers to zero
    if (do_reflux && level < finest_level) {
        FluxRegister& fr_adv  = getAdvFluxReg(level+1);
        fr_adv.setVal(0.0);
        FluxRegister& fr_visc = getViscFluxReg(level+1);
        fr_visc.setVal(0.0);
    }
    
    // alloc space for edge velocities (normal comp only)
    if (u_mac == 0) {
        u_mac = new MultiFab[BL_SPACEDIM];
        for (int dir = 0; dir < BL_SPACEDIM; dir++) {
            BoxArray edge_grids(grids);
            edge_grids.surroundingNodes(dir);
            u_mac[dir].define(edge_grids,1,0,Fab_allocate);
        }
    }

    // alloc multifab to hold advective tendencies
    assert( aofs == 0 );
    if ( do_diffusion ) {
        aofs = new MultiFab(grids,NUM_STATE,0,Fab_allocate);
    }
    
    // set rho_avg
    if ( !initial_step && (level > 0) && (iteration == 1) ) {
        Real fratio = (Real) ncycle;
        Real alpha = 0.5/fratio;
        initRhoAvg(alpha);
    }

    // set up state multifabs for the advance
    for ( int k = 0; k < num_state_type; k++) {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }
    MultiFab &temp = get_new_data(State_Type);
    temp.setVal(bogus_value);

    if( (level>0 || geom.isAnyPeriodic()) && do_diffusion ) {
        // this is neccessary so that diffusion works properly during the first
        // time step
        for ( int k = 0; k < num_state_type; k++) {
            if(k!=Press_Type) {
                MultiFab &new_state = get_new_data(k);
                MultiFab &old_state = get_old_data(k);
                FillStateBndry(time,k,0,old_state.nComp());
                FillStateBndry(time+dt,k,0,new_state.nComp());
            }
        }
        // --------------------- DS --------------
        // SaveNewBoundary( time+dt );
        // SaveOldBoundary( time    );
        // --------------------- DS --------------
    }

}


//clean up after the advance function
void NavierStokes::advance_cleanup(Real dt, int iteration, int ncycle)
{
    
    int finest_level = parent->finestLevel();
    
    // delete mac edge velocities
    if (level == finest_level) {
        delete [] u_mac;
        u_mac = 0;
    }
    
    // delete advective tendencies
    delete aofs;
    aofs = 0;
}


// -------------------------------------------------------------
// compute a timestep at a level
// return largest safe timestep that can be used
// -------------------------------------------------------------
Real NavierStokes::advance(Real time, Real dt, int iteration, int ncycle)
{
    // ------------------ Advance setup
    if (verbose) {
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
    
    SaveOldBoundary( time ); 
    MultiFab &Snew = get_new_data(State_Type);
    MultiFab &Sold = get_old_data(State_Type);

    RunStats vel_pred_stats("vel_predict", level);
    RunStats  vel_adv_stats("vel_advect" , level);
    RunStats scal_adv_stats("scal_advect", level);
    RunStats  vel_upd_stats("vel_update" , level);
    RunStats scal_upd_stats("scal_update", level);
    RunStats      mac_stats("mac_project", level);

    // ------------------ Advance starts here
    // compute traced states for normal comp of velocity
    // at half time level
    vel_pred_stats.start();
    Real dt_test = 0.0, dummy = 0.0;
    dt_test = predict_velocity(dt,dummy);
    vel_pred_stats.end();

    // do MAC projection and update edge velocities
    mac_stats.start();
    if (verbose) {
        cout << "... mac_projection\n";
    }

    MultiFab divu(grids,1,0,Fab_allocate);
    MultiFab dsdt(grids,1,0,Fab_allocate);
#if (USEOLDFILLPATCH == 1)
    for(MultiFabIterator divumfi(divu); divumfi.isValid(); ++divumfi) {
      DependentMultiFabIterator dsdtmfi(divumfi, dsdt);
        int k = divumfi.index();
        getDivCond(divumfi(),k,0,time);
        getDsdt(dsdtmfi(),k,0,time);
        dsdtmfi().mult(.5*dt);
        divumfi().plus(dsdtmfi());
    }
    
#else

    int destComp = 0;
    int srcComp  = 0;
    int nGrow    = 0;
    int nComp    = 1;
    if(have_divu && have_dsdt) {
      FillPatchIterator divufpi(*this, divu, nGrow, destComp, time, Divu_Type,
                                srcComp, nComp);
      FillPatchIterator dsdtfpi(*this, dsdt, nGrow, destComp, time, Dsdt_Type,
                                srcComp, nComp);
      for(;
          divufpi.isValid() &&
          dsdtfpi.isValid();
          ++dsdtfpi,
          ++divufpi)
      {
        DependentMultiFabIterator divumfi(divufpi, divu);
        DependentMultiFabIterator dsdtmfi(divufpi, dsdt);
          //int k = divumfi.index();
          //getDivCond(divumfi(),k,0,time);
              divumfi().copy(divufpi());  // copy from the fill patch iterator

          // getDsdt(dsdtmfi(),k,0,time);
              // dont need to resize dsdt (grow == 0)
              dsdtmfi().copy(dsdtfpi());  // copy from the fill patch iterator

          dsdtmfi().mult(.5*dt);
          divumfi().plus(dsdtmfi());
      }
    
    } else if(have_divu &&  ! have_dsdt) {

      FillPatchIterator divufpi(*this, divu, nGrow, destComp, time, Divu_Type,
                                srcComp, nComp);
      for(; divufpi.isValid(); ++divufpi) {
        DependentMultiFabIterator divumfi(divufpi, divu);
        DependentMultiFabIterator dsdtmfi(divufpi, dsdt);
          //int k = divumfi.index();
          //getDivCond(divumfi(),k,0,time);
              divumfi().copy(divufpi());  // copy from the fill patch iterator
      }
      dsdt.setVal(0.0);
    
    } else {  // dont have either
      divu.setVal(0.0);
      dsdt.setVal(0.0);
    }

#endif

    //------------------- compute mac velocities and maximum clf number
    if (do_mac_proj)
        mac_projector->mac_project(level,u_mac,Sold,dt,time,divu,have_divu);
    mac_stats.end();
    
    //------------------- advect velocities
    vel_adv_stats.start();
    velocity_advection(dt);
    vel_adv_stats.end();

    //------------------- advect scalars
    int first_scalar = Density;
    int last_scalar = first_scalar + NUM_SCALARS - 1;
    scal_adv_stats.start();
    scalar_advection(dt,first_scalar,last_scalar);
    scal_adv_stats.end();

    // add the advective and other terms to get scalars at t^{n+1}
    scal_upd_stats.start();
    scalar_update(dt,first_scalar,first_scalar);

    // compute rho at half time level including 1-zone bndry values
    makerhonph(dt);

    // add the advective and other terms to get scalars at t^{n+1}
    scalar_update(dt,first_scalar+1,last_scalar);
    scal_upd_stats.end();

    // add the advective and other terms to get velocity at t^{n+1}
    vel_upd_stats.start();
    velocity_update(dt);
    vel_upd_stats.end();

    if(have_divu) {
      MultiFab &Divu_new = get_new_data(Divu_Type);
      calc_divu(time+dt, dt, Divu_new);
      if(have_dsdt) {
        MultiFab &Dsdt_new = get_new_data(Dsdt_Type);
        calc_dsdt(time,dt,Dsdt_new);
      }
    }

    // clean up after the predicted value at t^n+1
    // estimate new timestep from umac cfl);
    advance_cleanup(dt,iteration,ncycle);
    
    // increment rho average
    if (!initial_step && level > 0) {
        Real fratio = (Real) ncycle;
        Real alpha = 1.0/fratio;
        if (iteration == ncycle) alpha = 0.5/fratio;
        incrRhoAvg(alpha);
    }
    
    // do a level project to update the pressure and velocity fields
    if (!initial_step) {
        if ( projector )
            level_projector(dt,time,iteration);
        if (level > 0) {
            Real alpha = 1.0/ (Real) ncycle;
            incrPAvg(iteration,alpha);
        }
    }

    // relax back to continuity constraint
    if (divu_relax_factor>0.0 && !initial_step) {
        MultiFab* delta_U = new MultiFab (grids,BL_SPACEDIM,0,Fab_allocate);
        compute_grad_divu_minus_s(time+dt, delta_U, 0);
        Snew.plus(*delta_U, 0, BL_SPACEDIM, 0);
        delete delta_U;
    }

    return dt_test;  // return estimate of best new timestep
}

// -------------------------------------------------------------
void
NavierStokes::level_projector(Real dt, Real time, int iteration)
{
   if (iteration > 0) {
    RunStats lp_stats("level_project",level);
    lp_stats.start();

    MultiFab &U_old = get_old_data(State_Type);
    MultiFab &U_new = get_new_data(State_Type);
    MultiFab &P_old = get_old_data(Press_Type);
    MultiFab &P_new = get_new_data(Press_Type);

    int finest_level = parent->finestLevel();
    SyncRegister * crse_ptr;
    if (level < finest_level && do_sync_proj) {
      crse_ptr = &(getLevel(level+1).getSyncReg()); 
    } else {
      crse_ptr = 0;
    }

    Real cur_pres_time = state[Press_Type].curTime();

    int** sync_bc =  new (int*[grids.length()]);

    int i;
    for (i = 0; i < grids.length(); i++) {
      sync_bc[i] = getBCArray( State_Type,i,Xvel,BL_SPACEDIM);
    }

    MultiFab dsdt(grids,1,1,Fab_allocate);
    MultiFab divuold(grids,1,1,Fab_allocate);
    //int k;
    if(have_divu) {
#if 1
      FillStateBndry(time,Divu_Type,0,1);
      FillStateBndry(time+dt,Divu_Type,0,1);
      int ngrids = grids.length();
      //for (k = 0; k < ngrids; k++)
      for(MultiFabIterator divuoldmfi(divuold); divuoldmfi.isValid(); ++divuoldmfi)
      {
        DependentMultiFabIterator dsdtmfi(divuoldmfi, dsdt);
        int k = divuoldmfi.index();
        getDivCond(divuoldmfi(),k,1,time);
        getDivCond(dsdtmfi(),k,1,time+dt);
      }
#else
      MultiFab &divu_old = get_old_data(Divu_Type);
      MultiFab &divu_new = get_new_data(Divu_Type);
      divu_new.FillBoundary();
      divu_old.FillBoundary();
      setPhysBoundaryValues(Divu_Type,0,1,time+dt);
      // use dsdt for divunp1
      int ngrids = grids.length();
      //for (k = 0; k < ngrids; k++)
      for(MultiFabIterator divuoldmfi(divuold); divuoldmfi.isValid(); ++divuolmfi) {
        DependentMultiFabIterator dsdtmfi(divuoldmfi, dsdt);
        int k = divuoldmfi.index();
        getDivCond(divuoldmfi(),k,1,time);
        getDivCond(dsdtmfi(),k,1,time+dt);
      }
#endif
      dsdt.minus(divuold,0,1,1);
      //for (k = 0; k < ngrids; k++) {
      for(MultiFabIterator dsdtmfi(dsdt); dsdtmfi.isValid(); ++dsdtmfi) {
        assert(dsdtmfi.validbox() == grids[dsdtmfi.index()]);
        Box grid = dsdtmfi.validbox();
        grid.grow(1);
        dsdtmfi().divide(dt,grid,0,1);
      }
    } else {
      dsdt.setVal(0.0);
    } 

    int crse_dt_ratio;
    if (level > 0) {
      crse_dt_ratio = parent->MaxRefRatio(level-1);
    } else {
      crse_dt_ratio = -1;
    }

    projector->level_project(level,dt,cur_pres_time,time,time+dt,geom,
                             U_old,U_new,P_old,P_new,rho_half,dsdt,
                             crse_ptr,sync_reg,crse_dt_ratio,sync_bc,iteration,
                             divu_minus_s_factor,divuold,have_divu);

    for (i = 0; i < grids.length(); i++) {
      delete sync_bc[i];
    }

    delete[] sync_bc;

    lp_stats.end();

   } else {

    Real cur_pres_time = state[Press_Type].curTime();
    MultiFab &P_old = get_old_data(Press_Type);

    cout << "NavierStokes::level_projector calling harmonic_project\n";
    ParallelDescriptor::Abort("NavierStokes::level_projector");
    projector->harmonic_project(level,dt,cur_pres_time,geom,P_old);

   } 
}

// -------------------------------------------------------------
void
NavierStokes::makerhonph(Real dt)
{
    Real prev_time = state[State_Type].prevTime();
    Real half_time = prev_time + 0.5*dt;

#if (USEOLDFILLPATCH == 1)
    assert(ParallelDescriptor::NProcs() == 1);
    for(MultiFabIterator rhohalfmfi(*rho_half); rhohalfmfi.isValid(); ++rhohalfmfi)
    {
        int k = rhohalfmfi.index();
        getState(rhohalfmfi(),k,1,Density,1,half_time);
    }
#else
    int nGrow    = 1;
    int srcComp  = Density;
    int destComp = 0;
    int nComp    = 1;
    for(FillPatchIterator rho_halffpi(*this, *rho_half, nGrow,
                                  destComp, half_time, State_Type,
                                  srcComp, nComp);
        rho_halffpi.isValid(); ++rho_halffpi)
    {
      DependentMultiFabIterator rho_halfmfi(rho_halffpi, *rho_half);

      // the old getState grows the fab it is given--test rho_half for
      // the correct size  vvvvvvvvvvvvvvvvvvvvvvv
      Box testGrowBox(grids[rho_halfmfi.index()]);
      testGrowBox.grow(nGrow);
      assert(rho_halfmfi.fabbox() == testGrowBox);
      // end box size test ^^^^^^^^^^^^^^^^^^^^^^^

      rho_halfmfi().copy(rho_halffpi());
    }
#endif
}




// -------------------------------------------------------------
// predict the edge velocities which go into forming u_mac.  This
// function also returns an estimate of dt for use in variable timesteping
// -------------------------------------------------------------

Real NavierStokes::predict_velocity( Real dt, Real &comp_cfl )
{
    //int i;
    FArrayBox U, Rho, tforces, Gp;


    if (verbose) {
        cout << "... predict edge velocities\n";
    }

    // get simulation parameters
    const Real *dx = geom.CellSize();
    Real prev_time = state[State_Type].prevTime();
    Real prev_pres_time = state[Press_Type].prevTime();
    MultiFab &S_old = get_old_data(State_Type);

    // compute viscous terms at level n
    MultiFab visc_terms(grids,BL_SPACEDIM,1,Fab_allocate);
    getViscTerms(visc_terms,Xvel,BL_SPACEDIM,prev_time);
    if(be_cn_theta==1.0)visc_terms.setVal(0.0,1);

    // set up the timestep estimation
    Real cflgrid,cflmax,u_max[3];
    cflmax    = 1.0e-10;
    comp_cfl  = ( level == 0 ? cflmax : comp_cfl );

    // for each grid, fill patch source terms and predict
    //for ( i = 0; i < grids.length(); i++)
#if (USEOLDFILLPATCH == 1)
    for(MultiFabIterator visc_termsmfi(visc_terms);
        visc_termsmfi.isValid(); ++visc_termsmfi)
    {
      DependentMultiFabIterator u_mac0mfi(visc_termsmfi, u_mac[0]);
      DependentMultiFabIterator u_mac1mfi(visc_termsmfi, u_mac[1]);
#if (BL_SPACEDIM == 3)                         
      DependentMultiFabIterator u_mac2mfi(visc_termsmfi, u_mac[2]);
#endif
      int i = visc_termsmfi.index();

        // get grid and boundary condition info
        const Box& grd = grids[i];

        // get the relevant state data
        getState(U,      i,HYP_GROW,Xvel,   BL_SPACEDIM,prev_time);
        getState(Rho,    i,1,       Density,1,          prev_time);
        getForce(tforces,i,1,       Xvel,   BL_SPACEDIM,prev_time);
        getGradP(Gp,     i,1,                           prev_pres_time);

        // test velocities, rho and cfl
        cflgrid  = godunov->test_u_rho( U, Rho, grd, dx, dt, u_max );
        cflmax   = Max(cflgrid,cflmax);
        comp_cfl = Max(cflgrid,comp_cfl);

        // compute the total forcing
        godunov->Sum_tf_gp_visc( tforces, visc_termsmfi(), Gp, Rho );

        // set up the Godunov box
        godunov->Setup( grd, dx, dt, 1,
                        u_mac0mfi(), getBCArray( State_Type,i,0,1),
                        u_mac1mfi(), getBCArray( State_Type,i,1,1),
#if (BL_SPACEDIM == 3)                         
                        u_mac2mfi(), getBCArray( State_Type,i,2,1),
#endif
                        U, Rho, tforces);

        // predict the mac velocities
        godunov->ComputeUmac( grd, dx, dt, 
                              u_mac0mfi(), getBCArray( State_Type,i,0,1), 
                              u_mac1mfi(), getBCArray( State_Type,i,1,1), 
#if (BL_SPACEDIM == 3)
                              u_mac2mfi(), getBCArray( State_Type,i,2,1),
#endif
                              U, tforces );
    }  // end for(MultiFabIterator...)

#else

    int destComp = 0;
    FillPatchIterator Ufpi(*this, S_old, HYP_GROW, destComp, prev_time,
                           State_Type, Xvel, BL_SPACEDIM);
    FillPatchIterator Rhofpi(*this, S_old, 1, destComp, prev_time,
                           State_Type, Density, 1);
    //FillPatchIterator tforcesfpi(*this, S_old, 1, destComp, prev_time,
                           //State_Type, Xvel, BL_SPACEDIM);
    MultiFab &P_old = get_old_data(Press_Type);
    FillPatchIterator pressurefpi(*this, P_old, 1, destComp, prev_pres_time,
                           Press_Type, 0, 1);
    for(;
        Ufpi.isValid()       &&
        Rhofpi.isValid()     &&
        //tforcesfpi.isValid() &&
        pressurefpi.isValid();
        ++Ufpi,
        ++Rhofpi,
        //++tforcesfpi,
        ++pressurefpi)
    {
      DependentMultiFabIterator visc_termsmfi(Ufpi, visc_terms);
      DependentMultiFabIterator u_mac0mfi(Ufpi, u_mac[0]);
      DependentMultiFabIterator u_mac1mfi(Ufpi, u_mac[1]);
#if (BL_SPACEDIM == 3)                         
      DependentMultiFabIterator u_mac2mfi(Ufpi, u_mac[2]);
#endif
      int i = visc_termsmfi.index();

        // get grid and boundary condition info
        const Box& grd = grids[i];

        // get the relevant state data
        //getState(U,      i,HYP_GROW,Xvel,   BL_SPACEDIM,prev_time);
        FArrayBox &U = Ufpi();
        //getState(Rho,    i,1,       Density,1,          prev_time);
        FArrayBox &Rho = Rhofpi();
        //getForce(tforces,i,1,       Xvel,   BL_SPACEDIM,prev_time);
        // from NS::getForce vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        Real grav = Abs(gravity);
        Box tfbox(grids[i]);
        tfbox.grow(1);
        tforces.resize(tfbox, BL_SPACEDIM);
        for(int dc = 0; dc < BL_SPACEDIM; dc++) {
          int sc = Xvel + dc;
#if (BL_SPACEDIM == 2)
          if (BL_SPACEDIM == 2 && sc == Yvel && grav > 0.001)
#endif
#if (BL_SPACEDIM == 3)
          if (BL_SPACEDIM == 3 && sc == Zvel && grav > 0.001)
#endif
          {
              // set force to -rho*g
            FArrayBox rho(Rho.box());
            //getState(rho,i,1,Density,1,prev_time);
            rho.copy(Rho);
            rho.mult(-grav);
            tforces.copy(rho,0,dc,1);
          } else {
            tforces.setVal(0.0,dc);
          }
        }  // end for(dc...)
        // from NS::getForce ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



        //getGradP(Gp,     i,1,                           prev_pres_time);
        // from NS::getGradP vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        {
        Box gpbx(grids[i]);
        gpbx.grow(1);
        Box p_box(surroundingNodes(gpbx));
        //FArrayBox p_fab(p_box,1);
        FArrayBox &p_fab = pressurefpi();
        assert(p_fab.box() == p_box);
        //FillPatch(p_fab,0,time,Press_Type,0,1);
        //-----------------------  size the pressure gradient storage
        //Box gpbx(grd);
        //gpbx.grow(ngrow);
        Gp.resize(gpbx,BL_SPACEDIM);

        //------------------------ test to see if p_fab contains gpbx
        Box test = p_fab.box();
        test = test.enclosedCells();
        assert( test.contains( gpbx ) == true );

        // ----------------------- set pointers
        const int *plo = p_fab.loVect();
        const int *phi = p_fab.hiVect();
        const int *glo = gpbx.loVect();
        const int *ghi = gpbx.hiVect();
        const Real *p_dat  = p_fab.dataPtr();
        const Real *gp_dat = Gp.dataPtr();
        const Real *dx     = geom.CellSize();

        // ----------------------- create the pressure gradient
        FORT_GRADP (p_dat,ARLIM(plo),ARLIM(phi),
                    gp_dat,ARLIM(glo),ARLIM(ghi),glo,ghi,dx);
        }

        // from NS::getGradP ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


        // test velocities, rho and cfl
        cflgrid  = godunov->test_u_rho( U, Rho, grd, dx, dt, u_max );
        cflmax   = Max(cflgrid,cflmax);
        comp_cfl = Max(cflgrid,comp_cfl);

        // compute the total forcing
        godunov->Sum_tf_gp_visc( tforces, visc_termsmfi(), Gp, Rho );

        // set up the Godunov box
        godunov->Setup( grd, dx, dt, 1,
                        u_mac0mfi(), getBCArray( State_Type,i,0,1),
                        u_mac1mfi(), getBCArray( State_Type,i,1,1),
#if (BL_SPACEDIM == 3)                         
                        u_mac2mfi(), getBCArray( State_Type,i,2,1),
#endif
                        U, Rho, tforces);

        // predict the mac velocities
        godunov->ComputeUmac( grd, dx, dt, 
                              u_mac0mfi(), getBCArray( State_Type,i,0,1), 
                              u_mac1mfi(), getBCArray( State_Type,i,1,1), 
#if (BL_SPACEDIM == 3)
                              u_mac2mfi(), getBCArray( State_Type,i,2,1),
#endif
                              U, tforces );
    }  // end for(FillPatchIterator...)

#endif


    // test for periodic u_mac
    if ( level == 0 && geom.isAnyPeriodic() ) {
        test_umac_periodic();
    }

    // compute estimate of the timestep
    Real tempdt = Min(change_max,cfl/cflmax);
    ParallelDescriptor::ReduceRealMin(tempdt);
    if(ParallelDescriptor::NProcs() > 1) {
      cout << "\n\nCheck reduction ops for zero grids on processor.\n\n\n";
    }
    return dt*tempdt;
    //return dt*Min(change_max,cfl/cflmax);
}  // end predict_velocity()


// test for periodic umac on a single level 0 grid
void NavierStokes::test_umac_periodic()
{
    // error block
    assert( ParallelDescriptor::NProcs() == 1);
    assert( level          == 0 );
    if ( grids.length() != 1 )
        return;

    // get the bounds and grid size
    const Box &grd = grids[0];
    const int *lo  = grd.loVect();
    const int *hi  = grd.hiVect();
    
    // get the velocities
    int xperiod        = ( geom.isPeriodic(0) ? 1 : 0 );
    const int *ux_lo   = u_mac[0][0].loVect();
    const int *ux_hi   = u_mac[0][0].hiVect();
    const Real *ux_dat = u_mac[0][0].dataPtr();

    int yperiod        = ( geom.isPeriodic(1) ? 1 : 0 );
    const int *uy_lo   = u_mac[1][0].loVect();
    const int *uy_hi   = u_mac[1][0].hiVect();
    const Real *uy_dat = u_mac[1][0].dataPtr();
#if (BL_SPACEDIM == 3)
    int zperiod        = ( geom.isPeriodic(2) ? 1 : 0 );
    const int *uz_lo   = u_mac[2][0].loVect();
    const int *uz_hi   = u_mac[2][0].hiVect();
    const Real *uz_dat = u_mac[2][0].dataPtr();
#endif

    // call the fortran
    Real udiff = 0.0;
    Real vdiff = 0.0;
#if (BL_SPACEDIM == 3)                                             
    Real wdiff = 0.0;
#endif
    
    FORT_TEST_UMAC_PERIODIC( lo, hi, 
                             ux_dat, ARLIM(ux_lo), ARLIM(ux_hi),
                             &xperiod, &udiff,
                             uy_dat, ARLIM(uy_lo), ARLIM(uy_hi),
                             &yperiod, &vdiff,
#if (BL_SPACEDIM == 3)                                             
                             uz_dat, ARLIM(uz_lo), ARLIM(uz_hi),
                             &zperiod, &wdiff,
#endif
                             &level );
    
    if ( level == 0 && udiff > 1.0e-10 ) {
        cout << "!!!!!!!!!!!!!!!!!!!!!!!! udiff = " << udiff << NL;
        ParallelDescriptor::Abort("Exiting.");
    }
    if ( level == 0 && vdiff > 1.0e-10 ) {
        cout << "!!!!!!!!!!!!!!!!!!!!!!!! vdiff = " << vdiff << NL;
        ParallelDescriptor::Abort("Exiting.");
    }
#if (BL_SPACEDIM == 3)                                             
    if ( level == 0 && vdiff > 1.0e-10 ) {
        cout << "!!!!!!!!!!!!!!!!!!!!!!!! vdiff = " << vdiff << NL;
        ParallelDescriptor::Abort("Exiting.");
    }
#endif
}


//
// This routine advects the velocities
//
void NavierStokes::velocity_advection( Real dt )
{

    if (verbose) {
        cout << "... advect velocities\n";
    }

    // get simulation parameters
    int finest_level    = parent->finestLevel();
    const Real *dx      = geom.CellSize();
    Real prev_time      = state[State_Type].prevTime();
    Real prev_pres_time = state[Press_Type].prevTime();

    // compute viscosity components
    MultiFab visc_terms(grids,BL_SPACEDIM,1,Fab_allocate);
    getViscTerms(visc_terms,Xvel,BL_SPACEDIM,prev_time);
    if(be_cn_theta==1.0)visc_terms.setVal(0.0,1);

    // set up the grid loop
    int comp;
    FArrayBox U, Rho, tforces, Gp;
    FArrayBox xflux, yflux, zflux, divu;

    // compute the advective forcing
    //for ( i = 0; i < grids.length(); i++)
#if (USEOLDFILLPATCH == 1)
    for(MultiFabIterator visc_termsmfi(visc_terms);
        visc_termsmfi.isValid(); ++visc_termsmfi)
    {
      DependentMultiFabIterator aofsmfi(visc_termsmfi, (*aofs));
      DependentMultiFabIterator volumemfi(visc_termsmfi, volume);
      DependentMultiFabIterator u_mac0mfi(visc_termsmfi, u_mac[0]);
      DependentMultiFabIterator u_mac1mfi(visc_termsmfi, u_mac[1]);
      DependentMultiFabIterator area0mfi(visc_termsmfi, area[0]);
      DependentMultiFabIterator area1mfi(visc_termsmfi, area[1]);
#if (BL_SPACEDIM == 3)
      DependentMultiFabIterator u_mac2mfi(visc_termsmfi, u_mac[2]);
      DependentMultiFabIterator area2mfi(visc_termsmfi, area[2]);
#endif
      int i = visc_termsmfi.index();

        assert(grids[visc_termsmfi.index()] == visc_termsmfi.validbox());
        const Box& grd = visc_termsmfi.validbox();

        // get needed data
        getState(U,      i,HYP_GROW,Xvel,   BL_SPACEDIM,prev_time);
        getState(Rho,    i,1,       Density,1,          prev_time);
        getForce(tforces,i,1,       Xvel,   BL_SPACEDIM,prev_time);
        getGradP(Gp,     i,1,                      prev_pres_time);

        // compute the total forcing
        godunov->Sum_tf_gp_visc( tforces, visc_termsmfi(), Gp, Rho );

        // set up the workspace for the godunov Box
        godunov->Setup( grd, dx, dt, 0,
                        xflux, getBCArray( State_Type,i,0,1),
                        yflux, getBCArray( State_Type,i,1,1),
#if (BL_SPACEDIM == 3)                          
                        zflux, getBCArray( State_Type,i,2,1),
#endif
                        U, Rho, tforces);
        
        // loop over the velocity components
        for ( comp = 0 ; comp < BL_SPACEDIM ; comp++ ) {
            godunov->AdvectState( grd, dx, dt, 
                                  area0mfi(), u_mac0mfi(), xflux,
                                  area1mfi(), u_mac1mfi(), yflux,
#if (BL_SPACEDIM == 3)                       
                                  area2mfi(), u_mac2mfi(), zflux,
#endif
                                  U, U, tforces, comp,
                                  aofsmfi(),    comp,
                                  is_conservative[comp],
                                  comp,
                                  getBCArray( State_Type,i,comp,1),
                                  volumemfi() );

            // get fluxes for diagnostics and refluxing
            pullFluxes( level, i, comp, 1, xflux, yflux, zflux, dt );
        } // end of velocity loop
    } // end of grid loop

#else

    MultiFab &S_old = get_old_data(State_Type);
    int destComp = 0;
    FillPatchIterator Ufpi(*this, S_old, HYP_GROW, destComp, prev_time,
                           State_Type, Xvel, BL_SPACEDIM);
    FillPatchIterator Rhofpi(*this, S_old, 1, destComp, prev_time,
                           State_Type, Density, 1);
    //FillPatchIterator tforcesfpi(*this, S_old, 1, destComp, prev_time,
                           //State_Type, Xvel, BL_SPACEDIM);
    MultiFab &P_old = get_old_data(Press_Type);
    FillPatchIterator pressurefpi(*this, P_old, 1, destComp, prev_pres_time,
                           Press_Type, 0, 1);
    for(;
        Ufpi.isValid()       &&
        Rhofpi.isValid()     &&
        //tforcesfpi.isValid() &&
        pressurefpi.isValid();
        ++Ufpi,
        ++Rhofpi,
        //++tforcesfpi,
        ++pressurefpi)
    {

      DependentMultiFabIterator visc_termsmfi(Ufpi, visc_terms);
      DependentMultiFabIterator aofsmfi(Ufpi, (*aofs));
      DependentMultiFabIterator volumemfi(Ufpi, volume);
      DependentMultiFabIterator u_mac0mfi(Ufpi, u_mac[0]);
      DependentMultiFabIterator u_mac1mfi(Ufpi, u_mac[1]);
      DependentMultiFabIterator area0mfi(Ufpi, area[0]);
      DependentMultiFabIterator area1mfi(Ufpi, area[1]);
#if (BL_SPACEDIM == 3)
      DependentMultiFabIterator u_mac2mfi(Ufpi, u_mac[2]);
      DependentMultiFabIterator area2mfi(Ufpi, area[2]);
#endif
      int i = visc_termsmfi.index();

        assert(grids[visc_termsmfi.index()] == visc_termsmfi.validbox());
        const Box& grd = visc_termsmfi.validbox();

        // get needed data
        //getState(U,      i,HYP_GROW,Xvel,   BL_SPACEDIM,prev_time);
        FArrayBox &U = Ufpi();
        //getState(Rho,    i,1,       Density,1,          prev_time);
        FArrayBox &Rho = Rhofpi();

        //getForce(tforces,i,1,       Xvel,   BL_SPACEDIM,prev_time);
        // from NS::getForce vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        Real grav = Abs(gravity);
        Box tfbox(grids[i]);
        tfbox.grow(1);
        tforces.resize(tfbox, BL_SPACEDIM);
        for(int dc = 0; dc < BL_SPACEDIM; dc++) {
          int sc = Xvel + dc;
#if (BL_SPACEDIM == 2)
          if (BL_SPACEDIM == 2 && sc == Yvel && grav > 0.001)
#endif
#if (BL_SPACEDIM == 3)
          if (BL_SPACEDIM == 3 && sc == Zvel && grav > 0.001)
#endif
          {
              // set force to -rho*g
            FArrayBox rho(Rho.box());
            //getState(rho,i,1,Density,1,prev_time);
            rho.copy(Rho);
            rho.mult(-grav);
            tforces.copy(rho,0,dc,1);
          } else {
            tforces.setVal(0.0,dc);
          }
        }  // end for(dc...)
        // from NS::getForce ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


        //getGradP(Gp,     i,1,                      prev_pres_time);
        // from NS::getGradP vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        {
        Box gpbx(grids[i]);
        gpbx.grow(1);
        Box p_box(surroundingNodes(gpbx));
        //FArrayBox p_fab(p_box,1);
        FArrayBox &p_fab = pressurefpi();
        assert(p_fab.box() == p_box);
        //FillPatch(p_fab,0,time,Press_Type,0,1);
        //-----------------------  size the pressure gradient storage
        //Box gpbx(grd);
        //gpbx.grow(ngrow);
        Gp.resize(gpbx,BL_SPACEDIM);

        //------------------------ test to see if p_fab contains gpbx
        Box test = p_fab.box();
        test = test.enclosedCells();
        assert( test.contains( gpbx ) == true );

        // ----------------------- set pointers
        const int *plo = p_fab.loVect();
        const int *phi = p_fab.hiVect();
        const int *glo = gpbx.loVect();
        const int *ghi = gpbx.hiVect();
        const Real *p_dat  = p_fab.dataPtr();
        const Real *gp_dat = Gp.dataPtr();
        const Real *dx     = geom.CellSize();

        // ----------------------- create the pressure gradient
        FORT_GRADP (p_dat,ARLIM(plo),ARLIM(phi),
                    gp_dat,ARLIM(glo),ARLIM(ghi),glo,ghi,dx);
        }
        // from NS::getGradP ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


        // compute the total forcing
        godunov->Sum_tf_gp_visc( tforces, visc_termsmfi(), Gp, Rho );

        // set up the workspace for the godunov Box
        godunov->Setup( grd, dx, dt, 0,
                        xflux, getBCArray( State_Type,i,0,1),
                        yflux, getBCArray( State_Type,i,1,1),
#if (BL_SPACEDIM == 3)                          
                        zflux, getBCArray( State_Type,i,2,1),
#endif
                        U, Rho, tforces);
        
        // loop over the velocity components
        for ( comp = 0 ; comp < BL_SPACEDIM ; comp++ ) {
            godunov->AdvectState( grd, dx, dt, 
                                  area0mfi(), u_mac0mfi(), xflux,
                                  area1mfi(), u_mac1mfi(), yflux,
#if (BL_SPACEDIM == 3)                       
                                  area2mfi(), u_mac2mfi(), zflux,
#endif
                                  U, U, tforces, comp,
                                  aofsmfi(),    comp,
                                  is_conservative[comp],
                                  comp,
                                  getBCArray( State_Type,i,comp,1),
                                  volumemfi() );

            // get fluxes for diagnostics and refluxing
            pullFluxes( level, i, comp, 1, xflux, yflux, zflux, dt );
        } // end of velocity loop
    } // end of grid loop
#endif
}




//
// This routine advects the scalars
//
void NavierStokes::scalar_advection( Real dt, int fscalar, int lscalar)
{

    if (verbose) {
        cout << "... advect scalars\n";
    }

    // get simulation parameters
    int num_scalars  = lscalar - fscalar + 1;
    int finest_level = parent->finestLevel();
    const Real *dx   = geom.CellSize();
    Real prev_time   = state[State_Type].prevTime();
    Real prev_pres_time   = state[Press_Type].prevTime();

    // get the viscous terms
    MultiFab visc_terms(grids,num_scalars,1,Fab_allocate);
    getViscTerms(visc_terms,fscalar,num_scalars,prev_time);
    if(be_cn_theta==1.0)visc_terms.setVal(0.0,1);

    // set up the grid loop
    int comp,state_ind;
    FArrayBox U, S, Rho, tforces;
    FArrayBox xflux, yflux, zflux, divu;

    int use_forces_in_trans = godunov->useForcesInTrans();

    FArrayBox tvelforces, Gp;
    MultiFab vel_visc_terms;
    if(use_forces_in_trans) {
      vel_visc_terms.define(grids,BL_SPACEDIM,1,Fab_allocate);
      getViscTerms(vel_visc_terms,Xvel,BL_SPACEDIM,prev_time);
      if(be_cn_theta==1.0)vel_visc_terms.setVal(0.0,1);
    }

    // compute the advective forcing
    //for ( i = 0; i < grids.length(); i++)
#if (USEOLDFILLPATCH == 1)
    for(MultiFabIterator visc_termsmfi(visc_terms);
        visc_termsmfi.isValid(); ++visc_termsmfi)
    {
      DependentMultiFabIterator vel_visc_termsmfi(visc_termsmfi, vel_visc_terms);
      DependentMultiFabIterator aofsmfi(visc_termsmfi, (*aofs));
      DependentMultiFabIterator volumemfi(visc_termsmfi, volume);
      DependentMultiFabIterator u_mac0mfi(visc_termsmfi, u_mac[0]);
      DependentMultiFabIterator u_mac1mfi(visc_termsmfi, u_mac[1]);
      DependentMultiFabIterator area0mfi(visc_termsmfi, area[0]);
      DependentMultiFabIterator area1mfi(visc_termsmfi, area[1]);
#if (BL_SPACEDIM == 3)
      DependentMultiFabIterator u_mac2mfi(visc_termsmfi, u_mac[2]);
      DependentMultiFabIterator area2mfi(visc_termsmfi, area[2]);
#endif
      int i = visc_termsmfi.index();

        assert(grids[visc_termsmfi.index()] == visc_termsmfi.validbox());
        const Box& grd = grids[i];

        // get needed data
        getState(U,      i,HYP_GROW,Xvel,   BL_SPACEDIM,prev_time);
        getState(S,      i,HYP_GROW,fscalar,num_scalars,prev_time);
        getState(Rho,    i,1,       Density,1,          prev_time);
        getForce(tforces,i,1,       fscalar,num_scalars,prev_time);
        getDivCond(divu, i,1,                           prev_time);

        if (use_forces_in_trans) {
          getForce(tvelforces,i,1,       Xvel,   BL_SPACEDIM,prev_time);
          getGradP(Gp,        i,1,                           prev_pres_time);
          godunov->Sum_tf_gp_visc( tvelforces, vel_visc_termsmfi(), Gp, Rho );
        }
        
        // set up the workspace for the godunov Box
        godunov->Setup( grd, dx, dt, 0,
                        xflux, getBCArray( State_Type,i,0,1),
                        yflux, getBCArray( State_Type,i,1,1),
#if (BL_SPACEDIM == 3)                         
                        zflux, getBCArray( State_Type,i,2,1),
#endif
                        U, Rho, tvelforces );
        
        // loop over the scalar components
        for ( comp = 0 ; comp < num_scalars ; comp++ ) {
            state_ind = fscalar + comp;

            // compute total forcing
            godunov->Sum_tf_divu_visc( S, tforces,    comp, 1,
                                       visc_termsmfi(), comp,
                                       divu, Rho,
                                       is_conservative[state_ind] );

            //  advect scalar
            godunov->AdvectState( grd, dx, dt, 
                                  area0mfi(), u_mac0mfi(), xflux,
                                  area1mfi(), u_mac1mfi(), yflux,
#if (BL_SPACEDIM == 3)                        
                                  area2mfi(), u_mac2mfi(), zflux,
#endif
                                  U, S, tforces, comp,
                                  aofsmfi(),    state_ind,
                                  is_conservative[state_ind],
                                  state_ind,
                                  getBCArray( State_Type,i,state_ind,1),
                                  volumemfi() );

            // get the fluxes for refluxing and diagnostic purposes
            pullFluxes( level, i, state_ind, 1, xflux, yflux, zflux, dt );
            
        } // end of scalar loop
    } // end of grid loop

#else
    MultiFab &S_old = get_old_data(State_Type);
    int destComp = 0;
    FillPatchIterator Ufpi(*this, S_old, HYP_GROW, destComp, prev_time,
                           State_Type, Xvel, BL_SPACEDIM);
    FillPatchIterator Sfpi(*this, S_old, HYP_GROW, destComp, prev_time,
                           State_Type, fscalar, num_scalars);
    FillPatchIterator Rhofpi(*this, S_old, 1, destComp, prev_time,
                           State_Type, Density, 1);
    //FillPatchIterator tforcesfpi(*this, S_old, 1, destComp, prev_time,
                           //State_Type, fscalar, num_scalars);
    FillPatchIterator tvelforcesfpi(*this, S_old, 1, destComp, prev_time,
                           State_Type, Xvel, BL_SPACEDIM);
    MultiFab &P_old = get_old_data(Press_Type);
    FillPatchIterator pressurefpi(*this, P_old, 1, destComp, prev_pres_time,
                           Press_Type, 0, 1);

    FillPatchIterator divufpi(*this, S_old);
    if(have_divu) {
      int boxGrow = 1;
      int srcComp = 0;
      int nComp   = 1;
      divufpi.Initialize(boxGrow, destComp, prev_time,
                           Divu_Type, srcComp, nComp);
    }

    bool bIteratorsValid;
    if(have_divu) {
      bIteratorsValid = 
        Ufpi.isValid()          &&
        Sfpi.isValid()          &&
        Rhofpi.isValid()        &&
        tvelforcesfpi.isValid() &&
        pressurefpi.isValid()   &&
        divufpi.isValid();
    } else {
      bIteratorsValid = 
        Ufpi.isValid()          &&
        Sfpi.isValid()          &&
        Rhofpi.isValid()        &&
        tvelforcesfpi.isValid() &&
        pressurefpi.isValid();
    }
    while(bIteratorsValid)
    {

      DependentMultiFabIterator visc_termsmfi(Ufpi, visc_terms);
      DependentMultiFabIterator vel_visc_termsmfi(Ufpi, vel_visc_terms);
      DependentMultiFabIterator aofsmfi(Ufpi, (*aofs));
      DependentMultiFabIterator volumemfi(Ufpi, volume);
      DependentMultiFabIterator u_mac0mfi(Ufpi, u_mac[0]);
      DependentMultiFabIterator u_mac1mfi(Ufpi, u_mac[1]);
      DependentMultiFabIterator area0mfi(Ufpi, area[0]);
      DependentMultiFabIterator area1mfi(Ufpi, area[1]);
#if (BL_SPACEDIM == 3)
      DependentMultiFabIterator u_mac2mfi(Ufpi, u_mac[2]);
      DependentMultiFabIterator area2mfi(Ufpi, area[2]);
#endif
      int i = visc_termsmfi.index();

        assert(grids[visc_termsmfi.index()] == visc_termsmfi.validbox());
        const Box& grd = grids[i];

        // get needed data
        //getState(U,      i,HYP_GROW,Xvel,   BL_SPACEDIM,prev_time);
        FArrayBox &U = Ufpi();
        //getState(S,      i,HYP_GROW,fscalar,num_scalars,prev_time);
        FArrayBox &S = Sfpi();
        //getState(Rho,    i,1,       Density,1,          prev_time);
        FArrayBox &Rho = Rhofpi();


        //getForce(tforces,i,1,       fscalar,num_scalars,prev_time);
        // from NS::getForce vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        Real grav = Abs(gravity);
        Box tfbox(grids[i]);
        tfbox.grow(1);
        tforces.resize(tfbox, BL_SPACEDIM);
        for(int dc = 0; dc < num_scalars; dc++) {
          int sc = fscalar + dc;
#if (BL_SPACEDIM == 2)
          if (BL_SPACEDIM == 2 && sc == Yvel && grav > 0.001)
#endif
#if (BL_SPACEDIM == 3)
          if (BL_SPACEDIM == 3 && sc == Zvel && grav > 0.001)
#endif
          {
              // set force to -rho*g
            FArrayBox rho(Rho.box());
            //getState(rho,i,1,Density,1,prev_time);
            rho.copy(Rho);
            rho.mult(-grav);
            tforces.copy(rho,0,dc,1);
          } else {
            tforces.setVal(0.0,dc);
          }
        }  // end for(dc...)
        // from NS::getForce ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        //getDivCond(divu, i,1,                           prev_time);
            Box bx(grids[i]);
            bx.grow(1);
            divu.resize(bx,1);
          if(have_divu) {
            if(divu.box() != divufpi().box()) {
              cerr << "NavierStokes::scalar_advection:  "
                   << "divu.box() != divufpi().box()\n";
              cerr << "divu.box()      = " << divu.box()      << NL;
              cerr << "divufpi().box() = " << divufpi().box() << NL;
              ParallelDescriptor::Abort("NavierStokes::scalar_advection(...)");
            }
            divu.copy(divufpi());
          } else {   // ! have_divu
            divu.setVal(0.0);  // can't filpatch what we don't have
          }

        if (use_forces_in_trans) {
          //getForce(tvelforces,i,1,       Xvel,   BL_SPACEDIM,prev_time);

        // from NS::getForce vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
          Real grav = Abs(gravity);
          for(int dc = 0; dc < BL_SPACEDIM; dc++) {
            int sc = Xvel + dc;
#if (BL_SPACEDIM == 2)
            if (BL_SPACEDIM == 2 && sc == Yvel && grav > 0.001)
#endif
#if (BL_SPACEDIM == 3)
            if (BL_SPACEDIM == 3 && sc == Zvel && grav > 0.001)
#endif
            {
                // set force to -rho*g
              FArrayBox rho(Rho.box());
              //getState(rho,i,1,Density,1,prev_time);
              rho.copy(Rho);
              rho.mult(-grav);
              tvelforces.copy(rho,0,dc,1);
            } else {
              tvelforces.setVal(0.0,dc);
            }
          }  // end for(dc...)
        // from NS::getForce ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


          //getGradP(Gp,        i,1,                           prev_pres_time);
        // from NS::getGradP vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        {
          Box gpbx(grids[i]);
          gpbx.grow(1);
          Box p_box(surroundingNodes(gpbx));
          //FArrayBox p_fab(p_box,1);
          FArrayBox &p_fab = pressurefpi();
          assert(p_fab.box() == p_box);
          //-----------------------  size the pressure gradient storage
          Gp.resize(gpbx,BL_SPACEDIM);

          //------------------------ test to see if p_fab contains gpbx
          Box test = p_fab.box();
          test = test.enclosedCells();
          assert( test.contains( gpbx ) == true );

          // ----------------------- set pointers
          const int *plo = p_fab.loVect();
          const int *phi = p_fab.hiVect();
          const int *glo = gpbx.loVect();
          const int *ghi = gpbx.hiVect();
          const Real *p_dat  = p_fab.dataPtr();
          const Real *gp_dat = Gp.dataPtr();
          const Real *dx     = geom.CellSize();

          // ----------------------- create the pressure gradient
          FORT_GRADP (p_dat,ARLIM(plo),ARLIM(phi),
                      gp_dat,ARLIM(glo),ARLIM(ghi),glo,ghi,dx);
        }
        // from NS::getGradP ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

          godunov->Sum_tf_gp_visc( tvelforces, vel_visc_termsmfi(), Gp, Rho );
        }
        
        // set up the workspace for the godunov Box
        godunov->Setup( grd, dx, dt, 0,
                        xflux, getBCArray( State_Type,i,0,1),
                        yflux, getBCArray( State_Type,i,1,1),
#if (BL_SPACEDIM == 3)                         
                        zflux, getBCArray( State_Type,i,2,1),
#endif
                        U, Rho, tvelforces );
        
        // loop over the scalar components
        for ( comp = 0 ; comp < num_scalars ; comp++ ) {
            state_ind = fscalar + comp;

            // compute total forcing
            godunov->Sum_tf_divu_visc( S, tforces,    comp, 1,
                                       visc_termsmfi(), comp,
                                       divu, Rho,
                                       is_conservative[state_ind] );

            //  advect scalar
            godunov->AdvectState( grd, dx, dt, 
                                  area0mfi(), u_mac0mfi(), xflux,
                                  area1mfi(), u_mac1mfi(), yflux,
#if (BL_SPACEDIM == 3)                        
                                  area2mfi(), u_mac2mfi(), zflux,
#endif
                                  U, S, tforces, comp,
                                  aofsmfi(),    state_ind,
                                  is_conservative[state_ind],
                                  state_ind,
                                  getBCArray( State_Type,i,state_ind,1),
                                  volumemfi() );

            // get the fluxes for refluxing and diagnostic purposes
            pullFluxes( level, i, state_ind, 1, xflux, yflux, zflux, dt );
            
        } // end of scalar loop

        ++Ufpi;
        ++Sfpi;
        ++Rhofpi;
        ++tvelforcesfpi;
        ++pressurefpi;
        if(have_divu) {
          ++divufpi;
        }

        if(have_divu) {
          bIteratorsValid = 
            Ufpi.isValid()          &&
            Sfpi.isValid()          &&
            Rhofpi.isValid()        &&
            tvelforcesfpi.isValid() &&
            pressurefpi.isValid()   &&
            divufpi.isValid();
        } else {
          bIteratorsValid = 
            Ufpi.isValid()          &&
            Sfpi.isValid()          &&
            Rhofpi.isValid()        &&
            tvelforcesfpi.isValid() &&
            pressurefpi.isValid();
        }

    } // end of while FillPatchIterator loop
#endif
}



//-------------------------------------------------------------
//  This subroutine updates the scalars, before the velocity update
//  and the level projection
//
//  AT this point in time, all we know is psi^n, rho^n+1/2, and the
//  general forcing terms at t^n, and after solving in this routine
//  viscous forcing at t^n+1/2.  Note, unless more complicated logic
//  is invoked earlier, we do not have any estimate of general forcing
//  terms at t^n+1/2.
//-------------------------------------------------------------

void NavierStokes::scalar_update(Real dt, int first_scalar, int last_scalar)
{

    if (verbose) {
        cout << "... update scalars\n";
    }

    scalar_advection_update(dt, first_scalar, last_scalar);
    scalar_diffusion_update(dt, first_scalar, last_scalar);

}

//-------------------------------------------------------------

void NavierStokes::scalar_advection_update(Real dt, int first_scalar, int last_scalar)
{

    // simulation parameters
    int finest_level = parent->finestLevel();
    MultiFab &S_old = get_old_data(State_Type);
    MultiFab &S_new = get_new_data(State_Type);
    MultiFab &Aofs  = *aofs;
    int ngrids        = grids.length();
    Real cur_time   = state[State_Type].curTime();
    Real prev_time  = state[State_Type].prevTime();
    Real half_time  = 0.5*(cur_time+prev_time);
    FArrayBox tforces;
    
    // loop over the desired scalars
    for (int sigma = first_scalar; sigma <= last_scalar; sigma++) {

#if (USEOLDFILLPATCH == 1)
        // compute inviscid estimate of scalars
        //for (int i = 0; i < ngrids; i++)
        for(MultiFabIterator S_oldmfi(S_old);
            S_oldmfi.isValid(); ++S_oldmfi)
        {
            DependentMultiFabIterator S_newmfi(S_oldmfi, S_new);
            DependentMultiFabIterator Aofsmfi(S_oldmfi, Aofs);
            assert(grids[S_oldmfi.index()] == S_oldmfi.validbox());
            const Box& grd = S_oldmfi.validbox();
            int i = S_oldmfi.index();
            getForce(tforces,i,0,sigma,1,half_time);
            godunov->Add_aofs_tf( S_oldmfi(),
                                  S_newmfi(), sigma, 1, 
                                  Aofsmfi(),  sigma,
                                  tforces,  0,
                                  grd,      dt );

            int do_minmax = (is_conservative[sigma] ? 0 : 1);
            if (do_minmax) {
                godunov->ScalMinMax( S_oldmfi(), S_newmfi(),
                                     sigma,
                                     getBCArray( State_Type,i,sigma,1),
                                     grd );
            }
        }

#else

        int destComp = 0;
        for(FillPatchIterator Rhofpi(*this, S_old, 0, destComp, half_time,
                                     State_Type, Density, 1);
            Rhofpi.isValid(); ++Rhofpi)
        {

            DependentMultiFabIterator S_oldmfi(Rhofpi, S_old);
            DependentMultiFabIterator S_newmfi(Rhofpi, S_new);
            DependentMultiFabIterator Aofsmfi(Rhofpi, Aofs);
            assert(grids[S_oldmfi.index()] == S_oldmfi.validbox());
            const Box& grd = S_oldmfi.validbox();
            int i = S_oldmfi.index();
            //getForce(tforces,i,0,sigma,1,half_time);
            // from NS::getForce vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
            Real grav = Abs(gravity);
            Box tfbox(grids[i]);
            tfbox.grow(0);
            tforces.resize(tfbox, BL_SPACEDIM);
            for(int dc = 0; dc < BL_SPACEDIM; dc++) {
              int sc = sigma + dc;
#if (BL_SPACEDIM == 2)
              if (BL_SPACEDIM == 2 && sc == Yvel && grav > 0.001)
#endif
#if (BL_SPACEDIM == 3)
              if (BL_SPACEDIM == 3 && sc == Zvel && grav > 0.001)
#endif
              {
                  // set force to -rho*g
                FArrayBox &Rho = Rhofpi();
                FArrayBox rho(Rho.box());
                //getState(rho,i,1,Density,1,half_time);
                rho.copy(Rhofpi());
                rho.mult(-grav);
                tforces.copy(rho,0,dc,1);
              } else {
                tforces.setVal(0.0,dc);
              }
            }  // end for(dc...)
            // from NS::getForce ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

            godunov->Add_aofs_tf( S_oldmfi(),
                                  S_newmfi(), sigma, 1, 
                                  Aofsmfi(),  sigma,
                                  tforces,  0,
                                  grd,      dt );

            int do_minmax = (is_conservative[sigma] ? 0 : 1);
            if (do_minmax) {
                godunov->ScalMinMax( S_oldmfi(), S_newmfi(),
                                     sigma,
                                     getBCArray( State_Type,i,sigma,1),
                                     grd );
            }
        }
#endif
    }
}

//-------------------------------------------------------------

void NavierStokes::scalar_diffusion_update(Real dt, int first_scalar, int last_scalar)
{

    // loop over the desired scalars
    for (int sigma = first_scalar; sigma <= last_scalar; sigma++) {
        // compute diffusion
        if (is_diffusive[sigma]) {
            int rho_flag = 0;
            diffuse_scalar_setup(sigma, &rho_flag);
            diffusion->diffuse_scalar(dt, sigma, be_cn_theta,
                                      rho_half, rho_flag);
        }  
    }    
}

// -------------------------------------------------------------
void NavierStokes::diffuse_scalar_setup(int sigma, int* rho_flag) 
{
    if(!is_conservative[sigma]) {
      (*rho_flag)=1;
    } else {
      (*rho_flag)=2;
    }
}

//-------------------------------------------------------------
//  This subroutine updates the velocity field before the level
//  projection
//
//  AT this point in time, all we know is u^n, rho^n+1/2, and the
//  general forcing terms at t^n, and after solving in this routine
//  viscous forcing at t^n+1/2.  Except for a simple buoyancy term,
//  b = -rho^n+1/2 g, it is usually not possible to estimate more
//  general forcing terms at t^n+1/2.  Since the default getForce, handles
//  this case automatically, F_new and F_old have been replaced by a single
//  tforces FArrayBox
//
//  we assume that if one component of velocity is viscous that
//  all must be.
//-------------------------------------------------------------

void NavierStokes::velocity_update(Real dt)
{
    if (verbose) {
        cout << "... update velocities\n";
    }
    velocity_advection_update(dt);
    if (!initial_iter) 
      velocity_diffusion_update(dt);
    else
      initial_velocity_diffusion_update(dt);

}

//-------------------------------------------------------------

void NavierStokes::velocity_advection_update(Real dt)
{
    FArrayBox Gp, tforces;

    // simulation parameters
    int finest_level    = parent->finestLevel();
    MultiFab &U_old     = get_old_data(State_Type);
    MultiFab &U_new     = get_new_data(State_Type);
    MultiFab &Aofs      = *aofs;
    Real cur_time       = state[State_Type].curTime();
    Real prev_time      = state[State_Type].prevTime();
    Real half_time      = 0.5*(prev_time+cur_time);
    Real pres_prev_time = state[Press_Type].prevTime();

    // estimate u^n+1 and put in U_new
    //for (i = 0; i < grids.length(); i++)

#if (USEOLDFILLPATCH == 1)

    for(MultiFabIterator U_oldmfi(U_old); U_oldmfi.isValid(); ++U_oldmfi) {
        DependentMultiFabIterator U_newmfi(U_oldmfi, U_new);
        DependentMultiFabIterator Aofsmfi(U_oldmfi, Aofs);
        DependentMultiFabIterator rho_halfmfi(U_oldmfi, (*rho_half));
        int i = U_oldmfi.index();
        
        // get the forcing terms
        getGradP(Gp,     i,0,            pres_prev_time);
        getForce(tforces,i,0,Xvel,BL_SPACEDIM,half_time);
        FArrayBox& Rh = rho_halfmfi();
        
#if 1   // do following only at initial iteration--per JBB
        if (initial_iter && is_diffusive[Xvel]) {
            tforces.setVal(0.); 
        }
#endif
        assert(grids[U_oldmfi.index()] == U_oldmfi.validbox());
        const Box& grd = U_oldmfi.validbox();
        godunov->Add_aofs_tf_gp( U_oldmfi(), U_newmfi(),
                                 Aofsmfi(),  tforces,
                                 Gp,       Rh,
                                 grd,      dt );
    }

#else

    int destComp = 0;
    int rhoNGrow = 0;
    int nComp    = 1;
    FillPatchIterator Rhofpi(*this, U_old, rhoNGrow, destComp, half_time,
                             State_Type, Density, nComp);
    MultiFab &P_old = get_old_data(Press_Type);
    int presNGrow   = 0;
    int presSrcComp = 0;
    FillPatchIterator pressurefpi(*this, P_old, presNGrow, destComp, pres_prev_time,
                           Press_Type, presSrcComp, nComp);
    for(;
        Rhofpi.isValid()     &&
        pressurefpi.isValid();
        ++Rhofpi,
        ++pressurefpi)
    {
        DependentMultiFabIterator U_oldmfi(Rhofpi, U_old);
        DependentMultiFabIterator U_newmfi(Rhofpi, U_new);
        DependentMultiFabIterator Aofsmfi(Rhofpi, Aofs);
        DependentMultiFabIterator rho_halfmfi(Rhofpi, (*rho_half));
        int i = U_oldmfi.index();
        
        // get the forcing terms
        //getGradP(Gp,     i,0,            pres_prev_time);
        // from NS::getGradP vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        {
        Box gpbx(grids[i]);
        gpbx.grow(presNGrow);
        Box p_box(surroundingNodes(gpbx));
        //FArrayBox p_fab(p_box,1);
        FArrayBox &p_fab = pressurefpi();
        assert(p_fab.box() == p_box);
        //FillPatch(p_fab,0,time,Press_Type,0,1);
        //-----------------------  size the pressure gradient storage
        //Box gpbx(grd);
        //gpbx.grow(ngrow);
        Gp.resize(gpbx,BL_SPACEDIM);

        //------------------------ test to see if p_fab contains gpbx
        Box test = p_fab.box();
        test = test.enclosedCells();
        assert( test.contains( gpbx ) == true );

        // ----------------------- set pointers
        const int *plo = p_fab.loVect();
        const int *phi = p_fab.hiVect();
        const int *glo = gpbx.loVect();
        const int *ghi = gpbx.hiVect();
        const Real *p_dat  = p_fab.dataPtr();
        const Real *gp_dat = Gp.dataPtr();
        const Real *dx     = geom.CellSize();

        // ----------------------- create the pressure gradient
        FORT_GRADP (p_dat,ARLIM(plo),ARLIM(phi),
                    gp_dat,ARLIM(glo),ARLIM(ghi),glo,ghi,dx);
        }
        // from NS::getGradP ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


        //getForce(tforces,i,0,Xvel,BL_SPACEDIM,half_time);
        // from NS::getForce vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        Real grav = Abs(gravity);
        Box tfbox(grids[i]);
        tfbox.grow(rhoNGrow);
        tforces.resize(tfbox, BL_SPACEDIM);
        for(int dc = 0; dc < BL_SPACEDIM; dc++) {
          int sc = Xvel + dc;
#if (BL_SPACEDIM == 2)
          if (BL_SPACEDIM == 2 && sc == Yvel && grav > 0.001)
#endif
#if (BL_SPACEDIM == 3)
          if (BL_SPACEDIM == 3 && sc == Zvel && grav > 0.001)
#endif
          {
              // set force to -rho*g
            FArrayBox rho(Rhofpi.fabbox());
            //getState(rho,i,1,Density,1,prev_time);
            rho.copy(Rhofpi());
            rho.mult(-grav);
            tforces.copy(rho,0,dc,1);
          } else {
            tforces.setVal(0.0,dc);
          }
        }  // end for(dc...)
        // from NS::getForce ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        FArrayBox& Rh = rho_halfmfi();
        
#if 1   // do following only at initial iteration--per JBB
        if (initial_iter && is_diffusive[Xvel]) {
            tforces.setVal(0.); 
        }
#endif
        assert(grids[U_oldmfi.index()] == U_oldmfi.validbox());
        const Box& grd = U_oldmfi.validbox();
        godunov->Add_aofs_tf_gp( U_oldmfi(), U_newmfi(),
                                 Aofsmfi(),  tforces,
                                 Gp,       Rh,
                                 grd,      dt );
    }
#endif

}

//-------------------------------------------------------------

void NavierStokes::velocity_diffusion_update(Real dt)
{
    // compute the viscous forcing
    // do following except at initial iteration--rbp, per jbb

    if (is_diffusive[Xvel]) {
        MultiFab *delta_rhs;
        diffuse_velocity_setup(dt, delta_rhs);
        diffusion->diffuse_velocity(dt, be_cn_theta, rho_half, 1);
        if (delta_rhs!=NULL) delete delta_rhs;
    }
}
    
// -------------------------------------------------------------

void
NavierStokes::diffuse_velocity_setup(Real dt,
              MultiFab*& delta_rhs)
{
    Real time = state[State_Type].prevTime();

    delta_rhs = NULL;

    if (S_in_vel_diffusion && have_divu) {
      Real time = state[State_Type].prevTime();
      delta_rhs = new MultiFab(grids,BL_SPACEDIM,0,Fab_allocate);
      delta_rhs->setVal(0.0);
      
      MultiFab divmusi(grids,BL_SPACEDIM,0,Fab_allocate);

      diffusion->compute_divmusi(time,visc_coef[Xvel],divmusi);
      divmusi.mult((1./3.)*(1.0-be_cn_theta),0,BL_SPACEDIM,0);
      (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);

      diffusion->compute_divmusi(time+dt,visc_coef[Xvel],divmusi);
      divmusi.mult((1./3.)*be_cn_theta,0,BL_SPACEDIM,0);
      (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);
    }
}

//-------------------------------------------------------------

void NavierStokes::initial_velocity_diffusion_update(Real dt)
{
#if 1
    //int i;
    FArrayBox Gp, tforces;

    // simulation parameters
    MultiFab &U_old     = get_old_data(State_Type);
    MultiFab &U_new     = get_new_data(State_Type);
    MultiFab &Aofs      = *aofs;
    Real prev_time      = state[State_Type].prevTime();
    Real pres_prev_time = state[Press_Type].prevTime();

    // do following only at initial iteration--rbp, per jbb
    if(is_diffusive[Xvel]) {
        
        // get viscous forcing terms
        MultiFab visc_terms(grids,BL_SPACEDIM,1,Fab_allocate);
        getViscTerms(visc_terms,Xvel,BL_SPACEDIM,prev_time);
        if(be_cn_theta==1.0)visc_terms.setVal(0.0,1);

        // update U_new with viscosity
        //for (i = 0; i < grids.length(); i++)
#if (USEOLDFILLPATCH == 1)
        for(MultiFabIterator U_oldmfi(U_old); U_oldmfi.isValid(); ++U_oldmfi) {
            DependentMultiFabIterator U_newmfi(U_oldmfi, U_new);
            DependentMultiFabIterator Aofsmfi(U_oldmfi, Aofs);
            DependentMultiFabIterator rho_halfmfi(U_oldmfi, (*rho_half));
            assert(grids[U_oldmfi.index()] == U_oldmfi.validbox());

            int i = U_oldmfi.index();

            FArrayBox& Rh = rho_halfmfi();
            const Box& grd = U_oldmfi.validbox();

            // get the forcing terms
            getGradP(Gp,     i,0,            pres_prev_time);
            getForce(tforces,i,0,Xvel,BL_SPACEDIM,prev_time);
            godunov->Sum_tf_gp_visc( tforces, visc_terms[i], Gp, Rh );
            
            // compute the inviscid update
            godunov->Add_aofs_tf( U_oldmfi(),
                                  U_newmfi(), 0, BL_SPACEDIM, 
                                  Aofsmfi(),  0,
                                  tforces,  0,
                                  grd,      dt );
        }  // end MultiFabIterator
#else

    int destComp = 0;
    int rhoNGrow = 0;
    int nComp    = 1;
    FillPatchIterator Rhofpi(*this, U_old, rhoNGrow, destComp, prev_time,
                             State_Type, Density, nComp);
    MultiFab &P_old = get_old_data(Press_Type);
    int presNGrow   = 0;
    int presSrcComp = 0;
    FillPatchIterator pressurefpi(*this, P_old, presNGrow, destComp, pres_prev_time,
                           Press_Type, presSrcComp, nComp);
    for(;
        Rhofpi.isValid()     &&
        pressurefpi.isValid();
        ++Rhofpi,
        ++pressurefpi)
    {
            DependentMultiFabIterator U_oldmfi(Rhofpi, U_old);
            DependentMultiFabIterator U_newmfi(Rhofpi, U_new);
            DependentMultiFabIterator Aofsmfi(Rhofpi, Aofs);
            DependentMultiFabIterator rho_halfmfi(Rhofpi, (*rho_half));
            DependentMultiFabIterator visc_termsmfi(Rhofpi, visc_terms);
            assert(grids[U_oldmfi.index()] == U_oldmfi.validbox());

            int i = U_oldmfi.index();

            FArrayBox& Rh = rho_halfmfi();
            const Box& grd = U_oldmfi.validbox();

            // get the forcing terms
            //getGradP(Gp,     i,0,            pres_prev_time);
        // from NS::getGradP vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        {
        Box gpbx(grids[i]);
        gpbx.grow(presNGrow);
        Box p_box(surroundingNodes(gpbx));
        //FArrayBox p_fab(p_box,1);
        FArrayBox &p_fab = pressurefpi();
        assert(p_fab.box() == p_box);
        //FillPatch(p_fab,0,time,Press_Type,0,1);
        //-----------------------  size the pressure gradient storage
        //Box gpbx(grd);
        //gpbx.grow(ngrow);
        Gp.resize(gpbx,BL_SPACEDIM);

        //------------------------ test to see if p_fab contains gpbx
        Box test = p_fab.box();
        test = test.enclosedCells();
        assert( test.contains( gpbx ) == true );

        // ----------------------- set pointers
        const int *plo = p_fab.loVect();
        const int *phi = p_fab.hiVect();
        const int *glo = gpbx.loVect();
        const int *ghi = gpbx.hiVect();
        const Real *p_dat  = p_fab.dataPtr();
        const Real *gp_dat = Gp.dataPtr();
        const Real *dx     = geom.CellSize();

        // ----------------------- create the pressure gradient
        FORT_GRADP (p_dat,ARLIM(plo),ARLIM(phi),
                    gp_dat,ARLIM(glo),ARLIM(ghi),glo,ghi,dx);
        }
        // from NS::getGradP ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


            //getForce(tforces,i,0,Xvel,BL_SPACEDIM,prev_time);
        // from NS::getForce vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        Real grav = Abs(gravity);
        Box tfbox(grids[i]);
        tfbox.grow(rhoNGrow);
        tforces.resize(tfbox, BL_SPACEDIM);
        for(int dc = 0; dc < BL_SPACEDIM; dc++) {
          int sc = Xvel + dc;
#if (BL_SPACEDIM == 2)
          if (BL_SPACEDIM == 2 && sc == Yvel && grav > 0.001)
#endif
#if (BL_SPACEDIM == 3)
          if (BL_SPACEDIM == 3 && sc == Zvel && grav > 0.001)
#endif
          {
              // set force to -rho*g
            FArrayBox rho(Rhofpi.fabbox());
            //getState(rho,i,1,Density,1,prev_time);
            rho.copy(Rhofpi());
            rho.mult(-grav);
            tforces.copy(rho,0,dc,1);
          } else {
            tforces.setVal(0.0,dc);
          }
        }  // end for(dc...)
        // from NS::getForce ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


            godunov->Sum_tf_gp_visc( tforces, visc_termsmfi(), Gp, Rh );
            
            // compute the inviscid update
            godunov->Add_aofs_tf( U_oldmfi(),
                                  U_newmfi(), 0, BL_SPACEDIM, 
                                  Aofsmfi(),  0,
                                  tforces,  0,
                                  grd,      dt );
        }  // end MultiFabIterator
#endif
    }  // end if(is_diffusive[Xvel])
#endif
}



//=================================================================
// Diagnostics and IO functions follow
//=================================================================

// -------------------------------------------------------------
void NavierStokes::errorEst(TagBoxArray &tags, int clearval, int tagval,
                            Real time)
{
    const Box &domain = geom.Domain();
    const int *domain_lo = domain.loVect();
    const int *domain_hi = domain.hiVect();
    const Real *dx = geom.CellSize();
    const Real *prob_lo = geom.ProbLo();

//  original code vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

/*
    for(FabArrayIterator<int, TagBox> tagsfai(tags); tagsfai.isValid(); ++tagsfai) {
        TagBox &tn = tagsfai();
        //assert(grids[tagsfai.index()] == tagsfai.validbox());
        int i = tagsfai.index();

        int *tptr = tn.dataPtr();
        const Box& tbox = tn.box();
        const int* tlo = tbox.loVect();
        const int* thi = tbox.hiVect();

        const int* lo = grids[i].loVect();
        const int* hi = grids[i].hiVect();
        const Real* xlo = grid_loc[i].lo();

          // loop over error estimation quantities, derive quantity
          // then call user supplied error tagging function
        for (int j = 0; j < err_list.length(); j++) {
            const ErrorRec *err = err_list[j];
            int ngrow = err->nGrow();
            Box bx(grow(grids[i],ngrow));
            FArrayBox *dfab = derive(bx,err->name(),time);
            Real* dat = dfab->dataPtr();
            const int* dat_lo = dfab->loVect();
            const int* dat_hi = dfab->hiVect();
            int ncomp = dfab->nComp();

            err->errFunc()(tptr,ARLIM(tlo),ARLIM(thi),&tagval,&clearval,
                           dat, ARLIM(dat_lo), ARLIM(dat_hi), 
                           lo, hi, &ncomp,
                           domain_lo, domain_hi,
                           dx, xlo, prob_lo, &time, &level);

            delete dfab;
        }
    }
//  end original code ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*/

//  new code vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    // loop over error estimation quantities, derive quantity
    // then call user supplied error tagging function
    for(int j = 0; j < err_list.length(); j++) {
      const ErrorRec *err = err_list[j];
      int ngrow = err->nGrow();
      int state_index, src_comp;
      if( ! isStateVariable(err->name(),state_index, src_comp)) {
        cerr << "Error in NavierStokes::errorEst:  " << err->name() << " is not"
             << " a state variable:  fix for parallel." << endl;
        ParallelDescriptor::Abort("Exiting.");
      }

      MultiFab &dataMF = get_new_data(state_index);
      assert(tags.length() == dataMF.length());

      int destComp = 0;
      int nComp = 1;
      for(FillPatchIterator fpi(*this, dataMF, ngrow, destComp, time,
                                state_index, src_comp, nComp);
          fpi.isValid();
          ++fpi)
      {
        int i = fpi.index();
        //DependentFabArrayIterator<int, TagBox> tagsfai(fpi, tags);
        //TagBox &tn = tagsfai();
        //assert(tags[i] != NULL);
        TagBox &tn = tags[i];

        int *tptr = tn.dataPtr();
        const Box &tbox = tn.box();
        const int *tlo = tbox.loVect();
        const int *thi = tbox.hiVect();

        const int *lo = grids[i].loVect();
        const int *hi = grids[i].hiVect();
        const Real *xlo = grid_loc[i].lo();

        FArrayBox &dfab = fpi();
        assert(dfab.box() == grow(grids[i],ngrow));
        Real *dat = dfab.dataPtr();
        const int *dat_lo = dfab.loVect();
        const int *dat_hi = dfab.hiVect();
        int ncomp = dfab.nComp();

        err->errFunc()(tptr,ARLIM(tlo),ARLIM(thi),&tagval,&clearval,
                       dat, ARLIM(dat_lo), ARLIM(dat_hi),
                       lo, hi, &ncomp,
                       domain_lo, domain_hi,
                       dx, xlo, prob_lo, &time, &level);
        //delete dfab;
      }
    }

}  // end errorEst(...)


// ---------------------------------------------------------------
Real NavierStokes::sumDerive(const aString& name, Real time)
{
if(ParallelDescriptor::NProcs() > 1) {
    ParallelDescriptor::Abort("NavierStokes::sumDerive(...) not implemented in parallel.");
} else {
    cerr << "NavierStokes::sumDerive(...) not implemented in parallel.\n";
}
    Real sum = 0.0;
    int finest_level = parent->finestLevel();
    int i;
    for (i = 0; i < grids.length(); i++) {
        FArrayBox* stuff = derive(grids[i],name,time);
        if (level < finest_level) {
            const BoxArray& f_box = parent->boxArray(level+1);
            for (int j = 0; j < f_box.length(); j++) {
                Box c_box = coarsen(f_box[j],fine_ratio);
                c_box &= grids[i];
                if (c_box.ok()) stuff->setVal(0.0,c_box,0);
            }
        }
        sum += stuff->sum(0);
        delete stuff;
    }
    ParallelDescriptor::ReduceRealSum(sum);
    return sum;
}

// ---------------------------------------------------------------
Real  NavierStokes::volWgtSum(const aString& name, Real time)
{
/*  original code vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    Real sum = 0.0;
    int finest_level = parent->finestLevel();
    int rz_flag = (CoordSys::IsRZ() ? 1 : 0);
    const Real* dx = geom.CellSize();
    int i;
    for (i = 0; i < grids.length(); i++) {
        FArrayBox* stuff = derive(grids[i],name,time);
        if (level < finest_level) {
            const BoxArray& f_box = parent->boxArray(level+1);
            for (int j = 0; j < f_box.length(); j++) {
                Box c_box = coarsen(f_box[j],fine_ratio);
                c_box &= grids[i];
                if (c_box.ok()) stuff->setVal(0.0,c_box,0);
            }
        }
        Real s;
        const Real* dat = stuff->dataPtr();
        const int* dlo = stuff->loVect();
        const int* dhi = stuff->hiVect();
        const int* lo = grids[i].loVect();
        const int* hi = grids[i].hiVect();
        int nz = hi[1]-lo[1]+1;
        Real *tmp = new Real[nz];

#if (BL_SPACEDIM == 2)
        int irlo = lo[0]-radius_grow;
        int irhi = hi[0]+radius_grow;
        Real *rad = &radius[i];
        FORT_SUMMASS(dat,ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),
                     dx,&s,rad,&irlo,&irhi,&rz_flag,tmp);
#endif
#if (BL_SPACEDIM == 3)
        FORT_SUMMASS(dat,ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),dx,&s,tmp);
#endif
        sum += s;
        delete stuff;
        delete tmp;
    }
    return sum;
*/


    int stateIndex, srcComp;
    if(isStateVariable(name,stateIndex,srcComp)) {
      // do nothing
    } else {
      ParallelDescriptor::Abort("Error in NS::volWgtSum:  bad variable.");
    }
    int boxGrow = 0;
    int destComp = 0;
    int nComp = 1;
    MultiFab &dataMF = get_data(stateIndex, time);

    Real sum = 0.0;
    int finest_level = parent->finestLevel();
    int rz_flag = (CoordSys::IsRZ() ? 1 : 0);
    const Real *dx = geom.CellSize();
    //for(int i = 0; i < grids.length(); i++)
    for(FillPatchIterator fpi(*this, dataMF, boxGrow, destComp, time,
                              stateIndex, srcComp, nComp); fpi.isValid(); ++fpi)
    {
      int i = fpi.index();
      assert(grids[i] == fpi.validbox());
        FArrayBox &derivedFabTemp = fpi();
        FArrayBox *derivedFab = new FArrayBox(grids[i], nComp);
        derivedFab->copy(derivedFabTemp);
        if(level < finest_level) {
            const BoxArray& f_box = parent->boxArray(level+1);
            for(int j = 0; j < f_box.length(); j++) {
                Box c_box = coarsen(f_box[j],fine_ratio);
                c_box &= grids[i];
                if(c_box.ok()) {
                  derivedFab->setVal(0.0,c_box,0);
                }
            }
        }
        Real s;
        const Real* dat = derivedFab->dataPtr();
        const int* dlo = derivedFab->loVect();
        const int* dhi = derivedFab->hiVect();
        const int* lo = grids[i].loVect();
        const int* hi = grids[i].hiVect();
        int nz = hi[1]-lo[1]+1;
        Real *tmp = new Real[nz];

#if (BL_SPACEDIM == 2)
        int irlo = lo[0]-radius_grow;
        int irhi = hi[0]+radius_grow;
        Real *rad = &radius[i];
        FORT_SUMMASS(dat,ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),
                     dx,&s,rad,&irlo,&irhi,&rz_flag,tmp);
#endif
#if (BL_SPACEDIM == 3)
        FORT_SUMMASS(dat,ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),dx,&s,tmp);
#endif
        sum += s;
        delete derivedFab;
        delete tmp;
    }
    ParallelDescriptor::ReduceRealSum(sum);
    return sum;
}

#ifdef BL_PARALLEL_IO
aString
NavierStokes::thePlotFileType () const
{
    //
    // Increment this whenever the writePlotFile() format changes.
    //
    static const aString the_plot_file_type("NavierStokes-V1.0");

    return the_plot_file_type;
}

void
NavierStokes::writePlotFile (const aString& dir,
                             ostream&       os,
                             VisMF::How     how)
{
    int i, n;
    List<DeriveRec*> dlist = derive_lst.dlist();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';
        //
        // Only write out velocity and scalar data.
        //
        int n_var = NUM_STATE;
        int n_data_items = n_var + dlist.length();

        if (have_dsdt)
            n_data_items += 2;
        else if (have_divu)
            n_data_items += 1;

        os << n_data_items << '\n';

        for (n = 0; n < NUM_STATE; n++)
            os << desc_lst[State_Type].name(n) << '\n';

        if (have_divu)
        {
            os << desc_lst[Divu_Type].name(0) << '\n';
            if (have_dsdt)
                os << desc_lst[Dsdt_Type].name(0) << '\n';
        }

        for (ListIterator<DeriveRec*> lidrp(dlist); lidrp; ++lidrp)
            os << dlist[lidrp]->name() << '\n';

        os << BL_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (i = 0; i < BL_SPACEDIM; i++) os << Geometry::ProbLo(i) << ' ';
        os << '\n';
        for (i = 0; i < BL_SPACEDIM; i++) os << Geometry::ProbHi(i) << ' ';
        os << '\n';
        for (i = 0; i < f_lev; i++) os << parent->refRatio(i)[0] << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++) os << parent->Geom(i).Domain() << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++) os << parent->levelSteps(i) << ' ';
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
    int ngrids = grids.length();
    Real cur_time = state[State_Type].curTime();
    MultiFab& cell_dat = state[State_Type].newData();
    MultiFab* divu_dat = 0;
    MultiFab* dsdt_dat = 0;

    if (have_divu)
    {
        divu_dat = &(state[Divu_Type].newData());
        if (have_dsdt)
            dsdt_dat = &(state[Dsdt_Type].newData());
    }

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << ngrids << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (i = 0; i < cell_dat.boxArray().length(); ++i)
        {
            for (n = 0; n < BL_SPACEDIM; n++)
                os << grid_loc[i].lo(n) << ' ' << grid_loc[i].hi(n) << '\n';
        }
    }
    //
    // There may be up to three MultiFab written out at each level,
    // not including any derived-type MultiFabs.  These are the ones with
    // hard-wired names.  The names of the derived-type MultiFabs will
    // be the names of the derived types themselves.
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

    if (!Utility::CreateDirectory(FullPath, 0755))
        Utility::CreateDirectoryFailed(FullPath);

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        aString PathNameInHeader = Level;
        PathNameInHeader += BaseName[0];
        os << PathNameInHeader << '\n';

        if (have_divu)
        {
            PathNameInHeader = Level;
            PathNameInHeader += BaseName[1];
            os << PathNameInHeader << '\n';

            if (have_dsdt)
            {
                PathNameInHeader = Level;
                PathNameInHeader += BaseName[2];
                os << PathNameInHeader << '\n';
            }
        }
        //
        // Don't forget the derived types.
        //
        for (ListIterator<DeriveRec*> lidrp(dlist); lidrp; ++lidrp)
        {
            PathNameInHeader = Level;
            //
            // The derived-type names aren't prefixed with a '/'
            //
            PathNameInHeader += '/';
            PathNameInHeader += dlist[lidrp]->name();
            os << PathNameInHeader << '\n';
        }
    }
    //
    // Use the Full pathname when naming the MultiFab.
    //
    aString TheFullPath = FullPath;
    TheFullPath += BaseName[0];
    RunStats::addBytes(VisMF::Write(cell_dat, TheFullPath, how));

    if (have_divu)
    {
        TheFullPath = FullPath;
        TheFullPath += BaseName[1];
        RunStats::addBytes(VisMF::Write(*divu_dat, TheFullPath, how));

        if (have_dsdt)
        {
            TheFullPath = FullPath;
            TheFullPath += BaseName[2];
            RunStats::addBytes(VisMF::Write(*dsdt_dat, TheFullPath, how));
        }
    }

    for (ListIterator<DeriveRec*> lidrp(dlist); lidrp; ++lidrp)
    {
        TheFullPath = FullPath;
        TheFullPath += dlist[lidrp]->name();

        MultiFab* mf = derive(dlist[lidrp]->name(), cur_time);

        RunStats::addBytes(VisMF::Write(*mf, TheFullPath, how));

        delete mf;
    }
}

#else

void
NavierStokes::writePlotFile (ostream& os)
{
    //
    // This is a hack to mimic the old Hyperbolic style plotfile so
    // we can use the graphics tools developed for that code.
    // The new version of this code will have a better interface.
    //
    int i, n;
    List<DeriveRec*> dlist = derive_lst.dlist();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // Only write out velocity and scalar data.
        //
        int n_var = NUM_STATE;
        int n_data_items = n_var + dlist.length();

        if (have_dsdt)
            n_data_items += 2;
        else if (have_divu)
            n_data_items += 1;

        os << n_data_items << '\n';

        for (n = 0; n < NUM_STATE; n++)
            os << desc_lst[State_Type].name(n) << '\n';

        if (have_divu)
        {
            os << desc_lst[Divu_Type].name(0) << '\n';
            if (have_dsdt)
                os << desc_lst[Dsdt_Type].name(0) << '\n';
        }

        for (ListIterator<DeriveRec*> lidrp(dlist); lidrp; ++lidrp)
            os << dlist[lidrp]->name() << '\n';

        os << BL_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (i = 0; i < BL_SPACEDIM; i++) os << Geometry::ProbLo(i) << ' ';
        os << '\n';
        for (i = 0; i < BL_SPACEDIM; i++) os << Geometry::ProbHi(i) << ' ';
        os << '\n';
        for (i = 0; i < f_lev; i++) os << parent->refRatio(i)[0] << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++) os << parent->Geom(i).Domain() << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++) os << parent->levelSteps(i) << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++)
                os << parent->Geom(i).CellSize()[k] << ' ';
            os << '\n';
        }
        os << (int) CoordSys::Coord() << '\n';
        //
        // Write bndry data.
        //
        os << "0\n";
    }
    //
    // Now write state data.
    //
    int ngrids = grids.length();
    Real cur_time = state[State_Type].curTime();
    MultiFab& cell_dat = state[State_Type].newData();
    MultiFab* divu_dat = 0;
    MultiFab* dsdt_dat = 0;
    if (have_divu)
    {
        divu_dat = &(state[Divu_Type].newData());
        if (have_dsdt)
            dsdt_dat = &(state[Dsdt_Type].newData());
    }

    if (ParallelDescriptor::IOProcessor())
        os << level << ' ' << ngrids << '\n';

    streampos filePosition;
    ParallelDescriptor::ShareVar(&filePosition, sizeof(streampos));
    ParallelDescriptor::Synchronize();
    int myproc = ParallelDescriptor::MyProc();

    for (i = 0; i < ngrids; i++)
    {
        int fabProc = cell_dat.DistributionMap()[i];
        if (fabProc == myproc)
        {
            const Box& grd = grids[i];
            os << grd << '\n'
               << level << '\n'
               << parent->levelSteps(level) << '\n'
               << cur_time << '\n';
            for (n = 0; n < BL_SPACEDIM; n++)
                os << grid_loc[i].lo(n) << ' ' << grid_loc[i].hi(n) << '\n';
            FArrayBox dat(grd,1);
            for (n = 0; n < NUM_STATE; n++)
            {
                dat.copy(cell_dat[i],n,0,1);
                dat.writeOn(os,0,1);
            }
            if (have_divu)
            {
                dat.copy((*divu_dat)[i],0,0,1);
                dat.writeOn(os,0,1);
                if (have_dsdt)
                {
                    dat.copy((*dsdt_dat)[i],0,0,1);
                    dat.writeOn(os,0,1);
                }
            }
            for (ListIterator<DeriveRec*> lidrp(dlist); lidrp; ++lidrp)
            {
                FArrayBox *dfab = derive(grd,dlist[lidrp]->name(),cur_time);
                dfab->writeOn(os,0,1);
                delete dfab;
            }
            filePosition = os.tellp();
        }
        ParallelDescriptor::Broadcast(fabProc, &filePosition, &filePosition);
        os.seekp(filePosition);
    }
    ParallelDescriptor::Synchronize();
    ParallelDescriptor::UnshareVar(&filePosition);
}
#endif /*BL_PARALLEL_IO*/

Real
NavierStokes::estTimeStep ()
{
    if (fixed_dt > 0.0) {
      Real factor;
      if (level == 0) {
        factor = 1.0;
      } else {
        int ratio = 1;
        for (int lev = 1; lev <= level; lev++) {
          ratio *= parent->MaxRefRatio(lev-1);
        }
        factor = 1.0/((double) ratio);
      }
      return factor*fixed_dt;
    }

    int n_grow = 0;
    Real cur_time       = state[State_Type].curTime();

    FArrayBox Gp, Rho, tforces;

    MultiFab &P_new = get_new_data(Press_Type);

    MultiFab &U_new = get_new_data(State_Type);
    const Real* dx = geom.CellSize();

    Real gr_max[BL_SPACEDIM], u_max[BL_SPACEDIM];
    int k;
    for (k = 0; k < BL_SPACEDIM; k++) {
        u_max[k] = 0.0;
    }
    Real estdt = 1.0e+20;

#if (USEOLDFILLPATCH == 1)

    for (int i = 0; i < grids.length(); i++) {

        // get the pressure
        Box p_box(surroundingNodes(grids[i]));
        FArrayBox p_fab(p_box,1);
        p_fab.copy(P_new[i],p_box);

        // get the density and and velocities
        getState(Rho,i,n_grow,Density,1,cur_time);
        FArrayBox& fab_vel = U_new[i];

        // get the velocity forcing
        // for some reason no viscous forcing
        getForce(tforces,i,n_grow,Xvel,BL_SPACEDIM,cur_time);
        getGradP( p_fab, Gp, grids[i], 0 );
        tforces.minus(Gp, 0,0,BL_SPACEDIM );

        // estimate the maximum allowable timestep from the Godunov box
        const Real *dx = geom.CellSize();
        Real dt = godunov->estdt( fab_vel, tforces, Rho,
                                  grids[i], dx, cfl, gr_max );

        for (k = 0; k < BL_SPACEDIM; k++) {
            u_max[k] = Max(u_max[k],gr_max[k]);
        }
        estdt = Min(estdt,dt);
    }

#else

    int destComp = 0;
    FillPatchIterator Ufpi(*this, U_new, n_grow, destComp, cur_time,
                           State_Type, Xvel, BL_SPACEDIM);
    FillPatchIterator Rhofpi(*this, U_new, n_grow, destComp, cur_time,
                           State_Type, Density, 1);
    //FillPatchIterator tforcesfpi(*this, U_new, 1, destComp, cur_time,
                           //State_Type, Xvel, BL_SPACEDIM);
    for(;
        Ufpi.isValid()       &&
        Rhofpi.isValid();
        //tforcesfpi.isValid();
        ++Ufpi,
        ++Rhofpi)
        //++tforcesfpi)
    {
        DependentMultiFabIterator P_newmfi(Ufpi, P_new);
        //DependentMultiFabIterator U_newmfi(Ufpi, U_new);
        //assert(grids[P_newmfi.index()] == P_newmfi.validbox());
        int i = Ufpi.index();

        // get the pressure
        assert(grids[i] == Ufpi.validbox());
        Box p_box(surroundingNodes(grids[i]));
        FArrayBox p_fab(p_box,1);
        p_fab.copy(P_newmfi(),p_box);

        // get the density and and velocities
        //getState(Rho,i,n_grow,Density,1,cur_time);
        FArrayBox &Rho     = Rhofpi();
        FArrayBox &fab_vel = Ufpi();

        // get the velocity forcing
        // for some reason no viscous forcing

        //getForce(tforces,i,n_grow,Xvel,BL_SPACEDIM,cur_time);
        // from NS::getForce vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        Real grav = Abs(gravity);
        Box tfbox(grids[i]);
        tfbox.grow(n_grow);
        tforces.resize(tfbox, BL_SPACEDIM);
        for(int dc = 0; dc < BL_SPACEDIM; dc++) {
          int sc = Xvel + dc;
#if (BL_SPACEDIM == 2)
          if (BL_SPACEDIM == 2 && sc == Yvel && grav > 0.001)
#endif
#if (BL_SPACEDIM == 3)
          if (BL_SPACEDIM == 3 && sc == Zvel && grav > 0.001)
#endif
          {
              // set force to -rho*g
            FArrayBox rho(Rho.box());
            //getState(rho,i,1,Density,1,cur_time);
            rho.copy(Rho);
            rho.mult(-grav);
            tforces.copy(rho,0,dc,1);
          } else {
            tforces.setVal(0.0,dc);
          }
        }  // end for(dc...)
        // from NS::getForce ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


        //getGradP( p_fab, Gp, grids[i], 0 );
        // from NS::getGradP vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        {
        Box gpbx(grids[i]);
        //gpbx.grow(0);
        Gp.resize(gpbx,BL_SPACEDIM);

        //------------------------ test to see if p_fab contains gpbx
        Box test = p_fab.box();
        test = test.enclosedCells();
        assert( test.contains( gpbx ) == true );

        // ----------------------- set pointers
        const int *plo = p_fab.loVect();
        const int *phi = p_fab.hiVect();
        const int *glo = gpbx.loVect();
        const int *ghi = gpbx.hiVect();
        const Real *p_dat  = p_fab.dataPtr();
        const Real *gp_dat = Gp.dataPtr();
        const Real *dx     = geom.CellSize();

        // ----------------------- create the pressure gradient
        FORT_GRADP (p_dat,ARLIM(plo),ARLIM(phi),
                    gp_dat,ARLIM(glo),ARLIM(ghi),glo,ghi,dx);
        }

        // from NS::getGradP ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



        tforces.minus(Gp, 0,0,BL_SPACEDIM );

        // estimate the maximum allowable timestep from the Godunov box
        const Real *dx = geom.CellSize();
        Real dt = godunov->estdt( fab_vel, tforces, Rho,
                                  grids[i], dx, cfl, gr_max );

        for (k = 0; k < BL_SPACEDIM; k++) {
            u_max[k] = Max(u_max[k],gr_max[k]);
        }
        estdt = Min(estdt,dt);
    }  // end for(MultiFabIterator ...)

#endif

    // parallel reductions
    for(k = 0; k < BL_SPACEDIM; k++) {
      ParallelDescriptor::ReduceRealMax(u_max[k]);
    }
    ParallelDescriptor::ReduceRealMin(estdt);

    if (verbose) {
        cout << "estTimeStep :: \n";
        cout << "LEV = " << level << " UMAX = ";
        for (k = 0; k < BL_SPACEDIM; k++) {
            cout << u_max[k] << "  ";
        }
        cout << NL;
    }
    Real vel_max = u_max[0];
    for (k = 1; k < BL_SPACEDIM; k++) {
        vel_max = Max(vel_max,u_max[k]);
    }
    if (vel_max < 1.0e-10) {
        Real grav = Abs(gravity);
        if (grav > 1.0e-5) {
            vel_max = grav;
        } else {
            vel_max = 1.0;
        }
        Real dx_min = dx[0];
        for (k = 1; k < BL_SPACEDIM; k++) {
            dx_min = Min(dx_min,dx[k]);
        }
        estdt = cfl*dx_min/vel_max;
    }

    return estdt;
}

// -------------------------------------------------------------
Real
NavierStokes::initialTimeStep()
{
    return init_shrink*estTimeStep();
}

// -------------------------------------------------------------
void NavierStokes::computeNewDt(int finest_level, int sub_cycle,
                           Array<int>& n_cycle,
                           const Array<IntVect>& ref_ratio,
                           Array<Real>& dt_min,
                           Array<Real>& dt_level, Real stop_time)
{
      // we are at the end of a coarse grid timecycle
      // compute the timesteps for the next iteration

    if (level > 0) return;

      // for Navier Stokes, we compute the new dt based on the
      // current velocity field.

    int max_level=parent->maxLevel();

    int i;
    n_cycle[0] = 1;
    for (i = 1; i <= max_level; i++) {
        if(sub_cycle) {
            n_cycle[i] = parent->MaxRefRatio(i-1);
        } else {
            n_cycle[i] = 1;
        }
    }

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (i = 0; i <= finest_level; i++) {
        NavierStokes& ns_level = getLevel(i);
        dt_min[i] = ns_level.estTimeStep();
        dt_min[i] = Min(dt_min[i],change_max*dt_level[i]);
        n_factor *= n_cycle[i];
        dt_0 = Min(dt_0,n_factor*dt_min[i]);
    }

    Real eps = 0.0001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if ( (cur_time + dt_0) > (stop_time - eps) ) dt_0 = stop_time - cur_time;

    // adjust the time step to be able to output checkpoints at specific times
    int a,b;
    Real check_per = parent->checkPer();
    if ( check_per > 0.0 ) {
        a = (cur_time + eps ) / check_per;
        b = (cur_time + dt_0) / check_per;
        if (a != b) dt_0 = b * check_per - cur_time;
    }

    // adjust the time step to be able to output plot files at specific times
    Real plot_per = parent->plotPer();
    if ( plot_per > 0.0 ) {
        a = (cur_time + eps ) / plot_per;
        b = (cur_time + dt_0) / plot_per;
        if (a != b) dt_0 = b * plot_per - cur_time;
    }

    n_factor = 1;
    for (i = 0; i <= max_level; i++) {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/( (float)n_factor );
    }

}

// -------------------------------------------------------------
void NavierStokes::computeInitialDt(int finest_level, int sub_cycle,
                               Array<int>& n_cycle,
                               const Array<IntVect>& ref_ratio,
                               Array<Real>& dt_level)
{
      // grids have been constructed, compute dt for all levels
    if (level > 0) return;

   // for Navier Stokes, we compute the new dt based on the
   // current velocity field.

    int max_level=parent->maxLevel();

    int i;
    n_cycle[0] = 1;
    for (i = 1; i <= max_level; i++) {
        if(sub_cycle) {
            n_cycle[i] = parent->MaxRefRatio(i-1);
        } else {
            n_cycle[i] = 1;
        }
    }

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (i = 0; i <= finest_level; i++) {
        NavierStokes& ns_level = getLevel(i);
        dt_level[i] = ns_level.initialTimeStep();
        n_factor *= n_cycle[i];
        dt_0 = Min(dt_0,n_factor*dt_level[i]);
    }

    n_factor = 1;
    for (i = 0; i <= max_level; i++) {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/( (float)n_factor );
    }

}

//-------------------------------------------------------------
// this function estimates the initial timesteping used by the model
//-------------------------------------------------------------
void NavierStokes::post_init_estDT( Real &dt_init,
                                    Array<int> &nc_save,
                                    Array<Real> &dt_save )
{
    int k;
    
    Real strt_time = state[State_Type].curTime();
    int finest_level = parent->finestLevel();
    const Array<IntVect>& ref_ratio = parent->refRatio();
    Array<Real> dt_level(finest_level+1);
    Array<int>  n_cycle(finest_level+1);

    dt_init = 1.0e+100;
    for (k = 0; k <= finest_level; k++) {
        nc_save[k] = parent->nCycle(k);
        dt_save[k] = getLevel(k).initialTimeStep();
        dt_init = Min(dt_init,dt_save[k]);
    }
    for (k = 0; k <= finest_level; k++) {
        dt_level[k] = dt_init;
        n_cycle[k] = 1;
    }
    Real dt0 = dt_save[0];
    int  n_factor = 1;
    for (k = 0; k <= finest_level; k++) {
        n_factor *= nc_save[k];
        dt0 = Min(dt0,n_factor*dt_save[k]);
    }

    n_factor = 1;
    for (k = 0; k <= finest_level; k++) {
        n_factor *= nc_save[k];
        dt_save[k] = dt0/( (float) n_factor);
    }

      // hack
    parent->setDtLevel(dt_level);
    parent->setNCycle(n_cycle);
    for (k = 0; k <= finest_level; k++) {
        getLevel(k).setTimeLevel(strt_time,dt_init,dt_init);
    }
}

// fills in amrLevel okToContinue
int NavierStokes::okToContinue()
{
    if (level > 0) return true;
    Real dt = parent->dtLevel(0);
    return (dt > dt_cutoff);
}


//=================================================================
// THE MAIN HOOKS INTO AMR AND AMRLEVEL
//=================================================================

// -------------------------------------------------------------
// integration cycle on fine level grids is complete 
// post_timestep is responsible for syncing levels together
//
// The registers used for level syncing are initialized in the
// coarse level advance and incremented in the fine level advance.
// These quantities are described in comments above advance_setup.
// -------------------------------------------------------------
void NavierStokes::post_timestep()
{
    int finest_level = parent->finestLevel();
    
    // reflux 
    if (do_reflux && (level < finest_level)) {
        reflux();
    }

    // average down
    if (level < finest_level) {
        avgDown();
        Real dt = parent->dtLevel(level);
        Real prev_time = state[State_Type].prevTime();
        Real half_time = prev_time + 0.5*dt;
#if (USEOLDFILLPATCH == 1)
        for (int k = 0; k < grids.length(); k++) {
          getState((*rho_half)[k],k,1,Density,1,half_time);
        }
#else
        int boxGrow  = 1;
        int destComp = 0;
        int srcComp  = Density;
        int nComp    = 1;
        for(FillPatchIterator rho_halffpi(*this, *rho_half, boxGrow, destComp,
                                          half_time, State_Type, srcComp, nComp
                                          /* , mapper = NULL */);
            rho_halffpi.isValid();
            ++rho_halffpi)
        {
          DependentMultiFabIterator rho_halfdest(rho_halffpi, (*rho_half));

          // copy from fillpatched fab to rho_half
          rho_halfdest().copy(rho_halffpi());
        }
#endif
    }

    // Mac Sync Correction
    if (do_mac_proj && (level < finest_level)) {
        RunStats mac_sync_stats("mac_sync", level);
        mac_sync_stats.start();
        mac_sync();
        mac_sync_stats.end();
    }
    
    // level Sync Correction
    if (do_sync_proj && (level < finest_level)) {
        level_sync();
    }

    // Test for conservation
    if (level == 0) {
        int nstep = parent->levelSteps(0);
        if ((sum_interval > 0) && (nstep%sum_interval == 0) ) {
            sum_integrated_quantities();
        }
    }
    
}



// build any additional data structures after restart
void NavierStokes::post_restart()
{
}


// build any additional data structures after regrid
void NavierStokes::post_regrid(int lbase, int new_finest)
{
    AmrLevel & amr_level = getLevel(level);
    int finest_level = parent->finestLevel();
    if ( projector )
        if (level == lbase) projector->setFinestLevel(new_finest);
}

// ensure state, and pressure are consistent
void NavierStokes::post_init()
{
    // nothing to sync up at level > 0
    if (level > 0)
        return;

    int finest_level = parent->finestLevel();
    Real dt_init     = 0.0;
    Array<Real> dt_save(finest_level+1);
    Array<int>  nc_save(finest_level+1);

    // ensure state is consistent, i.e. velocity field is non-divergent,
    // coarse levels are fine level averages, pressure is zero
    post_init_state();

    // estimate the initial timestepping
    post_init_estDT( dt_init, nc_save, dt_save );

    // initialize the pressure by iterating the initial timestep
    post_init_press( dt_init, nc_save, dt_save );
    
    // compute the initial estimate of conservation
    if (sum_interval > 0) {
        sum_integrated_quantities();
    }

}



//=================================================================
// MULTILEVEL SYNC FUNCTIONS
//=================================================================

// -------------------------------------------------------------
void NavierStokes::initRhoAvg(Real alpha)
{
    MultiFab &S_new = get_new_data(State_Type);
    rho_avg->setVal(0.);
    //int i;
    //for (i = 0; i < grids.length(); i++) {
    for(MultiFabIterator rho_avgmfi(*rho_avg); rho_avgmfi.isValid(); ++rho_avgmfi) {
      DependentMultiFabIterator S_newmfi(rho_avgmfi, S_new);
      assert(grids[S_newmfi.index()] == S_newmfi.validbox());
       rho_avgmfi().copy(S_newmfi(),S_newmfi.validbox(),Density,
                         S_newmfi.validbox(),0,1);
       rho_avgmfi().mult(alpha);
    }
}

// -------------------------------------------------------------
void NavierStokes::incrRhoAvg(Real alpha)
{
    MultiFab &S_new = get_new_data(State_Type);
    //int i;
    //for (i = 0; i < grids.length(); i++)
    for(MultiFabIterator S_newmfi(S_new); S_newmfi.isValid(); ++S_newmfi) {
        DependentMultiFabIterator rho_avgmfi(S_newmfi, (*rho_avg));
        assert(grids[S_newmfi.index()] == S_newmfi.validbox());
        const int* lo = S_newmfi.validbox().loVect();
        const int* hi = S_newmfi.validbox().hiVect();
        const FArrayBox &rho = S_newmfi();
        const int* rlo = rho.loVect();
        const int* rhi = rho.hiVect();
        const Real* rhodat = rho.dataPtr(Density);
        FArrayBox& avg = rho_avgmfi();
        const int* alo = avg.loVect();
        const int* ahi = avg.hiVect();
        Real* avgdat = avg.dataPtr();
        FORT_INCRMULT(avgdat,ARLIM(alo),ARLIM(ahi),
                      rhodat,ARLIM(rlo),ARLIM(rhi),
                      lo,hi,&alpha);
    }
}

// -------------------------------------------------------------
void NavierStokes::incrPAvg( int iteration, Real alpha)
{
    // set p_avg to zero for the initial iteration
    if (iteration == 1) p_avg->setVal(0.);

    // Then increment p_avg with fine grid pressure
    const BoxArray& P_grids = state[Press_Type].boxArray();
    MultiFab &P_new = get_new_data(Press_Type);

    //int i;
    //for (i = 0; i < P_grids.length(); i++) {
    for(MultiFabIterator P_newmfi(P_new); P_newmfi.isValid(); ++P_newmfi) {
        DependentMultiFabIterator p_avgmfi(P_newmfi, (*p_avg));
        assert(P_grids[P_newmfi.index()] == P_newmfi.validbox());

        const int* lo = P_newmfi.validbox().loVect();
        const int* hi = P_newmfi.validbox().hiVect();
        const FArrayBox &p = P_newmfi();
        const int* p_lo = p.loVect();
        const int* p_hi = p.hiVect();
        const Real* pdat = p.dataPtr();
        FArrayBox& avg = p_avgmfi();
        const int* alo = avg.loVect();
        const int* ahi = avg.hiVect();
        Real* avgdat = avg.dataPtr();
        FORT_INCRMULT(avgdat,ARLIM(alo),ARLIM(ahi),
                      pdat,ARLIM(p_lo),ARLIM(p_hi),
                      lo,hi,&alpha);
    }
}



// -------------------------------------------------------------
// This function ensures that the state is initially consitent
// with respect to the divergence condition and fields are initially consistent
// -------------------------------------------------------------
void NavierStokes::post_init_state()
{
    int k;
    int finest_level = parent->finestLevel();
    
    // do sync project to define divergence free velocity field
    Real divu_time;
    Real pres_time = state[Press_Type].curTime();
    if (have_divu) {
        divu_time = state[Divu_Type].curTime();
    } else {
        divu_time = pres_time;
    }
    if ( projector )
        projector->initialVelocityProject(0,divu_time,have_divu);
    NavierStokes::initial_step = true;
    
    // average velocity and scalar data down from finer levels
    // so that conserved data is consistant between levels
    for (k = finest_level-1; k>= 0; k--) getLevel(k).avgDown();
    
    // zero pressure field
    for (k = 0; k <= finest_level; k++) {
        NavierStokes& ns_level = getLevel(k);
        ns_level.zeroNewPress();
        ns_level.zeroOldPress();
    }
}



// -------------------------------------------------------------
// Initialize the pressure by iterating the initial timestep
// -------------------------------------------------------------
void NavierStokes::post_init_press( Real &dt_init,
                                    Array<int> &nc_save,
                                    Array<Real> &dt_save )
{
    Real strt_time = state[State_Type].curTime();
    int finest_level = parent->finestLevel();
    int iter,k; 
    NavierStokes::initial_iter = true;

    // iterate over the advance function
    for ( iter = 0; iter < init_iter; iter++ ) {

        for ( k = 0; k <= finest_level; k++ ) {
            getLevel(k).advance(strt_time,dt_init,1,1);
        }
          // this constructs a guess at P, also sets
          // p_old == p_new

        MultiFab **sig = new MultiFab*[finest_level+1];
        for (k = 0; k <= finest_level; k++) {
          sig[k] = getLevel(k).rho_half;
        }

        Real dt_crse = parent->dtLevel(0);
        if ( projector )
            projector->initialSyncProject(0,sig,dt_crse,strt_time,
                                      dt_init,have_divu);

        delete sig;

        for (k = finest_level-1; k>= 0; k--) {
            getLevel(k).avgDown();
        }
        for (k = 0; k <= finest_level; k++) {
              // reset state variables to initial time, do not reset
              // pressure variable
            NavierStokes& ns_level = getLevel(k);
            ns_level.resetState(strt_time, dt_save[k], dt_save[k]);
        }
        NavierStokes::initial_iter = false;
    }

    if (init_iter <= 0) NavierStokes::initial_iter = false; // just being 
                                                           // compulsive--rbp 
    NavierStokes::initial_step = false;

      // re-instate timestep
    for (k = 0; k <= finest_level; k++) {
        getLevel(k).setTimeLevel(strt_time,dt_save[k],dt_save[k]);
    }
    parent->setDtLevel(dt_save);
    parent->setNCycle(nc_save);
}




// ==================================================
// interpolate A cell centered Sync correction from a
// coarse level (c_lev) to a fine level (f_lev)
// ==================================================
void NavierStokes::SyncInterp( MultiFab &CrseSync, int c_lev,
                               MultiFab &FineSync, int f_lev, IntVect& ratio,
                               int src_comp, int dest_comp, int num_comp,
                               int increment , Real dt_clev, 
                               int** bc_orig_qty, int which_interp)
{
    int i,n,dir;

    if(ParallelDescriptor::NProcs() > 1) {  // punt
      ParallelDescriptor::Abort("NavierStokes::SyncInterp not implemented in parallel");
    } else {
      cerr << "NavierStokes::SyncInterp not implemented in parallel\n";
    }

    // This routine interpolates the num_comp components of CrseSync
    // (starting at src_comp) and either increments or puts the result into
    // the num_comp components of FineSync (starting at dest_comp)
    // Note: the components of bc_orig_qty corespond to the quantities
    //       of CrseSync
    
    assert (which_interp>=0 && which_interp<=2);
    Interpolater *interpolater=NULL;
    interpolater = (which_interp==0) ? &cell_cons_interp : interpolater;
    interpolater = (which_interp==1) ? &pc_interp : interpolater;
    interpolater = (which_interp==2) ? &unlimited_cc_interp : interpolater;

    FArrayBox cdata,fdata;
    
    // get fine parameters
    NavierStokes& fine_level = getLevel(f_lev);
    const BoxArray& fgrids   = fine_level.boxArray();
    const Geometry& fgeom    = parent->Geom(f_lev);
    int nfine                = fgrids.length();
    const Box& domain        = fgeom.Domain();

    // get coarse parameters
    NavierStokes& crse_level = getLevel(c_lev);
    const BoxArray& cgrids   = crse_level.boxArray();
    const Geometry& cgeom    = parent->Geom(c_lev);
    const Real* dx_crse      = cgeom.CellSize();

    const Box cdomain(coarsen(domain,ratio));
    const int* cdomlo = cdomain.loVect();
    const int* cdomhi = cdomain.hiVect();

    int* bc_new =  new int[2*BL_SPACEDIM*(src_comp+num_comp)];

    // loop over fine grids
    for (i = 0; i < nfine; i++) {

        // create storage for interpolation
        const Box& grd = fgrids[i];
        Box cgrd = interpolater->CoarseBox(grd,ratio);

        fdata.resize(grd,num_comp);
        cdata.resize(cgrd,num_comp);
        cdata.setVal(0.);
        CrseSync.copy(cdata,src_comp,0,num_comp);

        const int* clo = cdata.loVect();
        const int* chi = cdata.hiVect();
        const Real* xlo = fine_level.grid_loc[i].lo();

        for ( n = 0; n < num_comp; n++) {
 
           for ( dir = 0; dir < BL_SPACEDIM; dir++) {
              int bc_index = (n+src_comp)*(2*BL_SPACEDIM) + dir;
              bc_new[bc_index]             = INT_DIR;
              bc_new[bc_index+BL_SPACEDIM] = INT_DIR;
 
              if (clo[dir] < cdomlo[dir] || chi[dir] > cdomhi[dir]) 
                for (int crse = 0; crse < cgrids.length(); crse++) {

                const int* crse_lo = cgrids[crse].loVect();
                const int* crse_hi = cgrids[crse].hiVect();

                if (clo[dir] < cdomlo[dir] && crse_lo[dir] == cdomlo[dir]) {
                  bc_new[bc_index] = 
                    bc_orig_qty[crse][bc_index];
                }
                if (chi[dir] > cdomhi[dir] && crse_hi[dir] == cdomhi[dir]) {
                  bc_new[bc_index+BL_SPACEDIM] = 
                    bc_orig_qty[crse][bc_index+BL_SPACEDIM]; 
                }
             }
           }

           FORT_FILCC(cdata.dataPtr(n), ARLIM(clo), ARLIM(chi),
                      cdomlo, cdomhi, dx_crse, xlo, &(bc_new[2*BL_SPACEDIM*(n+src_comp)]));
        }

        // fill in periodic images
        if ( cgeom.isAnyPeriodic() ) {
            const Box& domain = cgeom.Domain();
            Array<IntVect> pshifts(27);
            cgeom.periodicShift(domain, cgrd, pshifts);

            for( int iiv=0; iiv<pshifts.length(); iiv++) {
                IntVect iv=pshifts[iiv];
                cdata.shift(iv);
                CrseSync.copy(cdata,src_comp,0,num_comp);
                cdata.shift(-iv);
            }
        }

        // set the boundary condition array for interpolation
        Array<BCRec> bc_interp(num_comp);
        for ( n = 0; n < num_comp; n++) {
          for (int dir = 0; dir < BL_SPACEDIM; dir++) {
            int bc_index = (n+src_comp)*(2*BL_SPACEDIM) + dir;
            bc_interp[n].setLo(dir,bc_new[bc_index]);
            bc_interp[n].setHi(dir,bc_new[bc_index+BL_SPACEDIM]);
          }
        }

        // scale coarse interpolant for anelastic
        ScaleCrseSyncInterp( cdata, c_lev, num_comp );
        
        // compute the interpolated correction

        interpolater->interp(cdata,0,
                             fdata,0,num_comp,grd,ratio,
                             cgeom,fgeom,bc_interp);

        // rescale Fine interpolant for anelastic
        reScaleFineSyncInterp( fdata, f_lev, num_comp );

        // set Fine Sync equal to the correction or add it in
        if ( increment ) {
            fdata.mult(dt_clev);
            FineSync[i].plus(fdata,0,dest_comp,num_comp);
        } else {
            FineSync[i].copy(fdata,0,dest_comp,num_comp);
        }
    }

    delete bc_new;
}

// interpolate sync pressure correction to a finer level
void NavierStokes::SyncProjInterp( MultiFab &phi,   int c_lev,
                                   MultiFab &P_new, int f_lev,
                                   IntVect& ratio )
{
    if(ParallelDescriptor::NProcs() > 1) {  // punt
      ParallelDescriptor::Abort("NavierStokes::SyncProjInterp not implemented in parallel");
      // because of phi.copy(crse_phi) within multifab loop below
    } else {
      cerr << "NavierStokes::SyncProjInterp not implemented in parallel\n";
    }
    // get fine parameters
    NavierStokes& fine_level = getLevel(f_lev);
    const BoxArray& fgrids   = fine_level.boxArray();
    const Geometry& fgeom    = parent->Geom(f_lev);
    const BoxArray& P_grids  = P_new.boxArray();
    int nfine                = fgrids.length();

    // get coarse parameters
    NavierStokes& crse_level = getLevel(c_lev);
    const Geometry& cgeom    = parent->Geom(c_lev);

    int i;
    FArrayBox fine_phi, crse_phi;
    Array<BCRec> bc(BL_SPACEDIM);
    
    // loop over the fine grids
    for ( i = 0; i < nfine; i++) {

        // get space for bi/trilinear interpolation
        const Box& fbox = P_grids[i];
        const Box  cbox(node_bilinear_interp.CoarseBox(fbox,ratio));
        fine_phi.resize(fbox,1);
        crse_phi.resize(cbox,1);
        fine_phi.setVal(1.e30);
        crse_phi.setVal(1.e30);
        phi.copy(crse_phi);

        // compute the interpolated pressure sync correction
        node_bilinear_interp.interp(crse_phi,0,
                                    fine_phi,0,1,
                                    fbox,ratio,cgeom,fgeom,bc);

        // add the fine corrections to the state
        P_new[i].plus(fine_phi);
    }
}




//------------------------------------------------------------
// This function averages a multifab of fine data down onto
// a multifab of coarse data
//
//  This should be an Amrlevel or Multifab function
//
//------------------------------------------------------------
void NavierStokes::avgDown( const BoxArray &grids,  const BoxArray &fgrids,
                            MultiFab &S_crse, MultiFab &S_fine,
                            MultiFab &volume, MultiFab &fvolume,
                            int c_level,            int f_level,
                            int strt_comp, int num_comp, IntVect& fratio )
{
    assert( S_crse.nComp() == S_fine.nComp() );
    int num_crse = grids.length();
    int num_fine = fgrids.length();

/*
    // loop over coarse grids, and intersect coarse with fine
    int crse, fine;
    for ( crse = 0; crse < num_crse; crse++) {
        const Box& cbox = grids[crse];
        for ( fine = 0; fine < num_fine; fine++) {
            const Box& fbox = fgrids[fine];
            Box ovlp(coarsen(fbox,fratio));
            ovlp &= cbox;
            if (ovlp.ok()) {
                avgDown( S_fine[fine],  S_crse[crse], 
                         fvolume[fine], volume[crse], 
                         f_level,       c_level,      
                         ovlp, strt_comp, num_comp, fratio );
            }
        }
    }
*/

//
//
//  this needs to be optimized to communicate only the
//  overlap part of the fab
//
//

    MultiFabCopyDescriptor mfcd(true);
    MultiFabId mfidS_fine  = mfcd.RegisterFabArray(&S_fine);
    MultiFabId mfidFineVol = mfcd.RegisterFabArray(&fvolume);

    List<FillBoxId> fillBoxIdList;
    List<FillBoxId> fillBoxIdListVol;

    // loop over coarse grids, and intersect coarse with fine
    for(ConstMultiFabIterator mfi(S_crse); mfi.isValid(); ++mfi) {
        assert(grids[mfi.index()] == mfi.validbox());
        const Box &cbox = mfi.validbox();
        for(int fine = 0; fine < num_fine; fine++) {
            const Box &fbox = fgrids[fine];
            Box ovlp(coarsen(fbox,fratio));
            ovlp &= cbox;
            if(ovlp.ok()) {
              BoxList unfilledBoxes(fbox.ixType());  // unused here
              FillBoxId fbidFine, fbidFineVol;
              const Box &sfineBox = S_fine.fabbox(fine);
              fbidFine = mfcd.AddBox(mfidS_fine, sfineBox, unfilledBoxes,
                                       0, 0, num_comp);

              const Box &fineVolBox = fvolume.fabbox(fine);
              fbidFineVol = mfcd.AddBox(mfidFineVol, fineVolBox, unfilledBoxes,
                                        0, 0, 1);
              fillBoxIdList.append(fbidFine);
              fillBoxIdListVol.append(fbidFineVol);
            }
        }
    }

    mfcd.CollectData();

    ListIterator<FillBoxId> fbidli(fillBoxIdList);
    ListIterator<FillBoxId> fbidlivol(fillBoxIdListVol);

    
    for(ConstMultiFabIterator mfi(S_crse); mfi.isValid(); ++mfi) {
        ConstDependentMultiFabIterator mfivol(mfi, volume);
        assert(grids[mfi.index()] == mfi.validbox());
        const Box &cbox = mfi.validbox();
        for(int fine = 0; fine < num_fine; fine++) {
            const Box &fbox = fgrids[fine];
            Box ovlp(coarsen(fbox,fratio));
            ovlp &= cbox;
            if(ovlp.ok()) {
                assert(fbidli);
                FillBoxId fbidFine = fbidli();
                ++fbidli;
                assert(fbidlivol);
                FillBoxId fbidFineVol = fbidlivol();
                ++fbidlivol;

                FArrayBox fine_fab(fbidFine.box(), num_comp);
                FArrayBox fine_vol(fbidFineVol.box(), 1);

                mfcd.FillFab(mfidS_fine,  fbidFine,    fine_fab);
                mfcd.FillFab(mfidFineVol, fbidFineVol, fine_vol);

                avgDown( fine_fab,  mfi(), 
                         fine_vol, mfivol(), 
                         f_level,       c_level,      
                         ovlp, strt_comp, num_comp, fratio );
            }
        }
    }
    
}  // end avgDown(...)


//------------------------------------------------------------
// average fine down to coarse in the ovlp intersection
void NavierStokes::avgDown( const FArrayBox& fine_fab, const FArrayBox& crse_fab, 
                            const FArrayBox& fine_vol, const FArrayBox& crse_vol,
                            int f_level,               int c_level,
                            const Box &ovlp,
                            int strt_comp, int num_comp, IntVect& fratio )
{
    // get the bounds
    const int* ovlo = ovlp.loVect();
    const int* ovhi = ovlp.hiVect();
    
    // get the fine data
    const int* flo            = fine_fab.loVect();
    const int* fhi            = fine_fab.hiVect();
    const Real* f_dat         = fine_fab.dataPtr(strt_comp);
    const int* fvlo           = fine_vol.loVect();
    const int* fvhi           = fine_vol.hiVect();
    const Real* fv_dat        = fine_vol.dataPtr();
    
    // get the coarse data
    const int* clo            = crse_fab.loVect();
    const int* chi            = crse_fab.hiVect();
    const Real* c_dat         = crse_fab.dataPtr(strt_comp);
    const int* cvlo           = crse_vol.loVect();
    const int* cvhi           = crse_vol.hiVect();
    const Real* cv_dat        = crse_vol.dataPtr();
    
    // call the fortran
    FORT_AVGDOWN (c_dat,ARLIM(clo),ARLIM(chi),&num_comp,
                  f_dat,ARLIM(flo),ARLIM(fhi),
                  cv_dat,ARLIM(cvlo),ARLIM(cvhi),
                  fv_dat,ARLIM(fvlo),ARLIM(fvhi),
                  ovlo,ovhi,fratio.getVect());
}


// -------------------------------------------------------------
// The Level Sync correction function
// -------------------------------------------------------------
void NavierStokes::level_sync()
{
    int i;

    // get parameters
    const Real* dx    = geom.CellSize();
    IntVect ratio     = parent->refRatio(level);
    int crse_dt_ratio = parent->MaxRefRatio(level);
    int finest_level  = parent->finestLevel();
    int ngrids        = grids.length();
    Real dt           = parent->dtLevel(level);
    Real prev_time    = state[State_Type].prevTime();
    Real half_time    = prev_time + 0.5*dt;

    // get objects
    MultiFab& pres              = get_new_data(Press_Type);
    MultiFab& vel               = get_new_data(State_Type);
    SyncRegister& rhs_sync_reg  = getLevel(level+1).getSyncReg();
    SyncRegister* crsr_sync_ptr = NULL;
    const Box& P_domain         = state[Press_Type].getDomain();
    NavierStokes& fine_level    = getLevel(level+1);
    MultiFab& pres_fine         = fine_level.get_new_data(Press_Type);

    // get rho at t^{n+1/2}

#if (USEOLDFILLPATCH == 1)
    for (int k = 0; k < ngrids; k++) {
        getState((*rho_half)[k],k,1,Density,1,half_time);
    }
#else
        int boxGrow  = 1;
        int destComp = 0;
        int srcComp  = Density;
        int nComp    = 1;
        for(FillPatchIterator rho_halffpi(*this, *rho_half, boxGrow, destComp,
                                          half_time, State_Type, srcComp, nComp
                                          /* , mapper = NULL */);
            rho_halffpi.isValid();
            ++rho_halffpi)
        {
          DependentMultiFabIterator rho_halfdest(rho_halffpi, (*rho_half));

          // copy from fillpatched fab to rho_half
          rho_halfdest().copy(rho_halffpi());
        }
#endif
    
    if (level > 0) crsr_sync_ptr = &(getLevel(level).getSyncReg());

    // get boundary conditions
    int** sync_bc =  new (int*[grids.length()]);
    for (i = 0; i < grids.length(); i++) {
      sync_bc[i] = getBCArray( State_Type,i,Xvel,BL_SPACEDIM);
    }
    
    
    // ================================ Multilevel Sync projection
    // ================================ or single level
    if (do_MLsync_proj)
    {
        
        MultiFab& vel_fine          = fine_level.get_new_data(State_Type);
        MultiFab& rho_fine          = *fine_level.rho_avg;
        const Geometry& fine_geom   = parent->Geom(level+1);
        const Geometry& crse_geom   = parent->Geom(level);
        const BoxArray&   finegrids = vel_fine.boxArray();
        const BoxArray& P_finegrids = pres_fine.boxArray();
        MultiFab    phi(P_finegrids,1,1);
        MultiFab V_corr(  finegrids,BL_SPACEDIM,1);
        V_corr.setVal(0.);
            
        // if periodic, enforce periodicity on Vsync
        projector->EnforcePeriodicity( *Vsync, BL_SPACEDIM, grids, crse_geom );
    
        // interpolate Vsync to fine grid correction in Vcorr
        SyncInterp( *Vsync, level, V_corr, level+1, ratio,
                    0, 0, BL_SPACEDIM, 0 , dt, sync_bc);
        
        // the multilevel projection
        // This computes the projection and adds in its contribution
        // to levels (level) and (level+1)
        projector->MLsyncProject(level,pres,vel,
                                 pres_fine,vel_fine,
                                 *rho_half, rho_fine, 
                                 Vsync,V_corr,phi,
                                 &rhs_sync_reg,crsr_sync_ptr,sync_bc,
                                 dx,dt,ratio,crse_dt_ratio,
                                 fine_geom,geom);
        
        // correct pressure and velocities after the projection
        ratio = IntVect::TheUnitVector();
        int ** fine_sync_bc = new (int*[finegrids.length()]);

        for (i = 0; i < finegrids.length(); i++) {
          fine_sync_bc[i] = 
            getLevel(level+1).getBCArray( State_Type,i,Xvel,BL_SPACEDIM);
        }

        for (int lev = level+2; lev <= finest_level; lev++) {
            ratio                 *= parent->refRatio(lev-1);
            NavierStokes& fine_lev = getLevel(lev);
            MultiFab &P_new        = fine_lev.get_new_data(Press_Type);
            MultiFab &U_new        = fine_lev.get_new_data(State_Type);

            SyncInterp( V_corr, level+1, U_new, lev, ratio,
                        0, 0, BL_SPACEDIM, 1 , dt, fine_sync_bc);
            SyncProjInterp( phi, level+1, P_new, lev, ratio );
        }

        for (int i = 0; i < finegrids.length(); i++) {
          delete fine_sync_bc[i];
        }
        delete[] fine_sync_bc;
        
    }
    else if ( do_sync_proj) 
    {
        
        const BoxArray& P_grids = pres.boxArray();
        BoxArray sync_boxes     = pres_fine.boxArray();
        MultiFab phi(P_grids,1,1);
        sync_boxes.coarsen(ratio);
        
        // the single level projection
        // This computes the projection and adds in its contribution
        // to level (level)
        projector->syncProject(level,pres,vel,rho_half,Vsync,phi,
                               &rhs_sync_reg,crsr_sync_ptr,sync_boxes,
                               sync_bc,geom,dx,dt,crse_dt_ratio);

        // correct pressure and velocities after the projection
        ratio = IntVect::TheUnitVector(); 
        for (int lev = level+1; lev <= finest_level; lev++) {
            ratio                 *= parent->refRatio(lev-1);
            NavierStokes& fine_lev = getLevel(lev);
            MultiFab &P_new        = fine_lev.get_new_data(Press_Type);
            MultiFab &U_new        = fine_lev.get_new_data(State_Type);
            
            SyncInterp( *Vsync, level, U_new, lev, ratio,
                        0, 0, BL_SPACEDIM, 1 , dt, sync_bc);
            SyncProjInterp( phi, level, P_new, lev, ratio );
        }
    }
    
    // garbage collection
    for ( i = 0; i < grids.length(); i++) {
        delete sync_bc[i];
    }
    delete[] sync_bc;
}


// -------------------------------------------------------------
// The Mac Sync correction function
// -------------------------------------------------------------
void NavierStokes::mac_sync()
{
    int i,lev,sigma;
    int finest_level = parent->finestLevel();
    int ngrids       = grids.length();
    int numscal      = NUM_STATE - BL_SPACEDIM;
    
    Real prev_time = state[State_Type].prevTime();
    Real prev_pres_time = state[Press_Type].prevTime();
    Real dt = parent->dtLevel(level);
    
    // compute the u_mac for the correction
    mac_projector->mac_sync_solve(level,u_mac,dt,rho_half,fine_ratio);
    
    // update coarse grid state by adding correction from mac_sync solve
    // the correction is the advective tendency of the new velocities
    if (do_reflux) {
        MultiFab& S_new = get_new_data(State_Type);
        mac_projector->mac_sync_compute(level,u_mac,Vsync,Ssync,
                                        rho_half,
                                        (level > 0) ? &getAdvFluxReg(level) : 0,
                                        is_conservative, prev_time,
                                        prev_pres_time,dt,
                                        NUM_STATE,be_cn_theta);
// the following used to be done in mac_sync_compute
        Ssync->mult(dt,Ssync->nGrow());

        // compute viscous sync
        if (is_diffusive[Xvel]) {
          diffusion->diffuse_Vsync(Vsync, dt, be_cn_theta, rho_half, 1);
        }
        for (sigma  = 0; sigma < numscal; sigma++) {
            int rho_flag = 0;
            int do_viscsyncflux = 1;
            if (!is_conservative[BL_SPACEDIM+sigma]) {
              rho_flag=1;
            } else {
              rho_flag=2;
            }
            if (is_diffusive[BL_SPACEDIM+sigma])
                diffusion->diffuse_Ssync(Ssync, sigma, dt, be_cn_theta,
                                         rho_half, rho_flag, do_viscsyncflux);
        }

        // add the sync correction to the state
        for (sigma  = 0; sigma < numscal; sigma++) {
            for(MultiFabIterator S_newmfi(S_new); S_newmfi.isValid(); ++S_newmfi) {
                DependentMultiFabIterator Ssyncmfi(S_newmfi, *Ssync);
                assert(grids[S_newmfi.index()] == S_newmfi.validbox());
                const Box& grd = S_newmfi.validbox();
                FArrayBox& s_sync = Ssyncmfi();
                S_newmfi().plus(s_sync,grd,sigma,BL_SPACEDIM+sigma,1);
            }
        }
    
        // get boundary conditions
        int** sync_bc =  new (int*[grids.length()]);
        for ( i = 0; i < ngrids; i++) {
          sync_bc[i] = getBCArray( State_Type,i,Density,numscal );
        }

        // interpolate the sync correction to the finer levels
        IntVect ratio = IntVect::TheUnitVector();
        Real mult = 1.0;
        for ( lev = level+1; lev <= finest_level; lev++) {
            ratio                 *= parent->refRatio(lev-1);
            NavierStokes& fine_lev = getLevel(lev);
            MultiFab &S_new        = fine_lev.get_new_data(State_Type);
            SyncInterp( *Ssync, level, S_new, lev, ratio,
                        0, BL_SPACEDIM, numscal, 1 , mult, sync_bc);
        }

        // garbage collection
        for ( i = 0; i < grids.length(); i++) {
            delete sync_bc[i];
        }
        delete[] sync_bc;
    }
}


// -------------------------------------------------------------
// The reflux function
// -------------------------------------------------------------
void NavierStokes::reflux()
{
    if (level == parent->finestLevel()) return;
    assert(do_reflux);

    MultiFab &S_crse = get_new_data(State_Type);

      // first do refluxing step
    FluxRegister& fr_adv  = getAdvFluxReg(level+1);
    FluxRegister& fr_visc = getViscFluxReg(level+1);

    Real dt_crse = parent->dtLevel(level);
    Real scale = 1.0/dt_crse;

    //  DONT ACTUALLY REFLUX HERE, JUST FILL VSYNC AND SSYNC
    //  IMPORTANT TO DO VISCOUS FIRST BECAUSE OF DENSITY-WEIGHTING OF V_SYNC
    fr_visc.Reflux(*Vsync,volume,scale,0,0,BL_SPACEDIM,geom);
    fr_visc.Reflux(*Ssync,volume,scale,
                   BL_SPACEDIM,0,NUM_STATE-BL_SPACEDIM,geom);

    for(MultiFabIterator rho_halfmfi(*rho_half); rho_halfmfi.isValid();
        ++rho_halfmfi)
    {
      DependentMultiFabIterator Vsyncmfi(rho_halfmfi, *Vsync);
      assert(grids[rho_halfmfi.index()] == rho_halfmfi.validbox());
      FArrayBox& Rh = rho_halfmfi();
      Vsyncmfi().divide(Rh,rho_halfmfi.validbox(),0,Xvel,1);
      Vsyncmfi().divide(Rh,rho_halfmfi.validbox(),0,Yvel,1);
#if (BL_SPACEDIM == 3)
      Vsyncmfi().divide(Rh,rho_halfmfi.validbox(),0,Zvel,1);
#endif
    }

    int istate;
    for (istate = BL_SPACEDIM; istate < NUM_STATE; istate++) {
      if (!is_conservative[istate]) {
        int sigma = istate -  BL_SPACEDIM;
        //for (i = 0; i < grids.length(); i++)
        for(MultiFabIterator rho_halfmfi(*rho_half); rho_halfmfi.isValid();
            ++rho_halfmfi)
        {
          DependentMultiFabIterator Ssyncmfi(rho_halfmfi, *Ssync);
          assert(grids[rho_halfmfi.index()] == rho_halfmfi.validbox());
          FArrayBox& Rh = rho_halfmfi();
          Ssyncmfi().divide(Rh,rho_halfmfi.validbox(),0,sigma,1);
        }
      }
    }

    fr_adv.Reflux(*Vsync,volume,scale,0,0,BL_SPACEDIM,geom);
    fr_adv.Reflux(*Ssync,volume,scale,
                  BL_SPACEDIM,0,NUM_STATE-BL_SPACEDIM,geom);

    const BoxArray& fine_boxes = getLevel(level+1).boxArray();
    int nfine = fine_boxes.length();
    int ngrids = grids.length();

    //  This is necessary in order to zero out the contribution to any
    //   coarse grid cells which underlie fine grid cells
    for (int kf = 0; kf < nfine; kf++) {
        Box bf(fine_boxes[kf]);
        bf.coarsen(fine_ratio);
        //for (int k = 0; k < ngrids; k++)
        for(MultiFabIterator Vsyncmfi(*Vsync); Vsyncmfi.isValid();
            ++Vsyncmfi)
        {
            DependentMultiFabIterator Ssyncmfi(Vsyncmfi, *Ssync);
            assert(grids[Vsyncmfi.index()] == Vsyncmfi.validbox());
            Box bx(Vsyncmfi.validbox());
            bx &= bf;
            if (bx.ok()) {
                FArrayBox& v_sync = Vsyncmfi();
                v_sync.setVal(0.0,bx,0,BL_SPACEDIM);
                FArrayBox& s_sync = Ssyncmfi();
                s_sync.setVal(0.0,bx,0,NUM_STATE-BL_SPACEDIM);
            }
        }
    }
}


// average down a single state component
void NavierStokes::avgDown(int comp)
{
    if (level == parent->finestLevel()) return;

    NavierStokes& fine_lev = getLevel(level+1);
    const BoxArray &fgrids = fine_lev.grids;
    MultiFab &fvolume      = fine_lev.volume;

    MultiFab &S_crse = get_new_data(State_Type);
    MultiFab &S_fine = fine_lev.get_new_data(State_Type);

    avgDown( grids, fgrids,
             S_crse, S_fine,
             volume, fvolume,
             level,  level+1,
             comp, 1, fine_ratio);
}


// inject fine pressure nodes down onto coarse nodes
void NavierStokes::injectDown( const Box &ovlp, FArrayBox &Pcrse,
                               const FArrayBox &Pfine, IntVect &fine_ratio )
{
    // get the coarse intersection bounds
    const int* ovlo = ovlp.loVect();
    const int* ovhi = ovlp.hiVect();

    // get the fine and coarse pressures
    Real* cpres       = Pcrse.dataPtr();
    const int* clo    = Pcrse.loVect();
    const int* chi    = Pcrse.hiVect();
    const Real* fpres = Pfine.dataPtr();
    const int* flo    = Pfine.loVect();
    const int* fhi    = Pfine.hiVect();

    // call the fortran
    FORT_PUTDOWN (cpres,ARLIM(clo),ARLIM(chi),
                  fpres,ARLIM(flo),ARLIM(fhi),
                  ovlo,ovhi,fine_ratio.getVect());
}



// test for consistency between fine and coarse nodes
void NavierStokes::testInject( const Box &ovlp, FArrayBox &Pcrse,
                               const FArrayBox &Pfine, IntVect &fine_ratio )
{
    // get the coarse intersection bounds
    const int* ovlo = ovlp.loVect();
    const int* ovhi = ovlp.hiVect();

    // get the fine and coarse pressures
    Real* cpres       = Pcrse.dataPtr();
    const int* clo    = Pcrse.loVect();
    const int* chi    = Pcrse.hiVect();
    const Real* fpres = Pfine.dataPtr();
    const int* flo    = Pfine.loVect();
    const int* fhi    = Pfine.hiVect();

    // call the fortran
    FORT_TESTINJECT( cpres,ARLIM(clo),ARLIM(chi),
                     fpres,ARLIM(flo),ARLIM(fhi),
                     ovlo,ovhi,fine_ratio.getVect());
}


//-------------------------------------------------------------
// average the fine information from the complete set of state types
// down to coarse
//-------------------------------------------------------------
void NavierStokes::avgDown()
{
    if (level == parent->finestLevel()) return;

    NavierStokes& fine_lev = getLevel(level+1);
    const BoxArray &fgrids = fine_lev.grids;
    MultiFab &fvolume      = fine_lev.volume;

    // ********************************************************** 
    // average down the sates at the new time
    // ********************************************************** 

    MultiFab &S_crse = get_new_data(State_Type);
    MultiFab &S_fine = fine_lev.get_new_data(State_Type);

    avgDown( grids, fgrids,
             S_crse, S_fine,
             volume, fvolume,
             level,  level+1,
             0, S_crse.nComp(), fine_ratio );

    // ********************************************************** 
    // Now average down pressure over time n-(n+1) interval      
    // ********************************************************** 

    // get objects
    MultiFab &P_crse      = get_new_data(Press_Type);
    MultiFab &P_fine_init = fine_lev.get_new_data(Press_Type);
    MultiFab &P_fine_avg  = *fine_lev.p_avg;
    MultiFab &P_fine      = ( initial_step ? P_fine_init : P_fine_avg );

    const BoxArray& P_cgrids = state[Press_Type].boxArray();
    const BoxArray& P_fgrids = fine_lev.state[Press_Type].boxArray();
    const Geometry& fgeom    = parent->Geom(level+1);
    Box domain               = fgeom.Domain();
    domain.surroundingNodes();

    // get looping parameters
    Box ovlp;
    Array<IntVect> pshifts(27);
    int ng_pres = P_crse.nGrow();
    int ncrse   = grids.length();
    int nfine   = fgrids.length();

/*
    int crse, fine;
    // inject fine pressure nodes down onto coarse nodes
    for (crse = 0; crse < ncrse; crse++) {
        // the coarse grid nodes
        const Box& cbox = P_cgrids[crse];
        
        // loop over fine grids and periodic extensions
        for (fine = 0; fine < nfine; fine++) {

            // inject fine down to coarse
            ovlp  = coarsen(P_fgrids[fine],fine_ratio);
            ovlp &= cbox;
            if (ovlp.ok()) {
                injectDown( ovlp, P_crse[crse], P_fine[fine], fine_ratio );
            }

            // // test for proper injecting of coarse to fine
            // if ( fgeom.isAnyPeriodic() ) {
            //     fgeom.periodicShift(domain, P_fgrids[fine], pshifts );
            //     for( int iiv=0; iiv<pshifts.length(); iiv++) {
            //         IntVect iv=pshifts[iiv];
            //         P_fine[fine].shift(iv);
            //         ovlp  = P_fgrids[fine];
            //         ovlp.shift(iv);
            //         ovlp.coarsen(fine_ratio);
            //         ovlp &= cbox;
            //         if (ovlp.ok()) {
            //           testInject( ovlp, P_crse[crse], P_fine[fine], fine_ratio );
            //         }
            //         P_fine[fine].shift(-iv);
            //     }
            // }
            
        } // end of fine grid loop
    } // end of coarse grid loop
*/

//
//
//  this needs to be optimized to communicate only the
//  overlap part of the fab
//
//


    MultiFabCopyDescriptor mfcd(true);
    MultiFabId mfidP_fine  = mfcd.RegisterFabArray(&P_fine);

    List<FillBoxId> fillBoxIdList;

    // inject fine pressure nodes down onto coarse nodes
    for(MultiFabIterator mfi(P_crse); mfi.isValid(); ++mfi) {
        assert(P_cgrids[mfi.index()] == mfi.validbox());
        const Box &cbox = mfi.validbox();
        
        // loop over fine grids and periodic extensions
        for(int fine = 0; fine < nfine; fine++) {
            ovlp  = coarsen(P_fgrids[fine],fine_ratio);
            ovlp &= cbox;
            if(ovlp.ok()) {         // inject fine down to coarse
              //injectDown( ovlp, P_crse[crse], P_fine[fine], fine_ratio );

              BoxList unfilledBoxes(ovlp.ixType());  // unused here
              FillBoxId fbidFine;
              const Box &sfineBox = P_fine.fabbox(fine);
              fbidFine = mfcd.AddBox(mfidP_fine, sfineBox, unfilledBoxes,
                                     0, 0, P_fine.nComp());
              fillBoxIdList.append(fbidFine);
            }
        } // end of fine grid loop
    } // end of coarse grid loop


    mfcd.CollectData();

    ListIterator<FillBoxId> fbidli(fillBoxIdList);


    for(MultiFabIterator mfi(P_crse); mfi.isValid(); ++mfi) {
        assert(P_cgrids[mfi.index()] == mfi.validbox());
        const Box &cbox = mfi.validbox();
        
        // loop over fine grids and periodic extensions
        for(int fine = 0; fine < nfine; fine++) {
            ovlp  = coarsen(P_fgrids[fine],fine_ratio);
            ovlp &= cbox;
            if(ovlp.ok()) {         // inject fine down to coarse
                assert(fbidli);
                FillBoxId fbidFine = fbidli();
                ++fbidli;

                FArrayBox fine_fab(fbidFine.box(), P_fine.nComp());
                mfcd.FillFab(mfidP_fine,  fbidFine,    fine_fab);

                injectDown( ovlp, mfi(), fine_fab, fine_ratio );
            }
        } // end of fine grid loop
    } // end of coarse grid loop



    // ******************************************************** 
    // Next average down divu and dSdT at new time                      
    // ******************************************************** 
    if(have_divu) {
        MultiFab &Divu_crse = get_new_data(Divu_Type);
        MultiFab &Divu_fine = fine_lev.get_new_data(Divu_Type);
        
        avgDown( grids,     fgrids,
                 Divu_crse, Divu_fine,
                 volume,    fvolume,
                 level,     level+1,
                 0, 1, fine_ratio );
    }
    if(have_dsdt) {
        MultiFab &Dsdt_crse = get_new_data(Dsdt_Type);
        MultiFab &Dsdt_fine = fine_lev.get_new_data(Dsdt_Type);
        
        avgDown( grids,     fgrids,
                 Dsdt_crse, Dsdt_fine,
                 volume,    fvolume,
                 level,     level+1,
                 0, 1, fine_ratio );
    }
}  // end avgDown()



//=================================================================
// ACCESS FUNCTIONS FOLLOW
//=================================================================

// ---------------------------------------------------------------
// virtual access function for getting the advective flux out of the
// advection routines for diagnostics and refluxing
// ---------------------------------------------------------------
void NavierStokes::pullFluxes( int level, int i, int start_ind, int ncomp,
                               FArrayBox &xflux,
                               FArrayBox &yflux,
                               FArrayBox &zflux, Real dt )
{
    int finest_level = parent->finestLevel();
    
    // add fluxes into the refluxing counters
    if (do_reflux) {
        if (level < finest_level) {
            FluxRegister& fr = getAdvFluxReg(level+1);
            fr.CrseInit(xflux,xflux.box(),0,0,start_ind,ncomp,-dt);
            fr.CrseInit(yflux,yflux.box(),1,0,start_ind,ncomp,-dt);
#if (BL_SPACEDIM == 3)                              
            fr.CrseInit(zflux,zflux.box(),2,0,start_ind,ncomp,-dt);
#endif
        }
        if (level > 0) {
            advflux_reg->FineAdd(xflux,0,i,0,start_ind,ncomp,dt);
            advflux_reg->FineAdd(yflux,1,i,0,start_ind,ncomp,dt);
#if (BL_SPACEDIM == 3)                                
            advflux_reg->FineAdd(zflux,2,i,0,start_ind,ncomp,dt);
#endif
        }
    }
}


// ---------------------------------------------------------------
// virtual access function for getting the forcing terms for the
// velocities and scalars.  The base version computes a buoyancy
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
// and velocity_advection routines
// ---------------------------------------------------------------
void NavierStokes::getForce( FArrayBox& force, int gridno, int ngrow,
                             int strt_comp, int num_comp, Real time)
{
//bool canUse_getForce = false;
//assert(canUse_getForce);
cerr << "\nError:  should not be in NavierStokes::getForce\n\n";

    Box bx(grids[gridno]);
    bx.grow(ngrow);
    force.resize(bx,num_comp);

    const int* lo = bx.loVect();
    const int* hi = bx.hiVect();

    Real grav = Abs(gravity);
    for (int dc = 0; dc < num_comp; dc++) {
        int sc = strt_comp + dc;
#if (BL_SPACEDIM == 2)
        if (BL_SPACEDIM == 2 && sc == Yvel && grav > 0.001) 
#endif
#if (BL_SPACEDIM == 3)
        if (BL_SPACEDIM == 3 && sc == Zvel && grav > 0.001) 
#endif
        {
              // set force to -rho*g
            FArrayBox rho;
            getState(rho,gridno,ngrow,Density,1,time);
            rho.mult(-grav);
            force.copy(rho,0,dc,1);
        } else {
            force.setVal(0.0,dc);
        }
    }
}



// fill patch and then call the other GradP function    
void NavierStokes::getGradP(FArrayBox& gp, int gridno, int ngrow, Real time)
{
//bool canUse_getGradP1 = false;
//assert(canUse_getGradP1);
cerr << "\nError:  should not be in NavierStokes::getGradP\n\n";

    Box gpbx(grids[gridno]);
    gpbx.grow(ngrow);
    Box p_box(surroundingNodes(gpbx));
    FArrayBox p_fab(p_box,1);
    FillPatch(p_fab,0,time,Press_Type,0,1);

    getGradP( p_fab, gp, grids[gridno], ngrow );
}




// given a nodal pressure P compute the pressure gradient at the
// contained cell centers
//
// This function should live in Projection, but it made life easier
// to have it here
//
void NavierStokes::getGradP(FArrayBox &p_fab, FArrayBox& gp,
                            const Box &grd, int ngrow )
{
//bool canUse_getGradP2 = false;
//assert(canUse_getGradP2);
cerr << "\nError:  should not be in NavierStokes::getGradP (2)\n\n";

    //-----------------------  size the pressure gradient storage
    Box gpbx(grd);
    gpbx.grow(ngrow);
    gp.resize(gpbx,BL_SPACEDIM);
    
    //------------------------ test to see if p_fab contains gpbx
    Box test = p_fab.box();
    test = test.enclosedCells();
    assert( test.contains( gpbx ) == true );
    
    // ----------------------- set pointers
    const int *plo = p_fab.loVect();
    const int *phi = p_fab.hiVect();
    const int *glo = gpbx.loVect();
    const int *ghi = gpbx.hiVect();
    const Real *p_dat  = p_fab.dataPtr();
    const Real *gp_dat = gp.dataPtr();
    const Real *dx     = geom.CellSize();

    // ----------------------- create the pressure gradient
    FORT_GRADP (p_dat,ARLIM(plo),ARLIM(phi),
                gp_dat,ARLIM(glo),ARLIM(ghi),glo,ghi,dx);
}



// fill patch divU
void NavierStokes::getDivCond(FArrayBox& fab, int gridno, int ngrow, 
                         Real time)
{
//bool canUse_getDivCond = false;
//assert(canUse_getDivCond);
cerr << "\nError:  should not be in NavierStokes::getDivCond (1)\n\n";

    getState(fab, gridno, ngrow, time, have_divu, Divu_Type,
             //divu_assoc, divu_unfilled);
             divu_unfilled);
}



// fill patch dSdt
void NavierStokes::getDsdt(FArrayBox& fab, int gridno, int ngrow, 
                         Real time)
{
//bool canUse_getDsDt = false;
//assert(canUse_getDsDt);
cerr << "\nError:  should not be in NavierStokes::getDsdt\n\n";

    getState(fab, gridno, ngrow, time, (have_dsdt && have_divu), Dsdt_Type,
             //dsdt_assoc, dsdt_unfilled);
             dsdt_unfilled);
}



// fill patch a state component
void NavierStokes::getState(FArrayBox& fab, int gridno, int ngrow,
                       int strt_comp, int num_comp, Real time)
{
//bool canUse_getState1 = false;
//assert(canUse_getState1);
    BoxLib::Error("Should not be in NavierStokes::getState (1)");

    getState(fab,gridno,ngrow,State_Type,strt_comp,num_comp,time);
}

// fill patch a state component
void NavierStokes::getState(FArrayBox& fab, int gridno, int ngrow,
                        int state_indx, int strt_comp, int num_comp, 
                        Real time)
{
//bool canUse_getState2 = false;
//assert(canUse_getState2);
    BoxLib::Error("Should not be in NavierStokes::getState (2)");

#if 0
    Box bx(grids[gridno]);
    bx.grow(ngrow);
    fab.resize(bx,num_comp);

    if (ngrow == 1 && !Geometry::isAnyPeriodic()) {
        //hyp_assoc.setCacheWidth(ngrow);
        FillPatch(fab,0,time,state_indx,strt_comp,num_comp,
                 //hyp_assoc,gridno,cc1_unfilled[gridno]);
                 cc1_unfilled[gridno]);

    } else if (ngrow == HYP_GROW && !Geometry::isAnyPeriodic()) {
        //hyp_assoc.setCacheWidth(ngrow);
        FillPatch(fab,0,time,state_indx,strt_comp,num_comp,
                 //hyp_assoc,gridno,hyp_unfilled[gridno]);
                 hyp_unfilled[gridno]);
    } else {
        FillPatch(fab,0,time,state_indx,strt_comp,num_comp);
    }
#endif
}

// fills ghost cells of state:
// For finer levels, the ghost cells that are interior to the
// problem domain, but exterior to the valid region of the state
// at that that level are not properly filled by either a multifab
// FillBoundary or a setPhysBoundaryValues. This routine takes care of those 
// cells. It is not particularly efficient, but probably not too
// inefficient either--rbp
void NavierStokes::FillStateBndry(Real time, int state_indx, int src_comp, 
                            int num_comp) 
{
  MultiFab& S = get_data(state_indx,time);
  int ngrow = S.nGrow();
  if (ngrow==0) {
    return;
  }
  MultiFab component(grids,1,ngrow,Fab_allocate);
  for (int istate = src_comp; istate < src_comp+num_comp; istate++) {
#if (USEOLDFILLPATCH == 1)
    for (int i = 0; i < grids.length(); i++) {
      Box grd = grids[i];
      grd.grow(ngrow);
      int need_ghost_cells_filled = !grids.contains(grd);
      if(need_ghost_cells_filled) {
        getState(component[i],i,ngrow,state_indx,istate,1,time);
        S[i].copy(component[i],0,istate,1);
      }
    }
#else
    int srcComp  = istate;
    int destComp = 0;
    int nComp    = 1;
    for(FillPatchIterator fpi(*this, component, ngrow, destComp, time, state_indx, 
                              srcComp, nComp);
        fpi.isValid();
        ++fpi)
    {
      DependentMultiFabIterator Smfi(fpi, S);
      assert(grids[fpi.index()] == fpi.UngrownBox());
      Box grd = fpi.UngrownBox();
      grd.grow(ngrow);
      int need_ghost_cells_filled = ! grids.contains(grd);
      if(need_ghost_cells_filled) {
        Smfi().copy(fpi(),0,istate,1);
      }
    }
#endif
  }  // end for(istate...)
  S.FillBoundary();
  setPhysBoundaryValues(state_indx,src_comp,num_comp,time);
}

//-----------------------------------------------------
// fill patch a state component.  This is a driver for calling
// the amrLevel filpatch function.  amrLevel filpatch works
// by first interpolating coarse data (time, then space), overwriting
// coarse data with fine data, and letting physical boundary conditions
// take precedence over interpolated ghost cells
//
// NOTE:: It might be a bug to resize the fab, since it could
// be part of a multifab with unforeseen consequences
//-----------------------------------------------------
void NavierStokes::getState(FArrayBox& fab, int gridno, int ngrow, 
                       Real time, int have_state, int state_indx,
                       //BoxAssoc& assoc, Array<Box>& unfilled, int
                       Array<Box>& unfilled, int
                       num_comp, int strt_comp)
{
//bool canUse_getState3 = false;
//assert(canUse_getState3);
    BoxLib::Error("Should not be in NavierStokes::getState (3)");

#if 0
    // create the storage
    Box bx(grids[gridno]);
    bx.grow(ngrow);
    fab.resize(bx,num_comp);
    if(fab.box()!=bx) {
      cout << "NavierStokes::getState : fab.box()!=bx\n";
      cout << "bx        = " << bx << NL;
      cout << "fab.box() = " << fab.box() << NL;
      ParallelDescriptor::Abort("NavierStokes::getState(...)");
    }
    
    if(!have_state) {
        
        fab.setVal(0.0);  // can't filpatch what we don't have
        
    } else {
        
        fab.setVal(1.0e30); // for debugging only
        
        if (ngrow == 1 && !Geometry::isAnyPeriodic()) {
            //assoc.setCacheWidth(ngrow);
            FillPatch(fab,0,time,state_indx,strt_comp,num_comp,
                     //assoc,gridno,unfilled[gridno]);
                     unfilled[gridno]);
        } else {
            FillPatch(fab,0,time,state_indx,strt_comp,num_comp);
        }
    }
#endif
}


// fill patch divU
void NavierStokes::getDivCond(FArrayBox& fab, int ngrow, Real time)
{
//bool canUse_getDivCond = false;
//assert(canUse_getDivCond);
cerr << "\nError:  should not be in NavierStokes::getDivCond (2)\n\n";

    // for NavierStokes, defaults divU = 0
    if(!have_divu) {
        fab.setVal(0.0);
    } else {
        fab.setVal(1.0e30); // for debugging only
        FillPatch(fab,0,time,Divu_Type,0,1);
    }
}

// default divU is set to zero
void NavierStokes::calc_divu(Real time, Real dt, MultiFab& divu)
{
    if (have_divu) {
      divu.setVal(0.0);
    }
}

// default dSdt is set to zero
void NavierStokes::calc_dsdt(Real time, Real dt, MultiFab& dsdt)
{
    if (have_divu && have_dsdt) {
      dsdt.setVal(0.0);
    }
}

// -------------------------------------------------------------
void
NavierStokes::getViscTerms(MultiFab& visc_terms, int src_comp, 
                           int num_comp, Real time)
{
    int rho_flag;
    for (int icomp = src_comp; icomp < src_comp+num_comp; icomp++) { 
      if(!is_conservative[icomp]) {
        rho_flag=1;
      } else {
        rho_flag=2;
      }
      diffusion->getViscTerms(visc_terms,src_comp,icomp,time,rho_flag);
    }
}

// -------------------------------------------------------------
void 
NavierStokes::compute_grad_divu_minus_s(Real time, MultiFab* grad_divu_minus_s,
                                        int scaleRhoDivDt)
{
    Real dt = parent->dtLevel(0);
    grad_divu_minus_s->setVal(0.0);
#if (BL_SPACEDIM == 2)
    MultiFab &S_new = get_new_data(State_Type);
    S_new.FillBoundary();
    setPhysBoundaryValues(State_Type,Xvel,BL_SPACEDIM,time);

    BoxArray node_grids = grids;
    node_grids.surroundingNodes();
    MultiFab S            (node_grids,1,0,Fab_allocate);
    MultiFab Dv           (node_grids,1,0,Fab_allocate);

    int is_rz = (CoordSys::IsRZ() ? 1 : 0);

    int nghost = 0;
    if(have_divu) {

      projector->put_divu_in_node_rhs(S, parent, level,
                  nghost,
                  node_grids, time);
    } else {
      S.setVal(0.0);
    }

    int nGrow = 1;
    MultiFab U (grids, BL_SPACEDIM, nGrow, Fab_allocate);
#if (USEOLDFILLPATCH == 1)
    FArrayBox Ufab;
    int ngrids = grids.length();
    for (int igrid = 0; igrid < ngrids; igrid++) {
      getState(Ufab, igrid, nGrow, 0, BL_SPACEDIM, time);
      U[igrid].copy(Ufab);
    }
#else
    int srcComp  = 0;
    int destComp = 0;
    int nComp    = BL_SPACEDIM;
    for(FillPatchIterator Ufpi(*this, U, nGrow,
                               destComp, time, State_Type,
                               srcComp, nComp);
        Ufpi.isValid(); ++Ufpi)
    {
      DependentMultiFabIterator Umfi(Ufpi, U);
      Umfi().copy(Ufpi());
    }
#endif

    if (is_rz) {
      for (int n = 0; n < BL_SPACEDIM; n++) 
        projector->radMult(level,U,n);
    }

    const Real* dx = geom.CellSize();
    projector->computeDV(Dv, U, Xvel, dx, is_rz);

    if (is_rz) {
      for (int n = 0; n < BL_SPACEDIM; n++) 
        projector->radDiv(level,U,n);
    }

// fix dv (mult by 2) at walls.
// Outflow, Inflow, and Symmetry should be ok already.

    const Box& domain = geom.Domain();
    Box node_domain = surroundingNodes(domain);
    const int* phys_lo = phys_bc.lo();
    const int* phys_hi = phys_bc.hi();
    const int* dlo = node_domain.loVect();
    const int* dhi = node_domain.hiVect();

    int dir;
    for (dir = 0; dir < BL_SPACEDIM; dir++) {
      //for (i = 0; i < node_grids.length(); i++)
      for(MultiFabIterator Dvmfi(Dv); Dvmfi.isValid(); ++Dvmfi) {
        Box domlo(node_domain), domhi(node_domain);
        domlo.setRange(dir,dlo[dir],1);
        domhi.setRange(dir,dhi[dir],1);
        Box blo(node_grids[Dvmfi.index()]);
        Box bhi(blo);
        blo &= domlo;
        bhi &= domhi;
        if(phys_lo[dir]==SlipWall || phys_lo[dir]==NoSlipWall) 
          if (blo.ok()) Dvmfi().mult(2.0,blo,0,1);
        if(phys_hi[dir]==SlipWall || phys_hi[dir]==NoSlipWall) 
          if (bhi.ok()) Dvmfi().mult(2.0,bhi,0,1);
      }
    }

// end, dv fix

#if 0
    for (i = 0; i < grids.length(); i++) {

        // compute nodal divU-S 
        FArrayBox divu_minus_s(node_grids[i],1);
        divu_minus_s.copy(Dv[i]);
        divu_minus_s.minus(S[i]);

        // compute gradient of this field
        getGradP( divu_minus_s, (*grad_divu_minus_s)[i], grids[i], 0 );
      
        if (scaleRhoDivDt) {
            (*grad_divu_minus_s)[i].mult(divu_relax_factor*dx[0]*dx[1]/parent->dtLevel(0));
            int n;
            for (n=0; n<BL_SPACEDIM; n++)
                (*grad_divu_minus_s)[i].mult((*rho_half)[i],grids[i],0,n,1);
        } else {
            (*grad_divu_minus_s)[i].mult(divu_relax_factor*dx[0]*dx[1]*dt/parent->dtLevel(0));
        }
    }

    if (is_rz) {
        for (int n = 0; n < BL_SPACEDIM; n++) 
            projector->radDiv(level,*grad_divu_minus_s,n);
    }
#endif
#if 1
    MultiFab divu_minus_s(node_grids,1,0,Fab_allocate);
    //for (i = 0; i < grids.length(); i++)
    for(MultiFabIterator divu_minus_smfi(divu_minus_s);
        divu_minus_smfi.isValid(); ++divu_minus_smfi)
    {
        DependentMultiFabIterator Dvmfi(divu_minus_smfi, Dv);
        DependentMultiFabIterator Smfi(divu_minus_smfi, S);

        // compute nodal divU-S 
        divu_minus_smfi().copy(Dvmfi());
        divu_minus_smfi().minus(Smfi());
    }

#if 1
    if(is_rz) {
      Box dom_lo(node_domain);
      dom_lo.setRange(0,dlo[0],1);
      //for (i = 0; i < grids.length(); i++)
      for(MultiFabIterator divu_minus_smfi(divu_minus_s);
          divu_minus_smfi.isValid(); ++divu_minus_smfi)
      {
        assert(grids[divu_minus_smfi.index()] == divu_minus_smfi.validbox());
        Box grid = divu_minus_smfi.validbox();
        grid.surroundingNodes();
        if(grid.intersects(dom_lo)) {
          Box left_strip = (grid & dom_lo);
          divu_minus_smfi().setVal(0.0,grid,0,1);
        }
      }
    }
#endif

    //for (i = 0; i < grids.length(); i++)
    for(MultiFabIterator divu_minus_smfi(divu_minus_s);
        divu_minus_smfi.isValid(); ++divu_minus_smfi)
    {
        DependentMultiFabIterator grad_divu_minus_smfi(divu_minus_smfi,
                                                       (*grad_divu_minus_s));
        assert(grids[divu_minus_smfi.index()] == divu_minus_smfi.validbox());
        // compute gradient of this field
        getGradP(divu_minus_smfi(), grad_divu_minus_smfi(), 
                 divu_minus_smfi.validbox(), 0 );
    }
#if 1
    if (is_rz) {
// what we just computed was grad (r * (dU - S))
// we now compute grad (dU - S) by 
//  grad (dU - S) = 1/r grad (r * (dU - S)) - 1/r grad (r) * (dU-S) (***)
        for (int n = 0; n < BL_SPACEDIM; n++) 
            projector->radDiv(level,*grad_divu_minus_s,n);

        MultiFab divu_minus_s_cc(grids,1,0,Fab_allocate);
        //for (i = 0; i < grids.length(); i++)
        for(MultiFabIterator divu_minus_smfi(divu_minus_s);
            divu_minus_smfi.isValid(); ++divu_minus_smfi)
        {
          DependentMultiFabIterator divu_minus_s_ccmfi(divu_minus_smfi,
                                                       divu_minus_s_cc);
          DEF_CLIMITS(divu_minus_smfi(),nodedat,nodelo,nodehi);
          DEF_LIMITS(divu_minus_s_ccmfi(),ccdat,cclo,cchi);
          int rweighted = 1;
          FORT_HGN2C(&is_rz,&rweighted, ARLIM(nodelo), ARLIM(nodehi), nodedat,
                     ARLIM(cclo), ARLIM(cchi), cclo, cchi, ccdat);
        }
// the first divide is because divu_minus_s_cc contains r*(dU-S)
// the second divide account for the 1/r in expression (***) above
        projector->radDiv(level,divu_minus_s_cc,0);
        projector->radDiv(level,divu_minus_s_cc,0);
// since grad r = (1,0), we subtract  1/r * (dU-S from the first component
        //for (i = 0; i < grids.length(); i++) 
        for(MultiFabIterator divu_minus_s_ccmfi(divu_minus_s_cc);
            divu_minus_s_ccmfi.isValid(); ++divu_minus_s_ccmfi)
        {
          DependentMultiFabIterator grad_divu_minus_smfi(divu_minus_s_ccmfi,
                                                         (*grad_divu_minus_s));
          grad_divu_minus_smfi().minus(divu_minus_s_ccmfi(),0,0,1);
        }
    }
#endif

    //for (i = 0; i < grids.length(); i++)
    for(MultiFabIterator grad_divu_minus_smfi(*grad_divu_minus_s);
        grad_divu_minus_smfi.isValid(); ++grad_divu_minus_smfi)
    {
        DependentMultiFabIterator rho_halfmfi(grad_divu_minus_smfi, (*rho_half));
        assert(grids[grad_divu_minus_smfi.index()] == grad_divu_minus_smfi.validbox());
        if (scaleRhoDivDt) {
            grad_divu_minus_smfi().mult(divu_relax_factor*dx[0]*dx[1]/parent->dtLevel(0));
            int n;
            for (n=0; n<BL_SPACEDIM; n++)
                grad_divu_minus_smfi().mult(rho_halfmfi(),
                                            grad_divu_minus_smfi.validbox(),0,n,1);
        } else {
            grad_divu_minus_smfi().mult(divu_relax_factor*dx[0]*dx[1]*dt/parent->dtLevel(0));
        }
    }

#endif

#endif
}

