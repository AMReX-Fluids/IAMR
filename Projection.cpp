
//
// $Id: Projection.cpp,v 1.39 1998-05-29 17:32:59 lijewski Exp $
//

#ifdef BL_T3E
#include <List.H>
#endif

#include <Misc.H>
#include <CoordSys.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <BoxDomain.H>
#include <MacOpProjDrivers.H>
#include <Projection.H>
#include <RunStats.H>
#include <SYNCREG_F.H>
#include <PROJECTION_F.H>
#include <NAVIERSTOKES_F.H>
#include <hg_projector.H>

const char NL = '\n';

#ifndef _NavierStokes_H_
enum StateType {State_Type=0, Press_Type};
#if (BL_SPACEDIM == 2)
enum StateNames  { Xvel=0, Yvel, Density};
#else
enum StateNames  { Xvel=0, Yvel, Zvel, Density};
#endif
#endif

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();
#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

#define SET_BOGUS_BNDRY 1
#define bogus_value 1.e20

#define MAX_LEV 10

// initialization of static members
int       Projection::verbose = 0;
int       Projection::P_code = 0;
Real      Projection::proj_tol = 1.0e-12;
Real      Projection::sync_tol = 1.0e-8;
Real      Projection::proj_abs_tol = 1.0e-16;
Real      Projection::filter_factor = 0.0;
int       Projection::filter_u = 0;
int       Projection::rho_wgt_vel_proj = 0;
int       Projection::do_outflow_bcs = 1;
int       Projection::make_sync_solvable = 0;

static RegType project_bc[] =
{interior, inflow, outflow, refWall, refWall, refWall};

// ==================================================
// Setup functions follow
// ==================================================

// the constructor
Projection::Projection(Amr * _parent,
                       BCRec * _phys_bc, 
                       int _do_sync_proj, int _finest_level, 
                       int _radius_grow )
   : parent(_parent), phys_bc(_phys_bc), 
     do_sync_proj(_do_sync_proj), finest_level(_finest_level),
     radius_grow(_radius_grow), 
     LevelData(_finest_level+1), 
     radius(_finest_level+1)
{

#ifdef BL_T3E
// force instantiation of List<int>::clear() for the T3E
List<int> tempIntList;
tempIntList.clear();
#endif

  read_params();

  if (verbose) 
  {
    cout << "Creating projector\n";
  }

#if BL_SPACEDIM == 2
  if ( CoordSys::IsRZ() ) amr_multigrid::SetRZ();
#endif
  setUpBcs();
  sync_proj = 0;
}




// the destructor
Projection::~Projection()
{
  if (verbose) 
  {
    cout << "Deleting projector\n";
  }
  delete sync_proj;
  delete projector_bndry;
}




// read parameters from input file and command line
void Projection::read_params()
{
  // read parameters from input file and command line
  ParmParse pp("proj");

  pp.query("v",verbose);
  pp.query("Pcode",P_code);

  pp.query("proj_tol",proj_tol);
  pp.query("sync_tol",sync_tol);
  pp.query("proj_abs_tol",proj_abs_tol);

  pp.query("make_sync_solvable",make_sync_solvable);

  pp.query("filter_factor",filter_factor);
  pp.query("filter_u",filter_u);


  pp.query("rho_wgt_vel_proj",rho_wgt_vel_proj);

  pp.query("do_outflow_bcs",do_outflow_bcs);
}




// set up bcs
void 
Projection::setUpBcs()
{
  const Geometry& geom = parent->Geom(0);

  // set up projector bndry
  const int* lo_bc = phys_bc->lo();
  const int* hi_bc = phys_bc->hi();
  RegType proj_bc[BL_SPACEDIM][2];
  for(int i = 0; i < BL_SPACEDIM; ++i)
  {
      proj_bc[i][0] = project_bc[lo_bc[i]];
      proj_bc[i][1] = project_bc[hi_bc[i]];
      if ( geom.isPeriodic(i) ) 
      {
          proj_bc[i][0] = periodic;
          proj_bc[i][1] = periodic;
      }
  }

  projector_bndry = new inviscid_fluid_boundary_class(proj_bc);
}

// install a level of the projection
void
Projection::install_level(int level, AmrLevel * level_data,
                          PArray<Real> * _radius )
{
  if (verbose) 
  {
    cout << "Installing projector level " << level << NL;
  }

  finest_level = parent->finestLevel();

  if (level > LevelData.length() - 1) 
  {
    LevelData.resize(finest_level+1);
    radius.resize(finest_level+1);
  }

  LevelData.clear(level);
  LevelData.set(level, level_data);
  radius.clear(level);
  radius.set(level, _radius);

  if (sync_proj != 0) 
  {
    delete sync_proj;
    sync_proj = 0;
  }
}

// Build the aliasLib projection object
void
Projection::bldSyncProject()
{
  const Box& fdomain = parent->Geom(finest_level).Domain();

  // destruct if it already exists
  if (sync_proj != 0) delete sync_proj;

  const Array<IntVect>& ref_ratio = parent->refRatio();

  // build mesh and ratio arrays for entire hierarchy
  Array<BoxArray> amesh(finest_level+1);
  Array<IntVect> gen_ratio(finest_level);
  for (int lev = 0; lev <= finest_level; lev++) 
  {
    amesh.set(lev, parent->boxArray(lev));
    if (lev > 0)
      gen_ratio.set(lev-1, ref_ratio[lev-1]);
  }

  if (verbose) 
  {
    cout << "bldSyncProject:: amr_mesh = \n";
    amr_multigrid::mesh_write(amesh, gen_ratio, fdomain, cout);
  }

  sync_proj = new holy_grail_amr_projector(amesh, gen_ratio, fdomain,
                                           0, finest_level, finest_level,
                                           *projector_bndry, false, true, false, P_code); 

#ifdef ATMOSPHERE
  // This is not the usual way of setting parameters
  sync_proj->line_solve_dim = BL_SPACEDIM - 1;  
#endif
  
  if (make_sync_solvable) 
  {
    sync_proj->make_it_so();
  }
}


// ==================================================
// Projection functions follow
// ==================================================



// ----------------------------------------------------
//  perform a level projection in the advance function
//  Explanation of arguments to the level projector:
//
//  rho_half  contains rho^{n+1/2}
//  U_old  contains the u^n velocities
//  U_new  starts as u^*, is converted to (u^* - u^n)/dt,
//         becomes (u^{n+1} - u^n)/dt in the solver,
//         and is converted back to u^n+1 at the end
//  P_old  contains p^{n-1/2}
//  P_new  gets cleared, initialized to an intial guess for p^{n+1/2}
//         using coarse grid data if available,
//         becomes pressure update phi in the solver,
//         and then converted into final prssure p^{n+1/2}
// ----------------------------------------------------

void
Projection::level_project(int level,
                          Real dt, Real cur_pres_time,
                          Real old_state_time, Real cur_state_time,
                          const Geometry& geom, 
                          MultiFab &U_old,
                          MultiFab &U_new,
                          MultiFab &P_old,
                          MultiFab &P_new,
                          MultiFab * rho_half, 
                          MultiFab& dsdt, 
                          SyncRegister * crse_sync_reg, 
                          SyncRegister * fine_sync_reg,  
                          int crse_dt_ratio, int ** bc, int iteration,
                          Real divu_minus_s_factor,
                          MultiFab &divuold, int have_divu) 
{
    if (ParallelDescriptor::IOProcessor() && verbose)
    {
	cout << "... level projector at level " << level << NL;
    }
    
    if (sync_proj == 0) bldSyncProject();
    
    int rzflag = CoordSys::IsRZ();
    
    // old time velocity has bndry values already
    // must gen valid bndry data for new time velocity.
    // must fill bndry cells in pressure with computable values
    // even though they are not used in calculation
    U_new.setBndry(bogus_value,Xvel,BL_SPACEDIM);
    P_old.setBndry(bogus_value);
    P_new.setBndry(bogus_value);
    LevelData[level].setPhysBoundaryValues(State_Type,Xvel,BL_SPACEDIM);

    const Real* dx = geom.CellSize();
    const BoxArray& grids = LevelData[level].boxArray();
    const BoxArray& P_grids = P_old.boxArray();

    if (level != 0)
    {
	LevelData[level].FillCoarsePatch(P_new,0,cur_pres_time,Press_Type,0,1);
    }

    // set up outflow bcs, BEFORE manipulating state, pressure data
#if (BL_SPACEDIM == 2)
    int outflow_at_top = (phys_bc->lo(0) != Outflow &&
			  phys_bc->lo(1) != Outflow &&
			  phys_bc->hi(0) != Outflow &&
			  phys_bc->hi(1) == Outflow);
    if (outflow_at_top && have_divu && do_outflow_bcs) 
    {
	set_level_projector_outflow_bcs(level,cur_state_time,P_new);
    }
#endif
    
    //
    // Convert Unew to acceleration and Pnew to an update.
    //
    if(level == 0) 
    {
	for(MultiFabIterator P_newmfi(P_new); P_newmfi.isValid(); ++P_newmfi) 
	{
	    DependentMultiFabIterator U_newmfi(P_newmfi, U_new);
	    DependentMultiFabIterator U_oldmfi(P_newmfi, U_old);
	    ConvertUnew( U_newmfi(), U_oldmfi(), dt, grids[P_newmfi.index()] );
	    P_newmfi().setVal(0.0);
	}
    } 
    else 
    {
	for(MultiFabIterator P_newmfi(P_new); P_newmfi.isValid(); ++P_newmfi) 
	{
	    DependentMultiFabIterator P_oldmfi(P_newmfi, P_old);
	    DependentMultiFabIterator U_newmfi(P_newmfi, U_new);
	    DependentMultiFabIterator U_oldmfi(P_newmfi, U_old);
	    
	    ConvertUnew(U_newmfi(), U_oldmfi(), dt, grids[P_newmfi.index()]);
	    P_newmfi().minus(P_oldmfi());
	    Box tempbox(P_newmfi().box());
	    tempbox.grow(-2);
	    P_newmfi().setVal(0.0,tempbox,0,1);
	}
    }
    
    if(divu_minus_s_factor>0.0 && divu_minus_s_factor<=1.0 && have_divu) 
  {
      // compute relaxation terms to account for approximate projection
      // add divu_old*divu...factor/dt to dsdt
      Real uoldfactor = divu_minus_s_factor*dt/parent->dtLevel(0);
      UpdateArg1( dsdt, uoldfactor/dt, divuold, 1, grids, 1 );
      
      // add U_old*divu...factor/dt to U_new
      UpdateArg1( U_new, uoldfactor/dt, U_old, BL_SPACEDIM, grids, 1 );
  }

  // scale the projection variables
  rho_half->setBndry(bogus_value);
  scaleVar( rho_half, 1, &U_new, grids, level );

  // application specific first guess
  ProjFirstGuess( U_new, P_new, level, grids );

  // Enforce periodicity of U_new and rho_half (i.e. coefficient of G phi)
  // *after* everything has been done to them
  EnforcePeriodicity( U_new,     BL_SPACEDIM, grids, geom );
  EnforcePeriodicity( *rho_half, 1,           grids, geom );
        
  // add the contribution from the un-projected V to syncregisters
  // ----------------------------------------------------
  int rz_flag = (CoordSys::IsRZ() ? 1 : 0);
  int ngrids = grids.length();
  if (do_sync_proj) 
  {

    int isrz = CoordSys::IsRZ();
    int bcxlo = phys_bc->lo(0);
    int bcxhi = phys_bc->hi(0);
    int lowfix = (isrz==1 && bcxlo!=EXTRAP && bcxlo!=HOEXTRAP);
    int  hifix = (isrz==1 && bcxhi!=EXTRAP && bcxhi!=HOEXTRAP);

    if (level < finest_level) 
    {   // init sync registers between level and level+1
      crse_sync_reg->CrseDVInit(U_new,geom,rz_flag,bc);
      if (have_divu) 
        crse_sync_reg->CrseDsdtAdd(dsdt, geom, rz_flag, bc, 
                                   lowfix, hifix);
    } 
    if (level > 0) 
    { // increment sync registers between level and level-1
      Real invrat = 1.0/(double)crse_dt_ratio;
      const Geometry& crse_geom = parent->Geom(level-1);
      fine_sync_reg->FineDVAdd(U_new,dx,crse_geom,rz_flag,bc,invrat);
      if (have_divu) 
        fine_sync_reg->FineDsdtAdd(dsdt,geom,rz_flag,bc,lowfix,
                                   hifix,invrat);
    }

  }

  // setup projection (note that u_real is a temporary copy)
  // ----------------------------------------------------

  // build amr_real data structures that projection wants
  //Array<BoxArray>& full_mesh = sync_proj->mesh();
  PArray<MultiFab> u_real[BL_SPACEDIM];
  PArray<MultiFab> p_real(level+1), s_real(level+1);
  int n;
  for (n = 0; n < BL_SPACEDIM; n++) 
  {
    u_real[n].resize(level+1);
    u_real[n].set(level, new MultiFab(grids, 1, 1));
    for(MultiFabIterator u_realmfi(u_real[n][level]); u_realmfi.isValid();
        ++u_realmfi)
    {
      DependentMultiFabIterator U_newmfi(u_realmfi, U_new);
      //u_real[n][level][i].copy(U_new[i], n, 0);
      u_realmfi().copy(U_newmfi(), n, 0);
    }
  }
  p_real.set(level, &P_new);
  s_real.set(level, rho_half);

  // project
  // ----------------------------------------------------

  if (!have_divu) 
  {
    // for divu=0 only
    sync_proj->project(u_real, p_real, null_amr_real, s_real, (Real*)dx,
                       proj_tol, level, level, proj_abs_tol);
  } 
  else 
  {
    // for general divu
    bool use_u = true;
    const int ngrids = grids.length();
    int nghost = 1; // required by aliaslib--rbp
    MultiFab rhs_cc(grids,1,nghost,Fab_allocate);
    rhs_cc.setVal(0.0);
    rhs_cc.copy(dsdt,0,0,1,nghost);
    radMult(level,rhs_cc,0);    
    for (int i=0;i<ngrids;i++) 
    {
#if (BL_SPACEDIM == 3)
      Real dsdtmin=dsdt[i].min();
      Real dsdtmax=dsdt[i].max();
      if(dsdtmin!=dsdtmax || dsdtmin!= 0.0) 
      {
        cout << "Projection::level_project: WARNING not yet " <<
          "implemented for 3-d, non-zero divu\n";
        ParallelDescriptor::Abort("Exiting.");
      }
#endif
    }
    rhs_cc.mult(-1.0,0,1,nghost);

    PArray<MultiFab> rhs_real(level+1);
    rhs_real.set(level, &rhs_cc);
    sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real, s_real,
                              use_u, (Real*)dx,
                              proj_tol, level, level, proj_abs_tol);
  }

  // copy and delete u_real
  // ----------------------------------------------------

  for (n = 0; n < BL_SPACEDIM; n++) 
  {
    //for (i = 0; i < grids.length(); i++) {
      //U_new[i].copy(u_real[n][level][i], 0, n);
    //}
    for(MultiFabIterator u_realmfi(u_real[n][level]); u_realmfi.isValid();
        ++u_realmfi)
    {
      DependentMultiFabIterator U_newmfi(u_realmfi, U_new);
      U_newmfi().copy(u_realmfi(), 0, n);
    }
    delete u_real[n].remove(level);
  }

  // increment sync registers by adding contribution from the full D(sig Gphi)
  // ----------------------------------------------------
  if (do_sync_proj) 
  {
    if (level < finest_level) 
    { // init sync registers between level and level+1
      crse_sync_reg->CrseLPhiAdd(P_new,*rho_half,geom,rz_flag);
    }
    if (level > 0) 
    { // increment sync registers between level and level-1
      Real invrat = 1.0/(double)crse_dt_ratio;
      const Geometry& crse_geom = parent->Geom(level-1);
      fine_sync_reg->FineLPhiAdd(P_new,*rho_half,dx,crse_geom,rz_flag,invrat);
    }
  }

  //----------------- reset state + pressure data ---------------------

  // unscale level projection variables
  rescaleVar( rho_half, 1, &U_new, grids, level );
  
  // put U_new back to "normal"
  // --> subtract U_old*divu...factor/dt from U_new
  if(divu_minus_s_factor>0.0 && divu_minus_s_factor<=1.0 && have_divu) 
  {
      Real uoldfactor = -divu_minus_s_factor*dt/parent->dtLevel(0);
      UpdateArg1( U_new, uoldfactor/dt, U_old, BL_SPACEDIM, grids, 1 );
  }

  // convert U back to a velocity, and phi into p^n+1/2
  UnConvertUnew( U_old, dt, U_new, grids );  // un = uo+dt*un
  AddPhi( P_new, P_old, grids );             // pn = pn + po

  // filter variables
  filterUandP( level, P_new, U_new, rho_half, grids, dx, dt );
}



// ----------------------------------------------------------------
//  This function is an attempt at improving the creation of fine
//  pressure data after regridding
// ----------------------------------------------------------------
void Projection::harmonic_project(int level, Real dt, Real cur_pres_time,
                                  const Geometry& geom, MultiFab &P_old)
{
  if (level == 0) 
  {
    cout << "NOT SUPPOSED TO BE IN HARMONIC AT LEVEL 0 \n";
    ParallelDescriptor::Abort("Exiting.");
  }

  if (verbose) 
  {
    cout << "... harmonic projector\n";
  }

  //----------------- manipulate state + pressure data ---------------------

  if (sync_proj == 0) bldSyncProject();

  int rzflag = CoordSys::IsRZ();

  const BoxArray& grids = LevelData[level].boxArray();
  const BoxArray& P_grids = P_old.boxArray();
  MultiFab* rhs      = new MultiFab(P_grids,1,1,Fab_allocate);
  MultiFab* harm_phi = new MultiFab(P_grids,1,1,Fab_allocate);
  MultiFab* temp_phi = new MultiFab(P_grids,1,1,Fab_allocate);
  MultiFab* rho      = new MultiFab(grids,1,1,Fab_allocate);
  MultiFab* harm_vel = new MultiFab(grids,BL_SPACEDIM,1,Fab_allocate);

  rhs->setVal(0.);
  harm_phi->setVal(0.);
  harm_vel->setVal(0.);
  rho->setVal(1.);

  harm_phi->setBndry(bogus_value);

  Real prev_pres_time = cur_pres_time - dt;

  LevelData[level].FillCoarsePatch(*temp_phi,0, prev_pres_time,Press_Type,0,1);
  LevelData[level].FillCoarsePatch(*harm_phi,0, cur_pres_time,Press_Type,0,1);

  for(MultiFabIterator temp_phimfi(*temp_phi); temp_phimfi.isValid(); ++temp_phimfi)
  {
    DependentMultiFabIterator harm_phimfi(temp_phimfi, *harm_phi);
    harm_phimfi().minus(temp_phimfi());
    Box tempbox(harm_phimfi().box());
    tempbox.grow(-2);
    harm_phimfi().setVal(0.,tempbox,0,1);
  }

  delete temp_phi;

  rho->setBndry(bogus_value);
  scaleVar(rho,1,0,grids,level);

  const Real* dx = geom.CellSize();

  // build alias lib structures
  //----------------------------------------------------------
  //Array<BoxArray>& full_mesh = sync_proj->mesh();
  PArray<MultiFab> u_real[BL_SPACEDIM];
  PArray<MultiFab> p_real(level+1), s_real(level+1), rhs_real(level+1);
  int n;
  for (n = 0; n < BL_SPACEDIM; n++) 
  {
    u_real[n].resize(level+1);
    u_real[n].set(level, new MultiFab(grids, 1, 1));
    for(MultiFabIterator u_realmfi(u_real[n][level]); u_realmfi.isValid();
        ++u_realmfi)
    {
      DependentMultiFabIterator harm_velmfi(u_realmfi, *harm_vel);
      //u_real[n][level][i].copy((*harm_vel)[i], n, 0);
      u_realmfi().copy(harm_velmfi(), n, 0);
    }
  }
  p_real.set(level, harm_phi);
  s_real.set(level, rho);
  rhs_real.set(level, rhs);

  // project
  //----------------------------------------------------------
  bool use_u = false;
  sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real, s_real,
                            use_u, (Real*)dx,
                            proj_tol, level, level, proj_abs_tol);

  // copy and delete u_real
  // ----------------------------------------------------

  for (n = 0; n < BL_SPACEDIM; n++) 
  {
    // copy unnecessary since harm_vel will be discarded anyway
    delete u_real[n].remove(level);
  }

  //----------------- reset state + pressure data ---------------------

  // unscale variables for harmonic projection
  rescaleVar(rho,1,0,grids,level);

  // update pressure
  AddPhi( P_old, *harm_phi, grids );

  delete rhs;
  delete harm_phi;
  delete harm_vel;
  delete rho;
}




//----------------------------------------------------------
//  SYNC_PROJECT
//----------------------------------------------------------

#define MAXLEV 10
void Projection::syncProject(int c_lev, MultiFab & pres, MultiFab & vel,
                             MultiFab * rho_half, MultiFab * Vsync,
                             MultiFab & phi,
                             SyncRegister * rhs_sync_reg,
                             SyncRegister * crsr_sync_reg,
                             const BoxArray& sync_boxes,  int ** sync_bc,
                             const Geometry& geom,
                             const Real * dx, Real dt_crse, int crse_dt_ratio)
{
  RunStats proj_stats("sync_project",c_lev);
  proj_stats.start();

  int rz_flag = (CoordSys::IsRZ() ? 1 : 0);

  if (verbose) 
  {
    cout << "SyncProject: level = " << c_lev
         << " correction to level " << finest_level << NL;
  }

  //----------------- manipulate state + pressure data ---------------------

  if (sync_proj == 0) bldSyncProject();

  int rzflag = CoordSys::IsRZ();

      // gather data
  const BoxArray& grids = LevelData[c_lev].boxArray();
  const BoxArray& P_grids = pres.boxArray();
  MultiFab  rhs(P_grids,1,1);
  MultiFab& sig = *rho_half;
  rhs_sync_reg->InitRHS(rhs,geom,phys_bc);

  phi.setVal(0.0);

#ifdef SET_BOGUS_BNDRY
  sig.setBndry(bogus_value);
#endif

  // scale sync projection variables
  scaleVar(&sig,1,Vsync,grids,c_lev);

  //  If this sync project is not at level 0 then we need to account for
  //  the changes made here in the level c_lev velocity in the sync registers
  //  going into the level (c_lev-1) sync project.
  //----------------------------------------------------------
  if (c_lev > 0) 
  {
    Real invrat = 1.0/(double)crse_dt_ratio;
    const Geometry& crsr_geom = parent->Geom(c_lev-1);
    crsr_sync_reg->CompDVAdd(*Vsync,sync_boxes,dx,geom,crsr_geom,
                             rz_flag,sync_bc,invrat);
  }

  // if periodic, copy into periodic translates of Vsync
  EnforcePeriodicity( *Vsync, BL_SPACEDIM, grids, geom );

  // build aliaslib structures
  //----------------------------------------------------------
  //Array<BoxArray>& full_mesh = sync_proj->mesh();
  int n;
  PArray<MultiFab> u_real[BL_SPACEDIM];
  PArray<MultiFab> p_real(c_lev+1), s_real(c_lev+1), rhs_real(c_lev+1);
  for (n = 0; n < BL_SPACEDIM; n++) 
  {
    u_real[n].resize(c_lev+1);
    u_real[n].set(c_lev, new MultiFab(grids, 1, 1));
    for(MultiFabIterator u_realmfi(u_real[n][c_lev]); u_realmfi.isValid();
        ++u_realmfi)
    {
      DependentMultiFabIterator Vsyncmfi(u_realmfi, *Vsync);
      u_realmfi().copy(Vsyncmfi(), n, 0);
    }
  }
  p_real.set(c_lev, &phi);
  s_real.set(c_lev, &sig);
  rhs_real.set(c_lev, &rhs);

  //  PROJECT
  //  if use_u = 0, then solves DGphi = RHS
  //  if use_u = 1, then solves DGphi = RHS + DV
  //  both return phi and (V-Gphi) as V
  //----------------------------------------------------------
  bool use_u = true;
  sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real, s_real,
                            use_u, (Real*)dx,
                            sync_tol, c_lev, c_lev, proj_abs_tol);

  // copy and delete u_real
  // ----------------------------------------------------

  for (n = 0; n < BL_SPACEDIM; n++) 
  {
    for(MultiFabIterator u_realmfi(u_real[n][c_lev]); u_realmfi.isValid();
        ++u_realmfi)
    {
      DependentMultiFabIterator Vsyncmfi(u_realmfi, *Vsync);
      Vsyncmfi().copy(u_realmfi(), 0, n);
    }
    delete u_real[n].remove(c_lev);
  }

  //  If this sync project is not at level 0 then we need to account for
  //  the changes made here in the level c_lev velocity in the sync registers
  //  going into the level (c_lev-1) sync project.  Note that this must be
  //  done before rho_half is scaled back.
  //----------------------------------------------------------
  if (c_lev > 0) 
  {
    Real invrat = 1.0/(double)crse_dt_ratio;
    const Geometry& crsr_geom = parent->Geom(c_lev-1);
    crsr_sync_reg->CompLPhiAdd(phi,sig,sync_boxes,
                               dx,geom,crsr_geom,rz_flag,invrat);
  }

  //----------------- reset state + pressure data ---------------------

  // unscale the sync projection variables for rz
  rescaleVar(&sig,1,Vsync,grids,c_lev);

  // add projected Vsync to new velocity at this level
  // and add phi to pressure
  AddPhi( pres, phi, grids );
  UpdateArg1( vel, dt_crse, *Vsync, BL_SPACEDIM, grids, 1 );

  proj_stats.end();
}




// ------------------------------------------------------------
//  MULTI-LEVEL SYNC_PROJECT
// ------------------------------------------------------------

void Projection::MLsyncProject(int c_lev,
                               MultiFab & pres_crse, MultiFab & vel_crse,
                               MultiFab & pres_fine, MultiFab & vel_fine,
                               MultiFab & rho_crse,  MultiFab & rho_fine,
                               MultiFab * Vsync, MultiFab & V_corr,
                               MultiFab & phi_fine,
                               SyncRegister * rhs_sync_reg,
                               SyncRegister * crsr_sync_reg,
                               int ** sync_bc,
                               const Real * dx, Real dt_crse, 
                               IntVect& ratio, int crse_dt_ratio,
                               const Geometry & fine_geom,
                               const Geometry & crse_geom)
{
  RunStats proj_stats("sync_project",c_lev);
  proj_stats.start();
    
  int lev;
  if (verbose) 
  {
    cout << "SyncProject: levels = " << c_lev << ", " << c_lev+1 << NL;
  }
    
  int rz_flag = (CoordSys::IsRZ() ? 1 : 0);
  if (sync_proj == 0) bldSyncProject();

    
  // Set up memory
  // -------------------------------------------------------------
  MultiFab *phi[MAXLEV];
    
  const BoxArray& Pgrids_crse = pres_crse.boxArray();
  phi[c_lev]    = new MultiFab(Pgrids_crse,1,1,Fab_allocate);
  phi[c_lev]->setVal(0.);
    
  MultiFab *crse_rhs = new MultiFab(Pgrids_crse,1,1,Fab_allocate);
    
  const BoxArray& Pgrids_fine = pres_fine.boxArray();
  phi[c_lev+1] = new MultiFab(Pgrids_fine,1,1,Fab_allocate);
  phi[c_lev+1]->setVal(0.);
    
  BoxArray sync_boxes = pres_fine.boxArray();
  sync_boxes.coarsen(ratio);
    
  // Set up RHS
  // -------------------------------------------------------------
  crse_rhs->setVal(0.);
  rhs_sync_reg->InitRHS(*crse_rhs,crse_geom,phys_bc);
    
  Box P_finedomain(surroundingNodes(crse_geom.Domain()));
  P_finedomain.refine(ratio);
  if (Pgrids_fine[0] == P_finedomain) crse_rhs->setVal(0.);
    
  const BoxArray& grids      = LevelData[c_lev].boxArray();
  const BoxArray& fine_grids = LevelData[c_lev+1].boxArray();
    
  // Do necessary scaling
  // -------------------------------------------------------------
  scaleVar( &rho_crse, 0, Vsync,   grids,      c_lev   );
  scaleVar( &rho_fine, 0, &V_corr, fine_grids, c_lev+1 );
    
  // -------------------------------------------------------------
  //  If this sync project is not at level 0 then we need to account for
  //  the changes made here in the level c_lev velocity in the sync registers
  //  going into the level (c_lev-1) sync project.
  // -------------------------------------------------------------
  if (c_lev > 0) 
  {
    Real invrat = 1.0/(double)crse_dt_ratio;
    const Geometry& crsr_geom = parent->Geom(c_lev-1);
    crsr_sync_reg->CompDVAdd(*Vsync,sync_boxes,dx,
                             crse_geom,crsr_geom,
                             rz_flag,sync_bc,invrat);
  }

  // set up alias lib
  // -------------------------------------------------------------
  //Array<BoxArray>& full_mesh = sync_proj->mesh();
  PArray<MultiFab> u_real[BL_SPACEDIM];
  PArray<MultiFab> p_real(c_lev+2), s_real(c_lev+2), crse_rhs_real(c_lev+1);
  int n;
  for (n = 0; n < BL_SPACEDIM; n++) 
  {
    u_real[n].resize(c_lev+2);
    u_real[n].set(c_lev,   new MultiFab(grids,      1, 1));
    u_real[n].set(c_lev+1, new MultiFab(fine_grids, 1, 1));
    //for (i = 0; i < grids.length(); i++) {
      //u_real[n][c_lev][i].copy((*Vsync)[i], n, 0);
    //}
    //for (i = 0; i < fine_grids.length(); i++) {
      //u_real[n][c_lev+1][i].copy(V_corr[i], n, 0);
    //}
    for(MultiFabIterator u_realmfi(u_real[n][c_lev]); u_realmfi.isValid();
        ++u_realmfi)
    {
      DependentMultiFabIterator Vsyncmfi(u_realmfi, *Vsync);
      u_realmfi().copy(Vsyncmfi(), n, 0);
    }
    for(MultiFabIterator u_realfinemfi(u_real[n][c_lev+1]); u_realfinemfi.isValid();
        ++u_realfinemfi)
    {
      DependentMultiFabIterator V_corrmfi(u_realfinemfi, V_corr);
      u_realfinemfi().copy(V_corrmfi(), n, 0);
    }

    restrict_level(u_real[n][c_lev], u_real[n][c_lev+1], ratio, default_restrictor(), level_interface(), 0);
  }

  s_real.set(c_lev,   &rho_crse);
  s_real.set(c_lev+1, &rho_fine);

  restrict_level(s_real[c_lev], s_real[c_lev+1], ratio, default_restrictor(), level_interface(), 0);

  p_real.set(c_lev,   phi[c_lev]);
  p_real.set(c_lev+1, phi[c_lev+1]);

  //  Note that crse_rhs_real is only built on the coarsest level
  crse_rhs_real.set(c_lev, crse_rhs);

  // -------------------------------------------------------------
  //  The Multilevel Projection
  //  if use_u = 0, then solves DGphi = RHS
  //  if use_u = 1, then solves DGphi = RHS + DV
  //  both return phi and (V-Gphi) as V
  // -------------------------------------------------------------
  bool use_u            = true;
  const Real* dx_fine  = parent->Geom(c_lev+1).CellSize();

  sync_proj->manual_project(u_real, p_real, null_amr_real,
                            crse_rhs_real, s_real,
                            use_u, (Real*)dx_fine,
                            sync_tol, c_lev, c_lev+1, proj_abs_tol);

  // copy and delete u_real
  // ----------------------------------------------------

  for (n = 0; n < BL_SPACEDIM; n++) 
  {
    //for (i = 0; i < grids.length(); i++) {
      //(*Vsync)[i].copy(u_real[n][c_lev][i], 0, n);
    //}
    //for (i = 0; i < fine_grids.length(); i++) {
      //V_corr[i].copy(u_real[n][c_lev+1][i], 0, n);
    //}
    for(MultiFabIterator u_realmfi(u_real[n][c_lev]); u_realmfi.isValid();
        ++u_realmfi)
    {
      DependentMultiFabIterator Vsyncmfi(u_realmfi, *Vsync);
      Vsyncmfi().copy(u_realmfi(), 0, n);
    }
    for(MultiFabIterator u_realfinemfi(u_real[n][c_lev+1]); u_realfinemfi.isValid();
        ++u_realfinemfi)
    {
      DependentMultiFabIterator V_corrmfi(u_realfinemfi, V_corr);
      V_corrmfi().copy(u_realfinemfi(), 0, n);
    }
    delete u_real[n].remove(c_lev);
    delete u_real[n].remove(c_lev+1);
  }

  // -------------------------------------------------------------
  //  If this sync project is not at levels 0-1 then we need to account for
  //  the changes made here in the level c_lev velocity in the sync registers
  //  going into the level (c_lev-1) sync project.  Note that this must be
  //  done before rho_half is scaled back.
  // -------------------------------------------------------------
  if (c_lev > 0) 
  {
    Real invrat = 1.0/(double)crse_dt_ratio;
    MultiFab * phi_crse = phi[c_lev];
    const Geometry& crsr_geom = parent->Geom(c_lev-1);
    crsr_sync_reg->CompLPhiAdd(*phi_crse,rho_crse,sync_boxes,
                               dx,crse_geom,crsr_geom,rz_flag,invrat);
  }


  // Do necessary un-scaling
  // -------------------------------------------------------------
  rescaleVar( &rho_crse, 0, Vsync,   grids,      c_lev   );
  rescaleVar( &rho_fine, 0, &V_corr, fine_grids, c_lev+1 );


  // Add projected vel to new velocity and add phi to pressure 
  // -------------------------------------------------------------
  AddPhi( pres_crse, *phi[c_lev],   grids      );
  AddPhi( pres_fine, *phi[c_lev+1], fine_grids );

  UpdateArg1( vel_crse, dt_crse, *Vsync, BL_SPACEDIM, grids,      1 );
  UpdateArg1( vel_fine, dt_crse, V_corr, BL_SPACEDIM, fine_grids, 1 );


  // Clean up memory
  // -------------------------------------------------------------
  //for ( i = 0; i < fine_grids.length(); i++ ) {
    //FArrayBox &e = (*phi[c_lev+1])[i];
    //phi_fine[i].copy(e,0,0,1);
  //}
  for(MultiFabIterator phimfi(*phi[c_lev+1]); phimfi.isValid(); ++phimfi) 
  {
    DependentMultiFabIterator phi_finemfi(phimfi, phi_fine);
    phi_finemfi().copy(phimfi(),0,0,1);
  }
  for (lev = c_lev; lev <= c_lev+1; lev++) 
  {
    delete phi[lev];
  }

  delete crse_rhs;

  proj_stats.end();
}




//------------------------------------------------------------
//  INITIAL_VELOCITY_PROJECT
//  the initial velocity projection in post_init.
//  this function ensures that the velocities are nondivergent
//------------------------------------------------------------

void Projection::initialVelocityProject(int c_lev,
                                        Real cur_divu_time, 
                                        int have_divu)
{
    //----------------- manipulate state + pressure data ---------------------
  RunStats proj_stats("sync_project",c_lev);
  proj_stats.start();

  int lev;
  int f_lev = finest_level;
  if (verbose) 
  {
    cout << "initialVelocityProject: levels = " << c_lev << "  "
         << f_lev << NL;
    if (rho_wgt_vel_proj) 
    {
      cout << "RHO WEIGHTED INITIAL VELOCITY PROJECTION\n";
    } 
    else 
    {
      cout << "CONSTANT DENSITY INITIAL VELOCITY PROJECTION\n";
    }
  }

  if (sync_proj == 0) bldSyncProject();

  int rzflag = CoordSys::IsRZ();
  MultiFab *vel[MAX_LEV];
  MultiFab *phi[MAX_LEV];
  MultiFab *sig[MAX_LEV];
  int rho_comp = 0;
  if (rho_wgt_vel_proj) rho_comp = Density;

  for (lev = c_lev; lev <= f_lev; lev++) 
  {
    MultiFab &P_old = LevelData[lev].get_old_data(Press_Type);
    P_old.setVal(0.0);
  }

  for (lev = c_lev; lev <= f_lev; lev++) 
  {
    vel[lev] = &LevelData[lev].get_new_data(State_Type);
    phi[lev] = &LevelData[lev].get_old_data(Press_Type);
    if (rho_wgt_vel_proj) 
    {
      sig[lev] = &LevelData[lev].get_new_data(State_Type);
    } 
    else 
    {
      const BoxArray& grids = LevelData[lev].boxArray();
      const int nghost = 1;
      sig[lev] = new MultiFab(grids,1,nghost);
      sig[lev]->setVal(1.0,nghost);
    }
  }

  // ------------------- Set boundary conditions + scale variables ----------
  // set up outflow bcs

#if (BL_SPACEDIM == 2)
  int outflow_at_top = (phys_bc->lo(0) != Outflow &&
                        phys_bc->lo(1) != Outflow && 
                        phys_bc->hi(0) != Outflow &&
                        phys_bc->hi(1) == Outflow);
  if (outflow_at_top && have_divu && do_outflow_bcs)
  {
      set_initial_projection_outflow_bcs(vel,sig,phi,c_lev,cur_divu_time);
  }
#endif

  for (lev = c_lev; lev <= f_lev; lev++) 
  {
      
#ifdef SET_BOGUS_BNDRY
      vel[lev]->setBndry(bogus_value,Xvel,BL_SPACEDIM);
      sig[lev]->setBndry(bogus_value,rho_comp,1);
#endif

      // set the physical boundary values
      AmrLevel& amr_level   = parent->getLevel(lev);
      const BoxArray& grids = amr_level.boxArray();
      amr_level.setPhysBoundaryValues(State_Type,Xvel,BL_SPACEDIM);
      if (rho_wgt_vel_proj) 
      {
          amr_level.setPhysBoundaryValues(State_Type,Density,1);
      }

      if (have_divu) 
      {
	  int Divu_Type, Divu;
	  if (!LevelData[lev].isStateVariable("divu", Divu_Type, Divu)) 
	  {
	      BoxLib::Error("Projection::initialVelocityProject(): Divu not found");
	  }
// make sure ghost cells are properly filled
	  MultiFab &divu_new = amr_level.get_new_data(Divu_Type);
	  divu_new.FillBoundary();
	  amr_level.setPhysBoundaryValues(Divu_Type,0,1,cur_divu_time);
      }
      
      
      // scale the projection variables
      scaleVar(sig[lev],1,vel[lev],grids,lev);
  }

  
  // setup alias lib
  //-------------------------------------------------------
  const Array<BoxArray>& full_mesh = sync_proj->mesh();
  PArray<MultiFab> u_real[BL_SPACEDIM];
  PArray<MultiFab> p_real(f_lev+1), s_real(f_lev+1);
  int n;
  for (n = 0; n < BL_SPACEDIM; n++) 
  {
    u_real[n].resize(f_lev+1);
  }
  for (lev = c_lev; lev <= f_lev; lev++) 
  {
    for (n = 0; n < BL_SPACEDIM; n++) 
    {
      u_real[n].set(lev, new MultiFab(full_mesh[lev], 1, 1));
      //for (i = 0; i < full_mesh[lev].length(); i++) {
        //u_real[n][lev][i].copy((*vel[lev])[i], Xvel+n, 0);
      //}
      for(MultiFabIterator u_realmfi(u_real[n][lev]); u_realmfi.isValid();
          ++u_realmfi)
      {
        DependentMultiFabIterator velmfi(u_realmfi, *vel[lev]);
        u_realmfi().copy(velmfi(), Xvel+n, 0);
      }
    }
    p_real.set(lev, phi[lev]);
    if (rho_comp == 0) 
    {
      s_real.set(lev, sig[lev]);
    }
    else 
    {
      s_real.set(lev, new MultiFab(full_mesh[lev], 1, 1));
      //for (i = 0; i < full_mesh[lev].length(); i++) {
        //s_real[lev][i].copy((*sig[lev])[i], rho_comp, 0);
      //}
      for(MultiFabIterator s_realmfi(s_real[lev]); s_realmfi.isValid();
          ++s_realmfi)
      {
        DependentMultiFabIterator sigmfi(s_realmfi, *sig[lev]);
        s_realmfi().copy(sigmfi(), rho_comp, 0);
      }
    }
  }

  // project
  //-------------------------------------------------------
  const Real* dx_lev = parent->Geom(f_lev).CellSize();
  if(!have_divu) 
  {
    // zero divu only
    sync_proj->project(u_real, p_real, null_amr_real, s_real, (Real*)dx_lev,
                       proj_tol, c_lev, f_lev, proj_abs_tol);
  } 
  else 
  {
    // general divu
    MultiFab *rhs_cc[MAX_LEV];
    int nghost = 1; 
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
      AmrLevel& amr_level   = parent->getLevel(lev);
      const BoxArray& grids = amr_level.boxArray();
      rhs_cc[lev]  = new MultiFab(grids,1,nghost,Fab_allocate);
      MultiFab* rhslev = rhs_cc[lev];
      put_divu_in_cc_rhs(*rhslev,lev,grids,cur_divu_time);
      rhslev->mult(-1.0,0,1,nghost);
      radMult(lev,*rhslev,0); 
    }
    bool use_u = true;
    PArray<MultiFab> rhs_real(f_lev+1);
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
      rhs_real.set(lev, rhs_cc[lev]);
    }
    sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real, s_real,
                              use_u, (Real*)dx_lev,
                              proj_tol, c_lev, f_lev, proj_abs_tol);
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
      delete rhs_cc[lev];
    }
  }

  // copy and delete u_real, delete s_real if appropriate
  // ----------------------------------------------------

  for (lev = c_lev; lev <= f_lev; lev++) 
  {
    for (n = 0; n < BL_SPACEDIM; n++) 
    {
      //for (i = 0; i < full_mesh[lev].length(); i++) {
        //(*vel[lev])[i].copy(u_real[n][lev][i], 0, Xvel+n);
      //}
      for(MultiFabIterator u_realmfi(u_real[n][lev]); u_realmfi.isValid();
          ++u_realmfi)
      {
        DependentMultiFabIterator velmfi(u_realmfi, *vel[lev]);
        velmfi().copy(u_realmfi(), 0, Xvel+n);
      }
      delete u_real[n].remove(lev);
    }
    if (rho_comp != 0)
      delete s_real.remove(lev);
  }

  //----------------- reset state + pressure data ---------------------
  // unscale initial projection variables
  for (lev = c_lev; lev <= f_lev; lev++) 
  {
      AmrLevel& amr_level   = parent->getLevel(lev);
      const BoxArray& grids = amr_level.boxArray();
      rescaleVar(sig[lev],1,vel[lev],grids,lev);
  }

  // delete sigma if not used later
  for (lev = c_lev; lev <= f_lev; lev++) 
  {
      if (!rho_wgt_vel_proj) 
      {
          delete sig[lev];
      }
  }

  proj_stats.end();
}

//
// INITIAL_SYNC_PROJECT
// the velocity projection in post_init, which computes the initial
// pressure used in the timestepping.
//

void
Projection::initialSyncProject (int       c_lev,
                                MultiFab* sig[],
                                Real      dt, 
                                Real      strt_time,
                                Real      dt_init,
                                int       have_divu)
{
    RunStats proj_stats("sync_project",c_lev);
    proj_stats.start();

    int lev;
    int f_lev = finest_level;
    if (verbose) 
    {
        cout << "SyncProject: levels = " << c_lev << "  " << f_lev << NL;
    }
    //
    //----------------- manipulate state + pressure data ---------------------
    //
    if (sync_proj == 0) bldSyncProject();
    //
    // Gather data.
    //
    int rzflag = CoordSys::IsRZ();
    MultiFab* vel[MAX_LEV];
    MultiFab* phi[MAX_LEV];
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        vel[lev] = &LevelData[lev].get_new_data(State_Type);
        phi[lev] = &LevelData[lev].get_old_data(Press_Type);
    }

    MultiFab* rhs[MAX_LEV];

    if (have_divu) 
    {
        //
        // Set up rhs for manual project.
        //
#if (BL_SPACEDIM == 3)
        Real mindsdt = 0.0;
        Real maxdsdt = 0.0;
#endif
        for (lev = c_lev; lev <= f_lev; lev++) 
        {
            AmrLevel& amr_level = parent->getLevel(lev);

            int Divu_Type, Divu;
            if (!LevelData[c_lev].isStateVariable("divu", Divu_Type, Divu)) 
            {
                BoxLib::Error("Projection::initialSyncProject(): Divu not found");
            }
            //
            // Make sure ghost cells are properly filled.
            //
            MultiFab& divu_new = amr_level.get_new_data(Divu_Type);
            divu_new.FillBoundary();
            MultiFab& divu_old = amr_level.get_new_data(Divu_Type);
            divu_old.FillBoundary();
            amr_level.setPhysBoundaryValues(Divu_Type,0,1,strt_time);
            amr_level.setPhysBoundaryValues(Divu_Type,0,1,strt_time+dt);

            const BoxArray& grids = amr_level.boxArray();
            const int nghost      = 1;
            rhs[lev] = new MultiFab(grids,1,nghost,Fab_allocate);
            MultiFab* rhslev = rhs[lev];
            rhslev->setVal(0.0);

            MultiFab* divu = getDivCond(lev,nghost,strt_time);
            MultiFab* dsdt = getDivCond(lev,nghost,strt_time+dt);

            for (MultiFabIterator mfi(*rhslev); mfi.isValid(false); ++mfi)
            {
                DependentMultiFabIterator divu_it(mfi,*divu);
                DependentMultiFabIterator dsdt_it(mfi,*dsdt);
                dsdt_it().minus(divu_it());
                dsdt_it().divide(dt);
                mfi().copy(dsdt_it());
#if (BL_SPACEDIM == 3)
                mindsdt = Min(mindsdt, dsdt_it().min());
                maxdsdt = Max(maxdsdt, dsdt_it().max());
#endif
            }

            delete divu;
            delete dsdt;

            rhslev->mult(-1.0,0,1);
            if (CoordSys::IsRZ()) 
                radMult(lev,*rhslev,0);    
        }
#if  (BL_SPACEDIM == 3)
        if (mindsdt != maxdsdt || mindsdt != 0.0) 
        {
            cout << "Projection::initialSyncProject: WARNING not yet " <<
                "implemented for 3-d, non-zero divu\n";
            ParallelDescriptor::Abort("Bye.");
        }
#endif
    }

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        MultiFab &P_old = LevelData[lev].get_old_data(Press_Type);
        P_old.setVal(0.0);
    }

#if (BL_SPACEDIM == 2)
    int outflow_at_top = phys_bc->lo(0) != Outflow && phys_bc->lo(1) != Outflow && 
        phys_bc->hi(0) != Outflow && phys_bc->hi(1) == Outflow; 
    if (outflow_at_top && have_divu && do_outflow_bcs) 
    {
        set_initial_syncproject_outflow_bcs(phi,c_lev,strt_time,dt);
    }
#endif

#ifdef SET_BOGUS_BNDRY
    //
    // Set velocity bndry values to bogus values.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        vel[lev]->setBndry(bogus_value,Xvel,BL_SPACEDIM);
        MultiFab &u_o = LevelData[lev].get_old_data(State_Type);
        u_o.setBndry(bogus_value,Xvel,BL_SPACEDIM);
        sig[lev]->setBndry(bogus_value);
    }
#endif
    //
    // Convert velocities to accelerations.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        LevelData[lev].setPhysBoundaryValues(State_Type,Xvel,BL_SPACEDIM,1);
        LevelData[lev].setPhysBoundaryValues(State_Type,Xvel,BL_SPACEDIM,0);
        MultiFab &u_o = LevelData[lev].get_old_data(State_Type);
        ConvertUnew(*vel[lev], u_o, dt, LevelData[lev].boxArray());
    }
    //
    // Scale initial sync projection variables.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        AmrLevel& amr_level = parent->getLevel(lev);
        scaleVar(sig[lev],1,vel[lev],amr_level.boxArray(),lev);
    }
    //
    // Build the aliaslib data structures
    //
    const Array<BoxArray>& full_mesh = sync_proj->mesh();
    PArray<MultiFab> u_real[BL_SPACEDIM];
    PArray<MultiFab> p_real(f_lev+1), s_real(f_lev+1);
    int n;
    for (n = 0; n < BL_SPACEDIM; n++) 
    {
        u_real[n].resize(f_lev+1);
    }
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        for (n = 0; n < BL_SPACEDIM; n++) 
        {
            u_real[n].set(lev, new MultiFab(full_mesh[lev], 1, 1));
            for (MultiFabIterator u_realmfi(u_real[n][lev]); u_realmfi.isValid();
                ++u_realmfi)
            {
                DependentMultiFabIterator velmfi(u_realmfi, *vel[lev]);
                u_realmfi->copy(*velmfi, Xvel+n, 0);
            }
        }
        p_real.set(lev, phi[lev]);
        s_real.set(lev, sig[lev]);
    }

    for (n = 0; n < BL_SPACEDIM; n++) 
    {
        for (lev = f_lev; lev >= c_lev+1; lev--) 
        {
            restrict_level(u_real[n][lev-1], u_real[n][lev], parent->refRatio(lev-1),
                           default_restrictor(), level_interface(), 0);
        }
    }

    const Real* dx_lev = parent->Geom(f_lev).CellSize();
    //
    // Project.
    //
    if (!have_divu) 
    {
        //
        // Zero divu only or debugging.
        //
        sync_proj->project(u_real, p_real, null_amr_real, s_real,
                           (Real*)dx_lev, proj_tol, c_lev, f_lev,proj_abs_tol);
    } 
    else 
    {
        //
        // General divu.
        //
        bool use_u = true;
        PArray<MultiFab> rhs_real(f_lev+1);
        for (lev = c_lev; lev <= f_lev; lev++) 
        {
            rhs_real.set(lev, rhs[lev]);
        }
        sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real, s_real,
                                  use_u, (Real*)dx_lev,
                                  proj_tol, c_lev, f_lev, proj_abs_tol);
    }
    //
    // Copy and delete u_real, delete s_real if appropriate.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        for (n = 0; n < BL_SPACEDIM; n++) 
        {
            for (MultiFabIterator u_realmfi(u_real[n][lev]); u_realmfi.isValid();
                ++u_realmfi)
            {
                DependentMultiFabIterator velmfi(u_realmfi, *vel[lev]);
                velmfi().copy(u_realmfi(), 0, Xvel+n);
            }
            delete u_real[n].remove(lev);
        }
    }
    //
    //----------------- reset state + pressure data ---------------------
    //
    // Unscale initial sync projection variables.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        AmrLevel& amr_level   = parent->getLevel(lev);
        rescaleVar(sig[lev],1,vel[lev],amr_level.boxArray(),lev);
    }
    //
    // Add correction at coarse and fine levels.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        incrPress(lev, 1.0);
    }

    proj_stats.end();
}

//
// Public nonprojector member functions follow.
//

//
// Filter the velocities and pressures after a projection.
//

void
Projection::filterUandP (int             level,
                         MultiFab&       P_new,
                         MultiFab&       U_new,
                         MultiFab*       rho_half,
                         const BoxArray& grids,
                         const Real*     dx,
                         const Real      dt )
{
    FArrayBox scratch;
    FArrayBox gradp;

    const Geometry& geom = parent->Geom(level);
    const BOX& domain    = geom.Domain();

    int wrap_around_x = 0;
    if (geom.isPeriodic(0) && grids[0] == domain)
        wrap_around_x = 1;

    int wrap_around_y = 0;
    if (geom.isPeriodic(1) && grids[0] == domain)
        wrap_around_y = 1;

#if (BL_SPACEDIM == 3)
    int wrap_around_z = 0;
    if (geom.isPeriodic(2) && grids[0] == domain)
        wrap_around_z = 1;
#endif

    for (MultiFabIterator P_newmfi(P_new); P_newmfi.isValid(); ++P_newmfi) 
    {
        DependentMultiFabIterator U_newmfi(P_newmfi, U_new);
        DependentMultiFabIterator rho_halfmfi(P_newmfi, *rho_half);

        const int* lo   = grids[P_newmfi.index()].loVect();
        const int* hi   = grids[P_newmfi.index()].hiVect();
        const int* p_lo = P_newmfi().loVect(); 
        const int* p_hi = P_newmfi().hiVect(); 
        
        scratch.resize(P_newmfi().box(),1);
        scratch.setVal(0.0);
        
        FORT_FILTERP(P_newmfi().dataPtr(), scratch.dataPtr(), 
                     ARLIM(p_lo), ARLIM(p_hi),
                     lo,hi,dx,&filter_factor,
                     &wrap_around_x, &wrap_around_y
#if (BL_SPACEDIM == 3)
                     &wrap_around_z
#endif
                     );

        if (filter_u) 
        {
            scratch.mult(dt);
            gradp.resize(grids[P_newmfi.index()],BL_SPACEDIM);
            const Real* gp_dat = gradp.dataPtr();
            FORT_GRADP(scratch.dataPtr(),ARLIM(p_lo),ARLIM(p_hi),
                       gp_dat,ARLIM(lo),ARLIM(hi),lo,hi,dx);
            for (int idim=0;idim<BL_SPACEDIM;idim++)
                gradp.divide(rho_halfmfi(),0,idim,1);
            U_newmfi().minus(gradp,Xvel,Xvel,BL_SPACEDIM);
        }
    }
}

//
// Compute the node-centered divergence of U.
//

void
Projection::computeDV (MultiFab&       DV,
                       const MultiFab& U,
                       int             src_comp,
                       const Real*     dx,
                       int             is_rz)
{
    const Real mult = 1.0;
    
    const BoxArray& U_boxes = U.boxArray();

    FArrayBox ufab;
    for (MultiFabIterator DVmfi(DV); DVmfi.isValid(); ++DVmfi) 
    {
        DependentMultiFabIterator Umfi(DVmfi, U);

        assert(U_boxes[Umfi.index()] == Umfi.validbox());

        ufab.resize(::grow(Umfi.validbox(),1),BL_SPACEDIM);
        ufab.setVal(0.0);
        ufab.copy(Umfi(),ufab.box(),src_comp,ufab.box(),0,BL_SPACEDIM);

        Box ndbox = ::surroundingNodes(Umfi.validbox());

        const int* ndlo = ndbox.loVect();
        const int* ndhi = ndbox.hiVect();
        const int* ulo  = ufab.box().loVect();
        const int* uhi  = ufab.box().hiVect();
        FORT_SRDIVU(ufab.dataPtr(),ARLIM(ulo),ARLIM(uhi),
                    DVmfi().dataPtr(),ARLIM(ndlo),ARLIM(ndhi),
                    ndlo,ndhi,dx,&mult,&is_rz);
    }
}

//
// Put S in the rhs of the projector--node based version.
//

void
Projection::put_divu_in_node_rhs (MultiFab&       rhs,
                                  int             level,
                                  const int&      nghost,
                                  const BoxArray& P_grids,
                                  Real            time,
                                  int             user_rz)
{
    assert(user_rz >= -1 && user_rz <= 1);

    rhs.setVal(0.0);

    const Geometry& geom = parent->Geom(level);

#if (BL_SPACEDIM == 2)
    int bcxlo      = phys_bc->lo(0);
    int bcxhi      = phys_bc->hi(0);
    int isrz       = (user_rz == -1) ? CoordSys::IsRZ() : user_rz;
    int lowfix     = (isrz==1 && bcxlo!=EXTRAP && bcxlo!=HOEXTRAP);
    int hifix      = (isrz==1 && bcxhi!=EXTRAP && bcxhi!=HOEXTRAP);
    const Real* dx = geom.CellSize();
    Real hx = dx[0];
#endif
    const Box& domain = geom.Domain();
    const int* domlo  = domain.loVect();
    const int* domhi  = domain.hiVect();

    FArrayBox divu;

    for (int i = 0; i < P_grids.length(); i++)
    {
        Box divubox = P_grids[i];
        divubox.convert(IntVect::TheCellVector());
        divubox.grow(1);
        divu.resize(divubox,1);
        getDivCond(level,divu,1,time);

#if (BL_SPACEDIM == 2)
        DEF_CLIMITS(divu,divudat,divulo,divuhi);
        DEF_LIMITS(rhs[i],rhsdat,rhslo,rhshi);
        Array<Real> rcen(divubox.length(0),1.0);
        if (isrz == 1) 
        {
            geom.GetCellLoc(rcen,divubox,0);
        } 
        int extrap_edges   = 0;
        int extrap_corners = 1;
        FORT_HGC2N(&nghost,ARLIM(divulo),ARLIM(divuhi),divudat,
                   rcen.dataPtr(), ARLIM(rhslo),ARLIM(rhshi),rhsdat,
                   domlo,domhi,lowfix,hifix,&hx,
                   &extrap_edges, &extrap_corners,&isrz);
#endif
#if (BL_SPACEDIM == 3)
        Real divumin = divu->min();
        Real divumax = divu->max();
        if (divumin != divumax || divumin != 0.0) 
        {
            cout << "Projection::put_divu_in_node_rhs: not yet "
                 << "implemented for 3-d, non-zero divu\n";
            ParallelDescriptor::Abort("Bye");
        }
#endif
    }
}

//
// Put S in the rhs of the projector--cell based version.
//

void
Projection::put_divu_in_cc_rhs (MultiFab&       rhs,
                                int             level,
                                const BoxArray& grids,
                                Real            time)
{
    rhs.setVal(0.0);

    MultiFab* divu = getDivCond(level,1,time);

    for (MultiFabIterator mfi(rhs); mfi.isValid(false); ++mfi)
    {
        DependentMultiFabIterator dmfi(mfi, *divu);

        mfi().copy(dmfi());

#if (BL_SPACEDIM == 3)
        Real divumin = divu.min();
        Real divumax = divu.max();
        if (divumin != divumax || divumin != 0.0) 
        {
            cout << "Projection::put_divu_in_cc_rhs: not yet "
                 << "implemented for 3-d, non-zero divu\n";
            BoxLib::Abort("Bye");
        }
#endif
    }

    delete divu;
}

//
// Fill patch to get the divU.
//

void
Projection::getDivCond (int        level,
                        FArrayBox& fab,
                        int        ngrow,
                        Real       time)
{
    int Divu_Type, Divu;

    if (!LevelData[level].isStateVariable("divu", Divu_Type, Divu)) 
    {
	BoxLib::Error("Projection::getDivCond(): Divu not found");
    }

#ifndef NDEBUG
    fab.setVal(bogus_value);
#endif
    
    BoxLib::Error("Projection::getDivCond(FAB): not implemented");

//  LevelData[level].FillPatch(fab,0,time,Divu_Type,0,1);
}

MultiFab*
Projection::getDivCond (int  level,
                        int  ngrow,
                        Real time)
{
    int Divu_Type, Divu;

    if (!LevelData[level].isStateVariable("divu", Divu_Type, Divu)) 
    {
	BoxLib::Error("Projection::getDivCond(): Divu not found");
    }

    MultiFab* divu = new MultiFab(LevelData[level].boxArray(),1,ngrow);

#ifndef NDEBUG
    divu->setVal(bogus_value);
#endif

    FillPatchIterator fpi(LevelData[level],*divu,ngrow,0,time,Divu_Type,0,1);

    for ( ; fpi.isValid(); ++fpi)
    {
        (*divu)[fpi.index()].copy(fpi());
    }

    return divu;
}

// ==================================================
// Simple Update function follow
// ==================================================


// enforce periodicity on a multifab
void Projection::EnforcePeriodicity( MultiFab &psi, int nvar,
                                     const BoxArray &grids,
                                     const Geometry &geom )
{
    if (!geom.isAnyPeriodic())
        return;

    assert(nvar <= psi.nComp());

    Array<IntVect> pshifts(27);
    FArrayBox temp;
    const Box& domain = geom.Domain();
    
    // loop over all of the grids
    for(MultiFabIterator psimfi(psi); psimfi.isValid(); ++psimfi) 
    {

        // compute available periodic shifts
        Box dbox(psimfi().box());
        temp.resize(dbox,nvar);
        temp.copy(psimfi(),0,0,nvar);
        geom.periodicShift( domain, dbox, pshifts);
        
        // loop over the periodic shifts
        for (int iiv = 0; iiv < pshifts.length(); iiv++) 
        {
            IntVect iv = pshifts[iiv];

            // copy from psi to temp
            temp.shift(-iv);
            if(ParallelDescriptor::NProcs() > 1) 
            {
              cerr << "Projection::EnforcePeriodicity not implemented in parallel." << NL;
              cerr << "Nested MultiFab loops.\n";
              ParallelDescriptor::Abort("Exiting.");
            } 
            psi.copy(temp,0,0,nvar);

            // shift and copy from temp back to psi
            temp.shift(iv);
            psimfi().copy(temp,0,0,nvar);
            
        }
    }
}


// convert U from an Accleration like quantity to a velocity
// Unew = Uold + alpha*Unew
void Projection::UnConvertUnew( MultiFab &Uold, Real alpha, MultiFab &Unew, 
                                const BoxArray &grids )
{
    for(MultiFabIterator Uoldmfi(Uold); Uoldmfi.isValid(); ++Uoldmfi) 
    {
        DependentMultiFabIterator Unewmfi(Uoldmfi, Unew);
        assert(grids[Uoldmfi.index()] == Uoldmfi.validbox());
        UnConvertUnew( Uoldmfi(), alpha, Unewmfi(), Uoldmfi.validbox() );
    }
}



// convert U from an Accleration like quantity to a velocity
// Unew = Uold + alpha*Unew
void Projection::UnConvertUnew( FArrayBox &Uold, Real alpha, FArrayBox &Unew,
                                const Box &grd )
{
    assert( Unew.nComp() >= BL_SPACEDIM );
    assert( Uold.nComp() >= BL_SPACEDIM );
    assert( Unew.contains(grd) == true );
    assert( Uold.contains(grd) == true );
    
    const int* lo   = grd.loVect();
    const int* hi   = grd.hiVect();

    const int* uo_lo = Uold.loVect(); 
    const int* uo_hi = Uold.hiVect(); 
    const Real *uold = Uold.dataPtr(0);

    const int* un_lo = Unew.loVect(); 
    const int* un_hi = Unew.hiVect(); 
    const Real *unew = Unew.dataPtr(0);
    
    FORT_ACCEL_TO_VEL( lo, hi,
                       uold, ARLIM(uo_lo), ARLIM(uo_hi),
                       &alpha,
                       unew, ARLIM(un_lo), ARLIM(un_hi) );
}




// convert U to an Accleration like quantity
// Unew = (Unew - Uold)/alpha
void Projection::ConvertUnew( MultiFab &Unew,  MultiFab &Uold, Real alpha,
                              const BoxArray &grids )
{
    for(MultiFabIterator Uoldmfi(Uold); Uoldmfi.isValid(); ++Uoldmfi) 
    {
        DependentMultiFabIterator Unewmfi(Uoldmfi, Unew);
        assert(grids[Uoldmfi.index()] == Uoldmfi.validbox());
        ConvertUnew( Unewmfi(), Uoldmfi(), alpha, Uoldmfi.validbox() );
    }
}



// convert U to an Accleration like quantity
// Unew = (Unew - Uold)/alpha
void Projection::ConvertUnew( FArrayBox &Unew, FArrayBox &Uold, Real alpha,
                              const Box &grd )
{
    assert( Unew.nComp() >= BL_SPACEDIM );
    assert( Uold.nComp() >= BL_SPACEDIM );
    assert( Unew.contains(grd) == true );
    assert( Uold.contains(grd) == true );
    
    const int* lo   = grd.loVect();
    const int* hi   = grd.hiVect();

    const int* uo_lo = Uold.loVect(); 
    const int* uo_hi = Uold.hiVect(); 
    const Real *uold = Uold.dataPtr(0);

    const int* un_lo = Unew.loVect(); 
    const int* un_hi = Unew.hiVect(); 
    const Real *unew = Unew.dataPtr(0);
                    
    FORT_VEL_TO_ACCEL( lo, hi, 
                       unew, ARLIM(un_lo), ARLIM(un_hi),
                       uold, ARLIM(uo_lo), ARLIM(uo_hi),
                       &alpha );
}


// update a quantity U using the formula
// Unew = Unew + alpha*Uold
void Projection::UpdateArg1( MultiFab &Unew, Real alpha, MultiFab &Uold,
                             int nvar, const BoxArray &grids, int ngrow )
{
    for(MultiFabIterator Uoldmfi(Uold); Uoldmfi.isValid(); ++Uoldmfi) 
    {
        DependentMultiFabIterator Unewmfi(Uoldmfi, Unew);
        assert(grids[Uoldmfi.index()] == Uoldmfi.validbox());
        UpdateArg1( Unewmfi(), alpha, Uoldmfi(), nvar, Uoldmfi.validbox(), ngrow );
    }
}

// update a quantity U using the formula
// currently only the velocity, but will do the pressure as well
// Unew = Unew + alpha*Uold
void Projection::UpdateArg1( FArrayBox &Unew, Real alpha, FArrayBox &Uold,
                             int nvar, const Box &grd, int ngrow )
{
    // test for enough variables
    assert( nvar <= Uold.nComp() );
    assert( nvar <= Unew.nComp() );

    // compute the bounds
    Box b(grow(grd,ngrow));
    const Box &bb = Unew.box();
    if ( bb.ixType() == IndexType::TheNodeType() )
        b.surroundingNodes();
    assert( Uold.contains(b) == true );
    assert( Unew.contains(b) == true );

    // call the FORTRAN
    const int* lo    = b.loVect();
    const int* hi    = b.hiVect();
    const int* uo_lo = Uold.loVect(); 
    const int* uo_hi = Uold.hiVect(); 
    const Real *uold = Uold.dataPtr(0);

    const int* un_lo = Unew.loVect(); 
    const int* un_hi = Unew.hiVect(); 
    const Real *unew = Unew.dataPtr(0);
                    
    FORT_PROJ_UPDATE( lo, hi, &nvar, &ngrow,
                      unew, ARLIM(un_lo), ARLIM(un_hi),
                      &alpha,
                      uold, ARLIM(uo_lo), ARLIM(uo_hi) );
}


// add phi to P
void Projection::AddPhi( MultiFab &p, MultiFab &phi, const BoxArray &grids )
{
    for(MultiFabIterator pmfi(p); pmfi.isValid(); ++pmfi) 
    {
        DependentMultiFabIterator phimfi(pmfi, phi);
        pmfi().plus(phimfi());
    }
}


// convert phi into p^n+1/2
void Projection::incrPress(int level, Real dt)
{
    MultiFab &P_old = LevelData[level].get_old_data(Press_Type);
    MultiFab &P_new = LevelData[level].get_new_data(Press_Type);
    const BoxArray& grids = LevelData[level].boxArray();
    for(MultiFabIterator P_newmfi(P_new); P_newmfi.isValid(); ++P_newmfi) 
    {
        DependentMultiFabIterator P_oldmfi(P_newmfi, P_old);
        //assert(grids[P_newmfi.index()] == P_newmfi.validbox());
        //UpdateArg1( P_newmfi(), 1.0/dt, P_oldmfi(), 1, P_newmfi.validbox(), 1 );
        UpdateArg1( P_newmfi(), 1.0/dt, P_oldmfi(), 1, grids[P_newmfi.index()], 1 );
        P_oldmfi().setVal(bogus_value);
    }
}


// ==================================================
// projection scaling functions follow
// ==================================================


// This function scales variables at the start of a projection
void Projection::scaleVar( MultiFab *sig, int sig_nghosts,
                           MultiFab *vel,
                           const BoxArray &grids, int level )
{
    if ( sig != 0 )
        assert( sig->nComp() == 1 );
    if ( vel != 0 )
        assert( vel->nComp() >= BL_SPACEDIM );

    // convert sigma from rho to 1/rho
    // nghosts info needed to avoid divide by zero
    if ( sig != 0 )
        sig->invert(1.0,sig_nghosts);
    
    // scale by radius for RZ
    if ( CoordSys::IsRZ() ) 
    {
        if ( sig != 0 )
            radMult(level,*sig,0);
        if ( vel != 0 )
            for (int n = 0; n < BL_SPACEDIM; n++)
                radMult(level,*vel,n);
    }
    
    // scale level projection variables for a particular projection
    proj_scale_var(sig,vel,grids,level);
}


// This function rescales variables at the end of a projection
void Projection::rescaleVar( MultiFab *sig, int sig_nghosts,
                             MultiFab *vel,
                             const BoxArray &grids, int level )
{
    if ( sig != 0 )
        assert( sig->nComp() == 1 );
    if ( vel != 0 )
        assert( vel->nComp() >= BL_SPACEDIM );
    
    // divide by radius to rescale for RZ coordinates
    if ( CoordSys::IsRZ() ) 
    {
        if ( sig != 0 )
            radDiv(level,*sig,0);
        if ( vel != 0 )
            for (int n = 0; n < BL_SPACEDIM; n++)
                radDiv(level,*vel,n);
    }

    // convert sigma from 1/rho to rho
    // NOTE: this must come after division by r to be correct,
    // nghosts info needed to avoid divide by zero
    if ( sig != 0 )
        sig->invert(1.0,sig_nghosts);
    
    // unscale level projection variables for a particular projection
    proj_unscale_var(sig,vel,grids,level);
}




// multiply by a radius for r-z coordinates
void Projection::radMult(int level, MultiFab &mf, int comp)
{
  int i;
  int ngrow = mf.nGrow();
  int ncomp = mf.nComp();
  assert( comp >= 0 && comp < ncomp );
  int nr = radius_grow;
  int n[BL_SPACEDIM];
  for (i = 0; i < BL_SPACEDIM; i++) n[i] = ngrow;
  assert( radius_grow >= ngrow );
  //int ngrid = mf.length();
  for(MultiFabIterator mfmfi(mf); mfmfi.isValid(); ++mfmfi) 
  {
    FArrayBox& fab = mfmfi();
    assert(mf.box(mfmfi.index()) == mfmfi.validbox());
    const Box& bx = mfmfi.validbox();
    const int* lo = bx.loVect();
    const int* hi = bx.hiVect();
    Real* dat = fab.dataPtr(comp);
    Real* rad = &radius[level][mfmfi.index()];
    FORT_RADMPY(dat,ARLIM(lo),ARLIM(hi),&ngrow,rad,&nr,n);
  }
}


// divide by a radius for r-z coordinates
void Projection::radDiv(int level, MultiFab &mf, int comp)
{
  int ngrow = mf.nGrow();
  int ncomp = mf.nComp();
  assert( comp >= 0 && comp < ncomp );
  int nr = radius_grow;
  int n[BL_SPACEDIM];
  int i;
  for (i = 0; i < BL_SPACEDIM; i++) n[i] = ngrow;
  assert( radius_grow >= ngrow );
  //int ngrid = mf.length();
  for(MultiFabIterator mfmfi(mf); mfmfi.isValid(); ++mfmfi) 
  {
    FArrayBox& fab = mfmfi();
    assert(mf.box(mfmfi.index()) == mfmfi.validbox());
    const Box& bx = mfmfi.validbox();
    const int* lo = bx.loVect();
    const int* hi = bx.hiVect();
    Real* dat = fab.dataPtr(comp);
    Real* rad = &radius[level][mfmfi.index()];
    FORT_RADDIV(dat,ARLIM(lo),ARLIM(hi),&ngrow,rad,&nr,n);
  }
}

// ==================================================
// Outflow Boundary Conditions follow
// ==================================================


void Projection::set_level_projector_outflow_bcs(int       level,
						 Real      cur_state_time,
						 MultiFab& phi)
{
    // Do as set_initial_projection_outflow_bcs, except this one for single specified level.
    // In this case, since the finest level has not been advanced to the new time, the data
    // there is not yet available.  We contruct the dirichlet phi values for the outflow
    // boundary at the requested level, fill-patching in coarser data, which presumably has
    // been advanced up to or past new_time.

    // FIXME: the three outflow boundary routines should be collapsed into one, since they
    //        all do roughly the same thing
    
#if (BL_SPACEDIM == 2)
    const int rzflag = CoordSys::IsRZ();

    const Real* dx = parent->Geom(level).CellSize();
    const Box& domain = parent->Geom(level).Domain();
    const int outDir = 1;
    const int ccStripWidth = 3;
    const int ncStripWidth = 1;
    const Box state_strip = Box(adjCellHi(domain, outDir, ccStripWidth)).shift(outDir,-ccStripWidth+1);
    const Box phi_strip = surroundingNodes(bdryHi(domain, outDir, ncStripWidth));
    
    const int nGrow = 0;
    const int dstComp = 0;
    const int srcCompRho = Density, nCompRho = 1;
    const int srcCompVel = Xvel,    nCompVel = BL_SPACEDIM;
    const int srcCompDivu = 0,      nCompDivu = 1;
    const BoxArray state_strip_ba(&state_strip,1);
    MultiFab cc_MultiFab(state_strip_ba, 1, nGrow, Fab_noallocate);
    
    const int nCompPhi = 1;
    const BoxArray phi_strip_ba(&phi_strip,1);
    MultiFab phi_fine_strip(phi_strip_ba, nCompPhi, nGrow, Fab_allocate);
    phi_fine_strip.setVal(0.0);
    
    // Make r_i needed in HGPHIBC (set = 1 if cartesian)
    const Box region = Box(adjCellHi(domain, outDir, 1)).shift(outDir, -1);
    Array<Real> rcen(region.length(0), 1.0);
    if (CoordSys::IsRZ() == 1) 
    {
	parent->Geom(level).GetCellLoc(rcen, region, 0);
    }
    
    int Divu_Type, Divu;
    if (!LevelData[level].isStateVariable("divu", Divu_Type, Divu)) 
    {
	BoxLib::Error("Projection::set_level_projector_outflow_bcs(): Divu not found");
    }
    
    FillPatchIterator rhoFpi(LevelData[level], cc_MultiFab, nGrow, dstComp, cur_state_time,
			     State_Type, srcCompRho, nCompRho);
    
    FillPatchIterator velFpi(LevelData[level], cc_MultiFab, nGrow, dstComp, cur_state_time,
			     State_Type, srcCompVel, nCompVel);
    
    FillPatchIterator divuFpi(LevelData[level], cc_MultiFab, nGrow, dstComp, cur_state_time,
			      Divu_Type, srcCompDivu, nCompDivu);
    
    for ( ;rhoFpi.isValid() && velFpi.isValid() && divuFpi.isValid(); ++rhoFpi, ++velFpi, ++divuFpi)
    {
	DependentMultiFabIterator phimfi(rhoFpi, phi_fine_strip);
	
	// Fill phi_fine_strip with boundary cell values for phi, then copy into arg data
	// Note: Though this looks like a distributed operation, the MultiFab is built
	// on a single box...this is necesary currently, since FORT_HGPHIBC requires
	// a slice across the entire domain
	const int isPeriodicInX = (int) parent->Geom(level).isPeriodic(0);
	FORT_HGPHIBC(ARLIM(velFpi().loVect()),  ARLIM(velFpi().hiVect()),  velFpi().dataPtr(),
		     ARLIM(divuFpi().loVect()), ARLIM(divuFpi().hiVect()), divuFpi().dataPtr(),
		     ARLIM(rhoFpi().loVect()),  ARLIM(rhoFpi().hiVect()),  rhoFpi().dataPtr(),
		     ARLIM(region.loVect()),    ARLIM(region.hiVect()),    rcen.dataPtr(),   &dx[0],
		     ARLIM(phimfi().loVect()),  ARLIM(phimfi().hiVect()),  phimfi().dataPtr(),
		     &isPeriodicInX);
    }

    // Copy into arg data
    phi.copy(phi_fine_strip);
#endif
}

void Projection::set_initial_projection_outflow_bcs(MultiFab** vel,
						    MultiFab** sig,
						    MultiFab** phi, int c_lev, Real cur_divu_time)
{
#if (BL_SPACEDIM == 2)
// what we do
//  at the finest level
//  1) we define 3 cell-wide strips for
//     rho, S, U at the top of the domain
//     The strips are three cells wide for convenience
//  2) we define a 1 node-wide for phi at the top of the domain
//  3) we fill the rho, S, and U strips
//  4) we compute phi on the phi strip with HGPHIBC
//  5) we copy the phi strip into the phi multifab
//  we then "putdown" onto coarser levels

    // FIXME: vel and sig unused here.  We needed a filpatch operation, so pulled data from
    //        the state (presumably is where vel, sig from anyway).  sig is poorly named, and
    //        we need other stuff as well, so we should probably just remove sig and vel from args
    const int f_lev = finest_level;
    const int rzflag = CoordSys::IsRZ();

    if(c_lev!=0) 
    {
	cout << "initialVelocityProject: clev!=0--something wrong?\n";
	ParallelDescriptor::Abort("Exiting.");
    } 
    else 
    {
	// Get 3-wide cc box, state_strip, along top, incl. 2 rows of int and 1 row of ghosts.
	// Get 1-wide nc box, phi_strip, along top, grown out by one tangential to face
	const Real* dx = parent->Geom(f_lev).CellSize();
	const Box& domain = parent->Geom(f_lev).Domain();
	const int outDir = 1;
	const int ccStripWidth = 3;
	const int ncStripWidth = 1;
	const Box state_strip = Box(adjCellHi(domain, outDir, ccStripWidth)).shift(outDir,-ccStripWidth+1);
	const Box phi_strip = surroundingNodes(bdryHi(domain, outDir, ncStripWidth));

	const int nGrow = 0;
	const int dstComp = 0;
	const int srcCompRho = Density, nCompRho = 1;
	const int srcCompVel = Xvel,    nCompVel = BL_SPACEDIM;
	const int srcCompDivu = 0,      nCompDivu = 1;
	const BoxArray state_strip_ba(&state_strip,1);
	MultiFab cc_MultiFab(state_strip_ba, 1, nGrow, Fab_noallocate);
	
	const int nCompPhi = 1;
	const BoxArray phi_strip_ba(&phi_strip,1);
	MultiFab phi_fine_strip(phi_strip_ba, nCompPhi, nGrow, Fab_allocate);
	phi_fine_strip.setVal(0.0);

	// Make r_i needed in HGPHIBC (set = 1 if cartesian)
	const Box region = Box(adjCellHi(domain, outDir, 1)).shift(outDir, -1);
	Array<Real> rcen(region.length(0), 1.0);
	if (CoordSys::IsRZ() == 1) 
	{
	    parent->Geom(f_lev).GetCellLoc(rcen, region, 0);
	} 
    
	int Divu_Type, Divu;
        FArrayBox rho;
	if (!LevelData[c_lev].isStateVariable("divu", Divu_Type, Divu)) 
	{
	    BoxLib::Error("Projection::set_initial_projection_outflow_bcs(): Divu not found");
	}
    
	FillPatchIterator rhoFpi(LevelData[f_lev], cc_MultiFab, nGrow, dstComp, cur_divu_time,
				 State_Type, srcCompRho, nCompRho);
	
	FillPatchIterator velFpi(LevelData[f_lev], cc_MultiFab, nGrow, dstComp, cur_divu_time,
				 State_Type, srcCompVel, nCompVel);
	
	FillPatchIterator divuFpi(LevelData[f_lev], cc_MultiFab, nGrow, dstComp, cur_divu_time,
				  Divu_Type, srcCompDivu, nCompDivu);
	
	for ( ; rhoFpi.isValid()  &&  velFpi.isValid()  &&  divuFpi.isValid(); ++rhoFpi, ++velFpi, ++divuFpi)
	{
	    DependentMultiFabIterator phimfi(rhoFpi, phi_fine_strip);

	    // Note: potentially wasteful FillPatch on Rho, but this seemed cleanest
            rho.resize(rhoFpi.validbox(), nCompRho);
	    if (rho_wgt_vel_proj)
	    {
		rho.copy(rhoFpi());
	    }
            else
            {
		rho.setVal(1.0);
	    }
	    
	    // Fill phi_fine_strip with boundary cell values for phi, then copy into arg data
	    // Note: Though this looks like a distributed operation, the MultiFab is built
	    // on a single box...this is necesary currently, since FORT_HGPHIBC requires
	    // a slice across the entire domain
	    const int isPeriodicInX = (int) parent->Geom(f_lev).isPeriodic(0);
	    FORT_HGPHIBC(ARLIM(velFpi().loVect()),  ARLIM(velFpi().hiVect()),  velFpi().dataPtr(),
			 ARLIM(divuFpi().loVect()), ARLIM(divuFpi().hiVect()), divuFpi().dataPtr(),
			 ARLIM(rho.loVect()),       ARLIM(rho.hiVect()),       rho.dataPtr(),
			 ARLIM(region.loVect()),    ARLIM(region.hiVect()),    rcen.dataPtr(),    &dx[0],
			 ARLIM(phimfi().loVect()),  ARLIM(phimfi().hiVect()),  phimfi().dataPtr(),
			 &isPeriodicInX);
	}

	// Copy into arg data
	phi[f_lev]->copy(phi_fine_strip);

	IntVect ratio = IntVect::TheUnitVector();
	for (int lev=f_lev-1;lev>=c_lev;lev--) 
	{
	    ratio *= parent->refRatio(lev);
	    const Box& domainC = parent->Geom(lev).Domain();
	    const Box top_phiC_strip = surroundingNodes(bdryHi(domainC, outDir, ncStripWidth));
	    
	    const BoxArray top_phiC_strip_ba(&top_phiC_strip,1);
	    MultiFab phi_crse_strip(top_phiC_strip_ba, nCompPhi, nGrow, Fab_allocate);
	    phi_crse_strip.setVal(0.0);
	    
	    for (MultiFabIterator finemfi(phi_fine_strip); finemfi.isValid(); ++finemfi)
	    {
		DependentMultiFabIterator crsemfi(finemfi, phi_crse_strip);
		Box ovlp = Box(finemfi.validbox()).coarsen(ratio) & crsemfi.validbox();
		FORT_PUTDOWN (crsemfi().dataPtr(),ARLIM(crsemfi().loVect()),ARLIM(crsemfi().hiVect()),
			      finemfi().dataPtr(),ARLIM(finemfi().loVect()),ARLIM(finemfi().hiVect()),
			      ovlp.loVect(),ovlp.hiVect(),ratio.getVect());
	    }
	    phi[lev]->copy(phi_crse_strip);
	}
    }
#endif
}

void Projection::set_initial_syncproject_outflow_bcs(MultiFab** phi, 
						     int c_lev, Real start_time, Real dt)
{
#if (BL_SPACEDIM == 2)
// what we do is similar to what we do in set_initial_projection_outflow_bcs
// except that we work with time derivatives of U and S,
    
    const int f_lev = finest_level;
    const int rzflag = CoordSys::IsRZ();

    if(c_lev!=0) 
    {
	cout << "initialSyncProject: clev!=0--something wrong?\n";
	ParallelDescriptor::Abort("Exiting.");
    } 
    else 
    {
	const Real* dx = parent->Geom(f_lev).CellSize();
	const Box& domain = parent->Geom(f_lev).Domain();
	const int outDir = 1;
	const int ccStripWidth = 3;
	const int ncStripWidth = 1;
	const Box state_strip = Box(adjCellHi(domain, outDir, ccStripWidth)).shift(outDir,-ccStripWidth+1);
	const Box phi_strip = surroundingNodes(bdryHi(domain, outDir, ncStripWidth));

	const int nGrow = 0;
	const int dstComp = 0;
	const int srcCompRho = Density, nCompRho = 1;
	const int srcCompVel = Xvel,    nCompVel = BL_SPACEDIM;
	const int srcCompDivu = 0,      nCompDivu = 1;
	const BoxArray state_strip_ba(&state_strip,1);
	MultiFab cc_MultiFab(state_strip_ba, 1, nGrow, Fab_noallocate);
	
	const int nCompPhi = 1;
	const BoxArray phi_strip_ba(&phi_strip,1);
	MultiFab phi_fine_strip(phi_strip_ba, nCompPhi, nGrow, Fab_allocate);
	phi_fine_strip.setVal(0.0);

	// Make r_i needed in HGPHIBC (set = 1 if cartesian)
	const Box region = Box(adjCellHi(domain, outDir, 1)).shift(outDir, -1);
	Array<Real> rcen(region.length(0), 1.0);
	if (CoordSys::IsRZ() == 1) 
	{
	    parent->Geom(f_lev).GetCellLoc(rcen, region, 0);
	}
    
	int Divu_Type, Divu;
	if (!LevelData[c_lev].isStateVariable("divu", Divu_Type, Divu)) 
	{
	    BoxLib::Error("Projection::set_syncproject_outflow_bcs(): Divu not found");
	}

        FArrayBox rhonph;
        FArrayBox dudt;
        FArrayBox dsdt;
    
	FillPatchIterator rhoOldFpi(LevelData[f_lev], cc_MultiFab, nGrow, dstComp, start_time,
				    State_Type, srcCompRho, nCompRho);
	
	FillPatchIterator rhoNewFpi(LevelData[f_lev], cc_MultiFab, nGrow, dstComp, start_time + dt,
				    State_Type, srcCompRho, nCompRho);
	
	FillPatchIterator velOldFpi(LevelData[f_lev], cc_MultiFab, nGrow, dstComp, start_time,
				    State_Type, srcCompVel, nCompVel);
	
	FillPatchIterator velNewFpi(LevelData[f_lev], cc_MultiFab, nGrow, dstComp, start_time + dt,
				    State_Type, srcCompVel, nCompVel);
	
	FillPatchIterator divuOldFpi(LevelData[f_lev], cc_MultiFab, nGrow, dstComp, start_time,
				     Divu_Type, srcCompDivu, nCompDivu);
	
	FillPatchIterator divuNewFpi(LevelData[f_lev], cc_MultiFab, nGrow, dstComp, start_time + dt,
				     Divu_Type, srcCompDivu, nCompDivu);
	
	for ( ;rhoOldFpi.isValid()  &&  velOldFpi.isValid()  &&  divuOldFpi.isValid() &&
		  rhoNewFpi.isValid()  &&  velNewFpi.isValid()  &&  divuNewFpi.isValid();
	      ++rhoOldFpi, ++velOldFpi, ++divuOldFpi, ++rhoNewFpi, ++velNewFpi, ++divuNewFpi)
	{
	    DependentMultiFabIterator phimfi(rhoOldFpi, phi_fine_strip);

	    // Make rhonph, du/dt, and dsdt
	    assert(rhoOldFpi.validbox() == rhoNewFpi.validbox());
            rhonph.resize(rhoOldFpi.validbox(),nCompRho);
	    rhonph.copy(rhoNewFpi());
	    rhonph.plus(rhoOldFpi());
	    rhonph.divide(2.0);

	    assert(divuOldFpi.validbox() == divuNewFpi.validbox());
            dudt.resize(velOldFpi.validbox(),nCompVel);
	    dudt.copy(velNewFpi());
	    dudt.minus(velOldFpi());
	    dudt.divide(dt);
	    
	    assert(velOldFpi.validbox() == velNewFpi.validbox());
            dsdt.resize(divuOldFpi.validbox(),nCompDivu);
	    dsdt.copy(divuNewFpi());
	    dsdt.minus(divuOldFpi());
	    dsdt.divide(dt);
            //
	    // Fill phi_fine_strip with boundary cell values for phi, then copy into arg data
	    // Note: Though this looks like a distributed operation, the MultiFab is built
	    // on a single box...this is necesary currently, since FORT_HGPHIBC requires
	    // a slice across the entire domain
            //
	    const int isPeriodicInX = (int) parent->Geom(f_lev).isPeriodic(0);
	    FORT_HGPHIBC(ARLIM(dudt.loVect()),     ARLIM(dudt.hiVect()),     dudt.dataPtr(),
			 ARLIM(dsdt.loVect()),     ARLIM(dsdt.hiVect()),     dsdt.dataPtr(),
			 ARLIM(rhonph.loVect()),   ARLIM(rhonph.hiVect()),   rhonph.dataPtr(),
			 ARLIM(region.loVect()),   ARLIM(region.hiVect()),   rcen.dataPtr(),     &dx[0],
			 ARLIM(phimfi().loVect()), ARLIM(phimfi().hiVect()), phimfi().dataPtr(),
			 &isPeriodicInX);
	}

	// Copy into arg data
	phi[f_lev]->copy(phi_fine_strip);

	// Put down to coarser levels
	IntVect ratio = IntVect::TheUnitVector();
	for (int lev=f_lev-1;lev>=c_lev;lev--) 
	{
	    ratio *= parent->refRatio(lev);
	    const Box& domainC = parent->Geom(lev).Domain();
	    const Box top_phiC_strip = surroundingNodes(bdryHi(domainC, outDir, ncStripWidth));
	    
	    const BoxArray top_phiC_strip_ba(&top_phiC_strip,1);
	    MultiFab phi_crse_strip(top_phiC_strip_ba, nCompPhi, nGrow, Fab_allocate);
	    phi_crse_strip.setVal(0.0);
	    
	    for (MultiFabIterator finemfi(phi_fine_strip); finemfi.isValid(); ++finemfi)
	    {
		DependentMultiFabIterator crsemfi(finemfi, phi_crse_strip);
		Box ovlp = Box(finemfi.validbox()).coarsen(ratio) & crsemfi.validbox();
		FORT_PUTDOWN (crsemfi().dataPtr(),ARLIM(crsemfi().loVect()),ARLIM(crsemfi().hiVect()),
			      finemfi().dataPtr(),ARLIM(finemfi().loVect()),ARLIM(finemfi().hiVect()),
			      ovlp.loVect(),ovlp.hiVect(),ratio.getVect());
	    }
	    phi[lev]->copy(phi_crse_strip);
	}
    }
#endif
}

