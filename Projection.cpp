
//
// $Id: Projection.cpp,v 1.20 1997-10-08 20:15:43 car Exp $
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
Real bogus_value = 1.0e+20;

#define MAX_LEV 10

// initialization of static members
int       Projection::verbose = false;
int       Projection::P_code = 0;
Real      Projection::proj_tol = 1.0e-10;
Real      Projection::sync_tol = 1.0e-8;
Real      Projection::proj_abs_error = 1.0e-15;
Real      Projection::filter_factor = 0.0;
int       Projection::filter_u = 0;
int       Projection::rho_wgt_vel_proj = 0;
int       Projection::Divu_Type = -1;
int       Projection::do_outflow_bcs = 1;
int       Projection::make_sync_solvable = 0;
int Projection::do_cache = 1;

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

#if BL_SPACEDIM != 2
  if (CoordSys::IsRZ() == 1) amr_multigrid::SetRZ();
#endif
  setUpBcs();
  sync_proj = NULL;
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

  pp.query("Pcode",P_code);

  if(ParallelDescriptor::IOProcessor()) 
  {
    verbose = pp.contains("v");
  } 
  else 
  {
    verbose = false;
  }
  pp.get("proj_tol",proj_tol);
  pp.get("sync_tol",sync_tol);
  pp.query("make_sync_solvable",make_sync_solvable);
  pp.query("bogus_value",bogus_value);

  pp.query("filter_factor",filter_factor);
  pp.query("filter_u",filter_u);

  pp.query("proj_abs_error",proj_abs_error);

  pp.query("rho_wgt_vel_proj",rho_wgt_vel_proj);

  pp.query("do_outflow_bcs",do_outflow_bcs);

  pp.query("do_cache", do_cache);
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

  projector_bndry = new inviscid_fluid_boundary(proj_bc);
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

  if (sync_proj != NULL) 
  {
    delete sync_proj;
    sync_proj = NULL;
  }
}

// Build the aliasLib projection object
void
Projection::bldSyncProject()
{
  const Box& fdomain = parent->Geom(finest_level).Domain();

  // destruct if it already exists
  if (sync_proj != NULL) delete sync_proj;

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
					   *projector_bndry, P_code, do_cache?true:false);

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
  if (verbose) 
  {
    cout << "... level projector at level " << level << NL;
  }

  //----------------- manipulate state + pressure data ---------------------

  if (sync_proj == NULL) bldSyncProject();

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
  //Real hx = dx[0];

  const BoxArray& grids = LevelData[level].boxArray();
  const BoxArray& P_grids = P_old.boxArray();
  
  // compute (Ustar - Un)/dt as input to projector
  int i;
  if(level == 0) 
  {
    for(MultiFabIterator P_newmfi(P_new); P_newmfi.isValid(); ++P_newmfi) 
    {
      DependentMultiFabIterator P_oldmfi(P_newmfi, P_old);
      DependentMultiFabIterator U_newmfi(P_newmfi, U_new);
      DependentMultiFabIterator U_oldmfi(P_newmfi, U_old);
      // convert Unew to acceleration and Pnew to an update
      //assert(grids[P_newmfi.index()] == P_newmfi.validbox());
      //ConvertUnew( U_newmfi(), U_oldmfi(), dt, P_newmfi.validbox() );
      ConvertUnew( U_newmfi(), U_oldmfi(), dt, grids[P_newmfi.index()] );
      P_newmfi().minus( P_oldmfi() );
      P_newmfi().setVal(0.0);			    // FIXME?????
    }
  } 
  else 
  {   // level > 0
    for(MultiFabIterator P_newmfi(P_new); P_newmfi.isValid(); ++P_newmfi) 
    {
      DependentMultiFabIterator P_oldmfi(P_newmfi, P_old);
      DependentMultiFabIterator U_newmfi(P_newmfi, U_new);
      DependentMultiFabIterator U_oldmfi(P_newmfi, U_old);
      // interpolate values for P_new from coarse grid
      if(ParallelDescriptor::NProcs() > 1) 
      {
	cerr << "Projection::level_project not implemented in parallel.\n";
	ParallelDescriptor::Abort("Exiting.");
      } 
      else 
      {
	cerr << "Projection::level_project not implemented in parallel.\n";
      }

      LevelData[level].FillCoarsePatch(P_newmfi(),0,cur_pres_time,Press_Type,0,1);
      //assert(grids[P_newmfi.index()] == P_newmfi.validbox());
      assert(P_newmfi().box() == P_newmfi.fabbox());
      // convert Unew to acceleration and Pnew to an update
      ConvertUnew( U_newmfi(), U_oldmfi(), dt, grids[P_newmfi.index()] );
      P_newmfi().minus( P_oldmfi() );
      Box tempbox(P_newmfi().box());
      tempbox.grow(-2);
      P_newmfi().setVal(0.,tempbox,0,1);
    }
  }

  // set up outflow bcs

#if (BL_SPACEDIM == 2)
  int outflow_at_top = (phys_bc->lo(0) != Outflow &&
			phys_bc->lo(1) != Outflow &&
			phys_bc->hi(0) != Outflow &&
			phys_bc->hi(1) == Outflow);
  if (outflow_at_top && have_divu && do_outflow_bcs) 
  {
    set_level_projector_outflow_bcs(level,cur_pres_time,old_state_time,
				    cur_state_time,geom,
				    U_new,P_old,P_new,rho_half,dsdt);
  }
#endif

  // end set up outflow bcs

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
      crse_sync_reg->CrseDsdtAdd(dsdt, geom, rz_flag, bc, 
				 lowfix, hifix);
    } 
    if (level > 0) 
    { // increment sync registers between level and level-1
      Real invrat = 1.0/(double)crse_dt_ratio;
      const Geometry& crse_geom = parent->Geom(level-1);
      fine_sync_reg->FineDVAdd(U_new,dx,crse_geom,rz_flag,bc,invrat);
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

  Real local_proj_tol  = proj_tol/pow(10.0,level);
  Real local_abs_error = proj_abs_error/pow(10.0,level);
  if (!have_divu) 
  {
    // for divu=0 only
    sync_proj->project(u_real, p_real, null_amr_real, s_real, (Real*)dx,
                       local_proj_tol, level, level, local_abs_error);
  } 
  else 
  {
    // for general divu
    int use_u = 1;
    const int ngrids = grids.length();
    int nghost = 1; // required by aliaslib--rbp
    MultiFab rhs_cc(grids,1,nghost,Fab_allocate);
    rhs_cc.setVal(0.0);
    rhs_cc.copy(dsdt,0,0,1,nghost);
    radMult(level,rhs_cc,0);    
    for (i=0;i<ngrids;i++) 
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
			      local_proj_tol, level, level, local_abs_error);
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
  filterUandP( P_new, U_new, rho_half, grids, dx, dt );
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

  if (sync_proj == NULL) bldSyncProject();

  int rzflag = CoordSys::IsRZ();

  const BoxArray& grids = LevelData[level].boxArray();
  const BoxArray& P_grids = P_old.boxArray();
  MultiFab * rhs      = new MultiFab(P_grids,1,1,Fab_allocate);
  MultiFab * harm_phi = new MultiFab(P_grids,1,1,Fab_allocate);
  MultiFab * temp_phi = new MultiFab(P_grids,1,1,Fab_allocate);

  MultiFab * rho      = new MultiFab(grids,1,1,Fab_allocate);
  MultiFab * harm_vel = new MultiFab(grids,BL_SPACEDIM,1,Fab_allocate);

  rhs->setVal(0.);
  harm_phi->setVal(0.);
  harm_vel->setVal(0.);
  rho->setVal(1.);

  harm_phi->setBndry(bogus_value);

  Real prev_pres_time = cur_pres_time - dt;

  //int i;
  //for (i = 0; i < grids.length(); i++)
  for(MultiFabIterator temp_phimfi(*temp_phi); temp_phimfi.isValid(); ++temp_phimfi)
  {
    DependentMultiFabIterator harm_phimfi(temp_phimfi, *harm_phi);
    if(ParallelDescriptor::NProcs() > 1) 
    {
      ParallelDescriptor::Abort("Projection::harmonic_project not implemented in parallel");
    } 
    else 
    {
      cerr << "Projection::harmonic_project not implemented in parallel\n";
    }
    LevelData[level].FillCoarsePatch(temp_phimfi(),0,
				  prev_pres_time,Press_Type,0,1);
    LevelData[level].FillCoarsePatch(harm_phimfi(),0,
				  cur_pres_time,Press_Type,0,1);
    harm_phimfi().minus(temp_phimfi());
    Box tempbox(harm_phimfi().box());
    tempbox.grow(-2);
    harm_phimfi().setVal(0.,tempbox,0,1);
  }

  delete temp_phi;

  rho->setBndry(bogus_value);
  scaleVar(rho,1,NULL,grids,level);

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
  int use_u = 0;
  sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real, s_real,
			    use_u, (Real*)dx,
			    proj_tol, level, level, proj_abs_error);

  // copy and delete u_real
  // ----------------------------------------------------

  for (n = 0; n < BL_SPACEDIM; n++) 
  {
    // copy unnecessary since harm_vel will be discarded anyway
    delete u_real[n].remove(level);
  }

  //----------------- reset state + pressure data ---------------------

  // unscale variables for harmonic projection
  rescaleVar(rho,1,NULL,grids,level);

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

  if (sync_proj == NULL) bldSyncProject();

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
  int use_u = 1;
  sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real, s_real,
			    use_u, (Real*)dx,
			    sync_tol, c_lev, c_lev, proj_abs_error);

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
  if (sync_proj == NULL) bldSyncProject();

    
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

    restrict_level(u_real[n][c_lev], u_real[n][c_lev+1], ratio);
  }

  s_real.set(c_lev,   &rho_crse);
  s_real.set(c_lev+1, &rho_fine);

  restrict_level(s_real[c_lev], s_real[c_lev+1], ratio);

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
  int use_u            = 1;
  Real local_sync_tol  = sync_tol/pow(10.0,c_lev);
  Real local_abs_error = proj_abs_error/pow(10.0,c_lev);
  const Real* dx_fine  = parent->Geom(c_lev+1).CellSize();

  sync_proj->manual_project(u_real, p_real, null_amr_real,
			    crse_rhs_real, s_real,
			    use_u, (Real*)dx_fine,
			    local_sync_tol, c_lev, c_lev+1, local_abs_error);

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

  if (sync_proj == NULL) bldSyncProject();

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
    set_initial_projection_outflow_bcs(vel,sig,phi,parent,c_lev,cur_divu_time);
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
        if (Divu_Type == -1) 
	{
          int Divu;
          if (!LevelData[lev].isStateVariable("divu", Divu_Type, Divu)) 
	  {
            BoxLib::Error("Projection::getDivCond(): Divu not found");
          }
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
  Array<BoxArray>& full_mesh = sync_proj->mesh();
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
		       proj_tol, c_lev, f_lev, proj_abs_error);
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
      put_divu_in_cc_rhs(*rhslev,parent,lev,
			 grids,cur_divu_time);
      rhslev->mult(-1.0,0,1,nghost);
      radMult(lev,*rhslev,0); 
    }
    int use_u = 1;
    PArray<MultiFab> rhs_real(f_lev+1);
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
      rhs_real.set(lev, rhs_cc[lev]);
    }
    sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real, s_real,
			      use_u, (Real*)dx_lev,
			      proj_tol, c_lev, f_lev, proj_abs_error);
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

//------------------------------------------------------------
// INITIAL_SYNC_PROJECT
// the velocity projection in post_init, which computes the initial
// pressure used in the timestepping.
//------------------------------------------------------------

void Projection::initialSyncProject(int c_lev, MultiFab *sig[], Real dt, 
                                    Real strt_time, Real dt_init,
				    int have_divu)
{
  RunStats proj_stats("sync_project",c_lev);
  proj_stats.start();

  int lev;
  int f_lev = finest_level;
  if (verbose) 
  {
    cout << "SyncProject: levels = " << c_lev << "  " << f_lev << NL;
  }

  //----------------- manipulate state + pressure data ---------------------

  if (sync_proj == NULL) bldSyncProject();

  // gather data
  int rzflag = CoordSys::IsRZ();
  MultiFab *vel[MAX_LEV];
  MultiFab *phi[MAX_LEV];
  for (lev = c_lev; lev <= f_lev; lev++) 
  {
    vel[lev] = &LevelData[lev].get_new_data(State_Type);
    phi[lev] = &LevelData[lev].get_old_data(Press_Type);
  }

  MultiFab *rhs[MAX_LEV];
  if (have_divu) 
  {
    // set up rhs for manual project
#if (BL_SPACEDIM == 3)
    Real mindsdt = 0.0;
    Real maxdsdt = 0.0;
#endif
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
      AmrLevel& amr_level   = parent->getLevel(lev);

// make sure ghost cells are properly filled
      MultiFab &divu_new = amr_level.get_new_data(Divu_Type);
      divu_new.FillBoundary();
      MultiFab &divu_old = amr_level.get_new_data(Divu_Type);
      divu_old.FillBoundary();
      amr_level.setPhysBoundaryValues(Divu_Type,0,1,strt_time);
      amr_level.setPhysBoundaryValues(Divu_Type,0,1,strt_time+dt);

      const BoxArray& grids = amr_level.boxArray();
      const int ngrids = grids.length();
      int nghost = 1; 
      rhs[lev]  = new MultiFab(grids,1,nghost,Fab_allocate);
      MultiFab* rhslev = rhs[lev];
      rhslev->setVal(0.0);
      if(ParallelDescriptor::NProcs() > 1) 
      {
	cerr << "Projection::initialSyncProject not implemented in parallel." << NL;
	cerr << "This loop contains a call to FillPatch (in getDivCond).\n";
	ParallelDescriptor::Abort("Exiting.");
      } 
      else 
      {
	cerr << "Projection::initialSyncProject not implemented in parallel." << NL;
      }
      for (int i=0;i<ngrids;i++) 
      {
	Box divubox = grids[i];
	divubox.grow(1);
	FArrayBox divu(divubox,1);
	getDivCond(lev,divu,1,strt_time);
	FArrayBox dsdt(divubox,1);
	getDivCond(lev,dsdt,1,strt_time+dt);
	dsdt.minus(divu);
	dsdt.divide(dt);
        (*rhslev)[i].copy(dsdt);
#if (BL_SPACEDIM == 3)
	mindsdt = Min(mindsdt, dsdt.min());
	maxdsdt = Max(maxdsdt, dsdt.max());
#endif
      }
      rhslev->mult(-1.0,0,1);
      if ( CoordSys::IsRZ() ) 
        radMult(lev,*rhslev,0);    
    }
#if  (BL_SPACEDIM == 3)
    if(mindsdt!=maxdsdt || mindsdt!= 0.0) 
    {
      cout << "Projection::initialSyncProject: WARNING not yet " <<
	"implemented for 3-d, non-zero divu\n";
      ParallelDescriptor::Abort("Exiting.");
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
    set_initial_syncproject_outflow_bcs(phi,parent,c_lev,strt_time,dt);
  }
#endif

#ifdef SET_BOGUS_BNDRY
  // set velocity bndry values to bogus values
  for (lev = c_lev; lev <= f_lev; lev++) 
  {
    vel[lev]->setBndry(bogus_value,Xvel,BL_SPACEDIM);
    MultiFab &u_o = LevelData[lev].get_old_data(State_Type);
    u_o.setBndry(bogus_value,Xvel,BL_SPACEDIM);
    sig[lev]->setBndry(bogus_value);
  }
#endif

  // convert velocities to accelerations
  for (lev = c_lev; lev <= f_lev; lev++) 
  {
    LevelData[lev].setPhysBoundaryValues(State_Type,Xvel,BL_SPACEDIM,1);
    LevelData[lev].setPhysBoundaryValues(State_Type,Xvel,BL_SPACEDIM,0);
    MultiFab &u_o = LevelData[lev].get_old_data(State_Type);
    const BoxArray& grids = LevelData[lev].boxArray();
    ConvertUnew( *vel[lev], u_o, dt, grids );
  }

  // scale initial sync projection variables
  for (lev = c_lev; lev <= f_lev; lev++) 
  {
      AmrLevel& amr_level   = parent->getLevel(lev);
      const BoxArray& grids = amr_level.boxArray();
      scaleVar(sig[lev],1,vel[lev],grids,lev);
  }

  // build the aliaslib data structures
  //-------------------------------------------------------
  Array<BoxArray>& full_mesh = sync_proj->mesh();
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
    s_real.set(lev, sig[lev]);
  }

  for ( n = 0; n < BL_SPACEDIM; n++) 
  {
    for (lev = f_lev; lev >= c_lev+1; lev--) 
    {
      restrict_level(u_real[n][lev-1], u_real[n][lev],
		     parent->refRatio(lev-1));
    }
  }

  const Real* dx_lev = parent->Geom(f_lev).CellSize();

  
  // project
  //-------------------------------------------------------
  if (!have_divu) 
  {
    // zero divu only or debugging
    sync_proj->project(u_real, p_real, null_amr_real, s_real, (Real*)dx_lev,
		       proj_tol, c_lev, f_lev, proj_abs_error);
  } 
  else 
  {
    // general divu

    int use_u = 1;
    PArray<MultiFab> rhs_real(f_lev+1);
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
      rhs_real.set(lev, rhs[lev]);
    }
    sync_proj->manual_project(u_real, p_real, rhs_real, null_amr_real, s_real,
			      use_u, (Real*)dx_lev,
			      proj_tol, c_lev, f_lev, proj_abs_error);
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
  }

  //----------------- reset state + pressure data ---------------------

  // unscale initial sync projection variables
  for (lev = c_lev; lev <= f_lev; lev++) 
  {
      AmrLevel& amr_level   = parent->getLevel(lev);
      const BoxArray& grids = amr_level.boxArray();
      rescaleVar(sig[lev],1,vel[lev],grids,lev);
  }

  // add correction at coarse and fine levels
  for (lev = c_lev; lev <= f_lev; lev++) 
  {
    incrPress(lev, 1.0);
  }

  proj_stats.end();
}


// ==================================================
// Public nonprojector member functions follow
// ==================================================



// filter the velocities and pressures after a projection
void Projection::filterUandP( MultiFab &P_new,
                              MultiFab &U_new,
                              MultiFab *rho_half,
                              const BoxArray &grids,
                              const Real *dx, const Real dt )
{
    //int i;
    FArrayBox scratch;
    FArrayBox gradp;
    for(MultiFabIterator P_newmfi(P_new); P_newmfi.isValid(); ++P_newmfi) 
    {
	DependentMultiFabIterator U_newmfi(P_newmfi, U_new);
	DependentMultiFabIterator rho_halfmfi(P_newmfi, *rho_half);
        
	//assert(grids[P_newmfi.index()] == P_newmfi.validbox());
	Box gridbox(grids[P_newmfi.index()]);
        const int* lo = gridbox.loVect();
        const int* hi = gridbox.hiVect();
        
        const int* p_lo = P_newmfi().loVect(); 
        const int* p_hi = P_newmfi().hiVect(); 
        
        scratch.resize(P_newmfi().box(),1);
        scratch.setVal(0.0);
        
        FORT_FILTERP(P_newmfi().dataPtr(), scratch.dataPtr(), 
                     ARLIM(p_lo), ARLIM(p_hi),
                     lo,hi,dx,&filter_factor);
        
        if(filter_u) 
	{
            scratch.mult(dt);
            gradp.resize(gridbox,BL_SPACEDIM);
            const Real *gp_dat = gradp.dataPtr();
            FORT_GRADP(scratch.dataPtr(),ARLIM(p_lo),ARLIM(p_hi),
                       gp_dat,ARLIM(lo),ARLIM(hi),lo,hi,dx);
            for (int idim=0;idim<BL_SPACEDIM;idim++)
                gradp.divide(rho_halfmfi(),0,idim,1);
            U_newmfi().minus(gradp,Xvel,Xvel,BL_SPACEDIM);
        }
        
    }
}


// compute the node-centered divergence of U
void Projection::computeDV(MultiFab& DV, const MultiFab& U,
			   int src_comp, const Real* dx,
			   int is_rz)
{
  Real mult = 1.0;
    
  const BoxArray& U_boxes = U.boxArray();
  int ngrids = U_boxes.length();

  FArrayBox ufab;
  for(MultiFabIterator DVmfi(DV); DVmfi.isValid(); ++DVmfi) 
  {
    DependentMultiFabIterator Umfi(DVmfi, U);
    assert(U_boxes[Umfi.index()] == Umfi.validbox());
    Box ubox(grow(Umfi.validbox(),1));
    ufab.resize(ubox,BL_SPACEDIM);
    ufab.setVal(0.0);
    ufab.copy(Umfi(),ubox,src_comp,ubox,0,BL_SPACEDIM);

    Box ndbox(surroundingNodes(Umfi.validbox()));

    const int* ndlo = ndbox.loVect();
    const int* ndhi = ndbox.hiVect();
    const int* ulo = ubox.loVect();
    const int* uhi = ubox.hiVect();
    FORT_SRDIVU(ufab.dataPtr(),ARLIM(ulo),ARLIM(uhi),
		DVmfi().dataPtr(),ARLIM(ndlo),ARLIM(ndhi),
		ndlo,ndhi,dx,&mult,&is_rz);
  }
}

// put S in the rhs of the projector--node based version
void Projection::put_divu_in_node_rhs(MultiFab& rhs, Amr* parent, int level,
				 const int& nghost, const BoxArray& P_grids,
				 Real time, int user_rz)
{
  assert(user_rz>=-1 && user_rz<=1);
  rhs.setVal(0.0);
  const Geometry& geom = parent->Geom(level);
  const int ngrids = P_grids.length();
  int bcxlo = phys_bc->lo(0);
  int bcxhi = phys_bc->hi(0);
#if (BL_SPACEDIM == 2)
  int isrz;
  if (user_rz==-1) 
  {
    isrz = CoordSys::IsRZ();
  } 
  else 
  {
    isrz = user_rz;
  }
  int lowfix = (isrz==1 && bcxlo!=EXTRAP && bcxlo!=HOEXTRAP);
  int hifix = (isrz==1 && bcxhi!=EXTRAP && bcxhi!=HOEXTRAP);
  const Real* dx = geom.CellSize();
  Real hx = dx[0];
#endif
  const Box& domain = geom.Domain();
  const int* domlo = domain.loVect();
  const int* domhi = domain.hiVect();
  int i;
  for (i=0;i<ngrids;i++) 
  {
    Box divubox = P_grids[i];
    divubox.convert(IntVect::TheCellVector());
    divubox.grow(1);
    FArrayBox* divu = new FArrayBox(divubox,1);
    getDivCond(level,*divu,1,time);
    if(ParallelDescriptor::NProcs() > 1) 
    {
      cerr << "Projection::put_divu_in_node_rhs not implemented in parallel." << NL;
      cerr << "Nested MultiFab loops.\n";
      ParallelDescriptor::Abort("Exiting.");
    } 
    else 
    {
      cerr << "Projection::put_divu_in_node_rhs not implemented in parallel." << NL;
    }

#if (BL_SPACEDIM == 2)
    DEF_CLIMITS((*divu),divudat,divulo,divuhi);
    DEF_LIMITS(rhs[i],rhsdat,rhslo,rhshi);
    int rlen  = divubox.length(0);
    Array<Real> rcen;
    rcen.resize(rlen);
    if (isrz == 1) 
    {
      geom.GetCellLoc(rcen,divubox, 0);
    } 
    else 
    {
      int ii;
      for (ii=0; ii<rlen; ii++) rcen[ii]=1.0;
    }
    int extrap_edges = 0;
    int extrap_corners = 1;
    FORT_HGC2N(&nghost,ARLIM(divulo),ARLIM(divuhi),divudat,
	       rcen.dataPtr(),
	       ARLIM(rhslo),ARLIM(rhshi),rhsdat,
	       domlo,domhi,lowfix,hifix,&hx,
	       &extrap_edges, &extrap_corners,&isrz);
    delete divu;
#endif
#if (BL_SPACEDIM == 3)
    Real divumin=divu->min();
    Real divumax=divu->max();
    if(divumin!=divumax || divumin!= 0.0) 
    {
      cout << "Projection::put_divu_in_node_rhs: not yet " <<
	"implemented for 3-d, non-zero divu\n";
      ParallelDescriptor::Abort("Exiting.");
    }
    delete divu;
#endif
  }
}

// put S in the rhs of the projector--cell based version
void Projection::put_divu_in_cc_rhs(MultiFab& rhs, Amr* parent, int level,
				 const BoxArray& grids, Real time)
{
  rhs.setVal(0.0);
  const int ngrids = grids.length();
  int i;
  for (i=0;i<ngrids;i++) 
  {
    Box divubox = grids[i];
    divubox.grow(1);
    FArrayBox divu(divubox,1);
    if(ParallelDescriptor::NProcs() > 1) 
    {
      cerr << "Projection::put_divu_in_cc_rhs not implemented in parallel" << NL;
      cerr << "Nested MultiFab loops.\n";
      ParallelDescriptor::Abort("Exiting.");
    } 
    else 
    {
      cerr << "Projection::put_divu_in_cc_rhs not implemented in parallel" << NL;
    }
    getDivCond(level,divu,1,time);
    rhs[i].copy(divu);
#if (BL_SPACEDIM == 3)
    Real divumin=divu.min();
    Real divumax=divu.max();
    if(divumin!=divumax || divumin!= 0.0) 
    {
      cout << "Projection::put_divu_in_cc_rhs: not yet " <<
	"implemented for 3-d, non-zero divu\n";
      exit(0);
    }
#endif
  }
}


// fill patch to get the divU
void Projection::getDivCond(int level, FArrayBox& fab, int ngrow, Real time)
{
  if (Divu_Type == -1) 
  {
    int Divu;
    if (!LevelData[level].isStateVariable("divu", Divu_Type, Divu)) 
    {
      BoxLib::Error("Projection::getDivCond(): Divu not found");
    }
  }
  fab.setVal(1.0e30); // for debugging only
  LevelData[level].FillPatch(fab,0,time,Divu_Type,0,1);
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

    assert ( nvar <= psi.nComp() );

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
	    else 
	    {
              cerr << "Projection::EnforcePeriodicity not implemented in parallel." << NL;
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
    if ( sig != NULL )
        assert( sig->nComp() == 1 );
    if ( vel != NULL )
        assert( vel->nComp() >= BL_SPACEDIM );

    // convert sigma from rho to 1/rho
    // nghosts info needed to avoid divide by zero
    if ( sig != NULL )
        sig->invert(1.0,sig_nghosts);
    
    // scale by radius for RZ
    if ( CoordSys::IsRZ() ) 
    {
        if ( sig != NULL )
            radMult(level,*sig,0);
        if ( vel != NULL )
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
    if ( sig != NULL )
        assert( sig->nComp() == 1 );
    if ( vel != NULL )
        assert( vel->nComp() >= BL_SPACEDIM );
    
    // divide by radius to rescale for RZ coordinates
    if ( CoordSys::IsRZ() ) 
    {
        if ( sig != NULL )
            radDiv(level,*sig,0);
        if ( vel != NULL )
            for (int n = 0; n < BL_SPACEDIM; n++)
                radDiv(level,*vel,n);
    }

    // convert sigma from 1/rho to rho
    // nghosts info needed to avoid divide by zero
    if ( sig != NULL )
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


void Projection::set_level_projector_outflow_bcs(int level,
				 Real cur_pres_time,
                                 Real old_time, Real new_time,
                                 const Geometry& geom,
                                 MultiFab &U_new,
                                 MultiFab &P_old,
                                 MultiFab &P_new,
                                 MultiFab * rho_half,
                                 MultiFab& dsdt)
{
#if (BL_SPACEDIM == 2)
  int rzflag = CoordSys::IsRZ();
  const Real* dx = geom.CellSize();
  Real hx = dx[0];

  const BoxArray& grids   = LevelData[level].boxArray();
  const BoxArray& P_grids = P_old.boxArray();

  // compute (Ustar - Un)/dt as input to projector

  if(level==0 && grids.length()==1) 
  {

    // FOLLOWING FOR SINGLE GRID, SINGLE LEVEL ONLY

    Box rbox = P_grids[0];
    rbox.convert(IntVect::TheCellVector());
    int rlen  = rbox.length(0);
    Array<Real> rcen;
    rcen.resize(rlen);
    if (CoordSys::IsRZ() == 1) 
    {
      geom.GetCellLoc(rcen,rbox, 0);
    } 
    else 
    {
      for (int i=0; i<rlen; i++) rcen[i]=1.0;
    }
    const int* rlo = rbox.loVect();
    const int* rhi = rbox.hiVect();
    DEF_CLIMITS(dsdt[0],sdat,slo,shi);
    DEF_CLIMITS(U_new[0],udat,ulo,uhi);
    DEF_CLIMITS((*rho_half)[0],rhodat,rholo,rhohi);
    DEF_LIMITS(P_new[0],pdat,p_lo,p_hi);
    FORT_HGPHIBC(ARLIM(ulo),ARLIM(uhi),udat,ARLIM(slo),ARLIM(shi),sdat,
		 ARLIM(rholo),ARLIM(rhohi),rhodat,
		 ARLIM(rlo),ARLIM(rhi),rcen.dataPtr(),&hx,
		 ARLIM(p_lo),ARLIM(p_hi),pdat);

  } 
  else if (level==0  && grids.length()!=1) 
  {
// what we do
//  1) we define 3 cell-wide strips for
//     rho, S, U at the top of the domain
//     The strips are three cells wide for convenience
//  2) we define a 2 cell-wide strips for
//     phi at the top of the domain
//     The strip is two cells wide for convenience.
//  3) we fill the rho, S, and U strips
//  4) we compute phi on the phi strip with HGPHIBC
//  5) we copy the phi strip into the phi multifab
    const Real* dx = parent->Geom(0).CellSize();
    Real hx = dx[0];
    const Box& domain = parent->Geom(0).Domain();
    IntVect bigend = domain.bigEnd();
    IntVect smallend = domain.smallEnd();
    int jhi = domain.bigEnd(1);
    smallend.setVal(1,jhi);
    Box top_strip(smallend,bigend,IntVect::TheCellVector());
    Box top_phi_strip = top_strip;
    Box rbox = top_strip;

    bigend.setVal(1,jhi+1);
    smallend.setVal(1,jhi-1);

    Box state_strip(smallend,bigend,IntVect::TheCellVector());

    FArrayBox rho_strip(state_strip,1);
    FArrayBox rho0_strip(state_strip,1);
    FArrayBox divu_strip(state_strip,1);
    FArrayBox vel_strip(state_strip,2);

    top_phi_strip.convert(IntVect::TheNodeVector());
    top_phi_strip.growLo(1,-1);
    top_phi_strip.growHi(1,1);
    top_phi_strip.grow(0,1);
    FArrayBox phi_strip(top_phi_strip,1);
    phi_strip.setVal(0.0);
    const BoxArray& P_grids = P_old.boxArray();
    const int ngrids = P_grids.length();
    int i;
    for(MultiFabIterator U_newmfi(U_new); U_newmfi.isValid(); ++U_newmfi) 
    {
      DependentMultiFabIterator dsdtmfi(U_newmfi, dsdt);
      assert(P_grids[U_newmfi.index()] == U_newmfi.validbox());
      Box destbox = U_newmfi.validbox();
      destbox.convert(IntVect::TheCellVector());
      destbox &= state_strip;
      Box srcbox = destbox;
      if(destbox.ok()) 
      {
	vel_strip.copy(U_newmfi(),srcbox,0,destbox,0,2);
	divu_strip.copy(dsdtmfi(),srcbox,0,destbox,0,1);
      }
    }

    LevelData[0].FillPatch(rho0_strip,0,old_time,State_Type,Density,1);
    LevelData[0].FillPatch(rho_strip,0,new_time,State_Type,Density,1);
    rho_strip.plus(rho0_strip);
    rho_strip.divide(2.0);
        
    int rlen  = rbox.length(0);
    Array<Real> rcen;
    rcen.resize(rlen);
    if (CoordSys::IsRZ() == 1) 
    {
      parent->Geom(0).GetCellLoc(rcen,rbox, 0);
    } 
    else 
    {
      for (i=0; i<rlen; i++) rcen[i]=1.0;
    }
    const int* rlo = rbox.loVect();
    const int* rhi = rbox.hiVect();
    DEF_CLIMITS(divu_strip,sdat,slo,shi);
    DEF_CLIMITS(vel_strip,udat,ulo,uhi);
    DEF_CLIMITS(rho_strip,rhodat,rholo,rhohi);
    DEF_LIMITS(phi_strip,pdat,p_lo,p_hi);
    FORT_HGPHIBC(ARLIM(ulo),ARLIM(uhi),udat,ARLIM(slo),ARLIM(shi),sdat,
		 ARLIM(rholo),ARLIM(rhohi),rhodat,
		 ARLIM(rlo),ARLIM(rhi),rcen.dataPtr(),&hx,
		 ARLIM(p_lo),ARLIM(p_hi),pdat);
    for(MultiFabIterator P_newmfi(P_new); P_newmfi.isValid(); ++P_newmfi) 
    {
      P_newmfi().copy(phi_strip);
    }
  } 
  else 
  {
// do nothing for now
  }
#endif
}

void Projection::set_initial_projection_outflow_bcs(MultiFab** vel,
   MultiFab** sig, MultiFab** phi, Amr* parent, int c_lev, Real cur_divu_time)
{
#if (BL_SPACEDIM == 2)
  int rho_comp = 0;
  if (rho_wgt_vel_proj) rho_comp = Density;

  int lev;
  int f_lev = finest_level;
  int rzflag = CoordSys::IsRZ();

  const BoxArray& P_grids = phi[0]->boxArray();
  const int ngrids = P_grids.length();
  if(c_lev!=0) 
  {
    cout << "initialVelocityProject: clev!=0--something wrong?\n";
    ParallelDescriptor::Abort("Exiting.");
  } 
  else if(c_lev==f_lev && ngrids==1) 
  {
    // FOLLOWING FOR SINGLE GRID, SINGLE LEVEL ONLY
    lev = c_lev;
    const Real* dx = parent->Geom(lev).CellSize();
    Real hx = dx[0];
    const BoxArray& P_grids = phi[lev]->boxArray();
    const int ngrids = P_grids.length();
    MultiFab* philev = phi[lev];
    int i=0;
    Box divubox = P_grids[i];
    divubox.convert(IntVect::TheCellVector());
    divubox.grow(1);
    FArrayBox* divu = new FArrayBox(divubox,1);
    getDivCond(lev,*divu,1,cur_divu_time);

    Box rbox = P_grids[i];
    rbox.convert(IntVect::TheCellVector());
    int rlen  = rbox.length(0);
    Array<Real> rcen;
    rcen.resize(rlen);
    if (CoordSys::IsRZ() == 1) 
    {
      parent->Geom(lev).GetCellLoc(rcen,rbox, 0);
    } 
    else 
    {
      int ii;
      for (ii=0; ii<rlen; ii++) rcen[ii]=1.0;
    }

// what we do
//  at the finest level
//  1) we define 3 cell-wide strips for
//     rho, S, U at the top of the domain
//     The strips are three cells wide for convenience
//  2) we define a 2 cell-wide strips for
//     phi at the top of the domain
//     The strip is two cells wide for convenience.
//  3) we fill the rho, S, and U strips
//  4) we compute phi on the phi strip with HGPHIBC
//  5) we copy the phi strip into the phi multifab
//  we then "putdown" onto coarser levels

    const int* rlo = rbox.loVect();
    const int* rhi = rbox.hiVect();
    DEF_CLIMITS((*divu),sdat,slo,shi);
    DEF_CLIMITS((((*vel)[lev])[i]),udat,ulo,uhi);
    const int* rholo = (((*sig)[lev])[i]).loVect();
    const int* rhohi = (((*sig)[lev])[i]).hiVect();
    const Real* rhodat = (((*sig)[lev])[i]).dataPtr(rho_comp);
    DEF_LIMITS(((*philev)[i]),pdat,p_lo,p_hi);

    FORT_HGPHIBC(ARLIM(ulo),ARLIM(uhi),udat,ARLIM(slo),ARLIM(shi),sdat,
		 ARLIM(rholo),ARLIM(rhohi),rhodat,
		 ARLIM(rlo),ARLIM(rhi),rcen.dataPtr(),&hx,
		 ARLIM(p_lo),ARLIM(p_hi),pdat);
    delete divu;
  } 
  else 
  {
    FArrayBox *phi_crse_strip;
    FArrayBox *phi_fine_strip;
    const Real* dx = parent->Geom(f_lev).CellSize();
    Real hx = dx[0];
    const Box& domain = parent->Geom(f_lev).Domain();
    IntVect bigend = domain.bigEnd();
    IntVect smallend = domain.smallEnd();
    int jhi = domain.bigEnd(1);
    smallend.setVal(1,jhi);
    Box top_strip(smallend,bigend,IntVect::TheCellVector());
    Box top_phi_strip = top_strip;
    Box rbox = top_strip;

    bigend.setVal(1,jhi+1);
    smallend.setVal(1,jhi-1);

    Box state_strip(smallend,bigend,IntVect::TheCellVector());

    FArrayBox rho_strip(state_strip,1);
    FArrayBox divu_strip(state_strip,1);
    FArrayBox vel_strip(state_strip,2);

    top_phi_strip.convert(IntVect::TheNodeVector());
    top_phi_strip.growLo(1,-1);
    top_phi_strip.growHi(1,1);
    top_phi_strip.grow(0,1);
    phi_fine_strip = new FArrayBox(top_phi_strip,1);
    phi_fine_strip->setVal(0.0);

    getDivCond(f_lev,divu_strip,0,cur_divu_time);
    LevelData[f_lev].FillPatch(vel_strip,0,cur_divu_time,State_Type,0,2);
    if (rho_wgt_vel_proj) 
    {
      LevelData[f_lev].FillPatch(rho_strip,0,cur_divu_time,
                                State_Type,Density,1);
    } 
    else 
    {
      rho_strip.setVal(1.0);
    }

    int rlen  = rbox.length(0);
    Array<Real> rcen;
    rcen.resize(rlen);
    if (CoordSys::IsRZ() == 1) 
    {
      parent->Geom(f_lev).GetCellLoc(rcen,rbox, 0);
    } 
    else 
    {
      for (int i=0; i<rlen; i++) rcen[i]=1.0;
    }
    const int* rlo = rbox.loVect();
    const int* rhi = rbox.hiVect();

    DEF_CLIMITS(divu_strip,sdat,slo,shi);
    DEF_CLIMITS(vel_strip,udat,ulo,uhi);
    DEF_CLIMITS(rho_strip,rhodat,rholo,rhohi);
    DEF_LIMITS((*phi_fine_strip),pdat,p_lo,p_hi);

    FORT_HGPHIBC(ARLIM(ulo),ARLIM(uhi),udat,ARLIM(slo),ARLIM(shi),sdat,
		 ARLIM(rholo),ARLIM(rhohi),rhodat,
		 ARLIM(rlo),ARLIM(rhi),rcen.dataPtr(),&hx,
		 ARLIM(p_lo),ARLIM(p_hi),pdat);

    const BoxArray& P_grids = phi[f_lev]->boxArray();
    const int ngrids = P_grids.length();
    MultiFab* philev = phi[f_lev];
    for(MultiFabIterator philevmfi(*philev); philevmfi.isValid(); ++philevmfi) 
    {
      philevmfi().copy((*phi_fine_strip));
    }
    if (c_lev!=f_lev) 
    {
      for (lev=f_lev-1;lev>=c_lev;lev--) 
      {
	IntVect ratio = parent->refRatio(lev);
	const Box& domain = parent->Geom(lev).Domain();
	IntVect bigend = domain.bigEnd();
	IntVect smallend = domain.smallEnd();
	int jhi = domain.bigEnd(1);
	smallend.setVal(1,jhi);
	Box top_strip(smallend,bigend,IntVect::TheCellVector());
	Box top_phi_strip = top_strip;
	top_phi_strip.convert(IntVect::TheNodeVector());
	top_phi_strip.growLo(1,-1);
	top_phi_strip.grow(0,1);
	phi_crse_strip = new FArrayBox(top_phi_strip,1);
	phi_crse_strip->setVal(0.0);
	//int ng_pres = 0;
	DEF_LIMITS((*phi_crse_strip),pcrse,clo,chi);
	DEF_CLIMITS((*phi_fine_strip),pfine,flo,fhi);
	int ovlo[2], ovhi[2];
	ovlo[0] = clo[0]+1;
	ovhi[0] = chi[0]-1;
	ovlo[1] = clo[1];
	ovhi[1] = chi[1];
	FORT_PUTDOWN (pcrse,ARLIM(clo),ARLIM(chi),
		      pfine,ARLIM(flo),ARLIM(fhi),
		      ovlo,ovhi,ratio.getVect());
	const BoxArray& P_grids = phi[lev]->boxArray();
	const int ncrsegrids = P_grids.length();
	MultiFab* philev = phi[lev];
	//for (int i=0;i<ncrsegrids;i++)
        for(MultiFabIterator philevmfi(*philev); philevmfi.isValid(); ++philevmfi) 
	{
	  philevmfi().copy(*phi_crse_strip);
	}
	delete phi_fine_strip;
	phi_fine_strip = phi_crse_strip;
      }
      delete phi_crse_strip;
    } 
    else 
    {
      delete phi_fine_strip;
    }
  } 
#endif
}

void Projection::set_initial_syncproject_outflow_bcs(MultiFab** phi, 
                 Amr* parent, int c_lev, Real start_time, Real dt)
{
#if (BL_SPACEDIM == 2)
  //int rho_comp = Density;
  int lev;
  int f_lev = finest_level;
  int rzflag = CoordSys::IsRZ();

  const BoxArray& P_grids = phi[0]->boxArray();
  const int ngrids = P_grids.length();
  if(c_lev!=0) 
  {
    cout << "initialSyncProject: clev!=0--something wrong?\n";
    ParallelDescriptor::Abort("Exiting.");
  } 
  else 
  {
// what we do is similar to what we do in set_initial_projection_outflow_bcs
// except that we work with derivatives of U and S,
    FArrayBox *phi_crse_strip;
    FArrayBox *phi_fine_strip;
    const Real* dx = parent->Geom(f_lev).CellSize();
    Real hx = dx[0];
    const Box& domain = parent->Geom(f_lev).Domain();
    IntVect bigend = domain.bigEnd();
    IntVect smallend = domain.smallEnd();
    int jhi = domain.bigEnd(1);
    smallend.setVal(1,jhi);
    Box top_strip(smallend,bigend,IntVect::TheCellVector());
    FArrayBox rho_strip(top_strip,1);
    FArrayBox rho0_strip(top_strip,1);
    FArrayBox divu_strip(top_strip,1);
    FArrayBox dsdt_strip(top_strip,1);
    FArrayBox vel_strip(top_strip,2);
    FArrayBox vel0_strip(top_strip,2);
    Box top_phi_strip = top_strip;
    top_phi_strip.convert(IntVect::TheNodeVector());
    top_phi_strip.growLo(1,-1);
    top_phi_strip.growHi(1,1);
    top_phi_strip.grow(0,1);
    phi_fine_strip = new FArrayBox(top_phi_strip,1);
    phi_fine_strip->setVal(0.0);
    getDivCond(f_lev,divu_strip,0,start_time);
    getDivCond(f_lev,dsdt_strip,0,start_time+dt);
    dsdt_strip.minus(divu_strip);
    dsdt_strip.divide(dt);
    LevelData[f_lev].FillPatch(vel0_strip,0,start_time,State_Type,0,2);
    LevelData[f_lev].FillPatch(vel_strip,0,start_time+dt,State_Type,0,2);
    vel_strip.minus(vel0_strip);
    vel_strip.divide(dt);
    LevelData[f_lev].FillPatch(rho0_strip,0,start_time,State_Type,Density,1);
    LevelData[f_lev].FillPatch(rho_strip,0,start_time+dt,State_Type,Density,1);
    rho_strip.plus(rho0_strip);
    rho_strip.divide(2.0);
    Box rbox = top_strip;
    int rlen  = rbox.length(0);
    Array<Real> rcen;
    rcen.resize(rlen);
    int i;
    if (CoordSys::IsRZ() == 1) 
    {
      parent->Geom(f_lev).GetCellLoc(rcen,rbox, 0);
    } 
    else 
    {
      for (i=0; i<rlen; i++) rcen[i]=1.0;
    }
    const int* rlo = rbox.loVect();
    const int* rhi = rbox.hiVect();
    DEF_CLIMITS(divu_strip,sdat,slo,shi);
    DEF_CLIMITS(vel_strip,udat,ulo,uhi);
    DEF_CLIMITS(rho_strip,rhodat,rholo,rhohi);
    DEF_LIMITS((*phi_fine_strip),pdat,p_lo,p_hi);

    FORT_HGPHIBC(ARLIM(ulo),ARLIM(uhi),udat,ARLIM(slo),ARLIM(shi),sdat,
		 ARLIM(rholo),ARLIM(rhohi),rhodat,
		 ARLIM(rlo),ARLIM(rhi),rcen.dataPtr(),&hx,
		 ARLIM(p_lo),ARLIM(p_hi),pdat);

    //const BoxArray& P_grids = phi[f_lev]->boxArray();
    //const int ngrids = P_grids.length();
    MultiFab* philev = phi[f_lev];
    for(MultiFabIterator philevmfi(*philev); philevmfi.isValid(); ++philevmfi) 
    {
      philevmfi().copy((*phi_fine_strip));
    } 
    if (c_lev!=f_lev) 
    {
      for (lev=f_lev-1;lev>=c_lev;lev--) 
      {
	IntVect ratio = parent->refRatio(lev);
	const Box& domain = parent->Geom(lev).Domain();
	IntVect bigend = domain.bigEnd();
	IntVect smallend = domain.smallEnd();
	int jhi = domain.bigEnd(1);
	smallend.setVal(1,jhi);
	Box top_strip(smallend,bigend,IntVect::TheCellVector());
	Box top_phi_strip = top_strip;
	top_phi_strip.convert(IntVect::TheNodeVector());
	top_phi_strip.growLo(1,-1);
	top_phi_strip.grow(0,1);
	phi_crse_strip = new FArrayBox(top_phi_strip,1);
	phi_crse_strip->setVal(0.0);
	//int ng_pres = 0;
	DEF_LIMITS((*phi_crse_strip),pcrse,clo,chi);
	DEF_CLIMITS((*phi_fine_strip),pfine,flo,fhi);
	int ovlo[2], ovhi[2];
	ovlo[0] = clo[0]+1;
	ovhi[0] = chi[0]-1;
	ovlo[1] = clo[1];
	ovhi[1] = chi[1];
	FORT_PUTDOWN (pcrse,ARLIM(clo),ARLIM(chi),
		      pfine,ARLIM(flo),ARLIM(fhi),
		      ovlo,ovhi,ratio.getVect());
	//const BoxArray& P_grids = phi[lev]->boxArray();
	//const int ncrsegrids = P_grids.length();
	MultiFab* philev = phi[lev];
	//int ii;
	//for (ii=0;ii<ncrsegrids;ii++) {
        for(MultiFabIterator philevmfi(*philev); philevmfi.isValid(); ++philevmfi) 
	{
	  philevmfi().copy(*phi_crse_strip);
	}
	delete phi_fine_strip;
	phi_fine_strip = phi_crse_strip;
      }
      delete phi_crse_strip;
    } 
    else 
    {
      delete phi_fine_strip;
    }
  } 
#endif
}

