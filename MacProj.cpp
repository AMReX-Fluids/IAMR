// $Id: MacProj.cpp,v 1.1 1997-07-08 23:08:12 vince Exp $

#include <stdio.h>
#include <Misc.H>
#include <LO_BCTYPES.H>
#include <CoordSys.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <BoxDomain.H>
#include <MacProj.H>
#include <RunStats.H>
#include <MacBndry.H>
#include <MacOpMacDrivers.H>
#include <NavierStokes.H>

#include <MACPROJ_F.H>

#ifndef _NavierStokes_H_
enum StateType {State_Type=0, Press_Type};
#  if (BL_SPACEDIM == 2)
enum StateNames  { Xvel=0, Yvel, Density};
#  else
enum StateNames  { Xvel=0, Yvel, Zvel, Density};
#  endif
#endif

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
REAL* fabdat = (fab).dataPtr();
#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const REAL* fabdat = (fab).dataPtr();

#define GEOM_GROW 1
#define HYP_GROW 3

// initialization of static members
//#ifdef ATMOSPHERE
//int  MacProj::use_viscosity    = 0;
//#else
//int  MacProj::use_viscosity    = 1;
//#endif

int  MacProj::use_viscosity    = 1;
int  MacProj::verbose          = false;
int  MacProj::use_cg_solve     = 0;
REAL MacProj::mac_tol          = 1.0e-8;
REAL MacProj::mac_abs_tol      = 1.0e-15;
REAL MacProj::mac_sync_tol     = 1.0e-8;
int  MacProj::do_outflow_bcs   = 1;
int  MacProj::fix_mac_sync_rhs = 0;

// -------------------------------------------------------------
// Setup functions follow
// -------------------------------------------------------------

MacProj::MacProj(Amr * _parent, int _finest_level,
                 BCRec * _phys_bc, int _radius_grow )
  : parent(_parent), finest_level(_finest_level),
    phys_bc(_phys_bc), 
    radius_grow(_radius_grow), 
    LevelData(_finest_level+1),
    phi_bcs(_finest_level+1),
    mac_phi_crse(_finest_level+1, PArrayManage),
    mac_reg(_finest_level+1, PArrayManage),
    volume(_finest_level+1), area(_finest_level+1),
    radius(_finest_level+1) 
{
  read_params();

  if (verbose) {
    cout << "Creating mac_projector" << endl;
  }

  finest_level_allocated = finest_level;
}



MacProj::~MacProj()
{
}



void MacProj::read_params()
{
  // read parameters from input file and command line
  ParmParse pp("mac");

  verbose = pp.contains("v");
  pp.get(   "mac_tol",          mac_tol          );
  pp.get(   "mac_sync_tol",     mac_sync_tol     );
  pp.query( "use_cg_solve",     use_cg_solve     );
  pp.query( "mac_abs_tol",      mac_abs_tol      );
  pp.query( "do_outflow_bcs",   do_outflow_bcs   );
  pp.query( "fix_mac_sync_rhs", fix_mac_sync_rhs );
  pp.query( "use_viscosity",    use_viscosity    ); 
}



void MacProj::install_level(int level, AmrLevel * level_data,
			    MultiFab &_volume, MultiFab *_area,
			    PArray<REAL> * _radius )
{
  if (verbose) {
    cout << "Installing MacProj level " << level << endl;
  }

  if (parent->finestLevel() < finest_level) {
    for (int lev = parent->finestLevel() + 1; lev <= finest_level; lev++) {
      mac_reg.clear(lev);
    }
  }

  finest_level = parent->finestLevel();

  if (level > finest_level_allocated) {
    finest_level_allocated = finest_level;
    LevelData.resize(finest_level+1);
    phi_bcs.resize(finest_level+1);
    mac_phi_crse.resize(finest_level+1);
    mac_reg.resize(finest_level+1);
    volume.resize(finest_level+1);
    area.resize(finest_level+1);
    radius.resize(finest_level+1);
  }

  LevelData.clear(level);
  LevelData.set(level, level_data);
  volume.clear(level);
  volume.set(level, &_volume);
  area.set(level, _area);
  radius.clear(level);
  radius.set(level, _radius);

  BuildPhiBC(level);

  if (level > 0) {
    const BoxArray& grids = LevelData[level].boxArray();
    mac_reg.clear(level);
    mac_reg.set(level, new FluxRegister(grids,parent->refRatio(level-1),
					level,1));
  }
}



void MacProj::BuildPhiBC(int level)
{
  const BoxArray& grids = LevelData[level].boxArray();
  const Geometry& geom  = parent->Geom(level);

  int ngrds = grids.length();
  phi_bcs[level].resize(ngrds);
  const BOX& domain = geom.Domain();
  const int* domlo = domain.loVect();
  const int* domhi = domain.hiVect();

  const int* phys_lo = phys_bc->lo();
  const int* phys_hi = phys_bc->hi();

  for (int i = 0; i < ngrds; i++) {
    BCRec &bc = phi_bcs[level][i];
    const int* lo = grids[i].loVect();
    const int* hi = grids[i].hiVect();

    for (int dir = 0; dir < BL_SPACEDIM; dir++) {
      if (lo[dir] == domlo[dir]) {
	if (phys_lo[dir] == Outflow) {
	  bc.setLo(dir,LO_DIRICHLET);
	} else {
	  bc.setLo(dir,LO_NEUMANN);
	}
      } else {
	bc.setLo(dir,LO_DIRICHLET);
      }
      if (hi[dir] == domhi[dir]) {
	if (phys_hi[dir] == Outflow) {
	  bc.setHi(dir,LO_DIRICHLET);
	} else {
	  bc.setHi(dir,LO_NEUMANN);
	}
      } else {
	bc.setHi(dir,LO_DIRICHLET);
      }
    }
  }
}



void MacProj::setup(int level)
{
  const BoxArray& grids = LevelData[level].boxArray();
  if (level < parent->maxLevel()) {
    if (!mac_phi_crse.defined(level)) {
      mac_phi_crse.set(level, new MultiFab(grids,1,1,Fab_allocate));
      mac_phi_crse[level].setVal(0.0);
    }
  }
}



void MacProj::cleanup(int level)
{
  if (level < parent->maxLevel()) {
    mac_phi_crse.clear(level);
  }
}

// -------------------------------------------------------------
// Projection functions follow
// -------------------------------------------------------------

// ==================================================
// compute the level advance mac projection
// ==================================================

void MacProj::mac_project(int level, MultiFab* u_mac, MultiFab & S,
			  REAL dt, REAL time,
			  const MultiFab& divu, int have_divu)
{
  if (verbose) {
    cout << "... mac_project at level " << level << endl;
  }

  const BoxArray& grids = LevelData[level].boxArray();
  const Geometry& geom  = parent->Geom(level);
  IntVect crse_ratio = (level > 0) ? 
       parent->refRatio(level-1) : IntVect::TheZeroVector();
  const REAL* dx = geom.CellSize();

  MultiFab *mac_phi;
  int max_level=parent->maxLevel();


  // if finest level possible no need to make permanent mac_phi for bcs
  if (level == max_level) {
      mac_phi = new MultiFab(grids,1,1,Fab_allocate);
  } else {
      mac_phi = &mac_phi_crse[level];
  }

  mac_phi->setVal(0.0);
  //mac_phi->setCacheWidth(1);

  REAL hx = dx[0];

#if (BL_SPACEDIM == 2)
  int outflow_at_top = phys_bc->lo(0) != Outflow && phys_bc->lo(1) != Outflow && 
    phys_bc->hi(0) != Outflow && phys_bc->hi(1) == Outflow; 
  if (outflow_at_top && have_divu && do_outflow_bcs) {
    set_outflow_bcs(level, mac_phi, parent, u_mac, S, divu);
  }
#endif

  // store the Dirichlet boundary condition for mac_phi in mac_bndry
  //-----------------------------------------------------------------
  MacBndry mac_bndry(grids,1,geom);
  int src_comp = 0;
  int dest_comp = 0;
  int num_comp = 1;
  if (level == 0) {
    mac_bndry.setBndryValues(*mac_phi,src_comp,dest_comp,num_comp,*phys_bc);
  } else {
    MultiFab& CPhi = mac_phi_crse[level-1];
    BoxArray crse_boxes(grids);
    crse_boxes.coarsen(crse_ratio);
    int in_rad = 0;
    int out_rad = 1;
    int extent_rad = 1;
    BndryRegister crse_br(crse_boxes,in_rad,out_rad,extent_rad,num_comp);
    crse_br.copyFrom(CPhi,extent_rad,src_comp,dest_comp,num_comp);

    mac_bndry.setBndryValues(crse_br,src_comp,*mac_phi,src_comp,
			     dest_comp,num_comp,crse_ratio,
			     *phys_bc);
  }

  // compute the nondivergent velocities, by creating the linop
  // and multigrid operator appropriate for the solved system
  //-----------------------------------------------------------------

  // initialize the rhs with divu
  REAL rhs_scale = 2.0/dt;
  MultiFab Rhs(grids,1,0,Fab_allocate);
  Rhs.copy(divu);

  mac_level_driver( mac_bndry, grids,
                    use_cg_solve, level, Density,
                    dx, dt, mac_tol, mac_abs_tol, rhs_scale, 
                    area[level], volume[level],
                    S, Rhs, u_mac, mac_phi );

  // test that u_mac is divergence free
  if (verbose) check_div_cond(level, u_mac);

  
  // store advection velocities in mac registers at crse/fine boundaries
  //-----------------------------------------------------------------

  // initialize advection velocity registers with coarse grid velocity
  if (level < finest_level) {
    FluxRegister& mr = mac_reg[level+1];
    mr.setVal(0.0);
    REAL mult = -1.0;
    for (int dir = 0; dir < BL_SPACEDIM; dir++) {
      mr.CrseInit(u_mac[dir],area[level][dir],dir,0,0,1,mult);
    }
    if (verbose) 
      cout << "LEVEL " << level << " MACREG: CrseInit sum = " 
	   << mr.SumReg(0) << endl;
  }

  // increment in fine grid velocity to velocity registers
  if (level > 0) {
    REAL mult = 1.0/( (double) parent->MaxRefRatio(level-1) );
    for (int dir = 0; dir < BL_SPACEDIM; dir++) {
      mac_reg[level].FineAdd(u_mac[dir],area[level][dir],dir,0,0,1,mult);
    }
    if (verbose) {
      cout << "LEVEL " << level << " MACREG: FineAdd sum = " 
	   << mac_reg[level].SumReg(0) << endl;
    }
  }

  // if finest level possible no need to keep phi for boundary conditions
  if (level == max_level) {
    delete mac_phi;
    mac_phi = 0;
  }
}

// ==================================================
// compute the corrective pressure used in the mac_sync
// ==================================================

void MacProj::mac_sync_solve(int level, MultiFab* u_mac, 
			     REAL dt, MultiFab * rho_half, IntVect& fine_ratio)
{
  if (verbose) {
    cout << "... mac_sync_solve at level " << level << endl;
  }

  assert(level < finest_level);
  const BoxArray& grids = LevelData[level].boxArray();
  const Geometry& geom  = parent->Geom(level);
  IntVect crse_ratio = (level > 0) ? 
       parent->refRatio(level-1) : IntVect::TheZeroVector();
  const REAL* dx        = geom.CellSize();
  int ngrds             = grids.length();
  const BoxArray& fine_boxes = LevelData[level+1].boxArray();

  // Reusing storage here, since there should be no more need for the
  // values in mac_phi at this level and mac_sync_phi only need to last
  // into the call to mac_sync_compute.  Hope this works...  (LHH)
  MultiFab *mac_sync_phi = &mac_phi_crse[level];

  // alloc and define RHS by doing a reflux-like operation in coarse
  // grid cells adjacent to fine grids.  The values in these
  // cells should be SUM{MR/VOL} where the sum is taken over
  // all edges of a cell that adjoin fine grids, MR = value in
  // MAC register, VOL = cell volume.  All other cells have a
  // value of zero (including crse cells under fine grids).
  MultiFab Rhs(grids,1,0,Fab_allocate);
  Rhs.setVal(0.0);
  
  // Reflux subtracts values at hi edge of coarse cell and
  // adds values at lo edge.  We want the opposite here so
  // set scale to -1
  // alloc space for Rhs
  FluxRegister& mr = mac_reg[level+1];
  REAL scale       = -1.0;
  mr.Reflux(Rhs,volume[level],scale,0,0,1,geom);

  int nfine = fine_boxes.length();
  for (int kf = 0; kf < nfine; kf++) {
      BOX bf(fine_boxes[kf]);
      bf.coarsen(fine_ratio);
      for (int k = 0; k < ngrds; k++) {
          BOX bx(grids[k]);
          bx &= bf;
          if (bx.ok()) {
              Rhs[k].setVal(0.0,bx,0);
          }
      }
  }

#if 1
  // remove constant null space component from the rhs of the solve
  // this code should go away when Marc makes this option
  // part of the multigrid code--rbp, 2/13/97
  if ( fix_mac_sync_rhs ) {
    REAL sum = 0;
    long size = 0;
	int i;
    for (i = 0; i < Rhs.length(); i++ ) {
      sum += Rhs[i].sum( 0, 1 );
      size += Rhs[i].box().numPts();
    }
    REAL fix = sum / size;
    cout << "Average correction = " << fix << endl;
    Rhs.plus( -fix, 0 );
    for (i = 0, sum = 0; i < Rhs.length(); i++ ) {
      sum += Rhs[i].sum( 0, 1 );
    }
    cout << "...new sum = " << sum << endl;
  }
#endif

  mac_sync_phi->setVal(0.0);
  //mac_sync_phi->setCacheWidth(1);

  // store the Dirichlet boundary condition for mac_sync_phi in mac_bndry
  //-----------------------------------------------------------------
  MacBndry mac_bndry(grids,1,geom);
  int src_comp = 0;
  int dest_comp = 0;
  int num_comp = 1;
  if (level == 0) {
    mac_bndry.setBndryValues(*mac_sync_phi,src_comp,dest_comp,num_comp,
			     *phys_bc);
  } else {
    BoxArray crse_boxes(grids);
    crse_boxes.coarsen(crse_ratio);
    int n_ghost = 0;
    BndryRegister crse_br(crse_boxes,n_ghost,1,1,num_comp);
    crse_br.setVal(0.);
    mac_bndry.setBndryValues(crse_br,src_comp,*mac_sync_phi,src_comp,
			     dest_comp,num_comp,crse_ratio,
			     *phys_bc);
  }

  // now define edge centered coefficients and adjust RHS
  // for MAC solve
  REAL rhs_scale = 2.0/dt;

  // solve the sync system
  //-----------------------------------------------------------------
  mac_sync_driver( mac_bndry,    grids,       
                   use_cg_solve, level, 
                   dx, dt, mac_sync_tol, mac_abs_tol, rhs_scale, 
                   area[level],  volume[level],
                   Rhs, rho_half, u_mac, mac_sync_phi );
}




// ================================================================
// after solving for mac_sync_phi in mac_sync_solve(), we
// can now do the sync advect step.  This consists of two steps
//
// 1.  compute u_corr as the gradient of mac_sync_phi
// 2.  compute advective tendency of u_corr and
//     add into Vsync or Ssync
//
// If increment_sync is non-null, the (i-BL_SPACEDIM)-th component 
// of (*Ssync) is incremented only when increment[i]==1
// This is useful if that component gets incrmnted in a non-standard
// way.
// ================================================================
void MacProj::mac_sync_compute(int level, MultiFab * u_mac, 
			       MultiFab * Vsync,  MultiFab * Ssync,
			       MultiFab * rho_half,
			       FluxRegister* adv_flux_reg,
			       Array<int> is_conservative, REAL prev_time, 
			       REAL pres_prev_time, REAL dt, 
                               int NUM_STATE, REAL be_cn_theta,
                               const int* increment_sync)
{
    // work variables
    int i,dir,comp;
    BOX b;
    FARRAYBOX S, Rho, tforces, Gp;
    FARRAYBOX xflux, yflux, zflux, divu;
    FARRAYBOX grad_phi[BL_SPACEDIM];
    
    // get parameters
    const BoxArray& grids  = LevelData[level].boxArray();
    const Geometry& geom   = parent->Geom(level);
    const REAL *dx         = geom.CellSize();
    int numscal            = NUM_STATE - BL_SPACEDIM;
    MultiFab *mac_sync_phi = &mac_phi_crse[level];
    NavierStokes& ns_level =  *(NavierStokes*) &(parent->getLevel(level));

    // get the godunov object
    Godunov *godunov = ns_level.godunov;

    int use_forces_in_trans = godunov->useForcesInTrans();
    FARRAYBOX tvelforces;
    MultiFab vel_visc_terms;
    if(use_forces_in_trans) {
      vel_visc_terms.define(grids,BL_SPACEDIM,1,Fab_allocate);
      ns_level.getViscTerms(vel_visc_terms,Xvel,BL_SPACEDIM,prev_time);
      if(be_cn_theta==1.0)vel_visc_terms.setVal(0.0,1);
    }

    // get viscous forcing
    MultiFab *visc_terms = NULL;
    if ( use_viscosity ) {
        visc_terms = new MultiFab(grids,NUM_STATE,1,Fab_allocate);
        if (be_cn_theta!=1.0) {
          ns_level.getViscTerms((*visc_terms),Xvel,NUM_STATE,prev_time);
        } else {
          visc_terms->setVal(0.0,1);
        }
    }

    // compute the mac sync correction
    for ( i = 0; i < grids.length(); i++) {

        // get the bounds
	const BOX &grd = grids[i];

        // Step 1: compute ucorr = grad(phi)/rhonph
        //---------------------------------------------------
        
        // create storage for corrective velocities
        grad_phi[0].resize(surroundingNodes(grd,0),1);
        grad_phi[1].resize(surroundingNodes(grd,1),1);
#if (BL_SPACEDIM == 3 )
        grad_phi[2].resize(surroundingNodes(grd,2),1);
#endif

        REAL mult = dt/2.0;
        mac_vel_update( 1,
                        grad_phi[0],
                        grad_phi[1],
#if (BL_SPACEDIM == 3 )
                        grad_phi[2],
#endif
                        (*mac_sync_phi)[i],
                        &(*rho_half)[i], 0,
                        grids[i], level, i, dx, mult );
        
        // Step 2: compute Mac correction by calling GODUNOV box
        //--------------------------------------------------------
        
        // get needed data
        ns_level.getState(S,      i,HYP_GROW,0,NUM_STATE,prev_time);
        ns_level.getForce(tforces,i,1,       0,NUM_STATE,prev_time);
        ns_level.getGradP(Gp,     i,1,              pres_prev_time);
        ns_level.getDivCond(divu, i,1,                   prev_time);
	Rho.resize(grow(grd,1),1);
	Rho.copy(S,Density,0,1);

        // compute total forcing terms
        if ( use_viscosity ) {
            godunov->Sum_tf_gp_visc( tforces, (*visc_terms)[i], Gp, Rho );
            godunov->Sum_tf_divu_visc( S, tforces,
                                       BL_SPACEDIM, numscal,
                                       (*visc_terms)[i],
                                       BL_SPACEDIM,
                                       divu, Rho, 1 );
        } else {
            godunov->Sum_tf_gp(   tforces, Gp, Rho );
            godunov->Sum_tf_divu( S, tforces,
                                  BL_SPACEDIM, numscal,
                                  divu, Rho, 1 );
        }

        if (use_forces_in_trans) {
	  ns_level.getForce(tvelforces,i,1,       Xvel,   BL_SPACEDIM,prev_time);
          godunov->Sum_tf_gp_visc( tvelforces, vel_visc_terms[i], Gp, Rho );
        }
        
        // set up the workspace for the godunov BOX
        godunov->Setup( grd, dx, dt, 0,
                        xflux, ns_level.getBCArray( State_Type,i,0,1),
                        yflux, ns_level.getBCArray( State_Type,i,1,1),
#if (BL_SPACEDIM == 3)
                        zflux, ns_level.getBCArray( State_Type,i,2,1),
#endif
                        S, Rho, tvelforces );

        // get the sync FABS
        FARRAYBOX& u_sync = (*Vsync)[i];
        FARRAYBOX& s_sync = (*Ssync)[i];
        
        // loop over the state components and
        // compute the sync advective component
        for ( comp = 0 ; comp < NUM_STATE ; comp++ ) {

            int do_comp = (increment_sync==NULL);
            if(!do_comp) {
              do_comp = comp < BL_SPACEDIM || (increment_sync[comp]==1);
            }
            if (!do_comp) continue; // skip to next value of comp

            int u_ind       = comp;
            int s_ind       = comp-BL_SPACEDIM;
            int sync_ind    = ( comp < BL_SPACEDIM ? u_ind  : s_ind  );
            FARRAYBOX &temp = ( comp < BL_SPACEDIM ? u_sync : s_sync );
            godunov->SyncAdvect( grd, dx, dt, level,
                                 area[level][0][i], u_mac[0][i],
                                 grad_phi[0],       xflux,
                                 
                                 area[level][1][i], u_mac[1][i],
                                 grad_phi[1],       yflux,
#if (BL_SPACEDIM == 3 )                            
                                 area[level][2][i], u_mac[2][i],
                                 grad_phi[2],       zflux,
#endif
                                 S, tforces, comp,
                                 temp,       sync_ind,
                                 is_conservative[comp],
                                 comp,
                                 ns_level.getBCArray( State_Type,i,comp,1),
                                 volume[level][i] );
            
            // NOTE: the signs here are opposite from VELGOD
            // NOTE: fluxes expected to be in extensive form
            if (level > 0) {
                
                adv_flux_reg->FineAdd(xflux,0,i,0,comp,1,-dt);
                adv_flux_reg->FineAdd(yflux,1,i,0,comp,1,-dt);
#if (BL_SPACEDIM == 3 )
                adv_flux_reg->FineAdd(zflux,2,i,0,comp,1,-dt);
#endif
            }
        } // end of advection loop

        // include grad_phi in the mac registers corresponding
        // to the next coarsest interface
        if (level > 0) {
            REAL mult =  -1.0/( (double) parent->MaxRefRatio(level-1));
            for (int dir = 0; dir < BL_SPACEDIM; dir++) {
                mac_reg[level].FineAdd(grad_phi[dir], area[level][dir][i],
                                       dir, i, 0, 0, 1, mult );
            }
        }

        // multiply the sync term by dt--this is now done in the
        // calling routine
        
    } // end of grid loop

    // garbage collection
    if ( visc_terms != NULL )
        delete visc_terms;
}

// ---------------------------------------------------------------

// ==================================================
// This routine does a sync advect step for a single 
// scalar component. Unlike the preceding routine, the
// half-time edge states are passed in from the calling routine.
// This routine is useful when the edge states are computed
// in a physics-class-specific manner. (For example, as they are
// in the calculation of div rho U h = div U sum_l (rho Y)_l h_l(T)).
// ==================================================

void MacProj::mac_sync_compute(int level, MultiFab * u_mac, 
			       MultiFab * Ssync, int comp,
                               MultiFab** sync_edges,
			       MultiFab * rho_half,
			       FluxRegister* adv_flux_reg,
			       Array<int> is_conservative, 
			       REAL dt)
{
    assert (comp>=BL_SPACEDIM);

    int s_ind = comp - BL_SPACEDIM;

    // work variables
    int i;
    FARRAYBOX xflux, yflux, zflux;
    FARRAYBOX grad_phi[BL_SPACEDIM];
    
    // get parameters
    const BoxArray& grids  = LevelData[level].boxArray();
    const Geometry& geom   = parent->Geom(level);
    const REAL *dx         = geom.CellSize();
    MultiFab *mac_sync_phi = &mac_phi_crse[level];
    NavierStokes& ns_level =  *(NavierStokes*) &(parent->getLevel(level));

    // create advection object
    Godunov godunov( 512 );

    // compute the mac sync correction
    for ( i = 0; i < grids.length(); i++) {

        // get the bounds
	const BOX &grd = grids[i];

        // Step 1: compute ucorr = grad(phi)/rhonph
        //---------------------------------------------------
        
        // create storage for corrective velocities
        grad_phi[0].resize(surroundingNodes(grd,0),1);
        grad_phi[1].resize(surroundingNodes(grd,1),1);
#if (BL_SPACEDIM == 3 )
        grad_phi[2].resize(surroundingNodes(grd,2),1);
#endif

        REAL mult = dt/2.0;
        mac_vel_update( 1,
                        grad_phi[0],
                        grad_phi[1],
#if (BL_SPACEDIM == 3 )
                        grad_phi[2],
#endif
                        (*mac_sync_phi)[i],
                        &(*rho_half)[i], 0,
                        grids[i], level, i, dx, mult );
        
        // Step 2: compute Mac correction by advecting the edge states
        //--------------------------------------------------------
        
        xflux.resize(surroundingNodes(grd,0),1);
        yflux.resize(surroundingNodes(grd,1),1);
        xflux.copy((*sync_edges[0])[i]);
        yflux.copy((*sync_edges[1])[i]);
#if (BL_SPACEDIM == 3 )
        zflux.resize(surroundingNodes(grd,2),1);
        zflux.copy((*sync_edges[2])[i]);
#endif

        // get the sync FABS
        FARRAYBOX& sync = (*Ssync)[i];

        godunov.ComputeSyncAofs( grd,
                          area[level][0][i],
                          grad_phi[0],       xflux,
                                
                          area[level][1][i],
                          grad_phi[1],       yflux,
#if (BL_SPACEDIM == 3 )                            
                          area[level][2][i],
                          grad_phi[2],       zflux,
#endif
                          volume[level][i], sync, s_ind, 
                          is_conservative[comp]);

            // NOTE: the signs here are opposite from VELGOD
            // NOTE: fluxes expected to be in extensive form
         if (level > 0) {
                
           adv_flux_reg->FineAdd(xflux,0,i,0,comp,1,-dt);
           adv_flux_reg->FineAdd(yflux,1,i,0,comp,1,-dt);
#if (BL_SPACEDIM == 3 )
           adv_flux_reg->FineAdd(zflux,2,i,0,comp,1,-dt);
#endif
         }

        // multiply the sync term by dt--this is done in the 
        // calling routine
        
    } // end of grid loop

}

// check the mac divergence 
void MacProj::check_div_cond(int level, MultiFab U_edge[]) const
{
  const BoxArray& grids = LevelData[level].boxArray();

  int ngrds = U_edge[0].length();
  REAL sum = 0.0;

  for (int k = 0; k < ngrds; k++) {
    FARRAYBOX dmac(grids[k],1);
    const FARRAYBOX& uxedge = U_edge[0][k];
    const FARRAYBOX& uyedge = U_edge[1][k];
    const FARRAYBOX& xarea = area[level][0][k];
    const FARRAYBOX& yarea = area[level][1][k];
    const FARRAYBOX& vol = volume[level][k];

    DEF_LIMITS(dmac,dmac_dat,dlo,dhi);
    DEF_CLIMITS(uxedge,ux_dat,uxlo,uxhi);
    DEF_CLIMITS(uyedge,uy_dat,uylo,uyhi);
    DEF_CLIMITS(xarea,ax_dat,axlo,axhi);
    DEF_CLIMITS(yarea,ay_dat,aylo,ayhi);
    DEF_CLIMITS(vol,vol_dat,vlo,vhi);
#if (BL_SPACEDIM == 2)
    FORT_MACDIV(dmac_dat,ARLIM(dlo),ARLIM(dhi),dlo,dhi,
		ux_dat,ARLIM(uxlo),ARLIM(uxhi),
		uy_dat,ARLIM(uylo),ARLIM(uyhi),
		ax_dat,ARLIM(axlo),ARLIM(axhi), 
		ay_dat,ARLIM(aylo),ARLIM(ayhi), 
		vol_dat,ARLIM(vlo),ARLIM(vhi));
#endif
#if (BL_SPACEDIM == 3)

    const FARRAYBOX& uzedge = U_edge[2][k];
    DEF_CLIMITS(uzedge,uz_dat,uzlo,uzhi);
    const FARRAYBOX& zarea = area[level][2][k];
    DEF_CLIMITS(zarea,az_dat,azlo,azhi);

    FORT_MACDIV(dmac_dat,ARLIM(dlo),ARLIM(dhi),dlo,dhi,
		ux_dat,ARLIM(uxlo),ARLIM(uxhi),
		uy_dat,ARLIM(uylo),ARLIM(uyhi),
		uz_dat,ARLIM(uzlo),ARLIM(uzhi),
		ax_dat,ARLIM(axlo),ARLIM(axhi),
		ay_dat,ARLIM(aylo),ARLIM(ayhi),
		az_dat,ARLIM(azlo),ARLIM(azhi),
		vol_dat,ARLIM(vlo),ARLIM(vhi));
#endif
    if (verbose) {
      REAL g_norm = dmac.norm(0);
      cout << "Max norm of div(U_edge) for grid  " << k << " = "
	   << g_norm << endl;
      sum += dmac.sum(0);
      cout << "SUM of DIV(U_edge) = " << sum << endl;
    }
  }
}



void MacProj::set_outflow_bcs(int level,
			      MultiFab* mac_phi, Amr* parent, MultiFab* u_mac, 
                              MultiFab & S, const  MultiFab& divu)
{
  const BoxArray& grids = LevelData[level].boxArray();
  const Geometry& geom  = parent->Geom(level);

  // so far only for 2-d
#if (BL_SPACEDIM == 2)
  int ngrds = grids.length();

  const REAL* dx = geom.CellSize();
  REAL hx = dx[0];
  if(level==0 && ngrds==1) {

    // WORKS FOR SINGLE GRID ONLY

    int i = 0;
    BOX rbox = grids[i];
    int rlen  = rbox.length(0);
    BOX redge_box(rbox);
    redge_box.growHi(0,1);
    Array<REAL> rcen;
    rcen.resize(rlen);
    Array<REAL> redge;
    redge.resize(rlen+1);
    if (CoordSys::IsRZ() == 1) {
      geom.GetCellLoc(rcen,rbox, 0);
      geom.GetEdgeLoc(redge,rbox, 0);
    } else {
      int ii;
      for(ii=0; ii<rlen;ii++) {
	rcen[ii] = 1.0;
	redge[ii] = 1.0;
      }
      redge[rlen] = 1.0;
    }
    const int* rcen_lo = rbox.loVect();
    const int* rcen_hi = rbox.hiVect();
    const int* redge_lo = rbox.loVect();
    const int* redge_hi = redge_box.hiVect();
    DEF_CLIMITS(divu[i],divudat,divu_lo,divu_hi);
    DEF_CLIMITS((*u_mac)[i],udat,u_lo,u_hi);
    REAL* uhalfx = u_mac[0][i].dataPtr();
    DEF_CLIMITS(S[i],rho,rho_lo,rho_hi);
    rho = S[i].dataPtr(Density);
    DEF_LIMITS((*mac_phi)[i],phi,phi_lo,phi_hi);
    FORT_MACPHIBC(ARLIM(u_lo),ARLIM(u_hi),uhalfx,
		  ARLIM(divu_lo),ARLIM(divu_hi),divudat,
		  ARLIM(rho_lo),ARLIM(rho_hi),rho,
		  ARLIM(rcen_lo),ARLIM(rcen_hi),rcen.dataPtr(),
		  ARLIM(redge_lo),ARLIM(redge_hi),redge.dataPtr(),
		  &hx,ARLIM(phi_lo),ARLIM(phi_hi),phi);
  } else if(level==0) {
    const REAL* dx = parent->Geom(0).CellSize();
    REAL hx = dx[0];
    const BOX& domain = parent->Geom(0).Domain();
    INTVECT bigend = domain.bigEnd();
    INTVECT smallend = domain.smallEnd();
    int jhi = domain.bigEnd(1);
    smallend.setVal(1,jhi);
    BOX top_strip(smallend,bigend,IntVect::TheCellVector());

    FARRAYBOX divu_strip(top_strip,1);
    BOX top_vel_strip = top_strip;
    top_vel_strip.growLo(0,1);
    top_vel_strip.shiftHalf(0,1);
    FARRAYBOX mac_vel_strip(top_vel_strip,1);
    BOX top_rho_strip = top_strip;
    top_rho_strip.grow(1);
    FARRAYBOX rho_strip(top_rho_strip,1);
    BOX top_phi_strip = top_strip;
    top_phi_strip.grow(0,1);
    top_phi_strip.growHi(1,1);
    FARRAYBOX mac_phi_strip(top_phi_strip,1);

    mac_phi_strip.setVal(0.0);
    int i;
    for (i=0;i<ngrds;i++) {
      BOX destbox = grids[i];
      destbox.grow(0,1);
      destbox &= top_rho_strip;
      BOX srcbox = destbox;
      if(destbox.ok()) {
	rho_strip.copy(S[i],srcbox,Density,destbox,0,1);
      }

      destbox = grids[i];
      destbox &= top_strip;
      srcbox = destbox;
      if(destbox.ok()) {
	divu_strip.copy(divu[i],srcbox,0,destbox,0,1);
      }

      destbox = grids[i];
      destbox.growLo(0,1);
      destbox.shiftHalf(0,1);
      destbox &= top_vel_strip;
      srcbox = destbox;
      if(destbox.ok()) {
	mac_vel_strip.copy(u_mac[0][i],srcbox,0,destbox,0,1);
      }
    }

    BOX rbox = top_strip;
    int rlen  = rbox.length(0);
    BOX redge_box(rbox);
    redge_box.growHi(0,1);
    Array<REAL> rcen;
    rcen.resize(rlen);
    Array<REAL> redge;
    redge.resize(rlen+1);
    if (CoordSys::IsRZ() == 1) {
      geom.GetCellLoc(rcen,rbox, 0);
      geom.GetEdgeLoc(redge,rbox, 0);
    } else {
      for(i=0; i<rlen;i++) {
	rcen[i] = 1.0;
	redge[i] = 1.0;
      }
      redge[rlen] = 1.0;
    }
    const int* rcen_lo  = rbox.loVect();
    const int* rcen_hi  = rbox.hiVect();
    const int* redge_lo = rbox.loVect();
    const int* redge_hi = redge_box.hiVect();
    DEF_CLIMITS(rho_strip,rhodat,rho_lo,rho_hi);
    DEF_CLIMITS(divu_strip,divudat,divu_lo,divu_hi);
    DEF_LIMITS(mac_phi_strip,phidat,phi_lo,phi_hi);
    DEF_CLIMITS(mac_vel_strip,udat,ulo,uhi);
    FORT_MACPHIBC(ARLIM(ulo),ARLIM(uhi),udat,
		  ARLIM(divu_lo),ARLIM(divu_hi),divudat,
		  ARLIM(rho_lo),ARLIM(rho_hi),rhodat,
		  ARLIM(rcen_lo),ARLIM(rcen_hi),rcen.dataPtr(),
		  ARLIM(redge_lo),ARLIM(redge_hi),redge.dataPtr(),
		  &hx,ARLIM(phi_lo),ARLIM(phi_hi),phidat);
    for (i=0;i<ngrds;i++) {
      (*mac_phi)[i].copy(mac_phi_strip);
    }
  }
#endif
}












