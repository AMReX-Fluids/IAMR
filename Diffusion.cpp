// $Id: Diffusion.cpp,v 1.1 1997-07-08 23:08:12 vince Exp $

// comment out this line to use diffusion class outside
// the context of NavierStokes and classes derived from it
#define USE_NAVIERSTOKES 1

#include <stdio.h>
#include <strstream.h>
#include <Box.H>
#include <BoxArray.H>
#include <Geometry.H>
#include <Misc.H>
#include <CoordSys.H>
#include <Geometry.H>
#include <BoxDomain.H>
#include <ParmParse.H>
#include <ErrorList.H>
#ifdef USE_NAVIERSTOKES
#include <NavierStokes.H>
#endif
#include <Diffusion.H>
#include <MultiGrid.H>
#include <CGSolver.H>

#include <DIFFUSION_F.H>
#include <VISCOPERATOR_F.H>

// include files for tensor solve
#if (BL_SPACEDIM==2) && defined (USE_TENSOR)
#include <DivVis.H>
#include <LO_BCTYPES.H>
#include <MCMultiGrid.H>
#include <MCCGSolver.H>
#include <ViscBndry2D.H>
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
#define bogus_value 1.e20

Array<int>  Diffusion::is_diffusive;
Array<REAL> Diffusion::visc_coef;
REAL      Diffusion::visc_tol = 1.0e-10;      // tolerance for viscous solve
REAL      Diffusion::visc_abs_tol = 1.0e-10;  // absolute tol. for visc solve

int  Diffusion::first = 1;
int  Diffusion::do_reflux = 1;
int  Diffusion::use_cg_solve = 0;
int  Diffusion::use_tensor_cg_solve = 0;
bool Diffusion::use_mg_precond_flag = false;
int  Diffusion::verbose = 0;
int  Diffusion::max_order = 2;
int  Diffusion::tensor_max_order = 2;
int  Diffusion::use_dv_constant_mu_def = 1;
int  Diffusion::scale_abec = 0;

Diffusion::Diffusion(Amr* Parent, AmrLevel* Caller, Diffusion* Coarser,
		     int num_state, FluxRegister *Viscflux_reg,
		     MultiFab& Volume, MultiFab* Area,
                     Array<int>  _is_diffusive,
                     Array<REAL> _visc_coef)
  : parent(Parent), caller(Caller), grids(caller->boxArray()),
    level(caller->Level()),
    coarser(Coarser), finer(NULL), NUM_STATE(num_state),
    viscflux_reg(Viscflux_reg), volume(Volume), area(Area)
{
  if (first) {
    first = 0;

    ParmParse ppdiff("diffuse");
    ppdiff.query("use_cg_solve",use_cg_solve);
    ppdiff.query("use_tensor_cg_solve",use_tensor_cg_solve);
    ppdiff.query("use_dv_constant_mu",use_dv_constant_mu_def);
    int use_mg_precond = 0;
    ppdiff.query("use_mg_precond",use_mg_precond);
    use_mg_precond_flag = (use_mg_precond ? true : false);
    ppdiff.query("max_order",max_order);
    ppdiff.query("tensor_max_order",tensor_max_order);
    ppdiff.query("scale_abec",scale_abec);
    use_mg_precond_flag = (use_mg_precond ? true : false);

    ParmParse pp("ns");

    pp.query("do_reflux",do_reflux);
    do_reflux = (do_reflux ? 1 : 0);

    pp.query("visc_tol",visc_tol);
    pp.query("visc_abs_tol",visc_abs_tol);

    int n_visc = _visc_coef.length();
    int n_diff = _is_diffusive.length();
    if (n_diff < NUM_STATE || n_visc < NUM_STATE) {
      cout << "Diffusion::Diffusion : is_diffusive and/or visc_coef arrays are " <<
              " not long enough" << endl;
      exit(0);
    }
    visc_coef.resize(NUM_STATE);
    is_diffusive.resize(NUM_STATE);
    int i;
    for (i = 0; i < NUM_STATE; i++) {
        is_diffusive[i] = _is_diffusive[i];
        visc_coef[i] = _visc_coef[i];
    }

    verbose = 1;
  }

  if (level > 0) {
    crse_ratio = parent->refRatio(level-1);
    coarser->finer = this;
  }
}

Diffusion::~Diffusion()
{
}

void Diffusion::diffuse_scalar(REAL dt, int sigma, REAL be_cn_theta,
			       MultiFab* rho_half, int rho_flag,
                               int do_viscreflux,
                               MultiFab* delta_rhs, 
                               MultiFab* alpha, 
                               MultiFab** betan, 
                               MultiFab** betanp1)
{

  if (verbose) {
    cout << "... diffuse_scalar " << sigma << endl;
  }

  int allnull, allthere;
  checkBetas(betan, betanp1, allthere, allnull);

  int finest_level = parent->finestLevel();

      // at this point, S_old has bndry at time N

  MultiFab &S_old = caller->get_old_data(State_Type);
  MultiFab &S_new = caller->get_new_data(State_Type);

  int ngrd = grids.length();
  int i;

  MultiFab* Rho_old;
  MultiFab* Rho_new;
  if(rho_flag==2) {
    Rho_old = new MultiFab(grids,1,1,Fab_allocate);
    Rho_old->copy(S_old,Density,0,1,1);
    Rho_new = new MultiFab(grids,1,1,Fab_allocate);
    Rho_new->copy(S_new,Density,0,1,1);

  }

  const BOX& domain = caller->Geom().Domain();
  const int* domlo = domain.loVect();
  const int* domhi = domain.hiVect();
  const REAL* dx = caller->Geom().CellSize();

  REAL cur_time  = caller->get_state_data(State_Type).curTime();
  REAL prev_time = caller->get_state_data(State_Type).prevTime();

  FARRAYBOX xflux,yflux,zflux;

  // We are now at this point in NS::scalar_update:
  //if (is_diffusive[sigma]) {

  // S_new now contains s(n) - dt*aofs.  We will use
  // this as part of the RHS for the viscous solve

  MultiFab Rhs(grids,1,0,Fab_allocate);
  MultiFab Soln(grids,1,1,Fab_allocate);

  REAL mf_norm = 0;

  { // set up Rhs

    REAL a = 0.0;
    REAL b;
    if(allnull) {
      b = -(1.0-be_cn_theta)*visc_coef[sigma]*dt;
    } else {
      b = -(1.0-be_cn_theta)*dt;
    }

    ViscBndry visc_bndry;
    ABecLaplacian *visc_op =   
      getViscOp(sigma,a,b,prev_time,visc_bndry,rho_half,rho_flag,NULL,betan);
    visc_op->maxOrder(max_order);

    if(rho_flag==2) {
// we are going to solve for S, not rho*S
      for (i = 0; i < ngrd; i++) {
        BOX box =grids[i];
        box.grow(1); 
        S_old[i].divide((*Rho_old)[i],box,0,sigma,1);
      }
    }

    // copy to single-component multifab
    // note: use Soln as a temporary here
    Soln.copy(S_old,sigma,0,1,1);

    visc_op->apply(Rhs,Soln);
    delete visc_op;

    // copy back to S_old for use in creating viscous fluxes 
    // NOTE: this requires that the visc_op->apply call returns
    //       Soln with the ghost cells correctly filled
    S_old.copy(Soln,0,sigma,1,1);

    // complete Rhs by adding body sources
    for (i = 0; i < ngrd; i++) {
      // scale inviscid part by volume
      S_new[i].mult(volume[i],grids[i],0,sigma,1);

      if(rho_flag==1) {
        // multiply by density at time nph
        FARRAYBOX& Rh = (*rho_half)[i];
        S_new[i].mult(Rh,grids[i],0,sigma,1);
      }

      if(alpha!=NULL) {
        S_new[i].mult((*alpha)[i],grids[i],0,sigma,1);
      }

      // add to rhs
      Rhs[i].plus(S_new[i],grids[i],sigma,0,1);
      REAL gr_norm = Rhs[i].norm(0);
      mf_norm = Max(gr_norm,mf_norm);
    }
    if (delta_rhs!=NULL) {
      mf_norm = 0.0;
      for (i = 0; i < ngrd; i++) {
        (*delta_rhs)[i].mult(dt);
        (*delta_rhs)[i].mult(volume[i],grids[i],0,0,1);
        Rhs[i].plus((*delta_rhs)[i],grids[i],0,0,1);
        REAL gr_norm = Rhs[i].norm(0);
        mf_norm = Max(gr_norm,mf_norm);
      }
    }
  }

  // construct viscous operator with bndry data at time N+1
  REAL a = 1.0;
  REAL b;
  if(allnull) {
    b = be_cn_theta*dt*visc_coef[sigma];
  } else {
    b = be_cn_theta*dt;
  }

  ViscBndry visc_bndry;
  REAL rhsscale = 1.0;
  ABecLaplacian *visc_op = getViscOp(sigma,a,b,cur_time,visc_bndry,
				     rho_half,rho_flag,&rhsscale,betanp1,alpha);
  Rhs.mult(rhsscale,0,1);
  visc_op->maxOrder(max_order);

  // construct solver and call it
  REAL S_tol = visc_tol;
  REAL S_tol_abs = visc_abs_tol*mf_norm;
  if(use_cg_solve) {
    CGSolver cg(*visc_op,use_mg_precond_flag);
    cg.solve(Soln,Rhs,S_tol,S_tol_abs);
  } else {
    MultiGrid mg(*visc_op);
    mg.solve(Soln,Rhs,S_tol,S_tol_abs);
  }

  int visc_op_lev = 0;
  visc_op->applyBC(Soln,visc_op_lev);

  delete visc_op;

  // copy into state variable at new time, with bcs hopefully

  S_new.copy(Soln,0,sigma,1,1);

  // create diffusive fluxes here

  if (do_reflux && do_viscreflux) {
    for (i = 0; i < ngrd; i++) {
      const BOX& grd = grids[i];
      const int* lo = grd.loVect();
      const int* hi = grd.hiVect();

      const int* slo = S_old[i].loVect();
      const int* shi = S_old[i].hiVect();

      BOX xflux_bx(grd);
      xflux_bx.surroundingNodes(0);
      xflux.resize(xflux_bx,1);
      DEF_LIMITS(xflux,xflux_dat,xflo,xfhi);

      BOX yflux_bx(grd);
      yflux_bx.surroundingNodes(1);
      yflux.resize(yflux_bx,1);
      DEF_LIMITS(yflux,yflux_dat,yflo,yfhi);

      FARRAYBOX& xarea = area[0][i];
      FARRAYBOX& yarea = area[1][i];

      DEF_CLIMITS(xarea,xarea_dat,axlo,axhi);
      DEF_CLIMITS(yarea,yarea_dat,aylo,ayhi);

      REAL mult;
      if(allnull) {
        mult = -visc_coef[sigma];
      } else {
        mult = -1.0;
      }

#if (BL_SPACEDIM == 2)
        FORT_VISCFLUX (S_old[i].dataPtr(sigma), 
		     S_new[i].dataPtr(sigma), 
                     ARLIM(slo), ARLIM(shi),
		     lo,hi,
                     xflux_dat,ARLIM(xflo),ARLIM(xfhi),
                     yflux_dat,ARLIM(yflo),ARLIM(yfhi),
                     xarea_dat,ARLIM(axlo),ARLIM(axhi),
                     yarea_dat,ARLIM(aylo),ARLIM(ayhi),
		     dx,&mult,&be_cn_theta);
#endif
#if (BL_SPACEDIM == 3)

        BOX zflux_bx(grd);
        zflux_bx.surroundingNodes(2);
        zflux.resize(zflux_bx,1);
        DEF_LIMITS(zflux,zflux_dat,zflo,zfhi);

        FARRAYBOX& zarea = area[2][i];
        DEF_CLIMITS(zarea,zarea_dat,azlo,azhi);

        FORT_VISCFLUX (S_old[i].dataPtr(sigma), 
		     S_new[i].dataPtr(sigma), 
                     ARLIM(slo), ARLIM(shi),
		     lo,hi,
                     xflux_dat,ARLIM(xflo),ARLIM(xfhi),
                     yflux_dat,ARLIM(yflo),ARLIM(yfhi),
                     zflux_dat,ARLIM(zflo),ARLIM(zfhi),
                     xarea_dat,ARLIM(axlo),ARLIM(axhi),
                     yarea_dat,ARLIM(aylo),ARLIM(ayhi),
                     zarea_dat,ARLIM(azlo),ARLIM(azhi),
		     dx,&mult,&be_cn_theta);
#endif
      if(allthere) {
        xflux.mult((*betanp1[0])[i]);
        yflux.mult((*betanp1[1])[i]);
#if (BL_SPACEDIM == 3)
        zflux.mult((*betanp1[2])[i]);
#endif
      }

      // NOTE: fluxes expected to be in extensive form
      if (level < finest_level) {
	FluxRegister& fr = *finer->viscflux_reg;
	fr.CrseInit(xflux,xflux_bx,0,0,sigma,1,-dt);
	fr.CrseInit(yflux,yflux_bx,1,0,sigma,1,-dt);
#if (BL_SPACEDIM == 3)
	fr.CrseInit(zflux,zflux_bx,2,0,sigma,1,-dt);
#endif
      }
      if (level > 0) {
	viscflux_reg->FineAdd(xflux,0,i,0,sigma,1,dt);
	viscflux_reg->FineAdd(yflux,1,i,0,sigma,1,dt);
#if (BL_SPACEDIM == 3)
	viscflux_reg->FineAdd(zflux,2,i,0,sigma,1,dt);
#endif
      }
    } 
  } 

  if(rho_flag==2) { 
// return S_old to hold rho*S
// also we solved for S, not rho*S, so fix S_new
    for (i = 0; i < ngrd; i++) { 
      BOX box = grids[i];
      box.grow(1); 
      S_new[i].mult((*Rho_new)[i],box,0,sigma,1);
      S_old[i].mult((*Rho_old)[i],box,0,sigma,1);
    }
    delete Rho_old;
    delete Rho_new;
  }

}

void Diffusion::diffuse_velocity(REAL dt, REAL be_cn_theta,
				 MultiFab* rho_half, int rho_flag,
                                 MultiFab* delta_rhs,
                                 MultiFab** betan, 
                                 MultiFab** betanp1)

{
  if (verbose) {
    cout << "... diffuse_velocity" << endl;
  }

  int allnull, allthere;
  checkBetas(betan, betanp1, allthere, allnull);

  int constant_viscosity = allnull;

  int use_dv_constant_mu = use_dv_constant_mu_def;
  for (int i=0; i<BL_SPACEDIM; i++) 
    use_dv_constant_mu = use_dv_constant_mu && visc_coef[Xvel+i] >= 0.0; 
  if(!use_dv_constant_mu && use_dv_constant_mu_def) {
    cout << "Diffusion::diffuse_velocity : must have velocity visc_coefs "
         << ">= 0.0 if use_dv_constant_mu == 1" << endl;
    exit(0);
  }

  if(constant_viscosity || use_dv_constant_mu)  {
    diffuse_velocity_constant_mu(dt, be_cn_theta, rho_half,
                                 delta_rhs);
  } else {
    diffuse_tensor_velocity(dt, be_cn_theta, rho_half,
                            delta_rhs, betan, betanp1);
  }
}

void Diffusion::diffuse_velocity_constant_mu(REAL dt, REAL be_cn_theta,
				 MultiFab* rho_half,
                                 MultiFab* delta_rhs)
{
  int finest_level = parent->finestLevel();

  // at this point, S_old has bndry at time N
  // S_new contains GRAD(SU)

  MultiFab &U_old = caller->get_old_data(State_Type);
  MultiFab &U_new = caller->get_new_data(State_Type);

  const BOX& domain = caller->Geom().Domain();
  const int* domlo = domain.loVect();
  const int* domhi = domain.hiVect();
  const REAL* dx = caller->Geom().CellSize();
  int ngrd = grids.length();

  REAL cur_time  = caller->get_state_data(State_Type).curTime();
  REAL prev_time = caller->get_state_data(State_Type).prevTime();

  FARRAYBOX xflux, yflux, zflux;

  // at this point in time we can only do decoupled scalar
  // so we loop over components
  int comp;
  for (comp = 0; comp < BL_SPACEDIM; comp++) {
    int sigma = Xvel + comp;
    // U_new now contains the inviscid update of U
    // this is part of the RHS for the viscous solve

    MultiFab Rhs(grids,1,0,Fab_allocate);
    MultiFab Soln(grids,1,1,Fab_allocate);

    REAL mf_norm = 0.;
    int i;
    int rho_flag = 1;

    { // set up Rhs
      REAL a = 0.0;
      REAL b = -(1.0-be_cn_theta)*visc_coef[sigma]*dt;

      ViscBndry visc_bndry;
      ABecLaplacian *visc_op = 
	getViscOp(sigma,a,b,prev_time,visc_bndry,rho_half,rho_flag);
      visc_op->maxOrder(max_order);

      // copy to single-component multifab
      // note: use Soln as a temporary here
      Soln.copy(U_old,sigma,0,1,1);
      visc_op->apply(Rhs,Soln);
      delete visc_op;

      // complete Rhs by adding body sources
      for (i = 0; i < ngrd; i++) {
	// scale inviscid part by volume
	U_new[i].mult(volume[i],grids[i],0,sigma,1);

	// multiply by density at time nph
	FARRAYBOX& Rh = (*rho_half)[i];
	U_new[i].mult(Rh,grids[i],0,sigma,1);

	// add to Rhs which contained mu/(2 dt) Lap(u)
	Rhs[i].plus(U_new[i],grids[i],sigma,0,1);
      }

      if (delta_rhs!=NULL) {
        mf_norm = 0.0;
        for (i = 0; i < ngrd; i++) {
          (*delta_rhs)[i].mult(dt,comp,1);
          (*delta_rhs)[i].mult(volume[i],grids[i],0,comp,1);
          Rhs[i].plus((*delta_rhs)[i],grids[i],comp,0,1);
          REAL gr_norm = Rhs[i].norm(0);
          mf_norm = Max(gr_norm,mf_norm);
        }
      }

//  Note: we have to add hoop stress explicitly because the hoop
//        stress which is added through the operator in getViscOp
//        is eliminated by setting a=0
#if (BL_SPACEDIM == 2) 
      if(sigma==Xvel && CoordSys::IsRZ()) {
	for (int k = 0; k < ngrd; k++) {
	  const BOX &bx = Rhs.box(k);
	  BOX sbx       = grow(U_old.box(k),U_old.nGrow());
	  int rlen      = bx.length(0);
	  Array<REAL> rcen;
	  rcen.resize(rlen);
	  parent->Geom(level).GetCellLoc(rcen, bx, 0);
	  const int *lo       = bx.loVect();
	  const int *hi       = bx.hiVect();
	  const int *slo      = sbx.loVect();
	  const int *shi      = sbx.hiVect();
	  REAL *rhs           = Rhs[k].dataPtr();
	  REAL *sdat          = U_old[k].dataPtr(sigma);
	  const REAL *rcendat = rcen.dataPtr();
	  const REAL coeff    = (1.0-be_cn_theta)*visc_coef[sigma]*dt;
	  const REAL *voli    = volume[k].dataPtr();
	  BOX vbox            = grow(volume.box(k),volume.nGrow());
	  const int *vlo      = vbox.loVect();
	  const int *vhi      = vbox.hiVect();
	  FORT_HOOPRHS(rhs, ARLIM(lo), ARLIM(hi), 
                       sdat, ARLIM(slo), ARLIM(shi),
		       rcendat, &coeff, 
                       voli, ARLIM(vlo), ARLIM(vhi));
	}
      }
#endif

      for (i = 0; i < ngrd; i++) {
	REAL gr_norm = Rhs[i].norm(0);
	mf_norm = Max(gr_norm,mf_norm);
      }
    }

    // compute guess of solution
    for (i = 0; i < grids.length(); i++) {
      if (level == 0) {
	Soln[i].copy(U_old[i],sigma,0,1);
      }
      else {
	// coarse grid data exists at this time
	// use interpolated crse grid data for guess
	caller->FillCoarsePatch(Soln[i],0,cur_time,State_Type,sigma,1);
      }
    }

    // copy guess into U_new ... needed for construction of
    // viscous operator.  We first set all of U_new to a bogus values 
    // so that its ghost cell values are not undefined.  These values
    // will be filled with better values later.
    int n_comp = 1;
    int n_ghost = 1;
    U_new.setVal(1.e30,sigma,n_comp,n_ghost);
    n_ghost = 0;
    U_new.copy(Soln,0,sigma,n_comp,0);

    // construct viscous operator with bndry data at time N+1
    REAL a = 1.0;
    REAL b = be_cn_theta*dt*visc_coef[sigma];
       
    ViscBndry visc_bndry;
    rho_flag = 1 ;
    REAL rhsscale = 1.0;
    ABecLaplacian *visc_op = getViscOp(sigma,a,b,cur_time,visc_bndry,
				       rho_half,rho_flag,&rhsscale);
    visc_op->maxOrder(max_order);
    Rhs.mult(rhsscale,0,1);

    // construct solver and call it
    REAL S_tol = visc_tol;
    REAL S_tol_abs = visc_abs_tol*mf_norm;
    if(use_cg_solve) {
      CGSolver cg(*visc_op,use_mg_precond_flag);
      cg.solve(Soln,Rhs,S_tol,S_tol_abs);
    } else {
      MultiGrid mg(*visc_op);
      mg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }

    int visc_op_lev = 0;
    visc_op->applyBC(Soln,visc_op_lev);

    delete visc_op;

    // copy into state variable at new time
    n_comp = 1;
    n_ghost = 1;
    U_new.copy(Soln,0,sigma,n_comp,n_ghost);

    // modify diffusive fluxes here
    if (do_reflux) {
      for (i = 0; i < ngrd; i++) {
	const BOX& grd = grids[i];
	const int* lo = grd.loVect();
	const int* hi = grd.hiVect();

	const int* ulo = U_old[i].loVect();
	const int* uhi = U_old[i].hiVect();

	BOX xflux_bx(grd);
	xflux_bx.surroundingNodes(0);
	xflux.resize(xflux_bx,1);
        DEF_LIMITS(xflux,xflux_dat,xflo,xfhi);

	BOX yflux_bx(grd);
	yflux_bx.surroundingNodes(1);
	yflux.resize(yflux_bx,1);
        DEF_LIMITS(yflux,yflux_dat,yflo,yfhi);

	FARRAYBOX& xarea = area[0][i];
	FARRAYBOX& yarea = area[1][i];

	DEF_CLIMITS(xarea,xarea_dat,axlo,axhi);
	DEF_CLIMITS(yarea,yarea_dat,aylo,ayhi);

	REAL mult = -visc_coef[sigma];

#if (BL_SPACEDIM == 2)
        FORT_VISCFLUX (U_old[i].dataPtr(sigma), 
  		       U_new[i].dataPtr(sigma), 
                       ARLIM(ulo), ARLIM(uhi),
		       lo,hi,
                       xflux_dat,ARLIM(xflo),ARLIM(xfhi),
                       yflux_dat,ARLIM(yflo),ARLIM(yfhi),
                       xarea_dat,ARLIM(axlo),ARLIM(axhi),
                       yarea_dat,ARLIM(aylo),ARLIM(ayhi),
		       dx,&mult,&be_cn_theta);
#endif
#if (BL_SPACEDIM == 3)

	BOX zflux_bx(grd);
	zflux_bx.surroundingNodes(2);
	zflux.resize(zflux_bx,1);
        DEF_LIMITS(zflux,zflux_dat,zflo,zfhi);

	FARRAYBOX& zarea = area[2][i];
	DEF_CLIMITS(zarea,zarea_dat,azlo,azhi);

        FORT_VISCFLUX (U_old[i].dataPtr(sigma), 
		       U_new[i].dataPtr(sigma), 
                       ARLIM(ulo), ARLIM(uhi),
		       lo,hi,
                       xflux_dat,ARLIM(xflo),ARLIM(xfhi),
                       yflux_dat,ARLIM(yflo),ARLIM(yfhi),
                       zflux_dat,ARLIM(zflo),ARLIM(zfhi),
                       xarea_dat,ARLIM(axlo),ARLIM(axhi),
                       yarea_dat,ARLIM(aylo),ARLIM(ayhi),
                       zarea_dat,ARLIM(azlo),ARLIM(azhi),
		       dx,&mult,&be_cn_theta);
#endif

	// NOTE: fluxes expected to be in extensive form
	if (level < finest_level) {
	  FluxRegister& fr = *finer->viscflux_reg;
	  fr.CrseInit(xflux,xflux_bx,0,0,sigma,1,-dt);
	  fr.CrseInit(yflux,yflux_bx,1,0,sigma,1,-dt);
#if (BL_SPACEDIM == 3)
	  fr.CrseInit(zflux,zflux_bx,2,0,sigma,1,-dt);
#endif
	}
	if (level > 0) {
	  viscflux_reg->FineAdd(xflux,0,i,0,sigma,1,dt);
	  viscflux_reg->FineAdd(yflux,1,i,0,sigma,1,dt);
#if (BL_SPACEDIM == 3)
	  viscflux_reg->FineAdd(zflux,2,i,0,sigma,1,dt);
#endif
	}
      }
    }
  }
}

void Diffusion::diffuse_tensor_velocity(REAL dt, REAL be_cn_theta,
				 MultiFab* rho_half,
                                 MultiFab* delta_rhs,
                                 MultiFab** betan, 
                                 MultiFab** betanp1)

{
#if (BL_SPACEDIM==3) 
    cout << "Diffusion::diffuse_tensor_velocity : " <<
            "not yet implemented for 3-D" << endl;
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0" <<
            " and rerun" << endl;
    exit(0);
#elif !defined(USE_TENSOR)
    cout << "Diffusion::diffuse_tensor_velocity :  " <<
    cout << "USE_TENSOR must be defined at compile time" << endl;
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0" <<
            " and rerun" << endl;
    exit(0);
#else

  int finest_level = parent->finestLevel();

  // at this point, S_old has bndry at time N
  // S_new contains GRAD(SU)

  MultiFab &U_old = caller->get_old_data(State_Type);
  MultiFab &U_new = caller->get_new_data(State_Type);

  const BOX& domain = caller->Geom().Domain();
  const int* domlo = domain.loVect();
  const int* domhi = domain.hiVect();
  const REAL* dx = caller->Geom().CellSize();
  int ngrd = grids.length();

  REAL cur_time  = caller->get_state_data(State_Type).curTime();
  REAL prev_time = caller->get_state_data(State_Type).prevTime();

  int comp;
  // U_new now contains the inviscid update of U
  // this is part of the RHS for the viscous solve

  MultiFab Rhs(grids,BL_SPACEDIM,0,Fab_allocate);
  Rhs.setVal(0.0);

  REAL mf_norm = 0.;
  int i, dim;

  MultiFab **tensorflux_old;
  { // set up Rhs
    int soln_old_grow = 1;
    MultiFab Soln_old(grids,BL_SPACEDIM,soln_old_grow,Fab_allocate);
    REAL a = 0.0;
    REAL b = -(1.0-be_cn_theta)*dt;

    ViscBndry2D visc_bndry;
    DivVis *tensor_op = 
      getTensorOp(a,b,prev_time,visc_bndry,rho_half,betan);
    tensor_op->maxOrder(tensor_max_order);

    // copy to single-component multifab
    // note: use Soln as a temporary here
    Soln_old.copy(U_old,Xvel,0,BL_SPACEDIM,soln_old_grow);
    tensor_op->apply(Rhs,Soln_old);
    if (do_reflux && (level<finest_level || level>0)) {
      allocFluxBoxesLevel(tensorflux_old,0,BL_SPACEDIM);
      tensor_op->compFlux(*(tensorflux_old[0]),*(tensorflux_old[1]),
#if(BL_SPACEDIM==3)
                          *(tensorflux_old[2]),
#endif
                          Soln_old);
      for (dim=0; dim<BL_SPACEDIM; dim++)
        tensorflux_old[dim]->mult(-(1.0-be_cn_theta),0);      
    }
    delete tensor_op;

    for (comp = 0; comp < BL_SPACEDIM; comp++) {
      int sigma = Xvel + comp;
      // complete Rhs by adding body sources
      for (i = 0; i < ngrd; i++) {
	// scale inviscid part by volume
	U_new[i].mult(volume[i],grids[i],0,sigma,1);

	// multiply by density at time nph
	FARRAYBOX& Rh = (*rho_half)[i];
	U_new[i].mult(Rh,grids[i],0,sigma,1);

	// add to Rhs which contained operator applied to U_old
	Rhs[i].plus(U_new[i],grids[i],sigma,comp,1);
      }

      if (delta_rhs!=NULL) {
        for (i = 0; i < ngrd; i++) {
          (*delta_rhs)[i].mult(dt,comp,1);
          (*delta_rhs)[i].mult(volume[i],grids[i],0,comp,1);
          Rhs[i].plus((*delta_rhs)[i],grids[i],comp,comp,1);
        }
      }
    }

#if (BL_SPACEDIM == 2) 
    if(CoordSys::IsRZ()) {
      int fort_xvel_comp = Xvel+1;
      for (int k = 0; k < ngrd; k++) {
	const BOX &bx = Rhs.box(k);
	BOX sbx       = grow(U_old.box(k),U_old.nGrow());
	int rlen      = bx.length(0);
	Array<REAL> rcen;
	rcen.resize(rlen);
	parent->Geom(level).GetCellLoc(rcen, bx, 0);
	const int *lo       = bx.loVect();
	const int *hi       = bx.hiVect();
	const int *slo      = sbx.loVect();
	const int *shi      = sbx.hiVect();
	REAL *rhs           = Rhs[k].dataPtr();
	REAL *sdat          = U_old[k].dataPtr(Xvel);
	const REAL *rcendat = rcen.dataPtr();
	const REAL coeff    = (1.0-be_cn_theta)*dt;
	const REAL *voli    = volume[k].dataPtr();
	BOX vbox            = grow(volume.box(k),volume.nGrow());
	const int *vlo      = vbox.loVect();
	const int *vhi      = vbox.hiVect();
    
        FARRAYBOX& betax = (*betanp1[0])[k];
        DEF_CLIMITS(betax,betax_dat,betax_lo,betax_hi);

        FARRAYBOX& betay = (*betanp1[1])[k];
        DEF_CLIMITS(betay,betay_dat,betay_lo,betay_hi);

 	FORT_TENSOR_HOOPRHS(&fort_xvel_comp, rhs, ARLIM(lo), ARLIM(hi), 
                       sdat, ARLIM(slo), ARLIM(shi),
		       rcendat, &coeff, 
                       voli, ARLIM(vlo), ARLIM(vhi),
                       betax_dat,ARLIM(betax_lo),ARLIM(betax_hi),
                       betay_dat,ARLIM(betay_lo),ARLIM(betay_hi));
      }
    }
#endif

    for (i = 0; i < ngrd; i++) {
      REAL gr_norm = Rhs[i].norm(0);
      mf_norm = Max(gr_norm,mf_norm);
    }
  }

// I am using a ghost cell in Soln even though Bill does not
    int soln_grow = 1;
    MultiFab Soln(grids,BL_SPACEDIM,soln_grow,Fab_allocate);

    Soln.setVal(0.0);
    // compute guess of solution
    for (i = 0; i < grids.length(); i++) {
      if (level == 0) {
	Soln[i].copy(U_old[i],Xvel,0,BL_SPACEDIM);
      }
      else {
	// coarse grid data exists at this time
	// use interpolated crse grid data for guess
	caller->FillCoarsePatch(Soln[i],0,cur_time,State_Type,Xvel,BL_SPACEDIM);
      }
    }

    // copy guess into U_new ... needed for construction of
    // viscous operator.  We first set all of U_new to a bogus values 
    // so that its ghost cell values are not undefined.  These values
    // will be filled with better values later.
    int n_comp = BL_SPACEDIM;
    int n_ghost = 1;
    U_new.setVal(1.e30,Xvel,n_comp,n_ghost);
    n_ghost = 0;
    U_new.copy(Soln,0,Xvel,n_comp,0);

    // construct viscous operator with bndry data at time N+1
    REAL a = 1.0;
    REAL b = be_cn_theta*dt;
       
    ViscBndry2D visc_bndry;
    DivVis *tensor_op = getTensorOp(a,b,cur_time,visc_bndry,
				       rho_half,betanp1);
    tensor_op->maxOrder(tensor_max_order);

    // construct solver and call it
    REAL S_tol = visc_tol;
    REAL S_tol_abs = visc_abs_tol*mf_norm;
    if(use_tensor_cg_solve) {
      int use_mg_pre = 0;
      MCCGSolver cg(*tensor_op,use_mg_pre);
      cg.solve(Soln,Rhs,S_tol,S_tol_abs);
    } else {
      MCMultiGrid mg(*tensor_op);
      mg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }

#if 1
    int visc_op_lev = 0;
    tensor_op->applyBC(Soln,visc_op_lev); // this may not be needed with
                                        // Bills stuff
#endif

    // copy into state variable at new time
    n_ghost = soln_grow;
    U_new.copy(Soln,0,Xvel,n_comp,n_ghost);

    // modify diffusive fluxes here
    if (do_reflux && (level<finest_level || level>0)) {
      MultiFab **tensorflux;
      allocFluxBoxesLevel(tensorflux,0,BL_SPACEDIM);
      tensor_op->compFlux(*(tensorflux[0]),*(tensorflux[1]),
#if(BL_SPACEDIM==3)
                          *(tensorflux[2]),
#endif
                          Soln);
      for (dim=0; dim<BL_SPACEDIM; dim++) {
        tensorflux[dim]->mult(-be_cn_theta,0);
        tensorflux[dim]->plus(*(tensorflux_old[dim]),0,BL_SPACEDIM,0);
      }       
      removeFluxBoxesLevel(tensorflux_old);
    
      FARRAYBOX xflux, yflux, zflux;

      for (int sigma = Xvel; sigma < BL_SPACEDIM+Xvel; sigma++) {

        for (i = 0; i < ngrd; i++) {
	  const BOX& grd = grids[i];
	  const int* lo = grd.loVect();
	  const int* hi = grd.hiVect();

	  BOX xflux_bx(grd);
	  xflux_bx.surroundingNodes(0);
	  xflux.resize(xflux_bx,1);
          xflux.copy((*(tensorflux[0]))[i],sigma,0,1);

	  BOX yflux_bx(grd);
	  yflux_bx.surroundingNodes(1);
	  yflux.resize(yflux_bx,1);
          yflux.copy((*(tensorflux[1]))[i],sigma,0,1);

#if (BL_SPACEDIM == 3)
	  BOX zflux_bx(grd);
	  zflux_bx.surroundingNodes(2);
	  zflux.resize(zflux_bx,1);
          zflux.copy((*(tensorflux[2]))[i],sigma,0,1);
#endif

	  if (level < finest_level) {
	    FluxRegister& fr = *finer->viscflux_reg;
	    fr.CrseInit(xflux,xflux_bx,0,0,sigma,1,-dt);
	    fr.CrseInit(yflux,yflux_bx,1,0,sigma,1,-dt);
#if (BL_SPACEDIM == 3)
	    fr.CrseInit(zflux,zflux_bx,2,0,sigma,1,-dt);
#endif
	  }
	  if (level > 0) {
	    viscflux_reg->FineAdd(xflux,0,i,0,sigma,1,dt);
	    viscflux_reg->FineAdd(yflux,1,i,0,sigma,1,dt);
#if (BL_SPACEDIM == 3)
	    viscflux_reg->FineAdd(zflux,2,i,0,sigma,1,dt);
#endif
	  }
        }
      }
      removeFluxBoxesLevel(tensorflux);
    }
    delete tensor_op;
#endif
    
}

void Diffusion::diffuse_Vsync(MultiFab *Vsync, REAL dt,
			      REAL be_cn_theta,
			      MultiFab* rho_half, int rho_flag,
                              MultiFab** beta)

{
  if (verbose) {
    cout << "Diffusion::diffuse_Vsync" << endl;
  }

  int allnull, allthere;
  checkBeta(beta, allthere, allnull);

  int constant_viscosity = allnull;

  int use_dv_constant_mu = use_dv_constant_mu_def;
  for (int i=0; i<BL_SPACEDIM; i++) 
    use_dv_constant_mu = use_dv_constant_mu && visc_coef[Xvel+i] >= 0.0; 
  if(!use_dv_constant_mu && use_dv_constant_mu_def) {
    cout << "Diffusion::diffuse_Vsync : must have velocity visc_coefs "
         << ">= 0.0 if use_dv_constant_mu == 1" << endl;
    exit(0);
  }

  if(constant_viscosity || use_dv_constant_mu)  {
    diffuse_Vsync_constant_mu(Vsync, dt, be_cn_theta, rho_half,
                              rho_flag);
  } else {
    diffuse_tensor_Vsync(Vsync, dt, be_cn_theta, rho_half,
                              rho_flag, beta);
  }

// applyBC has put "incorrect" values in the ghost cells
// outside external Dirichlet boundaries. Reset these to zero
// so that syncproject and conservative interpolation works correctly.

   BOX domain = caller->Geom().Domain();
   domain.grow(1);

   for (int n = Xvel; n < Xvel+BL_SPACEDIM; n++) {
     const BCRec& velbc = caller->get_desc_lst()[State_Type].getBC(n);

     for (int k = 0; k < BL_SPACEDIM; k++) {

       if (velbc.hi(k) == EXT_DIR) {
         INTVECT bigend   = domain.bigEnd();
         INTVECT smallend = domain.smallEnd();
         int hi = domain.bigEnd(k);

         smallend.setVal(k,hi);
         BOX top_strip(smallend,bigend,IntVect::TheCellVector());
         Vsync->setVal(0.0,top_strip,n-Xvel,1,1);
       }

       if (velbc.lo(k) == EXT_DIR) {
         INTVECT smallend = domain.smallEnd();
         INTVECT bigend   = domain.bigEnd();
         int lo = domain.smallEnd(k);

         bigend.setVal(k,lo);
         BOX bottom_strip(smallend,bigend,IntVect::TheCellVector());
         Vsync->setVal(0.0,bottom_strip,n-Xvel,1,1);
       }

     }
  }
}

void Diffusion::diffuse_Vsync_constant_mu(MultiFab *Vsync, REAL dt,
			      REAL be_cn_theta,
			      MultiFab* rho_half, int rho_flag)
{
  int i;

  if (verbose)
    cout << "Diffusion::diffuse_Vsync" << endl;

  int ngrds = grids.length();
  const REAL* dx = caller->Geom().CellSize();

  // at this point in time we can only do decoupled scalar
  // so we loop over components
  for (int comp = 0; comp < BL_SPACEDIM; comp++) {

    MultiFab Soln(grids,1,1,Fab_allocate);
    Soln.setVal(0.);

    MultiFab Rhs(grids,1,0,Fab_allocate);
    Rhs.copy(*Vsync,comp,0,1);

    REAL r_norm = Rhs[0].norm(0);
    cout << "Original max of Vsync " << r_norm << endl;

    //  Compute norm of RHS after multiplication by volume and density
    REAL mf_norm = 0.0;
    for (i = 0; i < ngrds; i++) {
      Rhs[i].mult(volume[i]); 
      Rhs[i].mult((*rho_half)[i]); 
      REAL gr_norm = Rhs[i].norm(0);
      mf_norm = Max(gr_norm,mf_norm);
    }

    //  SET UP COEFFICIENTS FOR VISCOUS SOLVER
    REAL a = 1.0;

    REAL b = be_cn_theta*dt*visc_coef[comp];

    REAL rhsscale = 1.0;
    ABecLaplacian * visc_op = getViscOp(comp,a,b,rho_half,rho_flag,&rhsscale);
    visc_op->maxOrder(max_order);
    Rhs.mult(rhsscale,0,1);

      // construct solver and call it
    REAL S_tol = visc_tol;
    REAL S_tol_abs = visc_abs_tol*mf_norm;
    if(use_cg_solve) {
      CGSolver cg(*visc_op,use_mg_precond_flag);
      cg.solve(Soln,Rhs,S_tol,S_tol_abs);
    } else {
      MultiGrid mg(*visc_op);
      mg.solve(Soln,Rhs,S_tol,S_tol_abs);
    }

    int visc_op_lev = 0;
    visc_op->applyBC(Soln,visc_op_lev);

    Vsync->copy(Soln,0,comp,1,1);

    REAL s_norm = Soln[0].norm(0);
    cout << "Final max of Vsync " << s_norm << endl;

    delete visc_op;

    FARRAYBOX xflux;
    FARRAYBOX yflux;

    if (level > 0) {

      for (i = 0; i < ngrds; i++) {
        const BOX& grd = grids[i];
        const int* lo = grd.loVect();
        const int* hi = grd.hiVect();

        FARRAYBOX& u_sync = (*Vsync)[i];

        const int* ulo = u_sync.loVect();
        const int* uhi = u_sync.hiVect();

        BOX xflux_bx(grd);
        xflux_bx.surroundingNodes(0);
        xflux.resize(xflux_bx,1);
        DEF_LIMITS(xflux,xflux_dat,xflux_lo,xflux_hi);

        BOX yflux_bx(grd);
        yflux_bx.surroundingNodes(1);
        yflux.resize(yflux_bx,1);
        DEF_LIMITS(yflux,yflux_dat,yflux_lo,yflux_hi);

        FARRAYBOX& xarea = area[0][i];
        FARRAYBOX& yarea = area[1][i];

        DEF_CLIMITS(xarea,xarea_dat,xarea_lo,xarea_hi);
        DEF_CLIMITS(yarea,yarea_dat,yarea_lo,yarea_hi);

	//      NOTE: the extra factor of dt comes from the fact that Vsync looks
	//            like dV/dt, not just an increment to V

        REAL mult = -be_cn_theta*dt*dt*visc_coef[comp];

#if (BL_SPACEDIM == 2)
        FORT_VISCSYNCFLUX (u_sync.dataPtr(comp), ARLIM(ulo), ARLIM(uhi),
                           lo,hi,
                           xflux_dat,ARLIM(xflux_lo),ARLIM(xflux_hi),
                           yflux_dat,ARLIM(yflux_lo),ARLIM(yflux_hi),
                           xarea_dat,ARLIM(xarea_lo),ARLIM(xarea_hi),
                           yarea_dat,ARLIM(yarea_lo),ARLIM(yarea_hi),
                           dx,&mult);
#endif
#if (BL_SPACEDIM == 3)

        FARRAYBOX zflux;
        BOX zflux_bx(grd);
        zflux_bx.surroundingNodes(2);
        zflux.resize(zflux_bx,1);
        DEF_LIMITS(zflux,zflux_dat,zflux_lo,zflux_hi);

        FARRAYBOX& zarea = area[2][i];
        DEF_CLIMITS(zarea,zarea_dat,zarea_lo,zarea_hi);

        FORT_VISCSYNCFLUX (u_sync.dataPtr(comp), ARLIM(ulo), ARLIM(uhi),
                           lo,hi,
                           xflux_dat,ARLIM(xflux_lo),ARLIM(xflux_hi),
                           yflux_dat,ARLIM(yflux_lo),ARLIM(yflux_hi),
                           zflux_dat,ARLIM(zflux_lo),ARLIM(zflux_hi),
                           xarea_dat,ARLIM(xarea_lo),ARLIM(xarea_hi),
                           yarea_dat,ARLIM(yarea_lo),ARLIM(yarea_hi),
                           zarea_dat,ARLIM(zarea_lo),ARLIM(zarea_hi),
                           dx,&mult);
#endif

        REAL one = 1.0;
        viscflux_reg->FineAdd(xflux,0,i,0,comp,1,one);
        viscflux_reg->FineAdd(yflux,1,i,0,comp,1,one);
#if (BL_SPACEDIM == 3)
        viscflux_reg->FineAdd(zflux,2,i,0,comp,1,one);
#endif
      }
    }
  }
}

void Diffusion::diffuse_tensor_Vsync(MultiFab *Vsync, REAL dt,
			      REAL be_cn_theta,
			      MultiFab* rho_half, int rho_flag,
                              MultiFab** beta)
{
#if (BL_SPACEDIM==3) 
    cout << "Diffusion::diffuse_tensor_Vsync : " <<
            "not yet implemented for 3-D" << endl;
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0" <<
            " and rerun" << endl;
    exit(0);
#elif !defined(USE_TENSOR)
    cout << "Diffusion::diffuse_tensor_Vsync :  " <<
    cout << "USE_TENSOR must be defined at compile time" << endl;
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0" <<
            " and rerun" << endl;
    exit(0);
#else

  int finest_level = parent->finestLevel();

  int i;

  int ngrds = grids.length();
  const REAL* dx = caller->Geom().CellSize();

  MultiFab Soln(grids,BL_SPACEDIM,1,Fab_allocate);
  Soln.setVal(0.);

  MultiFab Rhs(grids,BL_SPACEDIM,0,Fab_allocate);
  Rhs.copy(*Vsync,0,0,BL_SPACEDIM);

  REAL r_norm = 0.0;
  for (i = 0; i < ngrds; i++)
    r_norm = Max(r_norm,Rhs[i].norm(0));
  cout << "Original max of Vsync " << r_norm << endl;

  //  Compute norm of RHS after multiplication by volume and density
  REAL mf_norm = 0.0;
  for (int comp = 0; comp < BL_SPACEDIM; comp++) {
    for (i = 0; i < ngrds; i++) {
      Rhs[i].mult(volume[i],0,comp,1); 
      Rhs[i].mult((*rho_half)[i],0,comp,1); 
      REAL gr_norm = Rhs[i].norm(0);
      mf_norm = Max(gr_norm,mf_norm);
    }
  }

  //  SET UP COEFFICIENTS FOR VISCOUS SOLVER
  REAL a = 1.0;

  REAL b = be_cn_theta*dt;

  DivVis * tensor_op = getTensorOp(a,b,rho_half,beta);
  tensor_op->maxOrder(tensor_max_order);

  // construct solver and call it
  REAL S_tol = visc_tol;
  REAL S_tol_abs = visc_abs_tol*mf_norm;
  if(use_tensor_cg_solve) {
    MCCGSolver cg(*tensor_op,use_mg_precond_flag);
    cg.solve(Soln,Rhs,S_tol,S_tol_abs);
  } else {
    MCMultiGrid mg(*tensor_op);
    mg.solve(Soln,Rhs,S_tol,S_tol_abs);
  }

  int visc_op_lev = 0;
  tensor_op->applyBC(Soln,visc_op_lev); 

  Vsync->copy(Soln,0,0,BL_SPACEDIM,1);

  REAL s_norm = 0.0;
  for (i = 0; i < ngrds; i++)
    s_norm = Max(s_norm,Soln[i].norm(0));
  cout << "Final max of Vsync " << s_norm << endl;

  FARRAYBOX xflux, yflux, zflux;

  if (level > 0) {

    MultiFab **tensorflux;
    allocFluxBoxesLevel(tensorflux,0,BL_SPACEDIM);
    tensor_op->compFlux(*(tensorflux[0]),*(tensorflux[1]),
#if(BL_SPACEDIM==3)
                        *(tensorflux[2]),
#endif
                        Soln);
//  NOTE: the extra factor of dt comes from the fact that Vsync looks
//        like dV/dt, not just an increment to V

    for (int dim=0; dim<BL_SPACEDIM; dim++)
      tensorflux[0]->mult(-be_cn_theta*dt*dt,0);

    for (int sigma = Xvel; sigma < BL_SPACEDIM+Xvel; sigma++) {
      for (i = 0; i < ngrds; i++) {
        const BOX& grd = grids[i];
        const int* lo = grd.loVect();
        const int* hi = grd.hiVect();

        BOX xflux_bx(grd);
        xflux_bx.surroundingNodes(0);
        xflux.resize(xflux_bx,1);
        xflux.copy((*(tensorflux[0]))[i],sigma,0,1);

        BOX yflux_bx(grd);
        yflux_bx.surroundingNodes(1);
        yflux.resize(yflux_bx,1);
        yflux.copy((*(tensorflux[1]))[i],sigma,0,1); 

#if (BL_SPACEDIM == 3)
	BOX zflux_bx(grd);
	zflux_bx.surroundingNodes(2);
	zflux.resize(zflux_bx,1);
        zflux.copy((*(tensorflux[2]))[i],sigma,0,1);
#endif

        REAL one = 1.0;
        viscflux_reg->FineAdd(xflux,0,i,0,sigma,1,one);
        viscflux_reg->FineAdd(yflux,1,i,0,sigma,1,one);
#if (BL_SPACEDIM == 3)
        viscflux_reg->FineAdd(zflux,2,i,0,sigma,1,one);
#endif
      }
    }
    removeFluxBoxesLevel(tensorflux);
  }
  delete tensor_op;
#endif
}

void Diffusion::diffuse_Ssync(MultiFab *Ssync, int sigma, REAL dt,
			      REAL be_cn_theta,
			      MultiFab* rho_half, int rho_flag,
                              int do_viscsyncflux, 
                              MultiFab** beta,
                              MultiFab* alpha)


{
  if (verbose)
    cout << "Diffusion::diffuse_Ssync for scalar " << sigma << endl;

  int allnull, allthere;
  checkBeta(beta, allthere, allnull);

  int ngrds = grids.length();
  const REAL* dx = caller->Geom().CellSize();

  MultiFab Soln(grids,1,1,Fab_allocate);
  Soln.setVal(0.);

  MultiFab Rhs(grids,1,0,Fab_allocate);
  Rhs.copy(*Ssync,sigma,0,1);

  REAL r_norm = 0.0;
  int i;
  for (i = 0; i < ngrds; i++)
    r_norm = Max(r_norm,Rhs[i].norm(0));
  cout << "Original max of Ssync " << r_norm << endl;

  //  Compute norm of RHS
  REAL mf_norm = 0.0;
  for (i = 0; i < ngrds; i++) {
    Rhs[i].mult(volume[i]); 
    if (rho_flag==1) 
      Rhs[i].mult((*rho_half)[i]);
    REAL gr_norm = Rhs[i].norm(0);
    mf_norm = Max(gr_norm,mf_norm);
  }

  //  SET UP COEFFICIENTS FOR VISCOUS SOLVER
  REAL a = 1.0;
  REAL b;
  if (allnull) {
    b = be_cn_theta*dt*visc_coef[BL_SPACEDIM+sigma];
  } else {
    b = be_cn_theta*dt;
  }

  REAL rhsscale = 1.0;
  ABecLaplacian * visc_op = getViscOp(BL_SPACEDIM+sigma,a,b,
                                      rho_half,rho_flag,&rhsscale,beta,alpha);
  visc_op->maxOrder(max_order);
  Rhs.mult(rhsscale,0,1);

      // construct solver and call it
  REAL S_tol = visc_tol;
  REAL S_tol_abs = visc_abs_tol*mf_norm;
  if(use_cg_solve) {
    CGSolver cg(*visc_op,use_mg_precond_flag);
    cg.solve(Soln,Rhs,S_tol,S_tol_abs);
  } else {
    MultiGrid mg(*visc_op);
    mg.solve(Soln,Rhs,S_tol,S_tol_abs);
  }

  int visc_op_lev = 0;
  visc_op->applyBC(Soln,visc_op_lev);

  Ssync->copy(Soln,0,sigma,1,1);

  REAL s_norm = 0.0;
  for (i = 0; i < ngrds; i++)
    s_norm = Max(s_norm,Soln[i].norm(0));
  cout << "Final max of Ssync " << s_norm << endl;

  delete visc_op;

  FARRAYBOX xflux;
  FARRAYBOX yflux;

  if (level > 0 && do_viscsyncflux == 1) {
 
    for (i = 0; i < ngrds; i++) {
      const BOX& grd = grids[i];
      const int* lo = grd.loVect();
      const int* hi = grd.hiVect();

      FARRAYBOX& s_sync = (*Ssync)[i];

      const int* slo = s_sync.loVect();
      const int* shi = s_sync.hiVect();

      BOX xflux_bx(grd);
      xflux_bx.surroundingNodes(0);
      xflux.resize(xflux_bx,1);
      DEF_LIMITS(xflux,xflux_dat,xflux_lo,xflux_hi);

      BOX yflux_bx(grd);
      yflux_bx.surroundingNodes(1);
      yflux.resize(yflux_bx,1);
      DEF_LIMITS(yflux,yflux_dat,yflux_lo,yflux_hi);

      FARRAYBOX& xarea = area[0][i];
      FARRAYBOX& yarea = area[1][i];

      DEF_CLIMITS(xarea,xarea_dat,xarea_lo,xarea_hi);
      DEF_CLIMITS(yarea,yarea_dat,yarea_lo,yarea_hi);

      REAL mult = -b;

#if (BL_SPACEDIM == 2)
      FORT_VISCSYNCFLUX (s_sync.dataPtr(sigma), ARLIM(slo), ARLIM(shi),
                         lo,hi,
                         xflux_dat,ARLIM(xflux_lo),ARLIM(xflux_hi),
                         yflux_dat,ARLIM(yflux_lo),ARLIM(yflux_hi),
                         xarea_dat,ARLIM(xarea_lo),ARLIM(xarea_hi),
                         yarea_dat,ARLIM(yarea_lo),ARLIM(yarea_hi),
                         dx,&mult);
#endif
#if (BL_SPACEDIM == 3)

      FARRAYBOX zflux;
      BOX zflux_bx(grd);
      zflux_bx.surroundingNodes(2);
      zflux.resize(zflux_bx,1);
      DEF_LIMITS(zflux,zflux_dat,zflux_lo,zflux_hi);

      FARRAYBOX& zarea = area[2][i];
      DEF_CLIMITS(zarea,zarea_dat,zarea_lo,zarea_hi);

      FORT_VISCSYNCFLUX (s_sync.dataPtr(sigma), ARLIM(slo), ARLIM(shi),
                         lo,hi,
                         xflux_dat,ARLIM(xflux_lo),ARLIM(xflux_hi),
                         yflux_dat,ARLIM(yflux_lo),ARLIM(yflux_hi),
                         zflux_dat,ARLIM(zflux_lo),ARLIM(zflux_hi),
                         xarea_dat,ARLIM(xarea_lo),ARLIM(xarea_hi),
                         yarea_dat,ARLIM(yarea_lo),ARLIM(yarea_hi),
                         zarea_dat,ARLIM(zarea_lo),ARLIM(zarea_hi),
                         dx,&mult);
#endif

      if (allthere) {
        xflux.mult((*beta[0])[i]);
        yflux.mult((*beta[1])[i]);
#if (BL_SPACEDIM == 3)
        zflux.mult((*beta[2])[i]);
#endif
      }

      REAL one = 1.0;
      viscflux_reg->FineAdd(xflux,0,i,0,BL_SPACEDIM+sigma,1,one);
      viscflux_reg->FineAdd(yflux,1,i,0,BL_SPACEDIM+sigma,1,one);
#if (BL_SPACEDIM == 3)
      viscflux_reg->FineAdd(zflux,2,i,0,BL_SPACEDIM+sigma,1,one);
#endif
    }
  }

  if(rho_flag==2) {
// we just solved for S -- what we want is rho*S
    MultiFab &S_new = caller->get_new_data(State_Type);
    MultiFab Rho_new(grids,1,1,Fab_allocate);
    Rho_new.copy(S_new,Density,0,1,1);
    for (i = 0; i < ngrds; i++) { 
      BOX box = grids[i];
      box.grow(1); 
      (*Ssync)[i].mult(Rho_new[i],box,0,sigma,1);
    }
  }

// applyBC has put "incorrect" values in the ghost cells
// outside external Dirichlet boundaries. Reset these to zero
// so that conservative interpolation works correctly.

   const BCRec& scalbc =
     caller->get_desc_lst()[State_Type].getBC(sigma+BL_SPACEDIM);
   BOX domain = caller->Geom().Domain();
   domain.grow(1);
   for (int k = 0; k<BL_SPACEDIM; k++) {

     INTVECT bigend   = domain.bigEnd();
     INTVECT smallend = domain.smallEnd();
     
     int hi = domain.bigEnd(k);
     smallend.setVal(k,hi);
     BOX top_strip(smallend,bigend,IntVect::TheCellVector());
     if(scalbc.hi(k)==EXT_DIR) Ssync->setVal(0.0,top_strip,sigma,1,1);

     smallend = domain.smallEnd();
     bigend   = domain.bigEnd();
     int lo = domain.smallEnd(k);
     bigend.setVal(k,lo);
     BOX bottom_strip(smallend,bigend,IntVect::TheCellVector());
     if(scalbc.lo(k)==EXT_DIR) Ssync->setVal(0.0,bottom_strip,sigma,1,1);
  }
}

#if defined(USE_TENSOR)
DivVis* Diffusion::getTensorOp(REAL a, REAL b,
				    REAL time,
#if (BL_SPACEDIM==2)  
                                    ViscBndry2D& visc_bndry,
#else
                                    ViscBndry3D& visc_bndry,
#endif
				    MultiFab* rho_half,
                                    MultiFab** beta)
{
#if (BL_SPACEDIM==3) 
    cout << "Diffusion::getTensorOp :  " <<
            "not yet implemented for 3-D" << endl;
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0" <<
            " and rerun" << endl;
    exit(0);
#elif !defined(USE_TENSOR)
    cout << "Diffusion::getTensorOp :  " <<
    cout << "USE_TENSOR must be defined at compile time" << endl;
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0" <<
            " and rerun" << endl;
    exit(0);
#else

  int allthere;
  checkBeta(beta, allthere);

  const REAL* dx = caller->Geom().CellSize();

  getTensorBndryData(visc_bndry,time);

  DivVis *tensor_op = new DivVis(grids,visc_bndry,dx);
  tensor_op->maxOrder(tensor_max_order);

  int isrz = CoordSys::IsRZ();

  int nghost = 1; // just like bill

  MultiFab alpha;  // alpha should be the same size as volume
  alpha.define(grids,BL_SPACEDIM,nghost,Fab_allocate);
  alpha.setVal(0.0,nghost);
  int k;
  if (a!=0.0){
    for (k = 0; k < grids.length(); k++) {
      const BOX &bx = alpha.box(k);
      int rlen      = bx.length(0);
      Array<REAL> rcen;
      rcen.resize(rlen);
      parent->Geom(level).GetCellLoc(rcen, bx, 0);
      const int *lo       = bx.loVect();
      const int *hi       = bx.hiVect();
      REAL *alpha_dat     = alpha[k].dataPtr();
      BOX abx             = grow(bx,alpha.nGrow());
      const int *alo      = abx.loVect();
      const int *ahi      = abx.hiVect();
      const REAL *rcendat = rcen.dataPtr();
      const REAL *voli    = volume[k].dataPtr();
      BOX vbox            = grow(volume.box(k),volume.nGrow());
      const int *vlo      = vbox.loVect();
      const int *vhi      = vbox.hiVect();

      FARRAYBOX& Rh = (*rho_half)[k];
      DEF_LIMITS(Rh,rho_dat,rlo,rhi);

      FARRAYBOX& betax = (*beta[0])[k];
      DEF_CLIMITS(betax,betax_dat,betax_lo,betax_hi);

      FARRAYBOX& betay = (*beta[1])[k];
      DEF_CLIMITS(betay,betay_dat,betay_lo,betay_hi);

      FORT_SET_TENSOR_ALPHA(alpha_dat, ARLIM(alo), ARLIM(ahi),
                  lo, hi, rcendat, ARLIM(lo), ARLIM(hi), &b,
                  voli, ARLIM(vlo), ARLIM(vhi),
                  rho_dat,ARLIM(rlo),ARLIM(rhi),
                  betax_dat,ARLIM(betax_lo),ARLIM(betax_hi),
                  betay_dat,ARLIM(betay_lo),ARLIM(betay_hi),&isrz);

    }
  }

  tensor_op->setScalars(a,b);

  tensor_op->aCoefficients(alpha);
  for (int n = 0; n < BL_SPACEDIM; n++) {
    MultiFab bcoeffs(area[n].boxArray(),1,nghost,Fab_allocate);
    bcoeffs.setVal(0.0);
    bcoeffs.copy(area[n]);
    for (k = 0; k < grids.length(); k++) {
      bcoeffs[k].mult(dx[n]);
      bcoeffs[k].mult((*beta[n])[k]);
    }
    tensor_op->bCoefficients(bcoeffs,n);
  }

  return tensor_op;
#endif
}

DivVis* Diffusion::getTensorOp(REAL a, REAL b,
				    MultiFab* rho_half,
                                    MultiFab** beta)
{
#if (BL_SPACEDIM==3) 
    cout << "Diffusion::getTensorOp :  " <<
            "not yet implemented for 3-D" << endl;
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0" <<
            " and rerun" << endl;
    exit(0);
#elif !defined(USE_TENSOR)
    cout << "Diffusion::getTensorOp :  " <<
    cout << "USE_TENSOR must be defined at compile time" << endl;
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0" <<
            " and rerun" << endl;
    exit(0);
#else

  int allthere = beta != NULL;
  if(allthere) {
    for (int dim=0; dim<BL_SPACEDIM; dim++) {
      allthere = allthere && beta[dim] != NULL;
    }
  }
  if(!allthere) {
    cout << "Diffusion::getTensorOp : all betas must allocated" << endl;
    cout << "  all NULL or all non-NULL" << endl;
    exit(0);
  }

  const REAL* dx = caller->Geom().CellSize();
  const BOX& domain = caller->Geom().Domain();
  Array<BCRec> bcarray(2*BL_SPACEDIM);
  for (int idim=0; idim < BL_SPACEDIM; idim++) {
    bcarray[idim] = 
      caller->get_desc_lst()[State_Type].getBC(Xvel+idim);
    bcarray[idim+BL_SPACEDIM] = 
      BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
		   D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
  }

  IntVect ref_ratio;
  if (level > 0) {
    ref_ratio = parent->refRatio(level-1);
  } else {
    ref_ratio = IntVect::TheUnitVector();
  }

  ViscBndry2D bndry;
  bndry.define(grids,2*BL_SPACEDIM,caller->Geom());
  bndry.setHomogValues(bcarray, ref_ratio[0]);
  DivVis *tensor_op = new DivVis(grids,bndry,dx);
  tensor_op->maxOrder(tensor_max_order);

  int isrz = CoordSys::IsRZ();

  int nghost = 1; // just like bill

  MultiFab alpha;  // alpha should be the same size as volume
  alpha.define(grids,BL_SPACEDIM,nghost,Fab_allocate);
  alpha.setVal(0.0);
  int k;
  if (a!=0.0) {
    for (k = 0; k < grids.length(); k++) {
      const BOX &bx = alpha.box(k);
      int rlen      = bx.length(0);
      Array<REAL> rcen;
      rcen.resize(rlen);
      parent->Geom(level).GetCellLoc(rcen, bx, 0);
      const int *lo       = bx.loVect();
      const int *hi       = bx.hiVect();
      REAL *alpha_dat     = alpha[k].dataPtr();
      BOX abx             = grow(bx,alpha.nGrow());
      const int *alo      = abx.loVect();
      const int *ahi      = abx.hiVect();
      const REAL *rcendat = rcen.dataPtr();
      const REAL *voli    = volume[k].dataPtr();
      BOX vbox            = grow(volume.box(k),volume.nGrow());
      const int *vlo      = vbox.loVect();
      const int *vhi      = vbox.hiVect();

      FARRAYBOX& Rh = (*rho_half)[k];
      DEF_LIMITS(Rh,rho_dat,rlo,rhi);

      FARRAYBOX& betax = (*beta[0])[k];
      DEF_CLIMITS(betax,betax_dat,betax_lo,betax_hi);

      FARRAYBOX& betay = (*beta[1])[k];
      DEF_CLIMITS(betay,betay_dat,betay_lo,betay_hi);

      FORT_SET_TENSOR_ALPHA(alpha_dat, ARLIM(alo), ARLIM(ahi),
                  lo, hi, rcendat, ARLIM(lo), ARLIM(hi), &b,
                  voli, ARLIM(vlo), ARLIM(vhi),
                  rho_dat,ARLIM(rlo),ARLIM(rhi),
                  betax_dat,ARLIM(betax_lo),ARLIM(betax_hi),
                  betay_dat,ARLIM(betay_lo),ARLIM(betay_hi),&isrz);
    }
  }

  tensor_op->setScalars(a,b);

  tensor_op->aCoefficients(alpha);
  for (int n = 0; n < BL_SPACEDIM; n++) {
    MultiFab bcoeffs(area[n].boxArray(),1,nghost,Fab_allocate);
    bcoeffs.setVal(0.0);
    bcoeffs.copy(area[n]);
    for (k = 0; k < grids.length(); k++) {
      bcoeffs[k].mult(dx[n]);
      bcoeffs[k].mult((*beta[n])[k]);
    }
    tensor_op->bCoefficients(bcoeffs,n);
  }

  return tensor_op;
#endif
}
#endif

ABecLaplacian* Diffusion::getViscOp(int comp, REAL a, REAL b,
				    REAL time, ViscBndry & visc_bndry,
				    MultiFab* rho_half, int rho_flag, REAL* rhsscale,
                                    MultiFab** beta, MultiFab* alpha_in)
{

  int allnull, allthere;
  checkBeta(beta, allthere, allnull);

  const REAL* dx = caller->Geom().CellSize();

  getBndryData(visc_bndry,comp,1,time,rho_flag);

  ABecLaplacian *visc_op = new ABecLaplacian(grids,visc_bndry,dx);
  visc_op->maxOrder(max_order);

  int usehoop=( (comp==Xvel) && (CoordSys::IsRZ())   );

  int useden = (rho_flag == 1);

  MultiFab alpha;  // alpha should be the same size as volume
  alpha.define(grids,1,GEOM_GROW,Fab_allocate);
  int k;
  for (k = 0; k < grids.length(); k++) {
    const BOX &bx = alpha.box(k);
    int rlen      = bx.length(0);
    Array<REAL> rcen;
    rcen.resize(rlen);
    parent->Geom(level).GetCellLoc(rcen, bx, 0);
    const int *lo       = bx.loVect();
    const int *hi       = bx.hiVect();
    REAL *dat           = alpha[k].dataPtr();
    BOX abx             = grow(bx,alpha.nGrow());
    const int *alo      = abx.loVect();
    const int *ahi      = abx.hiVect();
    const REAL *rcendat = rcen.dataPtr();
    const REAL *voli    = volume[k].dataPtr();
    BOX vbox            = grow(volume.box(k),volume.nGrow());
    const int *vlo      = vbox.loVect();
    const int *vhi      = vbox.hiVect();

    FARRAYBOX& Rh = (*rho_half)[k];
    DEF_LIMITS(Rh,rho_dat,rlo,rhi);

    FORT_SETALPHA(dat, ARLIM(alo), ARLIM(ahi),
                  lo, hi, rcendat, ARLIM(lo), ARLIM(hi), &b,
                  voli, ARLIM(vlo), ARLIM(vhi),
                  rho_dat,ARLIM(rlo),ARLIM(rhi),&usehoop,&useden);


  }

//  visc_op->setScalars(a,dx[0]*b);
  if (rho_flag == 2) {
    // note: using conservative diffing for rho*T
    MultiFab &S_new = caller->get_new_data(State_Type);
    for (k = 0; k < grids.length(); k++) {
      const BOX& grd = grids[k];
      alpha[k].mult(S_new[k],grd,Density,0,1);
    }
  }

 if(alpha_in!=NULL) {
    for (k = 0; k < grids.length(); k++) {
      const BOX& grd = grids[k];
      alpha[k].mult((*alpha_in)[k],grd,0,0,1);
    }
  }

  if(rhsscale!=NULL) {
    if(scale_abec) {
      *rhsscale = 1.0/alpha.max(0);
    } else {
      *rhsscale = 1.0;
    }
    visc_op->setScalars(a*(*rhsscale),b*(*rhsscale));
  } else {
    visc_op->setScalars(a,b);
  }

  if(allnull) {
    visc_op->aCoefficients(alpha);
    for (int n = 0; n < BL_SPACEDIM; n++) {
      MultiFab bcoeffs(area[n].boxArray(),1,0);
      bcoeffs.copy(area[n]);
      for (k = 0; k < grids.length(); k++) {
        bcoeffs[k].mult(dx[n]);
      }
      visc_op->bCoefficients(bcoeffs,n);
    }
  } else {
    visc_op->aCoefficients(alpha);
    for (int n = 0; n < BL_SPACEDIM; n++) {
      MultiFab bcoeffs(area[n].boxArray(),1,0);
      bcoeffs.copy(area[n]);
      for (k = 0; k < grids.length(); k++) {
        bcoeffs[k].mult((*beta[n])[k]);
        bcoeffs[k].mult(dx[n]);
      }
      visc_op->bCoefficients(bcoeffs,n);
    }
  }

  return visc_op;
}

ABecLaplacian* Diffusion::getViscOp(int comp, REAL a, REAL b,
				    MultiFab* rho_half, int rho_flag, REAL* rhsscale,
                                    MultiFab** beta, MultiFab* alpha_in)
{
  int allnull, allthere;
  checkBeta(beta, allthere, allnull);

  const REAL* dx = caller->Geom().CellSize();
  const BOX& domain = caller->Geom().Domain();
  const BCRec& bc = caller->get_desc_lst()[State_Type].getBC(comp);

  IntVect ref_ratio;
  if (level > 0) {
    ref_ratio = parent->refRatio(level-1);
  } else {
    ref_ratio = IntVect::TheUnitVector();
  }

  ViscBndry bndry(grids,1,domain);
  bndry.setHomogValues(bc, ref_ratio);

  ABecLaplacian *visc_op = new ABecLaplacian(grids,bndry,dx);
  visc_op->maxOrder(max_order);

  int usehoop=( (comp==Xvel) && (CoordSys::IsRZ())   );

  int useden = (rho_flag == 1);

  MultiFab alpha;  // alpha should be the same size as volume
  alpha.define(grids,1,GEOM_GROW,Fab_allocate);
  int k;
  for (k = 0; k < grids.length(); k++) {
    const BOX &bx = alpha.box(k);
    int rlen      = bx.length(0);
    Array<REAL> rcen;
    rcen.resize(rlen);
    parent->Geom(level).GetCellLoc(rcen, bx, 0);
    const int *lo       = bx.loVect();
    const int *hi       = bx.hiVect();
    REAL *dat           = alpha[k].dataPtr();
    BOX abx             = grow(bx,alpha.nGrow());
    const int *alo      = abx.loVect();
    const int *ahi      = abx.hiVect();
    const REAL *rcendat = rcen.dataPtr();
    const REAL *voli    = volume[k].dataPtr();
    BOX vbox            = grow(volume.box(k),volume.nGrow());
    const int *vlo      = vbox.loVect();
    const int *vhi      = vbox.hiVect();

    FARRAYBOX& Rh = (*rho_half)[k];
    DEF_LIMITS(Rh,rho_dat,rlo,rhi);

    FORT_SETALPHA(dat, ARLIM(alo), ARLIM(ahi),
                  lo, hi, rcendat, ARLIM(lo), ARLIM(hi), &b,
                  voli, ARLIM(vlo), ARLIM(vhi),
                  rho_dat,ARLIM(rlo),ARLIM(rhi),&usehoop,&useden);


  }

  if (rho_flag == 2) {
    // note: using conservative diffing for rho*S
    MultiFab &S_new = caller->get_new_data(State_Type);
    for (k = 0; k < grids.length(); k++) {
      const BOX& grd = grids[k];
      alpha[k].mult(S_new[k],grd,Density,0,1);
    }
  }

  if (alpha_in!=NULL) {
    for (k = 0; k < grids.length(); k++) {
      const BOX& grd = grids[k];
      alpha[k].mult((*alpha_in)[k],grd,0,0,1);
    }
  }

  if(rhsscale!=NULL) {
    if(scale_abec) {
      *rhsscale = 1.0/alpha.max(0);
    } else {
      *rhsscale = 1.0;
    }
    visc_op->setScalars(a*(*rhsscale),b*(*rhsscale));
  } else {
    visc_op->setScalars(a,b);
  }

  visc_op->aCoefficients(alpha);
  if(allnull) {
    for (int n = 0; n < BL_SPACEDIM; n++) {
      MultiFab bcoeffs(area[n].boxArray(),1,0);
      bcoeffs.copy(area[n]);
      bcoeffs.mult(dx[n],0,1,0);
      visc_op->bCoefficients(bcoeffs,n);
    }
  } else {
    for (int n = 0; n < BL_SPACEDIM; n++) {
      MultiFab bcoeffs(area[n].boxArray(),1,0);
      bcoeffs.copy(area[n]);
      for (k = 0; k < grids.length(); k++) { 
        bcoeffs[k].mult((*beta[n])[k]);
      }
      bcoeffs.mult(dx[n],0,1,0);
      visc_op->bCoefficients(bcoeffs,n);
    }
  }

  return visc_op;
}

void Diffusion::getViscTerms(MultiFab& visc_terms, int src_comp, int comp, 
			     REAL time, int rho_flag, MultiFab** beta)
{
  int allnull, allthere;
  checkBeta(beta, allthere, allnull);

  // Before computing the godunov predicitors we may have to
  // precompute the viscous source terms.  To do this we must
  // construct a Laplacian operator, set the coeficients and apply
  // it to the time N data.  First, however, we must precompute the
  // tine N bndry values.  We will do this for each scalar that diffuses.
  const REAL* dx = caller->Geom().CellSize();
  MultiFab& S = caller->get_data(State_Type,time);
  visc_terms.setVal(0.0,comp-src_comp,1,1);
  int ngrd = grids.length();

  // FIXME
  // LinOp classes cannot handle multcomponent MultiFabs yet,
  // construct the components one at a time and copy to visc_terms.
  MultiFab visc_tmp(grids,1,1,Fab_allocate);
  MultiFab s_tmp(grids,1,1,Fab_allocate);

  int k,n;

  if (is_diffusive[comp]) {
    ViscBndry visc_bndry;
    getBndryData(visc_bndry,comp,1,time,rho_flag);

    // set up operator and apply to compute viscous terms
    REAL a = 0.0;
    REAL b;
    if (allnull) {
      b = -visc_coef[comp];
    } else {
      b = -1.0;
    }
    ABecLaplacian visc_op(grids,visc_bndry,dx);
//    visc_op.setScalars(a,b*dx[0]);
    visc_op.setScalars(a,b);
    visc_op.maxOrder(max_order);

    if(allnull) {
      for (n = 0; n < BL_SPACEDIM; n++) {
        MultiFab bcoeffs(area[n].boxArray(),1,0);
        bcoeffs.copy(area[n]);
        for (k = 0; k < ngrd; k++) { 
          bcoeffs[k].mult(dx[n]);
        }
        visc_op.bCoefficients(bcoeffs,n);
      }
    } else {
      for (n = 0; n < BL_SPACEDIM; n++) {
        MultiFab bcoeffs(area[n].boxArray(),1,0);
        bcoeffs.copy(area[n]);
        for (k = 0; k < ngrd; k++) { 
          bcoeffs[k].mult((*beta[n])[k]);
          bcoeffs[k].mult(dx[n]);
        }
        visc_op.bCoefficients(bcoeffs,n);
      }
    }

    // FIXME
    // copy to single component multifab for operator classes
    s_tmp.copy(S,comp,0,1,1);

    if(rho_flag==2) {
// we want to evaluate (div beta grad) S, not rho*S
      for (k=0;k<ngrd;k++) { 
        assert (S[k].min(Density)>0.0);
        s_tmp[k].divide(S[k],Density,0,1);
      }
    }

    visc_op.apply(visc_tmp,s_tmp);

    // must divide by volume
    for (k = 0; k < ngrd; k++) {
      visc_tmp[k].divide(volume[k],grids[k],0,0,1);
    }

#if (BL_SPACEDIM == 2)
    if(comp==Xvel && CoordSys::IsRZ()) {
      for (k = 0; k < ngrd; k++) {
	//visc_tmp[k] += -mu * u / r^2
	const BOX& bx  = visc_tmp.box(k);
	BOX  vbx       = grow(bx,visc_tmp.nGrow());
	BOX  sbx       = grow(s_tmp.box(k),s_tmp.nGrow());
	int rlen       = bx.length(0);
	Array<REAL> rcen;
	rcen.resize(rlen);
	parent->Geom(level).GetCellLoc(rcen, bx, 0);
	const int *lo       = bx.loVect();
	const int *hi       = bx.hiVect();
	const int *vlo      = vbx.loVect();
	const int *vhi      = vbx.hiVect();
	const int *slo      = sbx.loVect();
	const int *shi      = sbx.hiVect();
	REAL *vdat          = visc_tmp[k].dataPtr();
	REAL *sdat          = s_tmp[k].dataPtr();
	const REAL *rcendat = rcen.dataPtr();
	const REAL mu       = visc_coef[comp];
        FORT_HOOPSRC(ARLIM(lo), ARLIM(hi),
                     vdat, ARLIM(vlo), ARLIM(vhi),
                     sdat, ARLIM(slo), ARLIM(shi),
                     rcendat, &mu);


      }
    }
#endif

    assert(visc_tmp.nGrow() > 0);
    const BoxArray& ba = visc_tmp.boxArray();
    int ngrd = ba.length();
//    NOTE: THIS IS JUST A HACK - DONT KNOW HOW TO FILL VISC TERMS
//          IN GHOST CELLS OUTSIDE FINE GRIDS
    for (k = 0; k < ngrd; k++) {
      const BOX& grd(ba[k]);
      const int* lo = grd.loVect();
      const int* hi = grd.hiVect();
      FArrayBox& visc = visc_tmp[k];
      int ncomp = visc.nComp();
      DEF_LIMITS(visc,vdat,vlo,vhi);
      FORT_VISCEXTRAP(vdat,ARLIM(vlo),ARLIM(vhi),lo,hi,&ncomp);
    }

      // if periodic, copy into periodic translates of visc_tmp
    if ( caller->Geom().isAnyPeriodic() ) {

      const BOX& domain = caller->Geom().Domain();
      Array<IntVect> pshifts(27);
 
      FARRAYBOX dest;
      for (k = 0; k < ngrd; k++) {
  
        const BOX& dbox(visc_tmp[k].box());
        dest.resize(dbox,1);
        dest.copy(visc_tmp[k],dbox,0,dbox,0,1);
        caller->Geom().periodicShift( domain, dest.box(), pshifts);
        for( int iiv=0; iiv<pshifts.length(); iiv++) {
          IntVect iv=pshifts[iiv];
          dest.shift(iv);
          if( dbox.intersects(domain) ) visc_tmp.copy(dest);
          dest.shift(-iv);
          visc_tmp[k].copy(dest,dbox,0,dbox,0,1);
        }
      }
    }

    // copy from valid regions of overlapping grids
    visc_tmp.FillBoundary();

    // FIXME
    visc_terms.copy(visc_tmp,0,comp-src_comp,1,1);
  }
}

void Diffusion::getTensorViscTerms(MultiFab& visc_terms, 
			     REAL time, MultiFab** beta)
{
#if (BL_SPACEDIM==3)
    cout << "Diffusion::getTensorViscTerms :  " <<
            "not yet implemented for 3-D" << endl;
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0" <<
            " and rerun" << endl;
    exit(0);
#elif  !defined (USE_TENSOR)
    cout << "Diffusion::getTensorViscTerms :  " <<
            "USE_TENSOR must be defined at compile time" << endl;
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0" <<
            " and rerun" << endl;
    exit(0);
#else

  int allthere;
  checkBeta(beta, allthere);

  int src_comp = Xvel;
  int ncomp = visc_terms.nComp();
  if(ncomp < BL_SPACEDIM) {
    cout << "Diffusion::getTensorViscTerms : visc_terms must have" << endl;
    cout << "  at least BL_SPACEDIM components" << endl;
    exit(0);
  }
  int vel_ncomp = BL_SPACEDIM;


  // Before computing the godunov predicitors we may have to
  // precompute the viscous source terms.  To do this we must
  // construct a Laplacian operator, set the coeficients and apply
  // it to the time N data.  First, however, we must precompute the
  // tine N bndry values.  We will do this for each scalar that diffuses.
  const REAL* dx = caller->Geom().CellSize();
  MultiFab& S = caller->get_data(State_Type,time);
  visc_terms.setVal(0.0,src_comp,BL_SPACEDIM,1);
  int ngrd = grids.length();

  // FIXME
  // LinOp classes cannot handle multcomponent MultiFabs yet,
  // construct the components one at a time and copy to visc_terms.
  MultiFab visc_tmp(grids,BL_SPACEDIM,1,Fab_allocate);
  MultiFab s_tmp(grids,BL_SPACEDIM,1,Fab_allocate);

  int k,n;

  if (is_diffusive[src_comp]) {
    ViscBndry2D visc_bndry;
    getTensorBndryData(visc_bndry,time);

    // set up operator and apply to compute viscous terms
    REAL a = 0.0;
    REAL b = -1.0;

    DivVis tensor_op(grids,visc_bndry,dx);
    tensor_op.maxOrder(tensor_max_order);
    tensor_op.setScalars(a,b);

    int nghost = 1; // like bill
    MultiFab alpha;  // alpha should be the same size as volume
    alpha.define(grids,BL_SPACEDIM,nghost,Fab_allocate);
    alpha.setVal(0.0);
    tensor_op.aCoefficients(alpha);
    for (n = 0; n < BL_SPACEDIM; n++) {
      MultiFab bcoeffs(area[n].boxArray(),1,nghost,Fab_allocate);
      bcoeffs.setVal(0.0);
      bcoeffs.copy(area[n]);
      for (k = 0; k < ngrd; k++) { 
        bcoeffs[k].mult(dx[n]);
        bcoeffs[k].mult((*beta[n])[k]);
      }
      tensor_op.bCoefficients(bcoeffs,n);
    }

    s_tmp.copy(S,Xvel,0,BL_SPACEDIM,1);

    tensor_op.apply(visc_tmp,s_tmp);

    // must divide by volume
    for (k = 0; k < ngrd; k++) {
      visc_tmp[k].divide(volume[k],grids[k],0,0,1);
    }

#if (BL_SPACEDIM == 2)
    if(CoordSys::IsRZ()) {
      int fort_xvel_comp = Xvel+1;
      for (k = 0; k < ngrd; k++) {
	//visc_tmp[k] += -mu * u / r^2
	const BOX& bx  = visc_tmp.box(k);
	BOX  vbx       = grow(bx,visc_tmp.nGrow());
	BOX  sbx       = grow(s_tmp.box(k),s_tmp.nGrow());
	int rlen       = bx.length(0);
	Array<REAL> rcen;
	rcen.resize(rlen);
	parent->Geom(level).GetCellLoc(rcen, bx, 0);
	const int *lo       = bx.loVect();
	const int *hi       = bx.hiVect();
	const int *vlo      = vbx.loVect();
	const int *vhi      = vbx.hiVect();
	const int *slo      = sbx.loVect();
	const int *shi      = sbx.hiVect();
	REAL *vdat          = visc_tmp[k].dataPtr();
	REAL *sdat          = s_tmp[k].dataPtr();
	const REAL *rcendat = rcen.dataPtr();

        FARRAYBOX& betax = (*beta[0])[k];
        DEF_CLIMITS(betax,betax_dat,betax_lo,betax_hi);

        FARRAYBOX& betay = (*beta[1])[k];
        DEF_CLIMITS(betay,betay_dat,betay_lo,betay_hi);

        FORT_TENSOR_HOOPSRC(&fort_xvel_comp,ARLIM(lo), ARLIM(hi),
                     vdat, ARLIM(vlo), ARLIM(vhi),
                     sdat, ARLIM(slo), ARLIM(shi),
                     rcendat, 
                     betax_dat,ARLIM(betax_lo),ARLIM(betax_hi),
                     betay_dat,ARLIM(betay_lo),ARLIM(betay_hi));


      }
    }
#endif

    assert(visc_tmp.nGrow() > 0);
    const BoxArray& ba = visc_tmp.boxArray();
    int ngrd = ba.length();
//    NOTE: THIS IS JUST A HACK - DONT KNOW HOW TO FILL VISC TERMS
//          IN GHOST CELLS OUTSIDE FINE GRIDS
    for (k = 0; k < ngrd; k++) {
      const BOX& grd(ba[k]);
      const int* lo = grd.loVect();
      const int* hi = grd.hiVect();
      FArrayBox& visc = visc_tmp[k];
      DEF_LIMITS(visc,vdat,vlo,vhi);
      FORT_VISCEXTRAP(vdat,ARLIM(vlo),ARLIM(vhi),lo,hi,&vel_ncomp);
    }

      // if periodic, copy into periodic translates of visc_tmp
    if ( caller->Geom().isAnyPeriodic() ) {

      const BOX& domain = caller->Geom().Domain();
      Array<IntVect> pshifts(27);
 
      FARRAYBOX dest;
      for (k = 0; k < ngrd; k++) {
  
        const BOX& dbox(visc_tmp[k].box());
        dest.resize(dbox,1);
        dest.copy(visc_tmp[k],dbox,0,dbox,0,1);
        caller->Geom().periodicShift( domain, dest.box(), pshifts);
        for( int iiv=0; iiv<pshifts.length(); iiv++) {
          IntVect iv=pshifts[iiv];
          dest.shift(iv);
          if( dbox.intersects(domain) ) visc_tmp.copy(dest);
          dest.shift(-iv);
          visc_tmp[k].copy(dest,dbox,0,dbox,0,vel_ncomp);
        }
      }
    }

    // copy from valid regions of overlapping grids
    visc_tmp.FillBoundary();

    // FIXME
    visc_terms.copy(visc_tmp,0,0,BL_SPACEDIM,1);

  }
#endif
}

void Diffusion::getBndryData(ViscBndry& bndry, int src_comp,
			     int num_comp, REAL time, int rho_flag)
{
    if (num_comp != 1) {
       cout << "NEED NUM_COMP = 1" << endl;
       abort();
    }
    const BCRec& bc = caller->get_desc_lst()[State_Type].getBC(src_comp);
    bndry.define(grids,num_comp,caller->Geom());

    MultiFab& S = caller->get_data(State_Type,time);
    S.FillBoundary(src_comp,num_comp);
    caller->setPhysBoundaryValues(State_Type,src_comp,num_comp,time);
    if (rho_flag == 2) { 
      S.FillBoundary(Density,1);
      caller->setPhysBoundaryValues(State_Type,Density,1,time);
    }

    int use_old_code = 0; // FOR DEBUGGING ONLY
    if(rho_flag!=2 || use_old_code) {
      if (level == 0) {
	bndry.setBndryValues(S,src_comp,0,num_comp,bc);
      } else {
	BoxArray cgrids(grids);
	cgrids.coarsen(crse_ratio);
	BndryRegister crse_br(cgrids,0,1,1,num_comp);
        crse_br.setVal(1.0e30);
	coarser->FillBoundary(crse_br,src_comp,0,num_comp,time);
	bndry.setBndryValues(crse_br,0,S,src_comp,0,num_comp,
			     crse_ratio,bc);
      }
    } else {
// we want to fill the bndry with S, not rho*S
      MultiFab Stmp(S.boxArray(),1,S.nGrow());
      Stmp.copy(S,src_comp,0,1,S.nGrow());
      int ngrd = S.boxArray().length();
      for (int i=0;i<ngrd;i++) { 
        assert (S[i].min(Density)>0.0);
        Stmp[i].divide(S[i],Density,0,1);
      }
      if (level == 0) {
	bndry.setBndryValues(Stmp,0,0,num_comp,bc);
      } else {
	BoxArray cgrids(grids);
	cgrids.coarsen(crse_ratio);
	BndryRegister crse_br(cgrids,0,1,1,num_comp);
        crse_br.setVal(1.0e30);
	coarser->FillBoundary(crse_br,src_comp,0,num_comp,time,rho_flag);
	bndry.setBndryValues(crse_br,0,Stmp,0,0,num_comp,
			     crse_ratio,bc);
      }      
      
    }

}

void Diffusion::FillBoundary(BndryRegister& bdry, int src_comp, int dest_comp,
			 int num_comp, REAL time, int rho_flag)
{
    REAL old_time = caller->get_state_data(State_Type).prevTime();
    REAL new_time = caller->get_state_data(State_Type).curTime();
    REAL eps = 0.001*(new_time - old_time);
    assert( (time > old_time-eps) && (time < new_time + eps));

    // FillBoundary can be called before old data is defined, so
    // "need_old_data" really should be "have_old_data".
    // "need_old_data" implies "have_old_data", however, so
    // "need_old_data" is suffcient for checking for the existence
    // of old data--rbp
    int need_old_data = !(time > new_time - eps);

    MultiFab& S_new = caller->get_new_data(State_Type);
    // the next line is OK even if S_old is not defined yet
    MultiFab& S_old = caller->get_old_data(State_Type);

    if (need_old_data) {
      assert (S_old.nGrow() == S_new.nGrow());
    }
    int n_grow = S_new.nGrow();

    MultiFab sold_tmp;
    MultiFab snew_tmp;

    if(need_old_data) {
      sold_tmp.define(S_old.boxArray(),num_comp,n_grow,Fab_allocate);
      sold_tmp.copy(S_old,src_comp,0,num_comp,n_grow);
      sold_tmp.FillBoundary();
    }

    snew_tmp.define(S_new.boxArray(),num_comp,n_grow,Fab_allocate);
    snew_tmp.copy(S_new,src_comp,0,num_comp,n_grow);
    snew_tmp.FillBoundary();

    MultiFab rho_old;
    MultiFab rho_new;

    if (rho_flag == 2) {
      if(need_old_data) {
        rho_old.define(S_new.boxArray(),1,n_grow,Fab_allocate);
        rho_old.copy(S_old,Density,0,1,n_grow);
        rho_old.FillBoundary();
      }
      rho_new.define(S_new.boxArray(),1,n_grow,Fab_allocate);
      rho_new.copy(S_new,Density,0,1,n_grow);
      rho_new.FillBoundary();
    }

    int ngrd = S_new.boxArray().length();

//  Impose periodic boundary conditions if necessary on the state
//    data before it is copied into the crse bndry registers
    const Geometry& crse_geom  = parent->Geom(level);
    if ( crse_geom.isAnyPeriodic() ) {

      const BOX& crse_domain = crse_geom.Domain();
      Array<IntVect> pshifts(27);
      for (int grd = 0; grd < ngrd; grd++) {
        
        BOX cgrd(snew_tmp[grd].box());
        crse_geom.periodicShift(crse_domain, cgrd, pshifts);

        for( int iiv=0; iiv<pshifts.length(); iiv++) {
            IntVect iv=pshifts[iiv];

            if (need_old_data) {
              sold_tmp[grd].shift(iv);
              S_old.copy(sold_tmp[grd],src_comp,0,num_comp);
              sold_tmp[grd].shift(-iv);
            }

            snew_tmp[grd].shift(iv);
            S_new.copy(snew_tmp[grd],src_comp,0,num_comp);
            snew_tmp[grd].shift(-iv);
        }

        if (rho_flag == 2) {

         for( int iiv=0; iiv<pshifts.length(); iiv++) {
            IntVect iv=pshifts[iiv];

            if (need_old_data) {
              rho_old[grd].shift(iv);
              S_old.copy(rho_old[grd],Density,0,1);
              rho_old[grd].shift(-iv);
            }

            rho_new[grd].shift(iv);
            S_new.copy(rho_new[grd],Density,0,1);
            rho_new[grd].shift(-iv);
         }
        }
      }
    }

    int n_ghost = 1;

    if (rho_flag==2) {
//    we want to fill the bndry with S, not rho*S
      for (int i=0;i<ngrd;i++) { 
        if (need_old_data) sold_tmp[i].divide(rho_old[i],0,0,num_comp);
        assert (rho_new[i].min()>0.0);
        snew_tmp[i].divide(rho_new[i],0,0,num_comp);
      }
    }

    if (time < old_time - eps || time > new_time + eps) {
      BoxLib::Error("FillBoundary: bad interp time");
    } else if (time < old_time + eps) {
      bdry.copyFrom(sold_tmp,n_ghost,0,dest_comp,num_comp);
    } else if (time > new_time - eps) {
      bdry.copyFrom(snew_tmp,n_ghost,0,dest_comp,num_comp);
    } else {
      REAL a = (new_time - time)/(new_time - old_time);
      REAL b = (time - old_time)/(new_time - old_time);
      bdry.linComb(a,sold_tmp,0,b,snew_tmp,0,
                   dest_comp,num_comp,n_ghost);
    }
}

#if defined (USE_TENSOR)
void Diffusion::getTensorBndryData(
#if (BL_SPACEDIM==2) 
                            ViscBndry2D& bndry, 
#else 
                            ViscBndry3D& bndry, 
#endif
			    REAL time)
{
#if (BL_SPACEDIM==3)
    cout << "Diffusion::getTensorBndryData :  " <<
            "not yet implemented for 3-D" << endl;
    cout << "set use_dv_constant_mu = 1 with velocity visc_coef >=0.0" <<
            " and rerun" << endl;
    exit(0);
#else
    int num_comp = BL_SPACEDIM;
    int src_comp = Xvel;

    // Create the BCRec's interpreted by ViscBndry objects

    Array<BCRec> bcarray(2*BL_SPACEDIM);
    for (int idim=0; idim < BL_SPACEDIM; idim++) {
      bcarray[idim] = 
        caller->get_desc_lst()[State_Type].getBC(src_comp+idim);
      bcarray[idim+BL_SPACEDIM] = 
        BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
		     D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    }

    bndry.define(grids,2*num_comp,caller->Geom());

    MultiFab& S = caller->get_data(State_Type,time);
    S.FillBoundary(src_comp,num_comp);
    caller->setPhysBoundaryValues(State_Type,src_comp,num_comp,time);

    if (level == 0) {
      bndry.setBndryValues(S,src_comp,0,num_comp,bcarray);
    } else {
      BoxArray cgrids(grids);
      cgrids.coarsen(crse_ratio);
      BndryRegister crse_br(cgrids,0,1,1,num_comp);
      crse_br.setVal(1.0e30);
      coarser->FillBoundary(crse_br,src_comp,0,num_comp,time);
      bndry.setBndryValues(crse_br,0,S,src_comp,0,num_comp,
			     crse_ratio[0],bcarray);
    }

}
#endif
#endif

//-----------------------------------------------------------------

void Diffusion::checkBetas(MultiFab** beta1, MultiFab** beta2,
                           int& allthere, int& allnull)
{
  int allnull1, allnull2, allthere1, allthere2;
  checkBeta(beta1,allthere1,allnull1);
  checkBeta(beta2,allthere2,allnull2);
  allnull  = allnull1 && allnull2;
  allthere = allthere1 && allthere2;
  if(!(allthere || allnull)) {
    cout << "Diffusion::checkBetas : all betas must either be" << endl;
    cout << "  all NULL or all non-NULL" << endl;
    exit(0);
  }
}

//-----------------------------------------------------------------

void Diffusion::checkBeta(MultiFab** beta,
                           int& allthere, int& allnull)
{
  allnull = 1;
  allthere = beta != NULL;
  if(allthere) {
    for (int dim=0; dim<BL_SPACEDIM; dim++) {
      allnull = allnull && beta[dim] == NULL;
      allthere = allthere && beta[dim] != NULL;
    }
  }
  if(!(allthere || allnull)) {
    cout << "Diffusion::checkBeta : all betas must either be" << endl;
    cout << "  all NULL or all non-NULL" << endl;
    exit(0);
  }
}

//-----------------------------------------------------------------

void Diffusion::checkBeta(MultiFab** beta,
                           int& allthere)
{
  allthere = beta != NULL;
  if(allthere) {
    for (int dim=0; dim<BL_SPACEDIM; dim++) {
      allthere = allthere && beta[dim] != NULL;
    }
  }
  if(!allthere) {
    cout << "Diffusion::checkBeta : all betas must be" << endl;
    cout << "  all non-NULL" << endl;
    exit(0);
  }
}

// -------------------------------------------------------------
void Diffusion::allocFluxBoxesLevel(MultiFab**& fluxbox, 
                           int nghost, int nvar)
{
    fluxbox = new MultiFab*[BL_SPACEDIM];
    for (int dir = 0; dir < BL_SPACEDIM; dir++) {
      BoxArray edge_boxes(grids);
      edge_boxes.surroundingNodes(dir);
      fluxbox[dir] = new MultiFab(edge_boxes,nvar,nghost,Fab_allocate);
    }
}

// -------------------------------------------------------------
void Diffusion::removeFluxBoxesLevel(MultiFab**& fluxbox) 
{
  if(fluxbox!=NULL) {
    for (int i = 0; i<BL_SPACEDIM; i++) 
      if (fluxbox[i] != NULL) delete fluxbox[i];
    delete [] fluxbox;
  }
}

#ifdef USE_NAVIERSTOKES
// -------------------------------------------------------------
// this routine computes the vector div mu SI, where I is the identity 
// tensor, S = div U, and mu is constant

void
Diffusion::compute_divmusi(REAL time, REAL mu,
              MultiFab& divmusi)
{

  FARRAYBOX divu;
  const REAL* dx = caller->Geom().CellSize();
  NavierStokes& ns_level =  *(NavierStokes*) &(parent->getLevel(level));
  int nghost = divmusi.nGrow();
  for (int i = 0; i < grids.length(); i++) {
    const BOX& bx = grids[i];
    const int *lo       = bx.loVect();
    const int *hi       = bx.hiVect();
    DEF_LIMITS(divmusi[i],divmusi_dat,divmusilo,divmusihi);

    ns_level.getDivCond(divu,i,1,time);
    DEF_CLIMITS(divu,divu_dat,divulo,divuhi);

    FORT_DIV_MU_SI(lo,hi,dx,&mu,
        ARLIM(divulo), ARLIM(divuhi), divu_dat,
        ARLIM(divmusilo), ARLIM(divmusihi), divmusi_dat);

    if(nghost>0) {
      int ncomp = BL_SPACEDIM;
      FORT_VISCEXTRAP(divmusi_dat,ARLIM(divmusilo),ARLIM(divmusihi),lo,hi,
                      &ncomp);
    }
  }

}

// -------------------------------------------------------------
// this routine computes the vector div beta SI, where I is the identity 
// tensor, S = div U, and beta is non-constant

void
Diffusion::compute_divmusi(REAL time,
              MultiFab** beta, MultiFab& divmusi)
{

  FARRAYBOX divu;
  const REAL* dx = caller->Geom().CellSize();
  NavierStokes& ns_level =  *(NavierStokes*) &(parent->getLevel(level));
  int nghost = divmusi.nGrow();
  for (int i = 0; i < grids.length(); i++) {
    const BOX& bx = grids[i];
    const int *lo       = bx.loVect();
    const int *hi       = bx.hiVect();
    DEF_LIMITS(divmusi[i],divmusi_dat,divmusilo,divmusihi);

    ns_level.getDivCond(divu,i,1,time);
    DEF_CLIMITS((*beta[0])[i],betax,betaxlo,betaxhi);
    DEF_CLIMITS((*beta[1])[i],betay,betaylo,betayhi);
#if (BL_SPACEDIM==3)
    DEF_CLIMITS((*beta[2])[i],betaz,betazlo,betazhi);
#endif
    DEF_CLIMITS(divu,divu_dat,divulo,divuhi);

    FORT_DIV_VARMU_SI(lo,hi,dx,
      ARLIM(divulo), ARLIM(divuhi), divu_dat,
      ARLIM(betaxlo), ARLIM(betaxhi), betax,
      ARLIM(betaylo), ARLIM(betayhi), betay,
#if (BL_SPACEDIM==3)
      ARLIM(betazlo), ARLIM(betazhi), betaz,
#endif
      ARLIM(divmusilo), ARLIM(divmusihi), divmusi_dat);

    if(nghost>0) {
      int ncomp = BL_SPACEDIM;
      FORT_VISCEXTRAP(divmusi_dat,ARLIM(divmusilo),ARLIM(divmusihi),lo,hi,
                      &ncomp);
    }
  }

}
#endif

