//
// $Id: Diffusion.H,v 1.2 1997-09-23 19:25:40 lijewski Exp $
//

#ifndef _Diffusion_H_
#define _Diffusion_H_

// comment out this line to use diffusion class outside
// the context of NavierStokes and classes derived from it
#define USE_NAVIERSTOKES 1

#include <Box.H>
#include <BoxArray.H>
#include <Geometry.H>
#include <LevelBld.H>
#include <BC_TYPES.H>
#include <AmrLevel.H>
#include <ErrorList.H>
#include <FluxRegister.H>
#include <ABecLaplacian.H>
#include <ViscBndry.H>

// include files for tensor solve
#if(BL_SPACEDIM==2) && defined (USE_TENSOR)
#include <ViscBndry2D.H>
#include <DivVis.H>
#include <LO_BCTYPES.H>
#include <MCMultiGrid.H>
#include <MCCGSolver.H>
#endif

#ifndef _NavierStokes_H_
enum StateType {State_Type=0, Press_Type}; 
#  if (BL_SPACEDIM == 2)
enum StateNames  { Xvel=0, Yvel, Density};
#  else
enum StateNames  { Xvel=0, Yvel, Zvel, Density};
#  endif
#endif

class Diffusion {

  Amr* parent;
  AmrLevel* caller;
  Diffusion *coarser, *finer;
  const BoxArray& grids;
  int level, NUM_STATE;
  IntVect crse_ratio;
  FluxRegister *viscflux_reg;

  // volume and area fractions
  MultiFab& volume;
  MultiFab* area;

  static int first;
  static int do_reflux;
  static int use_cg_solve;
  static int use_tensor_cg_solve;
  static int max_order;
  static int tensor_max_order;
  static bool use_mg_precond_flag;
  static int use_dv_constant_mu_def;
  static int scale_abec;
  static Array<int> is_diffusive;     // does variable diffuse?
  static Array<REAL> visc_coef;       // const coef viscosity terms
  static REAL visc_tol;
  static REAL visc_abs_tol;
  static int verbose;

public:

  Diffusion();
  Diffusion(Amr* Parent, AmrLevel* Caller, Diffusion* coarser,
	    int num_state, FluxRegister* Viscflux_reg,
	    MultiFab& Volume, MultiFab* Area,
            Array<int>  _is_diffusive,
            Array<REAL> _visc_coef);

  virtual ~Diffusion();

  void diffuse_scalar(REAL dt, int sigma, REAL be_cn_theta,
		      MultiFab* rho_half, int rho_flag,
                      int do_viscreflux=1,
                      MultiFab* delta_rhs=NULL, MultiFab* alpha=NULL, 
                      MultiFab** betan = NULL, 
                      MultiFab** betanp1=NULL);

  void diffuse_velocity(REAL dt, REAL be_cn_theta,
			MultiFab* rho_half, int rho_flag,
                        MultiFab*  delta_rhs = NULL, 
                        MultiFab** betan = NULL, 
                        MultiFab** betanp1 = NULL);

  void diffuse_velocity_constant_mu(REAL dt, REAL be_cn_theta,
			MultiFab* rho_half,
                        MultiFab* delta_rhs); 

  void diffuse_tensor_velocity(REAL dt, REAL be_cn_theta,
			MultiFab* rho_half,
                        MultiFab* delta_rhs, 
                        MultiFab** betan, 
                        MultiFab** betanp1);

  void diffuse_Vsync(MultiFab *Vsync, REAL dt, REAL be_cn_theta,
		     MultiFab* rho_half, int rho_flag,
                     MultiFab** beta=NULL);

  void diffuse_Vsync_constant_mu(MultiFab *Vsync, REAL dt, REAL be_cn_theta,
		     MultiFab* rho_half, int rho_flag);

  void diffuse_tensor_Vsync(MultiFab *Vsync, REAL dt, REAL be_cn_theta,
		     MultiFab* rho_half, int rho_flag,
                     MultiFab** beta);

  void diffuse_Ssync(MultiFab *Ssync, int sigma, REAL dt, REAL be_cn_theta,
		     MultiFab* rho_half, int rho_flag,
                     int do_viscsyncflux, 
                     MultiFab** beta=NULL, MultiFab* alpha=NULL);

  ABecLaplacian* getViscOp(int src_comp, REAL a, REAL b, 
			   Real time, ViscBndry & visc_bndry,
			   MultiFab* rho_half, int rho_flag, REAL* rhsscale=NULL,
                           MultiFab** beta=NULL, MultiFab* alpha_in=NULL);

  ABecLaplacian* getViscOp(int src_comp, REAL a, REAL b,
			   MultiFab* rho_half, int rho_flag, REAL* rhsscale=NULL,
                           MultiFab** beta=NULL, MultiFab* alpha_in=NULL);

#if defined(USE_TENSOR)
  DivVis* getTensorOp(REAL a, REAL b, 
			   REAL time, 
#if (BL_SPACEDIM==2) 
                           ViscBndry2D& visc_bndry,
#else 
                           ViscBndry3D& visc_bndry,
#endif
			   MultiFab* rho_half,
                           MultiFab** beta);
  void getTensorBndryData(
#if (BL_SPACEDIM==2) 
                          ViscBndry2D& bndry, 
#else 
                          ViscBndry3D& bndry, 
#endif
                          REAL time);
  DivVis* getTensorOp(REAL a, REAL b, 
			   MultiFab* rho_half,
                           MultiFab** beta);
#endif

  void getViscTerms(MultiFab& visc_terms, int src_comp, int comp,
		    REAL time, int rho_flag=0, MultiFab**beta=NULL);

  void getTensorViscTerms(MultiFab& visc_terms, 
			     REAL time, MultiFab** beta);

  void getBndryData(ViscBndry& bndry, int src_comp,
		    int num_comp, REAL time, int rho_flag=0);

  void FillBoundary(BndryRegister& bdry, int src_comp, int dest_comp,
		int num_comp, REAL time, int rho_flag=0);

  int useDvConstantMu() { return use_dv_constant_mu_def; }

  void checkBetas(MultiFab** beta1, MultiFab** beta2,
                  int& allthere, int& allnull);

  void checkBeta(MultiFab** beta,
                 int& allthere, int& allnull);

  void checkBeta(MultiFab** beta,
                 int& allthere);

  void allocFluxBoxesLevel(MultiFab**& fluxbox, 
                           int nghost=0, int nvar=1);

  void removeFluxBoxesLevel(MultiFab**& fluxbox);

#ifdef USE_NAVIERSTOKES
  void compute_divmusi(REAL time, MultiFab** beta, MultiFab& divmusi);

  void compute_divmusi(REAL time, REAL mu, MultiFab& divmusi);
#endif
};

#endif

