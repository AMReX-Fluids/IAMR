// BL_COPYRIGHT_NOTICE

#ifndef _PROJOUTFLOWBC_H_
#define _PROJOUTFLOWBC_H_

//
// $Id: ProjOutFlowBC.H,v 1.3 1999-08-20 17:58:55 propp Exp $
//

#include "FArrayBox.H"
#include "Geometry.H"
#include "Orientation.H"

enum Proj_Solver_Type {
  /// back substitution
  HG_BACK= 0, 
  /// multigrid
  HG_MG};

//
// This class solves the BL_SPACEDIM-1 dimensional problem for the holy-grail
// projection at the outflow face.  It first extrapolates the interior values
// to the outflow edge.  Then, it solves using either back substitution or
// a multigrid solver.  If a multigrid solver is used, the problem is shifted
// to so that the last dimension is 1.
//

class ProjOutFlowBC 
{
protected:
  Proj_Solver_Type proj_solver;  
  static Real tol;
  static Real abs_tol;

#if (BL_SPACEDIM == 2)
  void solveBackSubstitution(FArrayBox& phi,
			     FArrayBox& divuExt,
			     FArrayBox& uExt,
			     FArrayBox& rhoExt,
			     Array<Real>& rcen,
			     int r_lo, int r_hi,
			     int* isPeriodicFiltered,
			     Real* dxFiltered,
			     Box& faceBox,
			     const Orientation& outFace);
#endif

  void computeCoefficients(FArrayBox& rhs,
			   FArrayBox& beta,
			   FArrayBox& uExt,
			   FArrayBox& divuExt,
			   FArrayBox& rhoExt,
			   Array<Real>& rcen,
			   int r_lo, int r_hi,
			   Box& faceBox,
			   Real* dxFiltered,
			   int* isPeriodic);

public:

  ProjOutFlowBC();

  ~ProjOutFlowBC();

  void computeProjBC(FArrayBox& velFab, FArrayBox& divuFab,
		     FArrayBox& rhoFab, FArrayBox& phiFab,
		     const Geometry& geom, const Orientation& outFace);

};

//
// This class is a BL_SPACEDIM-1 multigrid solver for the holy grail 
// projection.
// 

class ProjOutFlowBC_MG 
{
protected:
  Box domain;

  FArrayBox* phi;
  FArrayBox* rhs;
  FArrayBox* resid;
  FArrayBox* cgwork;
  FArrayBox* beta;

  Real h[BL_SPACEDIM];
  int isPeriodic[BL_SPACEDIM];

  ProjOutFlowBC_MG* next;

  static int verbose;
  static bool useCGbottomSolver;
  static Real cg_tol;
  static Real cg_abs_tol;
  static Real cg_max_jump;
  static int cg_maxiter;
  static int maxIters;

  void Restrict();
  void interpolate();
  Real residual();
  void step(int nSteps);
  void gsrb(int nSteps);

public:
  ProjOutFlowBC_MG(const Box& Domain,
		 FArrayBox* Phi,
		 FArrayBox* Rhs,
		 FArrayBox* Resid,
		 FArrayBox* Beta,
		 Real* H,
		 int* isPeriodic);
  ~ProjOutFlowBC_MG();

  void solve(Real tol, Real abs_tol,int i1, int i2);
  Real vcycle(int i1, int i2);
};

#endif


