#ifndef _PROJOUTFLOWBC_H_
#define _PROJOUTFLOWBC_H_

//
// $Id: ProjOutFlowBC.H,v 1.7 2003-02-19 19:15:46 almgren Exp $
//

#include "FArrayBox.H"
#include "Geometry.H"
#include "Orientation.H"
#include "OutFlowBC.H"

//
// This class solves the BL_SPACEDIM-1 dimensional problem for the holy-grail
// projection at the outflow face.  It first extrapolates the interior values
// to the outflow edge.  Then, it solves using either back substitution or
// a multigrid solver.  If a multigrid solver is used, the problem is shifted
// so that the last dimension is 1.
//

class ProjOutFlowBC 
    :
    public OutFlowBC
{
public:

    ProjOutFlowBC ();

    void computeBC (MultiFab*         velFab,
                    MultiFab&         divuFab,
                    MultiFab&         rhoFab,
                    MultiFab&         phiFab,
                    const Geometry&   geom,
                    Orientation*      outFaces,
                    int               numOutFlowFaces,
                    Real              gravity = 0.0);

    void computeRhoG (MultiFab&          rhoFab,
                      MultiFab&          phiFab,
                      const Geometry&    geom,
                      Orientation*       outFaces,
                      int                numOutFlowFaces,
                      Real               gravity);

protected:

#if (BL_SPACEDIM == 3)
    void computeCoefficients (FArrayBox&   rhs,
                              FArrayBox&   beta,
                              FArrayBox&   ccExt,
                              Box&         faceBox,
                              Real*        dxFiltered,
                              int*         isPeriodic);
    //
    // Static data.
    //
    static Real tol;
    static Real abs_tol;
#endif
};

#if (BL_SPACEDIM == 3)
//
// This class is a BL_SPACEDIM-1 multigrid solver for the holy grail 
// projection.
// 

class ProjOutFlowBC_MG 
    :
    public OutFlowBC_MG
{
public:

    ProjOutFlowBC_MG (const Box& Domain,
                      FArrayBox* Phi,
                      FArrayBox* Rhs,
                      FArrayBox* Resid,
                      FArrayBox* Beta,
                      Real*      H,
                      int*       isPeriodic);

    ~ProjOutFlowBC_MG ();

    static int Verbose () { return verbose; }

    static int MaxIters () { return maxIters; }

protected:

    void Restrict ();
    void interpolate ();
    Real residual ();
    void step (int nSteps);
    void gsrb (int nSteps);
    //
    // Static data.
    //
    static int  verbose;
    static bool useCGbottomSolver;
    static Real cg_tol;
    static Real cg_abs_tol;
    static Real cg_max_jump;
    static int  cg_maxiter;
    static int  maxIters;
};
#endif

#endif


