#ifndef _MACOUTFLOWBC_H_
#define _MACOUTFLOWBC_H_

//
// $Id: MacOutFlowBC.H,v 1.5 2000-11-02 18:08:01 lijewski Exp $
//

#include "FArrayBox.H"
#include "Geometry.H"
#include "Orientation.H"
#include "OutFlowBC.H"

//
// This class solves the BL_SPACEDIM-1 dimensional problem for the mac
// projection at the outflow face.  It first extrapolates the interior values
// to the outflow edge.  Then, it solves using either back substitution or
// a multigrid solver.  If a multigrid solver is used, the problem is shifted
// so that the last dimension is 1.
//

class MacOutFlowBC
    :
    public OutFlowBC
{
public:

    MacOutFlowBC ();
  
    void computeBC (FArrayBox*         velFab,
                    FArrayBox&         divuFab,
                    FArrayBox&         rhoFab,
                    FArrayBox&         phiFab,
                    const Geometry&    geom,
                    const Orientation& outFace);
protected:

#if (BL_SPACEDIM == 2)
    void solveBackSubstitution (FArrayBox&         phi,
                                FArrayBox&         divuExt,
                                FArrayBox*         uExt,
                                FArrayBox&         rhoExt,
                                Array<Real>&       rcen,
                                int                r_lo,
                                int                r_hi,
                                int*               isPeriodicFiltered,
                                Real*              dxFiltered,
                                Box&               faceBox,
                                const Orientation& outFace);
#endif

    void computeCoefficients (FArrayBox&   rhs,
                              FArrayBox*   beta,
                              FArrayBox*   uExt,
                              FArrayBox&   divuExt,
                              FArrayBox&   rhoExt,
                              Array<Real>& rcen,
                              int          r_lo,
                              int          r_hi,
                              Box&         faceBox,
                              Real*        dxFiltered,
                              int*         isPeriodic);
    //
    // Static data.
    //
    static Real tol;
    static Real abs_tol;
};

//
// This class is a BL_SPACEDIM-1 multigrid solver for the mac projection.
// 

class MacOutFlowBC_MG
    :
    public OutFlowBC_MG
{
public:

    MacOutFlowBC_MG (Box& Domain,
                     FArrayBox* Phi,
                     FArrayBox* Rhs,
                     FArrayBox* Resid,
                     FArrayBox* Beta,
                     Real*      H,
                     int*       IsPeriodic);

    ~MacOutFlowBC_MG ();

    static int Verbose () { return verbose; }

    static int MaxIters () { return maxIters; }

protected:

    void restrict ();
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

