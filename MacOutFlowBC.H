#ifndef _MACOUTFLOWBC_H_
#define _MACOUTFLOWBC_H_

//
// $Id: MacOutFlowBC.H,v 1.9 2003-02-19 21:42:35 almgren Exp $
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

    void computeBC (MultiFab*         velMF,
                    MultiFab&         divuMF,
                    MultiFab&         rhoMF,
                    MultiFab&         phiMF,
                    const Geometry&   geom,
                    Orientation*      outFaces,
                    int               numOutFlowFaces,
                    Real              gravity = 0.0);

protected:

#if (BL_SPACEDIM == 3)
    void computeCoefficients (FArrayBox&   rhs,
                              FArrayBox*   beta,
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

    static int Verbose ();

    static int MaxIters ();

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

