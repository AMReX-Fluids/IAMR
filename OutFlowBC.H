#ifndef _OUTFLOWBC_H_
#define _OUTFLOWBC_H_

//
// $Id: OutFlowBC.H,v 1.1 2000-11-02 18:08:02 lijewski Exp $
//

#include "FArrayBox.H"
#include "Geometry.H"
#include "Orientation.H"
#include "BCRec.H"

enum OutFlow_Solver_Type
{
    BC_BACK, // back substitution
    BC_MG    // multigrid
};

//
// The baseclass for MacOutFlowBC and ProjOutFlowBC.
//

class OutFlowBC
{
public:

    virtual void computeBC (FArrayBox*         velFab,
                            FArrayBox&         divuFab,
                            FArrayBox&         rhoFab,
                            FArrayBox&         phiFab,
                            const Geometry&    geom,
                            const Orientation& outFace) = 0;

    static Box SemiGrow (const Box& baseBox,
                         int        nGrow,
                         int        direction);

    static Box SemiCoarsen (const Box& baseBox,
                            int        ref_factor,
                            int        direction);

    static void GetOutFlowFace (bool&        haveOutFlow,
                                Orientation& outFace,
                                BCRec*       _phys_bc);

    static bool HasOutFlowBC (BCRec* _phys_bc);

protected:
    //
    // Per-object data.
    //
    OutFlow_Solver_Type solver;
};

//
// The baseclass for MacOutFlowBC_MG and ProjOutFlowBC_MG
//

class OutFlowBC_MG 
{
public:

    OutFlowBC_MG (const Box& Domain,
                  FArrayBox* Phi,
                  FArrayBox* Rhs,
                  FArrayBox* Resid,
                  FArrayBox* Beta,
                  Real*      H,
                  int*       IsPeriodic,
                  bool       is_scalar);

    ~OutFlowBC_MG ();

    void solve (Real tol,
                Real abs_tol,
                int  i1,
                int  i2,
                int  maxIters,
                int  verbose);

    Real vcycle (int i1, int i2);

    const Box& theDomain () const { return domain; }

    FArrayBox* theRhs () { return rhs; }

    FArrayBox* thePhi () { return phi; }

protected:

    virtual void restrict ()       = 0;
    virtual void interpolate ()    = 0;
    virtual Real residual ()       = 0;
    virtual void step (int nSteps) = 0;
    virtual void gsrb (int nSteps) = 0;
    //
    // The per-object data
    //
    Box           domain;
    FArrayBox*    phi;
    FArrayBox*    rhs;
    FArrayBox*    resid;
    FArrayBox*    cgwork;
    FArrayBox*    beta;
    Real          h[BL_SPACEDIM];
    int           isPeriodic[BL_SPACEDIM];
    OutFlowBC_MG* next;
    bool          beta_is_scalar;
};

#endif /*_OUTFLOWBC_H_*/

