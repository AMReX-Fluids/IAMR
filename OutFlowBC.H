#ifndef _OUTFLOWBC_H_
#define _OUTFLOWBC_H_

//
// $Id: OutFlowBC.H,v 1.7 2003-09-11 21:10:14 almgren Exp $
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

    OutFlowBC ();

    virtual ~OutFlowBC ();

    virtual void computeBC (FArrayBox         velMF[][2*BL_SPACEDIM],
                            FArrayBox         divuMF[2*BL_SPACEDIM],
                            FArrayBox         rhoMF[2*BL_SPACEDIM],
                            FArrayBox         phiMF[2*BL_SPACEDIM],
                            const Geometry&    geom,
                            Orientation*       outFaces,
                            int                numOutFlowFaces,
                            const int*        lo_bc,
                            const int*        hi_bc,
                            Real               gravity = 0.0) = 0;

    static Box SemiGrow (const Box& baseBox,
                         int        nGrow,
                         int        direction);

    static Box SemiCoarsen (const Box& baseBox,
                            int        ref_factor,
                            int        direction);

    static void GetOutFlowFaces (bool&        haveOutFlow,
                                 Orientation* outFaces,
                                 BCRec*       _phys_bc,
                                 int&         numOutFlowFaces);

    static bool HasOutFlowBC (BCRec* _phys_bc);

protected:
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

    virtual ~OutFlowBC_MG ();

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

    virtual void Restrict ()       = 0;
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

