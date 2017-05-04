
#ifndef _OUTFLOWBC_H_
#define _OUTFLOWBC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Orientation.H>
#include <AMReX_BCRec.H>

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

    virtual void computeBC (amrex::FArrayBox         velMF[][2*BL_SPACEDIM],
                            amrex::FArrayBox         divuMF[2*BL_SPACEDIM],
                            amrex::FArrayBox         rhoMF[2*BL_SPACEDIM],
                            amrex::FArrayBox         phiMF[2*BL_SPACEDIM],
                            const amrex::Geometry&   geom,
                            amrex::Orientation*      outFaces,
                            int               numOutFlowFaces,
                            const int*        lo_bc,
                            const int*        hi_bc,
                            amrex::Real              gravity) = 0;

    virtual void computeBC (amrex::FArrayBox         velMF[][2*BL_SPACEDIM],
                            amrex::FArrayBox         divuMF[2*BL_SPACEDIM],
                            amrex::FArrayBox         rhoMF[2*BL_SPACEDIM],
                            amrex::FArrayBox         phiMF[2*BL_SPACEDIM],
                            const amrex::Geometry&    geom,
                            amrex::Orientation*       outFaces,
                            int                numOutFlowFaces,
                            const int*        lo_bc,
                            const int*        hi_bc,
                            amrex::Real              small_udiff, 
                            amrex::Real              gravity) = 0;

    static amrex::Box SemiGrow (const amrex::Box& baseBox,
                         int        nGrow,
                         int        direction);

    static amrex::Box SemiCoarsen (const amrex::Box& baseBox,
                            int        ref_factor,
                            int        direction);

    static void GetOutFlowFaces (bool&        haveOutFlow,
                                 amrex::Orientation* outFaces,
                                 amrex::BCRec*       _phys_bc,
                                 int&         numOutFlowFaces);

    static bool HasOutFlowBC (amrex::BCRec* _phys_bc);

protected:
};

//
// The baseclass for MacOutFlowBC_MG and ProjOutFlowBC_MG
//

class OutFlowBC_MG 
{
public:

    OutFlowBC_MG (const amrex::Box& Domain,
                  amrex::FArrayBox* Phi,
                  amrex::FArrayBox* Rhs,
                  amrex::FArrayBox* Resid,
                  amrex::FArrayBox* Beta,
                  amrex::Real*      H,
                  int*       IsPeriodic,
                  bool       is_scalar);

    virtual ~OutFlowBC_MG ();

    void solve (amrex::Real tol,
                amrex::Real abs_tol,
                int  i1,
                int  i2,
                int  maxIters,
                int  verbose);

    amrex::Real vcycle (int i1, int i2);

    const amrex::Box& theDomain () const { return domain; }

    amrex::FArrayBox* theRhs () { return rhs; }

    amrex::FArrayBox* thePhi () { return phi; }

protected:

    virtual void Restrict ()       = 0;
    virtual void interpolate ()    = 0;
    virtual amrex::Real residual ()       = 0;
    virtual void step (int nSteps) = 0;
    virtual void gsrb (int nSteps) = 0;
    //
    // The per-object data
    //
    amrex::Box           domain;
    amrex::FArrayBox*    phi;
    amrex::FArrayBox*    rhs;
    amrex::FArrayBox*    resid;
    amrex::FArrayBox     cgwork;
    amrex::FArrayBox*    beta;
    amrex::Real          h[BL_SPACEDIM];
    int           isPeriodic[BL_SPACEDIM];
    OutFlowBC_MG* next;
    bool          beta_is_scalar;
};

#endif /*_OUTFLOWBC_H_*/

