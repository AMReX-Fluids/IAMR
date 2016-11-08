
#ifndef _MACOUTFLOWBC_H_
#define _MACOUTFLOWBC_H_

#include <FArrayBox.H>
#include <Geometry.H>
#include <Orientation.H>
#include <OutFlowBC.H>

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

    virtual void computeBC (FArrayBox         velMF[][2*BL_SPACEDIM],
			    FArrayBox         divuMF[2*BL_SPACEDIM],
			    FArrayBox         rhoMF[2*BL_SPACEDIM],
			    FArrayBox         phiMF[2*BL_SPACEDIM],
			    const Geometry&   geom,
			    Orientation*      outFaces,
			    int               numOutFlowFaces,
			    const int*        lo_bc,
			    const int*        hi_bc,
			    Real              gravity) override;

    virtual void computeBC (FArrayBox         velMF[][2*BL_SPACEDIM],
			    FArrayBox         divuMF[2*BL_SPACEDIM],
			    FArrayBox         rhoMF[2*BL_SPACEDIM],
			    FArrayBox         phiMF[2*BL_SPACEDIM],
			    const Geometry&   geom,
			    Orientation*      outFaces,
			    int               numOutFlowFaces,
			    const int*        lo_bc,
			    const int*        hi_bc,
			    Real              small_udiff,
			    Real              gravity) override;

    static void Initialize ();
    static void Finalize ();

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

    virtual ~MacOutFlowBC_MG ()  override;

    static int Verbose ();

    static int MaxIters ();

    static void Initialize ();
    static void Finalize ();

protected:

    virtual void Restrict () override;
    virtual void interpolate () override;
    virtual Real residual () override;
    virtual void step (int nSteps) override;
    virtual void gsrb (int nSteps) override;
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

