
#ifndef _PROJOUTFLOWBC_H_
#define _PROJOUTFLOWBC_H_

#include <FArrayBox.H>
#include <Geometry.H>
#include <Orientation.H>
#include <OutFlowBC.H>

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

    virtual void computeBC (FArrayBox         velFab[][2*BL_SPACEDIM],
			    FArrayBox         divuFab[2*BL_SPACEDIM],
			    FArrayBox         rhoFab[2*BL_SPACEDIM],
			    FArrayBox         phiFab[2*BL_SPACEDIM],
			    const Geometry&   geom,
			    Orientation*      outFaces,
			    int               numOutFlowFaces,
			    const int*        lo_bc,
			    const int*        hi_bc,
			    Real              gravity) override;

    virtual void computeBC (FArrayBox         velFab[][2*BL_SPACEDIM],
			    FArrayBox         divuFab[2*BL_SPACEDIM],
			    FArrayBox         rhoFab[2*BL_SPACEDIM],
			    FArrayBox         phiFab[2*BL_SPACEDIM],
			    const Geometry&   geom,
			    Orientation*      outFaces,
			    int               numOutFlowFaces,
			    const int*        lo_bc,
			    const int*        hi_bc,
			    Real              small_udiff,
			    Real              gravity) override;

    void computeRhoG (FArrayBox*         rhoFab,
                      FArrayBox*         phiFab,
                      const Geometry&    geom,
                      Orientation*       outFaces,
                      int                numOutFlowFaces,
                      Real               gravity,
                      const int*         lo_bc,
                      const int*         hi_bc);
protected:

    static void Initialize ();
    static void Finalize ();

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

    virtual ~ProjOutFlowBC_MG () override;

    static int Verbose () { return verbose; }

    static int MaxIters () { return maxIters; }

protected:

    static void Initialize ();
    static void Finalize ();

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


