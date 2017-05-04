
#ifndef _PROJOUTFLOWBC_H_
#define _PROJOUTFLOWBC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Orientation.H>
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

    virtual void computeBC (amrex::FArrayBox         velFab[][2*BL_SPACEDIM],
			    amrex::FArrayBox         divuFab[2*BL_SPACEDIM],
			    amrex::FArrayBox         rhoFab[2*BL_SPACEDIM],
			    amrex::FArrayBox         phiFab[2*BL_SPACEDIM],
			    const amrex::Geometry&   geom,
			    amrex::Orientation*      outFaces,
			    int               numOutFlowFaces,
			    const int*        lo_bc,
			    const int*        hi_bc,
			    amrex::Real              gravity) override;

    virtual void computeBC (amrex::FArrayBox         velFab[][2*BL_SPACEDIM],
			    amrex::FArrayBox         divuFab[2*BL_SPACEDIM],
			    amrex::FArrayBox         rhoFab[2*BL_SPACEDIM],
			    amrex::FArrayBox         phiFab[2*BL_SPACEDIM],
			    const amrex::Geometry&   geom,
			    amrex::Orientation*      outFaces,
			    int               numOutFlowFaces,
			    const int*        lo_bc,
			    const int*        hi_bc,
			    amrex::Real              small_udiff,
			    amrex::Real              gravity) override;

    void computeRhoG (amrex::FArrayBox*         rhoFab,
                      amrex::FArrayBox*         phiFab,
                      const amrex::Geometry&    geom,
                      amrex::Orientation*       outFaces,
                      int                numOutFlowFaces,
                      amrex::Real               gravity,
                      const int*         lo_bc,
                      const int*         hi_bc);
protected:

    static void Initialize ();
    static void Finalize ();

#if (BL_SPACEDIM == 3)
    void computeCoefficients (amrex::FArrayBox&   rhs,
                              amrex::FArrayBox&   beta,
                              amrex::FArrayBox&   ccExt,
                              amrex::Box&         faceBox,
                              amrex::Real*        dxFiltered,
                              int*         isPeriodic);
    //
    // Static data.
    //
    static amrex::Real tol;
    static amrex::Real abs_tol;
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

    ProjOutFlowBC_MG (const amrex::Box& Domain,
                      amrex::FArrayBox* Phi,
                      amrex::FArrayBox* Rhs,
                      amrex::FArrayBox* Resid,
                      amrex::FArrayBox* Beta,
                      amrex::Real*      H,
                      int*       isPeriodic);

    virtual ~ProjOutFlowBC_MG () override;

    static int Verbose () { return verbose; }

    static int MaxIters () { return maxIters; }

protected:

    static void Initialize ();
    static void Finalize ();

    virtual void Restrict () override;
    virtual void interpolate () override;
    virtual amrex::Real residual () override;
    virtual void step (int nSteps) override;
    virtual void gsrb (int nSteps) override;
    //
    // Static data.
    //
    static int  verbose;
    static bool useCGbottomSolver;
    static amrex::Real cg_tol;
    static amrex::Real cg_abs_tol;
    static amrex::Real cg_max_jump;
    static int  cg_maxiter;
    static int  maxIters;
};
#endif

#endif


