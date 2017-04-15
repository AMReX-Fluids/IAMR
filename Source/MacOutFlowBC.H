
#ifndef _MACOUTFLOWBC_H_
#define _MACOUTFLOWBC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Orientation.H>
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

    virtual void computeBC (amrex::FArrayBox         velMF[][2*BL_SPACEDIM],
			    amrex::FArrayBox         divuMF[2*BL_SPACEDIM],
			    amrex::FArrayBox         rhoMF[2*BL_SPACEDIM],
			    amrex::FArrayBox         phiMF[2*BL_SPACEDIM],
			    const amrex::Geometry&   geom,
			    amrex::Orientation*      outFaces,
			    int               numOutFlowFaces,
			    const int*        lo_bc,
			    const int*        hi_bc,
			    amrex::Real              gravity) override;

    virtual void computeBC (amrex::FArrayBox         velMF[][2*BL_SPACEDIM],
			    amrex::FArrayBox         divuMF[2*BL_SPACEDIM],
			    amrex::FArrayBox         rhoMF[2*BL_SPACEDIM],
			    amrex::FArrayBox         phiMF[2*BL_SPACEDIM],
			    const amrex::Geometry&   geom,
			    amrex::Orientation*      outFaces,
			    int               numOutFlowFaces,
			    const int*        lo_bc,
			    const int*        hi_bc,
			    amrex::Real              small_udiff,
			    amrex::Real              gravity) override;

    static void Initialize ();
    static void Finalize ();

protected:

#if (BL_SPACEDIM == 3)
    void computeCoefficients (amrex::FArrayBox&   rhs,
                              amrex::FArrayBox*   beta,
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
// This class is a BL_SPACEDIM-1 multigrid solver for the mac projection.
// 

class MacOutFlowBC_MG
    :
    public OutFlowBC_MG
{
public:

    MacOutFlowBC_MG (amrex::Box& Domain,
                     amrex::FArrayBox* Phi,
                     amrex::FArrayBox* Rhs,
                     amrex::FArrayBox* Resid,
                     amrex::FArrayBox* Beta,
                     amrex::Real*      H,
                     int*       IsPeriodic);

    virtual ~MacOutFlowBC_MG ()  override;

    static int Verbose ();

    static int MaxIters ();

    static void Initialize ();
    static void Finalize ();

protected:

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

