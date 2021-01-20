
#ifndef _PROJOUTFLOWBC_H_
#define _PROJOUTFLOWBC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Orientation.H>
#include <OutFlowBC.H>

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

};

#endif


