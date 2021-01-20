
#ifndef _MACOUTFLOWBC_H_
#define _MACOUTFLOWBC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Orientation.H>
#include <OutFlowBC.H>

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

};

#endif

