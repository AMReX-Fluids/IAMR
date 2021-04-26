
#ifndef _OUTFLOWBC_H_
#define _OUTFLOWBC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Orientation.H>
#include <AMReX_BCRec.H>

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

#endif /*_OUTFLOWBC_H_*/

