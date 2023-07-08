
#ifndef IAMR_OUTFLOWBC_H_
#define IAMR_OUTFLOWBC_H_

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
    static amrex::Box SemiGrow (const amrex::Box& baseBox,
                                int               nGrow,
                                int               direction);

    static amrex::Box SemiCoarsen (const amrex::Box& baseBox,
                                   int               ref_factor,
                                   int               direction);

    static void GetOutFlowFaces (bool&               haveOutFlow,
                                 amrex::Orientation* outFaces,
                                 amrex::BCRec*       _phys_bc,
                                 int&                numOutFlowBC);

    static bool HasOutFlowBC (amrex::BCRec* _phys_bc);
};

#endif /*_OUTFLOWBC_H_*/
