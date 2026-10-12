
#ifndef IAMR_OUTFLOWBC_H_
#define IAMR_OUTFLOWBC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Orientation.H>
#include <AMReX_BCRec.H>

//
// Queries on the outflow faces of the domain.
//

class OutFlowBC
{
public:
    static void GetOutFlowFaces (bool&               haveOutFlow,
                                 amrex::Orientation* outFaces,
                                 amrex::BCRec*       _phys_bc,
                                 int&                numOutFlowBC);

    static bool HasOutFlowBC (amrex::BCRec* _phys_bc);
};

#endif /*_OUTFLOWBC_H_*/
