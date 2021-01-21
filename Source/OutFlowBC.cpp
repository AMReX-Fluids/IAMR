

#include <iostream>
#include <algorithm>

#include <AMReX_Print.H>
#include <OutFlowBC.H>

using namespace amrex;

OutFlowBC::OutFlowBC () {}

OutFlowBC::~OutFlowBC () {}

Box
OutFlowBC::SemiGrow (const Box& baseBox,
                         int        nGrow,
                         int        direction)
{
    IntVect grow_factor(D_DECL(nGrow,nGrow,nGrow));
    grow_factor[direction] = 0;
    return amrex::grow(baseBox,grow_factor);
}

Box
OutFlowBC::SemiCoarsen (const Box& baseBox,
                        int        ref_factor,
                        int        direction)
{
    IntVect ref_ratio(D_DECL(ref_factor,ref_factor,ref_factor));
    ref_ratio[direction] = 1;
    return amrex::coarsen(baseBox,ref_ratio);
}

void
OutFlowBC::GetOutFlowFaces (bool&        haveOutFlow,
                            Orientation* outFaces,
                            BCRec*       _phys_bc,
                            int&        numOutFlowBC)
{
    haveOutFlow = false;

    numOutFlowBC = 0;

    for (int idir = 0; idir < BL_SPACEDIM; idir++)
    {
        if (_phys_bc->lo(idir) == Outflow)
        {
            haveOutFlow = true;
            outFaces[numOutFlowBC] = Orientation(idir,Orientation::low);
            numOutFlowBC++;
        }

        if (_phys_bc->hi(idir) == Outflow)
        {
            haveOutFlow = true;
            outFaces[numOutFlowBC] = Orientation(idir,Orientation::high);
            numOutFlowBC++;
        }
    }
}

bool
OutFlowBC::HasOutFlowBC (BCRec* _phys_bc)
{
    bool has_out_flow = false;
    int  numOutFlowBC = 0;

    for (int idir = 0; idir < BL_SPACEDIM; idir++)
    {
        if (_phys_bc->lo(idir) == Outflow)
        {
            has_out_flow = true;
            numOutFlowBC++;
        }

        if (_phys_bc->hi(idir) == Outflow)
        {
            has_out_flow = true;
            numOutFlowBC++;
        }
    }

    return has_out_flow;
}
