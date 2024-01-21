

#include <iostream>
#include <algorithm>

#include <AMReX_Print.H>
#include <OutFlowBC.H>

using namespace amrex;

Box
OutFlowBC::SemiGrow (const Box& baseBox,
                         int        nGrow,
                         int        direction)
{
    IntVect grow_factor(AMREX_D_DECL(nGrow,nGrow,nGrow));
    grow_factor[direction] = 0;
    return amrex::grow(baseBox,grow_factor);
}

Box
OutFlowBC::SemiCoarsen (const Box& baseBox,
                        int        ref_factor,
                        int        direction)
{
    IntVect ref_ratio(AMREX_D_DECL(ref_factor,ref_factor,ref_factor));
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

    for (int idir = 0; idir < AMREX_SPACEDIM; idir++)
    {
        if (_phys_bc->lo(idir) == PhysBCType::outflow)
        {
            haveOutFlow = true;
            outFaces[numOutFlowBC] = Orientation(idir,Orientation::low);
            numOutFlowBC++;
        }

        if (_phys_bc->hi(idir) == PhysBCType::outflow)
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

    for (int idir = 0; idir < AMREX_SPACEDIM; idir++)
    {
        if (_phys_bc->lo(idir) == PhysBCType::outflow)
        {
            has_out_flow = true;
        }

        if (_phys_bc->hi(idir) == PhysBCType::outflow)
        {
            has_out_flow = true;
        }
    }

    return has_out_flow;
}
