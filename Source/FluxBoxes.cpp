#include <FluxBoxes.H>

using namespace amrex;

MultiFab**
FluxBoxes::define (const AmrLevel* amr_level, int nvar, int nghost)
{
    AMREX_ASSERT(data == nullptr);
    data = new MultiFab*[AMREX_SPACEDIM];
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        const BoxArray& ba = amr_level->getEdgeBoxArray(dir);
        const DistributionMapping& dm = amr_level->DistributionMap();
        data[dir] = new MultiFab(ba,dm,nvar,nghost,MFInfo(),amr_level->Factory());
    }
    return data;
}

void
FluxBoxes::clear ()
{
    if (data != nullptr)
    {
        for (int i = 0; i<AMREX_SPACEDIM; i++) {
            delete data[i];
        }
        delete [] data;
        data = nullptr;
    }
}
