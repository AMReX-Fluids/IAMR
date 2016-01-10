#include <FluxBoxes.H>

MultiFab**
FluxBoxes::define (const AmrLevel* amr_level, int nvar, int nghost)
{
    BL_ASSERT(data == 0);
    data = new MultiFab*[BL_SPACEDIM];
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
	const BoxArray& ba = amr_level->getEdgeBoxArray(dir);
        data[dir] = new MultiFab(ba,nvar,nghost);
    }
    return data;
}

void
FluxBoxes::clear ()
{
    if (data != 0)
    {
        for (int i = 0; i<BL_SPACEDIM; i++)
            delete data[i];
        delete [] data;
        data = 0;
    }
}
