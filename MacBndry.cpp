//
// $Id: MacBndry.cpp,v 1.15 2001-08-22 16:42:00 car Exp $
//
#include <winstd.H>

#include <LO_BCTYPES.H>
#include <MacBndry.H>

MacBndry::MacBndry ()
    :
    InterpBndryData()
{}

MacBndry::MacBndry (const BoxArray& _grids,
                    int             _ncomp,
                    const Geometry& _geom)
    :
    InterpBndryData(_grids,_ncomp,_geom)
{}

void
MacBndry::setBndryConds (const BCRec& phys_bc,
                         IntVect&     ratio,
			 int          comp)
{
    //
    // ALL BCLOC VALUES ARE NOW DEFINED AS A LENGTH IN PHYSICAL
    // DIMENSIONS *RELATIVE* TO THE FACE, NOT IN ABSOLUTE PHYSICAL SPACE
    //
    const BoxArray& grids      = boxes();
    const int ngrds            = grids.size();
    const Real* dx             = geom.CellSize();
    const Box& domain          = geom.Domain();

    for (OrientationIter fi; fi; ++fi)
    {
        Array<Real>& bloc                = bcloc[fi()];
        Array< Array<BoundCond> >& bctag = bcond[fi()];

        int dir    = fi().coordDir();
        Real delta = dx[dir]*ratio[dir];
        int p_bc   = (fi().isLow() ? phys_bc.lo(dir) : phys_bc.hi(dir));

        for (int i = 0; i < ngrds; i++)
        {
            const Box& grd = grids[i];

            if (domain[fi()] == grd[fi()] && !geom.isPeriodic(dir))
            {
                //
                // All physical bc values are located on face.
                //
                bctag[i][comp] = (p_bc == Outflow) ? LO_DIRICHLET : LO_NEUMANN;
                bloc[i]        = 0;
            }
            else
            {
                //
                // Internal bndry.
                //
                bctag[i][comp] = LO_DIRICHLET;
		bloc[i] = 0.5*delta;
            }
        }
    }
}
