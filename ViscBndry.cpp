//BL_COPYRIGHT_NOTICE

//
// $Id: ViscBndry.cpp,v 1.14 1998-10-12 23:19:01 sstanley Exp $
//

#include <LO_BCTYPES.H>
#include <ViscBndry.H>

void
ViscBndry::setBndryConds (const BCRec&   bc,
                          /*const*/ IntVect& ratio)
{
    //
    //  NOTE: ALL BCLOC VALUES ARE NOW DEFINED AS A LENGTH IN PHYSICAL
    //        DIMENSIONS *RELATIVE* TO THE FACE, NOT IN ABSOLUTE PHYSICAL SPACE
    //
    const Real* dx    = geom.CellSize();
    const Box& domain = geom.Domain();

    int comp = 0;

    for (OrientationIter fi; fi; ++fi)
    {
        Array<Real>& bloc = bcloc[fi()];
        Array< Array<BoundCond> >& bctag = bcond[fi()];

        int dir          = fi().coordDir();
        const Real delta = dx[dir]*ratio[dir];
        int p_bc         = (fi().isLow() ? bc.lo(dir) : bc.hi(dir));

        for (int i = 0; i < boxes().length(); i++)
        {
            if (domain[fi()] == boxes()[i][fi()] && !geom.isPeriodic(dir))
            {
                //
                // All physical bc values are located on face.
                //
                if (p_bc == EXT_DIR)
                {
                    bctag[i][comp] = LO_DIRICHLET;
                    bloc[i] = 0.;
                }
                else if (p_bc == FOEXTRAP      ||
                         p_bc == HOEXTRAP      || 
                         p_bc == REFLECT_EVEN)
                {
                    bctag[i][comp] = LO_NEUMANN;
                    bloc[i] = 0.;
                }
                else if (p_bc == REFLECT_ODD)
                {
                    bctag[i][comp] = LO_REFLECT_ODD;
                    bloc[i] = 0.;
                }
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

void
ViscBndry::setHomogValues (const BCRec& bc,
                           /*const*/ IntVect& ratio)
{
    setBndryConds(bc, ratio);

    for (OrientationIter fi; fi; ++fi)
    {
        for (FabSetIterator fsi(bndry[fi()]); fsi.isValid(); ++fsi)
        {
            fsi().setVal(0.);
        }
    }
}
