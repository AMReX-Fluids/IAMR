
#include <winstd.H>

#include <LO_BCTYPES.H>
#include <ViscBndry.H>

void
ViscBndry::setBndryConds (const BCRec& phys_bc,
                          int          ratio)
{
    IntVect ratio_vect = ratio * IntVect::TheUnitVector();
    setBndryConds(phys_bc, ratio_vect);
}

void
ViscBndry::setBndryConds (const BCRec&   bc,
                          /*const*/ IntVect& ratio,
			  int comp)
{
    //
    //  NOTE: ALL BCLOC VALUES ARE NOW DEFINED AS A LENGTH IN PHYSICAL
    //        DIMENSIONS *RELATIVE* TO THE FACE, NOT IN ABSOLUTE PHYSICAL SPACE
    //
    const Real* dx    = geom.CellSize();
    const Box& domain = geom.Domain();

    for (OrientationIter fi; fi; ++fi)
    {
        const int  dir   = fi().coordDir();
        const Real delta = dx[dir]*ratio[dir];
        const int  p_bc  = (fi().isLow() ? bc.lo(dir) : bc.hi(dir));

        for (FabSetIter fsi(bndry[fi()]); fsi.isValid(); ++fsi)
        {
            const int i = fsi.index();

            Real& bloc = bcloc[fi()][i];

            Array<BoundCond>& bctag = bcond[fi()][i];

            if (domain[fi()] == boxes()[i][fi()] && !geom.isPeriodic(dir))
            {
                //
                // All physical bc values are located on face.
                //
                if (p_bc == EXT_DIR)
                {
                    bctag[comp] = LO_DIRICHLET;
                    bloc        = 0;
                }
                else if (p_bc == FOEXTRAP      ||
                         p_bc == HOEXTRAP      || 
                         p_bc == REFLECT_EVEN)
                {
                    bctag[comp] = LO_NEUMANN;
                    bloc        = 0;
                }
                else if (p_bc == REFLECT_ODD)
                {
                    bctag[comp] = LO_REFLECT_ODD;
                    bloc        = 0;
                }
            }
            else
            {
                //
                // Internal bndry.
                //
                bctag[comp] = LO_DIRICHLET;
                bloc        = 0.5*delta;
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
        for (FabSetIter fsi(bndry[fi()]); fsi.isValid(); ++fsi)
        {
            bndry[fi()][fsi].setVal(0);
        }
    }
}
