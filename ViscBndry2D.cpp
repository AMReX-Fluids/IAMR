#if (BL_SPACEDIM==2) && defined (USE_TENSOR)

//
// $Id: ViscBndry2D.cpp,v 1.14 1999-02-03 17:59:36 lijewski Exp $
//

#include <LO_BCTYPES.H>
#include <ViscBndry2D.H>

void
ViscBndry2D::setBndryConds (const BCRec& bc,
                            int          ratio,
                            int          comp)
{
#if BL_SPACEDIM == 2
    assert(comp < 2*2);     // u and v, plus derivs of same
#elif BL_SPACEDIM == 3
    assert(comp < 3*(3+1)); // u and v, plus derivs of same
#endif

    const REAL* dx     = geom.CellSize();
    const BOX&  domain = geom.Domain();

    for (OrientationIter fi; fi; ++fi)
    {
        Array<REAL>&               bloc  = bcloc[fi()];
        Array< Array<BoundCond> >& bctag = bcond[fi()];

        int dir = fi().coordDir();
#if 0
	Real delta = dx[dir]*ratio[dir];
#else
	Real delta = dx[dir]*ratio;
#endif
        int p_bc = fi().isLow() ? bc.lo(dir): bc.hi(dir);

        for (int i = 0; i < boxes().length(); i++)
        {
            if (domain[fi()] == boxes()[i][fi()] && !geom.isPeriodic(dir))
            {
                //
                // All physical bc values are located on face.
                //
                if (p_bc == EXT_DIR )
                {
                    bctag[i][comp] = LO_DIRICHLET;
                    bloc[i] = 0.0;
                }
                else if (p_bc == FOEXTRAP || p_bc == HOEXTRAP || p_bc == REFLECT_EVEN)
                {
                    bctag[i][comp] = LO_NEUMANN;
                    bloc[i] = 0.0;
                }
                else if (p_bc == REFLECT_ODD)
                {
                    bctag[i][comp] = LO_REFLECT_ODD;
                    bloc[i] = 0.0;
                }
            }
            else
            {
                //
                // Internal bndry, distance is half of crse.
                //
                bctag[i][comp] = LO_DIRICHLET;
                bloc[i] = 0.5*delta;
            }
        }
    }
}

void
ViscBndry2D::setHomogValues (const Array<BCRec>& bc,
                             int                 ratio)
{
    for (int n = 0; n < bc.length(); ++n)
        setBndryConds(bc[n], ratio, n);

    for (OrientationIter fi; fi; ++fi)
    {
        for (FabSetIterator fsi(bndry[fi()]); fsi.isValid(); ++fsi)
        {
            fsi().setVal(0);
        }
    }
}

#endif /*(BL_SPACEDIM==2) && defined (USE_TENSOR)*/
