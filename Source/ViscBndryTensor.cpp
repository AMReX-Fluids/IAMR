//
// $Id: ViscBndryTensor.cpp,v 1.7 2001-08-09 22:42:00 marc Exp $
//
#include <winstd.H>

#include <LO_BCTYPES.H>
#include <ViscBndryTensor.H>

void
ViscBndryTensor::setBndryConds (const BCRec& bc,
                                int          ratio,
                                int          comp)
{
    BL_ASSERT (comp < MCLinOp::bcComponentsNeeded());

    const Real* dx     = geom.CellSize();
    const Box&  domain = geom.Domain();

    for (OrientationIter fi; fi; ++fi)
    {
        Array<Real>&               bloc  = bcloc[fi()];
        Array< Array<BoundCond> >& bctag = bcond[fi()];

        int dir = fi().coordDir();
	Real delta = dx[dir]*ratio;
        int p_bc = fi().isLow() ? bc.lo(dir): bc.hi(dir);

        for (int i = 0; i < boxes().size(); i++)
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
                else if (p_bc == FOEXTRAP || p_bc == HOEXTRAP 
                                          || p_bc == REFLECT_EVEN)
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
ViscBndryTensor::setHomogValues (const Array<BCRec>& bc,
                                 int                 ratio)
{
    for (int n = 0; n < bc.size(); ++n)
        setBndryConds(bc[n], ratio, n);

    for (OrientationIter fi; fi; ++fi)
    {
        for (FabSetIter fsi(bndry[fi()]); fsi.isValid(); ++fsi)
        {
            bndry[fi()][fsi].setVal(0);
        }
    }
}
