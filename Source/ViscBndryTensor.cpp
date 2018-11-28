

#include <AMReX_LO_BCTYPES.H>
#include <ViscBndryTensor.H>

using namespace amrex;

void
ViscBndryTensor::setBndryConds (const BCRec& bc,
                                int          ratio,
                                int          comp)
{
    BL_ASSERT (comp < MCLinOp::bcComponentsNeeded());

    const Real* dx     = geom.CellSize();
    const Box&  domain = geom.Domain();
    //
    // We note that all orientations of the FabSets have the same distribution.
    // We'll use the low 0 side as the model.
    //
    for (FabSetIter fsi(bndry[Orientation(0,Orientation::low)]); fsi.isValid(); ++fsi)
    {
        const int                  i     = fsi.index();
        RealTuple&                 bloc  = bcloc[i];
        Vector< Vector<BoundCond> >& bctag = bcond[i];

        for (OrientationIter fi; fi; ++fi)
        {
            const Orientation face  = fi();
            const int         dir   = face.coordDir();

            if (domain[face] == boxes()[i][face] && !geom.isPeriodic(dir))
            {
                //
                // All physical bc values are located on face.
                //
                const int p_bc  = face.isLow() ? bc.lo(dir): bc.hi(dir);

                if (p_bc == EXT_DIR )
                {
                    bctag[face][comp] = LO_DIRICHLET;
                    bloc[face]        = 0;
                }
                else if (p_bc == FOEXTRAP || p_bc == HOEXTRAP || p_bc == REFLECT_EVEN)
                {
                    bctag[face][comp] = LO_NEUMANN;
                    bloc[face]        = 0;
                }
                else if (p_bc == REFLECT_ODD)
                {
                    bctag[face][comp] = LO_REFLECT_ODD;
                    bloc[face]        = 0;
                }
            }
            else
            {
                //
                // Internal bndry, distance is half of crse.
                //
                const Real delta = dx[dir]*ratio;

                bctag[face][comp] = LO_DIRICHLET;
                bloc[face]        = 0.5*delta;
            }
        }
    }
}

void
ViscBndryTensor::setHomogValues (const Vector<BCRec>& bc,
                                 int                 ratio)
{
    for (int n = 0; n < bc.size(); ++n)
        setBndryConds(bc[n], ratio, n);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (OrientationIter fi; fi; ++fi)
    {
        const Orientation face  = fi();

        for (FabSetIter fsi(bndry[face]); fsi.isValid(); ++fsi)
        {
            bndry[face][fsi].setVal(0);
        }
    }
}
