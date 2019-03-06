

#include <AMReX_LO_BCTYPES.H>
#include <ViscBndry.H>

using namespace amrex;

ViscBndry::ViscBndry ()
    :
    InterpBndryData()
{}

ViscBndry::ViscBndry (const BoxArray& _grids,
                      const DistributionMapping& _dmap,
                      int             _ncomp,
                      const Geometry& _geom)
    :
    InterpBndryData(_grids,_dmap,_ncomp,_geom)
{}

ViscBndry::~ViscBndry () {}

void
ViscBndry::setBndryConds (const BCRec& phys_bc,
                          int          ratio)
{
    const IntVect& ratio_vect = ratio * IntVect::TheUnitVector();
    setBndryConds(phys_bc, ratio_vect);
}

void
ViscBndry::setBndryConds (const BCRec&   bc,
                          const IntVect& ratio,
			  int comp)
{
    //
    //  NOTE: ALL BCLOC VALUES ARE NOW DEFINED AS A LENGTH IN PHYSICAL
    //        DIMENSIONS *RELATIVE* TO THE FACE, NOT IN ABSOLUTE PHYSICAL SPACE
    //
    const Real* dx    = geom.CellSize();
    const Box& domain = geom.Domain();
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
                const int p_bc  = (face.isLow() ? bc.lo(dir) : bc.hi(dir));

                if (p_bc == EXT_DIR)
                {
                    bctag[face][comp] = LO_DIRICHLET;
                    bloc[face]        = 0;
                }
                else if (p_bc == FOEXTRAP      ||
                         p_bc == HOEXTRAP      || 
                         p_bc == REFLECT_EVEN)
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
                // Internal bndry.
                //
                const Real delta = dx[dir]*ratio[dir];

                bctag[face][comp] = LO_DIRICHLET;
                bloc[face]        = 0.5*delta;
            }
        }
    }
}

void
ViscBndry::setHomogValues (const BCRec&   bc,
                           const IntVect& ratio)
{
    setBndryConds(bc, ratio);

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
