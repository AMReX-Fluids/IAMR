//
// $Id: ViscBndry.cpp,v 1.4 1997-10-01 01:03:17 car Exp $
//

#include <LO_BCTYPES.H>
#include <ViscBndry.H>

void ViscBndry::setBndryConds(const BCRec& bc,
			      const Geometry& geom, IntVect& ratio)
{

//  NOTE: ALL BCLOC VALUES ARE NOW DEFINED AS A LENGTH IN PHYSICAL DIMENSIONS
//        *RELATIVE* TO THE FACE, NOT IN ABSOLUTE PHYSICAL SPACE

    const BoxArray& grids = boxes();
    int ngrds = grids.length();
    const Real* dx = geom.CellSize();
    const BOX& domain = geom.Domain();
    const RealBox& prob_domain = geom.ProbDomain();

    for (OrientationIter fi; fi; ++fi) {
	Orientation face(fi());
	Array<Real> &bloc = bcloc[face];
	Array<BoundCond> &bctag = bcond[face];

	int dir = face.coordDir();
	Real delta = dx[dir]*ratio[dir];
	int p_bc = (face.isLow() ? bc.lo(dir) : bc.hi(dir));

	for (int i = 0; i < ngrds; i++) {
	    const BOX& grd = grids[i];

	    if (domain[face] == grd[face] && !geom.isPeriodic(dir)) {
		  // All physical bc values are located on face
		if (p_bc == EXT_DIR) {
		    bctag[i] = LO_DIRICHLET;
		    bloc[i] = 0.;
		} else if (p_bc == EXTRAP || p_bc == HOEXTRAP || 
                           p_bc == REFLECT_EVEN) {
		    bctag[i] = LO_NEUMANN;
		    bloc[i] = 0.;
		} else if (p_bc == REFLECT_ODD) {
		    bctag[i] = LO_REFLECT_ODD;
		    bloc[i] = 0.;
		}
	    } else {
		  // internal bndry
		bctag[i] = LO_DIRICHLET;
		bloc[i] = 0.5*delta;
	    }
	}
    }
}

// *************************************************************************

void
ViscBndry::setHomogValues(const BCRec& bc, IntVect& ratio)
{

    setBndryConds(bc, geom, ratio);

/*
    int ngrd = grids.length();
    for (int grd = 0; grd < ngrd; grd++) {
        const BOX& bx = grids[grd];
        for (OrientationIter fi; fi; ++fi) {
            Orientation face(fi());
            FArrayBox& bnd_fab = bndry[face][grd];
            bnd_fab.setVal(0.);
        }
    }
*/
    for(OrientationIter fi; fi; ++fi) {
      Orientation face(fi());
      for(FabSetIterator fsi(bndry[face]); fsi.isValid(); ++fsi) {
        fsi().setVal(0.);
      }
    }
}
