//
// $Id: MacBndry.cpp,v 1.3 1997-10-01 01:03:08 car Exp $
//

#include <LO_BCTYPES.H>
#include <MacBndry.H>

void MacBndry::setBndryConds(const BCRec& phys_bc,
			     const Geometry& geom, IntVect& ratio)
{

//  NOTE: ALL BCLOC VALUES ARE NOW DEFINED AS A LENGTH IN PHYSICAL DIMENSIONS
//        *RELATIVE* TO THE FACE, NOT IN ABSOLUTE PHYSICAL SPACE

    const BoxArray& grids = boxes();
    int ngrds = grids.length();
    const REAL* dx = geom.CellSize();
    const BOX& domain = geom.Domain();
    const RealBox& prob_domain = geom.ProbDomain();

    for (OrientationIter fi; fi; ++fi) {
	Orientation face(fi());
	Array<REAL> &bloc = bcloc[face];
	Array<BoundCond> &bctag = bcond[face];

	int dir = face.coordDir();
	REAL delta = dx[dir]*ratio[dir];
	int p_bc = (face.isLow() ? phys_bc.lo(dir) : phys_bc.hi(dir));

	for (int i = 0; i < ngrds; i++) {
	    const BOX& grd = grids[i];
	    int faceindx = grd[face] + (face.isLow() ? 0 : 1);

	    if (domain[face] == grd[face] && !geom.isPeriodic(dir)) {
		  // All physical bc values are located on face
		if (p_bc == Outflow) {
		    bctag[i] = LO_DIRICHLET;
		    bloc[i] = 0.;
		} else {
		    bctag[i] = LO_NEUMANN;
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

