//
// $Id: TestIBData.cpp,v 1.1 1997-07-08 23:08:02 vince Exp $
//
#include <LO_BCTYPES.H>
#include <TestIBData.H>

void TestIBData::setBndryConds(const BCRec& bc,
			       const Geometry& geom, IntVect& ratio)
{
    const BoxArray& grids = boxes();
    int ngrds = grids.length();
    const REAL* dx = geom.CellSize();
    const BOX& domain = geom.Domain();
    const REALBOX& prob_domain = geom.ProbDomain();

    for (OrientationIter fi; fi; ++fi) {
	Orientation face(fi());
	Array<REAL> &bloc = bcloc[face];
	Array<BoundCond> &bctag = bcond[face];

	int dir = face.coordDir();
	REAL delta = dx[dir]*ratio[dir];
	int problo = prob_domain.lo(dir);
	int domlo = domain.smallEnd(dir);
	int p_bc = (face.isLow() ? bc.lo(dir) : bc.hi(dir));

	for (int i = 0; i < ngrds; i++) {
	    const BOX& grd = grids[i];

	    if (domain[face] == grd[face]) {
		  // All physical bc values are located on face
		if (p_bc == EXT_DIR || p_bc == REFLECT_ODD) {
		    bctag[i] = LO_DIRICHLET;
		    bloc[i] = 0.0; // on face, distance to face = 0
		} else if (p_bc == EXTRAP || p_bc == HOEXTRAP || 
                           p_bc == REFLECT_EVEN) {
		    bctag[i] = LO_NEUMANN;
		    bloc[i] = 0.0; // on face, distance to face = 0
		}
	    } else {
		  // internal bndry
		bctag[i] = LO_DIRICHLET;
		bloc[i] = 0.5*delta; // internal, distance is half of crse
	    }
	}
    }
}
