#if (BL_SPACEDIM==2) && defined (USE_TENSOR)

//
// $Id: ViscBndry2D.cpp,v 1.8 1998-06-09 21:42:51 lijewski Exp $
//

#include <LO_BCTYPES.H>
#include <ViscBndry2D.H>

void ViscBndry2D::setBndryConds(const Array<BCRec>& bcarray, int ratio)
{
    int ncomp = bcond[0][0].length();
    assert(ncomp==2*2); // u and v, plus derivs of same
    assert(bcarray.length()==2*2);

    const BoxArray& grids = boxes();
    int ngrds = grids.length();
    const REAL* dx = geom.CellSize();
    const BOX& domain = geom.Domain();
    const RealBox& prob_domain = geom.ProbDomain();

    for (OrientationIter fi; fi; ++fi) {
        Orientation face(fi());
        Array<REAL> &bloc = bcloc[face];

        int dir = face.coordDir();
#if 0
        REAL delta = dx[dir]*ratio[dir];
#else
        REAL delta = dx[dir]*ratio;
#endif
        for( int icomp=0; icomp<ncomp; icomp++){
          int p_bc = (face.isLow() ? 
                      bcarray[icomp].lo(dir):bcarray[icomp].hi(dir));

          for (int i = 0; i < ngrds; i++) {
            const BOX& grd = grids[i];
            // bctag is bc type (with array info on orientation,grid,comp)gone
            BoundCond  &bctag = bcond[face][i][icomp];

            if (domain[face] == grd[face] && !geom.isPeriodic(dir)) {
              // All physical bc values are located on face
              if (p_bc == EXT_DIR ) {
                bctag = LO_DIRICHLET;
                bloc[i] = 0.0; // on face, distance to face = 0
              } else if (p_bc == EXTRAP || p_bc == HOEXTRAP || 
                         p_bc == REFLECT_EVEN) {
                bctag = LO_NEUMANN;
                bloc[i] = 0.0; // on face, distance to face = 0
              } else if( p_bc == REFLECT_ODD ){
                bctag = LO_REFLECT_ODD;
                bloc[i] = 0.0; // on face, distance to face = 0
              }
            } else {
              // internal bndry
              bctag = LO_DIRICHLET;
              bloc[i] = 0.5*delta; // internal, distance is half of crse
            }
          }
        }
    }
}


// *************************************************************************

void
ViscBndry2D::setHomogValues(const Array<BCRec>& bc, int ratio)
{

    setBndryConds(bc, ratio);

    for (OrientationIter fi; fi; ++fi)
    {
        Orientation face(fi());
        for (FabSetIterator fsi(bndry[face]); fsi.isValid(false); ++fsi)
        {
            fsi().setVal(0.);
        }
    }
}
#endif
