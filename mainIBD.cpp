#include <iostream.h>
#include <fstream.h>

#include <ParmParse.H>
#include <Boolean.H>

#include <LO_BCTYPES.H>
#include <TestIBData.H>
#include <CGSolver.H>
#include <MultiGrid.H>
#include <Laplacian.H>
#include <ABecLaplacian.H>
#include <ParallelDescriptor.H>
#include <WriteMultiFab.H>

#include <mainIBD_F.H>
#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
REAL* fabdat = (fab).dataPtr();



BoxList
readBoxList(const aString &file, BOX& domain );



int
main(int argc, char **argv)
{
    ParmParse pp(argc-2,argv+2,NULL,argv[1]);
    
      // Obtain prob domain and box-list, set H per phys domain [0:1]Xn
    BOX domain;
#if (BL_SPACEDIM == 2)
    aString boxfile("grids/gr.2_19boxes") ; pp.query("boxes", boxfile);
#endif
#if (BL_SPACEDIM == 3)
    aString boxfile("grids/gr.3_small_a") ; pp.query("boxes", boxfile);
#endif
    BoxArray bs(readBoxList(boxfile,domain));
    int Ngrids=bs.length();
    int Ncomp=1;
    int ratio=2; pp.query("ratio", ratio);

    // get ParmParse data out for Geometry
    Geometry::Setup();

      // Create the geometries
    Geometry geom(domain);
    Geometry geom_crse( Box(domain).coarsen(ratio) );
    const REAL* h = geom.CellSize();

      // Create the boundary object
    TestIBData vbd(bs,Ncomp,geom);

      // Create the BCRec's interpreted by TestIBData objects
    BCRec pbc(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
	      D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));

      // Create "foreground fine data" (internal IC's + BC's)
    int Nghost=1;
    MultiFab fine(bs,Ncomp,Nghost,Fab_allocate);
    //for (int i=0; i<Ngrids; i++) {
    for(MultiFabIterator finemfi(fine); finemfi.isValid(); ++finemfi) {
	DEF_LIMITS(finemfi(),fdat,flo,fhi);
	FORT_FILLFINE(fdat,ARLIM(flo),ARLIM(fhi),h,&Ncomp);
    }

      // Create "background coarse data"
    BOX crse_bx = geom_crse.Domain();
    crse_bx.grow(1);
    const REAL* h_crse = geom_crse.CellSize();
    FARRAYBOX crse_fab(crse_bx,Ncomp);
    DEF_LIMITS(crse_fab,cdat,clo,chi);
    FORT_FILLCRSE(cdat,ARLIM(clo),ARLIM(chi),h_crse,&Ncomp);

      // Create coarse boundary register, fill w/data from coarse FAB
    int bndry_InRad=0;
    int bndry_OutRad=1;
    int bndry_Extent=1;
    BoxArray cbs(bs);
    cbs.coarsen(ratio);
    BndryRegister cbr(cbs,bndry_InRad,bndry_OutRad,bndry_Extent,Ncomp);
    for (OrientationIter face; face; ++face) {
	Orientation f = face();
	FabSet& bnd_fs(cbr[f]);
	bnd_fs.copyFrom(crse_fab);
    }
    
      // Interpolate crse data to fine boundary, where applicable
    int cbr_Nstart=0;
    int fine_Nstart=0;
    int bndry_Nstart=0;
    vbd.setBndryValues(cbr,cbr_Nstart,fine,fine_Nstart,
		       bndry_Nstart,Ncomp,ratio,pbc);

      // Check results
      /*
    for (OrientationIter f_i; f_i; ++f_i) {
	Orientation face(f_i());
	const Array<BoundCond>& bcs(vbd.bndryConds(face));
	const Array<REAL>& bls(vbd.bndryLocs(face));
	const FabSet& bvs(vbd.bndryValues(face));
	int Nbcs=bls.length();
	for (int i=0; i<Nbcs; i++) {
	  if(ParallelDescriptor::IOProcessor()) {
	    cout
		<< "face= " << face
		<< " i= " << i
		<< " bl = " << bls[i]
		<< " bc = " << bcs[i]
		<< endl;
	    cout << bvs[i] << endl;
	  }
	}
    }
    */
    MultiFab soln(bs, Ncomp, Nghost, Fab_allocate); soln.setVal(0.0);
    MultiFab  rhs(bs, Ncomp, Nghost, Fab_allocate);  rhs.setVal(0.0);

    REAL alpha=1.0; pp.query("alpha",alpha);
    REAL beta=1.0; pp.query("beta",beta);
    MultiFab  acoefs;
    acoefs.define(bs, Ncomp, Nghost, Fab_allocate);
    REAL a=0.0; pp.query("a",  a);
    acoefs.setVal(a);

    MultiFab bcoefs[BL_SPACEDIM];
    Tuple<REAL, BL_SPACEDIM> b;
    b[0]=1.0; pp.query("b0", b[0]);
    b[1]=1.0; pp.query("b1", b[1]);
#if (BL_SPACEDIM > 2)
    b[2]=1.0; pp.query("b2", b[2]);
#endif
    for (int n=0; n<BL_SPACEDIM; ++n) {
	BoxArray bsC(bs);
	bcoefs[n].define(bsC.surroundingNodes(n),Ncomp,
			 Nghost,Fab_allocate);
	bcoefs[n].setVal(b[n]);
    } // -->> over dimension

    Laplacian lp(bs, vbd, h[0]);
    bool use_mg_precond=true;
    REAL rtol=1.e-10;
    REAL atol=1.e-10;
    soln.setVal(0.0);
    CGSolver cg(lp,use_mg_precond);
    cg.solve(soln, rhs, rtol, atol);
      //MultiGrid mg(lp);
      //mg.solve(soln, rhs, rtol, atol);
      /*
    ABecLaplacian lp(bs, vbd, h);
    lp.maxOrder(2);
    lp.setScalars(alpha,beta);
    lp.setCoefficients(acoefs,bcoefs);
    MultiGrid mg(lp);
    REAL tolerance = 1.0e-10; pp.query("tolerance", tolerance);
    mg.solve(soln, rhs, vbd, tolerance);
    */
    ofstream os("pltfile");
    int comp = 0;
    REAL s_min=soln.min(comp);
    REAL s_max=soln.max(comp);
    REAL bgValue = s_min - 0.1*(s_max - s_min);
    FArrayBox::setFormat(FABio::FAB_NATIVE);
    WriteMultiFab(os,soln,h[0],bs,domain,ratio,bgValue);

    return 0;
} // -->> main fnc


BoxList
readBoxList(const aString& file, BOX& domain )
{
    BoxList retval;
    ifstream boxspec(file.c_str());
    if( !boxspec ) {
	BoxLib::Error("readBoxList: unable to open " + *file.c_str());
    }
    boxspec >> domain;
    
    int numbox;
    boxspec >> numbox;

    for(int i=0; i<numbox; i++) {
	BOX tmpbox;
	boxspec >> tmpbox;
	if( ! domain.contains(tmpbox)) {
	    cerr << "readBoxList: bogus box " << tmpbox << endl;
	    exit(1);
        }
	retval.append(tmpbox);
    }
    boxspec.close();
    return retval;
}

