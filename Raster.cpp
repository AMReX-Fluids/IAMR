#if (BL_SPACEDIM == 2) || (BL_SPACEDIM == 3)

#include <Raster.H>
#include <Utility.H>
#include <ParmParse.H>
#include <Amr.H>
#include <AmrLevel.H>

typedef unsigned char  byte;
static  const int min_color = 2;
static  const int max_color = 254;
static  const int color_rng = max_color - min_color + 1;

int Raster::verbose = false;

#if (BL_SPACEDIM == 2)
enum outputType {HDF_SDS = 0, HDF_8, NO_FORM};
static const char* typ_names[NO_FORM] = {
    "hdf_sds",
    "hdf_8"
};
#endif
#if (BL_SPACEDIM == 3)
enum outputType {HDF_SDS = 0, AVS_FIELD, AVS_VOLUME, NO_FORM};
static const char* typ_names[NO_FORM] = {
    "hdf_sds"
};
#endif

// -------------------------------------------------------------
// External linkage to HDF modules
extern "C" {
    extern int DFSDsetmaxmin( float, float );
    extern int DFSDsetdims( const int rank, const int* dimsizes );
    extern int DFSDputdata( char* filename, const int rank, 
			    const int* dimsizes, float* data );
    extern int DFSDsetdatastrs( char*,char*,char*,char* );
    extern int DFSDsetdimstrs( int dim, char*label,char*unit,char*format);
    extern int DFSDadddata( char* filename, const int rank, 
			    const int* dimsizes, float* data );
    extern int DFSDgetdims( char* filename, int& rank, int* dimsizes,
			    int maxrank);
    extern int DFSDgetdata( char* filename, int rank, int* dimsizes,
			    float* data );
    extern int DFR8putimage( char* filename, char*image, int width,
			     int height, int compress );
    extern int DFR8addimage( char* filename, char*image, int width,
			     int height, int compress );
    extern int DFSDrestart();
}
#define DFTAG_RLE   (11) /* run length encoding */
#define DFTAG_IMC   (12) /* IMCOMP compression */
#define DFTAG_IMCOMP   (12) /* IMCOMP compression */


// -------------------------------------------------------------
RastFrame::RastFrame(aString& file_nm, aString& var_nm, int otyp, int freq,
                     int lev, BOX& graph_box, 
		     Amr& /*amrsys*/) {
    var_name = var_nm;
    file_name = file_nm;
    interval = freq;
    o_type   = otyp;
    level    = lev;
    box      = graph_box;
}

// -------------------------------------------------------------
RastFrame::~RastFrame() 
{}

//   ############################################################
//   ##### Raster members
//   ############################################################

// -------------------------------------------------------------
Raster::Raster(Amr &sys) 
{
    amrsys = &sys;

      // parse input file
    ParmParse pp("raster");
    dump_locations = pp.contains("grid_loc");
    verbose = pp.contains("verbose");
    int graph_lev = -1;
    pp.query("level",graph_lev);
    Array<int> g_box(2*BL_SPACEDIM);
    BOX graph_box;
    if (pp.queryarr("region",g_box,0,2*BL_SPACEDIM)) {
	INTVECT lo,hi;
        int k;
	for (k = 0; k < BL_SPACEDIM; k++) {
	    lo.setVal(k,g_box[2*k]);
	    hi.setVal(k,g_box[1+2*k]);
	}
	graph_box = BOX(lo,hi);
    }

    time_plot = 0;
    plot_cnt = 0;
    plot_dt = -100.0;
    if (pp.query("plot_dt",plot_dt)) {
	time_plot = 1;
    }
	
    int n = pp.countname("plot");
    int k;
    for (k = 0; k < n; k++) {
	int freq = -1;
	pp.getkth("plot",k,freq,2);
	if (freq > 0) {
	    aString fname;
	    aString vname;
	    aString typ;
	    pp.getkth("plot",k,fname,0);
	    pp.getkth("plot",k,vname,1);
	    pp.getkth("plot",k,typ,3);
	    int found = false;
            int otyp;
	    for (otyp = 0; otyp < NO_FORM; otyp++) {
		if (typ == typ_names[otyp]) {
		    found = true;
		    break;
		}
	    }
	    if (!found) {
		cerr << "RASTER:: invalid raster type = " << typ
		     << " ... will skip" << endl;
	    } else {
		RastFrame* rf = new RastFrame(fname,vname,otyp,freq,
					      graph_lev,graph_box,*amrsys);
		frames.append(rf);
	    }
	}
    }
}

// -------------------------------------------------------------
Raster::~Raster()
{
    clear();
}

// -------------------------------------------------------------
void Raster::clear()
{
    ListIterator<RastFrame*> li(frames);
    for (; li; ++li) {
	RastFrame *r = li();
	delete r;
	frames[li] = 0;
    }
    frames.clear();
}


// -------------------------------------------------------------
ostream& operator << (ostream &os, Raster &rast)
{
    ListIterator<RastFrame*> li(rast.frames);
    while (li) {
	RastFrame& r = *(li());
	os << "rastplot  " << typ_names[r.o_type]
           << "  " << r.var_name << " level = ";
        if (r.level < 0) {
	    os << "default";
        } else {
	    os << r.level;
        }
        os << " box = ";
        if (r.box.ok()) {
	    os << r.box;
        } else {
	    os << "default";
        }
        os << endl;
	++li;
    }
    return os;
}

#if (BL_SPACEDIM == 2)
// -------------------------------------------------------------
static void draw_hdf_8(FARRAYBOX& data, REAL u_min, REAL u_max, char* fname)
{
    BOX workbox(data.box());
    REAL u_rng = u_max - u_min;
    if (u_rng < 0.000001*(fabs(u_max)+fabs(u_min)) ) u_rng = 1.0;
    int  size = workbox.numPts();
    byte *rast = new byte[size];
    REAL *dat = data.dataPtr();
    const int* len = workbox.length();
    int nx = len[0];
    int ny = len[1];
    int j;
    for (j = 0; j < ny; j++) {
	int jflip = ny - j - 1;
        int i;
	for (i = 0; i < nx; i++) {
	    REAL v = dat[i + nx*jflip];
	    REAL v2 = (v-u_min)/u_rng;
	    int  icolor = min_color + int(v2*color_rng);
	    rast[i+nx*j] = byte(icolor);
	}
    }

      // now write out raster
    int err = DFR8putimage(fname, (char *)rast,nx,ny,DFTAG_RLE);
    if( err == -1 ){
        cerr << "Raster::dumpRast:HDF error" << endl;
        abort();
    }
    delete rast;
}
#endif



// -------------------------------------------------------------
static void draw_hdf_sds(FARRAYBOX &data, REAL /*u_min*/, REAL /*u_max*/,
                         char* fname, char* var_name)
{
    BOX workbox(data.box());
    int  size = workbox.numPts();

      // since REAL may not be 32 bit, need to make work array
    float *scidata = new float[size];
    if (scidata == NULL) {
	cerr << "draw_hdf_sds:: out of memory" << endl;
	abort();
    }
    REAL *dat = data.dataPtr();
    const int* len = workbox.length();
    int revlen[BL_SPACEDIM];
      // copy into work space
    int nx = len[0];
    int ny = len[1];
#  if (BL_SPACEDIM==2)
    revlen[0] = len[1];
    revlen[1] = len[0];
    int j;
    for (j = 0; j < ny; j++) {
	int jflip = ny - j - 1;
        int i;
	for (i = 0; i < nx; i++) {
	    scidata[i+nx*j] = dat[i+nx*jflip];
	}
    }
#  endif
#  if (BL_SPACEDIM==3)
    revlen[0] = len[0];
    revlen[1] = len[1];
    revlen[2] = len[2];
    int nz = len[2];
    int k;
    for(k = 0; k < nz; k++ ){
        int j;
	for (j = 0; j < ny; j++) {
	    int jflip = ny - j - 1;
            int i;
	    for (i = 0; i < nx; i++) {
		scidata[k+nz*(j+i*ny)] = dat[i+nx*(j+k*ny)];
	    }
	}
    }
#   endif

      // now write out scientific data set

      // prepare the data set
    if( DFSDsetdims( BL_SPACEDIM, revlen ) ){
     	cerr << "Raster::dumpData unable to set dims" << endl;
     	abort();
    }
    if( DFSDsetdimstrs(1,"x","nat","%f") ) {
     	cerr << "Raster::dumpData unable set dim str 1" << endl;
     	abort();
    }
    if( DFSDsetdimstrs(2,"y","nat","%f") ) {
     	cerr << "Raster::dumpData unable set dim str 2" << endl;
     	abort();
    }
#   if (BL_SPACEDIM==3)
    if( DFSDsetdimstrs(3,"z","nat","%f") ) {
	cerr << "Raster::dumpData unable  set dim str 3" << endl;
	abort();
    }
#   endif
#ifdef OLD_HDF
    if( DFSDsetmaxmin( u_max, u_min) ) {
    	cerr << "Raster::dumpData unable to setmaxmin " << endl;
     	abort();
    }
#endif    
    if( DFSDsetdatastrs(var_name,"nat","%e","") ){
     	cerr << "Raster::dumpData unable to set data str" << endl;
     	abort();
    }
    if( DFSDputdata( fname, BL_SPACEDIM, revlen,scidata ) ){
     	cerr << "Raster::dumpData unable to adddata" << endl;
     	abort();
    }

    delete scidata;
}


// -------------------------------------------------------------
int Raster::draw(REAL time, int nstep, int force_draw)
{
      // loop through frames, creating and writing images
    int num_frames = 0;
    int doit = false;
    if (force_draw) {
	doit = true;
    } else if (time_plot) {
	if (time < plot_cnt*plot_dt) return 0;
	doit = (time >= plot_cnt*plot_dt);
	if (doit) {
	    while (time >= plot_cnt*plot_dt) plot_cnt++;
	}
    }

    for(ListIterator<RastFrame*> li(frames); li; ++li) {
	RastFrame& r = *(li());
	int intval = r.interval;
	if (!time_plot) {
	    doit = doit || ((intval > 0) && (nstep%intval == 0));
	}
	if (doit && (intval > 0)) {
	    num_frames++;
	    outputType out = (outputType) r.o_type;

	    int rlev = r.level;
	    int amrlev = amrsys->finestLevel();
	    if (rlev < 0) rlev = amrlev;
	    rlev = Min(rlev,amrlev);

	    BOX workbox(r.box);
	    if (!workbox.ok()) {
		workbox = amrsys->Geom(0).Domain();
	    }
	    if (rlev > 0) {
		IntVect lrat = IntVect::TheUnitVector();
                int i;
		for (i = 0; i < rlev; i++) lrat *= amrsys->refRatio(i);
		workbox.refine(lrat);
	    }
	    workbox &= (amrsys->Geom(rlev).Domain());
	    char* v_name = (char*) (r.var_name).c_str();
	    char* f_name = (char*) (r.file_name).c_str();

	    FARRAYBOX* data = amrsys->derive(workbox,v_name,time,rlev);

	    char fname[80];
	    switch (out) {
	    case HDF_SDS:
		sprintf(fname,"%s%04d.hdf",f_name,nstep);
		break;
#if (BL_SPACEDIM == 2)
	    case HDF_8:
		sprintf(fname,"%s%04d.hdf",f_name,nstep);
		break;
#endif
	    case NO_FORM:
		cerr << "RASTER::draw(): OOPS no file format for this raster"
		     << endl;
		abort();
	    }
	    REAL u_min = data->min(0);
	    REAL u_max = data->max(0);
	    if (verbose) {
		cout << "RASTER " << v_name <<
	            "  min = " << u_min << " max = " << u_max << endl;
	    }

	    switch (out) {
	    case HDF_SDS:
		draw_hdf_sds(*data,u_min,u_max,fname,v_name);
		break;
	       
#if (BL_SPACEDIM == 2)
	    case HDF_8:
		draw_hdf_8(*data,u_min,u_max,fname);
		break;
#endif
	    }

	    delete data;
	    delete v_name;
	    delete f_name;
	}
    }

    if (num_frames > 0 && dump_locations) gridLocations(nstep);

    return num_frames;
}

// ###################################################################
// ##### GRID LOCATIONS
// ###################################################################
void Raster::gridLocations(int nstep) {
    char c_name[80];
    sprintf(c_name,"Grids%04d.ps",nstep);
    ofstream psfile(c_name,ios::out);
    psfile << "% Grids from at step " << nstep << '\n';
    psfile << "1 setlinecap\n";
    psfile << "1 setlinejoin\n\n";
    int finest_level = amrsys->finestLevel();
    int lev;
    for (lev = 0; lev <= finest_level; lev++) {
	AmrLevel& amr_lev = amrsys->getLevel(lev);
	int ngrids = amr_lev.numGrids();
	const Array<REALBOX>& grds = amr_lev.gridLocations();
	psfile << "%  " << ngrids << " Grids at level " << lev << '\n';
	switch (lev) {
	case 0 : psfile << "0 0 1 setrgbcolor\n"; break;
	case 1 : psfile << "0 1 0 setrgbcolor\n"; break;
	case 2 : psfile << "1 0 0 setrgbcolor\n"; break;
	default : psfile << "0 0 0 setrgbcolor\n"; break;
	}
        int i;
	for (i = 0; i < ngrids; i++) {
	    REAL x1 = grds[i].lo(0);
	    REAL y1 = grds[i].lo(1);
	    REAL x2 = grds[i].hi(0);
	    REAL y2 = grds[i].hi(1);
	    psfile << "% Grid " << i << '\n';
	    psfile << x1 << " inch " << y1 << " inch moveto\n";
	    psfile << x2 << " inch " << y1 << " inch lineto\n";
	    psfile << x2 << " inch " << y2 << " inch lineto\n";
	    psfile << x1 << " inch " << y2 << " inch lineto\n";
	    psfile << x1 << " inch " << y1 << " inch lineto\n";
	}
	psfile << "stroke" << endl;
    }
    psfile.close();
}

#endif
