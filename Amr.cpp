
#include <Array.H>
#include <CoordSys.H>
#include <ParmParse.H>
#include <BoxDomain.H>
#include <Cluster.H>
#include <RunStats.H>
#include <LevelBld.H>
#include <AmrLevel.H>
#include <TagBox.H>

#include <PROB_AMR_F.H>
#include <Amr.H>
#include <Amr_auxil.H>
#include <ParallelDescriptor.H>
#include <Utility.H>

#ifdef BL_USE_NEW_HFILES
using std::ifstream;
using std::ios;
#endif

const char NL = '\n';

// ###################################################################
// #####   INLINE-able functions
// ###################################################################
// ###################################################################
int Amr::maxGridSize() const		{ return max_grid_size; }

int Amr::maxLevel() const		{ return max_level;}

int Amr::finestLevel() const		{ return finest_level; }

IntVect Amr::refRatio(int level) const { return ref_ratio[level]; }

// This function should be deleted when an IntVect member function to
// return the maximum component is added to BoxLib.
int Amr::MaxRefRatio( int level ) const {
                                  int maxval=0;
                                  for( int n=0; n<BL_SPACEDIM; n++ ) 
				  maxval = Max(maxval,ref_ratio[level][n]);
				  return maxval;
                                  }

int Amr::nCycle(int level) const 	{ return n_cycle[level]; }

const Array<IntVect>& Amr::refRatio() const 	{ return ref_ratio; }

REAL Amr::dtLevel(int level) const  	{ return dt_level[level];  }

const Array<REAL>& Amr::dtLevel() const         { return dt_level; }

const Geometry&
Amr::Geom(int level) const
{
    return geom[level];
}

AmrLevel&
Amr::getLevel(int lev)
{
    return amr_level[lev];
}

PArray<AmrLevel>& Amr::getAmrLevels() { return amr_level; }

int Amr::levelSteps(int i) const	{ return level_steps[i]; }

REAL Amr::cumTime() const		{ return cumtime; }

int Amr::regridInt(int lev) const	{ return regrid_int[lev]; }

int Amr::nErrorBuf(int lev) const	{ return n_error_buf[lev]; }

REAL Amr::gridEff() const	{ return grid_eff; }

int Amr::subCycle() const     { return sub_cycle; }


// constructors and destructors

// ###################################################################
// ##### DEFAULT AMR CONSTRUCTOR
// ###################################################################
Amr::Amr()
    :amr_level(PArrayManage),
     grids_file(),
     check_file_root(),
     restart_file(),
     plot_file_root()
{
  // setup Geometry from ParmParse file.  May be needed for variableSetup or
  // even getLevelBld.
  Geometry::Setup();

      // determine physics class
#ifndef TESTINGMG
    levelbld = getLevelBld();
#endif

      // global function that define state variables
    levelbld->variableSetUp();

      // set default values
    max_level     = -1;
    record_run_info   = false;
    record_grid_info  = false;

    grid_eff	      = 0.7;
    blocking_factor   = 1;
    last_checkpoint   = 0;
    last_plotfile     = 0;
    plot_int          = -1;
    n_proper          = 1;
#if (BL_SPACEDIM == 2)
    max_grid_size     = 120;
#else
    max_grid_size     = 32;
#endif

    int i;
    for (i = 0; i < BL_SPACEDIM; i++) isPeriodic[i] = false;

        // construct parameter parsing object that will link into
	// global table of input and command line options
    ParmParse pp("amr");

      // check for command line flags
    trace = pp.contains("trace");
    debug = pp.contains("debug");
    silent = pp.contains("silent");
    verbose = pp.contains("verbose");
    sub_cycle = true;
    if (pp.contains("nosub")) sub_cycle = false;

    pp.query("regrid_file",grids_file);
    if (pp.contains("run_log")) {
	aString log_file_name;
	pp.get("run_log",log_file_name);
	setRecordRunInfo(log_file_name);
    }
    if (pp.contains("grid_log")) {
	aString grid_file_name;
	pp.get("grid_log",grid_file_name);
	setRecordGridInfo(grid_file_name);
    }

      // restart or run from scratch?
    pp.query("restart",restart_file);

      // read max_level and alloc memory for container objects
    pp.get("max_level",max_level);
    int nlev     = max_level+1;
    geom.resize(nlev);
    dt_level.resize(nlev);
    level_steps.resize(nlev);
    level_count.resize(nlev);
    regrid_int.resize(nlev);
    n_cycle.resize(nlev);
    dt_min.resize(nlev);
    n_error_buf.resize(nlev);
    amr_level.resize(nlev);
      // set bogus values
    for (i = 0; i < nlev; i++) {
	dt_level[i] = 0.0;
	level_steps[i] = 0;
	level_count[i] = 0;
	regrid_int[i] = 0;
	n_cycle[i] = 0;
	dt_min[i] = 0.0;
	n_error_buf[i] = 1;
    }
    ref_ratio.resize(max_level);
    for (i = 0; i < max_level; i++) ref_ratio[i] = IntVect::TheZeroVector();

      // read other amr specific values

    check_file_root = "chk";
    pp.query("check_file",check_file_root);

    check_int = -1;
    int got_check_int = pp.query("check_int",check_int);

    check_per = -1.0;
    int got_check_per = pp.query("check_per",check_per);

    if (got_check_int == 1 && got_check_per == 1) {
      BoxLib::Error("MUST ONLY SPECIFY amr.check_int OR amr.check_per");
    } else if (got_check_per == 1) {
      BoxLib::Warning("Specifying amr.check_per will change the time step.");
    }

    plot_file_root = "plt";
    pp.query("plot_file",plot_file_root);

    plot_int = -1;
    int got_plot_int = pp.query("plot_int",plot_int);

    plot_per = -1.0;
    int got_plot_per = pp.query("plot_per",plot_per);

    if (got_plot_int == 1 && got_plot_per == 1) {
      BoxLib::Error("MUST ONLY SPECIFY amr.plot_int OR amr.plot_per");
    } else if (got_plot_per == 1) {
      BoxLib::Warning("Specifying amr.plot_per will change the time step.");
    }

    pp.query("max_grid_size",max_grid_size);
    pp.query("n_proper",n_proper);
    pp.get("blocking_factor",blocking_factor);
    pp.get("grid_eff",grid_eff);

    pp.getarr("n_error_buf",n_error_buf,0,max_level);

    // read in the refinement ratio IntVects as integer BL_SPACEDIM-tuples
    if (max_level > 0) {

      int nratios_vect = max_level*BL_SPACEDIM;
      Array<int> ratios_vect(nratios_vect);

      int got_vect = pp.queryarr("ref_ratio_vect",ratios_vect,0,nratios_vect);

      Array<int> ratios(max_level);

      int got_int = pp.queryarr("ref_ratio"     ,ratios,0,max_level);
   
      if (got_int == 1 && got_vect == 1) {
        BoxLib::Error("Only input *either* ref_ratio or ref_ratio_vect");
      } else if (got_vect == 1) {
        int k = 0;
        for( i = 0; i < max_level; i++ ) {
    	  for( int n=0; n<BL_SPACEDIM; n++,k++ ) {
	    ref_ratio[i][n] = ratios_vect[k];
	  }
        }

      } else if (got_int == 1) {
        for( i = 0; i < max_level; i++ ) {
    	  for( int n=0; n<BL_SPACEDIM; n++) {
	    ref_ratio[i][n] = ratios[i];
	  }
        }

     } else {
       BoxLib::Error("Must input *either* ref_ratio or ref_ratio_vect");
     }
    }

      // read computational domain and set geometry
    Array<int> n_cell(BL_SPACEDIM);
    pp.getarr("n_cell",n_cell,0,BL_SPACEDIM);
    assert(n_cell.length() == BL_SPACEDIM);
    INTVECT lo(IntVect::TheZeroVector()), hi(n_cell);
    hi -= IntVect::TheUnitVector();
    BOX index_domain(lo,hi);
    for(i = 0; i <= max_level; i++) {
	geom[i].define(index_domain);
	if (i < max_level) index_domain.refine(ref_ratio[i]);
    }

      // Now define offset for CoordSys
    REAL offset[BL_SPACEDIM];
    for (i = 0; i < BL_SPACEDIM; i++) {
	REAL delta = Geometry::ProbLength(i)/(REAL)n_cell[i];
	offset[i] = Geometry::ProbLo(i) + delta*lo[i];
    }
    CoordSys::SetOffset(offset);

        // set regrid interval
    int ri;
    pp.get("regrid_int",ri);
    int k;
    for (k = 0; k <= max_level; k++) regrid_int[k] = ri;
    if (!sub_cycle) {
	  // must adjust regridding trigger
	int factor = 1;
	for (k = max_level-1; k >= 0; k--) {
	    factor *= MaxRefRatio(k);
	    regrid_int[k] = ri*factor;
	}
    }
}

// ###################################################################
// ##### DESTRUCTOR
// ###################################################################
Amr::~Amr(){
      // check to see if we should write a final checkpoint file
    if (level_steps[0] > last_checkpoint) checkPoint();
    if (level_steps[0] > last_plotfile) writePlotFile();

    levelbld->variableCleanUp();
}

// ###################################################################
// ##### SET_RECORD_GRID_INFO
// ###################################################################
void Amr::setRecordGridInfo( const aString& filename ){
    record_grid_info= true;
    gridlog.open(filename.c_str(),ios::out);
    if (!gridlog) {
	cerr << "can't open file: " << filename << '\n';
	BoxLib::Abort();
    }
}

// ###################################################################
// ##### SET_RECORD_RUN_INFO
// ###################################################################
void Amr::setRecordRunInfo( const aString& filename ){
    record_run_info= true;
    runlog.open(filename.c_str(),ios::out);
    if (!runlog) {
	cerr << "can't open file: " << filename << '\n';
	BoxLib::Abort();
    }
}

// ###################################################################
// ##### SET_DT_LEVEL
// ###################################################################
void Amr::setDtLevel(const Array<REAL>& dt_lev) {
    int i;
    for (i = 0; i <= finest_level; i++) {
	dt_level[i] = dt_lev[i];
    }
}

// ###################################################################
// ##### SET_N_CYCLE
// ###################################################################
void Amr::setNCycle(const Array<int>& ns) {
    int i;
    for (i = 0; i <= finest_level; i++) {
	n_cycle[i] = ns[i];
    }
}

// ###################################################################
// ##### CELL_COUNT (total)
// ###################################################################
long Amr::cellCount() {
    long cnt = 0;
    int i;
    for (i = 0; i <= finest_level; i++) {
	cnt += amr_level[i].countCells();
    }
    return cnt;
}

// ###################################################################
// ##### CELL_COUNT (on level)
// ###################################################################
long Amr::cellCount(int lev) {
    return amr_level[lev].countCells();
}

// ###################################################################
// ##### NUM_GRIDS (on level)
// ###################################################################
int Amr::numGrids(int lev) {
    return amr_level[lev].numGrids();
}

// ###################################################################
// ##### NUM_GRIDS (total)
// ###################################################################
int Amr::numGrids() {
    int cnt = 0;
    int i;
    for (i = 0; i <= finest_level; i++) {
	cnt += amr_level[i].numGrids();
    }
    return cnt;
}

// ###################################################################
// ##### BLOCKING_FACTOR
// ###################################################################
int Amr::blockingFactor() const {
    return blocking_factor;
}

const BoxArray& Amr::boxArray(int lev) const
{
    return amr_level[lev].boxArray();
}

int Amr::okToContinue() {
    int ok = true;
    int i;
    for (i = 0; ok && (i <= finest_level); i++) {
	ok = ok && amr_level[i].okToContinue();
    }
    return ok;
}

void Amr::writePlotFile()
{
    writePlotFile(plot_file_root,level_steps[0]);
}

static
aString
Concatenate (const aString& root,
             int            num)
{
    aString result = root;

    char buf[sizeof(int) + 1];

    sprintf(buf, "%04d", num);

    result += buf;

    return result;
}

void
Amr::writePlotFile (const aString& root,
                    int            num)
{
    aString pltfile = Concatenate(root, num);

    if (trace)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            cout << "PLOTFILE: file = " << pltfile << NL;
        }
    }
    if (record_run_info)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            runlog << "PLOTFILE: file = " << pltfile << NL;
        }
    }

    ofstream os;

#ifdef BL_T3E
    // ---------------- new buffer code
    int buffersize = 40960 * 32;
    signed char *iobuff = new signed char[buffersize];
    os.setbuf(iobuff, buffersize);
    // ------------ end new buffer code
#endif

    os.open(pltfile.c_str(), ios::out);

    if (os.fail())
    {
        cerr << "Error in Amr::writePlotFile:  os failed." << NL;
        return;
    }

    for (int k = 0; k <= finest_level; k++)
	amr_level[k].writePlotFile(os);

    os.close();

#ifdef BL_T3E
    // ---------------- new buffer code
    delete [] iobuff;
    // ------------ end new buffer code
#endif
}

void Amr::checkInput() {
    if( max_level < 0 ) {
	BoxLib::Error("checkInput: max_level not set");
    }
      // check that multigrid factor is a power of 2
    int k = blocking_factor;
    while ( k > 0 && (k%2 == 0) ) k = k/2;
    if (k != 1) {
	BoxLib::Error("Amr::checkInputs: multiGrid factor not a power of 2");
    }
      // check level dependent values
    int i;
    for (i = 0; i < max_level; i++){
	if( MaxRefRatio(i) < 2 || MaxRefRatio(i) > 12 ){
	    BoxLib::Error("checkInput bad ref_ratios");
	}
    }

    const BOX& domain = geom[0].Domain();
    if( !domain.ok() ){
	BoxLib::Error("level 0 domain bad or not set");
    }
      // check that domain size has a factor of blocking_factor
    for (i = 0; i < BL_SPACEDIM; i++) {
	int len = domain.length(i);
	if (len%blocking_factor != 0) {
	    BoxLib::Error("checkInputs: domain size not divisible by blocking_factor");
	}
    }
      // check that max_grid_size has a factor of blocking_factor
    for (i = 0; i < max_level; i++) {
      for( int n=0; n<BL_SPACEDIM; n++ ) {
	int lratio = blocking_factor*ref_ratio[i][n];
	if (max_grid_size%lratio != 0) {
	  BoxLib::Error("checkInputs: max_grid_size not divisible by blocking_factor*ref_ratio");
	}
      }
    }
    if (!Geometry::ProbDomain().ok()) {
	BoxLib::Error("checkInput: bad physical problem size");
    }
    if( regrid_int[0] <= 0 ){
	BoxLib::Error("checkinput: regrid_int not defined");
    }
}


// ###################################################################
// ##### INIT
// ###################################################################
void Amr::init(){

      // either read data from restart file or build initial amr structure
    if (restart_file != "") {
	restart(restart_file);

    } else {
	initialInit();
	checkPoint();
	if (plot_int > 0) writePlotFile();
    }
}

// ###################################################################
// ##### INITIAL_INIT
// ###################################################################
// init checks values and internally generates values
void Amr::initialInit() {

      // check inputs
    checkInput();

      // generate internal values from user-supplied values
    finest_level = 0;

      // init problem dependent data
    REAL strt_time = 0;
    int  init = true;

#ifndef TESTINGMG
    FORT_PROBINIT (&init);
#endif
    cumtime = strt_time;

      // define base level grids
    defBaseLevel(strt_time);

    if (max_level > 0) bldFineLevels(strt_time);

      // build any additional data structures after grid generation
    int lev;
    for (lev = 0; lev <= finest_level; lev++) {
	amr_level[lev].post_regrid(0,finest_level);
    }

      // compute dt and set time levels of all grid data
    amr_level[0].computeInitialDt(finest_level,sub_cycle,
				  n_cycle,ref_ratio,dt_level);
    for (lev = 0; lev <= finest_level; lev++) {
	amr_level[lev].setTimeLevel(strt_time,dt_level[lev],dt_level[lev]);
    }

      // perform any special post_initialization operations
    for (lev = 0; lev <= finest_level; lev++) {
	amr_level[lev].post_init();
    }

    for (lev = 0; lev <= finest_level; lev++) {
	level_count[lev] = 0;
	level_steps[lev] = 0;
    }

    if (record_grid_info) {
	gridlog << "INITIAL GRIDS \n";
	printGridInfo(gridlog,0,finest_level);
    }

}

// ###################################################################
// ##### DERIVE (all grids on a level)
// ###################################################################
PArray<FARRAYBOX>* Amr::derive(const aString& name, REAL time, int lev)
{
    return amr_level[lev].derive(name,time);
}


// ###################################################################
// ##### DERIVE (on a given subrange of a level)
// ###################################################################
FARRAYBOX*
Amr::derive(const BOX& b, const aString& name, REAL time, int lev)
{
    return amr_level[lev].derive(b,name,time);
}


// ###################################################################
// ##### RESTART
// ###################################################################
void Amr::restart(const aString& filename)
{
    int i;

    if (trace)
    {
      if(ParallelDescriptor::IOProcessor())
      {
	cout << "restarting calculation from file " << filename << NL;
      }
    }

      // init problem dependent data
    int  init = false;
#ifndef TESTINGMG
    FORT_PROBINIT (&init);
#endif

      // start calculation from given restart file
    if (record_run_info) {
      if(ParallelDescriptor::IOProcessor()) {
	runlog << "RESTART from file = " << filename << NL;
      }
    }

      // open stream
      // for parallel, everyone opens the stream and reads everything
      // but non local grids

      // NEED io buffer here for the t3e

    ifstream is(filename.c_str(),ios::in);

      // read global data
    int   spdim;
    is >> spdim;
    if (spdim != BL_SPACEDIM)
    {
	cerr << "restart: bad spacedim = " << spdim << '\n';
	BoxLib::Abort();
    }
    is >> cumtime;
    int mx_lev;
    is >> mx_lev;
    if (max_level < mx_lev) {
	BoxLib::Error("restart: different max_level");
    }
    is >> finest_level;

    for (i = 0; i <= mx_lev; i++) is >> geom[i];
    for (i = 0; i <  mx_lev; i++) is >> ref_ratio[i];

    for (i = 0; i <= mx_lev; i++) is >> dt_level[i];
    for (i = 0; i <= mx_lev; i++) is >> n_cycle[i];
    for (i = 0; i <= mx_lev; i++) is >> level_steps[i];
    for (i = 0; i <= mx_lev; i++) is >> level_count[i];

      // set bndry conditions
    if (max_level > mx_lev) {
	for (i = mx_lev+1; i <= max_level; i++) {
	    int rat = MaxRefRatio(i-1);
	    dt_level[i] = dt_level[i-1]/float(rat);
              // NEED SUB_CYCLE
	    if (sub_cycle) {
		n_cycle[i] = rat;
		level_steps[i] = rat*level_steps[i-1];
	    } else {
		n_cycle[i] = 1;
		level_steps[i] = level_steps[i-1];
	    }
	    level_count[i] = 0;
	}
	if (!sub_cycle) {
	    for (i = 0; i <= max_level; i++) {
		dt_level[i] = dt_level[max_level];
	    }
	}
    }

    checkInput();

    for (int lev = 0; lev <= finest_level; lev++)
    {
	amr_level.set(lev,(*levelbld)());
	amr_level[lev].restart(*this,is);
    }

    RunStats::readStats(is);

    for (int lev = 0; lev <= finest_level; lev++)
	amr_level[lev].post_restart();
}

void
Amr::checkPoint ()
{
    aString ckfile = Concatenate(check_file_root, level_steps[0]);

    if (trace)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            cout << "CHECKPOINT: file = " << ckfile << NL;
        }
    }
    if (record_run_info)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            runlog << "CHECKPOINT: file = " << ckfile << NL;
        }
    }

    ofstream os;

#ifdef BL_T3E
    // ---------------- new buffer code
    int buffersize = 40960 * 32;
    signed char *iobuff = new signed char[buffersize];
    os.setbuf(iobuff, buffersize);
    // ------------ end new buffer code
#endif

    os.open(ckfile.c_str(), ios::out);

    os.precision(15);

    int i;
    if (ParallelDescriptor::IOProcessor())
    {
        os << BL_SPACEDIM  << '\n';
        os << cumtime      << '\n';
        os << max_level    << '\n';
        os << finest_level << '\n';

        // write out problem domain
        for (i = 0; i <= max_level; i++)
            os << geom[i] << "  ";
        os << '\n';
        for (i = 0; i < max_level; i++)
            os << ref_ratio[i] << "  ";
        os << '\n';

        for (i = 0; i <= max_level; i++)
            os << dt_level[i] << "  ";
        os << '\n';
        for (i = 0; i <= max_level; i++)
            os << n_cycle[i] << "  ";
        os << '\n';
        for (i = 0; i <= max_level; i++)
            os << level_steps[i] << "  ";
        os << '\n';
        for (i = 0; i <= max_level; i++)
            os << level_count[i] << "  ";
        os << '\n';
    }

    for (i = 0; i <= finest_level; i++)
	amr_level[i].checkPoint(os);

    RunStats::dumpStats(os);

    os.close();

#ifdef BL_T3E
    // ---------------- new buffer code
    delete [] iobuff;
    // ------------ end new buffer code
#endif
}


// ###################################################################
// ##### TIME_STEP
// ###################################################################
void Amr::timeStep(int level, REAL time, int iteration, int niter)
{
      // time to regrid?
    int lev_top = Min(finest_level, max_level-1);
    int i;
    for (i = level; i <= lev_top; i++) {
	if (level_count[i] >= regrid_int[i]) {
	    int old_finest = finest_level;
	    regrid(i,time);
            int k;
	    for (k = i; k <= finest_level; k++) {
		level_count[k] = 0;
	    }
	    if( old_finest < finest_level ){
		  // the new levels will not have valid time steps
		  // and iteration counts
		for (k = old_finest+1; k <= finest_level; k++) {
// Changed by Dan Wake 4/16/97
                    if (sub_cycle) {
		       dt_level[k] = dt_level[k-1]/REAL(MaxRefRatio(k-1));
		       n_cycle[k] = MaxRefRatio(k-1);
                    }
                    else {
                       dt_level[k] = dt_level[k-1] ;
                       n_cycle[k] = 1 ;
                    }
		}
	    }
	}
    }

      // advance grids at this level
    if (trace) {
      if(ParallelDescriptor::IOProcessor()) {
	cout << "ADVANCE grids at level " << level
	     << " with dt = " << dt_level[level] << NL;
      }
    }
    REAL dt_new = amr_level[level].advance(time,dt_level[level],
					   iteration,niter);
    if (iteration == 1) {
	dt_min[level] = dt_new;
    } else {
	dt_min[level] = Min(dt_min[level],dt_new);
    }
    level_steps[level]++;
    level_count[level]++;
    RunStats::addCells(level,amr_level[level].countCells());

   // advance grids at higher level
    if (level < finest_level) {
	if (sub_cycle) {
	    int lev_fine = level+1;
	    int ncycle = n_cycle[lev_fine];
	    for (i = 1; i <= ncycle; i++) {
		timeStep(lev_fine,time + (i-1)*dt_level[lev_fine],i,ncycle);
	    }
	} else {
	    int lev_fine = level+1;
	    timeStep(lev_fine,time,1,1);
	}
    }

      // perform any post_integration actions
    amr_level[level].post_timestep();
}


// ###################################################################
// ##### COARSE_TIME_STEP
// ###################################################################
// make a time step on the coarsest grid.
void Amr::coarseTimeStep(REAL stop_time){

      // compute new dt
    if (level_steps[0] > 0) {
	amr_level[0].computeNewDt(finest_level,sub_cycle,
				  n_cycle,ref_ratio,
				  dt_min,dt_level,stop_time);
    }

    timeStep(0,cumtime,1,1);
    cumtime += dt_level[0];

    amr_level[0].postCoarseTimeStep(cumtime);


    if (!silent) {
      if(ParallelDescriptor::IOProcessor()) {
	cout << "\nSTEP = " << level_steps[0] << " TIME = "
	     << cumtime << " DT = " << dt_level[0] << NL << NL;
      }
    }
    if (record_run_info) {
      if(ParallelDescriptor::IOProcessor()) {
	runlog << "STEP = " << level_steps[0] << " TIME = "
	       << cumtime << " DT = " << dt_level[0] << NL;
      }
    }

    int check_test = 0;
    if (check_per > 0.0) {
      int num_per = int((cumtime+.001*dt_level[0]) / check_per);
      REAL resid = cumtime - num_per * check_per;
      if (resid < .001*dt_level[0]) check_test = 1;
    }

    if ((check_int > 0 && level_steps[0] % check_int == 0) || 
        (check_test == 1)) {
	last_checkpoint = level_steps[0];
	checkPoint();
    }

    int plot_test = 0;
    if (plot_per > 0.0) {
      int num_per = int((cumtime+.001*dt_level[0]) / plot_per);
      REAL resid = cumtime - num_per * plot_per;
      if (resid < .001*dt_level[0]) plot_test = 1;
    }

    if ((plot_int > 0 && level_steps[0] % plot_int == 0) || 
        (plot_test == 1)) {
	last_plotfile = level_steps[0];
	writePlotFile();
    }

}

// ###################################################################
// ##### DEF_BASE_LEVEL
// ###################################################################
void
Amr::defBaseLevel(REAL strt_time) {
      // check that base domain has even number of zones in all
      // directions
    const BOX& domain = geom[0].Domain();
    const int* d_len = domain.length();
    int idir;
    for (idir = 0; idir < BL_SPACEDIM; idir++) {
	if (d_len[idir]%2 != 0) {
	    BoxLib::Error("defBaseLevel: domain not have even number of cells");
	}
    }

      // coarsening before we split the grids ensures that
      // each resulting grid will have an even number of
      // cells in each direction
    BoxArray  lev0(1);
    lev0.set(0, ::coarsen(domain,2));

   // now split up into list of grids within max_grid_size limit
    lev0.maxSize(max_grid_size/2);

      // now refine these boxes back to level 0
    lev0.refine(2);

      // now build level 0 grids
    amr_level.set(0,(*levelbld)(*this,0,geom[0],lev0,strt_time));

    lev0.clear();
      // now init level 0 grids with data

    amr_level[0].initData();
}

// ###################################################################
// ##### REGRID
// ###################################################################
void
Amr::regrid(int lbase, REAL time)
{
    int       new_finest;

    if (!silent) {
      if(ParallelDescriptor::IOProcessor()) {
	cout << "REGRID: at level lbase = " << lbase << NL;
      }
    }
    if (record_run_info) {
      if(ParallelDescriptor::IOProcessor()) {
	runlog << "REGRID: at level lbase = " << lbase << NL;
      }
    }

      // remove old-time grid space at highest level
    if (finest_level == max_level) {
	amr_level[finest_level].removeOldData();
    }

      // compute positions of new grids
    Array<BoxArray> new_grid_places(max_level+1);
    assert(new_grid_places.ready());
      // compute positions of new grids
    if(lbase<=Min(finest_level,max_level-1)) {
      grid_places(lbase,time,new_finest, new_grid_places);
    }

      // reclaim old-time grid space for all remain levels > lbase
    int lev;
    for (lev = lbase+1; lev <= finest_level; lev++) {
	amr_level[lev].removeOldData();
    }

      // reclaim all remaining storage for levels > new_finest
    for (lev = new_finest+1; lev <= finest_level; lev++) {
	amr_level.clear(lev);
    }

    finest_level = new_finest;

      // define the new grids from level lbase+1 up to new_finest
    for (lev = lbase+1; lev <= new_finest; lev++) {
	  // construct skeleton of new level
	AmrLevel *a = (*levelbld)(*this,lev,geom[lev],
				  new_grid_places[lev],cumtime);
	assert(a);

	  // init with data from old structure then remove old structure
	if (!amr_level.defined(lev)) {
	    a->init();
	} else {
	    a->init(amr_level[lev]);
	    amr_level.clear(lev);
	}

	  // install new structure
	amr_level.set(lev,a);
    }

      // build any additional data structures after grid generation
    for (lev = 0; lev <= new_finest; lev++) {
	amr_level[lev].post_regrid(lbase,new_finest);
    }

      // report creation of new grids
    if (!silent || record_run_info) {
        int lev;
	for (lev = lbase+1; lev <= finest_level; lev++) {
	    int numgrids= amr_level[lev].numGrids();
	    long ncells = amr_level[lev].countCells();
	    long ntot = geom[lev].Domain().numPts();
	    float frac = 100.0F*(float(ncells) / float(ntot));
	    if (!silent) {
	      if(ParallelDescriptor::IOProcessor()) {
		cout   << "   level " << lev << ": "
		       << numgrids << " grids, "
		       << ncells << " cells  = " << frac
		       << " % of domain" << NL;
	      }
	    }
	    if (record_run_info) {
	      if(ParallelDescriptor::IOProcessor()) {
		runlog << "   level " << lev << ": "
		       << numgrids << " grids, "
		       << ncells << " cells  = " << frac
		       << " % of domain" << NL;
	      }
	    }
	}
    }

    if (record_grid_info) {
      if(ParallelDescriptor::IOProcessor()) {
        if (lbase == 0) {
          gridlog << "STEP = " << level_steps[0] <<
                     " TIME = " << time <<
                     " : REGRID  with lbase = " << lbase << '\n';
        } else {
          gridlog << "TIME = " << time <<
                     " : REGRID  with lbase = " << lbase << '\n';
        }
	printGridInfo(gridlog,lbase+1,finest_level);
      }
    }

}

// ###################################################################
// ##### PRINT_GRID_INFO
// ###################################################################
// ------------------------------------------------------------------
void Amr::printGridInfo(ostream &os, int min_lev, int max_lev) {
    int lev;
    for (lev = min_lev; lev <= max_lev; lev++) {
	const BoxArray& bs = amr_level[lev].boxArray();
	int numgrid = bs.length();
	long ncells = amr_level[lev].countCells();
	long ntot = geom[lev].Domain().numPts();
	float frac = 100.0F*(float(ncells) / float(ntot));
	os << "  Level " << lev << "   " << numgrid << " grids  "
	   << ncells << " cells  " << frac << " % of domain" << '\n';
        int k;
	for (k = 0; k < numgrid; k++) {
	    const BOX &b = bs[k];
	    os << "  " << lev << ": " << b << "   ";
            int i;
	    for (i = 0; i < BL_SPACEDIM; i++) {
		os << b.length(i) << ' ';
	    }
	    os << '\n';
	}
    }
    os << NL;
}

// ##################################################################
// ## prof_periodic -- a helper function for computing periodic translates
// ## of BoxArray's.  ba will contain only box's in domain, plus their 
// ## periodic translates
// ##################################################################
void
proj_periodic( BoxDomain& bd, const Geometry& geom)
{
  Box domain( geom.Domain() );
  // blout will initially contain all of bd, periodic translates
  // will be added to it
  BoxList blout;  
  BoxDomainIterator bdi(bd);
  for( ; bdi; ++bdi){
    blout.add(bdi());
  }
  blout.intersect(domain);
  // blorig will be the original bd intersected with domain
  BoxList blorig(blout);

  int nist,njst,nkst;
  int niend,njend,nkend;
  nist = njst = nkst = 0;
  niend = njend = nkend = 0;
  D_TERM( nist , =njst , =nkst ) = -1;
  D_TERM( niend , =njend , =nkend ) = +1;

  int ri,rj,rk;
  for( ri=nist; ri<=niend; ri++ ){
    if( ri!=0 && !geom.isPeriodic(0) ) continue;
    if( ri!=0 && geom.isPeriodic(0) ) blorig.shift(0,ri*domain.length(0));
    for( rj=njst; rj<=njend; rj++ ){
      if( rj!=0 && !geom.isPeriodic(1) ) continue;
      if( rj!=0 && geom.isPeriodic(1) ) blorig.shift(1,rj*domain.length(1));
      for( rk=nkst; rk<=nkend; rk++ ){
	if( rk!=0 && !geom.isPeriodic(2) ) continue;
	if( rk!=0 && geom.isPeriodic(2) ) blorig.shift(2,rk*domain.length(2));

	BoxList tmp(blorig);
	tmp.intersect(domain);
	blout.join(tmp);
	

	if( rk!=0 && geom.isPeriodic(2) ) blorig.shift(2,-rk*domain.length(2));
      }
      if( rj!=0 && geom.isPeriodic(1) ) blorig.shift(1,-rj*domain.length(1));
    }
    if( ri!=0 && geom.isPeriodic(0) ) blorig.shift(0,-ri*domain.length(0));
  }
  bd.clear();
  bd.add(blout);
}


// ###################################################################
// ##### GRID_PLACES
// ###################################################################
// ------------------------------------------------------------------
void
Amr::grid_places(
    int         lbase,        //  => finest level that does not change
    REAL        time,
    int         &new_finest,  // <=> new finest level
    Array<BoxArray> &new_grids   // <=> new grid structure
    )

{
    int i;
    int  max_crse = Min(finest_level,max_level-1);

    if (grids_file != "") {
#define STRIP while( is.get() != '\n' )
	ifstream is(grids_file.c_str(),ios::in);
	if (is.fail()) {
	    cerr << "can't open file: " << grids_file << '\n';
	    BoxLib::Abort();
	}
	new_finest = Min(max_level,(finest_level+1));
	int in_finest;
	is >> in_finest;
	STRIP;
	new_finest = Min(new_finest,in_finest);
	int ngrid;
        int lev;
	for (lev = 1; lev <= new_finest; lev++) {
	    BoxList bl;
	    is >> ngrid;
	    STRIP;
	    for (i = 0; i < ngrid; i++) {
		BOX bx;
		is >> bx;
		STRIP;
		if (lev > lbase) {
		    bx.refine(ref_ratio[lev-1]);
		    if (bx.longside() > max_grid_size) {
			cerr << "Grid " << bx << " too large" << '\n';
			BoxLib::Abort();
		    }
		    bl.append(bx);
		}
	    }
	    if (lev > lbase) new_grids[lev].define(bl);
	}
	is.close();
	return;
#undef STRIP
    }

      // construct problem domain at each level
    Array<IntVect> bf_lev(max_level);    // blocking factor at each level
    Array<IntVect> rr_lev(max_level);
    Array<BOX> pc_domain(max_level); // prob domain coarsened by blocking_factor
    for (i = lbase; i <= max_crse; i++) {
      for( int n=0; n<BL_SPACEDIM; n++ ) {
	bf_lev[i][n] = Max(1,blocking_factor/ref_ratio[i][n]);
      }
    }
    for (i = lbase; i < max_crse; i++) {
      for( int n=0; n<BL_SPACEDIM; n++ ) {
	rr_lev[i][n] = (ref_ratio[i][n]*bf_lev[i][n])/bf_lev[i+1][n];
      }
    }
    for (i = lbase; i <= max_crse; i++) {
	pc_domain[i] = ::coarsen(geom[i].Domain(),bf_lev[i]);
    }

      // construct proper nesting domains
    Array<BoxDomain> p_n(max_level);   // proper nesting domain
    Array<BoxDomain> p_n_comp(max_level);   // complement of proper nesting domain

    BoxDomain fd1;
    const BoxArray &bbase = amr_level[lbase].boxArray();
    for ( i = 0; i < bbase.length(); i++) {
	BOX tmp(::coarsen(bbase[i],bf_lev[lbase]));
	fd1.add(tmp);
    }

    p_n_comp[lbase].complementIn(pc_domain[lbase],fd1);
    p_n_comp[lbase].accrete(n_proper);
    Geometry tmp2(pc_domain[lbase]);
    proj_periodic( p_n_comp[lbase], tmp2 );
    p_n_comp[lbase].minimize();
    p_n[lbase].complementIn(pc_domain[lbase],p_n_comp[lbase]);
    p_n[lbase].minimize();
    fd1.clear();

    for (i = lbase+1; i <= max_crse;  i++) {

      /* The following boxlist "iterate and add" loop replaces the 
	 original two lines

	p_n_comp[i] = p_n_comp[i-1];
	p_n_comp[i].refine(rr_lev[i-1]);

	which should be restored after a refine(IntVect) member function
	is added to the BoxDomain class
	*/

        BoxList bl;
	BoxDomainIterator bdi(p_n_comp[i-1]);
	while(bdi) {
	  bl.add(refine(bdi(),rr_lev[i-1]));
	  bdi++;
	}
	p_n_comp[i].clear();
	p_n_comp[i].add(bl);

	p_n_comp[i].accrete(n_proper);
	Geometry tmp3(pc_domain[i]);
	proj_periodic( p_n_comp[i], tmp3 );
	p_n[i].complementIn(pc_domain[i],p_n_comp[i]);
	p_n[i].minimize();
    }

      // now generate grids from finest level down
    new_finest = lbase;
    int levc;
    for (levc = max_crse; levc >= lbase; levc--) {
	int levf = levc+1;

	  // construct TagBoxArray with sufficient grow factor to contain
	  // new levels projected down to this level
	const BoxArray& old_grids = amr_level[levc].boxArray();
	int ngrow = 0;
	BoxArray ba_proj;
	if (levf < new_finest) {
	    BoxList blst(old_grids);
	    ba_proj.define(new_grids[levf+1]);
	    ba_proj.coarsen(ref_ratio[levf]);
	    ba_proj.grow(n_proper);
	    ba_proj.coarsen(ref_ratio[levc]);
	    while (!blst.contains(ba_proj)) {
		blst.accrete(1);
		ngrow++;
	    }
	}

	TagBoxArray tags(old_grids,n_error_buf[levc]+ngrow);
	amr_level[levc].errorEst(tags,TAGBOX::CLEAR,TAGBOX::SET,time);

	  // if new grids have been constructed above this level, project
	  // those grids down and tag cells on intersections to ensure
	  // proper nesting
	if (levf < new_finest) {
	    tags.setVal(ba_proj,TAGBOX::SET);
	}

	  // buffer error cells
	tags.buffer(n_error_buf[levc]);

	  // coarsen the taglist by blocking_factor
	int bl_max = 0;
	for( int n=0; n<BL_SPACEDIM; n++ ) {
	  bl_max = Max(bl_max,bf_lev[levc][n]);
	}
	if (bl_max > 1) tags.coarsen(bf_lev[levc]);

	// map tagged points through periodic boundaries, if any
	Geometry tmpgeom(pc_domain[levc]);
	tags.mapPeriodic(tmpgeom);

	  // merge tagged points on overlap, remove redundant tags
	tags.mergeUnique();

	  // remove cells outside proper nesting domain for this level
	tags.setVal(p_n_comp[levc],TAGBOX::CLEAR);

	  // create initial cluster containing all tagged points
	Array<INTVECT> *pts = tags.colate();

	tags.clear();

	if (pts->length() == 0) {
	      // no tagged points, no new level
	      // do cleanup
	    delete pts;
	} else {

	      // created new level, now generate efficient grids
	    new_finest = Max(new_finest,levf);

	      // construct initial cluster
	    ClusterList clist(pts);
	    clist.chop(grid_eff);
	    clist.intersect(p_n[levc]);

	      // efficient properly nested CLUSTERs have been constructed
	      // now generate list of grids at level levf
	    BoxList new_bx;
	    clist.boxList(new_bx);

	      // MLW: (7/23/91) merge adjacent boxes
	    int nmerged = new_bx.minimize();

	    IntVect lratio = ref_ratio[levc]*bf_lev[levc];

	    IntVect largest_grid_size;
	    for (int n = 0; n < BL_SPACEDIM; n++ ) {
	      largest_grid_size[n] = max_grid_size / lratio[n];
	    }

	      // ensure new grid boxes are at most max_grid_size in
	      // each index direction
//	    new_bx.maxSize(largest_grid_size);
	    maxSizeBox(new_bx,largest_grid_size);

	    // refine up to levf

//	    new_bx.refine(lratio);
            refineBoxList(new_bx,lratio);

	    if (!new_bx.isDisjoint()) {
		cout << "WARNING: new grids at level "<<levf
		     << " not disjoint" << NL
		     << new_bx << NL;
	    }

	    new_grids[levf].define(new_bx);

	}  // new level?

    }  // loop over levels
}  // end grid_places()


// ###################################################################
// ##### BLD_FINE_LEVELS
// ###################################################################
void  Amr::bldFineLevels(REAL strt_time)
{

    finest_level = 0;
    int more_levels = true;

    int num_levs = max_level+1;
    Array<BoxArray> new_grid_places(num_levs);
    while (more_levels) {

	  // get new grid placement
	int     new_finest;
	grid_places(finest_level,strt_time,new_finest,new_grid_places);

	if (new_finest <= finest_level) {
	      // no new level
	    more_levels = false;

	} else {

	      // create a new level and link with others
	    int levf = new_finest;
	    finest_level = new_finest;

	    amr_level.set(levf,(*levelbld)(*this,levf,geom[levf],
					   new_grid_places[levf],
					   strt_time));

	      // init data on new level
	    amr_level[levf].initData();

	    more_levels = (finest_level < max_level);

	}  // if we created another level

    }  // more_levels

}

extern "C" void maxSizeBox (BoxList& bx_list, IntVect& block_size)
{
    for (BoxListIterator bli(bx_list); bli; ++bli)
    {
        const IntVect& ivlen = bli().length();
        const int* len       = ivlen.getVect();
        for (int i = 0; i < SpaceDim; i++)
        {
            if (len[i] > block_size[i])
            {
                //
                // Reduce by powers of 2.
                //
                int ratio = 1;
                int bs    = block_size[i];
                int nlen  = len[i];
                while ((bs%2==0) && (nlen%2==0))
                {
                    ratio *= 2;
                    bs /=2;
                    nlen /=2;
                }
                //
                // Determine number and size of (coarsened) cuts.
                //
                int numblk = nlen/bs + (nlen%bs ? 1 : 0);
                int size   = nlen/numblk;
                int extra  = nlen%numblk;
                //
                // Number of cuts = number of blocks - 1.
                //
                for (int k = 0; k < numblk-1; k++)
                {
                    //
                    // Compute size of this chunk, expand by power of 2.
                    //
                    int ksize = (k < extra) ? size+1 : size;
                    ksize *= ratio;
                    //
                    // Chop from high end.
                    //
                    IntVect iv = bli().bigEnd();
                    int pos    = iv[i] - ksize + 1;

                    bx_list.append(bx_list[bli].chop(i,pos));
                }
            }
        }
        //
        // b has been chopped down to size and pieces split off
        // have been added to the end of the list so that they
        // can be checked for splitting (in other directions) later.
        //
    }
}

extern "C" void refineBoxList (BoxList& bx_list, IntVect& lratio)
{
    BoxList bl;
    BoxListIterator bli(bx_list);
    while(bli) {
      bl.add(refine(bli(),lratio));
      bli++;
    }
    bx_list.clear();
    bx_list = bl;
}
