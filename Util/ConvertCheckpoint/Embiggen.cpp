
// This is the version that reads a Checkpoint file
// and writes it out again.
// ---------------------------------------------------------------
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <string>

#ifndef WIN32
#include <unistd.h>
#endif

#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>
#include <AMReX_Geometry.H>
#include <AMReX_StateDescriptor.H>
#include <AMReX_StateData.H>
#include <AMReX_BCRec.H>
#include <AMReX_LevelBld.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_Interpolater.H>


using namespace amrex;

#define VSHOWVAL(verbose, val) { if(verbose && \
                                   ParallelDescriptor::IOProcessor()) { \
                                   cout << #val << " = " << val << endl; } }
LevelBld *getLevelBld() {
  return 0;
}

#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif

using std::cout;
using std::cerr;
using std::endl;

using namespace amrex;

std::string CheckFileIn;
std::string CheckFileOut;
int nFiles(64);
bool verbose(true);
int      user_ratio(1);
int   max_grid_size(4096);
const std::string CheckPointVersion = "CheckPointVersion_1.0";

Vector<int> nsets_save(1);

VisMF::How how = VisMF::OneFilePerCPU;


// ---------------------------------------------------------------
struct FakeStateData {
    struct TimeInterval {
        Real start, stop;
    };
    const StateDescriptor *desc;
    Box domain;
    BoxArray grids;
    TimeInterval new_time;
    TimeInterval old_time;
    MultiFab *new_data;
    MultiFab *old_data;
    Vector< Vector<BCRec> > bc;
};


struct FakeAmrLevel {
    int level;                        // AMR level (0 is coarsest).
    Geometry geom;                    // Geom at this level.
    BoxArray grids;                   // Cell-centered locations of grids.
    IntVect crse_ratio;               // Refinement ratio to coarser level.
    IntVect fine_ratio;               // Refinement ratio to finer level.
    Vector<FakeStateData> state;       // Array of state data.
    Vector<FakeStateData> new_state;   // Array of new state data.
};


struct FakeAmr {
  int                  finest_level;
  Real                 cumtime;
  Vector<Real>          dt_level;
  Vector<int>           level_steps;
  Vector<int>           level_count;
  Vector<int>           n_cycle;
  Vector<Real>          dt_min;
  Vector<IntVect>       ref_ratio;
  Vector<Geometry>      geom;
  Vector<FakeAmrLevel> fakeAmrLevels;
};

FakeAmr fakeAmr;
FakeAmr fakeAmr_fine;

// ---------------------------------------------------------------
static void ScanArguments() {
    ParmParse pp;

    if(pp.contains("checkin")) {
      pp.get("checkin", CheckFileIn);
    }
    if(pp.contains("checkout")) {
      pp.get("checkout", CheckFileOut);
    }
    if(pp.contains("nfiles")) {
      pp.get("nfiles", nFiles);
    }
    if(pp.contains("verbose")) {
      pp.get("verbose", verbose);
    }
    if(pp.contains("user_ratio")) {
      pp.get("user_ratio", user_ratio);
    }

    if (user_ratio != 2 && user_ratio != 4)
       amrex::Abort("user_ratio must be 2 or 4");

}

// ---------------------------------------------------------------
static void PrintUsage (char *progName) {
    cout << "Usage: " << progName << " checkin=filename "
         << "checkout=outfilename "  
         << "user_ratio= 2 or 4 "   
         << "[nfiles=nfilesout] "
         << "[verbose=trueorfalse]" << endl;
    exit(1);
}

// ---------------------------------------------------------------
static void
set_bcrec_new (Vector<BCRec>  &bcrec,
               int             ncomp)

{
   for (int n = 0; n < ncomp; n++) {
      for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
      {
         bcrec[n].setLo(dir,INT_DIR);
         bcrec[n].setHi(dir,INT_DIR);
      }
   }
}

static void ReadCheckpointFile(const std::string& fileName) {
    int i;
    std::string File = fileName;

    File += '/';
    File += "Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ifstream is;
    is.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    is.open(File.c_str(), std::ios::in);
    if( ! is.good()) {
        amrex::FileOpenFailed(File);
    }

    //
    // Read global data.
    //
    // Attempt to differentiate between old and new CheckPointFiles.
    //
    int         spdim;
    bool        new_checkpoint_format = false;
    std::string first_line;

    std::getline(is,first_line);

    if (first_line == CheckPointVersion) {
        new_checkpoint_format = true;
        is >> spdim;
    } else {
        spdim = atoi(first_line.c_str());
    }

    if(spdim != BL_SPACEDIM) {
        cerr << "Amr::restart(): bad spacedim = " << spdim << '\n';
        amrex::Abort();
    }

    is >> fakeAmr.cumtime;
    int mx_lev;
    is >> mx_lev;
    is >> fakeAmr.finest_level;

    if(ParallelDescriptor::IOProcessor()) 
       std::cout << "previous finest_lev is " << fakeAmr.finest_level <<  std::endl;
std::cout << "DEBUG mx_lev is " << mx_lev << std::endl;


//    // ADDING LEVELS 
//    int n = num_new_levels;
//    mx_lev = mx_lev + n;
//    fakeAmr.finest_level = fakeAmr.finest_level + n;
//
//    if(ParallelDescriptor::IOProcessor()) {
//       std::cout << "     new finest_lev is " << fakeAmr.finest_level << std::endl;
//       std::cout << "previous     mx_lev is " << mx_lev-n << std::endl;
//       std::cout << "     new     mx_lev is " << mx_lev << std::endl;
//    }



    fakeAmr.geom.resize(mx_lev + 1);
    fakeAmr.ref_ratio.resize(mx_lev);
    fakeAmr.dt_level.resize(mx_lev + 1);
    fakeAmr.dt_min.resize(mx_lev + 1);
    fakeAmr.n_cycle.resize(mx_lev + 1);
    fakeAmr.level_steps.resize(mx_lev + 1);
    fakeAmr.level_count.resize(mx_lev + 1);
    fakeAmr.fakeAmrLevels.resize(mx_lev + 1);

    // READING GEOM, REF_RATIO, DT_LEVEL
    for (i = 0; i <= mx_lev; i++) {
      is >> fakeAmr.geom[i];
    }

    if(ParallelDescriptor::IOProcessor()) {
       std::cout << " " << std::endl;
       for (i = 0; i <= mx_lev; i++) {
          std::cout << "Old checkpoint level    " << i << std::endl;
          std::cout << " ... domain is       " << fakeAmr.geom[i].Domain() << std::endl;
          std::cout << " ...     dx is       " << fakeAmr.geom[i].CellSize()[0] << std::endl;
          std::cout << "  " << std::endl;
       }
    }

//    // Make sure current domain is divisible by 2*ref_ratio so length of coarsened domain is even
//    Box dom0(fakeAmr.geom[0].Domain());
//    for (int d = 0; d < BL_SPACEDIM; d++)
//    {
//      int dlen = dom0.size()[d];
//      int scaled = dlen / (2*ref_ratio);
//      if ( (scaled * 2 * ref_ratio) != dlen )
//        amrex::Abort("must have domain divisible by 2*ref_ratio");
//    }
//
//    if (grown_factor <= 1)  
//        amrex::Abort("must have grown_factor > 1");

    for (i = 1; i <=  mx_lev; i++) {
      is >> fakeAmr.ref_ratio[i-1];
std::cout << "DEBUG i=    " << i << std::endl;
std::cout << " DEBUG ref_ratio is       " << fakeAmr.ref_ratio[i-1] << std::endl;
    }
    for (i = 0; i <= mx_lev; i++) {
      is >> fakeAmr.dt_level[i];
std::cout << "DEBUG i=    " << i << std::endl;
std::cout << " DEBUG dt_level is       " << fakeAmr.dt_level[i] << std::endl;
    }

    if(ParallelDescriptor::IOProcessor()) {
       std::cout << " " << std::endl;
       for (i = 0; i <= mx_lev; i++) {
          std::cout << "Old checkpoint level    " << i << std::endl;
          if (i > 0){
            std::cout << " ... ref_ratio is       " << fakeAmr.ref_ratio[i-1] << std::endl;
          }
          std::cout << " ...  dt_level is       " << fakeAmr.dt_level[i] << std::endl;
          std::cout << "  " << std::endl;
       }
    }

// UNTIL HERE IT LOOKS OK FOR 2AMR LEVELS

/*
    Box          domain(fakeAmr.geom[1].Domain());
    RealBox prob_domain(fakeAmr.geom[1].ProbDomain());
    coord = fakeAmr.geom[1].Coord();

    // Define domain for new levels
    domain.coarsen(ref_ratio);
    fakeAmr.geom[0].define(domain,&prob_domain,coord);

    // Define ref_ratio for new levels
    fakeAmr.ref_ratio[0] = ref_ratio * IntVect::TheUnitVector();

    // Define dt_level for new levels
    fakeAmr.dt_level[0] = fakeAmr.dt_level[1] * ref_ratio;
*/


    if (new_checkpoint_format) {
      for (i = 0; i <= mx_lev; i++) is >> fakeAmr.dt_min[i];
//      fakeAmr.dt_min[0] = fakeAmr.dt_min[1] * ref_ratio;
    } else {
      for (i = 0; i <= mx_lev; i++) fakeAmr.dt_min[i] = fakeAmr.dt_level[i];
    }

for (i = 0; i <= mx_lev; i++) {
std::cout << "DEBUG i=    " << i << std::endl;
std::cout << " DEBUG fakeAmr.dt_min is       " << fakeAmr.dt_min[i] << std::endl;
    }

/*
Real test;
is >> test;
std::cout << " DEBUG  DEBUG      " << test << std::endl;
is >> test;
std::cout << " DEBUG  DEBUG      " << test << std::endl;
is >> test;
std::cout << " DEBUG  DEBUG      " << test << std::endl;
is >> test;
std::cout << " DEBUG  DEBUG      " << test << std::endl;
is >> test;
std::cout << " DEBUG  DEBUG      " << test << std::endl;
*/
    // READING N_CYCLE, LEVEL_STEPS, LEVEL_COUNT
    for (i = 0; i <= mx_lev; i++) {
is >> fakeAmr.n_cycle[i];
//std::cout << "DEBUG i=    " << i << std::endl;
//std::cout << " DEBUG n_cycle is       " << fakeAmr.n_cycle[i] << std::endl;
    }

    for (i = 0; i <= mx_lev; i++) {
      is >> fakeAmr.level_steps[i];
//std::cout << "DEBUG i=    " << i << std::endl;
//std::cout << " DEBUG n_cycle is       " << fakeAmr.level_steps[i] << std::endl;
    }
    for (i = 0; i <= mx_lev; i++) {
      is >> fakeAmr.level_count[i];
//std::cout << "DEBUG i=    " << i << std::endl;
//std::cout << " DEBUG n_cycle is       " << fakeAmr.level_count[i] << std::endl;
    }



    if(ParallelDescriptor::IOProcessor()) {
       std::cout << " " << std::endl;
       for (i = 0; i <= mx_lev; i++) {
          std::cout << "Old checkpoint level    " << i << std::endl;
            std::cout << " ... n_cycle is       " << fakeAmr.n_cycle[i] << std::endl;
          std::cout << " ...  level_steps is       " << fakeAmr.level_steps[i] << std::endl;
          std::cout << " ...  level_count is       " << fakeAmr.level_count[i] << std::endl;
          std::cout << "  " << std::endl;
       }
    }


// UNTIL HERE IT LOOKS OK FOR 2AMR LEVELS


/*


    // ADDING LEVELS

    // n_cycle is always equal to 1 at the coarsest level 
    fakeAmr.n_cycle[0] = 1;

    // At the old coarsest level, which is now level 1, we set n_cycle to ref_ratio 
    fakeAmr.n_cycle[1] = ref_ratio;

    fakeAmr.level_steps[0] = fakeAmr.level_steps[1] / ref_ratio;  
    if ( (fakeAmr.level_steps[0]*ref_ratio) != fakeAmr.level_steps[1] ) 
       amrex::Abort("Number of steps in original checkpoint must be divisible by ref_ratio");

    // level_count is how many steps we've taken at this level since the last regrid
    if (fakeAmr.level_count[1] == fakeAmr.level_steps[1]) 
    {
       fakeAmr.level_count[0] = fakeAmr.level_steps[0];

    // this is actually wrong but should work for now
    } else {
       fakeAmr.level_count[0] = std::min(fakeAmr.level_count[1],fakeAmr.level_steps[0]);;
    }
*/


    int ndesc_save;

    // READ LEVEL DATA
    for(int lev(0); lev <= fakeAmr.finest_level; ++lev) {
     
std::cout << "DEBUG lev=    " << lev << std::endl;
 
      FakeAmrLevel &falRef = fakeAmr.fakeAmrLevels[lev];

      is >> falRef.level;
//      falRef.level = falRef.level + n;

      is >> falRef.geom;

      falRef.fine_ratio = IntVect::TheUnitVector();
      falRef.fine_ratio.scale(-1);
      falRef.crse_ratio = IntVect::TheUnitVector();
      falRef.crse_ratio.scale(-1);

      if(falRef.level > 0) {
        falRef.crse_ratio = fakeAmr.ref_ratio[falRef.level-1];
      }
      if(falRef.level < mx_lev) {
        falRef.fine_ratio = fakeAmr.ref_ratio[falRef.level];
      }

      falRef.grids.readFrom(is);

      int nstate;
      is >> nstate;
      int ndesc = nstate;

      // This should be the same at all levels
      ndesc_save = ndesc;

      // ndesc depends on which descriptor so we store a value for each
      if (lev == 1) nsets_save.resize(ndesc_save);

      falRef.state.resize(ndesc);
      falRef.new_state.resize(ndesc);

      for(int ii = 0; ii < ndesc; ii++) {
        // ******* StateDescriptor::restart

        is >> falRef.state[ii].domain;

        falRef.state[ii].grids.readFrom(is);

        is >> falRef.state[ii].old_time.start;
        is >> falRef.state[ii].old_time.stop;
        is >> falRef.state[ii].new_time.start;
        is >> falRef.state[ii].new_time.stop;

        int nsets;
        is >> nsets;

        nsets_save[ii] = nsets;

        falRef.state[ii].old_data = 0;
        falRef.state[ii].new_data = 0;

        std::string mf_name;
        std::string FullPathName;

        // This reads the "new" data, if it's there
        if (nsets >= 1) {
           falRef.state[ii].new_data = new MultiFab;
           is >> mf_name;
           // Note that mf_name is relative to the Header file.
           // We need to prepend the name of the fileName directory.
           FullPathName = fileName;
           if( ! fileName.empty() && fileName[fileName.length()-1] != '/') {
             FullPathName += '/';
           }
           FullPathName += mf_name;
           VisMF::Read(*(falRef.state[ii].new_data), FullPathName);
        }

        // This reads the "old" data, if it's there
        if (nsets == 2) {
          falRef.state[ii].old_data = new MultiFab;
          is >> mf_name;
          // Note that mf_name is relative to the Header file.
          // We need to prepend the name of the fileName directory.
          FullPathName = fileName;
          if( ! fileName.empty() && fileName[fileName.length()-1] != '/') {
            FullPathName += '/';
	  }
          FullPathName += mf_name;
          VisMF::Read(*(falRef.state[ii].old_data), FullPathName);
        }



      }


    }

//VisMF::Write(*(fakeAmr.fakeAmrLevels[0].state[0].new_data), "initial_new_data_lev0");
//VisMF::Write(*(fakeAmr.fakeAmrLevels[1].state[0].new_data), "initial_new_data_lev1");
//VisMF::Write(*(fakeAmr.fakeAmrLevels[0].state[0].old_data), "initial_old_data_lev0");
//VisMF::Write(*(fakeAmr.fakeAmrLevels[1].state[0].old_data), "initial_old_data_lev1");



// UNTIL HERE IT LOOKS OK FOR 2AMR LEVELS




}

// ---------------------------------------------------------------

static void WriteCheckpointFile(const std::string& inFileName, const std::string &outFileName) {
    VisMF::SetNOutFiles(nFiles);
    // In checkpoint files always write out FABs in NATIVE format.
    FABio::Format thePrevFormat = FArrayBox::getFormat();
    FArrayBox::setFormat(FABio::FAB_NATIVE);

    const std::string ckfile = outFileName;

    // Only the I/O processor makes the directory if it doesn't already exist.
    if(ParallelDescriptor::IOProcessor()) {
      if( ! amrex::UtilCreateDirectory(ckfile, 0755)) {
        amrex::CreateDirectoryFailed(ckfile);
      }
    }
    // Force other processors to wait till directory is built.
    ParallelDescriptor::Barrier();

    // Write the main header file.

    std::string HeaderFileName = ckfile + "/Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    int old_prec(0), i;

    if(ParallelDescriptor::IOProcessor()) {
        // Only the IOProcessor() writes to the header file.
        HeaderFile.open(HeaderFileName.c_str(),
	                std::ios::out|std::ios::trunc|std::ios::binary);

        if( ! HeaderFile.good()) {
          amrex::FileOpenFailed(HeaderFileName);
	}

        old_prec = HeaderFile.precision(15);

	int max_level(fakeAmr_fine.finest_level);
        HeaderFile << CheckPointVersion << '\n'
                   << BL_SPACEDIM       << '\n'
                   << fakeAmr_fine.cumtime           << '\n'
                   << max_level                 << '\n'
                   << fakeAmr_fine.finest_level      << '\n';
        //
        // Write out problem domain.
        //
        for (i = 0; i <= max_level; i++) HeaderFile << fakeAmr_fine.geom[i]        << ' ';
        HeaderFile << '\n';
        for (i = 0; i < max_level; i++)  HeaderFile << fakeAmr_fine.ref_ratio[i]   << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << fakeAmr_fine.dt_level[i]    << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << fakeAmr_fine.dt_min[i]      << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << fakeAmr_fine.n_cycle[i]     << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << fakeAmr_fine.level_steps[i] << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << fakeAmr_fine.level_count[i] << ' ';
        HeaderFile << '\n';
    }

    for(int lev(0); lev <= fakeAmr_fine.finest_level; ++lev) {

      std::ostream &os = HeaderFile;
      FakeAmrLevel &falRef = fakeAmr_fine.fakeAmrLevels[lev];
      int ndesc = falRef.state.size();

      // Build directory to hold the MultiFabs in the StateData at this level.
      char buf[64];
      sprintf(buf, "Level_%d", lev);
      std::string Level = buf;

      // Now for the full pathname of that directory.
      std::string FullPath = ckfile;
      if( ! FullPath.empty() && FullPath[FullPath.length()-1] != '/') {
        FullPath += '/';
      }
      FullPath += Level;

      // Only the I/O processor makes the directory if it doesn't already exist.
      if(ParallelDescriptor::IOProcessor()) { if( ! amrex::UtilCreateDirectory(FullPath, 0755)) {
          amrex::CreateDirectoryFailed(FullPath);
        }
      }
      // Force other processors to wait till directory is built.
      ParallelDescriptor::Barrier();

      if(ParallelDescriptor::IOProcessor()) {
        os << lev << '\n' << falRef.geom  << '\n';
        falRef.grids.writeOn(os);
        os << ndesc << '\n';
      }
      //
      // Output state data.
      //
      for(int i(0); i < ndesc; ++i) {
        //
        // Now build the full relative pathname of the StateData.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        // There is only one MultiFab written out at each level in HyperCLaw.
        //
        std::string PathNameInHeader = Level;
        sprintf(buf, "/SD_%d", i);
        PathNameInHeader += buf;
        std::string FullPathName = FullPath;
        FullPathName += buf;
        // ++++++++++++ state[i].checkPoint(PathNameInHeader, FullPathName, os, how);
          static const std::string NewSuffix("_New_MF");
          static const std::string OldSuffix("_Old_MF");
          const std::string name(PathNameInHeader);
          const std::string fullpathname(FullPathName);

          bool dump_old(true);
          if(dump_old == true && falRef.state[i].old_data == 0) {
            dump_old = false;
          }

          if(ParallelDescriptor::IOProcessor()) {
            // The relative name gets written to the Header file.
            std::string mf_name_old = name;
	    mf_name_old += OldSuffix;
            std::string mf_name_new = name;
	    mf_name_new += NewSuffix;

            os << falRef.state[i].domain << '\n';

            falRef.state[i].grids.writeOn(os);

            os << falRef.state[i].old_time.start << '\n'
               << falRef.state[i].old_time.stop  << '\n'
               << falRef.state[i].new_time.start << '\n'
               << falRef.state[i].new_time.stop  << '\n';

            if (nsets_save[i] > 0) {
               if(dump_old) {
                 os << 2 << '\n' << mf_name_new << '\n' << mf_name_old << '\n';
               } else {
                 os << 1 << '\n' << mf_name_new << '\n';
               }
            } else {
              os << 0 << '\n';
            }

          }

          if (nsets_save[i] > 0) {
             BL_ASSERT(falRef.state[i].new_data);
             std::string mf_fullpath_new = fullpathname;
             mf_fullpath_new += NewSuffix;
             VisMF::Write(*(falRef.state[i].new_data),mf_fullpath_new,how);
          }

          if (nsets_save[i] > 1) {
            BL_ASSERT(dump_old);
            BL_ASSERT(falRef.state[i].old_data);
            std::string mf_fullpath_old = fullpathname;
	    mf_fullpath_old += OldSuffix;
            VisMF::Write(*(falRef.state[i].old_data),mf_fullpath_old,how);
          }
          // ++++++++++++
      }
      // ========================

       if (ParallelDescriptor::IOProcessor()) {
          if (lev == 0) {
             std::cout << " " << std::endl;
             std::cout << " **************************************** " << std::endl;
             std::cout << " " << std::endl;
          }
          std::cout << "New checkpoint level    " << lev << std::endl;
          std::cout << " ... domain is       " << fakeAmr_fine.geom[lev].Domain() << std::endl;
          std::cout << " ...     dx is       " << fakeAmr_fine.geom[lev].CellSize()[0] << std::endl;
          std::cout << "  " << std::endl;
       }

    }

    if(ParallelDescriptor::IOProcessor()) {
        HeaderFile.precision(old_prec);

        if( ! HeaderFile.good()) {
          amrex::Error("Amr::checkpoint() failed");
	}
    }

    FArrayBox::setFormat(thePrevFormat);
}



// ---------------------------------------------------------------
static void ConvertData() {


 fakeAmr_fine = fakeAmr;

int mx_lev = fakeAmr_fine.finest_level;

if(ParallelDescriptor::IOProcessor()) {
       std::cout << " " << std::endl;
       std::cout << " DEBUG fine AMR finest_lev is    " << fakeAmr_fine.finest_level <<  std::endl;
       std::cout << " DEBUG fine AMR cumtime is       " << fakeAmr_fine.cumtime << std::endl;
       std::cout << "  " << std::endl;

       for (int i = 0; i <= mx_lev; i++) {
          std::cout << "DEBUG fine AMR  level       " << i << std::endl;
          std::cout << " DEBUG dt_level is          " << fakeAmr_fine.dt_level[i] << std::endl;
          std::cout << " DEBUG level_steps is       " << fakeAmr_fine.level_steps[i] << std::endl;
          std::cout << " DEBUG level_count is       " << fakeAmr_fine.level_count[i] << std::endl;
          std::cout << " DEBUG n_cycle is           " << fakeAmr_fine.n_cycle[i] << std::endl;
          std::cout << " DEBUG dt_min is            " << fakeAmr_fine.dt_min[i] << std::endl;
          if (i > 0 ) std::cout << " DEBUG ref_ratio is         " << fakeAmr_fine.ref_ratio[i-1] << std::endl;
          std::cout << " DEBUG domain is       " << fakeAmr_fine.geom[i].Domain() << std::endl;
          std::cout << " DEBUG     dx is       " << fakeAmr_fine.geom[i].CellSize()[0] << std::endl;
          std::cout << "  " << std::endl;
       }
    }



for (int lev = 0; lev <= mx_lev; lev++)
   {  


   FakeAmrLevel &falRef_coarse = fakeAmr.fakeAmrLevels[lev];
   Box bx = fakeAmr.geom[lev].Domain();
   BoxArray ba = falRef_coarse.grids;
   DistributionMapping dm{ba};


// HERE WE REFINE THE LEVEL IN FINE STRUCTURE
    Box          domain_fine(fakeAmr_fine.geom[lev].Domain());
    RealBox prob_domain_fine(fakeAmr_fine.geom[lev].ProbDomain());
    int coord_fine = fakeAmr_fine.geom[lev].Coord();
    domain_fine.refine(user_ratio);
    fakeAmr_fine.geom[lev].define(domain_fine,&prob_domain_fine,coord_fine);
    fakeAmr_fine.dt_level[lev] = fakeAmr_fine.dt_level[lev] / user_ratio;
    fakeAmr_fine.dt_min[lev] = fakeAmr_fine.dt_min[lev] / user_ratio;

 std::cout << "  " << std::endl;
 std::cout << " DEBUG REFINEMENT " << std::endl;
if(ParallelDescriptor::IOProcessor()) {
       std::cout << " " << std::endl;
       std::cout << " DEBUG fine AMR finest_lev is    " << fakeAmr_fine.finest_level <<  std::endl;
       std::cout << " DEBUG fine AMR cumtime is       " << fakeAmr_fine.cumtime << std::endl;
       std::cout << "  " << std::endl;

       for (int i = 0; i <= mx_lev; i++) {
          std::cout << "DEBUG fine AMR  level       " << i << std::endl;
          std::cout << " DEBUG dt_level is          " << fakeAmr_fine.dt_level[i] << std::endl;
          std::cout << " DEBUG level_steps is       " << fakeAmr_fine.level_steps[i] << std::endl;
          std::cout << " DEBUG level_count is       " << fakeAmr_fine.level_count[i] << std::endl;
          std::cout << " DEBUG n_cycle is           " << fakeAmr_fine.n_cycle[i] << std::endl;
          std::cout << " DEBUG dt_min is            " << fakeAmr_fine.dt_min[i] << std::endl;
          if (i > 0 ) std::cout << " DEBUG ref_ratio is         " << fakeAmr_fine.ref_ratio[i-1] << std::endl;
          std::cout << " DEBUG domain is       " << fakeAmr_fine.geom[i].Domain() << std::endl;
          std::cout << " DEBUG     dx is       " << fakeAmr_fine.geom[i].CellSize()[0] << std::endl;
          std::cout << "  " << std::endl;
       }
    }

// NOW WE WORK ON DATA THAT ARE IN THE FakeAmrLevel structure
FakeAmrLevel &falRef_fine = fakeAmr_fine.fakeAmrLevels[lev];
BoxArray new_grids = falRef_fine.grids;
new_grids.refine(user_ratio);
falRef_fine.grids = new_grids;
falRef_fine.geom.define(domain_fine,&prob_domain_fine,coord_fine);


std::cout << " DEBUG fine_level  domain is         " << falRef_fine.geom.Domain() << std::endl;
std::cout << " DEBUG fine_level     dx is      " << falRef_fine.geom.CellSize()[0] << std::endl;
std::cout << " DEBUG fine_level     grids is      " << falRef_fine.grids << std::endl;
std::cout << " DEBUG nstate      is      " << falRef_fine.state.size() << std::endl;


DistributionMapping dm_fine{new_grids};


for (int n = 0; n < falRef_coarse.state.size(); n++){

// Assuming that OldState and NewState have the same number of components
   int ncomps = (falRef_coarse.state[n].old_data)->nComp();

   MultiFab * NewData_coarse = new MultiFab(ba,dm,ncomps,0);
   MultiFab * OldData_coarse = new MultiFab(ba,dm,ncomps,0);
   NewData_coarse -> setVal(0.); 
   OldData_coarse -> setVal(0.);

//   falRef_coarse.state[n].new_data = newNewData;
   
    NewData_coarse -> copy(*(falRef_coarse.state[n].new_data),0,0,ncomps);
    OldData_coarse -> copy(*(falRef_coarse.state[n].old_data),0,0,ncomps);

   falRef_fine.state[n].domain = domain_fine;
   falRef_fine.state[n].grids = falRef_fine.grids;

   MultiFab * NewData_fine = new MultiFab(new_grids,dm_fine,ncomps,0);
   MultiFab * OldData_fine = new MultiFab(new_grids,dm_fine,ncomps,0);
   NewData_fine->setVal(0.);
   OldData_fine->setVal(0.);

Interpolater*  interpolater = &pc_interp;


const Geometry& fgeom = falRef_fine.geom;
const Geometry& cgeom = falRef_coarse.geom;
IntVect new_ratio(AMREX_D_DECL(user_ratio,user_ratio,user_ratio));
const IntVect& rr = new_ratio; 

for (MFIter mfi(*NewData_fine); mfi.isValid(); ++mfi)
      {
         FArrayBox& ffab = (*NewData_fine)[mfi];
         const FArrayBox& cfab = (*NewData_coarse)[mfi];
         const Box&  bx   = mfi.tilebox();
 Vector<BCRec> bx_bcrec(ncomps);

         interpolater->interp(cfab,0,ffab,0,ncomps,bx,rr,
                              cgeom,fgeom,bx_bcrec,0,0,RunOn::Host);

}

for (MFIter mfi(*OldData_fine); mfi.isValid(); ++mfi)
      {
         FArrayBox& ffab = (*OldData_fine)[mfi];
         const FArrayBox& cfab = (*OldData_coarse)[mfi];
         const Box&  bx   = mfi.tilebox();
 Vector<BCRec> bx_bcrec(ncomps);

         interpolater->interp(cfab,0,ffab,0,ncomps,bx,rr,
                              cgeom,fgeom,bx_bcrec,0,0,RunOn::Host);

}

falRef_fine.state[n].new_data = NewData_fine;
falRef_fine.state[n].old_data = OldData_fine;

VisMF::Write(*NewData_fine,"pouet_fine");

}

}

VisMF::Write(*(fakeAmr_fine.fakeAmrLevels[0].state[0].new_data), "fine_new_data_lev0");


}





// ---------------------------------------------------------------
int main(int argc, char *argv[]) {
    amrex::Initialize(argc,argv);

    if(argc < 3) {
      PrintUsage(argv[0]);
    }

    ScanArguments();

    if(verbose && ParallelDescriptor::IOProcessor()) {
      cout << " " << std::endl;
      cout << "Reading from old checkpoint file: " <<  CheckFileIn << endl;
      cout << " " << std::endl;
    }

    // Read in the original checkpoint directory and add a coarser level covering the same domain 
    ReadCheckpointFile(CheckFileIn);

    // Enlarge the new level 0
    ConvertData();

    // Write out the new checkpoint directory
    WriteCheckpointFile(CheckFileIn, CheckFileOut);

    if(verbose && ParallelDescriptor::IOProcessor()) {
      cout << " " << std::endl;
      cout << "Finished writing to new checkpoint file: " <<  CheckFileOut << endl;
      cout << " " << std::endl;
    }

    amrex::Finalize();
}
// ---------------------------------------------------------------
// ---------------------------------------------------------------
