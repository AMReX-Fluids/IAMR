// ---------------------------------------------------------------
// This tool reads a Checkpoint file, interpolate to a finer grid
// and writes it out again.
// It is inspired from the Embiggening tool in the Castro distribution
// Validated for IAMR 2D and multi-levels.
// Please report bugs and problems to Emmanuel Motheau (emotheau@lbl.gov)
// November 14 2020
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
#include <AMReX_Extrapolater.H>

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
std::string interp_kind;


Vector<int> nsets_save(1);

VisMF::How how = VisMF::OneFilePerCPU;

int ngrow = 1; // WARNING: this works for IAMR with no EB, but this number may be different for EB

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

FakeAmr fakeAmr_src;
FakeAmr fakeAmr_trgt;

// ---------------------------------------------------------------
static void ScanArguments() {
    ParmParse pp;

    if(pp.contains("checkin")) {
      pp.get("checkin", CheckFileIn);
    }
    if(pp.contains("checkout")) {
      pp.get("checkout", CheckFileOut);
    }
    if(pp.contains("verbose")) {
      pp.get("verbose", verbose);
    }
    if(pp.contains("user_ratio")) {
      pp.get("user_ratio", user_ratio);
    }
    if(pp.contains("interp_kind")) {
      pp.get("interp_kind", interp_kind);
    }

    if (user_ratio != 2 && user_ratio != 4)
       amrex::Abort("user_ratio must be 2 or 4");

    if (interp_kind != "refine" && interp_kind != "coarsen" )
       amrex::Abort("interp_kind must be set to `refine` or `coarsen`");

}

// ---------------------------------------------------------------
static void PrintUsage (char *progName) {
    cout << "Usage: " << progName << " checkin=filename "
         << "checkout=outfilename "  
         << "user_ratio= 2 or 4 "  
         << "interp_kind= refine or coarsen " 
         << "[verbose=trueorfalse]" << endl;
    exit(1);
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

    is >> fakeAmr_src.cumtime;
    int mx_lev;
    is >> mx_lev;
    is >> fakeAmr_src.finest_level;

    fakeAmr_src.geom.resize(mx_lev + 1);
    fakeAmr_src.ref_ratio.resize(mx_lev);
    fakeAmr_src.dt_level.resize(mx_lev + 1);
    fakeAmr_src.dt_min.resize(mx_lev + 1);
    fakeAmr_src.n_cycle.resize(mx_lev + 1);
    fakeAmr_src.level_steps.resize(mx_lev + 1);
    fakeAmr_src.level_count.resize(mx_lev + 1);
    fakeAmr_src.fakeAmrLevels.resize(mx_lev + 1);

    // READING GEOM, REF_RATIO, DT_LEVEL
    for (i = 0; i <= mx_lev; i++) {
      is >> fakeAmr_src.geom[i];
    }

    for (i = 1; i <=  mx_lev; i++) {
      is >> fakeAmr_src.ref_ratio[i-1];
    }
    for (i = 0; i <= mx_lev; i++) {
      is >> fakeAmr_src.dt_level[i];
    }

    if (new_checkpoint_format) {
      for (i = 0; i <= mx_lev; i++) is >> fakeAmr_src.dt_min[i];
    } else {
      for (i = 0; i <= mx_lev; i++) fakeAmr_src.dt_min[i] = fakeAmr_src.dt_level[i];
    }

    // READING N_CYCLE, LEVEL_STEPS, LEVEL_COUNT
    for (i = 0; i <= mx_lev; i++) {
      is >> fakeAmr_src.n_cycle[i];
    }
    for (i = 0; i <= mx_lev; i++) {
      is >> fakeAmr_src.level_steps[i];
    }
    for (i = 0; i <= mx_lev; i++) {
      is >> fakeAmr_src.level_count[i];
    }

    int ndesc_save;

    // READ LEVEL DATA
    for(int lev(0); lev <= fakeAmr_src.finest_level; ++lev) {
    
      if (ParallelDescriptor::IOProcessor()) {
        if (lev == 0) {    
           std::cout << " " << std::endl;
           std::cout << " **************************************** " << std::endl;
           std::cout << " " << std::endl;
        }
           std::cout << "Old checkpoint level    " << lev << std::endl;
           std::cout << " ... domain is       " << fakeAmr_src.geom[lev].Domain() << std::endl;
           std::cout << " ...     dx is       " << fakeAmr_src.geom[lev].CellSize()[0] << std::endl;
           std::cout << "  " << std::endl;
      }
 
      FakeAmrLevel &falRef = fakeAmr_src.fakeAmrLevels[lev];

      is >> falRef.level;
      is >> falRef.geom;

      falRef.fine_ratio = IntVect::TheUnitVector();
      falRef.fine_ratio.scale(-1);
      falRef.crse_ratio = IntVect::TheUnitVector();
      falRef.crse_ratio.scale(-1);

      if(falRef.level > 0) {
        falRef.crse_ratio = fakeAmr_src.ref_ratio[falRef.level-1];
      }
      if(falRef.level < mx_lev) {
        falRef.fine_ratio = fakeAmr_src.ref_ratio[falRef.level];
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

	int max_level(fakeAmr_trgt.finest_level);
        HeaderFile << CheckPointVersion << '\n'
                   << BL_SPACEDIM       << '\n'
                   << fakeAmr_trgt.cumtime           << '\n'
                   << max_level                 << '\n'
                   << fakeAmr_trgt.finest_level      << '\n';

        //
        // Write out problem domain.
        //
        for (i = 0; i <= max_level; i++) HeaderFile << fakeAmr_trgt.geom[i]        << ' ';
        HeaderFile << '\n';
        for (i = 0; i < max_level; i++)  HeaderFile << fakeAmr_trgt.ref_ratio[i]   << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << fakeAmr_trgt.dt_level[i]    << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << fakeAmr_trgt.dt_min[i]      << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << fakeAmr_trgt.n_cycle[i]     << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << fakeAmr_trgt.level_steps[i] << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << fakeAmr_trgt.level_count[i] << ' ';
        HeaderFile << '\n';
    }

    for(int lev(0); lev <= fakeAmr_trgt.finest_level; ++lev) {

      std::ostream &os = HeaderFile;
      FakeAmrLevel &falRef = fakeAmr_trgt.fakeAmrLevels[lev];
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
          std::cout << " ... domain is       " << fakeAmr_trgt.geom[lev].Domain() << std::endl;
          std::cout << " ...     dx is       " << fakeAmr_trgt.geom[lev].CellSize()[0] << std::endl;
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

  fakeAmr_trgt = fakeAmr_src;

  int mx_lev = fakeAmr_trgt.finest_level;

  for (int lev = 0; lev <= mx_lev; lev++)
  {  
    FakeAmrLevel &falRef_src = fakeAmr_src.fakeAmrLevels[lev];
    Box bx = fakeAmr_src.geom[lev].Domain();
    BoxArray ba = falRef_src.grids;
    DistributionMapping dm{ba};

// HERE WE REFINE THE LEVEL IN FINE STRUCTURE
    Box          domain_trgt(fakeAmr_trgt.geom[lev].Domain());
    RealBox prob_domain_trgt(fakeAmr_trgt.geom[lev].ProbDomain());
    int coord_trgt = fakeAmr_trgt.geom[lev].Coord();

    if (interp_kind == "refine"){
       domain_trgt.refine(user_ratio);
    } else {
       domain_trgt.coarsen(user_ratio);
    }

    fakeAmr_trgt.geom[lev].define(domain_trgt,&prob_domain_trgt,coord_trgt);

    if (interp_kind == "refine"){
      fakeAmr_trgt.dt_level[lev] = fakeAmr_trgt.dt_level[lev] / user_ratio;
      fakeAmr_trgt.dt_min[lev] = fakeAmr_trgt.dt_min[lev] / user_ratio;
    } else {
      fakeAmr_trgt.dt_level[lev] = fakeAmr_trgt.dt_level[lev] * user_ratio;
      fakeAmr_trgt.dt_min[lev] = fakeAmr_trgt.dt_min[lev] * user_ratio;
    }

    const GpuArray<int,AMREX_SPACEDIM>& is_periodic_array = fakeAmr_src.geom[lev].isPeriodicArray();
    fakeAmr_trgt.geom[lev].setPeriodicity({{AMREX_D_DECL(is_periodic_array[0],is_periodic_array[1],is_periodic_array[2])}});

// NOW WE WORK ON DATA THAT ARE IN THE fine FakeAmrLevel structure
    FakeAmrLevel &falRef_trgt = fakeAmr_trgt.fakeAmrLevels[lev];
    BoxArray new_grids = falRef_trgt.grids;

    if (interp_kind == "refine"){
      new_grids.refine(user_ratio);
    } else {
      new_grids.coarsen(user_ratio);
    }


    falRef_trgt.grids = new_grids;
    falRef_trgt.geom.define(domain_trgt,&prob_domain_trgt,coord_trgt);

    falRef_trgt.geom.setPeriodicity({{AMREX_D_DECL(is_periodic_array[0],is_periodic_array[1],is_periodic_array[2])}});

    DistributionMapping dm_trgt{new_grids};

    for (int n = 0; n < falRef_src.state.size(); n++){

// Assuming that OldState and NewState have the same number of components
      int ncomps = (falRef_src.state[n].old_data)->nComp();

      BoxArray new_grids_state = falRef_trgt.state[n].grids;
      BoxArray save_grids_state = falRef_trgt.state[n].grids;

      if (interp_kind == "refine"){
        new_grids_state.refine(user_ratio);
      } else {
        new_grids_state.coarsen(user_ratio);
      }

      falRef_trgt.state[n].grids = new_grids_state;

       if (interp_kind == "refine"){
         falRef_trgt.state[n].domain.refine(user_ratio);
       } else {      
         falRef_trgt.state[n].domain.coarsen(user_ratio);
       }

      MultiFab * NewData_src = new MultiFab(save_grids_state,dm,ncomps,ngrow);
      MultiFab * OldData_src = new MultiFab(save_grids_state,dm,ncomps,ngrow);
      NewData_src -> setVal(10.); 
      OldData_src -> setVal(10.);

      NewData_src -> copy(*(falRef_src.state[n].new_data),0,0,ncomps,0,ngrow);
      OldData_src -> copy(*(falRef_src.state[n].old_data),0,0,ncomps,0,ngrow);

      MultiFab * NewData_trgt = new MultiFab(new_grids_state,dm_trgt,ncomps,ngrow);
      MultiFab * OldData_trgt = new MultiFab(new_grids_state,dm_trgt,ncomps,ngrow);
      NewData_trgt->setVal(10.);
      OldData_trgt->setVal(10.);

      IntVect new_ratio(AMREX_D_DECL(user_ratio,user_ratio,user_ratio));
      const IntVect& rr = new_ratio;

      if (interp_kind == "refine"){

        Interpolater*  interpolater = &cell_cons_interp;
        if ( n == 1)  interpolater = &node_bilinear_interp;

        const Geometry& fgeom = falRef_trgt.geom;
        const Geometry& cgeom = falRef_src.geom;

        for (MFIter mfi(*NewData_trgt); mfi.isValid(); ++mfi)
        {
          FArrayBox& ffab = (*NewData_trgt)[mfi];
          const FArrayBox& cfab = (*NewData_src)[mfi];
          const Box&  bx   = mfi.tilebox();
          Vector<BCRec> bx_bcrec(ncomps);

          interpolater->interp(cfab,0,ffab,0,ncomps,bx,rr,
                               cgeom,fgeom,bx_bcrec,0,0,RunOn::Host);
        }

        for (MFIter mfi(*OldData_trgt); mfi.isValid(); ++mfi)
        {
          FArrayBox& ffab = (*OldData_trgt)[mfi];
          const FArrayBox& cfab = (*OldData_src)[mfi];
          const Box&  bx   = mfi.tilebox();
          Vector<BCRec> bx_bcrec(ncomps);

          interpolater->interp(cfab,0,ffab,0,ncomps,bx,rr,
                               cgeom,fgeom,bx_bcrec,0,0,RunOn::Host);

        }

        NewData_trgt->FillBoundary(fgeom.periodicity());
        OldData_trgt->FillBoundary(fgeom.periodicity());

      } else {
      
        amrex::average_down (*NewData_src, *NewData_trgt,
                             0,  ncomps, new_ratio);

        amrex::average_down (*OldData_src, *OldData_trgt,
                             0,  ncomps, new_ratio);

      }

      falRef_trgt.state[n].new_data = NewData_trgt;
      falRef_trgt.state[n].old_data = OldData_trgt;

    }
  }
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

    // Read in the original checkpoint directory  
    ReadCheckpointFile(CheckFileIn);

    // Interpolate to a finest grid
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
