#define  _Main_C_ "%W% %G%"

#ifdef BL_ARCH_CRAY
#  ifdef BL_USE_DOUBLE
DOUBLE PRECISION NOT ALLOWED ON CRAY
#  endif
#endif

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <string.h>
#include <fstream.h>
#include <math.h>
#include <REAL.H>
#include <Misc.H>
#include <Utility.H>
#include <Tracer.H>
#include <IntVect.H>
#include <Box.H>
#include <Amr.H>
#include <RunStats.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#ifndef	WIN32
#include <unistd.h>
#endif

#include <AmrLevel.H>

// stupid thing to fix template instantiation bug
#include <AliasedDPtr.H>

void
boing()
{
    AliasedDPtr<REAL> adp(1);
}


#if BL_USE_WINDOWS&&(BL_SPACEDIM==2)
#   define HAS_CONTOUR  1
#   include <Contour.H>
#endif
#if HAS_RASTER && (BL_SPACEDIM==2 || BL_SPACEDIM == 3)
#   include <Raster.H>
#endif

// functions called
void   print_usage(int, char **);

// global (extern) objects
// RunStats stats;
// RunStats global_stats;

const int NPROCS = 1;

// ###################################################################
// ##### MAIN PROGRAM
// ###################################################################
main(int argc, char *argv[])
{
    TRACER("amr");
    cout << setprecision(10);

    // -------------------------------------------------
    // -----   parse command line and input file
    // -------------------------------------------------
    // infile must be first
    if (argc < 2) print_usage(argc,argv);
    if (argv[1][0] == '-') {
        cerr << "input file must be first argument\n";
        print_usage(argc, argv);
    }
    int  max_step;
    REAL stop_time;
    char* in_file = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,in_file); 


    int nprocs = NPROCS;
    pp.query("nprocs", nprocs);
#ifndef BL_USE_BSP
    if(nprocs > 1) {
      cerr << "Error in main:  multiple processors specified with "
	   << "code compiled without a parallel library." << endl;
      exit(-1);
    }
#endif
    StartParallel(nprocs);

#ifndef	WIN32
    int sleeptime = 0;
    pp.query("sleep", sleeptime);
    sleep(sleeptime);
#endif

    max_step = 0;    
    pp.query("max_step",max_step);
    stop_time =0.0;  pp.query("stop_time",stop_time);
    FARRAYBOX::init();

    // -------------------------------------------------
    // -----   construct objects
    // -------------------------------------------------
    Amr  *amrptr = new Amr;
#   ifdef HAS_CONTOUR
       Contour *contourptr;
       if(ParallelDescriptor::IOProcessor()) {
         contourptr = new Contour(*amrptr);
       }
       int n_contour = 0;
#   endif
#   ifdef HAS_RASTER
       Raster *rasterptr;
       if(ParallelDescriptor::IOProcessor()) {
         rasterptr = new Raster(*amrptr);
       }
       int n_raster = 0;
#   endif

    // -------------------------------------------------
    // -----   initialization section
    // -------------------------------------------------
      // stats.init();
    RunStats::init();
    amrptr->init();

//    cout << setprecision(10);

#   ifdef HAS_RASTER
       // dump raster plots
       n_raster = rasterptr->draw(amrptr->cumTime(),amrptr->levelSteps(0));
#   endif
#   ifdef HAS_CONTOUR
       // draw contour graphics on the fly
       if(ParallelDescriptor::IOProcessor()) {
         n_contour = contourptr->draw(amrptr->cumTime(),amrptr->levelSteps(0));
       }
#   endif

    // -------------------------------------------------
    // -----   loop until finished
    // -------------------------------------------------
    while (amrptr->okToContinue() &&
           amrptr->levelSteps(0) < max_step &&
	   amrptr->cumTime() < stop_time) {

        // do a timestep
        amrptr->coarseTimeStep(stop_time);

#       ifdef HAS_CONTOUR
           // draw contour graphics on the fly
         if(ParallelDescriptor::IOProcessor()) {
           n_contour = contourptr->draw(amrptr->cumTime(),
	                                amrptr->levelSteps(0));
         }
#       endif
#       ifdef HAS_RASTER
           // dump raster plots
         if(ParallelDescriptor::IOProcessor()) {
           n_raster = rasterptr->draw(amrptr->cumTime(),
	                              amrptr->levelSteps(0));
	 }
#       endif

    }

    // -------------------------------------------------
    // -----   final business
    // -------------------------------------------------
#   ifdef HAS_CONTOUR
       // dump final contours
       if (n_contour == 0) {
         if(ParallelDescriptor::IOProcessor()) {
          contourptr->draw(amrptr->cumTime(),amrptr->levelSteps(0), 1);
	 }
       }
#   endif
#   ifdef HAS_RASTER
       // dump final raster plots
       if (n_raster == 0) {
         if(ParallelDescriptor::IOProcessor()) {
          rasterptr->draw(amrptr->cumTime(),amrptr->levelSteps(0), 1);
	 }
       }
#   endif

      // cout << "Local timing stats (since restart)" << endl;
      // stats.report(cout);
      // cout << "Global timing stats (entire run)" << endl;
      // global_stats += stats;
      // global_stats.report(cout);
      RunStats::report(cout);

    // -------------------------------------------------
    // -----   delete memory
    // -------------------------------------------------
#   ifdef HAS_CONTOUR
     if(ParallelDescriptor::IOProcessor()) {
       delete contourptr;
     }
#   endif
#   ifdef HAS_RASTER
     if(ParallelDescriptor::IOProcessor()) {
       delete rasterptr;
     }
#   endif
    delete amrptr;

    EndParallel();

    // stop execution
    return(0);
}  // end main


// ###################################################################
// ##### PRINT_USAGE
// ###################################################################
void 
print_usage(int, char *argv[])
{
    cerr << "usage:\n";
    cerr << argv[0] << " infile [options] \n\tOptions:\n";
    cerr << "\t     [<root>.]<var>  = <val_list>\n";
    cerr << "\tor  -[<root>.]<var>\n";
    cerr << "\t where:\n";
    cerr << "\t    <root>     =  class name of variable\n";
    cerr << "\t    <var>      =  variable name\n";
    cerr << "\t    <val_list> =  list of values\n";

    ParallelDescriptor::Abort("Exiting.");
}
