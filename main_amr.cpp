
//
// $Id: main_amr.cpp,v 1.11 1998-04-01 17:01:54 car Exp $
//

#ifdef BL_ARCH_CRAY
#  ifdef BL_USE_DOUBLE
DOUBLE PRECISION NOT ALLOWED ON CRAY
#  endif
#endif

#ifdef BL_USE_NEW_HFILES
#include <cstdlib>
#include <new>
using std::setprecision;
using std::set_new_handler;
#else
#include <new.h>
#include <stdlib.h>
#endif

#include <REAL.H>
#include <Misc.H>
#include <Utility.H>
#include <IntVect.H>
#include <Box.H>
#include <Amr.H>
#include <RunStats.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <AmrLevel.H>

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
    //
    // Make sure to catch new failures.
    //
    set_new_handler(Utility::OutOfMemory);

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
    Real stop_time;
    char* in_file = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,in_file); 


    int nprocs = NPROCS;
    pp.query("nprocs", nprocs);

    ParallelDescriptor::StartParallel(nprocs, &argc, &argv);

    int sleeptime = 0;
    pp.query("sleep", sleeptime);
    sleep(sleeptime);

    max_step = 0;    
    pp.query("max_step",max_step);
    stop_time =0.0;  pp.query("stop_time",stop_time);
    FArrayBox::init();

    // -------------------------------------------------
    // -----   construct objects
    // -------------------------------------------------
    Amr  *amrptr = new Amr;
    // -------------------------------------------------
    // -----   initialization section
    // -------------------------------------------------
      // stats.init();
    RunStats::init();
    amrptr->init();

//    cout << setprecision(10);

    // -------------------------------------------------
    // -----   loop until finished
    // -------------------------------------------------
    while (amrptr->okToContinue() &&
           amrptr->levelSteps(0) < max_step &&
           amrptr->cumTime() < stop_time)
    {
        amrptr->coarseTimeStep(stop_time);
    }

    // -------------------------------------------------
    // -----   final business
    // -------------------------------------------------

      // cout << "Local timing stats (since restart)" << '\n';
      // stats.report(cout);
      // cout << "Global timing stats (entire run)" << '\n';
      // global_stats += stats;
      // global_stats.report(cout);
      RunStats::report(cout);

    // -------------------------------------------------
    // -----   delete memory
    // -------------------------------------------------
    delete amrptr;

    ParallelDescriptor::EndParallel();

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
