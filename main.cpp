//
// $Id: main.cpp,v 1.15 1997-12-11 23:30:35 lijewski Exp $
//

#ifdef BL_ARCH_CRAY
#ifdef BL_USE_DOUBLE
#error "DOUBLE PRECISION NOT ALLOWED ON CRAY"
#endif
#endif

#ifndef        WIN32
#include <unistd.h>
#endif

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
#include <AmrLevel.H>

#ifdef BL_USE_NEW_HFILES
#include <new>
using std::setprecision;
using std::set_new_handler;
#else
#include <new.h>
#endif

const int NPROCS = 1;

static
void 
print_usage (int,
             char* argv[])
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

static
void
OutOfMemory ()
{
    BoxLib::Error("Sorry, out of memory, bye ...");
}

int
main (int   argc,
      char* argv[])
{
    //
    // Make sure to catch new failures.
    //
    set_new_handler(OutOfMemory);

    TRACER("amr");
    cout << setprecision(10);

    if (argc < 2) print_usage(argc,argv);
    if (argv[1][0] == '-')
    {
        cerr << "input file must be first argument\n";
        print_usage(argc, argv);
    }
    int  max_step;
    Real stop_time;

    ParmParse pp(argc-2,argv+2,NULL,argv[1]); 

    int nprocs = NPROCS; pp.query("nprocs", nprocs);
#ifndef BL_USE_BSP
    if (nprocs > 1)
    {
      cerr << "Error in main:  multiple processors specified with "
           << "code compiled without a parallel library.\n";
      exit(-1);
    }
#endif
    StartParallel(nprocs);

#ifndef        WIN32
    int sleeptime = 0; pp.query("sleep", sleeptime);
    sleep(sleeptime);
#endif

    max_step  = 0;    pp.query("max_step",max_step);
    stop_time = 0.0;  pp.query("stop_time",stop_time);

    Amr *amrptr = new Amr;

    if (amrptr == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);

    amrptr->init();

    while (amrptr->okToContinue()           &&
           amrptr->levelSteps(0) < max_step &&
           amrptr->cumTime() < stop_time)
    {
        amrptr->coarseTimeStep(stop_time);
    }

    delete amrptr;
    //
    // This MUST follow the above delete as ~Amr() may dump files to disk.
    //
    RunStats::report(cout);

    EndParallel();

    return 0;
}
