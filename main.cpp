//
// $Id: main.cpp,v 1.21 1998-04-27 16:52:44 lijewski Exp $
//

#ifdef BL_ARCH_CRAY
#ifdef BL_USE_DOUBLE
#error "DOUBLE PRECISION NOT ALLOWED ON CRAY"
#endif
#endif

#ifndef WIN32
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
#ifndef WIN32
using std::set_new_handler;
#endif
#else
#include <new.h>
#endif

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

int
main (int   argc,
      char* argv[])
{
    //
    // Make sure to catch new failures.
    //
#ifndef WIN32
    set_new_handler(Utility::OutOfMemory);
#endif

#ifdef BL_USE_MPI
    int nprocs = 1;
    ParallelDescriptor::StartParallel(nprocs,&argc,&argv);
#endif

    cout << setprecision(10);

    if (argc < 2)
        print_usage(argc,argv);

    if (argv[1][0] == '-')
    {
        cerr << "input file must be first argument\n";
        print_usage(argc, argv);
    }
    int  max_step;
    Real stop_time;

    ParmParse pp(argc-2,argv+2,NULL,argv[1]); 

#ifndef BL_USE_MPI
    int nprocs = 1; pp.query("nprocs", nprocs);
    ParallelDescriptor::StartParallel(nprocs,&argc,&argv);
#endif
    //
    // Initialize random seed after we're running in parallel.
    //
    Utility::InitRandom(ParallelDescriptor::MyProc() + 1);
    //
    // Instantiate after we're running in Parallel.
    //
    TRACER("amr");

#ifndef WIN32
    int sleeptime = 0; pp.query("sleep", sleeptime);
    sleep(sleeptime);
#endif

    max_step  = 0;    pp.query("max_step",max_step);
    stop_time = 0.0;  pp.query("stop_time",stop_time);

    Amr* amrptr = new Amr;

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

    ParallelDescriptor::EndParallel();

    return 0;
}
