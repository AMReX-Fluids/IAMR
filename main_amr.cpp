//BL_COPYRIGHT_NOTICE

//
// $Id: main_amr.cpp,v 1.16 1998-11-18 21:06:14 lijewski Exp $
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

    BoxLib::Abort("Exiting.");
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

    ParallelDescriptor::StartParallel(&argc,&argv);

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
    char* in_file = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,in_file); 
    //
    // Initialize random seed after we're running in parallel.
    //
    Utility::InitRandom(ParallelDescriptor::MyProc() + 1);

#ifndef WIN32
    int sleeptime = 0; pp.query("sleep", sleeptime);
    sleep(sleeptime);
#endif

    max_step = 0;    pp.query("max_step",max_step);
    stop_time =0.0;  pp.query("stop_time",stop_time);

    Amr* amrptr = new Amr;

    amrptr->init();

    while (amrptr->okToContinue()           &&
           amrptr->levelSteps(0) < max_step &&
           amrptr->cumTime() < stop_time)
    {
        amrptr->coarseTimeStep(stop_time);
    }

    RunStats::report(cout);

    delete amrptr;

    ParallelDescriptor::EndParallel();

    return 0;
}
