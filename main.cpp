//BL_COPYRIGHT_NOTICE

//
// $Id: main.cpp,v 1.30 1999-04-01 18:15:17 lijewski Exp $
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
#include <IntVect.H>
#include <Box.H>
#include <Amr.H>
#include <RunStats.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <AmrLevel.H>
#include <hgparallel.h>

#ifdef BL_USE_NEW_HFILES
#include <new>
using std::setprecision;
#ifndef WIN32
using std::set_new_handler;
#endif
#else
#include <new.h>
#endif

#ifndef NDEBUG

extern "C"
{
    void PrintBoxArray (const BoxArray& ba);
    void PrintBoxDomain (const BoxDomain& bd);
    void PrintTagBox (const TagBox& tb);
    void PrintTagBoxArray (const TagBoxArray& tba);
    void TagBoxCount (const TagBox& tb);
    void TagBoxArrayCount (const TagBoxArray& tba);
}

void PrintBoxArray (const BoxArray& ba) { cout << ba << endl; }

void PrintBoxDomain (const BoxDomain& bd) { cout << bd << endl; }

void
PrintTagBox (const TagBox& tb)
{
    const Box& bx = tb.box();

    long count = 0;

    cout << "TagBox: box = " << bx << ":\n";

    for (IntVect p = bx.smallEnd(); p <= bx.bigEnd(); bx.next(p))
    {
        if (!(tb(p) == TagBox::CLEAR))
        {
            count++;
            cout << p << ' ';
        }
    }

    cout << "Total tagged cells = " << count << endl;
}

void
TagBoxCount (const TagBox& tb)
{
    const Box& bx = tb.box();

    long count = 0;

    for (IntVect p = bx.smallEnd(); p <= bx.bigEnd(); bx.next(p))
    {
        if (!(tb(p) == TagBox::CLEAR))
        {
            count++;
        }
    }

    cout << "Total tagged cells = " << count << endl;
}

void
PrintTagBoxArray (const TagBoxArray& tba)
{
    long count = 0;

    cout << "TagBoxArray:\n";

    for (int i = 0; i < tba.length(); i++)
    {
        const Box& bx = tba[i].box();

        cout << "\ti = " << i << ", box = " << bx << ":\n";

        for (IntVect p = bx.smallEnd(); p <= bx.bigEnd(); bx.next(p))
        {
            if (!(tba[i](p) == TagBox::CLEAR))
            {
                count++;
                cout << p << ' ';
            }
        }

        cout << '\n';
    }

    cout << "Total tagged cells = " << count << endl;
}

void
TagBoxArrayCount (const TagBoxArray& tba)
{
    long count = 0;

    for (int i = 0; i < tba.length(); i++)
    {
        const Box& bx = tba[i].box();

        for (IntVect p = bx.smallEnd(); p <= bx.bigEnd(); bx.next(p))
        {
            if (!(tba[i](p) == TagBox::CLEAR))
            {
                count++;
            }
        }
    }

    cout << "Total tagged cells = " << count << endl;
}
#endif /*NDEBUG*/

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
    Real strt_time;
    Real stop_time;

    ParmParse pp(argc-2,argv+2,NULL,argv[1]); 
    //
    // Initialize random seed after we're running in parallel.
    //
    Utility::InitRandom(ParallelDescriptor::MyProc() + 1);
    //
    // Initialize some Holy Grail junk.
    //
    HG::MPI_init();
    //
    // Instantiate after we're running in Parallel.
    //
#ifndef WIN32
    int sleeptime = 0; pp.query("sleep", sleeptime);
    sleep(sleeptime);
#endif

    max_step  = -1;    
    strt_time =  0.0;  
    stop_time = -1.0;  

    pp.query("max_step",max_step);
    pp.query("strt_time",strt_time);
    pp.query("stop_time",stop_time);

    if ( strt_time < 0.0 ) {
      BoxLib::Abort(
        "MUST SPECIFY a non-negative strt_time");
    }

    if ( max_step < 0 && stop_time < 0.0 ) {
      BoxLib::Abort(
       "Exiting because neither max_step nor stop_time is non-negative.");
    }

    Amr* amrptr = new Amr;

    amrptr->init(strt_time,stop_time);

    while ( amrptr->okToContinue()           &&
           (amrptr->levelSteps(0) < max_step || max_step < 0) &&
           (amrptr->cumTime() < stop_time || stop_time < 0.0) )
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
