//
// $Id: main.cpp,v 1.39 2001-08-01 21:51:01 lijewski Exp $
//

#include <cstdio>

#include <CArena.H>
#include <REAL.H>
#include <Utility.H>
#include <IntVect.H>
#include <Box.H>
#include <Amr.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <AmrLevel.H>
#include <TagBox.H>
#ifdef BL_USE_MPI
#include <hgparallel.H>
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

void PrintBoxArray (const BoxArray& ba) { std::cout << ba << std::endl; }

void PrintBoxDomain (const BoxDomain& bd) { std::cout << bd << std::endl; }

void
PrintTagBox (const TagBox& tb)
{
    const Box& bx = tb.box();

    long count = 0;

    std::cout << "TagBox: box = " << bx << ":\n";

    for (IntVect p = bx.smallEnd(); p <= bx.bigEnd(); bx.next(p))
    {
        if (!(tb(p) == TagBox::CLEAR))
        {
            count++;
            std::cout << p << ' ';
        }
    }

    std::cout << "Total tagged cells = " << count << std::endl;
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

    std::cout << "Total tagged cells = " << count << std::endl;
}

void
PrintTagBoxArray (const TagBoxArray& tba)
{
    long count = 0;

    std::cout << "TagBoxArray:\n";

    for (int i = 0; i < tba.size(); i++)
    {
        const Box& bx = tba[i].box();

        std::cout << "\ti = " << i << ", box = " << bx << ":\n";

        for (IntVect p = bx.smallEnd(); p <= bx.bigEnd(); bx.next(p))
        {
            if (!(tba[i](p) == TagBox::CLEAR))
            {
                count++;
                std::cout << p << ' ';
            }
        }

        std::cout << '\n';
    }

    std::cout << "Total tagged cells = " << count << std::endl;
}

void
TagBoxArrayCount (const TagBoxArray& tba)
{
    long count = 0;

    for (int i = 0; i < tba.size(); i++)
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

    std::cout << "Total tagged cells = " << count << std::endl;
}
#endif /*NDEBUG*/

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);

    const Real run_strt = ParallelDescriptor::second();

    int  max_step;
    Real strt_time;
    Real stop_time;

    ParmParse pp;
#ifdef BL_USE_MPI
    //
    // Initialize some Holy Grail junk.
    //
    HG::MPI_init();
#endif

    max_step  = -1;    
    strt_time =  0.0;  
    stop_time = -1.0;  

    pp.query("max_step",max_step);
    pp.query("strt_time",strt_time);
    pp.query("stop_time",stop_time);

    if (strt_time < 0.0)
        BoxLib::Abort("MUST SPECIFY a non-negative strt_time");

    if (max_step < 0 && stop_time < 0.0)
    {
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
    // Close down the Holy Grail junk.
    //
#ifdef BL_USE_MPI
    HG::MPI_finish();
#endif

    if (CArena* arena = dynamic_cast<CArena*>(The_FAB_Arena))
    {
        //
        // We're using a CArena -- output some FAB memory stats.
        // This'll output total # of bytes of heap space in the Arena.
        // It's actually the high water mark of heap space required by FABs.
        //
        char buf[256];

        sprintf(buf,
                "CPU(%d): Heap Space (bytes) used by Coalescing FAB Arena: %ld",
                ParallelDescriptor::MyProc(),
                arena->heap_space_used());

        std::cout << buf << std::endl;
    }

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_stop = ParallelDescriptor::second() - run_strt;

    ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "Run time = " << run_stop << std::endl;

    BoxLib::Finalize();

    return 0;
}
