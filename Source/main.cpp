
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

#ifndef MG_USE_F90_SOLVERS
#include <hgparallel.H>
#endif

#ifndef NDEBUG

extern "C"
{
    void PrintBoxArray (const BoxArray& ba);
    void PrintTagBox (const TagBox& tb);
    void PrintTagBoxArray (const TagBoxArray& tba);
    void TagBoxCount (const TagBox& tb);
    void TagBoxArrayCount (const TagBoxArray& tba);
}

void PrintBoxArray (const BoxArray& ba) { std::cout << ba << std::endl; }

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
#ifndef MG_USE_F90_SOLVERS
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

#ifndef MG_USE_F90_SOLVERS
    //
    // Close down the Holy Grail junk.
    //
    HG::MPI_finish();
#endif

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_stop = ParallelDescriptor::second() - run_strt;

    ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "Run time = " << run_stop << std::endl;

    BoxLib::Finalize();

    return 0;
}
