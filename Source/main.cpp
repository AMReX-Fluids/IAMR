
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
#include <FabSet.H>
#include <Profiler.H>

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

    BL_PROFILE_VAR("main()", pmain);

    const Real run_strt = ParallelDescriptor::second();

    int  max_step;
    Real strt_time;
    Real stop_time;

    ParmParse pp;

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

    //
    // If we set the regrid_on_restart flag and if we are *not* going to take
    // a time step then we want to go ahead and regrid here.
    //
    if (amrptr->RegridOnRestart())
    {
        if (    (amrptr->levelSteps(0) >= max_step ) ||
                ( (stop_time >= 0.0) &&
                  (amrptr->cumTime() >= stop_time)  )    )
        {
            //
            // Regrid only!
            //
            amrptr->RegridOnly(amrptr->cumTime());
        }
    }

    while ( amrptr->okToContinue()           &&
           (amrptr->levelSteps(0) < max_step || max_step < 0) &&
           (amrptr->cumTime() < stop_time || stop_time < 0.0) )
    {
        amrptr->coarseTimeStep(stop_time);
    }

    // Write final checkpoint and plotfile
    if (amrptr->stepOfLastCheckPoint() < amrptr->levelSteps(0)) {
        amrptr->checkPoint();
    }

    if (amrptr->stepOfLastPlotFile() < amrptr->levelSteps(0)) {
        amrptr->writePlotFile();
    }

    delete amrptr;

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_stop = ParallelDescriptor::second() - run_strt;

    ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "Run time = " << run_stop << std::endl;

    BL_PROFILE_VAR_STOP(pmain);
    BL_PROFILE_SET_RUN_TIME(run_stop);
    BL_PROFILE_FINALIZE();


    BoxLib::Finalize();

    return 0;
}
