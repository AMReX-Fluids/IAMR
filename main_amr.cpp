//
// $Id: main_amr.cpp,v 1.18 2001-08-01 21:51:01 lijewski Exp $
//

#include <cstdlib>

#include <REAL.H>
#include <Utility.H>
#include <IntVect.H>
#include <Box.H>
#include <Amr.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <AmrLevel.H>

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);

    int  max_step;
    Real stop_time;

    ParmParse pp;

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

    delete amrptr;

    BoxLib::Finalize();

    return 0;
}
