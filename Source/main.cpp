
#include <cstdio>

#include <AMReX_CArena.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfiler.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_AmrLevel.H>
#endif

using namespace amrex;

#ifdef AMREX_USE_EB
//skipping header file and just declaring eb2 init fn here as in CNS for now
void initialize_EB2 (const Geometry& geom, const int required_level,
		     const int max_level);
#endif

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    BL_PROFILE_REGION_START("main()");
    BL_PROFILE_VAR("main()", pmain);

    const Real run_strt = ParallelDescriptor::second();

    int  max_step;
    int  num_steps;
    Real strt_time;
    Real stop_time;
    Real stop_interval;

    ParmParse pp;

    max_step  = -1; 
    num_steps = -1; 
    strt_time =  0.0;
    stop_time = -1.0;
    stop_interval = 0.;

    pp.query("max_step",  max_step);
    pp.query("num_steps", num_steps);
    pp.query("strt_time", strt_time);
    pp.query("stop_time", stop_time);
    pp.query("stop_interval", stop_interval);

    if (strt_time < 0.0)
    {
        amrex::Abort("MUST SPECIFY a non-negative strt_time");
    }

    if (max_step < 0 && stop_time < 0)
    {
        amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    Amr* amrptr = new Amr;
    //    Amr amr;
#ifdef AMREX_USE_EB
    // fixme? not sure what level of support should be default
    // levels explianed in user guide
    // Ann suggested we need vol and area frac, and face and area centriod => full
    AmrLevel::SetEBSupportLevel(EBSupport::full);
    // set grow cells for basic, volume, full
    // fixme? not sure what these numbers should be
    // AmrLevel.cpp defaults 5, 4, 2
    // CNS::numGrow()= 5
    //AmrLevel::SetEBMaxGrowCells(CNS::numGrow(),4,2);
    // NavierStokesBase GEOM_GROW=1 currently. Change it? Make new var?
    // Using incflo values here
    AmrLevel::SetEBMaxGrowCells(4,4,4);

    //decide who should own max_coasening_level later
    int max_coarsening_level = 100;
    pp.query("max_coarsening_level", max_coarsening_level);
    initialize_EB2(amrptr->Geom(amrptr->maxLevel()), amrptr->maxLevel(),
		   max_coarsening_level);
#endif
		   
    amrptr->init(strt_time,stop_time);

    // This feature stop the simulation at a specfic time 
    // after the physical time of the checkpoint file 
    if (stop_interval > 0.) stop_time = amrptr->cumTime() + stop_interval;

    if (num_steps > 0)
    {
        if (max_step < 0)
        {
            max_step = num_steps + amrptr->levelSteps(0);
        }
        else
        {
            max_step = std::min(max_step, num_steps + amrptr->levelSteps(0));
        }

	amrex::Print() << "Using effective max_step = " << max_step << '\n';
    }
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

    while ( amrptr->okToContinue()                            &&
           (amrptr->levelSteps(0) < max_step || max_step < 0) &&
           (amrptr->cumTime() < stop_time || stop_time < 0.0) )
    {
        amrptr->coarseTimeStep(stop_time);
    }
    //
    // Write final checkpoint and plotfile.
    //
    if (amrptr->stepOfLastCheckPoint() < amrptr->levelSteps(0))
    {
        amrptr->checkPoint();
    }

    if (amrptr->stepOfLastPlotFile() < amrptr->levelSteps(0))
    {
        amrptr->writePlotFile();
    }

    delete amrptr;

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_stop = ParallelDescriptor::second() - run_strt;

    ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

    amrex::Print() << "Run time = " << run_stop << std::endl;

    BL_PROFILE_VAR_STOP(pmain);
    BL_PROFILE_REGION_STOP("main()");
    BL_PROFILE_SET_RUN_TIME(run_stop);
    BL_PROFILE_FINALIZE();


    amrex::Finalize();

    return 0;
}
