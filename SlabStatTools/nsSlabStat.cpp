#ifdef BL_ARCH_CRAY
#ifdef BL_USE_DOUBLE
#error "DOUBLE PRECISION NOT ALLOWED ON CRAY"
#endif
#endif

#ifndef        WIN32
#include <unistd.h>
#endif

#include "MultiFab.H"
#include "ParmParse.H"
#include "Utility.H"
#include "ParallelDescriptor.H"

#include <new>
using std::setprecision;
#include <iostream>
#ifndef WIN32
using std::set_new_handler;
#endif

#include "ArrayLim.H"

#include "SlabStatTools.H"
#include "StatTypes.H"
#include "ValidTypes.H"
#include "PackArrayArray.H"


static
void 
PrintUsage(const char* progName)
{
    if (ParallelDescriptor::IOProcessor())
    {

        cout << "Usage: " << endl;
        cout << progName << " [options]" << endl << endl;
        cout << "\t Options:" << endl;
        cout << "\t   statType = Type of statistics to calculate" << endl;
        cout << "\t   begItime = Beginning integer time" << endl;
        cout << "\t   endItime = Ending integer time" << endl;
        cout << "\t   axialDir = Streamwise direction.  (1->x, 2->y)" << endl;
        cout << endl;
        cout << "\t Available Statistic Types:" << endl;
        for (int nclc = 0; nclc < NumStatsTypes; ++nclc)
        {
            cout << "\t    " << StatsTypes[nclc].name << endl;
        }
        exit(0);
    }
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

    ParallelDescriptor::StartParallel(&argc, &argv);


    //-------------//
    // Print Usage //
    //-------------//

    if (argc == 1)
        PrintUsage(argv[0]);


//
//  Parse the command line
//
    ParmParse pp(argc-1,argv+1);

    if (pp.contains("help"))
        PrintUsage(argv[0]);

    aString statType;
    pp.query("statType", statType);
    if (statType.isNull() && ParallelDescriptor::IOProcessor())
        BoxLib::Abort("You must specify `statType'");

    int statIndx = ValidType(statType);


    int begItime = -1;
    pp.query("begItime", begItime);
    if (begItime == -1 && ParallelDescriptor::IOProcessor())
        BoxLib::Abort("You must specify `begItime'");

    int endItime = -1;
    pp.query("endItime", endItime);
    if (endItime == -1 && ParallelDescriptor::IOProcessor())
        BoxLib::Abort("You must specify `endItime'");

    int axialDir = -1;
    pp.query("axialDir", axialDir);
    if (axialDir == -1 && ParallelDescriptor::IOProcessor())
        BoxLib::Abort("You must specify `axialDir'");



    //---------------------//
    // Read SlabStats Data //
    //---------------------//

    const aString baseDir("slabstats");
    const aString slabStatName("basicStats");
    MultiFab ssData;
    int statLevel;
    Array<Real> dxLevel(BL_SPACEDIM), probLo(BL_SPACEDIM), probHi(BL_SPACEDIM);

    accumulateStats(baseDir, begItime, endItime, slabStatName, statLevel,
                    dxLevel, probLo, probHi, ssData);

    cout << statLevel << endl;
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        cout << dir << " " << probLo[dir] << " " 
                           << probHi[dir] << " " << dxLevel[dir] << endl;
    }


    //------------------------------//
    // Calculate and Write Profiles //
    //------------------------------//

    for (MultiFabIterator mfi(ssData); mfi.isValid(); ++mfi)
    {
        //
        // Calculate Statistics
        //
        const int nComp = mfi().nComp();

        const int *vblo = mfi.validbox().loVect();
        const int *vbhi = mfi.validbox().hiVect();
        const int *fblo = mfi.fabbox().loVect();
        const int *fbhi = mfi.fabbox().hiVect();

        const int nStats = StatsTypes[statIndx].nResStats;

        //
        // Note: This logic assumes that the profiles are always in the 
        //       direction having the maximum dimension on the box.  This also
        //       assumes that the constant and profile direction are always
        //       in the X or Y directions.  So, this essentially assumes that
        //       you never want statistics on constant Z slices or profiles in
        //       the Z-direction.  This works for flows which are either 
        //       periodic in the Z-direction or are round about the centerline
        //       in either the X and Z or the Y and Z directions.  The round
        //       flowfields works since the value of nStns is multiplied by
        //       the Z-dimension if the box is 3-d.
        //
        int cnstDir, profDir;
        cnstDir = (mfi.validbox().length(0) < mfi.validbox().length(1) ? 0 : 1);
        profDir = 1 - cnstDir;

        int cnstIndx = mfi.validbox().smallEnd(cnstDir);
        int nStns = mfi.validbox().length(profDir);
        if (BL_SPACEDIM > 2) nStns = nStns * mfi.validbox().length(2);

        Array<Real> physStn(nStns);
        Array< Array<Real> > stats(nStns);
        for (int n = 0; n < nStns; n++)
            stats[n].resize(nStats);

        Array<Real> passStats(nStns * nStats);
        passStats = PackArrayArray(nStns, nStats, stats);

        int nActualStations;

        StatsTypes[statIndx].func(mfi().dataPtr(), &nComp, 
                                    ARLIM(fblo), ARLIM(fbhi),
                                  &nStats, &nStns, 
                                  passStats.dataPtr(), physStn.dataPtr(),
                                  vblo, vbhi, dxLevel.dataPtr(),
                                  probLo.dataPtr(), probHi.dataPtr(),
                                  &axialDir, &nActualStations);

        stats = UnPackArrayArray(nStns, nStats, passStats);


        //
        // Write Stats
        //
        aString oFile = StatsTypes[statIndx].oFileBase;

        char* cnstDirChar = (cnstDir == 0 ? "i" : "j");
        oFile = oFile + cnstDirChar;
        oFile = Utility::Concatenate(oFile, cnstIndx);
        oFile = oFile + ".dat";

        ofstream ostr;
        ostr.open(oFile.c_str());

        for (int stn = 0; stn < nActualStations; stn++)
        {
            ostr << physStn[stn] << " ";
            for (int stat = 0; stat < nStats; stat++)
                ostr << stats[stn][stat] << " ";

            ostr << endl;
        }

        ostr.close();
    }
}
