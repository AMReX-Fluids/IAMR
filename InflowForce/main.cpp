//
// Merge swirl-type turbulence files together over specified range.
//
// Input file format:
//
//   N
//   time1 fabname1
//   time2 fabname2
//   time3 fabname3
//   time4 fabname4
//   ...
//   ...
//
// For a total of N lines.
//

#include <list>
#include <limits>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include <REAL.H>
#include <BLassert.H>
#include <BoxLib.H>
#include <FArrayBox.H>
#include <ParmParse.H>
#include <Utility.H>

struct Line
{
    Line (Real               time,
          const std::string& name)
        :
        m_time(time),
        m_name(name)
        {}

    bool operator<  (const Line& rhs) const { return m_time < rhs.m_time;  }

    bool operator== (const Line& rhs) const { return m_time == rhs.m_time; }

    Real        m_time;
    std::string m_name;
};

static
void
usage ()
{
    std::cout << "\nusage: turbMerge"
              << " ifile=inputfile"
              << " ofile=outputfile"
              << " DX=\"dx dy dz\"]\n"
              << " [range=\"tstart tend\"]\n"
              << std::endl;
    exit(1);
}

static
void
GetFab (const std::string& name, FArrayBox& fab)
{
    std::ifstream ifs(name.c_str(), std::ios::in);

    std::cout << "Attempting to read FAB from file: " << name << " ... " << std::flush;

    fab.readFrom(ifs);

    std::cout << "done\n\n" << std::flush;

    if (fab.box().length(2) != 1)
    {
        std::cout << "Z dimension of FAB box is not one!\n";
        usage();
    }

    if (fab.nComp() < 3)
    {
        std::cout << "FAB has fewer than three components!\n";
        usage();
    }
}

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);

    ParmParse pp;

    std::string       ofile;
    std::vector<Real> range;
    std::vector<Real> DX;
    std::string       ifile;

    pp.queryarr("range", range);

    if (range.empty())
    {
        //
        // Use full range of times in the input file.
        //
        range.push_back(-1);
        range.push_back( std::numeric_limits<Real>::max());
    }
    else
    {
        if (range.size() != 2) usage();
    }

    if (range[0] >= range[1])
    {
        std::cout << "range[0] MUST be less than range[1]!!!\n";
        usage();
    }

    std::cout << "\nrange: " << range[0] << ' ' << range[1] << '\n';

    pp.queryarr("DX", DX);

    if (DX.size() != 3)
    {
        std::cout << "DX must have three components!\n";
        usage();
    }

    std::cout << "\nDX: " << DX[0] << ' ' << DX[1] << ' ' << DX[2] << '\n';

    pp.query("ifile",ifile);

    if (ifile.empty()) usage();

    pp.query("ofile",ofile);

    if (ofile.empty()) usage();

    std::cout << "\nifile: " << ifile << '\n';
    std::cout << "\nofile: " << ofile << '\n';

    std::ifstream ifs(ifile.c_str(), std::ios::in);

    int N;

    ifs >> N;

    std::cout << "\nN: " << N << std::endl;

    if (N  < 2)
    {
        std::cout << "ifile needs to contain at least two lines!\n";
        usage();
    }

    std::list<Line> LL;

    Real        tm;
    std::string nm;

    for (int i = 0; i < N; i++)
    {
        ifs >> tm;
        ifs >> nm;

        if (tm < 0)
        {
            std::cout << "Got a negative time!\n";
            usage();
        }

        if (nm.empty())
        {
            std::cout << "Got an empty FAB file name!\n";
            usage();
        }

        LL.push_back(Line(tm,nm));
    }

    LL.sort();

    LL.unique();

    if (LL.empty())
    {
        std::cout << "No data falls within the specified range!\n";
        usage();
    }

    std::cout << '\n' << ofile
              << " will have "
              << LL.size()
              << " platters of data!\n\n";

    std::cout << "Time       File\n";
    std::cout << "---------------\n";

    for (std::list<Line>::const_iterator it = LL.begin(); it != LL.end(); ++it)
    {
        std::cout << it->m_time << '\t' << it->m_name << '\n';
    }
    std::cout << '\n';

    if (!BoxLib::UtilCreateDirectory(ofile, 0755))
        BoxLib::CreateDirectoryFailed(ofile);

    std::string Hdr = ofile; Hdr += "/HDR";
    std::string Dat = ofile; Dat += "/DAT";

    std::ofstream ohdr, odat;

    ohdr.open(Hdr.c_str(), std::ios::out|std::ios::trunc);
    if (!ohdr.good())
        BoxLib::FileOpenFailed(Hdr);

    odat.open(Dat.c_str(), std::ios::out|std::ios::trunc);
    if (!odat.good())
        BoxLib::FileOpenFailed(Dat);
    //
    // Gotta open a FAB to get NX and NY.
    //
    FArrayBox fab;

    GetFab(LL.front().m_name, fab);
    
    const int NX = fab.box().length(0);
    const int NY = fab.box().length(1);
    const int NZ = LL.size();

    ohdr << NX << ' '
         << NY << ' '
         << NZ << '\n';

    ohdr << DX[0] << ' '
         << DX[1] << ' '
         << DX[2] << '\n';

    ohdr << 1 << ' ' << 1 << ' ' << 1 << '\n';

    for (std::list<Line>::const_iterator it = LL.begin(); it != LL.end(); ++it)
    {
        ohdr << it->m_time << '\n';
    }
    //
    // Write out the first FABs info.
    //
    for (int i = 0; i < 3; i++)
    {
        ohdr << odat.tellp() << '\n';
        //
        // Write the FAB to the data file.
        //
        fab.writeOn(odat,i,1);
    }
    //
    // Now get the rest of the FABs;
    //
    std::list<Line>::const_iterator it = LL.begin();

    ++it; // Skip past first FAB.

    for ( ; it != LL.end(); ++it)
    {
        GetFab(it->m_name, fab);

        if (fab.box().length(0) != NX)
        {
            std::cout << "FAB does not have correct X dimension!\n";
            usage();
        }
        if (fab.box().length(1) != NY)
        {
            std::cout << "FAB does not have correct Y dimension!\n";
            usage();
        }
        if (fab.box().length(2) != 1)
        {
            std::cout << "FAB does not have correct Z dimension!\n";
            usage();
        }

        for (int i = 0; i < 3; i++)
        {
            ohdr << odat.tellp() << '\n';
            //
            // Write the FAB to the data file.
            //
            fab.writeOn(odat,i,1);
        }
    }

    BoxLib::Finalize();

    return 0;
}

