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

#if BL_SPACEDIM != 3
#error "This code only works for BL_SPACEDIM==3"
#endif

#include <list>
#include <limits>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include <AMReX_REAL.H>
#include <AMReX_BLassert.H>
#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

using namespace amrex;

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

static bool verbose = false;

static
void
usage ()
{
    std::cout << "\nusage: turbMerge"
              << " ifile=inputfile"
              << " ofile=outputfile"
              << " probsize=\"dx dy dz\"\n"
              << " [range=\"tstart tend\"]\n"
              << " [verbose=[0|1]]\n"
              << std::endl;
    exit(1);
}

static
void
GetFab (const std::string& name, FArrayBox& fab)
{
    std::ifstream ifs(name.c_str(), std::ios::in);

    if (verbose)
        std::cout << "Attempting to read FAB from file: " << name << " ... " << std::flush;

    fab.readFrom(ifs);

    if (verbose)
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

static
void
Extend (FArrayBox& xfab,
        FArrayBox& vfab,
        const Box& dm)
{
    Box tbx = vfab.box();

    tbx.setBig(0, dm.bigEnd(0) + 3);

    const int ygrow = BL_SPACEDIM==3 ? 3 : 1;

    tbx.setBig(1, dm.bigEnd(1) + ygrow);

    xfab.resize(tbx,vfab.nComp());

    xfab.copy(vfab);
    vfab.shift(0, dm.length(0));
    xfab.copy(vfab);
    vfab.shift(1, dm.length(1));
    xfab.copy(vfab);
    vfab.shift(0, -dm.length(0));
    xfab.copy(vfab);
}

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    ParmParse pp;

    std::string       ofile;
    std::vector<Real> range;
    std::vector<Real> probsize;
    std::string       ifile;

    pp.query("verbose",verbose);

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

    pp.queryarr("probsize", probsize);

    if (probsize.size() != 3)
    {
        std::cout << "probsize must have three components!\n";
        usage();
    }

    std::cout << "\nprobsize: " << probsize[0] << ' ' << probsize[1] << ' ' << probsize[2] << '\n';

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

    if (!amrex::UtilCreateDirectory(ofile, 0755))
        amrex::CreateDirectoryFailed(ofile);

    std::string Hdr = ofile; Hdr += "/HDR";
    std::string Dat = ofile; Dat += "/DAT";

    std::ofstream ohdr, odat;

    ohdr.open(Hdr.c_str(), std::ios::out|std::ios::trunc);
    if (!ohdr.good())
        amrex::FileOpenFailed(Hdr);

    odat.open(Dat.c_str(), std::ios::out|std::ios::trunc);
    if (!odat.good())
        amrex::FileOpenFailed(Dat);
    //
    // Gotta open a FAB to get initial NX and NY.
    //
    FArrayBox fab, xfab;

    GetFab(LL.front().m_name, fab);
    
    const int NX = fab.box().length(0);
    const int NY = fab.box().length(1);
    const int NZ = LL.size();

    const Box domain(IntVect(0,0,0), IntVect(NX-1,NY-1,0));

    const Real DX[BL_SPACEDIM] = {probsize[0]/NX, probsize[1]/NY, 1};

    ohdr << NX + 3 << ' '
         << NY + 3 << ' '
         << NZ     << '\n';

    ohdr << probsize[0] + 2*DX[0] << ' '
         << probsize[1] + 2*DX[1] << ' '
         << probsize[2]           << '\n';

    ohdr << 1 << ' ' << 1 << ' ' << 1 << '\n';

    for (std::list<Line>::const_iterator it = LL.begin(); it != LL.end(); ++it)
    {
        ohdr << it->m_time << '\n';
    }
    //
    // Write out the FABS one slab at a time.
    //
    // All X's, then all the Y,s, then all the Z's.
    //
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        int len;

        switch(i)
        {
        case 0: len = NX; break;
        case 1: len = NY; break;
        case 2: len = 1;  break;
        default:
            amrex::Abort("turbMerge: how did this happen?");
        }

        for (std::list<Line>::const_iterator it = LL.begin(); it != LL.end(); ++it)
        {
            GetFab(it->m_name, fab);

            if (fab.box().length(i) != len)
            {
                std::cout << "FAB does not have correct length in dim = " << i << "!\n";
                usage();
            }

            Extend(xfab, fab, domain);

            ohdr << odat.tellp() << '\n';
            //
            // Write the extended FAB to the data file.
            //
            xfab.writeOn(odat,i,1);
        }
    }

    amrex::Finalize();

    return 0;
}

