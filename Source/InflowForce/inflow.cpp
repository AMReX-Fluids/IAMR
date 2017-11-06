
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>

#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;

#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_GETPLANE    GETPLANE
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_GETPLANE    getplane
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_GETPLANE    getplane_
#  endif

extern "C" void FORT_GETPLANE(int* filename, int* len, Real* data, int* plane, int* ncomp, int* isswirltype);

void
FORT_GETPLANE (int* filename, int* len, Real* data, int* plane, int* ncomp, int* isswirltype)
{
    static int         kmax;
    static bool        first = true;
    static Vector<long> offset;
    std::string        flctfile;

    for (int i = 0; i < *len; i++)
    {
        char c = filename[i];

        flctfile += c;
    }

    if (first)
    {
        //
        // Read and save all the seekp() offsets in the inflow header file.
        //
        first = false;

        std::string hdr = flctfile; hdr += "/HDR";

        std::ifstream ifs;

        ifs.open(hdr.c_str(), std::ios::in);

        if (!ifs.good())
            amrex::FileOpenFailed(hdr);

        int  idummy;
        Real rdummy;
        //
        // Hardwire loop max to 3 regardless of spacedim.
        //
        for (int i = 0; i < 3; i++)
            ifs >> kmax;

        ifs >> rdummy >> rdummy >> rdummy;
        ifs >> idummy >> idummy >> idummy;

        if (*isswirltype)
        {
            //
            // Skip over fluct_times array.
            //
            for (int i = 0; i < kmax; i++)
                ifs >> rdummy;
        }

        offset.resize(kmax * BL_SPACEDIM);

        for (int i = 0; i < offset.size(); i++)
            ifs >> offset[i];
    }

    std::string dat = flctfile; dat += "/DAT";

    std::ifstream ifs;

    ifs.open(dat.c_str(), std::ios::in);

    if (!ifs.good())
        amrex::FileOpenFailed(dat);
    //
    // There are BL_SPACEDIM * kmax planes of FABs.
    // The first component are in the first kmax planes,
    // the second component in the next kmax planes, ....
    // Note also that both (*plane) and (*ncomp) start from
    // 1 not 0 since they're passed from Fortran.
    //
    const long start = offset[((*plane) - 1) + (((*ncomp) - 1) * kmax)] ;

    ifs.seekg(start, std::ios::beg);

    if (!ifs.good())
        amrex::Abort("getplane(): seekg() failed");

    FArrayBox fab;

    fab.readFrom(ifs);

    memcpy(data, fab.dataPtr(), fab.box().numPts()*sizeof(Real));
}
