//
// Merge swirl-type turbulence files together over specified range.
//

#include <iostream>

#include <list>
#include <REAL.H>
#include <vector>
#include <BLassert.H>
#include <BoxLib.H>
#include <ParmParse.H>

#include "main_F.H"

static
void
usage ()
{
    std::cout << "usage: turbMerge"
              << " range=\"tstart tend\""
              << " ifiles=\"ifile1 ... ifileN\""
              << " ofile=ofname"
              << std::endl;
    exit(1);
}

static
std::vector<int>
EncodeStringForFortran(const std::string& astr)
{
    int length = astr.size();
    std::vector<int> result(length);
    for (int i = 0; i < length; ++i)
        result[i] = astr[i];
    return result;
}

struct Triplet
{
    Triplet ()
        :
        first(0), second(0), third(0) {}

    Triplet (Real r ,int i, int j)
        :
        first(r), second(i), third(j) {}

    bool operator< (const Triplet& rhs) const { return first < rhs.first; }

    bool operator== (const Triplet& rhs) const { return first == rhs.first; }

    Real first;  // Time corresponding to the platter of data.
    int  second; // Index of the file containing that platter of data.
    int  third;  // Fortran index in that file of the platter of data.
};

typedef Triplet Platter;

typedef std::list<Platter> PList;

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);

    ParmParse pp;

    std::string              ofile;
    std::vector<Real>        range;
    std::vector<std::string> ifiles;

    pp.queryarr("range",range);

    if (range.size() != 2) usage();

    if (range[0] >= range[1])
    {
        std::cout << "range[0] MUST be less than range[1]!!!" << std::endl;
        usage();
    }

    switch(pp.countval("ifiles"))
    {
    case 0:  usage();                                      break;
    case 1:  ifiles.resize(1); pp.get("ifiles",ifiles[0]); break;
    default: pp.queryarr("ifiles",ifiles);                 break;
    }

    pp.query("ofile",ofile);

    if (ofile.empty()) usage();

    const int NN = 1024;
    int  imax, jmax, kmax, ncomp;
    Real scalex, scaley, scalez;
    std::vector<Real> times(NN);
    std::vector< std::vector<Real> > ftimes(ifiles.size());

    std::vector<int> fname = EncodeStringForFortran(ifiles[0]);
    int              flen  = fname.size();

    std::cout << "Attempting to read header from: " << ifiles[0] << std::endl;
    FORT_READ_SWIRL_HEADER(&fname[0],&flen,&imax,&jmax,&kmax,&ncomp,
                           &scalex,&scaley,&scalez,&times[0],&NN);

    ftimes[0].resize(kmax);
    for (int i = 0; i < kmax; i++)
        ftimes[0][i] = times[i];

    std::cout << "imax, jmax, ncomp: "
              << imax << ", "
              << jmax << ", "
              << ncomp << std::endl;

    std::cout << "scalex, scaley, scalez: "
              << scalex << ", "
              << scaley << ", "
              << scalez << std::endl;

    for (int i = 1; i < ifiles.size(); i++)
    {
        const int NN = 1024;
        int  imaxt, jmaxt, kmaxt, ncompt;
        Real scalext, scaleyt, scalezt;

        std::vector<int> fname = EncodeStringForFortran(ifiles[i]);
        int              flen  = fname.size();

        std::cout << "Attempting to read header from: " << ifiles[i] << std::endl;

        FORT_READ_SWIRL_HEADER(&fname[0],&flen,&imaxt,&jmaxt,&kmaxt,&ncompt,
                               &scalext,&scaleyt,&scalezt,&times[0],&NN);

        ftimes[i].resize(kmaxt);
        for (int j = 0; j < kmaxt; j++)
            ftimes[i][j] = times[j];

        if (imax != imaxt) {
            std::cout << "imax != imaxt" << std::endl; exit(1);
        }
        if (jmax != jmaxt) {
            std::cout << "jmax != jmaxt" << std::endl; exit(1);
        }
        if (scalex != scalext) {
            std::cout << "scalex != scalext" << std::endl; exit(1);
        }
        if (scaley != scaleyt) {
            std::cout << "scaley != scaleyt" << std::endl; exit(1);
        }
    }

    PList LP;

    for (int i = 0; i < ifiles.size(); i++)
        for (int j = 0; j < ftimes[i].size(); j++)
            if (ftimes[i][j] >= range[0] && ftimes[i][j] <= range[1])
                LP.push_back(Platter(ftimes[i][j],i,j+1));

    LP.sort();

    LP.unique();

    if (LP.empty())
    {
        std::cout << "No data falls within the specified range!" << std::endl;
        exit(1);
    }

    std::cout << ofile
              << " will have "
              << LP.size()
              << " platters of data!" << std::endl;

    std::cout << "Time       File      Index" << std::endl;
    std::cout << "--------------------------" << std::endl;

    for (PList::const_iterator it = LP.begin(); it != LP.end(); ++it)
    {
        std::cout << it->first
                  << "\t"
                  << ifiles[it->second]
                  << "\t"
                  << it->third
                  << std::endl;
    }

    kmax = LP.size();

    Real* fluct = new Real[kmax];

    Real* block = new Real[imax*jmax*kmax*ncomp];

    int putidx = 1;

    for (PList::const_iterator it = LP.begin(); it != LP.end(); )
    {
        std::vector<int> fname  = EncodeStringForFortran(ifiles[it->second]);
        int              flen   = fname.size();
        int              getidx = it->third;
        int              count  = 1;  // How many consecutive platters?

        fluct[putidx-1] = it->first;

        PList::const_iterator cit = it; ++cit;

        for ( ; cit != LP.end(); ++cit)
        {
            if (ifiles[it->second] == ifiles[cit->second])
            {
                BL_ASSERT((it->third + count) == cit->third);

                fluct[putidx+count-1] = cit->first;

                count++;
            }
            else
                break;
        }

        std::cout << "Reading "
                  << count
                  << " platters starting at "
                  << getidx
                  << " from "
                  << ifiles[it->second] << std::endl;

        FORT_FILL(block,
                  &imax,
                  &jmax,
                  &kmax,
                  &ncomp,
                  &fluct[putidx-1],
                  &putidx,
                  &getidx,
                  &fname[0],
                  &flen,
                  &count);

        putidx += count;
        for (int i = 0; i < count; i++) ++it;
    }

    std::cout << "Building swirl-type turbulence file: " << ofile << std::endl;    

    std::vector<int> oname  = EncodeStringForFortran(ofile);
    int              olen   = oname.size();

    FORT_WRITE(block,
               &imax,
               &jmax,
               &kmax,
               &ncomp,
               fluct,
               &scalex,
               &scaley,
               &scalez,
               &oname[0],
               &olen);

    BoxLib::Finalize();

    return 0;
}

