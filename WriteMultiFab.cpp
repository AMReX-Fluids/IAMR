#include <WriteMultiFab.H>
#include <Tuple.H>
#include <MultiFab.H>
#include <aString.H>
#include <ParallelDescriptor.H>
#ifdef BL_USE_BSP
#include <bsp.h>
#endif

#ifdef BL_USE_NEW_HFILES
#include <ostream>
#include <fstream>
using std::ofstream;
#else
#include <iostream.h>
#include <fstream.h>
#endif

/*
struct TagStruct {
  Box box;
  int level;
  int steps;
  Real time;
  Real probLo[BL_SPACEDIM];
  Real probHi[BL_SPACEDIM];
  int startVar;
  int nVar;
};
*/


void
WriteMultiFab(char *filename, const MultiFab &multifab,
	      Real _H, const BoxList &bs, const Box &container,
	      int ratio, Real bgValue)
{
    Real H[BL_SPACEDIM];
    int n;
    for(n = 0; n < BL_SPACEDIM; ++n) {
	H[n] = _H;
    }
    WriteMultiFab(filename, multifab, H, bs, container, ratio, bgValue);
}


void
WriteMultiFab(char *filename, const MultiFab &multifab,
	      const Real *H, const BoxList &bs, const Box &container,
	      int ratio, Real bgValue)
{
    int i, n, myproc = ParallelDescriptor::MyProc();
    Tuple< Real, BL_SPACEDIM > probLo;
    Tuple< Real, BL_SPACEDIM > probHi;
    const int *lo = container.loVect();
    const int *hi = container.hiVect();
    Real HBG[BL_SPACEDIM];
    for(n = 0; n < BL_SPACEDIM; n++) {
	probLo[n] = H[n] * lo[n];
	probHi[n] = H[n] * (hi[n]+1.0);
	HBG[n] = H[n] * ratio;
    }
    Box bxbg(container);
    bxbg.coarsen(ratio);
    FArrayBox fab(bxbg, multifab.nComp());

    fab.setVal(bgValue);

    ofstream os;
    if(ParallelDescriptor::IOProcessor()) {

      os.open(filename);
      os << multifab.nComp() << "\n"; // Number of states dumped
      for(int nv = 0; nv < multifab.nComp(); nv++) {
        os << "var" << nv << "\n"; // Name of state(s)
      }
      os << BL_SPACEDIM     << "\n"; // Dimension of data
      os << "0.0"        << "\n"; // Simulation time of dump
      os << "1"          << "\n"; // Finest level dumped
    
      for(n = 0; n < BL_SPACEDIM; n++)  os << probLo[n]  << " ";
      os << "\n";                                  // Position of lo-sides
      for(n = 0; n < BL_SPACEDIM; n++)  os << probHi[n]  << " "; 
      os << "\n";                                  // Position of hi-sides
      os << ratio << "\n";                         // 0:f_lev-2 ratio
      os << bxbg << " " << container << "\n";      // 0:f_lev-1 container
      os << "0 0\n";                               // 0:f_lev-1 nSteps
      for(n = 0; n < BL_SPACEDIM; n++) os << HBG[n] << " ";
      os << "\n";                                  // Grid spacing on background
      for(n = 0; n < BL_SPACEDIM; n++) os << H[n]   << " "; 
      os << "\n";                                  // Grid spacing of data
      os << "0\n";                                 // Coord sys flag (0=cart)
      os << "0\n";                                 // BC flag (0=no BC info)

      // dump base grid
      os << "0 1\n";                      // level and number of grids
      os << bxbg << "\n";                 // For each grid, dump box
      os << "0\n";                        //                level
      os << "0\n";                        //                steps
      os << "0.0\n";                      //                time
      for(n = 0; n < BL_SPACEDIM; n++) { 
          os << probLo[n] << " "          //                probLo
             << probHi[n] << "\n";        //                probHi
      }                                   
      fab.writeOn(os, 0, multifab.nComp());   //             fab dump valid data

      // dump actual data
      os << "1 " << bs.length() << "\n";        // level and number of grids

    }  // end if(ParallelDescriptor::IOProcessor())

    ParallelDescriptor::Synchronize();

    for(i = 0; i < multifab.length(); i++) {
	FArrayBox fab(multifab.boxArray()[i], multifab.nComp());
	multifab.copy(fab);

      if(ParallelDescriptor::IOProcessor()) {

        os << fab.box() << endl;    // For each grid, dump box
	os << "1\n";                    //                level
	os << "0\n";                    //                steps
	os << "0.0\n";                  //                time
	for(n = 0; n < BL_SPACEDIM; n++) {
	    os << probLo[n] << " "      //                probLo
               << probHi[n] << "\n";    //                probHi
	} // endfor: dimension

	//multifab.copy(fab);
	fab.writeOn(os);


      }  // end if(ParallelDescriptor::IOProcessor())
    }  // end for(i...)

    if(ParallelDescriptor::IOProcessor()) {
      os.close();
    }
    ParallelDescriptor::Synchronize();
}
