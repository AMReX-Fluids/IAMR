
//
// $Id: ViscBench3d.cpp,v 1.6 2001-10-17 17:53:33 lijewski Exp $
//

#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
using std::ios;
using std::set_new_handler;

#include <unistd.h>

#include "WritePlotFile.H"
#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "VisMF.H"
#include "ArrayLim.H"

#include "EXACT_F.H"

#define GARBAGE 666.e+40

static
void
PrintUsage (const char* progName)
{
    cout << '\n';
    cout << "This program reads a plot file and compares the results\n"
         << "against a 3-d analog to the 2-d viscous benchmark of\n"
         << "G.I. Taylor.  Essentially, this is a 2-d solution that\n"
         << "is uniform in the third-dimension.  But, the uniform\n"
         << "direction can be in any of the three coordinate directions\n";
    cout << "Usage:" << '\n';
    cout << progName << '\n';
    cout << "     infile = inputFileName" << '\n';
    cout << "    unifdir = uniformDirection" << '\n';
    cout << "    errfile = ErrorFileOutputFileName" << '\n';
    cout << "     exfile = ExactSolnOutputFileName" << '\n';
    cout << "         mu = viscosity" << '\n';
    cout << "       norm = integer norm (Ie. default is 2 for Ln norm)" << '\n';
    cout << "   [-help]" << '\n';
    cout << "   [-verbose]" << '\n';
    cout << '\n';
    exit(1);
}

bool
amrDatasHaveSameDerives(const AmrData& amrd1,
			const AmrData& amrd2);
int
main (int   argc,
      char* argv[])
{
    if (argc == 1)
        PrintUsage(argv[0]);
    //
    // Make sure to catch new failures.
    //
    set_new_handler(Utility::OutOfMemory);

    ParallelDescriptor::StartParallel(&argc, &argv);

    ParmParse pp(argc-1,argv+1);

    if (pp.contains("help"))
        PrintUsage(argv[0]);

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    //
    // Scan the arguments.
    //
    aString iFile, errFile, exFile;

    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
        AmrData::SetVerbose(true);
    }
    pp.query("infile", iFile);
    if (iFile.isNull())
        BoxLib::Abort("You must specify `infile'");

    pp.query("exfile", exFile);
    pp.query("errfile", errFile);

    int unifDir = -1;
    pp.query("unifdir", unifDir);
    if (unifDir < 0)
        BoxLib::Abort("You must specify `unifdir'");

    Real mu = -1.0;
    pp.query("mu", mu);
    if (mu < 0.0) 
        BoxLib::Abort("You must specify `mu'");

    int norm = 2;
    pp.query("norm", norm);

    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);
    
    DataServices dataServicesC(iFile, fileType);

    if (!dataServicesC.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);

    //
    // Generate AmrData Objects 
    //
    AmrData& amrDataI = dataServicesC.AmrDataRef();

    const int nComp          = amrDataI.NComp();
    const Real time = amrDataI.Time();
    const int finestLevel = amrDataI.FinestLevel();
    const Array<aString>& derives = amrDataI.PlotVarNames();
    

    //
    // Compute the error
    //
    Array<MultiFab*> error(finestLevel+1);
    Array<MultiFab*> dataE(finestLevel+1);
    
    cout << "Level Delta L"<< norm << " norm of Error in Each Component" << endl
         << "-----------------------------------------------" << endl;

    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        const BoxArray& baI = amrDataI.boxArray(iLevel);
        Array<Real> delI = amrDataI.DxLevel()[iLevel];

	error[iLevel] = new MultiFab(baI, nComp, 0);
	error[iLevel]->setVal(GARBAGE);

	dataE[iLevel] = new MultiFab(baI, nComp, 0);
	dataE[iLevel]->setVal(GARBAGE);

        MultiFab dataI(baI, nComp, 0);

	for (int iGrid=0; iGrid<baI.length(); ++iGrid)
	{
            const Box& dataGrid = baI[iGrid];

	    for (int iComp=0; iComp<nComp; ++iComp)
	    {
                FArrayBox tmpFab(dataGrid,1);

                amrDataI.FillVar(&tmpFab, dataGrid,
                                 iLevel, derives[iComp], 0);
                dataI[iGrid].copy(tmpFab,0,iComp,1);
            }

            //
            // Fill exact solution
            //
            const int *lo = dataGrid.loVect();
            const int *hi = dataGrid.hiVect();
            Array<Real> xlo = amrDataI.GridLocLo()[iLevel][iGrid];
            Array<Real> xhi = amrDataI.GridLocHi()[iLevel][iGrid];

            FORT_VISCBENCH(&time, &mu, &unifDir,
                           lo, hi, &nComp,
                           ((*dataE[iLevel])[iGrid]).dataPtr(), 
                           ARLIM(lo), ARLIM(hi),
                           delI.dataPtr(),
                           xlo.dataPtr(), xhi.dataPtr());
	}

        (*error[iLevel]).copy(dataI);
        (*error[iLevel]).minus((*dataE[iLevel]), 0, nComp, 0);

   
        //
        // Output Statistics
        //
        Real cellvol = 1.0;
        for (int i=0; i<BL_SPACEDIM; i++)
           cellvol = cellvol * delI[i];

        Real delAvg = pow(cellvol, (1.0/BL_SPACEDIM));

        cout << "  " << iLevel << " " << delAvg << "    ";
        for (int iComp=0; iComp<nComp; ++iComp)
        {
            Real Ln = 0.0;
            for (int iGrid=0; iGrid<baI.length(); ++iGrid)
            {
                Real grdLn = (*error[iLevel])[iGrid].norm(norm,iComp,1);
                Ln = Ln + pow(grdLn, norm) * cellvol;
            }
            Ln = pow(Ln, (1.0/norm));

            cout << Ln << "  ";
        }
        cout << endl;
    }

    //
    // Write Plot Files
    //
    if (!errFile.isNull())
        WritePlotFile(error, amrDataI, errFile, verbose);

    if (!exFile.isNull())
        WritePlotFile(dataE, amrDataI, exFile, verbose);


    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
	delete error[iLevel];

    //
    // This calls ParallelDescriptor::EndParallel() and exit()
    //
    DataServices::Dispatch(DataServices::ExitRequest, NULL);
}


bool
amrDatasHaveSameDerives(const AmrData& amrd1,
			const AmrData& amrd2)
{
    const Array<aString>& derives1 = amrd1.PlotVarNames();
    const Array<aString>& derives2 = amrd2.PlotVarNames();
    int length = derives1.length();
    if (length != derives2.length())
	return false;
    for (int i=0; i<length; ++i)
	if (derives1[i] != derives2[i])
	    return false;
    return true;
}
    
