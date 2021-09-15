
#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using std::ios;

#include <unistd.h>

#include <WritePlotFile.H>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>
#include <AMReX_ArrayLim.H>

#include <EXACT_F.H>

#define GARBAGE 666.e+40

using namespace amrex;

static
void
PrintUsage (const char* progName)
{
    std::cout << '\n';
    std::cout << "This program reads a plot file and compares the results\n"
         << "against a 3-d analog to the 2-d viscous benchmark of\n"
         << "G.I. Taylor.  Essentially, this is a 2-d solution that\n"
         << "is uniform in the third-dimension.  But, the uniform\n"
         << "direction can be in any of the three coordinate directions\n";
    std::cout << "Usage:" << '\n';
    std::cout << progName << '\n';
    std::cout << "     infile = inputFileName" << '\n';
    std::cout << "    unifdir = uniformDirection" << '\n';
    std::cout << "    errfile = ErrorFileOutputFileName" << '\n';
    std::cout << "     exfile = ExactSolnOutputFileName" << '\n';
    std::cout << "         mu = viscosity" << '\n';
    std::cout << "       norm = integer norm (Ie. default is 2 for Ln norm)" << '\n';
    std::cout << "   [-help]" << '\n';
    std::cout << "   [-verbose]" << '\n';
    std::cout << '\n';
    exit(1);
}

bool
amrDatasHaveSameDerives(const AmrData& amrd1,
			const AmrData& amrd2);
int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    if (argc == 1)
        PrintUsage(argv[0]);

    ParmParse pp;

    if (pp.contains("help"))
        PrintUsage(argv[0]);

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    //
    // Scan the arguments.
    //
    std::string iFile, errFile, exFile;

    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
        AmrData::SetVerbose(true);
    }
    pp.query("infile", iFile);
    if (iFile.empty())
        amrex::Abort("You must specify `infile'");

    pp.query("exfile", exFile);
    pp.query("errfile", errFile);

    int unifDir = -1;
    pp.query("unifdir", unifDir);
    if (unifDir < 0)
        amrex::Abort("You must specify `unifdir'");

    Real mu = -1.0;
    pp.query("mu", mu);
    if (mu < 0.0) 
        amrex::Abort("You must specify `mu'");

    int norm = 2;
    pp.query("norm", norm);

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    
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
    const Vector<std::string>& derives = amrDataI.PlotVarNames();
    

    //
    // Compute the error
    //
    Vector<MultiFab*> error(finestLevel+1);
    Vector<MultiFab*> dataE(finestLevel+1);
    
    std::cout << "Level Delta L"<< norm << " norm of Error in Each Component" << std::endl
         << "-----------------------------------------------" << std::endl;

    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        const BoxArray& baI = amrDataI.boxArray(iLevel);
        Vector<Real> delI = amrDataI.DxLevel()[iLevel];
	DistributionMapping dmap(baI);

	error[iLevel] = new MultiFab(baI, dmap, nComp, 0);
	error[iLevel]->setVal(GARBAGE);

	dataE[iLevel] = new MultiFab(baI, dmap, nComp, 0);
	dataE[iLevel]->setVal(GARBAGE);

        MultiFab dataI(baI, dmap, nComp, 0);

	for (int iGrid=0; iGrid<baI.size(); ++iGrid)
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
            Vector<Real> xlo = amrDataI.GridLocLo()[iLevel][iGrid];
            Vector<Real> xhi = amrDataI.GridLocHi()[iLevel][iGrid];

            FORT_VISCBENCH(&time, &mu, &unifDir,
                           lo, hi, &nComp,
                           ((*dataE[iLevel])[iGrid]).dataPtr(), 
                           ARLIM(lo), ARLIM(hi),
                           delI.dataPtr(),
                           xlo.dataPtr(), xhi.dataPtr());
	}

        (*error[iLevel]).ParallelCopy(dataI);
        (*error[iLevel]).minus((*dataE[iLevel]), 0, nComp, 0);

   
        //
        // Output Statistics
        //
        Real cellvol = 1.0;
        for (int i=0; i<BL_SPACEDIM; i++)
           cellvol = cellvol * delI[i];

        Real delAvg = pow(cellvol, (1.0/BL_SPACEDIM));

        std::cout << "  " << iLevel << " " << delAvg << "    ";
        for (int iComp=0; iComp<nComp; ++iComp)
        {
            Real Ln = 0.0;
            for (int iGrid=0; iGrid<baI.size(); ++iGrid)
            {
                Real grdLn = (*error[iLevel])[iGrid].norm(norm,iComp,1);
                Ln = Ln + pow(grdLn, norm) * cellvol;
            }
            Ln = pow(Ln, (1.0/norm));

            std::cout << Ln << "  ";
        }
        std::cout << std::endl;
    }

    //
    // Write Plot Files
    //
    if (!errFile.empty())
        WritePlotFile(error, amrDataI, errFile, verbose);

    if (!exFile.empty())
        WritePlotFile(dataE, amrDataI, exFile, verbose);


    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
	delete error[iLevel];

    amrex::Finalize();
}


bool
amrDatasHaveSameDerives(const AmrData& amrd1,
			const AmrData& amrd2)
{
    const Vector<std::string>& derives1 = amrd1.PlotVarNames();
    const Vector<std::string>& derives2 = amrd2.PlotVarNames();
    int length = derives1.size();
    if (length != derives2.size())
	return false;
    for (int i=0; i<length; ++i)
	if (derives1[i] != derives2[i])
	    return false;
    return true;
}
    
