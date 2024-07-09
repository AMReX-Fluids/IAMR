
#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using std::ios;

#include <unistd.h>

#include <AMReX_WritePlotFile.H>
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

// -----------------------------------------------------------
// This case is an unsteady viscous benchmark for which the
// exact solution in 2D is
//     u(x,y,t) =   Sin(2 Pi x) Cos(2 Pi y) Exp(-2 (2Pi)^2 Nu t)
//     v(x,y,t) = - Cos(2 Pi x) Sin(2 Pi y) Exp(-2 (2Pi)^2 Nu t)
//     p(x,y,t) = - {Cos(4 Pi x) + Cos(4 Pi y)} Exp(-4 (2Pi)^2 Nu t) / 4
// This tool reads a plotfile and compares it against this exact solution.
// This benchmark was originally derived by G.I. Taylor (Phil. Mag.,
// Vol. 46, No. 274, pp. 671-674, 1923) and Ethier and Steinman
// (Intl. J. Num. Meth. Fluids, Vol. 19, pp. 369-375, 1994) give
// the pressure field.
//
// For 3D, the 2D solution is used in conjunction with a uniform 3rd dimension.
// When generating 3D data, check IAMR/Source/prob/prob_init.cpp to ensure
// that the problem setup matches the solution here with t=0.
//
// Analytic, fully 3D solutions (not implemented here) are discussed in,
// for example,
// Antuono, M. (2020). Tri-periodic fully three-dimensional analytic solutions for the Navier–Stokes equations. Journal of Fluid Mechanics, 890, A23. doi:10.1017/jfm.2020.126
//
static
void
PrintUsage (const char* progName)
{
    std::cout << '\n';
    std::cout << "This program reads a plot file and compares the results\n"
              << "against the 2-d viscous benchmark of G.I. Taylor (for 2d build)\n"
              << "or a 3-d analog, which takes the 2-d solution and then is\n"
              << "is uniform in the third-dimension.  The uniform\n"
              << "direction can be in any of the three coordinate directions\n";
    std::cout << "Usage:" << '\n';
    std::cout << progName << '\n';
    std::cout << "     infile = inputFileName" << '\n';
    std::cout << "         mu = viscosity" << '\n';
    std::cout << "    unifdir = uniformDirection (optional, default is z-dir)" << '\n';
    std::cout << "    errfile = ErrorFileOutputFileName (optional)" << '\n';
    std::cout << "     exfile = ExactSolnOutputFileName (optional)" << '\n';
    std::cout << "       norm = integer norm (optional, default is 2 for L2 norm, 1 for L1 norm, 0 for L_{inf} norm)" << '\n';
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

    // Direction with uniform velocity. Default to z-direction
    int unifDir = 2;
    pp.query("unifdir", unifDir);
    if (unifDir < 0 || unifDir > 2)
        amrex::Abort("You must specify a valid direction of uniform velocity, `unifdir = ...'");

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

    const int nComp          = AMREX_SPACEDIM+1; //amrDataI.NComp();
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

#if (AMREX_SPACEDIM == 3)
            FORT_VISCBENCH(&time, &mu, &unifDir,
                           lo, hi, &nComp,
                           ((*dataE[iLevel])[iGrid]).dataPtr(),
                           AMREX_ARLIM(lo), AMREX_ARLIM(hi),
                           delI.dataPtr(),
                           xlo.dataPtr(), xhi.dataPtr());
#else
            FORT_VISCBENCH(&time, &mu,
                           lo, hi, &nComp,
                           ((*dataE[iLevel])[iGrid]).dataPtr(),
                           AMREX_ARLIM(lo), AMREX_ARLIM(hi),
                           delI.dataPtr(),
                           xlo.dataPtr(), xhi.dataPtr());
#endif
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
        Real Ln = 0.0;

        for (int iComp=0; iComp<nComp; ++iComp)
        {
            if(norm == 0){
                // L_{inf} norm
                Ln = (*error[iLevel]).norm0(iComp);
            }
            else if(norm == 1){
                // L_1 norm
                Ln = (*error[iLevel]).norm1(iComp)*cellvol;
            }
            else{
                // L_2 norm
                Ln = (*error[iLevel]).norm2(iComp)*sqrt(cellvol);
            }
            std::cout << Ln << "  ";
        }
        std::cout << std::endl;
    }

    //
    // Write Plot Files
    Vector<std::string> varNames{ AMREX_D_DECL("x_velocity", "y_velocity", "z_velocity"), "density"};
    //
    if (!errFile.empty())
        WritePlotFile(error, amrDataI, errFile, verbose, varNames);

    if (!exFile.empty())
        WritePlotFile(dataE, amrDataI, exFile, verbose, varNames);


    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel) {
        delete error[iLevel];
        delete dataE[iLevel];
    }

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
