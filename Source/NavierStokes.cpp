//
// "Divu_Type" means S, where divergence U = S
// "Dsdt_Type" means pd S/pd t, where S is as above
//
#include <unistd.h>

#include <algorithm>
#include <vector>
#include <cmath>

#include <AMReX_Geometry.H>
#include <AMReX_Extrapolater.H>
#include <AMReX_ParmParse.H>
#include <NavierStokes.H>
#include <NAVIERSTOKES_F.H>
#include <AMReX_BLProfiler.H>
#include <PROB_NS_F.H>
#include <NS_util.H>
#include <AMReX_BCUtil.H>
#include <AMReX_MultiFabUtil.H>
#ifdef BL_USE_VELOCITY
#include <AMReX_DataServices.H>
#include <AMReX_AmrData.H>
#endif

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#include <iamr_mol.H>
#endif

#include <AMReX_buildInfo.H>

#include <iamr_godunov.H>

using namespace amrex;

namespace
{
    bool initialized = false;
}

void
NavierStokes::variableCleanUp ()
{
    NavierStokesBase::variableCleanUp ();
}

void
NavierStokes::Initialize ()
{
    if (initialized) return;

    NavierStokesBase::Initialize();

    NavierStokesBase::Initialize_specific();

    amrex::ExecOnFinalize(NavierStokes::Finalize);

    initialized = true;
}

void
NavierStokes::Finalize ()
{
    initialized = false;
}

NavierStokes::NavierStokes () {}

NavierStokes::NavierStokes (Amr&            papa,
                            int             lev,
                            const Geometry& level_geom,
                            const BoxArray& bl,
                            const DistributionMapping& dm,
                            Real            time)
    :
    NavierStokesBase(papa,lev,level_geom,bl,dm,time)
{ }

NavierStokes::~NavierStokes () { }

//
// This function initializes the State and Pressure with data.
//
void
NavierStokes::initData ()
{
    //
    // Initialize the state and the pressure.
    //
    int testNumber;
    ParmParse pp("ns");
    pp.query("testNumber", testNumber);
    int         ns       = NUM_STATE - BL_SPACEDIM;
    const Real* dx       = geom.CellSize();
    MultiFab&   S_new    = get_new_data(State_Type);
    MultiFab&   P_new    = get_new_data(Press_Type);
    const Real  cur_time = state[State_Type].curTime();
    Vector<Vector<double>> solitonInit;
    Vector<double> cellData(3);

    solitonInit = readDataFile();
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter snewmfi(S_new,true); snewmfi.isValid(); ++snewmfi)
    {
        const Box& vbx = snewmfi.tilebox();

        FArrayBox& Sfab = S_new[snewmfi];
        FArrayBox& Pfab = P_new[snewmfi];

        Sfab.setVal<RunOn::Host>(0.0,snewmfi.growntilebox(),0,S_new.nComp());
        Pfab.setVal<RunOn::Host>(0.0,snewmfi.grownnodaltilebox(-1,P_new.nGrow()));

        RealBox    gridloc = RealBox(vbx,geom.CellSize(),geom.ProbLo());
        const int* lo      = vbx.loVect();
        const int* hi      = vbx.hiVect();
        const int* s_lo    = Sfab.loVect();
        const int* s_hi    = Sfab.hiVect();
        const int* p_lo    = Pfab.loVect();
        const int* p_hi    = Pfab.hiVect();

        // FORT_INITDATA (&level,&cur_time,lo,hi,&ns,
        //                Sfab.dataPtr(Xvel),
        //                Sfab.dataPtr(BL_SPACEDIM),
        //                ARLIM(s_lo), ARLIM(s_hi),
        //                Pfab.dataPtr(),
        //                ARLIM(p_lo), ARLIM(p_hi),
        //                dx,gridloc.lo(),gridloc.hi() );
        if(testNumber == 16)
        {
            Array4<Real> const& stateArray = Sfab.array();
            Dim3 lo = lbound(vbx);
            Dim3 hi = ubound(vbx);
            const Real* xlo = gridloc.lo();
            const Real* xhi = gridloc.hi();
            Real x, y;
            for(int i = lo.x; i <= hi.x; i++)
            {
                for(int j = lo.y; j <= hi.y; j++)
                {
                    for(int k = lo.z; k <= hi.z; k++)
                    {
                        x = xlo[0] + dx[0]*((i-lo.x) + 0.5);
                        y = xlo[1] + dx[1]*((j-lo.y) + 0.5);
                        // if(y < 0)
                        // {
                        //     stateArray(i,j,k,Density) = 999.2;
                        //
                        // }
                        // else
                        // {
                        //     stateArray(i,j,k,Density) = 1.225;
                        //
                        // }
                        lookUpData(solitonInit, x, y, cellData);
                        stateArray(i,j,k,Density) = cellData[0];
                        stateArray(i,j,k,Xvel) = cellData[1];
                        stateArray(i,j,k,Yvel) = cellData[2];
                        if(cellData[0] == 1000.0)
                        {
                            stateArray(i,j,k,Tracer) = 0.0;
                        }
                        else
                        {
                            stateArray(i,j,k,Tracer) = 1.0;

                        }
                    }
                }
            }
        }
        else
        {
            initialState(Sfab, vbx, gridloc, dx, testNumber);
        }

    }

#ifdef AMREX_USE_EB
    set_body_state(S_new);
#endif

#ifdef BL_USE_VELOCITY
    //
    // We want to add the velocity from the supplied plotfile
    // to what we already put into the velocity field via FORT_INITDATA.
    //
    // This code has a few drawbacks.  It assumes that the physical
    // domain size of the current problem is the same as that of the
    // one that generated the pltfile.  It also assumes that the pltfile
    // has at least as many levels (with the same refinement ratios) as does
    // the current problem.  If either of these are false this code is
    // likely to core dump.
    //
    ParmParse pp("ns");

    std::string velocity_plotfile;
    pp.query("velocity_plotfile", velocity_plotfile);

    std::string velocity_plotfile_xvel_name = "x_velocity";
    pp.query("velocity_plotfile_xvel_name", velocity_plotfile_xvel_name);

    Real velocity_plotfile_scale(1.0);
    pp.query("velocity_plotfile_scale",velocity_plotfile_scale);

    if (!velocity_plotfile.empty())
    {
        Print() << "initData: reading data from: " << velocity_plotfile << " ("
                << velocity_plotfile_xvel_name << ")" << '\n';

        DataServices::SetBatchMode();
        Amrvis::FileType fileType(Amrvis::NEWPLT);
        DataServices dataServices(velocity_plotfile, fileType);

        if (!dataServices.AmrDataOk())
            //
            // This calls ParallelDescriptor::EndParallel() and exit()
            //
            DataServices::Dispatch(DataServices::ExitRequest, NULL);

        AmrData&           amrData   = dataServices.AmrDataRef();
        Vector<std::string> plotnames = amrData.PlotVarNames();

        int idX = -1;
        for (int i = 0; i < plotnames.size(); ++i)
            if (plotnames[i] == velocity_plotfile_xvel_name) idX = i;

        if (idX == -1)
	       Abort("Could not find velocity fields in supplied velocity_plotfile");
	      else
	       Print() << "Found " << velocity_plotfile_xvel_name << ", idX = " << idX << '\n';

        MultiFab tmp(S_new.boxArray(), S_new.DistributionMap(), 1, 0);
        for (int i = 0; i < BL_SPACEDIM; i++)
        {
	    amrData.FillVar(tmp, level, plotnames[idX+i], 0);

	    MultiFab::Saxpy(S_new, velocity_plotfile_scale, tmp, 0, Xvel+i, 1, 0)

	    amrData.FlushGrids(idX+i);
        }

	Print() << "initData: finished init from velocity_plotfile" << '\n';
    }
#endif /*BL_USE_VELOCITY*/

    make_rho_prev_time();
    make_rho_curr_time();
    //
    // Initialize divU and dSdt.
    //
    if (have_divu)
    {
        const Real dt       = 1.0;
        const Real dtin     = -1.0; // Dummy value denotes initialization.
        const Real curTime = state[Divu_Type].curTime();
        MultiFab&  Divu_new = get_new_data(Divu_Type);

        state[State_Type].setTimeLevel(curTime,dt,dt);

        //Make sure something reasonable is in diffn_ec
        calcDiffusivity(cur_time);

        calc_divu(cur_time,dtin,Divu_new);

        if (have_dsdt)
            get_new_data(Dsdt_Type).setVal(0);
    }

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
    {
        get_new_data(Dpdt_Type).setVal(0);
    }

    is_first_step_after_regrid = false;
    old_intersect_new          = grids;

#ifdef AMREX_PARTICLES
    initParticleData ();
#endif
}


Vector<Vector<double>>
NavierStokes::readDataFile()
{
    // For the specific case of the soltion wave, where
    // initial data is generated from a different script file
    std::ifstream solitonInitDat;
    double iRead = 0.0;
    double xRead = 0.0;
    double yRead = 0.0; double rhoRead = 0.0; double uRead = 0.0;
    double vRead =0.0; double pRead = 0.0;
    Vector<double> x, y, rho, u, v;
    int Np, Nv;
    Vector<Vector<double>> solitonData;
    solitonInitDat.open("plotuv.dat");

    while(solitonInitDat >> iRead >> xRead >> yRead >> rhoRead
        >> uRead >> vRead >> pRead)
    {
        x.push_back(xRead);
        y.push_back(yRead);
        rho.push_back(rhoRead);
        u.push_back(uRead);
        v.push_back(vRead);
    }

    solitonInitDat.close();
    Nv = 5;
    Np = x.size();

    solitonData.resize(Nv, Vector<double>(Np));

    for(int i = 0; i < Np; i++)
    {
        solitonData[0][i] = x[i];
        solitonData[1][i] = y[i];
        solitonData[2][i] = rho[i];
        solitonData[3][i] = u[i];
        solitonData[4][i] = v[i];
        // std::cout << x[i] << " ";
    }

    return solitonData;
}


void
NavierStokes::lookUpData(Vector<Vector<double>>& solitonData, double& x, double& y,
    Vector<double>& cellData)
{

//std::cout << "Looking up data for soliton test case..." << std::endl;
    int Np = solitonData[0].size();
    int xInd = 0;
    int yInd = 0;
    int ind = 0;
    double yDat;
    double xDat;
    Vector<Real> diffY(Np);
    Vector<Real> diffX(Np);
    Vector<Real> diff(Np);
    for(int i = 0; i < Np; i++)
    {
        xDat = solitonData[0][i];
        yDat = solitonData[1][i];
        diffX[i] = fabs(xDat - x);
        diffY[i] = fabs(yDat - y);
        diff[i] = sqrt((xDat-x)*(xDat-x) + (yDat-y)*(yDat-y));
    }

    xInd = std::distance(std::begin(diffX),
        std::min_element(std::begin(diffX), std::end(diffX)));
    yInd = std::distance(std::begin(diffY),
        std::min_element(std::begin(diffY), std::end(diffY)));
    ind = std::distance(std::begin(diff),
        std::min_element(std::begin(diff), std::end(diff)));
    // if(xInd == yInd)
    // {
        cellData[0] = solitonData[2][ind];
        cellData[1] = solitonData[3][ind];
        cellData[2] = solitonData[4][ind];
        // std::cout << x << " " << y << " " <<
        //     solitonData[0][ind] << " " << solitonData[1][ind] << std::endl;

    // }


    //std::cout << "Successfully read in soliton data." << std::endl;

}

void
NavierStokes::initialState(FArrayBox& statein, const Box& bx,
    RealBox gridloc, const Real* dx, int testNumber)
{
    Array4<Real> const& stateArray = statein.array();
    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    const Real* xlo = gridloc.lo();
    const Real* xhi = gridloc.hi();
    Real x, y;

    // std::cout << "Break ";
    for(int i = lo.x; i <= hi.x; i++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int k = lo.z; k <= hi.z; k++)
            {
                // x = dx[0]*(i+0.5);
                // y = dx[1]*(j+0.5);
                x = xlo[0] + dx[0]*((i-lo.x) + 0.5);
                y = xlo[1] + dx[1]*((j-lo.y) + 0.5);
                if(testNumber == 0)
                {
                    // Rayleigh Problem or Stokes First Problem:
                    stateArray(i,j,Xvel) = 0.0;
                    stateArray(i,j,Density) = 1.0;
                    if(x < 0.5)
                    {
                        stateArray(i,j,Yvel) = -1.0;
                    }
                    else
                    {
                        // std::cout << "Break" << std::endl;
                        stateArray(i,j,Yvel) = 1.0;
                    }
                }

                else if(testNumber == 1)
                {
                    // Poiseuille flow
                    stateArray(i,j,Xvel) = 0.0;
                    stateArray(i,j,Yvel) = 0.0;
                    stateArray(i,j,Density) = 1.0;
                }

                else if(testNumber == 2)
                {
                    // Periodic shear layer from Bell, Colella and Glaz 1989
                    stateArray(i,j,k, Xvel) = tanh(30.0*(0.25-fabs(y-0.5)));
                    stateArray(i,j,k, Yvel) = 0.05*sin(2*M_PI*x);
                    stateArray(i,j,k, Density) = 1.0;

                }
                else if(testNumber == 3)
                {
                    // Two-dimensional shear layer from Chien et al., 1995
                    double rho1 = 1.0;
                    double rho2 = 7.0;
                    double f = 1.0;
                    double fm = 0.5;
                    double rhom, um;
                    double u1 = 1.45;
                    double u2 = 0.55;
                    double lRho = (rho1-rho2)/(rho1+rho2);
                    double l = (u1-u2)/(u1+u2);
                    double delta = 6.0;

                    rhom = (rho1+rho2)/2.0;
                    um = (u1+u2)/2.0;

                    stateArray(i,j,Density) = rhom*(1+lRho*tanh(2*y/delta));
                    stateArray(i,j,Xvel) = um*(1+l*tanh(2*y/delta));
                    stateArray(i,j,Yvel) = 0.0;
                    stateArray(i,j,Tracer) = 0.0;//0.5*(1+tanh(2*y/delta));


                }

                else if(testNumber == 4)
                {
                    // same test as above but for water and air.
                    ParmParse pp("ns");

                    double rhoL, rhoG, uL, uG;
                    pp.query("rhoL", rhoL);
                    pp.query("rhoG", rhoG);
                    pp.query("uG", uG);
                    pp.query("uL", uL);
                    double mu1 = 1.77625e-5;
                    double mu2 = 1.1377e-3;
                    double rhom, um, mum;
                    double lRho = (rhoG-rhoL)/(rhoG+rhoL);
                    double l = (uG-uL)/(uG+uL);
                    double lMu = (mu1-mu2)/(mu1+mu2);
                    double delta = 6.0;

                    rhom = (rhoL+rhoG)/2.0;
                    um = (uL+uG)/2.0;
                    mum = (mu1+mu2)/2.0;

                    stateArray(i,j,k,Density) = rhom*(1+lRho*tanh(2*y/delta));
                    stateArray(i,j,k, Xvel) = um*(1+l*tanh(2*y/delta));
                    stateArray(i,j,Tracer) = 0.5*(1+tanh(2*y/delta));

                    stateArray(i,j,Yvel) = 0.0;

                }
                else if(testNumber == 5)
                {
                    // same test as above but for water and air.
                    ParmParse pp("ns");

                    double rho2;
                    double rho1;
                    pp.query("rho1", rho1);
                    pp.query("rho2", rho2);

                    double mu1 = 1.77625e-5;
                    double mu2 = 1.1377e-3;
                    double r = 0.0015;
                    double x0 = 0.0035;
                    double y0 = 0.004;
                    double dist = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));

                    stateArray(i,j, k, Xvel) = 0.0;
                    stateArray(i,j, k, Yvel) = 0.0;
                    stateArray(i,j, k, Density) = rho1 + 0.5*(rho2-rho1)*(1-tanh(10000.0*(dist-r)));
                    stateArray(i,j, k, Tracer) = 0.5*(1-tanh(10000.0*(dist-r)));
                    // if(dist < r)
                    // {
                    //     stateArray(i,j, k, Tracer) = 1.0;
                    //     // stateArray(i,j, k, Viscosity) = mu1;
                    //     stateArray(i,j, k, Density) = rho2;
                    // }
                    // else
                    // {
                    //     stateArray(i,j, k, Tracer) = 0.0;
                    //     // stateArray(i,j, k, Viscosity) = mu2;
                    //     stateArray(i,j, k, Density) = rho1;
                    // }

                }
                else if(testNumber == 6)
                {
                    double rho1 = 1.14;
                    double rho2 = 1000.0;
                    double mu1 = 1.77625e-5;
                    double mu2 = 1.1377e-3;
                    double u1 = 6.532;
                    double u2 = 0.532;

                    if(y > 0)
                    {
                        stateArray(i,j, Density) = rho1;
                        // stateArray(i,j, Viscosity) = mu1;
                        stateArray(i,j, Tracer) = 1.0;
                        stateArray(i,j, Xvel) = u1;
                    }
                    else
                    {
                        stateArray(i,j, Density) = rho2;
                        // stateArray(i,j,Viscosity) = mu2;
                        stateArray(i,j, Tracer) = 0.0;
                        stateArray(i,j, Xvel) = u2;
                    }

                    stateArray(i,j, Yvel) = 0.0;
                }
                else if(testNumber == 7)
                {
                    double rho1 = 1.0;
                    double rho2 = 1000.0;
                    double mu1 = 1.77625e-5;
                    double mu2 = 1.1377e-3;
                    double u1 = 0.1;
                    double u2 = 0.01;
                    double epsilon = 0.025;
                    double ySurf = 0.25 + epsilon*cos(2.0*M_PI*x);
                    if(y < ySurf)
                    {
                        stateArray(i, j, Density) = rho2;
                        stateArray(i, j, Tracer) = 1.0;
                    }
                    else
                    {
                       stateArray(i, j, Density) = rho1;
                       stateArray(i, j, Tracer) = 0.0;
                    }

                    stateArray(i, j, Xvel) = 0.0;
                    stateArray(i, j, Yvel) = 0.0;
                    // stateArray(i, j, Viscosity) = 0.0;


                }

                else if(testNumber == 8)
                {
                    double Lx = 0.01;
                    ParmParse pp("ns");

                    double rhoG, muG;
                    double rhoL, muL;
                    pp.query("rhoG", rhoG);
                    pp.query("rhoL", rhoL);
                    pp.query("muG", muG);
                    pp.query("muL", muL);

                    double pertHeight = 2.0*Lx + 0.05*Lx*cos(2.0*M_PI*x/Lx);

                    stateArray(i, j, k, Density) = rhoG + 0.5*(rhoL-rhoG)*(1+tanh((y-pertHeight)/(0.0005)));
                    stateArray(i,j, k,  Tracer) = 1.0 - 0.5*(1+tanh((y-pertHeight)/(0.0005)));
                    // stateArray(i,j, k,  Tracer2) = muG + 0.5*(muL-muG)*(1+tanh((y-pertHeight)/0.05));
                    // stateArray(i, j, Viscosity) = mu1 + 0.5*(mu2-mu1)*(1+tanh((y-pertHeight)/0.001));
                    // if(y < pertHeight)
                    // {
                    //     stateArray(i,j, Tracer) = 1.0;
                    //     // stateArray(i,j,Density) = rhoG;
                    // }
                    // else
                    // {
                    //     stateArray(i, j, Tracer) = 0.0;
                    //     // stateArray(i,j, Density) = rhoL;
                    // }


                }

                else if(testNumber == 9)
                {
                    ParmParse pp("ns");

                    double rho2;
                    double rho1;
                    int do_variable_den = 0;
                    pp.query("rho1", rho1);
                    pp.query("rho2", rho2);
                    double x0 = 0.0;
                    double y0 = 0.0;
                    double f = 0.5;
                    double r = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
                    double R = 0.4;
                    double psi = atan2((y-y0),(x-x0));
                    double uPsi;
                    if(r < 0.2 && r >= 0)
                    {
                        uPsi = 5.0*r;
                    }
                    else if(r < 0.4 && r >= 0.2)
                    {
                        uPsi = 2.0-5.0*r;
                    }
                    else
                    {
                        uPsi = 0.0;
                    }

                    stateArray(i,j,k, Xvel) = -1.0*uPsi*sin(psi);
                    stateArray(i,j,k, Yvel) = uPsi*cos(psi);

                    // stateArray(i,j,k, Density) = f*rho1 + (1-f)*rho2;
                    // stateArray(i,j,k, Tracer) = 0.5;
                    stateArray(i,j, k, Density) = rho1 + 0.5*(rho2-rho1)*(1-tanh(50.0*(r-R)));
                    stateArray(i,j, k, Tracer) = 0.5*(1-tanh(50.0*(r-R)));

                }
                else if(testNumber == 10)
                {
                    ParmParse pp("ns");

                    double rho2;
                    double rho1;
                    int do_variable_den = 0;
                    pp.query("rho1", rho1);
                    pp.query("rho2", rho2);
                    double x0 = 0.0;
                    double y0 = 0.0;
                    double f = 0.5;
                    double r = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
                    double R = 0.4;
                    double psi = atan2((y-y0),(x-x0));
                    double uPsi;
                    if(r < 0.2 && r >= 0)
                    {
                        uPsi = 5.0*r;
                    }
                    else if(r < 0.4 && r >= 0.2)
                    {
                        uPsi = 2.0-5.0*r;
                    }
                    else
                    {
                        uPsi = 0.0;
                    }

                    stateArray(i,j,k, Xvel) = -1.0*uPsi*sin(psi) + 0.3;
                    stateArray(i,j,k, Yvel) = uPsi*cos(psi);

                    // stateArray(i,j,k, Density) = f*rho1 + (1-f)*rho2;
                    // stateArray(i,j,k, Tracer) = f;
                    stateArray(i,j, k, Density) = rho1 + 0.5*(rho2-rho1)*(1-tanh(50.0*(r-R)));
                    stateArray(i,j, k, Tracer) = 0.5*(1-tanh(50.0*(r-R)));

                }

                else if(testNumber == 11)
                {
                    // Boeck et al. (2007) Sect 4.3, balanced density
                    // ligament formation
                    double rhoG, rhoL, uG, uL, ampl, lbda;
                    Vector<double> domhi(BL_SPACEDIM);
                    Vector<double> domlo(BL_SPACEDIM);
                    ParmParse pp("ns");
                    ParmParse ppGeom("geometry");
                    pp.query("rhoG", rhoG);
                    pp.query("rhoL", rhoL);
                    pp.query("uG", uG);
                    pp.query("uL", uL);
                    pp.query("ampl", ampl);
                    pp.query("lbda", lbda);
                    ppGeom.getarr("prob_hi", domhi, 0, BL_SPACEDIM);
                    ppGeom.getarr("prob_lo", domlo, 0, BL_SPACEDIM);
                    double Lx = domhi[0]-domlo[0];
                    double deltaG = 0.2*Lx;
                    double lRho = (rhoG-rhoL)/(rhoG+rhoL);
                    double l = (uG-uL)/(uG+uL);
                    double rhom = (rhoL+rhoG)/2.0;
                    double um = (uL+uG)/2.0;
                    double delta = 0.01;
                    double pertHeight = 0.05*Lx*cos(2.0*M_PI*x/Lx);

                    // stateArray(i,j,k, Density) = rhom*(1+lRho*tanh(2*y/delta));
                    // stateArray(i,j,k,Xvel) = um*(1+l*tanh(2*y/delta));
                    // stateArray(i,j,k,Tracer) = 0.5*(1+tanh(2*y/delta));
                    // stateArray(i,j,k,Density) = rhoG;
                    if(y < 0)
                    {
                        stateArray(i,j,k,Xvel) = uL;
                        stateArray(i,j,k,Tracer) = 0.0;
                        stateArray(i,j,k,Density) = rhoL;

                    }
                    else if(y >= 0 && y < deltaG)
                    {
                        stateArray(i,j,k,Xvel) = uG;//(uG/deltaG)*y;
                        stateArray(i,j,k,Density) = rhoG;
                        stateArray(i,j,k, Tracer) = 1.0;
                    }
                    else
                    {
                        stateArray(i,j,k,Xvel) = uG;
                        stateArray(i,j,k, Tracer) = 1.0;
                        stateArray(i,j,k,Density) = rhoG;

                    }
                    stateArray(i,j,k,Yvel) = ampl*sin(2.0*M_PI*(x+0.5)/lbda);


                }

                else if(testNumber == 12)
                {
                    //Boeck case E
                    double rhoG, rhoL, uG, uL, circ;
                    Vector<double> domhi(BL_SPACEDIM);
                    Vector<double> domlo(BL_SPACEDIM);
                    ParmParse pp("ns");
                    ParmParse ppGeom("geometry");
                    pp.query("rhoG", rhoG);
                    pp.query("rhoL", rhoL);
                    pp.query("uG", uG);
                    pp.query("uL", uL);
                    pp.query("circ", circ);
                    ppGeom.getarr("prob_hi", domhi, 0, BL_SPACEDIM);
                    ppGeom.getarr("prob_lo", domlo, 0, BL_SPACEDIM);
                    double Lx = domhi[0]-domlo[0];
                    double Rc = Lx*4e-2;
                    double uPsi, y0;
                    double lbda = 0.5;
                    Vector<double> xVort(4), yVort(4);
                    Vector<double> dist(4);
                    Vector<double> psi(4);
                    double omega = ((uG-uL)*Lx*circ)/(Rc*Rc);
                    y0 =Lx/4.0;
                    //double omega = 0.90625;

                    if(y < Lx/3)// && r0 >= Rc && r1 >= Rc && r2 >= Rc && r3 >= Rc)
                    {
                        stateArray(i,j,k,Xvel) = uL*erf((Lx/3-y)/(0.2*Lx));
                        stateArray(i,j,k,Tracer) = 0.0;
                        stateArray(i,j,k,Density) = rhoL;

                    }

                    // else if(y >= Lx/3 && y < Lx/3 + 0.2*Lx)
                    // {
                    //     stateArray(i,j,k,Xvel) = uG*erf((y-Lx/3)/(0.1*Lx));
                    //     stateArray(i,j,k, Tracer) = 1.0;
                    //     stateArray(i,j,k,Density) = rhoG;
                    // }
                    else
                    {
                        stateArray(i,j,k,Xvel) = uG*erf((y-Lx/3)/(0.1*Lx));
                        stateArray(i,j,k, Tracer) = 1.0;
                        stateArray(i,j,k,Density) = rhoG;
                    }

                    // stateArray(i,j,k,Yvel) = 0.025*sin(2.0*M_PI*(x+0.5)/lbda);

                    //
                    xVort[0] = 0.15*Lx;
                    xVort[1] = 0.35*Lx;
                    xVort[2] = 0.65*Lx;
                    xVort[3] = 0.85*Lx;
                    // yVort[0] = 0.2*Lx;
                    // yVort[1] = 0.1*Lx;
                    // yVort[2] = 0.2*Lx;
                    // yVort[3] = 0.1*Lx;
                    for(int n = 0; n < 4; n++)
                    {
                        // xVort[n] = (n+1)*0.2*Lx;
                        yVort[n] = y0;
                        // yVort[n] = (n == 0 || n == 2) ? 0.2*Lx : 0.1*Lx;
                        // xVort[n] = (n == 1 || n == 3) ? 0.25*Lx : 0.75*Lx;
                        dist[n] = sqrt((x-xVort[n])*(x-xVort[n]) + (y-yVort[n])*(y-yVort[n]));
                        psi[n] = atan2((y-yVort[n]), (x-xVort[n]));

                        if(dist[n] <= Rc)
                        {
                           uPsi =  pow(-1.0,n)*omega*dist[n];
                            stateArray(i,j,k,Xvel) = -1.0*uPsi*sin(psi[n]);
                            stateArray(i,j,k,Yvel) = uPsi*cos(psi[n]);
                        }


                    }




                }


                else if(testNumber == 13)
                {
                    double rhoG, rhoL, uG, uL, circ;
                    Vector<double> domhi(BL_SPACEDIM);
                    Vector<double> domlo(BL_SPACEDIM);
                    ParmParse pp("ns");
                    ParmParse ppGeom("geometry");
                    pp.query("rhoG", rhoG);
                    pp.query("rhoL", rhoL);
                    pp.query("uG", uG);
                    pp.query("uL", uL);
                    pp.query("circ", circ);
                    ppGeom.getarr("prob_hi", domhi, 0, BL_SPACEDIM);
                    ppGeom.getarr("prob_lo", domlo, 0, BL_SPACEDIM);
                    double Lx = domhi[0]-domlo[0];
                    double Rc = Lx*4e-2;
                    double uPsi, y0;
                    double lbda = 0.5;
                    Vector<double> xVort(4), yVort(4);
                    Vector<double> dist(4);
                    Vector<double> psi(4);
                    double omega = ((uG-uL)*Lx*circ)/(Rc*Rc);
                    y0 =Lx/4.0;
                    //double omega = 0.90625;

                    if(y < Lx/3)// && r0 >= Rc && r1 >= Rc && r2 >= Rc && r3 >= Rc)
                    {
                        stateArray(i,j,k,Xvel) = uL;//*erf((Lx/3-y)/(0.2*Lx));
                        stateArray(i,j,k,Tracer) = 0.0;
                        stateArray(i,j,k,Density) = rhoL;

                    }

                    else if(y >= Lx/3 && y < Lx/3 + 0.2*Lx)
                    {
                        stateArray(i,j,k,Xvel) = uG*(y-Lx/3)/(0.2*Lx);//*erf((y-Lx/3)/(0.1*Lx));
                        stateArray(i,j,k, Tracer) = 1.0;
                        stateArray(i,j,k,Density) = rhoG;
                    }
                    else
                    {
                        stateArray(i,j,k,Xvel) = uG;//*(y-Lx/3)/(0.2*Lx);//*erf((y-Lx/3)/(0.1*Lx));
                        stateArray(i,j,k, Tracer) = 1.0;
                        stateArray(i,j,k,Density) = rhoG;
                    }

                    // stateArray(i,j,k,Yvel) = 0.025*sin(2.0*M_PI*(x+0.5)/lbda);

                    //
                    xVort[0] = 0.15*Lx;
                    xVort[1] = 0.35*Lx;
                    xVort[2] = 0.65*Lx;
                    xVort[3] = 0.85*Lx;
                    // yVort[0] = 0.2*Lx;
                    // yVort[1] = 0.1*Lx;
                    // yVort[2] = 0.2*Lx;
                    // yVort[3] = 0.1*Lx;
                    for(int n = 0; n < 4; n++)
                    {
                        // xVort[n] = (n+1)*0.2*Lx;
                        yVort[n] = y0;
                        // yVort[n] = (n == 0 || n == 2) ? 0.2*Lx : 0.1*Lx;
                        // xVort[n] = (n == 1 || n == 3) ? 0.25*Lx : 0.75*Lx;
                        dist[n] = sqrt((x-xVort[n])*(x-xVort[n]) + (y-yVort[n])*(y-yVort[n]));
                        psi[n] = atan2((y-yVort[n]), (x-xVort[n]));

                        if(dist[n] <= Rc)
                        {
                            uPsi = (n == 1 || n == 2) ? omega*dist[n] : -1.0*omega*dist[n];
                            // uPsi =  pow(-1.0,n+1)*omega*dist[n];
                            stateArray(i,j,k,Xvel) = -1.0*uPsi*sin(psi[n]);
                            stateArray(i,j,k,Yvel) = uPsi*cos(psi[n]);
                        }


                    }
                }

                else if(testNumber == 14)
                {
                    // Dam break
                    double rhoG, rhoL, uG, uL, circ, x0, y0;
                    Vector<double> domhi(BL_SPACEDIM);
                    Vector<double> domlo(BL_SPACEDIM);
                    ParmParse pp("ns");
                    ParmParse ppGeom("geometry");
                    pp.query("rhoG", rhoG);
                    pp.query("rhoL", rhoL);
                    pp.query("x0", x0);
                    pp.query("y0", y0);

                    if(x <= x0 && y <= y0)
                    {
                        // stateArray(i, j, k, Density) = rhoL + 0.5*(rhoG-rhoL)*(1+tanh((y-0.6)/(0.0005)));
                        stateArray(i,j,k, Density) = rhoL;
                        stateArray(i,j,k, Tracer) = 0.0;
                        if(do_trac3)
                        {
                            if(y <= 0.5*y0)
                            {
                                stateArray(i,j,k, Tracer3)  = 1.0;
                            }
                            else
                            {
                                stateArray(i,j,k, Tracer2)  = 1.0;

                            }
                        }

                    }
                    else
                    {
                        stateArray(i,j,k, Density) = rhoG;
                        stateArray(i,j,k, Tracer)  = 1.0;
                    }
                }

                else if(testNumber == 15)
                {
                    double rhoG, rhoL, hL;
                    ParmParse pp("ns");
                    ParmParse ppGeom("geometry");
                    pp.query("rhoG", rhoG);
                    pp.query("rhoL", rhoL);
                    pp.query("hL", hL);
                    if(y <= hL)
                    {
                        stateArray(i,j,k, Density) = rhoL;
                        stateArray(i,j,k, Tracer) = 0.0;
                    }
                    else
                    {
                        stateArray(i,j,k, Density) = rhoG;
                        stateArray(i,j,k, Tracer) = 1.0;
                    }
                    stateArray(i,j,k, Xvel) = 1.0;
                }

                else if(testNumber == 17)
                {
                    double rhoG, rhoL, H, h, L0;
                    ParmParse pp("ns");
                    ParmParse ppGeom("geometry");
                    pp.query("rhoG", rhoG);
                    pp.query("rhoL", rhoL);
                    pp.query("H", H);
                    pp.query("hB", h);
                    pp.query("L0", L0);

                    if(x < L0 && y < H)
                    {
                        stateArray(i,j,k, Density) = rhoL;
                        stateArray(i,j,k, Tracer) = 0.0;
                        stateArray(i,j,k, Tracer2) = 1.0;
                        stateArray(i,j,k, Tracer3) = 0.0;
                    }

                    else if(x >= L0 && y < h)
                    {
                        stateArray(i,j,k, Density) = rhoL;
                        stateArray(i,j,k, Tracer) = 0.0;
                        stateArray(i,j,k, Tracer2) = 0.0;
                        stateArray(i,j,k, Tracer3) = 1.0;
                    }
                    else
                    {
                        stateArray(i,j,k, Density) = rhoG;
                        stateArray(i,j,k, Tracer) = 1.0;
                        stateArray(i,j,k, Tracer2) = 0.0;
                        stateArray(i,j,k, Tracer3) = 0.0;
                    }

                }

                else if(testNumber == 18)
                {
                    //Equilibrium rod, surface tension validation
                    double rho1, rho2, sigma, R, x0, y0;
                    ParmParse pp("ns");
                    ParmParse ppGeom("geometry");
                    pp.query("rho1",    rho1);
                    pp.query("rho2",    rho2);
                    pp.query("R",          R);
                    pp.query("sigma",   sigma);
                    pp.query("x0",         x0);
                    pp.query("y0",         y0);
                    double dist = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));

                    stateArray(i,j, k, Density) = rho1 + 0.5*(rho2-rho1)*(1-tanh(1000.0*(dist-R)));
                    stateArray(i,j, k, Tracer) = 0.5*(1-tanh(1000.0*(dist-R)));

                    stateArray(i,j, k, Xvel) = 0.0;
                    stateArray(i,j, k, Yvel) = 0.0;
                    // if(dist < R)
                    // {
                    //     stateArray(i,j, k, Tracer) = 1.0;
                    //     // stateArray(i,j, k, Viscosity) = mu1;
                    //     stateArray(i,j, k, Density) = rho2;
                    // }
                    // else
                    // {
                    //     stateArray(i,j, k, Tracer) = 0.0;
                    //     // stateArray(i,j, k, Viscosity) = mu2;
                    //     stateArray(i,j, k, Density) = rho1;
                    // }
                }
                else if(testNumber == 19)
                {
                    //Equilibrium rod, surface tension validation
                    double rho1, rho2, sigma, R, x0, y0;
                    ParmParse pp("ns");
                    ParmParse ppGeom("geometry");
                    pp.query("rho1",    rho1);
                    pp.query("rho2",    rho2);
                    pp.query("R",          R);
                    pp.query("sigma",   sigma);
                    pp.query("x0",         x0);
                    pp.query("y0",         y0);
                    double dist = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
                    if(do_trac3)
                    {
                        if(y < y0)
                        {
                            // stateArray(i,j, k, Tracer) = 0.0;
                            stateArray(i,j, k, Tracer2) = 1.0 - 0.5*(1-tanh(2500.0*(dist-R)));
                            stateArray(i,j, k, Tracer3) = 0.0;
                            // stateArray(i,j, k, Density) = rho1;
                        }
                        else
                        {
                            stateArray(i,j, k, Tracer3) = 1.0 - 0.5*(1-tanh(2500.0*(dist-R)));
                            stateArray(i,j, k, Tracer2) = 0.0;
                            // stateArray(i,j, k, Tracer) = 0.0;
                            // stateArray(i,j, k, Density) = rho1;
                        }
                        stateArray(i,j, k, Density) = rho1 + 0.5*(rho2-rho1)*(1-tanh(2500.0*(dist-R)));
                        // if(nTrac == 3)
                        // {
                        stateArray(i,j, k, Tracer) = 0.5*(1-tanh(2500.0*(dist-R)));
                    }
                    else
                    {
                        stateArray(i,j, k, Density) = rho1 + 0.5*(rho2-rho1)*(1-tanh(2500.0*(dist-R)));
                        // if(nTrac == 3)
                        // {
                        stateArray(i,j, k, Tracer) = 0.5*(1-tanh(2500.0*(dist-R)));
                    }

                }

                else if(testNumber == 20)
                {
                    double rhow, rhoo, rhog;
                    ParmParse pp("ns");
                    ParmParse ppGeom("geometry");
                    pp.query("rhow",    rhow);
                    pp.query("rhoo",    rhoo);
                    pp.query("rhog",    rhog);
                    double wallthck = 0.01;
                    double w1 = 0.50;
                    double h1 = 0.1; double h2 = 0.25;
                    if(x <= w1)
                    {
                        stateArray(i,j, k, Tracer) = 1.0;
                        stateArray(i,j, k, Density) = rhow;
                    }
                    else
                    {
                        stateArray(i,j, k, Tracer2) =  1.0;
                        stateArray(i,j, k, Density) = rhoo;

                    }

                }
                else if(testNumber == 21)
                {
                    double rhow, rhoo, rhog;
                    ParmParse pp("ns");
                    ParmParse ppGeom("geometry");
                    pp.query("rhow",    rhow);
                    pp.query("rhoo",    rhoo);
                    pp.query("rhog",    rhog);
                    double wallthck = 0.01;
                    double w1 = 1.50;
                    double h1 = 0.1; double h2 = 0.25;
                    if(y <= h1 && x <= w1)
                    {
                        stateArray(i,j, k, Tracer2) = 1.0;
                        stateArray(i,j, k, Density) = rhow;
                    }
                    else if(x  > w1 && y <= h2)
                    {
                        stateArray(i,j, k, Tracer3) = 1.0;
                        stateArray(i,j, k, Density) = rhoo;
                    }
                    else
                    {
                        stateArray(i,j, k, Tracer) = 1.0;
                        stateArray(i,j, k, Density) = rhog;

                    }
                }
                else if(testNumber == 22)
                {
                    double rhow, rhoo, rhog, h;
                    ParmParse pp("ns");
                    ParmParse ppGeom("geometry");
                    pp.query("rhow",    rhow);
                    pp.query("rhoo",    rhoo);
                    pp.query("rhog",    rhog);
                    pp.query("h", h);
                    if(y <= h)
                    {
                        stateArray(i,j, k, Density) = rhow;
                        stateArray(i,j, k, Tracer2) = 1.0;
                    }
                    else if(y <= 2.0*h && y > h)
                    {
                        stateArray(i,j, k, Density) = rhoo;
                        stateArray(i,j, k, Tracer3) = 1.0;
                    }
                    else
                    {
                        stateArray(i,j, k, Density) = rhog;
                        stateArray(i,j, k, Tracer) = 1.0;
                    }
                    // stateArray(i,j, k, Yvel) = 0.1;

                }

                else if(testNumber == 23)
                {
                    //Oscillating bubble
                    double rho1, rho2, sigma, R, x0, y0;
                    ParmParse pp("ns");
                    ParmParse ppGeom("geometry");
                    pp.query("rho1",    rho1);
                    pp.query("rho2",    rho2);
                    pp.query("R",          R);
                    pp.query("sigma",   sigma);
                    pp.query("x0",         x0);
                    pp.query("y0",         y0);
                    double dist = sqrt(((x-x0)*(x-x0))/(0.2*0.2)+((y-y0)*(y-y0))/(0.12*0.12));

                        stateArray(i,j, k, Density) = rho1 + 0.5*(rho2-rho1)*(1-tanh(50.0*(dist-1.0)));

                            stateArray(i,j, k, Tracer) = 1.0-0.5*(1-tanh(50.0*(dist-1.0)));
                            stateArray(i,j, k, Tracer2) = 0.5*(1-tanh(50.0*(dist-1.0)));

                }

                else if(testNumber == 24)
                {
                    double rho1, rho2, sigma, R, x0, y0;
                    ParmParse pp("ns");
                    ParmParse ppGeom("geometry");
                    pp.query("rho1",    rho1);
                    pp.query("rho2",    rho2);
                    pp.query("R",          R);
                    pp.query("sigma",   sigma);
                    pp.query("x0",         x0);
                    pp.query("y0",         y0);

                    if(fabs(x-x0) <= R && fabs(y-y0) <= R)
                    {
                        stateArray(i,j, k, Density) = rho1;
                        stateArray(i,j,k,Tracer2)   = 1.0;
                        stateArray(i,j, k, Tracer) = 0.0;

                    }
                    else
                    {
                        stateArray(i,j, k, Density) = rho2;
                        stateArray(i,j, k, Tracer) = 1.0;
                        stateArray(i,j,k,Tracer2)   = 0.0;

                    }

                }
                else if(testNumber == 25)
                {
                    ParmParse pp("ns");
                    double rhow = 0.0;
                    pp.query("rhow",    rhow);
                    stateArray(i,j,k, Density) = rhow;
                    stateArray(i,j,k, Tracer) = 1.0;
                    stateArray(i,j,k, Xvel) = 0.0;
                    stateArray(i,j,k, Yvel) = 0.0;
                }
            }
        }
    }
    // std::cout << "Break" << std::endl;
}

//
// ADVANCE FUNCTIONS
//

//
// This function ensures that the multifab registers and boundary
// flux registers needed for syncing the composite grid
//
//     u_mac, Vsync, Ssync, rhoavg, fr_adv, fr_visc
//
// are initialized to zero.  In general these quantities
// along with the pressure sync registers (sync_reg) and
// advective velocity registers (mac_reg) are compiled by first
// setting them to the coarse value acquired during a coarse timestep
// and then incrementing in the fine values acquired during the
// subcycled fine timesteps.  This compilation procedure occurs in
// different parts for different quantities
//
// * u_mac is set in predict_velocity and mac_project.
// * fr_adv, fr_visc are set in velocity_advect and scalar_advect
// * Vsync, Ssync are set in subcycled calls to post_timestep
// * mac_reg is set in mac_project
// * sync_reg is set in level_project
// * rhoavg, pavg are set in advance_setup and advance
//
// After these quantities have been compiled during a coarse
// timestep and subcycled fine timesteps.  The post_timestep function
// uses them to sync the fine and coarse levels.  If the coarse level
// is not the base level, post_timestep modifies the next coarsest levels
// registers appropriately.
//
// Note :: There is a little ambiguity as to which level owns the
// boundary flux registers.  The Multifab registers are quantities
// sized by the coarse level BoxArray and belong to the coarse level.
// The fine levels own the boundary registers, since they are sized by
// the boundaries of the fine level BoxArray.
//

//
// Compute a timestep at a level. Return largest safe timestep.
//




Real
NavierStokes::advance (Real time,
                       Real dt,
                       int  iteration,
                       int  ncycle)
{
    BL_PROFILE("NavierStokes::advance()");

    if (verbose)
    {
        Print() << "Advancing grids at level " << level
                << " : starting time = "       << time
                << " with dt = "               << dt
                << std::endl;
    }

    advance_setup(time,dt,iteration,ncycle);

    //
    // Calculate the time N viscosity and diffusivity
    //   Note: The viscosity and diffusivity at time N+1 are
    //         initialized here to the time N values just to
    //         have something reasonable.
    //
    const Real prev_time = state[State_Type].prevTime();
    const int num_diff = NUM_STATE-BL_SPACEDIM-1;

    calcViscosity(prev_time,dt,iteration,ncycle);
    calcDiffusivity(prev_time);
    MultiFab::Copy(*viscnp1_cc, *viscn_cc, 0, 0, 1, viscn_cc->nGrow());
    MultiFab::Copy(*diffnp1_cc, *diffn_cc, 0, 0, num_diff, diffn_cc->nGrow());

    // Add this AFTER advance_setup()
    if (verbose)
    {
        Print() << "NavierStokes::advance(): before velocity update:"
                << std::endl;
        printMaxValues(false);
    }

    //
    // Compute traced states for normal comp of velocity at half time level.
    //
    Real dt_test = predict_velocity(dt);
    //
    // Do MAC projection and update edge velocities.
    //
    if (do_mac_proj)
    {
	// FIXME? rhs composed of divu and dSdt terms, which are FillPatch'ed
	// from the stored state
	// orig IAMR ng=0. mfix uses ng=4. Create NSBase variable???
	//
#ifdef AMREX_USE_EB
	int ng_rhs = 4;
#else
	int ng_rhs = 0;
#endif
	MultiFab mac_rhs(grids,dmap,1,ng_rhs,MFInfo(),Factory());
	create_mac_rhs(mac_rhs,ng_rhs,time,dt);
        MultiFab& S_old = get_old_data(State_Type);
	// NOTE have_divu is now a static var in NSBase
        mac_project(time,dt,S_old,&mac_rhs,umac_n_grow,true);
    } else {
	create_umac_grown(umac_n_grow);
    }
    //
    // Advect velocities.
    //
    if (do_mom_diff == 0)
        velocity_advection(dt);
    //
    // Advect scalars.
    //
    const int first_scalar = Density;
    const int last_scalar  = first_scalar + NUM_SCALARS - 1;
    scalar_advection(dt,first_scalar,last_scalar);
    //
    // Update Rho.
    //
    scalar_update(dt,first_scalar,first_scalar);
    make_rho_curr_time();
    //
    // Advect momenta after rho^(n+1) has been created.
    //
    if (do_mom_diff == 1)
        velocity_advection(dt);
    //
    // Add the advective and other terms to get scalars at t^{n+1}.
    //
    if (do_scalar_update_in_order)
    {
	for (int iComp=0; iComp<NUM_SCALARS-1; iComp++)
        {
	    int iScal = first_scalar+scalarUpdateOrder[iComp];
	    Print() << "... ... updating " << desc_lst[0].name(iScal) << '\n';
	    scalar_update(dt,iScal,iScal);
	}
    }
    else
    {
	scalar_update(dt,first_scalar+1,last_scalar);
    }
//
//     MultiFab& S_new = get_new_data(State_Type);
//     const Real* dx       = geom.CellSize();
//     MultiFab&   P_new    = get_new_data(Press_Type);
//     const Real  cur_time = state[State_Type].curTime();
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
//     for (MFIter snewmfi(S_new,true); snewmfi.isValid(); ++snewmfi)
//     {
//         const Box& vbx = snewmfi.tilebox();
//
//         FArrayBox& Sfab = S_new[snewmfi];
//         FArrayBox& Pfab = P_new[snewmfi];
//
//         Sfab.setVal<RunOn::Host>(0.0,snewmfi.growntilebox(),0,S_new.nComp());
//         Pfab.setVal<RunOn::Host>(0.0,snewmfi.grownnodaltilebox(-1,P_new.nGrow()));
//
//         RealBox    gridloc = RealBox(vbx,geom.CellSize(),geom.ProbLo());
//         const int* lo      = vbx.loVect();
//         const int* hi      = vbx.hiVect();
//         const int* s_lo    = Sfab.loVect();
//         const int* s_hi    = Sfab.hiVect();
//         const int* p_lo    = Pfab.loVect();
//         const int* p_hi    = Pfab.hiVect();
//
//         // FORT_INITDATA (&level,&cur_time,lo,hi,&ns,
//         //                Sfab.dataPtr(Xvel),
//         //                Sfab.dataPtr(BL_SPACEDIM),
//         //                ARLIM(s_lo), ARLIM(s_hi),
//         //                Pfab.dataPtr(),
//         //                ARLIM(p_lo), ARLIM(p_hi),
//         //                dx,gridloc.lo(),gridloc.hi() );
//         computeIntHeight(Sfab, vbx, gridloc, dx);
//
//     }

    //
    // S appears in rhs of the velocity update, so we better do it now.
    //
    if (have_divu)
    {
        calc_divu(time+dt,dt,get_new_data(Divu_Type));
        if (have_dsdt)
        {
            calc_dsdt(time,dt,get_new_data(Dsdt_Type));
            if (initial_step)
                MultiFab::Copy(get_old_data(Dsdt_Type),
                               get_new_data(Dsdt_Type),0,0,1,0);
        }
    }
    //
    // Add the advective and other terms to get velocity at t^{n+1}.
    //
    velocity_update(dt);

    //
    // Increment rho average.
    //
    if (!initial_step)
    {
        if (level > 0)
            incrRhoAvg((iteration==ncycle ? 0.5 : 1.0) / Real(ncycle));

        if (verbose)
        {
            Print() << "NavierStokes::advance(): before nodal projection " << std::endl;
            printMaxValues();
        }

        //
        // Do a level project to update the pressure and velocity fields.
        //
        if (projector)
            level_projector(dt,time,iteration);
        if (level > 0 && iteration == 1)
           p_avg.setVal(0);
    }

#ifdef AMREX_PARTICLES
    if (theNSPC() != 0 and NavierStokes::initial_iter != true)
    {
        theNSPC()->AdvectWithUmac(u_mac, level, dt);
    }
#endif
    //
    // Clean up after the predicted value at t^n+1.
    // Estimate new timestep from umac cfl.
    //
    advance_cleanup(iteration,ncycle);

    if (verbose)
    {
        Print() << "NavierStokes::advance(): after velocity update" << std::endl;
        printMaxValues();
    }

    return dt_test;  // Return estimate of best new timestep.
}

//
void
NavierStokes::computeIntHeight(FArrayBox& statein, const Box& bx,
    RealBox gridloc, const Real* dx)
{
    Array4<Real> const& stateArray = statein.array();
    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    const Real* xlo = gridloc.lo();
    const Real* xhi = gridloc.hi();
    Real x, y;
    for(int i = lo.x; i <= hi.x; i++)
    {
        for(int j = lo.y; j < hi.y; j++)
        {
            for(int k = lo.z; k <= hi.z; k++)
            {
                // x = dx[0]*(i+0.5);
                // y = dx[1]*(j+0.5);
                x = xlo[0] + dx[0]*((i-lo.x) + 0.5);
                y = xlo[1] + dx[1]*((j-lo.y) + 0.5);

                if(stateArray(i,j,k,Tracer) <= 0.5 && stateArray(i,j+1,k,Tracer) > 0.5)
                {
                    stateArray(lo.x,j,k,Tracer+1) = y;
                }
            }
        }
    }
}
//
// Predict the edge velocities which go into forming u_mac.  This
// function also returns an estimate of dt for use in variable timesteping.
//
void
NavierStokes::floor(MultiFab& mf){

  int ncomp = mf.nComp();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box gbx=mfi.growntilebox(godunov_hyp_grow);
        auto const& fab_a = mf.array(mfi);
        AMREX_PARALLEL_FOR_4D ( gbx, ncomp, i, j, k, n,
        {
            auto& val = fab_a(i,j,k,n);
            val = amrex::Math::abs(val) > 1.e-20 ? val : 0;
        });
    }
}


Real
NavierStokes::predict_velocity (Real  dt)
{
    BL_PROFILE("NavierStokes::predict_velocity()");

    if (verbose) Print() << "... predict edge velocities\n";
    //
    // Get simulation parameters.
    //
    const int   nComp          = BL_SPACEDIM;
    const Real* dx             = geom.CellSize();
    const Real  prev_time      = state[State_Type].prevTime();
    const Real  prev_pres_time = state[Press_Type].prevTime();
    //
    // Compute viscous terms at level n.
    // Ensure reasonable values in 1 grow cell.  Here, do extrap for
    // c-f/phys boundary, since we have no interpolator fn, also,
    // preserve extrap for corners at periodic/non-periodic intersections.
    //
    MultiFab visc_terms(grids,dmap,nComp,1,MFInfo(), Factory());
    if (be_cn_theta != 1.0)
    {
	  getViscTerms(visc_terms,Xvel,nComp,prev_time);
    }
    else
    {
	  visc_terms.setVal(0);
    }

    FillPatchIterator U_fpi(*this,visc_terms,godunov_hyp_grow,prev_time,State_Type,Xvel,BL_SPACEDIM);
    MultiFab& Umf=U_fpi.get_mf();

    // Floor small values of states to be extrapolated
    floor(Umf);

    FillPatchIterator S_fpi(*this,visc_terms,6,prev_time,State_Type,Density,NUM_SCALARS);
    MultiFab& Smf=S_fpi.get_mf();
    // MultiFab SAlias(Smf, amrex::make_alias,0,Smf.nComp());
    // const BoxArray& ba = Smf.boxArray();
    // const DistributionMapping& dm = Smf.DistributionMap();
    // MultiFab SBorder(ba, dm, Smf.nComp(),12);

    MultiFab scalGradFab(grids, dmap, 3*nTrac, 4);

    //
    // Compute "grid cfl number" based on cell-centered time-n velocities
    //
    auto umax = VectorMaxAbs({&Umf},FabArrayBase::mfiter_tile_size,0,BL_SPACEDIM,Umf.nGrow());

    Real cflmax = dt*umax[0]/dx[0];
    for (int d=1; d<BL_SPACEDIM; ++d)
    {
        cflmax = std::max(cflmax,dt*umax[d]/dx[d]);
    }
    Real tempdt = cflmax==0 ? change_max : std::min(change_max,cfl/cflmax);


#if AMREX_USE_EB

    Vector<BCRec> math_bcs(AMREX_SPACEDIM);
    math_bcs = fetchBCArray(State_Type,Xvel,AMREX_SPACEDIM);

    MOL::ExtrapVelToFaces( Umf,
                           D_DECL(u_mac[0], u_mac[1], u_mac[2]),
                           geom, math_bcs );
#else
    //
    // Non-EB version
    //
    const int ngrow = 1;
    MultiFab Gp(grids, dmap, AMREX_SPACEDIM,ngrow);
    getGradP(Gp, prev_pres_time);

    MultiFab forcing_term( grids, dmap, AMREX_SPACEDIM, ngrow );

    //
    // Compute forcing
    //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox tforces;
        Vector<int> bndry[BL_SPACEDIM];
        for(MFIter SGmfi(scalGradFab,true); SGmfi.isValid(); ++SGmfi)
        {
            const Box& bx = SGmfi.tilebox();
            computeScalarGradient(Smf[SGmfi], scalGradFab[SGmfi], bx, dx, Smf.nComp());
        }

        for (MFIter U_mfi(Umf,TilingIfNotGPU()); U_mfi.isValid(); ++U_mfi)
        {
            Box bx=U_mfi.tilebox();
            auto const gbx=U_mfi.growntilebox(ngrow);

            FArrayBox& Ufab = Umf[U_mfi];
            // tforces.resize(gbx,BL_SPACEDIM);
            if (getForceVerbose) {
                Print() << "---\nA - Predict velocity:\n Calling getForce...\n";
            }

	    // const Box& forcebx = grow(bx,1);
        // 	    tforces.resize(forcebx,AMREX_SPACEDIM);
        getForce(forcing_term[U_mfi],gbx,1,Xvel,BL_SPACEDIM,prev_time,
            Ufab,Smf[U_mfi],scalGradFab[U_mfi],0);

            //
            // Compute the total forcing.
            //
            auto const& tf   = forcing_term.array(U_mfi,Xvel);
            auto const& visc = visc_terms.const_array(U_mfi,Xvel);
            auto const& gp   = Gp.const_array(U_mfi);
            auto const& rho  = rho_ptime.const_array(U_mfi);

            amrex::ParallelFor(gbx, AMREX_SPACEDIM, [tf, visc, gp, rho]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                tf(i,j,k,n) = ( tf(i,j,k,n) + visc(i,j,k,n) - gp(i,j,k,n) ) / rho(i,j,k);
            });
        }
    }

    Vector<BCRec> math_bcs(AMREX_SPACEDIM);
    math_bcs = fetchBCArray(State_Type,Xvel,AMREX_SPACEDIM);

    //velpred=1 only, use_minion=1, ppm_type, slope_order
    Godunov::ExtrapVelToFaces( Umf, forcing_term, AMREX_D_DECL(u_mac[0], u_mac[1], u_mac[2]),
                               math_bcs, geom, dt, godunov_use_ppm,
                               godunov_use_forces_in_trans );

#endif

    return dt*tempdt;
}




// This routine advects the scalars
//

void
NavierStokes::scalar_advection (Real dt,
                                int  fscalar,
                                int  lscalar)
{
    BL_PROFILE("NavierStokes::scalar_advection()");

    if (verbose) Print() << "... advect scalars\n";
    //
    // Get simulation parameters.
    //
    const int   num_scalars    = lscalar - fscalar + 1;
    const Real* dx             = geom.CellSize();
    const Real  prev_time      = state[State_Type].prevTime();

    //
    // Get the viscous terms.
    //
    MultiFab visc_terms(grids,dmap,num_scalars,1,MFInfo(),Factory());

    if (be_cn_theta != 1.0)
    {
        getViscTerms(visc_terms,fscalar,num_scalars,prev_time);
    }
    else
    {
        visc_terms.setVal(0.0,1);
    }

    int nGrowF = 1;
    MultiFab* divu_fp = getDivCond(nGrowF,prev_time);
    MultiFab* dsdt    = getDsdt(nGrowF,prev_time);
    MultiFab::Saxpy(*divu_fp, 0.5*dt, *dsdt, 0, 0, 1, nGrowF);
    delete dsdt;

    MultiFab fluxes[BL_SPACEDIM];

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        const BoxArray& ba = getEdgeBoxArray(i);
        fluxes[i].define(ba, dmap, num_scalars, 0, MFInfo(), Factory());
    }

    //
    // Compute the advective forcing.
    //
    {
        FillPatchIterator S_fpi(*this,visc_terms,godunov_hyp_grow+6,prev_time,State_Type,fscalar,num_scalars);
        MultiFab& Smf=S_fpi.get_mf();

        // Floor small values of states to be extrapolated
	floor(Smf);
    MultiFab SBorder(grids, dmap, NUM_STATE,godunov_hyp_grow+3,MFInfo(), Factory());
    MultiFab scalGradFab(grids, dmap, 3*nTrac, godunov_hyp_grow+4);
    FillPatch(*this, SBorder, godunov_hyp_grow+3, prev_time, State_Type, 0, NUM_STATE);

        FillPatchIterator U_fpi(*this,visc_terms,godunov_hyp_grow,prev_time,State_Type,Xvel,BL_SPACEDIM);
        const MultiFab& Umf=U_fpi.get_mf();



#ifdef AMREX_USE_EB
        //////////////////////////////////////////////////////////////////////////////
        //  EB ALGORITHM
        //////////////////////////////////////////////////////////////////////////////

        const Box& domain = geom.Domain();

        Vector<BCRec> math_bc(num_scalars);
        math_bc = fetchBCArray(State_Type,fscalar,num_scalars);


        MultiFab cfluxes[AMREX_SPACEDIM];
        MultiFab edgstate[AMREX_SPACEDIM];
        int nghost = 2;

        for (int i(0); i < AMREX_SPACEDIM; i++)
        {
            const BoxArray& ba = getEdgeBoxArray(i);
            cfluxes[i].define(ba, dmap, num_scalars, nghost, MFInfo(), Factory());
            edgstate[i].define(ba, dmap, num_scalars, nghost, MFInfo(), Factory());
        }

        Vector<BCRec> math_bcs(num_scalars);
        math_bcs = fetchBCArray(State_Type, fscalar, num_scalars);

        MOL::ComputeAofs(*aofs, fscalar, num_scalars, Smf, 0,
                         D_DECL(u_mac[0],u_mac[1],u_mac[2]),
                         D_DECL(edgstate[0],edgstate[1],edgstate[2]), 0, false,
                         D_DECL(cfluxes[0],cfluxes[1],cfluxes[2]), 0,
                         math_bcs, geom  );

        if (do_reflux)
        {
            for (int d(0); d < AMREX_SPACEDIM; d++)
                MultiFab::Copy(fluxes[d], cfluxes[d], 0, 0, num_scalars, 0 );

        }

#else
        //////////////////////////////////////////////////////////////////////////////
        //  NON-EB ALGORITHM
        //////////////////////////////////////////////////////////////////////////////
        MultiFab cfluxes[AMREX_SPACEDIM];
        MultiFab edgestate[AMREX_SPACEDIM];
	MultiFab forcing_term( grids, dmap, num_scalars, nGrowF );

        // NO Gghost nodes???
        int nghost = 0;
        for (int i(0); i < AMREX_SPACEDIM; i++)
        {
            const BoxArray& ba = getEdgeBoxArray(i);
            cfluxes[i].define(ba, dmap, num_scalars, nghost);
            cfluxes[i].setVal(0.0);
            edgestate[i].define(ba, dmap, num_scalars, nghost);
        }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            for(MFIter SGmfi(scalGradFab,true); SGmfi.isValid(); ++SGmfi)
            {
                const Box& bx = SGmfi.tilebox();
                computeScalarGradient(Smf[SGmfi], scalGradFab[SGmfi], bx, dx, Smf.nComp());
            }
            Vector<int> state_bc;
            FArrayBox tforces;

            for (MFIter S_mfi(Smf,TilingIfNotGPU()); S_mfi.isValid(); ++S_mfi)
            {
                if (getForceVerbose)
                {
                    Print() << "---" << '\n' << "C - scalar advection:" << '\n'
                            << " Calling getForce..." << '\n';
                }
                const Box& gbx = S_mfi.growntilebox(nGrowF);
                		// tforces.resize(forcebx,num_scalars);
        // getForce(tforces,bx,nGrowF,fscalar,num_scalars,prev_time,Umf[S_mfi],Smf[S_mfi], scalGradFab[S_mfi],0);

                getForce(forcing_term[S_mfi],gbx,nGrowF,fscalar,num_scalars,
			 prev_time,Umf[S_mfi],Smf[S_mfi],scalGradFab[S_mfi],0);

                for (int n=0; n<num_scalars; ++n)
                {
                    // FIXME: Loop rqd b/c function does not take array conserv_diff
		    auto const& tf    = forcing_term.array(S_mfi,n);
                    auto const& visc  = visc_terms.const_array(S_mfi,n);


                    if (advectionType[fscalar+n] == Conservative)
                    {
                        auto const& divu  = divu_fp -> const_array(S_mfi);
                        auto const& S     = Smf.array(S_mfi);

                        amrex::ParallelFor(gbx, [tf, visc, S, divu]
                        AMREX_GPU_DEVICE (int i, int j, int k ) noexcept
                        {
			    tf(i,j,k) += visc(i,j,k) - S(i,j,k) * divu(i,j,k);
                        });
		    }
                    else
                    {
                        auto const& rho   = rho_ptime.const_array(S_mfi);

                        amrex::ParallelFor(gbx, [tf, visc, rho]
                        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            tf(i,j,k) = ( tf(i,j,k) + visc(i,j,k) ) / rho(i,j,k);
                        });
		    }

                }
            }
        }


        Vector<BCRec> math_bcs(num_scalars);
        math_bcs = fetchBCArray(State_Type, fscalar, num_scalars);

        amrex::Gpu::DeviceVector<int> iconserv;
        iconserv.resize(num_scalars, 0);
        for (int comp = 0; comp < num_scalars; ++comp)
        {
            iconserv[comp] = (advectionType[fscalar+comp] == Conservative) ? 1 : 0;
	}

        Godunov::ComputeAofs(*aofs, fscalar, num_scalars,
                             Smf, 0,
                             AMREX_D_DECL( u_mac[0], u_mac[1], u_mac[2] ),
                             AMREX_D_DECL( edgestate[0], edgestate[1], edgestate[2] ), 0, false,
                             AMREX_D_DECL( cfluxes[0], cfluxes[1], cfluxes[2] ), 0,
                             forcing_term, 0, *divu_fp, math_bcs, geom, iconserv,
                             dt, godunov_use_ppm, godunov_use_forces_in_trans, false );

        if (do_reflux)
        {
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                MultiFab::Copy(fluxes[d], cfluxes[d], 0, 0, num_scalars, 0 );
	}
#endif
    } // FillPathIterator

    delete divu_fp;

    if (do_reflux)
    {
        if (level > 0 )
        {
            for (int d = 0; d < BL_SPACEDIM; d++)
                advflux_reg->FineAdd(fluxes[d],d,0,fscalar,num_scalars,dt);
        }
        if (level < parent->finestLevel())
        {
            for (int i = 0; i < BL_SPACEDIM; i++)
                getAdvFluxReg(level+1).CrseInit(fluxes[i],i,0,fscalar,num_scalars,-dt);
        }
    }
}

//
// This subroutine updates the scalars, before the velocity update
// and the level projection
//
// AT this point in time, all we know is psi^n, rho^n+1/2, and the
// general forcing terms at t^n, and after solving in this routine
// viscous forcing at t^n+1/2.  Note, unless more complicated logic
// is invoked earlier, we do not have any estimate of general forcing
// terms at t^n+1/2.
//

void
NavierStokes::scalar_update (Real dt,
                             int  first_scalar,
                             int  last_scalar)
{
    BL_PROFILE("NavierStokes::scalar_update()");

    if (verbose) Print() << "... update scalars\n";

    scalar_advection_update(dt, first_scalar, last_scalar);

    bool do_any_diffuse = false;
    for (int sigma = first_scalar; sigma <= last_scalar; sigma++)
        if (is_diffusive[sigma]) do_any_diffuse = true;

    if (do_any_diffuse)
      scalar_diffusion_update(dt, first_scalar, last_scalar);

    MultiFab&  S_new     = get_new_data(State_Type);
//#ifdef AMREX_USE_EB
//  set_body_state(S_new);
//#endif
    for (int sigma = first_scalar; sigma <= last_scalar; sigma++)
    {
       if (S_new.contains_nan(sigma,1,0))
       {
	 Print() << "New scalar " << sigma << " contains Nans" << '\n';
	 exit(0);
       }
    }
}

void
NavierStokes::scalar_diffusion_update (Real dt,
                                       int  first_scalar,
                                       int  last_scalar)
{
    BL_PROFILE("NavierStokes::scalar_diffusion_update()");

    const MultiFab& Rh = get_rho_half_time();

    int ng=1;
    const Real prev_time = state[State_Type].prevTime();
    const Real curr_time = state[State_Type].curTime();

    //fixme? why fillpatch all of state when only doing scalars?
    FillPatch(*this,get_old_data(State_Type),ng,prev_time,State_Type,0,NUM_STATE);
    FillPatch(*this,get_new_data(State_Type),ng,curr_time,State_Type,0,NUM_STATE);

    auto Snc = std::unique_ptr<MultiFab>(new MultiFab());
    auto Snp1c = std::unique_ptr<MultiFab>(new MultiFab());

    if (level > 0) {
      auto& crselev = getLevel(level-1);
      Snc->define(crselev.boxArray(), crselev.DistributionMap(), NUM_STATE, ng, MFInfo(), crselev.Factory());
      FillPatch(crselev,*Snc  ,ng,prev_time,State_Type,0,NUM_STATE);

      Snp1c->define(crselev.boxArray(), crselev.DistributionMap(), NUM_STATE, ng, MFInfo(), crselev.Factory());
      FillPatch(crselev,*Snp1c,ng,curr_time,State_Type,0,NUM_STATE);
    }

    const int nlev = (level ==0 ? 1 : 2);
    Vector<MultiFab*> Sn(nlev,0), Snp1(nlev,0);
    Sn[0]   = &(get_old_data(State_Type));
    Snp1[0] = &(get_new_data(State_Type));

    if (nlev>1) {
      Sn[1]   =  Snc.get() ;
      Snp1[1] =  Snp1c.get() ;
    }

    const Vector<BCRec>& theBCs = AmrLevel::desc_lst[State_Type].getBCs();

    FluxBoxes fb_diffn, fb_diffnp1;
    MultiFab **cmp_diffn = 0, **cmp_diffnp1 = 0;

    MultiFab *delta_rhs = 0;
    MultiFab *alpha = 0;
    const int rhsComp = 0, alphaComp = 0, fluxComp  = 0;

    FluxBoxes fb_fluxn  (this);
    FluxBoxes fb_fluxnp1(this);
    MultiFab** fluxn   = fb_fluxn.get();
    MultiFab** fluxnp1 = fb_fluxnp1.get();

    Vector<int> diffuse_comp(1);

    for (int sigma = first_scalar; sigma <= last_scalar; sigma++)
    {
      if (verbose)
	Print()<<"scalar_diffusion_update "<<sigma<<" of "<<last_scalar<<"\n";

      if (is_diffusive[sigma])
      {
        if (be_cn_theta != 1)
        {
          cmp_diffn = fb_diffn.define(this);
          getDiffusivity(cmp_diffn, prev_time, sigma, 0, 1);
        }

        cmp_diffnp1 = fb_diffnp1.define(this);
        getDiffusivity(cmp_diffnp1, curr_time, sigma, 0, 1);

        diffuse_comp[0] = is_diffusive[sigma];
        const int rho_flag = Diffusion::set_rho_flag(diffusionType[sigma]);

        const Diffusion::SolveMode& solve_mode = Diffusion::ONEPASS;
        const bool add_old_time_divFlux = true;

        const int betaComp = 0;
        const int Rho_comp = Density;
	const int bc_comp  = sigma;

        diffusion->diffuse_scalar (Sn, Sn, Snp1, Snp1, sigma, 1, Rho_comp,
                                   prev_time,curr_time,be_cn_theta,Rh,rho_flag,
                                   fluxn,fluxnp1,fluxComp,delta_rhs,rhsComp,
                                   alpha,alphaComp,
                                   cmp_diffn,cmp_diffnp1,betaComp,
                                   crse_ratio,theBCs[bc_comp],geom,
                                   solve_mode,add_old_time_divFlux,
                                   diffuse_comp);

        if(alpha!=0) delete alpha;

        //
        // Increment the viscous flux registers
        //
        if (do_reflux)
        {

	  FArrayBox fluxtot;
	  for (int d = 0; d < BL_SPACEDIM; d++)
          {
	    MultiFab fluxes;

	    if (level < parent->finestLevel())
	    {
	      fluxes.define(fluxn[d]->boxArray(), fluxn[d]->DistributionMap(), 1, 0, MFInfo(), Factory());
	    }

	    for (MFIter fmfi(*fluxn[d]); fmfi.isValid(); ++fmfi)
	    {
	      const Box& ebox = (*fluxn[d])[fmfi].box();//fmfi.tilebox();

	      fluxtot.resize(ebox,1);
	      fluxtot.copy<RunOn::Host>((*fluxn[d])[fmfi],ebox,0,ebox,0,1);
	      fluxtot.plus<RunOn::Host>((*fluxnp1[d])[fmfi],ebox,0,0,1);

	      if (level < parent->finestLevel())
		fluxes[fmfi].copy<RunOn::Host>(fluxtot);

	      if (level > 0)
		getViscFluxReg().FineAdd(fluxtot,d,fmfi.index(),0,sigma,1,dt,RunOn::Host);
	    }

	    if (level < parent->finestLevel())
	      getLevel(level+1).getViscFluxReg().CrseInit(fluxes,d,0,sigma,1,-dt);

	  }
	}

	if (be_cn_theta != 1)
	  fb_diffn.clear();
	fb_diffnp1.clear();

      }//end if(is_diffusive)
    }
}

void
NavierStokes::velocity_diffusion_update (Real dt)
{
    BL_PROFILE("NavierStokes::velocity_diffusion_update()");

    const Real strt_time = ParallelDescriptor::second();
    //
    // Compute the viscous forcing.
    // Do following except at initial iteration.
    //
    if (is_diffusive[Xvel])
    {
        int rho_flag = (do_mom_diff == 0) ? 1 : 3;

        MultiFab* delta_rhs = 0;
        if (S_in_vel_diffusion && have_divu)
        {
            delta_rhs = new MultiFab(grids,dmap,BL_SPACEDIM,0, MFInfo(),Factory());
            delta_rhs->setVal(0);
        }

	FluxBoxes fb_viscn, fb_viscnp1;
        MultiFab** loc_viscn   = 0;
        MultiFab** loc_viscnp1 = 0;

        Real viscTime = state[State_Type].prevTime();
        loc_viscn = fb_viscn.define(this);
        getViscosity(loc_viscn, viscTime);

        viscTime = state[State_Type].curTime();
        loc_viscnp1 = fb_viscnp1.define(this);
        getViscosity(loc_viscnp1, viscTime);

        diffuse_velocity_setup(dt, delta_rhs, loc_viscn, loc_viscnp1);

        diffusion->diffuse_velocity(dt,be_cn_theta,get_rho_half_time(),rho_flag,
                                    delta_rhs,loc_viscn,viscn_cc,loc_viscnp1,viscnp1_cc);

        delete delta_rhs;
    }

    if (verbose)
    {
        Real run_time    = ParallelDescriptor::second() - strt_time;
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	Print() << "NavierStokes:velocity_diffusion_update(): lev: " << level
		       << ", time: " << run_time << '\n';
    }
}

void
NavierStokes::diffuse_velocity_setup (Real       dt,
                                      MultiFab*& delta_rhs,
                                      MultiFab**& viscn,
                                      MultiFab**& viscnp1)
{
    if (S_in_vel_diffusion && have_divu)
    {
        //
        // Include div mu S*I terms in rhs
        //  (i.e. make nonzero delta_rhs to add into RHS):
        //
        // The scalar and tensor solvers incorporate the relevant pieces of
        //  of Div(tau), provided the flow is divergence-free.  However, if
        //  Div(U) =/= 0, there is an additional piece not accounted for,
        //  which is of the form A.Div(U).
        //
        // Now we only use the tensor solver.
        // For history, before for constant viscosity, Div(tau)_i
        //  = Lapacian(U_i) + mu/3 d[Div(U)]/dx_i.
        // Now because  mu not constant,
        //  Div(tau)_i = d[ mu(du_i/dx_j + du_j/dx_i) ]/dx_i - 2mu/3 d[Div(U)]/dx_i
        //
        // As a convenience, we treat this additional term as a "source" in
        // the diffusive solve, computing Div(U) in the "normal" way we
        // always do--via a call to calc_divu.  This routine computes delta_rhs
        // if necessary, and stores it as an auxilliary rhs to the viscous solves.
        // This is a little strange, but probably not bad.
        //
        const Real time = state[State_Type].prevTime();

        MultiFab divmusi(grids,dmap,BL_SPACEDIM,0,MFInfo(),Factory());

        diffusion->compute_divmusi(time,viscn,divmusi);
        divmusi.mult((-2./3.)*(1.0-be_cn_theta),0,BL_SPACEDIM,0);
                      (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);

        diffusion->compute_divmusi(time+dt,viscnp1,divmusi);
        divmusi.mult((-2./3.)*be_cn_theta,0,BL_SPACEDIM,0);
                (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);
    }
}

//fixme? is there now an amrex fn for this?
Real
NavierStokes::MaxVal (const std::string& name,
                      Real           time)
{
    Real        mxval = 0.0;
    auto        mf = derive(name,time,0);
    BoxArray    baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    //Add and test this OMP
    //#ifdef _OPENMP
    //#pragma omp parallel if (!system::regtest_reduction) reduction(max:mxval,s)
    //#endif
    //{

    std::vector< std::pair<int,Box> > isects;

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        const int  i   = mfi.index();
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            baf.intersections(grids[i],isects);

            for (int ii = 0, N = isects.size(); ii < N; ii++)
              fab.setVal<RunOn::Host>(0,isects[ii].second,0,fab.nComp());
        }
        Real        s;
        const Real* dat = fab.dataPtr();
        const int*  dlo = fab.loVect();
        const int*  dhi = fab.hiVect();
	const Box&  bx  = grids[i];
        const int*  lo  = bx.loVect();
        const int*  hi  = bx.hiVect();

        fort_maxval(dat,ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),&s);

        mxval = std::max(mxval, s);
    }
    //} end OMP parallel

    ParallelDescriptor::ReduceRealMax(mxval);

    return mxval;
}

void
NavierStokes::sum_integrated_quantities ()
{
    const int finest_level = parent->finestLevel();

    Real time = state[State_Type].curTime();
    Real mass = 0.0;
    Real trac = 0.0;
    Real energy = 0.0;
    Real mgvort = 0.0;
#if defined(DO_IAMR_FORCE)
    Real forcing = 0.0;
#endif
#if (BL_SPACEDIM==3)
    Real udotlapu = 0.0;
#endif

    for (int lev = 0; lev <= finest_level; lev++)
    {
        NavierStokes& ns_level = getLevel(lev);
	mass += ns_level.volWgtSum("density",time);
	trac += ns_level.volWgtSum("tracer",time);
        energy += ns_level.volWgtSum("energy",time);
        mgvort = std::max(mgvort,ns_level.MaxVal("mag_vort",time));
#if defined(DO_IAMR_FORCE)
        forcing += ns_level.volWgtSum("forcing",time);
#endif
#if (BL_SPACEDIM==3)
        udotlapu += ns_level.volWgtSum("udotlapu",time);
#endif
    }

    Print() << '\n';
    Print().SetPrecision(12) << "TIME= " << time << " MASS= " << mass << '\n';
    Print().SetPrecision(12) << "TIME= " << time << " TRAC= " << trac << '\n';
    Print().SetPrecision(12) << "TIME= " << time << " KENG= " << energy << '\n';
    Print().SetPrecision(12) << "TIME= " << time << " MAGVORT= " << mgvort << '\n';
    Print().SetPrecision(12) << "TIME= " << time << " ENERGY= " << energy << '\n';
#if defined(DO_IAMR_FORCE)
    //NOTE: FORCING_T gives only the energy being injected by the forcing
    //      term used for generating turbulence in probtype 14, 15.
    //      Defaults to 0 for other probtypes.
    Print().SetPrecision(12) << "TIME= " << time << " FORCING_T= " << forcing << '\n';
#endif
#if (BL_SPACEDIM==3)
    Print().SetPrecision(12) << "TIME= " << time << " UDOTLAPU= " << udotlapu << '\n';
#endif
}

void
NavierStokes::setPlotVariables()
{
    AmrLevel::setPlotVariables();
}

void
NavierStokes::writePlotFile (const std::string& dir,
                             std::ostream&  os,
                             VisMF::How     how)
{
    if ( ! Amr::Plot_Files_Output() ) return;

    BL_PROFILE("NavierStokes::writePlotFile()");

    int i, n;
    //
    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    //
    std::vector<std::pair<int,int> > plot_var_map;
    for (int typ = 0; typ < desc_lst.size(); typ++)
        for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
            if (parent->isStatePlotVar(desc_lst[typ].name(comp)) &&
                desc_lst[typ].getType() == IndexType::TheCellType())
                plot_var_map.push_back(std::pair<int,int>(typ,comp));

    int num_derive = 0;
    std::list<std::string> derive_names;
    const std::list<DeriveRec>& dlist = derive_lst.dlist();

    for (std::list<DeriveRec>::const_iterator it = dlist.begin(), end = dlist.end();
         it != end;
         ++it)
    {
        if (parent->isDerivePlotVar(it->name()))
	{
            derive_names.push_back(it->name());
            num_derive += it->numDerive();
	}
    }

    int n_data_items = plot_var_map.size() + num_derive;
#ifdef AMREX_USE_EB
    // add in vol frac
    n_data_items++;
#endif
    Real cur_time = state[State_Type].curTime();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0)
            Error("Must specify at least one valid data item to plot");

        os << n_data_items << '\n';

	//
	// Names of variables -- first state, then derived
	//
	for (i =0; i < plot_var_map.size(); i++)
        {
	    int typ  = plot_var_map[i].first;
	    int comp = plot_var_map[i].second;
	    os << desc_lst[typ].name(comp) << '\n';
        }

	for (std::list<std::string>::const_iterator it = derive_names.begin(), end = derive_names.end();
             it != end;
             ++it)
        {
	    const DeriveRec* rec = derive_lst.get(*it);
	    for (i = 0; i < rec->numDerive(); i++)
                os << rec->variableName(i) << '\n';
        }
#ifdef AMREX_USE_EB
	//add in vol frac
	os << "volFrac\n";
#endif

        os << BL_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geom().ProbLo(i) << ' ';
        os << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geom().ProbHi(i) << ' ';
        os << '\n';
        for (i = 0; i < f_lev; i++)
            os << parent->refRatio(i)[0] << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->Geom(i).Domain() << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->levelSteps(i) << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++)
                os << parent->Geom(i).CellSize()[k] << ' ';
            os << '\n';
        }
        os << (int) Geom().Coord() << '\n';
        os << "0\n"; // Write bndry data.


        // job_info file with details about the run
	std::ofstream jobInfoFile;
	std::string FullPathJobInfoFile = dir;
	FullPathJobInfoFile += "/job_info";
	jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

	std::string PrettyLine = "===============================================================================\n";
	std::string OtherLine = "--------------------------------------------------------------------------------\n";
	std::string SkipSpace = "        ";


	// job information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Job Information\n";
	jobInfoFile << PrettyLine;

	jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << '\n';
#ifdef _OPENMP
	jobInfoFile << "number of threads:       " << omp_get_max_threads() << '\n';
#endif
	jobInfoFile << "\n\n";

        // plotfile information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Plotfile Information\n";
	jobInfoFile << PrettyLine;

	time_t now = time(0);

	// Convert now to tm struct for local timezone
	tm* localtm = localtime(&now);
	jobInfoFile   << "output data / time: " << asctime(localtm);

	char currentDir[FILENAME_MAX];
	if (getcwd(currentDir, FILENAME_MAX)) {
	  jobInfoFile << "output dir:         " << currentDir << '\n';
	}

	jobInfoFile << "\n\n";


        // build information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Build Information\n";
	jobInfoFile << PrettyLine;

	jobInfoFile << "build date:    " << buildInfoGetBuildDate() << '\n';
	jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << '\n';
	jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << '\n';
	jobInfoFile << "BoxLib dir:    " << buildInfoGetAMReXDir() << '\n';

        jobInfoFile << '\n';

        jobInfoFile << "COMP:          " << buildInfoGetComp() << '\n';
	jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";
        jobInfoFile << "FCOMP:         " << buildInfoGetFcomp() << '\n';
	jobInfoFile << "FCOMP version: " << buildInfoGetFcompVersion() << "\n";

	jobInfoFile << "\n";

	const char* githash1 = buildInfoGetGitHash(1);
	const char* githash2 = buildInfoGetGitHash(2);
	if (strlen(githash1) > 0) {
	  jobInfoFile << "IAMR   git hash: " << githash1 << "\n";
	}
	if (strlen(githash2) > 0) {
	  jobInfoFile << "BoxLib git hash: " << githash2 << "\n";
	}

	jobInfoFile << "\n\n";


	// runtime parameters
	jobInfoFile << PrettyLine;
	jobInfoFile << " Inputs File Parameters\n";
	jobInfoFile << PrettyLine;

	ParmParse::dumpTable(jobInfoFile, true);

	jobInfoFile.close();

    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";

    std::string LevelStr = Concatenate("Level_", level, 1);
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
        FullPath += '/';
    FullPath += LevelStr;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!UtilCreateDirectory(FullPath, 0755))
            CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
            for (n = 0; n < BL_SPACEDIM; n++)
                os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (n_data_items > 0)
        {
            std::string PathNameInHeader = LevelStr;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }

#ifdef AMREX_USE_EB
	// volfrac threshhold for amrvis
	// fixme? pulled directly from CNS, might need adjustment
        if (level == parent->finestLevel()) {
            for (int lev = 0; lev <= parent->finestLevel(); ++lev) {
                os << "1.0e-6\n";
            }
        }
#endif
    }
    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    // NOTE: we are assuming that each state variable has one component,
    // but a derived variable is allowed to have multiple components.
    int       cnt   = 0;
    int       ncomp = 1;
    const int nGrow = 0;
    MultiFab  plotMF(grids,dmap,n_data_items,nGrow,MFInfo(),Factory());
    MultiFab* this_dat = 0;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (i = 0; i < plot_var_map.size(); i++)
    {
	int typ  = plot_var_map[i].first;
	int comp = plot_var_map[i].second;
	this_dat = &state[typ].newData();
	MultiFab::Copy(plotMF,*this_dat,comp,cnt,ncomp,nGrow);
	cnt+= ncomp;
    }
    //
    // Cull data from derived variables.
    //
    Real plot_time;

    if (derive_names.size() > 0)
    {
	for (std::list<std::string>::const_iterator it = derive_names.begin(), end = derive_names.end();
             it != end;
             ++it)
	{
            if (*it == "avg_pressure" ||
                *it == "gradpx"       ||
                *it == "gradpy"       ||
                *it == "gradpz")
            {
                if (state[Press_Type].descriptor()->timeType() ==
                    StateDescriptor::Interval)
                {
                    plot_time = cur_time;
                }
                else
                {
                    int f_lev = parent->finestLevel();
                    plot_time = getLevel(f_lev).state[Press_Type].curTime();
                }
            }
            else
            {
                plot_time = cur_time;
            }
	    const DeriveRec* rec = derive_lst.get(*it);
	    ncomp = rec->numDerive();
	    auto derive_dat = derive(*it,plot_time,nGrow);
	    MultiFab::Copy(plotMF,*derive_dat,0,cnt,ncomp,nGrow);
	    cnt += ncomp;
	}
    }

#ifdef AMREX_USE_EB
    // add volume fraction to plotfile
    plotMF.setVal(0.0, cnt, 1, nGrow);
    MultiFab::Copy(plotMF,*volfrac,0,cnt,1,nGrow);

    // set covered values for ease of viewing
    EB_set_covered(plotMF, 0.0);
#endif

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(plotMF,TheFullPath,how,true);
}

std::unique_ptr<MultiFab>
NavierStokes::derive (const std::string& name,
                      Real               time,
                      int                ngrow)
{
#ifdef AMREX_PARTICLES
    return ParticleDerive(name, time, ngrow);
#else
    return AmrLevel::derive(name, time, ngrow);
#endif
}

void
NavierStokes::derive (const std::string& name,
                      Real               time,
                      MultiFab&          mf,
                      int                dcomp)
{
#ifdef AMREX_PARTICLES
        ParticleDerive(name,time,mf,dcomp);
#else
        AmrLevel::derive(name,time,mf,dcomp);
#endif
}

//
// Ensure state, and pressure are consistent.
//
void
NavierStokes::post_init (Real stop_time)
{
    if (level > 0)
        //
        // Nothing to sync up at level > 0.
        //
        return;

    const int   finest_level = parent->finestLevel();
    Real        dt_init      = 0.0;
    Vector<Real> dt_save(finest_level+1);
    Vector<int>  nc_save(finest_level+1);
    //
    // Ensure state is consistent, i.e. velocity field is non-divergent,
    // Coarse levels are fine level averages, pressure is zero.
    //
    post_init_state();
    //
    // Estimate the initial timestepping.
    //
    post_init_estDT(dt_init, nc_save, dt_save, stop_time);
    //
    // Initialize the pressure by iterating the initial timestep.
    //
    post_init_press(dt_init, nc_save, dt_save);
    //
    // Compute the initial estimate of conservation.
    //
    if (sum_interval > 0)
        sum_integrated_quantities();
#if (BL_SPACEDIM==3)
    //
    // Derive turbulent statistics
    //
    if (turb_interval > 0)
        sum_turbulent_quantities();
#ifdef SUMJET
    //
    // Derive turbulent statistics for the round jet
    //
    if (jet_interval > 0)
        sum_jet_quantities();
#endif
#endif
}

//
// Initialize the pressure by iterating the initial timestep
//

void
NavierStokes::post_init_press (Real&        dt_init,
                               Vector<int>&  nc_save,
                               Vector<Real>& dt_save)
{
    if ( init_iter <= 0 ){
      // make sure there's not NANs in old pressure field
      // end up with P_old = P_new as is the case when doing initial iters
      MultiFab& p_old=get_old_data(Press_Type);
      MultiFab& p_new=get_new_data(Press_Type);
      MultiFab::Copy(p_old, p_new, 0, 0, 1, p_new.nGrow());

      parent->setDtLevel(dt_save);
      parent->setNCycle(nc_save);

      NavierStokes::initial_step = false;

      Print()<< "WARNING! post_init_press(): exiting without doing inital iterations because init_iter == "<<init_iter<<std::endl;

      return;
    }

    const Real strt_time       = state[State_Type].curTime();
    const int  finest_level    = parent->finestLevel();
    NavierStokes::initial_iter = true;

    if (verbose)
    {
        Print() << std::endl
                << "post_init_press(): "
                << "doing initial pressure iterations with dt = "
                << dt_init
                << std::endl;
    }

    //
    // Iterate over the advance function.
    //
    for (int iter = 0; iter < init_iter; iter++)
    {

        if (verbose)
        {
            Print() << std::endl
                    << "post_init_press(): iter = " << iter
                    << std::endl;
        }

        for (int k = 0; k <= finest_level; k++ )
        {
            getLevel(k).advance(strt_time,dt_init,1,1);
        }

        //
        // This constructs a guess at P, also sets p_old == p_new.
        //
        Vector<MultiFab*> sig(finest_level+1, nullptr);

        for (int k = 0; k <= finest_level; k++)
        {
            sig[k] = &(getLevel(k).get_rho_half_time());
        }

        if (projector)
        {
            projector->initialSyncProject(0,sig,dt_init,
                                          strt_time,have_divu);
        }

        for (int k = finest_level-1; k >= 0; k--)
        {
            getLevel(k).avgDown();
        }

        if (verbose)
        {
            // initSyncProject project d(u)/dt, so new velocity
            // is actually the projected acceleration
            // We don't actually care because initial velocity state will be
            // recovered at the end of each iteration.
            // However, we need to recover u_new from d(u)/dt if we want to print
            // correct diagnostics
            MultiFab& S_new = get_new_data(State_Type);
            MultiFab& S_old = get_old_data(State_Type);
            MultiFab::Xpay(S_new, dt_init, S_old, Xvel, Xvel, AMREX_SPACEDIM, 0);

            Print() << "After nodal projection:" << std::endl;
            printMaxValues();
        }

        for (int k = 0; k <= finest_level; k++)
        {
            //
            // Reset state variables to initial time, but
            // do not reset pressure variable, only pressure time.
            //
            getLevel(k).resetState(strt_time, dt_init, dt_init);
        }

        make_rho_curr_time();

        NavierStokes::initial_iter = false;
    }

    NavierStokes::initial_step = false;
    //
    // Re-instate timestep.
    //
    for (int k = 0; k <= finest_level; k++)
    {
        getLevel(k).setTimeLevel(strt_time,dt_save[k],dt_save[k]);
        if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
        {
            getLevel(k).state[Press_Type].setNewTimeLevel(.5*dt_init);
            getLevel(k).get_old_data(Dpdt_Type).setVal(0);
        }
    }

    parent->setDtLevel(dt_save);
    parent->setNCycle(nc_save);

    // Add space to output if verbose
    if (verbose)
    {
        Print() << std::endl
                << "post_init_press(): exiting after " << init_iter << " iterations"
                << std::endl
                << "After initial iterations: "
                << std::endl;
        printMaxValues();
        Print() << std::endl << std::endl;
    }

}

//
// The Mac Sync correction function
//
void
NavierStokes::mac_sync ()
{
    BL_PROFILE_REGION_START("R::NavierStokes::mac_sync()");
    BL_PROFILE("NavierStokes::mac_sync()");

    if (!do_reflux) return;

    if (verbose)
    {
        Print() << std::endl
                << "mac_sync() on level "<<level
                << std::endl;
    }

    const int  numscal        = NUM_STATE - BL_SPACEDIM;
    const Real prev_time      = state[State_Type].prevTime();
    const Real prev_pres_time = state[Press_Type].prevTime();
    const Real dt             = parent->dtLevel(level);
    MultiFab*  DeltaSsync     = 0;// hold (Delta rho)*q for conserved quantities
    // does this have ghosts filled?
    MultiFab&  Rh             = get_rho_half_time();

#ifdef AMREX_USE_EB
    // fixme? unsure how many ghost cells...
    // for umac, inflo uses: use_godunov ? 4 : 3;
    // for now, match umac which uses 4
    const int nghost = umac_n_grow; // ==4; For redistribution ... We may not need 4 but for now we play safe
#else
    const int nghost = 0;
#endif


    Array<MultiFab*,AMREX_SPACEDIM> Ucorr;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
      const BoxArray& edgeba = getEdgeBoxArray(idim);

      Ucorr[idim]= new MultiFab(edgeba,dmap,1,nghost,MFInfo(),Factory());
    }

    sync_setup(DeltaSsync);
    //
    // Compute the u_mac for the correction.
    //
    Vector<BCRec> rho_math_bc = fetchBCArray(State_Type,Density,1);
    mac_projector->mac_sync_solve(level,dt,Rh,rho_math_bc[0],fine_ratio,Ucorr);
    //
    // Update coarse grid state by adding correction from mac_sync solve
    // the correction is the advective tendency of the new velocities.
    //
    MultiFab& S_new = get_new_data(State_Type);
    mac_projector->mac_sync_compute(level,Ucorr,u_mac,Vsync,Ssync,Rh,
				    level > 0 ? &getAdvFluxReg(level) : 0,
				    advectionType, prev_time,
				    prev_pres_time,dt,
				    NUM_STATE,be_cn_theta,
				    modify_reflux_normal_vel,
				    do_mom_diff);
    //
    // Delete Ucorr; we're done with it.
    //
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
      delete Ucorr[idim];


    //
    // For all conservative variables Q (other than density)
    // express Q as rho*q and increment sync by -(sync_for_rho)*q
    // (See Pember, et. al., LBNL-41339, Jan. 1989)
    //
    int iconserved = -1;
    for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
    {
      if (istate != Density && advectionType[istate] == Conservative)
      {
	iconserved++;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
	for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
	    const Box&  bx       = mfi.tilebox();
	    auto const& rho      = S_new.array(mfi,Density);
	    auto const& Snew     = S_new.array(mfi,istate);
	    auto const& dSsync   = DeltaSsync->array(mfi);
	    auto const& drhosync = Ssync.array(mfi,Density-AMREX_SPACEDIM);
	    auto const& ssync    = Ssync.array(mfi,istate-AMREX_SPACEDIM);

	    amrex::ParallelFor(bx, [rho, Snew, dSsync, drhosync, ssync, iconserved ]
	    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	    {
	      dSsync(i,j,k,iconserved) = Snew(i,j,k) * drhosync(i,j,k) / rho(i,j,k);
	      ssync(i,j,k) -= dSsync(i,j,k);
	    });
	}
      }
    }

    if (do_mom_diff == 1)
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(rho_ctime, TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
	const Box& bx = mfi.tilebox();
	auto const& rho_c    = rho_ctime.array(mfi);
	auto const& vsync    = Vsync.array(mfi,Xvel);
	amrex::ParallelFor(bx, [rho_c, vsync]
	AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
	  for (int n = 0; n < AMREX_SPACEDIM; n++) {
	    vsync(i,j,k,n) /= rho_c(i,j,k);
	  }
	});
      }
    }
    //
    // Compute viscous sync.
    //
    if (is_diffusive[Xvel])
    {
      int rho_flag = (do_mom_diff == 0) ? 1 : 3;

      MultiFab** loc_viscn = 0;
      FluxBoxes fb_viscn;

      Real viscTime = state[State_Type].prevTime();
      loc_viscn = fb_viscn.define(this);
      getViscosity(loc_viscn, viscTime);

      diffusion->diffuse_Vsync(Vsync,dt,be_cn_theta,Rh,rho_flag,loc_viscn,0);
    }

    FluxBoxes fb_SC;
    MultiFab** fluxSC        = 0;
    bool       any_diffusive = false;
    for (int sigma  = 0; sigma < numscal; sigma++)
      if (is_diffusive[BL_SPACEDIM+sigma])
	any_diffusive = true;

    if (any_diffusive) {
      fluxSC = fb_SC.define(this);
    }

    Vector<int> diffuse_comp(1);
    int ng=1;
    const Real curr_time = state[State_Type].curTime();

    // Diffusion solver switches
    // together implies that Diff solve does NOT need Sold
    const Diffusion::SolveMode& solve_mode = Diffusion::ONEPASS;
    const bool add_old_time_divFlux = false;

    const int nlev = 1;
    Vector<MultiFab*> Snp1(nlev,0);

    MultiFab dSsync(grids,dmap,NUM_STATE,1,MFInfo(),Factory());
    Snp1[0] = &dSsync;

    Vector<MultiFab*> Rhonp1(nlev,0);
    Rhonp1[0] = &(get_new_data(State_Type));
    int Rho_comp = Density;

    FluxBoxes fb_fluxn  (this);
    MultiFab** fluxn   = fb_fluxn.get();

    const Vector<BCRec>& theBCs = AmrLevel::desc_lst[State_Type].getBCs();

    for (int sigma = 0; sigma<numscal; sigma++)
    {
      const int state_ind = BL_SPACEDIM + sigma;
      const int rho_flag  = Diffusion::set_rho_flag(diffusionType[state_ind]);

      if (is_diffusive[state_ind])
      {
        Snp1[0]->setVal(0.,state_ind,1,ng);

        FluxBoxes fb_diffnp1;
        MultiFab** cmp_diffnp1=0, **cmp_diffn=0;

        Real diffTime = state[State_Type].curTime();
        cmp_diffnp1 = fb_diffnp1.define(this);
        getDiffusivity(cmp_diffnp1, diffTime, BL_SPACEDIM+sigma,0,1);

        int S_comp = state_ind;
 	const int num_comp = 1;
	const int fluxComp  = 0;
        MultiFab *delta_rhs = &Ssync;
        int rhsComp = sigma;
        MultiFab *alpha_in = 0;
        const int alphaComp = 0;
        int betaComp = 0;

        diffuse_comp[0] = is_diffusive[BL_SPACEDIM+sigma];

        diffusion->diffuse_scalar ({},{},Snp1,Rhonp1,
 	                           S_comp,num_comp,Rho_comp,
                                   prev_time,curr_time,be_cn_theta,
                                   Rh,rho_flag,
                                   fluxn,fluxSC,fluxComp,
                                   delta_rhs,rhsComp,
                                   alpha_in,alphaComp,
                                   cmp_diffn,cmp_diffnp1,betaComp,
                                   crse_ratio,theBCs[state_ind],geom,
                                   solve_mode,
                                   add_old_time_divFlux,diffuse_comp);

        if (alpha_in!=0) delete alpha_in;

        MultiFab::Copy(Ssync,*Snp1[0],state_ind,sigma,1,0);

        //
        // Increment the viscous flux registers
        //
        if (level > 0)
        {
          for (int d = 0; d < BL_SPACEDIM; d++)
          {
             getViscFluxReg().FineAdd(*fluxSC[d],d,0,state_ind,1,dt);
          }
        }
      }
      else // state component not diffusive
      {
      //
      // The following used to be done in mac_sync_compute.  Ssync is
      // the source for a rate of change to S over the time step, so
      // Ssync*dt is the source to the actual sync amount.
      //
        Ssync.mult(dt,sigma,1,Ssync.nGrow());
      }
    }

    //
    // For all conservative variables Q (other than density)
    // increment sync by (sync_for_rho)*q_presync.
    // (See Pember, et. al., LBNL-41339, Jan. 1989)
    //
    iconserved = -1;
    for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
    {
      if (istate != Density && advectionType[istate] == Conservative)
      {
	iconserved++;

	MultiFab::Add(Ssync,*DeltaSsync,iconserved,istate-AMREX_SPACEDIM,1,0);
      }
    }
    //
    // Add the sync correction to the state.
    //
    MultiFab::Add(S_new,Ssync,0,AMREX_SPACEDIM,numscal,0);
    //
    // Update rho_ctime after rho is updated with Ssync.
    //
    make_rho_curr_time();

    if (level > 0) incrRhoAvg(Ssync,Density-BL_SPACEDIM,1.0);
    //
    // Get boundary conditions.
    //
    const int N = grids.size();

    Vector<int*>         sync_bc(N);
    Vector< Vector<int> > sync_bc_array(N);

    for (int i = 0; i < N; i++)
    {
      sync_bc_array[i] = getBCArray(State_Type,i,Density,numscal);
      sync_bc[i]       = sync_bc_array[i].dataPtr();
    }
    //
    // Interpolate the sync correction to the finer levels,
    //  and update rho_ctime, rhoAvg at those levels.
    //
    IntVect    ratio = IntVect::TheUnitVector();
    const Real mult  = 1.0;
    for (int lev = level+1; lev <= parent->finestLevel(); lev++)
    {
      ratio                     *= parent->refRatio(lev-1);
      NavierStokes&     fine_lev = getLevel(lev);
      const BoxArray& fine_grids = fine_lev.boxArray();
      MultiFab sync_incr(fine_grids,fine_lev.DistributionMap(),numscal,0,MFInfo(),fine_lev.Factory());
      sync_incr.setVal(0.0);

      SyncInterp(Ssync,level,sync_incr,lev,ratio,0,0,
		 numscal,1,mult,sync_bc.dataPtr());

      MultiFab& Sf_new = fine_lev.get_new_data(State_Type);
      MultiFab::Add(Sf_new,sync_incr,0,Density,numscal,0);

      fine_lev.make_rho_curr_time();
      fine_lev.incrRhoAvg(sync_incr,Density-BL_SPACEDIM,1.0);
    }

    sync_cleanup(DeltaSsync);

    BL_PROFILE_REGION_STOP("R::NavierStokes::mac_sync()");
}

//
// The reflux function
//
void
NavierStokes::reflux ()
{
    if (level == parent->finestLevel())
        return;

    BL_PROFILE("NavierStokes::reflux()");

    BL_ASSERT(do_reflux);
    //
    // First do refluxing step.
    //
    FluxRegister& fr_adv  = getAdvFluxReg(level+1);
    FluxRegister& fr_visc = getViscFluxReg(level+1);
    const Real    dt_crse = parent->dtLevel(level);
    const Real    scale   = 1.0/dt_crse;
    //
    // It is important, for do_mom_diff == 0, to do the viscous
    //   refluxing first, since this will be divided by rho_half
    //   before the advective refluxing is added.  In the case of
    //   do_mom_diff == 1, both components of the refluxing will
    //   be divided by rho^(n+1) in level_sync.
    //

    fr_visc.Reflux(Vsync,volume,scale,0,0,BL_SPACEDIM,geom);
    fr_visc.Reflux(Ssync,volume,scale,BL_SPACEDIM,0,NUM_STATE-BL_SPACEDIM,geom);

    const MultiFab& Rh = get_rho_half_time();

    if (do_mom_diff == 0)
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(Vsync,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box&  bx      = mfi.tilebox();
         auto const& vsync   = Vsync.array(mfi);
         auto const& rhohalf = Rh.array(mfi);

         amrex::ParallelFor(bx, AMREX_SPACEDIM, [vsync, rhohalf]
         AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
            vsync(i,j,k,n) /= rhohalf(i,j,k);
         });
      }
    }

    for (int istate = AMREX_SPACEDIM; istate < NUM_STATE; istate++)
    {
      if (advectionType[istate] == NonConservative)
      {
	MultiFab::Divide(Ssync,Rh,0,istate-AMREX_SPACEDIM,1,0);
      }
    }

    fr_adv.Reflux(Vsync,volume,scale,0,0,BL_SPACEDIM,geom);
    fr_adv.Reflux(Ssync,volume,scale,BL_SPACEDIM,0,NUM_STATE-BL_SPACEDIM,geom);

    const BoxArray& fine_boxes = getLevel(level+1).boxArray();
    //
    // Zero out coarse grid cells which underlie fine grid cells.
    //
    BoxArray baf = fine_boxes;

    baf.coarsen(fine_ratio);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   {
      std::vector< std::pair<int,Box> > isects;
      for (MFIter mfi(Vsync,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.growntilebox();
         auto const& vsync   = Vsync.array(mfi);
         auto const& ssync   = Ssync.array(mfi);
         int nstate          = NUM_STATE;

         baf.intersections(bx,isects);

         for (int it = 0, N = isects.size(); it < N; it++) {
            amrex::ParallelFor(isects[it].second, [vsync, ssync, nstate]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               for (int n = 0; n < AMREX_SPACEDIM; n++) {
                  vsync(i,j,k,n) = 0.0;
               }
               for (int n = 0; n < nstate-AMREX_SPACEDIM; n++) {
                  ssync(i,j,k,n) = 0.0;
               }
            });
        }
      }
   }
}

//
// Average fine information from the complete set of state types to coarse.
//

void
NavierStokes::avgDown ()
{
    if (level == parent->finestLevel())
        return;

    NavierStokes&   fine_lev = getLevel(level+1);
    //
    // Average down the states at the new time.
    //
    MultiFab& S_crse = get_new_data(State_Type);
    MultiFab& S_fine = fine_lev.get_new_data(State_Type);

    average_down(S_fine, S_crse, 0, S_crse.nComp());

    //
    // Now average down pressure over time n-(n+1) interval.
    //
    MultiFab&       P_crse      = get_new_data(Press_Type);
    MultiFab&       P_fine_init = fine_lev.get_new_data(Press_Type);
    MultiFab&       P_fine_avg  = fine_lev.p_avg;
    MultiFab&       P_fine      = initial_step ? P_fine_init : P_fine_avg;

    // NOTE: this fills ghost cells, but amrex::average_down does not.
    amrex::average_down_nodal(P_fine,P_crse,fine_ratio);
    //
    // Next average down divu and dSdT at new time.
    //
    if (have_divu)
    {
        MultiFab& Divu_crse = get_new_data(Divu_Type);
        MultiFab& Divu_fine = fine_lev.get_new_data(Divu_Type);

	average_down(Divu_fine, Divu_crse, 0, 1);
    }
    if (have_dsdt)
    {
        MultiFab& Dsdt_crse = get_new_data(Dsdt_Type);
        MultiFab& Dsdt_fine = fine_lev.get_new_data(Dsdt_Type);

	average_down(Dsdt_fine, Dsdt_crse, 0, 1);
    }
    //
    // Fill rho_ctime at the current and finer levels with the correct data.
    //
    for (int lev = level; lev <= parent->finestLevel(); lev++)
    {
        getLevel(lev).make_rho_curr_time();
    }
}

//
// Default divU is set to zero.
//

void
NavierStokes::calc_divu (Real      time,
                         Real      /*dt*/,
                         MultiFab& divu)
{
    BL_PROFILE("NavierStokes::calc_divu()");

    if (have_divu)
    {
      // Don't think we need this here, but then ghost cells are uninitialized
      // divu.setVal(0);

        if (do_temp && visc_coef[Temp] > 0.0)
        {
            //
            // Compute Div(U) = Div(visc_cond_coef * Grad(T))/(c_p*rho*T)
            //
            getViscTerms(divu,Temp,1,time);

            const MultiFab&   rhotime = get_rho(time);

            FillPatchIterator temp_fpi(*this,divu,0,time,State_Type,Temp,1);
	    MultiFab& tmf = temp_fpi.get_mf();

	    Real THERMO_cp = 1004.6;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for ( MFIter rho_mfi(rhotime,TilingIfNotGPU()); rho_mfi.isValid(); ++rho_mfi)
            {
	        const Box&  bx  = rho_mfi.tilebox();
		auto const& div = divu.array(rho_mfi);
		auto const& rho = rhotime.array(rho_mfi);
		auto const& temp = tmf.array(rho_mfi);
#ifdef AMREX_USE_EB
		auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());
		auto const& flagfab = ebfactory.getMultiEBCellFlagFab()[rho_mfi];
		auto const& flag    = flagfab.const_array();

		if (flagfab.getType(bx) == FabType::covered)
		{
		  amrex::ParallelFor(bx, [div]
		  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
		  {
		    div( i, j, k ) = COVERED_VAL;
		  });
		}
		else if (flagfab.getType(bx) != FabType::regular)
		{
		  auto vfrac = ebfactory.getVolFrac().const_array(rho_mfi);

		  amrex::ParallelFor(bx, [div, rho, temp, vfrac, THERMO_cp]
		  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
		  {
		    if ( vfrac(i,j,k) > 0.0 )
		    {
		      div(i,j,k) /= ( rho(i,j,k)*temp(i,j,k)*THERMO_cp );
		    }
		    else
		    {
		      div(i,j,k) = COVERED_VAL;
		    }
		  });
		}
		else
#endif
		{
		  amrex::ParallelFor(bx, [div, rho, temp, THERMO_cp]
		  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
		  {
		    div(i,j,k) /= ( rho(i,j,k)*temp(i,j,k)*THERMO_cp );
		  });
		}
	    }
        }
	else
	{
	  divu.setVal(0);
	}
    }
}

void
NavierStokes::getViscTerms (MultiFab& visc_terms,
                            int       src_comp,
                            int       ncomp,
                            Real      time)
{
    BL_PROFILE("NavierStokes::getViscTerms()");
    //
    // The logic below for selecting between scalar or tensor solves does
    // not allow for calling NavierStokes::getViscTerms with src_comp=Yvel
    // or Zvel
    //
#ifdef AMREX_DEBUG
    if (src_comp<BL_SPACEDIM && (src_comp!=Xvel || ncomp<BL_SPACEDIM))
    {
      Print() << "src_comp=" << src_comp << "   ncomp=" << ncomp << '\n';
      Error("must call NavierStokes::getViscTerms with all three velocity components");
    }
#endif
    //
    // Initialize all viscous terms to zero
    //
    const int nGrow = visc_terms.nGrow();

    bool diffusive = false;
    //
    // Get Velocity Viscous Terms
    //
    if (src_comp == Xvel && !is_diffusive[Xvel])
    {
	visc_terms.setVal(0.0,0,ncomp,nGrow);
    }
    else if (src_comp == Xvel && is_diffusive[Xvel])
    {
	diffusive = true;

	FluxBoxes fb;
        MultiFab** viscosity = 0;

        viscosity = fb.define(this);
        getViscosity(viscosity, time);

	auto whichTime = which_time(State_Type,time);
	BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

	auto viscosityCC = (whichTime == AmrOldTime ? viscn_cc : viscnp1_cc);

        diffusion->getTensorViscTerms(visc_terms,time,viscosity,viscosityCC,0);

        //
        // Add Div(u) term if desired, if this is velocity, and if Div(u)
        // is nonzero.  If const-visc, term is mu.Div(u)/3, else
        // it's -Div(mu.Div(u).I)*2/3
        //
        if (have_divu && S_in_vel_diffusion)
        {
            MultiFab divmusi(grids,dmap,BL_SPACEDIM,1,MFInfo(),Factory());

            diffusion->compute_divmusi(time,viscosity,divmusi);
            divmusi.mult((-2./3.),0,BL_SPACEDIM,0);

            visc_terms.plus(divmusi,Xvel,BL_SPACEDIM,0);
        }
    }
    //
    // Get Scalar Diffusive Terms
    //
    const int first_scal = (src_comp==Xvel) ? BL_SPACEDIM : src_comp;
    const int num_scal   = (src_comp==Xvel) ? ncomp-BL_SPACEDIM : ncomp;

    if (num_scal > 0)
    {
        for (int icomp = first_scal; icomp < first_scal+num_scal; icomp++)
        {
            if (is_diffusive[icomp])
            {
		diffusive = true;

                int rho_flag = Diffusion::set_rho_flag(diffusionType[icomp]);

		FluxBoxes fb;
                MultiFab** cmp_diffn = 0;

                cmp_diffn = fb.define(this);
                getDiffusivity(cmp_diffn, time, icomp, 0, 1);

                diffusion->getViscTerms(visc_terms,src_comp,icomp,
                                        time,rho_flag,cmp_diffn,0);
            }
	    else {
	        visc_terms.setVal(0.0,icomp-src_comp,1,nGrow);
	    }

        }
    }
    //
    // Ensure consistent grow cells
    //
    if (diffusive && nGrow > 0)
    {
	    visc_terms.FillBoundary(0, ncomp, geom.periodicity());
	    Extrapolater::FirstOrderExtrap(visc_terms, geom, 0, ncomp);
    }
}

//
// Functions calcViscosity/Diffusivity and getViscosity/Diffusivity are
// for calculating variable viscosity and diffusivity. Here we default to
// constant visc/diff and set the variable viscosity and diffusivity arrays
// to the values in visc_coef and diff_coef.
// For variable viscosity/diffusivity, (per MSD) calcViscosity/Diffusivity
// should compute the transport coefficients at cell centers (or cell centroids
// for EB) and getViscosity/Diffusivity should interpolate those to faces (or
// face-centroids for EB).
//




void
NavierStokes::calcViscosity (const Real time,
                             const Real /*dt*/,
                             const int  /*iteration*/,
                             const int  /*ncycle*/)
{
    double muG = 0.0;
    double muL = 0.0;
    double muL2 = 0.0;
    int do_variable_visc = 0;
    int do_three_fluid = 0;
    ParmParse pp("ns");
    pp.query("muG", muG);
    pp.query("muL", muL);
    pp.query("muL2", muL2);
    pp.query("do_variable_visc", do_variable_visc);
    if (is_diffusive[Xvel])
    {
        if (visc_coef[Xvel] >= 0.0)
        {
            auto whichTime = which_time(State_Type,time);
            BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

            MultiFab& visc = (whichTime == AmrOldTime ? *viscn_cc : *viscnp1_cc);
               MultiFab& S = (whichTime == AmrOldTime) ? get_old_data(State_Type) : get_new_data(State_Type);

            FillPatchIterator fpi(*this,S,1,time,State_Type,0,S.nComp());
            MultiFab& S_cc = fpi.get_mf();
            MultiFab viscL(visc, amrex::make_alias, 0, visc.nComp());
            MultiFab trac_fab;
            trac_fab.define(grids, dmap, 1, S.nGrow());
            MultiFab::Copy(trac_fab, S, Tracer, 0, 1, S.nGrow());
            trac_fab.FillBoundary(geom.periodicity());
            trac_fab.mult((muG-muL), 0, trac_fab.nComp(), trac_fab.nGrow());
            viscL.setVal(muL);

            if(do_variable_visc)
            {
                // visc.plus(viscL,0, 1, visc.nGrow());
                // visc.plus(trac_fab,0, 1, visc.nGrow());
                #ifdef _OPENMP
                #pragma omp parallel if (Gpu::notInLaunchRegion())
                #endif
                   for (MFIter mfi(S_cc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
                   {
                      const Box& gbx = mfi.growntilebox();
                      auto lo = lbound(gbx);
                      auto hi = ubound(gbx);
                      auto const& rho    = S_cc.array(mfi,Density);
                      auto const& Trac       = S_cc.array(mfi,Tracer);
                      auto const& Trac2       = S_cc.array(mfi,Tracer2);
                      // Array4<Real>  Trac3;
                      // if(do_trac3)
                      // {
                          auto const& Trac3       = S_cc.array(mfi,Tracer3);
                      // }
                      auto const& mu      = visc.array(mfi);
                      for(int i = lo.x; i <= hi.x; i++)
                      {
                          for(int j = lo.y; j <= hi.y; j++)
                          {
                              for(int k = lo.z; k <= hi.z; k++)
                              {
                                  if(do_trac3)
                                  {
                                      mu(i,j,k) = muG*Trac(i,j,k,0) + muL*Trac2(i,j,k,0)
                                        + muL2*Trac3(i,j,k,0);
                                        // if(mu(i,j,k) < fmin(muL, fmin(muL2, muG))){ mu(i,j,k) = fmin(muL, fmin(muL2, muG)); }
                                        // else if(mu(i,j,k) > fmax(muL, fmax(muL2, muG))) { mu(i,j,k) = fmax(muL, fmax(muL2, muG)); }

                                  }
                                  else
                                  {
                                      mu(i,j,k) = muL + (muG-muL)*Trac(i,j,k,0);
                                      if(mu(i,j,k) < fmin(muL, muG)){ mu(i,j,k) = fmin(muL, muG); }
                                      else if(mu(i,j,k) > fmax(muL, muG)) { mu(i,j,k) = fmax(muL, muG); }

                                  }
                              }

                          }
                      }
                    }
              // MultiFab::LinComb(*visc, 1.0, viscL, 0, (muG-muL), face_trac_fabs[dir], 0, 0, visc[dir]->nComp(), visc[dir]->nGrow());

            }
            else
            {
                visc.setVal(visc_coef[Xvel], 0, visc.nComp(), visc.nGrow());

            }

        }
        else
	{
            Abort("NavierStokes::calcViscosity() : must have velocity visc_coef >= 0.0");
        }
    }
}


void
NavierStokes::calcDiffusivity (const Real time)
{
    //
    // NOTE:  In the diffusivity
    //        arrays, there is an offset since no diffusivity array
    //        is kept for the velocities or the density.  So, the scalar
    //        component Density+1 in the state corresponds to component
    //        0 in the arrays diffn and diffnp1.
    //
    int src_comp = Density+1;
    int ncomp = NUM_STATE - BL_SPACEDIM -1;

    const TimeLevel whichTime = which_time(State_Type,time);
    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    MultiFab* diff = (whichTime == AmrOldTime ? diffn_cc : diffnp1_cc);
    for (int comp=src_comp; comp<src_comp+ncomp; comp++)
    {
        int diff_comp = comp - Density - 1;

        if (is_diffusive[comp])
        {
            if (visc_coef[comp] >= 0.0)
            {
	      diff->setVal(visc_coef[comp], diff_comp, 1, diff->nGrow());
            }
            else
            {
                Abort("NavierStokes::calcDiffusivity() : must have scalar diff_coefs >= 0.0");
            }
        }
    }
}




void
NavierStokes::getViscosity (MultiFab* viscosity[BL_SPACEDIM],
                            const Real time)
{
    ParmParse pp("ns");
    int do_variable_visc = 0;
    pp.query("do_variable_visc", do_variable_visc);

    // //
    // // Select time level to work with (N or N+1)
    // //
    if(do_variable_visc)
    {

    const TimeLevel whichTime = which_time(State_Type,time);
    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    MultiFab *visc = (whichTime == AmrOldTime ? viscn_cc : viscnp1_cc);
   const Box& domain = geom.Domain();

   #ifdef AMREX_USE_EB
      // EB : use EB CCentroid -> FCentroid
      auto math_bc = fetchBCArray(State_Type,0,1);
      EB_interp_CellCentroid_to_FaceCentroid(*visc, D_DECL(*viscosity[0],*viscosity[1],*viscosity[2]),
                                             0, 0, 1, geom, math_bc);
      EB_set_covered_faces({D_DECL(viscosity[0],viscosity[1],viscosity[2])},0.0);
   #else
    #ifdef _OPENMP
    #pragma omp parallel if (Gpu::notInLaunchRegion())
    #endif
       for (MFIter mfi(*visc,TilingIfNotGPU()); mfi.isValid();++mfi)
       {
          for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
          {
             const Box ebx = mfi.nodaltilebox(dir);
             const Box& edomain = amrex::surroundingNodes(domain,dir);
             const auto& visc_c  = (*visc)[mfi];
             auto& visc_ed = (*viscosity[dir])[mfi];

             center_to_edge_plain(visc_c, visc_ed, ebx,0,0, viscosity[dir]->nComp());
         }
        }

    #endif
    }
    else
    {
        for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
          viscosity[dir]->setVal(visc_coef[Xvel], 0, viscosity[dir]->nComp(), viscosity[dir]->nGrow());
        }
    }
    // For non-const viscosity, uncomment above and add interp from
    // cell-center/centroid to faces.
    // But here we simply do constant viscosity.



    if (do_LES)
    {
      FluxBoxes mu_LES(this,1,0);
      MultiFab** mu_LES_mf = mu_LES.get();
      for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
	mu_LES_mf[dir]->setVal(0., 0, mu_LES_mf[dir]->nComp(), mu_LES_mf[dir]->nGrow());
      }

      NavierStokesBase::calc_mut_LES(mu_LES_mf,time);

      for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
	MultiFab::Add(*viscosity[dir], *mu_LES_mf[dir], 0, 0, 1, 0);
      }
    }
}

void
NavierStokes::getDiffusivity (MultiFab* diffusivity[BL_SPACEDIM],
                              const Real time,
                              const int state_comp,
                              const int dst_comp,
                              const int ncomp)
{
    BL_ASSERT(state_comp > Density);
    // //
    // // Pick correct component in the diffn/diffnp1 array
    // //
    // int diff_comp = state_comp - Density - 1;
    // //
    // // Select time level to work with (N or N+1)
    // //
    // const TimeLevel whichTime = which_time(State_Type,time);
    // BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    // MultiFab *diff = (whichTime == AmrOldTime ? diffn_cc : diffnp1_cc);

    // For non-const diffusivity, uncomment above and add interp from
    // cell-center/centroid to faces.
    // But here we simply do constant diffusivity.

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      diffusivity[dir]->setVal(visc_coef[state_comp], dst_comp, ncomp, diffusivity[dir]->nGrow());
    }
}


void
NavierStokes::center_to_edge_plain (const FArrayBox& ccfab,
                                    FArrayBox&       ecfab,
				    const Box&       bx,
                                    int              sComp,
                                    int              dComp,
                                    int              nComp)
{
    //
    // This routine fills an edge-centered FAB from a cell-centered FAB.
    // It assumes that the data in all cells of the cell-centered FAB is
    // valid and totally ignores any concept of boundary conditions.
    // It is assummed that the cell-centered FAB fully contains the
    // edge-centered FAB.  If anything special needs to be done at boundaries,
    // a varient of this routine needs to be written.  See
    // HeatTransfer::center_to_edge_fancy().
    //
    const Box&      ccbox = ccfab.box();
    const IndexType ixt   = ecfab.box().ixType();
    //
    // Get direction for interpolation to edges
    //
    int dir = -1;
    for (int d = 0; d < BL_SPACEDIM; d++)
        if (ixt.test(d))
            dir = d;
    //
    // Miscellanious checks
    //
    BL_ASSERT(!(ixt.cellCentered()) && !(ixt.nodeCentered()));
    BL_ASSERT(grow(ccbox,-BASISV(dir)).contains(enclosedCells(bx)));
    BL_ASSERT(sComp+nComp <= ccfab.nComp() && dComp+nComp <= ecfab.nComp());

    //
    // Shift cell-centered data to edges
    //
    const int isharm = def_harm_avg_cen2edge;

    cen2edg(bx.loVect(), bx.hiVect(),
	    ARLIM(ccfab.loVect()), ARLIM(ccfab.hiVect()),
	    ccfab.dataPtr(sComp),
	    ARLIM(ecfab.loVect()), ARLIM(ecfab.hiVect()),
	    ecfab.dataPtr(dComp),
	    &nComp, &dir, &isharm);
}
