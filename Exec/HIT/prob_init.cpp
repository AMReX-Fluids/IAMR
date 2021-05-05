#include <NavierStokes.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

int NavierStokes::probtype = -1;

// For now, define pi here, but maybe later make iamr_constants.H
namespace {
  constexpr Real Pi    = 3.141592653589793238462643383279502884197;
  constexpr Real TwoPi = 2.0 * 3.141592653589793238462643383279502884197;
}

//
// Initialize state and pressure with problem-specific data
//
void NavierStokes::prob_initData ()
{

    //
    // Create struct to hold initial conditions parameters
    //
    InitialConditions IC;

    //
    // Read problem parameters from inputs file
    //
    ParmParse pp("prob");

    // for HIT
    pp.query("inres",IC.inres);
    if (IC.inres == 0){
      amrex::Abort("for HIT, inres cannot be 0 !");
    }
    pp.query("iname",IC.iname);
    pp.query("binfmt",IC.binfmt);
    pp.query("urms0",IC.urms0);

    //
    // Fill state and, optionally, pressure
    //
    MultiFab& P_new = get_new_data(Press_Type);
    MultiFab& S_new = get_new_data(State_Type);
    const int nscal = NUM_STATE-Density;

    S_new.setVal(0.0);
    P_new.setVal(0.0);

    // Integer indices of the lower left and upper right corners of the
    // valid region of the entire domain.
    Box const&  domain = geom.Domain();
    auto const&     dx = geom.CellSizeArray();
    // Physical coordinates of the lower left corner of the domain
    auto const& problo = geom.ProbLoArray();
    // Physical coordinates of the upper right corner of the domain
    auto const& probhi = geom.ProbHiArray();

    amrex::Real lambda0 = 0.5;
    amrex::Real tau  = lambda0 / IC.urms0;
    // Output IC
    std::ofstream ofs("ic.txt", std::ofstream::out);
    amrex::Print(ofs)
      << "lambda0, urms0, tau "
      << std::endl;
    amrex::Print(ofs).SetPrecision(17)
      << lambda0 << "," << IC.urms0 << "," << tau << std::endl;
    ofs.close();


    // Load velocity fields from file. Assume data set ordered in Fortran
    // format and reshape the data accordingly. One thing to keep in mind
    // is that this contains the entire input data. We will interpolate
    // this data later to just match our box. Another assumption is that
    // the input data is a periodic cube. If the input cube is smaller
    // than our domain size, the cube will be repeated throughout the
    // domain (hence the mod operations in the interpolation).
    const size_t nx = IC.inres;
    const size_t ny = IC.inres;
    const size_t nz = IC.inres;
    amrex::Vector<amrex::Real> data(
      nx * ny * nz * 6); /* this needs to be double */
    if (IC.binfmt) {
      read_binary(IC.iname, nx, ny, nz, 6, data);
    } else {
      read_csv(IC.iname, nx, ny, nz, data);
    }

    // Extract position and velocities
    amrex::Vector<amrex::Real> xinput;
    amrex::Vector<amrex::Real> uinput;
    amrex::Vector<amrex::Real> vinput;
    amrex::Vector<amrex::Real> winput;
    amrex::Vector<amrex::Real> xdiff;
    amrex::Vector<amrex::Real> xarray;

    xinput.resize(nx * ny * nz);
    uinput.resize(nx * ny * nz);
    vinput.resize(nx * ny * nz);
    winput.resize(nx * ny * nz);

    for (long i = 0; i < xinput.size(); i++) {
      xinput[i] = data[0 + i * 6];
      uinput[i] = data[3 + i * 6] * IC.urms0 / IC.uin_norm;
      vinput[i] = data[4 + i * 6] * IC.urms0 / IC.uin_norm;
      winput[i] = data[5 + i * 6] * IC.urms0 / IC.uin_norm;
    }

    // Get the xarray table and the differences.
    xarray.resize(nx);
    for (long i = 0; i < xarray.size(); i++) {
      xarray[i] = xinput[i];
    }
    xdiff.resize(nx);
    std::adjacent_difference(
      xarray.begin(),
      xarray.end(),
      xdiff.begin());
    xdiff[0] = xdiff[1];

    // Make sure the search array is increasing
    if (not std::is_sorted(
          xarray.begin(),
          xarray.end())) {
      amrex::Abort("Error: non ascending x-coordinate array.");
    }

    // Dimensions of the input box.
    IC.Linput = xarray[nx - 1] + 0.5 * xdiff[nx - 1];

    IC.d_xarray = xarray.data();
    IC.d_xdiff  = xdiff.data();
    IC.d_uinput = uinput.data(); 
    IC.d_vinput = vinput.data();
    IC.d_winput = winput.data();


#ifdef _OPENMP
#pragma omp parallel  if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.tilebox();

          init_HIT(vbx, P_new.array(mfi), S_new.array(mfi, Xvel),
                           S_new.array(mfi, Density), nscal,
                           domain, dx, problo, probhi, IC);
    }
}


void NavierStokes::init_HIT (Box const& vbx,
                               Array4<Real> const& press,
                               Array4<Real> const& vel,
                               Array4<Real> const& scal,
                               const int nscal,
                               Box const& domain,
                               GpuArray<Real, AMREX_SPACEDIM> const& dx,
                               GpuArray<Real, AMREX_SPACEDIM> const& problo,
                               GpuArray<Real, AMREX_SPACEDIM> const& probhi,
                               InitialConditions IC)
{
  const auto domlo = amrex::lbound(domain);

#if (AMREX_SPACEDIM != 3)
    amrex::Abort("NavierStokes::init_HIT: This is a 3D problem, please recompile with DIM=3 in makefile");
#endif

  constexpr Real eps_input=0.05, rho_input=0.15;

  amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    amrex::Real x[3] = {
      problo[0] + static_cast<amrex::Real>(i + 0.5) * dx[0],
      problo[1] + static_cast<amrex::Real>(j + 0.5) * dx[1],
      problo[2] + static_cast<amrex::Real>(k + 0.5) * dx[2]};

    // Fill in the velocities and energy.
    amrex::Real u[3] = {0.0};
    amrex::Real uinterp[3] = {0.0};

    // Interpolation factors
    amrex::Real mod[3] = {0.0};
    int idx[3] = {0};
    int idxp1[3] = {0};
    amrex::Real slp[3] = {0.0};
    for (int cnt = 0; cnt < 3; cnt++) {
      mod[cnt] = std::fmod(x[cnt], IC.Linput);
      locate(IC.d_xarray, IC.inres, mod[cnt], idx[cnt]);
      idxp1[cnt] = (idx[cnt] + 1) % IC.inres;
      slp[cnt] =
        (mod[cnt] - IC.d_xarray[idx[cnt]]) / IC.d_xdiff[idx[cnt]];
    }

    const amrex::Real f0 = (1 - slp[0]) * (1 - slp[1]) * (1 - slp[2]);
    const amrex::Real f1 = slp[0] * (1 - slp[1]) * (1 - slp[2]);
    const amrex::Real f2 = (1 - slp[0]) * slp[1] * (1 - slp[2]);
    const amrex::Real f3 = (1 - slp[0]) * (1 - slp[1]) * slp[2];
    const amrex::Real f4 = slp[0] * (1 - slp[1]) * slp[2];
    const amrex::Real f5 = (1 - slp[0]) * slp[1] * slp[2];
    const amrex::Real f6 = slp[0] * slp[1] * (1 - slp[2]);
    const amrex::Real f7 = slp[0] * slp[1] * slp[2];

    uinterp[0] =
      IC.d_uinput
          [idx[0] + IC.inres * (idx[1] + IC.inres * idx[2])] *
        f0 +
      IC.d_uinput
          [idxp1[0] + IC.inres * (idx[1] + IC.inres * idx[2])] *
        f1 +
      IC.d_uinput
          [idx[0] + IC.inres * (idxp1[1] + IC.inres * idx[2])] *
        f2 +
      IC.d_uinput
          [idx[0] + IC.inres * (idx[1] + IC.inres * idxp1[2])] *
        f3 +
      IC.d_uinput
          [idxp1[0] + IC.inres * (idx[1] + IC.inres * idxp1[2])] *
        f4 +
      IC.d_uinput
          [idx[0] + IC.inres * (idxp1[1] + IC.inres * idxp1[2])] *
        f5 +
      IC.d_uinput
          [idxp1[0] + IC.inres * (idxp1[1] + IC.inres * idx[2])] *
        f6 +
      IC.d_uinput
          [idxp1[0] + IC.inres * (idxp1[1] + IC.inres * idxp1[2])] *
        f7;

    uinterp[1] =
      IC.d_vinput
          [idx[0] + IC.inres * (idx[1] + IC.inres * idx[2])] *
        f0 +
      IC.d_vinput
          [idxp1[0] + IC.inres * (idx[1] + IC.inres * idx[2])] *
        f1 +
      IC.d_vinput
          [idx[0] + IC.inres * (idxp1[1] + IC.inres * idx[2])] *
        f2 +
      IC.d_vinput
          [idx[0] + IC.inres * (idx[1] + IC.inres * idxp1[2])] *
        f3 +
      IC.d_vinput
          [idxp1[0] + IC.inres * (idx[1] + IC.inres * idxp1[2])] *
        f4 +
      IC.d_vinput
          [idx[0] + IC.inres * (idxp1[1] + IC.inres * idxp1[2])] *
        f5 +
      IC.d_vinput
          [idxp1[0] + IC.inres * (idxp1[1] + IC.inres * idx[2])] *
        f6 +
      IC.d_vinput
          [idxp1[0] + IC.inres * (idxp1[1] + IC.inres * idxp1[2])] *
        f7;

    uinterp[2] =
      IC.d_winput
          [idx[0] + IC.inres * (idx[1] + IC.inres * idx[2])] *
        f0 +
      IC.d_winput
          [idxp1[0] + IC.inres * (idx[1] + IC.inres * idx[2])] *
        f1 +
      IC.d_winput
          [idx[0] + IC.inres * (idxp1[1] + IC.inres * idx[2])] *
        f2 +
      IC.d_winput
          [idx[0] + IC.inres * (idx[1] + IC.inres * idxp1[2])] *
        f3 +
      IC.d_winput
          [idxp1[0] + IC.inres * (idx[1] + IC.inres * idxp1[2])] *
        f4 +
      IC.d_winput
          [idx[0] + IC.inres * (idxp1[1] + IC.inres * idxp1[2])] *
        f5 +
      IC.d_winput
          [idxp1[0] + IC.inres * (idxp1[1] + IC.inres * idx[2])] *
        f6 +
      IC.d_winput
          [idxp1[0] + IC.inres * (idxp1[1] + IC.inres * idxp1[2])] *
        f7;

    //
    // Fill Velocity
    //
    vel(i,j,k,0) = uinterp[0];
    vel(i,j,k,1) = uinterp[1];
    vel(i,j,k,2) = uinterp[2];

    //
    // Scalars, ordered as Density, Tracer(s)
    //
    scal(i,j,k,0) = 1.0;
    scal(i,j,k,1) = 0.0;

    // Additional Tracer, if using
    for ( int nt=2; nt<nscal; nt++)
    {
      scal(i,j,k,nt) = 1.0;
    }
  });
}

