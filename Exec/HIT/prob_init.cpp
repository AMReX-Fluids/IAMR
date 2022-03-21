#include <NavierStokes.H>
#include <AMReX_ParmParse.H>
#include <iamr_constants.H>

#ifdef AMREX_USE_TURBULENT_FORCING
#include <TurbulentForcing_def.H>
#endif

using namespace amrex;

int NavierStokes::probtype = -1;


//
// Initialize state and pressure with problem-specific data
//
void NavierStokes::prob_initData ()
{
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
    auto const& probhi = geom.ProbHiArray();

#ifdef AMREX_USE_TURBULENT_FORCING
    //
    // Initialize data structures used for homogenous isentropic forced turbulence.
    //
    TurbulentForcing::init_turbulent_forcing(problo,probhi);
#endif

    //
    // Create struct to hold initial conditions parameters
    //
    InitialConditions IC;

    //
    // Read problem parameters from inputs file
    //
    ParmParse pp("prob");

    pp.query("probtype",probtype);

    if ( probtype == 100 )
    {
	//
	// Random combination of cosine waves to be used with forced turbulence,
	// where ICs are less important as the forcing takes over with time.
	//

	pp.query("turb_scale",IC.turb_scale);
	pp.query("density_ic",IC.density);

#ifdef _OPENMP
#pragma omp parallel  if (Gpu::notInLaunchRegion())
#endif
	for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
	    const Box& vbx = mfi.tilebox();

	    init_forced(vbx, /*P_new.array(mfi),*/ S_new.array(mfi, Xvel),
			S_new.array(mfi, Density), nscal,
			domain, dx, problo, probhi, IC);
	}
    }
    else if ( probtype == 101 )
    {
	//
	// Initialize homogenous isentropic turbulence from file.
	// Can be run without forcing for a decaying HIT problem.

	// Some additional parameters for creating IC's
	std::string iname;
	bool binfmt = false;
	amrex::Real urms0 = 1.0;
	amrex::Real uin_norm = 1.0;

	pp.query("inres",IC.inres);
	if (IC.inres == 0){
	    amrex::Abort("for HIT, inres cannot be 0 !");
	}
	pp.query("iname",iname);
	pp.query("binfmt",binfmt);
	pp.query("urms0",urms0);


	//
	// Output info about initial turbulence
	//
	amrex::Real lambda0 = 0.5;
	amrex::Real tau  = lambda0 / urms0;
	// Output IC
	std::ofstream ofs("ic.txt", std::ofstream::out);
	amrex::Print(ofs)
	    << "lambda0, urms0, tau "
	    << std::endl;
	amrex::Print(ofs).SetPrecision(17)
	    << lambda0 << "," << urms0 << "," << tau << std::endl;
	ofs.close();

	//
	// Load velocity fields from file.
	//
	// Assumes data set ordered in Fortran format.
	// We will interpolate this data onto our domain box.
	// Assumes the input data is a periodic cube. If the input cube is smaller
	// than our domain size, the cube will be repeated throughout the
	// domain (hence the mod operations in the interpolation).
	const size_t nx = IC.inres;
	const size_t ny = IC.inres;
	const size_t nz = IC.inres;
	amrex::Vector<amrex::Real> data(
	    nx * ny * nz * 6); /* this needs to be double */
	if (binfmt) {
	    read_binary(iname, nx, ny, nz, 6, data);
	} else {
	    read_csv(iname, nx, ny, nz, data);
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
	    uinput[i] = data[3 + i * 6] * urms0 / uin_norm;
	    vinput[i] = data[4 + i * 6] * urms0 / uin_norm;
	    winput[i] = data[5 + i * 6] * urms0 / uin_norm;
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

	// fixme? this is not the most efficient way to get the data on the gpu,
	// since AsyncArray also makes it's own copy of the host data.
	amrex::Gpu::AsyncArray<amrex::Real> xarray_aa(xarray.data(),xarray.size());
	amrex::Gpu::AsyncArray<amrex::Real> xdiff_aa(xdiff.data(),xdiff.size());
	amrex::Gpu::AsyncArray<amrex::Real> uinput_aa(uinput.data(),uinput.size());
	amrex::Gpu::AsyncArray<amrex::Real> vinput_aa(vinput.data(),vinput.size());
	amrex::Gpu::AsyncArray<amrex::Real> winput_aa(winput.data(),winput.size());

	IC.d_xarray = xarray_aa.data();
	IC.d_xdiff  = xdiff_aa.data();
	IC.d_uinput = uinput_aa.data();
	IC.d_vinput = vinput_aa.data();
	IC.d_winput = winput_aa.data();



#ifdef _OPENMP
#pragma omp parallel  if (Gpu::notInLaunchRegion())
#endif
	for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
	    const Box& vbx = mfi.tilebox();

	    init_HIT(vbx, /*P_new.array(mfi),*/ S_new.array(mfi, Xvel),
		     S_new.array(mfi, Density), nscal,
		     /*domain,*/ dx, problo, /*probhi,*/ IC);
	}
    }
    else
    {
	amrex::Abort("NavierStokes::prob_init: unknown probtype");
    }
}

void NavierStokes::init_forced (Box const& vbx,
				/* Array4<Real> const& press, */
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

  amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    Real x = problo[0] + (i - domlo.x + 0.5)*dx[0];
    Real y = problo[1] + (j - domlo.y + 0.5)*dx[1];
#if (AMREX_SPACEDIM == 3)
    Real z = problo[2] + (k - domlo.z + 0.5)*dx[2];
#else
    constexpr Real z = 0.0;
#endif

    const Real Lx    = (probhi[0] - problo[0]);
    const Real Ly    = (probhi[1] - problo[1]);
#if (AMREX_SPACEDIM == 3)
    const Real Lz    = (probhi[2] - problo[1]);
#else
    const Real Lz    = 1.0;
#endif
    //
    // Fill Velocity
    //
    AMREX_D_TERM(vel(i,j,k,0) =  IC.turb_scale * std::cos(TwoPi*y/Ly) * std::cos(TwoPi*z/Lz);,
		 vel(i,j,k,1) =  IC.turb_scale * std::cos(TwoPi*x/Lx) * std::cos(TwoPi*z/Lz);,
		 vel(i,j,k,2) =  IC.turb_scale * std::cos(TwoPi*x/Lx) * std::cos(TwoPi*y/Ly););

    //
    // Scalars, ordered as Density, Tracer(s)
    //
    scal(i,j,k,0) = IC.density;

    // Tracers
    for ( int nt=1; nt<nscal; nt++)
    {
      scal(i,j,k,nt) = 1.0;
    }
  });
}

void NavierStokes::init_HIT (Box const& vbx,
			     /* Array4<Real> const& press, */
			     Array4<Real> const& vel,
			     Array4<Real> const& scal,
			     const int nscal,
			     /*Box const& domain,*/
			     GpuArray<Real, AMREX_SPACEDIM> const& dx,
			     GpuArray<Real, AMREX_SPACEDIM> const& problo,
			     /* GpuArray<Real, AMREX_SPACEDIM> const& probhi,*/
			     InitialConditions IC)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(parent->maxLevel()==0, "Decaying turbulence is single level only. Set amr.max_level = 0");
#if (AMREX_SPACEDIM != 3)
    amrex::Abort("NavierStokes::init_HIT: This is a 3D problem, please recompile with DIM=3 in makefile");
#endif

    amrex::ParallelFor(vbx, [vel, scal, nscal, problo, dx, IC]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
	amrex::Real x[3] = {
	    problo[0] + static_cast<amrex::Real>(i + 0.5) * dx[0],
	    problo[1] + static_cast<amrex::Real>(j + 0.5) * dx[1],
	    problo[2] + static_cast<amrex::Real>(k + 0.5) * dx[2]};

	// Fill in the velocities and energy.
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
