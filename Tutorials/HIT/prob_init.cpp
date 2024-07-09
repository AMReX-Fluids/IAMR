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

// FIXME - should remove ifdef and use runtime parameter instead...
#ifdef AMREX_USE_TURBULENT_FORCING
    //
    // Initialize data structures used for homogeneous isentropic forced turbulence.
    // Only need to do it once.
    if (level == 0)
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
