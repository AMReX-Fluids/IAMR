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

    pp.query("probtype",probtype);
    pp.query("density_ic",IC.density);
    pp.query("direction",IC.direction);
    pp.query("interface_width",IC.interface_width);

    Vector<Real> velocity(AMREX_SPACEDIM, 0.);
    pp.queryarr("velocity_ic",velocity,0,AMREX_SPACEDIM);
    AMREX_D_TERM(IC.v_x = velocity[0];,
		 IC.v_y = velocity[1];,
		 IC.v_z = velocity[2];);

    pp.query("blob_radius",IC.blob_radius);
    Vector<Real> blob_center(AMREX_SPACEDIM, 0.);
    pp.queryarr("blob_center",blob_center,0,AMREX_SPACEDIM);
    AMREX_D_TERM(IC.blob_x = blob_center[0];,
		 IC.blob_y = blob_center[1];,
		 IC.blob_z = blob_center[2];);

    // For Rayleigh-Taylor problem
    pp.query("rho_1",IC.rho_1);
    pp.query("rho_2",IC.rho_2);
    pp.query("tra_1",IC.tra_1);
    pp.query("tra_2",IC.tra_2);
    pp.query("perturbation_amplitude",IC.pertamp);

    // for Taylor-Green
    pp.query("velocity_factor",IC.v_x);

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

#ifdef _OPENMP
#pragma omp parallel  if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.tilebox();

	if ( 1 == probtype )
        {
	  //
	  // Start from rest, constant density of 1
	  // Introduced for LidDrivenCavity problem
	  //
	  S_new[mfi].setVal<RunOn::Gpu>(1.0,Density);
	}
        else if ( 2 == probtype || 6 == probtype )
	{
	  init_bubble(vbx, P_new.array(mfi), S_new.array(mfi, Xvel),
		      S_new.array(mfi, Density), nscal,
		      domain, dx, problo, probhi, IC);
	}
        else if ( 4 == probtype )
	{
	  init_constant_vel_rho(vbx, P_new.array(mfi), S_new.array(mfi, Xvel),
				S_new.array(mfi, Density), nscal,
				domain, dx, problo, probhi, IC);
	}
        else if ( 5 == probtype )
	{
	  init_DoubleShearLayer(vbx, P_new.array(mfi), S_new.array(mfi, Xvel),
				  S_new.array(mfi, Density), nscal,
				  domain, dx, problo, probhi, IC);
	}
        else if ( 7 == probtype )
	{
	  init_Euler(vbx, P_new.array(mfi), S_new.array(mfi, Xvel),
		     S_new.array(mfi, Density), nscal,
		     domain, dx, problo, probhi, IC);
	}
        else if ( 10 == probtype )
	{
	  init_RayleighTaylor(vbx, P_new.array(mfi), S_new.array(mfi, Xvel),
			      S_new.array(mfi, Density), nscal,
			      domain, dx, problo, probhi, IC);
	}
        else if ( 11 == probtype )
	{
	  init_TaylorGreen(vbx, P_new.array(mfi), S_new.array(mfi, Xvel),
			   S_new.array(mfi, Density), nscal,
			   domain, dx, problo, probhi, IC);
	}
	else
        {
            amrex::Abort("NavierStokes::prob_init: unknown probtype");
        };
    }
}

void NavierStokes::init_bubble (Box const& vbx,
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

  bool rise = probtype==6; //HotSpot problem

  amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    Real x = problo[0] + (i - domlo.x + 0.5)*dx[0];
    Real y = problo[1] + (j - domlo.y + 0.5)*dx[1];

    //
    // Fill Velocity
    //
    vel(i,j,k,0) = IC.v_x;
    vel(i,j,k,1) = IC.v_y;

#if (AMREX_SPACEDIM == 3)
    Real z = problo[2] + (k - domlo.z + 0.5)*dx[2];

    vel(i,j,k,2) = IC.v_z;
#endif

    Real dist = std::sqrt( (x-IC.blob_x)*(x-IC.blob_x)
			  + (y-IC.blob_y)*(y-IC.blob_y)
#if (AMREX_SPACEDIM == 3)
			  + (z-IC.blob_z)*(z-IC.blob_z)
#endif
			  );
    //
    // Scalars, ordered as Density, Tracer(s), Temp (if using)
    //

    // Tracers
    scal(i,j,k,1) = dist < IC.blob_radius ? 1.0 : 0.0;
    for ( int nt=2; nt<nscal; nt++)
    {
      scal(i,j,k,nt) = 1.0;
    }

    if ( rise )
    {
      // Density for Hot/less dense bubble rising
      scal(i,j,k,0) = 1.0/IC.density + 0.5*(1.0 - 1.0/IC.density)*(1.0 + std::tanh(40.*(dist - IC.blob_radius)/IC.interface_width));
      //Temp
      scal(i,j,k,nscal-1) = 1/scal(i,j,k,0);
    }
    else
    {
      // Density for dense bubble falling.
      scal(i,j,k,0) = 1.0 + 0.5*(IC.density-1.0)*(1.0-std::tanh(30.*(dist-IC.blob_radius)/IC.interface_width));
    }

  });
}

void NavierStokes::init_constant_vel_rho (Box const& vbx,
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

  amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    Real x = problo[0] + (i - domlo.x + 0.5)*dx[0];
    Real y = problo[1] + (j - domlo.y + 0.5)*dx[1];

    //
    // Fill Velocity
    //
    vel(i,j,k,0) = IC.v_x;
    vel(i,j,k,1) = IC.v_y;

#if (AMREX_SPACEDIM == 3)
    Real z = problo[2] + (k - domlo.z + 0.5)*dx[2];

    vel(i,j,k,2) = IC.v_z;
#endif

    Real dist = std::sqrt( (x-IC.blob_x)*(x-IC.blob_x)
			  + (y-IC.blob_y)*(y-IC.blob_y)
#if (AMREX_SPACEDIM == 3)
			  + (z-IC.blob_z)*(z-IC.blob_z)
#endif
			  );
    //
    // Scalars, ordered as Density, Tracer(s), Temp (if using)
    //
    scal(i,j,k,0) = IC.density;

    // Tracers
    //scal(i,j,k,1) = dist < IC.blob_radius ? 1.0 : 0.0;
    scal(i,j,k,1) = 0.5*(1.0-std::tanh(25.*(dist-IC.blob_radius)/IC.interface_width));
    for ( int nt=2; nt<nscal; nt++)
    {
      scal(i,j,k,nt) = 1.0;
    }
  });
}

void NavierStokes::init_DoubleShearLayer (Box const& vbx,
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

  if ( !(IC.direction != 0 || IC.direction != 1) )
    amrex::Abort("\n    init_DoubleShearLayer: Must set a direction with prob.direction = 0 or 1\n    in the inputs file.  Shear layer along the z-direction not yet written");

  amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    Real x = problo[0] + (i - domlo.x + 0.5)*dx[0];
    Real y = problo[1] + (j - domlo.y + 0.5)*dx[1];
#if (AMREX_SPACEDIM == 3)
    Real z = problo[2] + (k - domlo.z + 0.5)*dx[2];
#endif

    //
    // Fill Velocity
    //
    if ( IC.direction == 1 )
    {
      // shear layer in y-dir
      vel(i,j,k,0) = -.05*std::sin(Pi*y);
      vel(i,j,k,1) = std::tanh(30.*(.5-amrex::Math::abs(x))/IC.interface_width);
    }
    else
    {
      // shear layer in x-dir
      vel(i,j,k,0) = std::tanh(30.*(.5-amrex::Math::abs(y))/IC.interface_width);
      vel(i,j,k,1) = .05*std::sin(Pi*x);
    }

    Real dist = std::sqrt( (x-IC.blob_x)*(x-IC.blob_x)
			  + (y-IC.blob_y)*(y-IC.blob_y)
#if (AMREX_SPACEDIM == 3)
			  + (z-IC.blob_z)*(z-IC.blob_z)
#endif
			  );

    //
    // Scalars, ordered as Density, Tracer(s)
    //
    scal(i,j,k,0) = IC.density;

    // Tracers
    scal(i,j,k,1) = dist < IC.blob_radius ? 1.0 : 0.0;
    for ( int nt=2; nt<nscal; nt++)
    {
      scal(i,j,k,nt) = 1.0;
    }
  });
}

void NavierStokes::init_RayleighTaylor (Box const& vbx,
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

  //
  // Velocity already initialized to 0
  //

  //
  // Scalars, ordered as Density, Tracer(s), Temp (if using)
  //
  const Real Lx    = (probhi[0] - problo[0]);

#if (AMREX_SPACEDIM == 2)
  amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    Real x = problo[0] + (i - domlo.x + 0.5)*dx[0];
    Real y = problo[1] + (j - domlo.y + 0.5)*dx[1];

    const Real pertheight = 0.5 + IC.pertamp*(std::cos(2.0*Pi*x/Lx)
					      + std::cos(2.0*Pi*(Lx-x)/Lx));

    scal(i,j,k,0) = IC.rho_1 + ((IC.rho_2-IC.rho_1)/2.0)*(1.0+std::tanh((y-pertheight)/IC.interface_width));
    scal(i,j,k,1) = IC.tra_1 + ((IC.tra_2-IC.tra_1)/2.0)*(1.0+std::tanh((y-pertheight)/IC.interface_width));
    for ( int nt=2; nt<nscal; nt++)
    {
      scal(i,j,k,nt) = 1.0;
    }

  });

#elif (AMREX_SPACEDIM == 3)

  const Real Ly    = (probhi[1] - problo[1]);
  const Real splitz = 0.5*(problo[2] + probhi[2]);

  Real rn;

  // Create random amplitudes and phases for the perturbation

  //This doens't work for OMP. Just hard-code results below.
  // amrex::InitRandom(111397);
  // rn = amrex::Random();
  // const Real ranampl = 2.*(rn-0.5);

  // rn = amrex::Random();
  // const Real ranphse1 = 2.*Pi*rn;

  // rn = amrex::Random();
  // const Real ranphse2 = 2.*Pi*rn;

  const Real ranampl = 2.*(0.6544437533747718 - 0.5);
  const Real ranphse1 = 2.*Pi*0.1556190326530211;
  const Real ranphse2 = 2.*Pi*0.4196144025537369;

  amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    Real x = problo[0] + (i - domlo.x + 0.5)*dx[0];
    Real y = problo[1] + (j - domlo.y + 0.5)*dx[1];
    Real z = problo[2] + (k - domlo.z + 0.5)*dx[2];

    Real pert = ranampl * std::sin(2.0*Pi*x/Lx + ranphse1 )
                        * std::sin(2.0*Pi*y/Ly + ranphse2 );

    Real pertheight = splitz - IC.pertamp*pert;

    scal(i,j,k,0) = IC.rho_1 + ((IC.rho_2-IC.rho_1)/2.0)*(1.0+std::tanh((z-pertheight)/IC.interface_width));
    scal(i,j,k,1) = IC.tra_1 + ((IC.tra_2-IC.tra_1)/2.0)*(1.0+std::tanh((z-pertheight)/IC.interface_width));
    for ( int nt=2; nt<nscal; nt++)
    {
      scal(i,j,k,nt) = 1.0;
    }
  });

#endif
}

// -----------------------------------------------------------
// This case is an unsteady viscous benchmark for which the
// exact solution in 2D is
//     u(x,y,t) =   Sin(2 Pi x) Cos(2 Pi y) Exp(-2 (2Pi)^2 Nu t)
//     v(x,y,t) = - Cos(2 Pi x) Sin(2 Pi y) Exp(-2 (2Pi)^2 Nu t)
//     p(x,y,t) = - {Cos(4 Pi x) + Cos(4 Pi y)} Exp(-4 (2Pi)^2 Nu t) / 4
// In Exec/benchmarks, there is a tool ViscBench2d.cpp that reads a
// plot file and compares the solution against this exact solution.
// This benchmark was originally derived by G.I. Taylor (Phil. Mag.,
// Vol. 46, No. 274, pp. 671-674, 1923) and Ethier and Steinman
// (Intl. J. Num. Meth. Fluids, Vol. 19, pp. 369-375, 1994) give
// the pressure field.
//
void NavierStokes::init_TaylorGreen (Box const& vbx,
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

  if ( IC.v_x == 0.0 )
    amrex::Abort("NavierStokes::init_TaylorGreen: Must provide prob.velocity_factor. If unsure, prob.velocity_factor = 1.0 is a good choice.");

  amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    Real x = problo[0] + (i - domlo.x + 0.5)*dx[0];
    Real y = problo[1] + (j - domlo.y + 0.5)*dx[1];
#if (AMREX_SPACEDIM == 3)
    Real z = problo[2] + (k - domlo.z + 0.5)*dx[2];
#else
    constexpr Real z = 0.0;
#endif

    //
    // Fill Velocity
    //
    AMREX_D_TERM(vel(i,j,k,0) =  IC.v_x*std::sin(TwoPi*x) * std::cos(TwoPi*y) * std::cos(TwoPi*z);,
		 vel(i,j,k,1) = -IC.v_x*std::cos(TwoPi*x) * std::sin(TwoPi*y) * std::cos(TwoPi*z);,
		 vel(i,j,k,2) = 0.0;);

    //
    // Scalars, ordered as Density, Tracer(s)
    //
    scal(i,j,k,0) = IC.density;

    // The theoretical pressure perturbation from p_0
#if ( AMREX_SPACEDIM == 2 )
    scal(i,j,k,1) = (IC.density*IC.v_x*IC.v_x/4.0)*(cos(2.0*TwoPi*x)+cos(2.0*TwoPi*y));
#else
    scal(i,j,k,1) = (IC.density*IC.v_x*IC.v_x/16.0)*(2.0+cos(2.0*TwoPi*z))*(cos(2.0*TwoPi*x)+cos(2.0*TwoPi*y));
#endif

    // Tracers
    for ( int nt=2; nt<nscal; nt++)
    {
      scal(i,j,k,nt) = 1.0;
    }
  });
}

void NavierStokes::init_Euler (Box const& vbx,
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
    amrex::Abort("NavierStokes::init_Euler: This is a 3D problem, please recompile with DIM=3 in makefile");
#endif

  constexpr Real eps_input=0.05, rho_input=0.15;
  constexpr Real beta_input=15.0, delta_input=0.0333;
  constexpr Real kappa_input=500.0;

  amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    Real x = problo[0] + (i - domlo.x + 0.5)*dx[0] - 0.5;
    Real y = problo[1] + (j - domlo.y + 0.5)*dx[1] - 0.5;
    Real z = problo[2] + (k - domlo.z + 0.5)*dx[2] - 0.5;

    Real r_yz = std::sqrt(y*y+z*z);

    //
    // Fill Velocity
    //
    vel(i,j,k,0) = tanh( (rho_input - r_yz) / delta_input);
    vel(i,j,k,1) = 0.0;
    vel(i,j,k,2) = eps_input * std::exp(-beta_input * (x*x + y*y) );

    //
    // Scalars, ordered as Density, Tracer(s)
    //
    scal(i,j,k,0) = IC.density;
    scal(i,j,k,1) = std::exp( -kappa_input * (rho_input - r_yz)*(rho_input - r_yz) );

    // Additional Tracer, if using
    for ( int nt=2; nt<nscal; nt++)
    {
      scal(i,j,k,nt) = 1.0;
    }
  });
}
