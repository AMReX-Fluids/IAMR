#include <NavierStokes.H>
#include <AMReX_ParmParse.H>

using namespace amrex;
 
int NavierStokes::probtype = -1;

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
    //
    // Fill state and, optionally, pressure
    //
    MultiFab& P_new = get_new_data(Press_Type);
    MultiFab& S_new = get_new_data(State_Type);
    const int nscal = NUM_STATE-Density;
    
    S_new.setVal(0.0);
    P_new.setVal(0.0);

    Box const&  domain = geom.Domain();
    auto const&     dx = geom.CellSizeArray();
    auto const& problo = geom.ProbLoArray();
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
	  init_channel(vbx, P_new.array(mfi), S_new.array(mfi, Xvel),
		       S_new.array(mfi, Density), nscal,
		       domain, dx, problo, probhi, IC);
	}
	else
        {
            amrex::Abort("prob_init: unknown probtype");
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
      scal(i,j,k,0) = 1.0/IC.density + 0.5*(1.0 - 1.0/IC.density) 
	                               *(1.0 + std::tanh(40.*(dist - IC.blob_radius)));
      //Temp
      scal(i,j,k,nscal-1) = 1/scal(i,j,k,0);
    }
    else
    {
      // Density for dense bubble falling. 
      scal(i,j,k,0) = 1.0 + 0.5*(IC.density-1.0)*(1.0-std::tanh(30.*(dist-IC.blob_radius)));
    }
    
  });
}

void NavierStokes::init_channel (Box const& vbx,
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
    scal(i,j,k,1) = dist < IC.blob_radius ? 1.0 : 0.0;
    for ( int nt=2; nt<nscal; nt++)
    {
      scal(i,j,k,nt) = 1.0;
    }
  });
}
