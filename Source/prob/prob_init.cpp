#include <NavierStokes.H>
#include <AMReX_ParmParse.H>

using namespace amrex;
 
int NavierStokes::probtype = -1;

//
// Initial Conditions parameters
//
namespace IC
{
  Real density_ic = 1.0;
  GpuArray<Real, AMREX_SPACEDIM> velocity_ic {AMREX_D_DECL( 0., 0., 0. )};
  
  Real radblob = 0.1;
  GpuArray<Real, AMREX_SPACEDIM> blob {AMREX_D_DECL( 0., 0., 0. )};
}

//
// Initialize state and pressure with problem-specific data
//
void NavierStokes::prob_initData ()
{
    //
    // Read problem parameters from inputs file
    //
    {
      ParmParse pp("prob");
      
      pp.query("probtype",probtype);
      pp.query("density_ic",IC::density_ic);

      Vector<Real> vel(AMREX_SPACEDIM, 0.);
      pp.queryarr("velocity_ic",vel,0,AMREX_SPACEDIM);
      if (!vel.empty())
      {
	IC::velocity_ic = {AMREX_D_DECL(vel[0], vel[1], vel[2])};
      }
      
      pp.query("blob_radius",IC::radblob);
      Vector<Real> cen(AMREX_SPACEDIM, 0.);
      pp.queryarr("blob_center",cen,0,AMREX_SPACEDIM);
      if (!cen.empty())
      {
	IC::blob = {AMREX_D_DECL(cen[0], cen[1], cen[2])};
      }
    }

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
		      domain, dx, problo, probhi);
	}
        else if ( 4 == probtype )
	{
	  init_channel(vbx, P_new.array(mfi), S_new.array(mfi, Xvel),
		       S_new.array(mfi, Density), nscal,
		       domain, dx, problo, probhi);
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
				GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{
  const auto domlo = amrex::lbound(domain);
  //
  // FIXME - don't think static vars get captured on GPU, so have to make local copies
  //
  bool rise = probtype==6; //HotSpot problem

  amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    Real x = problo[0] + (i - domlo.x + 0.5)*dx[0];
    Real y = problo[1] + (j - domlo.y + 0.5)*dx[1];

    //
    // Fill Velocity
    //
    vel(i,j,k,0) = IC::velocity_ic[0];
    vel(i,j,k,1) = IC::velocity_ic[1];
    
#if (AMREX_SPACEDIM == 3)
    Real z = problo[2] + (k - domlo.z + 0.5)*dx[2];

    vel(i,j,k,2) = IC::velocity_ic[2];
#endif

    Real dist = std::sqrt( (x-IC::blob[0])*(x-IC::blob[0])
			  + (y-IC::blob[1])*(y-IC::blob[1])
#if (AMREX_SPACEDIM == 3)
			  + (z-IC::blob[2])*(z-IC::blob[2])
#endif
			  );
    //
    // Scalars, ordered as Density, Tracer(s), Temp (if using)
    //

    // Tracers
    scal(i,j,k,1) = dist < IC::radblob ? 1.0 : 0.0;
    for ( int nt=2; nt<nscal; nt++)
    {
      scal(i,j,k,nt) = 1.0;
    }

    if ( rise )
    {
      // Density for Hot/less dense bubble rising
      scal(i,j,k,0) = 1.0/IC::density_ic + 0.5*(1.0 - 1.0/IC::density_ic) 
	                               *(1.0 + std::tanh(40.*(dist - IC::radblob)));
      //Temp
      scal(i,j,k,nscal-1) = 1/scal(i,j,k,0);
    }
    else
    {
      // Density for dense bubble falling. 
      scal(i,j,k,0) = 1.0 + 0.5*(IC::density_ic-1.0)*(1.0-std::tanh(30.*(dist-IC::radblob)));
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
				 GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{  
  const auto domlo = amrex::lbound(domain);

  amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    Real x = problo[0] + (i - domlo.x + 0.5)*dx[0];
    Real y = problo[1] + (j - domlo.y + 0.5)*dx[1];

    //
    // Fill Velocity
    //
    vel(i,j,k,0) = IC::velocity_ic[0];
    vel(i,j,k,1) = IC::velocity_ic[1];
    
#if (AMREX_SPACEDIM == 3)
    Real z = problo[2] + (k - domlo.z + 0.5)*dx[2];

    vel(i,j,k,2) = IC::velocity_ic[2];
#endif

    Real dist = std::sqrt( (x-IC::blob[0])*(x-IC::blob[0])
			  + (y-IC::blob[1])*(y-IC::blob[1])
#if (AMREX_SPACEDIM == 3)
			  + (z-IC::blob[2])*(z-IC::blob[2])
#endif
			  );
    //
    // Scalars, ordered as Density, Tracer(s), Temp (if using)
    //
    scal(i,j,k,0) = IC::density_ic;
      
    // Tracers
    scal(i,j,k,1) = dist < IC::radblob ? 1.0 : 0.0;
    for ( int nt=2; nt<nscal; nt++)
    {
      scal(i,j,k,nt) = 1.0;
    }
  });
}
