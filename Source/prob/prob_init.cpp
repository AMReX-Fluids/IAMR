#include <NavierStokes.H>
#include <AMReX_ParmParse.H>

using namespace amrex;
 

bool NavierStokes::prob_initialized = false;
int NavierStokes::probtype = -1;
Real NavierStokes::adv_vel = 0.;
int NavierStokes::adv_dir = 0;
Real NavierStokes::denfact = 1.0;

Vector<Real> NavierStokes::blob(AMREX_SPACEDIM, 0.);
Real NavierStokes::radblob = 0.1;

//
// Initialize state and pressure with problem-specific data
//
void NavierStokes::prob_initData ()
{
    //
    // Read problem parameters from inputs file
    //
    if ( !prob_initialized )
    {
      ParmParse pp("prob");

      pp.query("probtype",probtype);
      pp.query("adv_vel",adv_vel);
      pp.query("adv_dir",adv_dir);
      pp.query("denfact",denfact);
      
      pp.queryarr("blob_center",blob,0,AMREX_SPACEDIM);
      pp.query("blob_radius",radblob);

      prob_initialized = true;
    }

    //
    // Fill state and, optionally, pressure
    //
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& P_new = get_new_data(Press_Type);

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

	if (1 == probtype )
        {
	  //
	  // Start from rest, constant density of 1
	  //
	  S_new[mfi].setVal<RunOn::Gpu>(1.0,Density);
	}
        else if (2 == probtype || 6 == probtype )
	{
	  init_bubble(vbx, P_new.array(mfi),
		      S_new.array(mfi, Xvel), S_new.array(mfi, Density),
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
				Box const& domain,
				GpuArray<Real, AMREX_SPACEDIM> const& dx,
				GpuArray<Real, AMREX_SPACEDIM> const& problo,
				GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{
  //
  // set up flow direction
  //
  Real xvel=0.0, yvel=0.0, zvel=0.0;

  if ( adv_dir == 0 )
  {
    xvel = adv_vel;
  }
  else if ( adv_dir == 1 )
  {
    yvel = adv_vel;
  }
#if (AMREX_SPACEDIM == 3)
  else if ( adv_dir == 2 )
  {
    zvel = adv_vel;
  }
#endif
  else
  {
    Abort("init_bubble: unknown adv_dir");
  }

  const auto domlo = amrex::lbound(domain);
  //
  // FIXME - don't think static vars get captured on GPU, so have to make local copies
  //
  const int ntrac = do_trac2 ? 2 : 1;
  bool rise = probtype==6; //HotSpot problem

  amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    Real x = problo[0] + (i - domlo.x + 0.5)*dx[0];
    Real y = problo[1] + (j - domlo.y + 0.5)*dx[1];

    //
    // Fill Velocity
    //
    vel(i,j,k,0) = xvel;
    vel(i,j,k,1) = yvel;
    
#if (AMREX_SPACEDIM == 3)
    Real z = problo[2] + (k - domlo.z + 0.5)*dx[2];

    vel(i,j,k,2) = zvel;
#endif

    Real dist = std::sqrt((x-blob[0])*(x-blob[0]) + (y-blob[1])*(y-blob[1])
#if (AMREX_SPACEDIM == 3)
			  + (z-blob[2])*(z-blob[2])
#endif
			  );
    //
    // Scalars, ordered as Density, Tracer(s), Temp (if using)
    //

    // Tracers
    scal(i,j,k,1) = dist < radblob ? 1.0 : 0.0;
    for ( int nt=2; nt<ntrac; nt++)
    {
      scal(i,j,k,nt) = 1.0;
    }

    if ( rise )
    {
      // Density for Hot/less dense bubble rising
      scal(i,j,k,0) = 1.0/denfact + 0.5*(1.0 - 1.0/denfact) 
	                               *(1.0 + std::tanh(40.*(dist - radblob)));
      //Temp
      scal(i,j,k,ntrac+1) = 1/scal(i,j,k,0);
    }
    else
    {
      // Density for dense bubble falling. 
      scal(i,j,k,0) = 1.0 + 0.5*(denfact-1.0)*(1.0-std::tanh(30.*(dist-radblob)));
    }
    
  });
}
