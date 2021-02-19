#ifndef PROB_INIT_H_
#define PROB_INIT_H_

// This header is included by NavierStokes.H.  These are members of NavierStokes

//
// struct to hold initial conditions parameters
//
struct InitialConditions
{
  amrex::Real density = 1.0;
  
  amrex::Real v_x = 0.;
  amrex::Real v_y = 0.;
  amrex::Real v_z = 0.;
  
  amrex::Real blob_radius = 0.1;
  
  amrex::Real blob_x = 0.;
  amrex::Real blob_y = 0.;
  amrex::Real blob_z = 0.;
};

//
// Problem initialization functions
//
void prob_initData();

void init_bubble (amrex::Box const& vbx, 
		  amrex::Array4<amrex::Real> const& press,
		  amrex::Array4<amrex::Real> const& vel,
		  amrex::Array4<amrex::Real> const& scal,
		  const int nscal,
		  amrex::Box const& domain,
		  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
		  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
		  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
		  InitialConditions IC);

void init_channel (amrex::Box const& vbx, 
		   amrex::Array4<amrex::Real> const& press,
		   amrex::Array4<amrex::Real> const& vel,
		   amrex::Array4<amrex::Real> const& scal,
		   const int nscal,
		   amrex::Box const& domain,
		   amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
		   amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
		   amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
		   InitialConditions IC);

//FIXME - maybe make use of an enum like this?
//enum struct Prob { HotSpot, bubble };

//
// Problems parameters, to be read from inputs file
//
static int probtype;

#endif
