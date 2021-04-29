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

  amrex::Real interface_width = 1.0;

  // for DoubleShearLayer
  int direction = 0;

  // for Rayleigh-Taylor
  // FIXME? perhaps this should be it's own struct
  amrex::Real rho_1 = 1.0;
  amrex::Real rho_2 = 2.0;
  amrex::Real tra_1 = 0.0;
  amrex::Real tra_2 = 1.0;
  amrex::Real pertamp = 1.0;

  // for HIT
  int inres = 0.0;
  std::string iname;
  bool binfmt = false;
  amrex::Real urms0 = 1.0;
  amrex::Real uin_norm = 1.0;

  amrex::Vector<amrex::Real> xdiff;
  amrex::Vector<amrex::Real> xarray;
  amrex::Real Linput = 0.0;

  amrex::Real* d_uinput = nullptr;
  amrex::Real* d_vinput = nullptr;
  amrex::Real* d_winput = nullptr;
  amrex::Real* d_xarray = nullptr;
  amrex::Real* d_xdiff = nullptr;

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

void init_constant_vel_rho (amrex::Box const& vbx, 
			    amrex::Array4<amrex::Real> const& press,
			    amrex::Array4<amrex::Real> const& vel,
			    amrex::Array4<amrex::Real> const& scal,
			    const int nscal,
			    amrex::Box const& domain,
			    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
			    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
			    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
			    InitialConditions IC);

void init_DoubleShearLayer (amrex::Box const& vbx, 
			    amrex::Array4<amrex::Real> const& press,
			    amrex::Array4<amrex::Real> const& vel,
			    amrex::Array4<amrex::Real> const& scal,
			    const int nscal,
			    amrex::Box const& domain,
			    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
			    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
			    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
			    InitialConditions IC);

void init_RayleighTaylor (amrex::Box const& vbx, 
			  amrex::Array4<amrex::Real> const& press,
			  amrex::Array4<amrex::Real> const& vel,
			  amrex::Array4<amrex::Real> const& scal,
			  const int nscal,
			  amrex::Box const& domain,
			  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
			  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
			  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
			  InitialConditions IC);

void init_TaylorGreen (amrex::Box const& vbx, 
		       amrex::Array4<amrex::Real> const& press,
		       amrex::Array4<amrex::Real> const& vel,
		       amrex::Array4<amrex::Real> const& scal,
		       const int nscal,
		       amrex::Box const& domain,
		       amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
		       amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
		       amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
		       InitialConditions IC);

void init_Euler (amrex::Box const& vbx, 
		 amrex::Array4<amrex::Real> const& press,
		 amrex::Array4<amrex::Real> const& vel,
		 amrex::Array4<amrex::Real> const& scal,
		 const int nscal,
		 amrex::Box const& domain,
		 amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
		 amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
		 amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
		 InitialConditions IC);

void init_HIT (amrex::Box const& vbx,
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
