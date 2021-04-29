#ifndef PROB_INIT_H_
#define PROB_INIT_H_


// This header is included by NavierStokes.H.  These are members of NavierStokes

//
// struct to hold initial conditions parameters
//
struct InitialConditions
{
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

//
// Problems parameters, to be read from inputs file
//
static int probtype;

#endif
