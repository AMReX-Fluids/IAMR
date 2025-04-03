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

  // Parameters for use in multiple setups
  amrex::Real a = 1.0;
  amrex::Real b = 1.0;
  amrex::Real c = 1.0;

  // for DoubleShearLayer
  int direction = 0;

  // for Rayleigh-Taylor
  amrex::Real rho_1 = 1.0;
  amrex::Real rho_2 = 2.0;
  amrex::Real tra_1 = 0.0;
  amrex::Real tra_2 = 1.0;
  amrex::Real pertamp = 1.0;

  // for convected vortex
  int meanFlowDir = 0;
  amrex::Real meanFlowMag = 0.0;
  amrex::Real forcevort = 6.0;

  //
  // ls related
  // total time for rsv problem
  amrex::Real totalTimeRsv = 4.0;
    
  //
  // ls related
  //
  int do_phi  = 0;
  int Density = 0;
  int phicomp = 0;
  amrex::Real rho_w = 3.0;
  amrex::Real rho_a = 1.0;

};

//
// Problem initialization functions
//
void prob_initData();

static void init_bubble (amrex::Box const& vbx,
          amrex::Array4<amrex::Real> const& press,
          amrex::Array4<amrex::Real> const& vel,
          amrex::Array4<amrex::Real> const& scal,
          int nscal,
          amrex::Box const& domain,
          amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
          amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
          amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
          InitialConditions IC);

//
// ls related
// initialize velocity for the rsv problem
static void init_rsv (amrex::Box const& vbx,
          amrex::Array4<amrex::Real> const& press,
          amrex::Array4<amrex::Real> const& vel,
          amrex::Array4<amrex::Real> const& scal,
          int nscal,
          amrex::Box const& domain,
          amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
          amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
          amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
          InitialConditions IC);

static void init_constant_vel_rho (amrex::Box const& vbx,
                amrex::Array4<amrex::Real> const& press,
                amrex::Array4<amrex::Real> const& vel,
                amrex::Array4<amrex::Real> const& scal,
                int nscal,
                amrex::Box const& domain,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
                InitialConditions IC);

static void init_channel (amrex::Box const& vbx,
                amrex::Array4<amrex::Real> const& press,
                amrex::Array4<amrex::Real> const& vel,
                amrex::Array4<amrex::Real> const& scal,
                int nscal,
                amrex::Box const& domain,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
                InitialConditions IC);

static void init_jump (amrex::Box const& vbx,
                amrex::Array4<amrex::Real> const& press,
                amrex::Array4<amrex::Real> const& vel,
                amrex::Array4<amrex::Real> const& scal,
                int nscal,
                amrex::Box const& domain,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
                InitialConditions IC);

static void init_DoubleShearLayer (amrex::Box const& vbx,
                amrex::Array4<amrex::Real> const& press,
                amrex::Array4<amrex::Real> const& vel,
                amrex::Array4<amrex::Real> const& scal,
                int nscal,
                amrex::Box const& domain,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
                InitialConditions IC);

static void init_RayleighTaylor (amrex::Box const& vbx,
              amrex::Array4<amrex::Real> const& press,
              amrex::Array4<amrex::Real> const& vel,
              amrex::Array4<amrex::Real> const& scal,
              int nscal,
              amrex::Box const& domain,
              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
              InitialConditions IC);

//
// ls related
//
static void init_RayleighTaylor_LS (amrex::Box const& vbx,
              amrex::Array4<amrex::Real> const& press,
              amrex::Array4<amrex::Real> const& vel,
              amrex::Array4<amrex::Real> const& scal,
              int nscal,
              amrex::Box const& domain,
              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
              InitialConditions IC);

static void init_BreakingWave (amrex::Box const& vbx,
              amrex::Array4<amrex::Real> const& press,
              amrex::Array4<amrex::Real> const& vel,
              amrex::Array4<amrex::Real> const& scal,
              int nscal,
              amrex::Box const& domain,
              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
              InitialConditions IC);

static void init_TaylorGreen (amrex::Box const& vbx,
               amrex::Array4<amrex::Real> const& press,
               amrex::Array4<amrex::Real> const& vel,
               amrex::Array4<amrex::Real> const& scal,
               int nscal,
               amrex::Box const& domain,
               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
               InitialConditions IC);

//
// diffused ib
//
static void init_particles (amrex::Box const& vbx,
               amrex::Array4<amrex::Real> const& press,
               amrex::Array4<amrex::Real> const& vel,
               amrex::Array4<amrex::Real> const& scal,
               int nscal,
               amrex::Box const& domain,
               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
               InitialConditions IC);

static void init_Euler (amrex::Box const& vbx,
         amrex::Array4<amrex::Real> const& press,
         amrex::Array4<amrex::Real> const& vel,
         amrex::Array4<amrex::Real> const& scal,
         int nscal,
         amrex::Box const& domain,
         amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
         amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
         amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
         InitialConditions IC);

static void init_ConvectedVortex (amrex::Box const& vbx,
              amrex::Array4<amrex::Real> const& press,
              amrex::Array4<amrex::Real> const& vel,
              amrex::Array4<amrex::Real> const& scal,
              int nscal,
              amrex::Box const& domain,
              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
              InitialConditions IC);

//
// Problems parameters, to be read from inputs file
//
static int probtype;
static amrex::Real ub;
static amrex::Real shearrate;

//
// ls related
// perform prescribed velocity
//
static void set_rsv_vel (amrex::Box const& vbx,
          amrex::Array4<amrex::Real> const& vel,
          amrex::Box const& domain,
          amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
          amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
          amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
          InitialConditions IC,
          amrex::Real time);

//
// diffused ib
// set initial phi nodal
//
static void set_initial_phi_nodal (amrex::Box const& bx,
          amrex::Array4<amrex::Real> const& phi_nodal,
          amrex::Box const& domain,
          amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
          amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
          amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
          InitialConditions IC,
          amrex::Real time);

// Sphere near the channel wall
static void SphereNearWall (amrex::Box const& vbx,
               amrex::Array4<amrex::Real> const& press,
               amrex::Array4<amrex::Real> const& vel,
               amrex::Array4<amrex::Real> const& scal,
               int nscal,
               amrex::Box const& domain,
               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
               InitialConditions IC);

// Falling Sphere
static void FallingSphere (amrex::Box const& vbx,
               amrex::Array4<amrex::Real> const& press,
               amrex::Array4<amrex::Real> const& vel,
               amrex::Array4<amrex::Real> const& scal,
               int nscal,
               amrex::Box const& domain,
               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi,
               InitialConditions IC);

#endif
