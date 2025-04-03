#include <unistd.h>
#include <AMReX_Geometry.H>
#include <AMReX_Extrapolater.H>
#include <AMReX_ParmParse.H>
#include <AMReX_buildInfo.H>
#include <AMReX_BLProfiler.H>
#include <NavierStokes.H>
#include <NS_util.H>
#include <iamr_constants.H>
#include <NS_LS.H>
#include <NS_kernels.H>

#ifdef AMREX_PARTICLES
#include <DiffusedIB.H>
#endif

#ifdef BL_USE_VELOCITY
#include <AMReX_DataServices.H>
#include <AMReX_AmrData.H>
#endif

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif


using namespace amrex;

int NavierStokes::set_plot_coveredCell_val = 1;

namespace
{
    bool initialized = false;
}

Vector<AMRErrorTag> NavierStokes::errtags;
GpuArray<GpuArray<Real, NavierStokes::NUM_STATE_MAX>, AMREX_SPACEDIM*2> NavierStokes::m_bc_values;

void
NavierStokes::Initialize ()
{
    if (initialized) return;

    NavierStokesBase::Initialize();

    //
    // Set number of state variables.
    //
    NUM_STATE = Density + 1;
    Tracer = NUM_STATE++;
    if (do_trac2)
        Tracer2 = NUM_STATE++;
    if (do_temp)
        Temp = NUM_STATE++;

    //
    // ls related
    //
    if (do_phi)
        phicomp = NUM_STATE++;
    if (verbose) 
        amrex::Print() << "do_phi, phicomp, NUM_STATE " << do_phi << " " << phicomp << " " << NUM_STATE << std::endl;

    //
    // ls related
    //
    // NUM_STATE_MAX is defined in NavierStokes.H
    // to be AMREX_SPACEDIM + 5 (for Density, 2 scalars, Temp, ls)
    AMREX_ALWAYS_ASSERT(NUM_STATE <= NUM_STATE_MAX);

    NUM_SCALARS = NUM_STATE - Density;

    NavierStokes::Initialize_bcs();

    NavierStokes::Initialize_diffusivities();

    amrex::ExecOnFinalize(NavierStokes::Finalize);

    initialized = true;
}

void
NavierStokes::Initialize_bcs ()
{
    //
    // Default BC values
    //
    int ntrac = do_trac2 ? 2 : 1;
    for (OrientationIter face; face; ++face)
    {
      int ori = int(face());
      AMREX_D_TERM(m_bc_values[ori][0] = 0.0;,
           m_bc_values[ori][1] = 0.0;,
           m_bc_values[ori][2] = 0.0;);
      m_bc_values[ori][Density] = 1.0;
      for ( int nc = 0; nc < ntrac; nc++ )
    m_bc_values[ori][Tracer+nc] = 0.0;
      if (do_temp)
	  m_bc_values[ori][Temp] = 1.0;
    //
    // ls related
    //
      if (do_phi)
    m_bc_values[ori][phicomp] = 0.0;
    }

    ParmParse pp("ns");

    //
    // Check for integer BC type specification in inputs file (older style)
    //
    if ( pp.contains("lo_bc") )
    {
      Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
      pp.getarr("lo_bc",lo_bc,0,AMREX_SPACEDIM);
      pp.getarr("hi_bc",hi_bc,0,AMREX_SPACEDIM);
      for (int i = 0; i < AMREX_SPACEDIM; i++)
      {
      phys_bc.setLo(i,lo_bc[i]);
      phys_bc.setHi(i,hi_bc[i]);
      }
    }

    //
    // Read string BC specifications and BC values
    //
    {
      int bc_tmp[2*AMREX_SPACEDIM];

      auto f = [&bc_tmp] (std::string const& bcid, Orientation ori)
      {
      ParmParse pbc(bcid);
      std::string bc_type_in = "null";
      pbc.query("type", bc_type_in);
      std::string bc_type = amrex::toLower(bc_type_in);

      if (bc_type == "no_slip_wall" or bc_type == "nsw"
          or phys_bc.data()[ori] == PhysBCType::noslipwall)
      {
          amrex::Print() << bcid <<" set to no-slip wall.\n";

          bc_tmp[ori] = PhysBCType::noslipwall;

          // Note that m_bc_velocity defaults to 0 above so we are ok if
          //      queryarr finds nothing
          std::vector<Real> v;
          if (pbc.queryarr("velocity", v, 0, AMREX_SPACEDIM)) {
        // Here we make sure that we only use the tangential components
        //      of a specified velocity field -- the wall is not allowed
        //      to move in the normal direction
                v[ori.coordDir()] = 0.0;
                for (int i=0; i<AMREX_SPACEDIM; i++){
          m_bc_values[ori][Xvel+i] = v[i];
                }
          }
      }
      else if (bc_type == "slip_wall" or bc_type == "sw")
      {
          amrex::Print() << bcid <<" set to slip wall.\n";

          bc_tmp[ori] = PhysBCType::slipwall;

          // These values are set by default above -
          //      note that we only actually use the zero value for the normal direction;
          //      the tangential components are set to be first order extrap
          // m_bc_velocity[ori] = {0.0, 0.0, 0.0};
      }
      else if (bc_type == "mass_inflow" or bc_type == "mi"
           or phys_bc.data()[ori] == PhysBCType::inflow)
      {
          amrex::Print() << bcid << " set to mass inflow.\n";

          bc_tmp[ori] = PhysBCType::inflow;

          std::vector<Real> v;
          if (pbc.queryarr("velocity", v, 0, AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; i++){
          m_bc_values[ori][Xvel+i] = v[i];
        }
          }

          pbc.query("density", m_bc_values[ori][Density]);
          pbc.query("tracer", m_bc_values[ori][Tracer]);
          if (do_trac2) {
        if ( pbc.countval("tracer") > 1 )
          amrex::Abort("NavierStokes::Initialize_specific: Please set tracer 2 inflow bc value with it's own entry in inputs file, e.g. xlo.tracer2 = bc_value");
        pbc.query("tracer2", m_bc_values[ori][Tracer2]);
          }
          if (do_temp)
        pbc.query("temp", m_bc_values[ori][Temp]);
      }
      else if (bc_type == "pressure_inflow" or bc_type == "pi")
      {
          amrex::Abort("NavierStokes::Initialize_specific: Pressure inflow boundary condition not yet implemented. If needed for your simulation, please contact us.");

          // amrex::Print() << bcid << " set to pressure inflow.\n";

          // bc_tmp[ori] = PhysBCType::pressure_inflow;

          // pbc.get("pressure", m_bc_pressure[ori]);
      }
      else if (bc_type == "pressure_outflow" or bc_type == "po"
           or phys_bc.data()[ori] == PhysBCType::outflow)
          {
          amrex::Print() << bcid << " set to pressure outflow.\n";

          bc_tmp[ori] = PhysBCType::outflow;

          //pbc.query("pressure", m_bc_pressure[ori]);
          Real tmp = 0.;
          pbc.query("pressure", tmp);

          if ( tmp != 0. )
        amrex::Abort("NavierStokes::Initialize_specific: Pressure outflow boundary condition != 0 not yet implemented. If needed for your simulation, please contact us.");
      }
      else if (bc_type == "symmetry" or bc_type == "sym")
      {
          amrex::Print() << bcid <<" set to symmetry.\n";

          bc_tmp[ori] = PhysBCType::symmetry;
      }
      else
      {
          bc_tmp[ori] = BCType::bogus;
      }

      if ( DefaultGeometry().isPeriodic(ori.coordDir()) ) {
            if (bc_tmp[ori] == BCType::bogus || phys_bc.data()[ori] == PhysBCType::interior) {
          bc_tmp[ori] = PhysBCType::interior;
            } else {
          std::cerr <<  " Wrong BC type for periodic boundary at "
                    << bcid << ". Please correct inputs file."<<std::endl;
          amrex::Abort();
            }
      }

      if ( bc_tmp[ori] == BCType::bogus && phys_bc.data()[ori] == BCType::bogus ) {
        std::cerr <<  " No valid BC type specified for "
                  << bcid  << ". Please correct inputs file."<<std::endl;
        amrex::Abort();
      }

      if ( bc_tmp[ori] != BCType::bogus && phys_bc.data()[ori] != BCType::bogus
           && bc_tmp[ori] != phys_bc.data()[ori] ) {
        std::cerr<<" Multiple conflicting BCs specified for "
                 << bcid << ": "<<bc_tmp[ori]<<", "<<phys_bc.data()[ori]
                 <<". Please correct inputs file."<<std::endl;
        amrex::Abort();
      }
      };

      f("xlo", Orientation(Direction::x,Orientation::low));
      f("xhi", Orientation(Direction::x,Orientation::high));
      f("ylo", Orientation(Direction::y,Orientation::low));
      f("yhi", Orientation(Direction::y,Orientation::high));
#if (AMREX_SPACEDIM == 3)
      f("zlo", Orientation(Direction::z,Orientation::low));
      f("zhi", Orientation(Direction::z,Orientation::high));
#endif

      if ( phys_bc.lo(0) == BCType::bogus )
      {
    //
    // Then valid BC types must be in bc_tmp.
    // Load them into phys_bc.
    //
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
      phys_bc.setLo(dir,bc_tmp[dir]);
      phys_bc.setHi(dir,bc_tmp[dir+AMREX_SPACEDIM]);
    }
      } // else, phys_bc already has valid BC types
    }

    //
    // This checks for RZ and makes sure phys_bc is consistent with that.
    //
    read_geometry();
}

void
NavierStokes::Initialize_diffusivities ()
{
    //
    // Read viscous/diffusive coeffs.
    //
    ParmParse pp("ns");

    const int n_vel_visc_coef   = pp.countval("vel_visc_coef");
    const int n_temp_cond_coef  = pp.countval("temp_cond_coef");
    const int n_scal_diff_coefs = pp.countval("scal_diff_coefs");

    if (n_vel_visc_coef != 1)
        amrex::Abort("NavierStokesBase::Initialize(): Only one vel_visc_coef allowed");

    if (do_temp && n_temp_cond_coef != 1)
        amrex::Abort("NavierStokesBase::Initialize(): Only one temp_cond_coef allowed");

    //
    // ls related
    //
    if (do_phi==0) {
        if (n_scal_diff_coefs+n_temp_cond_coef != NUM_SCALARS-1)
            amrex::Abort("NavierStokesBase::Initialize(): One scal_diff_coef required for each tracer");
    }

    visc_coef.resize(NUM_STATE);
    is_diffusive.resize(NUM_STATE);

    pp.get("vel_visc_coef",visc_coef[0]);
    for (int i = 1; i < AMREX_SPACEDIM; i++)
      visc_coef[i] = visc_coef[0];

    //
    // Here we set the coefficient for density, which does not diffuse.
    //
    visc_coef[Density] = -1;
    //
    // Set the coefficients for the scalars, temperature.
    //
    Vector<Real> scal_diff_coefs(n_scal_diff_coefs);
    pp.getarr("scal_diff_coefs",scal_diff_coefs,0,n_scal_diff_coefs);

    int scalId = Density;

    for (int i = 0; i < n_scal_diff_coefs; i++)
    {
        visc_coef[++scalId] = scal_diff_coefs[i];
    }
    //
    // Set the coefficient for temperature.
    //
    if (do_temp)
    {
        pp.get("temp_cond_coef",visc_coef[++scalId]);
    }
    //
    // ls related
    //
    if (do_phi)
    {
        visc_coef[phicomp] = -1;
    }    
}

void
NavierStokes::Finalize ()
{
    initialized = false;
}

NavierStokes::NavierStokes () = default;

NavierStokes::NavierStokes (Amr&            papa,
                            int             lev,
                            const Geometry& level_geom,
                            const BoxArray& bl,
                            const DistributionMapping& dm,
                            Real            time)
    :
    NavierStokesBase(papa,lev,level_geom,bl,dm,time)
{ }

//
// This function initializes the State and Pressure with data.
//
void
NavierStokes::initData ()
{
    //
    // Initialize the state and the pressure.
    //
    prob_initData();

    //
    // Initialize GradP
    //
    computeGradP(state[Press_Type].curTime());

    //
    // Initialize averaging, if using
    //
    if (avg_interval > 0){
      MultiFab&   Save_new    = get_new_data(Average_Type);
      Save_new.setVal(0.);
    }


#ifdef BL_USE_VELOCITY
    {
      //
      // We want to add the velocity from the supplied plotfile
      // to what we already put into the velocity field via FORT_INITDATA.
      //
      // This code has a few drawbacks.  It assumes that the physical
      // domain size of the current problem is the same as that of the
      // one that generated the pltfile.  It also assumes that the pltfile
      // has at least as many levels (with the same refinement ratios) as does
      // the current problem.  If either of these are false this code is
      // likely to core dump.
      //

      ParmParse pp("ns");

      std::string velocity_plotfile;
      pp.query("velocity_plotfile", velocity_plotfile);

      std::string velocity_plotfile_xvel_name = "x_velocity";
      pp.query("velocity_plotfile_xvel_name", velocity_plotfile_xvel_name);

      Real velocity_plotfile_scale(1.0);
      pp.query("velocity_plotfile_scale",velocity_plotfile_scale);

      if (!velocity_plotfile.empty())
      {
        Print() << "initData: reading data from: " << velocity_plotfile << " ("
                << velocity_plotfile_xvel_name << ")" << '\n';

        DataServices::SetBatchMode();
        Amrvis::FileType fileType(Amrvis::NEWPLT);
        DataServices dataServices(velocity_plotfile, fileType);

        if (!dataServices.AmrDataOk())
            //
            // This calls ParallelDescriptor::EndParallel() and exit()
            //
            DataServices::Dispatch(DataServices::ExitRequest, NULL);

        AmrData&           amrData   = dataServices.AmrDataRef();
        Vector<std::string> plotnames = amrData.PlotVarNames();

        int idX = -1;
        for (int i = 0; i < plotnames.size(); ++i)
            if (plotnames[i] == velocity_plotfile_xvel_name) idX = i;

        if (idX == -1)
           Abort("Could not find velocity fields in supplied velocity_plotfile");
          else
           Print() << "Found " << velocity_plotfile_xvel_name << ", idX = " << idX << '\n';

    MultiFab& S_new = get_new_data(State_Type);
        MultiFab  tmp(S_new.boxArray(), S_new.DistributionMap(), 1, 0);
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        {
        amrData.FillVar(tmp, level, plotnames[idX+i], 0);

        MultiFab::Saxpy(S_new, velocity_plotfile_scale, tmp, 0, Xvel+i, 1, 0);

        amrData.FlushGrids(idX+i);
        }

    Print() << "initData: finished init from velocity_plotfile" << '\n';
      }
    }
#endif /*BL_USE_VELOCITY*/

#ifdef AMREX_USE_EB
    //
    // Perform redistribution on initial fields
    // This changes the input velocity fields
    //
    if (redistribution_type == "StateRedist") {
        InitialRedistribution();
    }

    //
    // Make sure EB covered cell are set, and that it's not zero, as we
    // sometimes divide by rho.
    //
    {
        MultiFab&   S_new    = get_new_data(State_Type);
        EB_set_covered(S_new, COVERED_VAL);

        //
        // In some cases, it may be necessary for the covered cells to
        // contain a value typical (or around the same order of magnitude)
        // to the uncovered cells (e.g. if code to compute variable
        // viscosity fails for COVERED_VAL).
        //
        // set_body_state(S_new);
    }
#endif

    //
    // Make rho MFs with filled ghost cells
    // Not really sure why these are needed as opposed to just filling the
    // the ghost cells in state and using that.
    //
    make_rho_prev_time();
    make_rho_curr_time();

    //
    // Initialize divU and dSdt.
    //
    if (have_divu)
    {
        const Real dt       = 1.0;
        const Real dtin     = -1.0; // Dummy value denotes initialization.
        const Real curTime = state[Divu_Type].curTime();
        MultiFab&  Divu_new = get_new_data(Divu_Type);

        state[State_Type].setTimeLevel(curTime,dt,dt);

        //Make sure something reasonable is in diffn_ec
        calcDiffusivity(curTime);

        calc_divu(curTime,dtin,Divu_new);

        if (have_dsdt)
            get_new_data(Dsdt_Type).setVal(0);
    }

#ifdef AMREX_PARTICLES
    if (do_nspc) {
        initParticleData ();
    }
#endif
}

//
// Build/fill any additional data structures after restart.
//
void
NavierStokes::post_restart ()
{
    NavierStokesBase::post_restart();

    //Get probtype, ub, and shearrate; NS_bcfill.H may expect them.
    ParmParse pp("prob");
    pp.query("probtype",probtype);
    pp.query("ub",ub);
    pp.query("shearrate",shearrate);
}

//
// ADVANCE FUNCTIONS
//

//
// This function ensures that the multifab registers and boundary
// flux registers needed for syncing the composite grid
//
//     u_mac, Vsync, Ssync, rhoavg, fr_adv, fr_visc
//
// are initialized to zero.  In general these quantities
// along with the pressure sync registers (sync_reg) and
// advective velocity registers (mac_reg) are compiled by first
// setting them to the coarse value acquired during a coarse timestep
// and then incrementing in the fine values acquired during the
// subcycled fine timesteps.  This compilation procedure occurs in
// different parts for different quantities
//
// * u_mac is set in predict_velocity and mac_project.
// * fr_adv, fr_visc are set in velocity_advect and scalar_advect
// * Vsync, Ssync are set in subcycled calls to post_timestep
// * mac_reg is set in mac_project
// * sync_reg is set in level_project
// * rhoavg, pavg are set in advance_setup and advance
//
// After these quantities have been compiled during a coarse
// timestep and subcycled fine timesteps.  The post_timestep function
// uses them to sync the fine and coarse levels.  If the coarse level
// is not the base level, post_timestep modifies the next coarsest levels
// registers appropriately.
//
// Note :: There is a little ambiguity as to which level owns the
// boundary flux registers.  The Multifab registers are quantities
// sized by the coarse level BoxArray and belong to the coarse level.
// The fine levels own the boundary registers, since they are sized by
// the boundaries of the fine level BoxArray.
//

//
// Compute a timestep at a level. Return largest safe timestep.
//

Real
NavierStokes::advance (Real time,
                       Real dt,
                       int  iteration,
                       int  ncycle)
{
    BL_PROFILE("NavierStokes::advance()");

    //if (verbose)
    //{
        Print() << "Advancing grids at level " << level
                << " : starting time = "       << time
                << " with dt = "               << dt
                << std::endl;
    //}

    advance_setup(time,dt,iteration,ncycle);

    amrex::Real dt_test = 0.0;
    if (isolver==0) {
        dt_test = advance_semistaggered_twophase_ls(time,dt,iteration,ncycle);
    }
    else if(isolver==1 && do_diffused_ib==1) {
#ifdef AMREX_PARTICLES
        dt_test = advance_semistaggered_fsi_diffusedib(time,dt,iteration,ncycle);
#endif
    }
    else if (isolver==2) { // To be implemented
        dt_test = advance_semistaggered_twophase_phasefield(time,dt,iteration,ncycle);
    }
    else{
        amrex::Abort("Wrong isolver");
    }    

    //
    // Clean up after the predicted value at t^n+1.
    // Estimate new timestep from umac cfl.
    //
    advance_cleanup(iteration,ncycle);

    //if (verbose)
    //{
        Print() << "NavierStokes::advance(): exiting." << std::endl;
        printMaxValues();
    //}

    return dt_test;  // Return estimate of best new timestep.
}

//
// This routine advects the scalars
//

void
NavierStokes::scalar_advection (Real dt,
                                int  fscalar,
                                int  lscalar)
{
    BL_PROFILE("NavierStokes::scalar_advection()");

    if (advect_and_update_scalar) {

        if (verbose) Print() << "... advect scalars\n";
        //
        // Get simulation parameters.
        //
        const int   num_scalars    = lscalar - fscalar + 1;
        const Real  prev_time      = state[State_Type].prevTime();

        // divu
        std::unique_ptr<MultiFab> divu_fp(getDivCond(nghost_force(),prev_time));

        //
        // Start FillPatchIterator block
        //
        MultiFab forcing_term( grids, dmap, num_scalars, nghost_force(),MFInfo(),Factory());

        FillPatchIterator S_fpi(*this,forcing_term,nghost_state(),prev_time,State_Type,fscalar,num_scalars);
        MultiFab& Smf=S_fpi.get_mf();

        // Floor small values of states to be extrapolated
        floor(Smf);

        if ( advection_scheme == "Godunov_PLM" || advection_scheme == "Godunov_PPM" || advection_scheme == "BDS")
        {
            MultiFab visc_terms(grids,dmap,num_scalars,nghost_force(),MFInfo(),Factory());
            FillPatchIterator U_fpi(*this,visc_terms,nghost_state(),prev_time,State_Type,Xvel,AMREX_SPACEDIM);
            const MultiFab& Umf=U_fpi.get_mf();

            {
                std::unique_ptr<MultiFab> dsdt(getDsdt(nghost_force(),prev_time));
                MultiFab::Saxpy(*divu_fp, 0.5*dt, *dsdt, 0, 0, 1, nghost_force());
            }

            // Compute viscous term
            if (be_cn_theta != 1.0)
                getViscTerms(visc_terms,fscalar,num_scalars,prev_time);
            else
                visc_terms.setVal(0.0,1);

    #ifdef _OPENMP
    #pragma omp parallel if (Gpu::notInLaunchRegion())
    #endif
            for (MFIter S_mfi(Smf,TilingIfNotGPU()); S_mfi.isValid(); ++S_mfi)
            {

                // Box for forcing terms
                auto const force_bx = S_mfi.growntilebox(nghost_force());

                if (getForceVerbose)
                {
                    Print() << "---" << '\n' << "C - scalar advection:" << '\n'
                            << " Calling getForce..." << '\n';
                }

                getForce(forcing_term[S_mfi],force_bx,fscalar,num_scalars,
                        prev_time,Umf[S_mfi],Smf[S_mfi],0,S_mfi);

                for (int n=0; n<num_scalars; ++n)
                {
                    auto const& tf    = forcing_term.array(S_mfi,n);
                    auto const& visc  = visc_terms.const_array(S_mfi,n);
                    auto const& rho = Smf.const_array(S_mfi); //Previous time, nghost_state() grow cells filled. It's always true that nghost_state > nghost_force.

            if ( do_temp && n+fscalar==Temp )
            {
            //
            // Solving
            //   dT/dt + U dot del T = ( del dot lambda grad T + H_T ) / (rho c_p)
            // with tforces = H_T/c_p (since it's always density-weighted), and
            // visc = del dot mu grad T, where mu = lambda/c_p
            //
            amrex::ParallelFor(force_bx, [tf, visc, rho]
                    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    { tf(i,j,k) = ( tf(i,j,k) + visc(i,j,k) ) / rho(i,j,k); });
            }
            else
            {
            if (advectionType[fscalar+n] == Conservative)
            {
                //
                // For tracers, Solving
                //   dS/dt + del dot (U S) = del dot beta grad (S/rho) + rho H_q
                // where S = rho q, q is a concentration
                // tforces = rho H_q (since it's always density-weighted)
                // visc = del dot beta grad (S/rho)
                //
                        amrex::ParallelFor(force_bx, [tf, visc]
                        AMREX_GPU_DEVICE (int i, int j, int k ) noexcept
                        { tf(i,j,k) += visc(i,j,k); });
            }
            else
            {
                //
                // Solving
                //   dS/dt + U dot del S = del dot beta grad S + H_q
                // where S = q, q is a concentration
                // tforces = rho H_q (since it's always density-weighted)
                // visc = del dot beta grad S
                //
                        amrex::ParallelFor(force_bx, [tf, visc, rho]
                        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        { tf(i,j,k) = tf(i,j,k) / rho(i,j,k) + visc(i,j,k); });
            }
            }
                }
            }
        }

        ComputeAofs(fscalar, num_scalars, Smf, 0, forcing_term, *divu_fp, false, dt);

    }
}

//
// This subroutine updates the scalars, before the velocity update
// and the level projection
//
// AT this point in time, all we know is psi^n, rho^n+1/2, and the
// general forcing terms at t^n, and after solving in this routine
// viscous forcing at t^n+1/2.  Note, unless more complicated logic
// is invoked earlier, we do not have any estimate of general forcing
// terms at t^n+1/2.
//

void
NavierStokes::scalar_update (Real dt,
                             int  first_scalar,
                             int  last_scalar)
{
    BL_PROFILE("NavierStokes::scalar_update()");

    if (advect_and_update_scalar) {

        if (verbose) Print() << "... update scalars\n";

        scalar_advection_update(dt, first_scalar, last_scalar);

        bool do_any_diffuse = false;
        for (int sigma = first_scalar; sigma <= last_scalar; sigma++)
            if (is_diffusive[sigma]) do_any_diffuse = true;

        if (do_any_diffuse)
        scalar_diffusion_update(dt, first_scalar, last_scalar);

        MultiFab&  S_new     = get_new_data(State_Type);
    //#ifdef AMREX_USE_EB
    //  set_body_state(S_new);
    //#endif
        for (int sigma = first_scalar; sigma <= last_scalar; sigma++)
        {
        if (S_new.contains_nan(sigma,1,0))
        {
            Print() << "New scalar " << sigma << " contains Nans" << '\n';
            exit(0);
        }
        }
    }
}

void
NavierStokes::scalar_diffusion_update (Real dt,
                                       int  first_scalar,
                                       int  last_scalar)
{
    BL_PROFILE("NavierStokes::scalar_diffusion_update()");

    const MultiFab& Rh = get_rho_half_time();

    int ng=1;
    const Real prev_time = state[State_Type].prevTime();
    const Real curr_time = state[State_Type].curTime();

    FillPatch(*this,get_old_data(State_Type),ng,prev_time,State_Type,Density,NUM_SCALARS,Density);
    FillPatch(*this,get_new_data(State_Type),ng,curr_time,State_Type,Density,NUM_SCALARS,Density);

    auto Snc = std::make_unique<MultiFab>();
    auto Snp1c = std::make_unique<MultiFab>();

    if (level > 0) {
      auto& crselev = getLevel(level-1);
      Snc->define(crselev.boxArray(), crselev.DistributionMap(), NUM_STATE, ng, MFInfo(), crselev.Factory());
      FillPatch(crselev,*Snc  ,ng,prev_time,State_Type,0,NUM_STATE);

      Snp1c->define(crselev.boxArray(), crselev.DistributionMap(), NUM_STATE, ng, MFInfo(), crselev.Factory());
      FillPatch(crselev,*Snp1c,ng,curr_time,State_Type,0,NUM_STATE);
    }

    const int nlev = (level ==0 ? 1 : 2);
    Vector<MultiFab*> Sn(nlev,nullptr), Snp1(nlev,nullptr);
    Sn[0]   = &(get_old_data(State_Type));
    Snp1[0] = &(get_new_data(State_Type));

    if (nlev>1) {
      Sn[1]   =  Snc.get() ;
      Snp1[1] =  Snp1c.get() ;
    }

    const Vector<BCRec>& theBCs = AmrLevel::desc_lst[State_Type].getBCs();

    FluxBoxes fb_diffn, fb_diffnp1;
    MultiFab **cmp_diffn = nullptr, **cmp_diffnp1 = nullptr;

    MultiFab *delta_rhs = nullptr;
    MultiFab *alpha = nullptr;
    const int rhsComp = 0, alphaComp = 0, fluxComp  = 0;

    FluxBoxes fb_fluxn  (this);
    FluxBoxes fb_fluxnp1(this);
    MultiFab** fluxn   = fb_fluxn.get();
    MultiFab** fluxnp1 = fb_fluxnp1.get();

    Vector<int> diffuse_comp(1);

    for (int sigma = first_scalar; sigma <= last_scalar; sigma++)
    {
      if (verbose) {
          Print()<<"scalar_diffusion_update "<<sigma<<" of "<<last_scalar<<"\n";
      }

      if (is_diffusive[sigma])
      {
        if (be_cn_theta != 1)
        {
          cmp_diffn = fb_diffn.define(this);
          getDiffusivity(cmp_diffn, prev_time, sigma, 0, 1);
        }

        cmp_diffnp1 = fb_diffnp1.define(this);
        getDiffusivity(cmp_diffnp1, curr_time, sigma, 0, 1);

        diffuse_comp[0] = is_diffusive[sigma];
        const int rho_flag = Diffusion::set_rho_flag(diffusionType[sigma]);

        const bool add_old_time_divFlux = true;

        const int betaComp = 0;
        const int Rho_comp = Density;
        const int bc_comp  = sigma;

        diffusion->diffuse_scalar (Sn, Sn, Snp1, Snp1, sigma, 1, Rho_comp,
                                   prev_time,curr_time,be_cn_theta,Rh,rho_flag,
                                   fluxn,fluxnp1,fluxComp,delta_rhs,rhsComp,
                                   alpha,alphaComp,
                                   cmp_diffn,cmp_diffnp1,betaComp,
                                   crse_ratio,theBCs[bc_comp],geom,
                                   add_old_time_divFlux,
                                   diffuse_comp);

        delete alpha;

        //
        // Increment the viscous flux registers
        //
        if (do_reflux)
        {
            FArrayBox fluxtot;
            for (int d = 0; d < AMREX_SPACEDIM; d++)
            {
                MultiFab fluxes;

                if (level < parent->finestLevel())
                {
                    fluxes.define(fluxn[d]->boxArray(), fluxn[d]->DistributionMap(), 1, 0, MFInfo(), Factory());
                }

                for (MFIter fmfi(*fluxn[d]); fmfi.isValid(); ++fmfi)
                {
                    const Box& ebox = (*fluxn[d])[fmfi].box();

                    fluxtot.resize(ebox,1);
                    Elixir fdata_i = fluxtot.elixir();

                    auto const& ftot = fluxtot.array();
                    auto const& fn   = fluxn[d]->array(fmfi);
                    auto const& fnp1 = fluxnp1[d]->array(fmfi);

                    amrex::ParallelFor(ebox, [ftot, fn, fnp1 ]
                    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        ftot(i,j,k) = fn(i,j,k) + fnp1(i,j,k);
                    });

                    if (level < parent->finestLevel()) {
                        fluxes[fmfi].copy<RunOn::Gpu>(fluxtot);
                    }

                    if (level > 0) {
                        getViscFluxReg().FineAdd(fluxtot,d,fmfi.index(),0,sigma,1,dt,RunOn::Gpu);
                    }
                  } // mfi

                  if (level < parent->finestLevel()) {
                    getLevel(level+1).getViscFluxReg().CrseInit(fluxes,d,0,sigma,1,-dt);
                  }

            } // d
        } // do_reflux

        if (be_cn_theta != 1) {
            fb_diffn.clear();
        }
        fb_diffnp1.clear();

    }//end if(is_diffusive)
    }
}
void
NavierStokes::velocity_diffusion_update (Real dt)
{
    BL_PROFILE("NavierStokes::velocity_diffusion_update()");

    const Real strt_time = ParallelDescriptor::second();
    //
    // Compute the viscous forcing.
    // Do following except at initial iteration.
    //
    if (is_diffusive[Xvel])
    {
        int rho_flag = (do_mom_diff == 0) ? 1 : 3;

        FluxBoxes fb_viscn, fb_viscnp1;
        MultiFab** loc_viscn   = nullptr;
        MultiFab** loc_viscnp1 = nullptr;

        Real viscTime = state[State_Type].prevTime();
        loc_viscn = fb_viscn.define(this);
        getViscosity(loc_viscn, viscTime);

        viscTime = state[State_Type].curTime();
        loc_viscnp1 = fb_viscnp1.define(this);
        getViscosity(loc_viscnp1, viscTime);

        diffusion->diffuse_velocity(dt,be_cn_theta,get_rho_half_time(),rho_flag,
                                    nullptr,loc_viscn,viscn_cc,loc_viscnp1,viscnp1_cc);
    }

    if (verbose)
    {
        Real run_time    = ParallelDescriptor::second() - strt_time;
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        Print() << "NavierStokes:velocity_diffusion_update(): lev: " << level
                << ", time: " << run_time << '\n';
    }
}

void
NavierStokes::sum_integrated_quantities ()
{
    const int finest_level = parent->finestLevel();
    const Real time = state[State_Type].curTime();

    Vector<const MultiFab*> mf(finest_level+1);

    // Container to hold derived data
    Vector<std::unique_ptr<MultiFab>> derived_mf(finest_level+1);
    // Proper container for call to volumeWeightedSum
    Vector<const MultiFab*> derived_mf_ptrs(finest_level+1);

    for (int lev = 0; lev <= finest_level; lev++)
    {
        NavierStokes& ns_level = getLevel(lev);

        mf[lev] = &ns_level.get_new_data(State_Type);
        derived_mf[lev] = ns_level.derive("energy",time,0);
        derived_mf_ptrs[lev] = derived_mf[lev].get();
    }

    Real mass   = volumeWeightedSum(mf, Density,
                                    parent->Geom(), parent->refRatio());
    Real trac   = volumeWeightedSum(mf, Tracer,
                                    parent->Geom(), parent->refRatio());
    Real energy = volumeWeightedSum(derived_mf_ptrs, 0,
                                    parent->Geom(), parent->refRatio());

    Print() << '\n';
    Print().SetPrecision(12) << "TIME= " << time << " MASS= " << mass << '\n';
    Print().SetPrecision(12) << "TIME= " << time << " TRAC= " << trac << '\n';
    Print().SetPrecision(12) << "TIME= " << time << " KINETIC ENERGY= " << energy << '\n';

    //
    // ls related
    //
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream ofs("mass.txt", std::ios::app); // append mode
        // std::ofstream ofs("mass.txt", std::ios::out); // override mode
        if (ofs.is_open())
        {
            amrex::Print(ofs).SetPrecision(12) << time << " " << mass << " " << trac << " " << energy << '\n';
            ofs.close();
        }
        else
        {
            amrex::Abort("Failed to open mass.txt for writing");
        }        
    }
}

void
NavierStokes::setPlotVariables()
{
    AmrLevel::setPlotVariables();
}

void
NavierStokes::writePlotFilePre (const std::string& /*dir*/,
                                std::ostream&  /*os*/)
{
#ifdef AMREX_USE_EB
    if ( set_plot_coveredCell_val )
    {
        for (amrex::Long i =0; i < state.size(); i++)
        {
            auto& sdata = state[i].newData();
        // only cell-centered state data goes into plotfile
            if ( sdata.ixType().cellCentered() ){
                // set covered values for ease of viewing
                EB_set_covered(sdata, 0.0);
            }
        }
    }
#endif
}

void
NavierStokes::writePlotFilePost (const std::string& dir,
                                 std::ostream&  /*os*/)
{
    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        // job_info file with details about the run
        std::ofstream jobInfoFile;
        std::string FullPathJobInfoFile = dir;
        FullPathJobInfoFile += "/job_info";
        jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

        std::string PrettyLine = "===============================================================================\n";
        std::string OtherLine = "--------------------------------------------------------------------------------\n";
        std::string SkipSpace = "        ";


        // job information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Job Information\n";
        jobInfoFile << PrettyLine;

        jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << '\n';
#ifdef _OPENMP
        jobInfoFile << "number of threads:       " << omp_get_max_threads() << '\n';
#endif
        jobInfoFile << "\n\n";

        // plotfile information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Plotfile Information\n";
        jobInfoFile << PrettyLine;

        time_t now = time(nullptr);

        // Convert now to tm struct for local timezone
        tm* localtm = localtime(&now);
        jobInfoFile   << "output data / time: " << asctime(localtm);

        char currentDir[FILENAME_MAX];
        if (getcwd(currentDir, FILENAME_MAX)) {
          jobInfoFile << "output dir:         " << currentDir << '\n';
        }

        jobInfoFile << "\n\n";


        // build information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Build Information\n";
        jobInfoFile << PrettyLine;

        jobInfoFile << "build date:    " << buildInfoGetBuildDate() << '\n';
        jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << '\n';
        jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << '\n';
        jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << '\n';

        jobInfoFile << '\n';

        jobInfoFile << "COMP:          " << buildInfoGetComp() << '\n';
        jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";
        jobInfoFile << "FCOMP:         " << buildInfoGetFcomp() << '\n';
        jobInfoFile << "FCOMP version: " << buildInfoGetFcompVersion() << "\n";

        jobInfoFile << "\n";

        const char* githash1 = buildInfoGetGitHash(1);
        const char* githash2 = buildInfoGetGitHash(2);
        if (strlen(githash1) > 0) {
          jobInfoFile << "IAMR   git hash: " << githash1 << "\n";
        }
        if (strlen(githash2) > 0) {
          jobInfoFile << "AMReX  git hash: " << githash2 << "\n";
        }

        jobInfoFile << "\n\n";


        // runtime parameters
        jobInfoFile << PrettyLine;
        jobInfoFile << " Inputs File Parameters\n";
        jobInfoFile << PrettyLine;

        ParmParse::dumpTable(jobInfoFile, true);

        jobInfoFile.close();

    }

    //
    // Put particles in plotfile?
    // Used in regression testing, but also useful if you want to use
    // AMReX's tool particle_compare without writing a full checkpoint
    //
#ifdef AMREX_PARTICLES
    if (level == 0 && theNSPC() != 0 && particles_in_plotfile)
    {
      theNSPC()->Checkpoint(dir,"Particles");
    }
#endif

#ifdef AMREX_PARTICLES
    if(level == parent->finestLevel()){
        Particles::get_particles()->mContainer->Checkpoint(dir, "particles");
    }
#endif

#ifdef AMREX_USE_EB
    if ( set_plot_coveredCell_val )
    {
        for (amrex::Long i =0; i < state.size(); i++)
        {
            auto& sdata = state[i].newData();
        // only cell-centered state data goes into plotfile
            if ( sdata.ixType().cellCentered() ){
                // put COVERED_VAL back or set_body_state
                EB_set_covered(sdata, COVERED_VAL);
            }
        }
    }
#endif

}

std::unique_ptr<MultiFab>
NavierStokes::derive (const std::string& name,
                      Real               time,
                      int                ngrow)
{
#ifdef AMREX_PARTICLES
    return ParticleDerive(name, time, ngrow);
#else
    return AmrLevel::derive(name, time, ngrow);
#endif
}

void
NavierStokes::derive (const std::string& name,
                      Real               time,
                      MultiFab&          mf,
                      int                dcomp)
{
#ifdef AMREX_PARTICLES
        ParticleDerive(name,time,mf,dcomp);
#else
        AmrLevel::derive(name,time,mf,dcomp);
#endif
}

//
// Ensure state, and pressure are consistent.
//
void
NavierStokes::post_init (Real stop_time)
{

    if (level > 0)
        //
        // Nothing to sync up at level > 0.
        //
        return;

    const int   finest_level = parent->finestLevel();
    Real        dt_init      = 0.0;
    Vector<Real> dt_save(finest_level+1);
    Vector<int>  nc_save(finest_level+1);
    //
    // Ensure state is consistent, i.e. velocity field is non-divergent,
    // pressure & Gradp have been initialized, Coarse levels are fine
    // level averages.
    //
    post_init_state();
    //
    // Estimate the initial timestepping.
    //
    post_init_estDT(dt_init, nc_save, dt_save, stop_time);
    //
    // Initialize the pressure by iterating the initial timestep.
    //
    post_init_press(dt_init, nc_save, dt_save);
    //
    // Compute the initial estimate of conservation.
    //
    if (sum_interval > 0)
        sum_integrated_quantities();

    if (NavierStokesBase::avg_interval > 0)
    {
      NavierStokesBase::time_avg.resize(finest_level+1);
      NavierStokesBase::time_avg_fluct.resize(finest_level+1);
      NavierStokesBase::dt_avg.resize(finest_level+1);
      NavierStokesBase::time_avg[level] = 0.;
      NavierStokesBase::time_avg_fluct[level] = 0.;
      NavierStokesBase::dt_avg[level] = 0.;
      const amrex::Real dt_level = parent->dtLevel(level);
      time_average(NavierStokesBase::time_avg[level], NavierStokesBase::time_avg_fluct[level], NavierStokesBase::dt_avg[level], dt_level);
    }

}

//
// Initialize the pressure by iterating the initial timestep
//

void
NavierStokes::post_init_press (Real&        dt_init,
                               Vector<int>&  nc_save,
                               Vector<Real>& dt_save)
{
    if ( init_iter <= 0 )
    {
      parent->setDtLevel(dt_save);
      parent->setNCycle(nc_save);

      NavierStokes::initial_step = false;

      Print()<< "WARNING! post_init_press(): exiting without doing initial iterations because init_iter <= 0."<<std::endl;

      return;
    }

    const Real strt_time       = state[State_Type].curTime();
    const int  finest_level    = parent->finestLevel();
    NavierStokes::initial_iter = true;

    if (verbose)
    {
        Print() << std::endl
                << "post_init_press(): "
                << "doing initial pressure iterations with dt = "
                << dt_init
                << std::endl;
    }

    //
    // Iterate over the advance function.
    //
    for (int iter = 0; iter < init_iter; iter++)
    {

        if (verbose)
        {
            Print() << std::endl
                    << "post_init_press(): iter = " << iter
                    << std::endl;
        }

        for (int k = 0; k <= finest_level; k++ )
        {
            getLevel(k).advance(strt_time,dt_init,1,1);
        }

        //
        // Constructs a guess at P.
        //
        Vector<MultiFab*> sig(finest_level+1, nullptr);

        for (int k = 0; k <= finest_level; k++)
        {
            sig[k] = &(getLevel(k).get_rho_half_time());
        }

        if (projector)
        {
            projector->initialSyncProject(0,sig,dt_init,
                                          strt_time,have_divu);
        }

        for (int k = finest_level-1; k >= 0; k--)
        {
            getLevel(k).avgDown();
        }

        if (verbose)
        {
            // initSyncProject project d(u)/dt, so new velocity
            // is actually the projected acceleration
            // We don't actually care because initial velocity state will be
            // recovered at the end of each iteration.
            // However, we need to recover u_new from d(u)/dt if we want to print
            // correct diagnostics
            MultiFab& S_new = get_new_data(State_Type);
            MultiFab& S_old = get_old_data(State_Type);
            MultiFab::Xpay(S_new, dt_init, S_old, Xvel, Xvel, AMREX_SPACEDIM, 0);

            Print() << "After sync projection and avgDown:" << std::endl;
            printMaxValues();
        }

        for (int k = 0; k <= finest_level; k++)
        {
            //
            // Reset state variables to initial time, but
            // do not reset pressure variable, only pressure time.
        // do not reset dsdt variable, only dsdt time.
            //
        // The reset of data is ultimately achieved via a swap and swap back to
        // preserve old_data. resetState ends up calling swap(old_data, new_data)
        // but then advance_setup also ends up calling swap(old_data, new_data).
            // Since pressure is not reset, after advance_setup, the next step will
            // progress with p_old==p_new.
            getLevel(k).resetState(strt_time, dt_init, dt_init);
        }

        NavierStokes::initial_iter = false;
    }

    NavierStokes::initial_step = false;
    //
    // Re-instate timestep.
    //
    for (int k = 0; k <= finest_level; k++)
    {
        getLevel(k).setTimeLevel(strt_time,dt_save[k],dt_save[k]);
    }

    parent->setDtLevel(dt_save);
    parent->setNCycle(nc_save);

    // Add space to output if verbose
    if (verbose)
    {
        Print() << std::endl
                << "post_init_press(): exiting after " << init_iter << " iterations"
                << std::endl
                << "After initial iterations: "
                << std::endl;
        printMaxValues();
        Print() << std::endl << std::endl;
    }

}

//
// The Mac Sync correction function
//
void
NavierStokes::mac_sync ()
{
    BL_PROFILE_REGION_START("R::NavierStokes::mac_sync()");
    BL_PROFILE("NavierStokes::mac_sync()");

    if (!do_reflux) return;

    if (verbose)
    {
        Print() << std::endl
                << "mac_sync() on level "<<level
                << std::endl;
    }

    const int  numscal        = NUM_STATE - AMREX_SPACEDIM;
    const Real prev_time      = state[State_Type].prevTime();
    const Real dt             = parent->dtLevel(level);
    MultiFab*  DeltaSsync     = nullptr;// hold (Delta rho)*q for conserved quantities
    // does this have ghosts filled?
    MultiFab&  Rh             = get_rho_half_time();

#ifdef AMREX_USE_EB
    const int nghost = umac_n_grow;
#else
    const int nghost = 0;
#endif

    Array<MultiFab*,AMREX_SPACEDIM> Ucorr;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
      const BoxArray& edgeba = getEdgeBoxArray(idim);

      Ucorr[idim]= new MultiFab(edgeba,dmap,1,nghost,MFInfo(),Factory());
    }

    sync_setup(DeltaSsync);
    //
    // Compute the u_mac for the correction.
    //
    Vector<BCRec> rho_math_bc = fetchBCArray(State_Type,Density,1);
    mac_projector->mac_sync_solve(level,dt,Rh,rho_math_bc[0],fine_ratio,Ucorr);
    //
    // Update coarse grid state by adding correction from mac_sync solve
    // the correction is the advective tendency of the new velocities.
    //
    MultiFab& S_new = get_new_data(State_Type);
    mac_projector->mac_sync_compute(level,Ucorr,Vsync,Ssync,
                    advectionType, prev_time,dt,
                    NUM_STATE,be_cn_theta,
                    do_mom_diff);
    //
    // Delete Ucorr; we're done with it.
    //
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
      delete Ucorr[idim];


    //
    // For all conservative variables Q (other than density)
    // Q is taken as rho*q and we increment sync by -(sync_for_rho)*q
    // (See Pember, et. al., LBNL-41339, Jan. 1989)
    // This increment will be added back after the diffusive sync.
    //
    int iconserved = -1;
    for (int istate = AMREX_SPACEDIM; istate < NUM_STATE; istate++)
    {
      if (istate != Density && advectionType[istate] == Conservative)
      {
    iconserved++;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box&  bx       = mfi.tilebox();
        auto const& rho      = S_new.array(mfi,Density);
        auto const& Snew     = S_new.array(mfi,istate);
        auto const& dSsync   = DeltaSsync->array(mfi);
        auto const& drhosync = Ssync.array(mfi,Density-AMREX_SPACEDIM);
        auto const& ssync    = Ssync.array(mfi,istate-AMREX_SPACEDIM);

        amrex::ParallelFor(bx, [rho, Snew, dSsync, drhosync, ssync, iconserved ]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          dSsync(i,j,k,iconserved) = Snew(i,j,k) * drhosync(i,j,k) / rho(i,j,k);
          ssync(i,j,k) -= dSsync(i,j,k,iconserved);
        });
    }
      }
    }

    if (do_mom_diff == 1)
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        const Box& bx = mfi.tilebox();
        auto const& rho_c    = S_new.array(mfi,Density);
        auto const& vsync    = Vsync.array(mfi,Xvel);
        amrex::ParallelFor(bx, [rho_c, vsync]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
      for (int n = 0; n < AMREX_SPACEDIM; n++) {
        vsync(i,j,k,n) /= rho_c(i,j,k);
      }
    });
      }
    }
    //
    // Compute viscous sync.
    //
    if (is_diffusive[Xvel])
    {
      int rho_flag = (do_mom_diff == 0) ? 1 : 3;

      MultiFab** loc_viscn = nullptr;
      FluxBoxes fb_viscn;

      Real viscTime = state[State_Type].prevTime();
      loc_viscn = fb_viscn.define(this);
      getViscosity(loc_viscn, viscTime);

      diffusion->diffuse_Vsync(Vsync,dt,be_cn_theta,Rh,rho_flag,loc_viscn,0);
    }

    FluxBoxes fb_SC;
    MultiFab** fluxSC        = nullptr;
    bool       any_diffusive = false;
    for (int sigma  = 0; sigma < numscal; sigma++)
      if (is_diffusive[AMREX_SPACEDIM+sigma])
    any_diffusive = true;

    if (any_diffusive) {
      fluxSC = fb_SC.define(this);
    }

    Vector<int> diffuse_comp(1);
    int ng=1;
    const Real curr_time = state[State_Type].curTime();

    // Diffusion solver switches
    // implies that Diff solve does NOT need Sold
    const bool add_old_time_divFlux = false;

    const int nlev = 1;
    Vector<MultiFab*> Snp1(nlev,nullptr);

    MultiFab dSsync(grids,dmap,NUM_STATE,1,MFInfo(),Factory());
    Snp1[0] = &dSsync;

    Vector<MultiFab*> Rhonp1(nlev,nullptr);
    Rhonp1[0] = &(get_new_data(State_Type));
    int Rho_comp = Density;

    FluxBoxes fb_fluxn  (this);
    MultiFab** fluxn   = fb_fluxn.get();

    const Vector<BCRec>& theBCs = AmrLevel::desc_lst[State_Type].getBCs();

    for (int sigma = 0; sigma<numscal; sigma++)
    {
      const int state_ind = AMREX_SPACEDIM + sigma;
      const int rho_flag  = Diffusion::set_rho_flag(diffusionType[state_ind]);

      if (is_diffusive[state_ind])
      {
        Snp1[0]->setVal(0.,state_ind,1,ng);

        FluxBoxes fb_diffnp1;
        MultiFab** cmp_diffnp1=nullptr, **cmp_diffn=nullptr;

        Real diffTime = state[State_Type].curTime();
        cmp_diffnp1 = fb_diffnp1.define(this);
        getDiffusivity(cmp_diffnp1, diffTime, AMREX_SPACEDIM+sigma,0,1);

        int S_comp = state_ind;
        const int num_comp = 1;
        const int fluxComp  = 0;
        MultiFab *delta_rhs = &Ssync;
        int rhsComp = sigma;
        MultiFab *alpha_in = nullptr;
        const int alphaComp = 0;
        int betaComp = 0;

        diffuse_comp[0] = is_diffusive[AMREX_SPACEDIM+sigma];

        diffusion->diffuse_scalar ({},{},Snp1,Rhonp1,
                                   S_comp,num_comp,Rho_comp,
                                   prev_time,curr_time,be_cn_theta,
                                   Rh,rho_flag,
                                   fluxn,fluxSC,fluxComp,
                                   delta_rhs,rhsComp,
                                   alpha_in,alphaComp,
                                   cmp_diffn,cmp_diffnp1,betaComp,
                                   crse_ratio,theBCs[state_ind],geom,
                                   add_old_time_divFlux,diffuse_comp);

        delete alpha_in;

        MultiFab::Copy(Ssync,*Snp1[0],state_ind,sigma,1,0);

        //
        // Increment the viscous flux registers
        //
        if (level > 0)
        {
          for (int d = 0; d < AMREX_SPACEDIM; d++)
          {
             getViscFluxReg().FineAdd(*fluxSC[d],d,0,state_ind,1,dt);
          }
        }
      }
      else // state component not diffusive
      {
      //
      // The following used to be done in mac_sync_compute.  Ssync is
      // the source for a rate of change to S over the time step, so
      // Ssync*dt is the source to the actual sync amount.
      //
        Ssync.mult(dt,sigma,1,Ssync.nGrow());
      }
    }

    //
    // For all conservative variables Q (other than density)
    // increment sync by (sync_for_rho)*q_presync.
    // (See Pember, et. al., LBNL-41339, Jan. 1989)
    //
    iconserved = -1;
    for (int istate = AMREX_SPACEDIM; istate < NUM_STATE; istate++)
    {
      if (istate != Density && advectionType[istate] == Conservative)
      {
    iconserved++;

    // Must multiply by dt to agree with the factor of dt just put
    // into Ssync above.
    MultiFab::Saxpy(Ssync,dt,*DeltaSsync,iconserved,istate-AMREX_SPACEDIM,1,0);
      }
    }
    //
    // Add the sync correction to the state.
    //
    MultiFab::Add(S_new,Ssync,0,AMREX_SPACEDIM,numscal,0);
    //
    // Update rho_ctime after rho is updated with Ssync.
    //
    make_rho_curr_time();

    if (level > 0) incrRhoAvg(Ssync,Density-AMREX_SPACEDIM,1.0);
    //
    // Get boundary conditions.
    //
    const auto N = int(grids.size());

    Vector<int*>         sync_bc(N);
    Vector< Vector<int> > sync_bc_array(N);

    for (int i = 0; i < N; i++)
    {
      sync_bc_array[i] = getBCArray(State_Type,i,Density,numscal);
      sync_bc[i]       = sync_bc_array[i].dataPtr();
    }
    //
    // Interpolate the sync correction to the finer levels,
    //  and update rho_ctime, rhoAvg at those levels.
    //
    IntVect    ratio = IntVect::TheUnitVector();
    const Real mult  = 1.0;
    for (int lev = level+1; lev <= parent->finestLevel(); lev++)
    {
      ratio                     *= parent->refRatio(lev-1);
      NavierStokes&     fine_lev = getLevel(lev);
      const BoxArray& fine_grids = fine_lev.boxArray();
      MultiFab sync_incr(fine_grids,fine_lev.DistributionMap(),numscal,0,MFInfo(),fine_lev.Factory());
      sync_incr.setVal(0.0);

      SyncInterp(Ssync,level,sync_incr,lev,ratio,0,0,
         numscal,1,mult,sync_bc.dataPtr());

      MultiFab& Sf_new = fine_lev.get_new_data(State_Type);
      MultiFab::Add(Sf_new,sync_incr,0,Density,numscal,0);

      fine_lev.make_rho_curr_time();
      fine_lev.incrRhoAvg(sync_incr,Density-AMREX_SPACEDIM,1.0);
    }

    sync_cleanup(DeltaSsync);

    BL_PROFILE_REGION_STOP("R::NavierStokes::mac_sync()");
}

//
// The reflux function
//
void
NavierStokes::reflux ()
{
    if (level == parent->finestLevel())
        return;

    BL_PROFILE("NavierStokes::reflux()");

    AMREX_ASSERT(do_reflux);
    //
    // First do refluxing step.
    //
    auto&         fr_adv  = getAdvFluxReg(level+1);
    FluxRegister& fr_visc = getViscFluxReg(level+1);
    const Real    dt_crse = parent->dtLevel(level);
    //
    // It is important, for do_mom_diff == 0, to do the viscous
    //   refluxing first, since this will be divided by rho_half
    //   before the advective refluxing is added.  In the case of
    //   do_mom_diff == 1, both components of the refluxing will
    //   be divided by rho^(n+1) in level_sync.
    //

    fr_visc.Reflux(Vsync,volume,1.0,0,0,AMREX_SPACEDIM,geom);
    fr_visc.Reflux(Ssync,volume,1.0,AMREX_SPACEDIM,0,NUM_STATE-AMREX_SPACEDIM,geom);

    const MultiFab& Rh = get_rho_half_time();

    if (do_mom_diff == 0)
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(Vsync,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box&  bx      = mfi.tilebox();
         auto const& vsync   = Vsync.array(mfi);
         auto const& rhohalf = Rh.array(mfi);

         amrex::ParallelFor(bx, AMREX_SPACEDIM, [vsync, rhohalf]
         AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
            vsync(i,j,k,n) /= rhohalf(i,j,k);
         });
      }
    }

    for (int istate = AMREX_SPACEDIM; istate < NUM_STATE; istate++)
    {
      if (advectionType[istate] == NonConservative)
      {
          MultiFab::Divide(Ssync,Rh,0,istate-AMREX_SPACEDIM,1,0);
      }
    }

#ifdef AMREX_USE_EB
    fr_adv.Reflux(Vsync,*volfrac,              0, 0,           AMREX_SPACEDIM);
    fr_adv.Reflux(Ssync,*volfrac, AMREX_SPACEDIM, 0, NUM_STATE-AMREX_SPACEDIM);
#else
    fr_adv.Reflux(Vsync,              0, 0,          AMREX_SPACEDIM);
    fr_adv.Reflux(Ssync, AMREX_SPACEDIM, 0,NUM_STATE-AMREX_SPACEDIM);
#endif
    const Real    scale   = 1.0/dt_crse;
    Vsync.mult(scale);
    Ssync.mult(scale);

    const BoxArray& fine_boxes = getLevel(level+1).boxArray();
    //
    // Zero out coarse grid cells which underlie fine grid cells.
    //
    BoxArray baf = fine_boxes;

    baf.coarsen(fine_ratio);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   {
      std::vector< std::pair<int,Box> > isects;
      for (MFIter mfi(Vsync,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.growntilebox();
         auto const& vsync   = Vsync.array(mfi);
         auto const& ssync   = Ssync.array(mfi);
         int nstate          = NUM_STATE;

         baf.intersections(bx,isects);

     for (const auto& is : isects)
     {
        amrex::ParallelFor(is.second, [vsync, ssync, nstate]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               for (int n = 0; n < AMREX_SPACEDIM; n++) {
                  vsync(i,j,k,n) = 0.0;
               }
               for (int n = 0; n < nstate-AMREX_SPACEDIM; n++) {
                  ssync(i,j,k,n) = 0.0;
               }
            });
        }
      }
   }
}

//
// Average fine information from the complete set of state types to coarse.
//

void
NavierStokes::avgDown ()
{
    if (level == parent->finestLevel())
        return;

    auto&   fine_lev = getLevel(level+1);
    //
    // Average down the State and Pressure at the new time.
    //
    avgDown_StatePress();

    //
    // Next average down divu and dSdT at new time.
    //
    if (have_divu)
    {
        MultiFab& Divu_crse = get_new_data(Divu_Type);
        MultiFab& Divu_fine = fine_lev.get_new_data(Divu_Type);

    average_down(Divu_fine, Divu_crse, 0, 1);
    }
    if (have_dsdt)
    {
        MultiFab& Dsdt_crse = get_new_data(Dsdt_Type);
        MultiFab& Dsdt_fine = fine_lev.get_new_data(Dsdt_Type);

    average_down(Dsdt_fine, Dsdt_crse, 0, 1);
    }
}

//
// Default divU is set to zero.
//

void
NavierStokes::calc_divu (Real      time,
                         Real      /*dt*/,
                         MultiFab& divu)
{
    BL_PROFILE("NavierStokes::calc_divu()");

    if (have_divu)
    {
        // Don't think we need this here, but then ghost cells are uninitialized
        // divu.setVal(0);

        if (do_temp && visc_coef[Temp] > 0.0)
        {
            // Compute Div(U) = Div(lambda * Grad(T))/(c_p*rho*T)
            //                = Div(temp_cond_coeff * Grad(T)) / (rho*T)
            // where temp_cond_coef = lambda/cp
            getViscTerms(divu,Temp,1,time);

            const MultiFab&   rhotime = get_rho(time);

            FillPatchIterator temp_fpi(*this,divu,0,time,State_Type,Temp,1);
            MultiFab& tmf = temp_fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for ( MFIter rho_mfi(rhotime,TilingIfNotGPU()); rho_mfi.isValid(); ++rho_mfi)
            {
                const Box&  bx  = rho_mfi.tilebox();
                auto const& div = divu.array(rho_mfi);
                auto const& rho = rhotime.array(rho_mfi);
                auto const& temp = tmf.array(rho_mfi);
#ifdef AMREX_USE_EB
                auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());
                auto const& flagfab = ebfactory.getMultiEBCellFlagFab()[rho_mfi];

                if (flagfab.getType(bx) == FabType::covered)
                {
                    amrex::ParallelFor(bx, [div]
                    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        div( i, j, k ) = COVERED_VAL;
                    });
                }
                else if (flagfab.getType(bx) != FabType::regular)
                {
                    auto vfrac = ebfactory.getVolFrac().const_array(rho_mfi);

                    amrex::ParallelFor(bx, [div, rho, temp, vfrac]
                    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        if ( vfrac(i,j,k) > 0.0 )
                        {
                            div(i,j,k) /= ( rho(i,j,k)*temp(i,j,k) );
                        }
                        else
                        {
                            div(i,j,k) = COVERED_VAL;
                        }
                    });
                }
                else
#endif
                {
                    amrex::ParallelFor(bx, [div, rho, temp]
                    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        div(i,j,k) /= ( rho(i,j,k)*temp(i,j,k) );
                    });
                }
            }
        }
        else
        {
            divu.setVal(0);
        }
    }
}

void
NavierStokes::getViscTerms (MultiFab& visc_terms,
                            int       src_comp,
                            int       ncomp,
                            Real      time)
{
    BL_PROFILE("NavierStokes::getViscTerms()");
    //
    // The logic below for selecting between scalar or tensor solves does
    // not allow for calling NavierStokes::getViscTerms with src_comp=Yvel
    // or Zvel
    //
#ifdef AMREX_DEBUG
    if (src_comp<AMREX_SPACEDIM && (src_comp!=Xvel || ncomp<AMREX_SPACEDIM))
    {
      Print() << "src_comp=" << src_comp << "   ncomp=" << ncomp << '\n';
      Error("must call NavierStokes::getViscTerms with all three velocity components");
    }
#endif
    //
    // Initialize boundary to bogus value so we know if we're using them
    // when we shouldn't
    //
    visc_terms.setBndry(1.e40);

    const int nGrow = visc_terms.nGrow();

    bool diffusive = false;
    //
    // Get Velocity Viscous Terms
    //
    if (src_comp == Xvel && !is_diffusive[Xvel])
    {
        visc_terms.setVal(0.0,0,ncomp,nGrow);
    }
    else if (src_comp == Xvel && is_diffusive[Xvel])
    {
        diffusive = true;

        FluxBoxes fb;
        MultiFab** viscosity = nullptr;

        viscosity = fb.define(this);
        getViscosity(viscosity, time);

        auto whichTime = which_time(State_Type,time);
        AMREX_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

        auto* viscosityCC = (whichTime == AmrOldTime ? viscn_cc : viscnp1_cc);

        diffusion->getTensorViscTerms(visc_terms,time,viscosity,viscosityCC,0);
    }
    //
    // Get Scalar Diffusive Terms
    //
    const int first_scal = (src_comp==Xvel) ? AMREX_SPACEDIM : src_comp;
    const int num_scal   = (src_comp==Xvel) ? ncomp-AMREX_SPACEDIM : ncomp;

    if (num_scal > 0)
    {
        for (int icomp = first_scal; icomp < first_scal+num_scal; icomp++)
        {
            if (is_diffusive[icomp])
            {
                diffusive = true;

                int rho_flag = Diffusion::set_rho_flag(diffusionType[icomp]);

                FluxBoxes fb;
                MultiFab** cmp_diffn = nullptr;

                cmp_diffn = fb.define(this);
                getDiffusivity(cmp_diffn, time, icomp, 0, 1);

                diffusion->getViscTerms(visc_terms,src_comp,icomp,
                                        time,rho_flag,cmp_diffn,0);
            }
            else {
                visc_terms.setVal(0.0,icomp-src_comp,1,nGrow);
            }
        }
    }
    //
    // Ensure consistent grow cells
    //
    if (diffusive && nGrow > 0)
    {
        visc_terms.FillBoundary(0, ncomp, geom.periodicity());
        Extrapolater::FirstOrderExtrap(visc_terms, geom, 0, ncomp);
    }
}

//
// Functions calcViscosity/Diffusivity and getViscosity/Diffusivity are
// for calculating variable viscosity and diffusivity. Here we default to
// constant visc/diff and set the variable viscosity and diffusivity arrays
// to the values in visc_coef and diff_coef.
// For variable viscosity/diffusivity, (per MSD) calcViscosity/Diffusivity
// should compute the transport coefficients at cell centers (or cell centroids
// for EB) and getViscosity/Diffusivity should interpolate those to faces (or
// face-centroids for EB).
//
void
NavierStokes::calcViscosity (const Real time,
                             const Real /*dt*/,
                             const int  /*iteration*/,
                             const int  /*ncycle*/)
{
    if (is_diffusive[Xvel])
    {
        if (visc_coef[Xvel] >= 0.0)
        {
            auto whichTime = which_time(State_Type,time);
            AMREX_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

            auto* visc = (whichTime == AmrOldTime ? viscn_cc : viscnp1_cc);
            visc->setVal(visc_coef[Xvel], 0, visc->nComp(), visc->nGrow());
        }
        else
        {
            Abort("NavierStokes::calcViscosity() : must have velocity visc_coef >= 0.0");
        }
    }
}

void
NavierStokes::calcDiffusivity (const Real time)
{
    //
    // NOTE:  In the diffusivity
    //        arrays, there is an offset since no diffusivity array
    //        is kept for the velocities or the density.  So, the scalar
    //        component Density+1 in the state corresponds to component
    //        0 in the arrays diffn and diffnp1.
    //
    int src_comp = Density+1;
    int ncomp = NUM_STATE - AMREX_SPACEDIM -1;

    const TimeLevel whichTime = which_time(State_Type,time);
    AMREX_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    MultiFab* diff = (whichTime == AmrOldTime ? diffn_cc : diffnp1_cc);
    for (int comp=src_comp; comp<src_comp+ncomp; comp++)
    {
        int diff_comp = comp - Density - 1;

        if (is_diffusive[comp])
        {
            if (visc_coef[comp] >= 0.0)
            {
          diff->setVal(visc_coef[comp], diff_comp, 1, diff->nGrow());
            }
            else
            {
                Abort("NavierStokes::calcDiffusivity() : must have scalar diff_coefs >= 0.0");
            }
        }
    }
}

void
NavierStokes::getViscosity (MultiFab* viscosity[AMREX_SPACEDIM],
                            const Real time)
{
    // //
    // // Select time level to work with (N or N+1)
    // //
    // const TimeLevel whichTime = which_time(State_Type,time);
    // AMREX_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    // MultiFab *visc = (whichTime == AmrOldTime ? viscn_cc : viscnp1_cc);

    // For non-const viscosity, uncomment above and add interp from
    // cell-center/centroid to faces.
    // But here we simply do constant viscosity.

    for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
      viscosity[dir]->setVal(visc_coef[Xvel], 0, viscosity[dir]->nComp(), viscosity[dir]->nGrow());
    }

    if (do_LES)
    {
      FluxBoxes mu_LES(this,1,0);
      MultiFab** mu_LES_mf = mu_LES.get();
      for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
    mu_LES_mf[dir]->setVal(0., 0, mu_LES_mf[dir]->nComp(), mu_LES_mf[dir]->nGrow());
      }

      NavierStokesBase::calc_mut_LES(mu_LES_mf,time);

      for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
    MultiFab::Add(*viscosity[dir], *mu_LES_mf[dir], 0, 0, 1, 0);
      }
    }
}

void
NavierStokes::getDiffusivity (MultiFab* diffusivity[AMREX_SPACEDIM],
                              const Real /*time*/,
                              const int state_comp,
                              const int dst_comp,
                              const int ncomp)
{
    AMREX_ASSERT(state_comp > Density);
    // //
    // // Pick correct component in the diffn/diffnp1 array
    // //
    // int diff_comp = state_comp - Density - 1;
    // //
    // // Select time level to work with (N or N+1)
    // //
    // const TimeLevel whichTime = which_time(State_Type,time);
    // AMREX_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    // MultiFab *diff = (whichTime == AmrOldTime ? diffn_cc : diffnp1_cc);

    // For non-const diffusivity, uncomment above and add interp from
    // cell-center/centroid to faces.
    // But here we simply do constant diffusivity.

    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        for (int n = 0; n < ncomp; n++)
        {
            // Each scalar has it's own diffusivity in the visc_coef array, so must do 1 at a time
            diffusivity[dir]->setVal(visc_coef[state_comp+n], dst_comp+n, 1, diffusivity[dir]->nGrow());
        }
    }
}


//
// ls related
//
Real
NavierStokes::advance_semistaggered_twophase_ls (Real time,
                       Real dt,
                       int  iteration,
                       int  ncycle)
{
    BL_PROFILE("NavierStokes::advance_semistaggered_twophase_ls()");

    //
    // Calculate the time N viscosity and diffusivity
    //   Note: The viscosity and diffusivity at time N+1 are
    //         initialized here to the time N values just to
    //         have something reasonable.
    //
    const Real prev_time = state[State_Type].prevTime();
    const int num_diff = NUM_STATE-AMREX_SPACEDIM-1;

    calcViscosity(prev_time,dt,iteration,ncycle);
    calcDiffusivity(prev_time);
    MultiFab::Copy(*viscnp1_cc, *viscn_cc, 0, 0, 1, viscn_cc->nGrow());
    MultiFab::Copy(*diffnp1_cc, *diffn_cc, 0, 0, num_diff, diffn_cc->nGrow());

    // Add this AFTER advance_setup()
    if (verbose)
    {
        Print() << "NavierStokes::advance_semistaggered_twophase_ls(): before velocity update:"
                << std::endl;
        printMaxValues(false);
    }
    //
    // Compute traced states for normal comp of velocity at half time level.
    // Returns best estimate for new timestep.
    //
    Real dt_test = predict_velocity(dt);
    //
    // Do MAC projection and update edge velocities.
    //
    if (do_mac_proj)
    {
        // To enforce div constraint on coarse-fine boundary, need 1 ghost cell
        int ng_rhs = 1;

        MultiFab mac_rhs(grids,dmap,1,ng_rhs,MFInfo(),Factory());
        create_mac_rhs(mac_rhs,ng_rhs,time,dt);
        MultiFab& S_old = get_old_data(State_Type);
        mac_project(time,dt,S_old,&mac_rhs,umac_n_grow,true);

    } else {
        // Use interpolation from coarse to fill grow cells.
        create_umac_grown(umac_n_grow, nullptr);
    }
    //
    // Advect velocities.
    //
    if (do_mom_diff == 0)
        velocity_advection(dt);
    //
    // Advect scalars.
    //
    const int first_scalar = Density;
    const int last_scalar  = first_scalar + NUM_SCALARS - 1;
    scalar_advection(dt,first_scalar,last_scalar);
    //
    // ls related
    // note: in the above scalar_advection function, we still advect rho.
    // 
    if (do_phi) {
        amrex::Print() << "After scalar_advection " << std::endl;
        // const Real  prev_time = state[State_Type].prevTime();
        // MultiFab& S_old = get_old_data(State_Type);
        // int nScomp = S_old.nComp();
        // fill_allgts(S_old,State_Type,phicomp,1,prev_time);
        // MultiFab::Copy(phi_ptime, S_old, phicomp, 0, 1, S_old.nGrow());

        amrex::Print()<< "scalar_update phi " << std::endl;
        amrex::Print()<< "phicomp " << phicomp << std::endl;
        scalar_update(dt,phicomp,phicomp);

        // amrex::Print()<< std::endl;
        // amrex::Print()<< "6 " << std::endl;

        const Real cur_time = state[State_Type].curTime();
        MultiFab& S_new = get_new_data(State_Type);
        fill_allgts(S_new,State_Type,phicomp,1,cur_time);
        MultiFab::Copy(phi_ctime, S_new, phicomp, 0, 1, S_new.nGrow());

        // amrex::Print()<< "7 " << std::endl;

        // reinitialization
        if (do_reinit == 1 && (parent->levelSteps(0)% lev0step_of_reinit == 0) ){
            amrex::Print() << "parent->levelSteps(0) " << parent->levelSteps(0) << std::endl;
            reinit();
        }

        if (do_mom_diff == 0) {
            // update the rho_ctime and density in S_new
            phi_to_heavi(geom, epsilon, phi_ctime, heaviside); 
            heavi_to_rhoormu(heaviside, rho_w, rho_a, rho_ctime);
            MultiFab::Copy(S_new, rho_ctime, 0, Density, 1, rho_ctime.nGrow());
            // update phi_half
            MultiFab& phi_half_temp = get_phi_half_time();
            // update rho_half
            phi_to_heavi(geom, epsilon, phi_half_temp, heaviside);
            heavi_to_rhoormu(heaviside, rho_w, rho_a, rho_half);

            // update mu_half
            MultiFab outmf_mu_half(grids, dmap, 1, 1, MFInfo(), Factory());
            heavi_to_rhoormu(heaviside, mu_w, mu_a, outmf_mu_half);
            MultiFab::Copy(*viscn_cc,   outmf_mu_half, 0, 0, 1, 1);
            MultiFab::Copy(*viscnp1_cc, outmf_mu_half, 0, 0, 1, 1);
        }
        else {
            //
            // Update Rho.
            //
            scalar_update(dt,first_scalar,first_scalar);
            make_rho_curr_time();

            // update phi_half
            MultiFab& phi_half_temp = get_phi_half_time();
            // update rho_half
            phi_to_heavi(geom, epsilon, phi_half_temp, heaviside);

            // update mu_half
            MultiFab outmf_mu_half(grids, dmap, 1, 1, MFInfo(), Factory());
            heavi_to_rhoormu(heaviside, mu_w, mu_a, outmf_mu_half);
            MultiFab::Copy(*viscn_cc,   outmf_mu_half, 0, 0, 1, 1);
            MultiFab::Copy(*viscnp1_cc, outmf_mu_half, 0, 0, 1, 1);
        }
    }
    else {
        //
        // Update Rho.
        //
        scalar_update(dt,first_scalar,first_scalar);
        make_rho_curr_time();
    }

    if (prescribed_vel)
    {
        BL_ASSERT(do_phi==1);
        dt_test = dt;

        //
        // Create struct to hold initial conditions parameters
        //
        InitialConditions IC;
        // Integer indices of the lower left and upper right corners of the
        // valid region of the entire domain.
        Box const&  domain = geom.Domain();
        auto const&     dx = geom.CellSizeArray();
        // Physical coordinates of the lower left corner of the domain
        auto const& problo = geom.ProbLoArray();
        // Physical coordinates of the upper right corner of the domain
        auto const& probhi = geom.ProbHiArray();

        // Step 1: do the reinitialization
        // which has been done before

        // Step 2: set vel of internal cells in S_new and fill the gts
        MultiFab&  S_new    = get_new_data(State_Type);

        int ncomp = S_new.nComp();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& vbx = mfi.tilebox();
            set_rsv_vel(vbx, S_new.array(mfi, Xvel), domain, dx, problo, probhi, IC, time);
        }

        // Step 3: copy vel in S_new back to S_old
        MultiFab&  S_old    = get_old_data(State_Type);
        MultiFab::Copy(S_old, S_new, 0, 0, ncomp, S_new.nGrow());

    }
    else {
        //
        // Advect momenta after rho^(n+1) has been created.
        //
        if (do_mom_diff == 1)
            velocity_advection(dt);
        //
        // ls related
        // 
        if (do_phi) {
            //
            // Add the advective and other terms to get scalars at t^{n+1} except 
            // the level set function.
            scalar_update(dt,first_scalar+1,phicomp-1);
        }
        else {
            //
            // Add the advective and other terms to get scalars at t^{n+1}.
            //
            scalar_update(dt,first_scalar+1,last_scalar);
        }
        //
        // S appears in rhs of the velocity update, so we better do it now.
        //
        if (have_divu)
        {
            calc_divu(time+dt,dt,get_new_data(Divu_Type));
            if (have_dsdt)
            {
                calc_dsdt(time,dt,get_new_data(Dsdt_Type));
                if (initial_step)
                    MultiFab::Copy(get_old_data(Dsdt_Type),
                                get_new_data(Dsdt_Type),0,0,1,0);
            }
        }
        //
        // Add the advective and other terms to get velocity at t^{n+1}.
        //
        velocity_update(dt);

        //
        // Increment rho average.
        //
        if (!initial_step)
        {
            if (level > 0)
                incrRhoAvg((iteration==ncycle ? 0.5 : 1.0) / Real(ncycle));

            if (verbose)
            {
                Print() << "NavierStokes::advance_semistaggered_twophase_ls(): before nodal projection " << std::endl;
                printMaxVel();
            // New P, Gp get updated in the projector (below). Check old here.
            printMaxGp(false);
            }

            //
            // Do a level project to update the pressure and velocity fields.
            //
            if (projector) {
                const int finest_level = parent->finestLevel();
                int solve_coarse_level = iteration % 2; 
                if (verbose)
                {
                    // Print() << "solve_coarse_level " << solve_coarse_level << std::endl;
                    // Print() << "skip_level_projector " << skip_level_projector << std::endl;
                    // Print() << "level " << level << std::endl;
                    // Print() << "finest_level " << finest_level << std::endl;
                }
                if (skip_level_projector==0 || level==finest_level || solve_coarse_level) {
                    level_projector(dt,time,iteration);
                }
                else {
                    MultiFab& P_old = get_old_data(Press_Type);
                    MultiFab& P_new = get_new_data(Press_Type);
                    // Set P_new to be P_old
                    MultiFab::Copy(P_new,P_old,0,0,1,P_old.nGrow());
                }
            }
            if (level > 0 && iteration == 1)
            p_avg.setVal(0);
        }

#ifdef AMREX_PARTICLES
        if (theNSPC() != 0 and NavierStokes::initial_step != true)
        {
            theNSPC()->AdvectWithUmac(u_mac, level, dt);
        }
#endif
    } // end prescribed_vel

    return dt_test;  // Return estimate of best new timestep.
}

#ifdef AMREX_PARTICLES
Real
NavierStokes::advance_semistaggered_fsi_diffusedib (Real time,
                       Real dt,
                       int  iteration,
                       int  ncycle)
{
    BL_PROFILE("NavierStokes::advance_semistaggered_fsi_diffusedib()");

    //
    // Calculate the time N viscosity and diffusivity
    //   Note: The viscosity and diffusivity at time N+1 are
    //         initialized here to the time N values just to
    //         have something reasonable.
    //
    const Real prev_time = state[State_Type].prevTime();
    const int num_diff = NUM_STATE-AMREX_SPACEDIM-1;

    calcViscosity(prev_time,dt,iteration,ncycle);
    calcDiffusivity(prev_time);
    MultiFab::Copy(*viscnp1_cc, *viscn_cc, 0, 0, 1, viscn_cc->nGrow());
    MultiFab::Copy(*diffnp1_cc, *diffn_cc, 0, 0, num_diff, diffn_cc->nGrow());

    // Add this AFTER advance_setup()
    if (verbose)
    {
        Print() << "NavierStokes::advance_semistaggered_fsi_diffusedib(): before velocity update:"
                << std::endl;
        printMaxValues(false);
    }
    //
    // Compute traced states for normal comp of velocity at half time level.
    // Returns best estimate for new timestep.
    //
    Real dt_test = predict_velocity(dt);
    //
    // Do MAC projection and update edge velocities.
    //
    if (do_mac_proj)
    {
        // To enforce div constraint on coarse-fine boundary, need 1 ghost cell
        int ng_rhs = 1;

        MultiFab mac_rhs(grids,dmap,1,ng_rhs,MFInfo(),Factory());
        create_mac_rhs(mac_rhs,ng_rhs,time,dt);
        MultiFab& S_old = get_old_data(State_Type);
        mac_project(time,dt,S_old,&mac_rhs,umac_n_grow,true);

    } else {
        // Use interpolation from coarse to fill grow cells.
        create_umac_grown(umac_n_grow, nullptr);
    }
    //
    // Advect velocities.
    //
    if (do_mom_diff == 0)
        velocity_advection(dt);

    // Copy fluid density from old to new and set old and new tracer to 0.0
    if (!advect_and_update_scalar)
    {
        MultiFab&  S_new    = get_new_data(State_Type);
        MultiFab&  S_old    = get_old_data(State_Type);
        S_old.setVal(fluid_rho, Density, 1, S_old.nGrow());
        S_new.setVal(fluid_rho, Density, 1, S_new.nGrow());
        S_old.setVal(0.0, Tracer, 1, S_old.nGrow());
        S_new.setVal(0.0, Tracer, 1, S_new.nGrow());
    }
    //
    // Advect scalars.
    //
    const int first_scalar = Density;
    const int last_scalar  = first_scalar + NUM_SCALARS - 1;
    scalar_advection(dt,first_scalar,last_scalar);
    //
    // ls related
    // note: in the above scalar_advection function, we still advect rho.
    // 
    if (do_phi) {
        amrex::Print() << "After scalar_advection " << std::endl;
        // const Real  prev_time = state[State_Type].prevTime();
        // MultiFab& S_old = get_old_data(State_Type);
        // int nScomp = S_old.nComp();
        // fill_allgts(S_old,State_Type,phicomp,1,prev_time);
        // MultiFab::Copy(phi_ptime, S_old, phicomp, 0, 1, S_old.nGrow());

        amrex::Print()<< "scalar_update phi " << std::endl;
        amrex::Print()<< "phicomp " << phicomp << std::endl;
        scalar_update(dt,phicomp,phicomp);

        // amrex::Print()<< std::endl;
        // amrex::Print()<< "6 " << std::endl;

        const Real cur_time = state[State_Type].curTime();
        MultiFab& S_new = get_new_data(State_Type);
        fill_allgts(S_new,State_Type,phicomp,1,cur_time);
        MultiFab::Copy(phi_ctime, S_new, phicomp, 0, 1, S_new.nGrow());

        // amrex::Print()<< "7 " << std::endl;

        // reinitialization
        if (do_reinit == 1 && (parent->levelSteps(0)% lev0step_of_reinit == 0) ){
            amrex::Print() << "parent->levelSteps(0) " << parent->levelSteps(0) << std::endl;
            reinit();
        }

        if (do_mom_diff == 0) {
            // update the rho_ctime and density in S_new
            phi_to_heavi(geom, epsilon, phi_ctime, heaviside); 
            heavi_to_rhoormu(heaviside, rho_w, rho_a, rho_ctime);
            MultiFab::Copy(S_new, rho_ctime, 0, Density, 1, rho_ctime.nGrow());
            // update phi_half
            MultiFab& phi_half_temp = get_phi_half_time();
            // update rho_half
            phi_to_heavi(geom, epsilon, phi_half_temp, heaviside);
            heavi_to_rhoormu(heaviside, rho_w, rho_a, rho_half);

            // update mu_half
            MultiFab outmf_mu_half(grids, dmap, 1, 1, MFInfo(), Factory());
            heavi_to_rhoormu(heaviside, mu_w, mu_a, outmf_mu_half);
            MultiFab::Copy(*viscn_cc,   outmf_mu_half, 0, 0, 1, 1);
            MultiFab::Copy(*viscnp1_cc, outmf_mu_half, 0, 0, 1, 1);
        }
        else {
            //
            // Update Rho.
            //
            scalar_update(dt,first_scalar,first_scalar);
            make_rho_curr_time();

            // update phi_half
            MultiFab& phi_half_temp = get_phi_half_time();
            // update rho_half
            phi_to_heavi(geom, epsilon, phi_half_temp, heaviside);

            // update mu_half
            MultiFab outmf_mu_half(grids, dmap, 1, 1, MFInfo(), Factory());
            heavi_to_rhoormu(heaviside, mu_w, mu_a, outmf_mu_half);
            MultiFab::Copy(*viscn_cc,   outmf_mu_half, 0, 0, 1, 1);
            MultiFab::Copy(*viscnp1_cc, outmf_mu_half, 0, 0, 1, 1);
        }
    }
    else {
        //
        // Update Rho.
        //
        scalar_update(dt,first_scalar,first_scalar);
        make_rho_curr_time();
    }

    if (prescribed_vel)
    {
        dt_test = dt;

        amrex::Print() << "Begin the PVF test " << std::endl;

        // Step 1: initialize the nodal level set function
        // Create struct to hold initial conditions parameters
        //
        InitialConditions IC;
        ParmParse pp("prob");
        pp.query("blob_radius",IC.blob_radius);
        Vector<Real> blob_center(AMREX_SPACEDIM, 0.);
        pp.queryarr("blob_center",blob_center,0,AMREX_SPACEDIM);
        AMREX_D_TERM(IC.blob_x = blob_center[0];,
            IC.blob_y = blob_center[1];,
            IC.blob_z = blob_center[2];);
        // amrex::Print() << "check " << IC.blob_radius << " "
        //                            << IC.blob_x << " "
        //                            << IC.blob_y << " "
        //                            << IC.blob_z << std::endl;

        // Integer indices of the lower left and upper right corners of the
        // valid region of the entire domain.
        Box const&  domain = geom.Domain();
        auto const&     dx = geom.CellSizeArray();
        // Physical coordinates of the lower left corner of the domain
        auto const& problo = geom.ProbLoArray();
        // Physical coordinates of the upper right corner of the domain
        auto const& probhi = geom.ProbHiArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(phi_nodal,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            set_initial_phi_nodal(bx, phi_nodal.array(mfi), domain, dx, problo, probhi, IC, time);
        }

        // Step 2: calculate pvf from the nodal level set function (only for internal cells)
        pvf.setVal(0.0);
        nodal_phi_to_pvf(pvf, phi_nodal);

        // Step 3 (optional): copy for visualization
        MultiFab&  S_new    = get_new_data(State_Type);
        MultiFab::Copy(S_new, pvf, 0, Tracer, 1, pvf.nGrow());

        // Step 4: calculate total volume
        Real vol = dx[0];
        for (int d=1; d<AMREX_SPACEDIM; ++d) {
            vol *= dx[d];
        }
        pvf.mult(vol);
        amrex::Print() << "Volume with numerical integration " << pvf.sum() << std::endl;

        const Real cur_time = state[State_Type].curTime();
        fill_allgts(S_new,State_Type,0,S_new.nComp(),cur_time);
    }
    else {
        //
        // Advect momenta after rho^(n+1) has been created.
        //
        if (do_mom_diff == 1)
            velocity_advection(dt);
        //
        // ls related
        // 
        if (do_phi) {
            //
            // Add the advective and other terms to get scalars at t^{n+1} except 
            // the level set function.
            scalar_update(dt,first_scalar+1,phicomp-1);
        }
        else {
            //
            // Add the advective and other terms to get scalars at t^{n+1}.
            //
            scalar_update(dt,first_scalar+1,last_scalar);
        }
        //
        // S appears in rhs of the velocity update, so we better do it now.
        //
        if (have_divu)
        {
            calc_divu(time+dt,dt,get_new_data(Divu_Type));
            if (have_dsdt)
            {
                calc_dsdt(time,dt,get_new_data(Dsdt_Type));
                if (initial_step)
                    MultiFab::Copy(get_old_data(Dsdt_Type),
                                get_new_data(Dsdt_Type),0,0,1,0);
            }
        }
        //
        // Add the advective and other terms to get velocity at t^{n+1}.
        //
        velocity_update(dt);
        
#ifdef AMREX_PARTICLES
        if (level == Particles::ParticleFinestLevel())//parent->finestLevel())
        {
            MultiFab&  S_new    = get_new_data(State_Type);
            // S_new.setVal(1.0, 0, 1, S_new.nGrow()); // u = 1
            // S_new.setVal(2.0, 1, 1, S_new.nGrow()); // v = 2
            // S_new.setVal(3.0, 2, 1, S_new.nGrow()); // w = 3
            MultiFab EulerForce(S_new.boxArray(), S_new.DistributionMap(), 3, S_new.nGrow());            
            Particles::get_particles()->InteractWithEuler(S_new, EulerForce, dt); // parent->levelSteps(0), time
        }
        //amrex::Abort("Stop here!");
#endif
        //
        // Increment rho average.
        //
        if (!initial_step)
        {
            if (level > 0)
                incrRhoAvg((iteration==ncycle ? 0.5 : 1.0) / Real(ncycle));

            if (verbose)
            {
                Print() << "NavierStokes::advance_semistaggered_fsi_diffusedib(): before nodal projection " << std::endl;
                printMaxVel();
            // New P, Gp get updated in the projector (below). Check old here.
            printMaxGp(false);
            }

            //
            // Do a level project to update the pressure and velocity fields.
            //
            if (projector) {
                const int finest_level = parent->finestLevel();
                int solve_coarse_level = iteration % 2; 
                if (verbose)
                {
                    Print() << "solve_coarse_level " << solve_coarse_level << std::endl;
                    // Print() << "skip_level_projector " << skip_level_projector << std::endl;
                    // Print() << "level " << level << std::endl;
                    // Print() << "finest_level " << finest_level << std::endl;
                }
                if (skip_level_projector==0 || level==finest_level || solve_coarse_level) {
                    level_projector(dt,time,iteration);
                }
                else {
                    MultiFab& P_old = get_old_data(Press_Type);
                    MultiFab& P_new = get_new_data(Press_Type);
                    // Set P_new to be P_old
                    MultiFab::Copy(P_new,P_old,0,0,1,P_old.nGrow());
                }
            }
            if (level > 0 && iteration == 1)
            p_avg.setVal(0);
        }
#ifdef AMREX_PARTICLES
        if (level == Particles::ParticleFinestLevel())//parent->finestLevel())
        {
            MultiFab&  S_new    = get_new_data(State_Type);
            MultiFab&  S_old    = get_old_data(State_Type);
            Particles::get_particles()->UpdateParticles(parent->levelSteps(0), time, S_old, S_new, phi_nodal, pvf, dt);
        }
#endif

#ifdef AMREX_PARTICLES
        if (theNSPC() != 0 and NavierStokes::initial_step != true)
        {
            theNSPC()->AdvectWithUmac(u_mac, level, dt);
        }
#endif
    } // end prescribed_vel

    return dt_test;  // Return estimate of best new timestep.
}
#endif

//
// phase field method
// 
Real
NavierStokes::advance_semistaggered_twophase_phasefield (Real time,
                       Real dt,
                       int  iteration,
                       int  ncycle)
{
    BL_PROFILE("NavierStokes::advance_semistaggered_twophase_phasefield()");

    //
    // Calculate the time N viscosity and diffusivity
    //   Note: The viscosity and diffusivity at time N+1 are
    //         initialized here to the time N values just to
    //         have something reasonable.
    //
    const Real prev_time = state[State_Type].prevTime();
    const int num_diff = NUM_STATE-AMREX_SPACEDIM-1;

    calcViscosity(prev_time,dt,iteration,ncycle);
    calcDiffusivity(prev_time);
    MultiFab::Copy(*viscnp1_cc, *viscn_cc, 0, 0, 1, viscn_cc->nGrow());
    MultiFab::Copy(*diffnp1_cc, *diffn_cc, 0, 0, num_diff, diffn_cc->nGrow());

    // Add this AFTER advance_setup()
    if (verbose)
    {
        Print() << "NavierStokes::advance_semistaggered_twophase_phasefield(): before velocity update:"
                << std::endl;
        printMaxValues(false);
    }
    //
    // Compute traced states for normal comp of velocity at half time level.
    // Returns best estimate for new timestep.
    //
    Real dt_test = predict_velocity(dt);
    //
    // Do MAC projection and update edge velocities.
    //
    if (do_mac_proj)
    {
        // To enforce div constraint on coarse-fine boundary, need 1 ghost cell
        int ng_rhs = 1;

        MultiFab mac_rhs(grids,dmap,1,ng_rhs,MFInfo(),Factory());
        create_mac_rhs(mac_rhs,ng_rhs,time,dt);
        MultiFab& S_old = get_old_data(State_Type);
        mac_project(time,dt,S_old,&mac_rhs,umac_n_grow,true);

    } else {
        // Use interpolation from coarse to fill grow cells.
        create_umac_grown(umac_n_grow, nullptr);
    }
    //
    // Advect velocities.
    //
    if (do_mom_diff == 0)
        velocity_advection(dt);
    //
    // Advect scalars.
    //
    const int first_scalar = Density;
    const int last_scalar  = first_scalar + NUM_SCALARS - 1;
    scalar_advection(dt,first_scalar,last_scalar);
    //
    // pm related
    // note: in the above scalar_advection function, we still advect rho.
    // 
    if (do_phi) {

        // SOLVE AND UPDATE THE PHASE FIELD EQUATION HERE 
        // BY REPLACING THE FOLLOWING LEVEL SET METHOD! - by ZDSJTU

        amrex::Print() << "After scalar_advection " << std::endl;
        // const Real  prev_time = state[State_Type].prevTime();
        // MultiFab& S_old = get_old_data(State_Type);
        // int nScomp = S_old.nComp();
        // fill_allgts(S_old,State_Type,phicomp,1,prev_time);
        // MultiFab::Copy(phi_ptime, S_old, phicomp, 0, 1, S_old.nGrow());

        amrex::Print()<< "scalar_update phi " << std::endl;
        amrex::Print()<< "phicomp " << phicomp << std::endl;
        scalar_update(dt,phicomp,phicomp);

        // amrex::Print()<< std::endl;
        // amrex::Print()<< "6 " << std::endl;

        const Real cur_time = state[State_Type].curTime();
        MultiFab& S_new = get_new_data(State_Type);
        fill_allgts(S_new,State_Type,phicomp,1,cur_time);
        MultiFab::Copy(phi_ctime, S_new, phicomp, 0, 1, S_new.nGrow());

        // amrex::Print()<< "7 " << std::endl;

        // reinitialization
        if (do_reinit == 1 && (parent->levelSteps(0)% lev0step_of_reinit == 0) ){
            amrex::Print() << "parent->levelSteps(0) " << parent->levelSteps(0) << std::endl;
            reinit();
        }

        if (do_mom_diff == 0) {
            // update the rho_ctime and density in S_new
            phi_to_heavi(geom, epsilon, phi_ctime, heaviside); 
            heavi_to_rhoormu(heaviside, rho_w, rho_a, rho_ctime);
            MultiFab::Copy(S_new, rho_ctime, 0, Density, 1, rho_ctime.nGrow());
            // update phi_half
            MultiFab& phi_half_temp = get_phi_half_time();
            // update rho_half
            phi_to_heavi(geom, epsilon, phi_half_temp, heaviside);
            heavi_to_rhoormu(heaviside, rho_w, rho_a, rho_half);

            // update mu_half
            MultiFab outmf_mu_half(grids, dmap, 1, 1, MFInfo(), Factory());
            heavi_to_rhoormu(heaviside, mu_w, mu_a, outmf_mu_half);
            MultiFab::Copy(*viscn_cc,   outmf_mu_half, 0, 0, 1, 1);
            MultiFab::Copy(*viscnp1_cc, outmf_mu_half, 0, 0, 1, 1);
        }
        else {
            amrex::Abort("Only do_mom_diff == 0 is considered now.");
        }
    }
    else {
        //
        // Update Rho.
        //
        amrex::Abort("Rho is not updated by themselves.");
    }

    //
    // Advect momenta after rho^(n+1) has been created.
    //
    if (do_mom_diff == 1)
        velocity_advection(dt);
    //
    // pm related
    // 
    if (do_phi) {
        //
        // Add the advective and other terms to get scalars at t^{n+1} except 
        // the level set function.
        scalar_update(dt,first_scalar+1,phicomp-1);
    }
    //
    // S appears in rhs of the velocity update, so we better do it now.
    //
    if (have_divu)
    {
        calc_divu(time+dt,dt,get_new_data(Divu_Type));
        if (have_dsdt)
        {
            calc_dsdt(time,dt,get_new_data(Dsdt_Type));
            if (initial_step)
                MultiFab::Copy(get_old_data(Dsdt_Type),
                            get_new_data(Dsdt_Type),0,0,1,0);
        }
    }
    //
    // Add the advective and other terms to get velocity at t^{n+1}.
    //
    velocity_update(dt);

    //
    // Increment rho average.
    //
    if (!initial_step)
    {
        if (level > 0)
            incrRhoAvg((iteration==ncycle ? 0.5 : 1.0) / Real(ncycle));

        if (verbose)
        {
            Print() << "NavierStokes::advance_semistaggered_twophase_ls(): before nodal projection " << std::endl;
            printMaxVel();
            // New P, Gp get updated in the projector (below). Check old here.
            printMaxGp(false);
        }

        //
        // Do a level project to update the pressure and velocity fields.
        //
        if (projector) {
            const int finest_level = parent->finestLevel();
            int solve_coarse_level = iteration % 2; 
            if (verbose)
            {
                // Print() << "solve_coarse_level " << solve_coarse_level << std::endl;
                // Print() << "skip_level_projector " << skip_level_projector << std::endl;
                // Print() << "level " << level << std::endl;
                // Print() << "finest_level " << finest_level << std::endl;
            }
            if (skip_level_projector==0 || level==finest_level || solve_coarse_level) {
                level_projector(dt,time,iteration);
            }
            else {
                MultiFab& P_old = get_old_data(Press_Type);
                MultiFab& P_new = get_new_data(Press_Type);
                // Set P_new to be P_old
                MultiFab::Copy(P_new,P_old,0,0,1,P_old.nGrow());
            }
        }
        if (level > 0 && iteration == 1)
        p_avg.setVal(0);
    }

#ifdef AMREX_PARTICLES
    if (theNSPC() != 0 and NavierStokes::initial_step != true)
    {
        theNSPC()->AdvectWithUmac(u_mac, level, dt);
    }
#endif

    return dt_test;  // Return estimate of best new timestep.
}