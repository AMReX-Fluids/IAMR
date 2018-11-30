//
// "Divu_Type" means S, where divergence U = S
// "Dsdt_Type" means pd S/pd t, where S is as above
//
#include <unistd.h>

#include <algorithm>
#include <vector>
#include <cmath>

#include <AMReX_Geometry.H>
#include <AMReX_Extrapolater.H>
#include <AMReX_ParmParse.H>
#include <NavierStokes.H>
#include <AMReX_MultiGrid.H>
#include <NAVIERSTOKES_F.H>
#include <AMReX_BLProfiler.H>
#include <PROB_NS_F.H>

#ifdef BL_USE_VELOCITY
#include <AMReX_DataServices.H>
#include <AMReX_AmrData.H>
#endif

#include <AMReX_buildInfo.H>

using namespace amrex;

namespace
{
    bool initialized = false;
}

void
NavierStokes::variableCleanUp ()
{
    NavierStokesBase::variableCleanUp ();
}

void
NavierStokes::Initialize ()
{
    if (initialized) return;

    NavierStokesBase::Initialize();

    amrex::ExecOnFinalize(NavierStokes::Finalize);

    initialized = true;
}

void
NavierStokes::Finalize ()
{
    initialized = false;
}

NavierStokes::NavierStokes () {}

NavierStokes::NavierStokes (Amr&            papa,
                            int             lev,
                            const Geometry& level_geom,
                            const BoxArray& bl,
                            const DistributionMapping& dm,
                            Real            time)
    :
    NavierStokesBase(papa,lev,level_geom,bl,dm,time)
{ }

NavierStokes::~NavierStokes () { }

//
// This function initializes the State and Pressure with data.
//
void
NavierStokes::initData ()
{
    //
    // Initialize the state and the pressure.
    //
    int         ns       = NUM_STATE - BL_SPACEDIM;
    const Real* dx       = geom.CellSize();
    MultiFab&   S_new    = get_new_data(State_Type);
    MultiFab&   P_new    = get_new_data(Press_Type);
    const Real  cur_time = state[State_Type].curTime();

    for (MFIter snewmfi(S_new); snewmfi.isValid(); ++snewmfi)
    {
        const Box& vbx = snewmfi.validbox();

        BL_ASSERT(grids[snewmfi.index()] == vbx);

        FArrayBox& Sfab = S_new[snewmfi];
        FArrayBox& Pfab = P_new[snewmfi];

	Sfab.setVal(0.0);
        Pfab.setVal(0.0);

        const int  i       = snewmfi.index();
        RealBox    gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
        const int* lo      = vbx.loVect();
        const int* hi      = vbx.hiVect();
        const int* s_lo    = Sfab.loVect();
        const int* s_hi    = Sfab.hiVect();
        const int* p_lo    = Pfab.loVect();
        const int* p_hi    = Pfab.hiVect();

        FORT_INITDATA (&level,&cur_time,lo,hi,&ns,
                       Sfab.dataPtr(Xvel),
                       Sfab.dataPtr(BL_SPACEDIM),
                       ARLIM(s_lo), ARLIM(s_hi),
                       Pfab.dataPtr(),
                       ARLIM(p_lo), ARLIM(p_hi),
                       dx,gridloc.lo(),gridloc.hi() );
    }

#ifdef BL_USE_VELOCITY
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
        amrex::Print() << "initData: reading data from: " << velocity_plotfile << " (" 
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
	  amrex::Abort("Could not find velocity fields in supplied velocity_plotfile");
	else
	  amrex::Print() << "Found " << velocity_plotfile_xvel_name << ", idX = " << idX << '\n';

        MultiFab tmp(S_new.boxArray(), S_new.DistributionMap(), 1, 0);
        for (int i = 0; i < BL_SPACEDIM; i++)
        {
            amrData.FillVar(tmp, level, plotnames[idX+i], 0);
            for (MFIter mfi(tmp); mfi.isValid(); ++mfi)
            {
                FArrayBox& tfab = tmp[mfi];
  	        tfab.mult(velocity_plotfile_scale, 0, 1);
                S_new[mfi].plus(tfab, tfab.box(), 0, Xvel+i, 1);
	    }
            amrData.FlushGrids(idX+i);
        }

	amrex::Print() << "initData: finished init from velocity_plotfile" << '\n';
    }
#endif /*BL_USE_VELOCITY*/

    make_rho_prev_time();
    make_rho_curr_time();
    //
    // Initialize divU and dSdt.
    //
    if (have_divu)
    {
        const Real dt       = 1.0;
        const Real dtin     = -1.0; // Dummy value denotes initialization.
        const Real cur_time = state[Divu_Type].curTime();
        MultiFab&  Divu_new = get_new_data(Divu_Type);

        state[State_Type].setTimeLevel(cur_time,dt,dt);

	if (variable_scal_diff)
	  //Make sure something reasonable is in diffn_cc
	  calcDiffusivity(cur_time);
	
        calc_divu(cur_time,dtin,Divu_new);

        if (have_dsdt)
            get_new_data(Dsdt_Type).setVal(0);
    }

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point) 
    {
        get_new_data(Dpdt_Type).setVal(0);
    }

    is_first_step_after_regrid = false;
    old_intersect_new          = grids;

#ifdef AMREX_PARTICLES
    initParticleData ();
#endif
}

//
// ADVANCE FUNCTIONS
//

//
// This function ensures that the multifab registers and boundary
// flux registers needed for syncing the composite grid
//
//     u_mac, umacG, Vsync, Ssync, rhoavg, fr_adv, fr_visc
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

    if (verbose) {
      amrex::Print() << "Advancing grids at level " << level
		     << " : starting time = "       << time
		     << " with dt = "               << dt << '\n';
    }
    advance_setup(time,dt,iteration,ncycle);
    //
    // Compute traced states for normal comp of velocity at half time level.
    //
    Real dummy   = 0.0;
    Real dt_test = predict_velocity(dt,dummy);
    //
    // Do MAC projection and update edge velocities.
    //
    if (do_mac_proj) 
    {
        MultiFab mac_rhs(grids,dmap,1,0);
        create_mac_rhs(mac_rhs,0,time,dt);
        MultiFab& S_old = get_old_data(State_Type);
        mac_project(time,dt,S_old,&mac_rhs,have_divu,umac_n_grow,true);
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
    // Update Rho.
    //
    scalar_update(dt,first_scalar,first_scalar);

    make_rho_curr_time();
    //
    // Advect momenta after rho^(n+1) has been created.
    //
    if (do_mom_diff == 1) 
        velocity_advection(dt);
    //
    // Add the advective and other terms to get scalars at t^{n+1}.
    //
#ifdef MOREGENGETFORCE
    if (do_scalar_update_in_order)
    {
	for (int iComp=0; iComp<NUM_SCALARS-1; iComp++)
        {
	    int iScal = first_scalar+scalarUpdateOrder[iComp];
	    amrex::Print() << "... ... updating " << desc_lst[0].name(iScal) << '\n';
	    scalar_update(dt,iScal,iScal);
	}
    }
    else
    {
	scalar_update(dt,first_scalar+1,last_scalar);
    }
#else
    scalar_update(dt,first_scalar+1,last_scalar);
#endif
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

        //
        // Do a level project to update the pressure and velocity fields.
        //
        if (projector)
            level_projector(dt,time,iteration);
        if (level > 0 && iteration == 1)
           p_avg.setVal(0);
    }

#ifdef AMREX_PARTICLES
    if (theNSPC() != 0 and NavierStokes::initial_iter != true)
    {
        theNSPC()->AdvectWithUmac(u_mac, level, dt);
    }
#endif
    //
    // Clean up after the predicted value at t^n+1.
    // Estimate new timestep from umac cfl.
    //
    advance_cleanup(iteration,ncycle);

    return dt_test;  // Return estimate of best new timestep.
}

//
// Predict the edge velocities which go into forming u_mac.  This
// function also returns an estimate of dt for use in variable timesteping.
//

Real
NavierStokes::predict_velocity (Real  dt,
                                Real& comp_cfl)
{
    BL_PROFILE("NavierStokes::predict_velocity()");

    if (verbose) amrex::Print() << "... predict edge velocities\n";
    //
    // Get simulation parameters.
    //
    const int   nComp          = BL_SPACEDIM;
    const Real* dx             = geom.CellSize();
    const Real  prev_time      = state[State_Type].prevTime();
    const Real  prev_pres_time = state[Press_Type].prevTime();
    //
    // Compute viscous terms at level n.
    // Ensure reasonable values in 1 grow cell.  Here, do extrap for
    // c-f/phys boundary, since we have no interpolator fn, also,
    // preserve extrap for corners at periodic/non-periodic intersections.
    //
    MultiFab visc_terms(grids,dmap,nComp,1);

    if (be_cn_theta != 1.0)
    {
	getViscTerms(visc_terms,Xvel,nComp,prev_time);
    }
    else
    {
	visc_terms.setVal(0);
    }
    //
    // Set up the timestep estimation.
    //
    Real cflgrid,u_max[3];
    Real cflmax = 1.0e-10;
    comp_cfl    = (level == 0) ? cflmax : comp_cfl;

    FArrayBox tforces;

    Vector<int> bndry[BL_SPACEDIM];

    MultiFab Gp(grids,dmap,BL_SPACEDIM,1);

    getGradP(Gp, prev_pres_time);
    
    FArrayBox null_fab;

    for (FillPatchIterator U_fpi(*this,visc_terms,Godunov::hypgrow(),
                                 prev_time,State_Type,Xvel,BL_SPACEDIM)
#ifdef BOUSSINESQ
             ,S_fpi(*this,visc_terms,1,prev_time,State_Type,Tracer,1);
	 S_fpi.isValid() && U_fpi.isValid();
	 ++S_fpi, ++U_fpi
#else
#ifdef MOREGENGETFORCE
	     , S_fpi(*this,visc_terms,1,prev_time,State_Type,Density,NUM_SCALARS);
	 S_fpi.isValid() && U_fpi.isValid();
	 ++S_fpi, ++U_fpi
#else
         ; U_fpi.isValid();
	 ++U_fpi
#endif
#endif
	)
    {
        const int i = U_fpi.index();

#ifdef BOUSSINESQ
        getForce(tforces,i,1,Xvel,BL_SPACEDIM,prev_time,S_fpi());
#else
#ifdef GENGETFORCE
        getForce(tforces,i,1,Xvel,BL_SPACEDIM,prev_time,rho_ptime[U_fpi]);
#elif MOREGENGETFORCE
	if (getForceVerbose) {
	  amrex::Print() << "---" << '\n' 
			 << "A - Predict velocity:" << '\n'
			 << " Calling getForce..." << '\n';
	}
        getForce(tforces,i,1,Xvel,BL_SPACEDIM,prev_time,U_fpi(),S_fpi(),0);
#else
	getForce(tforces,i,1,Xvel,BL_SPACEDIM,rho_ptime[U_fpi]);
#endif		 
#endif
        //
        // Test velocities, rho and cfl.
        //
        cflgrid  = godunov->test_u_rho(U_fpi(),rho_ptime[U_fpi],grids[i],dx,dt,u_max);
        cflmax   = std::max(cflgrid,cflmax);
        comp_cfl = std::max(cflgrid,comp_cfl);
        //
        // Compute the total forcing.
        //
        godunov->Sum_tf_gp_visc(tforces,visc_terms[U_fpi],Gp[U_fpi],rho_ptime[U_fpi]);

        D_TERM(bndry[0] = getBCArray(State_Type,i,0,1);,
               bndry[1] = getBCArray(State_Type,i,1,1);,
               bndry[2] = getBCArray(State_Type,i,2,1);)

        godunov->Setup(grids[i], dx, dt, 1,
                       null_fab, bndry[0].dataPtr(),
                       null_fab, bndry[1].dataPtr(),
#if (BL_SPACEDIM == 3)                         
                       null_fab, bndry[2].dataPtr(),
#endif
                       U_fpi(), rho_ptime[U_fpi], tforces);

        godunov->ComputeUmac(grids[i], dx, dt, 
                             u_mac[0][U_fpi], bndry[0].dataPtr(),
                             u_mac[1][U_fpi], bndry[1].dataPtr(),
#if (BL_SPACEDIM == 3)
                             u_mac[2][U_fpi], bndry[2].dataPtr(),
#endif
                             U_fpi(), tforces);
    }

    Real tempdt = std::min(change_max,cfl/cflmax);

    ParallelDescriptor::ReduceRealMin(tempdt);

    return dt*tempdt;
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

    if (verbose) amrex::Print() << "... advect scalars\n";
    //
    // Get simulation parameters.
    //
    const int   num_scalars    = lscalar - fscalar + 1;
    const Real* dx             = geom.CellSize();
    const Real  prev_time      = state[State_Type].prevTime();
    const Real  prev_pres_time = state[Press_Type].prevTime();
    //
    // Get the viscous terms.
    //
    MultiFab visc_terms(grids,dmap,num_scalars,1);

    if (be_cn_theta != 1.0) {
        getViscTerms(visc_terms,fscalar,num_scalars,prev_time);
    } else {
        visc_terms.setVal(0,1);
    }
    //
    // Set up the grid loop.
    //
    FArrayBox flux[BL_SPACEDIM], tforces, tvelforces;

    MultiFab Gp, vel_visc_terms, fluxes[BL_SPACEDIM];

    MultiFab* divu_fp = getDivCond(1,prev_time);
    MultiFab* dsdt    = getDsdt(1,prev_time);
    for (MFIter dsdtmfi(*dsdt); dsdtmfi.isValid(); ++dsdtmfi)
    {
        FArrayBox& dsdtfab = (*dsdt)[dsdtmfi];
        dsdtfab.mult(.5*dt);
        (*divu_fp)[dsdtmfi].plus(dsdtfab);
    }
    delete dsdt;

    if (do_reflux && level < parent->finestLevel())
    {
        for (int i = 0; i < BL_SPACEDIM; i++)
        {
            const BoxArray& ba = getEdgeBoxArray(i);
            fluxes[i].define(ba, dmap, num_scalars, 0);
        }
    }

    const int use_forces_in_trans = godunov->useForcesInTrans();

    if (use_forces_in_trans)
    {
        Gp.define(grids,dmap,BL_SPACEDIM,1);

        vel_visc_terms.define(grids,dmap,BL_SPACEDIM,1);

        getGradP(Gp, prev_pres_time);

        if (be_cn_theta != 1.0)
            getViscTerms(vel_visc_terms,Xvel,BL_SPACEDIM,prev_time);
        else
            vel_visc_terms.setVal(0,1);
    }
    Vector<int> state_bc, bndry[BL_SPACEDIM];
    //
    // Compute the advective forcing.
    //
    for (FillPatchIterator U_fpi(*this,visc_terms,Godunov::hypgrow(),prev_time,State_Type,Xvel,BL_SPACEDIM),
#ifdef BOUSSINESQ
             Scal_fpi(*this,visc_terms,Godunov::hypgrow(),prev_time,State_Type,Tracer,1),
#endif
             S_fpi(*this,visc_terms,Godunov::hypgrow(),prev_time,State_Type,fscalar,num_scalars);
         U_fpi.isValid() && S_fpi.isValid();
         ++U_fpi, ++S_fpi)
    {
        const int i = U_fpi.index();

#ifdef BOUSSINESQ
        getForce(tforces,i,1,fscalar,num_scalars,prev_time,Scal_fpi());
#else
#ifdef GENGETFORCE
        getForce(tforces,i,1,fscalar,num_scalars,prev_time,rho_ptime[U_fpi]);
#elif MOREGENGETFORCE
	if (getForceVerbose) {
	  amrex::Print() << "---" << '\n' << "C - scalar advection:" << '\n' 
			 << " Calling getForce..." << '\n';
	}
        getForce(tforces,i,1,fscalar,num_scalars,prev_time,U_fpi(),S_fpi(),0);
#else
        getForce(tforces,i,1,fscalar,num_scalars,rho_ptime[U_fpi]);
#endif		 
#endif		 
        
        if (use_forces_in_trans)
        {
#ifdef BOUSSINESQ
        getForce(tvelforces,i,1,Xvel,BL_SPACEDIM,prev_time,S_fpi());
#else
#ifdef GENGETFORCE
            getForce(tvelforces,i,1,Xvel,BL_SPACEDIM,prev_time,rho_ptime[U_fpi]);
#elif MOREGENGETFORCE
	    if (getForceVerbose) {
	      amrex::Print() << "---" << '\n' 
			     << "D - scalar advection (use_forces_in_trans):" << '\n' 
			     << " Calling getForce..." << '\n';
	    }
            getForce(tvelforces,i,1,Xvel,BL_SPACEDIM,prev_time,U_fpi(),S_fpi(),0);
#else
            getForce(tvelforces,i,1,Xvel,BL_SPACEDIM,rho_ptime[U_fpi]);
#endif		 
#endif
            godunov->Sum_tf_gp_visc(tvelforces,vel_visc_terms[U_fpi],Gp[U_fpi],rho_ptime[U_fpi]);
        }

        D_TERM(bndry[0] = getBCArray(State_Type,i,0,1);,
               bndry[1] = getBCArray(State_Type,i,1,1);,
               bndry[2] = getBCArray(State_Type,i,2,1);)

        godunov->Setup(grids[i], dx, dt, 0,
                       flux[0], bndry[0].dataPtr(),
                       flux[1], bndry[1].dataPtr(),
#if (BL_SPACEDIM == 3)                         
                       flux[2], bndry[2].dataPtr(),
#endif
                       U_fpi(),rho_ptime[U_fpi],tvelforces);

        //
        // Loop over the scalar components.
        //
        FArrayBox& divufab = (*divu_fp)[U_fpi];
        FArrayBox& aofsfab = (*aofs)[U_fpi];
        FArrayBox& rhopfab = rho_ptime[U_fpi];
        FArrayBox& Ufab    = U_fpi();
        FArrayBox& Sfab    = S_fpi();

        D_TERM(FArrayBox& u_mac_fab0 = u_mac[0][U_fpi];,
               FArrayBox& u_mac_fab1 = u_mac[1][U_fpi];,
               FArrayBox& u_mac_fab2 = u_mac[2][U_fpi];);

        for (int comp = 0 ; comp < num_scalars ; comp++)
        {
            int state_ind = fscalar + comp;
            //
            // Compute total forcing.
            //
            int use_conserv_diff = (advectionType[state_ind] == Conservative) ? true : false;

            //AdvectionScheme adv_scheme = FPU;
            AdvectionScheme adv_scheme = PRE_MAC;

            if (adv_scheme == PRE_MAC)
            {
                godunov->Sum_tf_divu_visc(Sfab,tforces,comp,1,visc_terms[U_fpi],
                                          comp,divufab,rhopfab,use_conserv_diff);
            }
            else
            {
                FArrayBox junkDivu(tforces.box(),1);
                junkDivu.setVal(0);
                godunov->Sum_tf_divu_visc(Sfab,tforces,comp,1,visc_terms[U_fpi],
                                          comp,junkDivu,rhopfab,use_conserv_diff);
            }
            //
            // Advect scalar.
            //
            state_bc = getBCArray(State_Type,i,state_ind,1);

            godunov->AdvectState(grids[i], dx, dt, 
                                 area[0][i], u_mac_fab0, flux[0],
                                 area[1][i], u_mac_fab1, flux[1],
#if (BL_SPACEDIM == 3)                        
                                 area[2][i], u_mac_fab2, flux[2],
#endif
                                 Ufab,Sfab,tforces,divufab,comp,
                                 aofsfab,state_ind,use_conserv_diff,
                                 state_ind,state_bc.dataPtr(),adv_scheme,volume[i]);
            if (do_reflux)
            {
                if (level < parent->finestLevel())
                {
                    for (int d = 0; d < BL_SPACEDIM; d++)
                        fluxes[d][U_fpi].copy(flux[d],0,comp,1);
                }
                if (level > 0)
                {
                    for (int d = 0; d < BL_SPACEDIM; d++)
                        advflux_reg->FineAdd(flux[d],d,i,0,state_ind,1,dt);
                }
            }
        }
    }

    delete divu_fp;

    if (do_reflux && level < parent->finestLevel())
    {
        for (int i = 0; i < BL_SPACEDIM; i++)
            getAdvFluxReg(level+1).CrseInit(fluxes[i],i,0,fscalar,num_scalars,-dt);
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

    if (verbose) amrex::Print() << "... update scalars\n";

    scalar_advection_update(dt, first_scalar, last_scalar);

    bool do_any_diffuse = false;
    for (int sigma = first_scalar; sigma <= last_scalar; sigma++)
        if (is_diffusive[sigma]) do_any_diffuse = true;

    if (do_any_diffuse)
      scalar_diffusion_update(dt, first_scalar, last_scalar);

    MultiFab&  S_new     = get_new_data(State_Type);
    for (int sigma = first_scalar; sigma <= last_scalar; sigma++)
    {
       if (S_new.contains_nan(sigma,1,0))
       {
	 amrex::Print() << "New scalar " << sigma << " contains Nans" << '\n';
	 exit(0);
       }
    }
}

void
NavierStokes::scalar_diffusion_update (Real dt,
                                       int  first_scalar,
                                       int  last_scalar)
{
    BL_PROFILE("NavierStokes::scalar_diffusion_update()");

    FluxBoxes fb_SCn  (this);
    FluxBoxes fb_SCnp1(this);

    MultiFab** fluxSCn   = fb_SCn.get();
    MultiFab** fluxSCnp1 = fb_SCnp1.get();

    const MultiFab& Rh = get_rho_half_time();

    for (int sigma = first_scalar; sigma <= last_scalar; sigma++)
    {
        if (is_diffusive[sigma])
        {
            int        rho_flag    = 0;
            MultiFab*  delta_rhs   = 0;
            MultiFab*  alpha       = 0;
            MultiFab** cmp_diffn   = 0;
            MultiFab** cmp_diffnp1 = 0;

	    FluxBoxes fb_diffn, fb_diffnp1;

            if (variable_scal_diff)
            {
                Real diffTime = state[State_Type].prevTime();
		cmp_diffn = fb_diffn.define(this);
                getDiffusivity(cmp_diffn, diffTime, sigma, 0, 1);

                diffTime = state[State_Type].curTime();
		cmp_diffnp1 = fb_diffnp1.define(this);
                getDiffusivity(cmp_diffnp1, diffTime, sigma, 0, 1);
            }

            diffuse_scalar_setup(sigma, rho_flag);

            const int betaComp = 0, rhsComp = 0, alphaComp = 0, fluxComp  = 0;

            diffusion->diffuse_scalar(dt,sigma,be_cn_theta,Rh,
                                      rho_flag,fluxSCn,fluxSCnp1,fluxComp,delta_rhs,
                                      rhsComp,alpha,alphaComp,cmp_diffn,cmp_diffnp1,betaComp);

            delete delta_rhs;
            delete alpha;
            //
            // Increment the viscous flux registers
            //
            if (do_reflux)
            {
                FArrayBox fluxtot;

                for (int d = 0; d < BL_SPACEDIM; d++)
                {
                    MultiFab fluxes;

                    if (level < parent->finestLevel()) {
                        fluxes.define(fluxSCn[d]->boxArray(), fluxSCn[d]->DistributionMap(), 1, 0);
                    }

                    for (MFIter fmfi(*fluxSCn[d]); fmfi.isValid(); ++fmfi)
                    {
                        const Box& ebox = (*fluxSCn[d])[fmfi].box();

                        fluxtot.resize(ebox,1);
                        fluxtot.copy((*fluxSCn[d])[fmfi],ebox,0,ebox,0,1);
                        fluxtot.plus((*fluxSCnp1[d])[fmfi],ebox,0,0,1);

                        if (level < parent->finestLevel())
                            fluxes[fmfi].copy(fluxtot);

                        if (level > 0)
                            getViscFluxReg().FineAdd(fluxtot,d,fmfi.index(),0,sigma,1,dt);
                    }

                    if (level < parent->finestLevel())
                        getLevel(level+1).getViscFluxReg().CrseInit(fluxes,d,0,sigma,1,-dt);
                }
            }
        }
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

        MultiFab* delta_rhs = 0;
        if (S_in_vel_diffusion && have_divu)
        {
            delta_rhs = new MultiFab(grids,dmap,BL_SPACEDIM,0);
            delta_rhs->setVal(0);
        }

        MultiFab** loc_viscn   = 0;
        MultiFab** loc_viscnp1 = 0;
	FluxBoxes fb_viscn, fb_viscnp1;

        if (variable_vel_visc)
        {
            Real viscTime = state[State_Type].prevTime();
	    loc_viscn = fb_viscn.define(this);
            getViscosity(loc_viscn, viscTime);

            viscTime = state[State_Type].curTime();
	    loc_viscnp1 = fb_viscnp1.define(this);
            getViscosity(loc_viscnp1, viscTime);
        }

        diffuse_velocity_setup(dt, delta_rhs, loc_viscn, loc_viscnp1);

        diffusion->diffuse_velocity(dt,be_cn_theta,get_rho_half_time(),rho_flag,
                                    delta_rhs,loc_viscn,loc_viscnp1);

        delete delta_rhs;
    }

    if (verbose)
    {
        Real run_time    = ParallelDescriptor::second() - strt_time;
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	amrex::Print() << "NavierStokes:velocity_diffusion_update(): lev: " << level
		       << ", time: " << run_time << '\n';
    }
}

void
NavierStokes::diffuse_velocity_setup (Real       dt,
                                      MultiFab*& delta_rhs,
                                      MultiFab**& viscn,
                                      MultiFab**& viscnp1)
{
    if (S_in_vel_diffusion && have_divu)
    {
        //
        // Include div mu S*I terms in rhs
        //  (i.e. make nonzero delta_rhs to add into RHS):
        //
        // The scalar and tensor solvers incorporate the relevant pieces of
        //  of Div(tau), provided the flow is divergence-free.  However, if
        //  Div(U) =/= 0, there is an additional piece not accounted for,
        //  which is of the form A.Div(U).  For constant viscosity, Div(tau)_i
        //  = Lapacian(U_i) + mu/3 d[Div(U)]/dx_i.  For mu not constant,
        //  Div(tau)_i = d[ mu(du_i/dx_j + du_j/dx_i) ]/dx_i - 2mu/3 d[Div(U)]/dx_i
        //
        // As a convenience, we treat this additional term as a "source" in
        // the diffusive solve, computing Div(U) in the "normal" way we
        // always do--via a call to calc_divu.  This routine computes delta_rhs
        // if necessary, and stores it as an auxilliary rhs to the viscous solves.
        // This is a little strange, but probably not bad.
        //
        const Real time = state[State_Type].prevTime();

        MultiFab divmusi(grids,dmap,BL_SPACEDIM,0);

        if (!variable_vel_visc)
        {
            diffusion->compute_divmusi(time,visc_coef[Xvel],divmusi);
            divmusi.mult((1./3.)*(1.0-be_cn_theta),0,BL_SPACEDIM,0);
            (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);

            diffusion->compute_divmusi(time+dt,visc_coef[Xvel],divmusi);
            divmusi.mult((1./3.)*be_cn_theta,0,BL_SPACEDIM,0);
            (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);
        }
        else
        {
            diffusion->compute_divmusi(time,viscn,divmusi);
            divmusi.mult((-2./3.)*(1.0-be_cn_theta),0,BL_SPACEDIM,0);
            (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);

            diffusion->compute_divmusi(time+dt,viscnp1,divmusi);
            divmusi.mult((-2./3.)*be_cn_theta,0,BL_SPACEDIM,0);
            (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);
        }
    }
}

Real 
NavierStokes::MaxVal (const std::string& name,
                      Real           time)
{
    Real        mxval = 0.0;
    auto        mf = derive(name,time,0);
    BoxArray    baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    std::vector< std::pair<int,Box> > isects;

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        const int  i   = mfi.index();
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            baf.intersections(grids[i],isects);

            for (int ii = 0, N = isects.size(); ii < N; ii++)
                fab.setVal(0,isects[ii].second,0,fab.nComp());
        }
        Real        s;
        const Real* dat = fab.dataPtr();
        const int*  dlo = fab.loVect();
        const int*  dhi = fab.hiVect();
	const Box&  bx  = grids[i];
        const int*  lo  = bx.loVect();
        const int*  hi  = bx.hiVect();

        FORT_MAXVAL(dat,ARLIM(dlo),ARLIM(dhi),ARLIM(lo),ARLIM(hi),&s);

        mxval = std::max(mxval, s);
    }

    ParallelDescriptor::ReduceRealMax(mxval);

    return mxval;
}

void
NavierStokes::sum_integrated_quantities ()
{
    const int finest_level = parent->finestLevel();

    Real time = state[State_Type].curTime();
    // Real mass = 0.0;
    // Real trac = 0.0;
    Real energy = 0.0;
    Real mgvort = 0.0;
    Real udotlapu = 0.0;
    Real forcing = 0.0;

    for (int lev = 0; lev <= finest_level; lev++)
    {
        NavierStokes& ns_level = getLevel(lev);
        // mass += ns_level.volWgtSum("density",time);
        // trac += ns_level.volWgtSum("tracer",time);
        energy += ns_level.volWgtSum("energy",time);
        mgvort = std::max(mgvort,ns_level.MaxVal("mag_vort",time));
#if (BL_SPACEDIM==3)
        udotlapu += ns_level.volWgtSum("udotlapu",time);
#if defined(GENGETFORCE) || defined(MOREGENGETFORCE)
	forcing += ns_level.volWgtSum("forcing",time);
#endif
#endif
    }

    amrex::Print() << '\n';
    // amrex::Print().SetPrecision(12) << "TIME= " << time << " MASS= " << mass << '\n';
    // amrex::Print().SetPrecision(12) << "TIME= " << time << " TRAC= " << trac << '\n';
    amrex::Print().SetPrecision(12) << "TIME= " << time << " KENG= " << energy << '\n';
    amrex::Print().SetPrecision(12) << "TIME= " << time << " MAGVORT= " << mgvort << '\n';
    amrex::Print().SetPrecision(12) << "TIME= " << time << " ENERGY= " << energy << '\n';
#if (BL_SPACEDIM==3)
    amrex::Print().SetPrecision(12) << "TIME= " << time << " UDOTLAPU= " << udotlapu << '\n';
#if defined(GENGETFORCE) || defined(MOREGENGETFORCE)
    amrex::Print().SetPrecision(12) << "TIME= " << time << " FORCING= " << forcing << '\n';
#endif
#endif
}

void
NavierStokes::setPlotVariables()
{
    AmrLevel::setPlotVariables();
}

void
NavierStokes::writePlotFile (const std::string& dir,
                             std::ostream&  os,
                             VisMF::How     how)
{
    if ( ! Amr::Plot_Files_Output() ) return;

    BL_PROFILE("NavierStokes::writePlotFile()");

    int i, n;
    //
    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    //
    std::vector<std::pair<int,int> > plot_var_map;
    for (int typ = 0; typ < desc_lst.size(); typ++)
        for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
            if (parent->isStatePlotVar(desc_lst[typ].name(comp)) &&
                desc_lst[typ].getType() == IndexType::TheCellType())
                plot_var_map.push_back(std::pair<int,int>(typ,comp));

    int num_derive = 0;
    std::list<std::string> derive_names;
    const std::list<DeriveRec>& dlist = derive_lst.dlist();

    for (std::list<DeriveRec>::const_iterator it = dlist.begin(), end = dlist.end();
         it != end;
         ++it)
    {
        if (parent->isDerivePlotVar(it->name()))
	{
            derive_names.push_back(it->name());
            num_derive += it->numDerive();
	}
    }

    int n_data_items = plot_var_map.size() + num_derive;
    Real cur_time = state[State_Type].curTime();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0)
            amrex::Error("Must specify at least one valid data item to plot");

        os << n_data_items << '\n';

	//
	// Names of variables -- first state, then derived
	//
	for (i =0; i < plot_var_map.size(); i++)
        {
	    int typ  = plot_var_map[i].first;
	    int comp = plot_var_map[i].second;
	    os << desc_lst[typ].name(comp) << '\n';
        }

	for (std::list<std::string>::const_iterator it = derive_names.begin(), end = derive_names.end();
             it != end;
             ++it)
        {
	    const DeriveRec* rec = derive_lst.get(*it);
	    for (i = 0; i < rec->numDerive(); i++)
                os << rec->variableName(i) << '\n';
        }
        os << BL_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbLo(i) << ' ';
        os << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbHi(i) << ' ';
        os << '\n';
        for (i = 0; i < f_lev; i++)
            os << parent->refRatio(i)[0] << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->Geom(i).Domain() << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->levelSteps(i) << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++)
                os << parent->Geom(i).CellSize()[k] << ' ';
            os << '\n';
        }
        os << (int) Geometry::Coord() << '\n';
        os << "0\n"; // Write bndry data.


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

	time_t now = time(0);

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
	jobInfoFile << "BoxLib dir:    " << buildInfoGetAMReXDir() << '\n';

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
	  jobInfoFile << "BoxLib git hash: " << githash2 << "\n";
	}

	jobInfoFile << "\n\n";


	// runtime parameters
	jobInfoFile << PrettyLine;
	jobInfoFile << " Inputs File Parameters\n";
	jobInfoFile << PrettyLine;
	
	ParmParse::dumpTable(jobInfoFile, true);

	jobInfoFile.close();
	
    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";

    std::string Level = amrex::Concatenate("Level_", level, 1);
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
        FullPath += '/';
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!amrex::UtilCreateDirectory(FullPath, 0755))
            amrex::CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
            for (n = 0; n < BL_SPACEDIM; n++)
                os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (n_data_items > 0)
        {
            std::string PathNameInHeader = Level;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }
    }
    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    // NOTE: we are assuming that each state variable has one component,
    // but a derived variable is allowed to have multiple components.
    int       cnt   = 0;
    int       ncomp = 1;
    const int nGrow = 0;
    MultiFab  plotMF(grids,dmap,n_data_items,nGrow);
    MultiFab* this_dat = 0;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (i = 0; i < plot_var_map.size(); i++)
    {
	int typ  = plot_var_map[i].first;
	int comp = plot_var_map[i].second;
	this_dat = &state[typ].newData();
	MultiFab::Copy(plotMF,*this_dat,comp,cnt,ncomp,nGrow);
	cnt+= ncomp;
    }
    //
    // Cull data from derived variables.
    // 
    Real plot_time;

    if (derive_names.size() > 0)
    {
	for (std::list<std::string>::const_iterator it = derive_names.begin(), end = derive_names.end();
             it != end;
             ++it) 
	{
            if (*it == "avg_pressure" ||
                *it == "gradpx"       ||
                *it == "gradpy"       ||
                *it == "gradpz") 
            {
                if (state[Press_Type].descriptor()->timeType() == 
                    StateDescriptor::Interval) 
                {
                    plot_time = cur_time;
                }
                else
                {
                    int f_lev = parent->finestLevel();
                    plot_time = getLevel(f_lev).state[Press_Type].curTime();
                }
            }
            else
            {
                plot_time = cur_time;
            } 
	    const DeriveRec* rec = derive_lst.get(*it);
	    ncomp = rec->numDerive();
	    auto derive_dat = derive(*it,plot_time,nGrow);
	    MultiFab::Copy(plotMF,*derive_dat,0,cnt,ncomp,nGrow);
	    cnt += ncomp;
	}
    }
    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(plotMF,TheFullPath,how,true);
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
    // Coarse levels are fine level averages, pressure is zero.
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
#if (BL_SPACEDIM==3)
    //
    // Derive turbulent statistics
    //
    if (turb_interval > 0)
        sum_turbulent_quantities();
#ifdef SUMJET
    //
    // Derive turbulent statistics for the round jet
    //
    if (jet_interval > 0)
        sum_jet_quantities();
#endif
#endif
}

//
// Initialize the pressure by iterating the initial timestep
//

void
NavierStokes::post_init_press (Real&        dt_init,
                               Vector<int>&  nc_save,
                               Vector<Real>& dt_save)
{
    const Real strt_time       = state[State_Type].curTime();
    const int  finest_level    = parent->finestLevel();
    NavierStokes::initial_iter = true;
    //
    // Iterate over the advance function.
    //
    for (int iter = 0; iter < init_iter; iter++)
    {
        for (int k = 0; k <= finest_level; k++ )
        {
            getLevel(k).advance(strt_time,dt_init,1,1);
        }
        //
        // This constructs a guess at P, also sets p_old == p_new.
        //
        Vector<MultiFab*> sig(finest_level+1, nullptr);

        for (int k = 0; k <= finest_level; k++)
        {
            sig[k] = &(getLevel(k).get_rho_half_time());
        }
        if (projector) {
            projector->initialSyncProject(0,sig,parent->dtLevel(0),
                                          strt_time,have_divu);
        }

        for (int k = finest_level-1; k >= 0; k--)
        {
            getLevel(k).avgDown();
        }
        for (int k = 0; k <= finest_level; k++)
        {
            //
            // Reset state variables to initial time, but 
            // do not reset pressure variable, only pressure time.
            //
            getLevel(k).resetState(strt_time, dt_init, dt_init);
        }

        make_rho_curr_time();

        NavierStokes::initial_iter = false;
    }

    if (init_iter <= 0)
        NavierStokes::initial_iter = false; // Just being compulsive -- rbp.

    NavierStokes::initial_step = false;
    //
    // Re-instate timestep.
    //
    for (int k = 0; k <= finest_level; k++)
    {
        getLevel(k).setTimeLevel(strt_time,dt_save[k],dt_save[k]);
        if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point)
        {
            getLevel(k).state[Press_Type].setNewTimeLevel(.5*dt_init);
            getLevel(k).get_old_data(Dpdt_Type).setVal(0);
        }
    }

    parent->setDtLevel(dt_save);
    parent->setNCycle(nc_save);
}

//
// The Mac Sync correction function
//

void
NavierStokes::mac_sync ()
{
    BL_PROFILE_REGION_START("R::NavierStokes::mac_sync()");
    BL_PROFILE("NavierStokes::mac_sync()");

    const int  numscal        = NUM_STATE - BL_SPACEDIM;
    const Real prev_time      = state[State_Type].prevTime();
    const Real prev_pres_time = state[Press_Type].prevTime();
    const Real dt             = parent->dtLevel(level);
    MultiFab*  DeltaSsync     = 0;// hold (Delta rho)*q for conserved quantities
    MultiFab&  Rh             = get_rho_half_time();

    sync_setup(DeltaSsync);
    //
    // Compute the u_mac for the correction.
    //
    mac_projector->mac_sync_solve(level,dt,Rh,fine_ratio);
    //
    // Update coarse grid state by adding correction from mac_sync solve
    // the correction is the advective tendency of the new velocities.
    //
    if (do_reflux)
    {
        MultiFab& S_new = get_new_data(State_Type);
        mac_projector->mac_sync_compute(level,u_mac,Vsync,Ssync,Rh,
                                        level > 0 ? &getAdvFluxReg(level) : 0,
                                        advectionType, prev_time,
                                        prev_pres_time,dt,
                                        NUM_STATE,be_cn_theta, 
                                        modify_reflux_normal_vel,
                                        do_mom_diff);
        //
        // The following used to be done in mac_sync_compute.  Ssync is
        // the source for a rate of change to S over the time step, so
        // Ssync*dt is the source to the actual sync amount.
        //
        Ssync.mult(dt,Ssync.nGrow());
        //
        // For all conservative variables Q (other than density)
        // express Q as rho*q and increment sync by -(sync_for_rho)*q
        // (See Pember, et. al., LBNL-41339, Jan. 1989)
        //
        FArrayBox delta_ssync;

        int iconserved = -1;
        for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
        {
            if (istate != Density && advectionType[istate] == Conservative)
            {
                iconserved++;

                for (MFIter Smfi(S_new); Smfi.isValid(); ++Smfi)
                {
                    const int  i   = Smfi.index();
                    const Box& grd = grids[i];

                    delta_ssync.resize(grd,1);
                    delta_ssync.copy(S_new[Smfi], grd, istate, grd, 0, 1);
                    delta_ssync.divide(S_new[Smfi], grd, Density, 0, 1);
                    delta_ssync.mult(Ssync[Smfi],grd,Density-BL_SPACEDIM,0,1);
                    (*DeltaSsync)[Smfi].copy(delta_ssync,grd,0,grd,iconserved,1);
                    Ssync[Smfi].minus(delta_ssync,grd,0,istate-BL_SPACEDIM,1);
                }
            }
        }

        delta_ssync.clear();

        if (do_mom_diff == 1)
        {
            for (MFIter Vsyncmfi(Vsync); Vsyncmfi.isValid(); ++Vsyncmfi)
            {
                const int        i      = Vsyncmfi.index();
                FArrayBox&       vfab   = Vsync[Vsyncmfi];
                const FArrayBox& rhofab = rho_ctime[Vsyncmfi];
                const Box&       bx     = rho_ctime.box(i);

                D_TERM(vfab.divide(rhofab,bx,0,Xvel,1);,
                       vfab.divide(rhofab,bx,0,Yvel,1);,
                       vfab.divide(rhofab,bx,0,Zvel,1););
            }
        }
        //
        // Compute viscous sync.
        //
        if (is_diffusive[Xvel])
        {
            int rho_flag = (do_mom_diff == 0) ? 1 : 3;

            MultiFab** loc_viscn = 0;
	    FluxBoxes fb_viscn;

            if (variable_vel_visc)
            {
                Real viscTime = state[State_Type].prevTime();
		loc_viscn = fb_viscn.define(this);
                getViscosity(loc_viscn, viscTime);
            }

            diffusion->diffuse_Vsync(Vsync,dt,be_cn_theta,Rh,rho_flag,loc_viscn,0);
        }

	FluxBoxes fb_SC;
        MultiFab** fluxSC        = 0;
        bool       any_diffusive = false;
        for (int sigma  = 0; sigma < numscal; sigma++)
            if (is_diffusive[BL_SPACEDIM+sigma])
                any_diffusive = true;

        if (any_diffusive) {
	    fluxSC = fb_SC.define(this);
	}

        for (int sigma = 0; sigma<numscal; sigma++)
        {
            const int state_ind = BL_SPACEDIM + sigma;
            const int rho_flag  = Diffusion::set_rho_flag(diffusionType[state_ind]);

            if (is_diffusive[state_ind])
            {
		FluxBoxes fb_diffn;
                MultiFab** cmp_diffn=0;

                if (variable_scal_diff)
                {
                    Real diffTime = state[State_Type].prevTime();
		    cmp_diffn = fb_diffn.define(this);
                    getDiffusivity(cmp_diffn, diffTime, BL_SPACEDIM+sigma,0,1);
                }

                diffusion->diffuse_Ssync(Ssync,sigma,dt,be_cn_theta,
                                         Rh,rho_flag,fluxSC,0,cmp_diffn,0,0,0);

                //
                // Increment the viscous flux registers
                //
                if (level > 0)
                {
                    for (int d = 0; d < BL_SPACEDIM; d++)
                    {
                        getViscFluxReg().FineAdd(*fluxSC[d],d,0,state_ind,1,dt);
                    }
                }
            }
        }

        //
        // For all conservative variables Q (other than density)
        // increment sync by (sync_for_rho)*q_presync.
        // (See Pember, et. al., LBNL-41339, Jan. 1989)
        //
        iconserved = -1;
        for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
        {
            if (istate != Density && advectionType[istate] == Conservative)
            {
                iconserved++;

                for (MFIter SsyncMfi(Ssync); SsyncMfi.isValid(); ++SsyncMfi)
                {
                    const int i = SsyncMfi.index();
                    Ssync[SsyncMfi].plus((*DeltaSsync)[SsyncMfi], grids[i],
                                     iconserved, istate-BL_SPACEDIM, 1);
                }
            }
        }
        //
        // Add the sync correction to the state.
        //
        for (int sigma  = 0; sigma < numscal; sigma++)
        {
            for (MFIter S_newmfi(S_new); S_newmfi.isValid(); ++S_newmfi)
            {
                S_new[S_newmfi].plus(Ssync[S_newmfi],S_newmfi.validbox(),
                                     sigma,BL_SPACEDIM+sigma,1);
            }
        }
        //
        // Update rho_ctime after rho is updated with Ssync.
        //
        make_rho_curr_time();

        if (level > 0) incrRhoAvg(Ssync,Density-BL_SPACEDIM,1.0);
        //
        // Get boundary conditions.
        //
        const int N = grids.size();

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
            MultiFab sync_incr(fine_grids,fine_lev.DistributionMap(),numscal,0);
            sync_incr.setVal(0.0);

            SyncInterp(Ssync,level,sync_incr,lev,ratio,0,0,
                       numscal,1,mult,sync_bc.dataPtr());

            MultiFab& S_new = fine_lev.get_new_data(State_Type);
            for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
                S_new[mfi].plus(sync_incr[mfi],fine_grids[mfi.index()],0,Density,numscal);

            fine_lev.make_rho_curr_time();
            fine_lev.incrRhoAvg(sync_incr,Density-BL_SPACEDIM,1.0);
        }
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

    BL_ASSERT(do_reflux);
    //
    // First do refluxing step.
    //
    FluxRegister& fr_adv  = getAdvFluxReg(level+1);
    FluxRegister& fr_visc = getViscFluxReg(level+1);
    const Real    dt_crse = parent->dtLevel(level);
    const Real    scale   = 1.0/dt_crse;
    //
    // It is important, for do_mom_diff == 0, to do the viscous
    //   refluxing first, since this will be divided by rho_half
    //   before the advective refluxing is added.  In the case of
    //   do_mom_diff == 1, both components of the refluxing will
    //   be divided by rho^(n+1) in level_sync.
    //

    fr_visc.Reflux(Vsync,volume,scale,0,0,BL_SPACEDIM,geom);
    fr_visc.Reflux(Ssync,volume,scale,BL_SPACEDIM,0,NUM_STATE-BL_SPACEDIM,geom);

    const MultiFab& Rh = get_rho_half_time();

    if (do_mom_diff == 0)
    {
        for (MFIter Vsyncmfi(Vsync); Vsyncmfi.isValid(); ++Vsyncmfi)
        {
            const int        i     = Vsyncmfi.index();
            FArrayBox&       vfab  = Vsync[Vsyncmfi];
            const FArrayBox& rhfab = Rh[Vsyncmfi];
            const Box&       bx    = Rh.box(i);

            D_TERM(vfab.divide(rhfab,bx,0,Xvel,1);,
                   vfab.divide(rhfab,bx,0,Yvel,1);,
                   vfab.divide(rhfab,bx,0,Zvel,1););
        }
    }

    for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
    {
        if (advectionType[istate] == NonConservative)
        {
            const int sigma = istate -  BL_SPACEDIM;

            for (MFIter Ssyncmfi(Ssync); Ssyncmfi.isValid(); ++Ssyncmfi)
            {
                const int i = Ssyncmfi.index();

                Ssync[Ssyncmfi].divide(Rh[Ssyncmfi],Rh.box(i),0,sigma,1);
            }
        }
    }

    fr_adv.Reflux(Vsync,volume,scale,0,0,BL_SPACEDIM,geom);
    fr_adv.Reflux(Ssync,volume,scale,BL_SPACEDIM,0,NUM_STATE-BL_SPACEDIM,geom);

    const BoxArray& fine_boxes = getLevel(level+1).boxArray();
    //
    // Zero out coarse grid cells which underlie fine grid cells.
    //
    BoxArray baf = fine_boxes;

    baf.coarsen(fine_ratio);

    std::vector< std::pair<int,Box> > isects;

    for (MFIter Vsyncmfi(Vsync); Vsyncmfi.isValid(); ++Vsyncmfi)
    {
        const int i     = Vsyncmfi.index();
        FArrayBox& vfab = Vsync[Vsyncmfi];
        FArrayBox& sfab = Ssync[Vsyncmfi];

        BL_ASSERT(grids[i] == Vsyncmfi.validbox());

        baf.intersections(Vsyncmfi.validbox(),isects);

        for (int ii = 0, N = isects.size(); ii < N; ii++)
        {
            vfab.setVal(0,isects[ii].second,0,BL_SPACEDIM);
            sfab.setVal(0,isects[ii].second,0,NUM_STATE-BL_SPACEDIM);
        }
    }
}

//
// Average down a single state component.
//

void
NavierStokes::avgDown (int comp)
{
    if (level == parent->finestLevel())
        return;

    NavierStokes&   crse_lev = getLevel(level  );
    NavierStokes&   fine_lev = getLevel(level+1);
    MultiFab&       S_crse   = get_new_data(State_Type);
    MultiFab&       S_fine   = fine_lev.get_new_data(State_Type);

    amrex::average_down(S_fine, S_crse, fine_lev.geom, crse_lev.geom, 
                         comp, 1, fine_ratio);

    if (comp == Density) 
    {
        //
        // Fill rho_ctime at current and finer levels with the correct data.
        //
        for (int lev = level; lev <= parent->finestLevel(); lev++)
            getLevel(lev).make_rho_curr_time();
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

    NavierStokes&   crse_lev = getLevel(level  );
    NavierStokes&   fine_lev = getLevel(level+1);
    //
    // Average down the states at the new time.
    //
    MultiFab& S_crse = get_new_data(State_Type);
    MultiFab& S_fine = fine_lev.get_new_data(State_Type);

    amrex::average_down(S_fine, S_crse, fine_lev.geom, crse_lev.geom, 
                         0, S_crse.nComp(), fine_ratio);

    //   
    // Now average down pressure over time n-(n+1) interval.
    //
    MultiFab&       P_crse      = get_new_data(Press_Type);
    MultiFab&       P_fine_init = fine_lev.get_new_data(Press_Type);
    MultiFab&       P_fine_avg  = fine_lev.p_avg;
    MultiFab&       P_fine      = initial_step ? P_fine_init : P_fine_avg;
    const BoxArray& P_fgrids    = fine_lev.state[Press_Type].boxArray();

    BoxArray crse_P_fine_BA = P_fgrids; crse_P_fine_BA.coarsen(fine_ratio);

    MultiFab crse_P_fine(crse_P_fine_BA,fine_lev.DistributionMap(),1,0);

    for (MFIter mfi(P_fine); mfi.isValid(); ++mfi)
    {
        const int i = mfi.index();

        injectDown(crse_P_fine_BA[i],crse_P_fine[mfi],P_fine[mfi],fine_ratio);
    }
    P_crse.copy(crse_P_fine, parent->Geom(level).periodicity());

    crse_P_fine.clear();
    //
    // Next average down divu and dSdT at new time.
    //
    if (have_divu)
    {
        MultiFab& Divu_crse = get_new_data(Divu_Type);
        MultiFab& Divu_fine = fine_lev.get_new_data(Divu_Type);
        
        amrex::average_down(Divu_fine, Divu_crse, fine_lev.geom, crse_lev.geom, 
                             0, 1, fine_ratio);
    }
    if (have_dsdt)
    {
        MultiFab& Dsdt_crse = get_new_data(Dsdt_Type);
        MultiFab& Dsdt_fine = fine_lev.get_new_data(Dsdt_Type);
        
        amrex::average_down(Dsdt_fine, Dsdt_crse, fine_lev.geom, crse_lev.geom, 
                             0, 1, fine_ratio);
    }
    //
    // Fill rho_ctime at the current and finer levels with the correct data.
    //
    for (int lev = level; lev <= parent->finestLevel(); lev++)
    {
        getLevel(lev).make_rho_curr_time();
    }
}

//
// Default divU is set to zero.
//

void
NavierStokes::calc_divu (Real      time,
                         Real      dt,
                         MultiFab& divu)
{
    BL_PROFILE("NavierStokes::calc_divu()");

    if (have_divu)
    {
        divu.setVal(0);

        if (do_temp && visc_coef[Temp] > 0.0)
        {
            //
            // Compute Div(U) = Div(visc_cond_coef * Grad(T))/(c_p*rho*T)
            //
            getViscTerms(divu,Temp,1,time);

            const MultiFab&   rhotime = get_rho(time);
            FillPatchIterator temp_fpi(*this,divu,0,time,State_Type,Temp,1);
            MFIter            rho_mfi(rhotime);

            for ( ;
                  rho_mfi.isValid() && temp_fpi.isValid();
                  ++rho_mfi, ++temp_fpi)
            {
                const int  i       = rho_mfi.index();
                FArrayBox& divufab = divu[rho_mfi];

                divufab.divide(rhotime[rho_mfi],divu.box(i),0,0,1);
                divufab.divide(temp_fpi(),divu.box(i),0,0,1);
            }
            Real THERMO_cp_inv = 1.0 / 1004.6;
            divu.mult(THERMO_cp_inv);
        }
    }
}

//
// Default dSdt is set to zero.
//

void
NavierStokes::calc_dsdt (Real      time,
                         Real      dt,
                         MultiFab& dsdt)
{
    if (have_divu && have_dsdt)
    {
        dsdt.setVal(0);

        if (do_temp)
        {
            MultiFab& Divu_new = get_new_data(Divu_Type);
            MultiFab& Divu_old = get_old_data(Divu_Type);

            for (MFIter mfi(dsdt); mfi.isValid(); ++mfi)
            {
                const Box& vbx     = mfi.validbox();
                FArrayBox& dsdtfab = dsdt[mfi];
                dsdtfab.copy(Divu_new[mfi],vbx,0,vbx,0,1);
                dsdtfab.minus(Divu_old[mfi],vbx,0,0,1);
                dsdtfab.divide(dt,vbx,0,1);
            }
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
    if (src_comp<BL_SPACEDIM && (src_comp!=Xvel || ncomp<BL_SPACEDIM))
    {
      amrex::Print() << "src_comp=" << src_comp << "   ncomp=" << ncomp << '\n';
      amrex::Error("must call NavierStokes::getViscTerms with all three velocity components");
    }
#endif
    // 
    // Initialize all viscous terms to zero
    //
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
        MultiFab** viscosity = 0;

        if (variable_vel_visc)
        {
	    viscosity = fb.define(this);
            getViscosity(viscosity, time);

            diffusion->getTensorViscTerms(visc_terms,time,viscosity,0);
        }
        else
        {
            for (int icomp = Xvel; icomp < BL_SPACEDIM; icomp++)
            {
                int rho_flag = Diffusion::set_rho_flag(diffusionType[icomp]);

                diffusion->getViscTerms(visc_terms,src_comp,icomp,time,rho_flag,0,0);
            }
        }
        //
        // Add Div(u) term if desired, if this is velocity, and if Div(u) 
        // is nonzero.  If const-visc, term is mu.Div(u)/3, else 
        // it's -Div(mu.Div(u).I)*2/3
        //
        if (have_divu && S_in_vel_diffusion)
        {
            MultiFab divmusi(grids,dmap,BL_SPACEDIM,1);

            if (variable_vel_visc)
            {
                diffusion->compute_divmusi(time,viscosity,divmusi);
                divmusi.mult((-2./3.),0,BL_SPACEDIM,0);
            }
            else
            {
                diffusion->compute_divmusi(time,visc_coef[Xvel],divmusi);
                divmusi.mult((1./3.),0,BL_SPACEDIM,0);
            }

            visc_terms.plus(divmusi,Xvel,BL_SPACEDIM,0);
        }
    }
    //
    // Get Scalar Diffusive Terms
    //
    const int first_scal = (src_comp==Xvel) ? BL_SPACEDIM : src_comp;
    const int num_scal   = (src_comp==Xvel) ? ncomp-BL_SPACEDIM : ncomp;

    if (num_scal > 0)
    {
        for (int icomp = first_scal; icomp < first_scal+num_scal; icomp++)
        {
            if (is_diffusive[icomp])
            {
		diffusive = true;

                int rho_flag = Diffusion::set_rho_flag(diffusionType[icomp]);

		FluxBoxes fb;
                MultiFab** cmp_diffn = 0;

                if (variable_scal_diff)
                {
		    cmp_diffn = fb.define(this);
                    getDiffusivity(cmp_diffn, time, icomp, 0, 1);
                }

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
// Functions for calculating the variable viscosity and diffusivity.
// These default to setting the variable viscosity and diffusivity arrays
// to the values in visc_coef and diff_coef.  These functions would
// need to be replaced in any class derived from NavierStokes that
// wants variable coefficients.
//
void 
NavierStokes::calcViscosity (const Real time, 
                             const Real dt,
                             const int  iteration,
                             const int  ncycle)
{
    //
    // Select time level to work with (N or N+1)
    //
    MultiFab* visc_cc = 0;

    const TimeLevel whichTime = which_time(State_Type,time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    if (whichTime == AmrOldTime)               // time N
    {
        visc_cc = viscn_cc;
    }
    else if (whichTime == AmrNewTime)          // time N+1
    {
        visc_cc = viscnp1_cc;
    }
    //
    // Calculate viscosity
    //
    const int nGrow = visc_cc->nGrow();

    if (is_diffusive[Xvel]) 
    {
        if (visc_coef[Xvel] >= 0.0)
        {
            visc_cc->setVal(visc_coef[Xvel], 0, 1, nGrow);
        }
        else
        {
            amrex::Abort("NavierStokes::calcViscosity() : must have velocity visc_coef >= 0.0");
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
    int ncomp = NUM_STATE - BL_SPACEDIM -1;
    //
    // Select time level to work with (N or N+1)
    //
    MultiFab* diff_cc = 0;

    const TimeLevel whichTime = which_time(State_Type,time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    if (whichTime == AmrOldTime)               // time N
    {
        diff_cc = diffn_cc;
    }
    else if (whichTime == AmrNewTime)          // time N+1
    {
        diff_cc = diffnp1_cc;
    }
    //
    // Calculate diffusivity
    //
    const int nGrow = diff_cc->nGrow();

    for (int comp=src_comp; comp<src_comp+ncomp; comp++)
    {
        int diff_comp = comp - Density - 1;

        if (is_diffusive[comp])
        {
            if (visc_coef[comp] >= 0.0)
            {
                diff_cc->setVal(visc_coef[comp], diff_comp, 1, nGrow);
            }
            else
            {
                amrex::Abort("NavierStokes::calcDiffusivity() : must have scalar diff_coefs >= 0.0");
            }
        }
    }
}

void 
NavierStokes::getViscosity (MultiFab* viscosity[BL_SPACEDIM],
                            const Real time)
{
    //
    // Select time level to work with (N or N+1)
    //
    MultiFab* visc_cc = 0;

    const TimeLevel whichTime = which_time(State_Type,time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    if (whichTime == AmrOldTime)               // time N
    {
        visc_cc = viscn_cc;
    }
    else if (whichTime == AmrNewTime)          // time N+1
    {
        visc_cc = viscnp1_cc;
    }
    //
    // Fill edge-centered viscosity
    //
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        for (MFIter ecMfi(*viscosity[dir]); ecMfi.isValid(); ++ecMfi)
        {
            center_to_edge_plain((*visc_cc)[ecMfi],(*viscosity[dir])[ecMfi],0,0,1);
        }
    }
}

void 
NavierStokes::getDiffusivity (MultiFab* diffusivity[BL_SPACEDIM],
                              const Real time,
                              const int state_comp,
                              const int dst_comp,
                              const int ncomp)
{
    BL_ASSERT(state_comp > Density);
    //
    // Pick correct diffusivity component
    //
    int diff_comp = state_comp - Density - 1;
    //
    // Select time level to work with (N or N+1)
    //
    MultiFab* diff_cc = 0;

    const TimeLevel whichTime = which_time(State_Type,time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    if (whichTime == AmrOldTime)               // time N
    {
        diff_cc = diffn_cc;
    }
    else if (whichTime == AmrNewTime)          // time N+1
    {
        diff_cc = diffnp1_cc;
    }
    //
    // Fill edge-centered diffusivities
    //
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        for (MFIter ecMfi(*diffusivity[dir]); ecMfi.isValid(); ++ecMfi)
        {
            center_to_edge_plain((*diff_cc)[ecMfi],(*diffusivity[dir])[ecMfi],
                                 diff_comp,dst_comp,ncomp);
        }
    }
}

void
NavierStokes::center_to_edge_plain (const FArrayBox& ccfab,
                                    FArrayBox&       ecfab,
                                    int              sComp,
                                    int              dComp,
                                    int              nComp)
{
    //
    // This routine fills an edge-centered FAB from a cell-centered FAB.
    // It assumes that the data in all cells of the cell-centered FAB is
    // valid and totally ignores any concept of boundary conditions.  
    // It is assummed that the cell-centered FAB fully contains the 
    // edge-centered FAB.  If anything special needs to be done at boundaries, 
    // a varient of this routine needs to be written.  See 
    // HeatTransfer::center_to_edge_fancy().
    //
    const Box&      ccbox = ccfab.box();
    const Box&      ecbox = ecfab.box();
    const IndexType ixt   = ecbox.ixType();
    //
    // Get direction for interpolation to edges
    //
    int dir = -1;
    for (int d = 0; d < BL_SPACEDIM; d++)
        if (ixt.test(d))
            dir = d;
    //
    // Miscellanious checks
    //
    BL_ASSERT(!(ixt.cellCentered()) && !(ixt.nodeCentered()));
    BL_ASSERT(amrex::grow(ccbox,-amrex::BASISV(dir)).contains(amrex::enclosedCells(ecbox)));
    BL_ASSERT(sComp+nComp <= ccfab.nComp() && dComp+nComp <= ecfab.nComp());
    //
    // Shift cell-centered data to edges
    //
    Box fillBox = ccbox; 
    for (int d = 0; d < BL_SPACEDIM; d++)
        if (d != dir)
            fillBox.setRange(d, ecbox.smallEnd(d), ecbox.length(d));
    
    const int isharm = def_harm_avg_cen2edge;
    FORT_CEN2EDG(fillBox.loVect(), fillBox.hiVect(),
                 ARLIM(ccfab.loVect()), ARLIM(ccfab.hiVect()),
                 ccfab.dataPtr(sComp),
                 ARLIM(ecfab.loVect()), ARLIM(ecfab.hiVect()),
                 ecfab.dataPtr(dComp),
                 &nComp, &dir, &isharm);
}

