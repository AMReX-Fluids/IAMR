

#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <NavierStokesBase.H>
#include <NS_BC.H>
#include <AMReX_BLProfiler.H>
#include <Projection.H>
#include <PROJECTION_F.H>
#include <NAVIERSTOKES_F.H>
#include <ProjOutFlowBC.H>

#include <AMReX_MGT_Solver.H>
#include <AMReX_stencil_types.H>
#include <mg_cpp_f.h>

#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>

using namespace amrex;

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

#define DEF_BOX_LIMITS(box,boxlo,boxhi)   \
const int* boxlo = (box).loVect();           \
const int* boxhi = (box).hiVect();

const Real Projection::BogusValue = 1.e200;

int  Projection::P_code              = 0;
int  Projection::proj_2              = 1;
int  Projection::verbose             = 0;
Real Projection::proj_tol            = 1.0e-12;
Real Projection::sync_tol            = 1.e-8;
Real Projection::proj_abs_tol        = 1.e-16;
int  Projection::add_vort_proj       = 0;
int  Projection::do_outflow_bcs      = 1;
int  Projection::rho_wgt_vel_proj    = 0;
int  Projection::make_sync_solvable  = 0;
Real Projection::divu_minus_s_factor = 0.0;

namespace
{
    bool initialized = false;

#if MG_USE_HYPRE
    bool use_hypre_solve = false;
#endif

    bool benchmarking = false;

    int hg_stencil = ND_DENSE_STENCIL;

    bool use_mlmg_solver = false;
    bool rz_correction = true;
    bool agglomeration = true;
    bool consolidation = true;
    int max_fmg_iter = 0;
    bool use_gauss_seidel = true;
    bool use_harmonic_average = false;
    bool test_mlmg_solver = false;
}


void
Projection::Initialize ()
{
    if (initialized) return;

    ParmParse pp("proj");

    pp.query("v",                   verbose);
    pp.query("Pcode",               P_code);
    pp.query("proj_2",              proj_2);
    pp.query("proj_tol",            proj_tol);
    pp.query("sync_tol",            sync_tol);
    pp.query("proj_abs_tol",        proj_abs_tol);
    pp.query("benchmarking",        benchmarking);
    pp.query("add_vort_proj",       add_vort_proj);
    pp.query("do_outflow_bcs",      do_outflow_bcs);
    pp.query("rho_wgt_vel_proj",    rho_wgt_vel_proj);
    pp.query("divu_minus_s_factor", divu_minus_s_factor);
    pp.query("make_sync_solvable",  make_sync_solvable);

    pp.query("use_mlmg_solver",     use_mlmg_solver);
    pp.query("rz_correction",       rz_correction);
    pp.query("agglomeration",       agglomeration);
    pp.query("consolidation",       consolidation);
    pp.query("max_fmg_iter",        max_fmg_iter);
    pp.query("use_gauss_seidel",    use_gauss_seidel);
    pp.query("use_harmonic_average", use_harmonic_average);
    pp.query("test_mlmg_solver",    test_mlmg_solver);

    if (!proj_2) 
	amrex::Error("With new gravity and outflow stuff, must use proj_2");

    std::string stencil;

    if ( pp.query("stencil", stencil) )
    {
        if ( stencil == "cross" )
        {
            hg_stencil = ND_CROSS_STENCIL;
        }
        else if ( stencil == "full" || stencil == "dense")
        {
            hg_stencil = ND_DENSE_STENCIL;
        }
        else
        {
            amrex::Error("Must set proj.stencil to be cross, full or dense");
        }
    }

    amrex::ExecOnFinalize(Projection::Finalize);

    initialized = true;
}

void
Projection::Finalize ()
{
    initialized = false;
}

Projection::Projection (Amr*   _parent,
                        BCRec* _phys_bc, 
                        int    _do_sync_proj,
                        int    /*_finest_level*/, 
                        int    _radius_grow )
   :
    parent(_parent),
    LevelData(_parent->finestLevel()+1),
    radius_grow(_radius_grow), 
    radius(_parent->finestLevel()+1),
    anel_coeff(_parent->finestLevel()+1),
    phys_bc(_phys_bc), 
    do_sync_proj(_do_sync_proj)
{

    BL_ASSERT ( parent->finestLevel()+1 <= maxlev );

    Initialize();

    if (verbose) amrex::Print() << "Creating projector\n";

    for (int lev = 0; lev <= parent->finestLevel(); lev++)
       anel_coeff[lev] = 0;
}

Projection::~Projection ()
{
  if (verbose) amrex::Print() << "Deleting projector\n";
}

//
// Install a level of the projection.
//

void
Projection::install_level (int                   level,
                           AmrLevel*             level_data,
                           Vector< Vector<Real> >* _radius)
{
    if (verbose) amrex::Print() << "Installing projector level " << level << '\n';

    int finest_level = parent->finestLevel();

    if (level > LevelData.size() - 1) 
    {
        LevelData.resize(finest_level+1);
        radius.resize(finest_level+1);
    }

    if (level > anel_coeff.size()-1) {
       anel_coeff.resize(level+1);
       anel_coeff[level] = 0;
    }

    LevelData[level] = level_data;
    radius[level] = _radius;
}

void
Projection::install_anelastic_coefficient (int                   level,
                                           Real                **_anel_coeff)
{
  if (verbose) {
    amrex::Print() << "Installing anel_coeff into projector level " << level << '\n';
  }
  if (level > anel_coeff.size()-1) 
    anel_coeff.resize(level+1);
  anel_coeff[level] =  _anel_coeff;
}

//
//  Perform a level projection in the advance function
//  Explanation of arguments to the level projector:
//
//  rho_half  contains rho^{n+1/2}
//  U_old  contains the u^n velocities
//  U_new  starts as u^*, is converted to (u^* - u^n)/dt,
//         becomes (u^{n+1} - u^n)/dt in the solver,
//         and is converted back to u^n+1 at the end
//  P_old  contains p^{n-1/2}
//  P_new  gets cleared, initialized to an intial guess for p^{n+1/2}
//         using coarse grid data if available,
//         becomes pressure update phi in the solver,
//         and then converted into final prssure p^{n+1/2}
//

void
Projection::level_project (int             level,
                           Real            time,
                           Real            dt,
                           Real            cur_pres_time,
                           Real            prev_pres_time,
                           const Geometry& geom, 
                           MultiFab&       U_old,
                           MultiFab&       U_new,
                           MultiFab&       P_old,
                           MultiFab&       P_new,
                           MultiFab&       rho_half, 
                           SyncRegister*   crse_sync_reg, 
                           SyncRegister*   fine_sync_reg,  
                           int             crse_dt_ratio,
                           int             iteration,
                           int             have_divu)
{
    BL_PROFILE("Projection::level_project()");

    BL_ASSERT(rho_half.nGrow() >= 1);
    BL_ASSERT(U_new.nGrow() >= 1);

    if (verbose) {
      amrex::Print() << "... level projector at level " << level << '\n';
    }

    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real strt_time = ParallelDescriptor::second();

    //
    // old time velocity has bndry values already
    // must gen valid bndry data for new time velocity.
    // must fill bndry cells in pressure with computable values
    // even though they are not used in calculation.
    //
    U_old.setBndry(BogusValue,Xvel,BL_SPACEDIM);
    U_new.setBndry(BogusValue,Xvel,BL_SPACEDIM);
    P_old.setBndry(BogusValue);
    P_new.setBndry(BogusValue);

    MultiFab& S_old = LevelData[level]->get_old_data(State_Type);
    MultiFab& S_new = LevelData[level]->get_new_data(State_Type);
    
    Real prev_time = LevelData[level]->get_state_data(State_Type).prevTime();
    Real curr_time = LevelData[level]->get_state_data(State_Type).curTime();
    
    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        LevelData[level]->setPhysBoundaryValues(S_old[mfi],State_Type,prev_time,
                                                Xvel,Xvel,BL_SPACEDIM);
        LevelData[level]->setPhysBoundaryValues(S_new[mfi],State_Type,curr_time,
                                                Xvel,Xvel,BL_SPACEDIM);
    }

    const BoxArray& grids   = LevelData[level]->boxArray();
    const DistributionMapping& dmap = LevelData[level]->DistributionMap();
    const BoxArray& P_grids = P_old.boxArray();
    const DistributionMapping& P_dmap = P_old.DistributionMap();

    NavierStokesBase* ns = dynamic_cast<NavierStokesBase*>(&parent->getLevel(level));
    BL_ASSERT(!(ns==0));

    //
    //  NOTE: IT IS IMPORTANT TO DO THE BOUNDARY CONDITIONS BEFORE
    //    MAKING UNEW HOLD U_t OR U/dt, BECAUSE UNEW IS USED IN
    //    CONSTRUCTING THE OUTFLOW BC'S.
    //
    // Set boundary values for P_new, to increment, if applicable
    // // Note: we don't need to worry here about using FillCoarsePatch because
    //       it will automatically use the "new dpdt" to interpolate,
    //       since once we get here at level > 0, we've already defined
    //       a new pressure at level-1.
    if (level != 0)
    {
	LevelData[level]->FillCoarsePatch(P_new,0,cur_pres_time,Press_Type,0,1);
        if (!proj_2) 
            P_new.minus(P_old,0,1,0); // Care about nodes on box boundary
    }

    const int nGrow = (level == 0  ?  0  :  -1);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(P_new,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.growntilebox(nGrow);
        P_new[mfi].setVal(0.0,bx,0,1);
    }

    //
    // Compute Ustar/dt + Gp                  for proj_2,
    //         (Ustar-Un)/dt for not proj_2 (ie the original).
    //
    // Compute DU/dt for proj_2,
    //         (DU-DU_old)/dt for not proj_2 (ie the original).
    //
    std::unique_ptr<MultiFab> divusource, divuold;

    if (have_divu)
    {
        divusource.reset(ns->getDivCond(1,time+dt));
        if (!proj_2) {
            divuold.reset(ns->getDivCond(1,time));
        }
    }

    const Real dt_inv = 1./dt;
    if (proj_2)
    {
        U_new.mult(dt_inv,0,BL_SPACEDIM,1);
        if (have_divu)
            divusource->mult(dt_inv,0,1,divusource->nGrow());
    }
    else
    {
        for (MFIter U_newmfi(U_new); U_newmfi.isValid(); ++U_newmfi) 
        {
            const int i = U_newmfi.index();

            ConvertUnew(U_new[U_newmfi],U_old[U_newmfi],dt,U_new.box(i));
        } 

        if (have_divu)
        {
            divusource->minus(*divuold,0,1,divusource->nGrow());
            divusource->mult(dt_inv,0,1,divusource->nGrow());

            if (divu_minus_s_factor>0.0 && divu_minus_s_factor<=1.0)
            {
                amrex::Error("Check this code....not recently tested");
                //
                // Compute relaxation terms to account for approximate projection
                // add divu_old*divu...factor/dt to divusource.
                //
                const Real uoldfactor = divu_minus_s_factor*dt/parent->dtLevel(0);
                UpdateArg1(*divusource, uoldfactor/dt, *divuold, 1, grids, 1);
                //
                // add U_old*divu...factor/dt to U_new
                //
                UpdateArg1(U_new, uoldfactor/dt, U_old, BL_SPACEDIM, grids, 1);
            }
        }
    }

    if (proj_2)
    {
        MultiFab Gp(grids,dmap,BL_SPACEDIM,1);
        ns->getGradP(Gp, prev_pres_time);

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(rho_half,true); mfi.isValid(); ++mfi) 
	{
	    const Box& bx = mfi.growntilebox(1);
	    FArrayBox& Gpfab = Gp[mfi];
	    const FArrayBox& rhofab = rho_half[mfi];
      
	    for (int i = 0; i < BL_SPACEDIM; i++) {
		Gpfab.divide(rhofab,bx,0,i,1);
	    }
      
	    U_new[mfi].plus(Gpfab,bx,0,0,BL_SPACEDIM);
	}
    }

    //
    // Outflow uses appropriately constructed "U_new" and "divusource"
    //   so make sure this call comes after those are set,
    //   but before fields are scaled by r or rho is set to 1/rho.
    //
    Real gravity = ns->getGravity();
    if (OutFlowBC::HasOutFlowBC(phys_bc) && (have_divu || std::fabs(gravity) > 0.0) 
                                         && do_outflow_bcs) 
    {
        Vector<MultiFab*> phi(maxlev, nullptr);
        phi[level] = &(LevelData[level]->get_new_data(Press_Type));

        Vector<MultiFab*> Vel_ML(maxlev, nullptr);
        Vel_ML[level] = &U_new;

        Vector<MultiFab*> Divu_ML(maxlev, nullptr);
        Divu_ML[level] = divusource.get();

        Vector<MultiFab*> Rho_ML(maxlev, nullptr);
        Rho_ML[level] = &rho_half;

        set_outflow_bcs(LEVEL_PROJ,phi,Vel_ML,Divu_ML,Rho_ML,level,level,have_divu);
    }

    //
    // Scale the projection variables.
    //
    rho_half.setBndry(BogusValue);
    scaleVar(LEVEL_PROJ,&rho_half, 1, &U_new, level);
    //
    // Enforce periodicity of U_new and rho_half (i.e. coefficient of G phi)
    // *after* everything has been done to them.
    //
    if (geom.isAnyPeriodic()) {
	U_new.FillBoundary(0, BL_SPACEDIM, geom.periodicity());
	rho_half.FillBoundary(0, 1, geom.periodicity());
    }
    //
    // Add the contribution from the un-projected V to syncregisters.
    //
    int is_rz = (Geometry::IsRZ() ? 1 : 0);

    Vector<MultiFab*> vel(maxlev, nullptr);
    Vector<MultiFab*> phi(maxlev, nullptr);
    Vector<MultiFab*> sig(maxlev, nullptr);

    vel[level] = &U_new;
    phi[level] = &P_new;

    BL_ASSERT( 1 == rho_half.nGrow());
    sig[level] = &rho_half;

    //
    // Project
    //
    std::unique_ptr<MultiFab> sync_resid_crse, sync_resid_fine;

    if (level < parent->finestLevel()) {
        sync_resid_crse.reset(new MultiFab(P_grids,P_dmap,1,1));
    }

    if (level > 0 && ((proj_2 && iteration == crse_dt_ratio) || !proj_2))
    {
        const int ngrow = parent->MaxRefRatio(level-1) - 1;
        sync_resid_fine.reset(new MultiFab(P_grids,P_dmap,1,ngrow));
    }

    if (!have_divu) 
    {
        Vector<MultiFab*> rhs(maxlev, nullptr);
        if (use_mlmg_solver) {
            doMLMGNodalProjection(level, 1, vel, phi, sig, rhs, {}, proj_tol, proj_abs_tol, 
                                  sync_resid_crse.get(), sync_resid_fine.get());
        } else {
            doNodalProjection(level, 1, vel, phi, sig, rhs, {}, proj_tol, proj_abs_tol, 
                              sync_resid_crse.get(), sync_resid_fine.get());
        }
    }
    else 
    {
        if (is_rz == 1) {
            radMultScal(level,*divusource);
        }
        const int nghost = 0;
        divusource->mult(-1.0,0,1,nghost);

	Vector<MultiFab*> rhs_cc(maxlev, nullptr);
	rhs_cc[level] = divusource.get();
        if (use_mlmg_solver) {
            doMLMGNodalProjection(level, 1, vel, phi, sig, rhs_cc, {}, proj_tol, proj_abs_tol,
                                  sync_resid_crse.get(), sync_resid_fine.get());
        } else {
            doNodalProjection(level, 1, vel, phi, sig, rhs_cc, {}, proj_tol, proj_abs_tol,
                              sync_resid_crse.get(), sync_resid_fine.get());
        }
    }

    //
    // Note: this must occur *after* the projection has been done
    //       (but before the modified velocity has been copied back)
    //       because the SyncRegister routines assume the projection
    //       has been set up.
    //
    if (do_sync_proj)
    {
       if (level < parent->finestLevel())
       {
          //
          // Init sync registers between level and level+1.
          //
          const Real mult = 1.0;
          crse_sync_reg->CrseInit(*sync_resid_crse,geom,mult);
       }
       if (level > 0 && ((proj_2 && iteration == crse_dt_ratio) || !proj_2))
       {
          //
          // Increment sync registers between level and level-1.
          //
	 // invrat is 1/crse_dt_ratio for both proj_2 and !proj_2, but for different reasons.
	 // For !proj_2, this is because the fine residue is added to the sync register 
	 //    for each fine step.
	 // For proj_2, this is because the level projection works on U/dt, not dU/dt, 
	 //    and dt on the fine level is crse_dt_ratio times smaller than dt one the 
	 //    coarse level.
	  const Real invrat = 1.0/(double)crse_dt_ratio;
          const Geometry& crse_geom = parent->Geom(level-1);
          fine_sync_reg->FineAdd(*sync_resid_fine,crse_geom,invrat);
       }
    }

    //
    // Reset state + pressure data.
    //
    // Unscale level projection variables.
    //
    rescaleVar(LEVEL_PROJ,&rho_half, 1, &U_new, level);
    //
    // Put U_new back to "normal"; subtract U_old*divu...factor/dt from U_new
    //
    if (!proj_2 && divu_minus_s_factor>0.0 && divu_minus_s_factor<=1.0 && have_divu) 
    {
        const Real uoldfactor = -divu_minus_s_factor*dt/parent->dtLevel(0);
        UpdateArg1(U_new, uoldfactor/dt, U_old, BL_SPACEDIM, grids, 1);
    }
    //
    // Convert U back to a velocity, and phi into p^n+1/2.
    //
    if (proj_2) 
    {
        //
        // un = dt*un
        //
        U_new.mult(dt,0,BL_SPACEDIM,1);
    }
    else
    {
        //
        // un = uo+dt*un
        //
        UnConvertUnew(U_old, dt, U_new, grids);
    }

    if (!proj_2) 
        AddPhi(P_new, P_old);             // pn = pn + po

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	amrex::Print() << "Projection::level_project(): lev: " << level
		       << ", time: " << run_time << '\n';
    }
}

//
// SYNC_PROJECT
//

void
Projection::syncProject (int             c_lev,
                         MultiFab&       pres,
                         MultiFab&       vel,
                         MultiFab&       rho_half,
                         MultiFab&       Vsync,
                         MultiFab&       phi,
                         SyncRegister*   rhs_sync_reg,
                         SyncRegister*   crsr_sync_reg,
                         const BoxArray& sync_boxes,
                         const Geometry& geom,
                         const Real*     dx,
                         Real            dt_crse,
                         int             crse_iteration,
                         int             crse_dt_ratio)
{
    BL_PROFILE("Projection::syncProject()");

    if (verbose)
    {
      amrex::Print() << "SyncProject: level = "
		     << c_lev
		     << " correction to level "
		     << parent->finestLevel() << '\n';
    }

    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real strt_time = ParallelDescriptor::second();

    //
    // Gather data.
    //
    const BoxArray& grids   = LevelData[c_lev]->boxArray();
    const BoxArray& P_grids = pres.boxArray();
    const DistributionMapping& P_dmap = pres.DistributionMap();
    MultiFab& sig = rho_half;

    MultiFab rhnd(P_grids,P_dmap,1,0);
    rhs_sync_reg->InitRHS(rhnd,geom,*phys_bc);

    phi.setVal(0);

    sig.setBndry(BogusValue);
    //
    // Scale sync projection variables.
    //
    scaleVar(SYNC_PROJ,&sig,1,&Vsync,c_lev);
    //
    // If periodic, copy into periodic translates of Vsync.
    //
    if (geom.isAnyPeriodic()) {
	Vsync.FillBoundary(0, BL_SPACEDIM, geom.periodicity());
    }

    Vector<MultiFab*> phis(maxlev, nullptr);
    Vector<MultiFab*> vels(maxlev, nullptr);
    Vector<MultiFab*> sigs(maxlev, nullptr);
    Vector<MultiFab*> rhss(maxlev, nullptr);
    phis[c_lev] = &phi;
    vels[c_lev] = &Vsync;
    sigs[c_lev] = &sig;

    //
    //  PROJECT
    //  if use_u = 0, then solves DGphi = RHS
    //  if use_u = 1, then solves DGphi = RHS + DV
    //  both return phi and (V-Gphi) as V
    //
    MultiFab* sync_resid_crse = 0;
    std::unique_ptr<MultiFab> sync_resid_fine;

    if (c_lev > 0 && (!proj_2 || crse_iteration == crse_dt_ratio))
    {
        const int ngrow = parent->MaxRefRatio(c_lev-1) - 1;
        sync_resid_fine.reset(new MultiFab(P_grids,P_dmap,1,ngrow));
    }

    if (use_mlmg_solver) {
        doMLMGNodalProjection(c_lev, 1, vels, phis, sigs, rhss, {&rhnd}, sync_tol, proj_abs_tol,
                              sync_resid_crse, sync_resid_fine.get());
    } else {
        doNodalProjection(c_lev, 1, vels, phis, sigs, rhss, {&rhnd}, sync_tol, proj_abs_tol,
                          sync_resid_crse, sync_resid_fine.get());
    }

    //
    // If this sync project is not at level 0 then we need to account for
    // the changes made here in the level c_lev velocity in the sync registers
    // going into the level (c_lev-1) sync project.  Note that this must be
    // done before rho_half is scaled back.
    //
    if (c_lev > 0 && (!proj_2 || crse_iteration == crse_dt_ratio))
    {
        const Real invrat         = 1.0/(double)crse_dt_ratio;
        const Geometry& crsr_geom = parent->Geom(c_lev-1);
        crsr_sync_reg->CompAdd(*sync_resid_fine,geom,crsr_geom,sync_boxes,invrat);
    }
    //
    // Reset state + pressure data ...
    //
    // Unscale the sync projection variables for rz.
    //
    rescaleVar(SYNC_PROJ,&sig,1,&Vsync,c_lev);
    //
    // Add projected Vsync to new velocity at this level & add phi to pressure.
    //
    AddPhi(pres, phi);
    UpdateArg1(vel, dt_crse, Vsync, BL_SPACEDIM, grids, 1);

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	amrex::Print() << "Projection:syncProject(): c_lev: " << c_lev
		       << ", time: " << run_time << '\n';
    }
}

//
//  MULTI-LEVEL SYNC_PROJECT
//

void
Projection::MLsyncProject (int             c_lev,
                           MultiFab&       pres_crse,
                           MultiFab&       vel_crse,
                           MultiFab&       cc_rhs_crse,
                           MultiFab&       pres_fine,
                           MultiFab&       vel_fine,
                           MultiFab&       cc_rhs_fine,
                           MultiFab&       rho_crse,
                           MultiFab&       rho_fine,
                           MultiFab&       Vsync,
                           MultiFab&       V_corr,
                           MultiFab&       phi_fine,
                           SyncRegister*   rhs_sync_reg,
                           SyncRegister*   crsr_sync_reg,
                           Real            dt_crse, 
                           IntVect&        ratio,
                           int             crse_iteration,
                           int             crse_dt_ratio,
                           const Geometry& crse_geom,
                           bool		   pressure_time_is_interval,
                           bool first_crse_step_after_initial_iters,
                           Real             cur_crse_pres_time,
                           Real            prev_crse_pres_time,
                           Real             cur_fine_pres_time,
                           Real            prev_fine_pres_time)
{
    BL_PROFILE("Projection::MLsyncProject()");

    if (verbose) {
      amrex::Print() << "SyncProject: levels = " << c_lev << ", " << c_lev+1 << '\n';
    }

    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real strt_time = ParallelDescriptor::second();

    //
    // Set up memory.
    //
    Vector<std::unique_ptr<MultiFab> > phi(maxlev);

    const BoxArray& grids      = LevelData[c_lev]->boxArray();
    const BoxArray& fine_grids = LevelData[c_lev+1]->boxArray();
    const BoxArray& Pgrids_crse = pres_crse.boxArray();
    const BoxArray& Pgrids_fine = pres_fine.boxArray();
    const DistributionMapping& Pdmap_crse = pres_crse.DistributionMap();
    const DistributionMapping& Pdmap_fine = pres_fine.DistributionMap();

    phi[c_lev].reset(new MultiFab(Pgrids_crse,Pdmap_crse,1,1));
    phi[c_lev]->setVal(0);
    
    phi[c_lev+1].reset(new MultiFab(Pgrids_fine,Pdmap_fine,1,1));
    phi[c_lev+1]->setVal(0);

    //
    // Set up crse RHS
    //
    MultiFab rhnd(Pgrids_crse,Pdmap_crse,1,0);
    rhs_sync_reg->InitRHS(rhnd,crse_geom,*phys_bc);

    Box P_finedomain(amrex::surroundingNodes(crse_geom.Domain()));
    P_finedomain.refine(ratio);
    if (Pgrids_fine[0] == P_finedomain) {
        rhnd.setVal(0);
    }
    //
    // Do necessary scaling
    //
    scaleVar(SYNC_PROJ,&rho_crse, 0, &Vsync,   c_lev  );
    scaleVar(SYNC_PROJ,&rho_fine, 0, &V_corr, c_lev+1);

    if (Geometry::IsRZ()) {
       radMultScal(c_lev  ,cc_rhs_crse);
       radMultScal(c_lev+1,cc_rhs_fine);
    }

    Vector<MultiFab*> vel(maxlev, nullptr);
    Vector<MultiFab*> sig(maxlev, nullptr);
    Vector<MultiFab*> rhs(maxlev, nullptr);

    vel[c_lev  ] = &Vsync;
    vel[c_lev+1] = &V_corr;
    sig[c_lev  ] = &rho_crse;
    sig[c_lev+1] = &rho_fine;
    rhs[c_lev  ] = &cc_rhs_crse; 
    rhs[c_lev+1] = &cc_rhs_fine;

    const Geometry& fine_geom = parent->Geom(c_lev+1);

    // restrict_level(v_crse, v_fine, ratio);
    amrex::average_down(*vel[c_lev+1],*vel[c_lev],fine_geom,crse_geom,
                         0, BL_SPACEDIM, ratio);

    // restrict_level(*sig[c_lev], *sig[c_lev+1], ratio);
    amrex::average_down(*sig[c_lev+1],*sig[c_lev],fine_geom,crse_geom,
                           0, sig[c_lev]->nComp(), ratio);

    MultiFab* sync_resid_crse = 0;
    std::unique_ptr<MultiFab> sync_resid_fine;

    if (c_lev > 0 && (!proj_2 || crse_iteration == crse_dt_ratio))
      //    if (c_lev > 0)
    {
        int ngrow = parent->MaxRefRatio(c_lev-1) - 1;
        sync_resid_fine.reset(new MultiFab(Pgrids_crse,Pdmap_crse,1,ngrow));
    }

    if (use_mlmg_solver) {
        doMLMGNodalProjection(c_lev, 2, vel,
                          amrex::GetVecOfPtrs(phi),
                          sig, rhs, {&rhnd}, sync_tol, proj_abs_tol,
                          sync_resid_crse, sync_resid_fine.get());
    } else {
        doNodalProjection(c_lev, 2, vel, 
                          amrex::GetVecOfPtrs(phi),
                          sig, rhs, {&rhnd}, sync_tol, proj_abs_tol,
                          sync_resid_crse, sync_resid_fine.get());
    }

    //
    // If this sync project is not at levels 0-1 then we need to account for
    // the changes made here in the level c_lev velocity in the sync registers
    // going into the level (c_lev-1) sync project.  Note that this must be
    // done before rho_half is scaled back.
    //
    if (c_lev > 0 && (!proj_2 || crse_iteration == crse_dt_ratio))
    {
        const Real invrat         = 1.0/(double)crse_dt_ratio;
        const Geometry& crsr_geom = parent->Geom(c_lev-1);
        BoxArray sync_boxes       = pres_fine.boxArray();
        sync_boxes.coarsen(ratio);
        crsr_sync_reg->CompAdd(*sync_resid_fine,crse_geom,crsr_geom,sync_boxes,invrat);
    }
    //
    // Do necessary un-scaling.
    //
    rescaleVar(SYNC_PROJ,&rho_crse, 0, &Vsync,   c_lev  );
    rescaleVar(SYNC_PROJ,&rho_fine, 0, &V_corr, c_lev+1);

    MultiFab::Copy(phi_fine, *phi[c_lev+1], 0, 0, 1, 1);

    //
    // Add phi to pressure.
    //
    AddPhi(pres_crse, *phi[c_lev]);

    if (pressure_time_is_interval) 
    {
        //
        // Only update the most recent pressure.
        //
        AddPhi(pres_fine, *phi[c_lev+1]);
    }
    else 
    {
        MultiFab& pres_fine_old = LevelData[c_lev+1]->get_old_data(Press_Type);
 
        if (first_crse_step_after_initial_iters)
        {
            Real time_since_zero =  cur_crse_pres_time - prev_crse_pres_time;
            Real dt_to_prev_time = prev_fine_pres_time - prev_crse_pres_time;
            Real dt_to_cur_time  =  cur_fine_pres_time - prev_crse_pres_time;

            Real cur_mult_factor = dt_to_cur_time / time_since_zero;
            (*phi[c_lev+1]).mult(cur_mult_factor);
            AddPhi(pres_fine, *phi[c_lev+1]);

            Real prev_mult_factor = dt_to_prev_time / dt_to_cur_time;
            (*phi[c_lev+1]).mult(prev_mult_factor);
            AddPhi(pres_fine_old, *phi[c_lev+1]);
        }
        else 
        {
            AddPhi(pres_fine    , *phi[c_lev+1]);
            AddPhi(pres_fine_old, *phi[c_lev+1]);
        }
    }
    //
    // Add projected vel to new velocity.
    //
    UpdateArg1(vel_crse, dt_crse, Vsync, BL_SPACEDIM, grids,      1);
    UpdateArg1(vel_fine, dt_crse, V_corr, BL_SPACEDIM, fine_grids, 1);

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	amrex::Print() << "Projection::MLsyncProject(): levels = "
                       << c_lev << ", " << c_lev+1
		       << ", time: " << run_time << '\n';
    }
}

//
// The initial velocity projection in post_init.
// this function ensures that the velocities are nondivergent
//

void
Projection::initialVelocityProject (int  c_lev,
                                    Real cur_divu_time, 
                                    int  have_divu)
{
    int lev;
    int f_lev = parent->finestLevel();

    if (verbose)
    {
      amrex::Print() << "initialVelocityProject: levels = " << c_lev
		     << "  " << f_lev << '\n';
      if (rho_wgt_vel_proj) 
	amrex::Print() << "RHO WEIGHTED INITIAL VELOCITY PROJECTION\n";
      else 
	amrex::Print() << "CONSTANT DENSITY INITIAL VELOCITY PROJECTION\n";
    }

    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real strt_time = ParallelDescriptor::second();

    Vector<MultiFab*> vel(maxlev, nullptr);
    Vector<MultiFab*> phi(maxlev, nullptr);
    Vector<std::unique_ptr<MultiFab> > sig(maxlev);

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        LevelData[lev]->get_old_data(Press_Type).setVal(0);
    }

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        vel[lev] = &(LevelData[lev]->get_new_data(State_Type));
        phi[lev] = &(LevelData[lev]->get_old_data(Press_Type));

        const int       nghost = 1;
        const BoxArray& grids  = LevelData[lev]->boxArray();
        const DistributionMapping& dmap = LevelData[lev]->DistributionMap();
        sig[lev].reset(new MultiFab(grids,dmap,1,nghost));

        if (rho_wgt_vel_proj) 
        {
            LevelData[lev]->get_new_data(State_Type).setBndry(BogusValue,Density,1);

	    AmrLevel& amr_level = parent->getLevel(lev);
	    
	    MultiFab& S_new = amr_level.get_new_data(State_Type);
	    
            Real curr_time = amr_level.get_state_data(State_Type).curTime();
	    
            for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
	      amr_level.setPhysBoundaryValues(S_new[mfi],State_Type,curr_time,
					      Density,Density,1);
            }

            MultiFab::Copy(*sig[lev],
                           LevelData[lev]->get_new_data(State_Type),
                           Density,
                           0,
                           1,
                           nghost);
        }
        else 
        {
            sig[lev]->setVal(1,nghost);
        }
    }

    Vector<std::unique_ptr<MultiFab> > rhs_cc(maxlev);
    const int nghost = 1; 

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        vel[lev]->setBndry(BogusValue,Xvel,BL_SPACEDIM);
        //
        // Set the physical boundary values.
        //
        AmrLevel& amr_level = parent->getLevel(lev);

        MultiFab& S_new = amr_level.get_new_data(State_Type);

        Real curr_time = amr_level.get_state_data(State_Type).curTime();

        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            amr_level.setPhysBoundaryValues(S_new[mfi],State_Type,curr_time,Xvel,Xvel,BL_SPACEDIM);
        }

        if (have_divu) 
        {
            int Divu_Type, Divu;
            if (!LevelData[lev]->isStateVariable("divu", Divu_Type, Divu)) 
                amrex::Error("Projection::initialVelocityProject(): Divu not found");

            // The FillBoundary seems unnecessary because
            // put_divu_in_cc_rhs will call FillPatch.  Moreover, why
            // does rhcc need to have ghost cells at all?  The solver
            // could create a temp MF with one ghost cell and it knows
            // how to properly fill ghost cell.  -- Weiqun

            //
            // Make sure ghost cells are properly filled.
            //
            MultiFab& divu_new = amr_level.get_new_data(Divu_Type);
            divu_new.FillBoundary();

            Real curr_divu_time_lev = amr_level.get_state_data(Divu_Type).curTime();

            for (MFIter mfi(divu_new); mfi.isValid(); ++mfi)
            {
                amr_level.setPhysBoundaryValues(divu_new[mfi],Divu_Type,curr_divu_time_lev,0,0,1);
            }
	    
            const BoxArray& grids     = amr_level.boxArray();
            const DistributionMapping& dmap = amr_level.DistributionMap();
            rhs_cc[lev].reset(new MultiFab(grids,dmap,1,nghost));
            put_divu_in_cc_rhs(*rhs_cc[lev],lev,cur_divu_time);
        }
    }

    if (OutFlowBC::HasOutFlowBC(phys_bc) && do_outflow_bcs && have_divu)
       set_outflow_bcs(INITIAL_VEL,phi,vel,
                       amrex::GetVecOfPtrs(rhs_cc),
                       amrex::GetVecOfPtrs(sig),
                       c_lev,f_lev,have_divu); 

     //
     // Scale the projection variables.
     //
    for (lev = c_lev; lev <= f_lev; lev++)  {
        scaleVar(INITIAL_VEL,sig[lev].get(),1,vel[lev],lev);
    }

    bool doing_initial_velproj = true;

    for (lev = f_lev-1; lev >= c_lev; --lev)
    {
        amrex::average_down(*vel[lev+1], *vel[lev], parent->Geom(lev+1), parent->Geom(lev),
                            0, BL_SPACEDIM, parent->refRatio(lev));
    }

    //
    // Project
    //
    if (!have_divu)
    {
        Vector<MultiFab*> rhs(maxlev, nullptr);
        if (use_mlmg_solver) {
            doMLMGNodalProjection(c_lev, f_lev+1, vel, phi,
                                  amrex::GetVecOfPtrs(sig),
                                  rhs, {}, 
                                  proj_tol, proj_abs_tol, 0, 0, doing_initial_velproj);
        } else {
            doNodalProjection(c_lev, f_lev+1, vel, phi,
                              amrex::GetVecOfPtrs(sig),
                              rhs, {}, 
                              proj_tol, proj_abs_tol, 0, 0, doing_initial_velproj);
        }
    } 
    else 
    {
        for (lev = c_lev; lev <= f_lev; lev++) 
        {
            if (Geometry::IsRZ()) radMultScal(lev,*rhs_cc[lev]); 
            rhs_cc[lev]->mult(-1.0,0,1,nghost);
        }

        if (use_mlmg_solver) {
            doMLMGNodalProjection(c_lev, f_lev+1, vel, phi,
                                  amrex::GetVecOfPtrs(sig),
                                  amrex::GetVecOfPtrs(rhs_cc),
                                  {},
                                  proj_tol, proj_abs_tol, 0, 0, doing_initial_velproj);
        } else {
            doNodalProjection(c_lev, f_lev+1, vel, phi,
                              amrex::GetVecOfPtrs(sig),
                              amrex::GetVecOfPtrs(rhs_cc),
                              {},
                              proj_tol, proj_abs_tol, 0, 0, doing_initial_velproj);
        }

    }

    //
    // Unscale initial projection variables.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
        rescaleVar(INITIAL_VEL,sig[lev].get(),1,vel[lev],lev);

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        LevelData[lev]->get_old_data(Press_Type).setVal(0);
        LevelData[lev]->get_new_data(Press_Type).setVal(0);
    }

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	amrex::Print() << "Projection::initialVelocityProject(): time: " << run_time << '\n';
    }
}

void
Projection::initialPressureProject (int  c_lev)
{
    int lev;
    int f_lev = parent->finestLevel();
    if (verbose) {
      amrex::Print() << "initialPressureProject: levels = " << c_lev
		     << "  " << f_lev << '\n';
    }

    Vector<MultiFab*> vel(maxlev, nullptr);
    Vector<MultiFab*> phi(maxlev, nullptr);
    Vector<std::unique_ptr<MultiFab> > sig(maxlev);

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        vel[lev] = &(LevelData[lev]->get_new_data(State_Type));
        phi[lev] = &(LevelData[lev]->get_old_data(Press_Type));

        const int       nghost = 1;
        const BoxArray& grids  = LevelData[lev]->boxArray();
        const DistributionMapping& dmap = LevelData[lev]->DistributionMap();
        sig[lev].reset(new MultiFab(grids,dmap,1,nghost));

        LevelData[lev]->get_new_data(State_Type).setBndry(BogusValue,Density,1);

        AmrLevel& amr_level = parent->getLevel(lev);

        MultiFab& S_new = amr_level.get_new_data(State_Type);

        S_new.setBndry(BogusValue,Density,1);

        Real curr_time = amr_level.get_state_data(State_Type).curTime();

        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            amr_level.setPhysBoundaryValues(S_new[mfi],State_Type,curr_time,Density,Density,1);
        }

        const Geometry& geom = parent->Geom(lev);

        MultiFab::Copy(*sig[lev],
                       LevelData[lev]->get_new_data(State_Type),
                       Density,
                       0,
                       1,
                       nghost);

	if (geom.isAnyPeriodic()) {
	    sig[lev]->FillBoundary(0,1,geom.periodicity());
	}
    }

    //
    // Set up outflow bcs.
    //
    NavierStokesBase* ns = dynamic_cast<NavierStokesBase*>(LevelData[c_lev]);
    Real gravity = ns->getGravity();

    if (OutFlowBC::HasOutFlowBC(phys_bc) && do_outflow_bcs)
    {
        int have_divu_dummy = 0;
        set_outflow_bcs(INITIAL_PRESS,phi,vel,
                        Vector<MultiFab*>(maxlev, nullptr),
                        amrex::GetVecOfPtrs(sig),
                        c_lev,f_lev,have_divu_dummy);
    }

    Vector<std::unique_ptr<MultiFab> > raii;
    for (lev = c_lev; lev <= f_lev; lev++) {
        const BoxArray& grids = vel[lev]->boxArray();
        const DistributionMapping& dmap = vel[lev]->DistributionMap();
        raii.push_back(std::unique_ptr<MultiFab>(new MultiFab(grids, dmap, BL_SPACEDIM, 1)));
        vel[lev] = raii.back().get();
        vel[lev]->setVal(0.0    , 0            , BL_SPACEDIM-1, 1);
        vel[lev]->setVal(gravity, BL_SPACEDIM-1, 1            , 1);
    }

    //
    // Scale the projection variables.
    //
    for (lev = c_lev; lev <= f_lev; lev++) {
        scaleVar(INITIAL_PRESS,sig[lev].get(),1,vel[lev],lev);
    }

    //
    // Project
    //
    Vector<MultiFab*> rhs(maxlev, nullptr);
    if (use_mlmg_solver) {
        doMLMGNodalProjection(c_lev, f_lev+1, vel, phi,
                              amrex::GetVecOfPtrs(sig),
                              rhs, {},
                              proj_tol, proj_abs_tol);
    } else {
        doNodalProjection(c_lev, f_lev+1, vel, phi,
                          amrex::GetVecOfPtrs(sig),
                          rhs, {},
                          proj_tol, proj_abs_tol);
    }

    //
    // Unscale initial projection variables.
    //
    for (lev = c_lev; lev <= f_lev; lev++) {
        rescaleVar(INITIAL_PRESS,sig[lev].get(),1,vel[lev],lev);
    }

    //
    // Copy "old" pressure just computed into "new" pressure as well.
    //
    for (lev = c_lev; lev <= f_lev; lev++) {
        MultiFab::Copy(LevelData[lev]->get_new_data(Press_Type),
                       LevelData[lev]->get_old_data(Press_Type),
                       0, 0, 1, 0);
    }
}

//
// The velocity projection in post_init, which computes the initial
// pressure used in the timestepping.
//

void
Projection::initialSyncProject (int       c_lev,
                                const Vector<MultiFab*> sig,
                                Real      dt, 
                                Real      strt_time,
                                int       have_divu)
{
    int lev;
    int f_lev = parent->finestLevel();

    if (verbose)
      amrex::Print() << "initialSyncProject: levels = " << c_lev << "  " << f_lev << '\n';

    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real stime = ParallelDescriptor::second();

    //
    // Gather data.
    //
    Vector<MultiFab*> vel(maxlev, nullptr);
    Vector<MultiFab*> phi(maxlev, nullptr);
    Vector<std::unique_ptr<MultiFab> > rhs(maxlev);

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        vel[lev] = &(LevelData[lev]->get_new_data(State_Type));
        phi[lev] = &(LevelData[lev]->get_old_data(Press_Type));
    }
  
    const Real dt_inv = 1./dt;

    if (have_divu) 
    {
        //
        // Set up rhs for manual project.
        //
        for (lev = c_lev; lev <= f_lev; lev++) 
        {
            AmrLevel& amr_level = parent->getLevel(lev);

            int Divu_Type, Divu;
            if (!LevelData[c_lev]->isStateVariable("divu", Divu_Type, Divu)) 
                amrex::Error("Projection::initialSyncProject(): Divu not found");
            //
            // Make sure ghost cells are properly filled.
            //
            MultiFab& divu_new = amr_level.get_new_data(Divu_Type);
            MultiFab& divu_old = amr_level.get_old_data(Divu_Type);
            divu_new.FillBoundary();
            divu_old.FillBoundary();

            Real prev_time = amr_level.get_state_data(Divu_Type).prevTime();
            Real curr_time = amr_level.get_state_data(Divu_Type).curTime();

            for (MFIter mfi(divu_new); mfi.isValid(); ++mfi)
            {
                amr_level.setPhysBoundaryValues(divu_old[mfi],Divu_Type,prev_time,0,0,1);
                amr_level.setPhysBoundaryValues(divu_new[mfi],Divu_Type,curr_time,0,0,1);
            }

            const int nghost = 1;
            rhs[lev].reset(new MultiFab(amr_level.boxArray(),
                                        amr_level.DistributionMap(),
                                        1,nghost));
            MultiFab* rhslev = rhs[lev].get();
            rhslev->setVal(0);

            NavierStokesBase* ns = dynamic_cast<NavierStokesBase*>(&parent->getLevel(lev));

            BL_ASSERT(!(ns == 0));

            std::unique_ptr<MultiFab> divu (ns->getDivCond(nghost,strt_time));
            std::unique_ptr<MultiFab> dsdt (ns->getDivCond(nghost,strt_time+dt));

            for (MFIter mfi(*rhslev); mfi.isValid(); ++mfi)
            {
                FArrayBox& dsdtfab = (*dsdt)[mfi];
                dsdtfab.minus((*divu)[mfi]);
                dsdtfab.mult(dt_inv);
                (*rhslev)[mfi].copy(dsdtfab);
            }
        }
    }

    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        MultiFab& P_old = LevelData[lev]->get_old_data(Press_Type);
        P_old.setVal(0);
    }
    //
    // Set velocity bndry values to bogus values.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        vel[lev]->setBndry(BogusValue,Xvel,BL_SPACEDIM);
        MultiFab &u_o = LevelData[lev]->get_old_data(State_Type);
        u_o.setBndry(BogusValue,Xvel,BL_SPACEDIM);
        sig[lev]->setBndry(BogusValue);
    }
    //
    // Convert velocities to accelerations (we always do this for the
    //  projections in these initial iterations).
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        MultiFab& S_old = LevelData[lev]->get_old_data(State_Type);
        MultiFab& S_new = LevelData[lev]->get_new_data(State_Type);

        Real prev_time = LevelData[lev]->get_state_data(State_Type).prevTime();
        Real curr_time = LevelData[lev]->get_state_data(State_Type).curTime();

        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            LevelData[lev]->setPhysBoundaryValues(S_old[mfi],State_Type,prev_time,Xvel,Xvel,BL_SPACEDIM);
            LevelData[lev]->setPhysBoundaryValues(S_new[mfi],State_Type,curr_time,Xvel,Xvel,BL_SPACEDIM);
        }

        MultiFab& u_o = LevelData[lev]->get_old_data(State_Type);
        ConvertUnew(*vel[lev], u_o, dt, LevelData[lev]->boxArray());
    }

    if (OutFlowBC::HasOutFlowBC(phys_bc) && have_divu && do_outflow_bcs) {
        set_outflow_bcs(INITIAL_SYNC,phi,vel,
                        amrex::GetVecOfPtrs(rhs),
                        sig,c_lev,f_lev,have_divu);
    }

    //
    // Scale initial sync projection variables.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
    {
        scaleVar(INITIAL_SYNC,sig[lev],1,vel[lev],lev);

        if (have_divu && Geometry::IsRZ()) 
          radMultScal(lev,*(rhs[lev]));
    }

    for (lev = f_lev; lev >= c_lev+1; lev--) {
      const BoxArray& crse_grids = vel[lev-1]->boxArray();
      const BoxArray& fine_grids = vel[lev  ]->boxArray();
      const DistributionMapping& crse_dmap = vel[lev-1]->DistributionMap();
      const DistributionMapping& fine_dmap = vel[lev  ]->DistributionMap();

      MultiFab v_crse(crse_grids, crse_dmap, BL_SPACEDIM, 1);
      MultiFab v_fine(fine_grids, fine_dmap, BL_SPACEDIM, 1);

      const Geometry& fine_geom = parent->Geom(lev  );
      const Geometry& crse_geom = parent->Geom(lev-1);

      MultiFab::Copy(v_crse, *vel[lev-1], 0, 0, BL_SPACEDIM, 1);
      MultiFab::Copy(v_fine, *vel[lev  ], 0, 0, BL_SPACEDIM, 1);

      // restrict_level(v_crse, v_fine, parent->refRatio(lev-1));
      amrex::average_down(v_fine,v_crse,fine_geom,crse_geom,
                           0, v_crse.nComp(), parent->refRatio(lev-1));
	
      MultiFab::Copy(*vel[lev-1], v_crse, 0, 0, BL_SPACEDIM, 1);
    }

    //
    // Project.
    //
    if (have_divu) {
        for (lev = c_lev; lev <= f_lev; lev++) 
        {
            rhs[lev]->mult(-1.0,0,1);
        }
    }

    if (use_mlmg_solver) {
        doMLMGNodalProjection(c_lev, f_lev+1, vel, phi, sig,
                              amrex::GetVecOfPtrs(rhs),
                              {}, proj_tol, proj_abs_tol);
    } else {
        doNodalProjection(c_lev, f_lev+1, vel, phi, sig,
                          amrex::GetVecOfPtrs(rhs),
                          {}, proj_tol, proj_abs_tol);
    }

    //
    // Unscale initial sync projection variables.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
        rescaleVar(INITIAL_SYNC,sig[lev],1,vel[lev],lev);
    //
    // Add correction at coarse and fine levels.
    //
    for (lev = c_lev; lev <= f_lev; lev++) 
        incrPress(lev, 1.0);

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - stime;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	amrex::Print() << "Projection::initialSyncProject(): time: " << run_time << '\n';
    }
}

//
// Put S in the rhs of the projector--cell based version.
//

void
Projection::put_divu_in_cc_rhs (MultiFab&       rhs,
                                int             level,
                                Real            time)
{
    rhs.setVal(0);

    NavierStokesBase* ns = dynamic_cast<NavierStokesBase*>(&parent->getLevel(level));

    BL_ASSERT(!(ns == 0));

    std::unique_ptr<MultiFab> divu (ns->getDivCond(1,time));

    for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
    {
        rhs[mfi].copy((*divu)[mfi]);
    }
}

//
// Convert U from an Accl-like quantity to a velocity: Unew = Uold + alpha*Unew
//

void
Projection::UnConvertUnew (MultiFab&       Uold,
                           Real            alpha,
                           MultiFab&       Unew, 
                           const BoxArray& grids)
{
    for (MFIter Uoldmfi(Uold); Uoldmfi.isValid(); ++Uoldmfi) 
    {
        BL_ASSERT(grids[Uoldmfi.index()] == Uoldmfi.validbox());

        UnConvertUnew(Uold[Uoldmfi],alpha,Unew[Uoldmfi],Uoldmfi.validbox());
    }
}

//
// Convert U from an Accleration like quantity to a velocity
// Unew = Uold + alpha*Unew.
//

void
Projection::UnConvertUnew (FArrayBox& Uold,
                           Real       alpha,
                           FArrayBox& Unew,
                           const Box& grd)
{
    BL_ASSERT(Unew.nComp() >= BL_SPACEDIM);
    BL_ASSERT(Uold.nComp() >= BL_SPACEDIM);
    BL_ASSERT(Unew.contains(grd) == true);
    BL_ASSERT(Uold.contains(grd) == true);
    
    const int*  lo    = grd.loVect();
    const int*  hi    = grd.hiVect();
    const int*  uo_lo = Uold.loVect(); 
    const int*  uo_hi = Uold.hiVect(); 
    const Real* uold  = Uold.dataPtr(0);
    const int*  un_lo = Unew.loVect(); 
    const int*  un_hi = Unew.hiVect(); 
    const Real* unew  = Unew.dataPtr(0);
    
    FORT_ACCEL_TO_VEL(lo, hi,
                      uold, ARLIM(uo_lo), ARLIM(uo_hi),
                      &alpha,
                      unew, ARLIM(un_lo), ARLIM(un_hi));
}

//
// Convert U to an Accleration like quantity: Unew = (Unew - Uold)/alpha
//

void
Projection::ConvertUnew (MultiFab&       Unew,
                         MultiFab&       Uold,
                         Real            alpha,
                         const BoxArray& grids)
{
    for (MFIter Uoldmfi(Uold); Uoldmfi.isValid(); ++Uoldmfi) 
    {
        BL_ASSERT(grids[Uoldmfi.index()] == Uoldmfi.validbox());

        ConvertUnew(Unew[Uoldmfi],Uold[Uoldmfi],alpha,Uoldmfi.validbox());
    }
}

//
// Convert U to an Accleration like quantity: Unew = (Unew - Uold)/alpha
//

void
Projection::ConvertUnew( FArrayBox &Unew, FArrayBox &Uold, Real alpha,
                              const Box &grd )
{
    BL_ASSERT(Unew.nComp() >= BL_SPACEDIM);
    BL_ASSERT(Uold.nComp() >= BL_SPACEDIM);
    BL_ASSERT(Unew.contains(grd) == true);
    BL_ASSERT(Uold.contains(grd) == true);
    
    const int*  lo    = grd.loVect();
    const int*  hi    = grd.hiVect();
    const int*  uo_lo = Uold.loVect(); 
    const int*  uo_hi = Uold.hiVect(); 
    const Real* uold  = Uold.dataPtr(0);
    const int*  un_lo = Unew.loVect(); 
    const int*  un_hi = Unew.hiVect(); 
    const Real* unew  = Unew.dataPtr(0);
                    
    FORT_VEL_TO_ACCEL(lo, hi, 
                      unew, ARLIM(un_lo), ARLIM(un_hi),
                      uold, ARLIM(uo_lo), ARLIM(uo_hi), &alpha );
}

//
// Update a quantity U using the formula: Unew = Unew + alpha*Uold
//

void
Projection::UpdateArg1 (MultiFab&       Unew,
                        Real            alpha,
                        MultiFab&       Uold,
                        int             nvar,
                        const BoxArray& grids,
                        int             ngrow)
{
    for (MFIter Uoldmfi(Uold); Uoldmfi.isValid(); ++Uoldmfi) 
    {
        BL_ASSERT(grids[Uoldmfi.index()] == Uoldmfi.validbox());

        UpdateArg1(Unew[Uoldmfi],alpha,Uold[Uoldmfi],nvar,Uoldmfi.validbox(),ngrow);
    }
}

//
// Update a quantity U using the formula
// currently only the velocity, but will do the pressure as well.
// Unew = Unew + alpha*Uold
//

void
Projection::UpdateArg1 (FArrayBox& Unew,
                        Real       alpha,
                        FArrayBox& Uold,
                        int        nvar,
                        const Box& grd,
                        int        ngrow)
{
    BL_ASSERT(nvar <= Uold.nComp());
    BL_ASSERT(nvar <= Unew.nComp());

    Box        b  = amrex::grow(grd,ngrow);
    const Box& bb = Unew.box();

    if (bb.ixType() == IndexType::TheNodeType())
        b.surroundingNodes();

    BL_ASSERT(Uold.contains(b) == true);
    BL_ASSERT(Unew.contains(b) == true);

    const int*  lo    = b.loVect();
    const int*  hi    = b.hiVect();
    const int*  uo_lo = Uold.loVect(); 
    const int*  uo_hi = Uold.hiVect(); 
    const Real* uold  = Uold.dataPtr(0);
    const int*  un_lo = Unew.loVect(); 
    const int*  un_hi = Unew.hiVect(); 
    const Real* unew  = Unew.dataPtr(0);
                    
    FORT_PROJ_UPDATE(lo,hi,&nvar,&ngrow,
                     unew, ARLIM(un_lo), ARLIM(un_hi),
                     &alpha,
                     uold, ARLIM(uo_lo), ARLIM(uo_hi) );
}

//
// Add phi to P.
//

void
Projection::AddPhi (MultiFab&        p,
                    MultiFab&       phi)
{
    for (MFIter pmfi(p); pmfi.isValid(); ++pmfi) 
    {
        p[pmfi].plus(phi[pmfi]);
    }
}

//
// Convert phi into p^n+1/2.
//

void
Projection::incrPress (int  level,
                       Real dt)
{
    MultiFab& P_old = LevelData[level]->get_old_data(Press_Type);
    MultiFab& P_new = LevelData[level]->get_new_data(Press_Type);

    const BoxArray& grids = LevelData[level]->boxArray();

    for (MFIter P_newmfi(P_new); P_newmfi.isValid(); ++P_newmfi)
    {
        const int i = P_newmfi.index();

        UpdateArg1(P_new[P_newmfi],1.0/dt,P_old[P_newmfi],1,grids[i],1);

        P_old[P_newmfi].setVal(BogusValue);
    }
}

//
// This function scales variables at the start of a projection.
//

void
Projection::scaleVar (int             which_call,
                      MultiFab*       sig,
                      int             sig_nghosts,
                      MultiFab*       vel,
                      int             level)
{
    BL_ASSERT((which_call == INITIAL_VEL  ) || 
              (which_call == INITIAL_PRESS) || 
              (which_call == INITIAL_SYNC ) ||
              (which_call == LEVEL_PROJ   ) ||
              (which_call == SYNC_PROJ    ) );

    if (sig != 0)
        BL_ASSERT(sig->nComp() == 1);
    if (vel != 0)
        BL_ASSERT(vel->nComp() >= BL_SPACEDIM);

    //
    // Convert sigma from rho to anel_coeff/rho if not INITIAL_PRESS.
    // nghosts info needed to avoid divide by zero.
    //
    if (sig != 0) {
      sig->invert(1.0,sig_nghosts);
      if (which_call  != INITIAL_PRESS &&
          anel_coeff[level] != 0) AnelCoeffMult(level,*sig,0);
    }
    

    //
    // Scale by radius for RZ.
    //
    if (Geometry::IsRZ()) 
    {
        if (sig != 0)
            radMultScal(level,*sig);
        if (vel != 0)
            radMultVel(level,*vel);
    }

    //
    // Scale velocity by anel_coeff if it exists
    //
    if (vel != 0 && anel_coeff[level] != 0)
      for (int n = 0; n < BL_SPACEDIM; n++) 
        AnelCoeffMult(level,*vel,n);
}

//
// This function rescales variables at the end of a projection.
//

void
Projection::rescaleVar (int             which_call,
                        MultiFab*       sig,
                        int             sig_nghosts,
                        MultiFab*       vel,
                        int             level)
{
    BL_ASSERT((which_call == INITIAL_VEL  ) || 
              (which_call == INITIAL_PRESS) || 
              (which_call == INITIAL_SYNC ) ||
              (which_call == LEVEL_PROJ   ) ||
              (which_call == SYNC_PROJ    ) );

    if (sig != 0)
        BL_ASSERT(sig->nComp() == 1);
    if (vel != 0)
        BL_ASSERT(vel->nComp() >= BL_SPACEDIM);

    if (which_call  != INITIAL_PRESS && sig != 0 &&
        anel_coeff[level] != 0) AnelCoeffDiv(level,*sig,0);
    //
    // Divide by radius to rescale for RZ coordinates.
    //
    if (Geometry::IsRZ()) 
    {
        if (sig != 0)
            radDiv(level,*sig,0);
        if (vel != 0)
        {
            for (int n = 0; n < BL_SPACEDIM; n++)
                radDiv(level,*vel,n);
        }
    }
    if (vel != 0 && anel_coeff[level] != 0) 
      for (int n = 0; n < BL_SPACEDIM; n++)
        AnelCoeffDiv(level,*vel,n);
    //
    // Convert sigma from 1/rho to rho
    // NOTE: this must come after division by r to be correct,
    // nghosts info needed to avoid divide by zero.
    //
    if (sig != 0)
        sig->invert(1.0,sig_nghosts);
}

//
// Multiply by a radius for r-z coordinates.
//

void
Projection::radMultScal (int       level,
                         MultiFab& mf)
{
#if (BL_SPACEDIM < 3)
    BL_ASSERT(radius_grow >= mf.nGrow());

    int ngrow = mf.nGrow();

    int nr = radius_grow;

    const Box& domain = parent->Geom(level).Domain();
    const int* domlo  = domain.loVect();
    const int* domhi  = domain.hiVect();

    Real bogus_value = BogusValue;

    for (MFIter mfmfi(mf); mfmfi.isValid(); ++mfmfi) 
    {
        BL_ASSERT(mf.box(mfmfi.index()) == mfmfi.validbox());

        const Box& bx = mfmfi.validbox();
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        Real* dat     = mf[mfmfi].dataPtr(0);
        Real* rad     = &(*radius[level])[mfmfi.index()][0];

        FORT_RADMPYSCAL(dat,ARLIM(lo),ARLIM(hi),domlo,domhi,&ngrow,
                        rad,&nr,&bogus_value);
    }
#endif
}

void
Projection::radMultVel (int       level,
                        MultiFab& mf)
{
#if (BL_SPACEDIM < 3)
    BL_ASSERT(radius_grow >= mf.nGrow());

    int ngrow = mf.nGrow();

    int nr = radius_grow;

    const Box& domain = parent->Geom(level).Domain();
    const int* domlo  = domain.loVect();
    const int* domhi  = domain.hiVect();

    Real bogus_value = BogusValue;

    for (int n = 0; n < BL_SPACEDIM; n++) 
    {
       for (MFIter mfmfi(mf); mfmfi.isValid(); ++mfmfi) 
       {
           BL_ASSERT(mf.box(mfmfi.index()) == mfmfi.validbox());
   
           const Box& bx = mfmfi.validbox();
           const int* lo = bx.loVect();
           const int* hi = bx.hiVect();
           Real* dat     = mf[mfmfi].dataPtr(n);
           Real* rad     = &(*radius[level])[mfmfi.index()][0];

           FORT_RADMPYVEL(dat,ARLIM(lo),ARLIM(hi),domlo,domhi,&ngrow,
                          rad,&nr,&n);
       }
    }
#endif
}

//
// Divide by a radius for r-z coordinates.
//

void
Projection::radDiv (int       level,
                    MultiFab& mf,
                    int       comp)
{
#if (BL_SPACEDIM < 3)
    BL_ASSERT(comp >= 0 && comp < mf.nComp());
    BL_ASSERT(radius_grow >= mf.nGrow());

    int ngrow = mf.nGrow();
    int nr    = radius_grow;

    const Box& domain = parent->Geom(level).Domain();
    const int* domlo  = domain.loVect();
    const int* domhi  = domain.hiVect();

    Real bogus_value = BogusValue;

    for (MFIter mfmfi(mf); mfmfi.isValid(); ++mfmfi) 
    {
        BL_ASSERT(mf.box(mfmfi.index()) == mfmfi.validbox());

        const Box& bx  = mfmfi.validbox();
        const int* lo  = bx.loVect();
        const int* hi  = bx.hiVect();
        Real*      dat = mf[mfmfi].dataPtr(comp);
        Real*      rad = &(*radius[level])[mfmfi.index()][0];

        FORT_RADDIV(dat,ARLIM(lo),ARLIM(hi),domlo,domhi,&ngrow,
                    rad,&nr,&bogus_value);
    }
#endif
}

//
// Multiply by anel_coeff if it is defined
//
void
Projection::AnelCoeffMult (int       level,
                           MultiFab& mf,
                           int       comp)
{
    BL_ASSERT(anel_coeff[level] != 0);
    BL_ASSERT(comp >= 0 && comp < mf.nComp());
    int ngrow = mf.nGrow();
    int nr    = 1;

    const Box& domain = parent->Geom(level).Domain();
    const int* domlo  = domain.loVect();
    const int* domhi  = domain.hiVect();

    Real bogus_value = BogusValue;

    int mult = 1;

    for (MFIter mfmfi(mf); mfmfi.isValid(); ++mfmfi) 
    {
        BL_ASSERT(mf.box(mfmfi.index()) == mfmfi.validbox());

        const Box& bx = mfmfi.validbox();
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        Real* dat     = mf[mfmfi].dataPtr(comp);

        FORT_ANELCOEFFMPY(dat,ARLIM(lo),ARLIM(hi),domlo,domhi,&ngrow,
                          anel_coeff[level][mfmfi.index()],&nr,&bogus_value,&mult);
    }
}

//
// Divide by anel_coeff if it is defined
//
void
Projection::AnelCoeffDiv (int       level,
                          MultiFab& mf,
                          int       comp)
{
    BL_ASSERT(comp >= 0 && comp < mf.nComp());
    BL_ASSERT(anel_coeff[level] != 0);
    int ngrow = mf.nGrow();
    int nr    = 1;

    const Box& domain = parent->Geom(level).Domain();
    const int* domlo  = domain.loVect();
    const int* domhi  = domain.hiVect();

    Real bogus_value = BogusValue;

    int mult = 0;

    for (MFIter mfmfi(mf); mfmfi.isValid(); ++mfmfi) 
    {
        BL_ASSERT(mf.box(mfmfi.index()) == mfmfi.validbox());

        const Box& bx = mfmfi.validbox();
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        Real* dat     = mf[mfmfi].dataPtr(comp);

        FORT_ANELCOEFFMPY(dat,ARLIM(lo),ARLIM(hi),domlo,domhi,&ngrow,
                          anel_coeff[level][mfmfi.index()],&nr,&bogus_value,&mult);
    }
}

//
// This projects the initial vorticity field (stored in pressure)
// to define an initial velocity field.
//
void
Projection::initialVorticityProject (int c_lev)
{
#if (BL_SPACEDIM == 2)
  int f_lev = parent->finestLevel();

  if (verbose)
  {
    amrex::Print() << "initialVorticityProject: levels = " << c_lev
		   << "  " << f_lev << std::endl;
  }
    //
    // Set up projector bndry just for this projection.
    //
    const Geometry& geom = parent->Geom(0);

    Vector<std::unique_ptr<MultiFab> > p_real(maxlev);
    Vector<std::unique_ptr<MultiFab> > s_real(maxlev);

    for (int lev = c_lev; lev <= f_lev; lev++)
    {
        MultiFab& P_new  = LevelData[lev]->get_new_data(Press_Type);
        const int nghost = 1;
        s_real[lev].reset(new MultiFab(LevelData[lev]->boxArray(),
                                       LevelData[lev]->DistributionMap(),
                                       1,nghost));
        s_real[lev]->setVal(1,nghost);
        p_real[lev].reset(new MultiFab(P_new.boxArray(),
                                       P_new.DistributionMap(),
                                       1,nghost));
        p_real[lev]->setVal(0,nghost);
    }
    //
    // Set up outflow bcs.
    //
    Vector<std::unique_ptr<MultiFab> > u_real(maxlev);
    Vector<std::unique_ptr<MultiFab> > rhnd(maxlev);

    for (int lev = c_lev; lev <= f_lev; lev++)
    {
        u_real[lev].reset(new MultiFab(parent->getLevel(lev).boxArray(),
                                       parent->getLevel(lev).DistributionMap(),
                                       BL_SPACEDIM, 1));
        u_real[lev]->setVal(0);
        //
        // The vorticity is stored in the new pressure variable for now.
        //
        MultiFab& P_new = LevelData[lev]->get_new_data(Press_Type);

        rhnd[lev].reset(new MultiFab(P_new.boxArray(), 
                                     P_new.DistributionMap(),
                                     1, 0));

        for (MFIter mfi(*rhnd[lev]); mfi.isValid(); ++mfi)
        {
            (*rhnd[lev])[mfi].setVal(0);
            (*rhnd[lev])[mfi].copy(P_new[mfi], 0, 0);
        }
    }

    //
    // Set BC for vorticity solve, save a copy of orig ones
    //
    BCRec phys_bc_save(phys_bc->lo(),phys_bc->hi());
    for (int i=0; i<BL_SPACEDIM; ++i) {
      phys_bc->setLo(i,Outflow);
      phys_bc->setHi(i,Outflow);
      if (geom.isPeriodic(i)) {
        phys_bc->setLo(i,Interior);
        phys_bc->setHi(i,Interior);
      }
    }
    //
    // Project.
    //
    if (use_mlmg_solver) {
        doMLMGNodalProjection(c_lev, f_lev+1,
                              amrex::GetVecOfPtrs(u_real), 
                              amrex::GetVecOfPtrs(p_real),
                              amrex::GetVecOfPtrs(s_real),
                              Vector<MultiFab*>(maxlev, nullptr),
                              amrex::GetVecOfPtrs(rhnd),
                              proj_tol, proj_abs_tol,
                              nullptr, nullptr, false, true);
    } else {
        doNodalProjection(c_lev, f_lev+1,
                          amrex::GetVecOfPtrs(u_real), 
                          amrex::GetVecOfPtrs(p_real),
                          amrex::GetVecOfPtrs(s_real),
                          Vector<MultiFab*>(maxlev, nullptr),
                          amrex::GetVecOfPtrs(rhnd),
                          proj_tol, proj_abs_tol);
    }

    //
    // Generate velocity field from potential
    //
    const int idx[2] = {1, 0};

    Vector<MultiFab*> vel(maxlev, nullptr);
    for (int lev = c_lev; lev <= f_lev; lev++)
    {
        vel[lev] = &(LevelData[lev]->get_new_data(State_Type));
        //
        // Note: Here u_real from projection is -grad(phi), but if
        //  phi is the stream function, u=dphi/dy, v=-dphi/dx
        //
        (*u_real[lev]).mult(-1,Yvel,1);

        for (int n = 0; n < BL_SPACEDIM; n++)
        {
            for (MFIter mfi(*vel[lev]); mfi.isValid(); ++mfi)
            {
                const Box& box = mfi.validbox();
                if (add_vort_proj)
                {
                  (*vel[lev])[mfi].plus((*u_real[lev])[mfi],box,Xvel+n,Xvel+idx[n], 1);
                }
                else
                {
                  (*vel[lev])[mfi].copy((*u_real[lev])[mfi],box,Xvel+n,box,Xvel+idx[n], 1);
                }
            }
        }
    }

    //
    // Restore bcs
    //
    for (int i=0; i<BL_SPACEDIM; ++i) {
      phys_bc->setLo(i,phys_bc_save.lo()[i]);
      phys_bc->setHi(i,phys_bc_save.hi()[i]);
    }


#else
    amrex::Error("Projection::initialVorticityProject(): not implented yet for 3D");
#endif
}

void 
Projection::putDown (const Vector<MultiFab*>& phi,
                     FArrayBox*         phi_fine_strip,
                     int                c_lev,
                     int                f_lev,
                     const Orientation* outFaces,
                     int                numOutFlowFaces,
                     int                ncStripWidth)
{
    BL_PROFILE("Projection::putDown()");
    //
    // Put down to coarser levels.
    //
    const int nCompPhi = 1; // phi_fine_strip.nComp();
    const int nGrow    = 0; // phi_fine_strip.nGrow();
    IntVect ratio      = IntVect::TheUnitVector();

    for (int lev = f_lev-1; lev >= c_lev; lev--)
    {
        ratio *= parent->refRatio(lev);
        const Box& domainC = parent->Geom(lev).Domain();

        for (int iface = 0; iface < numOutFlowFaces; iface++) 
        {
            Box phiC_strip = 
                amrex::surroundingNodes(amrex::bdryNode(domainC, outFaces[iface], ncStripWidth));
            phiC_strip.grow(nGrow);
            BoxArray ba(phiC_strip);
            DistributionMapping dm{ba};
            MultiFab phi_crse_strip(ba, dm, nCompPhi, 0);
            phi_crse_strip.setVal(0);

            for (MFIter mfi(phi_crse_strip); mfi.isValid(); ++mfi)
            {
                Box ovlp = amrex::coarsen(phi_fine_strip[iface].box(),ratio) & mfi.validbox();

                if (ovlp.ok())
                {
                    FArrayBox& cfab = phi_crse_strip[mfi];
                    FORT_PUTDOWN (BL_TO_FORTRAN(cfab),
                                  BL_TO_FORTRAN(phi_fine_strip[iface]),
                                  ovlp.loVect(), ovlp.hiVect(), ratio.getVect());
                }
            }

            phi[lev]->copy(phi_crse_strip);
        }
    }
}

void
Projection::getStreamFunction (Vector<std::unique_ptr<MultiFab> >& phi)
{
  amrex::Abort("Projection::getStreamFunction not implemented");
}

//
// Given a nodal pressure P compute the pressure gradient at the
// contained cell centers.

void
Projection::getGradP (FArrayBox& p_fab,
                      FArrayBox& gp,
                      const Box& gpbox_to_fill,
                      const Real* dx)
{
    BL_PROFILE("Projection::getGradP()");
    //
    // Test to see if p_fab contains gpbox_to_fill
    //
    BL_ASSERT(amrex::enclosedCells(p_fab.box()).contains(gpbox_to_fill));

    const int*  plo    = p_fab.loVect();
    const int*  phi    = p_fab.hiVect();
    const int*  glo    = gp.box().loVect();
    const int*  ghi    = gp.box().hiVect();
    const int*   lo    = gpbox_to_fill.loVect();
    const int*   hi    = gpbox_to_fill.hiVect();
    const Real* p_dat  = p_fab.dataPtr();
    const Real* gp_dat = gp.dataPtr();

#if (BL_SPACEDIM == 2)
    int is_full = 0;
    FORT_GRADP(p_dat,ARLIM(plo),ARLIM(phi),gp_dat,ARLIM(glo),ARLIM(ghi),lo,hi,dx,
               &is_full);
#elif (BL_SPACEDIM == 3)
    FORT_GRADP(p_dat,ARLIM(plo),ARLIM(phi),gp_dat,ARLIM(glo),ARLIM(ghi),lo,hi,dx);
#endif
}

void
Projection::set_outflow_bcs (int        which_call,
                             const Vector<MultiFab*>& phi,
                             const Vector<MultiFab*>& Vel_in,
                             const Vector<MultiFab*>& Divu_in,
                             const Vector<MultiFab*>& Sig_in,
                             int        c_lev,
                             int        f_lev,
                             int        have_divu)
{
    BL_ASSERT((which_call == INITIAL_VEL  ) || 
              (which_call == INITIAL_PRESS) || 
              (which_call == INITIAL_SYNC ) ||
              (which_call == LEVEL_PROJ   ) );

    if (which_call != LEVEL_PROJ)
      BL_ASSERT(c_lev == 0);

    if (verbose)
      amrex::Print() << "...setting outflow bcs for the nodal projection ... " << '\n';

    bool        hasOutFlow;
    Orientation outFaces[2*BL_SPACEDIM];
    Orientation outFacesAtThisLevel[maxlev][2*BL_SPACEDIM];

    int fine_level[2*BL_SPACEDIM];

    int numOutFlowFacesAtAllLevels;
    int numOutFlowFaces[maxlev];
    OutFlowBC::GetOutFlowFaces(hasOutFlow,outFaces,phys_bc,numOutFlowFacesAtAllLevels);

    //
    // Get 2-wide cc box, state_strip, along interior of top. 
    // Get 1-wide nc box, phi_strip  , along top.
    //
    const int ccStripWidth = 2;

//    const int nCompPhi    = 1;
//    const int srcCompVel  = Xvel;
//    const int srcCompDivu = 0;
//    const int   nCompVel  = BL_SPACEDIM;
//    const int   nCompDivu = 1;

    //
    // Determine the finest level such that the entire outflow face is covered
    // by boxes at this level (skip if doesnt touch, and bomb if only partially
    // covered).
    //
    Box state_strip[maxlev][2*BL_SPACEDIM];

    int icount[maxlev];
    for (int i=0; i < maxlev; i++) icount[i] = 0;

    //
    // This loop is only to define the number of outflow faces at each level.
    //
    Box temp_state_strip;
    for (int iface = 0; iface < numOutFlowFacesAtAllLevels; iface++) 
    {
      const int outDir    = outFaces[iface].coordDir();

      fine_level[iface] = -1;
      for (int lev = f_lev; lev >= c_lev; lev--)
      {
        Box domain = parent->Geom(lev).Domain();

        if (outFaces[iface].faceDir() == Orientation::high)
        {
            temp_state_strip = amrex::adjCellHi(domain,outDir,ccStripWidth);
            temp_state_strip.shift(outDir,-ccStripWidth);
        }
        else
        {
            temp_state_strip = amrex::adjCellLo(domain,outDir,ccStripWidth);
            temp_state_strip.shift(outDir,ccStripWidth);
        }
        // Grow the box by one tangentially in order to get velocity bc's.
        for (int dir = 0; dir < BL_SPACEDIM; dir++) 
          if (dir != outDir) temp_state_strip.grow(dir,1);

        const BoxArray& Lgrids               = parent->getLevel(lev).boxArray();
        const Box&      valid_state_strip    = temp_state_strip & domain;
        const BoxArray  uncovered_outflow_ba = amrex::complementIn(valid_state_strip,Lgrids);

        BL_ASSERT( !(uncovered_outflow_ba.size() &&
                     amrex::intersect(Lgrids,valid_state_strip).size()) );

        if ( !(uncovered_outflow_ba.size()) && fine_level[iface] == -1) {
            int ii = icount[lev];
            outFacesAtThisLevel[lev][ii] = outFaces[iface];
            state_strip[lev][ii] = temp_state_strip;
            fine_level[iface] = lev;
            icount[lev]++;
        }
      }
    }

    for (int lev = f_lev; lev >= c_lev; lev--) {
      numOutFlowFaces[lev] = icount[lev];
    }

    NavierStokesBase* ns0 = dynamic_cast<NavierStokesBase*>(LevelData[c_lev]);
    BL_ASSERT(!(ns0 == 0));
   
    int Divu_Type, Divu;
    Real gravity = 0;

    if (which_call == INITIAL_SYNC || which_call == INITIAL_VEL)
    {
      gravity = 0;
      if (!LevelData[c_lev]->isStateVariable("divu", Divu_Type, Divu))
        amrex::Error("Projection::set_outflow_bcs: No divu.");
    }

    if (which_call == INITIAL_PRESS || which_call == LEVEL_PROJ)
    {
      gravity = ns0->getGravity();
      if (!LevelData[c_lev]->isStateVariable("divu", Divu_Type, Divu) &&
          (gravity == 0) )
        amrex::Error("Projection::set_outflow_bcs: No divu or gravity.");
    }

    for (int lev = c_lev; lev <= f_lev; lev++) 
    {
      if (numOutFlowFaces[lev] > 0) 
        set_outflow_bcs_at_level (which_call,lev,c_lev,
                                  state_strip[lev],
                                  outFacesAtThisLevel[lev],
                                  numOutFlowFaces[lev],
                                  phi,
                                  Vel_in[lev],
                                  Divu_in[lev],
                                  Sig_in[lev],
                                  have_divu,
                                  gravity);
                                  
    }

}

void
Projection::set_outflow_bcs_at_level (int          which_call,
                                      int          lev,
                                      int          c_lev,
                                      Box*         state_strip,
                                      Orientation* outFacesAtThisLevel,
                                      int          numOutFlowFaces,
                                      const Vector<MultiFab*>&  phi, 
                                      MultiFab*    Vel_in,
                                      MultiFab*    Divu_in,
                                      MultiFab*    Sig_in,
                                      int          have_divu,
                                      Real         gravity)
{
    BL_ASSERT(dynamic_cast<NavierStokesBase*>(LevelData[lev]) != nullptr);

    Box domain = parent->Geom(lev).Domain();

    const int ncStripWidth = 1;

    FArrayBox  rho[2*BL_SPACEDIM];
    FArrayBox dsdt[2*BL_SPACEDIM];
    FArrayBox dudt[1][2*BL_SPACEDIM];
    FArrayBox phi_fine_strip[2*BL_SPACEDIM];

    const int ngrow = 1;

    for (int iface = 0; iface < numOutFlowFaces; iface++)
    {
        dsdt[iface].resize(state_strip[iface],1);
        dudt[0][iface].resize(state_strip[iface],BL_SPACEDIM);

        rho[iface].resize(state_strip[iface],1);

        (*Sig_in).copyTo(rho[iface],0,0,1,ngrow);

        Box phi_strip = 
            amrex::surroundingNodes(amrex::bdryNode(domain,
                                                      outFacesAtThisLevel[iface],
                                                      ncStripWidth));
        phi_fine_strip[iface].resize(phi_strip,1);
        phi_fine_strip[iface].setVal(0.);
    }

    ProjOutFlowBC projBC;
    if (which_call == INITIAL_PRESS) 
    {

        const int*      lo_bc = phys_bc->lo();
        const int*      hi_bc = phys_bc->hi();
        projBC.computeRhoG(rho,phi_fine_strip,
                           parent->Geom(lev),
                           outFacesAtThisLevel,numOutFlowFaces,gravity,
                           lo_bc,hi_bc);
    }
    else
    {
        Vel_in->FillBoundary();

	for (int iface = 0; iface < numOutFlowFaces; iface++) 
	    (*Vel_in).copyTo(dudt[0][iface],0,0,BL_SPACEDIM,1);

	if (have_divu) {
            for (int iface = 0; iface < numOutFlowFaces; iface++) 
                (*Divu_in).copyTo(dsdt[iface],0,0,1,1);
	} else {
            for (int iface = 0; iface < numOutFlowFaces; iface++) 
                dsdt[iface].setVal(0);
	}

        const int*      lo_bc = phys_bc->lo();
        const int*      hi_bc = phys_bc->hi();
        projBC.computeBC(dudt, dsdt, rho, phi_fine_strip,
                         parent->Geom(lev),
                         outFacesAtThisLevel,
                         numOutFlowFaces, lo_bc, hi_bc, gravity);
    }

    for (int i = 0; i < 2*BL_SPACEDIM; i++)
    {
        rho[i].clear();
        dsdt[i].clear();
        dudt[0][i].clear();
    }

    for ( int iface = 0; iface < numOutFlowFaces; iface++)
    {
        BoxArray phi_fine_strip_ba(phi_fine_strip[iface].box());
        DistributionMapping dm {phi_fine_strip_ba};
        MultiFab phi_fine_strip_mf(phi_fine_strip_ba,dm,1,0);

        for (MFIter mfi(phi_fine_strip_mf); mfi.isValid(); ++mfi) {
            phi_fine_strip_mf[mfi].copy(phi_fine_strip[iface]);
        }

        phi[lev]->copy(phi_fine_strip_mf);
    }

    if (lev > c_lev) 
    {
        putDown(phi, phi_fine_strip, c_lev, lev, outFacesAtThisLevel,
                numOutFlowFaces, ncStripWidth);
    }
}


//
// Given vel, rhs & sig, this solves Div (sig * Grad phi) = Div vel + rhs.
// On return, vel becomes vel  - sig * Grad phi.
//
void Projection::doMLMGNodalProjection (int c_lev, int nlevel, 
                                        const Vector<MultiFab*>& vel, 
                                        const Vector<MultiFab*>& phi,
                                        const Vector<MultiFab*>& sig,
                                        const Vector<MultiFab*>& rhs_cc,
                                        const Vector<MultiFab*>& rhnd,
                                        Real rel_tol, Real abs_tol,
                                        MultiFab* sync_resid_crse,
                                        MultiFab* sync_resid_fine,
                                        bool doing_initial_velproj,
                                        bool doing_initial_vortproj)
{
    BL_PROFILE("Projection:::doMLMGNodalProjection()");

    int f_lev = c_lev + nlevel - 1;

    Vector<MultiFab> vel_test(nlevel);
    Vector<MultiFab> phi_test(nlevel);
    if (test_mlmg_solver) {
        for (int i = 0; i < nlevel; ++i) {
            vel_test[i].define(vel[c_lev+i]->boxArray(), vel[c_lev+i]->DistributionMap(),
                              vel[c_lev+i]->nComp(), vel[c_lev+i]->nGrow());
            MultiFab::Copy(vel_test[i], *vel[c_lev+i], 0, 0,
                           vel[c_lev+i]->nComp(), vel[c_lev+i]->nGrow());

            phi_test[i].define(phi[c_lev+i]->boxArray(), phi[c_lev+i]->DistributionMap(),
                              phi[c_lev+i]->nComp(), phi[c_lev+i]->nGrow());

//            MultiFab::Copy(phi_test[i], *phi[c_lev+i], 0, 0,
//                           phi[c_lev+i]->nComp(), phi[c_lev+i]->nGrow());
        }
    }
    
    BL_ASSERT(vel[c_lev]->nGrow() == 1);
    BL_ASSERT(vel[f_lev]->nGrow() == 1);
    BL_ASSERT(phi[c_lev]->nGrow() == 1);
    BL_ASSERT(phi[f_lev]->nGrow() == 1);
    BL_ASSERT(sig[c_lev]->nGrow() == 1);
    BL_ASSERT(sig[f_lev]->nGrow() == 1);
    
    BL_ASSERT(sig[c_lev]->nComp() == 1);
    BL_ASSERT(sig[f_lev]->nComp() == 1);
    
    if (sync_resid_crse != 0) { 
        BL_ASSERT(nlevel == 1);
        BL_ASSERT(c_lev < parent->finestLevel());
    }
    
    if (sync_resid_fine != 0) { 
        BL_ASSERT((nlevel == 1 || nlevel == 2));
        BL_ASSERT(c_lev > 0);
    }
    
    if (rhs_cc[c_lev]) {
        AMREX_ALWAYS_ASSERT(rhs_cc[c_lev]->boxArray().ixType().cellCentered());
        BL_ASSERT(rhs_cc[c_lev]->nGrow() == 1);
        BL_ASSERT(rhs_cc[f_lev]->nGrow() == 1);
    }

    set_boundary_velocity(c_lev, nlevel, vel, doing_initial_velproj, true);

    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (Geometry::isPeriodic(idim))
        {
            mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Periodic;
        }
        else if (doing_initial_vortproj)
        {
            mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Dirichlet;
        }
        else
        {
            if (phys_bc->lo(idim) == Outflow) {
                mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            } else if (phys_bc->lo(idim) == Inflow) {
                mlmg_lobc[idim] = LinOpBCType::inflow;
            } else {
                mlmg_lobc[idim] = LinOpBCType::Neumann;
            }

            if (phys_bc->hi(idim) == Outflow) {
                mlmg_hibc[idim] = LinOpBCType::Dirichlet;
            } else if (phys_bc->hi(idim) == Inflow) {
                mlmg_hibc[idim] = LinOpBCType::inflow;
            } else {
                mlmg_hibc[idim] = LinOpBCType::Neumann;
            }
        }
    }

    Vector<Geometry> mg_geom(nlevel);
    for (int lev = 0; lev < nlevel; lev++) {
        mg_geom[lev] = parent->Geom(lev+c_lev);
    }  

    Vector<BoxArray> mg_grids(nlevel);
    for (int lev = 0; lev < nlevel; lev++) {
        mg_grids[lev] = parent->boxArray(lev+c_lev);
    }

    Vector<DistributionMapping> mg_dmap(nlevel);
    for (int lev=0; lev < nlevel; lev++ ) {
        mg_dmap[lev] = LevelData[lev+c_lev]->get_new_data(State_Type).DistributionMap();
    }

    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setMetricTerm(false);

    MLNodeLaplacian mlndlap(mg_geom, mg_grids, mg_dmap, info);
#if (AMREX_SPACEDIM == 2)
    if (rz_correction) {
        mlndlap.setRZCorrection(Geometry::IsRZ());
    }
#endif
    mlndlap.setGaussSeidel(use_gauss_seidel);
    mlndlap.setHarmonicAverage(use_harmonic_average);

    mlndlap.setDomainBC(mlmg_lobc, mlmg_hibc);
  
    for (int ilev = 0; ilev < nlevel; ++ilev) {
        mlndlap.setSigma(ilev, *sig[c_lev+ilev]);
    }

    Vector<MultiFab> rhs(nlevel);
    for (int ilev = 0; ilev < nlevel; ++ilev)
    {
        const auto& ba = amrex::convert(mg_grids[ilev], IntVect::TheNodeVector());
        rhs[ilev].define(ba, mg_dmap[ilev], 1, 0);
    }

    Vector<MultiFab*> vel_rebase{vel.begin()+c_lev, vel.begin()+c_lev+nlevel};
    Vector<const MultiFab*> rhnd_rebase{rhnd.begin(), rhnd.end()};
    rhnd_rebase.resize(nlevel,nullptr);
    Vector<MultiFab*> rhcc_rebase{rhs_cc.begin()+c_lev, rhs_cc.begin()+c_lev+nlevel};
    mlndlap.compRHS(amrex::GetVecOfPtrs(rhs), vel_rebase, rhnd_rebase, rhcc_rebase);

    MLMG mlmg(mlndlap);
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setVerbose(P_code);

    Vector<MultiFab*> phi_rebase(phi.begin()+c_lev, phi.begin()+c_lev+nlevel);
    Real mlmg_err = mlmg.solve(phi_rebase, amrex::GetVecOfConstPtrs(rhs), rel_tol, abs_tol);

    if (test_mlmg_solver) {
        Vector<MultiFab*> vel_ptmp(f_lev+1);
        Vector<MultiFab*> phi_ptmp(f_lev+1);
        for (int i = 0; i < nlevel; ++i) {
            vel_ptmp[c_lev+i] = &vel_test[i];
            phi_ptmp[c_lev+i] = &phi_test[i];

            MultiFab::Copy(phi_test[i], *phi[c_lev+i], 0, 0,
                           phi[c_lev+i]->nComp(), phi[c_lev+i]->nGrow());
        }
        doNodalProjection(c_lev, nlevel, vel_ptmp, phi_ptmp, sig, rhs_cc, rhnd, rel_tol, mlmg_err*1.01,
                          sync_resid_crse, sync_resid_fine, doing_initial_velproj);

        int niters = amrex_f90mg_get_niters();
        if (niters > 0) {
            amrex::Print() << "WARNING!!! F90MG # iters: " << niters << "\n";
        }
    }

    if (sync_resid_fine != 0 or sync_resid_crse != 0)
    {
        set_boundary_velocity(c_lev, 1, vel, doing_initial_velproj, false);
    }

    if (sync_resid_fine != 0)
    {
        MultiFab resid_save;
        Real rmin, rmax;
        if (test_mlmg_solver)
        {
            resid_save.define(sync_resid_fine->boxArray(),
                              sync_resid_fine->DistributionMap(), 1, 0);
            MultiFab::Copy(resid_save, *sync_resid_fine, 0, 0, 1, 0);
            rmin = resid_save.min(0);
            rmax = resid_save.max(0);
        }

        mlndlap.compSyncResidualFine(*sync_resid_fine, *phi[c_lev], *vel[c_lev], rhs_cc[c_lev]);

        if (test_mlmg_solver)
        {
            MultiFab::Subtract(resid_save, *sync_resid_fine, 0, 0, 1, 0);
            Real dmin = resid_save.min(0);
            Real dmax = resid_save.max(0);
            amrex::Print() << "TEST MLMG: fine resid diff "
                           << "(" << dmin <<", " << dmax << ") / ("
                           << rmin << ", " << rmax << ")\n"
                           << "                           "
                           << "(" << dmin/(rmin+1.e-50) << ", "
                           << dmax/(rmax+1.e-50) << ")\n";
        }
    }

    if (sync_resid_crse != 0) {  // only level solve will come to here

        MultiFab resid_save;
        Real rmin, rmax;
        if (test_mlmg_solver)
        {
            resid_save.define(sync_resid_crse->boxArray(),
                              sync_resid_crse->DistributionMap(), 1, 0);
            MultiFab::Copy(resid_save, *sync_resid_crse, 0, 0, 1, 0);
            rmin = resid_save.min(0);
            rmax = resid_save.max(0);
        }

        const BoxArray& fineGrids = parent->boxArray(c_lev+1);
        const IntVect& ref_ratio = parent->refRatio(c_lev);
        mlndlap.compSyncResidualCoarse(*sync_resid_crse, *phi[c_lev], *vel[c_lev], rhs_cc[c_lev],
                                       fineGrids, ref_ratio);

        if (test_mlmg_solver)
        {
            MultiFab::Subtract(resid_save, *sync_resid_crse, 0, 0, 1, 0);
            Real dmin = resid_save.min(0);
            Real dmax = resid_save.max(0);
            amrex::Print() << "TEST MLMG: crse resid diff "
                           << "(" << dmin <<", " << dmax << ") / ("
                           << rmin << ", " << rmax << ")\n"
                           << "                           "
                           << "(" << dmin/(rmin+1.e-50) << ", "
                           << dmax/(rmax+1.e-50) << ")\n";
        }
    }

    mlndlap.updateVelocity(vel_rebase, amrex::GetVecOfConstPtrs(phi_rebase));
}


//
// Given vel, rhs & sig, this solves Div (sig * Grad phi) = Div vel + rhs.
// On return, vel becomes vel  - sig * Grad phi.
//
void Projection::doNodalProjection (int c_lev, int nlevel, 
                                    const Vector<MultiFab*>& vel, 
                                    const Vector<MultiFab*>& phi,
                                    const Vector<MultiFab*>& sig,
				    const Vector<MultiFab*>& rhs_cc, 
                                    const Vector<MultiFab*>& rhnd, 
				    Real rel_tol, Real abs_tol,
				    MultiFab* sync_resid_crse,
				    MultiFab* sync_resid_fine,
                                    bool doing_initial_velproj) 
{
  BL_PROFILE("Projection:::doNodalProjection()");

  int f_lev = c_lev + nlevel - 1;

  BL_ASSERT(vel[c_lev]->nGrow() == 1);
  BL_ASSERT(vel[f_lev]->nGrow() == 1);
  BL_ASSERT(phi[c_lev]->nGrow() == 1);
  BL_ASSERT(phi[f_lev]->nGrow() == 1);
  BL_ASSERT(sig[c_lev]->nGrow() == 1);
  BL_ASSERT(sig[f_lev]->nGrow() == 1);

  BL_ASSERT(sig[c_lev]->nComp() == 1);
  BL_ASSERT(sig[f_lev]->nComp() == 1);

  if (sync_resid_crse != 0) { 
    BL_ASSERT(nlevel == 1);
    BL_ASSERT(c_lev < parent->finestLevel());
  }

  if (sync_resid_fine != 0) { 
    BL_ASSERT((nlevel == 1 || nlevel == 2));
    BL_ASSERT(c_lev > 0);
  }

  if (rhs_cc[c_lev]) {
    if (rhs_cc[c_lev]->box(0).type() == IntVect::TheNodeVector()) {
      amrex::Abort("Projection::doNodalProjection: rhs_cc cannot be nodal type");
    }
    BL_ASSERT(rhs_cc[c_lev]->nGrow() == 1);
    BL_ASSERT(rhs_cc[f_lev]->nGrow() == 1);
  }

  Vector<std::unique_ptr<MultiFab> > vold(maxlev);
  if (sync_resid_fine !=0 || sync_resid_crse != 0) {
    vold[c_lev].reset(new MultiFab(parent->boxArray(c_lev), 
                                   parent->DistributionMap(c_lev),
                                   BL_SPACEDIM, 1));
    MultiFab::Copy(*vold[c_lev], *vel[c_lev], 0, 0, BL_SPACEDIM, 1);

    set_boundary_velocity(c_lev, 1,
                          amrex::GetVecOfPtrs(vold),
                          doing_initial_velproj, false);
  }

  set_boundary_velocity(c_lev, nlevel, vel, doing_initial_velproj, true);

  int lo_inflow[3] = {0};
  int hi_inflow[3] = {0};

  const int* lo_bc = phys_bc->lo();
  const int* hi_bc = phys_bc->hi();
  
  for (int idir=0; idir<BL_SPACEDIM; idir++) {
    if (lo_bc[idir] == Inflow) {
      lo_inflow[idir] = 1;
    }
    if (hi_bc[idir] == Inflow) {
      hi_inflow[idir] = 1;
    }
  }

  Vector<Geometry> mg_geom(nlevel);
  for (int lev = 0; lev < nlevel; lev++) {
    mg_geom[lev] = parent->Geom(lev+c_lev);
  }  

  int mg_bc[2*BL_SPACEDIM];
  for ( int i = 0; i < BL_SPACEDIM; ++i ) {
    if ( mg_geom[0].isPeriodic(i) ) {
      mg_bc[i*2 + 0] = 0;
      mg_bc[i*2 + 1] = 0;
    }
    else {
      mg_bc[i*2 + 0] = phys_bc->lo(i)==Outflow? MGT_BC_DIR : MGT_BC_NEU;
      mg_bc[i*2 + 1] = phys_bc->hi(i)==Outflow? MGT_BC_DIR : MGT_BC_NEU;
    }
  }

  Vector<BoxArray> mg_grids(nlevel);
  for (int lev = 0; lev < nlevel; lev++) {
    mg_grids[lev] = parent->boxArray(lev+c_lev);
  }

  Vector<DistributionMapping> dmap(nlevel);
  for (int lev=0; lev < nlevel; lev++ ) {
    dmap[lev] = LevelData[lev+c_lev]->get_new_data(State_Type).DistributionMap();
  }

  bool nodal = true;

  bool have_rhcc;
  if (rhs_cc[c_lev] == nullptr) {
    have_rhcc = false;
  }
  else {
    have_rhcc = false;
    for (int lev=c_lev; lev<=f_lev; lev++) {
      if (rhs_cc[lev]->norm0() != 0.0) {
	have_rhcc = true;
	break;
      }
    }
  }

  MGT_Solver mgt_solver(mg_geom, mg_bc, mg_grids, dmap, nodal, hg_stencil, have_rhcc,
                        0, 1, P_code);

  mgt_solver.set_nodal_coefficients({sig.begin()+c_lev, sig.end()});

  mgt_solver.nodal_project({phi.begin()+c_lev, phi.end()}, 
                           {vel.begin()+c_lev, vel.end()},
                           {rhs_cc.begin()+c_lev, rhs_cc.end()}, 
                           rhnd, rel_tol, abs_tol, &lo_inflow[0], &hi_inflow[0]);

  // Must fill sync_resid_fine before sync_resid_crse because of the side effecs in the calls.

  if (sync_resid_fine != 0) {
    const BoxArray& levelGrids = mg_grids[0];
    const DistributionMapping& levelDmap = dmap[0];
    const Geometry& levelGeom = mg_geom[0];

    MultiFab msk(levelGrids, levelDmap, 1, 1); 

    mask_grids(msk, levelGeom);

    sync_resid_fine->setVal(0.0, sync_resid_fine->nGrow());

    int isCoarse = 0;
    mgt_solver.fill_sync_resid(*sync_resid_fine, msk, *vold[c_lev], isCoarse);
  }

  if (sync_resid_crse != 0) {  // only level solve will come to here
    const BoxArray& fineGrids = parent->boxArray(c_lev+1);
    const BoxArray& levelGrids = mg_grids[0];
    const DistributionMapping& levelDmap = dmap[0];
    const Geometry& levelGeom = mg_geom[0];
    IntVect ref_ratio = parent->refRatio(c_lev);

    MultiFab msk(levelGrids, levelDmap, 1, 1); 

    mask_grids(msk, levelGrids, levelGeom, fineGrids, ref_ratio);

    sync_resid_crse->setVal(0.0, sync_resid_crse->nGrow());

    int isCoarse = 1;
    mgt_solver.fill_sync_resid(*sync_resid_crse, msk, *vold[c_lev], isCoarse);
  }

  if (verbose >= 1)
    MGT_Solver::FlushFortranOutput();
}


void
Projection::mask_grids (MultiFab& msk, const BoxArray& grids, const Geometry& geom,
                        const BoxArray& fineGrids, const IntVect& ref_ratio)
{
  BL_PROFILE("Projection::mask_grids(1)");

  BoxArray localfine = fineGrids;
  localfine.coarsen(ref_ratio);

  msk.setBndry(BogusValue);

  const int* lo_bc = phys_bc->lo();
  const int* hi_bc = phys_bc->hi();
  const Box& domainBox = geom.Domain();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(msk); mfi.isValid(); ++mfi) {
    int i = mfi.index();

    FArrayBox& msk_fab = msk[mfi];

    const Box& reg  = grids[i]; 
    msk_fab.setVal(1.0, reg, 0); 

    for (int idir=0; idir<BL_SPACEDIM; idir++) {
      if (lo_bc[idir] == Inflow) {
	if (reg.smallEnd(idir) == domainBox.smallEnd(idir)) {
	  Box bx = amrex::adjCellLo(reg, idir);
	  msk_fab.setVal(1.0, bx, 0);
	}
      }
      if (hi_bc[idir] == Inflow) {
	if (reg.bigEnd(idir) == domainBox.bigEnd(idir)) {
	  Box bx = amrex::adjCellHi(reg, idir);
	  msk_fab.setVal(1.0, bx, 0);
	}
      }
    }

    std::vector< std::pair<int,Box> > isects;
    localfine.intersections(reg,isects);

    for (int ii = 0; ii < isects.size(); ii++) {
      const Box& fbox = isects[ii].second;
      msk_fab.setVal(0.0, fbox, 0);

      for (int idir=0; idir<BL_SPACEDIM; idir++) {
	if (lo_bc[idir] == Inflow) {
	  if (fbox.smallEnd(idir) == domainBox.smallEnd(idir)) {
	    Box bx = amrex::adjCellLo(fbox, idir);
	    msk_fab.setVal(0.0, bx, 0);
	  }
	}
	if (hi_bc[idir] == Inflow) {
	  if (fbox.bigEnd(idir) == domainBox.bigEnd(idir)) {
	    Box bx = amrex::adjCellHi(fbox, idir);
	    msk_fab.setVal(0.0, bx, 0);
	  }
	}
      }      
    }
  }

  msk.FillBoundary();
  msk.EnforcePeriodicity(geom.periodicity());
}

void Projection::mask_grids (MultiFab& msk, const Geometry& geom)
{
    BL_PROFILE("Projection::mask_grids(2)");
    
    Real fineghost = 1.0;
    Real interior  = 1.0;
    Real crseghost = 0.0;
    Real physbnd   = BogusValue;
    
    const Box& domainBox = geom.Domain();
    const Periodicity& period = geom.periodicity();
    
    msk.BuildMask(domainBox, period, fineghost, crseghost, physbnd, interior);
    
    const int* lo_bc = phys_bc->lo();
    const int* hi_bc = phys_bc->hi();
    
    bool has_inflow = false;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	if (lo_bc[i] == Inflow || hi_bc[i] == Inflow) {
	    has_inflow = true;
	    break;
	}
    }

    if (has_inflow) {
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(msk); mfi.isValid(); ++mfi)
	{
	    FArrayBox& msk_fab = msk[mfi];
	    const Box& regBox = mfi.validbox();
	    
	    for (int idir=0; idir<BL_SPACEDIM; idir++) {
		if (lo_bc[idir] == Inflow) {
		    if (regBox.smallEnd(idir) == domainBox.smallEnd(idir)) {
			const Box& bx = amrex::adjCellLo(regBox, idir);
			msk_fab.setVal(1.0, bx, 0);
		    }
		}
		if (hi_bc[idir] == Inflow) {
		    if (regBox.bigEnd(idir) == domainBox.bigEnd(idir)) {
			const Box& bx = amrex::adjCellHi(regBox, idir);
			msk_fab.setVal(1.0, bx, 0);
		    }
		}
	    }
	}
	
	msk.EnforcePeriodicity(period);
    }
}

// Set velocity in ghost cells to zero except for inflow
void Projection::set_boundary_velocity(int c_lev, int nlevel, const Vector<MultiFab*>& vel, 
                                       bool doing_initial_velproj, bool inflowCorner)
{
  const int* lo_bc = phys_bc->lo();
  const int* hi_bc = phys_bc->hi();

  // 1) At non-inflow faces, the normal component of velocity will be completely zero'd 
  // 2) If a face is an inflow face, then
  //     i) if inflowCorner = false then the normal velocity at corners -- even periodic corners -- 
  //                                just outside inflow faces will be zero'd
  //    ii) if inflowCorner =  true then the normal velocity at corners just outside inflow faces 
  //                                will be zero'd outside of Neumann boundaries 
  //                                (slipWall, noSlipWall, Symmetry) 
  //                                will retain non-zero values at periodic corners

  for (int lev=c_lev; lev < c_lev+nlevel; lev++) {
    const BoxArray& grids = parent->boxArray(lev);
    const Box& domainBox = parent->Geom(lev).Domain();

    const Geometry& geom = parent->Geom(lev);

    for (int idir=0; idir<BL_SPACEDIM; idir++) {

      if (lo_bc[idir] != Inflow && hi_bc[idir] != Inflow) {
	vel[lev]->setBndry(0.0, Xvel+idir, 1);
      }
      else {

	for (MFIter mfi(*vel[lev]); mfi.isValid(); ++mfi) {
	  int i = mfi.index();

	  FArrayBox& v_fab = (*vel[lev])[mfi];

	  const Box& reg = grids[i];
	  const Box& bxg1 = amrex::grow(reg, 1);

	  BoxList bxlist(reg);

	  if (lo_bc[idir] == Inflow && reg.smallEnd(idir) == domainBox.smallEnd(idir)) {
	    Box bx;                // bx is the region we *protect* from zero'ing

	    if (inflowCorner) {

              bx = amrex::adjCellLo(reg, idir);

              for (int odir = 0; odir < BL_SPACEDIM; odir++) {
                 if (odir != idir)
                 {
                    if (geom.isPeriodic(odir)) bx.grow(odir,1);
                    if (reg.bigEnd  (odir) != domainBox.bigEnd  (odir) ) bx.growHi(odir,1);
                    if (reg.smallEnd(odir) != domainBox.smallEnd(odir) ) bx.growLo(odir,1);
                 }
              }

	      // This is the old code
              // bx = amrex::adjCellLo(bxg1, idir);
	      // bx.shift(idir, +1);

	    } else {
	      bx = amrex::adjCellLo(reg, idir);
	    }
	    bxlist.push_back(bx);
	  }

	  if (hi_bc[idir] == Inflow && reg.bigEnd(idir) == domainBox.bigEnd(idir)) {
	    Box bx;                // bx is the region we *protect* from zero'ing

	    if (inflowCorner) {

	      bx = amrex::adjCellHi(reg, idir);

              for (int odir = 0; odir < BL_SPACEDIM; odir++) {
                 if (odir != idir)
                 {
                    if (geom.isPeriodic(odir)) bx.grow(odir,1);
                    if (reg.bigEnd  (odir) != domainBox.bigEnd  (odir) ) bx.growHi(odir,1);
                    if (reg.smallEnd(odir) != domainBox.smallEnd(odir) ) bx.growLo(odir,1);
                 }
              }

	      // This is the old code
              // bx = amrex::adjCellHi(bxg1, idir);
              // bx.shift(idir, -1);

	    } else {
	      bx = amrex::adjCellHi(reg, idir);
	    }

	    bxlist.push_back(bx);
	  }

	  BoxList bxlist2 = amrex::complementIn(bxg1, bxlist); 
 
	  for (BoxList::iterator it=bxlist2.begin(); it != bxlist2.end(); ++it) {
            Box ovlp = *it & v_fab.box();
            if (ovlp.ok()) {
              v_fab.setVal(0.0, ovlp, Xvel+idir, 1);
            }
	  }
	}
      }
    }

  }
}
