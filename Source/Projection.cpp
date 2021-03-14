

#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <NavierStokesBase.H>
#include <NS_BC.H>
#include <AMReX_BLProfiler.H>
#include <Projection.H>
#include <PROJECTION_F.H>
#include <OutFlowBC.H>
#include <NSB_K.H>

#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>

#include <AMReX_NodalProjector.H>


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

int  Projection::P_code              = -1;
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

    bool benchmarking = false;

    bool rz_correction = false;
    bool agglomeration = true;
    bool consolidation = true;
    int max_fmg_iter = 0;
    bool use_gauss_seidel = true;
    bool use_harmonic_average = false;
}


void
Projection::Initialize ()
{
    if (initialized) return;

    ParmParse pp("proj");

    pp.query("v",                   verbose);
    pp.query("proj_tol",            proj_tol);
    pp.query("sync_tol",            sync_tol);
    pp.query("proj_abs_tol",        proj_abs_tol);
    pp.query("benchmarking",        benchmarking);
    pp.query("add_vort_proj",       add_vort_proj);
    pp.query("do_outflow_bcs",      do_outflow_bcs);
    pp.query("rho_wgt_vel_proj",    rho_wgt_vel_proj);
    pp.query("divu_minus_s_factor", divu_minus_s_factor);
    pp.query("make_sync_solvable",  make_sync_solvable);

    pp.query("rz_correction",       rz_correction);
    pp.query("agglomeration",       agglomeration);
    pp.query("consolidation",       consolidation);
    pp.query("max_fmg_iter",        max_fmg_iter);
    pp.query("use_gauss_seidel",    use_gauss_seidel);
    pp.query("use_harmonic_average", use_harmonic_average);

    pp.query("proj_2",              proj_2);
    if (!proj_2)
      amrex::Abort("Must use proj_2==1 due to new gravity and outflow stuff. proj_2!=1 no longer supported.\n");

    pp.query("Pcode",               P_code);
    if (P_code >=0 )
      amrex::Abort("proj.Pcode is no more. Use nodal_proj.verbose.\n");

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
    phys_bc(_phys_bc),
    do_sync_proj(_do_sync_proj)
{

    AMREX_ASSERT ( parent->finestLevel()+1 <= maxlev );

    Initialize();

    if (verbose) amrex::Print() << "Creating projector\n";

#ifdef AMREX_USE_EB
    // size the EB factory array
    ebfactory.resize(maxlev + 1);
#endif
}

Projection::~Projection ()
{
  if (verbose) amrex::Print() << "Deleting projector\n";
}

//
// Install a level of the projection.
//

void
Projection::install_level (int                     level,
                           AmrLevel*               level_data,
                           Vector< Vector<Real> >* _radius)
{
    if (verbose) amrex::Print() << "Installing projector level " << level << '\n';

    int finest_level = parent->finestLevel();

    if (level > LevelData.size() - 1)
    {
        LevelData.resize(finest_level+1);
        radius.resize(finest_level+1);
    }

    LevelData[level] = level_data;
    radius[level] = _radius;

#ifdef AMREX_USE_EB
    const auto& _ebfactory =
      dynamic_cast<EBFArrayBoxFactory const&>(LevelData[level]->Factory());
    ebfactory[level] = &_ebfactory;
#endif
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

    AMREX_ASSERT(rho_half.nGrow() >= 1);
    AMREX_ASSERT(U_new.nGrow() >= 1);

    if (verbose) {
      amrex::Print() << "... Projection::level_project() at level " << level << '\n';
    }

    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real strt_time = ParallelDescriptor::second();

    //
    // old time velocity has bndry values already
    // must gen valid bndry data for new time velocity.
    // must fill bndry cells in pressure with computable values
    // even though they are not used in calculation.
    //
    U_old.setBndry(BogusValue,Xvel,AMREX_SPACEDIM);
    U_new.setBndry(BogusValue,Xvel,AMREX_SPACEDIM);
    P_old.setBndry(BogusValue);
    P_new.setBndry(BogusValue);

    MultiFab& S_old = LevelData[level]->get_old_data(State_Type);
    MultiFab& S_new = LevelData[level]->get_new_data(State_Type);

    Real prev_time = LevelData[level]->get_state_data(State_Type).prevTime();
    Real curr_time = LevelData[level]->get_state_data(State_Type).curTime();

    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        LevelData[level]->setPhysBoundaryValues(S_old[mfi],State_Type,prev_time,
                                                Xvel,Xvel,AMREX_SPACEDIM);
        LevelData[level]->setPhysBoundaryValues(S_new[mfi],State_Type,curr_time,
                                                Xvel,Xvel,AMREX_SPACEDIM);
    }

    const BoxArray& P_grids = P_old.boxArray();
    const DistributionMapping& P_dmap = P_old.DistributionMap();

    NavierStokesBase* ns = dynamic_cast<NavierStokesBase*>(&parent->getLevel(level));
    AMREX_ASSERT(!(ns==0));

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
    }

    const int nGrow = (level == 0  ?  0  :  -1);
    //
    // MultiFab::setVal() won't work here bc nGrow could be <0
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(P_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.growntilebox(nGrow);
       auto const& pnew = P_new.array(mfi);  
       amrex::ParallelFor(bx, [pnew]
       AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
          pnew(i,j,k) = 0.0;
       });
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
    }

    const Real dt_inv = 1./dt;
    U_new.mult(dt_inv,0,AMREX_SPACEDIM,1);
    if (have_divu)
      divusource->mult(dt_inv,0,1,divusource->nGrow());

    MultiFab& Gp = ns->get_old_data(Gradp_Type);
    
#ifndef NDEBUG
#ifdef AMREX_USE_EB
    // fixme - deal with case where covered cells are set to zero
    //   there's probably a better way to handle this..
    EB_set_covered(rho_half,0,1,1,1.2345e40);
#endif
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rho_half,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       // No ghost cells needed here. Outflow BCs extrapolate from interior.
       // Velocity ghost cells are filled in doMLMGNodalProjection().
       const Box& bx = mfi.tilebox();
       const auto& rho_h = rho_half.array(mfi);
       const auto& gradp = Gp.array(mfi);
       const auto& u_new = U_new.array(mfi);
       amrex::ParallelFor(bx, AMREX_SPACEDIM, [rho_h,gradp,u_new]
       AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
       {
	   u_new(i,j,k,n) += gradp(i,j,k,n) / rho_h(i,j,k);
       });
    }

    //
    // Outflow uses appropriately constructed "U_new" and "divusource" to
    //   compute BC for phi, so make sure this call comes after those are set,
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
       U_new.FillBoundary(0, AMREX_SPACEDIM, geom.periodicity());
       rho_half.FillBoundary(0, 1, geom.periodicity());
    }
    //
    // Add the contribution from the un-projected V to syncregisters.
    //
    int is_rz = (geom.IsRZ() ? 1 : 0);

    Vector<MultiFab*> vel(maxlev, nullptr);
    Vector<MultiFab*> phi(maxlev, nullptr);
    Vector<MultiFab*> sig(maxlev, nullptr);

    vel[level] = &U_new;
    phi[level] = &P_new;

    AMREX_ASSERT( 1 == rho_half.nGrow());
    sig[level] = &rho_half;

    //
    // Project
    //
    std::unique_ptr<MultiFab> sync_resid_crse, sync_resid_fine;
#ifdef AMREX_USE_EB
    const auto& factory =
      dynamic_cast<EBFArrayBoxFactory const&>(LevelData[level]->Factory());
#else
    const auto& factory = LevelData[level]->Factory();
#endif

    if (level < parent->finestLevel()) {
      sync_resid_crse.reset(new MultiFab(P_grids,P_dmap,1,1, MFInfo(), factory));
      sync_resid_crse->setVal(0.);
    }

    if (level > 0 && iteration == crse_dt_ratio)
    {
      const int ngrow = parent->MaxRefRatio(level-1) - 1;
      sync_resid_fine.reset(new MultiFab(P_grids,P_dmap,1,ngrow, MFInfo(), factory));
      sync_resid_fine->setVal(0.);
    }

    Vector<MultiFab*> rhcc;
    if (have_divu)
    {
       rhcc.resize(maxlev);
       if (is_rz == 1) {
          radMultScal(level,*divusource);
       }
       const int nghost = 0;
       divusource->mult(-1.0,0,1,nghost);
       rhcc[level] = divusource.get();
    }

    bool increment_gp = false;
    doMLMGNodalProjection(level, 1, vel, phi, sig, rhcc, {}, proj_tol,
			  proj_abs_tol, increment_gp,
                          sync_resid_crse.get(), sync_resid_fine.get());

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
       if (level > 0 && iteration == crse_dt_ratio)
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
    // un = dt*un
    //
    U_new.mult(dt,0,AMREX_SPACEDIM,1);

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
                           const Geometry& crse_geom)
{
    BL_PROFILE("Projection::MLsyncProject()");

    if (verbose) {
      amrex::Print() << "Projection::MLsyncProject(): levels = " << c_lev << ", " << c_lev+1 << '\n';
    }

    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real strt_time = ParallelDescriptor::second();

    //
    // Set up memory.
    //
    Vector<std::unique_ptr<MultiFab> > phi(maxlev);

    const BoxArray& Pgrids_crse = pres_crse.boxArray();
    const BoxArray& Pgrids_fine = pres_fine.boxArray();
    const DistributionMapping& Pdmap_crse = pres_crse.DistributionMap();
    const DistributionMapping& Pdmap_fine = pres_fine.DistributionMap();

    phi[c_lev].reset(new MultiFab(Pgrids_crse,Pdmap_crse,1,1,MFInfo(),LevelData[c_lev]->Factory()));
    phi[c_lev]->setVal(0);

    phi[c_lev+1].reset(new MultiFab(Pgrids_fine,Pdmap_fine,1,1,MFInfo(),LevelData[c_lev+1]->Factory()));
    phi[c_lev+1]->setVal(0);

    //
    // Set up crse RHS
    //
    MultiFab rhnd(Pgrids_crse,Pdmap_crse,1,0,MFInfo(),LevelData[c_lev]->Factory());
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

    if (crse_geom.IsRZ()) {
       radMultScal(c_lev  ,cc_rhs_crse);
       radMultScal(c_lev+1,cc_rhs_fine);
    }

    Vector<MultiFab*> vel (maxlev, nullptr);
    Vector<MultiFab*> sig (maxlev, nullptr);
    Vector<MultiFab*> rhcc(maxlev, nullptr);
    Vector<MultiFab*> rhnd_vec(maxlev, nullptr);

    vel[c_lev  ] = &Vsync;
    vel[c_lev+1] = &V_corr;
    sig[c_lev  ] = &rho_crse;
    sig[c_lev+1] = &rho_fine;
    rhcc[c_lev  ] = &cc_rhs_crse;
    rhcc[c_lev+1] = &cc_rhs_fine;
    rhnd_vec[c_lev] = &rhnd;

    NavierStokesBase* ns = dynamic_cast<NavierStokesBase*>(LevelData[c_lev]);
    ns->average_down(*vel[c_lev+1],*vel[c_lev],0,AMREX_SPACEDIM);

    ns->average_down(*sig[c_lev+1],*sig[c_lev],0,sig[c_lev]->nComp());

    MultiFab* sync_resid_crse = 0;
    std::unique_ptr<MultiFab> sync_resid_fine;

    if (c_lev > 0 &&  crse_iteration == crse_dt_ratio)
    {
        int ngrow = parent->MaxRefRatio(c_lev-1) - 1;
        sync_resid_fine.reset(new MultiFab(Pgrids_crse,Pdmap_crse,1,ngrow));
        sync_resid_fine->setVal(0.);
    }

    bool increment_gp = true;
    doMLMGNodalProjection(c_lev, 2, vel,
                          amrex::GetVecOfPtrs(phi),
                          sig, rhcc, rhnd_vec, sync_tol, proj_abs_tol, increment_gp,
                          sync_resid_crse, sync_resid_fine.get());

    //
    // If this sync project is not at levels 0-1 then we need to account for
    // the changes made here in the level c_lev velocity in the sync registers
    // going into the level (c_lev-1) sync project.  Note that this must be
    // done before rho_half is scaled back.
    //
    if (c_lev > 0 && crse_iteration == crse_dt_ratio)
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
    // Only update the most recent pressure.
    //
    AddPhi(pres_crse, *phi[c_lev]);
    AddPhi(pres_fine, *phi[c_lev+1]);

    //
    // Grad(P_new) incremented in doMLMGNodalProjection
    //

    //
    // Add projected vel to new velocity.
    //
    MultiFab::Saxpy(vel_crse,dt_crse,Vsync,0,0,AMREX_SPACEDIM,1);
    MultiFab::Saxpy(vel_fine,dt_crse,V_corr,0,0,AMREX_SPACEDIM,1);

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	amrex::Print() << "Projection::MLsyncProject(): levels = " << c_lev
		       << ", " << c_lev+1 << ", time: " << run_time << '\n';
    }
}

//
// The initial velocity projection in post_init.
// this function ensures that the velocities are nondivergent
//
//
void
Projection::initialVelocityProject (int  c_lev,
                                    Real cur_divu_time,
                                    int  have_divu,
                                    int init_vel_iter )
{
    int lev;
    int f_lev = parent->finestLevel();

    if ( init_vel_iter <= 0 ){
      if ( verbose ) Print()<<"Returning from initalVelocityProject() without projecting because init_vel_iter<=0\n";
      return;
    }
    
    if (verbose)
    {
        amrex::Print() << "Projection::initialVelocityProject(): levels = " << c_lev
                       << "-" << f_lev << '\n';
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

    for (int iter = 0; iter < init_vel_iter; ++iter)
    {
        if (verbose)
        {
            amrex::Print() << std::endl
                           << "Projection::initialVelocityProject(): iteration "
                           << iter
                           << std::endl;
        }

        for (lev = c_lev; lev <= f_lev; lev++)
        {
            LevelData[lev]->get_old_data(Press_Type).setVal(0.0);
        }

	// MLNodeLaplacian does not take any ghost cells from rhcc or sig.
	// copies only valid cells and fills ghosts internally.
	// However, vel and phi are assumed to have 1 ghost cell in MLMG.
	// set outflow bcs fills vel and phi using 1 ghost cell from sig (holding rho)
	// and rhcc (holding divu)
	const int nghost =
	  (OutFlowBC::HasOutFlowBC(phys_bc) && do_outflow_bcs && have_divu) ? 1 : 0;

        for (lev = c_lev; lev <= f_lev; lev++)
        {
            vel[lev] = &(LevelData[lev]->get_new_data(State_Type));
            phi[lev] = &(LevelData[lev]->get_old_data(Press_Type));

            const BoxArray& grids  = LevelData[lev]->boxArray();
            const DistributionMapping& dmap = LevelData[lev]->DistributionMap();
            sig[lev].reset(new MultiFab(grids,dmap,1,nghost,MFInfo(),LevelData[lev]->Factory()));

            if (rho_wgt_vel_proj)
            {
	      if ( nghost > 0 ){
                LevelData[lev]->get_new_data(State_Type).setBndry(BogusValue,Density,1);

                AmrLevel& amr_level = parent->getLevel(lev);

                MultiFab& S_new = amr_level.get_new_data(State_Type);

                Real curr_time = amr_level.get_state_data(State_Type).curTime();

                for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
                {
                    amr_level.setPhysBoundaryValues(S_new[mfi],State_Type,curr_time,
                                                    Density,Density,1);
                }
	      }
	      
	      MultiFab::Copy(*sig[lev],
			     LevelData[lev]->get_new_data(State_Type),
			     Density, 0, 1, nghost);
            }
            else
            {
                sig[lev]->setVal(1.,nghost);
            }
        }

        Vector<std::unique_ptr<MultiFab> > rhcc(maxlev);

        for (lev = c_lev; lev <= f_lev; lev++)
        {
            vel[lev]->setBndry(BogusValue,Xvel,AMREX_SPACEDIM);
            //
            // Set the physical boundary values.
            //
            AmrLevel& amr_level = parent->getLevel(lev);

            MultiFab& S_new = amr_level.get_new_data(State_Type);

            Real curr_time = amr_level.get_state_data(State_Type).curTime();

            for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
            {        
                amr_level.setPhysBoundaryValues(S_new[mfi],State_Type,curr_time,Xvel,Xvel,AMREX_SPACEDIM);
            }

            if (have_divu)
            {
                int Divu_Type, Divu;
                if (!LevelData[lev]->isStateVariable("divu", Divu_Type, Divu))
		  amrex::Error("Projection::initialVelocityProject(): Divu not found");

		NavierStokesBase* ns = dynamic_cast<NavierStokesBase*>(LevelData[lev]);
		AMREX_ASSERT(!(ns == 0));
		
		rhcc[lev].reset(ns->getDivCond(nghost,cur_divu_time));
            }
        }
  
        if (OutFlowBC::HasOutFlowBC(phys_bc) && do_outflow_bcs && have_divu)
        {
            set_outflow_bcs(INITIAL_VEL,phi,vel,
                            amrex::GetVecOfPtrs(rhcc),
                            amrex::GetVecOfPtrs(sig),
                            c_lev,f_lev,have_divu);
        }

        //
        // Scale the projection variables.
        //
        for (lev = c_lev; lev <= f_lev; lev++)
        {
            scaleVar(INITIAL_VEL,sig[lev].get(),nghost,vel[lev],lev);
        }

        //
        // Project
        //
        //
                
        bool increment_gp = false;
        if (!have_divu)
        {
            doMLMGNodalProjection(c_lev, f_lev-c_lev+1, vel, phi,
                                  amrex::GetVecOfPtrs(sig),
                                  {},
                                  {},
                                  proj_tol, proj_abs_tol, increment_gp, 0, 0);
        }
        else
        {
            for (lev = c_lev; lev <= f_lev; lev++)
            {
                if (parent->Geom(0).IsRZ()) radMultScal(lev,*rhcc[lev]);
                rhcc[lev]->mult(-1.0,0,1,nghost);
            }

            doMLMGNodalProjection(c_lev, f_lev-c_lev+1, vel, phi,
                                  amrex::GetVecOfPtrs(sig),
                                  amrex::GetVecOfPtrs(rhcc),
                                  {},
                                  proj_tol, proj_abs_tol, increment_gp, 0, 0);
        }
        
        //
        // Unscale initial projection variables.
        //
        for (lev = c_lev; lev <= f_lev; lev++)
            rescaleVar(INITIAL_VEL,sig[lev].get(),nghost,vel[lev],lev);

        for (lev = c_lev; lev <= f_lev; lev++)
        {
            LevelData[lev]->get_old_data(Press_Type).setVal(0.);
            LevelData[lev]->get_new_data(Press_Type).setVal(0.);

            LevelData[lev]->get_old_data(Gradp_Type).setVal(0.);
            LevelData[lev]->get_new_data(Gradp_Type).setVal(0.);
        }

        if (verbose)
        {
            amrex::Print() << "After nodal projection:" << std::endl;
            for (lev = c_lev; lev <= f_lev; ++lev)
            {
                amrex::Print() << "  lev " << lev << ": "
#if (AMREX_SPACEDIM==3)
                               << "max(abs(u,v,w)) = "
#else
                               << "max(abs(u,v)) = "
#endif
                               << vel[lev]->norm0(0,0,false,true) << " "
                               << vel[lev]->norm0(1,0,false,true) << " "
#if (AMREX_SPACEDIM==3)
                               << vel[lev]->norm0(2,0,false,true)
#endif
                               << std::endl;
            }
        }
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
      amrex::Print() << "Projection::initialPressureProject(): levels = " << c_lev
           << "  " << f_lev << '\n';
    }

    const Real strt_time = ParallelDescriptor::second();

    Vector<MultiFab*> vel(maxlev, nullptr);
    Vector<MultiFab*> phi(maxlev, nullptr);
    Vector<std::unique_ptr<MultiFab> > sig(maxlev);

    for (lev = c_lev; lev <= f_lev; lev++)
    {
        vel[lev] = &(LevelData[lev]->get_new_data(State_Type));
        phi[lev] = &(LevelData[lev]->get_new_data(Press_Type));

        const int       nghost = 1;
        const BoxArray& grids  = LevelData[lev]->boxArray();
        const DistributionMapping& dmap = LevelData[lev]->DistributionMap();
	sig[lev].reset(new MultiFab(grids,dmap,1,nghost,MFInfo(),LevelData[lev]->Factory()));

        AmrLevel& amr_level = parent->getLevel(lev);

        MultiFab& S_new = amr_level.get_new_data(State_Type);

        S_new.setBndry(BogusValue,Density,1);

        Real curr_time = amr_level.get_state_data(State_Type).curTime();

        const Geometry& geom = parent->Geom(lev);
	// fill ghost cells... call FillBoundary (fills interior bndry)
	// first to get reasonable data in corner cells
	S_new.FillBoundary(Density,1,geom.periodicity());
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            amr_level.setPhysBoundaryValues(S_new[mfi],State_Type,curr_time,Density,Density,1);
        }

        MultiFab::Copy(*sig[lev],
                       LevelData[lev]->get_new_data(State_Type),
                       Density,
                       0,
                       1,
                       nghost);
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
	raii.push_back(std::unique_ptr<MultiFab>(new MultiFab(grids, dmap, AMREX_SPACEDIM, 1,MFInfo(),LevelData[lev]->Factory())));
        vel[lev] = raii.back().get();
        vel[lev]->setVal(0.0    , 0            , AMREX_SPACEDIM-1, 1);
        vel[lev]->setVal(gravity, AMREX_SPACEDIM-1, 1            , 1);
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
    bool increment_gp = false;
    Vector<MultiFab*> rhcc(0);
    doMLMGNodalProjection(c_lev, f_lev-c_lev+1, vel, phi,
                          amrex::GetVecOfPtrs(sig),
                          rhcc, {},
                          proj_tol, proj_abs_tol, increment_gp, 0, 0);

    //
    // Unscale initial projection variables.
    //
    for (lev = c_lev; lev <= f_lev; lev++) {
        rescaleVar(INITIAL_PRESS,sig[lev].get(),1,vel[lev],lev);
    }

    for (lev = c_lev; lev <= f_lev; lev++) {
        //
        // Copy "new" pressure & gradp just computed into "old" as well.
        //
        MultiFab::Copy(LevelData[lev]->get_old_data(Press_Type),
                       LevelData[lev]->get_new_data(Press_Type),
                       0, 0, 1, 0);

	int ng = (LevelData[lev]->get_new_data(Gradp_Type)).nGrow();
	MultiFab::Copy(LevelData[lev]->get_old_data(Gradp_Type),
                       LevelData[lev]->get_new_data(Gradp_Type),
                       0, 0, AMREX_SPACEDIM, ng);
    }


    if (verbose) {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	amrex::Print() << "Projection::initialPressureProject(): time: " << run_time << '\n';
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
      amrex::Print() << "Projection::initialSyncProject(): levels = " << c_lev << " - " << f_lev << '\n';

    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real stime = ParallelDescriptor::second();

    //
    // Gather data.
    //
    Vector<MultiFab*> vel(maxlev, nullptr);
    Vector<MultiFab*> phi(maxlev, nullptr);

    for (lev = c_lev; lev <= f_lev; lev++)
    {
        vel[lev] = &(LevelData[lev]->get_new_data(State_Type));
        // phi points to P_old and we use it as a temporary to
        // store phi resulting from the nodal projection, AKA
        // the pressure increment
        phi[lev] = &(LevelData[lev]->get_old_data(Press_Type));
        phi[lev] -> setVal(0.0);
    }

    const Real dt_inv = 1./dt;

    //
    // Source term
    //
    Vector<std::unique_ptr<MultiFab> > rhcc;
    if (have_divu)
    {
        rhcc.resize(maxlev);

        //
        // Set up rhcc for manual project.
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

	    // Needed for set_outflow_bcs(). MLMG ignores rh ghost cells.
            const int nghost = 1;
            rhcc[lev].reset(new MultiFab(amr_level.boxArray(),
                            amr_level.DistributionMap(),
                            1,nghost,MFInfo(),amr_level.Factory()));
            MultiFab* rhcclev = rhcc[lev].get();
            rhcclev->setVal(0);

            NavierStokesBase* ns = dynamic_cast<NavierStokesBase*>(&parent->getLevel(lev));

            AMREX_ASSERT(!(ns == 0));

            std::unique_ptr<MultiFab> divu (ns->getDivCond(nghost,strt_time));
            std::unique_ptr<MultiFab> dsdt (ns->getDivCond(nghost,strt_time+dt));

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*rhcclev,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
               const Box& bx = mfi.growntilebox();
               const auto& du       = divu->array(mfi);
               const auto& dsdt_arr = dsdt->array(mfi);
               const auto& rhcc_arr = rhcclev->array(mfi);
               amrex::ParallelFor(bx, [du,dsdt_arr,rhcc_arr,dt_inv]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                  rhcc_arr(i,j,k) = ( dsdt_arr(i,j,k) - du(i,j,k) ) * dt_inv;
               });
            }
        }
    }

    //
    // Set velocity bndry values to bogus values.
    //
    for (lev = c_lev; lev <= f_lev; lev++)
    {
        vel[lev]->setBndry(BogusValue,Xvel,AMREX_SPACEDIM);
        MultiFab &u_o = LevelData[lev]->get_old_data(State_Type);
        u_o.setBndry(BogusValue,Xvel,AMREX_SPACEDIM);
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
            LevelData[lev]->setPhysBoundaryValues(S_old[mfi],State_Type,prev_time,Xvel,Xvel,AMREX_SPACEDIM);
            LevelData[lev]->setPhysBoundaryValues(S_new[mfi],State_Type,curr_time,Xvel,Xvel,AMREX_SPACEDIM);
        }

        MultiFab& u_o = LevelData[lev]->get_old_data(State_Type);
        ConvertUnew(*vel[lev], u_o, dt, LevelData[lev]->boxArray());
    }

    if (OutFlowBC::HasOutFlowBC(phys_bc) && have_divu && do_outflow_bcs) {
        set_outflow_bcs(INITIAL_SYNC,phi,vel,
                        amrex::GetVecOfPtrs(rhcc),
                        sig,c_lev,f_lev,have_divu);
    }

    //
    // Scale initial sync projection variables.
    //
    for (lev = c_lev; lev <= f_lev; lev++)
    {
        scaleVar(INITIAL_SYNC,sig[lev],1,vel[lev],lev);

        if (have_divu && parent->Geom(0).IsRZ())
          radMultScal(lev,*(rhcc[lev]));
    }

    for (lev = f_lev; lev >= c_lev+1; lev--) {
      const BoxArray& crse_grids = vel[lev-1]->boxArray();
      const BoxArray& fine_grids = vel[lev  ]->boxArray();
      const DistributionMapping& crse_dmap = vel[lev-1]->DistributionMap();
      const DistributionMapping& fine_dmap = vel[lev  ]->DistributionMap();

      MultiFab v_crse(crse_grids, crse_dmap, AMREX_SPACEDIM, 1, MFInfo(), LevelData[lev-1]->Factory());
      MultiFab v_fine(fine_grids, fine_dmap, AMREX_SPACEDIM, 1, MFInfo(), LevelData[lev]->Factory());

      MultiFab::Copy(v_crse, *vel[lev-1], 0, 0, AMREX_SPACEDIM, 1);
      MultiFab::Copy(v_fine, *vel[lev  ], 0, 0, AMREX_SPACEDIM, 1);

      NavierStokesBase* ns = dynamic_cast<NavierStokesBase*>(LevelData[lev-1]);
      ns->average_down(v_fine, v_crse, 0, v_crse.nComp());

      MultiFab::Copy(*vel[lev-1], v_crse, 0, 0, AMREX_SPACEDIM, 1);
    }

    //
    // Project.
    //
    if (have_divu) {
        for (lev = c_lev; lev <= f_lev; lev++)
        {
            rhcc[lev]->mult(-1.0,0,1);
        }
    }

    bool increment_gp = true;
    doMLMGNodalProjection(c_lev, f_lev-c_lev+1, vel, phi, sig,
                          amrex::GetVecOfPtrs(rhcc),
                          {}, proj_tol, proj_abs_tol, increment_gp, 0, 0);

    //
    // Unscale initial sync projection variables.
    //
    for (lev = c_lev; lev <= f_lev; lev++)
        rescaleVar(INITIAL_SYNC,sig[lev],1,vel[lev],lev);

    //
    // Add correction at coarse and fine levels.
    // Only update new. NSB::resetState will take care of setting
    // old = new
    //
    for (lev = c_lev; lev <= f_lev; lev++)
    {
        MultiFab& P_new = LevelData[lev]->get_new_data(Press_Type);
        MultiFab::Add(P_new, *phi[lev], 0, 0, 1, 1);
    }

    //
    // Grad(P_new) incremented in doMLMGNodalProjection.
    //
    
    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - stime;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        amrex::Print() << "Projection::initialSyncProject(): time: " << run_time << '\n';
    }
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
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter Uoldmfi(Uold,true); Uoldmfi.isValid(); ++Uoldmfi)
  {
        const Box& bx=Uoldmfi.growntilebox(1);
        AMREX_ASSERT(grids[Uoldmfi.index()].contains(Uoldmfi.tilebox())==true);

        ConvertUnew(Unew[Uoldmfi],Uold[Uoldmfi],alpha,bx);
  }
}

//
// Convert U to an Accleration like quantity: Unew = (Unew - Uold)/alpha
//

void
Projection::ConvertUnew( FArrayBox &Unew, FArrayBox &Uold, Real alpha,
                              const Box &grd )
{
    AMREX_ASSERT(Unew.nComp() >= AMREX_SPACEDIM);
    AMREX_ASSERT(Uold.nComp() >= AMREX_SPACEDIM);
    AMREX_ASSERT(Unew.contains(grd) == true);
    AMREX_ASSERT(Uold.contains(grd) == true);

    const auto& unew = Unew.array();
    const auto& uold = Uold.array();
    const Real dt_inv = 1.0/alpha;
    amrex::ParallelFor(grd, AMREX_SPACEDIM, [unew,uold,dt_inv]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
       unew(i,j,k,n) -= uold(i,j,k,n);
       unew(i,j,k,n) *= dt_inv;
    });
}

//
// Add phi to P.
//

void
Projection::AddPhi (MultiFab&        p,
                    MultiFab&       phi)
{

  MultiFab::Add(p,phi,0,0,1,p.nGrow());
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
    AMREX_ASSERT((which_call == INITIAL_VEL  ) ||
              (which_call == INITIAL_PRESS) ||
              (which_call == INITIAL_SYNC ) ||
              (which_call == LEVEL_PROJ   ) ||
              (which_call == SYNC_PROJ    ) );

    if (sig != 0)
        AMREX_ASSERT(sig->nComp() == 1);
    if (vel != 0)
        AMREX_ASSERT(vel->nComp() >= AMREX_SPACEDIM);

    //
    // Convert sigma from rho to 1/rho if not INITIAL_PRESS.
    // nghosts info needed to avoid divide by zero.
    //
    if (sig != 0) {
#ifndef NDEBUG
#ifdef AMREX_USE_EB
      // fixme - deal with case where covered cells are set to zero
      //   there's probably a better way to handle this..
      EB_set_covered(*sig,0,1,sig->nGrow(),1.2345e40);
#endif
#endif
      sig->invert(1.0,sig_nghosts);
    }

    //
    // Scale by radius for RZ.
    //
    if (parent->Geom(0).IsRZ())
    {
        if (sig != 0)
            radMultScal(level,*sig);
        if (vel != 0)
            radMultVel(level,*vel);
    }
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
    AMREX_ASSERT((which_call == INITIAL_VEL  ) ||
              (which_call == INITIAL_PRESS) ||
              (which_call == INITIAL_SYNC ) ||
              (which_call == LEVEL_PROJ   ) ||
              (which_call == SYNC_PROJ    ) );

    if (sig != 0)
        AMREX_ASSERT(sig->nComp() == 1);
    if (vel != 0)
        AMREX_ASSERT(vel->nComp() >= AMREX_SPACEDIM);
    //
    // Divide by radius to rescale for RZ coordinates.
    //
    if (parent->Geom(0).IsRZ())
    {
        if (sig != 0)
            radDiv(level,*sig,0);
        if (vel != 0)
        {
            for (int n = 0; n < AMREX_SPACEDIM; n++)
                radDiv(level,*vel,n);
        }
    }
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
#if (AMREX_SPACEDIM < 3)
    AMREX_ASSERT(radius_grow >= mf.nGrow());

    const Box& domain = parent->Geom(level).Domain();
    const int* domlo  = domain.loVect();
    const int* domhi  = domain.hiVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfmfi(mf,true); mfmfi.isValid(); ++mfmfi)
    {
      AMREX_ASSERT(mf.box(mfmfi.index()) == mfmfi.validbox());

      const Box& bx = mfmfi.growntilebox();
      const int* lo = bx.loVect();
      const int* hi = bx.hiVect();
      Real* dat        = mf[mfmfi].dataPtr(0);
      const int* datlo = mf[mfmfi].loVect();
      const int* dathi = mf[mfmfi].hiVect();
      Real* rad        = &(*radius[level])[mfmfi.index()][0];
      const Box& gbx = mfmfi.validbox();
      int rlo   = (gbx.loVect())[0]-radius_grow;
      int rhi   = (gbx.hiVect())[0]+radius_grow;
      
      radmpyscal(lo,hi,dat,ARLIM(datlo),ARLIM(dathi),
		 domlo,domhi,rad,&rlo,&rhi);
    }
#endif
}

void
Projection::radMultVel (int       level,
                        MultiFab& mf)
{
#if (AMREX_SPACEDIM < 3)
    AMREX_ASSERT(radius_grow >= mf.nGrow());

    const Box& domain = parent->Geom(level).Domain();
    const int* domlo  = domain.loVect();
    const int* domhi  = domain.hiVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (int n = 0; n < AMREX_SPACEDIM; n++)
    {
       for (MFIter mfmfi(mf,true); mfmfi.isValid(); ++mfmfi)
       {
           AMREX_ASSERT(mf.box(mfmfi.index()) == mfmfi.validbox());
	   
	   const Box& bx = mfmfi.growntilebox();
	   const int* lo = bx.loVect();
	   const int* hi = bx.hiVect();
	   Real* dat        = mf[mfmfi].dataPtr(n);
	   const int* datlo = mf[mfmfi].loVect();
	   const int* dathi = mf[mfmfi].hiVect();
	   Real* rad        = &(*radius[level])[mfmfi.index()][0];
	   const Box& gbx = mfmfi.validbox();
	   int rlo   = (gbx.loVect())[0]-radius_grow;
	   int rhi   = (gbx.hiVect())[0]+radius_grow;
	   
	   radmpyvel(lo,hi,dat,ARLIM(datlo),ARLIM(dathi),
		     domlo,domhi,rad,&rlo,&rhi,&n);
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
#if (AMREX_SPACEDIM < 3)
    AMREX_ASSERT(comp >= 0 && comp < mf.nComp());
    AMREX_ASSERT(radius_grow >= mf.nGrow());

    const Box& domain = parent->Geom(level).Domain();
    const int* domlo  = domain.loVect();
    const int* domhi  = domain.hiVect();

    Real bogus_value = BogusValue;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfmfi(mf,true); mfmfi.isValid(); ++mfmfi)
    {
        AMREX_ASSERT(mf.box(mfmfi.index()) == mfmfi.validbox());

	const Box& bx = mfmfi.growntilebox();
	const int* lo = bx.loVect();
	const int* hi = bx.hiVect();
	Real* dat        = mf[mfmfi].dataPtr(comp);
	const int* datlo = mf[mfmfi].loVect();
	const int* dathi = mf[mfmfi].hiVect();
	Real* rad        = &(*radius[level])[mfmfi.index()][0];
	const Box& gbx = mfmfi.validbox();
	int rlo   = (gbx.loVect())[0]-radius_grow;
	int rhi   = (gbx.hiVect())[0]+radius_grow;
	
        fort_raddiv(lo,hi,dat,ARLIM(datlo),ARLIM(dathi),
		    domlo,domhi,rad,&rlo,&rhi,&bogus_value);

    }
#endif
}

//
// This projects the initial vorticity field (stored in pressure)
// to define an initial velocity field.
//
void
Projection::initialVorticityProject (int c_lev)
{
#if (AMREX_SPACEDIM == 2)
    int f_lev = parent->finestLevel();

    if (verbose) {
      amrex::Print() << "Projection::initialVorticityProject(): levels = " << c_lev
                     << "  " << f_lev << std::endl;
    }
    const Real strt_time = ParallelDescriptor::second();

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
                                       AMREX_SPACEDIM, 1));
        u_real[lev]->setVal(0);
        //
        // The vorticity is stored in the new pressure variable for now.
        //
        MultiFab& P_new = LevelData[lev]->get_new_data(Press_Type);

        rhnd[lev].reset(new MultiFab(P_new.boxArray(),
                                     P_new.DistributionMap(),
                                     1, 0));
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*rhnd[lev],true); mfi.isValid(); ++mfi)
	{
	  // rhnd has ng=0 as declared above
	  const Box& bx = mfi.tilebox();
	  (*rhnd[lev])[mfi].setVal<RunOn::Gpu>(0,bx);
	  (*rhnd[lev])[mfi].copy<RunOn::Gpu>(P_new[mfi], bx, 0, bx, 0, 1);
	}
    }

    //
    // Set BC for vorticity solve, save a copy of orig ones
    //
    BCRec phys_bc_save(phys_bc->lo(),phys_bc->hi());
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
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
    // FIXME -- need to think about what proj2 should really be. Don't
    // think we actually want to update Gradp here at all. And subsequent
    // initialVelocityProject will set P=Gp=0 anyway, right? 
    bool proj2 = !add_vort_proj;
    doMLMGNodalProjection(c_lev, f_lev-c_lev+1,
                          amrex::GetVecOfPtrs(u_real),
                          amrex::GetVecOfPtrs(p_real),
                          amrex::GetVecOfPtrs(s_real),
                          Vector<MultiFab*>(maxlev, nullptr),
                          amrex::GetVecOfPtrs(rhnd),
                          proj_tol, proj_abs_tol, proj2,
                          nullptr, nullptr, true);

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

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (int n = 0; n < AMREX_SPACEDIM; n++)
	{
	    for (MFIter mfi(*vel[lev],true); mfi.isValid(); ++mfi)
            {
                const Box& box = mfi.tilebox();
                if (add_vort_proj)
                {
                  (*vel[lev])[mfi].plus<RunOn::Gpu>((*u_real[lev])[mfi],box,Xvel+n,Xvel+idx[n], 1);
                }
                else
                {
                  (*vel[lev])[mfi].copy<RunOn::Gpu>((*u_real[lev])[mfi],box,Xvel+n,box,Xvel+idx[n], 1);
                }
            }
        }
    }

    //
    // Restore bcs
    //
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
      phys_bc->setLo(i,phys_bc_save.lo()[i]);
      phys_bc->setHi(i,phys_bc_save.hi()[i]);
    }

    if (verbose) {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	amrex::Print() << "Projection::initialVorticityProject(): time: " << run_time << '\n';
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

            // FIXME: this size may need adjusting
            ba.maxSize(32);

            DistributionMapping dm{ba};
            MultiFab phi_crse_strip(ba, dm, nCompPhi, 0);
            phi_crse_strip.setVal(0);
            const auto& phi_f_arr = phi_fine_strip[iface].array();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(phi_crse_strip); mfi.isValid(); ++mfi)
            {
                Box ovlp = amrex::coarsen(phi_fine_strip[iface].box(),ratio) & mfi.validbox();
                if (ovlp.ok())
                {
                    const auto& phi_c_arr = phi_crse_strip.array(mfi);
                    ParallelFor(ovlp, [phi_c_arr,phi_f_arr,ratio]
                    AMREX_GPU_DEVICE (int i, int j, int k) noexcept 
                    {
                       phi_c_arr(i,j,k) = phi_f_arr(i*ratio[0],j*ratio[1],k*ratio[2]);
                    });
                }
            }
            phi[lev]->copy(phi_crse_strip);
        }
    }
}

void
Projection::getStreamFunction (Vector<std::unique_ptr<MultiFab> >& /*phi*/)
{
  amrex::Abort("Projection::getStreamFunction not implemented");
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
    AMREX_ASSERT((which_call == INITIAL_VEL  ) ||
              (which_call == INITIAL_PRESS) ||
              (which_call == INITIAL_SYNC ) ||
              (which_call == LEVEL_PROJ   ) );

    if (which_call != LEVEL_PROJ)
      AMREX_ASSERT(c_lev == 0);

    if (verbose)
      amrex::Print() << "...setting outflow bcs for the nodal projection ... " << '\n';

    bool        hasOutFlow;
    Orientation outFaces[2*AMREX_SPACEDIM];
    Orientation outFacesAtThisLevel[maxlev][2*AMREX_SPACEDIM];

    int fine_level[2*AMREX_SPACEDIM];

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
//    const int   nCompVel  = AMREX_SPACEDIM;
//    const int   nCompDivu = 1;

    //
    // Determine the finest level such that the entire outflow face is covered
    // by boxes at this level (skip if doesnt touch, and bomb if only partially
    // covered).
    //
    Box state_strip[maxlev][2*AMREX_SPACEDIM];

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
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
          if (dir != outDir) temp_state_strip.grow(dir,1);

        const BoxArray& Lgrids               = parent->getLevel(lev).boxArray();
        const Box&      valid_state_strip    = temp_state_strip & domain;
        const BoxArray  uncovered_outflow_ba = amrex::complementIn(valid_state_strip,Lgrids);

        AMREX_ASSERT( !(uncovered_outflow_ba.size() &&
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
    AMREX_ASSERT(!(ns0 == 0));

    int Divu_Type, Divu;
    Real gravity = 0;

    if (which_call == INITIAL_SYNC || which_call == INITIAL_VEL)
    {
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
    AMREX_ASSERT(dynamic_cast<NavierStokesBase*>(LevelData[lev]) != nullptr);

    Box domain = parent->Geom(lev).Domain();

    const int ncStripWidth = 1;

    FArrayBox phi_fine_strip[2*AMREX_SPACEDIM];
    FArrayBox            rho[2*AMREX_SPACEDIM];
    
    const int ngrow = 1;

    for (int iface = 0; iface < numOutFlowFaces; iface++)
    {
        Box phi_strip = amrex::surroundingNodes(amrex::bdryNode(domain,
                                                outFacesAtThisLevel[iface],
                                                ncStripWidth));
        phi_fine_strip[iface].resize(phi_strip,1);
	phi_fine_strip[iface].setVal<RunOn::Gpu>(0.);

	rho[iface].resize(state_strip[iface],1);
	(*Sig_in).copyTo(rho[iface],0,0,1,ngrow);
    }

    if (std::fabs(gravity) > 0.)
      computeRhoG(rho,phi_fine_strip,
		  parent->Geom(lev),
		  outFacesAtThisLevel,numOutFlowFaces,gravity);

    // fixme - there's a cleaner way to do this
    for ( int iface = 0; iface < numOutFlowFaces; iface++)
    {
        BoxArray phi_fine_strip_ba(phi_fine_strip[iface].box());
        // FIXME: this size may need adjusting
        phi_fine_strip_ba.maxSize(32);
        DistributionMapping dm {phi_fine_strip_ba};
        MultiFab phi_fine_strip_mf(phi_fine_strip_ba,dm,1,0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(phi_fine_strip_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.tilebox();
          AMREX_ASSERT((phi_fine_strip[iface].box()).contains(bx));
          const auto& phi_f_mf = phi_fine_strip_mf.array(mfi); 
          const auto& phi_f    = phi_fine_strip[iface].array();
          amrex::ParallelFor(bx, [phi_f_mf,phi_f]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
             phi_f_mf(i,j,k) = phi_f(i,j,k);
          });
        }

        phi[lev]->copy(phi_fine_strip_mf);
    }

    if (lev > c_lev)
    {
      putDown(phi, phi_fine_strip, c_lev, lev, outFacesAtThisLevel,
                numOutFlowFaces, ncStripWidth);
    }
}


void
Projection::computeRhoG(FArrayBox*         rhoFab,
			FArrayBox*         phiFab,
			const Geometry&    geom, 
			Orientation*       outFaces,
			int                numOutFlowFaces,
			Real               gravity)
{
    AMREX_ASSERT(std::fabs(gravity) > 0.);
  
    for (int iface = 0; iface < numOutFlowFaces; iface++)
    {
      int outDir             = outFaces[iface].coordDir();
      Orientation::Side side = outFaces[iface].faceDir();

      if (outDir == (AMREX_SPACEDIM-1))
      {
	  if (side == Orientation::high) {
	    //
	    // Hydrostatic pressure == 0 here, given IAMR definition of gravity.
	    // Do nothing, since phi already initialized to zero
	    //
	  } else {
	    amrex::Abort("Projection::computeRhoG : Simulation box has outflow boundary condition on the bottom and gravity != 0. If this is really the desired configuration, just comment out this Abort");
	  }
      }
      else // integrate rho * g* dh
      {
	  const auto   lo = amrex::lbound(phiFab[iface].box());
	  const auto   hi = amrex::ubound(phiFab[iface].box());
	  const auto& phi = phiFab[iface].array();
	  const auto& rho = rhoFab[iface].array();
	  const Real dh = geom.CellSize(AMREX_SPACEDIM-1);
	  
	  auto add_rhog = [gravity, dh] ( Real rho1, Real rho2,
					  Real& rhog_i, Real& phi_i )
	  {
	    Real rhoExt = 0.5*(3.*rho1-rho2);
	    rhog_i -= gravity * rhoExt * dh;
	    phi_i  += rhog_i;
	  };


#if (AMREX_SPACEDIM == 2)
	  //
	  // Only possibilities are XLO face or XHI
	  //
	  // Ok to only use low index of phi because phi is only one
	  // node wide in i-direction.
	  //
	  AMREX_ASSERT( lo.x==hi.x );
	  int i = lo.x;
	
	  Real rhog = 0.;
	  //
	  // Note that the integral here prevents parallelization
	  //
	  if (side == Orientation::low)
	  {
	    for (int j = hi.y-1; j >= lo.y; j--) {
	      add_rhog(rho(i,j,0),rho(i+1,j,0),rhog,phi(i,j,0));
	    }
	  }
	  else
	  {
	    for (int j = hi.y-1; j >= lo.y; j--) {
	      add_rhog(rho(i-1,j,0),rho(i-2,j,0),rhog,phi(i,j,0));
	    }
	  }
#else
	  const Box& domain = geom.Domain();
	  const auto domlo = amrex::lbound(domain);
	  const auto domhi = amrex::ubound(domain);
	  
	  // fixme? Could make use of NSB::m_bcrec_scalars here.
	  int        lo_bc[AMREX_SPACEDIM];
	  int        hi_bc[AMREX_SPACEDIM];
	  //
	  // change from phys_bcs of Inflow, SlipWall, etc.
	  // to mathematical bcs of EXT_DIR, FOEXTRAP, etc.
	  //
	  for (int i = 0; i < AMREX_SPACEDIM; i++)
	  {
	      const int* lbc = phys_bc->lo();
	      const int* hbc = phys_bc->hi();
	      
	      lo_bc[i]=scalar_bc[lbc[i]];
	      hi_bc[i]=scalar_bc[hbc[i]];
	  }

	  //
	  // fixme? - TODO: Could parallelize here by dividing the loop over i (or j)
	  // only and thus the k integration stays intact. However, would want to move
	  // this declaration of rho_i, rho_ii to ensure each k integration has it's own 
	  // copy.
	  //
	  Real rho_i, rho_ii;
	  
	  if ( outDir == int(Direction::x) )
	  {
	      // 
	      // Ok to only use low index of phi because phi is only one
	      // node wide in i-direction.
	      //
	      AMREX_ASSERT( lo.x==hi.x );
	      int i = lo.x;
	      
	      bool has_extdir_lo = (lo.y==domlo.y   && lo_bc[1]==BCType::ext_dir);
	      bool has_extdir_hi = (hi.y==domhi.y+1 && hi_bc[1]==BCType::ext_dir);
	      bool has_hoextrap_lo = (lo.y==domlo.y   && lo_bc[1]==BCType::hoextrap);
	      bool has_hoextrap_hi = (hi.y==domhi.y+1 && hi_bc[1]==BCType::hoextrap);
	      bool has_foextrap_lo = (lo.y==domlo.y   && lo_bc[1]==BCType::foextrap);
	      bool has_foextrap_hi = (hi.y==domhi.y+1 && hi_bc[1]==BCType::foextrap);

	      //
	      // If there are any of the above mentioned bcs, then we'll need to handle
	      // edges separately in accordance with those conditions, so set bounds for
	      // loop over j accordingly.
	      //
	      int jlo, jhi;
	      if ( has_extdir_lo || has_hoextrap_lo || has_foextrap_lo ) {
		jlo = lo.y+1;
	      } else {
		jlo = lo.y;
	      }
	      if ( has_extdir_hi || has_hoextrap_hi || has_foextrap_hi ) {
		jhi = hi.y-1;
	      } else {
		jhi = hi.y;
	      }

	      //
	      // xlo face
	      //
	      if (side == Orientation::low)
	      {
		for (int j = jlo; j <= jhi; j++)
		{
		  Real rhog = 0.;
		  
		  for (int k = hi.z-1; k >= lo.z; k--) {
		    rho_i  = 0.5 * (rho(i  ,j,k) + rho(i  ,j-1,k));
		    rho_ii = 0.5 * (rho(i+1,j,k) + rho(i+1,j-1,k));
		    add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		  }
		}
		//
		// Now compute y-edges if needed
		//
		if ( has_extdir_lo || has_hoextrap_lo || has_foextrap_lo )
		{
		  int j = lo.y; 
		  Real rhog = 0.;
		  
		  if ( has_extdir_lo ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i  = rho(i  ,j-1,k);
		      rho_ii = rho(i+1,j-1,k);
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } else if ( has_hoextrap_lo ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i  = 0.5*(3.*rho(i  ,j,k) - rho(i  ,j+1,k));
		      rho_ii = 0.5*(3.*rho(i+1,j,k) - rho(i+1,j+1,k));
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } else if ( has_foextrap_lo ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i  = rho(i  ,j,k);
		      rho_ii = rho(i+1,j,k);
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  }
		}

		if ( has_extdir_hi || has_hoextrap_hi || has_foextrap_hi )
		{
		  int j = hi.y;
		  Real rhog = 0;
		  
		  if ( has_extdir_hi ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i  = rho(i  ,j,k);
		      rho_ii = rho(i+1,j,k);
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } else if ( has_hoextrap_hi ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i  = 0.5*(3.*rho(i  ,j-1,k) - rho(i  ,j-2,k));
		      rho_ii = 0.5*(3.*rho(i+1,j-1,k) - rho(i+1,j-2,k));
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } else if ( has_foextrap_hi ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i  = rho(i  ,j-1,k);
		      rho_ii = rho(i+1,j-1,k);
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } 
		}
	      }
	      else // xhi face 
	      {
		for (int j = jlo; j <= jhi; j++)
		{
		  Real rhog = 0.;
		  
		  for (int k = hi.z-1; k >= lo.z; k--) {
		    rho_i   = 0.5 * (rho(i-1,j,k) + rho(i-1,j-1,k));
		    rho_ii = 0.5 * (rho(i-2,j,k) + rho(i-2,j-1,k));
		    add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		  }
		}
		//
		// Now compute y-edges if needed
		//
		if ( has_extdir_lo || has_hoextrap_lo || has_foextrap_lo )
		{
		  int j = lo.y; 
		  Real rhog = 0.;
		  
		  if ( has_extdir_lo ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i   = rho(i-1,j-1,k);
		      rho_ii = rho(i-2,j-1,k);
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } else if ( has_hoextrap_lo ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i   = 0.5*(3.*rho(i-1,j,k) - rho(i-1,j+1,k));
		      rho_ii = 0.5*(3.*rho(i-2,j,k) - rho(i-2,j+1,k));
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } else if ( has_foextrap_lo ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i   = rho(i-1,j,k);
		      rho_ii = rho(i-2,j,k);
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } 
		}

		if ( has_extdir_hi || has_hoextrap_hi || has_foextrap_hi )
		{
		  int j = hi.y;
		  Real rhog = 0;
		  
		  if ( has_extdir_hi ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i   = rho(i-1,j,k);
		      rho_ii = rho(i-2,j,k);
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } else if ( has_hoextrap_hi ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i   = 0.5*(3.*rho(i-1,j-1,k) - rho(i-1,j-2,k));
		      rho_ii = 0.5*(3.*rho(i-2,j-1,k) - rho(i-2,j-2,k));
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } else if ( has_foextrap_hi ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i   = rho(i-1,j-1,k);
		      rho_ii = rho(i-2,j-1,k);
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } 
		}
	      }
	  }
	  else // ( outDir == Direction::y ) 
	  {
	      //
	      // Ok to only use low index of phi because phi is only one
	      // node wide in i-direction.
	      //
	      AMREX_ASSERT( lo.y==hi.y );
	      int j = lo.y;

	      bool has_extdir_lo = (lo.x==domlo.x   && lo_bc[0]==BCType::ext_dir);
	      bool has_extdir_hi = (hi.x==domhi.x+1 && hi_bc[0]==BCType::ext_dir);
	      bool has_hoextrap_lo = (lo.x==domlo.x   && lo_bc[0]==BCType::hoextrap);
	      bool has_hoextrap_hi = (hi.x==domhi.x+1 && hi_bc[0]==BCType::hoextrap);
	      bool has_foextrap_lo = (lo.x==domlo.x   && lo_bc[0]==BCType::foextrap);
	      bool has_foextrap_hi = (hi.x==domhi.x+1 && hi_bc[0]==BCType::foextrap);
	      //
	      // If there are any of the above mentioned bcs, then we'll need to handle
	      // edges separately in accordance with those conditions, so set bounds for
	      // loop over i accordingly.
	      //
	      int ilo, ihi;
	      if ( has_extdir_lo || has_hoextrap_lo || has_foextrap_lo ) {
		ilo = lo.x+1;
	      } else {
		ilo = lo.x;
	      }
	      if ( has_extdir_hi || has_hoextrap_hi || has_foextrap_hi ) {
		ihi = hi.x-1;
	      } else {
		ihi = hi.x;
	      }
	      //
	      // ylo face
	      //
	      if (side == Orientation::low)
	      {
		for (int i = ilo; i <= ihi; i++)
		{
		  Real rhog = 0.;
		  
		  for (int k = hi.z-1; k >= lo.z; k--) {
		    rho_i   = 0.5 * (rho(i,j  ,k) + rho(i-1,j ,k));
		    rho_ii = 0.5 * (rho(i,j+1,k) + rho(i-1,j+1,k));
		    add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		  }
		}
		//
		// Now compute x-edges if needed
		//
		if ( has_extdir_lo || has_hoextrap_lo || has_foextrap_lo )
		{
		  int i = lo.x; 
		  Real rhog = 0.;
		  
		  if ( has_extdir_lo ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i   = rho(i-1,j  ,k);
		      rho_ii = rho(i-1,j+1,k);
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } else if ( has_hoextrap_lo ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i   = 0.5*(3.*rho(i,j  ,k) - rho(i+1,j  ,k));
		      rho_ii = 0.5*(3.*rho(i,j+1,k) - rho(i+1,j+1,k));
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } else if ( has_foextrap_lo ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i   = rho(i,j  ,k);
		      rho_ii = rho(i,j+1,k);
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  }
		}

		if ( has_extdir_hi || has_hoextrap_hi || has_foextrap_hi )
		{
		  int i = hi.x;
		  Real rhog = 0;
		  
		  if ( has_extdir_hi ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i   = rho(i,j  ,k);
		      rho_ii = rho(i,j+1,k);
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } else if ( has_hoextrap_hi ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i   = 0.5*(3.*rho(i-1,j  ,k) - rho(i-2,j  ,k));
		      rho_ii = 0.5*(3.*rho(i-1,j+1,k) - rho(i-2,j+1,k));
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } else if ( has_foextrap_hi ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i   = rho(i-1,j  ,k);
		      rho_ii = rho(i-1,j+1,k);
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } 
		}
	      }
	      else // yhi face 
	      {
		for (int i = ilo; i <= ihi; i++)
		{
		  Real rhog = 0.;
		  
		  for (int k = hi.z-1; k >= lo.z; k--) {
		    rho_i   = 0.5 * (rho(i,j-1,k) + rho(i-1,j-1,k));
		    rho_ii = 0.5 * (rho(i,j-1,k) + rho(i-1,j-2,k));
		    add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		  }
		}
		//
		// Now compute y-edges if needed
		//
		if ( has_extdir_lo || has_hoextrap_lo || has_foextrap_lo )
		{
		  int i = lo.x; 
		  Real rhog = 0.;
		  
		  if ( has_extdir_lo ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i   = rho(i-1,j-1,k);
		      rho_ii = rho(i-2,j-1,k);
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } else if ( has_hoextrap_lo ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i   = 0.5*(3.*rho(i,j-1,k) - rho(i+1,j-1,k));
		      rho_ii = 0.5*(3.*rho(i,j-2,k) - rho(i+1,j-2,k));
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } else if ( has_foextrap_lo ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i   = rho(i,j-1,k);
		      rho_ii = rho(i,j-2,k);
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } 
		}

		if ( has_extdir_hi || has_hoextrap_hi || has_foextrap_hi )
		{
		  int i = hi.x;
		  Real rhog = 0;
		  
		  if ( has_extdir_hi ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i   = rho(i,j-1,k);
		      rho_ii = rho(i,j-2,k);
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } else if ( has_hoextrap_hi ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i   = 0.5*(3.*rho(i-1,j-1,k) - rho(i-2,j-1,k));
		      rho_ii = 0.5*(3.*rho(i-1,j-2,k) - rho(i-2,j-2,k));
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } else if ( has_foextrap_hi ) {
		    for (int k = hi.z-1; k >= lo.z; k--) {
		      rho_i   = rho(i-1,j-1,k);
		      rho_ii = rho(i-1,j-2,k);
		      add_rhog(rho_i, rho_ii, rhog, phi(i,j,k));
		    }
		  } 
		}
	      } // endif over hi/low sides 
	  } // endif over directions
#endif
      } // endif integrate rho*g*dh
    } // end loop over outflow faces
}

//
// Given vel, rhcc, rhnd, & sig, this solves Div (sig * Grad phi) = Div vel + (rhcc + rhnd).
// On return, vel becomes vel  - sig * Grad phi.
//
void Projection::doMLMGNodalProjection (int c_lev, int nlevel,
                                        const Vector<MultiFab*>& vel,
                                        const Vector<MultiFab*>& phi,
                                        const Vector<MultiFab*>& sig,
                                        const Vector<MultiFab*>& rhcc,
                                        const Vector<MultiFab*>& rhnd,
                                        Real rel_tol, Real abs_tol,
                                        bool increment_gp,
                                        MultiFab* sync_resid_crse,
                                        MultiFab* sync_resid_fine,
                                        bool doing_initial_vortproj)
{
    BL_PROFILE("Projection:::doMLMGNodalProjection()");

    int f_lev = c_lev + nlevel - 1;

    Vector<MultiFab> vel_test(nlevel);
    Vector<MultiFab> phi_test(nlevel);

    AMREX_ASSERT(vel[c_lev]->nGrow() >= 1);
    AMREX_ASSERT(vel[f_lev]->nGrow() >= 1);
    AMREX_ASSERT(phi[c_lev]->nGrow() == 1);
    AMREX_ASSERT(phi[f_lev]->nGrow() == 1);
    // MLMG does not copy any ghost cells from sig, rhcc or rhnd; fills ghost cells internally 
    // AMREX_ASSERT(sig[c_lev]->nGrow() == 1);
    // AMREX_ASSERT(sig[f_lev]->nGrow() == 1);

    AMREX_ASSERT(sig[c_lev]->nComp() == 1);
    AMREX_ASSERT(sig[f_lev]->nComp() == 1);

    if (sync_resid_crse != 0) {
        AMREX_ASSERT(nlevel == 1);
        AMREX_ASSERT(c_lev < parent->finestLevel());
    }

    if (sync_resid_fine != 0) {
        AMREX_ASSERT((nlevel == 1 || nlevel == 2));
        AMREX_ASSERT(c_lev > 0);
    }

    if (!rhcc.empty() )
        AMREX_ALWAYS_ASSERT(rhcc[c_lev]->boxArray().ixType().cellCentered());

    if (!rhnd.empty() )
        AMREX_ALWAYS_ASSERT(rhnd[c_lev]->boxArray().ixType().nodeCentered());

    set_boundary_velocity(c_lev, nlevel, vel, true);

    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (parent->Geom(0).isPeriodic(idim))
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

    //
    // Setup objects for projection
    //
    Vector<Geometry>            mg_geom(nlevel);
    Vector<BoxArray>            mg_grids(nlevel);
    Vector<DistributionMapping> mg_dmap(nlevel);

    for (int lev(0); lev < nlevel; lev++)
    {
        mg_geom[lev] = parent->Geom(lev+c_lev);
        mg_grids[lev] = parent->boxArray(lev+c_lev);
        mg_dmap[lev] = LevelData[lev+c_lev]->get_new_data(State_Type).DistributionMap();
    }

    // Setup infos to pass to linear operator
    LPInfo info;
    //Fixme
    // does this max_coarsening level need to match the one in main.cpp????
    int max_coarsening_level(30);
    info.setMaxCoarseningLevel(max_coarsening_level);
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    // metric term stuff doesn't get used at all for nodal
    info.setMetricTerm(false);

    //
    // Setup variables to use in projection
    //
    Vector<MultiFab*> phi_rebase(phi.begin()+c_lev, phi.begin()+c_lev+nlevel);
    Vector<MultiFab*> vel_rebase{vel.begin()+c_lev, vel.begin()+c_lev+nlevel};
    Vector<MultiFab*> sigma_rebase(sig.begin()+c_lev, sig.begin()+c_lev+nlevel);

    Vector<const MultiFab*> rhnd_rebase;
    Vector<MultiFab*> rhcc_rebase;

    if (!rhnd.empty())
    {
        rhnd_rebase.assign(rhnd.begin()+c_lev, rhnd.begin()+c_lev+nlevel);
    }

    if (!rhcc.empty())
    {
        rhcc_rebase.assign(rhcc.begin()+c_lev, rhcc.begin()+c_lev+nlevel);
    }

    // Setup nodal projector object
    NodalProjector  nodal_projector(vel_rebase, GetVecOfConstPtrs(sigma_rebase), mg_geom, info, rhcc_rebase, rhnd_rebase);
    nodal_projector.setDomainBC(mlmg_lobc, mlmg_hibc);

// WARNING: we set the strategy to Sigma to get exactly the same results as the no EB code
// when we don't have interior geometry
//  nodal_projector.getLinOp().setCoarseningStrategy(MLNodeLaplacian::CoarseningStrategy::Sigma);
    
// MLNodeLaplacian.define() will set is_rz based on geom. Do we really need this and
    // the ability to set is_rz separately from inputs file?
    // Also, what of LPInfo::has_metric_term? why not just use that instead of is_rz??
#if (AMREX_SPACEDIM == 2)
    if (rz_correction)
    {
        nodal_projector.getLinOp().setRZCorrection(parent->Geom(0).IsRZ());
    }
#endif
    nodal_projector.getLinOp().setGaussSeidel(use_gauss_seidel);
    nodal_projector.getLinOp().setHarmonicAverage(use_harmonic_average);
    nodal_projector.getMLMG().setMaxFmgIter(max_fmg_iter);

    if (sync_resid_fine != 0)
    {
        nodal_projector.setSyncResidualFine(sync_resid_fine);
    }
    if (sync_resid_crse != 0)
    {
        nodal_projector.setSyncResidualCrse(sync_resid_crse, parent->refRatio(c_lev), parent->boxArray(c_lev+1));
    }

    //
    // Project to get new P and update velocity
    //
    nodal_projector.project(phi_rebase,rel_tol,abs_tol);
    
    //
    // Update gradP
    //
    const auto gradphi = nodal_projector.getGradPhi();

    for (int lev = 0; lev < nlevel; lev++)
    {
      NavierStokesBase& ns = *dynamic_cast<NavierStokesBase*>(LevelData[lev+c_lev]);
      MultiFab& Gp = ns.get_new_data(Gradp_Type);
      
      if ( increment_gp )
      {
	//
	// Add a correction to Gradp
	//
        MultiFab::Add(Gp, *gradphi[lev], 0, 0, AMREX_SPACEDIM, 0);
      }
      else
      {
	//
	// Replace Gradp with the gradient(P) computed in MLMG
	//
	MultiFab::Copy(Gp, *gradphi[lev], 0, 0, AMREX_SPACEDIM, 0);

      }
      //
      // FIXME - could we get away with only FillPatching in predict_velocity
      // and initialPressureProject? I think this would depend on the definition
      // of properly nested... For now, be safe and just fill them.
      // Fill ghost cells
      //
      const Real& time = (ns.state)[Gradp_Type].curTime();
      NavierStokesBase::FillPatch(ns, Gp, Gp.nGrow(), time, Gradp_Type, 0, AMREX_SPACEDIM);
    }
}

// Set velocity in ghost cells to zero except for inflow
void Projection::set_boundary_velocity(int c_lev, int nlevel, const Vector<MultiFab*>& vel,
                                       bool inflowCorner)
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

    for (int idir=0; idir<AMREX_SPACEDIM; idir++) {

      if (lo_bc[idir] != Inflow && hi_bc[idir] != Inflow) {
	vel[lev]->setBndry(0.0, Xvel+idir, 1);
      }
      else {
	//fixme: is it worth the overhead to have threads here?
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(*vel[lev]); mfi.isValid(); ++mfi) {
	  int i = mfi.index();

	  FArrayBox& v_fab = (*vel[lev])[mfi];

	  const Box& reg = grids[i];
	  const Box& bxg1 = amrex::grow(reg,1);
	  BoxList bxlist(reg);

	  //If tiling only need to redefine these (all the rest can stay the same):
	  // const Box& bxg1 = mfi.growntilebox(1);
	  // const Box& tile = mfi.tilebox();
	  // BoxList bxlist(tile);

	  if (lo_bc[idir] == Inflow && reg.smallEnd(idir) == domainBox.smallEnd(idir)) {
	    Box bx;                // bx is the region we *protect* from zero'ing
	    bx = amrex::adjCellLo(reg, idir);

	    if (inflowCorner) {

              for (int odir = 0; odir < AMREX_SPACEDIM; odir++) {
		if (odir != idir) {
                    if (geom.isPeriodic(odir)) bx.grow(odir,1);
                    if (reg.bigEnd  (odir) != domainBox.bigEnd  (odir) ) bx.growHi(odir,1);
                    if (reg.smallEnd(odir) != domainBox.smallEnd(odir) ) bx.growLo(odir,1);
		}
              }
	    }
	    bxlist.push_back(bx);
	  }

	  if (hi_bc[idir] == Inflow && reg.bigEnd(idir) == domainBox.bigEnd(idir)) {
	    Box bx;                // bx is the region we *protect* from zero'ing
	    bx = amrex::adjCellHi(reg, idir);

	    if (inflowCorner) {

              for (int odir = 0; odir < AMREX_SPACEDIM; odir++) {
		if (odir != idir)
		  {
                    if (geom.isPeriodic(odir)) bx.grow(odir,1);
                    if (reg.bigEnd  (odir) != domainBox.bigEnd  (odir) ) bx.growHi(odir,1);
                    if (reg.smallEnd(odir) != domainBox.smallEnd(odir) ) bx.growLo(odir,1);
		  }
              }
	    }

	    bxlist.push_back(bx);
	  }

	  BoxList bxlist2 = amrex::complementIn(bxg1, bxlist);

	  for (BoxList::iterator it=bxlist2.begin(); it != bxlist2.end(); ++it) {
            Box ovlp = *it & v_fab.box();
            if (ovlp.ok()) {
		v_fab.setVal<RunOn::Gpu>(0.0, ovlp, Xvel+idir, 1);
            }
	  }
	}
      }
    }
  }
}
