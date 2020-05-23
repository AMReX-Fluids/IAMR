

#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <NavierStokesBase.H>
#include <NS_BC.H>
#include <AMReX_BLProfiler.H>
#include <Projection.H>
#include <PROJECTION_F.H>
#include <NAVIERSTOKES_F.H>
#include <ProjOutFlowBC.H>

#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>

#include <AMReX_NodalProjector.H>

//fixme, for writesingle level plotfile
#include<AMReX_PlotFileUtil.H>

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
int  Projection::anel_grow           = 1;

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
	amrex::Abort("Must use proj_2==1 due to new gravity and outflow stuff. proj_2!=1 no longer supported.");

    pp.query("Pcode",               P_code);
    if (P_code >=0 )
      amrex::Abort("proj.Pcode is no more. Use nodal_proj.verbose.");

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

#ifdef AMREX_USE_EB
    const auto& _ebfactory =
      dynamic_cast<EBFArrayBoxFactory const&>(LevelData[level]->Factory());
    ebfactory[level] = &_ebfactory;
#endif
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


void
Projection::build_anelastic_coefficient (int      level,
					 Real**& _anel_coeff)
{
  const BoxArray& grids = parent->getLevel(level).boxArray();
  const int N = grids.size();
  _anel_coeff = new Real*[N];
  for (int i = 0; i < grids.size(); i++)
  {
    const int jlo = grids[i].smallEnd(BL_SPACEDIM-1)-anel_grow;
    const int jhi = grids[i].bigEnd(BL_SPACEDIM-1)+anel_grow;
    const int len = jhi - jlo + 1;

    _anel_coeff[i] = new Real[len];

    // FIXME!
    // This is just a placeholder for testing. Should create (problem
    // dependent) build_coefficient function in problem directory
    // ...Perhaps also need to worry about deleting anel_coeff
    // Also not sure why Projection and MacProj have separate anel_coeff
    // arrays, since they both appear to be cell centered.
    for (int j=0; j<len; j++)
      _anel_coeff[i][j] = 0.05*(jlo+j);
  }
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
    }

    const int nGrow = (level == 0  ?  0  :  -1);
    //
    // MultiFab::setVal() won't work here bc nGrow could be <0
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(P_new,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.growntilebox(nGrow);
        P_new[mfi].setVal<RunOn::Host>(0.0,bx,0,1);
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
    U_new.mult(dt_inv,0,BL_SPACEDIM,1);
    if (have_divu)
      divusource->mult(dt_inv,0,1,divusource->nGrow());


#ifdef AMREX_USE_EB
    MultiFab& Gp = ns->getGradP();
    Gp.FillBoundary(geom.periodicity());
#else
    MultiFab Gp(grids,dmap,BL_SPACEDIM,1);
    ns->getGradP(Gp, prev_pres_time);
#endif


#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(rho_half,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.growntilebox(1);
	const FArrayBox& rhofab = rho_half[mfi];

#ifdef AMREX_USE_EB
	FArrayBox Gpfab(bx,BL_SPACEDIM);
	Gpfab.copy<RunOn::Host>((Gp)[mfi],0,0,BL_SPACEDIM);
#else
	FArrayBox& Gpfab = Gp[mfi];
#endif

	for (int i = 0; i < BL_SPACEDIM; i++) {
	  Gpfab.divide<RunOn::Host>(rhofab,bx,0,i,1);
	}

	U_new[mfi].plus<RunOn::Host>(Gpfab,bx,0,0,BL_SPACEDIM);
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
    int is_rz = (geom.IsRZ() ? 1 : 0);

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

    bool proj2 = true;
    doMLMGNodalProjection(level, 1, vel, phi, sig, rhcc, {}, proj_tol,
			  proj_abs_tol, proj2,
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
    U_new.mult(dt,0,BL_SPACEDIM,1);

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
      amrex::Print() << "Projection::MLsyncProject(): levels = " << c_lev << ", " << c_lev+1 << '\n';
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

    const Geometry& fine_geom = parent->Geom(c_lev+1);

    NavierStokesBase* ns = dynamic_cast<NavierStokesBase*>(LevelData[c_lev]);
    ns->average_down(*vel[c_lev+1],*vel[c_lev],0,AMREX_SPACEDIM);
    //
    // // restrict_level(v_crse, v_fine, ratio);
    // amrex::average_down(*vel[c_lev+1],*vel[c_lev],fine_geom,crse_geom,
    //                      0, BL_SPACEDIM, ratio);

    ns->average_down(*sig[c_lev+1],*sig[c_lev],0,sig[c_lev]->nComp());
    // // restrict_level(*sig[c_lev], *sig[c_lev+1], ratio);
    // amrex::average_down(*sig[c_lev+1],*sig[c_lev],fine_geom,crse_geom,
    //                        0, sig[c_lev]->nComp(), ratio);

    MultiFab* sync_resid_crse = 0;
    std::unique_ptr<MultiFab> sync_resid_fine;

    if (c_lev > 0 &&  crse_iteration == crse_dt_ratio)
    {
        int ngrow = parent->MaxRefRatio(c_lev-1) - 1;
        sync_resid_fine.reset(new MultiFab(Pgrids_crse,Pdmap_crse,1,ngrow));
        sync_resid_fine->setVal(0.);
    }

    bool proj2 = true;
    doMLMGNodalProjection(c_lev, 2, vel,
                          amrex::GetVecOfPtrs(phi),
                          sig, rhcc, rhnd_vec, sync_tol, proj_abs_tol, proj2,
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
    MultiFab::Saxpy(vel_crse,dt_crse,Vsync,0,0,AMREX_SPACEDIM,1);
    MultiFab::Saxpy(vel_fine,dt_crse,V_corr,0,0,AMREX_SPACEDIM,1);

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	amrex::Print() << "Projection::MLsyncProject(): levels = " << c_lev << ", " << c_lev+1
		       << ", time: " << run_time << '\n';
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
	const int nghost = (OutFlowBC::HasOutFlowBC(phys_bc) && do_outflow_bcs && have_divu) ? 1 : 0;

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
			     Density,
			     0,
			     1,
			     nghost);
            }
            else
            {
                sig[lev]->setVal(1.,nghost);
            }
        }

        Vector<std::unique_ptr<MultiFab> > rhcc(maxlev);

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

		NavierStokesBase* ns = dynamic_cast<NavierStokesBase*>(LevelData[lev]);
		BL_ASSERT(!(ns == 0));

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
                
        bool proj2 = true;
        if (!have_divu)
        {
            doMLMGNodalProjection(c_lev, f_lev-c_lev+1, vel, phi,
                                  amrex::GetVecOfPtrs(sig),
                                  {},
                                  {},
                                  proj_tol, proj_abs_tol, proj2, 0, 0);
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
                                  proj_tol, proj_abs_tol, proj2, 0, 0);
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
#ifdef AMREX_USE_EB
            // gradP updated in MLMGNodalProjection so need to reset to zero here
            NavierStokesBase* ns = dynamic_cast<NavierStokesBase*>(LevelData[lev]);
            MultiFab& Gp = ns->getGradP();
            Gp.setVal(0.);
#endif
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
        phi[lev] = &(LevelData[lev]->get_old_data(Press_Type));

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
	raii.push_back(std::unique_ptr<MultiFab>(new MultiFab(grids, dmap, BL_SPACEDIM, 1,MFInfo(),LevelData[lev]->Factory())));
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
    bool proj2 = true;
    Vector<MultiFab*> rhcc(0);
    doMLMGNodalProjection(c_lev, f_lev-c_lev+1, vel, phi,
                          amrex::GetVecOfPtrs(sig),
                          rhcc, {},
                          proj_tol, proj_abs_tol, proj2, 0, 0);

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

            BL_ASSERT(!(ns == 0));

            std::unique_ptr<MultiFab> divu (ns->getDivCond(nghost,strt_time));
            std::unique_ptr<MultiFab> dsdt (ns->getDivCond(nghost,strt_time+dt));

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(*rhcclev,true); mfi.isValid(); ++mfi)
            {
       	        const Box& bx = mfi.growntilebox();
	        FArrayBox& dsdtfab = (*dsdt)[mfi];

                dsdtfab.minus<RunOn::Host>((*divu)[mfi],bx,bx,0,0,1);
                dsdtfab.mult<RunOn::Host>(dt_inv,bx);
                (*rhcclev)[mfi].copy<RunOn::Host>(dsdtfab,bx,0,bx,0,1);
            }
        }
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

      MultiFab v_crse(crse_grids, crse_dmap, BL_SPACEDIM, 1, MFInfo(), LevelData[lev-1]->Factory());
      MultiFab v_fine(fine_grids, fine_dmap, BL_SPACEDIM, 1, MFInfo(), LevelData[lev]->Factory());

      const Geometry& fine_geom = parent->Geom(lev  );
      const Geometry& crse_geom = parent->Geom(lev-1);

      MultiFab::Copy(v_crse, *vel[lev-1], 0, 0, BL_SPACEDIM, 1);
      MultiFab::Copy(v_fine, *vel[lev  ], 0, 0, BL_SPACEDIM, 1);

      NavierStokesBase* ns = dynamic_cast<NavierStokesBase*>(LevelData[lev-1]);
      ns->average_down(v_fine, v_crse, 0, v_crse.nComp());
      // amrex::average_down(v_fine,v_crse,fine_geom,crse_geom,
      //                      0, v_crse.nComp(), parent->refRatio(lev-1));

      MultiFab::Copy(*vel[lev-1], v_crse, 0, 0, BL_SPACEDIM, 1);
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

    bool proj2 = false;
    doMLMGNodalProjection(c_lev, f_lev-c_lev+1, vel, phi, sig,
                          amrex::GetVecOfPtrs(rhcc),
                          {}, proj_tol, proj_abs_tol, proj2, 0, 0);

    //
    // Unscale initial sync projection variables.
    //
    for (lev = c_lev; lev <= f_lev; lev++)
        rescaleVar(INITIAL_SYNC,sig[lev],1,vel[lev],lev);

    //
    // Add correction at coarse and fine levels.
    //
    for (lev = c_lev; lev <= f_lev; lev++)
    {
        MultiFab& P_new = LevelData[lev]->get_new_data(Press_Type);
        MultiFab::Add(P_new, *phi[lev], 0, 0, 1, 1);
    }

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
        BL_ASSERT(grids[Uoldmfi.index()].contains(Uoldmfi.tilebox())==true);

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

    vel_to_accel(lo, hi,
                      unew, ARLIM(un_lo), ARLIM(un_hi),
                      uold, ARLIM(uo_lo), ARLIM(uo_hi), &alpha );
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
    if (parent->Geom(0).IsRZ())
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
    if (parent->Geom(0).IsRZ())
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

    const Box& domain = parent->Geom(level).Domain();
    const int* domlo  = domain.loVect();
    const int* domhi  = domain.hiVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfmfi(mf,true); mfmfi.isValid(); ++mfmfi)
    {
      BL_ASSERT(mf.box(mfmfi.index()) == mfmfi.validbox());

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
#if (BL_SPACEDIM < 3)
    BL_ASSERT(radius_grow >= mf.nGrow());

    const Box& domain = parent->Geom(level).Domain();
    const int* domlo  = domain.loVect();
    const int* domhi  = domain.hiVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
      for (MFIter mfmfi(mf,true); mfmfi.isValid(); ++mfmfi)
       {
           BL_ASSERT(mf.box(mfmfi.index()) == mfmfi.validbox());

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
#if (BL_SPACEDIM < 3)
    BL_ASSERT(comp >= 0 && comp < mf.nComp());
    BL_ASSERT(radius_grow >= mf.nGrow());

    const Box& domain = parent->Geom(level).Domain();
    const int* domlo  = domain.loVect();
    const int* domhi  = domain.hiVect();

    Real bogus_value = BogusValue;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfmfi(mf,true); mfmfi.isValid(); ++mfmfi)
    {
        BL_ASSERT(mf.box(mfmfi.index()) == mfmfi.validbox());

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
// Multiply by anel_coeff if it is defined
//
void
Projection::AnelCoeffMult (int       level,
                           MultiFab& mf,
                           int       comp)
{
    BL_ASSERT(anel_coeff[level] != 0);
    BL_ASSERT(comp >= 0 && comp < mf.nComp());

    const Box& domain = parent->Geom(level).Domain();
    const int* domlo  = domain.loVect();
    const int* domhi  = domain.hiVect();

    Real bogus_value = BogusValue;

    int mult = 1;

    for (MFIter mfmfi(mf,true); mfmfi.isValid(); ++mfmfi)
    {
        const Box& bx = mfmfi.growntilebox();
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        Real* dat     = mf[mfmfi].dataPtr(comp);
	const int* datlo = mf[mfmfi].loVect();
	const int* dathi = mf[mfmfi].hiVect();
	//const int anel_len = std::size(anel_coeff[level][mfmfi.index()]);
	const Box& gbx = mfmfi.validbox();
	int anel_lo   = (gbx.loVect())[BL_SPACEDIM-1]-anel_grow;
	int anel_hi   = (gbx.hiVect())[BL_SPACEDIM-1]+anel_grow;

	anelcoeffmpy(lo,hi,dat,ARLIM(datlo),ARLIM(dathi),domlo,domhi,
		     anel_coeff[level][mfmfi.index()],&anel_lo,&anel_hi,
		     &bogus_value,&mult);
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

    const Box& domain = parent->Geom(level).Domain();
    const int* domlo  = domain.loVect();
    const int* domhi  = domain.hiVect();

    Real bogus_value = BogusValue;

    int mult = 0;

    for (MFIter mfmfi(mf,true); mfmfi.isValid(); ++mfmfi)
    {
        const Box& bx = mfmfi.growntilebox();
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        Real* dat     = mf[mfmfi].dataPtr(comp);
	const int* datlo = mf[mfmfi].loVect();
	const int* dathi = mf[mfmfi].hiVect();
	const Box& gbx = mfmfi.validbox();
	int anel_lo   = (gbx.loVect())[BL_SPACEDIM-1]-anel_grow;
	int anel_hi   = (gbx.hiVect())[BL_SPACEDIM-1]+anel_grow;

	anelcoeffmpy(lo,hi,dat,ARLIM(datlo),ARLIM(dathi),domlo,domhi,
		     anel_coeff[level][mfmfi.index()],&anel_lo,&anel_hi,
		     &bogus_value,&mult);

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
                                       BL_SPACEDIM, 1));
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

	  (*rhnd[lev])[mfi].setVal<RunOn::Host>(0,bx);
	  (*rhnd[lev])[mfi].copy<RunOn::Host>(P_new[mfi], bx, 0, bx, 0, 1);
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
        for (int n = 0; n < BL_SPACEDIM; n++)
	{
	  for (MFIter mfi(*vel[lev],true); mfi.isValid(); ++mfi)
            {
                const Box& box = mfi.tilebox();
                if (add_vort_proj)
                {
                  (*vel[lev])[mfi].plus<RunOn::Host>((*u_real[lev])[mfi],box,Xvel+n,Xvel+idx[n], 1);
                }
                else
                {
                  (*vel[lev])[mfi].copy<RunOn::Host>((*u_real[lev])[mfi],box,Xvel+n,box,Xvel+idx[n], 1);
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

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(phi_crse_strip); mfi.isValid(); ++mfi)
            {
                Box ovlp = amrex::coarsen(phi_fine_strip[iface].box(),ratio) & mfi.validbox();

                if (ovlp.ok())
                {
                    FArrayBox& cfab = phi_crse_strip[mfi];
                    fort_putdown (BL_TO_FORTRAN(cfab),
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
    gradp(p_dat,ARLIM(plo),ARLIM(phi),gp_dat,ARLIM(glo),ARLIM(ghi),lo,hi,dx,
               &is_full);
#elif (BL_SPACEDIM == 3)
    gradp(p_dat,ARLIM(plo),ARLIM(phi),gp_dat,ARLIM(glo),ARLIM(ghi),lo,hi,dx);
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
        phi_fine_strip[iface].setVal<RunOn::Host>(0.);
    }

    ProjOutFlowBC projBC;
    // These bcs just get passed into rhogbc() for all vals of which_call
    int        lo_bc[BL_SPACEDIM];
    int        hi_bc[BL_SPACEDIM];
    // change from phys_bcs of Inflow, SlipWall, etc.
    // to mathematical bcs of EXT_DIR, FOEXTRAP, etc.
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
      const int* lbc = phys_bc->lo();
      const int* hbc = phys_bc->hi();

      lo_bc[i]=scalar_bc[lbc[i]];
      hi_bc[i]=scalar_bc[hbc[i]];
    }
    if (which_call == INITIAL_PRESS)
    {
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
                dsdt[iface].setVal<RunOn::Host>(0.);
	}

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
	// FIXME: this size may need adjusting
	phi_fine_strip_ba.maxSize(32);
        DistributionMapping dm {phi_fine_strip_ba};
        MultiFab phi_fine_strip_mf(phi_fine_strip_ba,dm,1,0);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(phi_fine_strip_mf); mfi.isValid(); ++mfi) {
	    const Box& bx = mfi.validbox();
	    BL_ASSERT((phi_fine_strip[iface].box()).contains(bx));
	      phi_fine_strip_mf[mfi].copy<RunOn::Host>(phi_fine_strip[iface],bx,0,bx,0,1);
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
					bool proj2,
                                        MultiFab* sync_resid_crse,
                                        MultiFab* sync_resid_fine,
                                        bool doing_initial_vortproj)
{
    BL_PROFILE("Projection:::doMLMGNodalProjection()");

    int f_lev = c_lev + nlevel - 1;

    Vector<MultiFab> vel_test(nlevel);
    Vector<MultiFab> phi_test(nlevel);

    BL_ASSERT(vel[c_lev]->nGrow() >= 1);
    BL_ASSERT(vel[f_lev]->nGrow() >= 1);
    BL_ASSERT(phi[c_lev]->nGrow() == 1);
    BL_ASSERT(phi[f_lev]->nGrow() == 1);
    // MLMG does not copy any ghost cells from sig
    // BL_ASSERT(sig[c_lev]->nGrow() == 1);
    // BL_ASSERT(sig[f_lev]->nGrow() == 1);

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

    if (!rhcc.empty() )
    {
        AMREX_ALWAYS_ASSERT(rhcc[c_lev]->boxArray().ixType().cellCentered());
	// MLNodeLaplacian only uses vaild cells from rhcc and rhnd; fills ghost cells internally
        // BL_ASSERT(rhcc[c_lev]->nGrow() == 1);
        // BL_ASSERT(rhcc[f_lev]->nGrow() == 1);
    }

    if (!rhnd.empty() )
    {
        AMREX_ALWAYS_ASSERT(rhnd[c_lev]->boxArray().ixType().nodeCentered());
        // Do we need these two checks ??? -- no, see above
        // BL_ASSERT(rhnd[c_lev]->nGrow() == 1);
        // BL_ASSERT(rhnd[f_lev]->nGrow() == 1);
    }

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

    nodal_projector.project(phi_rebase,rel_tol,abs_tol);
    
#ifdef AMREX_USE_EB
        Vector< NavierStokesBase* > ns(nlevel);
        Vector< MultiFab* > Gp(nlevel);
        const auto gradphi = nodal_projector.getGradPhi();

        for (int lev = 0; lev < nlevel; lev++)
        {
            ns[lev] = dynamic_cast<NavierStokesBase*>(LevelData[lev+c_lev]);
            //fixme is this assert needed?
            BL_ASSERT(!(ns[lev]==0));
            Gp[lev] = &(ns[lev]->getGradP());

            // Do we need ghost cells here?
            if ( proj2 )
            {
                MultiFab::Copy(*Gp[lev],*gradphi[lev], 0, 0, AMREX_SPACEDIM,
                               gradphi[lev]->nGrow());
            }
            else
            {
                MultiFab::Add(*Gp[lev],*gradphi[lev], 0, 0, AMREX_SPACEDIM,
                              gradphi[lev]->nGrow());
            }
        }
#endif

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

    for (int idir=0; idir<BL_SPACEDIM; idir++) {

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

              for (int odir = 0; odir < BL_SPACEDIM; odir++) {
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

	  if (hi_bc[idir] == Inflow && reg.bigEnd(idir) == domainBox.bigEnd(idir)) {
	    Box bx;                // bx is the region we *protect* from zero'ing
	    bx = amrex::adjCellHi(reg, idir);

	    if (inflowCorner) {

              for (int odir = 0; odir < BL_SPACEDIM; odir++) {
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
              v_fab.setVal<RunOn::Host>(0.0, ovlp, Xvel+idir, 1);
            }
	  }
	}
      }
    }

  }
}
