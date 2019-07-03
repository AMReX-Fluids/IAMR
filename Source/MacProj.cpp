

#include <AMReX_LO_BCTYPES.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <MacProj.H>
#include <AMReX_MacBndry.H>
#include <MacOpMacDrivers.H>
#include <NavierStokesBase.H>
#include <MACPROJ_F.H>
#include <MacOutFlowBC.H>

using namespace amrex;

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

#define DEF_BOX_LIMITS(box,boxlo,boxhi)      \
const int* boxlo = (box).loVect();           \
const int* boxhi = (box).hiVect();

namespace
{
    bool initialized = false;
}
//
// Set defaults for these in Initialize()!!!
//
int  MacProj::verbose;
Real MacProj::mac_tol;
Real MacProj::mac_abs_tol;
Real MacProj::mac_sync_tol;
bool MacProj::use_cg_solve;
int  MacProj::do_outflow_bcs;
int  MacProj::fix_mac_sync_rhs;
int  MacProj::check_umac_periodicity;
int  MacProj::anel_grow       = 1;

namespace
{
    bool benchmarking;
    Real umac_periodic_test_Tol;
}

void
MacProj::Initialize ()
{
    if (initialized) return;
    //
    // Set defaults here!!!
    //
    benchmarking                    = false;
    umac_periodic_test_Tol          = 1.e-10;
    MacProj::verbose                = 0;
    MacProj::mac_tol                = 1.0e-12;
    MacProj::mac_abs_tol            = 1.0e-16;
    MacProj::mac_sync_tol           = 1.0e-8;
    MacProj::use_cg_solve           = false;
    MacProj::do_outflow_bcs         = 1;
    MacProj::fix_mac_sync_rhs       = 0;
    //
    // Only check umac periodicity when debugging.  Can be overridden on input.
    //
#ifndef AMREX_DEBUG
    MacProj::check_umac_periodicity = 0;
#else
    MacProj::check_umac_periodicity = 1;
#endif

    ParmParse pp("mac");

    pp.query("v",                      verbose);
    pp.query("mac_tol",                mac_tol);
    pp.query("mac_abs_tol",            mac_abs_tol);
    pp.query("mac_sync_tol",           mac_sync_tol);
    pp.query("use_cg_solve",           use_cg_solve);
    pp.query("benchmarking",           benchmarking);
    pp.query("do_outflow_bcs",         do_outflow_bcs);
    pp.query("fix_mac_sync_rhs",       fix_mac_sync_rhs);
    pp.query("check_umac_periodicity", check_umac_periodicity);
    pp.query("umac_periodic_test_Tol", umac_periodic_test_Tol);

    amrex::ExecOnFinalize(MacProj::Finalize);

    initialized = true;
}

void
MacProj::Finalize ()
{
    initialized = false;
}

//
// Setup functions follow
//

MacProj::MacProj (Amr*   _parent,
                  int    _finest_level,
                  BCRec* _phys_bc,
                  int    /*not used*/)
  :
    parent(_parent),
    LevelData(_finest_level+1),
    phys_bc(_phys_bc), 
    phi_bcs(_finest_level+1),
    mac_phi_crse(_finest_level+1),
    mac_reg(_finest_level+1),
    anel_coeff(_finest_level+1),
    finest_level(_finest_level)
{
    Initialize();

    if (verbose) amrex::Print() << "Creating mac_projector\n";

    finest_level_allocated = finest_level;

    for (int lev = 0; lev <= finest_level; lev++)
       anel_coeff[lev] = 0;
}

MacProj::~MacProj () {}


void
MacProj::install_level (int       level,
                        AmrLevel* level_data)
{
    if (verbose) amrex::Print() << "Installing MacProj level " << level << '\n';
    if (parent->finestLevel() < finest_level)
        for (int lev = parent->finestLevel() + 1; lev <= finest_level; lev++)
            mac_reg[lev].reset();

    finest_level = parent->finestLevel();

    if (level > finest_level_allocated)
    {
        finest_level_allocated = finest_level;
        LevelData.resize(finest_level+1);
        phi_bcs.resize(finest_level+1);
        mac_phi_crse.resize(finest_level+1);
        mac_reg.resize(finest_level+1);
    }

    LevelData[level] = level_data;

    BuildPhiBC(level);

    if (level > 0)
    {
        mac_reg[level].reset(new FluxRegister(LevelData[level]->boxArray(),
                                              LevelData[level]->DistributionMap(),
                                              parent->refRatio(level-1),level,1));
    }

    if (level > anel_coeff.size()-1) {
       anel_coeff.resize(level+1);
       anel_coeff[level] = 0;
    }
}
void
MacProj::install_anelastic_coefficient (int               level,
                                        Real**            _anel_coeff)
{
  if (verbose) amrex::Print() << "Installing anel_coeff into MacProj level " << level << '\n';

    if (level > anel_coeff.size()-1) anel_coeff.resize(level+1);
    anel_coeff[level] = _anel_coeff;
}
void
MacProj::build_anelastic_coefficient (int      level,
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
    // dependent) build_coefficient() function in problem directory 
    // ...Perhaps also need to worry about deleting anel_coeff
    // Also not sure why Projection and MacProj have separate anel_coeff
    // arrays, since they both appear to be cell centered.
    for (int j=0; j<len; j++)
      _anel_coeff[i][j] = 0.05*(jlo+j);    
  }
}

// xxxxx Can we skip this if using mlmg?
void
MacProj::BuildPhiBC (int level)
{
    const BoxArray& grids   = LevelData[level]->boxArray();
    const Geometry& geom    = parent->Geom(level);
    const int       ngrds   = grids.size();
    phi_bcs[level].resize(ngrds);
    const Box&      domain  = geom.Domain();
    const int*      domlo   = domain.loVect();
    const int*      domhi   = domain.hiVect();
    const int*      phys_lo = phys_bc->lo();
    const int*      phys_hi = phys_bc->hi();

    for (int i = 0; i < ngrds; i++)
    {
        BCRec&     bc = phi_bcs[level][i];
	const Box& grdbx = grids[i];
        const int* lo = grdbx.loVect();
        const int* hi = grdbx.hiVect();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (lo[dir] == domlo[dir])
            {
                bc.setLo(dir,phys_lo[dir]==Outflow ? LO_DIRICHLET : LO_NEUMANN);
            }
            else
            {
                bc.setLo(dir,LO_DIRICHLET);
            }
            if (hi[dir] == domhi[dir])
            {
                bc.setHi(dir,phys_hi[dir]==Outflow ? LO_DIRICHLET : LO_NEUMANN);
            }
            else
            {
                bc.setHi(dir,LO_DIRICHLET);
            }
        }
    }
}

void
MacProj::setup (int level)
{
    if (level < parent->maxLevel())
    {
        if (mac_phi_crse[level] == nullptr)
        {
            const BoxArray& grids = LevelData[level]->boxArray();
            const DistributionMapping& dmap = LevelData[level]->DistributionMap();
            mac_phi_crse[level].reset(new MultiFab(grids,dmap,1,1));
            mac_phi_crse[level]->setVal(0.0);
        }
    }
}

void
MacProj::cleanup (int level)
{
    if (level < parent->maxLevel())
        mac_phi_crse[level].reset();
}

//
// Projection functions follow ...
//
static
bool
grids_on_side_of_domain (const BoxArray&    grids,
                         const Box&         domain,
                         const Orientation& outFace)
{
    const int idir = outFace.coordDir();

    if (outFace.isLow())
    {
        for (int igrid = 0; igrid < grids.size(); igrid++)
        { 
            if (grids[igrid].smallEnd(idir) == domain.smallEnd(idir))
            { 
                return true;
            }
        }
    }
  
    if (outFace.isHigh())
    {
        for (int igrid = 0; igrid < grids.size(); igrid++)
        {
            if (grids[igrid].bigEnd(idir) == domain.bigEnd(idir))
            {
                return true;
            }
        }
    }

    return false;
}

//
// Compute the level advance mac projection.
//

void
MacProj::mac_project (int             level,
                      MultiFab*       u_mac,
                      MultiFab&       S,
                      Real            dt,
                      Real            time,
                      const MultiFab& divu,
                      int             have_divu,
                      bool            increment_vel_register)
{
    BL_PROFILE("MacProj::mac_project()");
    if (verbose) amrex::Print() << "... mac_project at level " << level << '\n';

    const BoxArray& grids      = LevelData[level]->boxArray();
    const DistributionMapping& dmap = LevelData[level]->DistributionMap();
    const Geometry& geom       = parent->Geom(level);
    const Real*     dx         = geom.CellSize();
    const int       max_level  = parent->maxLevel();
    MultiFab*       mac_phi    = 0;
    NavierStokesBase&   ns         = *(NavierStokesBase*) &(parent->getLevel(level));
    const MultiFab& volume     = ns.Volume();
    const MultiFab* area_level = ns.Area();
    IntVect         crse_ratio = level > 0 ? parent->refRatio(level-1)
                                           : IntVect::TheZeroVector();
    //
    // If finest level possible no need to make permanent mac_phi for bcs.
    //
    std::unique_ptr<MultiFab> raii;
    if (level == max_level) {
        raii.reset(new MultiFab(grids,dmap,1,1));
        mac_phi = raii.get();
    } else {
        mac_phi = mac_phi_crse[level].get();
    }

    mac_phi->setVal(0.0);
    //
    // HACK!!!
    //
    // Some of the routines we call assume that density has one valid
    // ghost cell.  We enforce that assumption by setting it here.
    //
    const MultiFab& rhotime = ns.get_rho(time);
    MultiFab::Copy(S, rhotime, 0, Density, 1, 1);

    if (OutFlowBC::HasOutFlowBC(phys_bc) && have_divu && do_outflow_bcs) {
        set_outflow_bcs(level, mac_phi, u_mac, S, divu);
    }

    //
    // Compute the nondivergent velocities, by creating the linop
    // and multigrid operator appropriate for the solved system.
    //
    // Initialize the rhs with divu.
    //
    const Real rhs_scale = 2.0/dt;
    MultiFab Rhs(grids,dmap,1,0);

    Rhs.copy(divu);

    MultiFab area_tmp[BL_SPACEDIM];
    if (anel_coeff[level] != 0) {
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    area_tmp[i].define(area_level[i].boxArray(), area_level[i].DistributionMap(), 1, 1);
	    MultiFab::Copy(area_tmp[i], area_level[i], 0, 0, 1, 1);
	}
        scaleArea(level,area_tmp,anel_coeff[level]);
    } 

    const MultiFab* area = (anel_coeff[level] != 0) ? area_tmp : area_level;

    MultiFab* cphi = (level == 0) ? nullptr : mac_phi_crse[level-1].get();
    mlmg_mac_level_solve(parent, cphi, *phys_bc, level, Density, mac_tol, mac_abs_tol,
                         rhs_scale, area, volume, S, Rhs, u_mac, mac_phi, verbose);

    Rhs.clear();
    //
    // Test that u_mac is divergence free
    //
    if (verbose)
        check_div_cond(level, u_mac);

    if (increment_vel_register)
    {
        //
        // Store advection velocities in mac registers at crse/fine boundaries.
        //
        // Initialize advection velocity registers with coarse grid velocity.
        //
        if (level < finest_level)
        {
            FluxRegister& mr = *mac_reg[level+1];

            mr.setVal(0.0);

            for (int dir = 0; dir < BL_SPACEDIM; dir++)
            {
                mr.CrseInit(u_mac[dir],area[dir],dir,0,0,1,-1.0);
            }

            if (verbose)
            {
                Real sumreg = mr.SumReg(0);

		amrex::Print() << "LEVEL " << level << " MACREG: CrseInit sum = " << sumreg << std::endl;
            }
        }
        //
        // Increment in fine grid velocity to velocity registers.
        //
        if (level > 0)
        {
            const Real mult = 1.0/parent->nCycle(level);

            for (int dir = 0; dir < BL_SPACEDIM; dir++)
            {
                mac_reg[level]->FineAdd(u_mac[dir],area[dir],dir,0,0,1,mult);
            }

            if (verbose)
            {
                Real sumreg = mac_reg[level]->SumReg(0);
		amrex::Print() << "LEVEL "                  << level
			       << " MACREG: FineAdd sum = " << sumreg << std::endl;
            }
        }
    }

    if (check_umac_periodicity)
        test_umac_periodic(level,u_mac);
}

//
// Compute the corrective pressure used in the mac_sync.
//
void
MacProj::mac_sync_solve (int       level,
                         Real      dt,
                         MultiFab& rho_half,
                         IntVect&  fine_ratio)
{
    BL_ASSERT(level < finest_level);

    if (verbose) amrex::Print() << "... mac_sync_solve at level " << level << '\n';

    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real      strt_time  = ParallelDescriptor::second();
    const BoxArray& grids      = LevelData[level]->boxArray();
    const DistributionMapping& dmap = LevelData[level]->DistributionMap();
    const Geometry& geom       = parent->Geom(level);
    const Real*     dx         = geom.CellSize();
    const BoxArray& fine_boxes = LevelData[level+1]->boxArray();
    IntVect         crse_ratio = level > 0 ? parent->refRatio(level-1)
                                           : IntVect::TheZeroVector();
    const NavierStokesBase& ns_level   = *(NavierStokesBase*) &(parent->getLevel(level));
    const MultiFab&     volume     = ns_level.Volume();
    const MultiFab*     area_level = ns_level.Area();
    //
    // Reusing storage here, since there should be no more need for the
    // values in mac_phi at this level and mac_sync_phi only need to last
    // into the call to mac_sync_compute.  Hope this works...  (LHH).
    //
    MultiFab* mac_sync_phi = mac_phi_crse[level].get();

    //
    // Solve the sync system.
    //
    //
    // Alloc and define RHS by doing a reflux-like operation in coarse
    // grid cells adjacent to fine grids.  The values in these
    // cells should be SUM{MR/VOL} where the sum is taken over
    // all edges of a cell that adjoin fine grids, MR = value in
    // MAC register, VOL = cell volume.  All other cells have a
    // value of zero (including crse cells under fine grids).
    //
    MultiFab Rhs(grids,dmap,1,0);
    Rhs.setVal(0.0);
    //
    // Reflux subtracts values at hi edge of coarse cell and
    // adds values at lo edge.  We want the opposite here so
    // set scale to -1 & alloc space for Rhs.
    //
    FluxRegister& mr = *mac_reg[level+1];
    const Real scale = -1.0;

    mr.Reflux(Rhs,volume,scale,0,0,1,geom);

    BoxArray baf = fine_boxes;

    baf.coarsen(fine_ratio);

#ifdef _OPENMP
#pragma omp parallel
#endif
    // fixme? Should do some real tests to see if tiling here is a win or not
    for (MFIter Rhsmfi(Rhs,true); Rhsmfi.isValid(); ++Rhsmfi)
    {
      BL_ASSERT(grids[Rhsmfi.index()].contains(Rhsmfi.tilebox()) );
      
	const std::vector< std::pair<int,Box> >& isects = baf.intersections(Rhsmfi.tilebox());

        for (int ii = 0, N = isects.size(); ii < N; ii++)
        {
            Rhs[Rhsmfi].setVal(0.0,isects[ii].second,0);
        }
    }

    mac_sync_phi->setVal(0.0);

    //
    // Now define edge centered coefficients and adjust RHS for MAC solve.
    //
    const Real rhs_scale = 2.0/dt;

    MultiFab area_tmp[BL_SPACEDIM];
    if (anel_coeff[level] != 0) {
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    area_tmp[i].define(area_level[i].boxArray(), area_level[i].DistributionMap(), 1, 1);
	    MultiFab::Copy(area_tmp[i], area_level[i], 0, 0, 1, 1);
	}
        scaleArea(level,area_tmp,anel_coeff[level]);
    } 

    const MultiFab* area = (anel_coeff[level] != 0) ? area_tmp : area_level;

    mlmg_mac_sync_solve(parent,*phys_bc, level, mac_sync_tol, mac_abs_tol,
                        rhs_scale, area, volume, rho_half, Rhs, mac_sync_phi, verbose);

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	amrex::Print() << "MacProj::mac_sync_solve(): time: " << run_time << std::endl;
    }
}

// this version is for the closed chamber LMC algorithm
void
MacProj::mac_sync_solve (int       level,
                         Real      dt,
                         MultiFab& rho_half,
                         IntVect&  fine_ratio,
			 MultiFab* Rhs_increment,
			 bool      subtract_avg,
			 Real&     offset)
{
    BL_ASSERT(level < finest_level);

    if (verbose) amrex::Print() << "... mac_sync_solve at level " << level << '\n';

    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real      strt_time  = ParallelDescriptor::second();
    const BoxArray& grids      = LevelData[level]->boxArray();
    const DistributionMapping& dmap = LevelData[level]->DistributionMap();
    const Geometry& geom       = parent->Geom(level);
    const Real*     dx         = geom.CellSize();
    const BoxArray& fine_boxes = LevelData[level+1]->boxArray();
    IntVect         crse_ratio = level > 0 ? parent->refRatio(level-1)
                                           : IntVect::TheZeroVector();
    const NavierStokesBase& ns_level   = *(NavierStokesBase*) &(parent->getLevel(level));
    const MultiFab&     volume     = ns_level.Volume();
    const MultiFab*     area_level = ns_level.Area();
    //
    // Reusing storage here, since there should be no more need for the
    // values in mac_phi at this level and mac_sync_phi only need to last
    // into the call to mac_sync_compute.  Hope this works...  (LHH).
    //
    MultiFab* mac_sync_phi = mac_phi_crse[level].get();

    //
    // Solve the sync system.
    //
    //
    // Alloc and define RHS by doing a reflux-like operation in coarse
    // grid cells adjacent to fine grids.  The values in these
    // cells should be SUM{MR/VOL} where the sum is taken over
    // all edges of a cell that adjoin fine grids, MR = value in
    // MAC register, VOL = cell volume.  All other cells have a
    // value of zero (including crse cells under fine grids).
    //
    MultiFab Rhs(grids,dmap,1,0);
    Rhs.setVal(0.0);
    //
    // Reflux subtracts values at hi edge of coarse cell and
    // adds values at lo edge.  We want the opposite here so
    // set scale to -1 & alloc space for Rhs.
    //
    FluxRegister& mr = *mac_reg[level+1];
    const Real scale = -1.0;

    mr.Reflux(Rhs,volume,scale,0,0,1,geom);

    BoxArray baf = fine_boxes;

    baf.coarsen(fine_ratio);

    // Use tiling here?
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter Rhsmfi(Rhs,true); Rhsmfi.isValid(); ++Rhsmfi)
    {
        BL_ASSERT(grids[Rhsmfi.index()] == Rhsmfi.validbox());

        const std::vector< std::pair<int,Box> >& isects = baf.intersections(Rhsmfi.validbox());

        FArrayBox& rhsfab = Rhs[Rhsmfi];

        for (int ii = 0, N = isects.size(); ii < N; ii++)
        {
            rhsfab.setVal(0.0,isects[ii].second,0);
        }
    }

    if (Rhs_increment)
    {
      MultiFab::Add(Rhs,*Rhs_increment,0,0,1,0);
    }


    mac_sync_phi->setVal(0.0);

    //
    // Now define edge centered coefficients and adjust RHS for MAC solve.
    //
    const Real rhs_scale = 2.0/dt;

    MultiFab area_tmp[BL_SPACEDIM];
    if (anel_coeff[level] != 0) {
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    area_tmp[i].define(area_level[i].boxArray(), area_level[i].DistributionMap(), 1, 1);
	    MultiFab::Copy(area_tmp[i], area_level[i], 0, 0, 1, 1);
	}
        scaleArea(level,area_tmp,anel_coeff[level]);
    } 

    const MultiFab* area = (anel_coeff[level] != 0) ? area_tmp : area_level;

    mlmg_mac_sync_solve(parent,*phys_bc, level, mac_sync_tol, mac_abs_tol,
                        rhs_scale, area, volume, rho_half, Rhs, mac_sync_phi, verbose);


    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	amrex::Print() << "MacProj::mac_sync_solve(): time: " << run_time << std::endl;
    }
}

//
// After solving for mac_sync_phi in mac_sync_solve(), we
// can now do the sync advect step.  This consists of two steps
//
// 1. compute u_corr as the gradient of mac_sync_phi
// 2. compute advective tendency of u_corr and
//    add into Vsync or Ssync
//
// If increment_sync is non-null, the (i-BL_SPACEDIM)-th component 
// of (*Ssync) is incremented only when increment[i]==1
// This is useful if that component gets incrmnted in a non-standard way.
//

void
MacProj::mac_sync_compute (int                   level,
                           MultiFab*             u_mac, 
                           MultiFab&             Vsync,
                           MultiFab&             Ssync,
                           MultiFab&             rho_half,
                           FluxRegister*         adv_flux_reg,
                           Vector<AdvectionForm>& advectionType,
                           Real                  prev_time, 
                           Real                  prev_pres_time,
                           Real                  dt, 
                           int                   NUM_STATE,
                           Real                  be_cn_theta,
                           bool                  modify_reflux_normal_vel,
                           int                   do_mom_diff,
                           const Vector<int>&     increment_sync,
			   bool                  update_fluxreg)
{
    if (modify_reflux_normal_vel)
        amrex::Abort("modify_reflux_normal_vel is no longer supported");
    //
    // Get parameters.
    //
    const BoxArray& grids               = LevelData[level]->boxArray();
    const DistributionMapping& dmap     = LevelData[level]->DistributionMap();
    const Geometry& geom                = parent->Geom(level);
    const Real*     dx                  = geom.CellSize();
    const int       numscal             = NUM_STATE - BL_SPACEDIM;
    MultiFab*       mac_sync_phi        = mac_phi_crse[level].get();
    NavierStokesBase&   ns_level        = *(NavierStokesBase*) &(parent->getLevel(level));
    const MultiFab& volume              = ns_level.Volume();
    const MultiFab* area                = ns_level.Area();
    Godunov*        godunov             = ns_level.godunov;
    bool            use_forces_in_trans = godunov->useForcesInTrans() ? true : false;

    MultiFab vel_visc_terms(grids,dmap,BL_SPACEDIM,1);
    MultiFab scal_visc_terms(grids,dmap,numscal,1);

    vel_visc_terms.setVal(0,1);  // Initialize to make calls below safe
    scal_visc_terms.setVal(0,1); // Initialize to make calls below safe
    //
    // Get viscous forcing.
    //
    if (be_cn_theta != 1.0) 
    {
        bool do_get_visc_terms = false;

        for (int i=0; i < BL_SPACEDIM; ++i)
            if (increment_sync.empty() || increment_sync[i]==1)
                do_get_visc_terms = true;

        if (do_get_visc_terms || use_forces_in_trans)
            ns_level.getViscTerms(vel_visc_terms,Xvel,BL_SPACEDIM,prev_time);

        do_get_visc_terms = false;
        for (int i=BL_SPACEDIM; i < increment_sync.size(); ++i)
            if (increment_sync.empty() || increment_sync[i]==1)
                do_get_visc_terms = true;

        if (do_get_visc_terms)
            ns_level.getViscTerms(scal_visc_terms,BL_SPACEDIM,numscal,prev_time);
    }

    MultiFab Gp(grids,dmap,BL_SPACEDIM,1);

    ns_level.getGradP(Gp, prev_pres_time);

    std::unique_ptr<MultiFab> divu_fp (ns_level.getDivCond(1,prev_time));
    //
    // Compute the mac sync correction.
    //
    MultiFab fluxes[BL_SPACEDIM];
    MultiFab mac_fluxes[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) {
      const BoxArray& ba = LevelData[level]->getEdgeBoxArray(i);
      fluxes[i].define(ba, dmap, NUM_STATE, 0);
      mac_fluxes[i].define(ba, dmap, 1, 0);
    }
        
    FillPatchIterator S_fpi(ns_level,vel_visc_terms,Godunov::hypgrow(),
                                 prev_time,State_Type,0,NUM_STATE);
    MultiFab& Smf = S_fpi.get_mf();
#ifdef _OPENMP
#pragma omp parallel 
#endif
    {
      Vector<int> ns_level_bc;
      FArrayBox tforces, tvelforces, U;
      FArrayBox grad_phi[BL_SPACEDIM], Rho;
      FArrayBox flux[BL_SPACEDIM];

      for (MFIter Smfi(Smf,true); Smfi.isValid(); ++Smfi)
      {
        const int i     = Smfi.index();
        FArrayBox& S    = Smf[Smfi];
        FArrayBox& divu = (*divu_fp)[Smfi];
        const Box bx = Smfi.tilebox();
        //
        // Step 1: compute ucorr = grad(phi)/rhonph
        //
        // Create storage for corrective velocities.
        //
        Rho.resize(amrex::grow(bx,1),1);

        D_TERM(grad_phi[0].resize(amrex::surroundingNodes(bx,0),1);,
               grad_phi[1].resize(amrex::surroundingNodes(bx,1),1);,
	       grad_phi[2].resize(amrex::surroundingNodes(bx,2),1););

        mac_vel_update(1,D_DECL(grad_phi[0],grad_phi[1],grad_phi[2]),
                       (*mac_sync_phi)[Smfi], rho_half[Smfi],
                       0, bx, bx, level, dx, dt/2.0);
        //
        // Step 2: compute Mac correction by calling GODUNOV box
        //
        // Get needed data.
        //
        Rho.copy(S,Density,0,1);

        ns_level.getForce(tforces,bx,1,0,NUM_STATE,prev_time,Smf[Smfi],Smf[Smfi],Density);

        //
        // Compute total forcing terms.
        //
        godunov->Sum_tf_gp_visc(tforces, 0, vel_visc_terms[Smfi], 0, Gp[Smfi], 0, Rho, 0);
        godunov->Sum_tf_divu_visc(S, BL_SPACEDIM, tforces, BL_SPACEDIM, numscal,
                                  scal_visc_terms[Smfi], 0, divu, 0, Rho, 0, 1);
        if (use_forces_in_trans)
        {
            ns_level.getForce(tvelforces,bx,1,Xvel,BL_SPACEDIM,prev_time,Smf[Smfi],Smf[Smfi],Density);
	          godunov->Sum_tf_gp_visc(tvelforces,0,vel_visc_terms[Smfi],0,Gp[Smfi],0,Rho,0);
        }
        //
        // Get the sync FABS.
        //
        FArrayBox& u_sync = Vsync[Smfi];
        FArrayBox& s_sync = Ssync[Smfi];

	//Why copy here and not just use S?
	// SyncAdvect doesn't even use U
        //U.resize(bx,BL_SPACEDIM);
        //U.copy(Smf[Smfi],0,0,BL_SPACEDIM);
        D_TERM(FArrayBox& u_mac_fab0 = u_mac[0][Smfi];,
               FArrayBox& u_mac_fab1 = u_mac[1][Smfi];,
               FArrayBox& u_mac_fab2 = u_mac[2][Smfi];);
        //
        // Loop over state components and compute the sync advective component.
        //
        FArrayBox* Sp;
	Box gbx = grow(Smfi.tilebox(),Smf.nGrow());
	FArrayBox rhoS(gbx,BL_SPACEDIM);
	
        for (int comp = 0; comp < NUM_STATE; comp++)
        {
            if (increment_sync.empty() || increment_sync[comp]==1)
            {
                const int  sync_ind = comp < BL_SPACEDIM ? comp  : comp-BL_SPACEDIM;
                FArrayBox& temp     = comp < BL_SPACEDIM ? u_sync : s_sync;
                ns_level_bc         = ns_level.fetchBCArray(State_Type,bx,comp,1);

                int use_conserv_diff = (advectionType[comp] == Conservative) ? true : false;

                if (do_mom_diff == 1 && comp < BL_SPACEDIM)
                {
		  rhoS.copy(Smf[Smfi],gbx,comp,gbx,comp,1);
		  rhoS.mult(Smf[Smfi],gbx,gbx,Density,comp,1);
		  Sp = &rhoS;
		  tforces.mult(Smf[Smfi],tforces.box(),tforces.box(),Density,comp,1);
                }
		else
		{
		  Sp = &Smf[Smfi];
		}
                
                for (int d=0; d<BL_SPACEDIM; ++d){
                  const Box& ebx = amrex::surroundingNodes(bx,d);
                  flux[d].resize(ebx,BL_SPACEDIM+1);
                }

                godunov->SyncAdvect(bx, dx, dt, level, 
                                    area[0][i], u_mac_fab0, grad_phi[0], flux[0], 
                                    area[1][i], u_mac_fab1, grad_phi[1], flux[1],
#if (BL_SPACEDIM == 3)                            
                                    area[2][i], u_mac_fab2, grad_phi[2], flux[2],
#endif
                                    S, *Sp, tforces, divu, comp, temp, sync_ind,
                                    use_conserv_diff, comp,
                                    ns_level_bc.dataPtr(), FPU, volume[i]);
                //
                // NOTE: the signs here are opposite from VELGOD.
                // NOTE: fluxes expected to be in extensive form.
                //
                if (level > 0 && update_fluxreg)
                { 
                  for (int d = 0; d < BL_SPACEDIM; d++){
		    const Box& ebx = Smfi.nodaltilebox(d);
		    fluxes[d][Smfi].copy(flux[d],ebx,0,ebx,comp,1);
		  } 
                }
            }
        }
        //
        // Include grad_phi in the mac registers corresponding
        // to the next coarsest interface.
        //
        if (level > 0 && update_fluxreg)
        {
	  for (int d = 0; d < BL_SPACEDIM; d++){
	    const Box& ebx = Smfi.nodaltilebox(d);
	    mac_fluxes[d][Smfi].copy(grad_phi[d],ebx,0,ebx,0,1);
	  }                 
        }
        //
        // Multiply the sync term by dt -- now done in the calling routine.
        //
      }
    } //end OMP parallel region
    if (level > 0 && update_fluxreg){
      const Real mlt =  -1.0/( (double) parent->nCycle(level));
      for (int d = 0; d < BL_SPACEDIM; d++){
	for (int comp = 0; comp < NUM_STATE; comp++){
	  if (increment_sync.empty() || increment_sync[comp]==1){
	    adv_flux_reg->FineAdd(fluxes[d],d,comp,comp,1,-dt);
	  }
	}
	mac_reg[level]->FineAdd(mac_fluxes[d],area[d],d,0,0,1,mlt);
      }
    }
}

//
// This routine does a sync advect step for a single 
// scalar component. Unlike the preceding routine, the
// half-time edge states are passed in from the calling routine.
// This routine is useful when the edge states are computed
// in a physics-class-specific manner. (For example, as they are
// in the calculation of div rho U h = div U sum_l (rho Y)_l h_l(T)).
//
void
MacProj::mac_sync_compute (int                    level,
                           MultiFab&              Sync,
                           int                    comp,
                           int                    s_ind,
                           const MultiFab* const* sync_edges,
			   int                    eComp,
                           MultiFab&              rho_half,
                           FluxRegister*          adv_flux_reg,
                           Vector<AdvectionForm>&  advectionType, 
			   bool                   modify_reflux_normal_vel,
                           Real                   dt,
			   bool                   update_fluxreg)
{
    if (modify_reflux_normal_vel)
        amrex::Abort("modify_reflux_normal_vel is no longer supported");

    const DistributionMapping& dmap     = LevelData[level]->DistributionMap();
    const Geometry& geom         = parent->Geom(level);
    MultiFab*       mac_sync_phi = mac_phi_crse[level].get();
    const NavierStokesBase& ns_level = *(NavierStokesBase*) &(parent->getLevel(level));
    const MultiFab& volume       = ns_level.Volume();
    const MultiFab* area         = ns_level.Area();

    Godunov godunov(512);
    MultiFab fluxes[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) {
      const BoxArray& ba = LevelData[level]->getEdgeBoxArray(i);
      fluxes[i].define(ba, dmap, 1, 0);
    }

    //
    // Compute the mac sync correction.
    //
#ifdef _OPENMP
#pragma omp parallel 
#endif
    {
      FArrayBox flux[BL_SPACEDIM], grad_phi[BL_SPACEDIM];
    
      for (MFIter Syncmfi(Sync,true); Syncmfi.isValid(); ++Syncmfi)
      {
	const Box& bx = Syncmfi.tilebox();

        //
        // Step 1: compute ucorr = grad(phi)/rhonph
        //
        D_TERM(grad_phi[0].resize(amrex::surroundingNodes(bx,0),1);,
               grad_phi[1].resize(amrex::surroundingNodes(bx,1),1);,
               grad_phi[2].resize(amrex::surroundingNodes(bx,2),1););

        mac_vel_update(1,
                       D_DECL(grad_phi[0],grad_phi[1],grad_phi[2]),
                       (*mac_sync_phi)[Syncmfi],
                       rho_half[Syncmfi], 0,
                       bx, bx, level,
                       geom.CellSize(), dt/2.0);
        //
        // Step 2: compute Mac correction by advecting the edge states.
        //
        D_TERM(flux[0].resize(amrex::surroundingNodes(bx,0),1);,
               flux[1].resize(amrex::surroundingNodes(bx,1),1);,
               flux[2].resize(amrex::surroundingNodes(bx,2),1););

        D_TERM(flux[0].copy((*sync_edges[0])[Syncmfi],eComp,0,1);,
               flux[1].copy((*sync_edges[1])[Syncmfi],eComp,0,1);,
               flux[2].copy((*sync_edges[2])[Syncmfi],eComp,0,1););

        int use_conserv_diff = (advectionType[comp] == Conservative) ? true : false;

        godunov.ComputeSyncAofs(bx,
                                area[0][Syncmfi],
                                grad_phi[0],       flux[0],
                                
                                area[1][Syncmfi],
                                grad_phi[1],       flux[1],
#if (BL_SPACEDIM == 3)                            
                                area[2][Syncmfi],
                                grad_phi[2],       flux[2],
#endif
                                volume[Syncmfi], Sync[Syncmfi],
                                s_ind, use_conserv_diff);

        D_TERM(grad_phi[0].clear();, grad_phi[1].clear();, grad_phi[2].clear(););
        //
        // NOTE: the signs here are opposite from VELGOD.
        // NOTE: fluxes expected to be in extensive form.
        //
        if (level > 0 && update_fluxreg)
        {
	  for (int d = 0; d < BL_SPACEDIM; d++){
	    const Box& ebx = Syncmfi.nodaltilebox(d);
	    fluxes[d][Syncmfi].copy(flux[d],ebx,0,ebx,0,1);
	  } 
        }
      }
    }//end OMP parallel region

    if (level > 0 && update_fluxreg){
      for (int d = 0; d < BL_SPACEDIM; d++){ 
	adv_flux_reg->FineAdd(fluxes[d],d,0,comp,1,-dt);
      }  
    }    
}

//
// Check the mac divergence.
//
void
MacProj::check_div_cond (int      level,
                         MultiFab U_edge[]) const
{
    const NavierStokesBase& ns_level = *(NavierStokesBase*) &(parent->getLevel(level));
    const MultiFab& volume       = ns_level.Volume();
    const MultiFab* area         = ns_level.Area();

    Real sum = 0.0;

#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sum)
#endif
    {
      FArrayBox dmac;
      
      for (MFIter U_edge0mfi(U_edge[0],true);U_edge0mfi.isValid();++U_edge0mfi)
      {
        const Box& bx = U_edge0mfi.tilebox(IntVect::Zero);

        dmac.resize(bx,1);

        const FArrayBox& uxedge = U_edge[0][U_edge0mfi];
        const FArrayBox& uyedge = U_edge[1][U_edge0mfi];
        const FArrayBox& xarea  = area[0][U_edge0mfi];
        const FArrayBox& yarea  = area[1][U_edge0mfi];
        const FArrayBox& vol    = volume[U_edge0mfi];

        DEF_LIMITS(dmac,dmac_dat,dlo,dhi);
        DEF_CLIMITS(uxedge,ux_dat,uxlo,uxhi);
        DEF_CLIMITS(uyedge,uy_dat,uylo,uyhi);
        DEF_CLIMITS(xarea,ax_dat,axlo,axhi);
        DEF_CLIMITS(yarea,ay_dat,aylo,ayhi);
        DEF_CLIMITS(vol,vol_dat,vlo,vhi);

#if (BL_SPACEDIM == 2)
        macdiv(dmac_dat,ARLIM(dlo),ARLIM(dhi),dlo,dhi,
                    ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                    uy_dat,ARLIM(uylo),ARLIM(uyhi),
                    ax_dat,ARLIM(axlo),ARLIM(axhi), 
                    ay_dat,ARLIM(aylo),ARLIM(ayhi), 
                    vol_dat,ARLIM(vlo),ARLIM(vhi));
#endif

#if (BL_SPACEDIM == 3)
        const FArrayBox& uzedge = U_edge[2][U_edge0mfi];
        DEF_CLIMITS(uzedge,uz_dat,uzlo,uzhi);
        const FArrayBox& zarea = area[2][U_edge0mfi];
        DEF_CLIMITS(zarea,az_dat,azlo,azhi);

        macdiv(dmac_dat,ARLIM(dlo),ARLIM(dhi),dlo,dhi,
                    ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                    uy_dat,ARLIM(uylo),ARLIM(uyhi),
                    uz_dat,ARLIM(uzlo),ARLIM(uzhi),
                    ax_dat,ARLIM(axlo),ARLIM(axhi),
                    ay_dat,ARLIM(aylo),ARLIM(ayhi),
                    az_dat,ARLIM(azlo),ARLIM(azhi),
                    vol_dat,ARLIM(vlo),ARLIM(vhi));
#endif
        sum += dmac.sum(0);
      }
    }
    
    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceRealSum(sum,IOProc);

	amrex::Print().SetPrecision(15) << "SUM of DIV(U_edge) = " << sum << '\n';
    }
}

void
MacProj::set_outflow_bcs (int             level,
                          MultiFab*       mac_phi,
                          const MultiFab* u_mac, 
                          const MultiFab& S,
                          const MultiFab& divu)
{
    //
    // This code is very similar to the outflow BC stuff in the Projection
    // class except that here the the phi to be solved for lives on the
    // out-directed faces.  The projection equation to satisfy is
    //
    //   (1/r)(d/dr)[r/rho dphi/dr] = dv/dr - S
    //
    bool hasOutFlow;
    Orientation outFaces[2*BL_SPACEDIM];
    int numOutFlowFaces;

    OutFlowBC::GetOutFlowFaces(hasOutFlow,outFaces,phys_bc,numOutFlowFaces);

    const BoxArray&   grids  = LevelData[level]->boxArray();
    const Geometry&   geom   = parent->Geom(level);
    const Box&        domain = parent->Geom(level).Domain();
    //
    // Create 1-wide cc box just outside boundary to hold phi.
    // Create 1-wide cc box just inside  boundary to hold rho,u,divu.
    //
    BoxList ccBoxList, phiBoxList;
    // numOutFlowFaces gives the number of outflow faces on the entire
    //   problem domain
    // nOutFlowTouched gives the number of outflow faces a level touches, so
    //   nOutFlowTouched = numOutFlowFaces for level 0, but
    //   nOutFlowTouched <= numOutFlowFaces for levels > 0, since
    //   finer levels may not span the entire problem domain
    int nOutFlowTouched = 0;
    for (int iface = 0; iface < numOutFlowFaces; iface++)
    {
        if (grids_on_side_of_domain(grids,geom.Domain(),outFaces[iface])) 
	{
	    nOutFlowTouched++;
            const int outDir    = outFaces[iface].coordDir();

            Box ccBndBox;
            if (outFaces[iface].faceDir() == Orientation::high)
	    {
                ccBndBox = amrex::adjCellHi(domain,outDir,2);
                ccBndBox.shift(outDir,-2);
	    } 
            else 
	    {
                ccBndBox = amrex::adjCellLo(domain,outDir,2);
                ccBndBox.shift(outDir,2);
	    }
            ccBoxList.push_back(ccBndBox);

            Box phiBox  = amrex::adjCell(domain,outFaces[iface],1);
            phiBoxList.push_back(phiBox);

            const Box&     valid_ccBndBox       = ccBndBox & domain;
            const BoxArray uncovered_outflow_ba = amrex::complementIn(valid_ccBndBox,grids);

            if (uncovered_outflow_ba.size() && 
                amrex::intersect(grids,valid_ccBndBox).size())
                amrex::Error("MacProj: Cannot yet handle partially refined outflow");
	}
    }
  
    if ( !ccBoxList.isEmpty() ) 
    {
        BoxArray  ccBoxArray( ccBoxList);
        BoxArray phiBoxArray(phiBoxList);
        ccBoxList.clear();
        phiBoxList.clear();

        FArrayBox rhodat[2*BL_SPACEDIM];
        FArrayBox divudat[2*BL_SPACEDIM];
        FArrayBox phidat[2*BL_SPACEDIM];

	for ( int iface = 0; iface < nOutFlowTouched; ++iface) 
	{
            rhodat[iface].resize(ccBoxArray[iface], 1);
            divudat[iface].resize(ccBoxArray[iface], 1);
            phidat[iface].resize(phiBoxArray[iface], 1);

            phidat[iface].setVal(0.0);
            divu.copyTo(divudat[iface]);
            S.copyTo(rhodat[iface], Density, 0, 1);
	}

        // rhodat.copy(S, Density, 0, 1);
        // divudat.copy(divu, 0, 0, 1);

        //
        // Load ec data.
        //

        FArrayBox uedat[BL_SPACEDIM][2*BL_SPACEDIM];
        for (int i = 0; i < BL_SPACEDIM; ++i)
	{
            BoxArray edgeArray(ccBoxArray);
            edgeArray.surroundingNodes(i);
	    for ( int iface = 0; iface < nOutFlowTouched; ++iface) 
	    {
                uedat[i][iface].resize(edgeArray[iface], 1);
                u_mac[i].copyTo(uedat[i][iface], 0, 0, 1);
	    }
	}
    
        MacOutFlowBC macBC;

        NavierStokesBase* ns_level = dynamic_cast<NavierStokesBase*>(&parent->getLevel(level));
        Real gravity = ns_level->getGravity();
        const int* lo_bc = phys_bc->lo();
        const int* hi_bc = phys_bc->hi();
        macBC.computeBC(uedat, divudat, rhodat, phidat,
			geom, outFaces, nOutFlowTouched, lo_bc, hi_bc,
			umac_periodic_test_Tol, gravity);
        //
        // Must do this kind of copy instead of mac_phi->copy(phidat);
        // because we're copying onto the ghost cells of the FABs,
        // not the valid regions.
        //
#ifdef _OPENMP
#pragma omp parallel
#endif
	for ( int iface = 0; iface < nOutFlowTouched; ++iface )
	{
            for (MFIter mfi(*mac_phi); mfi.isValid(); ++mfi)
	    {
                Box ovlp = (*mac_phi)[mfi].box() & phidat[iface].box();
                if (ovlp.ok())
		  (*mac_phi)[mfi].copy(phidat[iface],ovlp,0,ovlp,0,1);
	    }
	}
    }
}

//
// Structure used by test_umac_periodic().
//

struct TURec
{
    TURec ()
        :
        m_idx(-1),
        m_dim(-1)
    {}

    TURec (int        idx,
           int        dim,
           const Box& srcBox,
           const Box& dstBox)
        :
        m_srcBox(srcBox),
        m_dstBox(dstBox),
        m_idx(idx),
        m_dim(dim)
    {}

    FillBoxId m_fbid;
    Box       m_srcBox;
    Box       m_dstBox;
    int       m_idx;
    int       m_dim;
};

//
// Test that edge-based values agree across periodic boundary.
//

void
MacProj::test_umac_periodic (int       level,
                             MultiFab* u_mac)
{
    const Geometry& geom = parent->Geom(level);

    if (!geom.isAnyPeriodic()) return;

    FArrayBox              diff;
    MultiFabCopyDescriptor mfcd;
    MultiFabId             mfid[BL_SPACEDIM];
    std::vector<TURec>     pirm;
    Vector<IntVect>         pshifts(27);
    std::vector< std::pair<int,Box> > isects;


    for (int dim = 0; dim < BL_SPACEDIM; dim++)
    {
        if (geom.isPeriodic(dim))
        {
            Box eDomain = amrex::surroundingNodes(geom.Domain(),dim);

            mfid[dim] = mfcd.RegisterMultiFab(&u_mac[dim]);

	    // How to combine pirm into one global pirm?
	    // don't think std::vector::push_back() is thread safe
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
// {
//             std::vector<TURec>     pirm;
//             std::vector< std::pair<int,Box> > isects;
//             Vector<IntVect>         pshifts(27);
    
            for (MFIter mfi(u_mac[dim]); mfi.isValid(); ++mfi)
            {
                Box eBox = u_mac[dim].boxArray()[mfi.index()];

                geom.periodicShift(eDomain, eBox, pshifts);

                for (int iiv = 0, M = pshifts.size(); iiv < M; iiv++)
                {
                    eBox += pshifts[iiv];

                    u_mac[dim].boxArray().intersections(eBox,isects);

                    for (int i = 0, N = isects.size(); i < N; i++)
                    {
                        const Box& srcBox = isects[i].second;
                        const Box& dstBox = srcBox - pshifts[iiv];

                        TURec r(mfi.index(),dim,srcBox,dstBox);

                        r.m_fbid = mfcd.AddBox(mfid[dim],
                                               srcBox,
                                               0,
                                               isects[i].first,
                                               0,
                                               0,
                                               1);
                        pirm.push_back(r);
                    }

                    eBox -= pshifts[iiv];
                }
            }
	    // }// end OMP region
        }
    }

    int nrecv = pirm.size();
    ParallelDescriptor::ReduceIntMax(nrecv);
    if (nrecv == 0)
        //
        // There's no parallel work to do.
        //
        return;

    mfcd.CollectData();

    for (int i = 0; i < pirm.size(); i++)
    {
        const int dim = pirm[i].m_dim;

        BL_ASSERT(pirm[i].m_fbid.box() == pirm[i].m_srcBox);
        BL_ASSERT(pirm[i].m_srcBox.sameSize(pirm[i].m_dstBox));
        BL_ASSERT(u_mac[dim].DistributionMap()[pirm[i].m_idx] == ParallelDescriptor::MyProc());

        diff.resize(pirm[i].m_srcBox, 1);

        mfcd.FillFab(mfid[dim], pirm[i].m_fbid, diff);

        diff.minus(u_mac[dim][pirm[i].m_idx],pirm[i].m_dstBox,diff.box(),0,0,1);

        const Real max_norm = diff.norm(0);

        if (max_norm > umac_periodic_test_Tol )
        {
	  amrex::Print() << "dir = "         << dim
			 << ", diff norm = " << max_norm
			 << " for region: "  << pirm[i].m_dstBox << std::endl;
	  amrex::Error("Periodic bust in u_mac");
        }
    }
}

void
MacProj::scaleArea (int level, MultiFab* area, Real** anel_coeff_loc)
{
    const BoxArray& grids = LevelData[level]->boxArray();

    int mult = 1;

    for (MFIter mfi(area[0],true); mfi.isValid(); ++mfi)
    {
        const FArrayBox& xarea  = area[0][mfi];
        const FArrayBox& yarea  = area[1][mfi];
        DEF_CLIMITS(xarea,ax_dat,axlo,axhi);
        DEF_CLIMITS(yarea,ay_dat,aylo,ayhi);

        const Box& bx = mfi.tilebox(IntVect::Zero);
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
	
        const int  i        = mfi.index();
	const Box& gbx      = grids[i];
        const int* gbxhi    = gbx.hiVect();
	int anel_coeff_lo   = (gbx.loVect())[BL_SPACEDIM-1]-anel_grow;
	int anel_coeff_hi   = (gbx.hiVect())[BL_SPACEDIM-1]+anel_grow;

	
#if (BL_SPACEDIM == 2)
        fort_scalearea(lo,hi,gbxhi,
		       ax_dat,ARLIM(axlo),ARLIM(axhi), 
                       ay_dat,ARLIM(aylo),ARLIM(ayhi), 
                       anel_coeff_loc[i],&anel_coeff_lo,&anel_coeff_hi,
		       &mult);

#elif (BL_SPACEDIM == 3)
        const FArrayBox& zarea = area[2][mfi];
        DEF_CLIMITS(zarea,az_dat,azlo,azhi);
        fort_scalearea(lo,hi,gbxhi,
		       ax_dat,ARLIM(axlo),ARLIM(axhi), 
                       ay_dat,ARLIM(aylo),ARLIM(ayhi), 
                       az_dat,ARLIM(azlo),ARLIM(azhi), 
                       anel_coeff_loc[i],&anel_coeff_lo,&anel_coeff_hi,
		       &mult);

#endif
    }
}

