
#include <winstd.H>

#include <LO_BCTYPES.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <MacProj.H>
#include <MacBndry.H>
#include <MacOpMacDrivers.H>
#include <NavierStokes.H>
#include <MACPROJ_F.H>
#include <MacOutFlowBC.H>

#ifndef _NavierStokes_H_
enum StateType {State_Type=0, Press_Type};
#if (BL_SPACEDIM == 2)
enum StateNames  { Xvel=0, Yvel, Density};
#else
enum StateNames  { Xvel=0, Yvel, Zvel, Density};
#endif
#endif

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

#define GEOM_GROW 1

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

namespace
{
    bool benchmarking;
    Real umac_periodic_test_Tol;

#if MG_USE_HYPRE
    bool use_hypre_solve;
#endif
    bool use_fboxlib_mg;
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
#ifdef NDEBUG
    MacProj::check_umac_periodicity = 0;
#else
    MacProj::check_umac_periodicity = 1;
#endif

#if MG_USE_HYPRE
    use_hypre_solve  = false;
#endif
    use_fboxlib_mg   = false;

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

#if MG_USE_HYPRE
    pp.query("use_hypre_solve", use_hypre_solve);
#endif
    pp.query("use_fboxlib_mg",  use_fboxlib_mg);

#if MG_USE_HYPRE
    if ( use_cg_solve && use_hypre_solve )
	BoxLib::Error("MacProj::read_params: cg_solve && .not. hypre_solve");
#endif
    if ( use_cg_solve && use_fboxlib_mg )
	BoxLib::Error("MacProj::read_params: cg_solve && .not. fboxlib_solve");

    BoxLib::ExecOnFinalize(MacProj::Finalize);

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
    mac_phi_crse(_finest_level+1, PArrayManage),
    mac_reg(_finest_level+1, PArrayManage),
    anel_coeff(_finest_level+1),
    finest_level(_finest_level)
{
    Initialize();

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Creating mac_projector\n";

    finest_level_allocated = finest_level;

    for (int lev = 0; lev <= finest_level; lev++)
       anel_coeff[lev] = 0;
}

MacProj::~MacProj () {}


void
MacProj::install_level (int       level,
                        AmrLevel* level_data)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Installing MacProj level " << level << '\n';

    if (parent->finestLevel() < finest_level)
        for (int lev = parent->finestLevel() + 1; lev <= finest_level; lev++)
            mac_reg.clear(lev);

    finest_level = parent->finestLevel();

    if (level > finest_level_allocated)
    {
        finest_level_allocated = finest_level;
        LevelData.resize(finest_level+1);
        phi_bcs.resize(finest_level+1);
        mac_phi_crse.resize(finest_level+1);
        mac_reg.resize(finest_level+1);
    }

    LevelData.clear(level);
    LevelData.set(level, level_data);

    BuildPhiBC(level);

    if (level > 0)
    {
        mac_reg.clear(level);
        mac_reg.set(level,new FluxRegister(LevelData[level].boxArray(),
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
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Installing anel_coeff into MacProj level " << level << '\n';

    if (level > anel_coeff.size()-1) anel_coeff.resize(level+1);
    anel_coeff.set(level, _anel_coeff);
}

void
MacProj::BuildPhiBC (int level)
{
    const BoxArray& grids   = LevelData[level].boxArray();
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
        const int* lo = grids[i].loVect();
        const int* hi = grids[i].hiVect();

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
        if (!mac_phi_crse.defined(level))
        {
            const BoxArray& grids = LevelData[level].boxArray();
            mac_phi_crse.set(level,new MultiFab(grids,1,1));
            mac_phi_crse[level].setVal(0.0);
        }
    }
}

void
MacProj::cleanup (int level)
{
    if (level < parent->maxLevel())
        mac_phi_crse.clear(level);
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
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... mac_project at level " << level << '\n';

    const BoxArray& grids      = LevelData[level].boxArray();
    const Geometry& geom       = parent->Geom(level);
    const Real*     dx         = geom.CellSize();
    const int       max_level  = parent->maxLevel();
    MultiFab*       mac_phi    = 0;
    NavierStokes&   ns         = *(NavierStokes*) &(parent->getLevel(level));
    IntVect         crse_ratio = level > 0 ? parent->refRatio(level-1)
                                           : IntVect::TheZeroVector();
    //
    // If finest level possible no need to make permanent mac_phi for bcs.
    //
    if (level == max_level)
        mac_phi = new MultiFab(grids,1,1);
    else
        mac_phi = &mac_phi_crse[level];

    mac_phi->setVal(0.0);
    //
    // HACK!!!
    //
    // Some of the routines we call assume that density has one valid
    // ghost cell.  We enforce that assumption by setting it here.
    //
    const MultiFab& rhotime = ns.get_rho(time);

    for (MFIter mfi(rhotime); mfi.isValid(); ++mfi)
        S[mfi].copy(rhotime[mfi],0,Density,1);

    if (OutFlowBC::HasOutFlowBC(phys_bc) && have_divu && do_outflow_bcs)
        set_outflow_bcs(level, mac_phi, u_mac, S, divu);
    //
    // Store the Dirichlet boundary condition for mac_phi in mac_bndry.
    //
    MacBndry mac_bndry(grids,1,geom);
    const int src_comp  = 0;
    const int dest_comp = 0;
    const int num_comp  = 1;
    if (level == 0)
    {
        mac_bndry.setBndryValues(*mac_phi,src_comp,dest_comp,num_comp,*phys_bc);
    }
    else
    {
        MultiFab& CPhi = mac_phi_crse[level-1];
	CPhi.FillBoundary();
	parent->Geom(level-1).FillPeriodicBoundary(CPhi,true);

        BoxArray crse_boxes(grids);
        crse_boxes.coarsen(crse_ratio);
        const int in_rad     = 0;
        const int out_rad    = 1;
        //const int extent_rad = 1;
        const int extent_rad = 2;
        BndryRegister crse_br(crse_boxes,in_rad,out_rad,extent_rad,num_comp);
        crse_br.copyFrom(CPhi,CPhi.nGrow(),src_comp,dest_comp,num_comp);

        mac_bndry.setBndryValues(crse_br,src_comp,*mac_phi,src_comp,
                                 dest_comp,num_comp,crse_ratio,*phys_bc);
    }
    //
    // Compute the nondivergent velocities, by creating the linop
    // and multigrid operator appropriate for the solved system.
    //
    // Initialize the rhs with divu.
    //
    const Real rhs_scale = 2.0/dt;
    MultiFab Rhs(grids,1,0), area[BL_SPACEDIM];

    Rhs.copy(divu);

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        geom.GetFaceArea(area[dir],grids,dir,GEOM_GROW);
    }

    int the_solver = 0;
    if (use_cg_solve)
    {
	the_solver = 1;
    }
#if MG_USE_HYPRE
    else if ( use_hypre_solve )
    {
	the_solver = 2;
    }
#endif
    else if ( use_fboxlib_mg )
    {
	the_solver = 3;
    }

    if (anel_coeff[level] != 0)
        scaleArea(level,area,anel_coeff[level]);

    {
        MultiFab volume;

        geom.GetVolume(volume,grids,GEOM_GROW);

        mac_level_driver(parent, mac_bndry, *phys_bc, grids, the_solver, level, Density,
                         dx, dt, mac_tol, mac_abs_tol, rhs_scale, 
                         area, volume, S, Rhs, u_mac, mac_phi, verbose);
    }

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
            FluxRegister& mr = mac_reg[level+1];

            mr.setVal(0.0);

            for (int dir = 0; dir < BL_SPACEDIM; dir++)
            {
                mr.CrseInit(u_mac[dir],area[dir],dir,0,0,1,-1.0);
            }

            if (verbose)
            {
                Real sumreg = mr.SumReg(0);

                if (ParallelDescriptor::IOProcessor())
                {
                    std::cout << "LEVEL " << level << " MACREG: CrseInit sum = " << sumreg << std::endl;
                }
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
                mac_reg[level].FineAdd(u_mac[dir],area[dir],dir,0,0,1,mult);
            }

            if (verbose)
            {
                Real sumreg = mac_reg[level].SumReg(0);

                if (ParallelDescriptor::IOProcessor())
                {
                    std::cout << "LEVEL "                  << level
                              << " MACREG: FineAdd sum = " << sumreg << std::endl;
                }
            }
        }
    }
    //
    // If finest level possible no need to keep phi for boundary conditions.
    //
    if (level == max_level)
    {
        delete mac_phi;
        mac_phi = 0;
    }

    if (check_umac_periodicity)
        test_umac_periodic(level,u_mac);

    if (anel_coeff[level] != 0) rescaleArea(level,area,anel_coeff[level]);
}

//
// Compute the corrective pressure used in the mac_sync.
//

void
MacProj::mac_sync_solve (int       level,
                         Real      dt,
                         MultiFab* rho_half,
                         IntVect&  fine_ratio)
{
    BL_ASSERT(level < finest_level);

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... mac_sync_solve at level " << level << '\n';

    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real      strt_time  = ParallelDescriptor::second();
    const BoxArray& grids      = LevelData[level].boxArray();
    const Geometry& geom       = parent->Geom(level);
    const Real*     dx         = geom.CellSize();
    const BoxArray& fine_boxes = LevelData[level+1].boxArray();
    IntVect         crse_ratio = level > 0 ? parent->refRatio(level-1)
                                           : IntVect::TheZeroVector();
    //
    // Reusing storage here, since there should be no more need for the
    // values in mac_phi at this level and mac_sync_phi only need to last
    // into the call to mac_sync_compute.  Hope this works...  (LHH).
    //
    MultiFab* mac_sync_phi = &mac_phi_crse[level];
    //
    // Alloc and define RHS by doing a reflux-like operation in coarse
    // grid cells adjacent to fine grids.  The values in these
    // cells should be SUM{MR/VOL} where the sum is taken over
    // all edges of a cell that adjoin fine grids, MR = value in
    // MAC register, VOL = cell volume.  All other cells have a
    // value of zero (including crse cells under fine grids).
    //
    MultiFab Rhs(grids,1,0), area[BL_SPACEDIM];
    Rhs.setVal(0.0);
    //
    // Reflux subtracts values at hi edge of coarse cell and
    // adds values at lo edge.  We want the opposite here so
    // set scale to -1 & alloc space for Rhs.
    //
    FluxRegister& mr = mac_reg[level+1];
    const Real scale = -1.0;

    MultiFab volume;

    geom.GetVolume(volume,grids,GEOM_GROW);

    mr.Reflux(Rhs,volume,scale,0,0,1,geom);

    BoxArray baf = fine_boxes;

    baf.coarsen(fine_ratio);

    std::vector< std::pair<int,Box> > isects;

    for (MFIter Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
    {
        BL_ASSERT(grids[Rhsmfi.index()] == Rhsmfi.validbox());

        baf.intersections(Rhsmfi.validbox(),isects);

        FArrayBox& rhsfab = Rhs[Rhsmfi];

        for (int ii = 0, N = isects.size(); ii < N; ii++)
        {
            rhsfab.setVal(0.0,isects[ii].second,0);
        }
    }

    baf.clear_hash_bin();
    //
    // Remove constant null space component from the rhs of the solve
    // when appropriate (i.e. when the grids span the whole domain AND
    // the boundary conditions are Neumann for phi on all sides of the domain.)
    // Note that Rhs does not yet have the radial scaling in it for r-z
    // problems so we must do explicit volume-weighting here.
    //
    if (fix_mac_sync_rhs)
    {
        int all_neumann = 1;
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (phys_bc->lo()[dir] == Outflow || phys_bc->hi()[dir] == Outflow)
                all_neumann = 0;
        }

        if (Rhs.boxArray().contains(geom.Domain()) && all_neumann == 1)
        {
            Real sum = 0.0;
            Real vol = 0.0;
            FArrayBox vol_wgted_rhs;
            for (MFIter Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
            {
                const FArrayBox& rhsfab = Rhs[Rhsmfi];
                vol_wgted_rhs.resize(rhsfab.box());
                vol_wgted_rhs.copy(rhsfab);
                vol_wgted_rhs.mult(volume[Rhsmfi]);
                sum += vol_wgted_rhs.sum(0,1);
                vol += volume[Rhsmfi].sum(rhsfab.box(),0,1);
            }

            Real vals[2] = {sum,vol};

            ParallelDescriptor::ReduceRealSum(&vals[0],2);

            sum = vals[0];
            vol = vals[1];

            const Real fix = sum / vol;

            if (verbose && ParallelDescriptor::IOProcessor())
                std::cout << "Average correction on mac sync RHS = " << fix << '\n';

            Rhs.plus(-fix, 0);
        }
    }

    mac_sync_phi->setVal(0.0);
    //
    // store the Dirichlet boundary condition for mac_sync_phi in mac_bndry
    //
    MacBndry mac_bndry(grids,1,geom);
    const int src_comp = 0;
    const int dest_comp = 0;
    const int num_comp = 1;
    if (level == 0)
    {
        mac_bndry.setBndryValues(*mac_sync_phi,src_comp,dest_comp,num_comp,
                                 *phys_bc);
    }
    else
    {
        BoxArray crse_boxes(grids);
        crse_boxes.coarsen(crse_ratio);
        const int in_rad     = 0;
        const int out_rad    = 1;
        //const int extent_rad = 1;
        const int extent_rad = 2;
        BndryRegister crse_br(crse_boxes,in_rad,out_rad,extent_rad,num_comp);
        crse_br.setVal(0);
        mac_bndry.setBndryValues(crse_br,src_comp,*mac_sync_phi,src_comp,
                                 dest_comp,num_comp,crse_ratio, *phys_bc);
    }
    //
    // Now define edge centered coefficients and adjust RHS for MAC solve.
    //
    const Real rhs_scale = 2.0/dt;
    //
    // Solve the sync system.
    //
    int the_solver = 0;
    if ( use_cg_solve )
    {
	the_solver = 1;
    }
#if MG_USE_HYPRE
    else if ( use_hypre_solve )
    {
	the_solver = 2;
    }
#endif
    else if ( use_fboxlib_mg )
    {
	the_solver = 3;
    }

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        geom.GetFaceArea(area[dir],grids,dir,GEOM_GROW);
    }
    if (anel_coeff[level] != 0)
        scaleArea(level,area,anel_coeff[level]);

    mac_sync_driver(parent, mac_bndry, *phys_bc, grids, the_solver, level, dx, dt,
                    mac_sync_tol, mac_abs_tol, rhs_scale, area,
                    volume, Rhs, rho_half, mac_sync_phi, verbose);

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "MacProj::mac_sync_solve(): time: " << run_time << std::endl;
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
                           MultiFab*             Vsync,
                           MultiFab*             Ssync,
                           MultiFab*             rho_half,
                           FluxRegister*         adv_flux_reg,
                           Array<AdvectionForm>& advectionType,
                           Real                  prev_time, 
                           Real                  prev_pres_time,
                           Real                  dt, 
                           int                   NUM_STATE,
                           Real                  be_cn_theta,
                           bool                  modify_reflux_normal_vel,
                           int                   do_mom_diff,
                           const Array<int>&     increment_sync)
{
    if (modify_reflux_normal_vel)
        BoxLib::Abort("modify_reflux_normal_vel is no longer supported");
    //
    // Get parameters.
    //
    const BoxArray& grids               = LevelData[level].boxArray();
    const Geometry& geom                = parent->Geom(level);
    const Real*     dx                  = geom.CellSize();
    const int       numscal             = NUM_STATE - BL_SPACEDIM;
    MultiFab*       mac_sync_phi        = &mac_phi_crse[level];
    NavierStokes&   ns_level            = *(NavierStokes*) &(parent->getLevel(level));
    Godunov*        godunov             = ns_level.godunov;
    bool            use_forces_in_trans = godunov->useForcesInTrans() ? true : false;

    MultiFab vel_visc_terms(grids,BL_SPACEDIM,1);
    MultiFab scal_visc_terms(grids,numscal,1);

    vel_visc_terms.setVal(0,1);  // Initialize to make calls below safe
    scal_visc_terms.setVal(0,1); // Initialize to make calls below safe
    //
    // Get viscous forcing.
    //
    if (be_cn_theta != 1.0) 
    {
        bool do_get_visc_terms = false;

        for (int i=0; i < BL_SPACEDIM; ++i)
            if (!increment_sync.size() || increment_sync[i]==1)
                do_get_visc_terms = true;

        if (do_get_visc_terms || use_forces_in_trans)
            ns_level.getViscTerms(vel_visc_terms,Xvel,BL_SPACEDIM,prev_time);

        do_get_visc_terms = false;
        for (int i=BL_SPACEDIM; i < increment_sync.size(); ++i)
            if (!increment_sync.size() || increment_sync[i]==1)
                do_get_visc_terms = true;

        if (do_get_visc_terms)
            ns_level.getViscTerms(scal_visc_terms,BL_SPACEDIM,numscal,prev_time);
    }

    Array<int> ns_level_bc, bndry[BL_SPACEDIM];

    MultiFab Gp(grids,BL_SPACEDIM,1);

    ns_level.getGradP(Gp, prev_pres_time);

    MultiFab* divu_fp = ns_level.getDivCond(1,prev_time);
    //
    // Compute the mac sync correction.
    //
    FArrayBox xflux, yflux, zflux, tforces, tvelforces, U, area[BL_SPACEDIM], volume;
    FArrayBox grad_phi[BL_SPACEDIM], Rho;

    for (FillPatchIterator S_fpi(ns_level,vel_visc_terms,Godunov::hypgrow(),
                                 prev_time,State_Type,0,NUM_STATE);
         S_fpi.isValid();
         ++S_fpi)
    {
        const int i     = S_fpi.index();
        FArrayBox& S    = S_fpi();
        FArrayBox& divu = (*divu_fp)[i];
        //
        // Step 1: compute ucorr = grad(phi)/rhonph
        //
        // Create storage for corrective velocities.
        //
        Rho.resize(BoxLib::grow(grids[i],1),1);

        D_TERM(grad_phi[0].resize(BoxLib::surroundingNodes(grids[i],0),1);,
               grad_phi[1].resize(BoxLib::surroundingNodes(grids[i],1),1);,
               grad_phi[2].resize(BoxLib::surroundingNodes(grids[i],2),1););

        mac_vel_update(1,D_DECL(grad_phi[0],grad_phi[1],grad_phi[2]),
                       (*mac_sync_phi)[S_fpi], &(*rho_half)[S_fpi],
                       0, grids[i], level, i, dx, dt/2.0);
        //
        // Step 2: compute Mac correction by calling GODUNOV box
        //
        // Get needed data.
        //
        Rho.copy(S,Density,0,1);

#ifdef BOUSSINESQ
        ns_level.getForce(tforces,i,1,0,NUM_STATE,prev_time,S_fpi());
#else
#ifdef GENGETFORCE
        ns_level.getForce(tforces,i,1,0,NUM_STATE,prev_time,Rho);
#elif MOREGENGETFORCE
        ns_level.getForce(tforces,i,1,0,NUM_STATE,prev_time,S_fpi(),S_fpi(),Density);
#else
        ns_level.getForce(tforces,i,1,0,NUM_STATE,Rho);
#endif		 
#endif		 
        //
        // Compute total forcing terms.
        //
        godunov->Sum_tf_gp_visc(tforces, 0, vel_visc_terms[S_fpi], 0, Gp[i], 0, Rho, 0);
        godunov->Sum_tf_divu_visc(S, BL_SPACEDIM, tforces, BL_SPACEDIM, numscal,
                                  scal_visc_terms[S_fpi], 0, divu, 0, Rho, 0, 1);
        if (use_forces_in_trans)
        {
#ifdef BOUSSINESQ
            ns_level.getForce(tvelforces,i,1,Xvel,BL_SPACEDIM,prev_time,S_fpi());
#else
#ifdef GENGETFORCE
            ns_level.getForce(tvelforces,i,1,Xvel,BL_SPACEDIM,prev_time,Rho);
#elif MOREGENGETFORCE
            ns_level.getForce(tvelforces,i,1,Xvel,BL_SPACEDIM,prev_time,S_fpi(),S_fpi(),Density);
#else
            ns_level.getForce(tvelforces,i,1,Xvel,BL_SPACEDIM,Rho);
#endif		 
#endif		 
	    godunov->Sum_tf_gp_visc(tvelforces,0,vel_visc_terms[S_fpi],0,Gp[i],0,Rho,0);
        }
        //
        // Set up the workspace for the godunov Box.
        //
        D_TERM(bndry[0] = ns_level.getBCArray(State_Type,i,0,1);,
               bndry[1] = ns_level.getBCArray(State_Type,i,1,1);,
               bndry[2] = ns_level.getBCArray(State_Type,i,2,1);)

        godunov->Setup(grids[i], dx, dt, 0,
                       D_DECL(xflux,yflux,zflux),
                       D_DECL(bndry[0].dataPtr(),bndry[1].dataPtr(),bndry[2].dataPtr()),
                       D_DECL(S,S,S), D_DECL(0,1,2), tvelforces, 0);

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
            geom.GetFaceArea(area[dir],grids,i,dir,GEOM_GROW);

        geom.GetVolume(volume,grids,i,GEOM_GROW);
        //
        // Get the sync FABS.
        //
        FArrayBox& u_sync = (*Vsync)[S_fpi];
        FArrayBox& s_sync = (*Ssync)[S_fpi];

        U.resize(S.box(),BL_SPACEDIM);
        U.copy(S_fpi(),0,0,BL_SPACEDIM);
        //
        // Loop over state components and compute the sync advective component.
        //
        D_TERM(FArrayBox& u_mac_fab0 = u_mac[0][S_fpi];,
               FArrayBox& u_mac_fab1 = u_mac[1][S_fpi];,
               FArrayBox& u_mac_fab2 = u_mac[2][S_fpi];);

        for (int comp = 0; comp < NUM_STATE; comp++)
        {
            if (!increment_sync.size() || increment_sync[comp]==1)
            {
                const int  sync_ind = comp < BL_SPACEDIM ? comp  : comp-BL_SPACEDIM;
                FArrayBox& temp     = comp < BL_SPACEDIM ? u_sync : s_sync;
                ns_level_bc         = ns_level.getBCArray(State_Type,i,comp,1);

                int use_conserv_diff = (advectionType[comp] == Conservative) ? true : false;

                if (do_mom_diff == 1 && comp < BL_SPACEDIM)
                {
                        S.mult(S,      S.box(),      S.box(),Density,comp,1);
                  tforces.mult(S,tforces.box(),tforces.box(),Density,comp,1);
                }

                godunov->SyncAdvect(grids[i], dx, dt, level, 
                                    area[0], u_mac_fab0, grad_phi[0], xflux, 
                                    area[1], u_mac_fab1, grad_phi[1], yflux,
#if (BL_SPACEDIM == 3)                            
                                    area[2], u_mac_fab2, grad_phi[2], zflux,
#endif
                                    U, S, tforces, divu, comp, temp, sync_ind,
                                    use_conserv_diff, comp,
                                    ns_level_bc.dataPtr(), PRE_MAC, volume);
                //
                // NOTE: the signs here are opposite from VELGOD.
                // NOTE: fluxes expected to be in extensive form.
                //
                if (level > 0)
                {
                    D_TERM(adv_flux_reg->FineAdd(xflux,0,i,0,comp,1,-dt);,
                           adv_flux_reg->FineAdd(yflux,1,i,0,comp,1,-dt);,
                           adv_flux_reg->FineAdd(zflux,2,i,0,comp,1,-dt););
                }
            }
        }
        //
        // Include grad_phi in the mac registers corresponding
        // to the next coarsest interface.
        //
        if (level > 0)
        {
            const Real mlt =  -1.0/( (double) parent->nCycle(level));
            D_TERM(mac_reg[level].FineAdd(grad_phi[0],area[0],0,i,0,0,1,mlt);,
                   mac_reg[level].FineAdd(grad_phi[1],area[1],1,i,0,0,1,mlt);,
                   mac_reg[level].FineAdd(grad_phi[2],area[2],2,i,0,0,1,mlt););
        }
        //
        // Multiply the sync term by dt -- now done in the calling routine.
        //
    }

    delete divu_fp;
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
                           MultiFab*              Sync,
                           int                    comp,
                           int                    s_ind,
                           const MultiFab* const* sync_edges,
			   int                    eComp,
                           MultiFab*              rho_half,
                           FluxRegister*          adv_flux_reg,
                           Array<AdvectionForm>&  advectionType, 
			   bool                   modify_reflux_normal_vel,
                           Real                   dt)
{
    if (modify_reflux_normal_vel)
        BoxLib::Abort("modify_reflux_normal_vel is no longer supported");

    FArrayBox xflux, yflux, zflux, grad_phi[BL_SPACEDIM], area[BL_SPACEDIM], volume;

    const BoxArray& grids        = LevelData[level].boxArray();
    const Geometry& geom         = parent->Geom(level);
    MultiFab*       mac_sync_phi = &mac_phi_crse[level];

    Godunov godunov(512);
    //
    // Compute the mac sync correction.
    //
    for (MFIter Syncmfi(*Sync); Syncmfi.isValid(); ++Syncmfi)
    {
        const int  i   = Syncmfi.index();
        const Box& grd = grids[i];
        //
        // Step 1: compute ucorr = grad(phi)/rhonph
        //
        D_TERM(grad_phi[0].resize(BoxLib::surroundingNodes(grd,0),1);,
               grad_phi[1].resize(BoxLib::surroundingNodes(grd,1),1);,
               grad_phi[2].resize(BoxLib::surroundingNodes(grd,2),1););

        mac_vel_update(1,
                       D_DECL(grad_phi[0],grad_phi[1],grad_phi[2]),
                       (*mac_sync_phi)[Syncmfi],
                       &(*rho_half)[Syncmfi], 0,
                       grd, level, i,
                       geom.CellSize(), dt/2.0);
        //
        // Step 2: compute Mac correction by advecting the edge states.
        //
        D_TERM(xflux.resize(BoxLib::surroundingNodes(grd,0),1);,
               yflux.resize(BoxLib::surroundingNodes(grd,1),1);,
               zflux.resize(BoxLib::surroundingNodes(grd,2),1););

        D_TERM(xflux.copy((*sync_edges[0])[Syncmfi],eComp,0,1);,
               yflux.copy((*sync_edges[1])[Syncmfi],eComp,0,1);,
               zflux.copy((*sync_edges[2])[Syncmfi],eComp,0,1););

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            geom.GetFaceArea(area[dir],grids,i,dir,GEOM_GROW);
        }
        geom.GetVolume(volume,grids,i,GEOM_GROW);

        int use_conserv_diff = (advectionType[comp] == Conservative) ? true : false;

        godunov.ComputeSyncAofs(grd,
                                area[0],
                                grad_phi[0],       xflux,
                                
                                area[1],
                                grad_phi[1],       yflux,
#if (BL_SPACEDIM == 3)                            
                                area[2],
                                grad_phi[2],       zflux,
#endif
                                volume, (*Sync)[Syncmfi],
                                s_ind, use_conserv_diff);

        D_TERM(area[0].clear();, area[1].clear();, area[2].clear(););
        D_TERM(grad_phi[0].clear();, grad_phi[1].clear();, grad_phi[2].clear(););
        //
        // NOTE: the signs here are opposite from VELGOD.
        // NOTE: fluxes expected to be in extensive form.
        //
        if (level > 0)
        {
            D_TERM(adv_flux_reg->FineAdd(xflux,0,i,0,comp,1,-dt);,
                   adv_flux_reg->FineAdd(yflux,1,i,0,comp,1,-dt);,
                   adv_flux_reg->FineAdd(zflux,2,i,0,comp,1,-dt););
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
    const BoxArray& grids = LevelData[level].boxArray();

    Real sum = 0.0;

    FArrayBox dmac, area[BL_SPACEDIM], volume;

    for (MFIter U_edge0mfi(U_edge[0]); U_edge0mfi.isValid(); ++U_edge0mfi)
    {
        dmac.resize(grids[U_edge0mfi.index()],1);

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            parent->Geom(level).GetFaceArea(area[dir],grids,U_edge0mfi.index(),dir,GEOM_GROW);
        }
        parent->Geom(level).GetVolume(volume,grids,U_edge0mfi.index(),GEOM_GROW);

        const FArrayBox& uxedge = U_edge[0][U_edge0mfi];
        const FArrayBox& uyedge = U_edge[1][U_edge0mfi];
        const FArrayBox& xarea  = area[0];
        const FArrayBox& yarea  = area[1];
        const FArrayBox& vol    = volume;

        DEF_LIMITS(dmac,dmac_dat,dlo,dhi);
        DEF_CLIMITS(uxedge,ux_dat,uxlo,uxhi);
        DEF_CLIMITS(uyedge,uy_dat,uylo,uyhi);
        DEF_CLIMITS(xarea,ax_dat,axlo,axhi);
        DEF_CLIMITS(yarea,ay_dat,aylo,ayhi);
        DEF_CLIMITS(vol,vol_dat,vlo,vhi);

#if (BL_SPACEDIM == 2)
        FORT_MACDIV(dmac_dat,ARLIM(dlo),ARLIM(dhi),dlo,dhi,
                    ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                    uy_dat,ARLIM(uylo),ARLIM(uyhi),
                    ax_dat,ARLIM(axlo),ARLIM(axhi), 
                    ay_dat,ARLIM(aylo),ARLIM(ayhi), 
                    vol_dat,ARLIM(vlo),ARLIM(vhi));
#endif

#if (BL_SPACEDIM == 3)
        const FArrayBox& uzedge = U_edge[2][U_edge0mfi];
        DEF_CLIMITS(uzedge,uz_dat,uzlo,uzhi);
        const FArrayBox& zarea = area[2];
        DEF_CLIMITS(zarea,az_dat,azlo,azhi);

        FORT_MACDIV(dmac_dat,ARLIM(dlo),ARLIM(dhi),dlo,dhi,
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

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceRealSum(sum,IOProc);
        
        if (ParallelDescriptor::IOProcessor())
            std::cout << "SUM of DIV(U_edge) = " << sum << '\n';
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

    const BoxArray&   grids  = LevelData[level].boxArray();
    const Geometry&   geom   = parent->Geom(level);
    const Box&        domain = parent->Geom(level).Domain();
    //
    // Create 1-wide cc box just outside boundary to hold phi.
    // Create 1-wide cc box just inside  boundary to hold rho,u,divu.
    //
    BoxList ccBoxList, phiBoxList;

    for (int iface = 0; iface < numOutFlowFaces; iface++)
    {
        if (grids_on_side_of_domain(grids,geom.Domain(),outFaces[iface])) 
	{
            const int outDir    = outFaces[iface].coordDir();

            Box ccBndBox;
            if (outFaces[iface].faceDir() == Orientation::high)
	    {
                ccBndBox = BoxLib::adjCellHi(domain,outDir,2);
                ccBndBox.shift(outDir,-2);
	    } 
            else 
	    {
                ccBndBox = BoxLib::adjCellLo(domain,outDir,2);
                ccBndBox.shift(outDir,2);
	    }
            ccBoxList.push_back(ccBndBox);

            Box phiBox  = BoxLib::adjCell(domain,outFaces[iface],1);
            phiBoxList.push_back(phiBox);

            const Box      valid_ccBndBox       = ccBndBox & domain;
            const BoxArray uncovered_outflow_ba = BoxLib::complementIn(valid_ccBndBox,grids);

            if (uncovered_outflow_ba.size() && 
                BoxLib::intersect(grids,valid_ccBndBox).size())
                BoxLib::Error("MacProj: Cannot yet handle partially refined outflow");
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
      
        for ( int iface = 0; iface < numOutFlowFaces; ++iface) 
	{
            rhodat[iface].resize(ccBoxArray[iface], 1);
            divudat[iface].resize(ccBoxArray[iface], 1);
            phidat[iface].resize(phiBoxArray[iface], 1);

            phidat[iface].setVal(0.0);
            divu.copy(divudat[iface]);
            S.copy(rhodat[iface], Density, 0, 1);
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
            for ( int iface = 0; iface < numOutFlowFaces; ++iface) 
	    {
                uedat[i][iface].resize(edgeArray[iface], 1);
                u_mac[i].copy(uedat[i][iface], 0, 0, 1);
	    }
	}
    
        MacOutFlowBC macBC;

        NavierStokes* ns_level = dynamic_cast<NavierStokes*>(&parent->getLevel(level));
        Real gravity = ns_level->getGravity();
        const int* lo_bc = phys_bc->lo();
        const int* hi_bc = phys_bc->hi();
        macBC.computeBC(uedat, divudat, rhodat, phidat,
                        geom, outFaces, numOutFlowFaces, lo_bc, hi_bc, umac_periodic_test_Tol, gravity);
        //
        // Must do this kind of copy instead of mac_phi->copy(phidat);
        // because we're copying onto the ghost cells of the FABs,
        // not the valid regions.
        //
        for ( int iface = 0; iface < numOutFlowFaces; ++iface )
	{
            for (MFIter mfi(*mac_phi); mfi.isValid(); ++mfi)
	    {
                Box ovlp = (*mac_phi)[mfi].box() & phidat[iface].box();
                if (ovlp.ok())
                    (*mac_phi)[mfi].copy(phidat[iface],ovlp);
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
    Array<IntVect>         pshifts(27);
    MultiFabCopyDescriptor mfcd;
    std::vector<TURec>     pirm;
    MultiFabId             mfid[BL_SPACEDIM];

    std::vector< std::pair<int,Box> > isects;

    for (int dim = 0; dim < BL_SPACEDIM; dim++)
    {
        if (geom.isPeriodic(dim))
        {
            Box eDomain = BoxLib::surroundingNodes(geom.Domain(),dim);

            mfid[dim] = mfcd.RegisterMultiFab(&u_mac[dim]);

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
                        const Box  dstBox = srcBox - pshifts[iiv];

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

            u_mac[dim].boxArray().clear_hash_bin();
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
            std::cout << "dir = "         << dim
                      << ", diff norm = " << max_norm
                      << " for region: "  << pirm[i].m_dstBox << std::endl;
            BoxLib::Error("Periodic bust in u_mac");
        }
    }
}

void
MacProj::scaleArea (int level, MultiFab* area, Real** anel_coeff)
{
    const BoxArray& grids = LevelData[level].boxArray();

    int mult = 1;

    for (MFIter mfi(*area); mfi.isValid(); ++mfi)
    {
        const int        i      = mfi.index();
        const FArrayBox& xarea  = area[0][mfi];
        const FArrayBox& yarea  = area[1][mfi];

        DEF_CLIMITS(xarea,ax_dat,axlo,axhi);
        DEF_CLIMITS(yarea,ay_dat,aylo,ayhi);

        const int* lo = grids[i].loVect();
        const int* hi = grids[i].hiVect();

        int anel_coeff_lo = lo[BL_SPACEDIM-1]-1;
        int anel_coeff_hi = hi[BL_SPACEDIM-1]+1;

#if (BL_SPACEDIM == 2)
        FORT_SCALEAREA(ax_dat,ARLIM(axlo),ARLIM(axhi), 
                       ay_dat,ARLIM(aylo),ARLIM(ayhi), 
                       anel_coeff[i],&anel_coeff_lo,&anel_coeff_hi,lo,hi,&mult);

#elif (BL_SPACEDIM == 3)
        const FArrayBox& zarea = area[2][mfi];
        DEF_CLIMITS(zarea,az_dat,azlo,azhi);
        FORT_SCALEAREA(ax_dat,ARLIM(axlo),ARLIM(axhi), 
                       ay_dat,ARLIM(aylo),ARLIM(ayhi), 
                       az_dat,ARLIM(azlo),ARLIM(azhi), 
                       anel_coeff[i],&anel_coeff_lo,&anel_coeff_hi,lo,hi,&mult);

#endif
    }
}

void
MacProj::rescaleArea (int level, MultiFab* area, Real** anel_coeff)
{
    const BoxArray& grids = LevelData[level].boxArray();

    int mult = -1;

    for (MFIter mfi(*area); mfi.isValid(); ++mfi)
    {
        const int        i      = mfi.index();
        const FArrayBox& xarea  = area[0][mfi];
        const FArrayBox& yarea  = area[1][mfi];

        DEF_CLIMITS(xarea,ax_dat,axlo,axhi);
        DEF_CLIMITS(yarea,ay_dat,aylo,ayhi);

        const int* lo        = grids[i].loVect();
        const int* hi        = grids[i].hiVect();

        int anel_coeff_lo = lo[BL_SPACEDIM-1]-1;
        int anel_coeff_hi = hi[BL_SPACEDIM-1]+1;

#if (BL_SPACEDIM == 2)
        FORT_SCALEAREA(ax_dat,ARLIM(axlo),ARLIM(axhi), 
                       ay_dat,ARLIM(aylo),ARLIM(ayhi), 
                       anel_coeff[i],&anel_coeff_lo,&anel_coeff_hi,lo,hi,&mult);

#elif (BL_SPACEDIM == 3)
        const FArrayBox& zarea = area[2][mfi];
        DEF_CLIMITS(zarea,az_dat,azlo,azhi);
        FORT_SCALEAREA(ax_dat,ARLIM(axlo),ARLIM(axhi), 
                       ay_dat,ARLIM(aylo),ARLIM(ayhi), 
                       az_dat,ARLIM(azlo),ARLIM(azhi), 
                       anel_coeff[i],&anel_coeff_lo,&anel_coeff_hi,lo,hi,&mult);

#endif
    }
}
