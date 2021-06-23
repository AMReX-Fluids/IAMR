

#include <AMReX_LO_BCTYPES.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <MacProj.H>
#include <AMReX_MacBndry.H>
#include <MacOpMacDrivers.H>
#include <NavierStokesBase.H>
#include <OutFlowBC.H>

#ifdef AMREX_USE_EB
#include <hydro_ebgodunov.H>
#include <hydro_ebmol.H>
#else
#include <hydro_godunov.H>
#include <hydro_mol.H>
#endif



//fixme, for writesingle level plotfile
#include<AMReX_PlotFileUtil.H>

using namespace amrex;

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)      \
    const int* fablo = (fab).loVect();          \
    const int* fabhi = (fab).hiVect();          \
    Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)     \
    const int* fablo = (fab).loVect();          \
    const int* fabhi = (fab).hiVect();          \
    const Real* fabdat = (fab).dataPtr();

#define DEF_BOX_LIMITS(box,boxlo,boxhi)         \
    const int* boxlo = (box).loVect();          \
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
int  MacProj::do_outflow_bcs;
int  MacProj::fix_mac_sync_rhs;
int  MacProj::check_umac_periodicity;

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
    finest_level(_finest_level)
{
    Initialize();

    if (verbose) amrex::Print() << "Creating mac_projector\n";

    finest_level_allocated = finest_level;
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
            mac_phi_crse[level].reset(new MultiFab(grids,dmap,1,1, MFInfo(), LevelData[level]->Factory()));
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
                      const BCRec&    density_math_bc,
                      bool            increment_vel_register )
{
    BL_PROFILE("MacProj::mac_project()");
    if (verbose) amrex::Print() << "... mac_project at level " << level << '\n';

    const BoxArray& grids      = LevelData[level]->boxArray();
    const DistributionMapping& dmap = LevelData[level]->DistributionMap();
    const int       max_level  = parent->maxLevel();
    MultiFab*       mac_phi    = 0;
    NavierStokesBase&   ns     = *(NavierStokesBase*) &(parent->getLevel(level));
    const MultiFab*     area   = ns.Area();
    //
    // If finest level possible no need to make permanent mac_phi for bcs.
    //
    std::unique_ptr<MultiFab> raii;
    if (level == max_level) {
        raii.reset(new MultiFab(grids,dmap,1,1, MFInfo(), LevelData[level]->Factory()));
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
    MultiFab Rhs(grids,dmap,1,0, MFInfo(), LevelData[level]->Factory());

    MultiFab::Copy(Rhs,divu,0,0,1,0); 

    MultiFab* cphi = (level == 0) ? nullptr : mac_phi_crse[level-1].get();
    mlmg_mac_level_solve(parent, cphi, *phys_bc, density_math_bc, level, Density, mac_tol, mac_abs_tol,
                         rhs_scale, S, Rhs, u_mac, mac_phi, verbose);


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

            for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
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
                         const BCRec& rho_math_bc,
                         IntVect&  fine_ratio,
                         Array<MultiFab*,AMREX_SPACEDIM>& Ucorr,
                         MultiFab* Rhs_increment )
{
    BL_ASSERT(level < finest_level);

    if (verbose) amrex::Print() << "... mac_sync_solve at level " << level << '\n';

    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real      strt_time  = ParallelDescriptor::second();
    const BoxArray& grids      = LevelData[level]->boxArray();
    const DistributionMapping& dmap = LevelData[level]->DistributionMap();
    const Geometry& geom       = parent->Geom(level);
    const BoxArray& fine_boxes = LevelData[level+1]->boxArray();
    const NavierStokesBase& ns_level   = *(NavierStokesBase*) &(parent->getLevel(level));
    const MultiFab&     volume     = ns_level.Volume();
    const MultiFab*     area       = ns_level.Area();
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
    MultiFab Rhs(grids,dmap,1,0, MFInfo(), LevelData[level]->Factory());
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
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    // fixme? Should do some real tests to see if tiling here is a win or not
    for (MFIter mfi(Rhs,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        BL_ASSERT(grids[mfi.index()].contains(mfi.tilebox()) );

        const std::vector< std::pair<int,Box> >& isects = baf.intersections(mfi.tilebox());

        auto const& rhs = Rhs.array(mfi);
	for (const auto& is : isects)
	{
	    amrex::ParallelFor(is.second, [rhs]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               rhs(i,j,k) = 0.0;
            });
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

    MultiFab area_tmp[AMREX_SPACEDIM];

    //
    // Compute Ucorr, including filling ghost cells for EB
    //
    mlmg_mac_sync_solve(parent,*phys_bc, rho_math_bc, level, mac_sync_tol, mac_abs_tol,
                        rhs_scale, area, volume, rho_half, Rhs, mac_sync_phi,
                        Ucorr, verbose);

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
// 1. compute u_corr as the gradient of mac_sync_phi --- Now part of mac_sync_solve
// 2. compute advective tendency of u_corr and
//    add into Vsync or Ssync
//
// If increment_sync is non-null, the (i-BL_SPACEDIM)-th component
// of (*Ssync) is incremented only when increment[i]==1
// This is useful if that component gets incrmnted in a non-standard way.
//
void
MacProj::mac_sync_compute (int                   level,
                           Array<MultiFab*,AMREX_SPACEDIM>& Ucorr,
                           MultiFab*             u_mac,
                           MultiFab&             Vsync,
                           MultiFab&             Ssync,
                           MultiFab&             /*rho_half*/,
                           FluxRegister*         adv_flux_reg,
                           Vector<AdvectionForm>& advectionType,
                           Real                  prev_time,
                           Real                  /*prev_pres_time*/,
                           Real                  dt,
                           int                   NUM_STATE,
                           Real                  be_cn_theta,
                           bool                  modify_reflux_normal_vel,
                           int                   do_mom_diff,
                           const Vector<int>&    increment_sync,
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
    const int       numscal             = NUM_STATE - AMREX_SPACEDIM;
    NavierStokesBase&   ns_level        = *(NavierStokesBase*) &(parent->getLevel(level));
    const MultiFab* area                = ns_level.Area();

    const int  ncomp = 1;         // Number of components to process at once

    const int  nghost  = 0;

    //
    // Prep MFs to store fluxes and edge states
    //
    MultiFab fluxes[AMREX_SPACEDIM];
    MultiFab edgestate[AMREX_SPACEDIM];

    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        const BoxArray& ba = LevelData[level]->getEdgeBoxArray(i);
        fluxes[i].define(ba, dmap, NUM_STATE, nghost, MFInfo(),ns_level.Factory());
        edgestate[i].define(ba, dmap, ncomp, nghost, MFInfo(), ns_level.Factory());
    }


    // Visc terms, is not used for MOL but we define it here anyways as base for the following
    // FillPatch operator
    MultiFab vel_visc_terms(grids,dmap,AMREX_SPACEDIM,
                            ns_level.nghost_force(),MFInfo(),ns_level.Factory());
    FillPatchIterator S_fpi(ns_level,vel_visc_terms,ns_level.nghost_state(),
                            prev_time,State_Type,0,NUM_STATE);
    MultiFab& Smf = S_fpi.get_mf();

    //
    // Compute the mac sync correction.
    //
    if (!ns_level.use_godunov)   // MOL ====================================================================
    {
        Vector<BCRec>  math_bcs(ncomp);

        for (int comp = 0; comp < NUM_STATE; ++comp)
        {
            if (increment_sync.empty() || increment_sync[comp]==1)
            {
                // Get BCs for this component
                math_bcs = ns_level.fetchBCArray(State_Type, comp, ncomp);

                // Select sync MF and its component for processing
                const int  sync_comp = comp < AMREX_SPACEDIM ? comp   : comp-AMREX_SPACEDIM;
                MultiFab*  sync_ptr  = comp < AMREX_SPACEDIM ? &Vsync : &Ssync;
                bool    is_velocity  = comp < AMREX_SPACEDIM ? true   : false;
		// Using fetchBCArray seems a better option than this.
		// Vector<BCRec> const& bcs  = comp < AMREX_SPACEDIM
		// 				   ? ns_level.get_bcrec_velocity()
		// 				   : ns_level.get_bcrec_scalars();
		// // create a single element amrex::Vector
		// Vector<BCRec> const& math_bcs= {bcs[sync_comp]};

		BCRec  const* d_bcrec_ptr = comp < AMREX_SPACEDIM
						   ? ns_level.get_bcrec_velocity_d_ptr()
						   : ns_level.get_bcrec_scalars_d_ptr();

#ifdef AMREX_USE_EB
                EBMOL::ComputeSyncAofs(*sync_ptr, sync_comp, ncomp, Smf, comp,
                                       D_DECL(u_mac[0],u_mac[1],u_mac[2]),
                                       D_DECL(*Ucorr[0],*Ucorr[1],*Ucorr[2]),
                                       D_DECL(edgestate[0],edgestate[1],edgestate[2]), 0, false,
                                       D_DECL(fluxes[0],fluxes[1],fluxes[2]), comp,
                                       math_bcs, &d_bcrec_ptr[sync_comp], geom, dt,
                                       is_velocity, ns_level.redistribution_type );
#else
                MOL::ComputeSyncAofs(*sync_ptr, sync_comp, ncomp, Smf, comp,
                                     D_DECL(u_mac[0],u_mac[1],u_mac[2]),
                                     D_DECL(*Ucorr[0],*Ucorr[1],*Ucorr[2]),
                                     D_DECL(edgestate[0],edgestate[1],edgestate[2]), 0, false,
                                     D_DECL(fluxes[0],fluxes[1],fluxes[2]), comp,
                                     math_bcs, &d_bcrec_ptr[sync_comp], geom, is_velocity );
#endif
            }
        }
    }
    else // GODUNOV ========================================================================================
    {
        std::unique_ptr<MultiFab> divu_fp (ns_level.getDivCond(ns_level.nghost_force(),prev_time));

        MultiFab& Gp = ns_level.get_old_data(Gradp_Type);

        MultiFab scal_visc_terms(grids,dmap,numscal,ns_level.nghost_force(),
                                 MFInfo(),ns_level.Factory());
        vel_visc_terms.setVal(0.0);  // Initialize to make calls below safe
        scal_visc_terms.setVal(0.0); // Initialize to make calls below safe

        bool use_forces_in_trans = ns_level.GodunovUseForcesInTrans(); // This should always return False for EB Godunov

        // Get viscous forcing.
        if (be_cn_theta != 1.0)
        {
            bool do_get_visc_terms = false;

            for (int i=0; i < AMREX_SPACEDIM; ++i)
                if (increment_sync.empty() || increment_sync[i]==1)
                    do_get_visc_terms = true;

            if (do_get_visc_terms || use_forces_in_trans)
                ns_level.getViscTerms(vel_visc_terms,Xvel,AMREX_SPACEDIM,prev_time);

            do_get_visc_terms = false;
            for (int i=AMREX_SPACEDIM; i < increment_sync.size(); ++i)
                if (increment_sync.empty() || increment_sync[i]==1)
                    do_get_visc_terms = true;

            if (do_get_visc_terms)
                ns_level.getViscTerms(scal_visc_terms,AMREX_SPACEDIM,numscal,prev_time);
        }


        MultiFab forcing_term(grids, dmap, NUM_STATE, ns_level.nghost_force());

        // Store momenta multifab if conservative approach is used,
        // i.e. rho* u.
        // We make it with AMREX_SPACEDIM components instead of only one
        // (loop below is done component by component) because ComputeSyncAofs will
        // need to know which component of velocity is being processed.
        MultiFab momenta;
        if  (do_mom_diff == 1)
        {
            momenta.define(grids,dmap, AMREX_SPACEDIM, Smf.nGrow(), MFInfo(), Smf.Factory());
            MultiFab::Copy(momenta,Smf,0,0,AMREX_SPACEDIM, Smf.nGrow());
            for (int d=0; d < AMREX_SPACEDIM; ++d )
                MultiFab::Multiply( momenta, Smf, Density, d, 1, Smf.nGrow());
        }


        //
        // Compute forcing terms for all component
        //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            FArrayBox tforces;
            for (MFIter Smfi(Smf,TilingIfNotGPU()); Smfi.isValid(); ++Smfi)
            {
                const auto gbx = Smfi.growntilebox(ns_level.nghost_force());

                ns_level.getForce(forcing_term[Smfi],gbx,0,NUM_STATE,
                                  prev_time,Smf[Smfi],Smf[Smfi],Density,Smfi);
            }
        }

        for (int comp = 0; comp < NUM_STATE; ++comp)
        {
            if (increment_sync.empty() || increment_sync[comp]==1)
            {

                //
                // Compute total forcing term
                //
                for (MFIter Smfi(Smf,TilingIfNotGPU()); Smfi.isValid(); ++Smfi)
                {
                    auto const gbx = Smfi.growntilebox(ns_level.nghost_force());

                    //
                    // Compute total forcing terms.
                    //
                    auto const& tf    = forcing_term.array(Smfi,comp);

                    if (comp < AMREX_SPACEDIM)  // Velocity/Momenta
                    {
                        auto const& visc = vel_visc_terms[Smfi].const_array(comp);
                        auto const& gp   = Gp[Smfi].const_array(comp);

                        if ( do_mom_diff == 0 )
                        {
                            auto const& rho   = Smf[Smfi].const_array(Density);

                            amrex::ParallelFor(gbx, [tf, visc, gp, rho]
                            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                            {
                                tf(i,j,k)  += visc(i,j,k) - gp(i,j,k);
                                tf(i,j,k)  /= rho(i,j,k);
                            });
                        }
                        else
                        {
                            amrex::ParallelFor(gbx, [tf, visc, gp]
                            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                            {
                                tf(i,j,k)  += visc(i,j,k) - gp(i,j,k);
                            });
                        }

                    }
                    else  // Scalars
                    {
                        auto const& visc = scal_visc_terms[Smfi].const_array(comp-AMREX_SPACEDIM);
                        auto const& S    = Smf.const_array(Smfi,comp);
                        auto const& divu = divu_fp -> const_array(Smfi);
                        amrex::ParallelFor(gbx, [tf, visc, S, divu]
                        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            tf(i,j,k) += visc(i,j,k) - S(i,j,k) * divu(i,j,k);
                        });
                    }
                }


                //
                // Perform sync
                //

                // Select sync MF and its component for processing
                const int  sync_comp   = comp < AMREX_SPACEDIM ? comp   : comp-AMREX_SPACEDIM;
                MultiFab*  sync_ptr    = comp < AMREX_SPACEDIM ? &Vsync : &Ssync;
                const bool is_velocity = comp < AMREX_SPACEDIM ? true   : false;
                BCRec  const* d_bcrec_ptr = comp < AMREX_SPACEDIM
                                                   ? &(ns_level.get_bcrec_velocity_d_ptr())[sync_comp]
                                                   : &(ns_level.get_bcrec_scalars_d_ptr())[sync_comp];

                const auto& Q = (do_mom_diff == 1 and comp < AMREX_SPACEDIM) ? momenta : Smf;

                amrex::Gpu::DeviceVector<int> iconserv;
                iconserv.resize(1, 0);
                iconserv[0] = (advectionType[comp] == Conservative) ? 1 : 0;

#ifdef AMREX_USE_EB
                Vector<BCRec> bcrec_ptr = comp < AMREX_SPACEDIM
                                                 ? ns_level.m_bcrec_velocity
                                                 : ns_level.m_bcrec_scalars;
                EBGodunov::ComputeSyncAofs(*sync_ptr, sync_comp, ncomp,
                                           Q, comp,
                                           AMREX_D_DECL(u_mac[0],u_mac[1],u_mac[2]),
                                           AMREX_D_DECL(*Ucorr[0],*Ucorr[1],*Ucorr[2]),
                                           AMREX_D_DECL(edgestate[0],edgestate[1],edgestate[2]), 0, false,
                                           AMREX_D_DECL(fluxes[0],fluxes[1],fluxes[2]), comp,
                                           forcing_term, comp, *divu_fp,
                                           bcrec_ptr, d_bcrec_ptr,
                                           geom, iconserv, dt, is_velocity,
                                           ns_level.redistribution_type);
#else
                Godunov::ComputeSyncAofs(*sync_ptr, sync_comp, ncomp,
                                         Q, comp,
                                         AMREX_D_DECL(u_mac[0],u_mac[1],u_mac[2]),
                                         AMREX_D_DECL(*Ucorr[0],*Ucorr[1],*Ucorr[2]),
                                         AMREX_D_DECL(edgestate[0],edgestate[1],edgestate[2]), 0, false,
                                         AMREX_D_DECL(fluxes[0],fluxes[1],fluxes[2]), comp,
                                         forcing_term, comp, *divu_fp,
                                         d_bcrec_ptr, geom, iconserv, dt,
                                         ns_level.GodunovUsePPM(), ns_level.GodunovUseForcesInTrans(),
                                         is_velocity );
#endif

            }
        }
    }


    if (level > 0 && update_fluxreg)
    {
        const Real mlt =  -1.0/( (double) parent->nCycle(level));
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            for (int comp = 0; comp < NUM_STATE; ++comp)
            {
                if (increment_sync.empty() || increment_sync[comp]==1)
                {
                    adv_flux_reg->FineAdd(fluxes[d],d,comp,comp,1,-dt);
                }
            }
            //
            // Include grad_phi(aka Ucorr) in the mac registers corresponding
            // to the next coarsest interface.
            //
            mac_reg[level]->FineAdd(*Ucorr[d],area[d],d,0,0,1,mlt);
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
                           Array<MultiFab*,AMREX_SPACEDIM>& Ucorr,
                           MultiFab&              Sync,
                           int                    comp,
                           int                    s_ind,
                           MultiFab* const*       sync_edges,
                           int                    eComp,
                           MultiFab&              /*rho_half*/,
                           FluxRegister*          adv_flux_reg,
                           Vector<AdvectionForm>& /*advectionType*/,
                           bool                   modify_reflux_normal_vel,
                           Real                   dt,
                           bool                   update_fluxreg)
{
    if (modify_reflux_normal_vel)
        amrex::Abort("modify_reflux_normal_vel is no longer supported");

    const DistributionMapping& dmap     = LevelData[level]->DistributionMap();
    const Geometry& geom         = parent->Geom(level);
    NavierStokesBase& ns_level   = *(NavierStokesBase*) &(parent->getLevel(level));

    const int  nghost  = 0;

    const int  ncomp   = 1;         // Number of components to process at once

    MultiFab fluxes[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        const BoxArray& ba = LevelData[level]->getEdgeBoxArray(i);
        fluxes[i].define(ba, dmap, 1, sync_edges[0]->nGrow(), MFInfo(),ns_level.Factory());
    }

    //
    // Compute the mac sync correction.
    //
    if (!ns_level.use_godunov)
    {
        //
        // MOL algorithm
        //

	// Bogus arguments -- they will not be used since we don't need to recompute the edge states
	Vector<BCRec>  bcs;
	BCRec  const* d_bcrec_ptr = NULL;

#ifdef AMREX_USE_EB
        EBMOL::ComputeSyncAofs(Sync, s_ind, ncomp,
                               Sync, s_ind, // this is not used when we pass edge states
                               D_DECL(*Ucorr[0],*Ucorr[1],*Ucorr[2]),  // this is not used when we pass edge states
                               D_DECL(*Ucorr[0],*Ucorr[1],*Ucorr[2]),
                               D_DECL(*sync_edges[0],*sync_edges[1],*sync_edges[2]), eComp, true,
                               D_DECL(fluxes[0],fluxes[1],fluxes[2]), 0,
                               bcs, d_bcrec_ptr, geom, dt,
                               false,  // not used when we pass edge states
                               ns_level.redistribution_type);
#else
        MOL::ComputeSyncAofs(Sync, s_ind, ncomp,
			     Sync, s_ind, // this is not used when we pass edge states
			     D_DECL(*Ucorr[0],*Ucorr[1],*Ucorr[2]),  // this is not used when we pass edge states
			     D_DECL(*Ucorr[0],*Ucorr[1],*Ucorr[2]),
			     D_DECL(*sync_edges[0],*sync_edges[1],*sync_edges[2]), eComp, true,
			     D_DECL(fluxes[0],fluxes[1],fluxes[2]), 0,
			     bcs, d_bcrec_ptr, geom,
                             false ); // not used when we pass edge states
#endif

    }
    else
    {
        //
        // Godunov algorithm
        //

        // Bogus arguments -- they will not be used since we don't need to recompute the edge states
        BCRec  const* d_bcrec_ptr = NULL;
        Gpu::DeviceVector<int> iconserv;

#ifdef AMREX_USE_EB
        EBGodunov::ComputeSyncAofs(Sync, s_ind, ncomp,
                                   MultiFab(), s_ind,                      // this is not used when known_edgestate = true
                                   AMREX_D_DECL(*Ucorr[0],*Ucorr[1],*Ucorr[2]),  // this is not used when we pass edge states
                                   AMREX_D_DECL(*Ucorr[0],*Ucorr[1],*Ucorr[2]),
                                   AMREX_D_DECL(*sync_edges[0],*sync_edges[1],*sync_edges[2]), eComp, true,
                                   AMREX_D_DECL(fluxes[0],fluxes[1],fluxes[2]), 0,
                                   MultiFab(), 0, MultiFab(),                        // this is not used when known_edgestate = true
                                   {}, d_bcrec_ptr,
                                   geom, iconserv, dt, true,
                                   ns_level.redistribution_type);
#else
        Godunov::ComputeSyncAofs(Sync, s_ind, ncomp,
                                 MultiFab(), s_ind,                      // this is not used when known_edgestate = true
                                 AMREX_D_DECL(*Ucorr[0],*Ucorr[1],*Ucorr[2]),  // this is not used when we pass edge states
                                 AMREX_D_DECL(*Ucorr[0],*Ucorr[1],*Ucorr[2]),
                                 AMREX_D_DECL(*sync_edges[0],*sync_edges[1],*sync_edges[2]), eComp, true,
                                 AMREX_D_DECL(fluxes[0],fluxes[1],fluxes[2]), 0,
                                 MultiFab(), 0, MultiFab(),                        // this is not used when known_edgestate = true
                                 d_bcrec_ptr, geom, iconserv, 0.0, false, false, false  ); // this is not used when known_edgestate = true
#endif

    }

    if (level > 0 && update_fluxreg)
    {
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
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

    MultiFab dmac(volume.boxArray(),volume.DistributionMap(),1,0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dmac,TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& cc_divu   = dmac.array(mfi);
        D_TERM(auto const& ux_e = U_edge[0].array(mfi);,
               auto const& uy_e = U_edge[1].array(mfi);,
               auto const& uz_e = U_edge[2].array(mfi););
        D_TERM(auto const& xarea  = area[0].array(mfi);,
               auto const& yarea  = area[1].array(mfi);,
               auto const& zarea  = area[2].array(mfi););
        auto const& vol       = volume.array(mfi);

        amrex::ParallelFor(bx, [cc_divu,D_DECL(ux_e,uy_e,uz_e),
                                        D_DECL(xarea,yarea,zarea), vol]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cc_divu(i,j,k) = D_TERM(  xarea(i+1,j,k)*ux_e(i+1,j,k) - xarea(i,j,k)*ux_e(i,j,k),
                                    + yarea(i,j+1,k)*uy_e(i,j+1,k) - yarea(i,j,k)*uy_e(i,j,k),
                                    + zarea(i,j,k+1)*uz_e(i,j,k+1) - zarea(i,j,k)*uz_e(i,j,k));
            cc_divu(i,j,k) /= vol(i,j,k);
        });
    }

    Real sm = amrex::ReduceSum(dmac, 0, []
    AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& dmac_arr) -> Real
    {
        Real tmp = 0.0;
        AMREX_LOOP_3D(bx, i, j, k,
        {
            tmp += dmac_arr(i,j,k);
        });
        return tmp;
    });

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceRealSum(sm,IOProc);

        amrex::Print().SetPrecision(15) << "SUM of DIV(U_edge) = " << sm << '\n';
    }
}

void
MacProj::set_outflow_bcs (int             level,
                          MultiFab*       mac_phi,
                          const MultiFab* /*u_mac*/,
                          const MultiFab& /*S*/,
                          const MultiFab& /*divu*/)
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
        BoxArray phiBoxArray(phiBoxList);
        phiBoxList.clear();

        //
        // Must do this kind of copy instead of mac_phi->copy(phidat);
        // because we're copying onto the ghost cells of the FABs,
        // not the valid regions.
        //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( int iface = 0; iface < nOutFlowTouched; ++iface )
        {
            for (MFIter mfi(*mac_phi); mfi.isValid(); ++mfi)
            {
                Box ovlp = (*mac_phi)[mfi].box() & phiBoxArray[iface];
                if (ovlp.ok()) {
                      (*mac_phi)[mfi].setVal<RunOn::Gpu>(0,ovlp,0,1);
                }
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

    for (long unsigned int i = 0; i < pirm.size(); i++)
    {
        const int dim = pirm[i].m_dim;

        BL_ASSERT(pirm[i].m_fbid.box() == pirm[i].m_srcBox);
        BL_ASSERT(pirm[i].m_srcBox.sameSize(pirm[i].m_dstBox));
        BL_ASSERT(u_mac[dim].DistributionMap()[pirm[i].m_idx] == ParallelDescriptor::MyProc());

        diff.resize(pirm[i].m_srcBox, 1);

        mfcd.FillFab(mfid[dim], pirm[i].m_fbid, diff);

        diff.minus<RunOn::Host>(u_mac[dim][pirm[i].m_idx],pirm[i].m_dstBox,diff.box(),0,0,1);

        const Real max_norm = diff.norm<RunOn::Host>(0);

        if (max_norm > umac_periodic_test_Tol )
        {
            amrex::Print() << "dir = "         << dim
                           << ", diff norm = " << max_norm
                           << " for region: "  << pirm[i].m_dstBox << std::endl;
            amrex::Error("Periodic bust in u_mac");
        }
    }
}
