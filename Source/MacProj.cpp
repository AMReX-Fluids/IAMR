
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <MacProj.H>
#include <NavierStokesBase.H>
#include <OutFlowBC.H>
#include <hydro_MacProjector.H>

#ifdef AMREX_USE_EB
#include <hydro_ebgodunov.H>
#endif
#include <hydro_godunov.H>
#include <hydro_bds.H>


//fixme, for writesingle level plotfile
//#include<AMReX_PlotFileUtil.H>

using namespace amrex;

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
int  MacProj::check_umac_periodicity;
int  MacProj::max_order = 4;
int  MacProj::agglomeration = 1;
int  MacProj::consolidation = 1;
int  MacProj::max_fmg_iter = -1;

namespace
{
    Real umac_periodic_test_Tol;
}

void
MacProj::Initialize ()
{
    if (initialized) return;
    //
    // Set defaults here!!!
    //
    umac_periodic_test_Tol          = 1.e-10;
    MacProj::verbose                = 0;
    MacProj::mac_tol                = 1.0e-12;
    MacProj::mac_abs_tol            = 1.0e-16;
    MacProj::mac_sync_tol           = 1.0e-10;
    MacProj::do_outflow_bcs         = 1;
    //
    // Only check umac periodicity when debugging.  Can be overridden on input.
    //
#ifndef AMREX_DEBUG
    MacProj::check_umac_periodicity = 0;
#else
    MacProj::check_umac_periodicity = 1;
#endif

    // NOTE: IAMR uses a different max_order default than hydro::MacProjector,
    // which uses a default of 3
    static int max_order = 4;
    static int agglomeration = 1;
    static int consolidation = 1;
    static int max_fmg_iter = -1;


    ParmParse pp("mac_proj");

    pp.query("verbose",                verbose);
    pp.query("mac_tol",                mac_tol);
    pp.query("mac_abs_tol",            mac_abs_tol);
    pp.query("mac_sync_tol",           mac_sync_tol);
    pp.query("do_outflow_bcs",         do_outflow_bcs);
    pp.query("check_umac_periodicity", check_umac_periodicity);
    pp.query("umac_periodic_test_Tol", umac_periodic_test_Tol);

    pp.query("agglomeration", agglomeration);
    pp.query("consolidation", consolidation);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query( "maxorder"      , max_order );
#ifdef AMREX_USE_HYPRE
    if ( pp.contains("use_hypre") )
      amrex::Abort("use_hypre is no more. To use Hypre set mac_proj.bottom_solver = hypre.");
    if ( pp.contains("hypre_verbose") )
      amrex::Abort("hypre_verbose is no more. To make the bottom solver verbose set mac_proj.bottom_verbose = 1.");
#endif

    // Abort if old verbose flag is found
    if ( pp.countname("v") > 0 ) {
	amrex::Abort("mac_proj.v found in inputs. To set verbosity use mac_proj.verbose");
    }
    // Abort if old "mac." prefix is used.
    std::set<std::string> old_mac = ParmParse::getEntries("mac");
    if (!old_mac.empty()){
	Print()<<"All runtime options related to the mac projection now use 'mac_proj'.\n"
	       <<"Found these depreciated entries in the parameters list: \n";
	for ( auto param : old_mac ) {
	    Print()<<"  "<<param<<"\n";
	}
	amrex::Abort("Replace 'mac' prefix with 'mac_proj' in inputs");
    }

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
        mac_phi_crse.resize(finest_level+1);
        mac_reg.resize(finest_level+1);
    }

    LevelData[level] = level_data;

    if (level > 0)
    {
        mac_reg[level].reset(new FluxRegister(LevelData[level]->boxArray(),
                                              LevelData[level]->DistributionMap(),
                                              parent->refRatio(level-1),level,1));
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
    //  Set up the mac projection
    //
    const Real rhs_scale = 2.0/dt;

    const MultiFab* cphi = (level == 0) ? nullptr : mac_phi_crse[level-1].get();

    // Set bcoefs to the average of Density at the faces
    // In the EB case, they will be defined at the Face Centroid
    MultiFab rho(S.boxArray(),S.DistributionMap(), 1, S.nGrow(),
                 MFInfo(), (parent->getLevel(level)).Factory());
    MultiFab::Copy(rho, S, Density, 0, 1, S.nGrow()); // Extract rho component from S

    Array<MultiFab*,AMREX_SPACEDIM>  umac;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        umac[idim]= &(u_mac[idim]);

    Array<MultiFab*,AMREX_SPACEDIM>  fluxes;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
      // We don't need the fluxes, so nullptr means we won't compute them
      fluxes[idim]= nullptr;

    //
    // Perform projection
    //
    mlmg_mac_solve(parent, cphi, *phys_bc, density_math_bc, level,
		   mac_tol, mac_abs_tol, rhs_scale,
		   rho, divu, umac, mac_phi, fluxes);

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
            const Real mult = 1.0/Real(parent->nCycle(level));

            for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
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
    AMREX_ASSERT(level < finest_level);

    if (verbose) amrex::Print() << "... mac_sync_solve at level " << level << '\n';

    const Real      strt_time  = ParallelDescriptor::second();
    const BoxArray& grids      = LevelData[level]->boxArray();
    const DistributionMapping& dmap = LevelData[level]->DistributionMap();
    const Geometry& geom       = parent->Geom(level);
    const BoxArray& fine_boxes = LevelData[level+1]->boxArray();
    const NavierStokesBase& ns_level   = *(NavierStokesBase*) &(parent->getLevel(level));
    const MultiFab&     volume     = ns_level.Volume();
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
        AMREX_ASSERT(grids[mfi.index()].contains(mfi.tilebox()) );

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

    //
    // Compute Ucorr, including filling ghost cells
    //
    // null umac will prevent MacProject from computing div(umac)
    // and adding it to RHS
    //
    Array<MultiFab*,AMREX_SPACEDIM>  umac;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
      umac[idim]= nullptr;

    Rhs.negate();

    //
    // Perform projection
    //
    mlmg_mac_solve(parent, nullptr, *phys_bc, rho_math_bc, level,
                   mac_sync_tol, mac_abs_tol, rhs_scale,
                   rho_half, Rhs, umac, mac_sync_phi, Ucorr);

    for ( int idim=0; idim<AMREX_SPACEDIM; idim++)
    {
      //Ucorr = fluxes = -B grad phi
      // Make sure Ucorr has correct sign
      Ucorr[idim]->negate();
    }

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
void
MacProj::mac_sync_compute (int                   level,
                           Array<MultiFab*,AMREX_SPACEDIM>& Ucorr,
                           MultiFab&             Vsync,
                           MultiFab&             Ssync,
                           FluxRegister*         adv_flux_reg,
                           Vector<AdvectionForm>& advectionType,
                           Real                  prev_time,
                           Real                  dt,
                           int                   num_state_comps,
                           Real                  be_cn_theta,
                           int                   do_mom_diff,
                           bool                  update_fluxreg)
{
    //
    // Get parameters.
    //
    const BoxArray& grids               = LevelData[level]->boxArray();
    const DistributionMapping& dmap     = LevelData[level]->DistributionMap();
    NavierStokesBase&   ns_level        = *(NavierStokesBase*) &(parent->getLevel(level));
    const MultiFab* area                = ns_level.Area();

    // Only two options: do all state components (including velocity) or just the
    // velocity components. Anything else will likely cause problems with the
    // mac register update.
    AMREX_ASSERT(num_state_comps == AMREX_SPACEDIM || num_state_comps == ns_level.NUM_STATE);

    //
    // Prep MFs to store fluxes and edge states
    //
    Array<MultiFab, AMREX_SPACEDIM> fluxes;
    Array<MultiFab, AMREX_SPACEDIM> edgestate;
    const int  nghost  = 0;
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        const BoxArray& ba = LevelData[level]->getEdgeBoxArray(i);
        fluxes[i].define(ba, dmap, num_state_comps, nghost, MFInfo(),ns_level.Factory());
        edgestate[i].define(ba, dmap, num_state_comps, nghost, MFInfo(), ns_level.Factory());
    }

    MultiFab visc_terms(grids,dmap,num_state_comps,ns_level.nghost_force(),
                        MFInfo(),ns_level.Factory());
    FillPatchIterator S_fpi(ns_level,visc_terms,ns_level.nghost_state(),
                            prev_time,State_Type,0,num_state_comps);
    MultiFab& Smf = S_fpi.get_mf();

    // Get density -- Smf might or might not contain rho.
    FillPatchIterator rho_fpi(ns_level,visc_terms,ns_level.nghost_state(),
                              prev_time,State_Type,Density,1);
    MultiFab& rhoMF = rho_fpi.get_mf();

    // Store momenta multifab if conservative approach is used.
    MultiFab momenta;
    if  (do_mom_diff == 1)
    {
        momenta.define(grids,dmap, AMREX_SPACEDIM, Smf.nGrow(), MFInfo(), Smf.Factory());
        MultiFab::Copy(momenta,Smf,0,0,AMREX_SPACEDIM, Smf.nGrow());
        for (int d=0; d < AMREX_SPACEDIM; ++d )
            MultiFab::Multiply( momenta, rhoMF, 0, d, 1, Smf.nGrow());
    }

    std::unique_ptr<MultiFab> forcing_term;
    std::unique_ptr<MultiFab> divu_fp;


    //
    // Compute the mac sync correction.
    //
    if ( ns_level.advection_scheme == "Godunov_PLM" ||
         ns_level.advection_scheme == "Godunov_PPM" ||
         ns_level.advection_scheme == "BDS" )
    {

        forcing_term.reset(new MultiFab(grids, dmap, num_state_comps, ns_level.nghost_force()));
        divu_fp.reset(ns_level.getDivCond(ns_level.nghost_force(),prev_time));

        MultiFab& Gp = ns_level.get_old_data(Gradp_Type);

        visc_terms.setVal(0.0); // Initialize to make calls below safe

        // Get viscous forcing.
        if (be_cn_theta != 1.0)
        {
            ns_level.getViscTerms(visc_terms,0,num_state_comps,prev_time);
        }

        //
        // Compute forcing terms
        //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter Smfi(Smf,TilingIfNotGPU()); Smfi.isValid(); ++Smfi)
        {
            auto const gbx = Smfi.growntilebox(ns_level.nghost_force());

            //
            // Compute total forcing terms.
            //
            ns_level.getForce(forcing_term->get(Smfi),gbx,0,num_state_comps,
                              prev_time,Smf[Smfi],rhoMF[Smfi],0,Smfi);

            for (int comp = 0; comp < num_state_comps; ++comp)
            {
                auto const& tf    = forcing_term->array(Smfi,comp);

                // This part really should be a function within NS/NSB, since it relies on how
                // things were compute there ...
                if (comp < AMREX_SPACEDIM)  // Velocity/Momenta
                {
                    auto const& visc = visc_terms[Smfi].const_array(comp);
                    auto const& gp   = Gp[Smfi].const_array(comp);

                    if ( do_mom_diff == 0 )
                    {
                        auto const& rho   = rhoMF[Smfi].const_array();

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
                else  // Scalars. Reconstruct forcing terms as in scalar_advection
                {
                    auto const& visc = visc_terms[Smfi].const_array(comp);
                    auto const& rho = rhoMF[Smfi].const_array();

                    if ( NavierStokesBase::do_temp && comp==NavierStokesBase::Temp )
                    {
                        //
                        // Solving
                        //   dT/dt + U dot del T = ( del dot lambda grad T + H_T ) / (rho c_p)
                        // with tforces = H_T/c_p (since it's always density-weighted), and
                        // visc = del dot mu grad T, where mu = lambda/c_p
                        //
                        amrex::ParallelFor(gbx, [tf, visc, rho]
                        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        { tf(i,j,k) = ( tf(i,j,k) + visc(i,j,k) ) / rho(i,j,k); });
                    }
                    else
                    {
                        if (advectionType[comp] == Conservative)
                        {
                            //
                            // For tracers, Solving
                            //   dS/dt + del dot (U S) = del dot beta grad (S/rho) + rho H_q
                            // where S = rho q, q is a concentration
                            // tforces = rho H_q (since it's always density-weighted)
                            // visc = del dot beta grad (S/rho)
                            //
                            amrex::ParallelFor(gbx, [tf, visc]
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
                            amrex::ParallelFor(gbx, [tf, visc, rho]
                            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                            { tf(i,j,k) = tf(i,j,k) / rho(i,j,k) + visc(i,j,k); });
                        }
                    }
                }
            }
        }

    }
    else
    {
        Abort("MacProj::mac_sync_compute: Unkown adveciton scheme");
    }

        //
        // Perform sync
        //

    //
    // Do velocity sync first
    //
    ns_level.ComputeAofs(Vsync, /*Vsync_comp*/ 0, /*state_indx*/ 0,
                         /*num comps*/ AMREX_SPACEDIM,
                         /*State*/ (do_mom_diff == 1) ? momenta : Smf, /*S_comp*/ 0,
                         forcing_term.get(), /*forcing_term_comp*/ 0,
                         divu_fp.get(),
                         fluxes, /*flux_comp*/ 0,
                         edgestate, /*edge_comp*/ 0, /*known_edgestate*/ false,
                         /*is_velocity*/ true, dt,
                         /*is_sync*/ true, Ucorr);

    if (num_state_comps > AMREX_SPACEDIM)
            {
        //
        // Also compute scalar sync. Recall that it must be all the scalars here.
        //
        ns_level.ComputeAofs(Ssync, /*Ssync_comp*/ 0, /* state_indx*/ Density,
                             /*num comps*/ num_state_comps - AMREX_SPACEDIM,
                             /*State*/ Smf, /*S_comp*/ Density,
                             forcing_term.get(), /*forcing_term_comp*/ Density,
                             divu_fp.get(),
                             fluxes, /*flux_comp*/ Density,
                             edgestate, /*edge_comp*/ Density, /*known_edgestate*/ false,
                             /*is_velocity*/ false, dt,
                             /*is_sync*/ true, Ucorr);
    }


    if (level > 0 && update_fluxreg)
    {
        const Real mlt =  -1.0/Real(parent->nCycle(level));
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            for (int comp = 0; comp < num_state_comps; ++comp)
            {
                    adv_flux_reg->FineAdd(fluxes[d],d,comp,comp,1,-dt);
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
                           int                    state_comp, // index in AmrLevel.state
                           int                    Sync_indx,
                           MultiFab* const*       edgestate,
                           int                    edge_comp,
                           FluxRegister*          adv_flux_reg,
                           Real                   dt,
                           bool                   update_fluxreg)
{
    const DistributionMapping& dmap     = LevelData[level]->DistributionMap();
    NavierStokesBase& ns_level   = *(NavierStokesBase*) &(parent->getLevel(level));

    int ncomp = 1;

    // Put edgestate in desired containter
    Array<MultiFab,AMREX_SPACEDIM> edges{
        AMREX_D_DECL(MultiFab(*edgestate[0], amrex::make_alias, edge_comp, ncomp),
                     MultiFab(*edgestate[1], amrex::make_alias, edge_comp, ncomp),
                     MultiFab(*edgestate[2], amrex::make_alias, edge_comp, ncomp) ) };

    Array<MultiFab, AMREX_SPACEDIM> fluxes;
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        const BoxArray& ba = LevelData[level]->getEdgeBoxArray(i);
        fluxes[i].define(ba, dmap, ncomp, edgestate[0]->nGrow(), MFInfo(),ns_level.Factory());
    }

    //
    // Compute the mac sync correction.
    //
    ns_level.ComputeAofs(Sync, /*Ssync_comp*/ Sync_indx, state_comp,
                         ncomp,
                         /*State*/ MultiFab(), /*S_comp*/ int(),//not used when known_edgestates
                         /*forcing*/ nullptr, /*f_comp*/ int(), //not used when known_edgestates
                         /*constraint divU*/ nullptr,           //not used when known_edgestates
                         fluxes, /*flux_comp*/ 0,
                         edges, 0, /*known_edgestate*/ true,
                         /*is_velocity*/ false, dt,
                         /*is_sync*/ true, Ucorr);

    if (level > 0 && update_fluxreg)
    {
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            adv_flux_reg->FineAdd(fluxes[d],d,0,state_comp,1,-dt);
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
    Orientation outFaces[2*AMREX_SPACEDIM];
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
    MultiFabId             mfid[AMREX_SPACEDIM];
    std::vector<TURec>     pirm;
    Vector<IntVect>         pshifts(27);
    std::vector< std::pair<int,Box> > isects;


    for (int dim = 0; dim < AMREX_SPACEDIM; dim++)
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

        AMREX_ASSERT(pirm[i].m_fbid.box() == pirm[i].m_srcBox);
        AMREX_ASSERT(pirm[i].m_srcBox.sameSize(pirm[i].m_dstBox));
        AMREX_ASSERT(u_mac[dim].DistributionMap()[pirm[i].m_idx] == ParallelDescriptor::MyProc());

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

// project
void
MacProj::mlmg_mac_solve (Amr* a_parent, const MultiFab* cphi, const BCRec& a_phys_bc,
			 const BCRec& density_math_bc,
			 int level, Real a_mac_tol, Real a_mac_abs_tol, Real rhs_scale,
			 const MultiFab &rho, const MultiFab &Rhs,
			 Array<MultiFab*,AMREX_SPACEDIM>& u_mac, MultiFab *mac_phi,
			 Array<MultiFab*,AMREX_SPACEDIM>& fluxes)
{
    const Geometry& geom = a_parent->Geom(level);
    const BoxArray& ba = Rhs.boxArray();
    const DistributionMapping& dm = Rhs.DistributionMap();

    //
    // Compute beta coefficients
    //
    Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM> bcoefs;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        BoxArray nba = amrex::convert(ba,IntVect::TheDimensionVector(idim));
        bcoefs[idim].reset(new  MultiFab(nba, dm, 1, 0, MFInfo(),
					 (a_parent->getLevel(level)).Factory()));
    }

    //
    // Set bcoefs to the average of Density at the faces
    // In the EB case, they will be defined at the Face Centroid
    //
#ifdef AMREX_USE_EB
    EB_interp_CellCentroid_to_FaceCentroid( rho, GetArrOfPtrs(bcoefs), 0, 0, 1,
					    geom, {density_math_bc});
#else
    amrex::ignore_unused(density_math_bc);
    average_cellcenter_to_face(GetArrOfPtrs(bcoefs), rho, geom);
#endif

    //
    // Now invert the coefficients and apply scale factor
    //
    int ng_for_invert(0);
    Real scale_factor(1.0/rhs_scale);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        bcoefs[idim]->invert(scale_factor,ng_for_invert);
        bcoefs[idim]->FillBoundary( geom.periodicity() );
    }

    //
    // Create MacProjector Object
    //
    LPInfo info;
    int max_coarsening_level(100);
    ParmParse pp("mac_proj");
    pp.query("mg_max_coarsening_level", max_coarsening_level);

    info.setMaxCoarseningLevel(max_coarsening_level);
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);

    //
    // To use phi on CellCentroids, must also call
    // macproj.get_linop().setEBDirichlet (int amrlev, const MultiFab& phi, const MultiFab& beta)
    // or
    // macproj.get_linop().setEBHomogDirichlet (int amrlev, const MultiFab& beta)
    //
    // Location information is not used for non-EB
    //
    Hydro::MacProjector macproj( {u_mac}, MLMG::Location::FaceCentroid, // Location of umac (face center vs centroid)
                                {GetArrOfConstPtrs(bcoefs)}, MLMG::Location::FaceCentroid,  // Location of beta (face center vs centroid)
                                MLMG::Location::CellCenter,           // Location of solution variable phi (cell center vs centroid)
                                {geom}, info,
                                {&Rhs}, MLMG::Location::CellCentroid);  // Location of RHS (cell center vs centroid)

    //
    // Set BCs
    //
    std::array<MLLinOp::BCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<MLLinOp::BCType,AMREX_SPACEDIM> mlmg_hibc;
    set_mac_solve_bc(mlmg_lobc, mlmg_hibc, a_phys_bc, geom);

    macproj.setDomainBC(mlmg_lobc, mlmg_hibc);
    if (level > 0 && cphi)
    {
        macproj.setCoarseFineBC(cphi, a_parent->refRatio(level-1)[0]);
    }
    macproj.setLevelBC(0, mac_phi);

    // MacProj default max order is 3. Here we use a default of 4, so must
    // call setMaxOrder to overwrite MacProj default.
    macproj.getLinOp().setMaxOrder(max_order);
    if ( max_fmg_iter > -1 )
      macproj.getMLMG().setMaxFmgIter(max_fmg_iter);

    //
    // Perform projection
    //
    macproj.project({mac_phi}, a_mac_tol, a_mac_abs_tol);

    if ( fluxes[0] )
      // fluxes = -B grad phi
      macproj.getFluxes({fluxes}, {mac_phi}, MLMG::Location::FaceCentroid);
}

void
MacProj::set_mac_solve_bc (Array<MLLinOp::BCType,AMREX_SPACEDIM>& mlmg_lobc,
			   Array<MLLinOp::BCType,AMREX_SPACEDIM>& mlmg_hibc,
			   const BCRec& a_phys_bc, const Geometry& geom)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom.isPeriodic(idim)) {
            mlmg_lobc[idim] = MLLinOp::BCType::Periodic;
            mlmg_hibc[idim] = MLLinOp::BCType::Periodic;
        } else {
            if (a_phys_bc.lo(idim) == Outflow) {
                mlmg_lobc[idim] = MLLinOp::BCType::Dirichlet;
            } else {
                mlmg_lobc[idim] = MLLinOp::BCType::Neumann;
            }
            if (a_phys_bc.hi(idim) == Outflow) {
                mlmg_hibc[idim] = MLLinOp::BCType::Dirichlet;
            } else {
                mlmg_hibc[idim] = MLLinOp::BCType::Neumann;
            }
        }
    }
}
