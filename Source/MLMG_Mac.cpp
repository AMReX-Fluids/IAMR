
#include <AMReX_MacProjector.H>
#include <AMReX_MLMG.H>
#include <AMReX_ParmParse.H>

#include <MacOpMacDrivers.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MLEBABecLap.H>
//fixme
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_VisMF.H>
#else
#include <AMReX_MLABecLaplacian.H>
#endif

using namespace amrex;

namespace {
    static int initialized = false;
    // NOTE: IAMR uses a different max_order default than amrex::MacProjector,
    // which uses a default of 3
    static int max_order = 4;
    static int agglomeration = 1;
    static int consolidation = 1;
    static int max_fmg_iter = -1;

    // Likely want these to match defaults in MacProjector
    static int          bottom_verbose(0);
    static int          maxiter(200);
    static int          bottom_maxiter(200);
    static Real         bottom_rtol(1.0e-4_rt);
    static Real         bottom_atol(-1.0_rt);
    static std::string  bottom_solver("bicg");
    static int num_pre_smooth(2);
    static int num_post_smooth(2);
}

namespace {

static void set_mac_solve_bc (std::array<MLLinOp::BCType,AMREX_SPACEDIM>& mlmg_lobc,
                              std::array<MLLinOp::BCType,AMREX_SPACEDIM>& mlmg_hibc,
                              const BCRec& phys_bc, const Geometry& geom)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom.isPeriodic(idim)) {
            mlmg_lobc[idim] = MLLinOp::BCType::Periodic;
            mlmg_hibc[idim] = MLLinOp::BCType::Periodic;
        } else {
            if (phys_bc.lo(idim) == Outflow) {
                mlmg_lobc[idim] = MLLinOp::BCType::Dirichlet;
            } else {
                mlmg_lobc[idim] = MLLinOp::BCType::Neumann;
            }
            if (phys_bc.hi(idim) == Outflow) {
                mlmg_hibc[idim] = MLLinOp::BCType::Dirichlet;
            } else {
                mlmg_hibc[idim] = MLLinOp::BCType::Neumann;
            }
        }
    }
}

}

void mlmg_mac_level_solve (Amr* parent, const MultiFab* cphi, const BCRec& phys_bc,
                           const BCRec& density_math_bc,
                           int level, int Density, Real mac_tol, Real mac_abs_tol, Real rhs_scale,
                           const MultiFab &S, MultiFab &Rhs,
                           MultiFab *u_mac, MultiFab *mac_phi, int verbose)
{
    if (!initialized)
    {
        ParmParse ppmacop("macop");
        if ( ppmacop.contains("max_order") )
	  amrex::Abort("macop.max_order is no more. Please use mac_proj.maxorder.");

        ParmParse ppmac("mac");
        ppmac.query("agglomeration", agglomeration);
        ppmac.query("consolidation", consolidation);
        ppmac.query("max_fmg_iter", max_fmg_iter);
#ifdef AMREX_USE_HYPRE
        if ( ppmac.contains("use_hypre") )
	  amrex::Abort("use_hypre is no more. To use Hypre set mac_proj.bottom_solver = hypre.");
        if ( ppmac.contains("hypre_verbose") )
	  amrex::Abort("hypre_verbose is no more. To make the bottom solver verbose set mac_proj.bottom_verbose = 1.");
#endif

	// Read the mac_proj options so we can set them here for use in the mac_sync
	// which does not go to the MacProjector.
	ParmParse pp("mac_proj");
	pp.query( "verbose"       , verbose );
	pp.query( "maxorder"      , max_order );
	pp.query( "bottom_verbose", bottom_verbose );
	pp.query( "maxiter"       , maxiter );
	pp.query( "bottom_maxiter", bottom_maxiter );
	pp.query( "bottom_rtol"   , bottom_rtol );
	pp.query( "bottom_atol"   , bottom_atol );
	pp.query( "bottom_solver" , bottom_solver );

	pp.query( "num_pre_smooth"  , num_pre_smooth );
	pp.query( "num_post_smooth" , num_post_smooth );

        initialized = true;
    }

    const Geometry& geom = parent->Geom(level);
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
                                            (parent->getLevel(level)).Factory()));
    }

    // Set bcoefs to the average of Density at the faces
    // In the EB case, they will be defined at the Face Centroid
    MultiFab rho(S.boxArray(),S.DistributionMap(), 1, S.nGrow(),
                 MFInfo(), (parent->getLevel(level)).Factory());
    MultiFab::Copy(rho, S, Density, 0, 1, S.nGrow()); // Extract rho component from S

#ifdef AMREX_USE_EB
        EB_interp_CellCentroid_to_FaceCentroid( rho, GetArrOfPtrs(bcoefs), 0, 0, 1,
                                                geom, {density_math_bc});
#else
        average_cellcenter_to_face(GetArrOfPtrs(bcoefs), rho, geom);
#endif

    // Now invert the coefficients and apply scale factor
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
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);

    Array<MultiFab*,AMREX_SPACEDIM>  umac;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        umac[idim]= &(u_mac[idim]);

    // phi is always considered at cell centers.
    // This won't be an issue because  we won't use phi directly
    // and EB is always, by construction, far from fine/coarse interfaces
    MacProjector macproj( {umac}, MLMG::Location::FaceCentroid, // Location of umac (face center vs centroid)
                          {GetArrOfConstPtrs(bcoefs)}, MLMG::Location::FaceCentroid,  // Location of beta (face center vs centroid)
                          MLMG::Location::CellCenter,           // Location of solution variable phi (cell center vs centroid)
                          {geom}, info,
                          {&Rhs}, MLMG::Location::CellCentroid);  // Location of RHS (cell center vs centroid)

    //
    // Set BCs
    //
    std::array<MLLinOp::BCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<MLLinOp::BCType,AMREX_SPACEDIM> mlmg_hibc;
    set_mac_solve_bc(mlmg_lobc, mlmg_hibc, phys_bc, geom);

    macproj.setDomainBC(mlmg_lobc, mlmg_hibc);
    if (level > 0)
    {
        macproj.getLinOp().setCoarseFineBC(cphi, parent->refRatio(level-1)[0]);
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
    macproj.project({mac_phi}, mac_tol, mac_abs_tol);
}



void mlmg_mac_sync_solve (Amr* parent, const BCRec& phys_bc,
                          const BCRec& rho_math_bc,
                          int level, Real mac_tol, Real mac_abs_tol, Real rhs_scale,
                          const MultiFab* area, const MultiFab& volume,
                          const MultiFab& rho, MultiFab& Rhs,
                          MultiFab* mac_phi, Array<MultiFab*,AMREX_SPACEDIM>& Ucorr,
			  int verbose)
{
    if (!initialized)
    {
        // Normal code execution will not get here. mlmg_mac_level_solve will get called
        // first and set these variables. Leave this here for debugging purposes.

        ParmParse ppmac("mac");
        ppmac.query("agglomeration", agglomeration);
        ppmac.query("consolidation", consolidation);
        ppmac.query("max_fmg_iter", max_fmg_iter);

	ParmParse pp("mac_proj");
	pp.query( "verbose"       , verbose );
	pp.query( "maxorder"      , max_order );
	pp.query( "bottom_verbose", bottom_verbose );
	pp.query( "maxiter"       , maxiter );
	pp.query( "bottom_maxiter", bottom_maxiter );
	pp.query( "bottom_rtol"   , bottom_rtol );
	pp.query( "bottom_atol"   , bottom_atol );
	pp.query( "bottom_solver" , bottom_solver );

	pp.query( "num_pre_smooth"  , num_pre_smooth );
	pp.query( "num_post_smooth" , num_post_smooth );

        initialized = true;
    }

    const Geometry& geom = parent->Geom(level);
    const BoxArray& ba = Rhs.boxArray();
    const DistributionMapping& dm = Rhs.DistributionMap();

    // no need to set A coef because it's zero
    Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM> bcoefs;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        BoxArray nba = amrex::convert(ba,IntVect::TheDimensionVector(idim));
        bcoefs[idim].reset(new  MultiFab(nba, dm, 1, 0, MFInfo(),
                                            (parent->getLevel(level)).Factory()));
    }

    //
    // Want to go back to ML(EB)ABecLap. MacProjector takes div and don't need that here
    // Will MLEBABecLap take care of bcoef like MacProj?
    // Not so relevant at the moment with restrriction EB not crossing coarse-fine bndry

    // Set bcoefs to the average of Density at the faces
    // In the EB case, they will be defined at the Face Centroid
#ifdef AMREX_USE_EB
        EB_interp_CellCentroid_to_FaceCentroid( rho, GetArrOfPtrs(bcoefs), 0, 0, 1,
                                                geom, {rho_math_bc});
#else
        average_cellcenter_to_face(GetArrOfPtrs(bcoefs), rho, geom);
#endif


    // Now invert the coefficients and apply scale factor
    int ng_for_invert(0);
    Real scale_factor(1.0/rhs_scale);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        bcoefs[idim]->invert(scale_factor,ng_for_invert);
        bcoefs[idim]->FillBoundary( geom.periodicity() );
    }

    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);

    // MacProjector is not what we want for sync
    // takes Face-centered velocity and interally computes div, but already have div

#ifdef AMREX_USE_EB
    const auto& ebf = &dynamic_cast<EBFArrayBoxFactory const&>((parent->getLevel(level)).Factory());
    MLEBABecLap mlabec({geom}, {ba}, {dm}, info, {ebf});
#else
    MLABecLaplacian mlabec({geom}, {ba}, {dm}, info);
#endif

    std::array<MLLinOp::BCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<MLLinOp::BCType,AMREX_SPACEDIM> mlmg_hibc;
    set_mac_solve_bc(mlmg_lobc, mlmg_hibc, phys_bc, geom);

    mlabec.setDomainBC(mlmg_lobc, mlmg_hibc);
    mlabec.setLevelBC(0, mac_phi);
    mlabec.setScalars(0.0, 1.0);
#ifdef AMREX_USE_EB
    mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(bcoefs), MLMG::Location::FaceCentroid);
#else
    mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(bcoefs));
#endif

    MLMG mlmg(mlabec);
    // Set options
    mlabec.setMaxOrder(max_order);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);
    mlmg.setMaxIter(maxiter);
    mlmg.setBottomMaxIter(bottom_maxiter);
    mlmg.setBottomTolerance(bottom_rtol);
    mlmg.setBottomToleranceAbs(bottom_atol);
    mlmg.setPreSmooth(num_pre_smooth);
    mlmg.setPostSmooth(num_post_smooth);
    if ( max_fmg_iter > -1 )
      mlmg.setMaxFmgIter(max_fmg_iter);

    if (bottom_solver == "smoother")
    {
        mlmg.setBottomSolver(MLMG::BottomSolver::smoother);
    }
    else if (bottom_solver == "bicg")
    {
        mlmg.setBottomSolver(MLMG::BottomSolver::bicgstab);
    }
    else if (bottom_solver == "cg")
    {
        mlmg.setBottomSolver(MLMG::BottomSolver::cg);
    }
    else if (bottom_solver == "bicgcg")
    {
        mlmg.setBottomSolver(MLMG::BottomSolver::bicgcg);
    }
    else if (bottom_solver == "cgbicg")
    {
        mlmg.setBottomSolver(MLMG::BottomSolver::cgbicg);
    }
    else if (bottom_solver == "hypre")
    {
#ifdef AMREX_USE_HYPRE
        mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
#else
        amrex::Abort("AMReX was not built with HYPRE support");
#endif
    }

    Rhs.negate();

    mlmg.setFinalFillBC(true);
    mlmg.solve({mac_phi}, {&Rhs}, mac_tol, mac_abs_tol);

    //Ucorr = fluxes = -B grad phi
    mlmg.getFluxes({Ucorr}, MLMG::Location::FaceCentroid);

    for ( int idim=0; idim<AMREX_SPACEDIM; idim++)
    {
      // Make sure Ucorr has correct sign
      Ucorr[idim]->negate();

#ifdef AMREX_USE_EB
      // Set all ghost nodes to zero and then swap ghost nodes with neighboring boxes
      // This will set correction to zero at a box boundary unless the boundary is
      // of periodic type or it is shared with another box
      Ucorr[idim] -> setBndry(0.0);
      Ucorr[idim] -> FillBoundary(geom.periodicity());
#endif
    }
}
