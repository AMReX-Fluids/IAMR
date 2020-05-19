
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
    static int max_order = 4;
    static int agglomeration = 1;
    static int consolidation = 1;
    static int max_fmg_iter = 0;
    static int use_hypre = 0;
    static int hypre_verbose = 0;
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
        ppmacop.query("max_order", max_order);

        ParmParse ppmac("mac");
        ppmac.query("agglomeration", agglomeration);
        ppmac.query("consolidation", consolidation);
        ppmac.query("max_fmg_iter", max_fmg_iter);
#ifdef AMREX_USE_HYPRE
        ppmac.query("use_hypre", use_hypre);
        ppmac.query("hypre_verbose", hypre_verbose);
#endif

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

    //
    // Other MacProjector settings
    //
    macproj.getLinOp().setMaxOrder(max_order);
    if (use_hypre)
    {
        macproj.getMLMG().setBottomSolver(MLMG::BottomSolver::hypre);
        macproj.getMLMG().setBottomVerbose(hypre_verbose);
    }
    // Ask about this one
    macproj.getMLMG().setMaxFmgIter(max_fmg_iter);
    macproj.setVerbose(verbose);

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
        ParmParse ppmacop("macop");
        ppmacop.query("max_order", max_order);

        ParmParse ppmac("mac");
        ppmac.query("agglomeration", agglomeration);
        ppmac.query("consolidation", consolidation);
        ppmac.query("max_fmg_iter", max_fmg_iter);
#ifdef AMREX_USE_HYPRE
        ppmac.query("use_hypre", use_hypre);
        ppmac.query("hypre_verbose", hypre_verbose);
#endif

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
    mlabec.setMaxOrder(max_order);
    mlabec.setVerbose(verbose);

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
    if (use_hypre) {
      mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
      mlmg.setBottomVerbose(hypre_verbose);
      mlmg.setVerbose(hypre_verbose);
    }
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setVerbose(verbose);

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
