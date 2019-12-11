
#include <AMReX_MacProjector.H>
#include <AMReX_MLMG.H>
#include <AMReX_ParmParse.H>

#include <MacOpMacDrivers.H>

#include <IAMR_MLMG_F.H>

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

static void compute_mac_coefficient (std::array<MultiFab,AMREX_SPACEDIM>& bcoefs,
                                     const MultiFab& rho, int rho_comp, Real scale)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(rho,true); mfi.isValid(); ++mfi)
    {
        AMREX_D_TERM(const Box& xbx = mfi.nodaltilebox(0);,
                     const Box& ybx = mfi.nodaltilebox(1);,
                     const Box& zbx = mfi.nodaltilebox(2););
        iamr_mac_coef(AMREX_D_DECL(BL_TO_FORTRAN_BOX(xbx),
                                   BL_TO_FORTRAN_BOX(ybx),
                                   BL_TO_FORTRAN_BOX(zbx)),
                      AMREX_D_DECL(BL_TO_FORTRAN_ANYD(bcoefs[0][mfi]),
                                   BL_TO_FORTRAN_ANYD(bcoefs[1][mfi]),
                                   BL_TO_FORTRAN_ANYD(bcoefs[2][mfi])),
                      BL_TO_FORTRAN_N_ANYD(rho[mfi],rho_comp),
                      &scale);
    }
}

#ifdef AMREX_USE_EB
  static void compute_mac_coefficient_eb (std::array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM>& bcoefs,
					  const MultiFab& rho, int rho_comp, Real scale)
{
  BL_PROFILE("MacProjection::compute_mac_coefficient_eb");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(rho,true); mfi.isValid(); ++mfi)
    {
      // Boxes for staggered components
      AMREX_D_TERM(const Box& xbx = mfi.nodaltilebox(0);,
		   const Box& ybx = mfi.nodaltilebox(1);,
		   const Box& zbx = mfi.nodaltilebox(2););

      // this is to check efficiently if this tile contains any eb stuff
      // need cc fab with ebfactory...
      const EBFArrayBox&  ffab = static_cast<EBFArrayBox const&>(rho[mfi]);
      const EBCellFlagFab&  flags = ffab.getEBCellFlagFab();

      if (flags.getType(amrex::grow(mfi.tilebox(),0)) == FabType::covered )
      {
	D_TERM(bcoefs[0]->setVal( 1.2345e300, xbx, 0, 1);,
	       bcoefs[1]->setVal( 1.2345e300, ybx, 0, 1);,
	       bcoefs[2]->setVal( 1.2345e300, zbx, 0, 1));
      }
      else
      {
          iamr_mac_coef(AMREX_D_DECL(BL_TO_FORTRAN_BOX(xbx),
				     BL_TO_FORTRAN_BOX(ybx),
				     BL_TO_FORTRAN_BOX(zbx)),
			AMREX_D_DECL(BL_TO_FORTRAN_ANYD((*(bcoefs[0]))[mfi]),
				     BL_TO_FORTRAN_ANYD((*(bcoefs[1]))[mfi]),
				     BL_TO_FORTRAN_ANYD((*(bcoefs[2]))[mfi])),
			BL_TO_FORTRAN_N_ANYD(rho[mfi],rho_comp),
			&scale);
      }
    }
    //fixme
    // VisMF::Write(*(bcoefs[0]),"bx");
    // VisMF::Write(*(bcoefs[1]),"by");
    // VisMF::Write(*(bcoefs[2]),"bz");
}
#endif

static void compute_mac_rhs (MultiFab& rhs, const MultiFab* umac,
                             const MultiFab* area, const MultiFab& volume, const Real* dxinv)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(rhs,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        iamr_mac_rhs(BL_TO_FORTRAN_BOX(bx),
                     BL_TO_FORTRAN_ANYD(rhs[mfi]),
                     AMREX_D_DECL(BL_TO_FORTRAN_ANYD(umac[0][mfi]),
                                  BL_TO_FORTRAN_ANYD(umac[1][mfi]),
                                  BL_TO_FORTRAN_ANYD(umac[2][mfi])),
#if (AMREX_SPACEDIM == 2)
                     BL_TO_FORTRAN_ANYD(area[0][mfi]),
                     BL_TO_FORTRAN_ANYD(area[1][mfi]),
                     BL_TO_FORTRAN_ANYD(volume[mfi]));
#else
                     dxinv);
#endif
    }
}

}

void mlmg_mac_level_solve (Amr* parent, const MultiFab* cphi, const BCRec& phys_bc,
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
    // This does NOT compute the average at the CENTROIDS, but at the faces
    // Operator inside MAc projector will take care of this.
    MultiFab rho(S.boxArray(),S.DistributionMap(), 1, S.nGrow());
    MultiFab::Copy(rho, S, Density, 0, 1, S.nGrow()); // Extract rho component from S
    average_cellcenter_to_face( GetArrOfPtrs(bcoefs), rho, geom);

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

    MacProjector macproj( {umac}, {GetArrOfConstPtrs(bcoefs)}, {geom}, info, {&Rhs} );

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
    macproj.project({mac_phi}, mac_tol, mac_abs_tol, MLMG::Location::FaceCentroid);
}



void mlmg_mac_sync_solve (Amr* parent, const BCRec& phys_bc,
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
    // This does NOT compute the average at the CENTROIDS, but at the faces
    // Operator inside MAc projector will take care of this.
    average_cellcenter_to_face( GetArrOfPtrs(bcoefs), rho, geom);

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
    //
    // Create right containers for MacProjector
    //Array<MultiFab*,AMREX_SPACEDIM>  a_umac({mac_phi});
    //
    // Create MacProjection object
    // make sure u_mac already has bndry properly filled
    //MacProjector macproj( {a_umac}, {GetArrOfConstPtrs(bcoefs)}, {geom}, info, {&Rhs} );

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
    // weiqun says not needed... I think default is the 0 we want
    // if (level > 0) {
    //     mlabec.setCoarseFineBC(nullptr, parent->refRatio(level-1)[0]);
    // }
    mlabec.setLevelBC(0, mac_phi);

    mlabec.setScalars(0.0, 1.0);
    mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(bcoefs));

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
    mlmg.getFluxes({Ucorr});

    // // Solve with initial guess of zero
    // macproj.project(mac_tol,mac_abs_tol,MLMG::Location::FaceCentroid);

    // // Enforce perodicity
    // for (int i=0; i<AMREX_SPACEDIM; ++i)
    //   mac_phi[i].FillBoundary(geom.periodicity());

    //if (verbose)  // Always print this for now
//     {
//       MultiFab divu(ba, dm, 1, 0, MFInfo(), (parent->getLevel(level)).Factory());
// #ifdef AMREX_USE_EB
//       bool already_on_centroid(true);
//       EB_computeDivergence(divu,GetArrOfConstPtrs(a_umac),geom,already_on_centroid);
// #else
//       computeDivergence(divu,GetArrOfConstPtrs(a_umac),geom);
// #endif

//       Print() << "  MAC SYNC solve on level "<< level
//     	      << " BEFORE projection: max(abs(divu)) = " << divu.norm0() << "\n";
//     }
}
