
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_ParmParse.H>

#include <MacOpMacDrivers.H>

#include <IAMR_MLMG_F.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MacProjector.H>
//fixme
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_VisMF.H>
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
                           const MultiFab *area, const MultiFab &volume,
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

    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);

    MLABecLaplacian mlabec({geom}, {ba}, {dm}, info);
    mlabec.setMaxOrder(max_order);

    std::array<MLLinOp::BCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<MLLinOp::BCType,AMREX_SPACEDIM> mlmg_hibc;
    set_mac_solve_bc(mlmg_lobc, mlmg_hibc, phys_bc, geom);

    mlabec.setDomainBC(mlmg_lobc, mlmg_hibc);
    if (level > 0) {
        mlabec.setCoarseFineBC(cphi, parent->refRatio(level-1)[0]);
    }
    mlabec.setLevelBC(0, mac_phi);

    mlabec.setScalars(0.0, 1.0);

    // no need to set A coef because it's zero
     
    std::array<MultiFab,AMREX_SPACEDIM> bcoefs;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        const BoxArray& nba = amrex::convert(ba, IntVect::TheDimensionVector(idim));
        bcoefs[idim].define(nba, dm, 1, 0);
    }
    compute_mac_coefficient(bcoefs, S, Density, 1.0/rhs_scale);
    mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(bcoefs));

    MLMG mlmg(mlabec);
    if (use_hypre) {
        mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
        mlmg.setBottomVerbose(hypre_verbose);
    }
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setVerbose(verbose);

    const Real* dxinv = geom.InvCellSize();
    compute_mac_rhs(Rhs, u_mac, area, volume, dxinv);

    mlmg.solve({mac_phi}, {&Rhs}, mac_tol, mac_abs_tol);

    auto& fluxes = bcoefs;
    mlmg.getFluxes({amrex::GetArrOfPtrs(fluxes)});
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        MultiFab::Add(u_mac[idim], fluxes[idim], 0, 0, 1, 0);
    }
}

void mlmg_mac_sync_solve (Amr* parent, const BCRec& phys_bc,
                          int level, Real mac_tol, Real mac_abs_tol, Real rhs_scale,
                          const MultiFab* area, const MultiFab& volume,
                          const MultiFab& rho, MultiFab& Rhs,
                          MultiFab* mac_phi, int verbose)
{
    if (!initialized) 
    {
        ParmParse ppmacop("macop");
        ppmacop.query("max_order", max_order);
        
        ParmParse ppmac("mac");
        ppmac.query("agglomeration", agglomeration);
        ppmac.query("consolidation", consolidation);
        ppmac.query("max_fmg_iter", max_fmg_iter);
        
        initialized = true;
    }

    const Geometry& geom = parent->Geom(level);
    const BoxArray& ba = Rhs.boxArray();
    const DistributionMapping& dm = Rhs.DistributionMap();

    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);

    MLABecLaplacian mlabec({geom}, {ba}, {dm}, info);
    mlabec.setMaxOrder(max_order);

    std::array<MLLinOp::BCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<MLLinOp::BCType,AMREX_SPACEDIM> mlmg_hibc;
    set_mac_solve_bc(mlmg_lobc, mlmg_hibc, phys_bc, geom);

    mlabec.setDomainBC(mlmg_lobc, mlmg_hibc);
    if (level > 0) {
        mlabec.setCoarseFineBC(nullptr, parent->refRatio(level-1)[0]);
    }
    mlabec.setLevelBC(0, mac_phi);

    mlabec.setScalars(0.0, 1.0);

    // no need to set A coef because it's zero
     
    std::array<MultiFab,AMREX_SPACEDIM> bcoefs;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        const BoxArray& nba = amrex::convert(ba, IntVect::TheDimensionVector(idim));
        bcoefs[idim].define(nba, dm, 1, 0);
    }
    compute_mac_coefficient(bcoefs, rho, 0, 1.0/rhs_scale);
    mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(bcoefs));

    MLMG mlmg(mlabec);
    if (use_hypre) {
        mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
        mlmg.setBottomVerbose(hypre_verbose);
    }
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setVerbose(verbose);

    Rhs.negate();

    mlmg.setFinalFillBC(true);
    mlmg.solve({mac_phi}, {&Rhs}, mac_tol, mac_abs_tol);
}

//
// EB functions
//
#ifdef AMREX_USE_EB
void eb_mac_level_solve (Amr* parent, const MultiFab* cphi,
				  const BCRec& phys_bc,
			 int level, int Density,
			 Real mac_tol, Real mac_abs_tol, Real rhs_scale,
			 const MultiFab *area, const MultiFab &volume,
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
    // Set up b coef
    // no need to set A coef because it's zero     
    //
    Vector<Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM>> bcoefs;
    bcoefs.resize(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        const BoxArray& nba =
	  amrex::convert(ba,IntVect::TheDimensionVector(idim));
        bcoefs[0][idim].reset(new  MultiFab(nba, dm, 1, 0, MFInfo(), (parent->getLevel(level)).Factory()));
    }
    // set to huge val if covered, do normal geometric average otherwise
    compute_mac_coefficient_eb(bcoefs[0], S, Density, 1.0/rhs_scale);
    // make sure BCs are filled
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
      bcoefs[0][idim]->FillBoundary( geom.periodicity() );
    
    // create right containers for MacProjector
    Vector<Array<MultiFab*,AMREX_SPACEDIM>>  a_umac;
    a_umac.resize(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
      a_umac[0][idim]= &(u_mac[idim]);

    //
    // Create MacProjection object
    // make sure u_mac already has bndry properly filled

    MacProjector macproj( a_umac, GetVecOfArrOfPtrsConst(bcoefs), {geom}, {}, {&Rhs}); 
    //
    // Need to figure out how to access some switches now buried in
    //   amrex::MacProjector
    //   MacProjector owns MLLinOp* m_linop
    //   linop owns LPInfo info
    //   MacProjector owns MLEBABecLap m_eb_abeclap
    
    // LPInfo info;
    // info.setAgglomeration(agglomeration);
    // info.setConsolidation(consolidation);

    // MLABecLaplacian mlabec({geom}, {ba}, {dm}, info);
    //mlabec.setMaxOrder(max_order);

    std::array<MLLinOp::BCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<MLLinOp::BCType,AMREX_SPACEDIM> mlmg_hibc;
    set_mac_solve_bc(mlmg_lobc, mlmg_hibc, phys_bc, geom);

    macproj.setDomainBC(mlmg_lobc, mlmg_hibc);
    // if (level > 0) {
    //     mlabec.setCoarseFineBC(cphi, parent->refRatio(level-1)[0]);
    // }
    // mlabec.setLevelBC(0, mac_phi);
    macproj.setVerbose(verbose);
    
    if (use_hypre) {
      macproj.setBottomSolver(MLMG::BottomSolver::hypre);
      //macproj.getMLMG().setBottomVerbose(hypre_verbose);
      macproj.setVerbose(hypre_verbose);
    }
    // macproj.getMLMG().setMaxFmgIter(max_fmg_iter);
    // macproj.getMLMG().setVerbose(verbose);

    // not 100% sure whether to give initial guess or not
    bool steady_state=false;
    if (steady_state)
    {
      // Solve using m_phi as an initial guess
      macproj.project({mac_phi}, mac_tol,mac_abs_tol);
    } 
    else 
    {
      // Solve with initial guess of zero
      macproj.project(mac_tol,mac_abs_tol);
    }

    // set mac velocity bcs here?
    for (int i=0; i<BL_SPACEDIM; ++i)
      u_mac[i].FillBoundary(geom.periodicity());

    // FIXME Get this working
    if (verbose)
    {
      MultiFab divu(ba, dm, 1, 0, MFInfo(), (parent->getLevel(level)).Factory());
      EB_computeDivergence(divu,GetArrOfConstPtrs(a_umac[0]),geom);
      
      MultiFab tmp(ba, dm, 1, 0, MFInfo(), (parent->getLevel(level)).Factory());
      MultiFab::Copy( tmp, divu, 0, 0, 1, 0 );
      EB_set_covered( tmp, 0.0 );

      Print() << "  * On level "<< level
    	      << " max(abs(divu)) = " << divu.norm0() << "\n";
    } 

}


#endif
