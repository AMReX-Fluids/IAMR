
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_ParmParse.H>

#include <MacOpMacDrivers.H>

#include <IAMR_MLMG_F.H>

using namespace amrex;

namespace {
    static int initialized = false;
    static int max_order = 4;
    static int agglomeration = 1;
    static int consolidation = 1;
    static int max_fmg_iter = 0;
}

namespace {
static void compute_mac_coefficient (std::array<MultiFab,AMREX_SPACEDIM>& bcoefs, 
                                     const MultiFab& rho, int rho_comp)
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
                      BL_TO_FORTRAN_N_ANYD(rho[mfi],rho_comp));
    }
}

static void compute_mac_rhs (MultiFab& rhs, const MultiFab* umac, Real rhs_scale,
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
                     BL_TO_FORTRAN_ANYD(volume[mfi]),
#else
                     dxinv,
#endif
                     &rhs_scale);
    }
}

}

void mlmg_mac_level_solve (Amr* parent, const MultiFab* cphi, const BCRec& phys_bc,
                           int level, int Density, Real dt,
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
        
        initialized = true;
    }
    
    const int old_agg = MLLinOp::setAgglomeration(agglomeration);
    const int old_con = MLLinOp::setConsolidation(consolidation);

    const Geometry& geom = parent->Geom(level);
    const BoxArray& ba = Rhs.boxArray();
    const DistributionMapping& dm = Rhs.DistributionMap();

    MLABecLaplacian mlabec({geom}, {ba}, {dm});
    mlabec.setMaxOrder(max_order);

    std::array<MLLinOp::BCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<MLLinOp::BCType,AMREX_SPACEDIM> mlmg_hibc;
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

    mlabec.setDomainBC(mlmg_lobc, mlmg_hibc);
    if (level > 0) {
        mlabec.setBCWithCoarseData(cphi, parent->refRatio(level-1)[0]);
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
    compute_mac_coefficient(bcoefs, S, Density);
    mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(bcoefs));

    MLMG mlmg(mlabec);
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setVerbose(verbose);

    const Real* dxinv = geom.InvCellSize();
    compute_mac_rhs(Rhs, u_mac, rhs_scale, area, volume, dxinv);
        
    mlmg.solve({mac_phi}, {&Rhs}, mac_tol, mac_abs_tol);

    VisMF::Write(*mac_phi, "phi");

    MLLinOp::setAgglomeration(old_agg);
    MLLinOp::setConsolidation(old_con);

    amrex::Abort("xxxxx");
}

