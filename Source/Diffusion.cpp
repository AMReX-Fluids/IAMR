//fixme, for writesingle level plotfile
#include<AMReX_PlotFileUtil.H>
//

#include <AMReX_ParmParse.H>

#include <Diffusion.H>
#include <NavierStokesBase.H>

#include <algorithm>
#include <cfloat>
#include <iomanip>
#include <array>

#include <iostream>

#include <AMReX_Utility.H>
#include <AMReX_MLMG.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLEBTensorOp.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EB_utils.H>
#else
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLTensorOp.H>
#endif

using namespace amrex;

#if defined(BL_OSF1)
#if defined(BL_USE_DOUBLE)
const Real BL_BOGUS      = DBL_QNAN;
#else
const Real BL_BOGUS      = FLT_QNAN;
#endif
#else
const Real BL_BOGUS      = 1.e200;
#endif

const Real BL_SAFE_BOGUS = -666.e200;

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

namespace
{
    bool initialized = false;
    static int agglomeration = 1;
    static int consolidation = 1;
    static int max_fmg_iter = 0;
    static int use_hypre = 0;
    static int hypre_verbose = 0;
}
//
// Set default values in !initialized section of code in constructor!!!
//
int         Diffusion::verbose;
Real        Diffusion::visc_tol;
int         Diffusion::do_reflux;
int         Diffusion::max_order;
int         Diffusion::scale_abec;
int         Diffusion::tensor_max_order;

Vector<Real> Diffusion::visc_coef;
Vector<int>  Diffusion::is_diffusive;

void
Diffusion::Finalize ()
{
    visc_coef.clear();
    is_diffusive.clear();

    initialized = false;
}

Diffusion::Diffusion (Amr*               Parent,
                      NavierStokesBase*  Caller,
                      Diffusion*         Coarser,
                      int                num_state,
                      FluxRegister*      Viscflux_reg,
                      const Vector<int>&  _is_diffusive,
                      const Vector<Real>& _visc_coef)
    :
    parent(Parent),
    navier_stokes(Caller),
    grids(navier_stokes->boxArray()),
    dmap(navier_stokes->DistributionMap()),
    level(navier_stokes->Level()),
    coarser(Coarser),
    finer(0),
    NUM_STATE(num_state),
    viscflux_reg(Viscflux_reg)

{
    if (!initialized)
    {
        //
        // Set defaults here!!!
        //
        Diffusion::verbose             = 0;
        Diffusion::visc_tol            = 1.0e-10;
        Diffusion::do_reflux           = 1;
        Diffusion::scale_abec          = 0;
        //
        // It is essential that we set max_order of the solver to 2
        // if we want to use the standard sol(i)-sol(i-1) approximation
        // for the gradient at Dirichlet boundaries.
        // The solver's default order is 3 and this uses three points for the
        // gradient at a Dirichlet boundary.
        //
        Diffusion::max_order           = 2;
        Diffusion::tensor_max_order    = 2;

        ParmParse ppdiff("diffuse");

        ppdiff.query("v",                   verbose);
        ppdiff.query("scale_abec",          scale_abec);
        ppdiff.query("max_order",           max_order);
        ppdiff.query("tensor_max_order",    tensor_max_order);

        ppdiff.query("agglomeration", agglomeration);
        ppdiff.query("consolidation", consolidation);
        ppdiff.query("max_fmg_iter", max_fmg_iter);
#ifdef AMREX_USE_HYPRE
        ppdiff.query("use_hypre", use_hypre);
        ppdiff.query("hypre_verbose", hypre_verbose);
#endif

        ParmParse pp("ns");

        pp.query("visc_tol",  visc_tol);
        pp.query("do_reflux", do_reflux);

        do_reflux = (do_reflux ? 1 : 0);

        const int n_visc = _visc_coef.size();
        const int n_diff = _is_diffusive.size();

        if (n_diff < NUM_STATE)
            amrex::Abort("Diffusion::Diffusion(): is_diffusive array is not long enough");

        if (n_visc < NUM_STATE)
            amrex::Abort("Diffusion::Diffusion(): visc_coef array is not long enough");

        if (n_visc > NUM_STATE)
            amrex::Abort("Diffusion::Diffusion(): TOO MANY diffusion coeffs were given!");

        visc_coef.resize(NUM_STATE);
        is_diffusive.resize(NUM_STATE);

        for (int i = 0; i < NUM_STATE; i++)
        {
            is_diffusive[i] = _is_diffusive[i];
            visc_coef[i] = _visc_coef[i];
        }

        echo_settings();

        amrex::ExecOnFinalize(Diffusion::Finalize);

        initialized = true;
    }

    if (level > 0)
    {
        crse_ratio = parent->refRatio(level-1);
        coarser->finer = this;
    }
}

Diffusion::~Diffusion () {}

FluxRegister*
Diffusion::viscFluxReg ()
{
    return viscflux_reg;
}

int
Diffusion::maxOrder() const
{
    return max_order;
}

int
Diffusion::tensorMaxOrder() const
{
    return tensor_max_order;
}

void
Diffusion::echo_settings () const
{
    //
    // Print out my settings.
    //
  if (verbose)
    {
        amrex::Print() << "Diffusion settings...\n";
        amrex::Print() << "  From diffuse:\n";
        amrex::Print() << "   max_order           = " << max_order           << '\n';
        amrex::Print() << "   tensor_max_order    = " << tensor_max_order    << '\n';
        amrex::Print() << "   scale_abec          = " << scale_abec          << '\n';

        amrex::Print() << "\n\n  From ns:\n";
        amrex::Print() << "   do_reflux           = " << do_reflux << '\n';
        amrex::Print() << "   visc_tol            = " << visc_tol  << '\n';

        amrex::Print() << "   is_diffusive =";
        for (int i =0; i < NUM_STATE; i++)
            amrex::Print() << "  " << is_diffusive[i];

        amrex::Print() << "\n   visc_coef =";
        for (int i = 0; i < NUM_STATE; i++)
            amrex::Print() << "  " << visc_coef[i];

        amrex::Print() << '\n';
    }
}

Real
Diffusion::get_scaled_abs_tol (const MultiFab& rhs,
                               Real            reduction) //const
{
    Real oncomp(1.0/rhs.nComp());
    Real rhs_avg_norm(0.0);
    // Let's take an average on all the components
    // This is the only way I can think of to prevent absolute tolerance to be
    // execessively small
    for (int comp(0); comp < rhs.nComp(); ++comp)
        rhs_avg_norm += oncomp * rhs.norm0(comp,0,false,true);
    return reduction * rhs_avg_norm;
}

void
Diffusion::diffuse_scalar (const Vector<MultiFab*>&  S_old,
                           const Vector<MultiFab*>&  Rho_old,
                           Vector<MultiFab*>&        S_new,
                           const Vector<MultiFab*>&  Rho_new,
                           int                       S_comp,
                           int                       num_comp,
                           int                       Rho_comp,
                           Real                      prev_time,
                           Real                      curr_time,
                           Real                      be_cn_theta,
                           const MultiFab&           rho_half,
                           int                       rho_flag,
                           MultiFab* const*          fluxn,
                           MultiFab* const*          fluxnp1,
                           int                       fluxComp,
                           MultiFab*                 delta_rhs,
                           int                       rhsComp,
                           const MultiFab*           alpha_in,
                           int                       alpha_in_comp,
                           const MultiFab* const*    betan,
                           const MultiFab* const*    betanp1,
                           int                       betaComp,
                           const IntVect&            cratio,
                           const BCRec&              bc,
                           const Geometry&           geom,
                           const SolveMode&          solve_mode,
                           bool                      add_old_time_divFlux,
                           const amrex::Vector<int>& is_diffusive)
{
    //
    // This routine expects that physical BC's have been loaded into
    // the grow cells of the old and new state at this level.  If rho_flag==2,
    // the values there are rho.phi, where phi is the quantity being diffused.
    // Values in these cells will be preserved.  Also, if there are any
    // explicit update terms, these have already incremented the new state
    // on the valid region (i.e., on the valid region the new state is the old
    // state + dt*Div(explicit_fluxes), e.g.)
    //

    if (verbose)
      amrex::Print() << "... Diffusion::diffuse_scalar(): \n"
                     << " lev: " << level << '\n';

    // Velocity components should go to tensor solver
    if (S_comp <= Xvel && Xvel <= S_comp+num_comp-1){
       amrex::Abort("Diffusion::diffuse_scalar(): velocity component(s) attemping to use scalar solver. Velocity must use tensor solver.\n");
    }

    bool has_coarse_data = S_new.size() > 1;

    const Real strt_time = ParallelDescriptor::second();

    int allthere, allnull;
    checkBeta(betan, allthere, allnull);
    checkBeta(betanp1, allthere);
    if (allnull && add_old_time_divFlux && be_cn_theta!=1)
      amrex::Abort("Diffusion::diffuse_scalar: Constant diffusivity case no longer supported separately. Must set non-zero beta.");

    //
    // No ghost cells are needed for MLMG in most cases. Except for
    // for cell-centered solver, you need to call setLevelBC, and that
    // needs to have 1 ghost cell if there is Dirichlet BC.
    //
    const int ng = 1;

    Real dt = curr_time - prev_time;
    BL_ASSERT(S_new[0]->nGrow()>0);
    if (S_old.size()>0) BL_ASSERT(S_old[0]->nGrow()>0);
    const BoxArray& ba = S_new[0]->boxArray();
    const DistributionMapping& dm = S_new[0]->DistributionMap();
    const DistributionMapping* dmc = (has_coarse_data ? &(S_new[1]->DistributionMap()) : 0);
    const BoxArray* bac = (has_coarse_data ? &(S_new[1]->boxArray()) : 0);

    BL_ASSERT(solve_mode==ONEPASS || (delta_rhs && delta_rhs->boxArray()==ba));

    const auto& factory = S_new[0]->Factory();

    MultiFab Rhs(ba,dm,1,0,MFInfo(),factory);
    MultiFab Soln(ba,dm,1,ng,MFInfo(),factory);

    auto Solnc = std::unique_ptr<MultiFab>(new MultiFab());
    if (has_coarse_data)
    {
      Solnc->define(*bac, *dmc, 1, 0, MFInfo(), S_new[1]->Factory());
    }

    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;

    //
    // Create operator at time n and n+1
    //
    LPInfo infon;
    infon.setAgglomeration(agglomeration);
    infon.setConsolidation(consolidation);
    infon.setMetricTerm(false);
    infon.setMaxCoarseningLevel(0);

#ifdef AMREX_USE_EB
    const auto& ebf = &(dynamic_cast<EBFArrayBoxFactory const&>(factory));
    MLEBABecLap opn({geom}, {ba}, {dm}, infon, {ebf});
#else
    MLABecLaplacian opn({geom}, {ba}, {dm}, infon);
#endif

    opn.setMaxOrder(max_order);
    MLMG mgn(opn);
    mgn.setVerbose(verbose);

    LPInfo infonp1;
    infonp1.setAgglomeration(agglomeration);
    infonp1.setConsolidation(consolidation);
    infonp1.setMetricTerm(false);

#ifdef AMREX_USE_EB
    MLEBABecLap opnp1({geom}, {ba}, {dm}, infonp1, {ebf});
#else
    MLABecLaplacian opnp1({geom}, {ba}, {dm}, infonp1);
#endif

    opnp1.setMaxOrder(max_order);

    MLMG mgnp1(opnp1);
    if (use_hypre)
    {
      mgnp1.setBottomSolver(MLMG::BottomSolver::hypre);
      mgnp1.setBottomVerbose(hypre_verbose);
    }
    mgnp1.setMaxFmgIter(max_fmg_iter);
    mgnp1.setVerbose(verbose);

    setDomainBC(mlmg_lobc, mlmg_hibc, bc); // Same for all comps, by assumption
    opn.setDomainBC(mlmg_lobc, mlmg_hibc);
    opnp1.setDomainBC(mlmg_lobc, mlmg_hibc);

    for (int icomp=0; icomp<num_comp; ++icomp)
    {
       if (verbose)
       {
          amrex::Print() << "diffusing scalar "<< icomp+1 << " of " << num_comp << "\n"
                         << "rho flag "<< rho_flag << "\n";
       }

       int sigma = S_comp + icomp;

       if (is_diffusive[icomp] == 0)
       {
          for (int n = 0; n < AMREX_SPACEDIM; n++)
          {
            if (fluxn[n]!=0 && fluxnp1[n]!=0)
            {
              fluxn[n]->setVal(0.0,fluxComp+icomp,1);
              fluxnp1[n]->setVal(0.0,fluxComp+icomp,1);
            }
          }
          break;
       }

       if (add_old_time_divFlux && be_cn_theta!=1)
       {
          Real a = 0.0;
          Real b = -(1.0-be_cn_theta)*dt;

          if(verbose)
             Print()<<"Adding old time diff ...\n";

          {
            if (has_coarse_data)
            {
              MultiFab::Copy(*Solnc,*S_old[1],sigma,0,1,0);
              if (rho_flag == 2)
              {
                MultiFab::Divide(*Solnc,*Rho_old[1],Rho_comp,0,1,0);
              }
              opn.setCoarseFineBC(Solnc.get(), cratio[0]);
            }
            MultiFab::Copy(Soln,*S_old[0],sigma,0,1,ng);
            if (rho_flag == 2)
            {
              MultiFab::Divide(Soln,*Rho_old[0],Rho_comp,0,1,ng);
            }
            opn.setLevelBC(0, &Soln);
          }

          opn.setScalars(a,b);
          setBeta(opn,betan,betaComp+icomp);

#ifdef AMREX_USE_EB
          MultiFab rhs_tmp(ba,dm,1,2,MFInfo(),factory);
          rhs_tmp.setVal(0.);
          mgn.apply({&rhs_tmp},{&Soln});

          const amrex::MultiFab* weights;
          const auto& ebf = &(dynamic_cast<EBFArrayBoxFactory const&>(factory));
          weights = &(ebf->getVolFrac());

          amrex::single_level_weighted_redistribute(rhs_tmp, Rhs, *weights, 0, 1, geom);
#else
          mgn.apply({&Rhs},{&Soln});
#endif

          computeExtensiveFluxes(mgn, Soln, fluxn, fluxComp+icomp, 1, &geom, -b/dt);
       }
       else
       {
          for (int n = 0; n < AMREX_SPACEDIM; n++)
          {
            fluxn[n]->setVal(0.0,fluxComp+icomp,1);
          }
          Rhs.setVal(0.0);
       }

       // Get local version of the switch
       int has_alpha     = (alpha_in  != 0) ? 1 : 0;
       int has_delta_rhs = (delta_rhs != 0) ? 1 : 0;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter Smfi(*S_new[0], TilingIfNotGPU()); Smfi.isValid(); ++Smfi)
       {
          const Box& bx = Smfi.tilebox();
          auto const& rhs      = Rhs.array(Smfi);
          auto const& solution = Soln.array(Smfi);
          auto const& snew     = S_new[0]->array(Smfi,sigma);
	  Array4<Real> dummy;
          auto const& sold     = (solve_mode == PREDICTOR) ? S_old[0]->array(Smfi,sigma) : dummy;
          auto const& rhoHalf  = (rho_flag == 1) ? rho_half.array(Smfi) : dummy;
          auto const& rho_old  = (rho_flag == 3) ? Rho_old[0]->array(Smfi,Rho_comp) : dummy;
          auto const& alpha    = (has_alpha) ? alpha_in->array(Smfi,alpha_in_comp+icomp) : Soln.array(Smfi);
          auto const& deltarhs = (has_delta_rhs) ? delta_rhs->array(Smfi,rhsComp+icomp) : Soln.array(Smfi);
          Real dtinv = 1.0/dt;

          amrex::ParallelFor(bx, [rhs, solution, snew, sold, rhoHalf, rho_old, alpha, deltarhs,
                                  has_alpha, has_delta_rhs, solve_mode, rho_flag, dtinv, dt ]
          AMREX_GPU_DEVICE(int i, int j, int k) noexcept
          {
             // If this is a predictor step, put "explicit" updates into Rhs after scaling by rho_half if reqd
             if (solve_mode == PREDICTOR) {
                amrex::Real tmp = snew(i,j,k) - sold(i,j,k);
                snew(i,j,k) -= tmp;
                tmp *= dtinv;
                if ( rho_flag == 1 ) {
                   tmp *= rhoHalf(i,j,k);
                }
                if ( has_alpha ) {
                   tmp *= alpha(i,j,k);
                }
                rhs(i,j,k) += tmp;
             }

             // Add body sources
             if ( has_delta_rhs ) {
                rhs(i,j,k) += deltarhs(i,j,k) * dt;
                if ( rho_flag == 1 ) {
                   rhs(i,j,k) *= rhoHalf(i,j,k);
                }
             }
             // Always do : copy S_new into Soln
             //             add Soln to Rhs
             solution(i,j,k) = snew(i,j,k);
             if ( rho_flag == 1) {
                solution(i,j,k) *= rhoHalf(i,j,k);
             } else if ( rho_flag == 3 ) {
                solution(i,j,k) *= rho_old(i,j,k);
             }
             if ( has_alpha ) {
                solution(i,j,k) *= alpha(i,j,k);
             }
             rhs(i,j,k) += solution(i,j,k);
          });
       }


       //
       // Construct viscous operator with bndry data at time N+1.
       //
       Real a = 1.0;
       Real b = be_cn_theta*dt;
       Real rhsscale = 1.0;
       {
          if (has_coarse_data)
          {
            MultiFab::Copy(*Solnc,*S_new[1],sigma,0,1,0);
            if (rho_flag == 2)
            {
              MultiFab::Divide(*Solnc,*Rho_new[1],Rho_comp,0,1,0);
            }
            opnp1.setCoarseFineBC(Solnc.get(), cratio[0]);
          }
          MultiFab::Copy(Soln,*S_new[0],sigma,0,1,ng);
          if (rho_flag == 2)
          {
            MultiFab::Divide(Soln,*Rho_new[0],Rho_comp,0,1,ng);
          }
          opnp1.setLevelBC(0, &Soln);
       }

       {
          std::pair<Real,Real> scalars;
          MultiFab alpha;
          const MultiFab* rho = (rho_flag == 1) ? &rho_half : Rho_new[0];
          int rhoComp = (rho_flag == 1 ) ? 0 : Rho_comp;

          computeAlpha(alpha, scalars, a, b,
                       &rhsscale, alpha_in, alpha_in_comp+icomp,
                       rho_flag, rho, rhoComp);
          opnp1.setScalars(scalars.first, scalars.second);
          opnp1.setACoeffs(0, alpha);
       }

       setBeta(opnp1,betanp1,betaComp+icomp);

       Rhs.mult(rhsscale,0,1);
       const Real S_tol     = visc_tol;
       const Real S_tol_abs = get_scaled_abs_tol(Rhs, visc_tol);

       mgnp1.solve({&Soln}, {&Rhs}, S_tol, S_tol_abs);

       computeExtensiveFluxes(mgnp1, Soln, fluxnp1, fluxComp+icomp, 1, &geom, b/dt);

       //
       // Copy into state variable at new time, without bc's and * rho if needed
       //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter Smfi(*S_new[0], TilingIfNotGPU()); Smfi.isValid(); ++Smfi)
       {
          const Box& bx = Smfi.tilebox();
          const auto& snew     = S_new[0]->array(Smfi, sigma);
          const auto& solution = Soln.array(Smfi);
          const auto& rhonew   = (rho_flag == 2) ? Rho_new[0]->array(Smfi,Rho_comp) : S_new[0]->array(Smfi, sigma);
          amrex::ParallelFor(bx, [snew, solution, rhonew, rho_flag]
          AMREX_GPU_DEVICE(int i, int j, int k) noexcept
          {
             snew(i,j,k) = solution(i,j,k);
             if (rho_flag == 2) {
                snew(i,j,k) *= rhonew(i,j,k);
             }
          });
       }
       if (verbose) amrex::Print() << "Done with diffuse_scalar" << "\n";
    }

    if (verbose) {
       const int IOProc   = ParallelDescriptor::IOProcessorNumber();
       Real      run_time = ParallelDescriptor::second() - strt_time;
       ParallelDescriptor::ReduceRealMax(run_time,IOProc);
       amrex::Print() << "Diffusion::diffuse_scalar() time: " << run_time << '\n';
    }
}

void
Diffusion::diffuse_velocity (Real                   dt,
                             Real                   be_cn_theta,
                             const MultiFab&        rho_half,
                             int                    rho_flag,
                             MultiFab*              delta_rhs,
                             const MultiFab* const* betan,
                             const MultiFab* const  betanCC,
                             const MultiFab* const* betanp1,
                             const MultiFab* const  betanp1CC)
{
    diffuse_velocity(dt, be_cn_theta, rho_half, rho_flag,
                     delta_rhs, 0, betan, betanCC, betanp1, betanp1CC, 0);
}

void
Diffusion::diffuse_velocity (Real                   dt,
                             Real                   be_cn_theta,
                             const MultiFab&        rho_half,
                             int                    rho_flag,
                             MultiFab*              delta_rhs,
                             int                    rhsComp,
                             const MultiFab* const* betan,
                             const MultiFab* const  betanCC,
                             const MultiFab* const* betanp1,
                             const MultiFab* const  betanp1CC,
                             int                    betaComp)
{
  if (verbose) amrex::Print() << "... Diffusion::diffuse_velocity() lev: " << level << std::endl;

    const Real strt_time = ParallelDescriptor::second();

    diffuse_tensor_velocity(dt,be_cn_theta,rho_half,rho_flag,
                            delta_rhs,rhsComp,betan,betanCC,betanp1,betanp1CC,betaComp);

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        amrex::Print() << "Diffusion::diffuse_velocity(): lev: " << level
                       << ", time: " << run_time << '\n';
    }
}

void
Diffusion::diffuse_tensor_velocity (Real                   dt,
                                    Real                   be_cn_theta,
                                    const MultiFab&        rho_half,
                                    int                    rho_flag,
                                    MultiFab*              delta_rhs,
                                    int                    rhsComp,
                                    const MultiFab* const* betan,
                                    const MultiFab* const  betanCC,
                                    const MultiFab* const* betanp1,
                                    const MultiFab* const  betanp1CC,
                                    int                    betaComp)
{
    int allthere, allnull;
    checkBeta(betan, allthere, allnull);
    checkBeta(betanp1, allthere);
    if (allnull && be_cn_theta!=1)
      amrex::Abort("Diffusion::diffuse_tensor_velocity: Constant viscosity case no longer supported separately. Must set non-zero beta.");

    BL_ASSERT( rho_flag == 1 || rho_flag == 3);

#ifdef AMREX_DEBUG
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
        BL_ASSERT( betan[d]->min(0,0) >= 0.0 );
#endif

    const int finest_level = parent->finestLevel();
    //
    // At this point, S_old has bndry at time N S_new contains GRAD(SU).
    //
    MultiFab&  U_new     = navier_stokes->get_new_data(State_Type);
    const Real cur_time  = navier_stokes->get_state_data(State_Type).curTime();
    const Real prev_time = navier_stokes->get_state_data(State_Type).prevTime();

    //
    // U_new now contains the inviscid update of U.
    // This is part of the RHS for the viscous solve.
    //
    const int soln_ng = 1;
    int flux_ng = 0;
    MultiFab Rhs(grids,dmap,AMREX_SPACEDIM,0, MFInfo(),navier_stokes->Factory());
    MultiFab Soln(grids,dmap,AMREX_SPACEDIM,soln_ng,MFInfo(),navier_stokes->Factory());
    MultiFab** tensorflux_old;
    FluxBoxes fb_old;

    const Geometry& geom   = navier_stokes->Geom();

    //
    // Set up Rhs.
    //
    if ( be_cn_theta != 1 )
    {
      //
      // Compute time n viscous terms
      //
      const Real a = 0.0;
      Real       b = -(1.0-be_cn_theta)*dt;

      {
         LPInfo info;
         info.setAgglomeration(agglomeration);
         info.setConsolidation(consolidation);
         info.setMaxCoarseningLevel(0);
         info.setMetricTerm(false);

#ifdef AMREX_USE_EB
         const auto& ebf = &dynamic_cast<EBFArrayBoxFactory const&>(navier_stokes->Factory());
         MLEBTensorOp tensorop({navier_stokes->Geom()}, {grids}, {dmap}, info, {ebf});
#else
         MLTensorOp tensorop({navier_stokes->Geom()}, {grids}, {dmap}, info);
#endif

         tensorop.setMaxOrder(tensor_max_order);

         // create right container
         Array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc[AMREX_SPACEDIM];
         Array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc[AMREX_SPACEDIM];
         // fill it
         for (int i=0; i<AMREX_SPACEDIM; i++)
           setDomainBC(mlmg_lobc[i], mlmg_hibc[i], Xvel+i);
         // pass to op
         tensorop.setDomainBC({AMREX_D_DECL(mlmg_lobc[0],mlmg_lobc[1],mlmg_lobc[2])},
                              {AMREX_D_DECL(mlmg_hibc[0],mlmg_hibc[1],mlmg_hibc[2])});

         // set coarse-fine BCs
         {
           MultiFab crsedata;
           int ng = 0;
           if (level > 0) {
             auto& crse_ns = *(coarser->navier_stokes);
             crsedata.define(crse_ns.boxArray(), crse_ns.DistributionMap(),
                             AMREX_SPACEDIM, ng, MFInfo(), crse_ns.Factory());
             AmrLevel::FillPatch(crse_ns, crsedata, ng, prev_time, State_Type, Xvel,
                                 AMREX_SPACEDIM);
             tensorop.setCoarseFineBC(&crsedata, crse_ratio[0]);
           }
           AmrLevel::FillPatch(*navier_stokes,Soln,soln_ng,prev_time,State_Type,Xvel,AMREX_SPACEDIM);
           tensorop.setLevelBC(0, &Soln);
         }

         tensorop.setScalars(a, b);

#ifdef AMREX_USE_EB
         setViscosity(tensorop, betan, betaComp, *betanCC);
#else
         setViscosity(tensorop, betan, betaComp);
#endif

         MLMG mlmg(tensorop);
         // FIXME -- consider making new parameters max_iter and bottom_verbose
         //mlmg.setMaxIter(max_iter);
         mlmg.setMaxFmgIter(max_fmg_iter);
         mlmg.setVerbose(10);
         mlmg.setBottomVerbose(10);
         //mlmg.setBottomVerbose(bottom_verbose);

         int nghost(2);
         MultiFab Rhs_tmp(grids,dmap,AMREX_SPACEDIM,nghost, MFInfo(),navier_stokes->Factory());
         Rhs_tmp.setVal(0.);
         mlmg.apply({&Rhs_tmp}, {&Soln});

#ifdef AMREX_USE_EB
         //
         // Redistribution may go here, at least as long as we're enforcing that
         // the coarse-fine boundary cannot intersect the EB. This garauntees
         // that none of the redistributed fluxes are used/needed since :
         //   reflux only uses the fluxes at coarse-fine boundaries.
         //   redistribution only alters cut cells and their nearest-neighbors.
         //   regridding algorithm buffers the cells flagged for refinement
         //
         const amrex::MultiFab* weights;
         const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(navier_stokes->Factory());
         weights = &(ebfactory.getVolFrac());
         amrex::single_level_weighted_redistribute(Rhs_tmp, Rhs, *weights, 0, AMREX_SPACEDIM, navier_stokes->Geom());
#else
         amrex::Copy(Rhs, Rhs_tmp, 0, 0, AMREX_SPACEDIM, 0);
#endif

         if (do_reflux && (level<finest_level || level>0))
         {
           tensorflux_old = fb_old.define(navier_stokes, AMREX_SPACEDIM);

           computeExtensiveFluxes(mlmg, Soln, tensorflux_old, 0, AMREX_SPACEDIM, &geom, -b/dt);
         }
      }
    }
    else
    {
       Rhs.setVal(0.);
    }

    //
    // Complete Rhs by adding body sources.
    //
    int has_delta_rhs = (delta_rhs != 0) ? 1 : 0;

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(Rhs,mfi_info); mfi.isValid(); ++mfi)
    {
       const Box& bx  = mfi.tilebox();
       auto const& rhs      = Rhs.array(mfi);
       auto const& unew     = U_new.array(mfi,Xvel);
       auto const& rho      = (rho_flag == 1) ? rho_half.array(mfi) : navier_stokes->rho_ptime.array(mfi);
       auto const& deltarhs = (has_delta_rhs) ? delta_rhs->array(mfi,rhsComp) : U_new.array(mfi);
       amrex::ParallelFor(bx, [rhs, unew, rho, deltarhs, has_delta_rhs, dt]
       AMREX_GPU_DEVICE(int i, int j, int k) noexcept
       {
          for (int n = 0; n < AMREX_SPACEDIM; n++) {
             unew(i,j,k,n) *= rho(i,j,k);
             rhs(i,j,k,n) += unew(i,j,k,n);
             if ( has_delta_rhs ) {
                rhs(i,j,k,n) += deltarhs(i,j,k,n) * dt;
             }
          }
       });
    }

    //
    // Construct viscous operator at time N+1.
    //
    const Real a = 1.0;
    Real       b = be_cn_theta*dt;

    // MLMG solution
    {
      // genaric tol suggestion for MLMG
      // const Real tol_rel = 1.e-11;
      // const Real tol_abs = 0.0;
      // cribbing from scalar
      const Real tol_rel = visc_tol;
      const Real tol_abs = get_scaled_abs_tol(Rhs, visc_tol);

      LPInfo info;
      info.setAgglomeration(agglomeration);
      info.setConsolidation(consolidation);
      info.setMetricTerm(false);
      info.setMaxCoarseningLevel(100);

#ifdef AMREX_USE_EB
      const auto& ebf = &dynamic_cast<EBFArrayBoxFactory const&>(navier_stokes->Factory());
      MLEBTensorOp tensorop({navier_stokes->Geom()}, {grids}, {dmap}, info, {ebf});
#else
      MLTensorOp tensorop({navier_stokes->Geom()}, {grids}, {dmap}, info);
#endif

      tensorop.setMaxOrder(tensor_max_order);

      // create right container
      Array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc[AMREX_SPACEDIM];
      Array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc[AMREX_SPACEDIM];
      // fill it
      for (int i=0; i<AMREX_SPACEDIM; i++)
         setDomainBC(mlmg_lobc[i], mlmg_hibc[i], Xvel+i);
      // pass to op
      tensorop.setDomainBC({AMREX_D_DECL(mlmg_lobc[0],mlmg_lobc[1],mlmg_lobc[2])},
                           {AMREX_D_DECL(mlmg_hibc[0],mlmg_hibc[1],mlmg_hibc[2])});

      // set up level BCs
      {
         MultiFab crsedata;
         int ng = 0;
         if (level > 0) {
            auto& crse_ns = *(coarser->navier_stokes);
            crsedata.define(crse_ns.boxArray(), crse_ns.DistributionMap(), AMREX_SPACEDIM,
                            ng, MFInfo(),crse_ns.Factory());
            AmrLevel::FillPatch(crse_ns, crsedata, ng, cur_time, State_Type, Xvel,
                                AMREX_SPACEDIM);
            tensorop.setCoarseFineBC(&crsedata, crse_ratio[0]);
         }
         AmrLevel::FillPatch(*navier_stokes,Soln,soln_ng,cur_time,State_Type,Xvel,AMREX_SPACEDIM);
         tensorop.setLevelBC(0, &Soln);
      }

      {
         MultiFab acoef;
         std::pair<Real,Real> scalars;
         Real rhsscale = 1.0;
         const MultiFab& rho = (rho_flag == 1) ? rho_half : navier_stokes->rho_ctime;
         computeAlpha(acoef, scalars, a, b,
                      &rhsscale, nullptr, 0,
                      rho_flag, &rho, 0);
         tensorop.setScalars(scalars.first, scalars.second);
         tensorop.setACoeffs(0, acoef);
      }

#ifdef AMREX_USE_EB
      setViscosity(tensorop, betanp1, betaComp, *betanp1CC);
#else
      setViscosity(tensorop, betanp1, betaComp);
#endif

      MLMG mlmg(tensorop);
      //fixme?
      //mlmg.setMaxIter(max_iter);
      mlmg.setMaxFmgIter(max_fmg_iter);
      mlmg.setVerbose(verbose);
      //mlmg.setBottomVerbose(bottom_verbose);

      // ensures ghost cells of sol are correctly filled when returned from solver
      mlmg.setFinalFillBC(true);

      //    solution.setVal(0.0);
      mlmg.solve({&Soln}, {&Rhs}, tol_rel, tol_abs);

      //
      // Copy into state variable at new time.
      //
      MultiFab::Copy(U_new,Soln,0,Xvel,AMREX_SPACEDIM,soln_ng);

      //
      // Modify diffusive fluxes here.
      //
      if (do_reflux && (level < finest_level || level > 0))
      {
         FluxBoxes fb(navier_stokes, AMREX_SPACEDIM);
         MultiFab** tensorflux = fb.get();

         computeExtensiveFluxes(mlmg, Soln, tensorflux, 0, AMREX_SPACEDIM, &geom, b/dt);
         if ( be_cn_theta!=1 ) {
            for ( int i = 0; i < AMREX_SPACEDIM; i++)
               MultiFab::Add(*tensorflux[i], *tensorflux_old[i], 0, 0,AMREX_SPACEDIM, flux_ng);
         }

         if (level > 0)
         {
            for (int k = 0; k < AMREX_SPACEDIM; k++)
               viscflux_reg->FineAdd(*(tensorflux[k]),k,Xvel,Xvel,AMREX_SPACEDIM,dt);
         }
         if (level < finest_level)
         {
            for (int d = 0; d < AMREX_SPACEDIM; d++)
               finer->viscflux_reg->CrseInit(*tensorflux[d],d,0,Xvel,AMREX_SPACEDIM,-dt);
         }
      }
   }
}

void
Diffusion::diffuse_Vsync (MultiFab&              Vsync,
                          Real                   dt,
                          Real                   be_cn_theta,
                          const MultiFab&        rho_half,
                          int                    rho_flag,
                          const MultiFab* const* beta,
                          int                    betaComp,
                          bool                   update_fluxreg)
{
    BL_ASSERT(rho_flag == 1|| rho_flag == 3);

    int allthere;
    checkBeta(beta, allthere);

#ifdef AMREX_DEBUG
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
        BL_ASSERT(beta[d]->min(0,0) >= 0.0);
#endif

    diffuse_tensor_Vsync(Vsync,dt,be_cn_theta,rho_half,rho_flag,beta,betaComp,update_fluxreg);
    //
    // applyBC has put "incorrect" values in the ghost cells
    // outside external Dirichlet boundaries. Reset these to zero
    // so that syncproject and conservative interpolation works correctly.
    //
    Box domain = amrex::grow(navier_stokes->Geom().Domain(),1);

    for (int n = Xvel; n < Xvel+AMREX_SPACEDIM; n++)
    {
        const BCRec& velbc = navier_stokes->get_desc_lst()[State_Type].getBC(n);

        for (int k = 0; k < AMREX_SPACEDIM; k++)
        {
            if (velbc.hi(k) == EXT_DIR)
            {
                IntVect smallend = domain.smallEnd();
                smallend.setVal(k,domain.bigEnd(k));
                Box top_strip(smallend,domain.bigEnd(),IntVect::TheCellVector());
                Vsync.setVal(0,top_strip,n-Xvel,1,1);
            }
            if (velbc.lo(k) == EXT_DIR)
            {
                IntVect bigend = domain.bigEnd();
                bigend.setVal(k,domain.smallEnd(k));
                Box bottom_strip(domain.smallEnd(),bigend,IntVect::TheCellVector());
                Vsync.setVal(0,bottom_strip,n-Xvel,1,1);
            }
        }
    }
}

void
Diffusion::diffuse_tensor_Vsync (MultiFab&              Vsync,
                                 Real                   dt,
                                 Real                   be_cn_theta,
                                 const MultiFab&        rho_half,
                                 int                    rho_flag,
                                 const MultiFab* const* beta,
                                 int                    betaComp,
                                 bool                   update_fluxreg)
{
    BL_ASSERT(rho_flag == 1 || rho_flag == 3);

    if (verbose) amrex::Print() << "Diffusion::diffuse_tensor_Vsync ...\n";

    MultiFab Rhs(grids,dmap,AMREX_SPACEDIM,0,MFInfo(),navier_stokes->Factory());

    MultiFab::Copy(Rhs,Vsync,0,0,AMREX_SPACEDIM,0);

    if (verbose > 1)
    {
        Real r_norm = Rhs.norm0();
        amrex::Print() << "Original max of Vsync " << r_norm << '\n';
    }

    //
    // Multiply RHS by density.
    //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(Rhs, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.tilebox();
       auto const& rhs = Rhs.array(mfi);
       auto const& rho = (rho_flag == 1) ? rho_half.array(mfi) : navier_stokes->rho_ptime.array(mfi);

       amrex::ParallelFor(bx, [rhs, rho]
       AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
          for (int n = 0; n < AMREX_SPACEDIM; n++) {
             rhs(i,j,k,n) *= rho(i,j,k);
          }
       });
    }

    //
    // SET UP COEFFICIENTS FOR VISCOUS SOLVER.
    //
    const Real      a         = 1.0;
    const Real      b         = be_cn_theta*dt;
    Real rhsscale = 1.0;

    int soln_ng = 1;
    MultiFab Soln(grids,dmap,AMREX_SPACEDIM,soln_ng, MFInfo(),navier_stokes->Factory());
    Soln.setVal(0.0);

    // MLMG
    const Real tol_rel = visc_tol;
    const Real tol_abs = -1;

    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setMetricTerm(false);
    //info.setMaxCoarseningLevel(100);

#ifdef AMREX_USE_EB
    const auto& ebf = &dynamic_cast<EBFArrayBoxFactory const&>(navier_stokes->Factory());
    MLEBTensorOp tensorop({navier_stokes->Geom()}, {grids}, {dmap}, info, {ebf});
#else
    MLTensorOp tensorop({navier_stokes->Geom()}, {grids}, {dmap}, info);
#endif

    tensorop.setMaxOrder(tensor_max_order);

    // create right container
    Array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc[AMREX_SPACEDIM];
    Array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc[AMREX_SPACEDIM];
    // fill it
    for (int i=0; i<AMREX_SPACEDIM; i++)
      setDomainBC(mlmg_lobc[i], mlmg_hibc[i], Xvel+i);
    // pass to op
    tensorop.setDomainBC({AMREX_D_DECL(mlmg_lobc[0],mlmg_lobc[1],mlmg_lobc[2])},
                         {AMREX_D_DECL(mlmg_hibc[0],mlmg_hibc[1],mlmg_hibc[2])});

    // set up level BCs
    if (level > 0) {
      tensorop.setCoarseFineBC(nullptr, crse_ratio[0]);
    }
    tensorop.setLevelBC(0, nullptr);

    {
      MultiFab acoef;
      std::pair<Real,Real> scalars;
      const MultiFab& rho = (rho_flag == 1) ? rho_half : navier_stokes->rho_ctime;
      computeAlpha(acoef, scalars, a, b,
                   &rhsscale, nullptr, 0,
                   rho_flag, &rho, 0);
      tensorop.setScalars(scalars.first, scalars.second);
      tensorop.setACoeffs(0, acoef);
    }

    {
      FluxBoxes  fb_bcoef;
      MultiFab** face_bcoef = 0;
      face_bcoef = fb_bcoef.define(navier_stokes);
      for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
         face_bcoef[dir]->setVal(1.0);
      }

#ifdef AMREX_USE_EB
      MultiFab bcoefCC(grids,dmap,1,0,MFInfo(),navier_stokes->Factory());
      bcoefCC.setVal(1.0);
      setViscosity(tensorop, face_bcoef, 0, bcoefCC);
#else
      setViscosity(tensorop, face_bcoef, 0);
#endif
    }

    MLMG mlmg(tensorop);
    //fixme?
    //mlmg.setMaxIter(max_iter);
    //mlmg.setBottomVerbose(bottom_verbose);
    if (use_hypre) {
       mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
       mlmg.setBottomVerbose(hypre_verbose);
    }
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setVerbose(verbose);

    Rhs.mult(rhsscale,0,1);

    mlmg.setFinalFillBC(true);
    mlmg.solve({&Soln}, {&Rhs}, tol_rel, tol_abs);

    //
    // Copy into state variable at new time.
    //
    MultiFab::Copy(Vsync,Soln,0,0,AMREX_SPACEDIM,soln_ng);

    if (verbose > 1)
    {
        Real s_norm = Soln.norm0(0,Soln.nGrow());
        amrex::Print() << "Final max of Vsync " << s_norm << '\n';
    }

    if (level > 0)
    {
        FluxBoxes fb(navier_stokes, AMREX_SPACEDIM);
        MultiFab** tensorflux = fb.get();
        const Geometry& geom   = navier_stokes->Geom();
        //
        // The extra factor of dt comes from the fact that Vsync looks
        // like dV/dt, not just an increment to V.
        //
        computeExtensiveFluxes(mlmg, Soln, tensorflux, 0, AMREX_SPACEDIM, &geom, b/dt);

        if (update_fluxreg)
        {
           for (int k = 0; k < AMREX_SPACEDIM; k++) {
              viscflux_reg->FineAdd(*(tensorflux[k]),k,Xvel,Xvel,AMREX_SPACEDIM,dt*dt);
           }
        }
    }
}

//
// Used by PeleLM to sync species
//
void
Diffusion::diffuse_Ssync (MultiFab&              Ssync,
                          int                    sigma,
                          Real                   dt,
                          Real                   be_cn_theta,
                          const MultiFab&        rho_half,
                          int                    rho_flag,
                          MultiFab* const*       flux,
                          int                    fluxComp,
                          const MultiFab* const* beta,
                          int                    betaComp,
                          const MultiFab*        alpha,
                          int                    alphaComp)
{
    const int state_ind    = sigma + AMREX_SPACEDIM;
    if (verbose)
    {
        amrex::Print() << "Diffusion::diffuse_Ssync lev: " << level << " "
                       << navier_stokes->get_desc_lst()[State_Type].name(state_ind) << '\n';
    }

    const Real strt_time = ParallelDescriptor::second();

    int allthere;
    checkBeta(beta, allthere);

    MultiFab  Rhs(grids,dmap,1,0,MFInfo(),navier_stokes->Factory());

    MultiFab::Copy(Rhs,Ssync,sigma,0,1,0);

    if (verbose > 1)
    {
        MultiFab junk(grids,dmap,1,0,MFInfo(),navier_stokes->Factory());
        MultiFab::Copy(junk,Rhs,0,0,1,0);
        if (rho_flag == 2)
        {
            MultiFab& S_new = navier_stokes->get_new_data(State_Type);
            MultiFab::Divide(junk, S_new, Density, 0, 1, 0);
        }
        Real r_norm = junk.norm0();
        amrex::Print() << "Original max of Ssync " << r_norm << '\n';
    }
    //
    // SET UP COEFFICIENTS FOR VISCOUS SOLVER.
    //
    const Real  a = 1.0;
    Real        b = be_cn_theta*dt;
    Real rhsscale = 1.0;

    const Real S_tol     = visc_tol;
    const Real S_tol_abs = -1;

    MultiFab Soln(grids,dmap,1,1,MFInfo(),navier_stokes->Factory());
    Soln.setVal(0);

    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setMetricTerm(false);

#ifdef AMREX_USE_EB
    const auto& ebf = &dynamic_cast<EBFArrayBoxFactory const&>(navier_stokes->Factory());
    MLEBABecLap mlabec({navier_stokes->Geom()}, {grids}, {dmap}, info, {ebf});
#else
    MLABecLaplacian mlabec({navier_stokes->Geom()}, {grids}, {dmap}, info);
#endif
    mlabec.setMaxOrder(max_order);

    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
    setDomainBC(mlmg_lobc, mlmg_hibc, state_ind);

    mlabec.setDomainBC(mlmg_lobc, mlmg_hibc);
    if (level > 0) {
      mlabec.setCoarseFineBC(nullptr, crse_ratio[0]);
    }
    mlabec.setLevelBC(0, nullptr);

    {
      MultiFab acoef;
      std::pair<Real,Real> scalars;
      const Real cur_time = navier_stokes->get_state_data(State_Type).curTime();
      MultiFab& S = navier_stokes->get_data(State_Type,cur_time);
      const MultiFab& rho = (rho_flag == 1) ? rho_half : S;
      const int Rho_comp = (rho_flag ==1) ? 0 : Density;

      computeAlpha(acoef, scalars, a, b,
                   &rhsscale, alpha, alphaComp,
                   rho_flag, &rho, Rho_comp);
      mlabec.setScalars(scalars.first, scalars.second);
      mlabec.setACoeffs(0, acoef);
    }

    setBeta(mlabec, beta, betaComp);

    MLMG mlmg(mlabec);
    if (use_hypre) {
      mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
      mlmg.setBottomVerbose(hypre_verbose);
    }
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setVerbose(verbose);

    if (rho_flag == 1) {
       MultiFab::Multiply(Rhs,rho_half,0,0,1,0);
    }
    Rhs.mult(rhsscale,0,1,0);

    mlmg.solve({&Soln}, {&Rhs}, S_tol, S_tol_abs);

    int flux_allthere, flux_allnull;
    checkBeta(flux, flux_allthere, flux_allnull);
    if (flux_allthere)
    {
      computeExtensiveFluxes(mlmg, Soln, flux, fluxComp, 1,
                             &(navier_stokes->Geom()), b/dt);
    }

    MultiFab::Copy(Ssync,Soln,0,sigma,1,0);

    if (verbose > 1)
    {
        Real s_norm = Soln.norm0(0,Soln.nGrow());
        amrex::Print() << "Final max of Ssync " << s_norm << '\n';
    }

    if (rho_flag == 2)
    {
        MultiFab& S_new = navier_stokes->get_new_data(State_Type);
        MultiFab::Multiply(Ssync,S_new,Density,sigma,1,0);
    }

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);
        amrex::Print() << "Diffusion::diffuse_Ssync(): lev: " << level
                       << ", time: " << run_time << '\n';
    }
}

void
Diffusion::computeAlpha (MultiFab&             alpha,
                         std::pair<Real,Real>& scalars,
                         Real                  a,
                         Real                  b,
                         Real*                 rhsscale,
                         const MultiFab*       alpha_in,
                         int                   alpha_in_comp,
                         int                   rho_flag,
                         const MultiFab*       rho,
                         int                   rho_comp)
{
    const BoxArray& ba = rho->boxArray();
    const DistributionMapping& dm = rho->DistributionMap();
    const auto& factory = rho->Factory();

    alpha.define(ba, dm, 1, 0, MFInfo(), factory);

    if (alpha_in != 0)
    {
        BL_ASSERT(alpha_in_comp >= 0 && alpha_in_comp < alpha.nComp());
        MultiFab::Copy(alpha,*alpha_in,alpha_in_comp,0,1,0);
    }
    else
      alpha.setVal(1.0);

    if ( rho_flag > 0 )
    {
        MultiFab::Multiply(alpha,*rho,rho_comp,0,1,0);
    }

    if (rhsscale != 0)
    {
        *rhsscale = scale_abec ? 1.0/alpha.max(0) : 1.0;

        scalars.first = a*(*rhsscale);
        scalars.second = b*(*rhsscale);
    }
    else
    {
        scalars.first = a;
        scalars.second = b;
    }
}

#ifdef AMREX_USE_EB
void
Diffusion::setBeta(MLEBABecLap&           op,
                   const MultiFab* const* beta,
                   int                    betaComp)
{
    Array<MultiFab,AMREX_SPACEDIM> bcoeffs{
    AMREX_D_DECL(MultiFab(*beta[0], amrex::make_alias, betaComp, 1),
                 MultiFab(*beta[1], amrex::make_alias, betaComp, 1),
                 MultiFab(*beta[2], amrex::make_alias, betaComp, 1) ) };

    op.setBCoeffs(0, amrex::GetArrOfConstPtrs(bcoeffs), MLLinOp::Location::FaceCentroid);
}

void
Diffusion::setViscosity(MLEBTensorOp&          tensorop,
                        const MultiFab* const* beta,
                        int                    betaComp,
                        const MultiFab&        beta_cc)
{
    Array<MultiFab,AMREX_SPACEDIM> face_bcoef{
    AMREX_D_DECL(MultiFab(*beta[0], amrex::make_alias, betaComp, 1),
                 MultiFab(*beta[1], amrex::make_alias, betaComp, 1),
                 MultiFab(*beta[2], amrex::make_alias, betaComp, 1) ) };
    tensorop.setShearViscosity(0, amrex::GetArrOfConstPtrs(face_bcoef), MLLinOp::Location::FaceCentroid);

    MultiFab cc_bcoef = MultiFab(beta_cc, amrex::make_alias, betaComp, 1);
    tensorop.setEBShearViscosity(0, cc_bcoef);
}

#else

void
Diffusion::setBeta(MLABecLaplacian&       op,
                   const MultiFab* const* beta,
                   int                    betaComp)
{
    Array<MultiFab,AMREX_SPACEDIM> bcoeffs{
    AMREX_D_DECL(MultiFab(*beta[0], amrex::make_alias, betaComp, 1),
                 MultiFab(*beta[1], amrex::make_alias, betaComp, 1),
                 MultiFab(*beta[2], amrex::make_alias, betaComp, 1) ) };

    op.setBCoeffs(0, amrex::GetArrOfConstPtrs(bcoeffs));
}

void
Diffusion::setViscosity(MLTensorOp&            tensorop,
                        const MultiFab* const* beta,
                        int                    betaComp)
{
    Array<MultiFab,AMREX_SPACEDIM> face_bcoef{
    AMREX_D_DECL(MultiFab(*beta[0], amrex::make_alias, betaComp, 1),
                 MultiFab(*beta[1], amrex::make_alias, betaComp, 1),
                 MultiFab(*beta[2], amrex::make_alias, betaComp, 1) ) };
    tensorop.setShearViscosity(0, amrex::GetArrOfConstPtrs(face_bcoef));
}
#endif

//
// return the fluxes multiplied by the area -- extensive fluxes
//
void
Diffusion::computeExtensiveFluxes(MLMG& a_mg, MultiFab& Soln,
                                  MultiFab* const* flux, const int fluxComp,
                                  const int ncomp,
                                  const Geometry* a_geom, const Real fac )
{
   BL_ASSERT(flux[0]->nGrow()==0);

   AMREX_D_TERM(MultiFab flxx(*flux[0], amrex::make_alias, fluxComp, ncomp);,
                MultiFab flxy(*flux[1], amrex::make_alias, fluxComp, ncomp);,
                MultiFab flxz(*flux[2], amrex::make_alias, fluxComp, ncomp););
   std::array<MultiFab*,AMREX_SPACEDIM> fp{AMREX_D_DECL(&flxx,&flxy,&flxz)};

   a_mg.getFluxes({fp},{&Soln},MLLinOp::Location::FaceCentroid);

    //
    // fixme? which way to define area?
    //
    // amrex::MultiFab area[BL_SPACEDIM];
    // DistributionMapping dmap = Soln.DistributionMap();
    // BoxArray ba = Soln.boxArray();
    // for (int dir = 0; dir < BL_SPACEDIM; ++dir)
    // {
    //     area[dir].clear();
    // 	area[dir].define(ba.surroundingNodes(dir),dmap,1,nghost);
    //     a_geom->GetFaceArea(area[dir],dir);
    // }

   const Real*  dx = a_geom->CellSize();
#if ( AMREX_SPACEDIM == 3 )
   Real areax = dx[1]*dx[2];
   Real areay = dx[0]*dx[2];
   Real areaz = dx[0]*dx[1];
#else
   Real areax = dx[1];
   Real areay = dx[0];
#endif

#ifdef AMREX_USE_EB
   auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Soln.Factory());
   Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
   areafrac  = ebfactory.getAreaFrac();
#endif

   MFItInfo mfi_info;
   if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(Soln, mfi_info); mfi.isValid(); ++mfi)
   {
      D_TERM( const auto& fx = flxx.array(mfi);,
              const auto& fy = flxy.array(mfi);,
              const auto& fz = flxz.array(mfi););

      D_TERM( const Box ubx = mfi.nodaltilebox(0);,
              const Box vbx = mfi.nodaltilebox(1);,
              const Box wbx = mfi.nodaltilebox(2););

#ifdef AMREX_USE_EB
      Box bx = mfi.tilebox();

      const EBFArrayBox&     cc_fab = static_cast<EBFArrayBox const&>(Soln[mfi]);
      const EBCellFlagFab&    flags = cc_fab.getEBCellFlagFab();

      if ( flags.getType(amrex::grow(bx,0)) == FabType::covered )
      {
         //
         // For now, set to very large num so we know if you accidentally use it
         // MLMG will set covered fluxes to zero
         //
         D_TERM(AMREX_PARALLEL_FOR_4D(ubx, ncomp, i, j, k, n, {fx(i,j,k,n) = COVERED_VAL;});,
                AMREX_PARALLEL_FOR_4D(vbx, ncomp, i, j, k, n, {fy(i,j,k,n) = COVERED_VAL;});,
                AMREX_PARALLEL_FOR_4D(wbx, ncomp, i, j, k, n, {fz(i,j,k,n) = COVERED_VAL;}););
      }
      else if ( flags.getType(amrex::grow(bx,0)) != FabType::regular )
      {
	 D_TERM( const auto& afrac_x = areafrac[0]->array(mfi);,
		 const auto& afrac_y = areafrac[1]->array(mfi);,
		 const auto& afrac_z = areafrac[2]->array(mfi););

	 D_TERM(AMREX_PARALLEL_FOR_4D(ubx, ncomp, i, j, k, n, {fx(i,j,k,n) *= fac*areax*afrac_x(i,j,k);});,
		AMREX_PARALLEL_FOR_4D(vbx, ncomp, i, j, k, n, {fy(i,j,k,n) *= fac*areay*afrac_y(i,j,k);});,
                AMREX_PARALLEL_FOR_4D(wbx, ncomp, i, j, k, n, {fz(i,j,k,n) *= fac*areaz*afrac_z(i,j,k);}););
      }
      else
#endif
      {
	D_TERM(AMREX_PARALLEL_FOR_4D(ubx, ncomp, i, j, k, n, {fx(i,j,k,n) *= fac*areax;});,
	       AMREX_PARALLEL_FOR_4D(vbx, ncomp, i, j, k, n, {fy(i,j,k,n) *= fac*areay;});,
	       AMREX_PARALLEL_FOR_4D(wbx, ncomp, i, j, k, n, {fz(i,j,k,n) *= fac*areaz;}););
      }
   }
}

void
Diffusion::getViscTerms (MultiFab&              visc_terms,
                         int                    src_comp,
                         int                    comp,
                         Real                   time,
                         int                    rho_flag,
                         const MultiFab* const* beta,
                         int                    betaComp)
{
    int allthere;
    checkBeta(beta, allthere);
    //
    // Before computing the godunov predictors we may have to
    // precompute the viscous source terms.  To do this we must
    // construct a Laplacian operator, set the coeficients and apply
    // it to the time N data.  First, however, we must precompute the
    // fine N bndry values.  We will do this for each scalar that diffuses.
    //
    // Note: This routine DOES NOT fill grow cells
    //

    //
    // FIXME
    // LinOp classes cannot handle multcomponent MultiFabs yet,
    // construct the components one at a time and copy to visc_terms.
    //

    if (is_diffusive[comp])
    {
        int ng = 1;
        int ng_visc(2);// needed for redistribution
        MultiFab visc_tmp(grids,dmap,1,ng_visc,MFInfo(),navier_stokes->Factory()),
                 s_tmp(grids,dmap,1,ng,MFInfo(),navier_stokes->Factory());
	// not sure this is needed...
        visc_tmp.setVal(0.);
        //
        // Set up operator and apply to compute viscous terms.
        //
        const Real a = 0.0;
        const Real b = -1.0;

	LPInfo info;
	info.setAgglomeration(agglomeration);
	info.setConsolidation(consolidation);
	info.setMaxCoarseningLevel(0);
	//
	// For now, assume velocity always goes to tensor sovler, so it will not get here
	// Otherwise, I *think* we would need to check component and only turn on metric
	// for Xvel
	info.setMetricTerm(false);

#ifdef AMREX_USE_EB
	const auto& ebf = &(dynamic_cast<EBFArrayBoxFactory const&>(navier_stokes->Factory()));
	MLEBABecLap mlabec({navier_stokes->Geom()}, {grids}, {dmap}, info, {ebf});
#else
	MLABecLaplacian mlabec({navier_stokes->Geom()},{grids},{dmap},info);
#endif

	mlabec.setMaxOrder(max_order);

	{
	  // set BCs

	  std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
	  std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
	  setDomainBC(mlmg_lobc, mlmg_hibc, comp);

	  mlabec.setDomainBC(mlmg_lobc, mlmg_hibc);

	  MultiFab crsedata;

	  if (level > 0) {
	    auto& crse_ns = *(coarser->navier_stokes);
	    crsedata.define(crse_ns.boxArray(), crse_ns.DistributionMap(), 1, ng,MFInfo(),crse_ns.Factory());
	    AmrLevel::FillPatch(crse_ns,crsedata,ng,time,State_Type,comp,1);
	    if (rho_flag == 2) {
	      // We want to evaluate (div beta grad) S, not rho*S.
	      const MultiFab& rhotime = crse_ns.get_rho(time);
	      MultiFab::Divide(crsedata,rhotime,0,0,1,ng);
	    }
	    mlabec.setCoarseFineBC(&crsedata, crse_ratio[0]);
	  }

	  AmrLevel::FillPatch(*navier_stokes,s_tmp,ng,time,State_Type,comp,1);
	  if (rho_flag == 2) {
	    const MultiFab& rhotime = navier_stokes->get_rho(time);
	    MultiFab::Divide(s_tmp,rhotime,0,0,1,ng);
	  }
	  mlabec.setLevelBC(0, &s_tmp);
	}

	mlabec.setScalars(a,b);
	// mlabec.setACoeffs() not needed since a = 0.0


	setBeta(mlabec, beta, betaComp);

	// Do we need something like this cribbed from mfix???
	// This sets the coefficient on the wall and defines it as a homogeneous
	// Dirichlet bc for the solve. mu_g is the viscosity at cc in mfix
	// matches what's in bcoeff
	//mlabec.setEBHomogDirichlet ( 0, (*mu_g[lev]) );

	MLMG mgn(mlabec);
	mgn.setVerbose(verbose);
	mgn.apply({&visc_tmp},{&s_tmp});

#ifdef AMREX_USE_EB
        const amrex::MultiFab* weights;
        const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(navier_stokes->Factory());
        weights = &(ebfactory.getVolFrac());

        amrex::single_level_weighted_redistribute(visc_tmp, visc_terms, *weights, comp-src_comp, 1,
                                         navier_stokes->Geom());
#else
	MultiFab::Copy(visc_terms,visc_tmp,0,comp-src_comp,1,0);
#endif

    }
    else
    {
      int ngrow = visc_terms.nGrow();
      visc_terms.setVal(0.0,comp-src_comp,1,ngrow);
    }
}

void
Diffusion::getTensorViscTerms (MultiFab&              visc_terms,
                               Real                   time,
                               const MultiFab* const* beta,
                               const MultiFab* const  betaCC,
                               int                    betaComp)
{
    int allthere;
    checkBeta(beta, allthere);

    const int src_comp = Xvel;
    const int ncomp    = visc_terms.nComp();

    if (ncomp < AMREX_SPACEDIM)
        amrex::Abort("Diffusion::getTensorViscTerms(): visc_terms needs at least BL_SPACEDIM components");
    //
    // Before computing the godunov predicitors we may have to
    // precompute the viscous source terms.  To do this we must
    // construct a Laplacian operator, set the coeficients and apply
    // it to the time N data.  First, however, we must precompute the
    // fine N bndry values.  We will do this for each scalar that diffuses.
    //
    // Note: This routine DOES NOT fill grow cells
    //
    MultiFab&   S    = navier_stokes->get_data(State_Type,time);

    if (is_diffusive[src_comp])
    {
        int ng = 1;
        MultiFab visc_tmp(grids,dmap,AMREX_SPACEDIM,2,MFInfo(),navier_stokes->Factory()),
                 s_tmp(grids,dmap,AMREX_SPACEDIM,ng,MFInfo(),navier_stokes->Factory());
        visc_tmp.setVal(0.);
        MultiFab::Copy(s_tmp,S,Xvel,0,AMREX_SPACEDIM,0);

        //
        // Set up operator and apply to compute viscous terms.
        //
        const Real a = 0.0;
        const Real b = -1.0;

        // MLMG tensor solver
        {
           LPInfo info;
           info.setAgglomeration(agglomeration);
           info.setConsolidation(consolidation);
           info.setMaxCoarseningLevel(0);
           info.setMetricTerm(false);

#ifdef AMREX_USE_EB
           const auto& ebf = &dynamic_cast<EBFArrayBoxFactory const&>(navier_stokes->Factory());
           MLEBTensorOp tensorop({navier_stokes->Geom()}, {grids}, {dmap}, info, {ebf});
#else
           MLTensorOp tensorop({navier_stokes->Geom()}, {grids}, {dmap}, info);
#endif
           tensorop.setMaxOrder(tensor_max_order);

           // create right container
           Array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc[AMREX_SPACEDIM];
           Array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc[AMREX_SPACEDIM];
           // fill it
           for (int i=0; i<AMREX_SPACEDIM; i++)
             setDomainBC(mlmg_lobc[i], mlmg_hibc[i], Xvel+i);
           // pass to op
           tensorop.setDomainBC({AMREX_D_DECL(mlmg_lobc[0],mlmg_lobc[1],mlmg_lobc[2])},
                                {AMREX_D_DECL(mlmg_hibc[0],mlmg_hibc[1],mlmg_hibc[2])});

           // set coarse-fine BCs
           {
              MultiFab crsedata;

              if (level > 0) {
                auto& crse_ns = *(coarser->navier_stokes);
                crsedata.define(crse_ns.boxArray(), crse_ns.DistributionMap(),
                                AMREX_SPACEDIM, ng, MFInfo(),crse_ns.Factory());
                AmrLevel::FillPatch(crse_ns, crsedata, ng, time, State_Type, Xvel,
                                    AMREX_SPACEDIM);
                tensorop.setCoarseFineBC(&crsedata, crse_ratio[0]);
              }
              AmrLevel::FillPatch(*navier_stokes,s_tmp,ng,time,State_Type,Xvel,AMREX_SPACEDIM);
              tensorop.setLevelBC(0, &s_tmp);

              // FIXME: check divergence of vel
              // MLNodeLaplacian mllap({navier_stokes->Geom()}, {grids}, {dmap}, info);
              // mllap.setDomainBC(mlmg_lobc[0], mlmg_hibc[0]);
              // Rhs2.setVal(0.);
              // mllap.compDivergence({&Rhs2}, {&s_tmp});
              // amrex::WriteSingleLevelPlotfile("div_"+std::to_string(count), Rhs2, {AMREX_D_DECL("x","y","z")},navier_stokes->Geom(), 0.0, 0);
              //
           }

           tensorop.setScalars(a, b);

#ifdef AMREX_USE_EB
           setViscosity(tensorop, beta, betaComp, *betaCC);
#else
           setViscosity(tensorop, beta, betaComp);
#endif

           MLMG mlmg(tensorop);
           // FIXME -- consider making new parameters max_iter and bottom_verbose
           //mlmg.setMaxIter(max_iter);
           mlmg.setMaxFmgIter(max_fmg_iter);
           mlmg.setVerbose(10);
           mlmg.setBottomVerbose(10);
           //mlmg.setBottomVerbose(bottom_verbose);

           mlmg.apply({&visc_tmp}, {&s_tmp});
        }

#ifdef AMREX_USE_EB
        const amrex::MultiFab* weights;
        const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(navier_stokes->Factory());
        weights = &(ebfactory.getVolFrac());
        amrex::single_level_weighted_redistribute(visc_tmp, visc_terms, *weights, 0, AMREX_SPACEDIM,
                                                  navier_stokes->Geom() );
#else
        MultiFab::Copy(visc_terms,visc_tmp,0,0,AMREX_SPACEDIM,0);
#endif
    }
    else
    {
        int ngrow = visc_terms.nGrow();
        visc_terms.setVal(0.0,src_comp,AMREX_SPACEDIM,ngrow);
    }
}

void
Diffusion::FillBoundary (BndryRegister& bdry,
                         int            state_ind,
                         int            dest_comp,
                         int            num_comp,
                         Real           time,
                         int            rho_flag)
{
    //
    // Need one grow cell filled for linear solvers.
    // We assume filPatch gets this right, where possible.
    //
    const int     nGrow = 1;

    const MultiFab& rhotime = navier_stokes->get_rho(time);

    MultiFab S(navier_stokes->boxArray(),
               navier_stokes->DistributionMap(),
               num_comp,nGrow,MFInfo(),navier_stokes->Factory());

    AmrLevel::FillPatch(*navier_stokes,S,nGrow,time,State_Type,state_ind,num_comp);

    if (rho_flag == 2) {
        for (int n = 0; n < num_comp; ++n) {
	    MultiFab::Divide(S,rhotime,0,n,1,nGrow);
	}
    }

    //
    // Copy into boundary register.
    //
    bdry.copyFrom(S,nGrow,0,dest_comp,num_comp);

}

void
Diffusion::checkBeta (const MultiFab* const* beta,
                      int&                   allthere,
                      int&                   allnull)
{
    allnull  = 1;
    allthere = beta != 0;

    if (allthere)
    {
        for (int d = 0; d < BL_SPACEDIM; d++)
        {
            allnull = allnull && beta[d] == 0;
            allthere = allthere && beta[d] != 0;
        }
    }

    if (!(allthere || allnull))
        amrex::Abort("Diffusion::checkBeta(): beta must be all 0 or all non-0");
}

void
Diffusion::checkBeta (const MultiFab* const* beta,
                      int&                   allthere) const
{
    allthere = beta != 0;

    if (allthere)
    {
        for (int d = 0; d < BL_SPACEDIM; d++)
            allthere = allthere && beta[d] != 0;
    }

    if (!allthere)
        amrex::Abort("Diffusion::checkBeta(): beta must be all non-0");
}

//
// SAS: The following function is a temporary fix in the migration from
//      using is_conservative and rho_flag over to using advectionType
//      and diffusionType.
//
int
Diffusion::set_rho_flag(const DiffusionForm compDiffusionType)
{
    int rho_flag = 0;

    switch (compDiffusionType)
    {
        case Laplacian_S:
            rho_flag = 0;
            break;

        case RhoInverse_Laplacian_S:
            rho_flag = 1;
            break;

        case Laplacian_SoverRho:
            rho_flag = 2;
            break;

        //NOTE: rho_flag = 3 is used in a different context for
        //      do_mom_diff==1

        default:
            amrex::Print() << "compDiffusionType = " << compDiffusionType << '\n';
            amrex::Abort("An unknown NavierStokesBase::DiffusionForm was used in set_rho_flag");
    }

    return rho_flag;
}

bool
Diffusion::are_any(const Vector<DiffusionForm>& diffusionType,
                   const DiffusionForm         testForm,
                   const int                   sComp,
                   const int                   nComp)
{
    for (int comp = sComp; comp < sComp + nComp; ++comp)
    {
        if (diffusionType[comp] == testForm)
            return true;
    }

    return false;
}

int
Diffusion::how_many(const Vector<DiffusionForm>& diffusionType,
                    const DiffusionForm         testForm,
                    const int                   sComp,
                    const int                   nComp)
{
    int counter = 0;

    for (int comp = sComp; comp < sComp + nComp; ++comp)
    {
        if (diffusionType[comp] == testForm)
            ++counter;
    }

    return counter;
}

void
Diffusion::setDomainBC (std::array<LinOpBCType,AMREX_SPACEDIM>& mlmg_lobc,
                        std::array<LinOpBCType,AMREX_SPACEDIM>& mlmg_hibc,
                        int src_comp)
{
    const BCRec& bc = navier_stokes->get_desc_lst()[State_Type].getBC(src_comp);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (parent->Geom(0).isPeriodic(idim))
        {
            mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Periodic;
        }
        else
        {
            int pbc = bc.lo(idim);
            if (pbc == EXT_DIR)
            {
                mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            }
            else if (pbc == FOEXTRAP      ||
                     pbc == HOEXTRAP      ||
                     pbc == REFLECT_EVEN)
            {
                mlmg_lobc[idim] = LinOpBCType::Neumann;
            }
            else if (pbc == REFLECT_ODD)
            {
                mlmg_lobc[idim] = LinOpBCType::reflect_odd;
            }
            else
            {
                mlmg_lobc[idim] = LinOpBCType::bogus;
            }

            pbc = bc.hi(idim);
            if (pbc == EXT_DIR)
            {
                mlmg_hibc[idim] = LinOpBCType::Dirichlet;
            }
            else if (pbc == FOEXTRAP      ||
                     pbc == HOEXTRAP      ||
                     pbc == REFLECT_EVEN)
            {
                mlmg_hibc[idim] = LinOpBCType::Neumann;
            }
            else if (pbc == REFLECT_ODD)
            {
                mlmg_hibc[idim] = LinOpBCType::reflect_odd;
            }
            else
            {
                mlmg_hibc[idim] = LinOpBCType::bogus;
            }
        }
    }
}

void
Diffusion::setDomainBC (std::array<LinOpBCType,AMREX_SPACEDIM>& mlmg_lobc,
                            std::array<LinOpBCType,AMREX_SPACEDIM>& mlmg_hibc,
                            const BCRec& bc)
{

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
      // DefaultGeometry() is same as parent->Geom(0)
      if (DefaultGeometry().isPeriodic(idim))
        {
            mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Periodic;
        }
        else
        {
            int pbc = bc.lo(idim);
            if (pbc == EXT_DIR)
            {
                mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            }
            else if (pbc == FOEXTRAP      ||
                     pbc == HOEXTRAP      ||
                     pbc == REFLECT_EVEN)
            {
                mlmg_lobc[idim] = LinOpBCType::Neumann;
            }
            else if (pbc == REFLECT_ODD)
            {
                mlmg_lobc[idim] = LinOpBCType::reflect_odd;
            }
            else
            {
                mlmg_lobc[idim] = LinOpBCType::bogus;
            }

            pbc = bc.hi(idim);
            if (pbc == EXT_DIR)
            {
                mlmg_hibc[idim] = LinOpBCType::Dirichlet;
            }
            else if (pbc == FOEXTRAP      ||
                     pbc == HOEXTRAP      ||
                     pbc == REFLECT_EVEN)
            {
                mlmg_hibc[idim] = LinOpBCType::Neumann;
            }
            else if (pbc == REFLECT_ODD)
            {
                mlmg_hibc[idim] = LinOpBCType::reflect_odd;
            }
            else
            {
                mlmg_hibc[idim] = LinOpBCType::bogus;
            }
        }
    }
}
